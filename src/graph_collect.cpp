#include "graph_collect.hpp"
#include "collect_phase_noisy.hpp"
#include "collect_pipeline.hpp"

#include "arg_parse.hpp"
#include "collect_output.hpp"
#include "collect_phase.hpp"
#include "collect_phase_pgbam.hpp"
#include "collect_types.hpp"
#include "collect_var.hpp"
#include "gbz_ffi.h"
#include "graph_bam_adapter.hpp"
#include "graph_query.hpp"
#include "graph_sites.hpp"
#include "noise_filter.hpp"

#include <algorithm>
#include <atomic>
#include <cstdio>
#include <cstdint>
#include <fstream>
#include <getopt.h>
#include <iostream>
#include <memory>
#include <mutex>
#include <sstream>
#include <stdexcept>
#include <string>
#include <thread>
#include <unordered_map>
#include <unordered_set>
#include <vector>
#include <unistd.h>

#include <htslib/faidx.h>
#include <htslib/sam.h>

namespace pgphase_collect {

size_t synthesize_meta_for_appended_candidates(GraphChunkBuildResult& result,
                                               const std::string& contig) {
    PhasingChunk& chunk = result.chunk;
    if (chunk.ref_seq.empty()) return 0;
    if (chunk.candidates.size() <= result.site_meta.size()) return 0;

    const auto ref_at = [&](hts_pos_t pos, int len) -> std::string {
        if (len <= 0) return std::string();
        const hts_pos_t off = pos - chunk.ref_beg;
        if (off < 0) return std::string();
        if (static_cast<size_t>(off) + static_cast<size_t>(len) > chunk.ref_seq.size())
            return std::string();
        return chunk.ref_seq.substr(static_cast<size_t>(off), static_cast<size_t>(len));
    };

    size_t built = 0;
    for (size_t ci = result.site_meta.size(); ci < chunk.candidates.size(); ++ci) {
        const CandidateVariant& cand = chunk.candidates[ci];
        GraphSiteMeta meta;
        meta.chrom = contig;
        const auto anchored = [&](const std::string& allele) -> std::string {
            if (cand.key.type == VariantType::Snp) return allele;
            if (cand.key.type == VariantType::Insertion && cand.key.ref_len == 0)
                return ref_at(cand.key.pos - 1, 1) + allele;
            if (cand.key.type == VariantType::Insertion) return allele;
            return ref_at(cand.key.pos - 1, 1) + allele;
        };
        std::string vcf_ref, vcf_alt;
        hts_pos_t vcf_pos = cand.key.pos;
        if (cand.key.type == VariantType::Snp) {
            vcf_ref = ref_at(cand.key.pos, std::max(1, cand.key.ref_len));
        } else if (cand.key.type == VariantType::Insertion && cand.key.ref_len == 0) {
            vcf_pos = cand.key.pos - 1;
            vcf_ref = ref_at(vcf_pos, 1);
        } else if (cand.key.type == VariantType::Insertion) {
            vcf_ref = ref_at(cand.key.pos, cand.key.ref_len);
        } else {
            vcf_pos = cand.key.pos - 1;
            vcf_ref = ref_at(vcf_pos, 1) + ref_at(cand.key.pos, cand.key.ref_len);
        }
        vcf_alt = anchored(cand.key.alt);

        bool usable = !vcf_ref.empty() && !vcf_alt.empty();
        if (usable) {
            const VariantKey round_trip =
                vcf_to_variant_key(cand.key.tid, vcf_pos, vcf_ref, vcf_alt);
            usable = round_trip.type == cand.key.type && round_trip.pos == cand.key.pos &&
                     round_trip.ref_len == cand.key.ref_len && round_trip.alt == cand.key.alt;
        }
        std::vector<std::string> alts;
        if (cand.msa_insertion_alts.size() >= 2) {
            for (const std::string& allele : cand.msa_insertion_alts) {
                const std::string a = anchored(allele);
                if (a.empty()) { alts.clear(); break; }
                alts.push_back(a);
            }
        }
        if (alts.empty() && usable) alts.push_back(vcf_alt);
        if (usable && !alts.empty()) {
            meta.pos = vcf_pos;
            meta.ref = vcf_ref;
            meta.alts = alts;
            ++built;
        }
        // An unusable site still gets an EMPTY entry: the lookup is by index,
        // so skipping a push would shift every later site's metadata.
        result.site_meta.push_back(std::move(meta));
        result.site_ids.push_back(std::string());
        std::vector<int> orig;
        for (int i = 0; i <= static_cast<int>(alts.size()); ++i) orig.push_back(i);
        if (orig.size() < 2) orig = std::vector<int>{0, 1};
        result.site_allele_orig_idx.push_back(std::move(orig));
    }
    return built;
}

void seed_graph_noisy_regions(GraphChunkBuildResult& result,
                              const std::string& ref_seq,
                              hts_pos_t ref_beg,
                              hts_pos_t ref_end) {
    PhasingChunk& chunk = result.chunk;
    if (ref_seq.empty()) return;

    // The MSA reads the reference off the chunk; a graph chunk never carried
    // one (build_graph_chunk leaves ref_seq empty), which is the first reason
    // stage 2 could not run here.
    chunk.ref_seq = ref_seq;
    chunk.ref_beg = ref_beg;
    chunk.ref_end = ref_end;
    populate_low_complexity_intervals(chunk);

    const std::vector<Interval> lc = find_low_complexity_intervals(ref_seq, ref_beg);
    std::vector<Interval> seeded;
    for (const CandidateVariant& cand : chunk.candidates) {
        if (cand.counts.category != VariantCategory::RepeatHetIndel) continue;
        hts_pos_t beg = 0, end = 0;
        variant_genomic_span(cand.key, beg, end);
        // Same widening as cr_add_var_to_noisy_cr: a locus inside a
        // low-complexity tract takes the whole tract, so the MSA sees the
        // repeat it has to resolve rather than one anchor inside it.
        for (const Interval& iv : lc) {
            if (iv.end < beg || iv.beg > end) continue;
            beg = std::min(beg, iv.beg);
            end = std::max(end, iv.end);
        }
        if (beg <= 0 || end < beg) continue;
        seeded.push_back(Interval{beg, end, 0});
    }
    if (seeded.empty()) return;

    std::sort(seeded.begin(), seeded.end(),
              [](const Interval& a, const Interval& b) { return a.beg < b.beg; });
    std::vector<Interval> merged;
    for (const Interval& iv : seeded) {
        if (!merged.empty() && iv.beg <= merged.back().end + 1)
            merged.back().end = std::max(merged.back().end, iv.end);
        else
            merged.push_back(iv);
    }
    for (const Interval& iv : merged) chunk.noisy_regions.push_back(iv);
    std::sort(chunk.noisy_regions.begin(), chunk.noisy_regions.end(),
              [](const Interval& a, const Interval& b) { return a.beg < b.beg; });
}


namespace {

static bam_hdr_t* build_synthetic_header(faidx_t* fai) {
    bam_hdr_t* hdr = sam_hdr_init();
    if (!hdr) throw std::runtime_error("failed to allocate synthetic BAM header");
    const int n_seq = faidx_nseq(fai);
    for (int i = 0; i < n_seq; ++i) {
        const char* name = faidx_iseq(fai, i);
        const hts_pos_t len = faidx_seq_len(fai, name);
        const std::string len_str = std::to_string(len);
        if (sam_hdr_add_line(hdr, "SQ", "SN", name, "LN", len_str.c_str(), nullptr) < 0) {
            sam_hdr_destroy(hdr);
            throw std::runtime_error(std::string("failed to add SQ line for ") + name);
        }
    }
    return hdr;
}

// Converts phased multi-allelic graph candidates to biallelic CandidateVariants compatible
// with the existing TSV/VCF writers. One output row per passing alt allele per site.
static CandidateTable graph_chunks_to_candidate_table(
    const std::vector<GraphChunkBuildResult>& graph_chunks,
    const std::unordered_map<std::string, int>& contig_to_tid,
    const Options& opts)
{
    CandidateTable result;

    for (const GraphChunkBuildResult& graph_chunk : graph_chunks) {
        const PhasingChunk& chunk = graph_chunk.chunk;
        for (size_t ci = 0; ci < chunk.candidates.size(); ++ci) {
            const CandidateVariant& mcand = chunk.candidates[ci];

            if (ci >= graph_chunk.site_meta.size()) continue;
            const GraphSiteMeta& meta = graph_chunk.site_meta[ci];

            // A site merged in by the in-chunk recovery (no catalog id) is
            // written only when it carries phase information. The alignment's
            // in-gap discovery also calls homozygous variants, and emitting
            // those would add ~12,850 hom records to chr20 -- a 23% larger VCF
            // whose extra content is alignment-discovered calls appearing ONLY
            // inside recovery windows, which is a biased subset of the genome.
            // The arm's output is the catalog's sites plus what recovery
            // phased, not a variant call set.
            if (ci < graph_chunk.site_ids.size() && graph_chunk.site_ids[ci].empty()) {
                const auto& h = mcand.hap_to_cons_alle;
                const int n_alleles = static_cast<int>(mcand.counts.alle_covs.size());
                // Both haplotype consensus alleles must be VALID INDICES for this
                // site's allele set as merged, and they must differ. A merged site
                // is collapsed to biallelic (alle_covs = {ref_cov, alt_cov}), so a
                // consensus index of 2 survives from a wider alignment
                // representation and has no allele here; such sites were emitted
                // 1|1 -- 186 of them on chr20, carrying hap_to_cons_alle (2,1) or
                // (1,2).
                if (h.size() < 3 || h[1] < 0 || h[2] < 0 || h[1] >= n_alleles ||
                    h[2] >= n_alleles || h[1] == h[2])
                    continue;
            }

            auto tid_it = contig_to_tid.find(meta.chrom);
            if (tid_it == contig_to_tid.end()) continue;
            const int fai_tid = tid_it->second;

            const std::vector<int>* orig_idx =
                (ci < graph_chunk.site_allele_orig_idx.size())
                    ? &graph_chunk.site_allele_orig_idx[ci]
                    : nullptr;

            const std::vector<int>& alle_covs = mcand.counts.alle_covs;
            const int n_new_alleles = static_cast<int>(alle_covs.size());

            for (int new_a = 1; new_a < n_new_alleles; ++new_a) {
                // Map surviving allele index back to original walk index in GraphSite.
                int orig_walk_idx = new_a;
                if (orig_idx != nullptr && new_a < static_cast<int>(orig_idx->size())) {
                    orig_walk_idx = (*orig_idx)[new_a];
                }

                // orig_walk_idx: 0 = ref walk, 1 = first alt walk → meta.alts[0], etc.
                const int alt_idx = orig_walk_idx - 1;
                if (alt_idx < 0 || alt_idx >= static_cast<int>(meta.alts.size())) continue;
                const std::string& raw_alt = meta.alts[static_cast<size_t>(alt_idx)];
                if (raw_alt.empty() || raw_alt == "*") continue;

                if (meta.ref.empty()) continue;

                // Normalize to minimal VCF form before deriving the key. Catalog
                // alleles carry flanking repeat context on both sides; an equal-
                // length SNP padded with repeat bases (e.g. the AGGG array at
                // chr20:49031440) otherwise falls through to the MNP branch with a
                // misleading multi-bp ref/alt and the wrong POS, and the same SNP
                // is emitted once per overlapping snarl. Trimming the shared suffix
                // then prefix yields the BAM pipeline's canonical representation so
                // keys match and duplicates collapse. (apply_graph_noise_filter
                // already trims for the noise check; this trims the key itself.)
                std::string ref_seq = meta.ref;
                std::string alt_seq = raw_alt;
                hts_pos_t var_pos = meta.pos;
                trim_to_minimal_vcf(var_pos, ref_seq, alt_seq);

                CandidateVariant cand;
                cand.key.tid = fai_tid;

                // Derive VariantKey from VCF-anchored alleles matching BAM-path convention:
                //   SNP : key.pos = site.pos, key.alt = alt base, ref_len = 1
                //   INS : key.pos = site.pos+1 (after anchor), key.alt = inserted bases
                //   DEL : key.pos = site.pos+1 (first deleted base), key.alt = deleted bases
                if (ref_seq.size() == 1 && alt_seq.size() == 1) {
                    // SNP
                    cand.key.type = VariantType::Snp;
                    cand.key.pos = var_pos;
                    cand.key.alt = alt_seq;
                    cand.key.ref_len = 1;
                } else if (alt_seq.size() > ref_seq.size() && ref_seq[0] == alt_seq[0]) {
                    // Left-anchored insertion: strip shared prefix; pos after anchor.
                    cand.key.type = VariantType::Insertion;
                    cand.key.pos = var_pos + 1;
                    cand.key.alt = alt_seq.substr(ref_seq.size());
                    cand.key.ref_len = 0;
                } else if (ref_seq.size() > alt_seq.size() && ref_seq[0] == alt_seq[0]) {
                    // Left-anchored deletion: empty alt matches BAM-path convention.
                    cand.key.type = VariantType::Deletion;
                    cand.key.pos = var_pos + static_cast<hts_pos_t>(alt_seq.size());
                    cand.key.alt = "";
                    cand.key.ref_len = static_cast<int>(ref_seq.size() - alt_seq.size());
                } else {
                    // Complex / MNP: no shared anchor base — classify by net length change.
                    cand.key.pos = var_pos;
                    cand.key.alt = alt_seq;
                    cand.key.ref_len = static_cast<int>(ref_seq.size());
                    if (alt_seq.size() > ref_seq.size()) {
                        cand.key.type = VariantType::Insertion;
                    } else if (ref_seq.size() > alt_seq.size()) {
                        cand.key.type = VariantType::Deletion;
                    } else {
                        cand.key.type = VariantType::Snp;  // MNP (equal length)
                    }
                }

                // Heterozygous between two ALTERNATE alleles: the reference is not
                // one of this site's haplotypes, so every quantity defined against
                // it is degenerate -- ref_cov is 0, and the allele fraction
                // alt/(ref+alt) is 1.0 for BOTH alleles, which reads as a
                // homozygous alt and then as LOW_AF. The site's own coverage is the
                // denominator that means something: at 4,785,719 that turns two
                // fractions of 1.00 into 0.34 and 0.66.
                const auto& hcons = mcand.hap_to_cons_alle;
                const bool het_by_consensus = hcons.size() > 2 && hcons[1] >= 0 &&
                                              hcons[2] >= 0 && hcons[1] != hcons[2];
                const int ref_cov = alle_covs.empty() ? 0 : alle_covs[0];
                const int alt_cov = alle_covs[static_cast<size_t>(new_a)];
                int total_cov = ref_cov + alt_cov;
                if (het_by_consensus) {
                    total_cov = 0;
                    for (const int c : alle_covs) total_cov += c;
                }
                cand.counts.ref_cov = ref_cov;
                cand.counts.alt_cov = alt_cov;
                cand.counts.total_cov = total_cov;
                // Preserve the strand split computed during biallelic decomposition
                // (graph_bam_adapter.cpp). Candidates are biallelic here, so the
                // chunk candidate's forward/reverse fields already correspond to
                // this ref/alt pair; copying them keeps REVERSE counts non-zero.
                cand.counts.forward_ref = mcand.counts.forward_ref;
                cand.counts.reverse_ref = mcand.counts.reverse_ref;
                cand.counts.forward_alt = mcand.counts.forward_alt;
                cand.counts.reverse_alt = mcand.counts.reverse_alt;
                cand.counts.allele_fraction =
                    total_cov > 0 ? static_cast<double>(alt_cov) / total_cov : 0.0;
                cand.counts.n_uniq_alles = 2;
                cand.counts.alle_covs = {ref_cov, alt_cov};

                // A site whose two haplotype consensus alleles DIFFER is
                // heterozygous even when no read carries the reference: a 1|2
                // site has ref_cov == 0 by construction. Deciding hom from depth
                // alone discards exactly the sites that bridge a gap -- measured
                // in chr20:4,766,928-4,792,960, where 4,785,719 ('ATTTT' at 22
                // reads against a pure 25 bp deletion at 43) and 4,791,668 (16 T
                // at 34 against 17 T at 24) are the two heterozygotes the gap
                // needs and both were called CleanHom here.
                const bool is_hom_alt = (ref_cov == 0 && alt_cov >= opts.min_alt_depth) &&
                                        !het_by_consensus;
                if (alt_cov < opts.min_alt_depth || total_cov < opts.min_depth) {
                    cand.counts.category = VariantCategory::LowCoverage;
                    cand.counts.candvarcate_initial = VariantCategory::LowCoverage;
                    cand.lcd_var_i_to_cate = kLongcalldLowCovVar;
                } else if (is_hom_alt) {
                    cand.counts.category = VariantCategory::CleanHom;
                    cand.counts.candvarcate_initial = VariantCategory::CleanHom;
                    cand.lcd_var_i_to_cate = kCandCleanHom;
                } else if (cand.counts.allele_fraction < opts.min_af ||
                           cand.counts.allele_fraction > opts.max_af) {
                    cand.counts.category = VariantCategory::LowAlleleFraction;
                    cand.counts.candvarcate_initial = VariantCategory::LowAlleleFraction;
                    cand.lcd_var_i_to_cate = kLongcalldLowAfVar;
                } else if (cand.key.type == VariantType::Snp) {
                    cand.counts.category = VariantCategory::CleanHetSnp;
                    cand.counts.candvarcate_initial = VariantCategory::CleanHetSnp;
                    cand.lcd_var_i_to_cate = kCandCleanHetSnp;
                } else if (mcand.counts.category == VariantCategory::RepeatHetIndel) {
                    // Preserve the noise-filter demotion (apply_graph_noise_filter).
                    // Re-classifying from scratch would re-promote homopolymer/STR
                    // indels to CleanHetIndel, letting them pass the germline output
                    // gate as phased het indels even though they were excluded from
                    // k-means (hence never properly phased).
                    cand.counts.category = VariantCategory::RepeatHetIndel;
                    cand.counts.candvarcate_initial = VariantCategory::RepeatHetIndel;
                    cand.lcd_var_i_to_cate = kLongcalldRepHetVar;
                } else {
                    cand.counts.category = VariantCategory::CleanHetIndel;
                    cand.counts.candvarcate_initial = VariantCategory::CleanHetIndel;
                    cand.lcd_var_i_to_cate = kCandCleanHetIndel;
                }

                // Translate multi-allelic hap_to_cons_alle to biallelic for this alt.
                // Homozygous alt: both haplotypes carry the alt allele, no phase set.
                const int hap1 = mcand.hap_to_cons_alle[1];
                const int hap2 = mcand.hap_to_cons_alle[2];
                cand.hap_to_cons_alle[0] = -1;
                if (is_hom_alt) {
                    cand.hap_to_cons_alle[1] = 1;
                    cand.hap_to_cons_alle[2] = 1;
                    cand.hap_alt = 1;
                    cand.hap_ref = 1;
                    cand.phase_set = -1;
                } else {
                    cand.hap_to_cons_alle[1] = (hap1 == new_a) ? 1 : 0;
                    cand.hap_to_cons_alle[2] = (hap2 == new_a) ? 1 : 0;
                    cand.hap_alt = cand.hap_to_cons_alle[1];
                    cand.hap_ref = cand.hap_to_cons_alle[2];
                    cand.phase_set = mcand.phase_set;
                }
                cand.alt_ref_base = 4;  // use FASTA anchor (BAM-path default)
                cand.lcd_make_variants_region_pass = true;

                // A site merged in by the in-chunk recovery is written only when
                // this writer's own classification calls it a het. Two reasons.
                // The reclassification above sets CleanHom whenever ref_cov == 0
                // (is_hom_alt, line ~196), and an alignment candidate merged from
                // a gap often carries ref_cov = 0 -- measured, ref_cov 0 against
                // alt_cov 57 -- so the depth synthesis manufactures hom calls:
                // 408 extra 1|1 records on chr20 against 125 in the default. And
                // a hom site carries no phase information in any case, so a VCF
                // whose contract is the catalog's sites plus what recovery phased
                // has no reason to gain hom calls that appear only inside
                // recovery windows.
                if (ci < graph_chunk.site_ids.size() && graph_chunk.site_ids[ci].empty()) {
                    const VariantCategory c = cand.counts.category;
                    if (c != VariantCategory::CleanHetSnp &&
                        c != VariantCategory::CleanHetIndel &&
                        c != VariantCategory::NoisyCandHet)
                        continue;
                }

                result.push_back(std::move(cand));
            }
        }
    }

    std::stable_sort(result.begin(), result.end(),
                     [](const CandidateVariant& a, const CandidateVariant& b) {
                         return exact_comp_cand_var(&a, &b) < 0;
                     });

    // Collapse duplicate variants. Overlapping/nested snarls in the catalog can
    // emit the same physical variant from several sites; after minimal-VCF
    // normalization these share an identity and sort adjacently. Keep the
    // highest-coverage copy so the output VCF/TSV carries each variant once.
    // (Phasing already ran upstream per snarl site; this only de-duplicates the
    // emitted table.) If two copies disagree on the phased haplotype, the
    // higher-coverage call wins and the conflict is reported under --verbose.
    if (!result.empty()) {
        size_t write = 0;
        int collapsed = 0;
        int hap_conflicts = 0;
        for (size_t read = 1; read < result.size(); ++read) {
            if (exact_comp_cand_var(&result[write], &result[read]) == 0) {
                ++collapsed;
                if (result[write].phase_set == result[read].phase_set &&
                    result[write].hap_alt != result[read].hap_alt) {
                    ++hap_conflicts;
                }
                if (result[read].counts.total_cov > result[write].counts.total_cov) {
                    result[write] = std::move(result[read]);
                }
            } else {
                ++write;
                if (write != read) result[write] = std::move(result[read]);
            }
        }
        result.resize(write + 1);
        if (opts.verbose && collapsed > 0) {
            std::cerr << "graph: collapsed " << collapsed
                      << " duplicate variant record(s) from overlapping snarls ("
                      << hap_conflicts << " with conflicting haplotype calls)\n";
        }
    }

    drop_conflicting_haplotype_alleles(result);

    return result;
}

// Processes one batch of graph chunks in parallel (one thread pool per reg_chunk_i batch,
// mirroring collect_chunk_batch_parallel in collect_pipeline.cpp).
// Each worker queries overlapping reads via the gbz-base FFI (one SQLite connection
// per thread), then build_graph_chunk + assign_hap k-means.
// Peak memory = threads × (reads_per_chunk + sites_per_chunk).
// After all workers join: populate_graph_chunk_overlaps + stitch_chunk_haps.
static std::vector<GraphChunkBuildResult> process_graph_chunk_batch(
    const std::string& sites_vcf,
    const std::vector<RegionChunk>& chunks,
    size_t batch_begin,
    size_t batch_end,
    const bam_hdr_t* header,
    const GraphQueryConfig& qconfig,
    const std::string& ref_sample,
    const std::unordered_map<std::string, std::string>& fai_full_to_suffix,
    const std::unordered_map<std::string, std::string>& chrom_remap,
    const Options& opts,
    const PgbamSidecarData* pgbam_sidecar)
{
    const size_t batch_size = batch_end - batch_begin;
    std::vector<GraphChunkBuildResult> graph_chunks(batch_size);

    const std::string batch_contig = header->target_name[chunks[batch_begin].tid];

    // Pre-compute the query contig (pangenome suffix) once for the batch.
    const std::string batch_query_contig = [&]() -> std::string {
        auto it = fai_full_to_suffix.find(batch_contig);
        return (it != fai_full_to_suffix.end()) ? it->second : batch_contig;
    }();

    // Query the genome-coordinate range of the reference path so we can
    // skip or clamp chunks that fall outside the GBZ subgraph.  This is
    // done once on the main thread before spawning workers.
    PathRange path_range;
    {
        char* err = nullptr;
        void* gbz_tmp = pgphase_gbz_open(qconfig.gbz_db.c_str(), &err);
        if (gbz_tmp) {
            path_range = query_gbz_path_range(gbz_tmp, ref_sample, batch_query_contig);
            pgphase_gbz_close(gbz_tmp);
            if (opts.verbose >= 2 && path_range.valid)
                std::cerr << "GBZ path range for " << ref_sample << "#" << batch_query_contig
                          << ": [" << path_range.start << ", " << path_range.end << ")\n";
        } else {
            if (err) pgphase_gbz_free_string(err);
        }
    }

    const size_t worker_count = std::min<size_t>(static_cast<size_t>(opts.threads), batch_size);
    std::atomic<size_t> next_offset{0};
    std::exception_ptr first_error;
    std::mutex error_mutex;
    std::vector<std::thread> workers;
    workers.reserve(worker_count);

    for (size_t w = 0; w < worker_count; ++w) {
        workers.emplace_back([&]() {
            char* err = nullptr;
            void* gbz_h = pgphase_gbz_open(qconfig.gbz_db.c_str(), &err);
            if (!gbz_h) {
                std::string msg = err ? std::string(err) : "unknown error";
                if (err) pgphase_gbz_free_string(err);
                std::lock_guard<std::mutex> lock(error_mutex);
                if (!first_error)
                    first_error = std::make_exception_ptr(
                        std::runtime_error("failed to open GBZ: " + msg));
                return;
            }
            void* gaf_h = pgphase_gaf_open(qconfig.gaf_db.c_str(), &err);
            if (!gaf_h) {
                std::string msg = err ? std::string(err) : "unknown error";
                if (err) pgphase_gbz_free_string(err);
                pgphase_gbz_close(gbz_h);
                std::lock_guard<std::mutex> lock(error_mutex);
                if (!first_error)
                    first_error = std::make_exception_ptr(
                        std::runtime_error("failed to open GAF-base: " + msg));
                return;
            }
            if (pgphase_gbz_gaf_validate(gbz_h, gaf_h, &err) != 0) {
                std::string msg = err ? std::string(err) : "unknown error";
                if (err) pgphase_gbz_free_string(err);
                pgphase_gaf_close(gaf_h);
                pgphase_gbz_close(gbz_h);
                std::lock_guard<std::mutex> lock(error_mutex);
                if (!first_error)
                    first_error = std::make_exception_ptr(
                        std::runtime_error("GBZ/GAF-base incompatible: " + msg));
                return;
            }
            struct HandleCleanup {
                void* gbz; void* gaf;
                ~HandleCleanup() {
                    if (gaf) pgphase_gaf_close(gaf);
                    if (gbz) pgphase_gbz_close(gbz);
                }
            } cleanup{gbz_h, gaf_h};

            try {
                // One sites VCF handle per thread.
                SitesVcfHandle sites_handle(sites_vcf);

                // Per-thread reference index for noise detection.
                std::unique_ptr<faidx_t, FaiDeleter> thread_fai(
                    load_reference_index(opts.ref_fasta));

                // Per-thread recovery context. htslib handles are not shareable,
                // which is why the parent pipeline opens a set per worker, and it
                // is the only reason recovery could not already live here. Built
                // once per thread rather than per chunk.
                std::unique_ptr<WorkerContext> thread_recovery_ctx;
                if (!opts.bam_files.empty())
                    thread_recovery_ctx = std::make_unique<WorkerContext>(opts);

                while (true) {
                    const size_t offset = next_offset.fetch_add(1);
                    if (offset >= batch_size) break;
                    const RegionChunk& region = chunks[batch_begin + offset];

                    // Load sites for this chunk's region via tabix.
                    GraphSiteCatalog chunk_catalog = load_sites_for_region(
                        sites_handle, batch_contig, region.beg, region.end);
                    // Normalize contig names to match FAI convention.
                    for (GraphSite& s : chunk_catalog.sites) {
                        auto it = chrom_remap.find(s.chrom);
                        if (it != chrom_remap.end()) s.chrom = it->second;
                        if (!s.ref_contig.empty()) {
                            auto it2 = chrom_remap.find(s.ref_contig);
                            if (it2 != chrom_remap.end()) s.ref_contig = it2->second;
                        }
                    }

                    GraphSiteCatalogView chunk_view = chunk_catalog.view_all();

                    std::vector<GraphReadAllele> chunk_rows;
                    if (!chunk_view.empty()) {
                        hts_pos_t q_beg = region.beg - 1;
                        hts_pos_t q_end = region.end;
                        if (path_range.valid) {
                            q_beg = std::max(q_beg, path_range.start);
                            q_end = std::min(q_end, path_range.end);
                        }
                        if (q_beg < q_end) {
                            chunk_rows = query_gbz_interval_gaf_ffi(
                                gbz_h, gaf_h, ref_sample, batch_query_contig,
                                q_beg, q_end,
                                chunk_view, qconfig.min_mapq);
                        }
                    }

                    graph_chunks[offset] = build_graph_chunk(
                        chunk_view,
                        chunk_rows,
                        batch_contig,
                        region.beg - 1,
                        region.end,
                        region.chunk_id,
                        opts);

                    // Noise filter: fetch reference slice and reclassify
                    // indels in homopolymer/repeat/low-complexity contexts.
                    {
                        hts_pos_t ref_len = 0;
                        char* ref_raw = faidx_fetch_seq64(
                            thread_fai.get(), batch_contig.c_str(),
                            region.beg - 1, region.end - 1, &ref_len);
                        if (ref_raw && ref_len > 0) {
                            std::string ref_slice(ref_raw, ref_raw + ref_len);
                            std::free(ref_raw);
                            apply_graph_noise_filter(
                                graph_chunks[offset], ref_slice,
                                region.beg, region.beg + ref_len - 1,
                                opts.noisy_reg_max_xgaps);
                            if (opts.graph_noisy_msa)
                                seed_graph_noisy_regions(graph_chunks[offset], ref_slice,
                                                         region.beg, region.beg + ref_len - 1);
                        } else {
                            std::free(ref_raw);
                        }
                    }

                    // A repeat-context indel may earn its way back in before the solve.
                    promote_link_supported_repeat_indels(graph_chunks[offset], opts);

                    assign_hap_based_on_germline_het_vars_kmeans(
                        graph_chunks[offset].chunk, opts, kCandGermlineClean);

                    // Stage 2, as the alignment arm runs it: the noisy-region MSA
                    // reconstructs the demoted loci, then the second round solves over
                    // clean plus the rebuilt noisy sites.
                    // The noisy-region MSA cannot run on the graph chunk itself:
                    // collect_noisy_reg_reads skips every read with no digars
                    // (collect_phase_noisy.cpp:1042) and graph-only reads have none, so
                    // a seeded region has zero reads and the MSA has nothing to align.
                    // The seeded regions go to the recovery instead, whose sub-solve
                    // re-reads the BAM through process_chunk and does build digars.

                    // Recovery, in the chunk, while this chunk is still the unit
                    // of work. Each chunk's windows are its own, so this needs no
                    // batching and no cross-chunk coordination -- the worker pool
                    // already provides the parallelism the post-hoc path had to
                    // rebuild for itself.
                    //
                    // Running it HERE, before populate_graph_chunk_overlaps and
                    // stitch_chunk_haps, is what removes the reconciliation
                    // layer: phase sets are not final yet, so merged sites are
                    // part of a better first solve rather than an answer grafted
                    // on with a vote and a relabel. The same merge run after the
                    // stitch fragmented the output into 4,624 blocks against 47.
                    if (thread_recovery_ctx != nullptr) {
                        const size_t merged = retry_unphased_windows_in_place(
                            graph_chunks[offset], opts, *thread_recovery_ctx,
                            batch_contig.c_str());
                        if (merged > 0) {
                            // The re-solve needs the recovery windows: allele_depths_call_het
                            // is scoped by them, and without it each haplotype takes its own
                            // majority allele and a genuine het collapses to one side.
                            Options solve_opts = opts;
                            // The depth-based het escape is keyed on the
                            // candidate's own provenance now
                            // (CandidateVariant::bam_injected), not on a window
                            // list, so nothing needs to be carried here. Window
                            // keying admitted every site in the window and cost
                            // 1.153% -> 2.067% read hamming on the default path.
                            // Two rounds, as the alignment pipeline solves: the
                            // clean sites set the gauge, then the merged in-gap
                            // sites -- NOISY_CAND_HET, which the clean mask does
                            // not admit -- join it. Anchored, so the noisy sites
                            // can extend a block without re-deciding the parity
                            // the catalog's clean sites already established.
                            // The merged sites belong to the FIRST solve, not to a
                            // correction applied after one. Both rounds run again
                            // over the union and neither is anchored: the chunk is
                            // solved once, with every site it will ever have.
                            //
                            // Anchoring here was measured and removed. On
                            // chr20:42,500,000-43,000,000, where recovery merges 6
                            // sites across 3 windows, pinning the pre-merge
                            // consensus left 840 of 2,008 reads misplaced inside a
                            // single block -- a switch, not an inversion: the pinned
                            // flank keeps one parity while the merged sites decide
                            // the other. Re-running both rounds unanchored places
                            // every read correctly (0 misplaced). Restricting the
                            // pin to the recovered intervals does NOT help (same 840),
                            // because the parity that has to change is the chunk's,
                            // not the gap's.
                            // Block import: the parent keeps the read labels it
                            // already has and only incorporates the imported
                            // sites. Re-solving from scratch is what every
                            // earlier admission mechanism did, and all of them
                            // degraded the result.
                            assign_hap_based_on_germline_het_vars_kmeans(
                                graph_chunks[offset].chunk, solve_opts, kCandGermlineClean,
                                opts.stitch_recovered);

                            // Stage 2, as the alignment arm runs it: the noisy-region MSA
                            // reconstructs the demoted loci, then the second round solves over
                            // clean plus the rebuilt noisy sites.
                            // The noisy-region MSA cannot run on the graph chunk itself:
                            // collect_noisy_reg_reads skips every read with no digars
                            // (collect_phase_noisy.cpp:1042) and graph-only reads have none, so
                            // a seeded region has zero reads and the MSA has nothing to align.
                            // The seeded regions go to the recovery instead, whose sub-solve
                            // re-reads the BAM through process_chunk and does build digars.
                            assign_hap_based_on_germline_het_vars_kmeans(
                                graph_chunks[offset].chunk, solve_opts, kCandGermlineVarCate,
                                false);
                        }
                    }
                }
            } catch (...) {
                std::lock_guard<std::mutex> lock(error_mutex);
                if (!first_error) first_error = std::current_exception();
            }
        });
    }
    for (std::thread& w : workers) w.join();
    if (first_error) std::rethrow_exception(first_error);

    populate_graph_chunk_overlaps(graph_chunks);

    std::vector<PhasingChunk> phasing_chunks;
    phasing_chunks.reserve(batch_size);
    for (GraphChunkBuildResult& gc : graph_chunks)
        phasing_chunks.push_back(std::move(gc.chunk));
    stitch_chunk_haps(phasing_chunks, &opts, pgbam_sidecar);
    for (size_t i = 0; i < batch_size; ++i)
        graph_chunks[i].chunk = std::move(phasing_chunks[i]);

    return graph_chunks;
}

// Per-chunk tabix queries on an indexed GAF file.  Each worker seeks directly
// to the overlapping region — only the relevant reads are decompressed and
// parsed, making this efficient even for very large (100+ GB) GAF files.
static std::vector<GraphChunkBuildResult> process_graph_chunk_batch_indexed_gaf(
    const std::string& sites_vcf,
    const std::vector<RegionChunk>& chunks,
    size_t batch_begin,
    size_t batch_end,
    const bam_hdr_t* header,
    const std::string& gaf_file,
    int min_mapq,
    const std::unordered_map<std::string, std::string>& fai_full_to_suffix,
    const std::unordered_map<std::string, std::string>& chrom_remap,
    const Options& opts,
    const PgbamSidecarData* pgbam_sidecar)
{
    const size_t batch_size = batch_end - batch_begin;
    std::vector<GraphChunkBuildResult> graph_chunks(batch_size);

    const std::string batch_contig_gaf = header->target_name[chunks[batch_begin].tid];
    const std::string batch_query_contig_gaf = [&]() -> std::string {
        auto it = fai_full_to_suffix.find(batch_contig_gaf);
        return (it != fai_full_to_suffix.end()) ? it->second : batch_contig_gaf;
    }();

    const size_t worker_count = std::min<size_t>(static_cast<size_t>(opts.threads), batch_size);
    std::atomic<size_t> next_offset{0};
    std::exception_ptr first_error;
    std::mutex error_mutex;
    std::vector<std::thread> workers;
    workers.reserve(worker_count);

    for (size_t w = 0; w < worker_count; ++w) {
        workers.emplace_back([&]() {
            try {
                IndexedGafHandle gaf_handle(gaf_file);
                SitesVcfHandle sites_handle(sites_vcf);

                // Per-thread reference index for noise detection.
                std::unique_ptr<faidx_t, FaiDeleter> thread_fai(
                    load_reference_index(opts.ref_fasta));

                // Per-thread recovery context. htslib handles are not shareable,
                // which is why the parent pipeline opens a set per worker, and it
                // is the only reason recovery could not already live here. Built
                // once per thread rather than per chunk.
                std::unique_ptr<WorkerContext> thread_recovery_ctx;
                if (!opts.bam_files.empty())
                    thread_recovery_ctx = std::make_unique<WorkerContext>(opts);

                while (true) {
                    const size_t offset = next_offset.fetch_add(1);
                    if (offset >= batch_size) break;
                    const RegionChunk& region = chunks[batch_begin + offset];

                    GraphSiteCatalog chunk_catalog = load_sites_for_region(
                        sites_handle, batch_contig_gaf, region.beg, region.end);
                    for (GraphSite& s : chunk_catalog.sites) {
                        auto it = chrom_remap.find(s.chrom);
                        if (it != chrom_remap.end()) s.chrom = it->second;
                        if (!s.ref_contig.empty()) {
                            auto it2 = chrom_remap.find(s.ref_contig);
                            if (it2 != chrom_remap.end()) s.ref_contig = it2->second;
                        }
                    }

                    GraphSiteCatalogView chunk_view = chunk_catalog.view_all();

                    std::vector<GraphReadAllele> chunk_rows;
                    if (!chunk_view.empty()) {
                        const hts_pos_t pad = static_cast<hts_pos_t>(opts.gaf_pad);
                        chunk_rows = scan_indexed_gaf_chunk(
                            gaf_handle, batch_query_contig_gaf,
                            std::max<hts_pos_t>(0, region.beg - 1 - pad), region.end + pad,
                            chunk_view, min_mapq);
                    }

                    graph_chunks[offset] = build_graph_chunk(
                        chunk_view,
                        chunk_rows,
                        batch_contig_gaf,
                        region.beg - 1,
                        region.end,
                        region.chunk_id,
                        opts);

                    // Noise filter: fetch reference slice and reclassify
                    // indels in homopolymer/repeat/low-complexity contexts.
                    {
                        hts_pos_t ref_len = 0;
                        char* ref_raw = faidx_fetch_seq64(
                            thread_fai.get(), batch_contig_gaf.c_str(),
                            region.beg - 1, region.end - 1, &ref_len);
                        if (ref_raw && ref_len > 0) {
                            std::string ref_slice(ref_raw, ref_raw + ref_len);
                            std::free(ref_raw);
                            apply_graph_noise_filter(
                                graph_chunks[offset], ref_slice,
                                region.beg, region.beg + ref_len - 1,
                                opts.noisy_reg_max_xgaps);
                            if (opts.graph_noisy_msa)
                                seed_graph_noisy_regions(graph_chunks[offset], ref_slice,
                                                         region.beg, region.beg + ref_len - 1);
                        } else {
                            std::free(ref_raw);
                        }
                    }

                    // A repeat-context indel may earn its way back in before the solve.
                    promote_link_supported_repeat_indels(graph_chunks[offset], opts);

                    assign_hap_based_on_germline_het_vars_kmeans(
                        graph_chunks[offset].chunk, opts, kCandGermlineClean);

                    // Stage 2, as the alignment arm runs it: the noisy-region MSA
                    // reconstructs the demoted loci, then the second round solves over
                    // clean plus the rebuilt noisy sites.
                    // The noisy-region MSA cannot run on the graph chunk itself:
                    // collect_noisy_reg_reads skips every read with no digars
                    // (collect_phase_noisy.cpp:1042) and graph-only reads have none, so
                    // a seeded region has zero reads and the MSA has nothing to align.
                    // The seeded regions go to the recovery instead, whose sub-solve
                    // re-reads the BAM through process_chunk and does build digars.

                    // Recovery, in the chunk, while this chunk is still the unit
                    // of work. Each chunk's windows are its own, so this needs no
                    // batching and no cross-chunk coordination -- the worker pool
                    // already provides the parallelism the post-hoc path had to
                    // rebuild for itself.
                    //
                    // Running it HERE, before populate_graph_chunk_overlaps and
                    // stitch_chunk_haps, is what removes the reconciliation
                    // layer: phase sets are not final yet, so merged sites are
                    // part of a better first solve rather than an answer grafted
                    // on with a vote and a relabel. The same merge run after the
                    // stitch fragmented the output into 4,624 blocks against 47.
                    if (thread_recovery_ctx != nullptr) {
                        const size_t merged = retry_unphased_windows_in_place(
                            graph_chunks[offset], opts, *thread_recovery_ctx,
                            batch_contig_gaf.c_str());
                        if (merged > 0) {
                            // The re-solve needs the recovery windows: allele_depths_call_het
                            // is scoped by them, and without it each haplotype takes its own
                            // majority allele and a genuine het collapses to one side.
                            Options solve_opts = opts;
                            // The depth-based het escape is keyed on the
                            // candidate's own provenance now
                            // (CandidateVariant::bam_injected), not on a window
                            // list, so nothing needs to be carried here. Window
                            // keying admitted every site in the window and cost
                            // 1.153% -> 2.067% read hamming on the default path.
                            // Two rounds, as the alignment pipeline solves: the
                            // clean sites set the gauge, then the merged in-gap
                            // sites -- NOISY_CAND_HET, which the clean mask does
                            // not admit -- join it. Anchored, so the noisy sites
                            // can extend a block without re-deciding the parity
                            // the catalog's clean sites already established.
                            // The merged sites belong to the FIRST solve, not to a
                            // correction applied after one. Both rounds run again
                            // over the union and neither is anchored: the chunk is
                            // solved once, with every site it will ever have.
                            //
                            // Anchoring here was measured and removed. On
                            // chr20:42,500,000-43,000,000, where recovery merges 6
                            // sites across 3 windows, pinning the pre-merge
                            // consensus left 840 of 2,008 reads misplaced inside a
                            // single block -- a switch, not an inversion: the pinned
                            // flank keeps one parity while the merged sites decide
                            // the other. Re-running both rounds unanchored places
                            // every read correctly (0 misplaced). Restricting the
                            // pin to the recovered intervals does NOT help (same 840),
                            // because the parity that has to change is the chunk's,
                            // not the gap's.
                            // Block import: the parent keeps the read labels it
                            // already has and only incorporates the imported
                            // sites. Re-solving from scratch is what every
                            // earlier admission mechanism did, and all of them
                            // degraded the result.
                            assign_hap_based_on_germline_het_vars_kmeans(
                                graph_chunks[offset].chunk, solve_opts, kCandGermlineClean,
                                opts.stitch_recovered);

                            // Stage 2, as the alignment arm runs it: the noisy-region MSA
                            // reconstructs the demoted loci, then the second round solves over
                            // clean plus the rebuilt noisy sites.
                            // The noisy-region MSA cannot run on the graph chunk itself:
                            // collect_noisy_reg_reads skips every read with no digars
                            // (collect_phase_noisy.cpp:1042) and graph-only reads have none, so
                            // a seeded region has zero reads and the MSA has nothing to align.
                            // The seeded regions go to the recovery instead, whose sub-solve
                            // re-reads the BAM through process_chunk and does build digars.
                            assign_hap_based_on_germline_het_vars_kmeans(
                                graph_chunks[offset].chunk, solve_opts, kCandGermlineVarCate,
                                false);
                        }
                    }
                }
            } catch (...) {
                std::lock_guard<std::mutex> lock(error_mutex);
                if (!first_error) first_error = std::current_exception();
            }
        });
    }
    for (std::thread& w : workers) w.join();
    if (first_error) std::rethrow_exception(first_error);

    if (opts.verbose >= 1) graph_query_report_match_stats();
    populate_graph_chunk_overlaps(graph_chunks);

    std::vector<PhasingChunk> phasing_chunks;
    phasing_chunks.reserve(batch_size);
    for (GraphChunkBuildResult& gc : graph_chunks)
        phasing_chunks.push_back(std::move(gc.chunk));
    stitch_chunk_haps(phasing_chunks, &opts, pgbam_sidecar);
    for (size_t i = 0; i < batch_size; ++i)
        graph_chunks[i].chunk = std::move(phasing_chunks[i]);

    return graph_chunks;
}

void run_collect_graph_variation(const Options& opts) {
    const bool use_indexed_gaf = !opts.gaf_file.empty();
    if (!use_indexed_gaf) {
        if (opts.gbz_db.empty())
            throw std::runtime_error("--gbz-db is required for collect-graph-variation without --gaf");
        if (opts.gaf_db.empty())
            throw std::runtime_error("--gaf-db is required for collect-graph-variation without --gaf");
    } else {
        // The --gaf path requires a tabix-indexed, bgzip-compressed GAF with
        // annotated coordinate columns so per-chunk region queries are efficient.
        require_indexed_gaf(opts.gaf_file);
        if (!opts.gbz_db.empty() || !opts.gaf_db.empty()) {
            std::cerr << "Warning: --gaf provided; ignoring --gbz-db/--gaf-db\n";
        }
    }

    // Load optional .pgbam sidecar for fallback chunk stitching.
    std::unique_ptr<PgbamSidecarData> pgbam_sidecar;
    if (!opts.pgbam_file.empty()) {
        pgbam_sidecar = std::make_unique<PgbamSidecarData>(load_pgbam_sidecar(opts.pgbam_file));
        if (opts.verbose >= 1)
            std::cerr << "Loaded pgbam sidecar with "
                      << pgbam_sidecar->set_to_threads.size() << " sets from "
                      << opts.pgbam_file << "\n";
    }

    // 1. Reference FASTA index first — needed to resolve autosome contig names for
    //    region filters, which are then passed to the VCF loader so only sites in the
    //    requested regions are parsed (tabix-assisted for bgzipped + indexed VCFs).
    std::unique_ptr<faidx_t, FaiDeleter> fai(load_reference_index(opts.ref_fasta));
    std::unique_ptr<bam_hdr_t, HeaderDeleter> header(build_synthetic_header(fai.get()));

    // 2. Build a chrom alias map between the FASTA and VCF naming conventions.
    //    Pangenome FASTAs use "SAMPLE#HAP#CHROM" (e.g. "CHM13#0#chr20") while graph VCFs
    //    use the plain reference name ("chr20"), or vice versa.  We inspect the FAI names
    //    and build a bidirectional suffix map so both directions resolve automatically.
    //
    //    fai_suffix_to_full : "chr20" → "CHM13#0#chr20"  (used when VCF is short, FAI is full)
    //    fai_full_to_suffix : "CHM13#0#chr20" → "chr20"  (used when VCF is full, FAI is short)
    std::unordered_map<std::string, std::string> fai_suffix_to_full;
    std::unordered_map<std::string, std::string> fai_full_to_suffix;
    {
        const int nseq = faidx_nseq(fai.get());
        fai_suffix_to_full.reserve(static_cast<size_t>(nseq));
        for (int i = 0; i < nseq; ++i) {
            const std::string full(faidx_iseq(fai.get(), i));
            const size_t h = full.rfind('#');
            if (h != std::string::npos) {
                fai_suffix_to_full.emplace(full.substr(h + 1), full);
                fai_full_to_suffix.emplace(full, full.substr(h + 1));
            }
        }
    }

    // Resolves a contig name to its canonical FAI form (the name that exists in the header).
    // Handles both "chr20" → "CHM13#0#chr20" and "CHM13#0#chr20" → "chr20".
    auto resolve_contig = [&](const std::string& name) -> std::string {
        if (faidx_has_seq(fai.get(), name.c_str())) return name;
        // Short name → full pangenome name ("chr20" → "CHM13#0#chr20")
        auto it = fai_suffix_to_full.find(name);
        if (it != fai_suffix_to_full.end()) return it->second;
        // Full pangenome name → short name ("CHM13#0#chr20" → "chr20")
        auto it2 = fai_full_to_suffix.find(name);
        if (it2 != fai_full_to_suffix.end() && faidx_has_seq(fai.get(), it2->second.c_str()))
            return it2->second;
        return name; // unchanged; add_filter_chunks will throw a clear error
    };

    std::vector<RegionFilter> region_filters;
    for (const std::string& r : opts.regions) {
        RegionFilter f = parse_region(r);
        f.chrom = resolve_contig(f.chrom);
        region_filters.push_back(std::move(f));
    }
    if (!opts.region_file.empty()) {
        auto bed = load_bed_regions(opts.region_file);
        for (RegionFilter& f : bed) f.chrom = resolve_contig(f.chrom);
        region_filters.insert(region_filters.end(), bed.begin(), bed.end());
    }
    if (opts.autosome) {
        // With pangenome FASTAs, "chr1"–"chr22" won't be found directly; resolve_contig
        // maps them to the full name (e.g. "CHM13#0#chr1").
        for (int i = 1; i <= 22; ++i) {
            for (const std::string& candidate :
                     {"chr" + std::to_string(i), std::to_string(i)}) {
                const std::string resolved = resolve_contig(candidate);
                if (faidx_has_seq(fai.get(), resolved.c_str())) {
                    region_filters.push_back(RegionFilter{true, resolved, 1, -1});
                    break;
                }
            }
        }
    }

    // When no region filters were specified, build whole-chromosome filters for every
    // FASTA contig so the VCF loader can use tabix instead of streaming the entire file.
    if (region_filters.empty()) {
        const int nseq = faidx_nseq(fai.get());
        region_filters.reserve(static_cast<size_t>(nseq));
        for (int i = 0; i < nseq; ++i)
            region_filters.push_back(
                RegionFilter{true, std::string(faidx_iseq(fai.get(), i)), 1, -1});
    }

    // 3. Build VCF-name → FAI-name chrom remap for per-chunk site normalization.
    //    Sites are loaded per-chunk via tabix, so no global catalog is needed.
    std::unordered_map<std::string, std::string> chrom_remap;
    {
        // VCF uses short name, FAI uses full ("chr20" → "CHM13#0#chr20")
        for (const auto& kv : fai_suffix_to_full) chrom_remap.emplace(kv.first, kv.second);
        // VCF uses full name, FAI uses short ("CHM13#0#chr20" → "chr20")
        for (const auto& kv : fai_full_to_suffix) {
            if (faidx_has_seq(fai.get(), kv.second.c_str()))
                chrom_remap.emplace(kv.first, kv.second);
        }
    }

    // 5. contig name → FAI-order tid (used when remapping candidate tids for output).
    std::unordered_map<std::string, int> contig_to_tid;
    contig_to_tid.reserve(static_cast<size_t>(header->n_targets));
    for (int32_t tid = 0; tid < header->n_targets; ++tid) {
        contig_to_tid[header->target_name[tid]] = tid;
    }

    // 6. Tile genome into RegionChunks using the resolved region_filters (which may have
    //    contig names like "CHM13#0#chr20" resolved from a user-supplied "chr20").
    // Alignment recovery, when a BAM is given. One context for the whole run:
    // it owns the BAM handles the targeted sub-solves reopen per window.
    std::unique_ptr<WorkerContext> bam_recovery_ctx;
    if (!opts.bam_files.empty()) {
        bam_recovery_ctx = std::make_unique<WorkerContext>(opts);
        std::cerr << "Alignment recovery enabled from " << opts.primary_bam_file() << "\n";
    }
    const std::vector<RegionChunk> chunks =
        build_region_chunks(opts, header.get(), fai.get(), region_filters);
    if (chunks.empty()) {
        std::cerr << "No region chunks to process\n";
        return;
    }
    if (opts.verbose >= 1) {
        std::cerr << "Tiled genome into " << chunks.size() << " chunks ("
                  << opts.chunk_size << " bp, " << opts.threads << " thread(s))\n";
    }

    // 7. Build per-chunk query config for the GBZ/GAF-base FFI path.
    GraphQueryConfig qconfig;
    qconfig.gbz_db   = opts.gbz_db;
    qconfig.gaf_db   = opts.gaf_db;
    qconfig.min_mapq  = opts.min_mapq;

    // Derive the reference sample name for GBZ interval queries.
    // For pangenome FASTAs ("CHM13#0#chr20") extract "CHM13".
    // For plain FASTAs ("chr20") leave empty so query uses its GENERIC_SAMPLE default.
    std::string ref_sample = opts.graph_sample;
    if (ref_sample.empty()) {
        const int nseq = faidx_nseq(fai.get());
        for (int i = 0; i < nseq && ref_sample.empty(); ++i) {
            const std::string name(faidx_iseq(fai.get(), i));
            const size_t h = name.find('#');
            if (h != std::string::npos && h > 0)
                ref_sample = name.substr(0, h);
        }
    }
    if (!use_indexed_gaf && opts.verbose >= 1 && !ref_sample.empty())
        std::cerr << "Using reference sample \"" << ref_sample << "\" for GBZ interval queries\n";

    // 8. Open output streams.
    std::ofstream variant_out(opts.output_tsv);
    if (!variant_out) throw std::runtime_error("failed to open output: " + opts.output_tsv);
    write_variants_tsv_header(variant_out);

    std::ofstream vcf_out;
    if (!opts.output_vcf.empty()) {
        vcf_out.open(opts.output_vcf);
        if (!vcf_out) throw std::runtime_error("failed to open VCF output: " + opts.output_vcf);
        write_variants_vcf_header(vcf_out, opts, header.get());
    }
    std::ofstream phased_vcf_out;
    if (!opts.output_phased_vcf.empty()) {
        phased_vcf_out.open(opts.output_phased_vcf);
        if (!phased_vcf_out)
            throw std::runtime_error("failed to open phased VCF output: " + opts.output_phased_vcf);
        write_phased_variants_vcf_header(phased_vcf_out, opts, header.get());
    }

    ReferenceCache ref(fai.get());

    // Phased BAM: open output file and write header.
    struct SamFileCloser { void operator()(samFile* fp) const { if (fp) hts_close(fp); } };
    std::unique_ptr<samFile, SamFileCloser> phased_bam_fp;
    std::unique_ptr<sam_hdr_t, decltype(&sam_hdr_destroy)> phased_bam_hdr(nullptr, &sam_hdr_destroy);
    std::unordered_map<std::string, PhaseReadOutputRow> phased_bam_rows;
    std::unordered_set<std::string> phased_bam_emitted;
    const bool emit_phased_bam = !opts.output_phased_bam.empty();
    if (emit_phased_bam) {
        samFile* raw_fp = hts_open(opts.output_phased_bam.c_str(), "wb");
        if (!raw_fp)
            throw std::runtime_error("failed to open phased BAM: " + opts.output_phased_bam);
        phased_bam_fp.reset(raw_fp);
        sam_hdr_t* hdr = sam_hdr_init();
        if (!hdr) throw std::runtime_error("failed to allocate phased BAM header");
        phased_bam_hdr.reset(hdr);
        if (sam_hdr_write(phased_bam_fp.get(), phased_bam_hdr.get()) < 0)
            throw std::runtime_error("failed to write phased BAM header");
    }

    // 10. Process chunks in reg_chunk_i batches (one contig per batch), streaming output.
    //     Mirrors run_collect_bam_variation's batch loop exactly.
    size_t n_variants = 0;
    size_t n_filtered = 0;
    // Diagnostic: why catalog sites never became candidates. Streamed alongside
    // the batch loop so a whole-chromosome run does not buffer millions of rows.
    std::unique_ptr<std::FILE, int (*)(std::FILE*)> filtered_out(nullptr, std::fclose);
    if (!opts.output_filtered_sites.empty()) {
        std::FILE* raw = std::fopen(opts.output_filtered_sites.c_str(), "w");
        if (raw == nullptr)
            throw std::runtime_error("failed to open filtered sites file: " +
                                     opts.output_filtered_sites);
        filtered_out.reset(raw);
        std::fprintf(filtered_out.get(),
                     "CHROM\tPOS\tSITE_ID\tREF_COV\tALT_COV\tTOTAL_COV\tAF\tREASON\n");
    }
    std::ofstream phase_sites_out;
    if (!opts.output_phase_sites.empty()) {
        phase_sites_out.open(opts.output_phase_sites);
        if (!phase_sites_out)
            throw std::runtime_error("failed to open phase sites file: " +
                                     opts.output_phase_sites);
        write_graph_phase_sites_tsv_header(phase_sites_out);
    }
    // Per-read phasing evidence, accumulated across chunks. Used to diagnose
    // which reads get a haplotype on thin or contradictory evidence: a read is
    // assigned by init_assign_read_hap with no minimum-observation or margin
    // requirement, so a single informative site is enough to commit it.
    struct PhaseReadDiag {
        int hap = 0;
        hts_pos_t phase_set = -1;
        int n_obs = 0;
        int agree = 0;
        int conflict = 0;
        int score_margin = 0;
        int n_scored = 0;
    };
    std::unordered_map<std::string, PhaseReadDiag> phase_read_diag;
    const bool emit_phase_reads = !opts.output_phase_reads.empty();

    size_t batch_begin = 0;
    while (batch_begin < chunks.size()) {
        size_t batch_end = batch_begin + 1;
        while (batch_end < chunks.size() &&
               chunks[batch_end].reg_chunk_i == chunks[batch_begin].reg_chunk_i) {
            ++batch_end;
        }

        std::vector<GraphChunkBuildResult> graph_chunks =
            use_indexed_gaf
                ? process_graph_chunk_batch_indexed_gaf(
                      opts.graph_sites_vcf, chunks, batch_begin, batch_end,
                      header.get(), opts.gaf_file, opts.min_mapq,
                      fai_full_to_suffix, chrom_remap, opts,
                      pgbam_sidecar.get())
                : process_graph_chunk_batch(
                      opts.graph_sites_vcf, chunks, batch_begin, batch_end,
                      header.get(), qconfig, ref_sample, fai_full_to_suffix,
                      chrom_remap, opts, pgbam_sidecar.get());


        if (emit_phase_reads) {
            for (const GraphChunkBuildResult& gc : graph_chunks) {
                const PhasingChunk& pc = gc.chunk;
                for (const ReadVariantProfile& profile : pc.read_var_profile) {
                    const size_t read_i = static_cast<size_t>(profile.read_id);
                    if (read_i >= pc.reads.size()) continue;
                    const ReadRecord& rr = pc.reads[read_i];
                    PhaseReadDiag& d = phase_read_diag[rr.qname];
                    int obs = 0;
                    for (int allele : profile.alleles)
                        if (allele >= 0) ++obs;
                    d.n_obs += obs;
                    d.agree += rr.n_clean_agree_snps;
                    d.conflict += rr.n_clean_conflict_snps;
                    if (rr.hap_score_margin > d.score_margin)
                        d.score_margin = rr.hap_score_margin;
                    if (rr.n_vars_scored > d.n_scored) d.n_scored = rr.n_vars_scored;
                    const int hap = read_i < pc.haps.size() ? pc.haps[read_i] : 0;
                    if (hap != 0) {
                        d.hap = hap;
                        d.phase_set = read_i < pc.phase_sets.size()
                                          ? pc.phase_sets[read_i]
                                          : static_cast<hts_pos_t>(-1);
                    }
                }
            }
        }

        for (const GraphChunkBuildResult& gc : graph_chunks) {
            n_filtered += gc.filtered_sites.size();
            if (filtered_out) {
                for (const FilteredGraphSite& fs : gc.filtered_sites) {
                    std::fprintf(filtered_out.get(), "%s\t%lld\t%s\t%d\t%d\t%d\t%.4f\t%s\n",
                                 fs.chrom.c_str(), static_cast<long long>(fs.pos),
                                 fs.site_id.c_str(), fs.ref_cov, fs.alt_cov,
                                 fs.total_cov, fs.allele_fraction,
                                 fs.filter_reason.c_str());
                }
            }
            if (phase_sites_out) {
                write_graph_phase_sites_tsv_rows(phase_sites_out, gc);
            }
        }

        CandidateTable variants =
            graph_chunks_to_candidate_table(graph_chunks, contig_to_tid, opts);

        n_variants += variants.size();
        write_variants_tsv_records(variant_out, header.get(), ref, variants);
        if (!opts.output_vcf.empty())
            write_variants_vcf_records(vcf_out, opts, header.get(), ref, variants);
        if (!opts.output_phased_vcf.empty())
            write_phased_variants_vcf_records(phased_vcf_out, opts, header.get(), ref, variants);

        // Phased BAM: accumulate post-stitch read assignments, then flush.
        // Batches are per-contig so cross-batch read overlap is impossible;
        // flush everything after each batch.
        if (emit_phased_bam) {
            for (const GraphChunkBuildResult& gc : graph_chunks)
                merge_graph_chunk_into_read_rows(phased_bam_rows, gc,
                                                 opts.min_read_hap_margin);
            flush_graph_phase_bam_after_merge(
                phased_bam_fp.get(), phased_bam_hdr.get(),
                phased_bam_rows, nullptr, phased_bam_emitted);
        }

        batch_begin = batch_end;
    }

    if (emit_phase_reads) {
        std::FILE* fp = std::fopen(opts.output_phase_reads.c_str(), "w");
        if (fp == nullptr)
            throw std::runtime_error("failed to open phase reads file: " +
                                     opts.output_phase_reads);
        std::fprintf(fp, "READ\tHAP\tPHASE_SET\tN_OBS\tCLEAN_AGREE\tCLEAN_CONFLICT\tSCORE_MARGIN\tN_SCORED\n");
        for (const auto& [qname, d] : phase_read_diag) {
            std::fprintf(fp, "%s\t%d\t%lld\t%d\t%d\t%d\t%d\t%d\n", qname.c_str(), d.hap,
                         static_cast<long long>(d.phase_set), d.n_obs, d.agree, d.conflict,
                         d.score_margin, d.n_scored);
        }
        std::fclose(fp);
        std::cerr << "Wrote per-read phasing evidence to " << opts.output_phase_reads << "\n";
    }

    std::cerr << "Processed " << chunks.size() << " region chunks with " << opts.threads
              << " worker thread(s)\n";
    std::cerr << "Collected " << n_variants << " candidate variant sites ("
              << n_filtered << " filtered) into " << opts.output_tsv << "\n";
    if (!opts.output_vcf.empty())
        std::cerr << "Wrote candidate VCF to " << opts.output_vcf << "\n";
    if (!opts.output_phased_vcf.empty())
        std::cerr << "Wrote phased candidate VCF to " << opts.output_phased_vcf << "\n";
    if (!opts.output_phase_sites.empty())
        std::cerr << "Wrote retained graph sites to " << opts.output_phase_sites << "\n";
    if (emit_phased_bam)
        std::cerr << "Wrote phased BAM to " << opts.output_phased_bam << "\n";
}

static void print_graph_collect_help() {
    std::cout
        << "Usage: pgphase collect-graph-variation [options]\n"
        << "\n"
        << "Required:\n"
        << "      --ref FILE                Reference FASTA (indexed)\n"
        << "      --sites FILE              Sites VCF from build-snarl-catalog (bgzipped + tabix-indexed)\n"
        << "\n"
        << "Read input:\n"
        << "      --gaf FILE                Coordinate-indexed GAF from pggaf (bgzipped + tabix-indexed)\n"
        << "      --gbz-db FILE             GBZ graph database (legacy --gaf-db path)\n"
        << "      --gaf-db FILE             GAF-base read alignment database (legacy path)\n"
        << "\n"
        << "Options:\n"
        << "  -o, --output FILE             Output TSV [output.tsv]\n"
        << "  -v, --vcf-output FILE         Candidate VCF output\n"
        << "      --phased-vcf-out FILE     Phased VCF with GT:DP:AD:VAF:GQ:PS\n"
        << "      --phased-bam-out FILE     Unaligned BAM with HP/PS tags per read\n"
        << "      --graph-noisy-msa         Run the alignment pipeline's noisy-region MSA\n"
        << "      --stitch-recovered        Import a recovery sub-solve as a block: orient it\n"
        << "                                against the parent on shared reads and keep its\n"
        << "                                per-site consensus instead of re-solving\n"
        << "      --recovery-audit-out FILE One row per candidate the recovery sub-solve\n"
        << "                                found, and what the merge did with it\n"
        << "                                over repeat-context loci (stage 2)\n"
        << "      --link-earned-repeat-indels  Re-admit a repeat-context het indel when it\n"
        << "                                agrees with a nearby clean het SNP on >= 15 reads\n"
        << "      --bam FILE                Indexed BAM used ONLY to recover what the graph\n"
        << "                                of grafting a sub-solve on afterwards\n"
        << "                                set, instead of refining stage 1 (pre-2026-09-18)\n"
        << "                                 sites could not phase: each unphased window and\n"
        << "                                 each seam between blocks is re-solved from the\n"
        << "                                 alignment as its own chunk and stitched in\n"
        << "      --filtered-sites-out FILE Diagnostic TSV of dropped catalog sites and why\n"
        << "      --phase-sites-out FILE    Diagnostic TSV of retained graph sites with SITE_ID\n"
        << "      --phase-reads-out FILE    Diagnostic TSV of per-read phasing evidence\n"
        << "      --graph-indel-af-margin F  Max |AF-0.5| for a het-indel k-means anchor [0.11]\n"
        << "      --graph-indel-min-alt INT  Min alt support for a het-indel k-means anchor [0]\n"
        << "      --min-read-margin INT     Min clean-SNP (agree-conflict) to phase a read [0=off]\n"
        << "      --stitch-min-margin INT   Abstain on chunk seams below this vote margin [0]\n"
        << "      --stitch-rule INT         0=net-margin 1=both-strands 2=literal 3=both+margin [0]\n"
        << "      --anchor-af-margin F      Max |AF-0.5| for a site to vote in k-means [0.5=off]\n"
        << "      --min-block-link-reads INT Spanning reads needed to carry a phase block [2]\n"
        << "      --block-link-window INT   Preceding het variants searched for that link [1]\n"
        << "      --link-by-alleles         Let untagged reads link variants by allele pattern\n"
        << "      --emit-nonanchor-hets     Emit/phase hets outside --anchor-af-margin (never anchor)\n"
        << "      --gaf-pad INT             Widen the per-chunk GAF read query by INT bp [0]\n"
        << "      --snarl-allele-phasing    Score multi-allelic snarls as alt-vs-other, not alt-vs-ref\n"
        << "      --snarl-keep-whole        Keep multi-allelic snarls as single n-allelic anchors\n"
        << "      --snarl-top2-frac FLOAT   Min read share on a snarl's top 2 alleles to anchor [0.9]\n"
        << "      --af-vs-site-depth        Score allele fraction against total site depth,\n"
        << "                                recovering hets between two non-reference alleles\n"
        << "  -t, --threads INT             Worker threads [1]\n"
        << "  -q, --min-mapq INT            Minimum read mapping quality [30]\n"
        << "  -D, --min-depth INT           Minimum total depth [5]\n"
        << "      --min-alt-depth INT       Minimum alt depth [2]\n"
        << "      --min-af FLOAT            Minimum allele fraction [0.20]\n"
        << "      --max-af FLOAT            Maximum allele fraction [0.80]\n"
        << "      --min-sv-len INT          Min SV length for SVTYPE/SVLEN tags [30]\n"
        << "      --chunk-size INT          Region chunk size in bp [500000]\n"
        << "  -r, --region STR              Restrict to region (may be repeated)\n"
        << "      --region-file FILE        BED file of regions\n"
        << "      --autosome                Process chr1-22 / 1-22 only\n"
        << "      --sample NAME             Reference sample name for GBZ interval queries\n"
        << "                                (auto-derived from FASTA if not provided)\n"

        << "      --hifi                    HiFi read mode [default]\n"
        << "      --ont                     ONT read mode (enables strand-bias filter)\n"
        << "      --strand-bias-pval FLOAT  Max p-value for ONT strand-bias filter [0.01]\n"
        << "\n"
        << "Pgbam stitching:\n"
        << "      --pgbam-file FILE         Optional .pgbam sidecar for fallback chunk stitching\n"
        << "      --pgbam-primary-margin INT         Thread polarity margin for primary stitching [2]\n"
        << "      --pgbam-primary-min-winning INT    Winning shared polarized threads for primary stitching [2]\n"
        << "      --no-pgbam-cleanup-pass            Disable final .pgbam cleanup pass\n"
        << "      --pgbam-cleanup-margin INT         Thread polarity margin for cleanup pass [2]\n"
        << "      --pgbam-cleanup-min-winning INT    Winning shared polarized threads for cleanup pass [1]\n"
        << "      --no-pgbam-relaxed-cleanup-pass    Disable relaxed .pgbam cleanup pass\n"
        << "      --pgbam-relaxed-cleanup-margin INT Thread polarity margin for relaxed cleanup [1]\n"
        << "      --pgbam-relaxed-cleanup-min-winning INT Winning threads for relaxed cleanup [1]\n"
        << "\n"
        << "  -V, --verbose INT             Verbosity level [0]\n"
        << "  -h, --help                    Print this help\n"
        << "\n"
        << "Examples:\n"
        << "  pgphase collect-graph-variation \\\n"
        << "      --ref ref.fa \\\n"
        << "      --sites sites.vcf.gz \\\n"
        << "      --gaf reads.gaf \\\n"
        << "      --phased-vcf-out phased.vcf \\\n"
        << "      --phased-bam-out phased.bam \\\n"
        << "      -t 8\n"
        << "\n"
        << "  pgphase collect-graph-variation \\\n"
        << "      --ref ref.fa \\\n"
        << "      --sites sites.vcf.gz \\\n"
        << "      --gbz-db reads.gaf.db \\\n"
        << "      --ont \\\n"
        << "      --phased-vcf-out phased.vcf \\\n"
        << "      -t 16\n";
}

enum GraphCollectOption {
    kGcRecoveryBam = 2000,
    kGcMinAltDepth = 1000,
    kGcMinAf,
    kGcMaxAf,
    kGcMinSvLen,
    kGcChunkSize,
    kGcPhasedVcf,
    kGcGbzDb,
    kGcGafFile,
    kGcGafDb,
    kGcRegionFile,
    kGcAutosome,
    kGcSample,
    kGcOnt,
    kGcHifi,
    kGcStrandBiasPval,
    kGcPhasedBam,
    kGcLinkEarnedRepeatIndels,
    kGcGraphNoisyMsa,
    kGcRecoveryAuditOut,
    kGcStitchRecovered,
    kGcRef,
    kGcSites,
    kGcPgbamFile,
    kGcPgbamPrimaryMargin,
    kGcPgbamPrimaryMinWinning,
    kGcNoPgbamCleanupPass,
    kGcPgbamCleanupMargin,
    kGcPgbamCleanupMinWinning,
    kGcNoPgbamRelaxedCleanupPass,
    kGcPgbamRelaxedCleanupMargin,
    kGcPgbamRelaxedCleanupMinWinning,
    kGcFilteredSitesOut,
    kGcPhaseSitesOut,
    kGcPhaseReadsOut,
    kGcGraphIndelAfMargin,
    kGcGraphIndelMinAlt,
    kGcMinReadHapMargin,
    kGcStitchMinMargin,
    kGcStitchRule,
    kGcAnchorAfMargin,
    kGcAfVsSiteDepth,
    kGcBlockLink,
    kGcBlockLinkWindow,
    kGcLinkByAlleles,
    kGcEmitNonAnchorHets,
    kGcGafPad,
    kGcSnarlAllelePhasing,
    kGcSnarlKeepWhole,
    kGcSnarlTop2Frac,
};

} // namespace

} // namespace pgphase_collect

int collect_graph_variation(int argc, char* argv[]) {
    using namespace pgphase_collect;
    Options opts;

    {
        std::ostringstream cmd;
        cmd << "pgphase collect-graph-variation";
        for (int i = 1; i < argc; ++i) cmd << ' ' << argv[i];
        opts.command_line = cmd.str();
    }

    optind = 1;
    const struct option long_options[] = {
        {"output",            required_argument, nullptr, 'o'},
        {"vcf-output",        required_argument, nullptr, 'v'},
        {"phased-vcf-out",    required_argument, nullptr, kGcPhasedVcf},
        {"phased-bam-out",   required_argument, nullptr, kGcPhasedBam},
        {"link-earned-repeat-indels", no_argument, nullptr, kGcLinkEarnedRepeatIndels},
        {"graph-noisy-msa", no_argument, nullptr, kGcGraphNoisyMsa},
        {"recovery-audit-out", required_argument, nullptr, kGcRecoveryAuditOut},
        {"stitch-recovered", no_argument, nullptr, kGcStitchRecovered},
        {"bam",              required_argument, nullptr, kGcRecoveryBam},
        {"filtered-sites-out", required_argument, nullptr, kGcFilteredSitesOut},
        {"phase-sites-out",   required_argument, nullptr, kGcPhaseSitesOut},
        {"phase-reads-out",   required_argument, nullptr, kGcPhaseReadsOut},
        {"graph-indel-af-margin", required_argument, nullptr, kGcGraphIndelAfMargin},
        {"graph-indel-min-alt",   required_argument, nullptr, kGcGraphIndelMinAlt},
        {"min-read-margin",   required_argument, nullptr, kGcMinReadHapMargin},
        {"stitch-min-margin", required_argument, nullptr, kGcStitchMinMargin},
        {"stitch-rule",       required_argument, nullptr, kGcStitchRule},
        {"anchor-af-margin",  required_argument, nullptr, kGcAnchorAfMargin},
        {"af-vs-site-depth",  no_argument,       nullptr, kGcAfVsSiteDepth},
        {"min-block-link-reads", required_argument, nullptr, kGcBlockLink},
        {"block-link-window",    required_argument, nullptr, kGcBlockLinkWindow},
        {"link-by-alleles",      no_argument,       nullptr, kGcLinkByAlleles},
        {"emit-nonanchor-hets",  no_argument,       nullptr, kGcEmitNonAnchorHets},
        {"gaf-pad",              required_argument, nullptr, kGcGafPad},
        {"snarl-allele-phasing", no_argument,       nullptr, kGcSnarlAllelePhasing},
        {"snarl-keep-whole",     no_argument,       nullptr, kGcSnarlKeepWhole},
        {"snarl-top2-frac",      required_argument, nullptr, kGcSnarlTop2Frac},
        {"threads",           required_argument, nullptr, 't'},
        {"min-mapq",          required_argument, nullptr, 'q'},
        {"min-depth",         required_argument, nullptr, 'D'},
        {"min-alt-depth",     required_argument, nullptr, kGcMinAltDepth},
        {"min-af",            required_argument, nullptr, kGcMinAf},
        {"max-af",            required_argument, nullptr, kGcMaxAf},
        {"min-sv-len",        required_argument, nullptr, kGcMinSvLen},
        {"chunk-size",        required_argument, nullptr, kGcChunkSize},
        {"region",            required_argument, nullptr, 'r'},
        {"region-file",       required_argument, nullptr, kGcRegionFile},
        {"autosome",          no_argument,       nullptr, kGcAutosome},
        {"gaf",               required_argument, nullptr, kGcGafFile},
        {"gaf-file",          required_argument, nullptr, kGcGafFile},
        {"gbz-db",            required_argument, nullptr, kGcGbzDb},
        {"gaf-db",            required_argument, nullptr, kGcGafDb},
        {"sample",            required_argument, nullptr, kGcSample},
        {"hifi",              no_argument,       nullptr, kGcHifi},
        {"ont",               no_argument,       nullptr, kGcOnt},
        {"strand-bias-pval",  required_argument, nullptr, kGcStrandBiasPval},
        {"ref",               required_argument, nullptr, kGcRef},
        {"sites",             required_argument, nullptr, kGcSites},
        {"pgbam-file",                required_argument, nullptr, kGcPgbamFile},
        {"pgbam-primary-margin",      required_argument, nullptr, kGcPgbamPrimaryMargin},
        {"pgbam-primary-min-winning", required_argument, nullptr, kGcPgbamPrimaryMinWinning},
        {"no-pgbam-cleanup-pass",     no_argument,       nullptr, kGcNoPgbamCleanupPass},
        {"pgbam-cleanup-margin",      required_argument, nullptr, kGcPgbamCleanupMargin},
        {"pgbam-cleanup-min-winning", required_argument, nullptr, kGcPgbamCleanupMinWinning},
        {"no-pgbam-relaxed-cleanup-pass", no_argument,   nullptr, kGcNoPgbamRelaxedCleanupPass},
        {"pgbam-relaxed-cleanup-margin", required_argument, nullptr, kGcPgbamRelaxedCleanupMargin},
        {"pgbam-relaxed-cleanup-min-winning", required_argument, nullptr, kGcPgbamRelaxedCleanupMinWinning},
        {"verbose",           required_argument, nullptr, 'V'},
        {"help",              no_argument,       nullptr, 'h'},
        {nullptr, 0, nullptr, 0}
    };

    int opt = 0;
    int long_index = 0;
    while ((opt = getopt_long(argc, argv, "o:v:t:q:D:r:V:h", long_options, &long_index)) != -1) {
        switch (opt) {
            case 'o': opts.output_tsv = optarg; break;
            case 'v': opts.output_vcf = optarg; break;
            case kGcPhasedVcf:    opts.output_phased_vcf = optarg; break;
            case kGcPhasedBam:    opts.output_phased_bam = optarg; break;
            case kGcLinkEarnedRepeatIndels: opts.link_earned_repeat_indels = true; break;
            case kGcGraphNoisyMsa: opts.graph_noisy_msa = true; break;
            case kGcRecoveryAuditOut: opts.recovery_audit_out = optarg; break;
            case kGcStitchRecovered: opts.stitch_recovered = true; break;
            // Recovery only. The graph pass never reads this BAM; it is used to
            // re-solve the intervals the catalog's sites could not phase.
            case kGcRecoveryBam:  opts.bam_files.push_back(optarg); break;

            case kGcFilteredSitesOut: opts.output_filtered_sites = optarg; break;
            case kGcPhaseSitesOut: opts.output_phase_sites = optarg; break;
            case kGcPhaseReadsOut: opts.output_phase_reads = optarg; break;
            case kGcGraphIndelAfMargin:
                opts.graph_indel_af_margin = parse_double_arg(optarg, "--graph-indel-af-margin");
                break;
            case kGcGraphIndelMinAlt:
                opts.graph_indel_min_alt = parse_int_arg(optarg, "--graph-indel-min-alt");
                break;
            case kGcMinReadHapMargin:
                opts.min_read_hap_margin = parse_int_arg(optarg, "--min-read-margin");
                break;
            case kGcStitchMinMargin:
                opts.stitch_min_margin = parse_int_arg(optarg, "--stitch-min-margin");
                break;
            case kGcStitchRule:
                opts.stitch_rule = parse_int_arg(optarg, "--stitch-rule");
                break;
            case kGcAnchorAfMargin:
                opts.anchor_af_margin = parse_double_arg(optarg, "--anchor-af-margin");
                break;
            case kGcAfVsSiteDepth: opts.af_vs_site_depth = true; break;
            case kGcBlockLink:
                opts.min_block_link_reads = parse_int_arg(optarg, "--min-block-link-reads");
                break;
            case kGcBlockLinkWindow:
                opts.block_link_window = parse_int_arg(optarg, "--block-link-window");
                break;
            case kGcLinkByAlleles: opts.link_by_alleles = true; break;
            case kGcEmitNonAnchorHets: opts.emit_nonanchor_hets = true; break;
            case kGcGafPad: opts.gaf_pad = parse_int_arg(optarg, "--gaf-pad"); break;
            case kGcSnarlAllelePhasing: opts.snarl_allele_phasing = true; break;
            case kGcSnarlKeepWhole: opts.snarl_keep_whole = true; opts.snarl_allele_phasing = true; break;
            case kGcSnarlTop2Frac: opts.snarl_top2_frac = parse_double_arg(optarg, "--snarl-top2-frac"); break;
            case 't': opts.threads = parse_int_arg(optarg, "--threads"); break;
            case 'q': opts.min_mapq = parse_int_arg(optarg, "--min-mapq"); break;
            case 'D': opts.min_depth = parse_int_arg(optarg, "--min-depth"); break;
            case kGcMinAltDepth:  opts.min_alt_depth = parse_int_arg(optarg, "--min-alt-depth"); break;
            case kGcMinAf:        opts.min_af = parse_double_arg(optarg, "--min-af"); break;
            case kGcMaxAf:        opts.max_af = parse_double_arg(optarg, "--max-af"); break;
            case kGcMinSvLen:     opts.min_sv_len = parse_int_arg(optarg, "--min-sv-len"); break;
            case kGcChunkSize:    opts.chunk_size = parse_ll_arg(optarg, "--chunk-size"); break;
            case 'r': opts.regions.push_back(optarg); break;
            case kGcRegionFile:   opts.region_file = optarg; break;
            case kGcAutosome:     opts.autosome = true; break;
            case kGcGafFile:      opts.gaf_file = optarg; break;
            case kGcGbzDb:        opts.gbz_db = optarg; break;
            case kGcGafDb:        opts.gaf_db = optarg; break;
            case kGcSample:       opts.graph_sample = optarg; break;
            case kGcHifi:         opts.read_technology = ReadTechnology::Hifi; break;
            case kGcOnt:          opts.read_technology = ReadTechnology::Ont; break;
            case kGcStrandBiasPval: opts.strand_bias_pval = parse_double_arg(optarg, "--strand-bias-pval"); break;
            case kGcRef:          opts.ref_fasta = optarg; break;
            case kGcSites:        opts.graph_sites_vcf = optarg; break;
            case kGcPgbamFile:      opts.pgbam_file = optarg; break;
            case kGcPgbamPrimaryMargin: opts.pgbam_primary_polarity_margin = parse_int_arg(optarg, "--pgbam-primary-margin"); break;
            case kGcPgbamPrimaryMinWinning: opts.pgbam_primary_min_winning_threads = parse_int_arg(optarg, "--pgbam-primary-min-winning"); break;
            case kGcNoPgbamCleanupPass: opts.pgbam_cleanup_pass = false; break;
            case kGcPgbamCleanupMargin: opts.pgbam_cleanup_polarity_margin = parse_int_arg(optarg, "--pgbam-cleanup-margin"); break;
            case kGcPgbamCleanupMinWinning: opts.pgbam_cleanup_min_winning_threads = parse_int_arg(optarg, "--pgbam-cleanup-min-winning"); break;
            case kGcNoPgbamRelaxedCleanupPass: opts.pgbam_relaxed_cleanup_pass = false; break;
            case kGcPgbamRelaxedCleanupMargin: opts.pgbam_relaxed_cleanup_polarity_margin = parse_int_arg(optarg, "--pgbam-relaxed-cleanup-margin"); break;
            case kGcPgbamRelaxedCleanupMinWinning: opts.pgbam_relaxed_cleanup_min_winning_threads = parse_int_arg(optarg, "--pgbam-relaxed-cleanup-min-winning"); break;
            case 'V': opts.verbose = parse_int_arg(optarg, "--verbose"); break;
            case 'h': print_graph_collect_help(); return 0;
            default:  print_graph_collect_help(); return 1;
        }
    }

    if (opts.ref_fasta.empty() || opts.graph_sites_vcf.empty()) {
        std::cerr << "Error: --ref and --sites are required\n";
        print_graph_collect_help();
        return 1;
    }

    require_graph_site_vcf_tabix_index(opts.graph_sites_vcf);

    if (opts.gaf_file.empty() && (opts.gbz_db.empty() || opts.gaf_db.empty())) {
        std::cerr << "Error: provide --gaf, or provide --gbz-db and --gaf-db\n";
        print_graph_collect_help();
        return 1;
    }

    if (opts.threads < 1 || opts.min_mapq < 0 || opts.min_depth < 0 ||
        opts.min_alt_depth < 0 || opts.min_af < 0.0 || opts.max_af < opts.min_af ||
        opts.min_sv_len < 0 || opts.chunk_size < 1 || opts.verbose < 0 ||
        opts.strand_bias_pval < 0.0 || opts.strand_bias_pval > 1.0) {
        std::cerr << "Error: invalid numeric threshold\n";
        return 1;
    }

    try {
        run_collect_graph_variation(opts);
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << "\n";
        return 1;
    }
    return 0;
}
