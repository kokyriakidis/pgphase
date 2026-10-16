#include "graph_collect.hpp"
#include "collect_pipeline.hpp"

#include "arg_parse.hpp"
#include "collect_output.hpp"
#include "collect_phase.hpp"
#include "collect_phase_noisy.hpp"
#include "collect_types.hpp"
#include "collect_var.hpp"
#include "graph_bam_adapter.hpp"
#include "graph_query.hpp"
#include "graph_sites.hpp"
#include "noise_filter.hpp"
#include "union_phase.hpp"

#include "edlib.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cctype>
#include <atomic>
#include <cstdlib>
#include <iomanip>
#include <chrono>
#include <cstdio>
#include <cstdint>
#include <fstream>
#include <functional>
#include <getopt.h>
#include <iostream>
#include <memory>
#include <limits>
#include <map>
#include <set>
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

namespace {

// Graph alignments encode path support directly, and the chr20 gap audit
// found useful linking reads between MAPQ 5 and 29. Keep this separate
// from longcallD/BAM, whose established default remains MAPQ 30.
constexpr int kDefaultGraphMinMapq = 5;
// One megabase gives the graph and its independent BAM solve enough shared
// context to stabilize read assignments while retaining bounded worker memory.
constexpr hts_pos_t kDefaultGraphBamChunkSize = 1000000;

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
            if (meta.bam_alt_deletion_no_ref) continue;

            // A site injected by union gap phasing (no catalog id) is
            // written only when it carries phase information. The alignment's
            // in-gap discovery also calls homozygous variants, and emitting
            // those would add ~12,850 hom records to chr20 -- a 23% larger VCF
            // whose extra content is alignment-discovered calls appearing ONLY
            // where the alignment solve found sites, which is a biased subset of
            // the genome. The output is the catalog's sites plus what union
            // phasing phased, not a variant call set.
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
                // Output reconstruction must retain the injection provenance used by
                // duplicate selection. Otherwise a verified noisy heterozygote
                // loses to a deeper unphased graph repeat at the same allele.
                cand.bam_injected = mcand.bam_injected;
                cand.msa_verified = mcand.msa_verified;
                cand.alignment_verified = mcand.alignment_verified;

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
                // Physical validation can establish HOM ALT despite graph REF
                // observations; retain that decision at biallelic emission.
                const bool validated_hom_alt =
                    mcand.counts.category == VariantCategory::CleanHom &&
                    hcons.size() > 2 && hcons[1] == new_a && hcons[2] == new_a;
                const bool is_hom_alt = !het_by_consensus &&
                    ((ref_cov == 0 && alt_cov >= opts.min_alt_depth) || validated_hom_alt);
                // Recovery already classified an MSA heterozygote from its
                // own reads and haplotype consensus. Its REF/ALT depth can be
                // very uneven; reapplying graph discovery thresholds here
                // turns a valid phased source row into LOW_AF at emission.
                const bool verified_injected_het = mcand.bam_injected &&
                    mcand.counts.category == VariantCategory::NoisyCandHet &&
                    mcand.msa_verified && mcand.alignment_verified &&
                    mcand.phase_set > 0 && het_by_consensus;
                if (verified_injected_het) {
                    cand.counts.category = mcand.counts.category;
                    cand.counts.candvarcate_initial =
                        mcand.counts.candvarcate_initial;
                    cand.lcd_var_i_to_cate = mcand.lcd_var_i_to_cate;
                } else if (alt_cov < opts.min_alt_depth || total_cov < opts.min_depth) {
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
                    // Homozygous candidates have no phase block. Keep the BAM
                    // candidate convention (0) rather than the read sentinel
                    // (-1), so both pipelines serialize the same state.
                    cand.phase_set = kUnsetCandidatePhaseSet;
                } else {
                    cand.hap_to_cons_alle[1] = (hap1 == new_a) ? 1 : 0;
                    cand.hap_to_cons_alle[2] = (hap2 == new_a) ? 1 : 0;
                    cand.hap_alt = cand.hap_to_cons_alle[1];
                    cand.hap_ref = cand.hap_to_cons_alle[2];
                    cand.phase_set = mcand.phase_set;
                }
                cand.alt_ref_base = 4;  // use FASTA anchor (BAM-path default)
                cand.lcd_make_variants_region_pass = true;

                // A site injected by union gap phasing is written only when
                // this writer's own classification calls it a het. Two reasons.
                // The reclassification above sets CleanHom whenever ref_cov == 0
                // (is_hom_alt, line ~196), and an alignment candidate merged from
                // a gap often carries ref_cov = 0 -- measured, ref_cov 0 against
                // alt_cov 57 -- so the depth synthesis manufactures hom calls:
                // 408 extra 1|1 records on chr20 against 125 in the default. And
                // a hom site carries no phase information in any case, so a VCF
                // whose contract is the catalog's sites plus what union phasing
                // phased has no reason to gain hom calls that appear only where
                // the alignment solve found sites.
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

    // The same physical variant can come from a catalog snarl and a verified
    // injected BAM row. Retain its phased heterozygote when the other copy is
    // an unphased repeat call; coverage alone would discard the injected
    // phase whenever the graph row had a few more observations.
    const auto phased_het = [](const CandidateVariant& candidate) {
        const auto category = candidate.counts.category;
        return candidate.phase_set > 0 &&
            candidate.hap_to_cons_alle[1] >= 0 &&
            candidate.hap_to_cons_alle[2] >= 0 &&
            candidate.hap_to_cons_alle[1] != candidate.hap_to_cons_alle[2] &&
            (category == VariantCategory::CleanHetSnp ||
             category == VariantCategory::CleanHetIndel ||
             (category == VariantCategory::NoisyCandHet &&
              candidate.msa_verified && candidate.alignment_verified));
    };
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
                const bool write_phased = phased_het(result[write]);
                const bool read_phased = phased_het(result[read]);
                if ((read_phased && !write_phased) ||
                    (read_phased == write_phased &&
                     result[read].counts.total_cov >
                         result[write].counts.total_cov)) {
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

    make_colocated_alleles_complementary(result, opts.min_alt_depth);
    drop_conflicting_haplotype_alleles(result);

    return result;
}

// A graph SNP can look heterozygous when the physical molecules carry its ALT
// base on one haplotype and a deletion over its REF base on the other. The
// graph's REF/ALT labels then give the wrong read orientation. Check only
// graph-clean, biallelic SNPs and require decisive physical evidence before
// reclassifying homozygous ALT or excluding an ALT/deletion contrast. This uses the existing per-thread BAM
// handles and scans each chunk once rather than seeking for every site.
static bool exclude_ref_absent_graph_snps(
        GraphSiteCatalog& catalog, const GraphChunkBuildResult& built,
        const RegionChunk& region, const std::string& contig,
        WorkerContext& context, const Options& opts) {
    struct SiteEvidence {
        hts_pos_t pos;
        size_t catalog_index;
        char ref;
        char alt;
        int ref_count = 0;
        int alt_count = 0;
        int deletion_count = 0;
        int other_count = 0;
        bool ref_absence_eligible = false;
    };
    std::unordered_map<std::string, size_t> catalog_index;
    catalog_index.reserve(catalog.sites.size());
    for (size_t i = 0; i < catalog.sites.size(); ++i)
        catalog_index.emplace(graph_site_key_str(catalog.sites[i]), i);

    std::unordered_map<std::string, size_t> retained_candidates;
    std::unordered_map<std::string, size_t> retained_snps;
    for (const CandidateVariant& candidate : built.chunk.candidates) {
        ++retained_candidates[candidate.key.alt];
        if (candidate.counts.category == VariantCategory::CleanHetSnp)
            ++retained_snps[candidate.key.alt];
    }
    std::vector<SiteEvidence> sites;
    for (size_t i = 0; i < built.chunk.candidates.size(); ++i) {
        if (built.chunk.candidates[i].counts.category != VariantCategory::CleanHetSnp)
            continue;
        const auto found = catalog_index.find(built.chunk.candidates[i].key.alt);
        if (found == catalog_index.end()) continue;

        const std::string* alt = selected_graph_candidate_alt(built, i);
        if (alt == nullptr) continue;
        const GraphSiteMeta& meta = built.site_meta[i];
        const VariantKey key = vcf_to_variant_key(
            built.chunk.candidates[i].key.tid, meta.pos, meta.ref, *alt);
        if (key.type != VariantType::Snp || key.ref_len != 1 || key.alt.size() != 1 ||
            key.pos < meta.pos || key.pos - meta.pos >= static_cast<hts_pos_t>(meta.ref.size())) continue;
        SiteEvidence evidence{key.pos, found->second, meta.ref[key.pos - meta.pos], key.alt[0]};
        evidence.ref_absence_eligible = built.site_ids[i] == graph_site_key_str(catalog.sites[found->second]) &&
            retained_candidates.at(built.chunk.candidates[i].key.alt) == 1;
        sites.push_back(evidence);
    }
    if (sites.empty()) return false;
    std::sort(sites.begin(), sites.end(),
              [](const SiteEvidence& a, const SiteEvidence& b) {
                  return a.pos < b.pos;
              });

    // Graph GAF phasing admits MAPQ 5, but rejecting a graph allele needs
    // the BAM pipeline's high-confidence alignment floor.
    constexpr int kMinPhysicalValidationMapq = 30;
    const int min_mapq = std::max(opts.min_mapq, kMinPhysicalValidationMapq);
    for (size_t input = 0; input < context.bams.size(); ++input) {
        const int tid = sam_hdr_name2tid(context.headers[input].get(), contig.c_str());
        if (tid < 0) continue;
        std::unique_ptr<hts_itr_t, decltype(&hts_itr_destroy)> itr(
            sam_itr_queryi(context.indexes[input].get(), tid,
                           region.beg - 1, region.end), hts_itr_destroy);
        if (!itr) throw std::runtime_error("failed to query BAM for graph SNP validation: " + contig);
        std::unique_ptr<bam1_t, decltype(&bam_destroy1)> read(bam_init1(), bam_destroy1);
        if (!read) throw std::runtime_error("failed to allocate BAM record for graph SNP validation");
        int status = 0;
        while ((status = sam_itr_next(context.bams[input]->get(), itr.get(), read.get())) >= 0) {
            const bam1_core_t& core = read->core;
            if ((core.flag & (BAM_FUNMAP | BAM_FSECONDARY | BAM_FSUPPLEMENTARY)) ||
                (!opts.include_filtered && (core.flag & (BAM_FQCFAIL | BAM_FDUP))) ||
                core.qual < min_mapq || core.qual == 255)
                continue;
            hts_pos_t ref_pos = core.pos + 1;
            int query_pos = 0;
            const uint32_t* cigar = bam_get_cigar(read.get());
            const uint8_t* sequence = bam_get_seq(read.get());
            const uint8_t* quality = bam_get_qual(read.get());
            for (uint32_t op_i = 0; op_i < core.n_cigar; ++op_i) {
                const int op = bam_cigar_op(cigar[op_i]);
                const int length = bam_cigar_oplen(cigar[op_i]);
                const int consumed = bam_cigar_type(op);
                if (consumed & 2) {
                    auto it = std::lower_bound(
                        sites.begin(), sites.end(), ref_pos,
                        [](const SiteEvidence& site, hts_pos_t pos) {
                            return site.pos < pos;
                        });
                    for (; it != sites.end() && it->pos < ref_pos + length; ++it) {
                        if (op == BAM_CDEL) {
                            ++it->deletion_count;
                        } else if (op == BAM_CMATCH || op == BAM_CEQUAL || op == BAM_CDIFF) {
                            const int qi = query_pos + static_cast<int>(it->pos - ref_pos);
                            if (qi < 0 || qi >= core.l_qseq || quality[qi] < opts.min_bq)
                                continue;
                            const char base = seq_nt16_str[bam_seqi(sequence, qi)];
                            if (base == it->ref) ++it->ref_count;
                            else if (base == it->alt) ++it->alt_count;
                            else ++it->other_count;
                        }
                    }
                    ref_pos += length;
                }
                if (consumed & 1) query_pos += length;
            }
        }
        if (status < -1)
            throw std::runtime_error("failed to read BAM for graph SNP validation: " + contig);
    }

    const size_t ref_absence_sites = static_cast<size_t>(std::count_if(sites.begin(), sites.end(),
        [](const SiteEvidence& site) { return site.ref_absence_eligible; }));
    bool excluded = false;
    for (const SiteEvidence& site : sites) {
        const bool padded = catalog.sites[site.catalog_index].ref.size() > 1;
        GraphSite& physical_site = catalog.sites[site.catalog_index];

        const bool mixed_indel = std::any_of(physical_site.alts.begin(), physical_site.alts.end(),
            [&physical_site](const std::string& alt) { return alt.size() != physical_site.ref.size(); });
        const std::string physical_key = graph_site_key_str(physical_site);
        // This binary physical test validates a two-ALT decomposition;
        // additional selected branches need their own joint validation.
        if (padded && mixed_indel && retained_candidates.at(physical_key) == 2 &&
            retained_snps.at(physical_key) == 1 &&
            graph_snp_low_alt_fraction_supported(
                site.ref_count, site.alt_count, site.deletion_count, site.other_count,
                opts.min_af, sites.size())) {
            physical_site.bam_low_fraction_snp = true;
            excluded = true;
            continue;
        }
        if (!site.ref_absence_eligible) continue;
        const bool significant_deletion = padded && graph_snp_padded_deletion_supported(
            site.ref_count, site.alt_count, site.deletion_count, site.other_count, ref_absence_sites);
        if (padded ? !significant_deletion :
            (catalog.sites[site.catalog_index].alts.size() != 1 ||
             !graph_snp_ref_absence_supported(site.ref_count, site.alt_count,
                site.deletion_count, site.other_count, ref_absence_sites))) continue;
        GraphSite& graph_site = catalog.sites[site.catalog_index];
        if (site.deletion_count == 0) {
            graph_site.bam_homozygous_alt = true;
        } else if (padded) {
            graph_site.bam_alt_deletion_no_ref = true;
        } else {
            graph_site.eligible = false;
            graph_site.skip_reason = "bam_alt_deletion_no_ref";
        }
        excluded = true;
        if (opts.verbose)
            std::cerr << "graph: validated REF-absent SNP " << contig << ':' << site.pos
                      << " (REF " << site.ref_count << ", ALT " << site.alt_count
                      << ", deletion " << site.deletion_count << ")\n";
    }
    return excluded;
}

// One chunk's GAF rows for the sites in view, from one worker thread's own
// handles (opened once per thread; see run_graph_chunk_batch).
using ChunkRowFetcher =
    std::function<std::vector<GraphReadAllele>(const RegionChunk&, const GraphSiteCatalogView&)>;

// One chunk, start to finish: its sites (tabix), its GAF rows, the graph
// chunk, the noise filter, union gap phasing's injected sites (with --bam),
// the clean solve and the union solve. Diagnostic, off by default:
// PGPHASE_TIMING=1 prints the chunk's seconds per stage to stderr.
static void process_graph_chunk(GraphChunkBuildResult& out,
                                const RegionChunk& region,
                                const std::string& contig,
                                SiteCatalogHandle& sites_handle,
                                faidx_t* fai,
                                WorkerContext* bam_ctx,
                                const std::unordered_map<std::string, std::string>& chrom_remap,
                                const Options& opts,
                                const ChunkRowFetcher& fetch_rows) {
    static const bool timing = std::getenv("PGPHASE_TIMING") != nullptr;
    using Clock = std::chrono::steady_clock;
    std::vector<double> stage_s;
    Clock::time_point stage_t = Clock::now();
    const auto lap = [&]() {
        if (!timing) return;
        const Clock::time_point now = Clock::now();
        stage_s.push_back(std::chrono::duration<double>(now - stage_t).count());
        stage_t = now;
    };

    GraphSiteCatalog chunk_catalog = load_sites_for_region(sites_handle, contig, region.beg, region.end);
    lap();
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
    if (!chunk_view.empty()) chunk_rows = fetch_rows(region, chunk_view);
    lap();

    out = build_graph_chunk(chunk_view, chunk_rows, contig, region.beg - 1, region.end, region.chunk_id, opts);
    if (bam_ctx != nullptr && exclude_ref_absent_graph_snps(chunk_catalog, out, region, contig, *bam_ctx, opts))
        out = build_graph_chunk(chunk_view, chunk_rows, contig, region.beg - 1, region.end, region.chunk_id, opts);

    std::string reference_slice;
    // Noise filter: fetch reference slice and reclassify indels in
    // homopolymer/repeat/low-complexity contexts.
    {
        hts_pos_t ref_len = 0;
        char* ref_raw = faidx_fetch_seq64(fai, contig.c_str(), region.beg - 1, region.end - 1, &ref_len);
        if (ref_raw && ref_len > 0) {
            reference_slice.assign(ref_raw, ref_raw + ref_len);
            std::free(ref_raw);
            apply_graph_noise_filter(out, reference_slice, region.beg, region.beg + ref_len - 1,
                                     opts.noisy_reg_max_xgaps);
        } else {
            std::free(ref_raw);
        }
    }
    lap();

    // With --bam, the alignment's sample-specific hets join the catalog's
    // before the clean solve (union gap phasing).
    const bool union_phasing = bam_ctx != nullptr;
    std::unique_ptr<PhasingChunk> union_bam;
    if (union_phasing) inject_alignment_private_sites(out, opts, *bam_ctx, contig.c_str(), &union_bam);
    lap();

    phase_joint_graph_candidates(out, opts, [&](hts_pos_t pos) {
        const hts_pos_t index = pos - region.beg;
        return index >= 0 && static_cast<size_t>(index) < reference_slice.size()
            ? reference_slice[static_cast<size_t>(index)] : 'N';
    });
    reclassify_physically_validated_graph_snps(chunk_view, out);
    if (bam_ctx != nullptr) supplement_phased_snp_branches(chunk_view, chunk_rows, out, opts);
    lap();
    if (union_phasing) phase_chunk_with_alignment_sites(out, union_bam.get());
    lap();
    if (timing) {
        std::ostringstream line;
        line << "[timing] chunk " << region.beg << " sites " << chunk_catalog.sites.size()
             << " gaf_rows " << chunk_rows.size() << " seconds";
        static const char* const names[] = {"catalog", "gaf", "build", "inject", "clean", "union"};
        double total = 0.0;
        for (size_t i = 0; i < stage_s.size() && i < 6; ++i) {
            line << ' ' << names[i] << ' ' << std::fixed << std::setprecision(2) << stage_s[i];
            total += stage_s[i];
        }
        line << " total " << total << '\n';
        std::cerr << line.str();
    }
}

// Processes one batch of graph chunks in parallel (one thread pool per
// reg_chunk_i batch, mirroring collect_chunk_batch_parallel in
// collect_pipeline.cpp): each worker opens its own row source (open_rows),
// sites VCF, reference index and, with --bam, BAM context, then takes chunks
// until none are left. Peak memory = threads x (reads_per_chunk +
// sites_per_chunk). After all workers join: populate_graph_chunk_overlaps,
// stitch_chunk_haps, the rescue and the island labels.
static std::vector<GraphChunkBuildResult> run_graph_chunk_batch(
    const std::string& sites,
    const std::unordered_map<std::string, hts_pos_t>& max_span,
    const std::vector<RegionChunk>& chunks,
    size_t batch_begin,
    size_t batch_end,
    const std::string& batch_contig,
    const std::unordered_map<std::string, std::string>& chrom_remap,
    const Options& opts,
    const std::function<ChunkRowFetcher()>& open_rows)
{
    const size_t batch_size = batch_end - batch_begin;
    std::vector<GraphChunkBuildResult> graph_chunks(batch_size);
    const size_t worker_count = std::min<size_t>(static_cast<size_t>(opts.threads), batch_size);
    std::atomic<size_t> next_offset{0};
    std::exception_ptr first_error;
    std::mutex error_mutex;
    std::vector<std::thread> workers;
    workers.reserve(worker_count);

    for (size_t w = 0; w < worker_count; ++w) {
        workers.emplace_back([&]() {
            try {
                const ChunkRowFetcher fetch_rows = open_rows();
                SiteCatalogHandle sites_handle(sites);
                sites_handle.max_span = max_span;
                std::unique_ptr<faidx_t, FaiDeleter> thread_fai(load_reference_index(opts.ref_fasta));
                // Per-thread BAM context for union gap phasing. htslib handles are
                // not shareable, so each worker opens its own set, once per thread
                // rather than per chunk.
                std::unique_ptr<WorkerContext> thread_bam_ctx;
                if (!opts.bam_files.empty()) thread_bam_ctx = std::make_unique<WorkerContext>(opts);
                while (true) {
                    const size_t offset = next_offset.fetch_add(1);
                    if (offset >= batch_size) break;
                    process_graph_chunk(graph_chunks[offset], chunks[batch_begin + offset], batch_contig,
                                        sites_handle, thread_fai.get(), thread_bam_ctx.get(), chrom_remap, opts,
                                        fetch_rows);
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
    stitch_chunk_haps(phasing_chunks, &opts);
    for (size_t i = 0; i < batch_size; ++i) {
        graph_chunks[i].chunk = std::move(phasing_chunks[i]);
        rescue_unphased_graph_reads(graph_chunks[i].chunk);
        apply_island_labels(graph_chunks[i].chunk);
    }
    return graph_chunks;
}

// The pangenome contig a FAI contig is queried under (its suffix in a PanSN
// name), or the contig itself.
static std::string query_contig_of(const std::string& contig,
                                   const std::unordered_map<std::string, std::string>& fai_full_to_suffix) {
    auto it = fai_full_to_suffix.find(contig);
    return it != fai_full_to_suffix.end() ? it->second : contig;
}

// Reads from a tabix-indexed GAF file: each worker seeks directly to the
// overlapping region -- only the relevant reads are decompressed and parsed,
// making this efficient even for very large (100+ GB) GAF files.
static std::vector<GraphChunkBuildResult> process_graph_chunk_batch_indexed_gaf(
    const std::string& sites,
    const std::unordered_map<std::string, hts_pos_t>& max_span,
    const std::vector<RegionChunk>& chunks,
    size_t batch_begin,
    size_t batch_end,
    const bam_hdr_t* header,
    const std::string& gaf_file,
    int min_mapq,
    const std::unordered_map<std::string, std::string>& fai_full_to_suffix,
    const std::unordered_map<std::string, std::string>& chrom_remap,
    const Options& opts)
{
    const std::string batch_contig = header->target_name[chunks[batch_begin].tid];
    const std::string query_contig = query_contig_of(batch_contig, fai_full_to_suffix);
    const auto open_rows = [&]() -> ChunkRowFetcher {
        auto handle = std::make_shared<IndexedGafHandle>(gaf_file);
        return [handle, &query_contig, min_mapq](const RegionChunk& region, const GraphSiteCatalogView& view) {
            return scan_indexed_gaf_chunk(*handle, query_contig, std::max<hts_pos_t>(0, region.beg - 1),
                                          region.end, view, min_mapq);
        };
    };
    return run_graph_chunk_batch(sites, max_span, chunks, batch_begin, batch_end, batch_contig, chrom_remap, opts,
                                 open_rows);
}

void run_collect_graph_variation(const Options& opts) {
    // The GAF must be bgzip-compressed and tabix-indexed on index-gaf's
    // coordinate columns, so each chunk reads only its own region.
    require_indexed_gaf(opts.gaf_file);
    // A site whose every allele spans more graph bases than any alignment on
    // its contig cannot be traversed whole by a read: it is not loaded.
    const std::unordered_map<std::string, hts_pos_t> max_span = read_gaf_max_spans(opts.gaf_file);

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
    // Union gap phasing, when a BAM is given. One context for the whole run
    // owns the BAM handles each worker's alignment solve reads.
    std::unique_ptr<WorkerContext> bam_ctx;
    if (!opts.bam_files.empty()) {
        bam_ctx = std::make_unique<WorkerContext>(opts);
        std::cerr << "Union gap phasing with alignments from " << opts.primary_bam_file() << "\n";
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
    // assigned by init_assign_read_hap_based_on_cons_alle with no minimum-observation or margin
    // requirement, so a single informative site is enough to commit it.
    struct PhaseReadDiag {
        int hap = 0;
        hts_pos_t phase_set = kUnphasedReadPhaseSet;
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

        std::vector<GraphChunkBuildResult> graph_chunks = process_graph_chunk_batch_indexed_gaf(
            opts.graph_sites, max_span, chunks, batch_begin, batch_end, header.get(), opts.gaf_file,
            opts.min_mapq, fai_full_to_suffix, chrom_remap, opts);


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
                                          : kUnphasedReadPhaseSet;
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
                merge_graph_chunk_into_read_rows(phased_bam_rows, gc);
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
        << "Phase long reads from their pangenome graph alignments (GAF); with --bam,\n"
        << "also from their linear alignments (union gap phasing).\n"
        << "\n"
        << "Required:\n"
        << "      --ref FILE                Reference FASTA (indexed)\n"
        << "      --sites PATH              Site catalogs from build-snarl-catalog: its output folder\n"
        << "                                (<contig>.pgsites), one .pgsites file, or a bgzipped,\n"
        << "                                tabix-indexed sites VCF\n"
        << "      --gaf FILE                The reads' GAF indexed by reference position, from\n"
        << "                                pgphase index-gaf (bgzipped + tabix-indexed)\n"
        << "\n"
        << "Input:\n"
        << "      --bam FILE                Indexed BAM of the same reads: alignment-only sites and\n"
        << "                                reads join each chunk's graph solve (recommended;\n"
        << "                                --min-mapq defaults to 1 with --bam)\n"
        << "\n"
        << "Output:\n"
        << "  -o, --output FILE             Candidate TSV [output.tsv]\n"
        << "  -v, --vcf-output FILE         Candidate VCF\n"
        << "      --phased-vcf-out FILE     Phased VCF with GT:DP:AD:VAF:GQ:PS\n"
        << "      --phased-bam-out FILE     Unaligned BAM with HP/PS tags per read\n"
        << "\n"
        << "Options:\n"
        << "  -t, --threads INT             Worker threads [all cores]\n"
        << "  -q, --min-mapq INT            Minimum read mapping quality [5]\n"
        << "  -D, --min-depth INT           Minimum total depth [5]\n"
        << "      --min-alt-depth INT       Minimum alt depth [2]\n"
        << "      --min-af FLOAT            Minimum allele fraction [0.20]\n"
        << "      --max-af FLOAT            Maximum allele fraction [0.80]\n"
        << "      --min-sv-len INT          Min SV length for SVTYPE/SVLEN tags [30]\n"
        << "      --chunk-size INT          Region chunk size [500000; 1000000 with --bam]\n"
        << "  -r, --region STR              Restrict to region (may be repeated)\n"
        << "      --region-file FILE        BED file of regions\n"
        << "      --autosome                Process chr1-22 / 1-22 only\n"
        << "      --hifi                    HiFi read mode [default]\n"
        << "      --ont                     ONT read mode (enables strand-bias filter)\n"
        << "      --strand-bias-pval FLOAT  Max p-value for ONT strand-bias filter [0.01]\n"
        << "  -V, --verbose INT             Verbosity level [0]\n"
        << "  -h, --help                    Print this help\n"
        << "\n"
        << "Diagnostics:\n"
        << "      --filtered-sites-out FILE Catalog sites that did not become candidates, and why\n"
        << "      --phase-sites-out FILE    Retained graph sites with SITE_ID\n"
        << "      --phase-reads-out FILE    Per-read phasing evidence\n"
        << "      --phase-matrix-dump PATH  Phasing inputs and incoming assignments\n"
        << "\n"
        << "Example:\n"
        << "  pgphase collect-graph-variation \\\n"
        << "      --ref ref.fa \\\n"
        << "      --sites chr20.sites.vcf.gz \\\n"
        << "      --gaf reads.coord.gaf.gz \\\n"
        << "      --bam reads.bam \\\n"
        << "      --phased-vcf-out phased.vcf \\\n"
        << "      --phased-bam-out phased.bam \\\n"
        << "      -t 16\n";
}

enum GraphCollectOption {
    kGcBam = 2000,
    kGcMinAltDepth = 1000,
    kGcMinAf,
    kGcMaxAf,
    kGcMinSvLen,
    kGcChunkSize,
    kGcPhasedVcf,
    kGcGafFile,
    kGcRegionFile,
    kGcAutosome,
    kGcOnt,
    kGcHifi,
    kGcStrandBiasPval,
    kGcPhasedBam,
    kGcRef,
    kGcSites,
    kGcFilteredSitesOut,
    kGcPhaseSitesOut,
    kGcPhaseReadsOut,
    kGcPhaseMatrixDump,
};

} // namespace

} // namespace pgphase_collect

int collect_graph_variation(int argc, char* argv[]) {
    using namespace pgphase_collect;
    Options opts;
    opts.threads = default_thread_count();
    opts.min_mapq = kDefaultGraphMinMapq;
    bool chunk_size_explicit = false;
    bool min_mapq_explicit = false;

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
        {"bam",              required_argument, nullptr, kGcBam},
        {"filtered-sites-out", required_argument, nullptr, kGcFilteredSitesOut},
        {"phase-sites-out",   required_argument, nullptr, kGcPhaseSitesOut},
        {"phase-reads-out",   required_argument, nullptr, kGcPhaseReadsOut},
        {"phase-matrix-dump", required_argument, nullptr, kGcPhaseMatrixDump},
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
        {"hifi",              no_argument,       nullptr, kGcHifi},
        {"ont",               no_argument,       nullptr, kGcOnt},
        {"strand-bias-pval",  required_argument, nullptr, kGcStrandBiasPval},
        {"ref",               required_argument, nullptr, kGcRef},
        {"sites",             required_argument, nullptr, kGcSites},
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
            // Union gap phasing: alignment-only sites and reads join each
            // chunk's graph solve.
            case kGcBam:  opts.bam_files.push_back(optarg); break;

            case kGcFilteredSitesOut: opts.output_filtered_sites = optarg; break;
            case kGcPhaseSitesOut: opts.output_phase_sites = optarg; break;
            case kGcPhaseReadsOut: opts.output_phase_reads = optarg; break;
            case kGcPhaseMatrixDump: opts.phase_matrix_dump_prefix = optarg; break;
            case 't': opts.threads = parse_int_arg(optarg, "--threads"); break;
            case 'q':
                opts.min_mapq = parse_int_arg(optarg, "--min-mapq");
                min_mapq_explicit = true;
                break;
            case 'D': opts.min_depth = parse_int_arg(optarg, "--min-depth"); break;
            case kGcMinAltDepth:  opts.min_alt_depth = parse_int_arg(optarg, "--min-alt-depth"); break;
            case kGcMinAf:        opts.min_af = parse_double_arg(optarg, "--min-af"); break;
            case kGcMaxAf:        opts.max_af = parse_double_arg(optarg, "--max-af"); break;
            case kGcMinSvLen:     opts.min_sv_len = parse_int_arg(optarg, "--min-sv-len"); break;
            case kGcChunkSize:
                opts.chunk_size = parse_ll_arg(optarg, "--chunk-size");
                chunk_size_explicit = true;
                break;
            case 'r': opts.regions.push_back(optarg); break;
            case kGcRegionFile:   opts.region_file = optarg; break;
            case kGcAutosome:     opts.autosome = true; break;
            case kGcGafFile:      opts.gaf_file = optarg; break;
            case kGcHifi:         opts.read_technology = ReadTechnology::Hifi; break;
            case kGcOnt:          opts.read_technology = ReadTechnology::Ont; break;
            case kGcStrandBiasPval: opts.strand_bias_pval = parse_double_arg(optarg, "--strand-bias-pval"); break;
            case kGcRef:          opts.ref_fasta = optarg; break;
            case kGcSites:        opts.graph_sites = optarg; break;
            case 'V': opts.verbose = parse_int_arg(optarg, "--verbose"); break;
            case 'h': print_graph_collect_help(); return 0;
            default:  print_graph_collect_help(); return 1;
        }
    }

    if (!chunk_size_explicit && !opts.primary_bam_file().empty())
        opts.chunk_size = kDefaultGraphBamChunkSize;

    if (opts.ref_fasta.empty() || opts.graph_sites.empty()) {
        std::cerr << "Error: --ref and --sites are required\n";
        print_graph_collect_help();
        return 1;
    }

    require_site_catalog(opts.graph_sites);

    if (opts.gaf_file.empty()) {
        std::cerr << "Error: --gaf is required\n";
        print_graph_collect_help();
        return 1;
    }

    // With --bam, every read the EM can place is labelled; low-MAPQ reads are
    // labelled but do not shape the model.
    if (!opts.bam_files.empty() && !min_mapq_explicit) opts.min_mapq = 1;
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
