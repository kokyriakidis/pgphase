// The graph command's records: each chunk's rows as VCF records for the shared writer (graph_records.hpp).

#include "graph_records.hpp"

#include "collect_output.hpp"
#include "collect_phase.hpp"
#include "collect_var.hpp"
#include "noise_filter.hpp"
#include "trace.hpp"

#include <algorithm>
#include <array>
#include <cctype>
#include <iostream>
#include <map>
#include <set>
#include <string>
#include <tuple>
#include <unordered_map>
#include <unordered_set>
#include <vector>

namespace pgphase_collect {

CandidateTable graph_chunks_to_candidate_table(
    const std::vector<GraphChunkBuildResult>& graph_chunks,
    const std::unordered_map<std::string, int>& contig_to_tid,
    const Options& opts,
    ReferenceCache& reference,
    const bam_hdr_t* header)
{
    CandidateTable result;

    for (const GraphChunkBuildResult& graph_chunk : graph_chunks) {
        const PhasingChunk& chunk = graph_chunk.chunk;
        const pgphase_trace::Scope trace_scope(chunk.region.chunk_id);
        std::unordered_map<size_t, std::string> rephase_votes;  // trace: the vote behind each placement
        const std::unordered_map<size_t, UnlinkedHetPlacement> rephase =
            place_unlinked_hets(chunk, pgphase_trace::on() ? &rephase_votes : nullptr);
        for (size_t ci = 0; ci < chunk.candidates.size(); ++ci) {
            CandidateVariant rephased;
            const CandidateVariant* mp = &chunk.candidates[ci];
            if (const auto it = rephase.find(ci); it != rephase.end()) {
                rephased = chunk.candidates[ci];
                rephased.phase_set = it->second.phase_set;
                rephased.hap_to_cons_alle[1] = it->second.hap1_allele;
                rephased.hap_to_cons_alle[2] = it->second.hap2_allele;
                mp = &rephased;
            }
            const CandidateVariant& mcand = *mp;

            if (ci >= graph_chunk.site_meta.size()) continue;
            const GraphSiteMeta& meta = graph_chunk.site_meta[ci];
            // A row the rebuild writes nothing for, by reason (stage trace).
            const auto skip = [&](std::string_view reason, std::string detail = {}) {
                if (!pgphase_trace::on()) return;
                const hts_pos_t beg = meta.pos > 0 ? meta.pos : mcand.key.pos;
                std::string key = meta.ref.empty() ? pgphase_trace::key_text(mcand.key) : meta.ref + ">";
                for (size_t i = 0; !meta.ref.empty() && i < meta.alts.size(); ++i)
                    key += (i ? "," : "") + meta.alts[i].substr(0, 60);
                pgphase_trace::emit("output", "skip", beg, beg + std::max<hts_pos_t>(1, meta.ref.size()) - 1, key,
                                    reason, detail + (detail.empty() ? "" : " ") + pgphase_trace::genotype_text(mcand) +
                                    (mcand.bam_injected ? " source=injected" : " source=catalog"));
            };
            if (rephase.count(ci) != 0)
                pgphase_trace::candidate("output", "rephase", chunk.candidates[ci], "lone_het_rephased",
                                         pgphase_trace::genotype_text(chunk.candidates[ci]) + " " +
                                         pgphase_trace::genotype_text(mcand, "to_") +
                                         (rephase_votes.count(ci) ? " " + rephase_votes[ci] : ""));
            if (meta.bam_alt_deletion_no_ref || meta.no_record) {
                skip(meta.no_record ? (meta.msa_hom ? "no_record" : "no_record") : "bam_alt_deletion_no_ref");
                continue;
            }

            // A site injected by union gap phasing (no catalog id) is
            // written only when it carries phase information. The alignment's
            // in-gap discovery also calls homozygous variants, and emitting
            // those would add ~12,850 hom records to chr20 -- a 23% larger VCF
            // whose extra content is alignment-discovered calls appearing ONLY
            // where the alignment solve found sites, which is a biased subset of
            // the genome. The output is the catalog's sites plus what union
            // phasing phased, not a variant call set. The exception is the MSA's
            // homozygous calls inside a region whose two consensus sequences
            // differ (GraphSiteMeta::msa_hom): without them the region's
            // haplotypes are spelled wrong (chr12:17,633,561, the insertion both
            // carry beside the het T run). With the catalog rows the MSA does not
            // carry now unwritten, records over HG002 HiFi chr20 99,721 -> 98,837,
            // chr12 167,368 -> 167,866, chr18 100,853 -> 100,592.
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
                    h[2] >= n_alleles || (h[1] == h[2] && !meta.msa_hom)) {
                    skip("injected_no_het_consensus");
                    continue;
                }
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
                if (alt_idx < 0 || alt_idx >= static_cast<int>(meta.alts.size())) {
                    skip(meta.alts.empty() ? "no_meta" : "allele_index_out_of_range", "allele=" + std::to_string(new_a));
                    continue;
                }
                const std::string& raw_alt = meta.alts[static_cast<size_t>(alt_idx)];
                if (raw_alt.empty() || raw_alt == "*") {
                    skip("allele_empty_or_star", "allele=" + std::to_string(new_a));
                    continue;
                }

                if (meta.ref.empty()) {
                    skip("no_ref", "allele=" + std::to_string(new_a));
                    continue;
                }

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
                cand.vcf_unphased = mcand.vcf_unphased;
                cand.window_phased = mcand.window_phased;

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
                cand.split_allele = n_new_alleles > 2;
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
                const bool is_hom_alt = meta.msa_hom || (!het_by_consensus &&
                    ((ref_cov == 0 && alt_cov >= opts.min_alt_depth) || validated_hom_alt));
                // Recovery already classified an MSA heterozygote from its
                // own reads and haplotype consensus. Its REF/ALT depth can be
                // very uneven; reapplying graph discovery thresholds here
                // turns a valid phased source row into LOW_AF at emission.
                const bool verified_injected_het = mcand.bam_injected &&
                    mcand.counts.category == VariantCategory::NoisyCandHet &&
                    mcand.msa_verified && mcand.alignment_verified &&
                    (mcand.phase_set > 0 || mcand.vcf_unphased) && het_by_consensus;
                if (verified_injected_het && !meta.msa_hom) {
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
                // the alignment solve found sites -- save the MSA's, which spell
                // the het region they sit in (msa_hom, above).
                if (ci < graph_chunk.site_ids.size() && graph_chunk.site_ids[ci].empty() && !meta.msa_hom) {
                    const VariantCategory c = cand.counts.category;
                    if (c != VariantCategory::CleanHetSnp &&
                        c != VariantCategory::CleanHetIndel &&
                        c != VariantCategory::NoisyCandHet) {
                        pgphase_trace::candidate("output", "skip", cand, "injected_not_het_category",
                                                 pgphase_trace::genotype_text(cand));
                        continue;
                    }
                }

                // An equal-length multi-base substitution (an MNP: a catalog walk
                // changing several adjacent bases) is written as its single-base
                // changes, each with the MNP's genotype and phase set -- all its
                // bases change on one haplotype. Written whole (GA>CC), the truth
                // and other callers' single-base records never matched it: over
                // chr12/18/20 +3,792 phased truth hets, reads unchanged. A copy
                // of one of those SNPs from another row collapses below.
                if (cand.key.type == VariantType::Snp && cand.key.ref_len > 1 &&
                    static_cast<int>(cand.key.alt.size()) == cand.key.ref_len &&
                    static_cast<int>(ref_seq.size()) == cand.key.ref_len) {
                    for (int i = 0; i < cand.key.ref_len; ++i) {
                        const auto up = [](char c) { return static_cast<char>(std::toupper(static_cast<unsigned char>(c))); };
                        if (up(ref_seq[static_cast<size_t>(i)]) == up(cand.key.alt[static_cast<size_t>(i)])) continue;
                        CandidateVariant snp = cand;
                        snp.key.pos = cand.key.pos + i;
                        snp.key.alt = std::string(1, cand.key.alt[static_cast<size_t>(i)]);
                        snp.key.ref_len = 1;
                        pgphase_trace::candidate("output", "record", snp, "mnp_split",
                                                 pgphase_trace::genotype_text(snp) + (meta.msa_hom ? " msa_hom" : ""));
                        result.push_back(std::move(snp));
                    }
                    continue;
                }
                pgphase_trace::candidate("output", "record", cand, cand.bam_injected ? "injected" : "catalog",
                                         pgphase_trace::genotype_text(cand) + (meta.msa_hom ? " msa_hom" : ""));
                result.push_back(std::move(cand));
            }
        }
    }

    // Indels are left-aligned (the VCF convention; truth sets and
    // alignment-based callers write them so), before duplicates collapse, so
    // two descriptions of one indel become one record. Catalog alleles sit
    // wherever their snarl put them: written so, over chr12/18/20 919 phased
    // truth hets did not match their own truth record; reads unchanged.
    if (header != nullptr) {
        const auto base = [&](int tid, hts_pos_t p) {
            return static_cast<char>(std::toupper(static_cast<unsigned char>(reference.base(tid, p, header))));
        };
        for (CandidateVariant& c : result) {
            const int tid = c.key.tid;
            if (c.key.type == VariantType::Insertion && c.key.ref_len == 0 && !c.key.alt.empty()) {
                std::string alt = c.key.alt;
                for (char& ch : alt) ch = static_cast<char>(std::toupper(static_cast<unsigned char>(ch)));
                hts_pos_t pos = c.key.pos;
                for (int n = 0; n < 10000 && pos > 2; ++n) {
                    const char prev = base(tid, pos - 1);
                    if (prev == 'N' || prev != alt.back()) break;
                    alt = prev + alt.substr(0, alt.size() - 1);
                    --pos;
                }
                c.key.pos = pos;
                c.key.alt = alt;
            } else if (c.key.type == VariantType::Deletion && c.key.ref_len > 0 && c.key.alt.empty()) {
                hts_pos_t pos = c.key.pos;
                const hts_pos_t len = c.key.ref_len;
                for (int n = 0; n < 10000 && pos > 2; ++n) {
                    const char prev = base(tid, pos - 1);
                    if (prev == 'N' || prev != base(tid, pos + len - 1)) break;
                    --pos;
                }
                c.key.pos = pos;
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
    // phase whenever the graph row had a few more observations. Next, a het
    // written unphased (the union EM's MSA-genotyped het its reads cannot
    // phase): a repeat call of the same deletion with more reads would
    // otherwise drop it from the VCF (chr12:5,924,102, 86 truth hets on
    // chr12 alone). Then any copy with an ALT genotype over one with none,
    // which writes no record: an uncalled catalog repeat row with more reads
    // dropped the MSA's 1/1 insertion at chr20:9,664,305.
    const auto het_rank = [](const CandidateVariant& candidate) {
        const auto category = candidate.counts.category;
        const auto& h = candidate.hap_to_cons_alle;
        const bool het = h[1] >= 0 && h[2] >= 0 && h[1] != h[2] &&
            (category == VariantCategory::CleanHetSnp ||
             category == VariantCategory::CleanHetIndel ||
             (category == VariantCategory::NoisyCandHet &&
              candidate.msa_verified && candidate.alignment_verified));
        if (het) return candidate.phase_set > 0 ? 3 : candidate.vcf_unphased ? 2 : 1;
        const bool written_category = category == VariantCategory::CleanHetSnp ||
            category == VariantCategory::CleanHetIndel || category == VariantCategory::CleanHom ||
            category == VariantCategory::NoisyCandHet || category == VariantCategory::NoisyCandHom;
        return written_category && (h[1] > 0 || h[2] > 0) ? 1 : 0;
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
                const int write_rank = het_rank(result[write]);
                const int read_rank = het_rank(result[read]);
                const bool replace = read_rank > write_rank ||
                    (read_rank == write_rank && result[read].counts.total_cov > result[write].counts.total_cov);
                const CandidateVariant& loser = replace ? result[write] : result[read];
                const CandidateVariant& winner = replace ? result[read] : result[write];
                if (pgphase_trace::on())
                    pgphase_trace::candidate("output", "drop", loser, "duplicate_collapsed",
                                             pgphase_trace::genotype_text(loser) + " rank=" +
                                             std::to_string(het_rank(loser)) + " " +
                                             pgphase_trace::genotype_text(winner, "kept_") + " kept_rank=" +
                                             std::to_string(het_rank(winner)));
                if (replace) result[write] = std::move(result[read]);
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

    merge_adjacent_deletion_insertion(result);
    make_colocated_alleles_complementary(result, opts.min_alt_depth);
    drop_conflicting_haplotype_alleles(result);
    drop_contradicted_walks(result, [&reference, header](int tid, hts_pos_t pos) {
        return static_cast<char>(std::toupper(static_cast<unsigned char>(reference.base(tid, pos, header))));
    });

    return result;
}

}  // namespace pgphase_collect
