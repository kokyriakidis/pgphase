/**
 * @file collect_phase.cpp
 * k-means read-haplotype clustering for diploid phasing.
 */

#include "collect_phase.hpp"

#include <stdexcept>
#include "collect_phase_pgbam.hpp"

#include <algorithm>
#include <array>
#include <cinttypes>
#include <climits>
#include <cmath>
#include <cstdio>
#include <functional>
#include <map>
#include <optional>
#include <set>
#include <string>
#include <vector>

extern "C" {
#include "cgranges.h"
}

namespace pgphase_collect {

// Read visit count / index (one slot per read).
static inline int read_visit_count(const PhasingChunk& chunk) {
    return chunk.ordered_read_ids.empty() ? static_cast<int>(chunk.reads.size())
                                         : static_cast<int>(chunk.ordered_read_ids.size());
}

static inline int read_at_visit_ord(const PhasingChunk& chunk, int ord) {
    return chunk.ordered_read_ids.empty() ? ord : chunk.ordered_read_ids[static_cast<size_t>(ord)];
}

static bool parse_debug_site_pos(const std::string& site, hts_pos_t& pos_out) {
    const size_t colon = site.find(':');
    if (colon == std::string::npos) return false;
    std::string pos_s = site.substr(colon + 1);
    pos_s.erase(std::remove(pos_s.begin(), pos_s.end(), ','), pos_s.end());
    if (pos_s.empty()) return false;
    try {
        const long long p = std::stoll(pos_s);
        if (p < 1) return false;
        pos_out = static_cast<hts_pos_t>(p);
        return true;
    } catch (...) {
        return false;
    }
}

// ════════════════════════════════════════════════════════════════════════════
// Category flag mapping
// ════════════════════════════════════════════════════════════════════════════

uint32_t category_to_flag(VariantCategory c) {
    // Map VariantCategory enum to its lcd_var_i_to_cate bitmask flag.
    switch (c) {
        case VariantCategory::LowCoverage:
            return kLongcalldLowCovVar;
        case VariantCategory::LowAlleleFraction:
            return kLongcalldLowAfVar;
        case VariantCategory::StrandBias:
            return kLongcalldStrandBiasVar;
        case VariantCategory::CleanHetSnp:
            return kCandCleanHetSnp;
        case VariantCategory::CleanHetIndel:
            return kCandCleanHetIndel;
        case VariantCategory::CleanHom:
            return kCandCleanHom;
        case VariantCategory::NoisyCandHet:
            return kCandNoisyCandHet;
        case VariantCategory::NoisyCandHom:
            return kCandNoisyCandHom;
        case VariantCategory::NoisyResolved:
            return 0; // unmapped category
        case VariantCategory::RepeatHetIndel:
            return kLongcalldRepHetVar;
        case VariantCategory::NonVariant:
            return kLongcalldNonVar;
    }
    return 0;
}

// ════════════════════════════════════════════════════════════════════════════
// Helper functions for k-means phasing.
// ════════════════════════════════════════════════════════════════════════════

// Read_init_hap_phase_set.
static void read_init_hap_phase_set(PhasingChunk& chunk) {
    for (size_t i = 0; i < chunk.reads.size(); ++i) {
        chunk.haps[i] = 0;
        chunk.phase_sets[i] = kUnphasedReadPhaseSet;
    }
}

// Width of hap_to_alle_profile[1/2] (n_uniq_alles, or 2 if alle_covs is empty).
static int variant_allele_slots(const CandidateVariant& var) {
    const VariantCounts& c = var.counts;
    if (c.n_uniq_alles > 0) return c.n_uniq_alles;
    if (!c.alle_covs.empty()) return static_cast<int>(c.alle_covs.size());
    return 2;
}

// Returns majority allele index or -1 if inconclusive.
// Strict > on coverage so ties keep the lower index (ref-first).
static int get_var_init_max_cov_allele(bool is_ont, const CandidateVariant& var) {
    if (is_ont && var.is_homopolymer_indel) return -1;
    const VariantCounts& c = var.counts;
    int max_cov = 0;
    int max_cov_alle_i = -1;
    if (!c.alle_covs.empty()) {
        const int n = c.n_uniq_alles > 0
                          ? c.n_uniq_alles
                          : static_cast<int>(c.alle_covs.size());
        for (int i = 0; i < n; ++i) {
            if (c.alle_covs[i] > max_cov) {
                max_cov = c.alle_covs[i];
                max_cov_alle_i = i;
            }
        }
        return max_cov_alle_i;
    }
    const int n = (c.n_uniq_alles > 0) ? c.n_uniq_alles : 2;
    for (int i = 0; i < n; ++i) {
        const int cov = (i == 0) ? c.ref_cov : (i == 1 ? c.alt_cov : 0);
        if (cov > max_cov) {
            max_cov = cov;
            max_cov_alle_i = i;
        }
    }
    return max_cov_alle_i;
}

// Initialize hap_to_alle_profile (zeroed) and hap_to_cons_alle for all valid vars.
// Hom vars get cons_alle[1]=cons_alle[2]=1; het vars get -1/-1.
// Var_init_hap_profile_cons_allele.
void var_init_hap_profile_cons_allele(bool is_ont,
                                              CandidateTable& variants,
                                              const std::vector<int>& valid_var_idx,
                                              bool preserve_decided) {
    for (int vi : valid_var_idx) {
        CandidateVariant& var = variants[vi];
        const int na = variant_allele_slots(var);
        const bool first_init = var.hap_to_alle_profile[1].empty() && var.hap_to_alle_profile[2].empty();
        if (first_init) {
            var.hap_to_alle_profile[0].assign(na, 0);
        }
        var.hap_to_alle_profile[1].assign(na, 0);
        var.hap_to_alle_profile[2].assign(na, 0);
        var.hap_to_cons_alle[0] = get_var_init_max_cov_allele(is_ont, var);
        // Anchored: a site the previous round already decided keeps its
        // consensus. The per-iteration vote tallies above are still zeroed --
        // those are scratch -- but hap_to_cons_alle IS the gauge, and clearing it
        // is what let the second round re-decide the first round's parity.
        if (preserve_decided &&
            (var.hap_to_cons_alle[1] != -1 || var.hap_to_cons_alle[2] != -1))
            continue;
        const uint32_t vic = var.lcd_var_i_to_cate;
        if (vic == kCandNoisyCandHom || vic == kCandCleanHom) {
            var.hap_to_cons_alle[1] = 1;
            var.hap_to_cons_alle[2] = 1;
        } else {
            var.hap_to_cons_alle[1] = -1;
            var.hap_to_cons_alle[2] = -1;
        }
    }
}

// Zero all hap allele profiles (0, 1, 2) for the valid vars; preserves cons_alle.
// Initialize per-haplotype allele count profiles to zero.
static void var_init_hap_to_alle_profile(CandidateTable& variants,
                                          const std::vector<int>& valid_var_idx) {
    for (int vi : valid_var_idx) {
        CandidateVariant& v = variants[vi];
        const int na = variant_allele_slots(v);
        v.hap_to_alle_profile[0].assign(na, 0);
        v.hap_to_alle_profile[1].assign(na, 0);
        v.hap_to_alle_profile[2].assign(na, 0);
    }
}

// Pick pivot: deepest CleanHetSnp > CleanHetIndel > NoisyHetSnp > NoisyHetIndel.
// Returns index into valid_var_idx, or -1.
// Select_init_var.
static int select_init_var(const CandidateTable& variants,
                            const std::vector<int>& valid_var_idx) {
    int snp_i = -1, indel_i = -1, noisy_snp_i = -1, noisy_indel_i = -1;
    int snp_dp = 0, indel_dp = 0, noisy_snp_dp = 0, noisy_indel_dp = 0;
    for (int _vi = 0; _vi < (int)valid_var_idx.size(); ++_vi) {
        const CandidateVariant& var = variants[valid_var_idx[_vi]];
        const uint32_t vic = var.lcd_var_i_to_cate;
        if (vic == kCandCleanHetSnp) {
            if (snp_i == -1 || var.counts.total_cov > snp_dp) {
                snp_i = _vi; snp_dp = var.counts.total_cov;
            }
        } else if (vic == kCandCleanHetIndel) {
            if (indel_i == -1 || var.counts.total_cov > indel_dp) {
                indel_i = _vi; indel_dp = var.counts.total_cov;
            }
        } else if (vic == kCandNoisyCandHet) {
            if (var.key.type == VariantType::Snp) {
                if (noisy_snp_i == -1 || var.counts.total_cov > noisy_snp_dp) {
                    noisy_snp_i = _vi; noisy_snp_dp = var.counts.total_cov;
                }
            } else if (!var.is_homopolymer_indel) {
                if (noisy_indel_i == -1 || var.counts.total_cov > noisy_indel_dp) {
                    noisy_indel_i = _vi; noisy_indel_dp = var.counts.total_cov;
                }
            }
        }
    }
    if (snp_i != -1) return snp_i;
    if (indel_i != -1) return indel_i;
    if (noisy_snp_i != -1) return noisy_snp_i;
    return noisy_indel_i;
}

// Update hap_to_cons_alle[hap] via argmax of the allele profile, with ONT 67% guard.
// Update_var_hap_to_cons_alle.
void update_var_hap_to_cons_alle(bool is_ont, CandidateVariant& var, int hap) {
    if (hap == 0) return;
    const auto& prof = var.hap_to_alle_profile[hap];
    int max_cov = 0, max_alle = -1, total = 0;
    for (size_t a = 0; a < prof.size(); ++a) {
        total += prof[a];
        if (prof[a] > max_cov) {
            max_cov = prof[a];
            max_alle = static_cast<int>(a);
        }
    }
    if (is_ont && var.is_homopolymer_indel && max_cov < total * 0.67) max_alle = -1;
    var.hap_to_cons_alle[hap] = max_alle;
}

// Score a read allele against the hap consensus; infers complement when one hap is unknown.
// Returns +var_score (match), −var_score (mismatch), or 0 (no info).
// A record with exactly two ALTs and reads behind each describes a locus whose
// haplotypes both differ from the reference. This decides only how such a record
// is ORIENTED once it is already being phased -- not whether it may enter the
// solve. Letting these records in generally was measured and rejected: with the
// noisy class in the solve it took the panel from 1 concordant-to-discordant
// read to 31 and one window from 100.00% to 90.94%, because a co-located
// deletion pair in a repeat tract is usually the uninformative class. Orienting
// the ones already in the solve is a different matter: independent haplotype
// majorities can select the same allele twice, which emitted four of the panel's
// twenty-seven multiallelic loci as 1|1 or 2|2 with reads on both alleles.
static bool two_allele_het(const CandidateVariant& var) {
    return var.msa_insertion_alts.size() == 2 && var.counts.alle_covs.size() == 3 &&
           var.counts.alle_covs[1] > 0 && var.counts.alle_covs[2] > 0;
}

// Side effect: may set hap_to_cons_alle[hap] or [3-hap] when one is -1.
// Score a read against consensus alleles: +1 for agreement, -1 for conflict.
int read_to_cons_allele_score(CandidateVariant& var, int hap, int allele_i,
                                     bool msa_sites_vote_without_gap_link,
                                     bool infer_complement_at_multiallelic,
                                     bool upstream_read_scoring) {
    const uint32_t var_i_to_cate = var.lcd_var_i_to_cate;
    // longcallD has no equivalent gate (assign_hap.c:127-147): every candidate
    // in the mask votes. See Options::msa_sites_vote_without_gap_link.
    if (!msa_sites_vote_without_gap_link &&
        !var.msa_insertion_alts.empty() && !var.gap_link_supported) return 0;
    int var_score = 1;
    const bool two_alle_ok = upstream_read_scoring || var.counts.n_uniq_alles <= 2;
    if (two_alle_ok && var_i_to_cate == kCandCleanHetSnp) var_score = 2;
    else if (two_alle_ok && var_i_to_cate == kCandCleanHetIndel) var_score = 2;
    if (var.hap_to_cons_alle[hap] == -1 && var.hap_to_cons_alle[3 - hap] == -1) return 0;
    // A multiallelic site has no unique complementary allele. Infer it only
    // for a biallelic site; otherwise let actual read observations resolve it.
    if (infer_complement_at_multiallelic || variant_allele_slots(var) == 2) {
        if (var.hap_to_cons_alle[hap] == -1) var.hap_to_cons_alle[hap] = 1 - var.hap_to_cons_alle[3 - hap];
        if (var.hap_to_cons_alle[3 - hap] == -1) var.hap_to_cons_alle[3 - hap] = 1 - var.hap_to_cons_alle[hap];
    }
    // A non-anchor het still gets a consensus and a phase set -- it is a real
    // call -- but it must not influence which haplotype a read is assigned to,
    // which is the decision --anchor-af-margin exists to protect.  The return
    // must come *after* the consensus fill-in above, or the site never acquires
    // a genotype and is dropped at output instead of merely not voting.
    if (var_i_to_cate == kCandNonAnchorHet) return 0;
    if (var.hap_to_cons_alle[hap] == allele_i) return var_score;
    if (var.hap_to_cons_alle[hap] == -1) return 0;
    return -var_score;
}

// Weight a candidate gets in k-means scoring (mirrors read_to_cons_allele_score):
// clean het SNP/indel count double, everything else single.
static int phase_matrix_var_weight(const CandidateVariant& var) {
    const uint32_t var_i_to_cate = var.lcd_var_i_to_cate;
    if (var_i_to_cate == kCandNonAnchorHet) return 0;
    if (var.counts.n_uniq_alles <= 2 && var_i_to_cate == kCandCleanHetSnp) return 2;
    if (var.counts.n_uniq_alles <= 2 && var_i_to_cate == kCandCleanHetIndel) return 2;
    return 1;
}

// Dump the per-read x per-variant allele matrix that k-means consumes, in long
// format (one row per read/variant observation). Debug-only; gated on
// opts.phase_matrix_dump_prefix. Must be called BEFORE any phasing mutates
// hap_to_cons_alle. Truth labels are NOT emitted (the pipeline has none) — join
// by read qname offline. Variants are emitted in candidate-table order so the
// offline reader can reconstruct read ordering along the contig.
static void dump_phase_matrix(const PhasingChunk& chunk,
                              const std::vector<int>& valid_var_idx,
                              const std::vector<bool>& var_is_valid,
                              const Options& opts, uint32_t flags,
                              const char* label) {
    if (opts.phase_matrix_dump_prefix.empty()) return;

    std::string path = opts.phase_matrix_dump_prefix + ".chunk" +
                       std::to_string(chunk.region.chunk_id);
    if (label != nullptr)
        path += "." + std::string(label);
    else
        path += ".flags" + std::to_string(flags);
    path += ".tsv";
    std::FILE* fp = std::fopen(path.c_str(), "w");
    if (fp == nullptr) {
        std::fprintf(stderr, "[dump-phase-matrix] cannot open %s\n", path.c_str());
        return;
    }

    // Variant header block: idx, pos, type, category, weight. Variant idx is the
    // position within valid_var_idx (0-based, contig order).
    std::fprintf(fp, "#tid\t%d\tchunk\t%d\tflags\t%u\tn_vars\t%zu\n",
                 chunk.region.tid, chunk.region.chunk_id, flags, valid_var_idx.size());
    std::fprintf(fp, "#VAR\tvar_idx\tpos\ttype\tcate\tweight"
                     "\tref_cov\talt_cov\tallele_fraction\n");
    std::vector<int> global_to_vidx(chunk.candidates.size(), -1);
    for (int vidx = 0; vidx < (int)valid_var_idx.size(); ++vidx) {
        const int gi = valid_var_idx[vidx];
        global_to_vidx[gi] = vidx;
        const CandidateVariant& var = chunk.candidates[gi];
        const char t = (var.key.type == VariantType::Snp ? 'X'
                     : (var.key.type == VariantType::Insertion ? 'I' : 'D'));
        std::fprintf(fp,
                     "VAR\t%d\t%" PRId64 "\t%c\t%u\t%d\t%d\t%s\t%" PRId64
                     "\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%.17g\n",
                     vidx, static_cast<int64_t>(var.key.pos), t,
                     var.lcd_var_i_to_cate, phase_matrix_var_weight(var),
                     var.key.ref_len, var.key.alt.c_str(),
                     static_cast<int64_t>(var.phase_set),
                     var.hap_to_cons_alle[1], var.hap_to_cons_alle[2],
                     var.bam_injected, var.alignment_verified,
                     var.gap_link_supported, var.counts.ref_cov,
                     var.counts.alt_cov, var.counts.allele_fraction);
        std::fprintf(fp, "#META\t%d\tmsa_verified=%d\thomopolymer=%d\tgap_link_supported=%d\n",
                     vidx, var.msa_verified, var.is_homopolymer_indel, var.gap_link_supported);
    }

    // Preserve the incoming labels as well as the allele matrix. For recovery
    // rounds these are the two neighboring blocks that the next solve is trying
    // to connect; recording them makes the boundary independently replayable.
    std::fprintf(fp,
                 "#READ\tqname\tbeg\tend\tmapq\tskipped\thap\tphase_set"
                 "\tagree\tconflict\n");
    for (size_t read_i = 0; read_i < chunk.reads.size(); ++read_i) {
        const ReadRecord& read = chunk.reads[read_i];
        const int hap = read_i < chunk.haps.size() ? chunk.haps[read_i] : 0;
        const hts_pos_t phase_set = read_i < chunk.phase_sets.size()
                                          ? chunk.phase_sets[read_i]
                                          : kUnphasedReadPhaseSet;
        std::fprintf(fp,
                     "READ\t%s\t%" PRId64 "\t%" PRId64
                     "\t%d\t%d\t%d\t%" PRId64 "\t%d\t%d\n",
                     read.qname.c_str(), static_cast<int64_t>(read.beg),
                     static_cast<int64_t>(read.end), read.mapq, read.is_skipped,
                     hap, static_cast<int64_t>(phase_set),
                     read.n_clean_agree_snps,
                     read.n_clean_conflict_snps);
    }

    // Observation rows: qname, var_idx, allele (0=ref,1=alt,-1=non-inf,-2=lowqual).
    std::fprintf(fp, "#OBS\tqname\tvar_idx\tallele\n");
    for (size_t read_i = 0; read_i < chunk.reads.size(); ++read_i) {
        const ReadRecord& read = chunk.reads[read_i];
        if (read.is_skipped) continue;
        const ReadVariantProfile& prof = chunk.read_var_profile[read_i];
        if (prof.start_var_idx < 0) continue;
        for (int vi = prof.start_var_idx; vi <= prof.end_var_idx; ++vi) {
            if (!var_is_valid[vi]) continue;
            const int allele = prof.alleles[vi - prof.start_var_idx];
            std::fprintf(fp, "OBS\t%s\t%d\t%d\n",
                         read.qname.c_str(), global_to_vidx[vi], allele);
        }
    }
    std::fclose(fp);
}

void dump_recovery_phase_state(const PhasingChunk& chunk, const Options& opts,
                               const char* label) {
    std::vector<int> candidate_indices(chunk.candidates.size());
    for (size_t i = 0; i < candidate_indices.size(); ++i)
        candidate_indices[i] = static_cast<int>(i);
    const std::vector<bool> included(chunk.candidates.size(), true);
    dump_phase_matrix(chunk, candidate_indices, included, opts, 0, label);
}

// Assign a read to hap 1, 2, 0 (tied), or -1 (no informative variants).
// Updates n_clean_agree_snps / n_clean_conflict_snps on the read (for max_hap only).
// CleanHom variants contribute to agree/conflict stats but not to hap_scores.
// Init_assign_read_hap_based_on_cons_alle.
int init_assign_read_hap_based_on_cons_alle(PhasingChunk& chunk, int read_i, uint32_t flags,
                                std::optional<hts_pos_t> phase_set, bool msa_sites_vote,
                                bool infer_complement, bool upstream_read_scoring) {
    ReadRecord& read = chunk.reads[read_i];
    read.n_clean_agree_snps = 0;
    read.n_clean_conflict_snps = 0;
    read.n_bridge_agree_snps = 0;
    read.n_bridge_conflict_snps = 0;
    read.n_hp_gap_agree = 0;
    read.n_hp_gap_conflict = 0;

    const ReadVariantProfile& prof = chunk.read_var_profile[read_i];
    if (prof.start_var_idx < 0) return -1;

    int hap_scores[3] = {0, 0, 0};
    int n_vars_used[3] = {0, 0, 0};
    int n_clean_agree[3] = {0, 0, 0};
    int n_clean_conflict[3] = {0, 0, 0};
    int n_bridge_agree[3] = {0, 0, 0};
    int n_bridge_conflict[3] = {0, 0, 0};
    int n_hp_gap_agree[3] = {0, 0, 0};
    int n_hp_gap_conflict[3] = {0, 0, 0};

    for (int vi = prof.start_var_idx; vi <= prof.end_var_idx; ++vi) {
        CandidateVariant& var = chunk.candidates[vi];
        const uint32_t vic = var.lcd_var_i_to_cate;
        if ((vic & flags) == 0 || (phase_set && var.phase_set != *phase_set)) continue;
        // Homopolymer indels and noisy homozygous sites do not contribute read
        // scores. The one exception is a site the homopolymer tier admitted as a
        // gap's last-resort evidence: excluding it there left the reads inside
        // the gap with no interior evidence at all, so the tier could earn link
        // support and still not phase the reads it was reached for. Measured on
        // chr20:36,247,421-36,268,291, the admitted homopolymer deletion
        // segregates at 0.923 against read truth while the non-homopolymer
        // verified indels beside it sit at 0.509-0.644.
        if ((var.is_homopolymer_indel && (upstream_read_scoring || !var.hp_gap_scorable)) ||
            vic == kCandNoisyCandHom) continue;

        const int aidx = prof.alleles[vi - prof.start_var_idx];
        if (aidx < 0) continue;

        for (int hap = 1; hap <= 2; ++hap) {
            const int score = read_to_cons_allele_score(var, hap, aidx, msa_sites_vote, infer_complement,
                                                      upstream_read_scoring);
            if (score != 0) {
                if (vic != kCandCleanHom) n_vars_used[hap]++;
                const bool clean_snp_upstream =
                    (vic & kCandGermlineClean) != 0 && var.key.type == VariantType::Snp;
                if (upstream_read_scoring ? clean_snp_upstream
                                          : (vic == kCandCleanHetSnp &&
                                             var.counts.n_uniq_alles <= 2)) {
                    if (score > 0) n_clean_agree[hap]++;
                    else n_clean_conflict[hap]++;
                } else if (vic == kCandNoisyCandHet && var.key.type == VariantType::Snp &&
                          var.counts.n_uniq_alles <= 2) {
                    if (score > 0) n_bridge_agree[hap]++;
                    else n_bridge_conflict[hap]++;
                } else if (var.hp_gap_scorable) {
                    if (score > 0) n_hp_gap_agree[hap]++;
                    else n_hp_gap_conflict[hap]++;
                }
            }
            if (vic != kCandCleanHom) hap_scores[hap] += score;
        }
    }

    int max_hap = 0, max_score = 0, min_hap = 0, min_score = 0;
    for (int hap = 1; hap <= 2; ++hap) {
        if (hap_scores[hap] > max_score) { max_hap = hap; max_score = hap_scores[hap]; }
        else if (hap_scores[hap] < min_score) { min_hap = hap; min_score = hap_scores[hap]; }
    }

    read.hap_score_margin = std::abs(hap_scores[1] - hap_scores[2]);

    if (n_vars_used[1] == 0 && n_vars_used[2] == 0) return -1;
    if (max_score == 0 && min_score == 0) return 0;
    if (max_score > 0) {
        read.n_clean_agree_snps = n_clean_agree[max_hap];
        read.n_clean_conflict_snps = n_clean_conflict[max_hap];
        read.n_bridge_agree_snps = n_bridge_agree[max_hap];
        read.n_bridge_conflict_snps = n_bridge_conflict[max_hap];
        read.n_hp_gap_agree = n_hp_gap_agree[max_hap];
        read.n_hp_gap_conflict = n_hp_gap_conflict[max_hap];
        read.n_vars_scored = n_vars_used[max_hap];
        return max_hap;
    }
    read.n_vars_scored = n_vars_used[3 - min_hap];
    return 3 - min_hap;
}

// Update allele profile + cons_alle for all vars a read covers (used in Phase 1).
// hap==0 means unassigned — both hap 1 and 2 are updated identically.
// Update_var_hap_profile_cons_alle_based_on_read_hap.
void update_var_hap_profile_cons_alle_based_on_read_hap(PhasingChunk& chunk, bool is_ont,
                                              int read_i, int hap, uint32_t flags) {
    const ReadVariantProfile& prof = chunk.read_var_profile[read_i];
    if (prof.start_var_idx < 0) return;
    for (int vi = prof.start_var_idx; vi <= prof.end_var_idx; ++vi) {
        const uint32_t vic = chunk.candidates[vi].lcd_var_i_to_cate;
        if ((vic & flags) == 0) continue;
        const int aidx = prof.alleles[vi - prof.start_var_idx];
        if (aidx < 0) continue;
        CandidateVariant& var = chunk.candidates[vi];
        if (hap == 0) {
            for (int h = 1; h <= 2; ++h) {
                var.hap_to_alle_profile[h][aidx]++;
                update_var_hap_to_cons_alle(is_ont, var, h);
            }
        } else {
            var.hap_to_alle_profile[hap][aidx]++;
            update_var_hap_to_cons_alle(is_ont, var, hap);
        }
    }
}

// Update allele profile only for all vars a read covers (used in Phase 2).
// Update_var_hap_profile_based_on_read_hap.
void update_var_hap_profile_based_on_read_hap(PhasingChunk& chunk, int read_i, int hap, uint32_t flags,
                                   std::optional<hts_pos_t> phase_set) {
    const ReadVariantProfile& prof = chunk.read_var_profile[read_i];
    if (prof.start_var_idx < 0) return;
    for (int vi = prof.start_var_idx; vi <= prof.end_var_idx; ++vi) {
        const uint32_t vic = chunk.candidates[vi].lcd_var_i_to_cate;
        if ((vic & flags) == 0 || (phase_set && chunk.candidates[vi].phase_set != *phase_set)) continue;
        const int aidx = prof.alleles[vi - prof.start_var_idx];
        if (aidx < 0) continue;
        CandidateVariant& cand = chunk.candidates[vi];
        if (hap == 0) {
            cand.hap_to_alle_profile[1][aidx]++;
            cand.hap_to_alle_profile[2][aidx]++;
        } else {
            cand.hap_to_alle_profile[hap][aidx]++;
        }
    }
}

// Returns 1=agree, 0=conflict, -1=uninformative.
// Check_agree_haps.
int check_agree_haps(const PhasingChunk& chunk, int read_i, int hap, int var1, int var2) {
    const ReadVariantProfile& prof = chunk.read_var_profile[read_i];
    if (var1 < prof.start_var_idx || var2 > prof.end_var_idx) return -1;
    if (hap == 0) return -1;
    const int a1 = prof.alleles[var1 - prof.start_var_idx];
    const int a2 = prof.alleles[var2 - prof.start_var_idx];
    if (a1 < 0 || a2 < 0) return -1;
    const CandidateVariant& v1 = chunk.candidates[var1];
    const CandidateVariant& v2 = chunk.candidates[var2];
    const bool agree   = (v1.hap_to_cons_alle[hap] == a1 && v2.hap_to_cons_alle[hap] == a2);
    const bool conflict = (v1.hap_to_cons_alle[hap] == a1 && v2.hap_to_cons_alle[3 - hap] == a2);
    if (agree) return 1;
    if (conflict) return 0;
    return -1;
}

// Link two variants by the allele pattern a read carries across them, without
// requiring the read to have been assigned a haplotype.
//
// check_agree_haps() consults chunk.haps[read_i] and so ignores every unassigned
// read.  That is a real loss: at the junctions where blocks break, most spanning
// reads are untagged precisely *because* the block broke there, so the evidence
// that would repair the break is discarded.  A read does not need a haplotype to
// say "these two variants' alleles travel together on me" -- which is the edge a
// read-backed phasing graph is built from.  Returns 1 if the read carries the
// same haplotype's consensus at both variants, 0 if opposite ones, -1 if it does
// not resolve both.
int check_agree_alleles(const PhasingChunk& chunk, int read_i, int var1, int var2) {
    const ReadVariantProfile& prof = chunk.read_var_profile[read_i];
    if (var1 < prof.start_var_idx || var2 > prof.end_var_idx) return -1;
    const int a1 = prof.alleles[var1 - prof.start_var_idx];
    const int a2 = prof.alleles[var2 - prof.start_var_idx];
    if (a1 < 0 || a2 < 0) return -1;
    const CandidateVariant& v1 = chunk.candidates[var1];
    const CandidateVariant& v2 = chunk.candidates[var2];
    // Both are het here, so hap_to_cons_alle[1] != [2] and the match is unambiguous.
    int h1 = 0, h2 = 0;
    if (v1.hap_to_cons_alle[1] == a1) h1 = 1;
    else if (v1.hap_to_cons_alle[2] == a1) h1 = 2;
    if (v2.hap_to_cons_alle[1] == a2) h2 = 1;
    else if (v2.hap_to_cons_alle[2] == a2) h2 = 2;
    if (h1 == 0 || h2 == 0) return -1;
    return (h1 == h2) ? 1 : 0;
}

static bool merge_phase_sets_in_place(PhasingChunk& chunk,
                                      hts_pos_t upstream_phase_set,
                                      hts_pos_t downstream_phase_set,
                                      bool flip) {
    if (upstream_phase_set <= 0 || downstream_phase_set <= 0 ||
        upstream_phase_set == downstream_phase_set)
        return false;

    const auto flip_hap = [](int hap) {
        return hap == 1 ? 2 : (hap == 2 ? 1 : hap);
    };
    for (CandidateVariant& candidate : chunk.candidates) {
        if (candidate.phase_set != downstream_phase_set) continue;
        if (flip) {
            std::swap(candidate.hap_to_cons_alle[1],
                      candidate.hap_to_cons_alle[2]);
            std::swap(candidate.hap_to_alle_profile[1],
                      candidate.hap_to_alle_profile[2]);
            candidate.hap_alt = flip_hap(candidate.hap_alt);
            candidate.hap_ref = flip_hap(candidate.hap_ref);
        }
        candidate.phase_set = upstream_phase_set;
    }
    for (size_t read_i = 0; read_i < chunk.reads.size(); ++read_i) {
        if (read_i >= chunk.phase_sets.size() ||
            chunk.phase_sets[read_i] != downstream_phase_set)
            continue;
        if (flip && read_i < chunk.haps.size() &&
            (chunk.haps[read_i] == 1 || chunk.haps[read_i] == 2))
            chunk.haps[read_i] = flip_hap(chunk.haps[read_i]);
        chunk.phase_sets[read_i] = upstream_phase_set;
    }
    return true;
}

// Score one read only on loci that earned an edge in the recovery chain. The
// scoring primitive is the longcallD port; this wrapper supplies the graph
// recovery's selected-site set so unrelated MSA candidates cannot cancel a
// decisive bridge allele. Returns -1 when the read observes no selected locus.
static int assign_read_hap_from_recovery_chain(PhasingChunk& chunk,
                                               size_t read_i,
                                               hts_pos_t phase_set) {
    if (read_i >= chunk.reads.size() ||
        read_i >= chunk.read_var_profile.size())
        return -1;
    const ReadVariantProfile& profile = chunk.read_var_profile[read_i];
    if (profile.start_var_idx < 0) return -1;

    int scores[3] = {0, 0, 0};
    int sites_scored = 0;
    constexpr bool kMsaSitesVote = true;
    constexpr bool kInferComplement = true;
    constexpr bool kUpstreamReadScoring = true;
    for (int candidate_i = profile.start_var_idx;
         candidate_i <= profile.end_var_idx; ++candidate_i) {
        if (candidate_i < 0 ||
            static_cast<size_t>(candidate_i) >= chunk.candidates.size())
            continue;
        const size_t offset =
            static_cast<size_t>(candidate_i - profile.start_var_idx);
        if (offset >= profile.alleles.size() || profile.alleles[offset] < 0)
            continue;
        CandidateVariant& candidate =
            chunk.candidates[static_cast<size_t>(candidate_i)];
        if (candidate.phase_set != phase_set ||
            !candidate.gap_link_supported)
            continue;

        const int hap1_score = read_to_cons_allele_score(
            candidate, 1, profile.alleles[offset], kMsaSitesVote,
            kInferComplement, kUpstreamReadScoring);
        const int hap2_score = read_to_cons_allele_score(
            candidate, 2, profile.alleles[offset], kMsaSitesVote,
            kInferComplement, kUpstreamReadScoring);
        if (hap1_score == 0 && hap2_score == 0) continue;
        scores[1] += hap1_score;
        scores[2] += hap2_score;
        ++sites_scored;
    }
    if (sites_scored == 0) return -1;

    ReadRecord& read = chunk.reads[read_i];
    read.hap_score_margin = std::abs(scores[1] - scores[2]);
    read.n_vars_scored = sites_scored;
    if (scores[1] == scores[2]) return 0;
    return scores[1] > scores[2] ? 1 : 2;
}

using DiploidLinkCounts = std::array<std::array<int, 2>, 2>;

struct SupportedAlleleEdge {
    int upstream_candidate = -1;
    int downstream_candidate = -1;
    int same = 0;
    int cross = 0;
    size_t comparisons = 0;
    bool flip = false;
};

// Recovery parity is a binary choice: keep the downstream block's gauge or
// flip it. Test the winning same/cross count against a 50:50 null so confidence
// scales with depth instead of depending on a fixed read-count difference.
constexpr double kRecoveryParityPValue = 0.01;
constexpr double kRecoveryCandidateAnchorPValue = 0.05;

static std::optional<bool> parity_flip_at_p(
        int same, int cross, size_t comparisons, double max_p_value) {
    const int total = same + cross;
    const int winner = std::max(same, cross);
    if (total <= 0 || same == cross || comparisons == 0) return std::nullopt;

    // Under the null, X~Binomial(total, 0.5). Because winner > total/2,
    // successive upper-tail terms decrease; the recurrence avoids factorials.
    const double log_term =
        std::lgamma(static_cast<double>(total + 1)) -
        std::lgamma(static_cast<double>(winner + 1)) -
        std::lgamma(static_cast<double>(total - winner + 1)) -
        static_cast<double>(total) * std::log(2.0);
    double term = std::exp(log_term);
    double tail = term;
    for (int k = winner; k < total; ++k) {
        term *= static_cast<double>(total - k) /
                static_cast<double>(k + 1);
        tail += term;
    }

    const double corrected_p = std::min(
        1.0, tail * static_cast<double>(comparisons));
    if (corrected_p > max_p_value) return std::nullopt;
    return cross > same;
}

static std::optional<bool> significant_parity_flip(
        int same, int cross, size_t comparisons = 1) {
    return parity_flip_at_p(
        same, cross, comparisons, kRecoveryParityPValue);
}

static std::optional<bool> significant_parity_flip(
        const DiploidLinkCounts& counts) {
    return significant_parity_flip(
        counts[0][0] + counts[1][1],
        counts[0][1] + counts[1][0]);
}

static std::optional<SupportedAlleleEdge> strongest_phase_set_edge(
        PhasingChunk& chunk, hts_pos_t upstream_phase_set,
        hts_pos_t downstream_phase_set, const Options& opts) {
    if (upstream_phase_set <= 0 || downstream_phase_set <= 0 ||
        upstream_phase_set == downstream_phase_set || chunk.read_var_cr == nullptr)
        return std::nullopt;

    const auto oriented_in = [](const CandidateVariant& candidate,
                                hts_pos_t phase_set) {
        return candidate.phase_set == phase_set &&
               candidate.hap_to_cons_alle[1] >= 0 &&
               candidate.hap_to_cons_alle[2] >= 0 &&
               candidate.hap_to_cons_alle[1] != candidate.hap_to_cons_alle[2];
    };

    std::vector<int> upstream;
    std::vector<int> downstream;
    const size_t window = static_cast<size_t>(std::max(1, opts.block_link_window));
    for (size_t i = 0; i < chunk.candidates.size(); ++i) {
        if (oriented_in(chunk.candidates[i], upstream_phase_set))
            upstream.push_back(static_cast<int>(i));
        if (oriented_in(chunk.candidates[i], downstream_phase_set) &&
            downstream.size() < window)
            downstream.push_back(static_cast<int>(i));
    }
    if (upstream.size() > window)
        upstream.erase(upstream.begin(), upstream.end() - window);
    if (upstream.empty() || downstream.empty()) return std::nullopt;

    int best_margin = 0;
    int best_total = 0;
    int best_agree = 0;
    int best_conflict = 0;
    int best_upstream = -1;
    int best_downstream = -1;
    bool parity_tied = false;
    size_t comparisons = 0;
    int64_t* overlaps = nullptr;
    int64_t overlap_capacity = 0;
    for (auto left = upstream.rbegin(); left != upstream.rend(); ++left) {
        for (const int right : downstream) {
            if (*left >= right) continue;
            ++comparisons;
            int agree = 0;
            int conflict = 0;
            const int64_t overlap_count = cr_overlap(
                chunk.read_var_cr.get(), "cr", *left, right + 1,
                &overlaps, &overlap_capacity);
            for (int64_t oi = 0; oi < overlap_count; ++oi) {
                const int read_i = static_cast<int>(
                    cr_label(chunk.read_var_cr.get(), overlaps[oi]));
                if (chunk.reads[read_i].is_skipped) continue;
                const int vote = check_agree_alleles(
                    chunk, read_i, *left, right);
                if (vote > 0) ++agree;
                else if (vote == 0) ++conflict;
            }
            const int margin = std::abs(agree - conflict);
            const int total = agree + conflict;
            const bool stronger = margin > best_margin ||
                                  (margin == best_margin && total > best_total);
            if (stronger) {
                best_margin = margin;
                best_total = total;
                best_agree = agree;
                best_conflict = conflict;
                best_upstream = *left;
                best_downstream = right;
                parity_tied = false;
            } else if (margin == best_margin && total == best_total && margin > 0 &&
                       (conflict > agree) != (best_conflict > best_agree)) {
                parity_tied = true;
            }
        }
    }
    free(overlaps);
    if (parity_tied ||
        best_margin < std::max(1, opts.min_block_link_reads))
        return std::nullopt;

    return SupportedAlleleEdge{
        best_upstream, best_downstream, best_agree, best_conflict,
        comparisons, best_conflict > best_agree};
}

static void mark_supported_locus(PhasingChunk& chunk, int candidate_i) {
    // Record the exact locus that earned this edge. Co-located rows remain
    // independent candidates, but both alleles must be eligible when reads are
    // rescored: a read carrying one deletion has no exact observation on the
    // complementary deletion row.
    if (candidate_i < 0 ||
        static_cast<size_t>(candidate_i) >= chunk.candidates.size())
        return;
    const CandidateVariant& selected =
        chunk.candidates[static_cast<size_t>(candidate_i)];
    for (CandidateVariant& candidate : chunk.candidates) {
        if (candidate.phase_set == selected.phase_set &&
            candidate.key.sort_pos() == selected.key.sort_pos() &&
            candidate.hap_to_cons_alle[1] >= 0 &&
            candidate.hap_to_cons_alle[2] >= 0 &&
            candidate.hap_to_cons_alle[1] !=
                candidate.hap_to_cons_alle[2]) {
            candidate.gap_link_supported = true;
        }
    }
}

bool stitch_phase_sets_by_alleles(PhasingChunk& chunk,
                                  hts_pos_t upstream_phase_set,
                                  hts_pos_t downstream_phase_set,
                                  const Options& opts) {
    const std::optional<SupportedAlleleEdge> edge = strongest_phase_set_edge(
        chunk, upstream_phase_set, downstream_phase_set, opts);
    if (!edge) return false;
    mark_supported_locus(chunk, edge->upstream_candidate);
    mark_supported_locus(chunk, edge->downstream_candidate);

    return merge_phase_sets_in_place(
        chunk, upstream_phase_set, downstream_phase_set,
        edge->flip);
}


// Combine two explicit candidate sets before voting so each read contributes
// at most one parity decision. Keeping the candidate identities separate from
// their mutable phase-set labels lets recovery test one imported BAM block
// after a preceding flank merge has renamed it.
static std::optional<bool> aggregate_candidate_set_flip(
        const PhasingChunk& chunk, std::vector<int> upstream,
        std::vector<int> downstream, const Options& opts) {
    if (chunk.read_var_cr == nullptr) return std::nullopt;
    const size_t window =
        static_cast<size_t>(std::max(1, opts.block_link_window));
    if (upstream.size() > window)
        upstream.erase(upstream.begin(), upstream.end() - window);
    if (downstream.size() > window)
        downstream.erase(
            downstream.begin() +
                static_cast<std::vector<int>::difference_type>(window),
            downstream.end());
    if (upstream.empty() || downstream.empty()) return std::nullopt;

    int same = 0;
    int cross = 0;
    int64_t* overlaps = nullptr;
    int64_t overlap_capacity = 0;
    const int64_t overlap_count = cr_overlap(
        chunk.read_var_cr.get(), "cr", upstream.front(),
        downstream.back() + 1, &overlaps, &overlap_capacity);
    for (int64_t overlap_i = 0; overlap_i < overlap_count; ++overlap_i) {
        const int64_t read_label = cr_label(
            chunk.read_var_cr.get(), overlaps[overlap_i]);
        if (read_label < 0 ||
            static_cast<size_t>(read_label) >= chunk.reads.size() ||
            static_cast<size_t>(read_label) >= chunk.read_var_profile.size()) {
            continue;
        }
        const size_t read_i = static_cast<size_t>(read_label);
        if (chunk.reads[read_i].is_skipped) continue;
        const ReadVariantProfile& profile = chunk.read_var_profile[read_i];
        if (profile.start_var_idx < 0) continue;
        std::array<std::array<int, 2>, 2> votes{};
        const auto score_block = [&](const std::vector<int>& candidates,
                                     size_t side) {
            for (const int candidate_i : candidates) {
                if (candidate_i < profile.start_var_idx ||
                    candidate_i > profile.end_var_idx) {
                    continue;
                }
                const size_t offset = static_cast<size_t>(
                    candidate_i - profile.start_var_idx);
                if (offset >= profile.alleles.size()) continue;
                const int allele = profile.alleles[offset];
                const CandidateVariant& candidate =
                    chunk.candidates[static_cast<size_t>(candidate_i)];
                if (allele == candidate.hap_to_cons_alle[1])
                    ++votes[side][0];
                else if (allele == candidate.hap_to_cons_alle[2])
                    ++votes[side][1];
            }
        };
        score_block(upstream, 0);
        score_block(downstream, 1);
        if (votes[0][0] == votes[0][1] ||
            votes[1][0] == votes[1][1]) {
            continue;
        }
        const bool upstream_hap2 = votes[0][1] > votes[0][0];
        const bool downstream_hap2 = votes[1][1] > votes[1][0];
        if (upstream_hap2 == downstream_hap2) ++same;
        else ++cross;
    }
    free(overlaps);

    return significant_parity_flip(same, cross);
}

// This graph-only wrapper follows current phase-set labels. Recovery boundaries
// that contain BAM blocks instead pass their immutable candidate memberships to
// aggregate_candidate_set_flip().
static std::optional<bool> aggregate_phase_set_flip(
        const PhasingChunk& chunk, hts_pos_t upstream_phase_set,
        hts_pos_t downstream_phase_set, const Options& opts) {
    const auto oriented_in = [](const CandidateVariant& candidate,
                                hts_pos_t phase_set) {
        const int hap1 = candidate.hap_to_cons_alle[1];
        const int hap2 = candidate.hap_to_cons_alle[2];
        return candidate.phase_set == phase_set && hap1 >= 0 && hap1 <= 1 &&
               hap2 >= 0 && hap2 <= 1 && hap1 != hap2;
    };

    std::vector<int> upstream;
    std::vector<int> downstream;
    for (size_t i = 0; i < chunk.candidates.size(); ++i) {
        if (oriented_in(chunk.candidates[i], upstream_phase_set))
            upstream.push_back(static_cast<int>(i));
        if (oriented_in(chunk.candidates[i], downstream_phase_set))
            downstream.push_back(static_cast<int>(i));
    }
    return aggregate_candidate_set_flip(
        chunk, std::move(upstream), std::move(downstream), opts);
}

constexpr int kRecoveryPathMinEdgeMargin = 10;
constexpr int kRecoveryPathMinOrientationMargin = 2;
constexpr hts_pos_t kRecoveryPathMaxStepBp = 20000;
constexpr size_t kRecoveryPathMinFlankAnchors = 2;

struct RecoveryPairEvidence {
    int same = 0;
    int cross = 0;
    bool right_ref = false;
    bool right_alt = false;
};

static RecoveryPairEvidence recovery_pair_evidence(
        const PhasingChunk& chunk, int left_i, int right_i) {
    RecoveryPairEvidence evidence;
    if (left_i < 0 || right_i < 0 || left_i >= right_i ||
        static_cast<size_t>(right_i) >= chunk.candidates.size() ||
        chunk.read_var_cr == nullptr) {
        return evidence;
    }
    int64_t* overlaps = nullptr;
    int64_t overlap_capacity = 0;
    const int64_t overlap_count = cr_overlap(
        chunk.read_var_cr.get(), "cr", left_i, right_i + 1,
        &overlaps, &overlap_capacity);
    for (int64_t oi = 0; oi < overlap_count; ++oi) {
        const size_t read_i = static_cast<size_t>(
            cr_label(chunk.read_var_cr.get(), overlaps[oi]));
        if (read_i >= chunk.reads.size() ||
            read_i >= chunk.read_var_profile.size() ||
            chunk.reads[read_i].is_skipped) {
            continue;
        }
        const ReadVariantProfile& profile = chunk.read_var_profile[read_i];
        if (profile.start_var_idx < 0 ||
            left_i < profile.start_var_idx || right_i > profile.end_var_idx) {
            continue;
        }
        const size_t left_offset = static_cast<size_t>(
            left_i - profile.start_var_idx);
        const size_t right_offset = static_cast<size_t>(
            right_i - profile.start_var_idx);
        if (left_offset >= profile.alleles.size() ||
            right_offset >= profile.alleles.size()) {
            continue;
        }
        const int left_allele = profile.alleles[left_offset];
        const int right_allele = profile.alleles[right_offset];
        if ((left_allele != 0 && left_allele != 1) ||
            (right_allele != 0 && right_allele != 1)) {
            continue;
        }
        evidence.right_ref |= right_allele == 0;
        evidence.right_alt |= right_allele == 1;
        if (left_allele == right_allele) ++evidence.same;
        else ++evidence.cross;
    }
    free(overlaps);
    return evidence;
}

struct RecoveryPathScore {
    bool valid = false;
    int bottleneck = 0;
    int total_margin = 0;
    int total_support = 0;
    int edges = 0;
};

static bool better_recovery_path(const RecoveryPathScore& lhs,
                                 const RecoveryPathScore& rhs) {
    if (!lhs.valid) return false;
    if (!rhs.valid) return true;
    if (lhs.bottleneck != rhs.bottleneck)
        return lhs.bottleneck > rhs.bottleneck;
    if (lhs.total_margin != rhs.total_margin)
        return lhs.total_margin > rhs.total_margin;
    if (lhs.total_support != rhs.total_support)
        return lhs.total_support > rhs.total_support;
    return lhs.edges < rhs.edges;
}

static RecoveryPathScore extend_recovery_path(
        const RecoveryPathScore& prefix, int support, int conflict) {
    const int margin = support - conflict;
    if (margin < kRecoveryPathMinEdgeMargin) return {};
    if (!prefix.valid) {
        return RecoveryPathScore{true, margin, margin, support, 1};
    }
    return RecoveryPathScore{
        true, std::min(prefix.bottleneck, margin),
        prefix.total_margin + margin,
        prefix.total_support + support, prefix.edges + 1};
}

struct RecoveryDiploidPath {
    bool flip_right = false;
    // A direct relation between two centered SNPs can validate a graph/BAM
    // edge even when the callers describe no sequence-identical site.
    bool direct_snp_bridge = false;
    // Candidate index and the allele assigned to haplotype 1 in the left
    // phase set's gauge.
    std::vector<std::pair<int, int>> sites;
};

// Solve an ordered two-state DAG. A node is one eligible injected BAM site in
// one of its two allele orientations. An edge exists when two sites are within
// the molecule-scale distance limit and their shared reads favor one parity by
// the required margin. Jumping over a site omits it from the chain.
//
// The recurrence is exact for this DAG and lexicographic objective: maximize
// the weakest edge first, then total margin, then supporting observations, and
// finally prefer fewer edges. Both right-flank orientations are solved; an
// ambiguous winner is rejected. This is deliberately a path solver, not a
// whole-matrix MEC claim.
static std::optional<RecoveryDiploidPath> optimal_recovery_diploid_path(
        const PhasingChunk& chunk, const RecoverySeam& window,
        hts_pos_t left_phase_set, hts_pos_t right_phase_set) {
    std::vector<int> left_anchors;
    std::vector<int> right_anchors;
    std::vector<int> internal;
    std::vector<bool> seen_ref(chunk.candidates.size(), false);
    std::vector<bool> seen_alt(chunk.candidates.size(), false);

    for (size_t read_i = 0; read_i < chunk.read_var_profile.size(); ++read_i) {
        if (read_i >= chunk.reads.size() || chunk.reads[read_i].is_skipped)
            continue;
        const ReadVariantProfile& profile = chunk.read_var_profile[read_i];
        if (profile.start_var_idx < 0) continue;
        for (size_t offset = 0; offset < profile.alleles.size(); ++offset) {
            const int candidate_i = profile.start_var_idx +
                                    static_cast<int>(offset);
            if (candidate_i < 0 ||
                static_cast<size_t>(candidate_i) >= chunk.candidates.size()) {
                break;
            }
            if (profile.alleles[offset] == 0)
                seen_ref[static_cast<size_t>(candidate_i)] = true;
            else if (profile.alleles[offset] == 1)
                seen_alt[static_cast<size_t>(candidate_i)] = true;
        }
    }

    for (size_t i = 0; i < chunk.candidates.size(); ++i) {
        const CandidateVariant& candidate = chunk.candidates[i];
        const int hap1 = candidate.hap_to_cons_alle[1];
        const int hap2 = candidate.hap_to_cons_alle[2];
        const bool oriented = hap1 >= 0 && hap1 <= 1 &&
                              hap2 >= 0 && hap2 <= 1 && hap1 != hap2;
        if (candidate.phase_set == left_phase_set && oriented)
            left_anchors.push_back(static_cast<int>(i));
        if (candidate.phase_set == right_phase_set && oriented)
            right_anchors.push_back(static_cast<int>(i));

        const hts_pos_t pos = candidate.key.sort_pos();
        const bool exact_bam_site =
            candidate.key.type == VariantType::Snp || candidate.alignment_verified;
        if (pos > window.beg && pos < window.end && candidate.bam_injected &&
            candidate.phase_set != left_phase_set &&
            candidate.phase_set != right_phase_set && exact_bam_site &&
            seen_ref[i] && seen_alt[i] &&
            (candidate.lcd_var_i_to_cate & kCandGermlineVarCate) != 0) {
            internal.push_back(static_cast<int>(i));
        }
    }
    // One anchor cannot corroborate the orientation of its phase set. The
    // ordinary stitcher may still join such blocks directly, but the recovery
    // fallback must not use a single noisy row to absorb an entire flank.
    if (left_anchors.size() < kRecoveryPathMinFlankAnchors ||
        right_anchors.size() < kRecoveryPathMinFlankAnchors ||
        internal.empty()) {
        return std::nullopt;
    }

    struct State {
        RecoveryPathScore score;
        int prev_node = -1;
        int prev_orientation = -1;
    };
    std::vector<std::array<State, 2>> dp(internal.size());
    std::map<std::pair<int, int>, RecoveryPairEvidence> evidence_cache;
    const auto evidence = [&](int left_i, int right_i) {
        const auto key = std::make_pair(left_i, right_i);
        const auto found = evidence_cache.find(key);
        if (found != evidence_cache.end()) return found->second;
        const RecoveryPairEvidence value =
            recovery_pair_evidence(chunk, left_i, right_i);
        evidence_cache.emplace(key, value);
        return value;
    };

    for (size_t j = 0; j < internal.size(); ++j) {
        const int right_i = internal[j];
        const hts_pos_t right_pos =
            chunk.candidates[static_cast<size_t>(right_i)].key.sort_pos();
        for (int right_orientation = 0; right_orientation <= 1;
             ++right_orientation) {
            State best;
            for (const int left_i : left_anchors) {
                const CandidateVariant& left =
                    chunk.candidates[static_cast<size_t>(left_i)];
                const hts_pos_t left_pos = left.key.sort_pos();
                if (left_pos >= right_pos ||
                    right_pos - left_pos > kRecoveryPathMaxStepBp) {
                    continue;
                }
                const RecoveryPairEvidence edge = evidence(left_i, right_i);
                if (!edge.right_ref || !edge.right_alt) continue;
                const bool same_orientation =
                    left.hap_to_cons_alle[1] == right_orientation;
                const int support = same_orientation ? edge.same : edge.cross;
                const int conflict = same_orientation ? edge.cross : edge.same;
                const RecoveryPathScore score =
                    extend_recovery_path({}, support, conflict);
                if (better_recovery_path(score, best.score)) {
                    best.score = score;
                    best.prev_node = -1;
                    best.prev_orientation = -1;
                }
            }
            for (size_t i = 0; i < j; ++i) {
                const int left_i = internal[i];
                const hts_pos_t left_pos =
                    chunk.candidates[static_cast<size_t>(left_i)].key.sort_pos();
                if (left_pos >= right_pos ||
                    right_pos - left_pos > kRecoveryPathMaxStepBp) {
                    continue;
                }
                const RecoveryPairEvidence edge = evidence(left_i, right_i);
                if (!edge.right_ref || !edge.right_alt) continue;
                for (int left_orientation = 0; left_orientation <= 1;
                     ++left_orientation) {
                    if (!dp[i][left_orientation].score.valid) continue;
                    const bool same_orientation =
                        left_orientation == right_orientation;
                    const int support = same_orientation ? edge.same : edge.cross;
                    const int conflict = same_orientation ? edge.cross : edge.same;
                    const RecoveryPathScore score = extend_recovery_path(
                        dp[i][left_orientation].score, support, conflict);
                    if (better_recovery_path(score, best.score)) {
                        best.score = score;
                        best.prev_node = static_cast<int>(i);
                        best.prev_orientation = left_orientation;
                    }
                }
            }
            dp[j][right_orientation] = best;
        }
    }

    struct Endpoint {
        RecoveryPathScore score;
        int node = -1;
        int orientation = -1;
    };
    std::array<Endpoint, 2> endpoints;
    for (int flip = 0; flip <= 1; ++flip) {
        Endpoint best;
        for (size_t i = 0; i < internal.size(); ++i) {
            const int left_i = internal[i];
            const hts_pos_t left_pos =
                chunk.candidates[static_cast<size_t>(left_i)].key.sort_pos();
            for (int left_orientation = 0; left_orientation <= 1;
                 ++left_orientation) {
                if (!dp[i][left_orientation].score.valid) continue;
                for (const int right_i : right_anchors) {
                    const CandidateVariant& right =
                        chunk.candidates[static_cast<size_t>(right_i)];
                    const hts_pos_t right_pos = right.key.sort_pos();
                    if (left_pos >= right_pos ||
                        right_pos - left_pos > kRecoveryPathMaxStepBp) {
                        continue;
                    }
                    const int right_orientation = flip
                        ? right.hap_to_cons_alle[2]
                        : right.hap_to_cons_alle[1];
                    const RecoveryPairEvidence edge = evidence(left_i, right_i);
                    const bool same_orientation =
                        left_orientation == right_orientation;
                    const int support = same_orientation ? edge.same : edge.cross;
                    const int conflict = same_orientation ? edge.cross : edge.same;
                    const RecoveryPathScore score = extend_recovery_path(
                        dp[i][left_orientation].score, support, conflict);
                    if (better_recovery_path(score, best.score)) {
                        best.score = score;
                        best.node = static_cast<int>(i);
                        best.orientation = left_orientation;
                    }
                }
            }
        }
        endpoints[flip] = best;
    }

    const int winner = better_recovery_path(
        endpoints[1].score, endpoints[0].score) ? 1 : 0;
    const int runner = 1 - winner;
    if (!endpoints[winner].score.valid) return std::nullopt;
    if (endpoints[runner].score.valid &&
        endpoints[winner].score.bottleneck -
                endpoints[runner].score.bottleneck <
            kRecoveryPathMinOrientationMargin) {
        return std::nullopt;
    }

    RecoveryDiploidPath path;
    path.flip_right = winner != 0;
    int node = endpoints[winner].node;
    int orientation = endpoints[winner].orientation;
    while (node >= 0) {
        path.sites.emplace_back(internal[static_cast<size_t>(node)], orientation);
        const State& state = dp[static_cast<size_t>(node)][orientation];
        orientation = state.prev_orientation;
        node = state.prev_node;
    }
    std::reverse(path.sites.begin(), path.sites.end());
    return path;
}

static bool commit_recovery_diploid_path(
        PhasingChunk& chunk, const RecoveryDiploidPath& path,
        hts_pos_t left_phase_set, hts_pos_t right_phase_set) {
    // A DP node may already belong to a locally phased BAM block. Determine
    // each block's required parity before mutating anything: two nodes from the
    // same block must agree on one atomic flip.
    std::map<hts_pos_t, bool> local_flips;
    for (const auto& [candidate_i, orientation] : path.sites) {
        const CandidateVariant& candidate =
            chunk.candidates[static_cast<size_t>(candidate_i)];
        if (candidate.phase_set <= 0 ||
            candidate.phase_set == left_phase_set ||
            candidate.phase_set == right_phase_set) {
            continue;
        }
        const int hap1 = candidate.hap_to_cons_alle[1];
        const int hap2 = candidate.hap_to_cons_alle[2];
        if (hap1 < 0 || hap1 > 1 || hap2 < 0 || hap2 > 1 || hap1 == hap2)
            continue;
        const bool flip = hap1 != orientation;
        const auto [it, inserted] =
            local_flips.emplace(candidate.phase_set, flip);
        if (!inserted && it->second != flip) return false;
    }

    for (const auto& [phase_set, flip] : local_flips) {
        if (!merge_phase_sets_in_place(
                chunk, left_phase_set, phase_set, flip)) {
            return false;
        }
    }
    if (!merge_phase_sets_in_place(
            chunk, left_phase_set, right_phase_set, path.flip_right)) {
        return false;
    }

    for (const auto& [candidate_i, orientation] : path.sites) {
        CandidateVariant& candidate =
            chunk.candidates[static_cast<size_t>(candidate_i)];
        const int hap1 = candidate.hap_to_cons_alle[1];
        const int hap2 = candidate.hap_to_cons_alle[2];
        const bool already_oriented =
            candidate.phase_set == left_phase_set &&
            hap1 >= 0 && hap1 <= 1 && hap2 >= 0 && hap2 <= 1 &&
            hap1 != hap2;
        if (!already_oriented) {
            candidate.hap_to_cons_alle[1] = orientation;
            candidate.hap_to_cons_alle[2] = 1 - orientation;
            candidate.hap_alt = orientation == 1 ? 1 : 2;
            candidate.hap_ref = orientation == 0 ? 1 : 2;
            candidate.phase_set = left_phase_set;
        }
        candidate.gap_link_supported = true;
    }
    return true;
}


struct RecoveryMecRow {
    int fixed_mismatches = 0;
    int fixed_observations = 0;
    std::map<size_t, std::array<int, 2>> variable_mismatches;
};

struct RecoveryMecEffect {
    size_t row = 0;
    int mismatch0 = 0;
    int mismatch1 = 0;
    int observations = 0;
};

struct RecoveryMecOptimum {
    int score = INT_MAX;
    std::vector<int> bits;
};

// Find the exact minimum-error diploid assignment for one fixed orientation of
// the right block. A read pays the smaller Hamming distance to haplotype 1 or
// its complement. Branch-and-bound is exact: its lower bound lets every
// unassigned observation choose any mismatch count, which can only
// underestimate the attainable cost and therefore cannot prune an optimum.
static RecoveryMecOptimum solve_recovery_mec(
        const std::vector<RecoveryMecRow>& rows,
        size_t variable_count) {
    std::vector<std::vector<RecoveryMecEffect>> effects(variable_count);
    std::vector<int> total_observations(rows.size(), 0);
    std::vector<int> mismatches(rows.size(), 0);
    std::vector<int> assigned_observations(rows.size(), 0);
    for (size_t row_i = 0; row_i < rows.size(); ++row_i) {
        const RecoveryMecRow& row = rows[row_i];
        mismatches[row_i] = row.fixed_mismatches;
        assigned_observations[row_i] = row.fixed_observations;
        total_observations[row_i] = row.fixed_observations;
        for (const auto& [variable, costs] : row.variable_mismatches) {
            const int observations = costs[0] + costs[1];
            if (variable >= variable_count || observations == 0) continue;
            effects[variable].push_back(
                RecoveryMecEffect{row_i, costs[0], costs[1], observations});
            total_observations[row_i] += observations;
        }
    }

    std::vector<size_t> order(variable_count);
    for (size_t i = 0; i < variable_count; ++i) order[i] = i;
    std::stable_sort(order.begin(), order.end(), [&](size_t a, size_t b) {
        const auto support = [&](size_t variable) {
            int n = 0;
            for (const RecoveryMecEffect& effect : effects[variable])
                n += effect.observations;
            return n;
        };
        return support(a) > support(b);
    });

    const auto lower_bound = [&]() {
        int bound = 0;
        for (size_t row_i = 0; row_i < rows.size(); ++row_i) {
            const int remaining = total_observations[row_i] -
                                  assigned_observations[row_i];
            const int low = mismatches[row_i];
            const int high = low + remaining;
            const int total = total_observations[row_i];
            bound += std::min(std::min(low, total - low),
                              std::min(high, total - high));
        }
        return bound;
    };
    const auto score_assignment = [&](const std::vector<int>& assignment) {
        int score = 0;
        for (const RecoveryMecRow& row : rows) {
            int mismatch = row.fixed_mismatches;
            int observations = row.fixed_observations;
            for (const auto& [variable, costs] : row.variable_mismatches) {
                if (variable >= assignment.size()) continue;
                mismatch += costs[static_cast<size_t>(assignment[variable])];
                observations += costs[0] + costs[1];
            }
            score += std::min(mismatch, observations - mismatch);
        }
        return score;
    };

    // A coordinate-descent seed supplies a tight feasible upper bound before
    // exact search. It affects only search order and pruning; branch-and-bound
    // still proves that no lower score exists.
    std::vector<int> seed(variable_count, 0);
    int seed_score = score_assignment(seed);
    bool improved = true;
    while (improved) {
        improved = false;
        for (size_t variable : order) {
            seed[variable] ^= 1;
            const int flipped_score = score_assignment(seed);
            if (flipped_score < seed_score) {
                seed_score = flipped_score;
                improved = true;
            } else {
                seed[variable] ^= 1;
            }
        }
    }

    RecoveryMecOptimum optimum;
    optimum.score = seed_score;
    optimum.bits = seed;
    std::vector<int> bits(variable_count, 0);
    std::function<void(size_t)> search = [&](size_t depth) {
        const int bound = lower_bound();
        // The incumbent is already feasible. Equal-score assignments cannot
        // change source-to-sink parity, so only a strict improvement matters.
        if (bound >= optimum.score) return;
        if (depth == order.size()) {
            optimum.score = bound;
            optimum.bits = bits;
            return;
        }

        const size_t variable = order[depth];
        for (int branch = 0; branch <= 1; ++branch) {
            const int bit = branch == 0 ? seed[variable]
                                        : 1 - seed[variable];
            bits[variable] = bit;
            for (const RecoveryMecEffect& effect : effects[variable]) {
                mismatches[effect.row] +=
                    bit == 0 ? effect.mismatch0 : effect.mismatch1;
                assigned_observations[effect.row] += effect.observations;
            }
            search(depth + 1);
            for (const RecoveryMecEffect& effect : effects[variable]) {
                mismatches[effect.row] -=
                    bit == 0 ? effect.mismatch0 : effect.mismatch1;
                assigned_observations[effect.row] -= effect.observations;
            }
        }
    };
    search(0);
    return optimum;
}

// Solve one adjacent recovery edge with the same binary MEC objective used by
// HiPhase, over observations pgphase already collected. Clean, allele-balanced
// SNPs define the primary problem. Centered indels enter only when SNPs alone
// do not connect the two blocks; one verified boundary row may substitute when
// a split multi-allelic event makes row-wise AF off center. The full read set
// and two deterministic, disjoint read halves must choose the same unique
// parity; this stability test replaces an arbitrary vote-margin threshold.
static std::optional<RecoveryDiploidPath> trusted_recovery_mec_path(
        const PhasingChunk& chunk, const RecoverySeam& edge,
        hts_pos_t left_phase_set, hts_pos_t right_phase_set,
        bool boundary_scope = false,
        bool* exceeded_variable_limit = nullptr) {
    constexpr double kTrustedAlleleFractionMargin = 0.12;
    constexpr size_t kTrustedMecMaxVariables = 20;
    if (exceeded_variable_limit != nullptr)
        *exceeded_variable_limit = false;
    const size_t candidate_count = chunk.candidates.size();

    std::vector<bool> seen_ref(candidate_count, false);
    std::vector<bool> seen_alt(candidate_count, false);
    for (size_t read_i = 0; read_i < chunk.read_var_profile.size(); ++read_i) {
        if (read_i >= chunk.reads.size() || chunk.reads[read_i].is_skipped)
            continue;
        const ReadVariantProfile& profile = chunk.read_var_profile[read_i];
        if (profile.start_var_idx < 0) continue;
        for (size_t offset = 0; offset < profile.alleles.size(); ++offset) {
            const int candidate_i = profile.start_var_idx +
                                    static_cast<int>(offset);
            if (candidate_i < 0 ||
                static_cast<size_t>(candidate_i) >= candidate_count) break;
            if (profile.alleles[offset] == 0)
                seen_ref[static_cast<size_t>(candidate_i)] = true;
            else if (profile.alleles[offset] == 1)
                seen_alt[static_cast<size_t>(candidate_i)] = true;
        }
    }

    const auto oriented = [](const CandidateVariant& candidate) {
        return candidate.phase_set > 0 &&
               candidate.hap_to_cons_alle[1] >= 0 &&
               candidate.hap_to_cons_alle[1] <= 1 &&
               candidate.hap_to_cons_alle[2] >= 0 &&
               candidate.hap_to_cons_alle[2] <= 1 &&
               candidate.hap_to_cons_alle[1] !=
                   candidate.hap_to_cons_alle[2];
    };
    const auto centered = [&](const CandidateVariant& candidate) {
        return candidate.counts.ref_cov > 0 && candidate.counts.alt_cov > 0 &&
               std::abs(candidate.counts.allele_fraction - 0.5) <=
                   kTrustedAlleleFractionMargin;
    };
    const auto closest_anchor = [&](hts_pos_t phase_set, hts_pos_t target,
                                    VariantType required_type,
                                    bool allow_verified_off_center)
            -> std::optional<int> {
        std::optional<int> best;
        hts_pos_t best_distance = std::numeric_limits<hts_pos_t>::max();
        for (size_t i = 0; i < candidate_count; ++i) {
            const CandidateVariant& candidate = chunk.candidates[i];
            const bool trusted_off_center_indel =
                allow_verified_off_center &&
                required_type != VariantType::Snp &&
                candidate.bam_injected && candidate.alignment_verified;
            if (candidate.phase_set != phase_set || !oriented(candidate) ||
                (!centered(candidate) && !trusted_off_center_indel) ||
                candidate.key.type != required_type ||
                !seen_ref[i] || !seen_alt[i]) {
                continue;
            }
            const hts_pos_t distance =
                std::llabs(candidate.key.sort_pos() - target);
            if (!best || distance < best_distance) {
                best = static_cast<int>(i);
                best_distance = distance;
            }
        }
        return best;
    };

    std::optional<int> left_anchor =
        closest_anchor(left_phase_set, edge.beg, VariantType::Snp, false);
    std::optional<int> right_anchor =
        closest_anchor(right_phase_set, edge.end, VariantType::Snp, false);
    bool allow_indels = false;
    // Prefer a centered indel of either kind. A verified injected boundary row
    // is the final fallback: splitting a diploid multi-allelic event into two
    // separate BAM rows can move both row-wise AFs away from 0.5 even though
    // their complementary read alleles determine the phase unambiguously.
    if (!left_anchor)
        left_anchor = closest_anchor(
            left_phase_set, edge.beg, VariantType::Insertion, false);
    if (!left_anchor)
        left_anchor = closest_anchor(
            left_phase_set, edge.beg, VariantType::Deletion, false);
    if (!left_anchor)
        left_anchor = closest_anchor(
            left_phase_set, edge.beg, VariantType::Insertion, true);
    if (!left_anchor)
        left_anchor = closest_anchor(
            left_phase_set, edge.beg, VariantType::Deletion, true);
    if (!right_anchor)
        right_anchor = closest_anchor(
            right_phase_set, edge.end, VariantType::Insertion, false);
    if (!right_anchor)
        right_anchor = closest_anchor(
            right_phase_set, edge.end, VariantType::Deletion, false);
    if (!right_anchor)
        right_anchor = closest_anchor(
            right_phase_set, edge.end, VariantType::Insertion, true);
    if (!right_anchor)
        right_anchor = closest_anchor(
            right_phase_set, edge.end, VariantType::Deletion, true);
    if (!left_anchor || !right_anchor) return std::nullopt;
    if (chunk.candidates[static_cast<size_t>(*left_anchor)].key.type !=
            VariantType::Snp ||
        chunk.candidates[static_cast<size_t>(*right_anchor)].key.type !=
            VariantType::Snp) {
        allow_indels = true;
    }

    hts_pos_t solve_beg = std::min(
        chunk.candidates[static_cast<size_t>(*left_anchor)].key.sort_pos(),
        chunk.candidates[static_cast<size_t>(*right_anchor)].key.sort_pos());
    hts_pos_t solve_end = std::max(
        chunk.candidates[static_cast<size_t>(*left_anchor)].key.sort_pos(),
        chunk.candidates[static_cast<size_t>(*right_anchor)].key.sort_pos());
    if (!boundary_scope) {
        // The normal solve validates the read-connected extent of both atomic
        // blocks. It can therefore expose an older polarity change away from
        // the immediate boundary, but unrelated unphased sites may make that
        // exact problem too large.
        for (size_t i = 0; i < candidate_count; ++i) {
            const CandidateVariant& candidate = chunk.candidates[i];
            if ((candidate.phase_set != left_phase_set &&
                 candidate.phase_set != right_phase_set) ||
                !oriented(candidate) || !centered(candidate) ||
                !seen_ref[i] || !seen_alt[i]) {
                continue;
            }
            const bool boundary_uses_indels =
                chunk.candidates[static_cast<size_t>(*left_anchor)].key.type !=
                    VariantType::Snp ||
                chunk.candidates[static_cast<size_t>(*right_anchor)].key.type !=
                    VariantType::Snp;
            if (!boundary_uses_indels &&
                candidate.key.type != VariantType::Snp) {
                continue;
            }
            solve_beg = std::min(solve_beg, candidate.key.sort_pos());
            solve_end = std::max(solve_end, candidate.key.sort_pos());
        }
    }
    // Boundary scope keeps only the selected-anchor interval. The caller may
    // request it solely after the complete atomic-block problem exceeds the
    // exact-search bound and independent candidate evidence validates the same
    // graph/BAM gauge. It cannot override a tie or contradictory full solve.

    const auto build_problem = [&](bool include_indels) {
        std::vector<int> nodes;
        std::vector<int> node_of(candidate_count, -1);
        std::set<std::tuple<hts_pos_t, int, int, std::string, hts_pos_t>> seen_keys;
        for (size_t i = 0; i < candidate_count; ++i) {
            const CandidateVariant& candidate = chunk.candidates[i];
            const hts_pos_t pos = candidate.key.sort_pos();
            const bool selected_boundary =
                static_cast<int>(i) == *left_anchor ||
                static_cast<int>(i) == *right_anchor;
            if (pos < solve_beg || pos > solve_end ||
                (!centered(candidate) && !selected_boundary) ||
                !seen_ref[i] || !seen_alt[i] ||
                (!include_indels && candidate.key.type != VariantType::Snp)) {
                continue;
            }
            const auto key = std::make_tuple(
                pos, static_cast<int>(candidate.key.type), candidate.key.ref_len,
                candidate.key.alt, candidate.phase_set);
            if (!seen_keys.insert(key).second) continue;
            node_of[i] = static_cast<int>(nodes.size());
            nodes.push_back(static_cast<int>(i));
        }

        std::vector<int> parent(nodes.size());
        for (size_t i = 0; i < parent.size(); ++i)
            parent[i] = static_cast<int>(i);
        const auto root_of = [&](int node) {
            int root = node;
            while (parent[static_cast<size_t>(root)] != root)
                root = parent[static_cast<size_t>(root)];
            return root;
        };
        const auto join = [&](int a, int b) {
            const int root_a = root_of(a);
            const int root_b = root_of(b);
            if (root_a != root_b)
                parent[static_cast<size_t>(root_b)] = root_a;
        };
        for (size_t read_i = 0; read_i < chunk.read_var_profile.size(); ++read_i) {
            if (read_i >= chunk.reads.size() || chunk.reads[read_i].is_skipped)
                continue;
            const ReadVariantProfile& profile = chunk.read_var_profile[read_i];
            if (profile.start_var_idx < 0) continue;
            int first = -1;
            for (size_t offset = 0; offset < profile.alleles.size(); ++offset) {
                if (profile.alleles[offset] != 0 &&
                    profile.alleles[offset] != 1) continue;
                const int candidate_i = profile.start_var_idx +
                                        static_cast<int>(offset);
                if (candidate_i < 0 ||
                    static_cast<size_t>(candidate_i) >= node_of.size()) break;
                const int node = node_of[static_cast<size_t>(candidate_i)];
                if (node < 0) continue;
                if (first < 0) first = node;
                else join(first, node);
            }
        }

        std::set<int> left_roots;
        std::set<int> right_roots;
        for (size_t node = 0; node < nodes.size(); ++node) {
            const CandidateVariant& candidate =
                chunk.candidates[static_cast<size_t>(nodes[node])];
            if (candidate.phase_set == left_phase_set && oriented(candidate))
                left_roots.insert(root_of(static_cast<int>(node)));
            if (candidate.phase_set == right_phase_set && oriented(candidate))
                right_roots.insert(root_of(static_cast<int>(node)));
        }
        std::set<int> bridge_roots;
        std::set_intersection(
            left_roots.begin(), left_roots.end(), right_roots.begin(),
            right_roots.end(), std::inserter(bridge_roots,
                                             bridge_roots.begin()));
        std::vector<bool> retained(candidate_count, false);
        for (size_t node = 0; node < nodes.size(); ++node) {
            if (bridge_roots.count(root_of(static_cast<int>(node))) != 0)
                retained[static_cast<size_t>(nodes[node])] = true;
        }
        return retained;
    };

    std::vector<bool> retained = build_problem(allow_indels);
    const auto has_both_flanks = [&]() {
        bool left = false;
        bool right = false;
        for (size_t i = 0; i < candidate_count; ++i) {
            if (!retained[i]) continue;
            left = left || chunk.candidates[i].phase_set == left_phase_set;
            right = right || chunk.candidates[i].phase_set == right_phase_set;
        }
        return left && right;
    };
    if (!has_both_flanks() && !allow_indels) {
        allow_indels = true;
        retained = build_problem(true);
    }
    if (!has_both_flanks()) return std::nullopt;

    std::vector<int> variable_of(candidate_count, -1);
    std::vector<int> base_orientation(candidate_count, 0);
    std::vector<int> variable_candidates;
    std::vector<bool> is_left(candidate_count, false);
    std::vector<bool> is_right(candidate_count, false);
    for (size_t i = 0; i < candidate_count; ++i) {
        if (!retained[i]) continue;
        const CandidateVariant& candidate = chunk.candidates[i];
        is_left[i] = candidate.phase_set == left_phase_set && oriented(candidate);
        is_right[i] = candidate.phase_set == right_phase_set && oriented(candidate);
        if (is_left[i] || is_right[i]) continue;
        variable_of[i] = static_cast<int>(variable_candidates.size());
        variable_candidates.push_back(static_cast<int>(i));
        base_orientation[i] = oriented(candidate)
                                  ? candidate.hap_to_cons_alle[1]
                                  : 0;
    }
    if (variable_candidates.size() > kTrustedMecMaxVariables) {
        if (exceeded_variable_limit != nullptr)
            *exceeded_variable_limit = true;
        return std::nullopt;
    }

    int indel_cost_bound = 0;
    if (allow_indels) {
        for (const ReadVariantProfile& profile : chunk.read_var_profile) {
            if (profile.start_var_idx < 0) continue;
            for (size_t offset = 0; offset < profile.alleles.size(); ++offset) {
                if (profile.alleles[offset] != 0 &&
                    profile.alleles[offset] != 1) continue;
                const int candidate_i = profile.start_var_idx +
                                        static_cast<int>(offset);
                if (candidate_i < 0 ||
                    static_cast<size_t>(candidate_i) >= retained.size()) break;
                if (!retained[static_cast<size_t>(candidate_i)]) continue;
                const CandidateVariant& candidate =
                    chunk.candidates[static_cast<size_t>(candidate_i)];
                if (candidate.key.type != VariantType::Snp)
                    indel_cost_bound += phase_matrix_var_weight(candidate);
            }
        }
    }
    const int snp_multiplier = indel_cost_bound + 1;
    const auto stable_hash = [](const std::string& value) {
        uint64_t hash = 14695981039346656037ULL;
        for (const unsigned char byte : value) {
            hash ^= byte;
            hash *= 1099511628211ULL;
        }
        return hash;
    };
    const auto solve_fold = [&](int fold) {
        std::array<RecoveryMecOptimum, 2> result;
        for (int right_flip = 0; right_flip <= 1; ++right_flip) {
            std::vector<RecoveryMecRow> rows;
            for (size_t read_i = 0; read_i < chunk.read_var_profile.size();
                 ++read_i) {
                if (read_i >= chunk.reads.size() ||
                    chunk.reads[read_i].is_skipped ||
                    (fold >= 0 &&
                     static_cast<int>(stable_hash(chunk.reads[read_i].qname) & 1ULL) !=
                         fold)) {
                    continue;
                }
                const ReadVariantProfile& profile = chunk.read_var_profile[read_i];
                if (profile.start_var_idx < 0) continue;
                RecoveryMecRow row;
                int observed_sites = 0;
                for (size_t offset = 0; offset < profile.alleles.size(); ++offset) {
                    const int allele = profile.alleles[offset];
                    if (allele != 0 && allele != 1) continue;
                    const int candidate_i = profile.start_var_idx +
                                            static_cast<int>(offset);
                    if (candidate_i < 0 ||
                        static_cast<size_t>(candidate_i) >= retained.size()) break;
                    const size_t ci = static_cast<size_t>(candidate_i);
                    if (!retained[ci]) continue;
                    const CandidateVariant& candidate = chunk.candidates[ci];
                    int weight = phase_matrix_var_weight(candidate);
                    if (candidate.key.type == VariantType::Snp)
                        weight *= snp_multiplier;
                    ++observed_sites;
                    if (is_left[ci] || is_right[ci]) {
                        int orientation = candidate.hap_to_cons_alle[1];
                        if (is_right[ci] && right_flip != 0)
                            orientation = 1 - orientation;
                        row.fixed_mismatches +=
                            weight * static_cast<int>(allele != orientation);
                        row.fixed_observations += weight;
                    } else {
                        const int variable = variable_of[ci];
                        if (variable < 0) continue;
                        std::array<int, 2>& costs =
                            row.variable_mismatches[static_cast<size_t>(variable)];
                        const int mismatch0 = allele != base_orientation[ci];
                        costs[0] += weight * mismatch0;
                        costs[1] += weight * (1 - mismatch0);
                    }
                }
                // One heterozygous observation can always match one of the
                // two haplotypes, so it contributes zero information to MEC
                // parity while increasing the branch-and-bound workload.
                if (observed_sites >= 2) rows.push_back(std::move(row));
            }
            result[static_cast<size_t>(right_flip)] =
                solve_recovery_mec(rows, variable_candidates.size());
        }
        if (result[0].score == result[1].score)
            return std::make_pair(std::optional<bool>{}, result);
        return std::make_pair(
            std::optional<bool>{result[1].score < result[0].score}, result);
    };

    const auto full = solve_fold(-1);
    const auto half0 = solve_fold(0);
    const auto half1 = solve_fold(1);
    if (!full.first || !half0.first || !half1.first ||
        *full.first != *half0.first || *full.first != *half1.first) {
        return std::nullopt;
    }

    const RecoveryMecOptimum& winner =
        full.second[static_cast<size_t>(*full.first)];
    if (winner.bits.size() != variable_candidates.size()) return std::nullopt;
    RecoveryDiploidPath path;
    path.flip_right = *full.first;

    // A pair of distinct SNPs is a valid representation-independent anchor.
    // Use the nearest selected boundary pair, then require the full reads and
    // both disjoint halves to select the MEC parity independently. Restricting
    // this check to one pair keeps recovery linear in the read count and avoids
    // selecting the strongest result from many correlated SNP comparisons.
    const size_t left_i = static_cast<size_t>(*left_anchor);
    const size_t right_i = static_cast<size_t>(*right_anchor);
    if (chunk.candidates[left_i].key.type == VariantType::Snp &&
        chunk.candidates[right_i].key.type == VariantType::Snp &&
        retained[left_i] && retained[right_i]) {
        const auto direct_snp_flip = [&](int fold) -> std::optional<bool> {
            int same = 0;
            int cross = 0;
            const int left_hap1 =
                chunk.candidates[left_i].hap_to_cons_alle[1];
            const int right_hap1 =
                chunk.candidates[right_i].hap_to_cons_alle[1];
            for (size_t read_i = 0; read_i < chunk.read_var_profile.size();
                 ++read_i) {
                if (read_i >= chunk.reads.size() ||
                    chunk.reads[read_i].is_skipped ||
                    (fold >= 0 &&
                     static_cast<int>(stable_hash(
                         chunk.reads[read_i].qname) & 1ULL) != fold)) {
                    continue;
                }
                const ReadVariantProfile& profile =
                    chunk.read_var_profile[read_i];
                if (profile.start_var_idx < 0 ||
                    static_cast<int>(left_i) < profile.start_var_idx ||
                    static_cast<int>(right_i) < profile.start_var_idx ||
                    static_cast<int>(left_i) > profile.end_var_idx ||
                    static_cast<int>(right_i) > profile.end_var_idx) {
                    continue;
                }
                const size_t left_offset =
                    left_i - static_cast<size_t>(profile.start_var_idx);
                const size_t right_offset =
                    right_i - static_cast<size_t>(profile.start_var_idx);
                if (left_offset >= profile.alleles.size() ||
                    right_offset >= profile.alleles.size()) {
                    continue;
                }
                const int left_allele = profile.alleles[left_offset];
                const int right_allele = profile.alleles[right_offset];
                if ((left_allele != 0 && left_allele != 1) ||
                    (right_allele != 0 && right_allele != 1)) {
                    continue;
                }
                const bool flip =
                    (left_allele == left_hap1) !=
                    (right_allele == right_hap1);
                if (flip) ++cross;
                else ++same;
            }
            return parity_flip_at_p(
                same, cross, 1, kRecoveryCandidateAnchorPValue);
        };
        const std::optional<bool> pair_full = direct_snp_flip(-1);
        const std::optional<bool> pair_half0 = direct_snp_flip(0);
        const std::optional<bool> pair_half1 = direct_snp_flip(1);
        path.direct_snp_bridge =
            pair_full && pair_half0 && pair_half1 &&
            *pair_full == path.flip_right &&
            *pair_half0 == path.flip_right &&
            *pair_half1 == path.flip_right;
    }

    for (size_t variable = 0; variable < variable_candidates.size(); ++variable) {
        const int candidate_i = variable_candidates[variable];
        path.sites.emplace_back(
            candidate_i,
            base_orientation[static_cast<size_t>(candidate_i)] ^
                winner.bits[variable]);
    }
    return path;
}

// Optimize every biallelic site in the read-connected source-to-sink component.
// Existing local phase sets are variables as atomic units, while an unphased
// site is one variable. This preserves the relative orientation already solved
// inside a BAM block and lets many individually weak observations establish one
// globally supported parity. Tied block parity is left unjoined; equal internal
// optima need not be enumerated because they cannot change that block parity.
static std::optional<RecoveryDiploidPath> optimal_recovery_mec_path(
        const PhasingChunk& chunk, const RecoverySeam& window,
        hts_pos_t left_phase_set, hts_pos_t right_phase_set) {
    const size_t candidate_count = chunk.candidates.size();
    std::vector<bool> seen_ref(candidate_count, false);
    std::vector<bool> seen_alt(candidate_count, false);
    for (size_t read_i = 0; read_i < chunk.read_var_profile.size(); ++read_i) {
        if (read_i >= chunk.reads.size() || chunk.reads[read_i].is_skipped)
            continue;
        const ReadVariantProfile& profile = chunk.read_var_profile[read_i];
        if (profile.start_var_idx < 0) continue;
        for (size_t offset = 0; offset < profile.alleles.size(); ++offset) {
            const int candidate_i = profile.start_var_idx +
                                    static_cast<int>(offset);
            if (candidate_i < 0 ||
                static_cast<size_t>(candidate_i) >= candidate_count) break;
            if (profile.alleles[offset] == 0)
                seen_ref[static_cast<size_t>(candidate_i)] = true;
            else if (profile.alleles[offset] == 1)
                seen_alt[static_cast<size_t>(candidate_i)] = true;
        }
    }

    std::vector<int> nodes;
    std::vector<int> node_of(candidate_count, -1);
    std::vector<bool> is_left(candidate_count, false);
    std::vector<bool> is_right(candidate_count, false);
    std::vector<bool> is_internal(candidate_count, false);
    for (size_t candidate_i = 0; candidate_i < candidate_count; ++candidate_i) {
        const CandidateVariant& candidate = chunk.candidates[candidate_i];
        const int hap1 = candidate.hap_to_cons_alle[1];
        const int hap2 = candidate.hap_to_cons_alle[2];
        const bool oriented = hap1 >= 0 && hap1 <= 1 &&
                              hap2 >= 0 && hap2 <= 1 && hap1 != hap2;
        is_left[candidate_i] =
            candidate.phase_set == left_phase_set && oriented;
        is_right[candidate_i] =
            candidate.phase_set == right_phase_set && oriented;
        const hts_pos_t pos = candidate.key.sort_pos();
        const bool exact_bam_site =
            candidate.bam_injected &&
            (candidate.key.type == VariantType::Snp ||
             candidate.alignment_verified);
        is_internal[candidate_i] =
            pos > window.beg && pos < window.end && exact_bam_site &&
            candidate.phase_set != left_phase_set &&
            candidate.phase_set != right_phase_set &&
            seen_ref[candidate_i] && seen_alt[candidate_i] &&
            (candidate.lcd_var_i_to_cate & kCandGermlineVarCate) != 0;
        if (is_left[candidate_i] || is_right[candidate_i] ||
            is_internal[candidate_i]) {
            node_of[candidate_i] = static_cast<int>(nodes.size());
            nodes.push_back(static_cast<int>(candidate_i));
        }
    }
    if (nodes.empty()) return std::nullopt;

    std::vector<int> parent(nodes.size());
    for (size_t i = 0; i < parent.size(); ++i) parent[i] = static_cast<int>(i);
    const auto root_of = [&](int node) {
        int root = node;
        while (parent[static_cast<size_t>(root)] != root)
            root = parent[static_cast<size_t>(root)];
        return root;
    };
    const auto join = [&](int a, int b) {
        const int root_a = root_of(a);
        const int root_b = root_of(b);
        if (root_a != root_b) parent[static_cast<size_t>(root_b)] = root_a;
    };
    for (size_t read_i = 0; read_i < chunk.read_var_profile.size(); ++read_i) {
        if (read_i >= chunk.reads.size() || chunk.reads[read_i].is_skipped)
            continue;
        const ReadVariantProfile& profile = chunk.read_var_profile[read_i];
        if (profile.start_var_idx < 0) continue;
        int first_node = -1;
        for (size_t offset = 0; offset < profile.alleles.size(); ++offset) {
            if (profile.alleles[offset] != 0 && profile.alleles[offset] != 1)
                continue;
            const int candidate_i = profile.start_var_idx +
                                    static_cast<int>(offset);
            if (candidate_i < 0 ||
                static_cast<size_t>(candidate_i) >= node_of.size()) break;
            const int node = node_of[static_cast<size_t>(candidate_i)];
            if (node < 0) continue;
            if (first_node < 0) first_node = node;
            else join(first_node, node);
        }
    }

    std::vector<bool> component_left(nodes.size(), false);
    std::vector<bool> component_right(nodes.size(), false);
    for (size_t node = 0; node < nodes.size(); ++node) {
        const size_t candidate_i = static_cast<size_t>(nodes[node]);
        const size_t root = static_cast<size_t>(root_of(static_cast<int>(node)));
        component_left[root] = component_left[root] || is_left[candidate_i];
        component_right[root] = component_right[root] || is_right[candidate_i];
    }
    std::vector<bool> retained(candidate_count, false);
    bool connected = false;
    for (size_t node = 0; node < nodes.size(); ++node) {
        const size_t root = static_cast<size_t>(root_of(static_cast<int>(node)));
        if (!component_left[root] || !component_right[root]) continue;
        retained[static_cast<size_t>(nodes[node])] = true;
        connected = true;
    }
    if (!connected) return std::nullopt;

    struct Variable {
        hts_pos_t phase_set = 0;
        std::vector<int> candidates;
    };
    std::vector<Variable> variables;
    std::map<hts_pos_t, size_t> variable_of_phase_set;
    std::vector<int> variable_of(candidate_count, -1);
    std::vector<int> base_orientation(candidate_count, 0);
    for (size_t candidate_i = 0; candidate_i < candidate_count; ++candidate_i) {
        if (!retained[candidate_i] || !is_internal[candidate_i]) continue;
        const CandidateVariant& candidate = chunk.candidates[candidate_i];
        const int hap1 = candidate.hap_to_cons_alle[1];
        const int hap2 = candidate.hap_to_cons_alle[2];
        const bool oriented = hap1 >= 0 && hap1 <= 1 &&
                              hap2 >= 0 && hap2 <= 1 && hap1 != hap2;
        size_t variable = variables.size();
        if (candidate.phase_set > 0 && oriented) {
            const auto [it, inserted] = variable_of_phase_set.emplace(
                candidate.phase_set, variables.size());
            if (inserted)
                variables.push_back(Variable{candidate.phase_set, {}});
            variable = it->second;
            base_orientation[candidate_i] = hap1;
        } else {
            variables.push_back(Variable{0, {}});
            base_orientation[candidate_i] = 0;
        }
        variable_of[candidate_i] = static_cast<int>(variable);
        variables[variable].candidates.push_back(static_cast<int>(candidate_i));
    }

    // With no internal BAM site this is a direct block edge, for which the
    // aggregate voter already enforces the validated support threshold. MEC is
    // reserved for the case it adds: optimizing a diploid chain through one or
    // more injected sites.
    if (variables.empty()) return std::nullopt;

    // Exact diploid MEC is NP-hard. Targeted recovery gaps are small (the
    // chr20 panel needs at most 13 injected variables), while an unrelated
    // 52.8 kb seam contains 131. Keep production latency bounded by abstaining
    // before search; every solution returned below is still the exact optimum.
    constexpr size_t kRecoveryMecMaxVariables = 20;
    if (variables.size() > kRecoveryMecMaxVariables) return std::nullopt;

    const auto build_rows = [&](bool flip_right) {
        std::vector<RecoveryMecRow> rows;
        for (size_t read_i = 0; read_i < chunk.read_var_profile.size(); ++read_i) {
            if (read_i >= chunk.reads.size() || chunk.reads[read_i].is_skipped)
                continue;
            const ReadVariantProfile& profile = chunk.read_var_profile[read_i];
            if (profile.start_var_idx < 0) continue;
            RecoveryMecRow row;
            int observations = 0;
            for (size_t offset = 0; offset < profile.alleles.size(); ++offset) {
                const int allele = profile.alleles[offset];
                if (allele != 0 && allele != 1) continue;
                const int candidate_i = profile.start_var_idx +
                                        static_cast<int>(offset);
                if (candidate_i < 0 ||
                    static_cast<size_t>(candidate_i) >= retained.size()) break;
                const size_t ci = static_cast<size_t>(candidate_i);
                if (!retained[ci]) continue;
                ++observations;
                if (is_left[ci] || is_right[ci]) {
                    int orientation = chunk.candidates[ci].hap_to_cons_alle[1];
                    if (is_right[ci] && flip_right) orientation = 1 - orientation;
                    row.fixed_mismatches += allele != orientation;
                    ++row.fixed_observations;
                } else {
                    const int variable = variable_of[ci];
                    if (variable < 0) continue;
                    std::array<int, 2>& costs =
                        row.variable_mismatches[static_cast<size_t>(variable)];
                    const int mismatch0 = allele != base_orientation[ci];
                    costs[0] += mismatch0;
                    costs[1] += 1 - mismatch0;
                }
            }
            if (observations >= 2) rows.push_back(std::move(row));
        }
        return rows;
    };

    const std::vector<RecoveryMecRow> same_rows = build_rows(false);
    const std::vector<RecoveryMecRow> flipped_rows = build_rows(true);
    const RecoveryMecOptimum same =
        solve_recovery_mec(same_rows, variables.size());
    const RecoveryMecOptimum flipped =
        solve_recovery_mec(flipped_rows, variables.size());
    if (same.score == flipped.score) return std::nullopt;
    const bool flip_right = flipped.score < same.score;
    const RecoveryMecOptimum& winner = flip_right ? flipped : same;
    if (winner.bits.size() != variables.size()) return std::nullopt;

    RecoveryDiploidPath path;
    path.flip_right = flip_right;
    for (size_t variable = 0; variable < variables.size(); ++variable) {
        for (const int candidate_i : variables[variable].candidates) {
            path.sites.emplace_back(
                candidate_i,
                base_orientation[static_cast<size_t>(candidate_i)] ^
                    winner.bits[variable]);
        }
    }
    return path;
}

// Visit the spatial chain of graph and recovered BAM phase sets from left to
// right. A BAM subsolve may contain several phase sets whose HP labels are
// independent. Keep those blocks separate until direct read evidence determines
// the parity of one adjacent pair. Graph-only seams retain the ordinary allele,
// aggregate, and exact-path fallbacks.
size_t stitch_recovery_phase_sets_left_to_right(
        PhasingChunk& chunk,
        const std::vector<RecoverySeam>& windows,
        const std::vector<RecoveryPhaseGauge>& gauges,
        const Options& opts) {
    const auto oriented = [](const CandidateVariant& candidate) {
        return candidate.phase_set > 0 &&
               candidate.hap_to_cons_alle[1] >= 0 &&
               candidate.hap_to_cons_alle[2] >= 0 &&
               candidate.hap_to_cons_alle[1] != candidate.hap_to_cons_alle[2];
    };
    const int min_support = std::max(1, opts.min_block_link_reads);
    const auto orientation_in_gauge = [min_support](
            const RecoveryPhaseGauge& gauge, hts_pos_t phase_set) {
        if (std::find(gauge.imported_phase_sets.begin(),
                      gauge.imported_phase_sets.end(), phase_set) !=
            gauge.imported_phase_sets.end()) {
            return 1;  // imported candidates already use the BAM solve's gauge
        }
        const auto vote = std::find_if(
            gauge.graph_votes.begin(), gauge.graph_votes.end(),
            [phase_set](const PhaseSetGaugeVote& v) {
                return v.phase_set == phase_set;
            });
        if (vote == gauge.graph_votes.end() ||
            std::abs(vote->same - vote->cross) < min_support)
            return 0;
        return vote->same > vote->cross ? 1 : -1;
    };
    const auto significant_graph_orientation = [](
            const RecoveryPhaseGauge& gauge, hts_pos_t phase_set) {
        const auto vote = std::find_if(
            gauge.graph_votes.begin(), gauge.graph_votes.end(),
            [phase_set](const PhaseSetGaugeVote& candidate) {
                return candidate.phase_set == phase_set;
            });
        if (vote == gauge.graph_votes.end()) return 0;
        const std::optional<bool> flip =
            significant_parity_flip(vote->same, vote->cross);
        if (!flip) return 0;
        return *flip ? -1 : 1;
    };
    const auto candidate_anchor_flip = [](
            const RecoveryPhaseGauge& gauge, hts_pos_t graph_phase_set,
            hts_pos_t bam_phase_set) -> std::optional<bool> {
        const auto vote = std::find_if(
            gauge.block_votes.begin(), gauge.block_votes.end(),
            [graph_phase_set, bam_phase_set](
                    const RecoveryBlockGaugeVote& candidate) {
                return (candidate.graph_phase_set == graph_phase_set &&
                        candidate.bam_phase_set == bam_phase_set) ||
                       (candidate.graph_phase_set == bam_phase_set &&
                        candidate.bam_phase_set == graph_phase_set);
            });
        if (vote == gauge.block_votes.end()) return std::nullopt;
        return parity_flip_at_p(
            vote->shared_candidate_same, vote->shared_candidate_cross, 1,
            kRecoveryCandidateAnchorPValue);
    };
    const auto find_block_gauge_vote = [](
            const RecoveryPhaseGauge& gauge, hts_pos_t first_phase_set,
            hts_pos_t second_phase_set) {
        return std::find_if(
            gauge.block_votes.begin(), gauge.block_votes.end(),
            [first_phase_set, second_phase_set](
                    const RecoveryBlockGaugeVote& candidate) {
                return (candidate.graph_phase_set == first_phase_set &&
                        candidate.bam_phase_set == second_phase_set) ||
                       (candidate.graph_phase_set == second_phase_set &&
                        candidate.bam_phase_set == first_phase_set);
            });
    };
    const auto block_read_gauge_flip = [&find_block_gauge_vote](
            const RecoveryPhaseGauge& gauge, hts_pos_t first_phase_set,
            hts_pos_t second_phase_set) -> std::optional<bool> {
        const auto vote = find_block_gauge_vote(
            gauge, first_phase_set, second_phase_set);
        if (vote == gauge.block_votes.end()) return std::nullopt;
        return significant_parity_flip(vote->counts);
    };
    const auto block_gauge_flip = [&find_block_gauge_vote](
            const RecoveryPhaseGauge& gauge, hts_pos_t first_phase_set,
            hts_pos_t second_phase_set) -> std::optional<bool> {
        const auto vote = find_block_gauge_vote(
            gauge, first_phase_set, second_phase_set);
        if (vote == gauge.block_votes.end()) return std::nullopt;
        const std::optional<bool> read_flip =
            significant_parity_flip(vote->counts);
        if (!read_flip) return std::nullopt;

        // Sequence-identical clean heterozygotes validate a decisive molecule
        // vote. They cannot establish a join alone: one BAM phase block may
        // contain an internal switch, so its two distant consensus anchors can
        // each be locally correct while their implied long-range join is not.
        const bool have_candidate_vote =
            vote->shared_candidate_same > 0 ||
            vote->shared_candidate_cross > 0;
        if (!have_candidate_vote) return std::nullopt;
        if (vote->shared_candidate_same > 0 &&
            vote->shared_candidate_cross == 0) {
            return *read_flip ? std::nullopt : read_flip;
        }
        if (vote->shared_candidate_cross > 0 &&
            vote->shared_candidate_same == 0) {
            return *read_flip ? read_flip : std::nullopt;
        }
        return std::nullopt;
    };

    // Block votes are measured before any recovery stitch mutates the chunk.
    // Keep one oriented candidate from each original phase set so a later vote
    // can be translated into the current gauge after an earlier flank merge.
    std::map<hts_pos_t, std::pair<size_t, int>> initial_orientations;
    std::map<hts_pos_t, std::vector<int>> initial_phase_set_candidates;
    std::vector<hts_pos_t> initial_phase_sets;
    initial_phase_sets.reserve(chunk.candidates.size());
    for (size_t candidate_i = 0; candidate_i < chunk.candidates.size();
         ++candidate_i) {
        const CandidateVariant& candidate = chunk.candidates[candidate_i];
        initial_phase_sets.push_back(candidate.phase_set);
        if (!oriented(candidate)) continue;
        initial_orientations.try_emplace(
            candidate.phase_set,
            std::make_pair(candidate_i, candidate.hap_to_cons_alle[1]));
        initial_phase_set_candidates[candidate.phase_set].push_back(
            static_cast<int>(candidate_i));
    }
    const auto flipped_from_initial = [&](hts_pos_t phase_set) {
        const auto initial = initial_orientations.find(phase_set);
        if (initial == initial_orientations.end()) return false;
        const CandidateVariant& candidate =
            chunk.candidates[initial->second.first];
        return candidate.hap_to_cons_alle[1] != initial->second.second;
    };

    // Earlier seams can absorb a phase set that is a later seam's left side.
    // Preserve that relabeling explicitly so the detector's original IDs stay
    // usable while the chunk is updated in place from left to right.
    std::map<hts_pos_t, hts_pos_t> phase_set_aliases;
    // A component first created beyond an unsupported edge may grow only on
    // aggregate/DP evidence or a statistically supported edge from an injected BAM site.
    // Preserve that rule when a following seam refers to the component by its
    // surviving phase-set label.
    std::set<hts_pos_t> strong_only_roots;
    // Exact MEC may otherwise bridge around an intentionally independent BAM
    // block after the first pass abstains at one of its boundaries.
    std::set<std::pair<hts_pos_t, hts_pos_t>> bam_block_windows;
    const auto resolve_phase_set = [&phase_set_aliases](hts_pos_t phase_set) {
        for (;;) {
            const auto alias = phase_set_aliases.find(phase_set);
            if (alias == phase_set_aliases.end() || alias->second == phase_set)
                return phase_set;
            phase_set = alias->second;
        }
    };

    size_t joined = 0;
    // A recovered block may appear in several overlapping seam windows. Keep
    // this set for the whole chunk so the trusted fallback cannot attach the
    // same original gauge twice through different windows.
    std::set<hts_pos_t> trusted_used_blocks;
    for (const RecoverySeam& window : windows) {
        const auto gauge_it = std::find_if(
            gauges.begin(), gauges.end(), [&](const RecoveryPhaseGauge& gauge) {
                return gauge.beg <= window.beg && gauge.end >= window.end;
            });
        const RecoveryPhaseGauge* recovery_gauge =
            gauge_it != gauges.end() ? &*gauge_it : nullptr;

        // Candidate order is reference order. Keep the first occurrence of
        // each local PS so the chain follows the blocks spatially.
        std::vector<std::pair<size_t, hts_pos_t>> local_blocks;
        for (size_t candidate_i = 0; candidate_i < chunk.candidates.size();
             ++candidate_i) {
            const CandidateVariant& candidate = chunk.candidates[candidate_i];
            const hts_pos_t pos = candidate.key.sort_pos();
            if (!candidate.bam_injected || !oriented(candidate) ||
                pos <= window.beg || pos >= window.end)
                continue;
            const hts_pos_t source_phase_set =
                initial_phase_sets[candidate_i];
            if (source_phase_set <= 0) continue;
            const bool seen = std::any_of(
                local_blocks.begin(), local_blocks.end(),
                [&](const auto& block) {
                    return block.second == source_phase_set;
                });
            if (!seen)
                local_blocks.emplace_back(candidate_i, source_phase_set);
        }
        if (!local_blocks.empty())
            bam_block_windows.emplace(window.beg, window.end);

        // The detector captured these identities while canonical graph and BAM
        // coordinates were both available. Reusing them avoids attempting to
        // rediscover an indel flank from a differently padded representation.
        const hts_pos_t left_phase_set =
            resolve_phase_set(window.left_phase_set);
        const hts_pos_t right_phase_set =
            resolve_phase_set(window.right_phase_set);
        if (left_phase_set <= 0 || right_phase_set <= 0) continue;

        struct ChainBlock {
            hts_pos_t current_phase_set;
            hts_pos_t gauge_phase_set;
        };
        std::vector<ChainBlock> chain;
        chain.reserve(local_blocks.size() + 2);
        chain.push_back(ChainBlock{left_phase_set, window.left_phase_set});
        for (const auto& block : local_blocks) {
            const hts_pos_t current = resolve_phase_set(block.second);
            const bool seen = std::any_of(
                chain.begin(), chain.end(),
                [current](const ChainBlock& item) {
                    return item.current_phase_set == current;
                });
            if (current > 0 && !seen)
                chain.push_back(ChainBlock{current, block.second});
        }
        const bool right_seen = std::any_of(
            chain.begin(), chain.end(),
            [right_phase_set](const ChainBlock& item) {
                return item.current_phase_set == right_phase_set;
            });
        if (!right_seen)
            chain.push_back(
                ChainBlock{right_phase_set, window.right_phase_set});

        const auto is_imported_phase_set = [&](hts_pos_t phase_set) {
            if (recovery_gauge != nullptr &&
                std::find(recovery_gauge->imported_phase_sets.begin(),
                          recovery_gauge->imported_phase_sets.end(),
                          phase_set) !=
                    recovery_gauge->imported_phase_sets.end()) {
                return true;
            }
            return std::any_of(
                local_blocks.begin(), local_blocks.end(),
                [phase_set](const auto& block) {
                    return block.second == phase_set;
                });
        };
        std::optional<bool> outer_relation;
        const auto left_candidates =
            initial_phase_set_candidates.find(window.left_phase_set);
        const auto right_candidates =
            initial_phase_set_candidates.find(window.right_phase_set);
        if (left_candidates != initial_phase_set_candidates.end() &&
            right_candidates != initial_phase_set_candidates.end()) {
            const std::optional<bool> current_relation =
                aggregate_candidate_set_flip(
                    chunk, left_candidates->second,
                    right_candidates->second, opts);
            if (current_relation) {
                // Convert the current candidate orientations back to the
                // pre-stitch gauges used by the saved phase-set identities.
                outer_relation =
                    *current_relation ^
                    flipped_from_initial(window.left_phase_set) ^
                    flipped_from_initial(window.right_phase_set);
            }
        }
        // Every recovered chain must preserve a statistically decisive
        // relation between the already phased graph flanks. Prefer direct
        // outer-candidate observations. A multi-block BAM solve can also
        // provide a whole-window graph gauge; when both exist they must agree.
        std::optional<bool> gauge_outer_relation;
        if (local_blocks.size() > 1 && recovery_gauge != nullptr) {
            const std::optional<bool> left_anchor = candidate_anchor_flip(
                *recovery_gauge, window.left_phase_set,
                local_blocks.front().second);
            const std::optional<bool> right_anchor = candidate_anchor_flip(
                *recovery_gauge, window.right_phase_set,
                local_blocks.back().second);
            if (left_anchor && right_anchor) {
                const int left_orientation = significant_graph_orientation(
                    *recovery_gauge, window.left_phase_set);
                const int right_orientation = significant_graph_orientation(
                    *recovery_gauge, window.right_phase_set);
                if (left_orientation != 0 && right_orientation != 0) {
                    gauge_outer_relation =
                        left_orientation != right_orientation;
                }
            }
        }
        if (outer_relation && gauge_outer_relation &&
            *outer_relation != *gauge_outer_relation) {
            continue;
        }
        if (!outer_relation) outer_relation = gauge_outer_relation;

        const auto merge_imported_boundary =
            [&](hts_pos_t upstream_phase_set,
                hts_pos_t downstream_phase_set, bool flip) {
                const hts_pos_t left_root =
                    resolve_phase_set(window.left_phase_set);
                const hts_pos_t right_root =
                    resolve_phase_set(window.right_phase_set);
                if (outer_relation && left_root != right_root) {
                    std::optional<bool> resulting_relation;
                    if (upstream_phase_set == left_root &&
                        downstream_phase_set == right_root) {
                        resulting_relation =
                            flipped_from_initial(window.left_phase_set) ^
                            flipped_from_initial(window.right_phase_set) ^
                            flip;
                    } else if (upstream_phase_set == right_root &&
                               downstream_phase_set == left_root) {
                        resulting_relation =
                            (flipped_from_initial(window.left_phase_set) ^
                             flip) ^
                            flipped_from_initial(window.right_phase_set);
                    }
                    if (resulting_relation &&
                        *resulting_relation != *outer_relation) {
                        return false;
                    }
                }
                return merge_phase_sets_in_place(
                    chunk, upstream_phase_set, downstream_phase_set, flip);
            };

        // Imported recovery is an atomic seam transaction. A prefix merge that
        // fails at a later edge adds no continuity and can only relabel reads
        // into a less accurate block. Snapshot the mutable phase state and keep
        // it only if the complete left-to-right chain closes.
        const CandidateTable candidates_before_window = chunk.candidates;
        const std::vector<int> haps_before_window = chunk.haps;
        const std::vector<hts_pos_t> phase_sets_before_window =
            chunk.phase_sets;
        const auto aliases_before_window = phase_set_aliases;
        const auto strong_roots_before_window = strong_only_roots;
        const size_t joined_before_window = joined;

        const hts_pos_t primary_upstream_phase_set =
            chain.front().current_phase_set;
        bool merged_into_upstream = false;
        bool component_contiguous =
            strong_only_roots.count(primary_upstream_phase_set) == 0;
        int upstream_orientation = recovery_gauge != nullptr
            ? orientation_in_gauge(*recovery_gauge,
                                   chain.front().gauge_phase_set)
            : 0;
        if (recovery_gauge != nullptr && upstream_orientation == 0 &&
            chain.front().current_phase_set !=
                chain.front().gauge_phase_set) {
            // A preceding seam may already have absorbed this boundary. The
            // targeted BAM solve then votes for the surviving graph label,
            // while the detector correctly retains the original local label.
            upstream_orientation = orientation_in_gauge(
                *recovery_gauge, chain.front().current_phase_set);
        }
        for (size_t i = 1; i < chain.size(); ++i) {
            // A failed edge separates two components; it must not hide a later
            // supported edge. Resolve the immediately preceding block after
            // every merge so each component still grows from left to right.
            const hts_pos_t upstream_phase_set =
                resolve_phase_set(chain[i - 1].current_phase_set);
            const hts_pos_t downstream_phase_set =
                resolve_phase_set(chain[i].current_phase_set);
            if (downstream_phase_set == upstream_phase_set) continue;

            bool stitched = false;
            bool used_strong_fallback = false;
            const bool upstream_imported =
                is_imported_phase_set(chain[i - 1].gauge_phase_set);
            const bool downstream_imported =
                is_imported_phase_set(chain[i].gauge_phase_set);
            const bool imported_boundary =
                upstream_imported || downstream_imported;

            // BAM blocks are complete local phase solutions. Their arbitrary
            // HP gauges carry no relationship until reads spanning an exact
            // adjacent pair establish one. Candidate membership was captured
            // before mutation, so neither this vote nor a graph-flank vote can
            // accidentally pool a different recovered block.
            if (imported_boundary) {
                if (upstream_imported && downstream_imported) {
                    const auto upstream_candidates =
                        initial_phase_set_candidates.find(
                            chain[i - 1].gauge_phase_set);
                    const auto downstream_candidates =
                        initial_phase_set_candidates.find(
                            chain[i].gauge_phase_set);
                    if (upstream_candidates !=
                            initial_phase_set_candidates.end() &&
                        downstream_candidates !=
                            initial_phase_set_candidates.end()) {
                        const std::optional<bool> current_flip =
                            aggregate_candidate_set_flip(
                                chunk, upstream_candidates->second,
                                downstream_candidates->second, opts);
                        if (current_flip) {
                            stitched = merge_imported_boundary(
                                upstream_phase_set, downstream_phase_set,
                                *current_flip);
                        }
                    }
                } else if (recovery_gauge != nullptr) {
                    const std::optional<bool> gauge_flip = block_gauge_flip(
                        *recovery_gauge, chain[i - 1].gauge_phase_set,
                        chain[i].gauge_phase_set);
                    std::optional<bool> current_flip;
                    if (gauge_flip) {
                        current_flip =
                            *gauge_flip ^
                            flipped_from_initial(
                                chain[i - 1].gauge_phase_set) ^
                            flipped_from_initial(chain[i].gauge_phase_set);
                    }
                    if (current_flip) {
                        stitched = merge_imported_boundary(
                            upstream_phase_set, downstream_phase_set,
                            *current_flip);
                    }
                }
                if (!stitched) {
                    component_contiguous = false;
                    continue;
                }
            }
            // One imported block shares one BAM gauge with both graph flanks.
            // Multiple imported blocks do not share an orientation with each
            // other, so only the entry edge may use that common gauge.
            const bool gauge_covers_edge =
                local_blocks.size() <= 1 || i == 1;
            if (!imported_boundary && !stitched && component_contiguous &&
                gauge_covers_edge &&
                recovery_gauge != nullptr && upstream_orientation != 0) {
                const int downstream_orientation =
                    orientation_in_gauge(*recovery_gauge,
                                         chain[i].gauge_phase_set);
                if (downstream_orientation != 0) {
                    const bool gauge_flip =
                        upstream_orientation != downstream_orientation;
                    // Alleles select which locus may later assign reads. The
                    // shared-read gauge still controls the block orientation;
                    // mark the locus only when both evidence paths agree.
                    const std::optional<SupportedAlleleEdge> edge =
                        strongest_phase_set_edge(
                            chunk, upstream_phase_set,
                            downstream_phase_set, opts);
                    if (edge && edge->flip == gauge_flip) {
                        mark_supported_locus(
                            chunk, edge->upstream_candidate);
                        mark_supported_locus(
                            chunk, edge->downstream_candidate);
                    }
                    stitched = merge_phase_sets_in_place(
                        chunk, upstream_phase_set, downstream_phase_set,
                        gauge_flip);
                }
            }
            if (!imported_boundary && !stitched && component_contiguous) {
                stitched = stitch_phase_sets_by_alleles(
                    chunk, upstream_phase_set, downstream_phase_set, opts);
            }
            if (!imported_boundary && !stitched && !component_contiguous) {
                // A failed edge begins a new component, but it must not hide a
                // later, independently decisive site pair. Apply the same
                // parity test used by the aggregate fallback before extending
                // that new component. The ordinary one-read recovery floor is
                // too weak once the outer gauge has disconnected.
                const std::optional<SupportedAlleleEdge> edge =
                    strongest_phase_set_edge(
                        chunk, upstream_phase_set, downstream_phase_set, opts);
                const std::optional<bool> tested_flip = edge
                    ? significant_parity_flip(
                          edge->same, edge->cross, edge->comparisons)
                    : std::nullopt;
                if (edge && tested_flip &&
                    chunk.candidates[static_cast<size_t>(
                        edge->upstream_candidate)].bam_injected) {
                    mark_supported_locus(chunk, edge->upstream_candidate);
                    mark_supported_locus(chunk, edge->downstream_candidate);
                    stitched = merge_phase_sets_in_place(
                        chunk, upstream_phase_set, downstream_phase_set,
                        *tested_flip);
                    used_strong_fallback = stitched;
                }
            }
            if (!imported_boundary && !stitched) {
                const std::optional<bool> aggregate_flip =
                    aggregate_phase_set_flip(
                        chunk, upstream_phase_set,
                        downstream_phase_set, opts);
                if (aggregate_flip) {
                    stitched = merge_phase_sets_in_place(
                        chunk, upstream_phase_set,
                        downstream_phase_set, *aggregate_flip);
                    used_strong_fallback = stitched;
                }
            }
            if (!imported_boundary && !stitched) {
                std::optional<hts_pos_t> right_pos;
                for (const CandidateVariant& candidate : chunk.candidates) {
                    if (candidate.phase_set != downstream_phase_set ||
                        !oriented(candidate)) {
                        continue;
                    }
                    const hts_pos_t pos = candidate.key.sort_pos();
                    if (pos < window.beg || pos > window.end) continue;
                    if (!right_pos || pos < *right_pos) right_pos = pos;
                }
                std::optional<hts_pos_t> left_pos;
                if (right_pos) {
                    for (const CandidateVariant& candidate : chunk.candidates) {
                        if (candidate.phase_set != upstream_phase_set ||
                            !oriented(candidate)) {
                            continue;
                        }
                        const hts_pos_t pos = candidate.key.sort_pos();
                        if (pos < window.beg || pos >= *right_pos) continue;
                        if (!left_pos || pos > *left_pos) left_pos = pos;
                    }
                }
                if (left_pos && right_pos) {
                    const RecoverySeam edge_window{
                        *left_pos, *right_pos, upstream_phase_set,
                        downstream_phase_set};
                    const std::optional<RecoveryDiploidPath> optimal_path =
                        optimal_recovery_diploid_path(
                            chunk, edge_window, upstream_phase_set,
                            downstream_phase_set);
                    if (optimal_path) {
                        stitched = commit_recovery_diploid_path(
                            chunk, *optimal_path, upstream_phase_set,
                            downstream_phase_set);
                        used_strong_fallback = stitched;
                    }
                }
            }
            if (!stitched) {
                component_contiguous = false;
                continue;
            }

            // Redirect every known alias of the absorbed label. Later seams
            // still refer to their original detector IDs.
            for (auto& [source, target] : phase_set_aliases) {
                (void)source;
                if (resolve_phase_set(target) == downstream_phase_set)
                    target = upstream_phase_set;
            }
            phase_set_aliases[downstream_phase_set] = upstream_phase_set;
            phase_set_aliases[chain[i].gauge_phase_set] = upstream_phase_set;
            const bool keep_strong_only = used_strong_fallback ||
                strong_only_roots.count(upstream_phase_set) != 0 ||
                strong_only_roots.count(downstream_phase_set) != 0;
            strong_only_roots.erase(downstream_phase_set);
            if (keep_strong_only)
                strong_only_roots.insert(upstream_phase_set);
            // merge_phase_sets_in_place already transfers and, when needed,
            // flips every read belonging to an imported block. Re-scoring the
            // whole seam would erase the BAM subsolve's independent labels.
            if (!imported_boundary) merged_into_upstream = true;
            ++joined;
        }

        if (!local_blocks.empty() &&
            (!outer_relation ||
             resolve_phase_set(window.left_phase_set) !=
                 resolve_phase_set(window.right_phase_set))) {
            chunk.candidates = candidates_before_window;
            chunk.haps = haps_before_window;
            chunk.phase_sets = phase_sets_before_window;
            phase_set_aliases = aliases_before_window;
            strong_only_roots = strong_roots_before_window;
            joined = joined_before_window;

            // The outer transaction could not prove every edge. Revisit each
            // graph/BAM neighbor with one exact decision flow: solve the full
            // atomic blocks, and use boundary scope only when that problem is
            // over the exact-search bound. Unsupported neighbors remain
            // separate blocks.
            // Attach recovered BAM blocks only to graph blocks. Independent
            // BAM phase blocks remain separate; joining two arbitrary subsolve
            // gauges here can alter reads without closing a graph seam. Never
            // use one original block in two fallback joins during the
            // same pass: chaining locally stable edges can conceal an older
            // polarity change inside an atomic block.
            for (size_t edge_i = 1; edge_i < chain.size(); ++edge_i) {
                const hts_pos_t upstream_source =
                    chain[edge_i - 1].gauge_phase_set;
                const hts_pos_t downstream_source =
                    chain[edge_i].gauge_phase_set;
                const bool upstream_imported =
                    is_imported_phase_set(upstream_source);
                const bool downstream_imported =
                    is_imported_phase_set(downstream_source);
                const bool graph_bam_edge =
                    upstream_imported != downstream_imported;
                if (!graph_bam_edge ||
                    trusted_used_blocks.count(upstream_source) != 0 ||
                    trusted_used_blocks.count(downstream_source) != 0) {
                    continue;
                }
                const hts_pos_t upstream_phase_set =
                    resolve_phase_set(chain[edge_i - 1].current_phase_set);
                const hts_pos_t downstream_phase_set =
                    resolve_phase_set(chain[edge_i].current_phase_set);
                if (upstream_phase_set <= 0 || downstream_phase_set <= 0 ||
                    upstream_phase_set == downstream_phase_set) {
                    continue;
                }
                if (!is_imported_phase_set(
                        chain[edge_i - 1].gauge_phase_set) &&
                    !is_imported_phase_set(
                        chain[edge_i].gauge_phase_set)) {
                    continue;
                }

                const auto upstream_candidates =
                    initial_phase_set_candidates.find(
                        chain[edge_i - 1].gauge_phase_set);
                const auto downstream_candidates =
                    initial_phase_set_candidates.find(
                        chain[edge_i].gauge_phase_set);
                if (upstream_candidates ==
                        initial_phase_set_candidates.end() ||
                    downstream_candidates ==
                        initial_phase_set_candidates.end()) {
                    continue;
                }
                std::optional<hts_pos_t> upstream_pos;
                for (const int candidate_i : upstream_candidates->second) {
                    const hts_pos_t pos = chunk.candidates[
                        static_cast<size_t>(candidate_i)].key.sort_pos();
                    if (!upstream_pos || pos > *upstream_pos)
                        upstream_pos = pos;
                }
                std::optional<hts_pos_t> downstream_pos;
                for (const int candidate_i : downstream_candidates->second) {
                    const hts_pos_t pos = chunk.candidates[
                        static_cast<size_t>(candidate_i)].key.sort_pos();
                    if (!downstream_pos || pos < *downstream_pos)
                        downstream_pos = pos;
                }
                if (!upstream_pos || !downstream_pos) continue;
                const RecoverySeam trusted_edge{
                    std::min(*upstream_pos, *downstream_pos),
                    std::max(*upstream_pos, *downstream_pos),
                    upstream_phase_set, downstream_phase_set};
                bool whole_block_too_large = false;
                std::optional<RecoveryDiploidPath> trusted_path =
                    trusted_recovery_mec_path(
                        chunk, trusted_edge, upstream_phase_set,
                        downstream_phase_set, false,
                        &whole_block_too_large);
                // The narrow scope is a resource fallback only. Never use it
                // to override a full-block tie or conflicting parity.
                if (!trusted_path && whole_block_too_large &&
                    graph_bam_edge && recovery_gauge != nullptr &&
                    candidate_anchor_flip(
                        *recovery_gauge, upstream_source,
                        downstream_source)) {
                    trusted_path = trusted_recovery_mec_path(
                        chunk, trusted_edge, upstream_phase_set,
                        downstream_phase_set, true);
                }
                if (!trusted_path) continue;
                if (graph_bam_edge) {
                    if (recovery_gauge == nullptr) continue;
                    std::optional<bool> gauge_flip = block_gauge_flip(
                        *recovery_gauge, upstream_source, downstream_source);
                    // Distinct centered SNPs with independently significant
                    // full/half evidence replace the exact-key anchor. Keep
                    // the source-specific read gauge as a separate check.
                    if (!gauge_flip && trusted_path->direct_snp_bridge) {
                        gauge_flip = block_read_gauge_flip(
                            *recovery_gauge, upstream_source,
                            downstream_source);
                    }
                    if (!gauge_flip) continue;
                    const bool current_gauge_flip =
                        *gauge_flip ^
                        flipped_from_initial(upstream_source) ^
                        flipped_from_initial(downstream_source);
                    if (current_gauge_flip != trusted_path->flip_right)
                        continue;
                }
                if (!commit_recovery_diploid_path(
                        chunk, *trusted_path, upstream_phase_set,
                        downstream_phase_set)) {
                    continue;
                }

                std::set<hts_pos_t> absorbed{downstream_phase_set};
                for (const auto& site : trusted_path->sites) {
                    const hts_pos_t phase_set = chunk.candidates[
                        static_cast<size_t>(site.first)].phase_set;
                    if (phase_set > 0 && phase_set != upstream_phase_set)
                        absorbed.insert(phase_set);
                }
                for (auto& [source, target] : phase_set_aliases) {
                    (void)source;
                    if (absorbed.count(resolve_phase_set(target)) != 0)
                        target = upstream_phase_set;
                }
                for (const hts_pos_t phase_set : absorbed)
                    phase_set_aliases[phase_set] = upstream_phase_set;
                phase_set_aliases[chain[edge_i].gauge_phase_set] =
                    upstream_phase_set;
                joined += absorbed.size();
                trusted_used_blocks.insert(upstream_source);
                trusted_used_blocks.insert(downstream_source);
            }
            continue;
        }

        if (!merged_into_upstream) continue;
        // Imported HP labels were assigned before the final chain existed and
        // are stale as soon as a later block is flipped. Reuse longcallD's read
        // scorer against the completed, oriented phase set so every retained
        // read derives its HP from the same candidate consensus emitted in VCF.
        // Include previously unphased reads: once the chain orients an exact MSA
        // allele, that allele supplies the phase evidence the subsolve lacked.
        for (size_t read_i = 0; read_i < chunk.reads.size(); ++read_i) {
            if (read_i >= chunk.phase_sets.size() ||
                read_i >= chunk.haps.size() || chunk.reads[read_i].is_skipped)
                continue;
            const hts_pos_t current_phase_set = chunk.phase_sets[read_i];
            if (current_phase_set > 0 &&
                current_phase_set != primary_upstream_phase_set)
                continue;
            const int hap = assign_read_hap_from_recovery_chain(
                chunk, read_i, primary_upstream_phase_set);
            if (hap < 0) continue;
            chunk.haps[read_i] = hap;
            if (hap > 0)
                chunk.phase_sets[read_i] = primary_upstream_phase_set;
        }
    }

    // Preserve the validated left-to-right stitch above as pass one. Exact MEC
    // revisits only seams that remain open after every ordinary decision has
    // settled, so an experimental join cannot change the evidence or policy of
    // a later baseline seam.
    for (const RecoverySeam& window : windows) {
        if (bam_block_windows.count({window.beg, window.end}) != 0) continue;
        const hts_pos_t left_phase_set =
            resolve_phase_set(window.left_phase_set);
        const hts_pos_t right_phase_set =
            resolve_phase_set(window.right_phase_set);
        if (left_phase_set <= 0 || right_phase_set <= 0 ||
            left_phase_set == right_phase_set) {
            continue;
        }

        const std::optional<RecoveryDiploidPath> mec_path =
            optimal_recovery_mec_path(
                chunk, window, left_phase_set, right_phase_set);
        if (!mec_path) continue;

        std::set<hts_pos_t> absorbed_phase_sets;
        absorbed_phase_sets.insert(right_phase_set);
        for (const auto& site : mec_path->sites) {
            const hts_pos_t phase_set =
                chunk.candidates[static_cast<size_t>(site.first)].phase_set;
            if (phase_set > 0 && phase_set != left_phase_set)
                absorbed_phase_sets.insert(phase_set);
        }
        if (!commit_recovery_diploid_path(
                chunk, *mec_path, left_phase_set, right_phase_set)) {
            continue;
        }

        for (auto& [source, target] : phase_set_aliases) {
            (void)source;
            if (absorbed_phase_sets.count(resolve_phase_set(target)) != 0)
                target = left_phase_set;
        }
        for (const hts_pos_t phase_set : absorbed_phase_sets)
            phase_set_aliases[phase_set] = left_phase_set;
        joined += absorbed_phase_sets.size();

        for (size_t read_i = 0; read_i < chunk.reads.size(); ++read_i) {
            if (read_i >= chunk.phase_sets.size() ||
                read_i >= chunk.haps.size() || chunk.reads[read_i].is_skipped)
                continue;
            const hts_pos_t current_phase_set = chunk.phase_sets[read_i];
            if (current_phase_set > 0 && current_phase_set != left_phase_set)
                continue;
            const int hap = assign_read_hap_from_recovery_chain(
                chunk, read_i, left_phase_set);
            if (hap < 0) continue;
            chunk.haps[read_i] = hap;
            if (hap > 0) chunk.phase_sets[read_i] = left_phase_set;
        }
    }
    return joined;
}

bool allele_depths_call_het(const CandidateVariant& var,
                                  const Options& opts);

// Phase-set assignment + flip for one k-means iteration.
// Returns 1 if any flip occurred (changed), 0 if converged.
// Iter_update_var_hap_cons_phase_set.
int iter_update_var_hap_cons_phase_set(PhasingChunk& chunk,
                                               const std::vector<int>& valid_var_idx,
                                               const Options& opts) {
    const int n = (int)valid_var_idx.size();
    // Collect het var positions (indices into valid_var_idx).
    std::vector<int> het_var_idx;
    std::vector<bool> is_het(n, false);
    int seeded_het = 0;
    for (int _vi = 0; _vi < n; ++_vi) {
        CandidateVariant& var = chunk.candidates[valid_var_idx[_vi]];
        // The same validation used for repeat/multiallelic bridges also applies
        // to ordinary MSA indels. Otherwise a rejected site can reconnect the
        // flanks through its abundant reference observations alone.
        // A homopolymer indel is kept out of the LINK list because its allele is
        // unreliable for linking, and the emit loop below hands it the running
        // phase set without applying `parity` -- the same shape as longcallD
        // (assign_hap.c:409 flips inside the `is_het` guard, :418 assigns the
        // phase set to every candidate).
        //
        // That asymmetry is NOT the cause of the switch described below, which
        // was tested directly: applying the running parity to every site that
        // only inherits the phase set moved 128 genotypes on chr20 and, scored
        // against each site's own block orientation from read truth, 31 of them
        // went from agreeing to disagreeing against 18 the other way. The solve
        // is iterative -- the caller re-derives read labels from the flipped
        // consensus and calls this again -- so by convergence an inherited site
        // is already in the block's gauge and applying the parity double-counts
        // it. Rejected and reverted; see
        // evaluations/2026-09-20-inherited-parity/. The mechanism behind the
        // site below is therefore still open. On
        // chr20:48,225,786 (CAAAA>C in an A run, segregation 1.000 against read
        // truth) that put a maternal-on-hap1 site inside a block whose body is
        // paternal-on-hap1 -- a switch invisible in read space, because the
        // reads covering it belong to the next block. A site whose own allele
        // depths call it a clear heterozygote is therefore admitted to the link
        // list, so its orientation is decided by spanning reads like any other
        // het rather than inherited.
        const bool hp_indel_blocks_link =
            var.is_homopolymer_indel &&
            !allele_depths_call_het(var, opts);
        // A record carrying both of the locus' alleles was admitted to the link
        // list only when a gap link had vouched for it, which made every merged
        // multiallelic record invisible to the linker outside gap recovery. On
        // chr20:24,121,714 that site has both alleles, 61 observations, and a
        // phase set, and its two alleles separate the haplotypes perfectly
        // against read truth (+4 on 13 reads all paternal, +8 on 16 all
        // maternal) -- yet the linker's het list in 24.10-24.15 Mb held only
        // clean het SNPs, so the right boundary linked back across 37 kb to the
        // gap's left edge and reported agree=0 conflict=0. Reads behind both
        // alleles are the warrant here, exactly as for any other het: its
        // orientation is then decided by spanning reads rather than inherited.
        if (var.hap_to_cons_alle[1] != -1 && var.hap_to_cons_alle[2] != -1 &&
            var.hap_to_cons_alle[1] != var.hap_to_cons_alle[2] &&
            (opts.upstream_assign_hap ||
             var.msa_insertion_alts.empty() || var.gap_link_supported ||
             two_allele_het(var)) &&
            (opts.upstream_assign_hap
                 ? (!var.is_homopolymer_indel || var.gap_link_supported ||
                    (var.bam_injected && var.alignment_verified))
                 : !hp_indel_blocks_link)) {
            is_het[_vi] = true;
            het_var_idx.push_back(_vi);
        }
    }

    const int n_het = (int)het_var_idx.size();
    // Strongest link found for each het var: which earlier het it attaches to,
    // and the agree/conflict counts of that link.
    std::vector<int> link_h(n_het, -1), link_agree(n_het, 0), link_conflict(n_het, 0);
    int64_t* ovlp_b = nullptr;
    int64_t max_b = 0;
    cgranges_t* cr = chunk.read_var_cr.get();

    // Link each het var to a preceding one by spanning-read evidence.
    //
    // Search the configured number of preceding heterozygotes and retain the
    // strongest supported edge. Ties keep the nearer edge because the scan is
    // nearest first. The ordinary BAM defaults restrict this to one predecessor.
    // The BAM defaults remain the upstream rule (one predecessor and two
    // reads). Recovery can widen this local search through its copied Options;
    // keeping the parameters here avoids a second link implementation.
    const int window = std::max(1, opts.block_link_window);
    const int min_link_reads = std::max(1, opts.min_block_link_reads);
    for (int hi = 1; hi < n_het; ++hi) {
        const int vi = valid_var_idx[het_var_idx[hi]];
        const int lo = std::max(0, hi - window);
        int best_h = -1, best_support = -1, best_a = 0, best_c = 0;
        for (int hj = hi - 1; hj >= lo; --hj) {
            const int vj = valid_var_idx[het_var_idx[hj]];
            int a = 0, c = 0;
            const int64_t ovlp_n = cr_overlap(cr, "cr", vj, vi + 1, &ovlp_b, &max_b);
            for (int64_t oi = 0; oi < ovlp_n; ++oi) {
                const int read_i = (int)cr_label(cr, ovlp_b[oi]);
                if (chunk.reads[read_i].is_skipped) continue;
                const int agree = opts.link_by_alleles
                        ? check_agree_alleles(chunk, read_i, vj, vi)
                        : check_agree_haps(chunk, read_i, chunk.haps[read_i], vj, vi);
                if (agree > 0) a++;
                else if (agree == 0) c++;
            }
            // An MSA bridge must segregate both haplotypes. A large pile of
            // reference-only observations cannot compensate for missing or
            // contradictory support on the other allele.
            const int support = opts.link_by_alleles
                                    ? (a == c ? 0 : std::max(a, c))
                                    : (opts.upstream_assign_hap
                                           ? std::max(a, c)
                                           : (a == c ? 0 : std::max(a, c)));
            // Require a net margin for additional repeat links. Read-level
            // disagreements alone do not establish a wrong block orientation.
            if (!opts.upstream_assign_hap &&
                (chunk.candidates[vi].is_homopolymer_indel ||
                 chunk.candidates[vj].is_homopolymer_indel) &&
                std::abs(a - c) < min_link_reads) continue;
            if (support > best_support) {
                best_support = support; best_h = hj; best_a = a; best_c = c;
            }
        }
        link_h[hi] = best_h;
        link_agree[hi] = best_a;
        link_conflict[hi] = best_c;
    }
    free(ovlp_b);

    // Resolve orientation and phase set per het var.  parity[hi] is the
    // cumulative flip to apply; it is chained through whichever earlier het the
    // variant actually linked to, not blindly through its predecessor.
    std::vector<int> parity(n_het, 0);
    std::vector<hts_pos_t> het_ps(n_het, -1);
    for (int hi = 0; hi < n_het; ++hi) {
        const CandidateVariant& var = chunk.candidates[valid_var_idx[het_var_idx[hi]]];
        if (hi == 0) {
            het_ps[hi] = var.key.sort_pos();
            continue;
        }
        const int hj = link_h[hi];
        const int support = opts.link_by_alleles
                                ? (link_agree[hi] == link_conflict[hi]
                                       ? 0 : std::max(link_agree[hi], link_conflict[hi]))
                                : (opts.upstream_assign_hap
                                       ? std::max(link_agree[hi], link_conflict[hi])
                                       : (link_agree[hi] == link_conflict[hi]
                                              ? 0 : std::max(link_agree[hi], link_conflict[hi])));
        if (hj < 0 || support < min_link_reads) {
            // No sufficiently supported link anywhere in the window -- break.
            // Carry the running parity unchanged, as orientation within a fresh
            // phase set is arbitrary anyway.
            parity[hi] = parity[hi - 1];
            het_ps[hi] = var.key.sort_pos();
        } else {
            parity[hi] = parity[hj] ^ ((link_conflict[hi] > link_agree[hi]) ? 1 : 0);
            het_ps[hi] = het_ps[hj];
        }
    }


    // Map variant index -> het rank, so the emit loop can look up its decision.
    std::vector<int> het_rank(n, -1);
    for (int hi = 0; hi < n_het; ++hi) het_rank[het_var_idx[hi]] = hi;

    int changed = seeded_het;
    hts_pos_t phase_set = -1;
    for (int _vi = 0; _vi < n; ++_vi) {
        const int vi = valid_var_idx[_vi];
        CandidateVariant& var = chunk.candidates[vi];
        if (_vi == 0) {
            phase_set = var.key.sort_pos();
            var.phase_set = phase_set;
            continue;
        }
        const int hi = het_rank[_vi];
        if (opts.verbose >= 2 && hi > 0) {
            std::fprintf(stderr, "%" PRIi64 " %d %d\n",
                         static_cast<int64_t>(var.key.pos),
                         link_agree[hi], link_conflict[hi]);
        }
        if (hi >= 0) {
            phase_set = het_ps[hi];
            if (parity[hi] == 1) {
                changed = 1;
                // assign_hap.c swaps once per hap (1 then 2), returning the
                // consensus to its original order while still reporting change.
                // Allele-link recovery resolves orientation directly from
                // the observed allele pair, so its parity must be applied even
                // while read assignment follows the upstream BAM rule.
                if (!opts.upstream_assign_hap || opts.link_by_alleles)
                    std::swap(var.hap_to_cons_alle[1], var.hap_to_cons_alle[2]);
            }
        }
        var.phase_set = phase_set;
    }
    return changed;
}

// Rebuild allele profiles and consensus; returns 1 if any cons_alle changed.
// Iter_update_var_hap_to_cons_alle.
/// Do this biallelic candidate's own allele depths call it heterozygous? Such a
/// site must not be collapsed to a homozygous consensus by provisional read
/// labels; see the call site. Applies inside a window the retry is re-solving,
/// and everywhere when --joint-het-orientation is set.
/// Does this candidate's OWN allele depth call it a heterozygote?
///
/// The one place `retry_windows` is read, and therefore the only thing the retry
/// scopes positionally. A site admitted here is put on the LINK list, so its
/// orientation is decided by spanning reads like any other het instead of by
/// whatever its own consensus produced in isolation.
///
/// Why that matters: a homopolymer indel is normally kept off the link list
/// because its allele is unreliable for linking, and such a site inherits the
/// running phase set without `parity` being applied. Applying that parity was
/// tested and rejected -- it regressed 31 site orientations against 18 improved
/// (evaluations/2026-09-20-inherited-parity/) -- so admission to the link list,
/// where spanning reads decide the orientation, remains the route that fixes
/// the site below. On
/// chr20:48,225,786 (CAAAA>C in an A run, segregating 1.000 against read truth)
/// that put a maternal-on-hap1 site inside a paternal-on-hap1 block: a switch
/// invisible in read space, because the reads covering it belong to the NEXT
/// block, and visible only at site level.
///
/// The tests, in order, and each is a real exclusion:
///   - off unless a retry window exists or --joint-het-orientation is set;
///   - category must be noisy het or clean het (SNP or indel);
///   - multiallelic records are excluded -- they are oriented jointly elsewhere;
///   - both haplotype profiles must hold at least two observations, so a site
///     with one side empty cannot claim to segregate;
///   - both reference and alternate depth at or above min_alt_depth;
///   - allele fraction inside [min_af, max_af].
/// With --joint-het-orientation the verdict applies chunk-wide; otherwise the
/// position must fall inside one of `retry_windows`, which is what keeps the
/// widened admission confined to the intervals the first solve failed on.
bool allele_depths_call_het(const CandidateVariant& var, const Options& opts) {
    if (opts.retry_windows.empty() && !opts.joint_het_orientation) return false;
    if (var.lcd_var_i_to_cate != kCandNoisyCandHet &&
        var.lcd_var_i_to_cate != kCandCleanHetSnp &&
        var.lcd_var_i_to_cate != kCandCleanHetIndel) return false;
    if (!var.msa_insertion_alts.empty()) return false;  // handled jointly above
    if (var.hap_to_alle_profile[1].size() < 2 || var.hap_to_alle_profile[2].size() < 2)
        return false;
    if (var.counts.ref_cov < opts.min_alt_depth || var.counts.alt_cov < opts.min_alt_depth)
        return false;
    if (var.counts.allele_fraction < opts.min_af || var.counts.allele_fraction > opts.max_af)
        return false;
    if (opts.joint_het_orientation) return true;
    if (!var.alignment_verified) return false;
    const hts_pos_t pos = var.key.sort_pos();
    bool inside_window = false;
    for (const auto& [beg, end] : opts.retry_windows)
        if (pos >= beg && pos < end) {
            inside_window = true;
            break;
        }
    if (!inside_window) return false;
    // Repeat-derived rows can affect connectivity only after their merged BAM
    // observations meet the same net-margin rule as the graph block linker.
    return (!var.is_homopolymer_indel && var.msa_insertion_alts.empty()) ||
           var.gap_link_supported;
}

static int iter_update_var_hap_to_cons_alle(PhasingChunk& chunk, bool is_ont,
                                             const std::vector<int>& valid_var_idx,
                                             uint32_t flags, const Options& opts,
                                             const std::vector<char>* pinned = nullptr) {
    const int n = (int)valid_var_idx.size();

    // Save current consensus for convergence check.
    std::vector<std::array<int, 3>> saved(n);
    for (int _vi = 0; _vi < n; ++_vi)
        saved[_vi] = chunk.candidates[valid_var_idx[_vi]].hap_to_cons_alle;

    var_init_hap_to_alle_profile(chunk.candidates, valid_var_idx);

    for (int oi = 0; oi < read_visit_count(chunk); ++oi) {
        const int read_i = read_at_visit_ord(chunk, oi);
        if (read_i < 0 || static_cast<size_t>(read_i) >= chunk.reads.size()) continue;
        if (chunk.reads[read_i].is_skipped) continue;
        // Apply phase-set-local updates to clean-candidate rounds first.
        // MSA-round scoping exposes unresolved repeat-link regressions;
        // retain its existing update path (see CHECKPOINT.md).
        if (opts.upstream_assign_hap || (flags & kCandNoisyCandHet) ||
            !opts.phase_set_scoped_clean_rounds) {
            int hap = init_assign_read_hap_based_on_cons_alle(chunk, read_i, flags, std::nullopt,
                                                    opts.msa_sites_vote_without_gap_link,
                                       opts.infer_complement_at_multiallelic,
                                       opts.upstream_read_scoring);
            if (hap == -1) hap = 0;
            chunk.haps[read_i] = hap;
            update_var_hap_profile_based_on_read_hap(chunk, read_i, hap, flags);
            continue;
        }
        const ReadVariantProfile& prof = chunk.read_var_profile[read_i];
        if (prof.start_var_idx < 0) continue;
        std::vector<hts_pos_t> phase_sets;
        for (int vi = prof.start_var_idx; vi <= prof.end_var_idx; ++vi) {
            const CandidateVariant& var = chunk.candidates[vi];
            if ((var.lcd_var_i_to_cate & flags) == 0 ||
                prof.alleles[vi - prof.start_var_idx] < 0) continue;
            phase_sets.push_back(var.phase_set);
        }
        std::sort(phase_sets.begin(), phase_sets.end());
        phase_sets.erase(std::unique(phase_sets.begin(), phase_sets.end()), phase_sets.end());
        // HP integers in disconnected blocks have independent orientations.
        // A spanning read must update each block using that block's evidence.
        for (const hts_pos_t phase_set : phase_sets) {
            const int hap = std::max(0, init_assign_read_hap_based_on_cons_alle(chunk, read_i, flags, phase_set,
                                       opts.msa_sites_vote_without_gap_link,
                                       opts.infer_complement_at_multiallelic,
                                       opts.upstream_read_scoring));
            update_var_hap_profile_based_on_read_hap(chunk, read_i, hap, flags, phase_set);
        }
        chunk.haps[read_i] = std::max(0, init_assign_read_hap_based_on_cons_alle(chunk, read_i, flags,
            phase_sets.empty() ? std::nullopt : std::optional<hts_pos_t>(phase_sets.front())));
    }

    for (int _vi = 0; _vi < n; ++_vi) {
        CandidateVariant& var = chunk.candidates[valid_var_idx[_vi]];
        // Anchored: a pinned site keeps the consensus the previous round gave it.
        // Pinning must SKIP the write rather than revert it afterwards: the
        // convergence check below compares against a snapshot taken at entry, so a
        // revert-after would keep reporting 'changed' and never converge.
        if (pinned != nullptr && (*pinned)[static_cast<size_t>(valid_var_idx[_vi])]) continue;
        if ((var.gap_link_supported || two_allele_het(var)) &&
            var.msa_insertion_alts.size() == 2) {
            const int same = var.hap_to_alle_profile[1][1] + var.hap_to_alle_profile[2][2];
            const int flip = var.hap_to_alle_profile[1][2] + var.hap_to_alle_profile[2][1];
            // Optimize the orientation of the verified allele pair jointly.
            // Independent haplotype majorities can select the same allele twice.
            //
            // The condition reads the CURRENT read evidence, not `saved`. It
            // previously also required saved[1] > 0, saved[2] > 0 and
            // saved[1] != saved[2] -- that is, the pair had to ALREADY be a
            // het in the previous iteration before this branch would orient it.
            // `saved` is the previous iteration's consensus, kept for the
            // convergence check below, and when that iteration had collapsed
            // both haplotypes onto one allele -- the exact failure this branch
            // exists to fix -- the guard refused and the collapse became
            // permanent. It left three of the panel's twenty-seven multiallelic
            // loci emitted 1|1 with reads on both alternates: 55,795,217
            // (AD 0,40,28), 55,815,775 (0,18,10) and 5,379,662 (0,45,31).
            // two_allele_het already establishes from the counts that reads sit
            // behind both alternates, which is the precondition that matters.
            if (same != flip) {
                var.hap_to_cons_alle[1] = same > flip ? 1 : 2;
                var.hap_to_cons_alle[2] = same > flip ? 2 : 1;
            } else {
                // No preference either way: seed a het, for the reason the
                // biallelic branch below gives -- an arbitrary orientation is
                // resolvable by the link votes, a collapsed one is not.
                var.hap_to_cons_alle[1] = 1;
                var.hap_to_cons_alle[2] = 2;
            }
        } else if (allele_depths_call_het(var, opts)) {
            // Same hazard the branch above guards against, for a plain biallelic
            // site: taking each haplotype's majority independently lets both
            // pick the same allele, which emits a genuine het as 1|1 and makes
            // it link nothing. Inside a window the first solve could not phase
            // the reads carry no labels, so both majorities are the deeper
            // allele and the collapse is certain -- and it is self-sustaining,
            // because the window then stays unphasable. Measured on
            // chr20:48,204,383 (AT>A, 30 ref / 41 alt, AF 0.577), the only
            // heterozygote between 48,183,976 and 48,225,786 and the site
            // hiphase bridges this gap with: emitted homozygous, the solve
            // jumps 41.8 kb with no spanning read instead.
            //
            // So orient the pair jointly, exactly as the verified-allele branch
            // does, and when the labels carry no preference at all seed a het
            // rather than a hom: an arbitrary orientation is resolvable by the
            // link votes, a collapsed one is not.
            const int same = var.hap_to_alle_profile[1][0] + var.hap_to_alle_profile[2][1];
            const int flip = var.hap_to_alle_profile[1][1] + var.hap_to_alle_profile[2][0];
            var.hap_to_cons_alle[1] = same >= flip ? 0 : 1;
            var.hap_to_cons_alle[2] = same >= flip ? 1 : 0;
        } else {
            for (int hap = 1; hap <= 2; ++hap)
                update_var_hap_to_cons_alle(is_ont, var, hap);
        }
    }

    int changed = 0;
    for (int _vi = 0; _vi < n && !changed; ++_vi) {
        const auto& cur = chunk.candidates[valid_var_idx[_vi]].hap_to_cons_alle;
        for (int hap = 1; hap <= 2; ++hap)
            if (cur[hap] != saved[_vi][hap]) { changed = 1; break; }
    }
    return changed;
}

// Assign per-read phase_sets from the first phased het var each read covers.
// Update_read_phase_set.
static void update_read_phase_set(PhasingChunk& chunk, const std::vector<bool>& var_is_valid,
                                 bool upstream_read_scoring) {
    for (int oi = 0; oi < read_visit_count(chunk); ++oi) {
        const int read_i = read_at_visit_ord(chunk, oi);
        if (read_i < 0 || static_cast<size_t>(read_i) >= chunk.reads.size()) continue;
        if (chunk.reads[read_i].is_skipped) continue;
        const ReadVariantProfile& prof = chunk.read_var_profile[read_i];
        if (prof.start_var_idx < 0) continue;
        hts_pos_t ps = -1;
        for (int vi = prof.start_var_idx; vi <= prof.end_var_idx; ++vi) {
            if (!var_is_valid[vi]) continue;
            const CandidateVariant& var = chunk.candidates[vi];
            // Use the same eligible evidence as init_assign_read_hap_based_on_cons_alle. An
            // excluded repeat can inherit a preceding PS without a link.
            if (!upstream_read_scoring) {
                if (var.is_homopolymer_indel || var.lcd_var_i_to_cate == kCandNoisyCandHom ||
                    (!var.msa_insertion_alts.empty() && !var.gap_link_supported)) continue;
                const int allele = prof.alleles[vi - prof.start_var_idx];
                if (allele < 0 || (allele != var.hap_to_cons_alle[1] &&
                                   allele != var.hap_to_cons_alle[2])) continue;
            }
            if (var.hap_to_cons_alle[1] != -1 && var.hap_to_cons_alle[2] != -1 &&
                var.hap_to_cons_alle[1] != var.hap_to_cons_alle[2]) {
                ps = var.phase_set;
            }
            if (ps != -1) break;
        }
        chunk.phase_sets[read_i] = ps;
    }
}

// ════════════════════════════════════════════════════════════════════════════
// Public entry point
// ════════════════════════════════════════════════════════════════════════════


bool read_carries_phase_tags(const int mapq, const Options& opts) {
    return mapq >= opts.min_assign_mapq;
}

/// The solver: cluster the chunk's reads into two haplotypes over the
/// candidates whose category matches `flags`, then assign each read a haplotype
/// and a phase set.
///
/// `flags` is the whole difference between the two solves in a retry:
///   kCandGermlineClean   -- clean het SNPs, clean het indels, clean hom. The
///                           first solve, and the only one the default runs.
///   kCandGermlineVarCate -- that mask PLUS noisy het and noisy hom. What
///                           collect_noisy_vars_step4 uses once
///                           skip_noisy_kmeans is clear.
///
/// It is a GLOBAL clustering with NO pairwise evidence test, which is the
/// property to keep in mind when widening `flags`: any admitted site can act as
/// a bridge between two components, so one unreliable site can join two blocks
/// on evidence no pairwise link check would accept. Measured at
/// chr20:24,105,188, where the chain evidence reads 53/0, 10/1, 9/5, 0/0: a
/// single near-even repeat-tract site supplied the path, and the entire right
/// flank -- 13 sites each segregating 1.000 against read truth -- was absorbed
/// INVERTED. The block linker is not at fault there and cannot be: it correctly
/// starts a new phase set for a het with no supported link, but by then the
/// k-means has already merged the components.
void assign_hap_based_on_germline_het_vars_kmeans(PhasingChunk& chunk,
                                                   const Options& opts,
                                                   uint32_t flags,
                                                   bool anchored) {
    const int n_cands = (int)chunk.candidates.size();
    std::vector<int> valid_var_idx;
    valid_var_idx.reserve(n_cands);
    std::vector<bool> var_is_valid(n_cands, false);
    for (int i = 0; i < n_cands; ++i) {
        if ((chunk.candidates[i].lcd_var_i_to_cate & flags) == 0) continue;
        valid_var_idx.push_back(i);
        var_is_valid[i] = true;
    }
    if (valid_var_idx.empty()) return;

    // Debug: dump the read x variant allele matrix before phasing mutates state.
    dump_phase_matrix(chunk, valid_var_idx, var_is_valid, opts, flags, nullptr);

    const bool is_ont = opts.is_ont();
    const size_t n_reads = chunk.reads.size();

    // Anchored: keep the incoming read labels and the consensuses already
    // decided, so this round refines the previous one instead of replacing it.
    // The pin list is taken BEFORE the init call, which is what decides the
    // gauge; preserving read labels alone anchors nothing, because Phase 3 below
    // recomputes every read's haplotype from the consensus.
    std::vector<char> pinned;
    bool any_pinned = false;
    if (anchored) {
        pinned.assign(chunk.candidates.size(), 0);
        for (int vi : valid_var_idx) {
            const auto& cons = chunk.candidates[vi].hap_to_cons_alle;
            if (cons[1] != -1 || cons[2] != -1) pinned[static_cast<size_t>(vi)] = 1;
        }
        any_pinned = true;
    } else {
        chunk.haps.assign(n_reads, 0);
        chunk.phase_sets.assign(n_reads, kUnphasedReadPhaseSet);
        read_init_hap_phase_set(chunk);
    }
    var_init_hap_profile_cons_allele(is_ont, chunk.candidates, valid_var_idx, anchored);

    // Phase 1: initial sweep from highest-confidence pivot variant outward.
    // Phase 1 sweeps outward from a pivot chosen over THIS round's site set, and
    // that is how a noisy site can carry a parity across the chunk. Anchored,
    // the incoming consensuses are the gauge, so there is nothing to seed.
    const int init_vi = anchored ? -1 : select_init_var(chunk.candidates, valid_var_idx);
    if (init_vi != -1) {
        const int nv = (int)valid_var_idx.size();

        // Build sweep order: [pivot, pivot-1, ..., 0, pivot+1, ..., nv-1].
        std::vector<int> var_ii(nv);
        var_ii[0] = init_vi;
        for (int vi = init_vi - 1; vi >= 0; --vi) var_ii[init_vi - vi] = vi;
        for (int vi = init_vi + 1; vi < nv; ++vi) var_ii[vi] = vi;

        int64_t* ovlp_b = nullptr;
        int64_t max_b = 0;
        cgranges_t* cr = chunk.read_var_cr.get();

        for (int idx = 0; idx < nv; ++idx) {
            const int vi = valid_var_idx[var_ii[idx]];
            const uint32_t vic = chunk.candidates[vi].lcd_var_i_to_cate;
            // HOM variants are not used in the first round.
            if (vic == kCandCleanHom || vic == kCandNoisyCandHom) continue;

            const int64_t ovlp_n = cr_overlap(cr, "cr", vi, vi + 1, &ovlp_b, &max_b);
            for (int64_t oi = 0; oi < ovlp_n; ++oi) {
                const int read_i = (int)cr_label(cr, ovlp_b[oi]);
                if (chunk.reads[read_i].is_skipped || chunk.haps[read_i] != 0) continue;
                int hap = init_assign_read_hap_based_on_cons_alle(chunk, read_i, flags, std::nullopt,
                                                    opts.msa_sites_vote_without_gap_link,
                                       opts.infer_complement_at_multiallelic,
                                       opts.upstream_read_scoring);
                if (hap == -1) hap = 1; // no informative vars yet — seed new phase set as hap1
                chunk.haps[read_i] = hap;
                update_var_hap_profile_cons_alle_based_on_read_hap(chunk, is_ont, read_i, hap, flags);
            }
        }
        free(ovlp_b);
    }

    // Phase 2: iterative k-means (up to 10 rounds, stop on convergence).
    for (int iter = 0; iter < 10; ++iter) {
        const int c1 = iter_update_var_hap_cons_phase_set(chunk, valid_var_idx, opts);
        const int c2 = iter_update_var_hap_to_cons_alle(chunk, is_ont, valid_var_idx, flags, opts,
                                        any_pinned && !pinned.empty() ? &pinned : nullptr);
        if (c1 == 0 && c2 == 0) break;
    }

    // Phase 3: report HP in its own phase set, using the final consensus.
    update_read_phase_set(chunk, var_is_valid, opts.upstream_read_scoring);
    // assign_hap.c returns the last iterative read labels after setting PS.
    if (!opts.upstream_assign_hap) {
        for (size_t ri = 0; ri < n_reads; ++ri) {
            if (chunk.reads[ri].is_skipped) continue;
            chunk.haps[ri] = chunk.phase_sets[ri] < 0 ? 0 :
                std::max(0, init_assign_read_hap_based_on_cons_alle(chunk, static_cast<int>(ri), flags,
                           chunk.phase_sets[ri], opts.msa_sites_vote_without_gap_link,
                                           opts.infer_complement_at_multiallelic,
                                           opts.upstream_read_scoring));
        }
    }

    // Phase 4: fill hap_alt / hap_ref from finalized hap_to_cons_alle.
    // Resolve hap_alt/hap_ref from finalized consensus alleles:
    //   both -1 -> fall back to hap_to_cons_alle[0] (hom_idx);
    //   one -1  -> treat as ref (0);
    //   any non-zero allele index is treated as ALT for GT projection.
    for (auto& var : chunk.candidates) {
        int c1 = var.hap_to_cons_alle[1];
        int c2 = var.hap_to_cons_alle[2];
        var.hap_alt = 0;
        var.hap_ref = 0;
        if (c1 == -1 && c2 == -1) {
            c1 = c2 = var.hap_to_cons_alle[0]; // hom_idx fallback
        }
        if (c1 == -1) c1 = 0; // unknown hap → ref
        if (c2 == -1) c2 = 0;
        const bool h1_alt = (c1 != 0);
        const bool h2_alt = (c2 != 0);
        if (h1_alt && h2_alt) {
            var.hap_alt = 3; var.hap_ref = 0;
        } else if (h1_alt && !h2_alt) {
            var.hap_alt = 1; var.hap_ref = 2;
        } else if (!h1_alt && h2_alt) {
            var.hap_alt = 2; var.hap_ref = 1;
        }
        // both ref/unresolved: leave hap_alt/hap_ref as 0/0

        if (opts.verbose >= 2 && !opts.debug_site.empty()) {
            hts_pos_t dbg = 0;
            if (parse_debug_site_pos(opts.debug_site, dbg)) {
                hts_pos_t vcf_pos = var.key.pos;
                if (var.key.type == VariantType::Insertion || var.key.type == VariantType::Deletion) {
                    vcf_pos = std::max<hts_pos_t>(1, var.key.pos - 1);
                }
                if (vcf_pos == dbg) {
                    const char t = (var.key.type == VariantType::Snp ? 'X'
                                 : (var.key.type == VariantType::Insertion ? 'I' : 'D'));
                    std::fprintf(stderr,
                                 "DebugPhase\tpos=%" PRId64 "\tkey_pos=%" PRId64 "\ttype=%c\talt=%s\tc1=%d\tc2=%d\thap_alt=%d\thap_ref=%d\thp=%d\n",
                                 static_cast<int64_t>(vcf_pos),
                                 static_cast<int64_t>(var.key.pos),
                                 t,
                                 var.key.alt.c_str(),
                                 var.hap_to_cons_alle[1],
                                 var.hap_to_cons_alle[2],
                                 var.hap_alt,
                                 var.hap_ref,
                                 var.is_homopolymer_indel ? 1 : 0);
                }
            }
        }
    }
}

// ════════════════════════════════════════════════════════════════════════════
// Chunk-boundary stitching: flip downstream hap labels when overlap reads disagree.
// ════════════════════════════════════════════════════════════════════════════

// Stitch one pair of adjacent chunks: optionally flip hap labels in the
// downstream chunk, then merge its phase set into the upstream chunk's PS.
// longcallD updates read HP/PS only when it writes a phased alignment.
static void apply_chunk_flip_and_merge(PhasingChunk& cur,
                                       bool do_flip,
                                       hts_pos_t max_pre_ps,
                                       hts_pos_t min_cur_ps,
                                       bool update_reads) {
    // Flip hap labels when overlap reads voted for a flip and a valid PS exists.
    if (do_flip && min_cur_ps != INT64_MAX && min_cur_ps != static_cast<hts_pos_t>(-1)) {
        for (CandidateVariant& v : cur.candidates) {
            if (v.phase_set != min_cur_ps) continue;
            std::swap(v.hap_to_cons_alle[1], v.hap_to_cons_alle[2]);
        }
        if (update_reads) {
            for (size_t read_i = 0; read_i < cur.reads.size(); ++read_i) {
                if (cur.reads[read_i].is_skipped || cur.haps[read_i] == 0) continue;
                if (cur.phase_sets[read_i] == min_cur_ps) cur.haps[read_i] = 3 - cur.haps[read_i];
            }
        }
    }
    if (max_pre_ps != -1 && min_cur_ps != INT64_MAX) {
        for (CandidateVariant& v : cur.candidates) {
            if (v.phase_set == -1) continue;
            if (v.phase_set == min_cur_ps) v.phase_set = max_pre_ps;
        }
        if (update_reads) {
            for (size_t read_i = 0; read_i < cur.reads.size(); ++read_i) {
                if (cur.phase_sets[read_i] == kUnphasedReadPhaseSet) continue;
                if (cur.phase_sets[read_i] == min_cur_ps) cur.phase_sets[read_i] = max_pre_ps;
            }
        }
    }
}

bool select_stitch_orientation(const std::array<int, 4>& votes,
                               const Options* opts, bool& do_flip) {
    const int n11 = votes[0], n12 = votes[1], n21 = votes[2], n22 = votes[3];
    const int flip_hap_score = n12 + n21 - n11 - n22;
    const int margin = (opts != nullptr) ? opts->stitch_min_margin : 0;
    const int rule = (opts != nullptr) ? opts->stitch_rule : kStitchRuleNetMargin;

    // Decide which orientation (if any) to merge under the selected rule.
    // do_flip is only meaningful when merge == true.
    bool merge = false;
    do_flip = false;
    switch (rule) {
        case kStitchRuleBothStrands: {
            // Reading A: merge only when the winning orientation has evidence on
            // BOTH of its haplotype links.  no-flip needs n11>=1 AND n22>=1;
            // flip needs n12>=1 AND n21>=1.  When both orientations qualify,
            // pick the one with more net votes (ties stay unmerged).
            const bool noflip_ok = (n11 >= 1 && n22 >= 1);
            const bool flip_ok = (n12 >= 1 && n21 >= 1);
            const int noflip_votes = n11 + n22;
            const int flip_votes = n12 + n21;
            if (noflip_ok && flip_ok) {
                if (flip_votes > noflip_votes) { merge = true; do_flip = true; }
                else if (noflip_votes > flip_votes) { merge = true; do_flip = false; }
            } else if (noflip_ok) {
                merge = true; do_flip = false;
            } else if (flip_ok) {
                merge = true; do_flip = true;
            }
            break;
        }
        case kStitchRuleLiteral: {
            // Reading B: merge whenever both orientations have >=1 supporting
            // read (>=2 reads total, at least one each way).  Diagnostic: this
            // deliberately stitches contested seams.  Orientation follows the
            // net vote.
            const int noflip_votes = n11 + n22;
            const int flip_votes = n12 + n21;
            if (noflip_votes >= 1 && flip_votes >= 1) {
                merge = true;
                do_flip = flip_hap_score > 0;
            }
            break;
        }
        case kStitchRuleBothStrandsMargin: {
            // Reading C: Reading A AND the net vote still exceeds the margin.
            const bool noflip_ok = (n11 >= 1 && n22 >= 1);
            const bool flip_ok = (n12 >= 1 && n21 >= 1);
            if (std::abs(flip_hap_score) > margin) {
                if (flip_hap_score > 0 && flip_ok) { merge = true; do_flip = true; }
                else if (flip_hap_score < 0 && noflip_ok) { merge = true; do_flip = false; }
            }
            break;
        }
        case kStitchRuleNetMargin:
        default: {
            // Abstain on weakly supported boundaries: require the net vote
            // magnitude to strictly exceed the configured margin.  margin 0
            // reproduces the original behavior (merge on any non-zero score).
            if (std::abs(flip_hap_score) > margin) {
                merge = true;
                do_flip = flip_hap_score > 0;
            }
            break;
        }
    }

    return merge;
}

// Overlap-read voting between adjacent chunks: count reads that agree vs
// disagree on hap assignment, then flip + merge if disagreement wins.
static bool flip_chunk_hap(PhasingChunk& pre, PhasingChunk& cur,
                           const Options* opts, bool update_reads) {
    if (pre.region.tid != cur.region.tid) return false;

    int n_cur_ovlp_reads = 0;
    int n_pre_ovlp_reads = 0;
    const size_t n_bams = cur.up_ovlp_read_i.size();
    for (size_t bi = 0; bi < n_bams; ++bi) {
        n_cur_ovlp_reads += static_cast<int>(cur.up_ovlp_read_i[bi].size());
        if (bi < pre.down_ovlp_read_i.size()) {
            n_pre_ovlp_reads += static_cast<int>(pre.down_ovlp_read_i[bi].size());
        }
    }
    if (n_cur_ovlp_reads != n_pre_ovlp_reads) {
        throw std::runtime_error("overlap read count mismatch between adjacent chunks");
    }
    if (n_cur_ovlp_reads <= 0) return false;
    if (pre.candidates.empty() || cur.candidates.empty()) return false;

    hts_pos_t max_pre_read_ps = -1;
    hts_pos_t min_cur_read_ps = INT64_MAX;

    // Per-haplotype-link evidence counts.  Each overlap read maps a phased
    // upstream hap (1 or 2) to a phased downstream hap (1 or 2).  The four
    // link types group into two orientations:
    //   no-flip: pre1->cur1 (n11) and pre2->cur2 (n22)
    //   flip:    pre1->cur2 (n12) and pre2->cur1 (n21)
    int n11 = 0, n12 = 0, n21 = 0, n22 = 0;

    for (size_t bi = 0; bi < n_bams; ++bi) {
        const std::vector<int>& cur_list = cur.up_ovlp_read_i[bi];
        const std::vector<int>& pre_list = pre.down_ovlp_read_i[bi];
        for (size_t j = 0; j < cur_list.size(); ++j) {
            if (j >= pre_list.size()) {
                throw std::runtime_error("overlap read pairing mismatch between adjacent chunks");
            }
            const int cur_read_i = cur_list[j];
            const int pre_read_i = pre_list[j];
            if (pre_read_i < 0 || static_cast<size_t>(pre_read_i) >= pre.reads.size() ||
                cur_read_i < 0 || static_cast<size_t>(cur_read_i) >= cur.reads.size()) {
                throw std::runtime_error("overlap read index out of bounds during chunk stitching");
            }
            if (pre.reads[static_cast<size_t>(pre_read_i)].is_skipped ||
                pre.haps[static_cast<size_t>(pre_read_i)] == 0 ||
                cur.reads[static_cast<size_t>(cur_read_i)].is_skipped ||
                cur.haps[static_cast<size_t>(cur_read_i)] == 0) {
                continue;
            }
            const int pre_read_hap = pre.haps[static_cast<size_t>(pre_read_i)];
            const hts_pos_t pre_read_ps = pre.phase_sets[static_cast<size_t>(pre_read_i)];
            const int cur_read_hap = cur.haps[static_cast<size_t>(cur_read_i)];
            const hts_pos_t cur_read_ps = cur.phase_sets[static_cast<size_t>(cur_read_i)];
            if (pre_read_hap == 1 && cur_read_hap == 1) ++n11;
            else if (pre_read_hap == 1 && cur_read_hap == 2) ++n12;
            else if (pre_read_hap == 2 && cur_read_hap == 1) ++n21;
            else ++n22;
            if (max_pre_read_ps < pre_read_ps) max_pre_read_ps = pre_read_ps;
            if (min_cur_read_ps > cur_read_ps) min_cur_read_ps = cur_read_ps;
        }
    }

    bool do_flip = false;
    const bool merge = select_stitch_orientation({n11, n12, n21, n22}, opts, do_flip);

    if (!merge) return false;

    apply_chunk_flip_and_merge(cur,
                               do_flip,
                               max_pre_read_ps,
                               min_cur_read_ps,
                               update_reads);
    return true;
}

// Copy downstream HP/PS onto unphased upstream overlap reads.
// Only safe when the chunk pair has been merged (flip_chunk_hap returned true),
// because the merge rewrites the downstream PS to match the upstream PS and
// flips hap labels if needed.  For unmerged pairs the relative phase is unknown
// and propagation could assign reads to the wrong haplotype.
static void propagate_overlap_read_phase_to_output_owner(PhasingChunk& pre, const PhasingChunk& cur) {
    if (pre.region.tid != cur.region.tid) return;
    const size_t n_bams = std::min(pre.down_ovlp_read_i.size(), cur.up_ovlp_read_i.size());
    for (size_t bi = 0; bi < n_bams; ++bi) {
        const std::vector<int>& pre_list = pre.down_ovlp_read_i[bi];
        const std::vector<int>& cur_list = cur.up_ovlp_read_i[bi];
        if (pre_list.size() != cur_list.size()) {
            throw std::runtime_error("overlap read count mismatch between adjacent chunks");
        }
        for (size_t j = 0; j < pre_list.size(); ++j) {
            const int pre_read_i = pre_list[j];
            const int cur_read_i = cur_list[j];
            if (pre_read_i < 0 || static_cast<size_t>(pre_read_i) >= pre.reads.size() ||
                cur_read_i < 0 || static_cast<size_t>(cur_read_i) >= cur.reads.size()) {
                throw std::runtime_error("overlap read index out of bounds during output phase propagation");
            }
            if (pre.reads[static_cast<size_t>(pre_read_i)].is_skipped ||
                cur.reads[static_cast<size_t>(cur_read_i)].is_skipped) {
                continue;
            }
            if (pre.haps[static_cast<size_t>(pre_read_i)] != 0) continue;
            const int cur_hap = cur.haps[static_cast<size_t>(cur_read_i)];
            const hts_pos_t cur_ps = cur.phase_sets[static_cast<size_t>(cur_read_i)];
            if (cur_hap == 0 || cur_ps <= 0) continue;
            pre.haps[static_cast<size_t>(pre_read_i)] = cur_hap;
            pre.phase_sets[static_cast<size_t>(pre_read_i)] = cur_ps;
        }
    }
}

void verify_chunk_invariants(const PhasingChunk& chunk,
                             size_t site_ids_size,
                             size_t site_meta_size,
                             size_t site_orig_size,
                             hts_pos_t region_lo,
                             hts_pos_t region_hi) {
    const size_t n_cand = chunk.candidates.size();
    const size_t n_reads = chunk.reads.size();
    auto fail = [](const std::string& what) {
        throw std::runtime_error("chunk invariant violated: " + what);
    };
    if (site_ids_size != n_cand)
        fail("site_ids has " + std::to_string(site_ids_size) + " entries for " +
             std::to_string(n_cand) + " candidates");
    if (site_meta_size != n_cand)
        fail("site_meta has " + std::to_string(site_meta_size) + " entries for " +
             std::to_string(n_cand) + " candidates");
    if (site_orig_size != n_cand)
        fail("site_allele_orig_idx has " + std::to_string(site_orig_size) +
             " entries for " + std::to_string(n_cand) + " candidates");
    if (chunk.read_var_profile.size() != n_reads)
        fail("read_var_profile has " + std::to_string(chunk.read_var_profile.size()) +
             " entries for " + std::to_string(n_reads) + " reads");
    if (!chunk.haps.empty() && chunk.haps.size() != n_reads)
        fail("haps has " + std::to_string(chunk.haps.size()) + " entries for " +
             std::to_string(n_reads) + " reads");
    if (!chunk.phase_sets.empty() && chunk.phase_sets.size() != n_reads)
        fail("phase_sets has " + std::to_string(chunk.phase_sets.size()) +
             " entries for " + std::to_string(n_reads) + " reads");
    for (size_t i = 1; i < n_cand; ++i)
        if (chunk.candidates[i - 1].key.sort_pos() > chunk.candidates[i].key.sort_pos())
            fail("candidates are not position-sorted at index " + std::to_string(i));
    for (size_t i = 1; i < n_reads; ++i)
        if (chunk.reads[i - 1].qname > chunk.reads[i].qname)
            fail("reads are not qname-sorted at index " + std::to_string(i) +
                 " (the cross-chunk stitch pairs them with a merge-join)");
    for (size_t i = 0; i < chunk.read_var_profile.size(); ++i)
        if (chunk.read_var_profile[i].read_id != static_cast<int>(i))
            fail("read_var_profile[" + std::to_string(i) + "].read_id is " +
                 std::to_string(chunk.read_var_profile[i].read_id));
    if (region_hi > region_lo) {
        if (region_lo < chunk.ref_beg || region_hi > chunk.ref_end)
            fail("re-solved region " + std::to_string(region_lo) + "-" +
                 std::to_string(region_hi) + " reaches outside the chunk " +
                 std::to_string(chunk.ref_beg) + "-" + std::to_string(chunk.ref_end));
    }
}

void stitch_chunk_haps(std::vector<PhasingChunk>& chunks,
                       const Options* opts,
                       const PgbamSidecarData* pgbam_sidecar) {
    const bool use_pgbam = opts != nullptr && pgbam_sidecar != nullptr && !opts->pgbam_file.empty();
    if (use_pgbam) {
        for (PhasingChunk& chunk : chunks) {
            stitch_phase_blocks_with_pgbam(chunk,
                                          *pgbam_sidecar,
                                          opts->pgbam_primary_min_winning_threads,
                                          opts->pgbam_primary_polarity_margin);
        }
    }
    // Track which adjacent pairs were successfully merged so we can safely
    // propagate overlap-read phase only for those pairs (see CHECKPOINT.md
    // "BAM Pipeline Parity").
    std::vector<bool> pair_stitched(chunks.size(), false);
    // collect_var.c updates stitched read HP/PS only when writing an alignment.
    const bool update_reads = opts == nullptr || !opts->upstream_assign_hap ||
                              !opts->output_aln.empty();
    for (size_t ii = 1; ii < chunks.size(); ++ii) {
        const bool stitched = flip_chunk_hap(chunks[ii - 1], chunks[ii], opts,
                                             update_reads);
        if (stitched) {
            pair_stitched[ii] = true;
        } else if (use_pgbam) {
            stitch_adjacent_chunks_with_pgbam(chunks[ii - 1],
                                             chunks[ii],
                                             *pgbam_sidecar,
                                             opts->pgbam_primary_min_winning_threads,
                                             opts->pgbam_primary_polarity_margin);
        }
    }
    if (use_pgbam && opts->pgbam_cleanup_pass) {
        stitch_phase_blocks_with_pgbam(chunks,
                                      *pgbam_sidecar,
                                      opts->pgbam_cleanup_min_winning_threads,
                                      opts->pgbam_cleanup_polarity_margin);
    }
    if (use_pgbam && opts->pgbam_relaxed_cleanup_pass) {
        stitch_phase_blocks_with_pgbam(chunks,
                                      *pgbam_sidecar,
                                      opts->pgbam_relaxed_cleanup_min_winning_threads,
                                      opts->pgbam_relaxed_cleanup_polarity_margin);
    }
    // Propagate overlap-read phase only when read tags are being written and
    // the pair was merged. For an unmerged pair the relative phase is unknown.
    for (size_t ii = chunks.size(); update_reads && ii > 1; --ii) {
        if (pair_stitched[ii - 1]) {
            propagate_overlap_read_phase_to_output_owner(chunks[ii - 2], chunks[ii - 1]);
        }
    }
}

} // namespace pgphase_collect
