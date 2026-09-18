/**
 * @file collect_phase.cpp
 * k-means read-haplotype clustering for diploid phasing.
 */

#include "collect_phase.hpp"
#include "collect_phase_pgbam.hpp"

#include <algorithm>
#include <array>
#include <cinttypes>
#include <climits>
#include <cstdio>
#include <cstdlib>
#include <map>
#include <optional>
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
        chunk.phase_sets[i] = -1;
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
static void var_init_hap_profile_cons_allele(bool is_ont,
                                              CandidateTable& variants,
                                              const std::vector<int>& valid_var_idx) {
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
static void update_var_hap_to_cons_alle(bool is_ont, CandidateVariant& var, int hap) {
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
static int read_to_cons_allele_score(CandidateVariant& var, int hap, int allele_i) {
    const uint32_t var_i_to_cate = var.lcd_var_i_to_cate;
    if (!var.msa_insertion_alts.empty() && !var.gap_link_supported) return 0;
    int var_score = 1;
    if (var.counts.n_uniq_alles <= 2 && var_i_to_cate == kCandCleanHetSnp) var_score = 2;
    else if (var.counts.n_uniq_alles <= 2 && var_i_to_cate == kCandCleanHetIndel) var_score = 2;
    if (var.hap_to_cons_alle[hap] == -1 && var.hap_to_cons_alle[3 - hap] == -1) return 0;
    // A multiallelic site has no unique complementary allele. Infer it only
    // for a biallelic site; otherwise let actual read observations resolve it.
    if (variant_allele_slots(var) == 2) {
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
                              const Options& opts, uint32_t flags) {
    if (opts.phase_matrix_dump_prefix.empty()) return;

    std::string path = opts.phase_matrix_dump_prefix + ".chunk" +
                       std::to_string(chunk.region.chunk_id) + ".flags" +
                       std::to_string(flags) + ".tsv";
    std::FILE* fp = std::fopen(path.c_str(), "w");
    if (fp == nullptr) {
        std::fprintf(stderr, "[dump-phase-matrix] cannot open %s\n", path.c_str());
        return;
    }

    // Variant header block: idx, pos, type, category, weight. Variant idx is the
    // position within valid_var_idx (0-based, contig order).
    std::fprintf(fp, "#tid\t%d\tchunk\t%d\tflags\t%u\tn_vars\t%zu\n",
                 chunk.region.tid, chunk.region.chunk_id, flags, valid_var_idx.size());
    std::fprintf(fp, "#VAR\tvar_idx\tpos\ttype\tcate\tweight\n");
    std::vector<int> global_to_vidx(chunk.candidates.size(), -1);
    for (int vidx = 0; vidx < (int)valid_var_idx.size(); ++vidx) {
        const int gi = valid_var_idx[vidx];
        global_to_vidx[gi] = vidx;
        const CandidateVariant& var = chunk.candidates[gi];
        const char t = (var.key.type == VariantType::Snp ? 'X'
                     : (var.key.type == VariantType::Insertion ? 'I' : 'D'));
        std::fprintf(fp, "VAR\t%d\t%" PRId64 "\t%c\t%u\t%d\n",
                     vidx, static_cast<int64_t>(var.key.pos), t,
                     var.lcd_var_i_to_cate,
                     phase_matrix_var_weight(var));
        std::fprintf(fp, "#META\t%d\tmsa_verified=%d\thomopolymer=%d\tgap_link_supported=%d\n",
                     vidx, var.msa_verified, var.is_homopolymer_indel, var.gap_link_supported);
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

// Assign a read to hap 1, 2, 0 (tied), or -1 (no informative variants).
// Updates n_clean_agree_snps / n_clean_conflict_snps on the read (for max_hap only).
// CleanHom variants contribute to agree/conflict stats but not to hap_scores.
// Init_assign_read_hap_based_on_cons_alle.
static int init_assign_read_hap(PhasingChunk& chunk, int read_i, uint32_t flags,
                                std::optional<hts_pos_t> phase_set = std::nullopt) {
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
        if ((var.is_homopolymer_indel && !var.hp_gap_scorable) ||
            vic == kCandNoisyCandHom) continue;

        const int aidx = prof.alleles[vi - prof.start_var_idx];
        if (aidx < 0) continue;

        for (int hap = 1; hap <= 2; ++hap) {
            const int score = read_to_cons_allele_score(var, hap, aidx);
            if (score != 0) {
                if (vic != kCandCleanHom) n_vars_used[hap]++;
                if (vic == kCandCleanHetSnp && var.counts.n_uniq_alles <= 2) {
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
static void update_var_hap_profile_cons_alle(PhasingChunk& chunk, bool is_ont,
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
static void update_var_hap_profile(PhasingChunk& chunk, int read_i, int hap, uint32_t flags,
                                   std::optional<hts_pos_t> phase_set = std::nullopt) {
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
static int check_agree_haps(const PhasingChunk& chunk, int read_i, int hap, int var1, int var2) {
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
static int check_agree_alleles(const PhasingChunk& chunk, int read_i, int var1, int var2) {
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
        // unreliable for linking. But the emit loop below still hands such a
        // site the running phase set while leaving `parity` unapplied, so it
        // joins a block carrying whatever orientation its own consensus
        // produced, never reconciled against that block. On
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
            (var.msa_insertion_alts.empty() || var.gap_link_supported ||
             two_allele_het(var)) &&
            !hp_indel_blocks_link) {
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
    // Ordinary rounds retain the nearest sufficient link within the preceding
    // `block_link_window` hets. Recovery rounds retain all sufficient edges,
    // so a later verified site can connect earlier components.
    const int window = std::max(1, opts.block_link_window);
    struct Edge { int left, right, agree, conflict; };
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
            const int support = a == c ? 0 : std::max(a, c);
            // Require a net margin for additional repeat links. Read-level
            // disagreements alone do not establish a wrong block orientation.
            if ((chunk.candidates[vi].is_homopolymer_indel ||
                 chunk.candidates[vj].is_homopolymer_indel) &&
                std::abs(a - c) < opts.min_block_link_reads) continue;
            if (support > best_support) {
                best_support = support; best_h = hj; best_a = a; best_c = c;
            }
            if (support >= opts.min_block_link_reads) break;
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
        const int support = link_agree[hi] == link_conflict[hi]
                                ? 0 : std::max(link_agree[hi], link_conflict[hi]);
        if (hj < 0 || support < opts.min_block_link_reads) {
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
/// because its allele is unreliable for linking, but the emit loop still hands
/// such a site the running phase set while leaving `parity` unapplied -- so it
/// joins a block carrying its own unreconciled orientation. On
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
    const hts_pos_t pos = var.key.sort_pos();
    for (const auto& [beg, end] : opts.retry_windows)
        if (pos >= beg && pos < end) return true;
    return false;
}

static int iter_update_var_hap_to_cons_alle(PhasingChunk& chunk, bool is_ont,
                                             const std::vector<int>& valid_var_idx,
                                             uint32_t flags, const Options& opts) {
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
        if (flags & kCandNoisyCandHet) {
            int hap = init_assign_read_hap(chunk, read_i, flags);
            if (hap == -1) hap = 0;
            chunk.haps[read_i] = hap;
            update_var_hap_profile(chunk, read_i, hap, flags);
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
            const int hap = std::max(0, init_assign_read_hap(chunk, read_i, flags, phase_set));
            update_var_hap_profile(chunk, read_i, hap, flags, phase_set);
        }
        chunk.haps[read_i] = std::max(0, init_assign_read_hap(chunk, read_i, flags,
            phase_sets.empty() ? std::nullopt : std::optional<hts_pos_t>(phase_sets.front())));
    }

    for (int _vi = 0; _vi < n; ++_vi) {
        CandidateVariant& var = chunk.candidates[valid_var_idx[_vi]];
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
static void update_read_phase_set(PhasingChunk& chunk, const std::vector<bool>& var_is_valid) {
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
            // Use the same eligible evidence as init_assign_read_hap. An
            // excluded repeat can inherit a preceding PS without a link.
            if (var.is_homopolymer_indel || var.lcd_var_i_to_cate == kCandNoisyCandHom ||
                (!var.msa_insertion_alts.empty() && !var.gap_link_supported)) continue;
            const int allele = prof.alleles[vi - prof.start_var_idx];
            if (allele < 0 || (allele != var.hap_to_cons_alle[1] &&
                               allele != var.hap_to_cons_alle[2])) continue;
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
                                                   uint32_t flags) {
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
    dump_phase_matrix(chunk, valid_var_idx, var_is_valid, opts, flags);

    const bool is_ont = opts.is_ont();
    const size_t n_reads = chunk.reads.size();

    chunk.haps.assign(n_reads, 0);
    chunk.phase_sets.assign(n_reads, -1);
    read_init_hap_phase_set(chunk);
    var_init_hap_profile_cons_allele(is_ont, chunk.candidates, valid_var_idx);

    // Phase 1: initial sweep from highest-confidence pivot variant outward.
    const int init_vi = select_init_var(chunk.candidates, valid_var_idx);
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
                int hap = init_assign_read_hap(chunk, read_i, flags);
                if (hap == -1) hap = 1; // no informative vars yet — seed new phase set as hap1
                chunk.haps[read_i] = hap;
                update_var_hap_profile_cons_alle(chunk, is_ont, read_i, hap, flags);
            }
        }
        free(ovlp_b);
    }

    // Phase 2: iterative k-means (up to 10 rounds, stop on convergence).
    for (int iter = 0; iter < 10; ++iter) {
        const int c1 = iter_update_var_hap_cons_phase_set(chunk, valid_var_idx, opts);
        const int c2 = iter_update_var_hap_to_cons_alle(chunk, is_ont, valid_var_idx, flags, opts);
        if (c1 == 0 && c2 == 0) break;
    }

    // Phase 3: report HP in its own phase set, using the final consensus.
    update_read_phase_set(chunk, var_is_valid);
    for (size_t ri = 0; ri < n_reads; ++ri) {
        if (chunk.reads[ri].is_skipped) continue;
        chunk.haps[ri] = chunk.phase_sets[ri] < 0 ? 0 :
            std::max(0, init_assign_read_hap(chunk, static_cast<int>(ri), flags, chunk.phase_sets[ri]));
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
// Keeps both candidate and read-level HP/PS in sync.
static void apply_chunk_flip_and_merge(PhasingChunk& cur,
                                       bool do_flip,
                                       hts_pos_t max_pre_ps,
                                       hts_pos_t min_cur_ps) {
    // Flip hap labels when overlap reads voted for a flip and a valid PS exists.
    if (do_flip && min_cur_ps != INT64_MAX && min_cur_ps != static_cast<hts_pos_t>(-1)) {
        for (CandidateVariant& v : cur.candidates) {
            if (v.phase_set != min_cur_ps) continue;
            std::swap(v.hap_to_cons_alle[1], v.hap_to_cons_alle[2]);
        }
        for (size_t read_i = 0; read_i < cur.reads.size(); ++read_i) {
            if (cur.reads[read_i].is_skipped || cur.haps[read_i] == 0) continue;
            if (cur.phase_sets[read_i] == min_cur_ps) cur.haps[read_i] = 3 - cur.haps[read_i];
        }
    }
    if (max_pre_ps != -1 && min_cur_ps != INT64_MAX) {
        for (CandidateVariant& v : cur.candidates) {
            if (v.phase_set == -1) continue;
            if (v.phase_set == min_cur_ps) v.phase_set = max_pre_ps;
        }
        for (size_t read_i = 0; read_i < cur.reads.size(); ++read_i) {
            if (cur.phase_sets[read_i] == -1) continue;
            if (cur.phase_sets[read_i] == min_cur_ps) cur.phase_sets[read_i] = max_pre_ps;
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
static bool flip_chunk_hap(PhasingChunk& pre, PhasingChunk& cur, const Options* opts) {
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
                               min_cur_read_ps);
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
    for (size_t ii = 1; ii < chunks.size(); ++ii) {
        const bool stitched = flip_chunk_hap(chunks[ii - 1], chunks[ii], opts);
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
    // Propagate overlap-read phase from downstream to upstream only for
    // pairs that were successfully merged.  For unmerged pairs the relative
    // phase is unknown and propagation could assign the wrong haplotype.
    for (size_t ii = chunks.size(); ii > 1; --ii) {
        if (pair_stitched[ii - 1]) {
            propagate_overlap_read_phase_to_output_owner(chunks[ii - 2], chunks[ii - 1]);
        }
    }
}

} // namespace pgphase_collect
