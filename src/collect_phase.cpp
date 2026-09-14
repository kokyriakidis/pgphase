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
static int init_assign_read_hap(PhasingChunk& chunk, int read_i, uint32_t flags) {
    ReadRecord& read = chunk.reads[read_i];
    read.n_clean_agree_snps = 0;
    read.n_clean_conflict_snps = 0;
    read.n_bridge_agree_snps = 0;
    read.n_bridge_conflict_snps = 0;

    const ReadVariantProfile& prof = chunk.read_var_profile[read_i];
    if (prof.start_var_idx < 0) return -1;

    int hap_scores[3] = {0, 0, 0};
    int n_vars_used[3] = {0, 0, 0};
    int n_clean_agree[3] = {0, 0, 0};
    int n_clean_conflict[3] = {0, 0, 0};
    int n_bridge_agree[3] = {0, 0, 0};
    int n_bridge_conflict[3] = {0, 0, 0};

    for (int vi = prof.start_var_idx; vi <= prof.end_var_idx; ++vi) {
        CandidateVariant& var = chunk.candidates[vi];
        const uint32_t vic = var.lcd_var_i_to_cate;
        if ((vic & flags) == 0) continue;
        // HP-indel or noisy-hom candidates: score only, do not use as pivot.
        if (var.is_homopolymer_indel || vic == kCandNoisyCandHom) continue;

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
static void update_var_hap_profile(PhasingChunk& chunk, int read_i, int hap, uint32_t flags) {
    const ReadVariantProfile& prof = chunk.read_var_profile[read_i];
    if (prof.start_var_idx < 0) return;
    for (int vi = prof.start_var_idx; vi <= prof.end_var_idx; ++vi) {
        const uint32_t vic = chunk.candidates[vi].lcd_var_i_to_cate;
        if ((vic & flags) == 0) continue;
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

static bool clean_snp_has_confident_bam_observation(const ReadRecord& read,
                                                    const CandidateVariant& var, int allele) {
    constexpr int kGapBridgeMinBaseQuality = 30;
    if (!read.alignment || var.key.ref_len != 1 || var.key.alt.size() != 1 ||
        var.ref_base > 3 || (allele != 0 && allele != 1)) return false;
    const bam1_t* bam = read.alignment.get();
    const auto* cigar = bam_get_cigar(bam);
    hts_pos_t ref_pos = bam->core.pos + 1;
    int query_pos = 0;
    for (uint32_t ci = 0; ci < bam->core.n_cigar; ++ci) {
        const int op = bam_cigar_op(cigar[ci]), len = bam_cigar_oplen(cigar[ci]);
        const int consumption = bam_cigar_type(op);
        if ((consumption & 2) && var.key.pos >= ref_pos && var.key.pos < ref_pos + len) {
            if (!(consumption & 1)) return false;
            const int qi = query_pos + static_cast<int>(var.key.pos - ref_pos);
            const int quality = bam_get_qual(bam)[qi];
            const int expected = allele == 0 ? 1 << var.ref_base :
                seq_nt16_table[static_cast<unsigned char>(var.key.alt[0])];
            return quality != 255 && quality >= kGapBridgeMinBaseQuality &&
                   bam_seqi(bam_get_seq(bam), qi) == expected;
        }
        if (consumption & 1) query_pos += len;
        if (consumption & 2) ref_pos += len;
    }
    return false;
}

// Phase-set assignment + flip for one k-means iteration.
// Returns 1 if any flip occurred (changed), 0 if converged.
// Iter_update_var_hap_cons_phase_set.
int iter_update_var_hap_cons_phase_set(PhasingChunk& chunk,
                                               const std::vector<int>& valid_var_idx,
                                               const Options& opts) {
    const int n = (int)valid_var_idx.size();
    const bool recovery_graph = opts.recover_gaps && opts.link_by_alleles &&
                                opts.private_msa_admit_all_in_region;

    // Collect het var positions (indices into valid_var_idx).
    std::vector<int> het_var_idx;
    std::vector<bool> is_het(n, false);
    int seeded_het = 0;
    for (int _vi = 0; _vi < n; ++_vi) {
        CandidateVariant& var = chunk.candidates[valid_var_idx[_vi]];
        if (recovery_graph && var.gap_link_supported && var.msa_insertion_alts.size() == 2 &&
            var.counts.alle_covs.size() == 3 && var.counts.total_cov > 0 &&
            var.counts.alle_covs[1] >= opts.min_alt_depth &&
            var.counts.alle_covs[2] >= opts.min_alt_depth &&
            static_cast<double>(var.counts.alle_covs[1]) / var.counts.total_cov >= opts.min_af &&
            static_cast<double>(var.counts.alle_covs[2]) / var.counts.total_cov >= opts.min_af &&
            (var.hap_to_cons_alle[1] < 0 || var.hap_to_cons_alle[2] < 0 ||
             var.hap_to_cons_alle[1] == var.hap_to_cons_alle[2])) {
            // Provisional read labels cannot turn two verified alternate
            // sequences into a homozygous site before their links are solved.
            var.hap_to_cons_alle[1] = 1;
            var.hap_to_cons_alle[2] = 2;
            seeded_het = 1;
        }
        const bool gap_hp_link = recovery_graph && opts.gap_hp_link_beg >= 0 &&
                var.key.sort_pos() >= opts.gap_hp_link_beg &&
                var.key.sort_pos() <= opts.gap_hp_link_end &&
                var.lcd_var_i_to_cate == kCandNoisyCandHet && var.gap_link_supported;
        if (gap_hp_link && var.is_homopolymer_indel &&
            var.counts.ref_cov >= opts.min_alt_depth && var.counts.alt_cov >= opts.min_alt_depth &&
            var.counts.allele_fraction >= opts.min_af && var.counts.allele_fraction <= opts.max_af &&
            (var.hap_to_cons_alle[1] < 0 || var.hap_to_cons_alle[2] < 0 ||
             var.hap_to_cons_alle[1] == var.hap_to_cons_alle[2])) {
            // Independent blocks have arbitrary HP labels. Their provisional
            // read assignments must not erase a verified het before its
            // allele edges can resolve the relative block orientation.
            var.hap_to_cons_alle[1] = 0;
            var.hap_to_cons_alle[2] = 1;
            seeded_het = 1;
        }
        if (var.hap_to_cons_alle[1] != -1 && var.hap_to_cons_alle[2] != -1 &&
            var.hap_to_cons_alle[1] != var.hap_to_cons_alle[2] &&
            (var.msa_insertion_alts.empty() || var.gap_link_supported) &&
            (!var.is_homopolymer_indel || gap_hp_link)) {
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
    std::vector<Edge> edges;
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
            const int support = a == c ? 0 : std::max(a, c);
            if (opts.verbose >= 2 && opts.gap_hp_link_beg >= 0 && a + c > 0 &&
                (chunk.candidates[vi].is_homopolymer_indel || chunk.candidates[vj].is_homopolymer_indel))
                std::fprintf(stderr, "GapHpEdge\t%lld\t%lld\t%d\t%d\n",
                             static_cast<long long>(chunk.candidates[vj].key.sort_pos()),
                             static_cast<long long>(chunk.candidates[vi].key.sort_pos()), a, c);
            // Require a net margin for additional repeat links. Read-level
            // disagreements alone do not establish a wrong block orientation.
            if ((chunk.candidates[vi].is_homopolymer_indel ||
                 chunk.candidates[vj].is_homopolymer_indel) &&
                std::abs(a - c) < opts.min_block_link_reads) continue;
            if (recovery_graph && support >= opts.min_block_link_reads)
                edges.push_back({hj, hi, a, c});
            if (support > best_support) {
                best_support = support; best_h = hj; best_a = a; best_c = c;
            }
            if (!recovery_graph && support >= opts.min_block_link_reads) break;
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

    if (recovery_graph) {
        // A later MSA site can connect two earlier components. Keeping just
        // its nearest edge loses that bridge. Process all supported edges,
        // strongest net evidence first, without overturning stronger paths.
        std::sort(edges.begin(), edges.end(), [](const Edge& a, const Edge& b) {
            const int an = std::abs(a.agree - a.conflict), bn = std::abs(b.agree - b.conflict);
            if (an != bn) return an > bn;
            const int as = std::max(a.agree, a.conflict), bs = std::max(b.agree, b.conflict);
            if (as != bs) return as > bs;
            if (a.right - a.left != b.right - b.left) return a.right - a.left < b.right - b.left;
            if (a.right != b.right) return a.right < b.right;
            return a.left < b.left;
        });
        std::vector<int> parent(n_het), flip(n_het, 0);
        for (int hi = 0; hi < n_het; ++hi) parent[hi] = hi;
        auto root = [&](auto&& self, int hi) -> int {
            if (parent[hi] != hi) {
                const int previous = parent[hi];
                parent[hi] = self(self, previous);
                flip[hi] ^= flip[previous];
            }
            return parent[hi];
        };
        for (const auto& edge : edges) {
            const int left = root(root, edge.left), right = root(root, edge.right);
            if (left == right) continue;
            const int child = std::max(left, right);
            parent[child] = std::min(left, right);
            flip[child] = flip[edge.left] ^ flip[edge.right] ^ (edge.conflict > edge.agree);
        }
        // A long read can identify both established blocks through several
        // clean SNPs even when no individual SNP pair has two spanning reads.
        // Evaluate those observations in the components' resolved orientations.
        constexpr int kGapBridgeMinMapq = 30;
        constexpr int kGapBridgeMinAnchorSnps = 2;
        std::vector<int> component(n_het), orientation(n_het);
        for (int hi = 0; hi < n_het; ++hi) {
            component[hi] = root(root, hi);
            orientation[hi] = flip[hi];
        }
        struct BlockVotes { std::array<int, 2> all{}, strong{}; };
        std::map<std::pair<int, int>, BlockVotes> block_votes;
        for (size_t ri = 0; ri < chunk.reads.size(); ++ri) {
            const auto& read = chunk.reads[ri];
            if (read.is_skipped || read.mapq < std::max(opts.min_mapq, kGapBridgeMinMapq) ||
                (read.alignment && (read.alignment->core.flag & (BAM_FSECONDARY | BAM_FSUPPLEMENTARY)))) continue;
            const auto& profile = chunk.read_var_profile[ri];
            std::map<int, std::array<int, 2>> observations;
            std::map<int, bool> confident_base;
            for (int hi = 0; hi < n_het; ++hi) {
                const int vi = valid_var_idx[het_var_idx[hi]];
                const auto& var = chunk.candidates[vi];
                if (vi < profile.start_var_idx || vi > profile.end_var_idx ||
                    var.key.type != VariantType::Snp || var.lcd_var_i_to_cate != kCandCleanHetSnp ||
                    variant_allele_slots(var) != 2) continue;
                const int allele = profile.alleles[vi - profile.start_var_idx];
                if (allele < 0) continue;
                const int hap = allele == var.hap_to_cons_alle[1] ? 0 :
                                allele == var.hap_to_cons_alle[2] ? 1 : -1;
                if (hap >= 0) {
                    ++observations[component[hi]][hap ^ orientation[hi]];
                    confident_base[component[hi]] |= clean_snp_has_confident_bam_observation(read, var, allele);
                }
            }
            struct Anchor { int block, hap; bool strong; };
            std::vector<Anchor> anchors;
            for (const auto& [block, counts] : observations) {
                if (std::min(counts[0], counts[1]) != 0) continue;
                anchors.push_back({block, counts[1] > counts[0],
                    std::max(counts[0], counts[1]) >= kGapBridgeMinAnchorSnps || confident_base[block]});
            }
            for (size_t i = 0; i < anchors.size(); ++i) {
                for (size_t j = i + 1; j < anchors.size(); ++j) {
                    auto& votes = block_votes[{anchors[i].block, anchors[j].block}];
                    const int direction = anchors[i].hap ^ anchors[j].hap;
                    ++votes.all[direction];
                    if (anchors[i].strong && anchors[j].strong) ++votes.strong[direction];
                }
            }
        }
        std::vector<Edge> block_edges;
        for (const auto& [blocks, votes] : block_votes) {
            // Even a sparse opposing read vetoes a single-read bridge.
            if ((votes.all[0] && votes.all[1]) || !(votes.strong[0] || votes.strong[1])) continue;
            block_edges.push_back({blocks.first, blocks.second, votes.all[0], votes.all[1]});
        }
        std::stable_sort(block_edges.begin(), block_edges.end(), [](const Edge& a, const Edge& b) {
            return a.agree + a.conflict > b.agree + b.conflict;
        });
        for (const auto& edge : block_edges) {
            const int left = root(root, edge.left), right = root(root, edge.right);
            if (left == right) continue;
            const int child = std::max(left, right);
            parent[child] = std::min(left, right);
            flip[child] = flip[edge.left] ^ flip[edge.right] ^ (edge.conflict > edge.agree);
            if (opts.verbose >= 2)
                std::fprintf(stderr, "GapCleanBlockLink\t%lld\t%lld\t%d\t%d\n",
                    static_cast<long long>(chunk.candidates[valid_var_idx[het_var_idx[edge.left]]].key.sort_pos()),
                    static_cast<long long>(chunk.candidates[valid_var_idx[het_var_idx[edge.right]]].key.sort_pos()),
                    edge.agree, edge.conflict);
        }
        for (int hi = 0; hi < n_het; ++hi) {
            const int anchor = root(root, hi);
            parity[hi] = flip[hi];
            het_ps[hi] = chunk.candidates[valid_var_idx[het_var_idx[anchor]]].key.sort_pos();
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
static int iter_update_var_hap_to_cons_alle(PhasingChunk& chunk, bool is_ont,
                                             const std::vector<int>& valid_var_idx,
                                             uint32_t flags) {
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
        int hap = init_assign_read_hap(chunk, read_i, flags);
        if (hap == -1) hap = 0;
        chunk.haps[read_i] = hap;
        update_var_hap_profile(chunk, read_i, hap, flags);
    }

    for (int _vi = 0; _vi < n; ++_vi) {
        CandidateVariant& var = chunk.candidates[valid_var_idx[_vi]];
        if (var.gap_link_supported && var.msa_insertion_alts.size() == 2 && saved[_vi][1] > 0 &&
            saved[_vi][2] > 0 && saved[_vi][1] != saved[_vi][2]) {
            const int same = var.hap_to_alle_profile[1][1] + var.hap_to_alle_profile[2][2];
            const int flip = var.hap_to_alle_profile[1][2] + var.hap_to_alle_profile[2][1];
            // Optimize the orientation of the verified allele pair jointly.
            // Independent haplotype majorities can select the same allele twice.
            if (same != flip) {
                var.hap_to_cons_alle[1] = same > flip ? 1 : 2;
                var.hap_to_cons_alle[2] = same > flip ? 2 : 1;
            }
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

static void select_gap_link_sites(PhasingChunk& chunk, const Options& opts,
                                      const std::vector<int>& valid_var_idx) {
    int64_t* overlaps = nullptr;
    int64_t capacity = 0;
    for (const int vi : valid_var_idx) {
        auto& var = chunk.candidates[vi];
        const bool multi = !var.msa_insertion_alts.empty();
        var.gap_link_supported = multi;
        if (!multi && (!var.is_homopolymer_indel || var.lcd_var_i_to_cate != kCandNoisyCandHet ||
            opts.gap_hp_link_beg < 0 || var.key.sort_pos() < opts.gap_hp_link_beg ||
            var.key.sort_pos() > opts.gap_hp_link_end)) continue;
        std::map<hts_pos_t, std::array<int, 4>> votes;
        const int64_t n = cr_overlap(chunk.read_var_cr.get(), "cr", vi, vi + 1, &overlaps, &capacity);
        for (int64_t oi = 0; oi < n; ++oi) {
            const int ri = static_cast<int>(cr_label(chunk.read_var_cr.get(), overlaps[oi]));
            if (chunk.reads[ri].is_skipped || static_cast<size_t>(ri) >= chunk.haps.size() ||
                static_cast<size_t>(ri) >= chunk.phase_sets.size() ||
                chunk.haps[ri] < 1 || chunk.haps[ri] > 2 || chunk.phase_sets[ri] < 0) continue;
            const auto& profile = chunk.read_var_profile[ri];
            const int allele = profile.alleles[vi - profile.start_var_idx] - (multi ? 1 : 0);
            if (allele != 0 && allele != 1) continue;
            ++votes[chunk.phase_sets[ri]][2 * (chunk.haps[ri] - 1) + allele];
        }
        int best_anchor_depth = -1, best_anchor_total = -1;
        bool direct_supported = false;
        // Evaluate clean-site observations independently of provisional read
        // tags: unassigned spanning reads can still establish allele linkage.
        for (const int anchor_i : valid_var_idx) {
            const auto& anchor = chunk.candidates[anchor_i];
            if (anchor_i == vi || anchor.is_homopolymer_indel ||
                (anchor.lcd_var_i_to_cate & kCandGermlineClean) == 0 ||
                anchor.hap_to_cons_alle[1] < 0 || anchor.hap_to_cons_alle[2] < 0 ||
                anchor.hap_to_cons_alle[1] == anchor.hap_to_cons_alle[2]) continue;
            std::array<int, 4> v{};
            for (int64_t oi = 0; oi < n; ++oi) {
                const int ri = static_cast<int>(cr_label(chunk.read_var_cr.get(), overlaps[oi]));
                if (chunk.reads[ri].is_skipped) continue;
                const auto& profile = chunk.read_var_profile[ri];
                if (anchor_i < profile.start_var_idx || anchor_i > profile.end_var_idx) continue;
                const int allele = profile.alleles[vi - profile.start_var_idx] - (multi ? 1 : 0);
                const int anchor_allele = profile.alleles[anchor_i - profile.start_var_idx];
                if ((allele != 0 && allele != 1) || (anchor_allele != 0 && anchor_allele != 1)) continue;
                ++v[2 * anchor_allele + allele];
            }
            const int first = v[0] - v[1], second = v[3] - v[2];
            const int margin = opts.min_block_link_reads;
            const bool supported = multi ?
                ((first > 0 && second > 0) || (first < 0 && second < 0)) &&
                    std::abs(first + second) >= margin :
                (first >= margin && second >= margin) || (first <= -margin && second <= -margin);
            const int depth = std::min(v[0] + v[1], v[2] + v[3]);
            const int total = v[0] + v[1] + v[2] + v[3];
            if (depth > best_anchor_depth || (depth == best_anchor_depth && total > best_anchor_total)) {
                best_anchor_depth = depth;
                best_anchor_total = total;
                direct_supported = supported;
            } else if (depth == best_anchor_depth && total == best_anchor_total) {
                direct_supported &= supported;
            }
            if (opts.verbose >= 2 && v[0] + v[1] + v[2] + v[3] > 0)
                std::fprintf(stderr, "GapSiteAnchor\t%lld\t%lld\t%d\t%d\t%d\t%d\t%d\n",
                    static_cast<long long>(var.key.sort_pos()), static_cast<long long>(anchor.key.sort_pos()),
                    v[0], v[1], v[2], v[3], supported ? 1 : 0);
        }
        for (const auto& [ps, v] : votes) {
            const int first = v[0] - v[1], second = v[3] - v[2];
            const int margin = opts.min_block_link_reads;
            const bool supported = (first >= margin && second >= margin) ||
                                   (first <= -margin && second <= -margin);
            var.gap_link_supported |= supported;
            if (opts.verbose >= 2)
                std::fprintf(stderr, "GapHetAnchor\t%lld\t%lld\t%d\t%d\t%d\t%d\t%d\n",
                    static_cast<long long>(var.key.sort_pos()), static_cast<long long>(ps),
                    v[0], v[1], v[2], v[3], supported ? 1 : 0);
        }
        // A sparse distant subset must not establish heterozygosity when a
        // better-covered clean anchor shows the same allele on both haplotypes.
        if (best_anchor_depth >= opts.min_block_link_reads)
            var.gap_link_supported = direct_supported;
    }
    free(overlaps);
}

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

    if (opts.gap_hp_link_beg >= 0 ||
        (opts.recover_gaps && opts.private_msa_admit_all_in_region))
        select_gap_link_sites(chunk, opts, valid_var_idx);
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
        const int c2 = iter_update_var_hap_to_cons_alle(chunk, is_ont, valid_var_idx, flags);
        if (c1 == 0 && c2 == 0) break;
    }

    // Phase 3: finalize per-read phase_sets.
    update_read_phase_set(chunk, var_is_valid);

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
