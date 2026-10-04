/**
 * @file collect_phase_noisy.cpp
 * @brief Step 4: iterative noisy-region MSA variant calling.
 *
 * Step 4 implementation: noisy-region MSA variant recall.
 *   `sort_noisy_regs`, `collect_noisy_vars1`, and the outer while-loop in
 *   `collect_var_main` that gates k-means re-runs on `kCandGermlineVarCate`.
 *
 * Alignment-heavy functions (`collect_noisy_reg_aln_strs` and helpers) live in
 * `align.cpp`; this file owns the MSA-to-candidate conversion, read profile
 * update, old/new candidate merge, and iterative k-means trigger.
 */

#include "collect_phase_noisy.hpp"
#include <map>
#include "collect_phase.hpp" // assign_hap_based_on_germline_het_vars_kmeans, kCandGermlineVarCate
#include "collect_var.hpp"   // exact_comp_cand_var, exact_comp_var_site

#include <algorithm>
#include <cassert>
#include <cstdio>
#include <numeric>

namespace pgphase_collect {

namespace {


constexpr uint8_t kGapBase = 5;
constexpr const char* kNt4Bases = "ACGTN";

char var_type_code_local(VariantType type) {
    switch (type) {
        case VariantType::Snp:
            return 'X';
        case VariantType::Insertion:
            return 'I';
        case VariantType::Deletion:
            return 'D';
    }
    return '?';
}

const char* category_name_local(VariantCategory category) {
    switch (category) {
        case VariantCategory::LowCoverage:         return "LOW_COV";
        case VariantCategory::LowAlleleFraction:   return "LOW_AF";
        case VariantCategory::StrandBias:          return "STRAND_BIAS";
        case VariantCategory::CleanHetSnp:         return "CLEAN_HET_SNP";
        case VariantCategory::CleanHetIndel:       return "CLEAN_HET_INDEL";
        case VariantCategory::CleanHom:            return "CLEAN_HOM";
        case VariantCategory::NoisyCandHet:        return "NOISY_CAND_HET";
        case VariantCategory::NoisyCandHom:        return "NOISY_CAND_HOM";
        case VariantCategory::NoisyResolved:       return "NOISY_RESOLVED";
        case VariantCategory::RepeatHetIndel:      return "REP_HET_INDEL";
        case VariantCategory::NonVariant:          return "NON_VAR";
    }
    return "UNKNOWN";
}

char nt4_to_base(uint8_t b) {
    return b < 4 ? kNt4Bases[b] : 'N';
}


int alt_len_for_key(const VariantKey& key) {
    return key.type == VariantType::Deletion ? 0 : static_cast<int>(key.alt.size());
}

void update_variant_depth_fields(CandidateVariant& var) {
    if (var.counts.alle_covs.empty())
        var.counts.alle_covs.assign(static_cast<size_t>(std::max(2, var.counts.n_uniq_alles)), 0);
    var.counts.n_uniq_alles = static_cast<int>(var.counts.alle_covs.size());
    var.counts.ref_cov = var.counts.alle_covs.empty() ? 0 : var.counts.alle_covs[0];
    var.counts.alt_cov = std::accumulate(var.counts.alle_covs.begin() + 1, var.counts.alle_covs.end(), 0);
    int alle_cov_sum = 0;
    for (int cov : var.counts.alle_covs) alle_cov_sum += cov;
    if (var.counts.total_cov <= 0) var.counts.total_cov = alle_cov_sum;
    var.counts.allele_fraction =
        var.counts.total_cov > 0 ? static_cast<double>(var.counts.alt_cov) / var.counts.total_cov : 0.0;
}

void set_noisy_category(CandidateVariant& var, VariantCategory category) {
    var.counts.category = category;
    var.lcd_var_i_to_cate = category_to_flag(category);
}

static void restore_stashed_initial_if_any(PhasingChunk& chunk, CandidateVariant& var) {
    auto& stash = chunk.erased_clean_signal_initial;
    for (size_t i = 0; i < stash.size(); ++i) {
        if (exact_comp_var_site(&stash[i].first, &var.key) != 0) continue;
        var.counts.candvarcate_initial = stash[i].second;
        if (i + 1 < stash.size()) stash[i] = std::move(stash.back());
        stash.pop_back();
        return;
    }
}

} // namespace

void update_read_var_profile_with_allele(int var_idx, int allele, int alt_qi, ReadVariantProfile& profile) {
    if (var_idx < 0) return;
    if (profile.start_var_idx < 0) {
        profile.start_var_idx = var_idx;
        profile.end_var_idx = var_idx;
        profile.alleles.assign(1, allele);
        profile.alt_qi.assign(1, alt_qi);
        if (!profile.graph_alleles.empty()) profile.graph_alleles.assign(1, -1);
        return;
    }
    // Recovery revisits already-built sparse profiles and can discover an
    // observation on either side of their current range. Grow every populated
    // index-parallel channel before writing so the candidate offset remains
    // valid and the original BAM/graph evidence stays aligned.
    if (var_idx < profile.start_var_idx) {
        const size_t grow =
            static_cast<size_t>(profile.start_var_idx - var_idx);
        profile.alleles.insert(profile.alleles.begin(), grow, -1);
        profile.alt_qi.insert(profile.alt_qi.begin(), grow, -1);
        if (!profile.graph_alleles.empty())
            profile.graph_alleles.insert(profile.graph_alleles.begin(), grow, -1);
        if (!profile.bam_alleles.empty())
            profile.bam_alleles.insert(profile.bam_alleles.begin(), grow, -1);
        if (!profile.bam_qi.empty())
            profile.bam_qi.insert(profile.bam_qi.begin(), grow, -1);
        if (!profile.bam_base_qualities.empty())
            profile.bam_base_qualities.insert(profile.bam_base_qualities.begin(), grow, 0);
        profile.start_var_idx = var_idx;
        profile.alleles.front() = allele;
        profile.alt_qi.front() = alt_qi;
        return;
    }
    if (var_idx > profile.end_var_idx) {
        const int gap = var_idx - profile.end_var_idx - 1;
        profile.alleles.insert(profile.alleles.end(), static_cast<size_t>(gap), -1);
        profile.alt_qi.insert(profile.alt_qi.end(), static_cast<size_t>(gap), -1);
        if (!profile.graph_alleles.empty())
            profile.graph_alleles.insert(
                profile.graph_alleles.end(), static_cast<size_t>(gap), -1);
        if (!profile.bam_alleles.empty())
            profile.bam_alleles.insert(
                profile.bam_alleles.end(), static_cast<size_t>(gap), -1);
        if (!profile.bam_qi.empty())
            profile.bam_qi.insert(
                profile.bam_qi.end(), static_cast<size_t>(gap), -1);
        // No BAM base was measured for the newly added primary observation.
        if (!profile.bam_base_qualities.empty())
            profile.bam_base_qualities.insert(
                profile.bam_base_qualities.end(), static_cast<size_t>(gap), 0);
        profile.end_var_idx = var_idx;
        profile.alleles.push_back(allele);
        profile.alt_qi.push_back(alt_qi);
        if (!profile.graph_alleles.empty()) profile.graph_alleles.push_back(-1);
        if (!profile.bam_alleles.empty()) profile.bam_alleles.push_back(-1);
        if (!profile.bam_qi.empty()) profile.bam_qi.push_back(-1);
        if (!profile.bam_base_qualities.empty()) profile.bam_base_qualities.push_back(0);
        return;
    }
    const int offset = var_idx - profile.start_var_idx;
    profile.alleles[static_cast<size_t>(offset)] = allele;
    profile.alt_qi[static_cast<size_t>(offset)] = alt_qi;
}

namespace {

std::vector<ReadVariantProfile> init_read_profiles(size_t n_reads) {
    std::vector<ReadVariantProfile> profiles(n_reads);
    for (size_t i = 0; i < n_reads; ++i)
        profiles[i].read_id = static_cast<int>(i);
    return profiles;
}


CandidateVariant make_noisy_candidate(const PhasingChunk& chunk,
                                      hts_pos_t ref_pos,
                                      VariantType type,
                                      int ref_len,
                                      uint8_t ref_base,
                                      std::string alt,
                                      uint8_t alt_ref_base,
                                      bool is_homopolymer_indel) {
    CandidateVariant var;
    var.key.tid = chunk.region.tid;
    var.key.pos = ref_pos;
    var.key.type = type;
    var.key.ref_len = ref_len;
    var.key.alt = std::move(alt);
    var.ref_base = type == VariantType::Snp && ref_base < 4 ? ref_base : 4;
    // INS/DEL: keep consensus anchor encoding (e.g. abPOA gap == 5) for VCF output;
    // SNP path does not use alt_ref_base (filled separately).
    var.alt_ref_base = type == VariantType::Snp ? static_cast<uint8_t>(4) : alt_ref_base;
    var.is_homopolymer_indel = is_homopolymer_indel;
    var.msa_verified = true;
    var.counts.n_uniq_alles = 2;
    var.counts.alle_covs.assign(2, 0);
    var.counts.category = VariantCategory::NoisyCandHom;
    var.lcd_var_i_to_cate = category_to_flag(VariantCategory::NoisyCandHom);
    return var;
}

std::vector<CandidateVariant> make_cand_vars_from_baln0(const Options& opts,
                                                        const PhasingChunk& chunk,
                                                        hts_pos_t noisy_reg_beg,
                                                        const std::vector<uint8_t>& ref_msa_seq,
                                                        const std::vector<uint8_t>& cons_msa_seq,
                                                        bool no_end_var) {
    std::vector<CandidateVariant> vars;
    vars.reserve(ref_msa_seq.size());

    hts_pos_t ref_pos = noisy_reg_beg;
    int i = 0;
    const int msa_len = static_cast<int>(ref_msa_seq.size());
    while (i < msa_len) {
        const uint8_t ref_base = ref_msa_seq[static_cast<size_t>(i)];
        const uint8_t cons_base = cons_msa_seq[static_cast<size_t>(i)];
        if (ref_base == kGapBase && cons_base == kGapBase) {
            ++i;
            continue;
        }
        if (ref_base == cons_base) {
            ++i;
            ++ref_pos;
            continue;
        }
        if (ref_base != kGapBase && cons_base != kGapBase) {
            const bool next_ref_non_gap = (i + 1 == msa_len) || ref_msa_seq[static_cast<size_t>(i + 1)] != kGapBase;
            const bool next_cons_non_gap = (i + 1 == msa_len) || cons_msa_seq[static_cast<size_t>(i + 1)] != kGapBase;
            if (next_ref_non_gap && next_cons_non_gap) {
                vars.push_back(make_noisy_candidate(chunk, ref_pos, VariantType::Snp, 1,
                                                    ref_base, std::string(1, nt4_to_base(cons_base)),
                                                    4, false));
            }
            ++i;
            ++ref_pos;
            continue;
        }
        if (ref_base == kGapBase) {
            int gap_len = 1;
            while (i + gap_len < msa_len &&
                   ref_msa_seq[static_cast<size_t>(i + gap_len)] == kGapBase &&
                   cons_msa_seq[static_cast<size_t>(i + gap_len)] != kGapBase) {
                ++gap_len;
            }
            if (no_end_var && (i - 1 < 0 || i + gap_len >= msa_len ||
                               ref_msa_seq[static_cast<size_t>(i - 1)] == kGapBase ||
                               ref_msa_seq[static_cast<size_t>(i + gap_len)] == kGapBase ||
                               cons_msa_seq[static_cast<size_t>(i - 1)] == kGapBase ||
                               cons_msa_seq[static_cast<size_t>(i + gap_len)] == kGapBase)) {
                i += gap_len;
                continue;
            }
            std::string alt;
            alt.reserve(static_cast<size_t>(gap_len));
            for (int k = 0; k < gap_len; ++k)
                alt.push_back(nt4_to_base(cons_msa_seq[static_cast<size_t>(i + k)]));

            const bool hp = gap_len < opts.min_sv_len &&
                            var_is_homopolymer_indel(chunk, ref_pos, VariantType::Insertion, 0, alt,
                                                      opts.upstream_msa_insertion_hp);
            const uint8_t alt_ref_base = i - 1 >= 0 ? cons_msa_seq[static_cast<size_t>(i - 1)] : 4;
            vars.push_back(make_noisy_candidate(chunk, ref_pos, VariantType::Insertion, 0,
                                                4, alt, alt_ref_base, hp));
            i += gap_len;
            continue;
        }

        int gap_len = 1;
        while (i + gap_len < msa_len &&
               ref_msa_seq[static_cast<size_t>(i + gap_len)] != kGapBase &&
               cons_msa_seq[static_cast<size_t>(i + gap_len)] == kGapBase) {
            ++gap_len;
        }
        if (no_end_var && (i - 1 < 0 || i + gap_len >= msa_len ||
                           ref_msa_seq[static_cast<size_t>(i - 1)] == kGapBase ||
                           ref_msa_seq[static_cast<size_t>(i + gap_len)] == kGapBase ||
                           cons_msa_seq[static_cast<size_t>(i - 1)] == kGapBase ||
                           cons_msa_seq[static_cast<size_t>(i + gap_len)] == kGapBase)) {
            i += gap_len;
            ref_pos += gap_len;
            continue;
        }
        const bool hp = gap_len < opts.min_sv_len &&
                        var_is_homopolymer_indel(chunk, ref_pos, VariantType::Deletion, gap_len, {});
        const uint8_t alt_ref_base = i - 1 >= 0 ? cons_msa_seq[static_cast<size_t>(i - 1)] : 4;
        vars.push_back(make_noisy_candidate(chunk, ref_pos, VariantType::Deletion, gap_len,
                                            4, {}, alt_ref_base, hp));
        i += gap_len;
        ref_pos += gap_len;
    }
    return vars;
}

std::vector<CandidateVariant> make_cand_vars_from_msa(const Options& opts,
                                                      const PhasingChunk& chunk,
                                                      hts_pos_t noisy_reg_beg,
                                                      const std::vector<uint8_t>& ref_msa_seq,
                                                      const std::vector<uint8_t>& cons_msa_seq,
                                                      int msa_len,
                                                      bool no_end_var) {
    std::vector<uint8_t> ref_filtered;
    std::vector<uint8_t> cons_filtered;
    ref_filtered.reserve(static_cast<size_t>(msa_len));
    cons_filtered.reserve(static_cast<size_t>(msa_len));
    for (int i = 0; i < msa_len; ++i) {
        if (ref_msa_seq[static_cast<size_t>(i)] != kGapBase ||
            cons_msa_seq[static_cast<size_t>(i)] != kGapBase) {
            ref_filtered.push_back(ref_msa_seq[static_cast<size_t>(i)]);
            cons_filtered.push_back(cons_msa_seq[static_cast<size_t>(i)]);
        }
    }
    return make_cand_vars_from_baln0(opts, chunk, noisy_reg_beg, ref_filtered, cons_filtered, no_end_var);
}

int is_match_aln_str(const AlnStr& aln_str,
                     int target_pos,
                     int len,
                     float cons_sim_thres,
                     int* full_cover) {
    int cur_pos = -1;
    int n_eq = 0;
    int n_xid = 0;
    int cover_start = 0;
    int cover_end = 0;
    const int start_pos = target_pos < 0 ? 0 : target_pos;
    const int end_pos = target_pos < 0 ? len - 1 : target_pos + len - 1;

    for (int i = 0; i < aln_str.aln_len; ++i) {
        if (aln_str.target_aln[static_cast<size_t>(i)] != kGapBase) ++cur_pos;
        if (cur_pos == target_pos + len) break;
        if (i < aln_str.query_beg || i < aln_str.target_beg) continue;
        if (i > aln_str.query_end || i > aln_str.target_end) break;

        if (cur_pos == start_pos) cover_start = 1;
        if (cur_pos == end_pos) cover_end = 1;
        if (cur_pos >= target_pos) {
            if (aln_str.query_aln[static_cast<size_t>(i)] == aln_str.target_aln[static_cast<size_t>(i)]) ++n_eq;
            else ++n_xid;
        }
    }
    *full_cover = cover_start && cover_end;
    if (len >= 10) {
        if (n_eq >= len * cons_sim_thres) return 1;
        return *full_cover ? 0 : -1;
    }
    if (n_eq == len && n_xid == 0) return 1;
    return *full_cover ? 0 : -1;
}

int is_match_aln_str_del(const AlnStr& aln_str,
                         int target_del_left,
                         int target_del_right,
                         int* full_cover) {
    int cur_pos = -1;
    int started_check_del = 0;
    int n_non_del = 0;
    int cover_start = 0;
    int cover_end = 0;
    const int start_pos = target_del_left < 0 ? 0 : target_del_left;
    const int end_pos = target_del_left < 0 ? target_del_right : target_del_right;

    for (int i = 0; i < aln_str.aln_len; ++i) {
        if (aln_str.target_aln[static_cast<size_t>(i)] != kGapBase) ++cur_pos;
        if (cur_pos > target_del_right) break;
        if (i < aln_str.query_beg || i < aln_str.target_beg) continue;
        if (i > aln_str.query_end || i > aln_str.target_end) break;

        if (cur_pos == start_pos) cover_start = 1;
        if (cur_pos == end_pos) cover_end = 1;
        if (cur_pos >= target_del_left && cur_pos < target_del_right) {
            if (started_check_del == 0) {
                started_check_del = 1;
            } else if (aln_str.query_aln[static_cast<size_t>(i)] != kGapBase) {
                ++n_non_del;
            }
        }
    }
    if (cover_start && cover_end) {
        *full_cover = 1;
        return n_non_del == 0 ? 1 : 0;
    }
    *full_cover = 0;
    return -1;
}

int get_var_allele_i_from_cons_aln_str(const AlnStr& cons_aln_str,
                                       const CandidateVariant& var,
                                       int alt_pos,
                                       float cons_sim_thres,
                                       int* full_cover) {
    *full_cover = 0;
    if (var.key.type == VariantType::Snp) {
        assert(alt_len_for_key(var.key) == 1);
        return is_match_aln_str(cons_aln_str, alt_pos, 1, cons_sim_thres, full_cover);
    }
    if (var.key.type == VariantType::Insertion) {
        return is_match_aln_str(cons_aln_str, alt_pos, alt_len_for_key(var.key),
                                cons_sim_thres, full_cover);
    }
    return is_match_aln_str_del(cons_aln_str, alt_pos - 1, alt_pos, full_cover);
}

int is_cover_aln_str(const AlnStr& aln_str, int target_pos, int len) {
    int cur_pos = -1;
    int cover_start = 0;
    int cover_end = 0;
    const int start_pos = target_pos < 0 ? 0 : target_pos;
    const int end_pos = target_pos < 0 ? len - 1 : target_pos + len - 1;

    for (int i = 0; i < aln_str.aln_len; ++i) {
        if (aln_str.target_aln[static_cast<size_t>(i)] != kGapBase) ++cur_pos;
        if (i < aln_str.query_beg || i < aln_str.target_beg) continue;
        if (i > aln_str.query_end || i > aln_str.target_end) break;
        if (cur_pos == start_pos) cover_start = 1;
        if (cur_pos == end_pos) cover_end = 1;
        if (cover_start && cover_end) return 1;
    }
    return 0;
}

int get_full_cover_from_cons_aln_str(const AlnStr& cons_aln_str,
                                     const CandidateVariant& var,
                                     int alt_pos) {
    if (var.key.type == VariantType::Snp) {
        assert(var.key.ref_len == 1);
        return is_cover_aln_str(cons_aln_str, alt_pos, 1);
    }
    if (var.key.type == VariantType::Insertion) {
        return is_cover_aln_str(cons_aln_str, alt_pos, var.key.ref_len + 1);
    }
    return is_cover_aln_str(cons_aln_str, alt_pos - 1, var.key.ref_len + 1);
}

int get_full_cover_from_ref_cons_aln_str(const AlnStr& cons_aln_str,
                                         const AlnStr& ref_cons_aln_str,
                                         int beg_in_ref,
                                         int end_in_ref) {
    int cur_ref_pos = -1;
    int cur_cons_pos = -1;
    int beg_in_cons = -1;
    int end_in_cons = -1;
    int reach_end = 0;
    for (int i = 0; i < ref_cons_aln_str.aln_len; ++i) {
        if (ref_cons_aln_str.target_aln[static_cast<size_t>(i)] != kGapBase) ++cur_ref_pos;
        if (ref_cons_aln_str.query_aln[static_cast<size_t>(i)] != kGapBase) ++cur_cons_pos;
        if (i < ref_cons_aln_str.query_beg || i < ref_cons_aln_str.target_beg) continue;
        if (i > ref_cons_aln_str.query_end || i > ref_cons_aln_str.target_end) break;

        if (cur_ref_pos == beg_in_ref && beg_in_cons == -1) beg_in_cons = cur_cons_pos;
        if (cur_ref_pos == end_in_ref) reach_end = 1;
        if (reach_end && ref_cons_aln_str.query_aln[static_cast<size_t>(i)] != kGapBase) {
            end_in_cons = cur_cons_pos;
            break;
        }
    }
    return is_cover_aln_str(cons_aln_str, beg_in_cons, end_in_cons - beg_in_cons + 1);
}

void update_cand_var_profile_from_cons_aln_str(const AlnStr& cons_aln_str,
                                               hts_pos_t ref_pos_beg,
                                               std::vector<CandidateVariant>& vars,
                                               ReadVariantProfile& profile) {
    constexpr float cons_sim_thres = 0.9f;
    int delta_ref_alt = 0;
    for (int i = 0; i < static_cast<int>(vars.size()); ++i) {
        CandidateVariant& var = vars[static_cast<size_t>(i)];
        const int var_ref_pos = static_cast<int>(var.key.pos - ref_pos_beg);
        int full_cover = 0;
        const int allele_i = get_var_allele_i_from_cons_aln_str(
            cons_aln_str, var, var_ref_pos - delta_ref_alt, cons_sim_thres, &full_cover);
        if (full_cover) {
            ++var.counts.total_cov;
            if (var.counts.alle_covs.empty()) var.counts.alle_covs.assign(2, 0);
            if (allele_i >= 0 && allele_i < static_cast<int>(var.counts.alle_covs.size()))
                ++var.counts.alle_covs[static_cast<size_t>(allele_i)];
            update_read_var_profile_with_allele(i, allele_i, -1, profile);
        }
        if (var.key.type == VariantType::Insertion) delta_ref_alt -= alt_len_for_key(var.key);
        else if (var.key.type == VariantType::Deletion) delta_ref_alt += var.key.ref_len;
    }
}

void update_cand_var_profile_from_cons_aln_str1(int clu_n_seqs,
                                                const std::vector<int>& clu_read_ids,
                                                const std::vector<AlnStr>& clu_aln_strs,
                                                hts_pos_t ref_pos_beg,
                                                std::vector<CandidateVariant>& noisy_vars,
                                                std::vector<ReadVariantProfile>& profiles) {
    for (int i = 0; i < clu_n_seqs; ++i) {
        const size_t cons_read_slot = static_cast<size_t>((i + 1) * 2 - 1);
        if (cons_read_slot >= clu_aln_strs.size()) continue;
        const int read_id = clu_read_ids[static_cast<size_t>(i)];
        if (read_id < 0 || static_cast<size_t>(read_id) >= profiles.size()) continue;
        update_cand_var_profile_from_cons_aln_str(
            clu_aln_strs[cons_read_slot], ref_pos_beg, noisy_vars, profiles[static_cast<size_t>(read_id)]);
    }
}

void update_cand_var_profile_from_cons_aln_str21(int clu_idx,
                                                 const AlnStr& cons_aln_str,
                                                 const AlnStr& ref_cons_aln_str,
                                                 hts_pos_t ref_pos_beg,
                                                 std::vector<CandidateVariant>& vars,
                                                 const std::vector<int>& var_from_cons_idx,
                                                 ReadVariantProfile& profile) {
    constexpr float cons_sim_thres = 0.9f;
    int delta_ref_alt = 0;
    for (int i = 0; i < static_cast<int>(vars.size()); ++i) {
        CandidateVariant& var = vars[static_cast<size_t>(i)];
        const int var_beg_in_ref_str = static_cast<int>(var.key.pos - ref_pos_beg);
        const int var_end_in_ref_str = var.key.type == VariantType::Insertion
                                           ? var_beg_in_ref_str
                                           : var_beg_in_ref_str + var.key.ref_len - 1;
        int full_cover = 0;
        int allele_i = -1;
        if (var_from_cons_idx[static_cast<size_t>(i)] & clu_idx) {
            allele_i = get_var_allele_i_from_cons_aln_str(
                cons_aln_str, var, var_beg_in_ref_str - delta_ref_alt,
                cons_sim_thres, &full_cover);
        } else {
            if (var.key.type != VariantType::Deletion) {
                full_cover = get_full_cover_from_cons_aln_str(
                    cons_aln_str, var, var_beg_in_ref_str - delta_ref_alt);
            } else {
                full_cover = get_full_cover_from_ref_cons_aln_str(
                    cons_aln_str, ref_cons_aln_str,
                    var_beg_in_ref_str - 1, var_end_in_ref_str + 1);
            }
            allele_i = 0;
        }
        if (full_cover) {
            ++var.counts.total_cov;
            if (var.counts.alle_covs.empty()) var.counts.alle_covs.assign(2, 0);
            if (allele_i >= 0 && allele_i < static_cast<int>(var.counts.alle_covs.size()))
                ++var.counts.alle_covs[static_cast<size_t>(allele_i)];
            update_read_var_profile_with_allele(i, allele_i, -1, profile);
        }
        if (var_from_cons_idx[static_cast<size_t>(i)] & clu_idx) {
            if (var.key.type == VariantType::Insertion) delta_ref_alt -= alt_len_for_key(var.key);
            else if (var.key.type == VariantType::Deletion) delta_ref_alt += var.key.ref_len;
        }
    }
}

int update_cand_var_profile_from_cons_aln_str2(const Options& opts,
                                               const PhasingChunk& chunk,
                                               const std::array<int, 2>& clu_n_seqs,
                                               const std::array<std::vector<int>, 2>& clu_read_ids,
                                               const std::array<std::vector<AlnStr>, 2>& aln_strs,
                                               hts_pos_t noisy_reg_beg,
                                               std::vector<CandidateVariant>& hap1_vars,
                                               std::vector<CandidateVariant>& hap2_vars,
                                               std::vector<CandidateVariant>& noisy_vars,
                                               std::vector<VariantCategory>& noisy_var_cate,
                                               std::vector<ReadVariantProfile>& profiles) {
    if (hap1_vars.empty() && hap2_vars.empty()) return 0;

    std::vector<int> var_from_cons_idx;
    noisy_vars.clear();
    noisy_var_cate.clear();
    var_from_cons_idx.reserve(hap1_vars.size() + hap2_vars.size());

    size_t i1 = 0, i2 = 0;
    while (i1 < hap1_vars.size() && i2 < hap2_vars.size()) {
        const int ret = exact_comp_var_site(&hap1_vars[i1].key, &hap2_vars[i2].key);
        if (ret < 0) {
            set_noisy_category(hap1_vars[i1], VariantCategory::NoisyCandHet);
            noisy_var_cate.push_back(VariantCategory::NoisyCandHet);
            var_from_cons_idx.push_back(1);
            noisy_vars.push_back(hap1_vars[i1++]);
        } else if (ret > 0) {
            set_noisy_category(hap2_vars[i2], VariantCategory::NoisyCandHet);
            noisy_var_cate.push_back(VariantCategory::NoisyCandHet);
            var_from_cons_idx.push_back(2);
            noisy_vars.push_back(hap2_vars[i2++]);
        } else {
            set_noisy_category(hap1_vars[i1], VariantCategory::NoisyCandHom);
            noisy_var_cate.push_back(VariantCategory::NoisyCandHom);
            var_from_cons_idx.push_back(3);
            noisy_vars.push_back(hap1_vars[i1++]);
            ++i2;
        }
    }
    while (i1 < hap1_vars.size()) {
        set_noisy_category(hap1_vars[i1], VariantCategory::NoisyCandHet);
        noisy_var_cate.push_back(VariantCategory::NoisyCandHet);
        var_from_cons_idx.push_back(1);
        noisy_vars.push_back(hap1_vars[i1++]);
    }
    while (i2 < hap2_vars.size()) {
        set_noisy_category(hap2_vars[i2], VariantCategory::NoisyCandHet);
        noisy_var_cate.push_back(VariantCategory::NoisyCandHet);
        var_from_cons_idx.push_back(2);
        noisy_vars.push_back(hap2_vars[i2++]);
    }

    profiles = init_read_profiles(chunk.reads.size());
    for (int ci = 0; ci < 2; ++ci) {
        const std::vector<AlnStr>& clu_aln_strs = aln_strs[static_cast<size_t>(ci)];
        if (clu_aln_strs.empty()) continue;
        const AlnStr& ref_cons_aln_str = clu_aln_strs[0];
        for (int j = 0; j < clu_n_seqs[static_cast<size_t>(ci)]; ++j) {
            const size_t cons_read_slot = static_cast<size_t>((j + 1) * 2 - 1);
            if (cons_read_slot >= clu_aln_strs.size()) continue;
            const int read_id = clu_read_ids[static_cast<size_t>(ci)][static_cast<size_t>(j)];
            if (read_id < 0 || static_cast<size_t>(read_id) >= profiles.size()) continue;
            update_cand_var_profile_from_cons_aln_str21(
                ci + 1, clu_aln_strs[cons_read_slot], ref_cons_aln_str,
                noisy_reg_beg, noisy_vars, var_from_cons_idx,
                profiles[static_cast<size_t>(read_id)]);
        }
    }
    for (CandidateVariant& var : noisy_vars) update_variant_depth_fields(var);

    if (opts.verbose >= 2 && !noisy_vars.empty()) {
        for (size_t i = 0; i < noisy_vars.size(); ++i) {
            const CandidateVariant& v = noisy_vars[i];
            const VariantCounts& c = v.counts;
            std::fprintf(stderr, "Var: %" PRId64 ", %d-%c-%d %d(%d,%d) %s\t%s\n",
                         static_cast<int64_t>(v.key.pos),
                         v.key.ref_len,
                         var_type_code_local(v.key.type),
                         alt_len_for_key(v.key),
                         c.total_cov,
                         c.ref_cov,
                         c.alt_cov,
                         category_name_local(c.category),
                         v.key.alt.c_str());
        }
        for (int hap = 1; hap <= 2; ++hap) {
            const int ci = hap - 1;
            for (int j = 0; j < clu_n_seqs[static_cast<size_t>(ci)]; ++j) {
                const int read_id = clu_read_ids[static_cast<size_t>(ci)][static_cast<size_t>(j)];
                if (read_id < 0 || static_cast<size_t>(read_id) >= chunk.reads.size()) continue;
                const ReadVariantProfile& p = profiles[static_cast<size_t>(read_id)];
                std::fprintf(stderr, "Hap%d-Read: %s start_var_i: %d, end_var_i: %d\n",
                             hap,
                             chunk.reads[static_cast<size_t>(read_id)].qname.c_str(),
                             p.start_var_idx,
                             p.end_var_idx);
                if (p.start_var_idx < 0) continue;
                for (int k = 0; k <= (p.end_var_idx - p.start_var_idx); ++k) {
                    const int var_i = p.start_var_idx + k;
                    if (var_i < 0 || static_cast<size_t>(var_i) >= noisy_vars.size()) continue;
                    const CandidateVariant& v = noisy_vars[static_cast<size_t>(var_i)];
                    const int allele = p.alleles[static_cast<size_t>(k)];
                    std::fprintf(stderr, "P\tVar: (%d) %" PRId64 " %d-%c-%d, allele: %d\n",
                                 k,
                                 static_cast<int64_t>(v.key.pos),
                                 v.key.ref_len,
                                 var_type_code_local(v.key.type),
                                 alt_len_for_key(v.key),
                                 allele);
                }
            }
        }
    }
    return static_cast<int>(noisy_vars.size());
}

void merge_read_var_profile_entries(const ReadVariantProfile* old_profile,
                                    const std::vector<int>& old_to_merged,
                                    const ReadVariantProfile* new_profile,
                                    const std::vector<int>& new_to_merged,
                                    ReadVariantProfile& merged_profile,
                                    int merged_var_limit) {
    int old_var_i = old_profile != nullptr ? old_profile->start_var_idx : 1;
    int old_end_var_i = old_profile != nullptr ? old_profile->end_var_idx : 0;
    int new_var_i = new_profile != nullptr ? new_profile->start_var_idx : 1;
    int new_end_var_i = new_profile != nullptr ? new_profile->end_var_idx : 0;

    if (old_profile == nullptr || old_to_merged.empty() || old_var_i < 0 || old_end_var_i < old_var_i) {
        old_var_i = 1;
        old_end_var_i = 0;
    }
    if (new_profile == nullptr || new_to_merged.empty() || new_var_i < 0 || new_end_var_i < new_var_i) {
        new_var_i = 1;
        new_end_var_i = 0;
    }

    while (true) {
        while (old_var_i <= old_end_var_i &&
               (old_var_i >= static_cast<int>(old_to_merged.size()) || old_to_merged[static_cast<size_t>(old_var_i)] < 0)) {
            ++old_var_i;
        }
        while (new_var_i <= new_end_var_i &&
               (new_var_i >= static_cast<int>(new_to_merged.size()) || new_to_merged[static_cast<size_t>(new_var_i)] < 0)) {
            ++new_var_i;
        }

        const int old_merged_i = old_var_i <= old_end_var_i
                                     ? old_to_merged[static_cast<size_t>(old_var_i)]
                                     : merged_var_limit;
        const int new_merged_i = new_var_i <= new_end_var_i
                                     ? new_to_merged[static_cast<size_t>(new_var_i)]
                                     : merged_var_limit;
        if (old_merged_i == merged_var_limit && new_merged_i == merged_var_limit) break;

        if (old_merged_i <= new_merged_i) {
            const int old_profile_i = old_var_i - old_profile->start_var_idx;
            update_read_var_profile_with_allele(
                old_merged_i,
                old_profile->alleles[static_cast<size_t>(old_profile_i)],
                old_profile->alt_qi[static_cast<size_t>(old_profile_i)],
                merged_profile);
            if (static_cast<size_t>(old_profile_i) <
                old_profile->graph_alleles.size()) {
                if (merged_profile.graph_alleles.size() <
                    merged_profile.alleles.size())
                    merged_profile.graph_alleles.resize(
                        merged_profile.alleles.size(), -1);
                merged_profile.graph_alleles[static_cast<size_t>(
                    old_merged_i - merged_profile.start_var_idx)] =
                        old_profile->graph_alleles[
                            static_cast<size_t>(old_profile_i)];
            }
            ++old_var_i;
        } else {
            const int new_profile_i = new_var_i - new_profile->start_var_idx;
            update_read_var_profile_with_allele(
                new_merged_i,
                new_profile->alleles[static_cast<size_t>(new_profile_i)],
                new_profile->alt_qi[static_cast<size_t>(new_profile_i)],
                merged_profile);
            if (static_cast<size_t>(new_profile_i) < new_profile->graph_alleles.size()) {
                merged_profile.graph_alleles.resize(merged_profile.alleles.size(), -1);
                merged_profile.graph_alleles[static_cast<size_t>(
                    new_merged_i - merged_profile.start_var_idx)] =
                    new_profile->graph_alleles[static_cast<size_t>(new_profile_i)];
            }
            ++new_var_i;
        }
    }
}

} // namespace

// ---------------------------------------------------------------------------
// Predicates with external linkage: declared in collect_phase_noisy.hpp so the
// Catch2 predicate tests can exercise them directly. Kept out of the anonymous
// namespace for that reason alone; they have no other callers outside this file.
// ---------------------------------------------------------------------------

uint8_t base_to_nt4(char base) {
    switch (base) {
        case 'A':
        case 'a':
            return 0;
        case 'C':
        case 'c':
            return 1;
        case 'G':
        case 'g':
            return 2;
        case 'T':
        case 't':
        case 'U':
        case 'u':
            return 3;
        default:
            return 4;
    }
}

bool var_is_homopolymer_indel(const PhasingChunk& chunk,
                              hts_pos_t ref_pos,
                              VariantType type,
                              int ref_len,
                              const std::string& alt,
                              bool upstream_reference_bytes) {
    // longcallD compares raw FASTA bytes with nt4-coded MSA insertion bases.
    if (type == VariantType::Snp) return false;
    const hts_pos_t off = ref_pos - chunk.ref_beg;
    if (off < 0) return false;
    const size_t idx0 = static_cast<size_t>(off);
    if (type == VariantType::Insertion) {
        if (alt.empty()) return false;
        if (idx0 + 5 > chunk.ref_seq.size()) return false;
        const uint8_t ins_base0 = base_to_nt4(alt[0]);
        if (ins_base0 > 3) return false;
        for (size_t i = 1; i < alt.size(); ++i) {
            if (base_to_nt4(alt[i]) != ins_base0) return false;
        }
        for (int i = 0; i < 5; ++i) {
            const uint8_t ref_base = upstream_reference_bytes
                                         ? static_cast<uint8_t>(chunk.ref_seq[idx0 + static_cast<size_t>(i)])
                                         : base_to_nt4(chunk.ref_seq[idx0 + static_cast<size_t>(i)]);
            if (ref_base != ins_base0) return false;
        }
        return true;
    }
    const size_t span = static_cast<size_t>(ref_len > 5 ? ref_len : 5);
    if (idx0 + span > chunk.ref_seq.size()) return false;
    const uint8_t ref_base0 = base_to_nt4(chunk.ref_seq[idx0]);
    if (ref_base0 > 3) return false;
    for (int i = 1; i < ref_len; ++i) {
        if (base_to_nt4(chunk.ref_seq[idx0 + static_cast<size_t>(i)]) != ref_base0) return false;
    }
    for (int i = 0; i < 5; ++i) {
        if (base_to_nt4(chunk.ref_seq[idx0 + static_cast<size_t>(i)]) != ref_base0) return false;
    }
    return true;
}

/// Fill in the strand tallies for candidates whose counts the MSA path built.
///
/// update_variant_depth_fields derives ref_cov and alt_cov from alle_covs and
/// leaves forward_ref/reverse_ref/forward_alt/reverse_alt untouched, so a
/// candidate discovered or re-counted by the MSA carries coverage with no
/// strand at all -- measured at 29 to 45 candidates per panel window, for
/// example chr20:48,147,225 (SNP C>G) at 44 reference and 7 alternate with
/// 0+0 on both. That is not merely a reporting gap: the ONT strand-bias screen
/// computes its expected value from forward_alt + reverse_alt and declines to
/// test when the sum is not positive, so those candidates are silently exempt
/// from a screen every initially-discovered candidate faces.
///
/// Derived in one sweep from the read profiles rather than accumulated at each
/// site that touches a count, which is the same contract the graph-only
/// backfill follows: one writer, and the result a pure function of the
/// profiles.
///
/// A record is filled only when the derivation reproduces the record's OWN
/// ref_cov and alt_cov. The counts pass through several stages after the
/// profiles are written, so a disagreement means the two no longer describe
/// the same set of observations, and apportioning a strand split across counts
/// it does not match would be inventing the tally rather than measuring it.
/// Such a record keeps its zeros and stays visible as an unfilled one.

void derive_msa_candidate_strand_counts(PhasingChunk& chunk) {
    const size_t n_cands = chunk.candidates.size();
    if (n_cands == 0 || chunk.read_var_profile.empty()) return;
    std::vector<std::array<int, 4>> tally(n_cands, {0, 0, 0, 0});  // fwdRef revRef fwdAlt revAlt

    for (const ReadVariantProfile& profile : chunk.read_var_profile) {
        if (profile.read_id < 0 ||
            static_cast<size_t>(profile.read_id) >= chunk.reads.size()) continue;
        if (profile.start_var_idx < 0) continue;
        const bool reverse = chunk.reads[static_cast<size_t>(profile.read_id)].reverse;
        for (size_t i = 0; i < profile.alleles.size(); ++i) {
            const int allele = profile.alleles[i];
            if (allele < 0) continue;
            const size_t vi = static_cast<size_t>(profile.start_var_idx) + i;
            if (vi >= n_cands) break;
            if (allele == 0) tally[vi][reverse ? 1 : 0] += 1;
            else tally[vi][reverse ? 3 : 2] += 1;
        }
    }

    for (size_t vi = 0; vi < n_cands; ++vi) {
        VariantCounts& c = chunk.candidates[vi].counts;
        const bool ref_missing = c.ref_cov > 0 && c.forward_ref + c.reverse_ref == 0;
        const bool alt_missing = c.alt_cov > 0 && c.forward_alt + c.reverse_alt == 0;
        if (!ref_missing && !alt_missing) continue;
        const auto& t = tally[vi];
        if (ref_missing) {
            if (t[0] + t[1] != c.ref_cov) continue;
        }
        if (alt_missing) {
            if (t[2] + t[3] != c.alt_cov) continue;
        }
        if (ref_missing) { c.forward_ref = t[0]; c.reverse_ref = t[1]; }
        if (alt_missing) { c.forward_alt = t[2]; c.reverse_alt = t[3]; }
    }
}


static bool bam_aligned_base_quality(const bam1_t* bam, hts_pos_t target, int min_bq,
                                     int* query_index = nullptr) {
    const auto* cigar = bam_get_cigar(bam);
    hts_pos_t ref_pos = bam->core.pos + 1;
    int query_pos = 0;
    for (uint32_t ci = 0; ci < bam->core.n_cigar; ++ci) {
        const int len = bam_cigar_oplen(cigar[ci]);
        const int consumption = bam_cigar_type(bam_cigar_op(cigar[ci]));
        if ((consumption & 3) == 3 && target >= ref_pos && target < ref_pos + len) {
            const int qi = query_pos + static_cast<int>(target - ref_pos);
            const int quality = bam_get_qual(bam)[qi];
            if (quality == 255 || quality < min_bq) return false;
            if (query_index != nullptr) *query_index = qi;
            return true;
        }
        if (consumption & 1) query_pos += len;
        if (consumption & 2) ref_pos += len;
    }
    return false;
}

int bam_exact_indel_allele(const bam1_t* bam, const CandidateVariant& var,
                                  int min_bq, int* alt_qi) {
    const auto* cigar = bam_get_cigar(bam);
    hts_pos_t ref_pos = bam->core.pos + 1;
    int query_pos = 0;
    for (uint32_t ci = 0; ci < bam->core.n_cigar; ++ci) {
        const int op = bam_cigar_op(cigar[ci]);
        const int len = bam_cigar_oplen(cigar[ci]);
        if (var.key.type == VariantType::Insertion && op == BAM_CINS && ref_pos == var.key.pos) {
            if (len != static_cast<int>(var.key.alt.size())) return -1;
            for (int j = 0; j < len; ++j) {
                const int quality = bam_get_qual(bam)[query_pos + j];
                const int expected = seq_nt16_table[static_cast<unsigned char>(var.key.alt[j])];
                if (quality == 255 || quality < min_bq ||
                    bam_seqi(bam_get_seq(bam), query_pos + j) != expected) return -1;
            }
            *alt_qi = query_pos;
            return 1;
        }
        if (var.key.type == VariantType::Deletion && op == BAM_CDEL && ref_pos == var.key.pos) {
            if (len != var.key.ref_len ||
                !bam_aligned_base_quality(bam, var.key.pos - 1, min_bq) ||
                !bam_aligned_base_quality(bam, var.key.pos + var.key.ref_len, min_bq)) return -1;
            *alt_qi = query_pos;
            return 1;
        }
        if (op == BAM_CINS && ref_pos == var.key.pos) return -1;
        if (op == BAM_CDEL && ref_pos < var.key.pos + std::max(1, var.key.ref_len) &&
            ref_pos + len > var.key.pos) return -1;
        const int consumption = bam_cigar_type(op);
        if (consumption & 1) query_pos += len;
        if (consumption & 2) ref_pos += len;
    }
    if (var.key.type == VariantType::Insertion) {
        int qi = -1;
        if (!bam_aligned_base_quality(bam, var.key.pos - 1, min_bq) ||
            !bam_aligned_base_quality(bam, var.key.pos, min_bq, &qi)) return -1;
        *alt_qi = qi;
        return 0;
    }
    int qi = -1;
    for (hts_pos_t pos = var.key.pos; pos < var.key.pos + var.key.ref_len; ++pos)
        if (!bam_aligned_base_quality(bam, pos, min_bq, &qi)) return -1;
    if (!bam_aligned_base_quality(bam, var.key.pos - 1, min_bq) ||
        !bam_aligned_base_quality(bam, var.key.pos + var.key.ref_len, min_bq)) return -1;
    *alt_qi = qi;
    return 0;
}

static constexpr hts_pos_t kMaxEquivalentDeletionLength = 64;

// Source recall keeps its CIGAR edit certificate. The final physical stitch
// can additionally certify the complete query allele across a shifted edit;
// this must not change the source solve's consensuses, genotypes or weak cuts.
template<class ReferenceBase>
static int equivalent_deletion_allele(
        const bam1_t* read, const CandidateVariant& deletion,
        const ReferenceBase& reference_base, int min_baseq, int* alt_qi,
        bool verify_query_sequence = false) {
    constexpr hts_pos_t kMaxEquivalentShift = 32;
    constexpr int kUnknownQuality = 255;
    const hts_pos_t target_pos = deletion.key.pos;
    const hts_pos_t length = deletion.key.ref_len;
    if (target_pos <= kMaxEquivalentShift || length <= 0 ||
        length > kMaxEquivalentDeletionLength) return -1;
    const hts_pos_t window_beg = target_pos - kMaxEquivalentShift;
    const hts_pos_t window_end = target_pos + length + kMaxEquivalentShift;
    struct IndelEvent {
        int op;
        hts_pos_t pos;
        hts_pos_t length;
    };
    std::vector<IndelEvent> indels;
    hts_pos_t ref_pos = read->core.pos + 1;
    const uint32_t* cigar = bam_get_cigar(read);
    for (uint32_t ci = 0; ci < read->core.n_cigar; ++ci) {
        const int op = bam_cigar_op(cigar[ci]);
        const hts_pos_t op_length = bam_cigar_oplen(cigar[ci]);
        if ((op == BAM_CINS && ref_pos >= window_beg &&
             ref_pos <= window_end) ||
            (op == BAM_CDEL && ref_pos < window_end &&
             ref_pos + op_length > window_beg))
            indels.push_back(IndelEvent{op, ref_pos, op_length});
        if (bam_cigar_type(op) & 2) ref_pos += op_length;
    }
    const auto equivalent = [&](const IndelEvent& event) {
        if (event.op != BAM_CDEL || event.length != length ||
            event.pos < target_pos - kMaxEquivalentShift ||
            event.pos > target_pos + kMaxEquivalentShift)
            return false;
        const hts_pos_t beg = std::min(target_pos, event.pos);
        const hts_pos_t end = std::max(target_pos, event.pos) + length;
        if (!verify_query_sequence) {
            std::string ref;
            ref.reserve(static_cast<size_t>(end - beg));
            for (hts_pos_t pos = beg; pos < end; ++pos) {
                const char base = reference_base(pos);
                if (base == 'N') return false;
                ref.push_back(base);
            }
            std::string expected = ref;
            expected.erase(static_cast<size_t>(target_pos - beg),
                           static_cast<size_t>(length));
            ref.erase(static_cast<size_t>(event.pos - beg),
                      static_cast<size_t>(length));
            return ref == expected;
        }
        // A shifted deletion plus CIGAR mismatches can encode the same allele.
        // Compare the complete query between surviving reference anchors;
        // substitutions are admissible only when the final sequence is exact.
        int query_beg = -1, query_end = -1;
        if (!bam_aligned_base_quality(read, beg - 1, min_baseq, &query_beg) ||
            !bam_aligned_base_quality(read, end, min_baseq, &query_end) ||
            query_beg > query_end)
            return false;
        std::string allele;
        allele.reserve(static_cast<size_t>(end - beg + 2 - length));
        for (hts_pos_t pos = beg - 1; pos <= end; ++pos) {
            if (pos >= target_pos && pos < target_pos + length) continue;
            const char base = reference_base(pos);
            if (base == 'N') return false;
            allele.push_back(base);
        }
        if (query_end - query_beg + 1 != static_cast<int>(allele.size())) return false;
        for (size_t offset = 0; offset < allele.size(); ++offset) {
            const int qi = query_beg + static_cast<int>(offset);
            const int quality = bam_get_qual(read)[qi];
            if (quality == kUnknownQuality || quality < min_baseq ||
                bam_seqi(bam_get_seq(read), qi) !=
                    seq_nt16_table[static_cast<unsigned char>(allele[offset])])
                return false;
        }
        return true;
    };
    std::optional<size_t> selected;
    for (size_t i = 0; i < indels.size(); ++i) {
        if (!equivalent(indels[i])) continue;
        if (selected) return -1;
        selected = i;
    }
    const hts_pos_t observed_pos = selected ?
        indels[*selected].pos : target_pos;
    const hts_pos_t check_beg = std::min(target_pos, observed_pos);
    const hts_pos_t check_end = std::max(target_pos, observed_pos) + length;
    // A second indel outside the verified reference span does not change the
    // allele. An indel within it makes either the ALT or REF call ambiguous.
    for (size_t i = 0; i < indels.size(); ++i) {
        if (selected && i == *selected) continue;
        const IndelEvent& event = indels[i];
        const bool overlaps = event.op == BAM_CINS ?
            event.pos >= check_beg - 1 && event.pos <= check_end :
            event.pos <= check_end &&
                event.pos + event.length > check_beg - 1;
        if (overlaps) return -1;
    }
    if (selected && verify_query_sequence) {
        if (alt_qi != nullptr) {
            int qi = -1;
            bam_aligned_base_quality(read, check_end, min_baseq, &qi);
            *alt_qi = qi;
        }
        return 1;
    }
    int last_query_index = -1;
    for (hts_pos_t pos = check_beg - 1; pos <= check_end; ++pos) {
        if (selected && pos >= observed_pos &&
            pos < observed_pos + length)
            continue;
        const char ref = reference_base(pos);
        int query_index = -1;
        if (ref == 'N' ||
            !bam_aligned_base_quality(read, pos, min_baseq, &query_index) ||
            bam_seqi(bam_get_seq(read), query_index) !=
                seq_nt16_table[static_cast<unsigned char>(ref)])
            return -1;
        last_query_index = query_index;
    }
    if (alt_qi != nullptr) *alt_qi = last_query_index;
    return selected ? 1 : 0;
}


int bam_equivalent_deletion_allele(
        const bam1_t* read, const CandidateVariant& deletion,
        ReferenceCache& reference, int tid, const bam_hdr_t* header,
        int min_baseq) {
    return equivalent_deletion_allele(read, deletion,
        [&reference, tid, header](hts_pos_t pos) {
            return reference.base(tid, pos, header);
        }, min_baseq, nullptr);
}

bool bam_matches_deletion_sequence(
        const bam1_t* read, const CandidateVariant& deletion,
        ReferenceCache& reference, int tid, const bam_hdr_t* header,
        int min_baseq) {
    return equivalent_deletion_allele(read, deletion,
        [&reference, tid, header](hts_pos_t pos) {
            return reference.base(tid, pos, header);
        }, min_baseq, nullptr, true) == 1;
}

// Edit equivalence certifies the allele, not the source block's orientation.
// Reuse the deletion repair's independent physical SNP gauge for insertions.
static std::optional<bool> physical_snp_source_alt_gauge(
        const bam1_t* read, const CandidateVariant& variant,
        const PhasingChunk& chunk, int min_baseq) {
    constexpr int kMinRepairBaseq = 30;
    constexpr hts_pos_t kMinCorroboratingSnpSpacing = 100;
    if (!is_phase_set_anchor(variant) || variant.hap_to_cons_alle[1] > 1 ||
        variant.hap_to_cons_alle[2] > 1)
        return std::nullopt;
    const int variant_hap = variant.hap_to_cons_alle[1] == 1 ? 1 : 2;
    bool corroborated = false;
    for (const CandidateVariant& snp : chunk.candidates) {
        if (snp.phase_set != variant.phase_set ||
            snp.counts.category != VariantCategory::CleanHetSnp ||
            snp.key.type != VariantType::Snp || snp.key.ref_len != 1 ||
            snp.key.alt.size() != 1 || snp.ref_base > 3 ||
            snp.hap_to_cons_alle[1] < 0 || snp.hap_to_cons_alle[1] > 1 ||
            snp.hap_to_cons_alle[2] != 1 - snp.hap_to_cons_alle[1] ||
            std::llabs(snp.key.pos - variant.key.pos) < kMinCorroboratingSnpSpacing)
            continue;
        int qi = -1;
        if (!bam_aligned_base_quality(read, snp.key.pos,
                                      std::max(min_baseq, kMinRepairBaseq), &qi))
            continue;
        const int base = bam_seqi(bam_get_seq(read), qi);
        const int allele = base == (1 << snp.ref_base) ? 0 :
            base == seq_nt16_table[static_cast<unsigned char>(snp.key.alt[0])] ? 1 : -1;
        if (allele < 0) continue;
        const int snp_hap = allele == snp.hap_to_cons_alle[1] ? 1 : 2;
        if (snp_hap != variant_hap) return false;
        corroborated = true;
    }
    return corroborated ? std::optional<bool>{true} : std::nullopt;
}

static int bam_recovery_deletion_allele(
        const bam1_t* read, const CandidateVariant& deletion,
        const PhasingChunk& chunk, int min_baseq, int* alt_qi) {
    const int exact = bam_exact_indel_allele(read, deletion, min_baseq, alt_qi);
    // Exact ALT already identifies this edit. Exact-position REF can still
    // contain the same deletion shifted beyond this row's footprint.
    if (exact == 1) return exact;
    // A one-versus-other MSA locus needs its original multi-allele gauge.
    // Certifying one shifted deletion alone cannot orient a new read there.
    if (std::any_of(chunk.candidates.begin(), chunk.candidates.end(),
                    [&deletion](const CandidateVariant& other) {
                        return other.msa_verified &&
                            other.key.type != VariantType::Snp &&
                            other.key.pos == deletion.key.pos &&
                            exact_comp_var_site(&other.key, &deletion.key) != 0;
                    }))
        return exact;
    constexpr int kMinRepairMapq = 30;
    constexpr int kMinRepairBaseq = 30;
    constexpr int kUnknownMapq = 255;
    if (read->core.qual < kMinRepairMapq || read->core.qual == kUnknownMapq)
        return exact;
    const int shifted = equivalent_deletion_allele(read, deletion,
        [&chunk](hts_pos_t pos) {
            const hts_pos_t offset = pos - chunk.ref_beg;
            if (offset < 0 || static_cast<size_t>(offset) >= chunk.ref_seq.size())
                return 'N';
            return kNt4Bases[base_to_nt4(chunk.ref_seq[static_cast<size_t>(offset)])];
        }, std::max(min_baseq, kMinRepairBaseq), alt_qi);
    if (shifted != 1 || deletion.phase_set <= 0 ||
        deletion.hap_to_cons_alle[1] < 0 || deletion.hap_to_cons_alle[1] > 1 ||
        deletion.hap_to_cons_alle[2] != 1 - deletion.hap_to_cons_alle[1])
        return exact;
    // Allele certainty alone cannot orient a newly observed repeat allele.
    // Require a separate clean SNP on this molecule in the same source gauge.
    const std::optional<bool> gauge =
        physical_snp_source_alt_gauge(read, deletion, chunk, min_baseq);
    // Without a physical SNP, zero retains the MSA ALT-absence contrast.
    // A contradictory SNP makes this new call unknown rather than literal REF.
    if (!gauge) return exact;
    return *gauge ? 1 : -1;
}

// Keep the candidate coordinates: equivalence is an allele observation,
// not another variant row. Single-base insertions are equivalent precisely
// when the intervening reference consists of the same inserted base.
static int shifted_single_base_insertion_query_index(
        const bam1_t* bam, const CandidateVariant& var,
        const PhasingChunk& chunk, int min_bq) {
    constexpr hts_pos_t kMaxEquivalentShift = 32;
    constexpr int kUnknownQuality = 255;
    const hts_pos_t target = var.key.pos;
    const int alt = seq_nt16_table[static_cast<unsigned char>(var.key.alt[0])];
    const auto* cigar = bam_get_cigar(bam);
    hts_pos_t ref_pos = bam->core.pos + 1;
    int query_pos = 0;
    hts_pos_t observed_pos = -1;
    int observed_qi = -1;
    for (uint32_t ci = 0; ci < bam->core.n_cigar; ++ci) {
        const int op = bam_cigar_op(cigar[ci]);
        const int len = bam_cigar_oplen(cigar[ci]);
        if (op == BAM_CINS && ref_pos >= target - kMaxEquivalentShift &&
            ref_pos <= target + kMaxEquivalentShift) {
            const int quality = bam_get_qual(bam)[query_pos];
            if (observed_pos >= 0 || len != 1 || quality == kUnknownQuality ||
                quality < min_bq || bam_seqi(bam_get_seq(bam), query_pos) != alt)
                return -1;
            observed_pos = ref_pos;
            observed_qi = query_pos;
        }
        // A compound allele is not certified by a single inserted base.
        if (op == BAM_CDEL && ref_pos < target + kMaxEquivalentShift &&
            ref_pos + len > target - kMaxEquivalentShift)
            return -1;
        const int consumption = bam_cigar_type(op);
        if (consumption & 1) query_pos += len;
        if (consumption & 2) ref_pos += len;
    }
    if (observed_pos < 0 || observed_pos == target) return -1;
    const hts_pos_t beg = std::min(target, observed_pos);
    const hts_pos_t end = std::max(target, observed_pos);
    for (hts_pos_t pos = beg - 1; pos <= end; ++pos) {
        if (pos < chunk.ref_beg) return -1;
        const size_t offset = static_cast<size_t>(pos - chunk.ref_beg);
        if (offset >= chunk.ref_seq.size()) return -1;
        const uint8_t ref_base = base_to_nt4(chunk.ref_seq[offset]);
        if (ref_base > 3) return -1;
        const int ref = 1 << ref_base;
        if (pos >= beg && pos < end && ref != alt) return -1;
        int qi = -1;
        if (!bam_aligned_base_quality(bam, pos, min_bq, &qi) ||
            bam_seqi(bam_get_seq(bam), qi) != ref)
            return -1;
    }
    return observed_qi;
}

static void rebuild_read_variant_index(PhasingChunk& chunk) {
    cgranges_t* cr = cr_init();
    for (size_t read_i = 0; read_i < chunk.read_var_profile.size(); ++read_i) {
        const auto& profile = chunk.read_var_profile[read_i];
        if (profile.start_var_idx < 0 || profile.end_var_idx < profile.start_var_idx) continue;
        cr_add(cr, "cr", profile.start_var_idx, profile.end_var_idx + 1,
               static_cast<int32_t>(read_i));
    }
    cr_index(cr);
    chunk.read_var_cr.reset(cr);
}

template <class ReferenceBase>
static std::pair<hts_pos_t, hts_pos_t> insertion_equivalent_positions_impl(
        hts_pos_t pos, const std::string& alt, const ReferenceBase& reference) {
    hts_pos_t left = pos, right = pos;
    if (alt.empty()) return {left, right};
    const auto equal_base = [](char ref, char inserted) {
        return base_to_nt4(ref) <= 3 && base_to_nt4(ref) == base_to_nt4(inserted);
    };
    // Moving the edit one base rotates the insertion by one base. The
    // crossed reference base must equal the base rotated out of the allele.
    while (left > 1 && equal_base(reference(left - 1), alt[alt.size() - 1 -
               static_cast<size_t>(pos - left) % alt.size()])) --left;
    while (equal_base(reference(right), alt[static_cast<size_t>(right - pos) % alt.size()]))
        ++right;
    return {left, right};
}

std::pair<hts_pos_t, hts_pos_t> insertion_equivalent_positions(
        hts_pos_t pos, const std::string& alt, const PhasingChunk& chunk) {
    return insertion_equivalent_positions_impl(pos, alt, [&chunk](hts_pos_t base) {
        const hts_pos_t offset = base - chunk.ref_beg;
        return offset >= 0 && offset < static_cast<hts_pos_t>(chunk.ref_seq.size())
            ? chunk.ref_seq[static_cast<size_t>(offset)] : 'N';
    });
}

std::pair<hts_pos_t, hts_pos_t> insertion_equivalent_positions(
        hts_pos_t pos, const std::string& alt, ReferenceCache& reference,
        int tid, const bam_hdr_t* header) {
    return insertion_equivalent_positions_impl(pos, alt, [&](hts_pos_t base) {
        return reference.base(tid, base, header);
    });
}

bool insertion_edits_are_equivalent(
        hts_pos_t left_pos, std::string_view left_alt,
        hts_pos_t right_pos, std::string_view right_alt,
        std::string_view reference_between) {
    if (left_pos < 1 || right_pos < left_pos || left_alt.empty() ||
        left_alt.size() != right_alt.size() ||
        reference_between.size() != static_cast<size_t>(right_pos - left_pos))
        return false;
    const size_t shift = reference_between.size();
    const size_t length = left_alt.size();
    // Compare ALT + crossed background with crossed background + ALT.
    // No alignment, motif search limit or temporary edited strings are needed.
    for (size_t i = 0; i < shift + length; ++i) {
        const char left = i < length ? left_alt[i] : reference_between[i - length];
        const char right = i < shift ? reference_between[i] : right_alt[i - shift];
        if (base_to_nt4(left) > 3 || base_to_nt4(right) > 3 ||
            base_to_nt4(left) != base_to_nt4(right))
            return false;
    }
    return true;
}

// Certify a missing ALT by its reference edit, not its CIGAR coordinate.
// Keep existing MSA calls authoritative and reject compound repeat events.
template <class ReferenceBase>
static int shifted_repeat_insertion_query_index(
        const bam1_t* read, const CandidateVariant& var,
        const ReferenceBase& reference, int min_bq) {
    const hts_pos_t target = var.key.pos;
    const std::string& alt = var.key.alt;
    if (var.key.type != VariantType::Insertion || alt.size() < 2 ||
        !std::all_of(alt.begin(), alt.end(),
            [](char base) { return base_to_nt4(base) <= 3; })) return -1;
    const auto [left, right] = insertion_equivalent_positions_impl(target, alt, reference);
    hts_pos_t ref_pos = read->core.pos + 1;
    int query_pos = 0, observed_qi = -1;
    hts_pos_t observed_pos = 0;
    const uint32_t* cigar = bam_get_cigar(read);
    constexpr int kUnknownQuality = 255;
    for (uint32_t ci = 0; ci < read->core.n_cigar; ++ci) {
        const int op = bam_cigar_op(cigar[ci]);
        const int length = bam_cigar_oplen(cigar[ci]);
        if (op == BAM_CINS && left <= ref_pos && ref_pos <= right) {
            if (observed_qi >= 0 || static_cast<size_t>(length) != alt.size()) return -1;
            const hts_pos_t beg = std::min(ref_pos, target);
            const hts_pos_t end = std::max(ref_pos, target);
            std::string expected;
            expected.reserve(static_cast<size_t>(end - beg) + alt.size());
            for (hts_pos_t pos = beg; pos < end; ++pos)
                expected.push_back(reference(pos));
            std::string actual = expected;
            expected.insert(static_cast<size_t>(target - beg), alt);
            std::string inserted;
            for (int qi = query_pos; qi < query_pos + length; ++qi) {
                const int quality = bam_get_qual(read)[qi];
                if (quality < min_bq || quality == kUnknownQuality) return -1;
                inserted.push_back(seq_nt16_str[bam_seqi(bam_get_seq(read), qi)]);
            }
            actual.insert(static_cast<size_t>(ref_pos - beg), inserted);
            if (expected != actual) return -1;
            observed_qi = query_pos;
            observed_pos = ref_pos;
        }
        if (op == BAM_CDEL && ref_pos <= right && ref_pos + length > left) return -1;
        if (bam_cigar_type(op) & 1) query_pos += length;
        if (bam_cigar_type(op) & 2) ref_pos += length;
    }
    if (observed_qi < 0 || observed_pos == target) return -1;
    const hts_pos_t beg = std::min(target, observed_pos);
    const hts_pos_t end = std::max(target, observed_pos);
    for (hts_pos_t pos = beg - 1; pos <= end; ++pos) {
        int qi = -1;
        const char ref = reference(pos);
        if (base_to_nt4(ref) > 3 || !bam_aligned_base_quality(read, pos, min_bq, &qi) ||
            bam_seqi(bam_get_seq(read), qi) != seq_nt16_table[static_cast<unsigned char>(ref)])
            return -1;
    }
    return observed_qi;
}

int bam_shifted_repeat_insertion_query_index(
        const bam1_t* read, const CandidateVariant& insertion,
        const PhasingChunk& chunk, int min_baseq) {
    return shifted_repeat_insertion_query_index(read, insertion,
        [&chunk](hts_pos_t pos) {
            const hts_pos_t offset = pos - chunk.ref_beg;
            return offset >= 0 && offset < static_cast<hts_pos_t>(chunk.ref_seq.size())
                ? kNt4Bases[base_to_nt4(chunk.ref_seq[static_cast<size_t>(offset)])] : 'N';
        }, min_baseq);
}

int bam_shifted_repeat_insertion_query_index(
        const bam1_t* read, const CandidateVariant& insertion,
        ReferenceCache& reference, int tid, const bam_hdr_t* header, int min_baseq) {
    return shifted_repeat_insertion_query_index(read, insertion,
        [&reference, tid, header](hts_pos_t pos) {
            return reference.base(tid, pos, header);
        }, min_baseq);
}

int backfill_shifted_msa_insertions(PhasingChunk& chunk, const Options& opts) {
    constexpr int kMinRepairMapq = 30;
    constexpr int kMinRepairBaseq = 30;
    constexpr int kUnknownMapq = 255;
    if (opts.retry_windows.empty()) return 0;
    const int min_bq = std::max(opts.min_bq, kMinRepairBaseq);
    const int min_mapq = std::max(opts.min_mapq, kMinRepairMapq);
    int added = 0;
    for (size_t vi = 0; vi < chunk.candidates.size(); ++vi) {
        const CandidateVariant& var = chunk.candidates[vi];
        if (var.key.type != VariantType::Insertion || var.key.alt.size() != 1 ||
            base_to_nt4(var.key.alt[0]) > 3 || !var.msa_verified ||
            !var.msa_insertion_alts.empty() || var.is_homopolymer_indel ||
            var.lcd_var_i_to_cate != kCandNoisyCandHet ||
            !std::any_of(opts.retry_windows.begin(), opts.retry_windows.end(),
                         [&var](const auto& window) {
                             // Seams use VCF anchors; the inserted base's
                             // internal key lies one base to their right.
                             const hts_pos_t pos = var.key.sort_pos();
                             return window.first <= pos && pos <= window.second;
                         }))
            continue;
        for (size_t ri = 0; ri < chunk.reads.size(); ++ri) {
            const ReadRecord& read = chunk.reads[ri];
            if (read.is_skipped || !read.alignment ||
                read.mapq < min_mapq ||
                read.mapq == kUnknownMapq ||
                read.beg > var.key.pos || read.end < var.key.pos)
                continue;
            ReadVariantProfile& profile = chunk.read_var_profile[ri];
            if (static_cast<int>(vi) >= profile.start_var_idx &&
                static_cast<int>(vi) <= profile.end_var_idx &&
                profile.alleles[vi - profile.start_var_idx] != -1)
                continue;
            const int qi = shifted_single_base_insertion_query_index(
                read.alignment.get(), var, chunk, min_bq);
            if (qi < 0) continue;
            update_read_var_profile_with_allele(static_cast<int>(vi), 1, qi, profile);
            ++added;
        }
    }
    if (added > 0) rebuild_read_variant_index(chunk);
    return added;
}

int backfill_msa_retry_deletions(PhasingChunk& chunk, const Options& opts,
                                 hts_pos_t beg, hts_pos_t end, bool include_isolated) {
    // A homopolymer MSA row is an ALT-versus-other contrast, so a literal
    // CIGAR REF call alone cannot restore its allele. Two complementary rows
    // in one source block define the complete diploid contrast: a verified
    // ALT of exactly one row also proves ALT absence at the other row.
    std::map<hts_pos_t, std::vector<int>> deletion_loci;
    for (size_t vi = 0; vi < chunk.candidates.size(); ++vi) {
        const CandidateVariant& var = chunk.candidates[vi];
        const hts_pos_t pos = var.key.sort_pos();
        if (pos >= beg && pos <= end && var.msa_verified &&
            var.key.type == VariantType::Deletion)
            deletion_loci[var.key.pos].push_back(static_cast<int>(vi));
    }
    std::vector<std::pair<int, int>> deletion_pairs;
    std::vector<int> isolated_deletions;
    const auto eligible = [](const CandidateVariant& var) {
        return var.key.ref_len > 0 && var.key.ref_len <= kMaxEquivalentDeletionLength &&
            var.key.alt.empty() && var.msa_insertion_alts.empty() &&
            var.counts.n_uniq_alles == 2 &&
            var.counts.category == VariantCategory::NoisyCandHet &&
            var.lcd_var_i_to_cate == kCandNoisyCandHet &&
            var.phase_set > 0 && var.hap_to_cons_alle[1] >= 0 &&
            var.hap_to_cons_alle[1] <= 1 &&
            var.hap_to_cons_alle[2] == 1 - var.hap_to_cons_alle[1];
    };
    for (const auto& [pos, locus] : deletion_loci) {
        if (include_isolated && locus.size() == 1) {
            const CandidateVariant& var = chunk.candidates[locus[0]];
            // A physical observation diagnoses a missing or reversed MSA edge;
            // it must never directly inherit the frozen source haplotype.
            if (eligible(var) && var.is_homopolymer_indel &&
                std::none_of(chunk.candidates.begin(), chunk.candidates.end(),
                    [&var, pos](const CandidateVariant& other) {
                        return other.msa_verified && other.key.type != VariantType::Snp &&
                            other.key.pos == pos && &other != &var;
                    }))
                isolated_deletions.push_back(locus[0]);
            continue;
        }
        if (locus.size() != 2) continue;
        const CandidateVariant& first = chunk.candidates[locus[0]];
        const CandidateVariant& second = chunk.candidates[locus[1]];
        if (eligible(first) && eligible(second) &&
            first.key.ref_len != second.key.ref_len &&
            first.phase_set == second.phase_set &&
            first.hap_to_cons_alle[1] == second.hap_to_cons_alle[2])
            deletion_pairs.emplace_back(locus[0], locus[1]);
    }
    if (deletion_pairs.empty() && isolated_deletions.empty()) return 0;
    int added = 0;
    for (size_t read_i = 0; read_i < chunk.reads.size(); ++read_i) {
        const auto& read = chunk.reads[read_i];
        if (read.is_skipped || !read.alignment || read.mapq < opts.min_mapq) continue;
        const bam1_t* bam = read.alignment.get();
        constexpr int kMinRetryMapq = 30;
        constexpr int kMinRetryBaseq = 30;
        constexpr int kUnknownMapq = 255;
        if (read.mapq >= kMinRetryMapq && read.mapq != kUnknownMapq) {
            ReadVariantProfile& profile = chunk.read_var_profile[read_i];
            const auto unknown = [&](int vi) {
                return vi < profile.start_var_idx || vi > profile.end_var_idx ||
                    profile.alleles[static_cast<size_t>(vi - profile.start_var_idx)] < 0;
            };
            const auto reference_base = [&chunk](hts_pos_t pos) {
                const hts_pos_t offset = pos - chunk.ref_beg;
                if (offset < 0 || static_cast<size_t>(offset) >= chunk.ref_seq.size())
                    return 'N';
                return kNt4Bases[base_to_nt4(chunk.ref_seq[static_cast<size_t>(offset)])];
            };
            for (const int vi : isolated_deletions) {
                if (!unknown(vi)) continue;
                int qi = -1;
                const int allele = equivalent_deletion_allele(
                    bam, chunk.candidates[vi], reference_base,
                    std::max(opts.min_bq, kMinRetryBaseq), &qi);
                if (allele < 0) continue;
                update_read_var_profile_with_allele(vi, allele, qi, profile);
                ++added;
            }
            for (const auto& [first, second] : deletion_pairs) {
                // Existing MSA projections own their gauge. Only restore a
                // jointly absent pair; never overwrite or complete one call
                // using a different physical observation path.
                if (!unknown(first) || !unknown(second)) continue;
                int first_qi = -1;
                int second_qi = -1;
                const int baseq = std::max(opts.min_bq, kMinRetryBaseq);
                const int first_call = equivalent_deletion_allele(
                    bam, chunk.candidates[first], reference_base, baseq, &first_qi);
                const int second_call = equivalent_deletion_allele(
                    bam, chunk.candidates[second], reference_base, baseq, &second_qi);
                if ((first_call == 1) == (second_call == 1)) continue;
                const bool first_alt = first_call == 1;
                const int qi = first_alt ? first_qi : second_qi;
                update_read_var_profile_with_allele(first, first_alt ? 1 : 0, qi, profile);
                update_read_var_profile_with_allele(second, first_alt ? 0 : 1, qi, profile);
                added += 2;
            }
        }
    }
    if (added > 0) rebuild_read_variant_index(chunk);
    return added;
}

int backfill_msa_observations(PhasingChunk& chunk, const Options& opts,
                              hts_pos_t beg, hts_pos_t end) {
    std::vector<int> sites;
    for (size_t vi = 0; vi < chunk.candidates.size(); ++vi) {
        const auto& var = chunk.candidates[vi];
        const bool supported_snp = var.key.type == VariantType::Snp && var.key.ref_len == 1 &&
                                   var.key.alt.size() == 1 && var.ref_base <= 3;
        const bool supported_indel = var.key.type != VariantType::Snp &&
                                     var.msa_insertion_alts.empty() &&
                                     !var.is_homopolymer_indel;
        // Recovery seams use inclusive VCF anchors. Use that coordinate for
        // membership, but keep the internal position for CIGAR allele calls.
        const hts_pos_t pos = var.key.sort_pos();
        if (pos >= beg && pos <= end && var.msa_verified &&
            (supported_snp || supported_indel) &&
            var.lcd_var_i_to_cate == kCandNoisyCandHet)
            sites.push_back(static_cast<int>(vi));
    }
    if (sites.empty()) return 0;
    int added = 0;
    for (size_t read_i = 0; read_i < chunk.reads.size(); ++read_i) {
        const auto& read = chunk.reads[read_i];
        if (read.is_skipped || !read.alignment || read.mapq < opts.min_mapq) continue;
        const bam1_t* bam = read.alignment.get();
        const auto* cigar = bam_get_cigar(bam);
        hts_pos_t ref_pos = bam->core.pos + 1;
        int query_pos = 0;
        for (uint32_t ci = 0; ci < bam->core.n_cigar; ++ci) {
            const int len = bam_cigar_oplen(cigar[ci]);
            const int consumption = bam_cigar_type(bam_cigar_op(cigar[ci]));
            if ((consumption & 3) == 3) {
                const hts_pos_t match_end = ref_pos + len;
                for (const int vi : sites) {
                    const auto& var = chunk.candidates[static_cast<size_t>(vi)];
                    if (var.key.pos < ref_pos || var.key.pos >= match_end) continue;
                    auto& profile = chunk.read_var_profile[read_i];
                    if (vi >= profile.start_var_idx && vi <= profile.end_var_idx &&
                        profile.alleles[static_cast<size_t>(vi - profile.start_var_idx)] != -1) continue;
                    if (var.key.type != VariantType::Snp) continue;
                    const int qi = query_pos + static_cast<int>(var.key.pos - ref_pos);
                    const int quality = bam_get_qual(bam)[qi];
                    if (quality == 255 || quality < opts.min_bq) continue;
                    const int base = bam_seqi(bam_get_seq(bam), qi);
                    const int ref_base = 1 << var.ref_base;
                    const int alt_base = seq_nt16_table[static_cast<unsigned char>(var.key.alt[0])];
                    const int allele = base == ref_base ? 0 : base == alt_base ? 1 : -1;
                    if (allele < 0) continue;
                    update_read_var_profile_with_allele(vi, allele, qi, profile);
                    ++added;
                }
            }
            if (consumption & 1) query_pos += len;
            if (consumption & 2) ref_pos += len;
        }
        for (const int vi : sites) {
            const auto& var = chunk.candidates[static_cast<size_t>(vi)];
            if (var.key.type == VariantType::Snp) continue;
            auto& profile = chunk.read_var_profile[read_i];
            if (vi >= profile.start_var_idx && vi <= profile.end_var_idx &&
                profile.alleles[static_cast<size_t>(vi - profile.start_var_idx)] != -1) continue;
            int qi = -1;
            // Exact CIGARs can miss a deletion at another repeat placement.
            // Repair only an independently certified missing ALT; callable
            // binary MSA contrasts retain their original source gauge.
            int allele = var.key.type == VariantType::Deletion &&
                               var.key.ref_len <= kMaxEquivalentDeletionLength
                ? bam_recovery_deletion_allele(bam, var, chunk, opts.min_bq, &qi)
                : bam_exact_indel_allele(bam, var, opts.min_bq, &qi);
            // An exact-coordinate REF call can be the same insertion shifted
            // elsewhere in its repeat. Certify the unchanged edit before
            // adding that missing observation; existing MSA contrasts above
            // remain authoritative, including co-located complementary rows.
            constexpr int kMinShiftedInsertionMapq = 30;
            constexpr int kMinShiftedInsertionBaseq = 30;
            constexpr int kUnknownMapq = 255;
            if (var.key.type == VariantType::Insertion && allele != 1 &&
                read.mapq >= kMinShiftedInsertionMapq && read.mapq != kUnknownMapq) {
                const int shifted_qi = bam_shifted_repeat_insertion_query_index(
                    bam, var, chunk, std::max(opts.min_bq, kMinShiftedInsertionBaseq));
                if (shifted_qi >= 0) {
                    const std::optional<bool> gauge =
                        physical_snp_source_alt_gauge(bam, var, chunk, opts.min_bq);
                    if (gauge.value_or(false)) {
                        allele = 1;
                    } else {
                        // Separate complementary rows use ALT absence rather
                        // than literal REF. Preserve that source contrast until
                        // a joint call certifies both rows' allele gauge.
                        const bool complementary = std::any_of(
                            sites.begin(), sites.end(), [&](int other_i) {
                                const CandidateVariant& other = chunk.candidates[other_i];
                                return other.key.type == VariantType::Insertion &&
                                    other.key.pos == var.key.pos && other.key.alt != var.key.alt &&
                                    other.phase_set == var.phase_set &&
                                    is_phase_set_anchor(var) && is_phase_set_anchor(other) &&
                                    var.hap_to_cons_alle[1] >= 0 && var.hap_to_cons_alle[1] <= 1 &&
                                    var.hap_to_cons_alle[2] == 1 - var.hap_to_cons_alle[1] &&
                                    other.hap_to_cons_alle[1] == var.hap_to_cons_alle[2] &&
                                    other.hap_to_cons_alle[2] == var.hap_to_cons_alle[1];
                            });
                        // For an unpaired binary row the exact-anchor REF
                        // lookup missed its ALT. Contradictory clean SNPs mean
                        // unknown, not evidence for the opposite allele.
                        if (!complementary) {
                            if (gauge) allele = -1;
                            else {
                                // Preserve the discovery/retry matrix here. An
                                // exact shifted edit without a per-read SNP
                                // gauge can be admitted after source selection
                                // only if its block has one coherent graph gauge.
                                chunk.pending_msa_observations.push_back(
                                    {var.key, static_cast<int>(read_i), 1, false});
                            }
                        }
                    }
                    qi = shifted_qi;
                }
            }
            if (allele < 0) continue;
            if (allele == 0 && var.key.type == VariantType::Insertion) {
                for (const int other_i : sites) {
                    const auto& other = chunk.candidates[other_i];
                    if (other.key.type != VariantType::Insertion ||
                        other.key.pos != var.key.pos || other.key.alt == var.key.alt ||
                        other.phase_set != var.phase_set || !is_phase_set_anchor(var) ||
                        !is_phase_set_anchor(other) || var.hap_to_cons_alle[1] < 0 ||
                        var.hap_to_cons_alle[1] > 1 ||
                        var.hap_to_cons_alle[2] != 1 - var.hap_to_cons_alle[1] ||
                        other.hap_to_cons_alle[1] != var.hap_to_cons_alle[2] ||
                        other.hap_to_cons_alle[2] != var.hap_to_cons_alle[1]) continue;
                    const auto& shorter = var.key.alt.size() < other.key.alt.size()
                        ? var.key.alt : other.key.alt;
                    const auto& longer = var.key.alt.size() < other.key.alt.size()
                        ? other.key.alt : var.key.alt;
                    if (longer.compare(0, shorter.size(), shorter) != 0) continue;
                    for (uint32_t ci = 0; ci < bam->core.n_cigar; ++ci) {
                        const size_t length = bam_cigar_oplen(cigar[ci]);
                        if (bam_cigar_op(cigar[ci]) != BAM_CINS ||
                            length <= shorter.size() || length >= longer.size()) continue;
                        CandidateVariant third = var;
                        third.key.alt = longer.substr(0, length);
                        if (bam_shifted_repeat_insertion_query_index(
                                bam, third, chunk, opts.min_bq) < 0) continue;
                        // A certified third repeat length matches neither
                        // source allele. Literal-anchor absence cannot vote
                        // for both sides of their complementary contrast.
                        chunk.pending_msa_observations.push_back(
                            {var.key, static_cast<int>(read_i), -1, false});
                        break;
                    }
                }
            }
            update_read_var_profile_with_allele(vi, allele, qi, profile);
            ++added;
        }
    }
    if (added == 0) return 0;
    // Post-solve CIGAR calls extend linkage evidence, not the discovery
    // genotype. Preserve that census: re-genotyping from this supplementary
    // projection can change retry admission and the existing block gauge.
    rebuild_read_variant_index(chunk);
    return added;
}

// ════════════════════════════════════════════════════════════════════════════
// sort_noisy_regs
// ════════════════════════════════════════════════════════════════════════════

std::vector<int> sort_noisy_regs(const PhasingChunk& chunk) {
    const int n = static_cast<int>(chunk.noisy_regions.size());
    std::vector<int> idx(static_cast<size_t>(n));
    std::iota(idx.begin(), idx.end(), 0);

    // Bubble sort: primary key = Interval::label (variant count)
    // ascending; secondary key = noisy_reg_lens[i] (= cr_end - cr_start = end - beg) ascending.
    for (int i = 0; i < n; ++i) {
        for (int j = i + 1; j < n; ++j) {
            const Interval& a = chunk.noisy_regions[static_cast<size_t>(idx[i])];
            const Interval& b = chunk.noisy_regions[static_cast<size_t>(idx[j])];
            bool do_swap = false;
            if (a.label > b.label) {
                do_swap = true;
            } else if (a.label == b.label) {
                if ((a.end - a.beg) > (b.end - b.beg)) do_swap = true;
            }
            if (do_swap) std::swap(idx[i], idx[j]);
        }
    }
    return idx;
}

// ════════════════════════════════════════════════════════════════════════════
// collect_noisy_reg_reads
// ════════════════════════════════════════════════════════════════════════════

std::vector<int> collect_noisy_reg_reads(const PhasingChunk& chunk,
                                         hts_pos_t beg, hts_pos_t end) {
    std::vector<int> result;
    result.reserve(chunk.reads.size());

    // Iterate ordered_read_ids, skip
    // skipped reads and reads with no overlap (beg > end || end <= beg).
    // ReadRecord::beg/end = aligned region boundaries.
    const auto visit = [&](int ri) {
        if (ri < 0 || static_cast<size_t>(ri) >= chunk.reads.size()) return;
        const ReadRecord& r = chunk.reads[static_cast<size_t>(ri)];
        if (r.is_skipped) return;
        if (r.digars.empty()) return;  // graph-only reads have no sequences
        if (r.beg > end || r.end <= beg) return;
        result.push_back(ri);
    };
    if (chunk.ordered_read_ids.empty()) {
        for (size_t read_i = 0; read_i < chunk.reads.size(); ++read_i)
            visit(static_cast<int>(read_i));
    } else {
        for (int ri : chunk.ordered_read_ids) visit(ri);
    }
    return result;
}

// ════════════════════════════════════════════════════════════════════════════
// collect_reg_ref_bseq
// ════════════════════════════════════════════════════════════════════════════

// nt4 encoding table (A=0, C=1, G=2, T/U=3, other=4).
static const uint8_t kNstNt4Table[256] = {
    4,4,4,4, 4,4,4,4, 4,4,4,4, 4,4,4,4, // 0x00-0x0F
    4,4,4,4, 4,4,4,4, 4,4,4,4, 4,4,4,4, // 0x10-0x1F
    4,4,4,4, 4,4,4,4, 4,4,4,4, 4,4,4,4, // 0x20-0x2F (space, !"#...)
    4,4,4,4, 4,4,4,4, 4,4,4,4, 4,4,4,4, // 0x30-0x3F (0-9, :;<=>?)
    4,0,4,1, 4,4,4,2, 4,4,4,4, 4,4,4,4, // 0x40-0x4F (@A-O): A=0,C=1,G=2
    4,4,4,4, 3,3,4,4, 4,4,4,4, 4,4,4,4, // 0x50-0x5F (P-_): T=3,U=3
    4,0,4,1, 4,4,4,2, 4,4,4,4, 4,4,4,4, // 0x60-0x6F (`a-o): a=0,c=1,g=2
    4,4,4,4, 3,3,4,4, 4,4,4,4, 4,4,4,4, // 0x70-0x7F (p-DEL): t=3,u=3
    4,4,4,4, 4,4,4,4, 4,4,4,4, 4,4,4,4, // 0x80-0x8F
    4,4,4,4, 4,4,4,4, 4,4,4,4, 4,4,4,4, // 0x90-0x9F
    4,4,4,4, 4,4,4,4, 4,4,4,4, 4,4,4,4, // 0xA0-0xAF
    4,4,4,4, 4,4,4,4, 4,4,4,4, 4,4,4,4, // 0xB0-0xBF
    4,4,4,4, 4,4,4,4, 4,4,4,4, 4,4,4,4, // 0xC0-0xCF
    4,4,4,4, 4,4,4,4, 4,4,4,4, 4,4,4,4, // 0xD0-0xDF
    4,4,4,4, 4,4,4,4, 4,4,4,4, 4,4,4,4, // 0xE0-0xEF
    4,4,4,4, 4,4,4,4, 4,4,4,4, 4,4,4,4, // 0xF0-0xFF
};

std::vector<uint8_t> collect_reg_ref_bseq(const PhasingChunk& chunk,
                                           hts_pos_t& beg, hts_pos_t& end) {
    // Clip beg/end to ref slice boundaries.
    if (beg < chunk.ref_beg) beg = chunk.ref_beg;
    if (end > chunk.ref_end) end = chunk.ref_end;

    const hts_pos_t len = end - beg + 1;
    if (len <= 0 || chunk.ref_seq.empty()) return {};

    std::vector<uint8_t> out(static_cast<size_t>(len));
    for (hts_pos_t i = beg; i <= end; ++i) {
        out[static_cast<size_t>(i - beg)] =
            kNstNt4Table[static_cast<unsigned char>(
                chunk.ref_seq[static_cast<size_t>(i - chunk.ref_beg)])];
    }
    return out;
}

// ════════════════════════════════════════════════════════════════════════════
// make_vars_from_msa_cons_aln
// ════════════════════════════════════════════════════════════════════════════

static void split_nested_msa_deletions(const Options& opts, const PhasingChunk& chunk,
                                       std::vector<CandidateVariant>& first,
                                       std::vector<CandidateVariant>& second) {
    const size_t first_size = first.size(), second_size = second.size();
    for (size_t i = 0; i < first_size; ++i) {
        for (size_t j = 0; j < second_size; ++j) {
            const auto a = first[i], b = second[j];
            if (a.key.type != VariantType::Deletion || b.key.type != VariantType::Deletion ||
                a.key.pos != b.key.pos || a.key.ref_len == b.key.ref_len) continue;
            // Both haplotypes delete the common prefix. Treating each full
            // deletion as an independent het lets a longer-deletion read
            // support both ALT rows, manufacturing a false phase bridge.
            const bool first_longer = a.key.ref_len > b.key.ref_len;
            const auto& common = first_longer ? b : a;
            auto residual = first_longer ? a : b;
            residual.key.pos += common.key.ref_len;
            residual.key.ref_len -= common.key.ref_len;
            residual.alt_ref_base = base_to_nt4(chunk.ref_seq[
                static_cast<size_t>(residual.key.pos - 1 - chunk.ref_beg)]);
            residual.is_homopolymer_indel = residual.key.ref_len < opts.min_sv_len &&
                var_is_homopolymer_indel(chunk, residual.key.pos, VariantType::Deletion,
                                          residual.key.ref_len, {});
            auto& longer = first_longer ? first : second;
            longer[first_longer ? i : j] = common;
            longer.push_back(std::move(residual));
        }
    }
    const auto less = [](const CandidateVariant& a, const CandidateVariant& b) {
        return exact_comp_cand_var(&a, &b) < 0;
    };
    std::stable_sort(first.begin(), first.end(), less);
    std::stable_sort(second.begin(), second.end(), less);
}

static void merge_msa_insertion_alleles(std::vector<CandidateVariant>& vars,
    std::vector<VariantCategory>& categories, std::vector<ReadVariantProfile>& profiles) {
    for (size_t i = 0; i + 1 < vars.size(); ++i) {
        auto& first = vars[i];
        const auto& second = vars[i + 1];
        if (first.key.type != VariantType::Insertion || second.key.type != VariantType::Insertion ||
            first.key.pos != second.key.pos || first.key.alt == second.key.alt ||
            first.counts.category != VariantCategory::NoisyCandHet ||
            second.counts.category != VariantCategory::NoisyCandHet) continue;
        first.msa_insertion_alts = {first.key.alt, second.key.alt};
        first.counts.alle_covs.assign(3, 0);
        first.counts.total_cov = 0;
        for (auto& profile : profiles) {
            if (profile.start_var_idx < 0) continue;
            const auto allele_at = [&](int index) {
                return index < profile.start_var_idx || index > profile.end_var_idx ? -1 :
                       profile.alleles[index - profile.start_var_idx];
            };
            const int a = allele_at(static_cast<int>(i)), b = allele_at(static_cast<int>(i + 1));
            const int allele = a == 1 && b != 1 ? 1 : b == 1 && a != 1 ? 2 :
                               a == 0 && b == 0 ? 0 : -1;
            ReadVariantProfile merged;
            merged.read_id = profile.read_id;
            for (int vi = profile.start_var_idx; vi <= profile.end_var_idx; ++vi) {
                if (vi == static_cast<int>(i + 1)) {
                    if (profile.start_var_idx == vi)
                        update_read_var_profile_with_allele(static_cast<int>(i), allele, -1, merged);
                    continue;
                }
                update_read_var_profile_with_allele(vi > static_cast<int>(i) ? vi - 1 : vi,
                    vi == static_cast<int>(i) ? allele : allele_at(vi),
                    vi == static_cast<int>(i) ? -1 : profile.alt_qi[vi - profile.start_var_idx], merged);
            }
            profile = std::move(merged);
            if (allele >= 0) ++first.counts.alle_covs[allele];
        }
        update_variant_depth_fields(first);
        vars.erase(vars.begin() + i + 1);
        categories.erase(categories.begin() + i + 1);
    }
}

static void merge_msa_colocated_deletions(const PhasingChunk& chunk,
    std::vector<CandidateVariant>& vars,
    std::vector<VariantCategory>& categories, std::vector<ReadVariantProfile>& profiles) {
    for (size_t i = 0; i + 1 < vars.size(); ++i) {
        auto& first = vars[i];
        const auto& second = vars[i + 1];
        // A homozygous verdict on one of a co-located pair is the symptom, not a
        // reason to leave the pair alone: that record scores the other
        // haplotype's reads against its own allele, so its allele fraction
        // approaches 1 and the classifier calls it homozygous. Requiring both to
        // be heterozygous refuses exactly the loci that need merging, so one
        // heterozygous verdict between them is enough.
        const auto noisy = [](const CandidateVariant& v) {
            return v.counts.category == VariantCategory::NoisyCandHet ||
                   v.counts.category == VariantCategory::NoisyCandHom;
        };
        if (first.key.type != VariantType::Deletion || second.key.type != VariantType::Deletion ||
            first.key.pos != second.key.pos || first.key.ref_len == second.key.ref_len ||
            !first.msa_insertion_alts.empty() || !second.msa_insertion_alts.empty() ||
            !noisy(first) || !noisy(second) ||
            (first.counts.category != VariantCategory::NoisyCandHet &&
             second.counts.category != VariantCategory::NoisyCandHet)) continue;
        // Two deletions of different length at one position are two alleles of
        // one site, not two sites. Emitted separately, each scores the other
        // haplotype's reads against its own allele, so neither can express the
        // locus: at chr20:55,896,396 the 3 bp and 16 bp deletions are the
        // maternal and paternal alleles and no read is reference, which left the
        // 16 bp record at 0 reference / 30 alt, allele fraction 1.000.
        //
        // Merge them the way merge_msa_insertion_alleles merges co-located
        // insertions: the longer deletion's span becomes REF, ALT 1 deletes all
        // of it and ALT 2 retains the bases the shorter deletion leaves, so the
        // writer emits anchor + retained sequence for each allele.
        const int long_len = std::max(first.key.ref_len, second.key.ref_len);
        const int short_len = std::min(first.key.ref_len, second.key.ref_len);
        const bool first_is_long = first.key.ref_len == long_len;
        const size_t tail_beg = static_cast<size_t>(first.key.pos - chunk.ref_beg) +
                                static_cast<size_t>(short_len);
        const size_t tail_len = static_cast<size_t>(long_len - short_len);
        if (tail_beg + tail_len > chunk.ref_seq.size()) continue;
        const std::string retained = chunk.ref_seq.substr(tail_beg, tail_len);
        // Allele order follows record order, so a read's existing vote at each
        // record maps onto the merged allele without re-deriving it.
        first.key.ref_len = long_len;
        first.msa_insertion_alts = first_is_long
            ? std::vector<std::string>{std::string(), retained}
            : std::vector<std::string>{retained, std::string()};
        first.counts.alle_covs.assign(3, 0);
        first.counts.total_cov = 0;
        for (auto& profile : profiles) {
            if (profile.start_var_idx < 0) continue;
            const auto allele_at = [&](int index) {
                return index < profile.start_var_idx || index > profile.end_var_idx ? -1 :
                       profile.alleles[index - profile.start_var_idx];
            };
            const int a = allele_at(static_cast<int>(i)), b = allele_at(static_cast<int>(i + 1));
            const int allele = a == 1 && b != 1 ? 1 : b == 1 && a != 1 ? 2 :
                               a == 0 && b == 0 ? 0 : -1;
            ReadVariantProfile merged;
            merged.read_id = profile.read_id;
            for (int vi = profile.start_var_idx; vi <= profile.end_var_idx; ++vi) {
                if (vi == static_cast<int>(i + 1)) {
                    if (profile.start_var_idx == vi)
                        update_read_var_profile_with_allele(static_cast<int>(i), allele, -1, merged);
                    continue;
                }
                update_read_var_profile_with_allele(vi > static_cast<int>(i) ? vi - 1 : vi,
                    vi == static_cast<int>(i) ? allele : allele_at(vi),
                    vi == static_cast<int>(i) ? -1 : profile.alt_qi[vi - profile.start_var_idx], merged);
            }
            profile = std::move(merged);
            if (allele >= 0) ++first.counts.alle_covs[allele];
        }
        update_variant_depth_fields(first);
        vars.erase(vars.begin() + i + 1);
        categories.erase(categories.begin() + i + 1);
        // A third co-located deletion can remain, and it duplicates an allele the
        // merged record now carries. Left in place it is classified from its own
        // double-counted evidence and emitted as a second record at the same
        // position: at chr20:55,883,019 the merged record reads
        // AATATAT -> AAT,A at 1|2, matching the competitor, while a leftover 4 bp
        // record reported 1 reference against 36 alt at allele fraction 0.973 and
        // was emitted 1|1 -- a homozygous call at a locus whose two haplotypes
        // are the alleles beside it. Demote the leftover so the prune pass drops
        // it; its allele is not lost, it is allele 1 or 2 of the merged record.
        for (size_t j = 0; j < vars.size(); ++j) {
            if (j == i) continue;
            CandidateVariant& other = vars[j];
            if (other.key.type != VariantType::Deletion ||
                other.key.pos != first.key.pos ||
                !other.msa_insertion_alts.empty()) continue;
            const int retained = first.key.ref_len - other.key.ref_len;
            bool duplicate = retained == 0;
            for (const std::string& allele : first.msa_insertion_alts)
                if (retained > 0 && static_cast<size_t>(retained) == allele.size())
                    duplicate = true;
            if (!duplicate) continue;
            other.counts.category = VariantCategory::LowCoverage;
            if (j < categories.size()) categories[j] = VariantCategory::LowCoverage;
        }
    }
}

static void refresh_assigned_msa_observations(const Options& opts,
    const std::array<int, 2>& clu_n_seqs,
    const std::array<std::vector<int>, 2>& clu_read_ids,
    const std::array<std::vector<AlnStr>, 2>& alignments, hts_pos_t ref_beg,
    std::vector<CandidateVariant>& vars, std::vector<ReadVariantProfile>& profiles);

int make_vars_from_msa_cons_aln(
    const Options& opts, PhasingChunk& chunk,
    int /*n_noisy_reads*/, const std::vector<int>& /*read_ids*/,
    hts_pos_t noisy_reg_beg,
    int n_cons,
    const std::array<int, 2>& clu_n_seqs,
    const std::array<std::vector<int>, 2>& clu_read_ids,
    const std::array<std::vector<AlnStr>, 2>& aln_strs,
    std::vector<CandidateVariant>& noisy_vars,
    std::vector<VariantCategory>& noisy_var_cate,
    std::vector<ReadVariantProfile>& noisy_rvp) {
    noisy_vars.clear();
    noisy_var_cate.clear();
    noisy_rvp.clear();

    // Branch on the number of consensus sequences from MSA:
    // returned by `collect_noisy_reg_aln_strs`, not on whether auxiliary vectors are empty.
    if (n_cons == 0) return 0;

    if (n_cons == 1) {
        const std::vector<AlnStr>& clu_aln_strs = aln_strs[0];
        if (clu_aln_strs.empty()) return 0;
        const AlnStr& ref_cons_aln_str = clu_aln_strs[0];
        noisy_vars = make_cand_vars_from_msa(opts, chunk, noisy_reg_beg,
                                             ref_cons_aln_str.target_aln,
                                             ref_cons_aln_str.query_aln,
                                             ref_cons_aln_str.aln_len,
                                             false);
        if (noisy_vars.empty()) return 0;
        noisy_var_cate.assign(noisy_vars.size(), VariantCategory::NoisyCandHom);
        for (CandidateVariant& var : noisy_vars) var.counts.category = VariantCategory::NoisyCandHom;

        noisy_rvp = init_read_profiles(chunk.reads.size());
        update_cand_var_profile_from_cons_aln_str1(
            clu_n_seqs[0], clu_read_ids[0], clu_aln_strs,
            noisy_reg_beg, noisy_vars, noisy_rvp);
        for (CandidateVariant& var : noisy_vars) update_variant_depth_fields(var);
        return static_cast<int>(noisy_vars.size());
    }

    std::vector<CandidateVariant> hap1_vars;
    std::vector<CandidateVariant> hap2_vars;
    if (!aln_strs[0].empty()) {
        const AlnStr& ref_cons_aln_str = aln_strs[0][0];
        hap1_vars = make_cand_vars_from_msa(opts, chunk, noisy_reg_beg,
                                            ref_cons_aln_str.target_aln,
                                            ref_cons_aln_str.query_aln,
                                            ref_cons_aln_str.aln_len,
                                            false);
    }
    if (!aln_strs[1].empty()) {
        const AlnStr& ref_cons_aln_str = aln_strs[1][0];
        hap2_vars = make_cand_vars_from_msa(opts, chunk, noisy_reg_beg,
                                            ref_cons_aln_str.target_aln,
                                            ref_cons_aln_str.query_aln,
                                            ref_cons_aln_str.aln_len,
                                            false);
    }
    // Asked for by name by the retry (retry_opts.force_noisy_msa). One tandem
    // repeat can be described by two nested deletion forms; without this split
    // both forms are emitted as independent heterozygotes at the same locus.
    // The retry sets this field rather than clearing a broader flag precisely so
    // that this guard stays on: an earlier version cleared the whole recovery
    // flag to reach the MSA and silently disabled this split as a side effect.
    if (opts.force_noisy_msa)
        split_nested_msa_deletions(opts, chunk, hap1_vars, hap2_vars);
    update_cand_var_profile_from_cons_aln_str2(
        opts, chunk, clu_n_seqs, clu_read_ids, aln_strs, noisy_reg_beg,
        hap1_vars, hap2_vars, noisy_vars, noisy_var_cate, noisy_rvp);
    // A locus whose two haplotypes both differ from the reference has to be
    // emitted once with both alleles, in every arm. While this merge was gated
    // on the recovery pass it never ran by default, so the default pipeline
    // described such a locus as two competing biallelic records, each scoring
    // the other haplotype's reads against its own allele.
    if (opts.merge_colocated_msa_alleles)
        merge_msa_insertion_alleles(noisy_vars, noisy_var_cate, noisy_rvp);
    if (opts.merge_colocated_msa_alleles)
        merge_msa_colocated_deletions(chunk, noisy_vars, noisy_var_cate, noisy_rvp);
    // An assigned read's allele at an MSA site is read from its own cluster
    // alignment, in every arm: these counts describe the reads. While this
    // refresh was gated on the recovery pass, the alignment-only channel was
    // left with the pre-refresh counts, which are wrong -- with two clusters of
    // two reads each carrying their own consensus the site is 2 ref / 2 alt, and
    // without the refresh it was counted 3 ref / 1 alt. It is also why eight of
    // the seventy-six candidates shared with the hybrid arm on
    // chr20:48,176,830-48,229,446 carried different counts.
    if (opts.refresh_msa_observations && !aln_strs[0].empty() && !aln_strs[1].empty())
        refresh_assigned_msa_observations(opts, clu_n_seqs, clu_read_ids, aln_strs,
                                           noisy_reg_beg, noisy_vars, noisy_rvp);
    return static_cast<int>(noisy_vars.size());
}

// ════════════════════════════════════════════════════════════════════════════
// merge_var_profile
// ════════════════════════════════════════════════════════════════════════════

int merge_var_profile(PhasingChunk& chunk,
                      const std::vector<CandidateVariant>& noisy_vars,
                      const std::vector<VariantCategory>& noisy_var_cate,
                      const std::vector<ReadVariantProfile>& noisy_rvp,
                      const VariantKeySet* site_whitelist,
                      bool admit_all_in_region,
                      const VariantKeySet* replace_sites) {
    if (noisy_vars.empty()) return 0;
    // Region-trust mode: collect_noisy_vars1 already restricted which noisy
    // regions ran MSA to those overlapping a whitelisted window (the exact
    // position inside the window is not knowable in advance -- that is the
    // whole reason MSA is being asked to find it). Requiring calls to also
    // land on an exact whitelist key defeats that purpose; here the containment
    // gate already did the trust decision, so admit every call.
    if (admit_all_in_region) site_whitelist = nullptr;

    std::vector<CandidateVariant> new_vars = noisy_vars;
    std::vector<VariantCategory> new_cats = noisy_var_cate;
    if (new_cats.size() != new_vars.size()) {
        new_cats.resize(new_vars.size());
        for (size_t i = 0; i < new_vars.size(); ++i)
            new_cats[i] = new_vars[i].counts.category;
    }

    std::vector<size_t> order(new_vars.size());
    std::iota(order.begin(), order.end(), 0);
    std::stable_sort(order.begin(), order.end(), [&](size_t a, size_t b) {
        return exact_comp_cand_var(&new_vars[a], &new_vars[b]) < 0;
    });
    std::vector<CandidateVariant> sorted_new;
    std::vector<VariantCategory> sorted_cats;
    sorted_new.reserve(new_vars.size());
    sorted_cats.reserve(new_cats.size());
    for (size_t idx : order) {
        sorted_new.push_back(std::move(new_vars[idx]));
        sorted_cats.push_back(new_cats[idx]);
    }
    new_vars = std::move(sorted_new);
    new_cats = std::move(sorted_cats);

    const CandidateTable old_vars = chunk.candidates;
    const std::vector<ReadVariantProfile> old_profiles = chunk.read_var_profile.empty()
                                                             ? init_read_profiles(chunk.reads.size())
                                                             : chunk.read_var_profile;

    CandidateTable merged_vars;
    merged_vars.reserve(old_vars.size() + new_vars.size());
    std::vector<int> old_to_merged(old_vars.size(), -1);
    std::vector<int> new_to_merged(new_vars.size(), -1);

    size_t old_i = 0;
    size_t new_i = 0;
    int admitted = 0;
    while (old_i < old_vars.size() && new_i < new_vars.size()) {
        const int ret = exact_comp_var_site(&old_vars[old_i].key, &new_vars[new_i].key);
        if (ret < 0) {
            old_to_merged[old_i] = static_cast<int>(merged_vars.size());
            merged_vars.push_back(old_vars[old_i++]);
        } else if (ret > 0) {
            if ((site_whitelist != nullptr &&
                 site_whitelist->find(new_vars[new_i].key) == site_whitelist->end())) {
                ++new_i;
                continue;
            }
            set_noisy_category(new_vars[new_i], new_cats[new_i]);
            restore_stashed_initial_if_any(chunk, new_vars[new_i]);
            new_to_merged[order[new_i]] = static_cast<int>(merged_vars.size());
            merged_vars.push_back(new_vars[new_i++]);
            ++admitted;
        } else {
            // A RepeatHetIndel is screened out of k-means by construction, so it
            // is never assigned a haplotype and never emitted: chr20:48,177,726
            // (INS, DP 74, 35/39, AF 0.527) and chr20:48,234,101 (INS, DP 57,
            // 23/34, AF 0.596) both sit in the candidate table with
            // HAP_ALT = HAP_REF = 0 and PS = 0, and hiphase phases both. Such a
            // candidate carries no phase information, so the MSA's own call at
            // the same key strictly dominates it -- and that call is the
            // verified one, which is the evidence this pipeline is supposed to
            // admit. Requiring region-trust mode or a whitelist hit before the
            // swap meant a plain run always kept the screened, unusable version.
            const bool replace_repeat =
                old_vars[old_i].counts.category == VariantCategory::RepeatHetIndel;
            const bool replace_selected = replace_sites != nullptr &&
                replace_sites->count(new_vars[new_i].key) != 0;
            // An old candidate in a pruned category is about to be deleted by
            // prune_not_candidate_variants, so keeping it in preference to the
            // MSA's own call at the same key discards the only version of the
            // site that will survive -- and the MSA's version is the verified
            // one, since `msa_verified` is set where the MSA CONSTRUCTS a
            // candidate. That is how chr20:48,173,317 was lost: the catalog
            // claims the locus, classify_graph_only_candidates judges it by the
            // graph window (AF 0.652 against [0.39, 0.61]) and marks it
            // LowCoverage, and the MSA call that collect-bam-variation emits as
            // a phased het was then thrown away in favour of it. Replacing a
            // pruned candidate can only add information: the alternative is no
            // site at all.
            const bool replace_pruned =
                (old_vars[old_i].counts.category == VariantCategory::LowCoverage ||
                 old_vars[old_i].counts.category == VariantCategory::LowAlleleFraction);
            // A record carrying both of the locus' alleles describes a locus
            // whose haplotypes are those two alleles. A single-allele record at
            // the same key cannot describe it: whichever allele it names, the
            // other haplotype's reads are scored against that one, so the depth
            // collapses and the allele fraction runs to 1. Measured at
            // chr20:55,896,396, a 16 bp homopolymer deletion the competitor
            // phases and read truth confirms at purity 1.000 over 28 reads: the
            // alignment channel merges it to REF=T(16) with alleles retaining 13
            // and 0 bases at DP 69, and the catalog claims the locus, so the
            // hybrid kept the graph-claimed single-allele version at DP 37, 2/35,
            // allele fraction 0.946, classified CleanHom -- the one candidate
            // injection lost across the whole panel. The claim is still honoured:
            // graph_site and the initial category carry over as for any other
            // replacement. Only the allele set is preserved.
            const bool replace_single_allele =
                new_vars[new_i].msa_insertion_alts.size() == 2 &&
                old_vars[old_i].msa_insertion_alts.empty();
            if (replace_repeat || replace_selected || replace_pruned ||
                replace_single_allele) {
                set_noisy_category(new_vars[new_i], new_cats[new_i]);
                new_vars[new_i].graph_site |= old_vars[old_i].graph_site;
                new_vars[new_i].counts.candvarcate_initial =
                    old_vars[old_i].counts.candvarcate_initial;
                new_to_merged[order[new_i]] = static_cast<int>(merged_vars.size());
                merged_vars.push_back(new_vars[new_i]);
                ++admitted;
            } else {
                old_to_merged[old_i] = static_cast<int>(merged_vars.size());
                merged_vars.push_back(old_vars[old_i]);
            }
            ++old_i;
            ++new_i;
        }
    }
    while (old_i < old_vars.size()) {
        old_to_merged[old_i] = static_cast<int>(merged_vars.size());
        merged_vars.push_back(old_vars[old_i++]);
    }
    while (new_i < new_vars.size()) {
        if ((site_whitelist != nullptr &&
             site_whitelist->find(new_vars[new_i].key) == site_whitelist->end())) {
            ++new_i;
            continue;
        }
        set_noisy_category(new_vars[new_i], new_cats[new_i]);
        restore_stashed_initial_if_any(chunk, new_vars[new_i]);
        new_to_merged[order[new_i]] = static_cast<int>(merged_vars.size());
        merged_vars.push_back(new_vars[new_i++]);
        ++admitted;
    }

    std::vector<ReadVariantProfile> merged_profiles = init_read_profiles(chunk.reads.size());
    const std::vector<ReadVariantProfile> empty_new_profiles =
        noisy_rvp.empty() ? init_read_profiles(chunk.reads.size()) : std::vector<ReadVariantProfile>{};
    const std::vector<ReadVariantProfile>& new_profiles = noisy_rvp.empty() ? empty_new_profiles : noisy_rvp;

    const int n_reads = static_cast<int>(chunk.reads.size());
    for (int ord = 0; ord < (chunk.ordered_read_ids.empty() ? n_reads : static_cast<int>(chunk.ordered_read_ids.size())); ++ord) {
        const int read_i = chunk.ordered_read_ids.empty()
                               ? ord
                               : chunk.ordered_read_ids[static_cast<size_t>(ord)];
        if (read_i < 0 || read_i >= n_reads) continue;
        if (chunk.reads[static_cast<size_t>(read_i)].is_skipped) continue;
        const ReadVariantProfile* old_p = static_cast<size_t>(read_i) < old_profiles.size()
                                              ? &old_profiles[static_cast<size_t>(read_i)]
                                              : nullptr;
        const ReadVariantProfile* new_p = static_cast<size_t>(read_i) < new_profiles.size()
                                              ? &new_profiles[static_cast<size_t>(read_i)]
                                              : nullptr;
        merge_read_var_profile_entries(old_p, old_to_merged,
                                       new_p, new_to_merged,
                                       merged_profiles[static_cast<size_t>(read_i)],
                                       static_cast<int>(merged_vars.size()));
    }

    // MSA replaces a working observation, not its original BAM/GAF history.
    for (size_t ri = 0; ri < merged_profiles.size() && ri < old_profiles.size(); ++ri) {
        auto& dest = merged_profiles[ri];
        const auto& source = old_profiles[ri];
        dest.bam_alleles.assign(dest.alleles.size(), -1);
        dest.bam_qi.assign(dest.alleles.size(), -1);
        for (size_t pi = 0; pi < source.alleles.size(); ++pi) {
            const int vi = source.start_var_idx + static_cast<int>(pi);
            if (vi < 0 || static_cast<size_t>(vi) >= old_vars.size()) continue;
            const auto found = std::lower_bound(merged_vars.begin(), merged_vars.end(), old_vars[vi].key,
                [](const CandidateVariant& v, const VariantKey& key) {
                    return exact_comp_var_site(&v.key, &key) < 0;
                });
            if (found == merged_vars.end() || exact_comp_var_site(&found->key, &old_vars[vi].key)) continue;
            const int di = static_cast<int>(found - merged_vars.begin()) - dest.start_var_idx;
            if (di < 0 || static_cast<size_t>(di) >= dest.alleles.size()) continue;
            if (pi < source.bam_alleles.size()) dest.bam_alleles[di] = source.bam_alleles[pi];
            if (pi < source.bam_qi.size()) dest.bam_qi[di] = source.bam_qi[pi];
            if (pi < source.graph_alleles.size()) {
                dest.graph_alleles.resize(dest.alleles.size(), -1);
                dest.graph_alleles[di] = source.graph_alleles[pi];
            }
        }
    }

    cgranges_t* merged_cr = cr_init();
    for (int ord = 0; ord < (chunk.ordered_read_ids.empty() ? n_reads : static_cast<int>(chunk.ordered_read_ids.size())); ++ord) {
        const int read_i = chunk.ordered_read_ids.empty()
                               ? ord
                               : chunk.ordered_read_ids[static_cast<size_t>(ord)];
        if (read_i < 0 || read_i >= n_reads) continue;
        if (chunk.reads[static_cast<size_t>(read_i)].is_skipped) continue;
        const ReadVariantProfile& p = merged_profiles[static_cast<size_t>(read_i)];
        if (p.start_var_idx < 0 || p.end_var_idx < p.start_var_idx) continue;
        cr_add(merged_cr, "cr", p.start_var_idx, p.end_var_idx + 1,
               static_cast<int32_t>(read_i));
    }
    cr_index(merged_cr);

    chunk.candidates = std::move(merged_vars);
    chunk.read_var_profile = std::move(merged_profiles);
    chunk.read_var_cr.reset(merged_cr);
    return admitted;
}

// ════════════════════════════════════════════════════════════════════════════
// collect_noisy_vars1
// ════════════════════════════════════════════════════════════════════════════

struct MsaSiteSlice {
    std::string ref, query;
    std::array<std::string, 2> flank_ref, flank_query;
    bool covered = false;
};

static MsaSiteSlice slice_msa_site(const AlnStr& aln, const VariantKey& key,
                                   hts_pos_t ref_beg) {
    constexpr int kSiteFlankBases = 3;
    const int beg = static_cast<int>(key.pos - ref_beg);
    const int end = beg + (key.type == VariantType::Insertion ? 0 : key.ref_len);
    int ref_pos = -1;
    MsaSiteSlice result;
    for (int i = 0; i < aln.aln_len; ++i) {
        const auto t = aln.target_aln[i], q = aln.query_aln[i];
        if (t != kGapBase) ++ref_pos;
        const int pos = ref_pos + (t == kGapBase ? 1 : 0);
        if (pos < beg - kSiteFlankBases || pos >= end + kSiteFlankBases) continue;
        if (i < aln.target_beg || i > aln.target_end ||
            i < aln.query_beg || i > aln.query_end ||
            (t > 3 && t != kGapBase) || (q > 3 && q != kGapBase)) return result;
        const bool at_site = (pos >= beg && pos < end) ||
            (key.type == VariantType::Insertion && pos == beg && t == kGapBase);
        const int side = pos < beg ? 0 : 1;
        auto& ref = at_site ? result.ref : result.flank_ref[side];
        auto& query = at_site ? result.query : result.flank_query[side];
        if (t != kGapBase) ref.push_back(nt4_to_base(t));
        if (q != kGapBase) query.push_back(nt4_to_base(q));
    }
    result.covered = result.flank_ref[0].size() == kSiteFlankBases &&
                     result.flank_ref[1].size() == kSiteFlankBases;
    return result;
}

static int msa_site_event_allele(const MsaSiteSlice& site, const VariantKey& key,
                                 const std::array<MsaSiteSlice, 2>& context,
                                 const std::vector<std::string>* insertion_alts = nullptr) {
    if (insertion_alts != nullptr && !insertion_alts->empty()) {
        // What the reference looks like at a multi-allele site depends on the
        // event: an insertion site carries no query sequence when the read has
        // no insertion, while a deletion site carries the whole footprint when
        // the read deletes nothing. Treating an empty query as reference is
        // right for insertions and wrong for deletions, where it is the allele
        // that removes the entire span -- at chr20:55,896,395 that scored 34
        // reads carrying the 16 bp deletion as reference and left the allele
        // itself at zero support.
        const bool deletion = key.type == VariantType::Deletion;
        if (deletion ? site.query == site.ref : site.query.empty()) return 0;
        for (size_t ai = 0; ai < insertion_alts->size(); ++ai)
            if (site.query == (*insertion_alts)[ai]) return static_cast<int>(ai + 1);
        return -1;
    }
    const std::string alt = key.type == VariantType::Deletion ? std::string() : key.alt;
    if (site.query == site.ref) return 0;
    if (site.query == alt) return 1;
    // Separate insertion rows use zero for absence of their own ALT. When
    // both haplotypes insert, literal REF cannot represent the other allele.
    // Project only an exact, complete diploid contrast with common flanks;
    // unrelated insertion sequences retain their unknown allele identity.
    if (key.type == VariantType::Insertion &&
        context[0].covered && context[1].covered &&
        context[0].ref == site.ref && context[1].ref == site.ref &&
        context[0].flank_query == context[1].flank_query &&
        ((context[0].query == alt) != (context[1].query == alt))) {
        const auto& other = context[context[0].query == alt ? 1 : 0];
        if (site.query == other.query) return 0;
    }
    // The non-deleted haplotype may carry a verified substitution within the
    // deletion footprint. Require its exact consensus sequence, not just length.
    if (key.type == VariantType::Deletion && site.query.size() == site.ref.size())
        for (const auto& consensus : context)
            if (consensus.covered && consensus.ref == site.ref &&
                consensus.query == site.query) return 0;
    // Zero is absence of this ALT, not necessarily literal reference. The
    // other fixed haplotype may insert bases inside this deletion footprint.
    // Require the target edit in exactly one complete consensus and an exact
    // match to the other; a third sequence or ambiguous contrast stays missing.
    if (key.type == VariantType::Deletion &&
        context[0].covered && context[1].covered &&
        context[0].ref == site.ref && context[1].ref == site.ref &&
        context[0].query.empty() != context[1].query.empty()) {
        const auto& other = context[context[0].query.empty() ? 1 : 0];
        if (other.query.size() > site.ref.size() && site.query == other.query)
            return 0;
    }
    return -1;
}

static int call_msa_site_with_context(const std::array<AlnStr, 2>& alignments,
                                      const VariantKey& key, hts_pos_t ref_beg,
                                      const std::array<MsaSiteSlice, 2>& context,
                                      const std::vector<std::string>* insertion_alts = nullptr) {
    int first = -1;
    for (int ci = 0; ci < 2; ++ci) {
        const auto site = slice_msa_site(alignments[ci], key, ref_beg);
        if (!site.covered) return -1;
        for (int side = 0; side < 2; ++side) {
            bool supported = site.flank_query[side] == site.flank_ref[side];
            for (const auto& consensus : context)
                supported |= consensus.covered &&
                    site.flank_ref[side] == consensus.flank_ref[side] &&
                    site.flank_query[side] == consensus.flank_query[side];
            if (!supported) return -1;
        }
        const int allele = msa_site_event_allele(site, key, context, insertion_alts);
        if (allele < 0 || (ci == 1 && allele != first)) return -1;
        first = allele;
    }
    return first;
}

int call_msa_site_allele(const std::array<AlnStr, 2>& alignments,
                         const VariantKey& key, hts_pos_t ref_beg,
                         const std::array<AlnStr, 2>* consensuses) {
    std::array<MsaSiteSlice, 2> context;
    if (consensuses != nullptr)
        for (int ci = 0; ci < 2; ++ci) context[ci] = slice_msa_site((*consensuses)[ci], key, ref_beg);
    return call_msa_site_with_context(alignments, key, ref_beg, context);
}

static int local_allele_distance(const std::string& first, const std::string& second) {
    constexpr int kBeyondOneEdit = 2;
    if (first == second) return 0;
    if (first.size() + 1 < second.size() || second.size() + 1 < first.size()) return kBeyondOneEdit;
    size_t i = 0, j = 0;
    int edits = 0;
    while (i < first.size() && j < second.size()) {
        if (first[i] == second[j]) { ++i; ++j; continue; }
        if (++edits > 1) return kBeyondOneEdit;
        if (first.size() >= second.size()) ++i;
        if (second.size() >= first.size()) ++j;
    }
    return edits + (i < first.size() || j < second.size());
}

static constexpr int kMaxLocalAlleleEdits = 1;

static int call_local_msa_allele(const AlnStr& read, const VariantKey& key,
                                     hts_pos_t ref_beg, const std::array<AlnStr, 2>& consensuses,
                                     const std::vector<std::string>* insertion_alts = nullptr) {
    const std::array<MsaSiteSlice, 2> context = {
        slice_msa_site(consensuses[0], key, ref_beg), slice_msa_site(consensuses[1], key, ref_beg)};
    const int exact = call_msa_site_with_context({read, read}, key, ref_beg, context, insertion_alts);
    if (exact >= 0) return exact;
    // A DEL/INS contrast contains two different edits. Recognizing the exact
    // alternate-absent insertion must not also enable the one-error fallback
    // that previously rejected this context. Its compound alleles need the
    // complete exact sequence and flanks checked above.
    if (key.type == VariantType::Deletion &&
        ((context[0].query.empty() && context[1].query.size() > context[1].ref.size()) ||
         (context[1].query.empty() && context[0].query.size() > context[0].ref.size())))
        return -1;
    const auto observed = slice_msa_site(read, key, ref_beg);
    if (!observed.covered || !context[0].covered || !context[1].covered ||
        context[0].flank_query != context[1].flank_query) return -1;
    std::array<int, 2> alleles, distances;
    const std::string query = observed.flank_query[0] + observed.query + observed.flank_query[1];
    for (int ci = 0; ci < 2; ++ci) {
        alleles[ci] = msa_site_event_allele(context[ci], key, context, insertion_alts);
        const std::string expected = context[ci].flank_query[0] + context[ci].query + context[ci].flank_query[1];
        distances[ci] = local_allele_distance(query, expected);
    }
    if (alleles[0] < 0 || alleles[1] < 0 || alleles[0] == alleles[1] ||
        distances[0] == distances[1]) return -1;
    const int best = distances[0] < distances[1] ? 0 : 1;
    // Permit one local sequencing error only when the two fixed haplotypes
    // differ solely at this site and one allele is strictly closer. Cluster
    // membership itself never supplies the answer.
    return distances[best] <= kMaxLocalAlleleEdits ? alleles[best] : -1;
}

static void refresh_assigned_msa_observations(const Options& opts,
    const std::array<int, 2>& clu_n_seqs,
    const std::array<std::vector<int>, 2>& clu_read_ids,
    const std::array<std::vector<AlnStr>, 2>& alignments, hts_pos_t ref_beg,
    std::vector<CandidateVariant>& vars, std::vector<ReadVariantProfile>& profiles) {
    constexpr int kLeftGapAlignment = 1;
    const std::array<AlnStr, 2> consensuses = {alignments[0][0], alignments[1][0]};

    // call_local_msa_allele requires a consistent, valid allele call from
    // BOTH haplotype consensuses before it will attribute any read at all --
    // true whenever the two haplotypes carry the same kind of event (the
    // ordinary case this refresh targets: fixing cross-cluster count
    // contamination), but not when a homopolymer is unstable enough that the
    // two consensuses disagree on event TYPE, not just allele (one an
    // insertion, the other a deletion, relative to reference, at the same
    // locus). When that happens every read fails classification for that
    // variant specifically.
    //
    // The refresh below still needs to run per-read (skipping only the reads
    // that individually fail, the way it always has) for every variant it CAN
    // say something about -- that per-read skip is what actually fixes the
    // contamination this function exists for, and doing it selectively here
    // instead would just reintroduce that contamination for every variant
    // with a single stray unclassifiable read, which is most of them. The
    // one thing a variant with NO valid classification anywhere is owed is
    // being left alone entirely: this pre-pass identifies those variants so
    // the write loop below can skip them outright, rather than resetting
    // their counts and read entries to nothing. Measured at chr20:411,654 (a
    // homopolymer T-insertion candidate whose sibling haplotype consensus is
    // a deletion at the same run): unconditional reset-and-rebuild left the
    // variant with zero reads at every position (real DP 51, matching
    // longcallD's own DP 51 at the same locus), and the site silently
    // vanished from the phased VCF. A per-read "skip on failure" instead of
    // this per-variant one measurably regresses read-phasing accuracy
    // chromosome-wide (evaluations/2026-09-20-refresh-observations-bug/):
    // it reintroduces the contamination on every partially-classifiable
    // variant to fix the rare fully-unclassifiable one.
    std::vector<bool> any_valid_call(vars.size(), false);
    for (int ci = 0; ci < 2; ++ci) {
        for (int ri = 0; ri < clu_n_seqs[ci]; ++ri) {
            const auto& cons_read = alignments[ci][2 * ri + 1];
            if (cons_read.target_beg != 0 || cons_read.query_beg != 0 ||
                cons_read.target_end != cons_read.aln_len - 1 ||
                cons_read.query_end != cons_read.aln_len - 1) continue;
            AlnStr ref_read;
            make_ref_read_aln_str(opts, consensuses[ci], cons_read, ref_read);
            ref_read.target_beg = ref_read.query_beg = 0;
            ref_read.target_end = ref_read.query_end = ref_read.aln_len - 1;
            if (opts.gap_aln == kLeftGapAlignment) left_normalize_msa_alignment(ref_read);
            for (size_t vi = 0; vi < vars.size(); ++vi) {
                if (vars[vi].counts.category != VariantCategory::NoisyCandHet) continue;
                if (any_valid_call[vi]) continue;
                if (call_local_msa_allele(ref_read, vars[vi].key, ref_beg, consensuses,
                                          &vars[vi].msa_insertion_alts) >= 0)
                    any_valid_call[vi] = true;
            }
        }
    }

    for (int ci = 0; ci < 2; ++ci) {
        for (int ri = 0; ri < clu_n_seqs[ci]; ++ri) {
            const auto& cons_read = alignments[ci][2 * ri + 1];
            // Preserve existing partial-read coverage handling. Full reads can
            // be checked directly against reference without guessing an allele
            // from their whole-consensus cluster membership.
            if (cons_read.target_beg != 0 || cons_read.query_beg != 0 ||
                cons_read.target_end != cons_read.aln_len - 1 ||
                cons_read.query_end != cons_read.aln_len - 1) continue;
            AlnStr ref_read;
            make_ref_read_aln_str(opts, consensuses[ci], cons_read, ref_read);
            ref_read.target_beg = ref_read.query_beg = 0;
            ref_read.target_end = ref_read.query_end = ref_read.aln_len - 1;
            if (opts.gap_aln == kLeftGapAlignment) left_normalize_msa_alignment(ref_read);
            const int rid_ = clu_read_ids[ci][ri];
            if (rid_ < 0 || static_cast<size_t>(rid_) >= profiles.size()) {
                if (getenv("PGPHASE_UBPROBE") != nullptr)
                    fprintf(stderr, "UB profiles.size=%zu index=%d clu=%d/%d\n",
                            profiles.size(), rid_, ci, ri);
                continue;
            }
            auto& profile = profiles[static_cast<size_t>(rid_)];
            for (size_t vi = 0; vi < vars.size(); ++vi) {
                if (vars[vi].counts.category != VariantCategory::NoisyCandHet) continue;
                if (!any_valid_call[vi]) continue;
                const int allele = call_local_msa_allele(ref_read, vars[vi].key, ref_beg, consensuses, &vars[vi].msa_insertion_alts);
                update_read_var_profile_with_allele(static_cast<int>(vi), allele, -1, profile);
            }
        }
    }
    for (size_t vi = 0; vi < vars.size(); ++vi) {
        auto& var = vars[vi];
        if (var.counts.category != VariantCategory::NoisyCandHet) continue;
        if (!any_valid_call[vi]) continue;
        var.counts.alle_covs.assign(var.msa_insertion_alts.empty() ? 2 : var.msa_insertion_alts.size() + 1, 0);
        for (const auto& profile : profiles) {
            if (static_cast<int>(vi) < profile.start_var_idx || static_cast<int>(vi) > profile.end_var_idx) continue;
            const int allele = profile.alleles[vi - profile.start_var_idx];
            if (allele >= 0 && static_cast<size_t>(allele) < var.counts.alle_covs.size())
                ++var.counts.alle_covs[allele];
        }
        var.counts.total_cov = std::accumulate(var.counts.alle_covs.begin(), var.counts.alle_covs.end(), 0);
        update_variant_depth_fields(var);
    }
}

void add_msa_site_observations(const Options& opts,
                                const std::vector<UnassignedMsaRead>& reads,
                                hts_pos_t ref_beg,
                                std::vector<CandidateVariant>& vars,
                                std::vector<ReadVariantProfile>& profiles,
                                const std::array<AlnStr, 2>* consensuses) {
    for (size_t vi = 0; vi < vars.size(); ++vi) {
        auto& var = vars[vi];
        // Missing reference observations can make a noisy het appear hom.
        // Recount both categories before applying the allele-fraction gate.
        const bool was_hom = var.counts.category == VariantCategory::NoisyCandHom;
        if (var.counts.category != VariantCategory::NoisyCandHet && !was_hom) continue;
        std::vector<std::pair<int, int>> observations;
        auto counts = var.counts.alle_covs;
        std::array<MsaSiteSlice, 2> context;
        if (consensuses != nullptr)
            for (int ci = 0; ci < 2; ++ci)
                context[ci] = slice_msa_site((*consensuses)[ci], var.key, ref_beg);
        for (const auto& read : reads) {
            int allele = call_msa_site_with_context(read.ref_read, var.key, ref_beg, context, &var.msa_insertion_alts);
            if (allele < 0 && consensuses != nullptr) {
                const int first = call_local_msa_allele(read.ref_read[0], var.key, ref_beg, *consensuses, &var.msa_insertion_alts);
                const int second = call_local_msa_allele(read.ref_read[1], var.key, ref_beg, *consensuses, &var.msa_insertion_alts);
                // Whole-window assignment is unnecessary, but both independently
                // composed paths must favor the same local allele.
                if (first >= 0 && first == second) allele = first;
            }
            if (allele < 0) continue;
            observations.emplace_back(read.read_id, allele);
            ++counts[allele];
        }
        if (observations.empty()) continue;
        // Two different consensuses alone do not establish heterozygosity.
        // Do not extend a site whose expanded evidence fails the existing AF gate.
        const int total = std::accumulate(counts.begin(), counts.end(), 0);
        const double af = static_cast<double>(counts[1]) / total;
        if (af < opts.min_af || af > opts.max_af) continue;
        if (!var.msa_insertion_alts.empty()) {
            const double second_af = static_cast<double>(counts[2]) / total;
            if (second_af < opts.min_af || second_af > opts.max_af) continue;
        }
        // Correct depth and profiles without changing the site's category.
        for (const auto& [read_id, allele] : observations)
            update_read_var_profile_with_allele(static_cast<int>(vi), allele, -1,
                                                profiles[read_id]);
        var.counts.total_cov += static_cast<int>(observations.size());
        var.counts.alle_covs = std::move(counts);
        update_variant_depth_fields(var);
    }
}

// A longer deletion also satisfies shorter deletion windows at the same locus.
// Keep its read on the longest supported candidate so one event votes once.
void make_colocated_deletions_exclusive(std::vector<CandidateVariant>& vars,
                                               std::vector<ReadVariantProfile>& profiles) {
    std::map<hts_pos_t, std::vector<int>> by_pos;
    for (int vi = 0; vi < static_cast<int>(vars.size()); ++vi) {
        const CandidateVariant& v = vars[static_cast<size_t>(vi)];
        if (v.key.type == VariantType::Deletion && v.key.ref_len > 0)
            by_pos[v.key.pos].push_back(vi);
    }
    for (auto& [pos, group] : by_pos) {
        (void)pos;
        if (group.size() < 2) continue;
        for (ReadVariantProfile& pr : profiles) {
            if (pr.start_var_idx < 0) continue;
            std::vector<int> alt_hits;
            for (const int vi : group) {
                if (vi < pr.start_var_idx || vi > pr.end_var_idx) continue;
                if (pr.alleles[static_cast<size_t>(vi - pr.start_var_idx)] == 1)
                    alt_hits.push_back(vi);
            }
            if (alt_hits.size() < 2) continue;
            // The read is alt at several records for one event. Keep the longest
            // record it supports -- the read's own event reached that far, so a
            // shorter record is a prefix of the same deletion, not a separate
            // allele -- and make it reference at the rest.
            int keep = alt_hits.front();
            for (const int vi : alt_hits)
                if (vars[static_cast<size_t>(vi)].key.ref_len >
                    vars[static_cast<size_t>(keep)].key.ref_len) keep = vi;
            for (const int vi : alt_hits) {
                if (vi == keep) continue;
                pr.alleles[static_cast<size_t>(vi - pr.start_var_idx)] = 0;
                CandidateVariant& v = vars[static_cast<size_t>(vi)];
                if (v.counts.alle_covs.size() > 1 && v.counts.alle_covs[1] > 0) {
                    --v.counts.alle_covs[1];
                    ++v.counts.alle_covs[0];
                }
            }
        }
        for (const int vi : group) update_variant_depth_fields(vars[static_cast<size_t>(vi)]);
    }
}

void restrict_msa_observations_to_read_coverage(
        const std::vector<ReadRecord>& reads,
        std::vector<CandidateVariant>& vars,
        std::vector<ReadVariantProfile>& profiles) {
    for (CandidateVariant& var : vars) {
        std::fill(var.counts.alle_covs.begin(), var.counts.alle_covs.end(), 0);
        var.counts.total_cov = 0;
    }
    for (size_t ri = 0; ri < profiles.size(); ++ri) {
        ReadVariantProfile& profile = profiles[ri];
        if (profile.start_var_idx < 0) continue;
        const ReadRecord& read = reads[ri];
        for (size_t offset = 0; offset < profile.alleles.size(); ++offset) {
            const size_t ci = static_cast<size_t>(profile.start_var_idx) + offset;
            if (ci >= vars.size()) break;
            CandidateVariant& var = vars[ci];
            const VariantKey& key = var.key;
            // BAM bounds and keys are one-based. SNPs need the base itself;
            // indels also need a surviving base on each side of their event.
            const hts_pos_t first = key.pos - (key.type != VariantType::Snp);
            const hts_pos_t last = key.pos + key.ref_len - (key.type == VariantType::Snp);
            int& allele = profile.alleles[offset];
            if (first < read.beg || last > read.end) {
                allele = -1;
                if (offset < profile.alt_qi.size()) profile.alt_qi[offset] = -1;
            }
            if (allele < 0 || static_cast<size_t>(allele) >= var.counts.alle_covs.size())
                continue;
            ++var.counts.alle_covs[allele];
            ++var.counts.total_cov;
        }
    }
    for (CandidateVariant& var : vars) update_variant_depth_fields(var);
}

static int collect_noisy_vars1(PhasingChunk& chunk, const Options& opts, int noisy_reg_i,
                        const VariantKeySet* site_whitelist,
                        std::vector<DeferredMsaObservation>* deferred) {
    const Interval& reg = chunk.noisy_regions[static_cast<size_t>(noisy_reg_i)];
    // Enter with noisy_reg_beg/end from the interval tree.
    hts_pos_t noisy_reg_beg = reg.beg;
    hts_pos_t noisy_reg_end = reg.end;

    if (site_whitelist != nullptr) {
        bool contains_whitelisted_site = false;
        for (const VariantKey& key : *site_whitelist) {
            if (key.tid == chunk.region.tid &&
                key.pos >= noisy_reg_beg && key.pos <= noisy_reg_end) {
                contains_whitelisted_site = true;
                break;
            }
        }
        if (!contains_whitelisted_site) return 0;
    }

    // Skip regions longer than max_noisy_reg_len (return 0 = done, no vars).
    if (noisy_reg_end - noisy_reg_beg + 1 > static_cast<hts_pos_t>(opts.max_noisy_reg_len))
        return 0;

    // First clip beg/end via collect_reg_ref_bseq, then collect reads on the
    // finalized region bounds.
    std::vector<uint8_t> ref_seq =
        collect_reg_ref_bseq(chunk, noisy_reg_beg, noisy_reg_end);
    if (ref_seq.empty())
        return 0;

    const std::vector<int> read_ids =
        collect_noisy_reg_reads(chunk, noisy_reg_beg, noisy_reg_end);
    const int n_noisy_reads = static_cast<int>(read_ids.size());

    // Skip regions with more than max_noisy_reg_cov reads.
    if (n_noisy_reads > opts.max_noisy_reg_cov)
        return 0;

    // MSA + consensus via WFA2-lib / abPOA (align.cpp).
    std::array<int, 2> clu_n_seqs{};
    std::array<std::vector<int>, 2> clu_read_ids{};
    std::array<std::vector<AlnStr>, 2> aln_strs{};
    std::vector<UnassignedMsaRead> unassigned;
    const bool recall_unplaced = opts.add_unplaced_msa_observations ||
        opts.recall_unplaced_msa_insertions;
    const int n_cons = collect_noisy_reg_aln_strs(
        opts, chunk, noisy_reg_beg, noisy_reg_end,
        read_ids, ref_seq,
        clu_n_seqs, clu_read_ids, aln_strs,
        // A broad retry may extend clusters. Local insertion recall instead
        // composes both fixed paths and defers read-local allele observations.
        recall_unplaced ? &unassigned : nullptr);

    // n_cons == 0 → MSA could not resolve; return -1 so the
    // outer loop leaves this region undone and may retry if another region makes progress.
    if (n_cons == 0)
        return -1;

    std::vector<CandidateVariant>    noisy_vars;
    std::vector<VariantCategory>     noisy_var_cate;
    std::vector<ReadVariantProfile>  noisy_rvp;

    make_vars_from_msa_cons_aln(
        opts, chunk,
        n_noisy_reads, read_ids, noisy_reg_beg,
        n_cons, clu_n_seqs, clu_read_ids, aln_strs,
        noisy_vars, noisy_var_cate, noisy_rvp);

    // A real neighboring variant is valid flank sequence, not alignment noise.
    std::array<AlnStr, 2> consensuses;
    if (n_cons == 2) consensuses = {aln_strs[0][0], aln_strs[1][0]};
    if (opts.recall_unplaced_msa_insertions && !opts.add_unplaced_msa_observations &&
        deferred != nullptr && !unassigned.empty()) {
        // Discovery must finish on the original BAM evidence. Feeding these
        // calls back between MSA regions can change later consensuses and
        // remove an already discovered allele representation.
        std::vector<bool> eligible(noisy_vars.size(), false);
        for (size_t vi = 0; vi < noisy_vars.size(); ++vi) {
            const CandidateVariant& site = noisy_vars[vi];
            if (site.key.type != VariantType::Insertion || site.key.ref_len != 0 ||
                site.counts.category != VariantCategory::NoisyCandHet ||
                !site.msa_insertion_alts.empty() || site.key.alt.empty())
                continue;
            // The full neighboring blocks supply context, not permission to
            // change their read evidence at unrelated repeat contrasts.
            const hts_pos_t anchor = site.key.pos - 1;
            if (!std::any_of(opts.retry_windows.begin(), opts.retry_windows.end(),
                    [anchor](const auto& window) {
                        return window.first <= anchor && anchor <= window.second;
                    })) continue;
            eligible[vi] = std::any_of(noisy_vars.begin(), noisy_vars.end(),
                [&site](const CandidateVariant& other) {
                    return other.key.type == VariantType::Insertion &&
                        other.key.pos == site.key.pos && other.key.ref_len == 0 &&
                        other.key.alt != site.key.alt && other.msa_insertion_alts.empty() &&
                        !other.key.alt.empty() &&
                        (other.key.alt.compare(0, site.key.alt.size(), site.key.alt) == 0 ||
                         site.key.alt.compare(0, other.key.alt.size(), other.key.alt) == 0) &&
                        other.counts.category == VariantCategory::NoisyCandHet &&
                        // One-base error neighborhoods of the two insertion
                        // lengths must not overlap when recalling new reads.
                        std::llabs(static_cast<long long>(other.key.alt.size()) -
                                   static_cast<long long>(site.key.alt.size())) >
                            2 * kMaxLocalAlleleEdits;
                });
        }
        auto recalled_vars = noisy_vars;
        auto recalled_profiles = init_read_profiles(chunk.reads.size());
        add_msa_site_observations(opts, unassigned, noisy_reg_beg,
            recalled_vars, recalled_profiles, n_cons == 2 ? &consensuses : nullptr);
        restrict_msa_observations_to_read_coverage(chunk.reads, recalled_vars, recalled_profiles);
        for (const UnassignedMsaRead& read : unassigned) {
            ReadVariantProfile& profile = recalled_profiles[read.read_id];
            for (size_t offset = 0; offset < profile.alleles.size(); ++offset) {
                if (profile.alleles[offset] < 0) continue;
                const size_t vi = static_cast<size_t>(profile.start_var_idx) + offset;
                if (!eligible[vi]) continue;
                // Two independent non-reference rows describe one diploid
                // contrast. A third allele can be REF to both; that cannot
                // certify either fixed haplotype. Require complementary calls.
                const bool complementary = std::any_of(
                    recalled_vars.begin(), recalled_vars.end(),
                    [&](const CandidateVariant& other) {
                        const size_t oi = static_cast<size_t>(&other - recalled_vars.data());
                        if (oi == vi || !eligible[oi] ||
                            other.key.pos != recalled_vars[vi].key.pos ||
                            static_cast<int>(oi) < profile.start_var_idx) return false;
                        const size_t other_offset = oi - static_cast<size_t>(profile.start_var_idx);
                        return other_offset < profile.alleles.size() &&
                            profile.alleles[other_offset] == 1 - profile.alleles[offset];
                    });
                if (!complementary) continue;
                deferred->push_back({recalled_vars[vi].key, read.read_id,
                                     profile.alleles[offset]});
            }
        }
    } else if (opts.add_unplaced_msa_observations) {
        add_msa_site_observations(opts, unassigned, noisy_reg_beg,
                                  noisy_vars, noisy_rvp, n_cons == 2 ? &consensuses : nullptr);
    }
    // Partial consensus alignments may inherit a cluster allele beyond their
    // physical BAM coverage. Remove those synthetic calls before rephasing.
    if (opts.add_unplaced_msa_observations)
        restrict_msa_observations_to_read_coverage(chunk.reads, noisy_vars, noisy_rvp);
    make_colocated_deletions_exclusive(noisy_vars, noisy_rvp);

    return merge_var_profile(
        chunk, noisy_vars, noisy_var_cate, noisy_rvp, site_whitelist,
        /*admit_all_in_region=*/false);
}

int collect_noisy_vars1(PhasingChunk& chunk, const Options& opts, int noisy_reg_i,
                        const VariantKeySet* site_whitelist) {
    return collect_noisy_vars1(chunk, opts, noisy_reg_i, site_whitelist, nullptr);
}

static bool apply_deferred_msa_observations(PhasingChunk& chunk,
                            std::vector<DeferredMsaObservation>& observations,
                            bool commit_pending = false) {
    std::sort(observations.begin(), observations.end(),
        [](const DeferredMsaObservation& a, const DeferredMsaObservation& b) {
            const int cmp = exact_comp_var_site(&a.key, &b.key);
            return cmp != 0 ? cmp < 0 : a.read_id < b.read_id;
        });
    bool changed = false;
    size_t ci = 0;
    for (size_t oi = 0; oi < observations.size();) {
        const DeferredMsaObservation& observation = observations[oi];
        size_t end = oi + 1;
        int allele = observation.allele;
        bool update_counts = observation.update_counts;
        bool conflicting = false;
        while (end < observations.size() &&
               exact_comp_var_site(&observation.key, &observations[end].key) == 0 &&
               observation.read_id == observations[end].read_id) {
            if (observations[end].allele != allele) conflicting = true;
            update_counts = update_counts || observations[end].update_counts;
            ++end;
        }
        while (ci < chunk.candidates.size() &&
               exact_comp_var_site(&chunk.candidates[ci].key, &observation.key) < 0) ++ci;
        if (!conflicting && allele >= -1 && ci < chunk.candidates.size() &&
            exact_comp_var_site(&chunk.candidates[ci].key, &observation.key) == 0) {
            CandidateVariant& site = chunk.candidates[ci];
            if (site.msa_verified && is_phase_set_anchor(site) &&
                (allele < 0 || static_cast<size_t>(allele) < site.counts.alle_covs.size())) {
                const int hap = chunk.haps[observation.read_id];
                // A fixed source label belongs to its own block gauge. Local
                // recall may extend that membership, not silently move the
                // read to another independently oriented source block.
                if (allele >= 0 && (hap == 1 || hap == 2) && chunk.phase_sets[observation.read_id] > 0 &&
                    (chunk.phase_sets[observation.read_id] != site.phase_set ||
                     site.hap_to_cons_alle[hap] != allele)) {
                    oi = end;
                    continue;
                }
                ReadVariantProfile& profile = chunk.read_var_profile[observation.read_id];
                const int offset = static_cast<int>(ci) - profile.start_var_idx;
                if (commit_pending || profile.start_var_idx < 0 || offset < 0 ||
                    static_cast<size_t>(offset) >= profile.alleles.size() ||
                    profile.alleles[static_cast<size_t>(offset)] < 0) {
                    if (!commit_pending && site.key.alt.find_first_not_of(
                            site.key.alt.front()) != std::string::npos) {
                        // Retry selection must see the original discovery
                        // matrix. Remember that this MSA call was missing now;
                        // later CIGAR REF backfill cannot erase that provenance.
                        chunk.pending_msa_observations.push_back(observation);
                        oi = end;
                        continue;
                    }
                    update_read_var_profile_with_allele(static_cast<int>(ci), allele, -1, profile);
                    if (commit_pending && allele < 0)
                        chunk.rejected_msa_observations.push_back(observation);
                    if (update_counts && allele >= 0) {
                        ++site.counts.alle_covs[static_cast<size_t>(allele)];
                        ++site.counts.total_cov;
                        update_variant_depth_fields(site);
                    }
                    changed = true;
                }
            }
        }
        oi = end;
    }
    return changed;
}

bool apply_pending_msa_observations(PhasingChunk& chunk) {
    std::vector<DeferredMsaObservation> observations;
    observations.swap(chunk.pending_msa_observations);
    const bool changed = apply_deferred_msa_observations(chunk, observations, true);
    if (changed) rebuild_read_variant_index(chunk);
    return changed;
}

// ════════════════════════════════════════════════════════════════════════════
// collect_noisy_vars_step4
// ════════════════════════════════════════════════════════════════════════════

// Attempt each region once per pass; newly admitted sites can help later regions.
static void run_noisy_pass(PhasingChunk& chunk, const Options& opts,
                           const VariantKeySet* site_whitelist,
                           const std::vector<int>& regions,
                           std::vector<bool>& done,
                           std::vector<DeferredMsaObservation>& deferred) {
    while (true) {
        bool any_done = false, any_new_var = false;
        for (int reg_idx : regions) {
            if (done[static_cast<size_t>(reg_idx)]) continue;
            const int ret = collect_noisy_vars1(chunk, opts, reg_idx, site_whitelist, &deferred);
            if (ret >= 0) {
                done[static_cast<size_t>(reg_idx)] = true;
                any_done = true;
                if (ret > 0) any_new_var = true;
            }
        }
        if (any_new_var && !opts.skip_noisy_kmeans) {
            // Recovery needs the repaired calls in the source solve itself.
            // Adding them only after phasing leaves its HP/PS gauge based on
            // the incomplete MSA matrix. Ordinary BAM runs have no retry windows.
            backfill_shifted_msa_insertions(chunk, opts);
            assign_hap_based_on_germline_het_vars_kmeans(chunk, opts, kCandGermlineVarCate,
                                                         opts.anchored_stage2);
        }
        if (!any_done) break;
    }
}

/// Recall candidates from noisy regions with MSA. Regions that cannot be
/// resolved stay pending; a newly called variant triggers a wider k-means pass
/// and may let a pending region resolve on the next iteration.
void collect_noisy_vars_step4(PhasingChunk& chunk, const Options& opts,
                              const VariantKeySet* site_whitelist) {
    if (chunk.noisy_regions.empty()) return;

    const std::vector<int> sorted = sort_noisy_regs(chunk);
    const int n_regs = static_cast<int>(chunk.noisy_regions.size());

    std::vector<bool> done(static_cast<size_t>(n_regs), false);
    std::vector<DeferredMsaObservation> deferred;
    run_noisy_pass(chunk, opts, site_whitelist, sorted, done, deferred);
    if (apply_deferred_msa_observations(chunk, deferred)) {
        // Keep the original BAM phase gauges and genotypes. Transfer and
        // stitching will evaluate the expanded matrix in the graph chunk.
        rebuild_read_variant_index(chunk);
    }
}

} // namespace pgphase_collect
