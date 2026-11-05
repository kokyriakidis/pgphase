/**
 * @file read_profiles.cpp
 * Per-read candidate allele profiling from alignment DIGARs.
 */

#include "collect_var.hpp"
#include "collect_phase.hpp"

#include <cstdlib>
#include <cstring>
#include <utility>
#include <vector>

namespace pgphase_collect {

static inline int read_order_visit_count(const PhasingChunk& chunk) {
    return chunk.ordered_read_ids.empty() ? static_cast<int>(chunk.reads.size())
                                         : static_cast<int>(chunk.ordered_read_ids.size());
}

// Compact variant-site representation for read-vs-candidate profile matching.
struct LcdProfileSite {
    hts_pos_t pos = 0;
    int var_type = 0;
    int ref_len = 0;
    int alt_len = 0;
    const uint8_t* alt_seq = nullptr;
};

static void make_profile_site_from_cand(const CandidateVariant& cv, LcdProfileSite* out) {
    out->pos = cv.key.pos;
    out->alt_seq = reinterpret_cast<const uint8_t*>(cv.key.alt.data());
    if (cv.key.type == VariantType::Snp) {
        out->var_type = BAM_CDIFF;
        out->ref_len = 1;
        out->alt_len = static_cast<int>(cv.key.alt.size());
    } else if (cv.key.type == VariantType::Insertion) {
        out->var_type = BAM_CINS;
        out->ref_len = 0;
        out->alt_len = static_cast<int>(cv.key.alt.size());
    } else {
        out->var_type = BAM_CDEL;
        out->ref_len = cv.key.ref_len;
        out->alt_len = 0;
    }
}

// Build a profile-site from a digar op (including clip rows).
static void make_profile_site_from_digar(const DigarOp& op, LcdProfileSite* out) {
    out->pos = op.pos;
    out->var_type = BAM_CEQUAL;
    out->ref_len = 1;
    out->alt_len = op.len;
    out->alt_seq = nullptr;
    switch (op.type) {
        case DigarType::Snp:
            out->var_type = BAM_CDIFF;
            out->ref_len = 1;
            out->alt_len = static_cast<int>(op.alt.size());
            out->alt_seq = reinterpret_cast<const uint8_t*>(op.alt.data());
            return;
        case DigarType::Insertion:
            out->var_type = BAM_CINS;
            out->ref_len = 0;
            out->alt_len = static_cast<int>(op.alt.size());
            out->alt_seq = reinterpret_cast<const uint8_t*>(op.alt.data());
            return;
        case DigarType::Deletion:
            out->var_type = BAM_CDEL;
            out->ref_len = op.len;
            out->alt_len = 0;
            return;
        case DigarType::SoftClip:
            out->var_type = BAM_CSOFT_CLIP;
            return;
        case DigarType::HardClip:
            out->var_type = BAM_CHARD_CLIP;
            return;
        case DigarType::Equal:
        default:
            out->var_type = BAM_CEQUAL;
            return;
    }
}

static int profile_ovlp_var_site(const LcdProfileSite* var1, const LcdProfileSite* var2) {
    const int beg1 = static_cast<int>(var1->pos);
    const int end1 = static_cast<int>(var1->pos + var1->ref_len);
    const int beg2 = static_cast<int>(var2->pos);
    const int end2 = static_cast<int>(var2->pos + var2->ref_len);
    if (var1->ref_len == 0 && var2->ref_len == 0) return (beg1 == beg2) ? 1 : 0;
    if (var1->ref_len == 0) return (beg1 > beg2 && end1 < end2) ? 1 : 0;
    if (var2->ref_len == 0) return (beg2 > beg1 && end2 < end1) ? 1 : 0;
    if (beg1 >= end2 || beg2 >= end1) return 0;
    return 1;
}

static int profile_exact_comp_var_site(const LcdProfileSite* var1, const LcdProfileSite* var2) {
    hts_pos_t var1_pos = (var1->var_type == BAM_CDIFF) ? var1->pos : var1->pos - 1;
    hts_pos_t var2_pos = (var2->var_type == BAM_CDIFF) ? var2->pos : var2->pos - 1;
    if (var1_pos < var2_pos) return -1;
    if (var1_pos > var2_pos) return 1;
    if (var1->var_type < var2->var_type) return -1;
    if (var1->var_type > var2->var_type) return 1;
    if (var1->ref_len < var2->ref_len) return -1;
    if (var1->ref_len > var2->ref_len) return 1;
    if (var1->alt_len < var2->alt_len) return -1;
    if (var1->alt_len > var2->alt_len) return 1;
    if (var1->var_type == BAM_CDIFF || var1->var_type == BAM_CINS) {
        return std::memcmp(var1->alt_seq, var2->alt_seq, static_cast<size_t>(var1->alt_len));
    }
    return 0;
}

static int profile_comp_ovlp_var_site(const LcdProfileSite* var1, const LcdProfileSite* var2, int* is_ovlp) {
    *is_ovlp = profile_ovlp_var_site(var1, var2);
    return profile_exact_comp_var_site(var1, var2);
}

/**
 * Per-read candidate allele sweep: match each read's digars against sorted candidates.
 *
 * Visits reads in ordered_read_ids order when populated, else 0..n-1.
 * Fills `profiles` (one per read, indices into `variants`) and, when given, `cr` ([start_var_idx, end_var_idx+1) →
 * read_i, not indexed).
 */
static void sweep_read_var_profiles(const Options& opts, const PhasingChunk& chunk, const CandidateTable& variants,
                                    std::vector<ReadVariantProfile>& profiles, cgranges_t* cr) {
    const int n = static_cast<int>(variants.size());
    profiles.assign(chunk.reads.size(), ReadVariantProfile{});
    for (size_t i = 0; i < chunk.reads.size(); ++i)
        profiles[i].read_id = static_cast<int>(i);
    if (n == 0) return;

    for (int ord = 0; ord < read_order_visit_count(chunk); ++ord) {
        const int read_i = chunk.ordered_read_ids.empty()
                               ? ord
                               : chunk.ordered_read_ids[static_cast<size_t>(ord)];
        if (read_i < 0 || static_cast<size_t>(read_i) >= chunk.reads.size()) continue;
        const ReadRecord& read = chunk.reads[static_cast<size_t>(read_i)];
        if (read.is_skipped || read.tid < 0) continue;

        ReadVariantProfile cur;
        cur.read_id = static_cast<int>(read_i);
        const hts_pos_t pos_start = read.beg;
        const hts_pos_t pos_end = read.end;
        int site_i = get_var_site_start(variants, pos_start, n);
        size_t digar_i = 0;
        const std::vector<DigarOp>& dig = read.digars;

        cgranges_t* noisy_cr = intervals_to_cr(read.noisy_regions);
        if (noisy_cr != nullptr) cr_index(noisy_cr);
        // Check if position falls inside a noisy region.
        // Buffer hoisted outside the lambda so cr_overlap reuses it across calls
        // instead of malloc/free per call.
        int64_t* noisy_ovlp_b = nullptr;
        int64_t noisy_max_b = 0;
        const auto is_in_noisy_reg = [&](hts_pos_t pos) -> bool {
            if (noisy_cr == nullptr) return false;
            const int64_t ovlp_n =
                cr_overlap(noisy_cr, "cr", static_cast<int32_t>(pos), static_cast<int32_t>(pos + 1), &noisy_ovlp_b, &noisy_max_b);
            return ovlp_n > 0;
        };

        const auto append_observation = [&](int idx, int allele, int alt_qi) {
            if (cur.start_var_idx < 0) cur.start_var_idx = idx;
            while (cur.start_var_idx + static_cast<int>(cur.alleles.size()) < idx) {
                cur.alleles.push_back(-1);
                cur.alt_qi.push_back(-1);
            }
            cur.end_var_idx = idx;
            cur.alleles.push_back(allele);
            cur.alt_qi.push_back(alt_qi);
        };

        // Skip NON_VAR candidates by bitmask.
        while (site_i < n && digar_i < dig.size()) {
            if (variants[static_cast<size_t>(site_i)].lcd_var_i_to_cate == kLongcalldNonVar) {
                ++site_i;
                continue;
            }
            if (dig[digar_i].type == DigarType::Equal) {
                ++digar_i;
                continue;
            }
            LcdProfileSite var_site0{};
            make_profile_site_from_cand(variants[static_cast<size_t>(site_i)], &var_site0);
            LcdProfileSite digar_site{};
            make_profile_site_from_digar(dig[digar_i], &digar_site);

            const DigarOp& d = dig[digar_i];
            const int ave_qual = get_digar_ave_qual(d, read.qual);
            const int var_read_pos = d.qi;

            int is_ovlp = 0;
            const int ret = profile_comp_ovlp_var_site(&var_site0, &digar_site, &is_ovlp);

            if (is_ovlp == 0) {
                if (ret < 0) {
                    append_observation(site_i, 0, -1);
                    ++site_i;
                } else if (ret > 0) {
                    ++digar_i;
                } else {
                    ++site_i;
                    ++digar_i;
                }
            } else {
                if (ret == 0) {
                    int allele_i = 1;
                    if (ave_qual < opts.min_bq) allele_i = -2;
                    append_observation(site_i, allele_i, var_read_pos);
                    ++site_i;
                } else {
                    append_observation(site_i, -1, -1);
                    ++site_i;
                }
            }
        }

        // Tail loop: mark remaining candidates as non-informative (-1) if in noisy region.
        for (; site_i < n; ++site_i) {
            if (variants[static_cast<size_t>(site_i)].key.pos > pos_end) break;
            if (is_in_noisy_reg(variants[static_cast<size_t>(site_i)].key.pos)) continue;
            append_observation(site_i, 0, -1);
        }

        free(noisy_ovlp_b);
        if (noisy_cr != nullptr) cr_destroy(noisy_cr);

        if (cur.start_var_idx >= 0) {
            if (cr != nullptr)
                cr_add(cr, "cr", cur.start_var_idx, cur.end_var_idx + 1,
                       static_cast<int32_t>(read_i));
            profiles[static_cast<size_t>(read_i)] = std::move(cur);
        }
    }
}

/// `chunk.read_var_profile` and `chunk.read_var_cr` over `chunk.candidates`.
static void collect_read_var_profile(const Options& opts, PhasingChunk& chunk) {
    cgranges_t* cr = cr_init();
    sweep_read_var_profiles(opts, chunk, chunk.candidates, chunk.read_var_profile, cr);
    cr_index(cr);
    chunk.read_var_cr.reset(cr);
}

std::vector<ReadVariantProfile> read_var_profiles_at(const Options& opts, const PhasingChunk& chunk,
                                                     const CandidateTable& variants) {
    std::vector<ReadVariantProfile> profiles;
    sweep_read_var_profiles(opts, chunk, variants, profiles, nullptr);
    return profiles;
}

/// Record each read's allele at every overlapping candidate. Profiles use 0
/// for reference, positive indices for alternate alleles, and negative values
/// when the read is uninformative or low quality. K-means consumes these vectors.
void collect_var_build_profiles(PhasingChunk& chunk, const Options& opts) {
    if (chunk.candidates.empty()) return;
    collect_read_var_profile(opts, chunk);
}

} // namespace pgphase_collect
