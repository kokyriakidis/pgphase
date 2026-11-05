/**
 * @file collect_phase_noisy.cpp
 * @brief Step 4: iterative noisy-region MSA variant calling and profile merge.
 *
 * Step 4 implementation: noisy-region sorting, read collection, reference slice
 * encoding, old/new candidate profile merge, and the outer iterative k-means loop.
 *
 * Alignment-heavy functions live in `align.cpp`; MSA consensus candidate
 * extraction and slice-based allele calling live in `msa_candidates.cpp`.
 */

#include "trace.hpp"
#include "collect_phase_noisy.hpp"
#include "plegma_msa.hpp"
#include "task_pool.hpp"
#include "collect_phase.hpp" // assign_hap_based_on_germline_het_vars_kmeans, kCandGermlineVarCate
#include "collect_var.hpp"   // exact_comp_cand_var, exact_comp_var_site

#include <algorithm>
#include <cassert>
#include <cmath>
#include <numeric>

namespace pgphase_collect {

namespace {

constexpr uint8_t kGapBase = 5;

void set_noisy_category(CandidateVariant& var, VariantCategory category) {
    var.counts.category = category;
    var.lcd_var_i_to_cate = category_to_flag(category);
}

std::vector<ReadVariantProfile> init_read_profiles(size_t n_reads) {
    std::vector<ReadVariantProfile> profiles(n_reads);
    for (size_t i = 0; i < n_reads; ++i) profiles[i].read_id = static_cast<int>(i);
    return profiles;
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
// merge_var_profile
// ════════════════════════════════════════════════════════════════════════════

int merge_var_profile(PhasingChunk& chunk,
                      const std::vector<CandidateVariant>& noisy_vars,
                      const std::vector<VariantCategory>& noisy_var_cate,
                      const std::vector<ReadVariantProfile>& noisy_rvp,
                      const VariantKeySet* site_whitelist,
                      bool admit_all_in_region,
                      const VariantKeySet* replace_sites,
                      bool trace_admitted) {
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

    // The chunk's tables are rebuilt below: take them rather than copy them.
    CandidateTable old_vars = std::move(chunk.candidates);
    std::vector<ReadVariantProfile> old_profiles = chunk.read_var_profile.empty()
                                                       ? init_read_profiles(chunk.reads.size())
                                                       : std::move(chunk.read_var_profile);

    CandidateTable merged_vars;
    merged_vars.reserve(old_vars.size() + new_vars.size());
    std::vector<int> old_to_merged(old_vars.size(), -1);
    std::vector<int> new_to_merged(new_vars.size(), -1);
    // Merged row carrying each old row's key: the row itself, or the new row
    // that replaced it at the same key.
    std::vector<int> old_key_to_merged(old_vars.size(), -1);

    size_t old_i = 0;
    size_t new_i = 0;
    int admitted = 0;
    while (old_i < old_vars.size() && new_i < new_vars.size()) {
        const int ret = exact_comp_var_site(&old_vars[old_i].key, &new_vars[new_i].key);
        if (ret < 0) {
            old_to_merged[old_i] = old_key_to_merged[old_i] = static_cast<int>(merged_vars.size());
            merged_vars.push_back(std::move(old_vars[old_i++]));
        } else if (ret > 0) {
            if ((site_whitelist != nullptr &&
                 site_whitelist->find(new_vars[new_i].key) == site_whitelist->end())) {
                pgphase_trace::candidate("msa", "drop", new_vars[new_i], "not_whitelisted", "");
                ++new_i;
                continue;
            }
            set_noisy_category(new_vars[new_i], new_cats[new_i]);
            if (trace_admitted)
                pgphase_trace::candidate("msa", "call", new_vars[new_i], "admitted",
                                         "cat=" + std::to_string(static_cast<int>(new_cats[new_i])));
            new_to_merged[order[new_i]] = static_cast<int>(merged_vars.size());
            merged_vars.push_back(std::move(new_vars[new_i++]));
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
                pgphase_trace::candidate("msa", "call", new_vars[new_i], "replaced_existing",
                                         "cat=" + std::to_string(static_cast<int>(new_cats[new_i])) + " old_cat=" +
                                         std::to_string(static_cast<int>(old_vars[old_i].counts.category)));
                new_vars[new_i].graph_site |= old_vars[old_i].graph_site;
                new_vars[new_i].counts.candvarcate_initial =
                    old_vars[old_i].counts.candvarcate_initial;
                new_to_merged[order[new_i]] = old_key_to_merged[old_i] = static_cast<int>(merged_vars.size());
                merged_vars.push_back(std::move(new_vars[new_i]));
                ++admitted;
            } else {
                pgphase_trace::candidate("msa", "drop", new_vars[new_i], "existing_candidate_kept",
                                         "cat=" + std::to_string(static_cast<int>(new_cats[new_i])) + " old_cat=" +
                                         std::to_string(static_cast<int>(old_vars[old_i].counts.category)));
                old_to_merged[old_i] = old_key_to_merged[old_i] = static_cast<int>(merged_vars.size());
                merged_vars.push_back(std::move(old_vars[old_i]));
            }
            ++old_i;
            ++new_i;
        }
    }
    while (old_i < old_vars.size()) {
        old_to_merged[old_i] = old_key_to_merged[old_i] = static_cast<int>(merged_vars.size());
        merged_vars.push_back(std::move(old_vars[old_i++]));
    }
    while (new_i < new_vars.size()) {
        if ((site_whitelist != nullptr &&
             site_whitelist->find(new_vars[new_i].key) == site_whitelist->end())) {
            pgphase_trace::candidate("msa", "drop", new_vars[new_i], "not_whitelisted", "");
            ++new_i;
            continue;
        }
        set_noisy_category(new_vars[new_i], new_cats[new_i]);
        new_to_merged[order[new_i]] = static_cast<int>(merged_vars.size());
        merged_vars.push_back(std::move(new_vars[new_i++]));
        ++admitted;
    }

    std::vector<ReadVariantProfile> merged_profiles = init_read_profiles(chunk.reads.size());
    const std::vector<ReadVariantProfile> empty_new_profiles =
        noisy_rvp.empty() ? init_read_profiles(chunk.reads.size()) : std::vector<ReadVariantProfile>{};
    const std::vector<ReadVariantProfile>& new_profiles = noisy_rvp.empty() ? empty_new_profiles : noisy_rvp;

    const int n_reads = static_cast<int>(chunk.reads.size());
    // A read with no new observation whose old rows all keep their order with
    // nothing inserted or dropped among them (most reads: the merge adds one
    // noisy region's calls) keeps its profile, shifted. Rebuilt entry by entry
    // it comes out the same: alleles and alt_qi as before; graph_alleles,
    // bam_alleles and bam_qi padded (or cut) to the span with -1; base
    // qualities and MAPQ not carried (as the full rebuild below).
    std::vector<char> shifted(static_cast<size_t>(n_reads), 0);
    for (int ord = 0; ord < (chunk.ordered_read_ids.empty() ? n_reads : static_cast<int>(chunk.ordered_read_ids.size())); ++ord) {
        const int read_i = chunk.ordered_read_ids.empty()
                               ? ord
                               : chunk.ordered_read_ids[static_cast<size_t>(ord)];
        if (read_i < 0 || read_i >= n_reads) continue;
        if (chunk.reads[static_cast<size_t>(read_i)].is_skipped) continue;
        if (static_cast<size_t>(read_i) >= old_profiles.size() || shifted[static_cast<size_t>(read_i)]) continue;
        ReadVariantProfile& old_p = old_profiles[static_cast<size_t>(read_i)];
        const ReadVariantProfile* new_p = static_cast<size_t>(read_i) < new_profiles.size()
                                              ? &new_profiles[static_cast<size_t>(read_i)]
                                              : nullptr;
        if (new_p != nullptr && new_p->start_var_idx >= 0 && new_p->end_var_idx >= new_p->start_var_idx) continue;
        const int s0 = old_p.start_var_idx, e0 = old_p.end_var_idx;
        if (s0 < 0 || e0 < s0 || static_cast<size_t>(e0) >= old_to_merged.size() ||
            old_p.alleles.size() != static_cast<size_t>(e0 - s0 + 1) || old_p.alt_qi.size() != old_p.alleles.size())
            continue;
        const int m0 = old_to_merged[static_cast<size_t>(s0)];
        bool contiguous = m0 >= 0;
        for (int v = s0; contiguous && v <= e0; ++v)
            contiguous = old_to_merged[static_cast<size_t>(v)] == m0 + (v - s0) &&
                         old_key_to_merged[static_cast<size_t>(v)] == m0 + (v - s0);
        if (!contiguous) continue;
        const size_t n = old_p.alleles.size();
        ReadVariantProfile& dest = merged_profiles[static_cast<size_t>(read_i)];
        dest.start_var_idx = m0;
        dest.end_var_idx = m0 + (e0 - s0);
        dest.alleles = std::move(old_p.alleles);
        dest.alt_qi = std::move(old_p.alt_qi);
        if (!old_p.graph_alleles.empty()) {
            dest.graph_alleles = std::move(old_p.graph_alleles);
            dest.graph_alleles.resize(n, -1);
        }
        dest.bam_alleles = std::move(old_p.bam_alleles);
        dest.bam_alleles.resize(n, -1);
        dest.bam_qi = std::move(old_p.bam_qi);
        dest.bam_qi.resize(n, -1);
        shifted[static_cast<size_t>(read_i)] = 1;
    }
    for (int ord = 0; ord < (chunk.ordered_read_ids.empty() ? n_reads : static_cast<int>(chunk.ordered_read_ids.size())); ++ord) {
        const int read_i = chunk.ordered_read_ids.empty()
                               ? ord
                               : chunk.ordered_read_ids[static_cast<size_t>(ord)];
        if (read_i < 0 || read_i >= n_reads) continue;
        if (chunk.reads[static_cast<size_t>(read_i)].is_skipped) continue;
        if (shifted[static_cast<size_t>(read_i)]) continue;
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
        if (ri < shifted.size() && shifted[ri]) continue;
        auto& dest = merged_profiles[ri];
        const auto& source = old_profiles[ri];
        dest.bam_alleles.assign(dest.alleles.size(), -1);
        dest.bam_qi.assign(dest.alleles.size(), -1);
        for (size_t pi = 0; pi < source.alleles.size(); ++pi) {
            const int vi = source.start_var_idx + static_cast<int>(pi);
            if (vi < 0 || static_cast<size_t>(vi) >= old_vars.size()) continue;
            const int mi = old_key_to_merged[static_cast<size_t>(vi)];
            if (mi < 0) continue;
            const int di = mi - dest.start_var_idx;
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

// query the consensus).
static void record_msa_region(PhasingChunk& chunk, hts_pos_t beg, hts_pos_t end, int n_cons,
                              const std::array<std::vector<AlnStr>, 2>& aln_strs, const std::array<int, 2>& clu_n_seqs,
                              const std::array<std::vector<int>, 2>& clu_read_ids) {
    PhasingChunk::MsaRegion region{beg, end, {}};
    for (int c = 0; c < n_cons && c < 2; ++c) {
        if (aln_strs[static_cast<size_t>(c)].empty()) return;
        const AlnStr& aln = aln_strs[static_cast<size_t>(c)][0];
        std::vector<std::string> spelled(static_cast<size_t>(end - beg + 1));
        const auto letter = [](uint8_t q) { return q < 4 ? "ACGT"[q] : 'N'; };
        hts_pos_t ref_pos = beg;
        for (int i = 0; i < aln.aln_len; ++i) {
            const uint8_t r = aln.target_aln[static_cast<size_t>(i)];
            const uint8_t q = aln.query_aln[static_cast<size_t>(i)];
            if (r == kGapBase && q == kGapBase) continue;
            if (r == kGapBase) {
                region.differs.push_back(Interval{ref_pos - 1, ref_pos, 0});
                if (ref_pos > beg && ref_pos - 1 - beg < static_cast<hts_pos_t>(spelled.size()))
                    spelled[static_cast<size_t>(ref_pos - 1 - beg)] += letter(q);
                continue;
            }
            if (q != r) region.differs.push_back(Interval{ref_pos, ref_pos, 0});
            if (q != kGapBase && ref_pos - beg < static_cast<hts_pos_t>(spelled.size()))
                spelled[static_cast<size_t>(ref_pos - beg)] += letter(q);
            ++ref_pos;
        }
        region.spelled.push_back(std::move(spelled));
        const auto& ids = clu_read_ids[static_cast<size_t>(c)];
        region.cluster_reads[static_cast<size_t>(c)].assign(
            ids.begin(), ids.begin() + std::min<long>(static_cast<long>(ids.size()), clu_n_seqs[static_cast<size_t>(c)]));
    }
    chunk.msa_regions.push_back(std::move(region));
}

NoisyRegionCall compute_noisy_region(const PhasingChunk& chunk, const Options& opts, int noisy_reg_i,
                                     const VariantKeySet* site_whitelist) {
    NoisyRegionCall call;
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
        if (!contains_whitelisted_site) return call;
    }

    // Skip regions longer than max_noisy_reg_len (return 0 = done, no vars).
    if (noisy_reg_end - noisy_reg_beg + 1 > static_cast<hts_pos_t>(opts.max_noisy_reg_len)) {
        pgphase_trace::emit("msa", "region", noisy_reg_beg, noisy_reg_end, "-", "too_long_for_msa");
        return call;
    }

    // First clip beg/end via collect_reg_ref_bseq, then collect reads on the
    // finalized region bounds.
    std::vector<uint8_t> ref_seq =
        collect_reg_ref_bseq(chunk, noisy_reg_beg, noisy_reg_end);
    if (ref_seq.empty()) {
        pgphase_trace::emit("msa", "region", noisy_reg_beg, noisy_reg_end, "-", "no_reference");
        return call;
    }

    std::vector<int> read_ids =
        collect_noisy_reg_reads(chunk, noisy_reg_beg, noisy_reg_end);
    // A region mostly of reads placed at low MAPQ (a collapsed repeat: a satellite array, a segmental duplication) holds
    // several copies' reads; their MSA clusters by copy, not haplotype. Its MSA takes the confidently placed reads
    // (kDefaultMinMapq) alone, those of this copy; the others still give their calls at the sites it makes. With fewer
    // than kMsaMinConfidentReads of them, a region takes all its reads, as before (skipped, chr20's centromere, 26-30
    // Mb, lost 750 truth hets of 2,550, with 2,198 false het loci and 1,339 wrong alleles fewer), and a piece of a split
    // region none (skipped whole before the split; all its reads took the centromere's four chunks from 41 s to over
    // 14 minutes). plegma only (abPOA keeps longcallD's reads).
    if (use_plegma_msa() && !read_ids.empty()) {
        std::vector<int> confident;
        for (int ri : read_ids)
            if (chunk.reads[static_cast<size_t>(ri)].mapq >= kDefaultMinMapq) confident.push_back(ri);
        const bool collapsed = confident.size() < kMsaMinConfidentShare * read_ids.size();
        if (collapsed && reg.piece && static_cast<int>(confident.size()) < kMsaMinConfidentReads) {
            pgphase_trace::emit("msa", "region", noisy_reg_beg, noisy_reg_end, "-", "few_confident_reads",
                                "reads=" + std::to_string(read_ids.size()) + " confident=" +
                                    std::to_string(confident.size()));
            return call;
        }
        if (collapsed && static_cast<int>(confident.size()) >= kMsaMinConfidentReads) {
            pgphase_trace::emit("msa", "region", noisy_reg_beg, noisy_reg_end, "-", "confident_reads_only",
                                "reads=" + std::to_string(read_ids.size()) + " confident=" +
                                    std::to_string(confident.size()));
            read_ids = std::move(confident);
        }
    }
    const int n_noisy_reads = static_cast<int>(read_ids.size());

    // Skip regions with more than max_noisy_reg_cov reads.
    if (n_noisy_reads > opts.max_noisy_reg_cov) {
        pgphase_trace::emit("msa", "region", noisy_reg_beg, noisy_reg_end, "-", "too_many_reads",
                            "reads=" + std::to_string(n_noisy_reads));
        return call;
    }

    // MSA + consensus via WFA2-lib / abPOA (align.cpp).
    std::array<int, 2>& clu_n_seqs = call.clu_n_seqs;
    std::array<std::vector<int>, 2>& clu_read_ids = call.clu_read_ids;
    std::array<std::vector<AlnStr>, 2>& aln_strs = call.aln_strs;
    std::vector<UnassignedMsaRead> unassigned;
    const bool recall_unplaced = opts.add_unplaced_msa_observations;
    const int n_cons = collect_noisy_reg_aln_strs(
        opts, chunk, noisy_reg_beg, noisy_reg_end,
        read_ids, ref_seq,
        clu_n_seqs, clu_read_ids, aln_strs,
        // A broad retry may extend clusters. Local insertion recall instead
        // composes both fixed paths and defers read-local allele observations.
        recall_unplaced ? &unassigned : nullptr, &call.spans);

    // n_cons == 0 → MSA could not resolve; return -1 so the
    // outer loop leaves this region undone and may retry if another region makes progress.
    const std::string why_none = take_msa_failure();
    if (n_cons == 0) {
        pgphase_trace::emit("msa", "region", noisy_reg_beg, noisy_reg_end, "-", "no_consensus",
                            "reads=" + std::to_string(n_noisy_reads) + (why_none.empty() ? "" : " why=" + why_none));
        call.status = -1;
        return call;
    }
    if (pgphase_trace::on())
        pgphase_trace::emit("msa", "region", noisy_reg_beg, noisy_reg_end, "-", "resolved",
                            "consensus=" + std::to_string(n_cons) + " clusters=" + std::to_string(clu_n_seqs[0]) + "/" +
                            std::to_string(clu_n_seqs[1]) + " reads=" + std::to_string(n_noisy_reads));

    call.status = 1;
    call.beg = noisy_reg_beg;
    call.end = noisy_reg_end;
    call.n_cons = n_cons;
    std::vector<CandidateVariant>&   noisy_vars = call.vars;
    std::vector<VariantCategory>&    noisy_var_cate = call.var_cate;
    std::vector<ReadVariantProfile>& noisy_rvp = call.rvp;

    make_vars_from_msa_cons_aln(
        opts, chunk,
        n_noisy_reads, read_ids, noisy_reg_beg,
        n_cons, clu_n_seqs, clu_read_ids, aln_strs,
        noisy_vars, noisy_var_cate, noisy_rvp);
    for (CandidateVariant& v : noisy_vars) {
        v.msa_region_reads = n_noisy_reads;
        v.msa_cluster_reads = clu_n_seqs[0] + (n_cons > 1 ? clu_n_seqs[1] : 0);
    }

    // A real neighboring variant is valid flank sequence, not alignment noise.
    std::array<AlnStr, 2> consensuses;
    if (n_cons == 2) consensuses = {aln_strs[0][0], aln_strs[1][0]};
    if (opts.add_unplaced_msa_observations) {
        add_msa_site_observations(opts, unassigned, noisy_reg_beg,
                                  noisy_vars, noisy_rvp, n_cons == 2 ? &consensuses : nullptr);
    }
    // Partial consensus alignments may inherit a cluster allele beyond their
    // physical BAM coverage. Remove those synthetic calls before rephasing.
    if (opts.add_unplaced_msa_observations)
        restrict_msa_observations_to_read_coverage(chunk.reads, noisy_vars, noisy_rvp);
    make_colocated_deletions_exclusive(noisy_vars, noisy_rvp);
    return call;
}

int apply_noisy_region(PhasingChunk& chunk, const Options& opts, NoisyRegionCall& call,
                       const VariantKeySet* site_whitelist) {
    if (call.status != 1) return call.status;
    if (opts.refine_msa_digars && !call.spans.full_covers.empty())
        update_digars_from_msa(chunk, call.beg, call.end, call.spans, call.aln_strs, call.n_cons, call.clu_n_seqs,
                               call.clu_read_ids);
    if (use_plegma_msa())
        record_msa_region(chunk, call.beg, call.end, call.n_cons, call.aln_strs, call.clu_n_seqs, call.clu_read_ids);
    return merge_var_profile(
        chunk, call.vars, call.var_cate, call.rvp, site_whitelist,
        /*admit_all_in_region=*/false);
}

int collect_noisy_vars1(PhasingChunk& chunk, const Options& opts, int noisy_reg_i,
                        const VariantKeySet* site_whitelist) {
    NoisyRegionCall call = compute_noisy_region(chunk, opts, noisy_reg_i, site_whitelist);
    return apply_noisy_region(chunk, opts, call, site_whitelist);
}

// ════════════════════════════════════════════════════════════════════════════
// collect_noisy_vars_step4
// ════════════════════════════════════════════════════════════════════════════

// A clean hom SNP held back for its noisy region's MSA to recall (PhasingChunk::held_noisy_homs) comes back as the
// pileup called it when no MSA resolved the region -- no consensus (too few reads span it: chr6:32,495,446-32,505,330,
// HLA-DRB, 0 of 13 reads full length), a split piece of few confident reads, too many reads, too long: nothing else
// calls the locus then. chr6:32,502,269 A>G, 65 of 65 reads ALT at MAPQ 20 or more, the truth's hom: the pileup's call
// pruned for the MSA, the MSA failed and the graph dropped its pair (pair_high_af), and nothing was written. Back as a
// CleanHom candidate (restored_noisy_hom), with its reads' alleles, it is written as the alignment's other homs are, if
// both haplotypes of the alignment solve carry it and the graph's reads call a catalog pair with its letter homozygous
// (select_sites: an SV breakpoint or a deletion leaves one haplotype's reads alone there -- chr14's IGH, at half
// depth, the truth .|1 -- and alignments through a repeat spell letters no haplotype has). Only regions no MSA
// resolved, and no key a candidate already has.
static void restore_failed_noisy_homs(PhasingChunk& chunk, const Options& opts, const std::vector<char>& resolved) {
    CandidateTable held = std::move(chunk.held_noisy_homs);
    chunk.held_noisy_homs.clear();
    if (held.empty()) return;
    std::vector<std::pair<hts_pos_t, hts_pos_t>> spans[2];  // regions not resolved, resolved
    for (size_t i = 0; i < chunk.noisy_regions.size() && i < resolved.size(); ++i)
        spans[resolved[i] ? 1 : 0].emplace_back(chunk.noisy_regions[i].beg, chunk.noisy_regions[i].end);
    if (spans[0].empty()) return;
    const auto covers = [](const std::vector<std::pair<hts_pos_t, hts_pos_t>>& s, hts_pos_t p) {
        for (const auto& [b, e] : s)
            if (b <= p && p <= e) return true;
        return false;
    };
    CandidateTable back;
    for (CandidateVariant& cv : held) {
        if (!covers(spans[0], cv.key.pos) || covers(spans[1], cv.key.pos)) continue;
        const auto it = std::lower_bound(chunk.candidates.begin(), chunk.candidates.end(), cv,
                                         [](const CandidateVariant& a, const CandidateVariant& b) {
                                             return exact_comp_var_site(&a.key, &b.key) < 0;
                                         });
        if (it != chunk.candidates.end() && exact_comp_var_site(&it->key, &cv.key) == 0) continue;
        cv.counts.category = VariantCategory::CleanHom;
        cv.lcd_var_i_to_cate = category_to_flag(VariantCategory::CleanHom);
        cv.restored_noisy_hom = true;
        pgphase_trace::candidate("bam", "restore", cv, "failed_msa_region",
                                 "alt=" + std::to_string(cv.counts.alt_cov) + "/" +
                                     std::to_string(cv.counts.total_cov));
        back.push_back(std::move(cv));
    }
    if (back.empty()) return;
    std::sort(back.begin(), back.end(), [](const CandidateVariant& a, const CandidateVariant& b) {
        return exact_comp_cand_var(&a, &b) < 0;
    });
    const std::vector<ReadVariantProfile> profiles = read_var_profiles_at(opts, chunk, back);
    const std::vector<VariantCategory> categories(back.size(), VariantCategory::CleanHom);
    merge_var_profile(chunk, back, categories, profiles, nullptr, false, nullptr, /*trace_admitted=*/false);
}

// Attempt each region once per pass; newly admitted sites can help later regions.
static void run_noisy_pass(PhasingChunk& chunk, const Options& opts,
                           const VariantKeySet* site_whitelist,
                           const std::vector<int>& regions,
                           std::vector<bool>& done,
                           std::vector<char>& resolved) {
    while (true) {
        bool any_done = false, any_new_var = false;
        const auto settle = [&](int reg_idx, int ret) {
            if (ret >= 0) {
                done[static_cast<size_t>(reg_idx)] = true;
                any_done = true;
                if (ret > 0) any_new_var = true;
            }
        };
        if (opts.tasks != nullptr && !opts.refine_msa_digars) {
            // The pass's regions computed at once on the batch's workers, then applied in order: compute_noisy_region
            // reads nothing apply_noisy_region writes, so the chunk comes out as one at a time. Not with
            // refine_msa_digars: a region's digar rewrite changes the reads the next one aligns.
            std::vector<int> pending;
            for (int reg_idx : regions)
                if (!done[static_cast<size_t>(reg_idx)]) pending.push_back(reg_idx);
            std::vector<NoisyRegionCall> calls(pending.size());
            const pgphase_trace::Context trace_context = pgphase_trace::context();
            opts.tasks->run_all(pending.size(), [&](size_t k) {
                const pgphase_trace::Scope scope(trace_context);
                calls[k] = compute_noisy_region(chunk, opts, pending[k], site_whitelist);
                resolved[static_cast<size_t>(pending[k])] |= calls[k].status == 1;
                // Held until applied: of the reads' rows only the digar rewrite reads any (record_msa_region takes
                // each consensus's own row).
                for (std::vector<AlnStr>& rows : calls[k].aln_strs)
                    if (rows.size() > 1) rows.resize(1);
            });
            for (size_t k = 0; k < pending.size(); ++k)
                settle(pending[k], apply_noisy_region(chunk, opts, calls[k], site_whitelist));
        } else {
            for (int reg_idx : regions) {
                if (done[static_cast<size_t>(reg_idx)]) continue;
                NoisyRegionCall call = compute_noisy_region(chunk, opts, reg_idx, site_whitelist);
                resolved[static_cast<size_t>(reg_idx)] |= call.status == 1;
                settle(reg_idx, apply_noisy_region(chunk, opts, call, site_whitelist));
            }
        }
        if (any_new_var) {
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
    std::vector<char> resolved(static_cast<size_t>(n_regs), 0);
    run_noisy_pass(chunk, opts, site_whitelist, sorted, done, resolved);
    if (opts.restore_failed_noisy_homs) restore_failed_noisy_homs(chunk, opts, resolved);
}

} // namespace pgphase_collect
