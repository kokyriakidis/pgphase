/**
 * @file test_phase_block_stitch.cpp
 * Unit tests for chunk stitching (stitch_chunk_haps / flip logic).
 *
 * Exercises boundary stitching, within-chunk allele links, and MSA profile admission.
 *
 * Build (from repo root):
 *   make collect_phase.o && g++ -O0 -g -std=c++17 -Wall -Wextra \
 *       src/test_phase_block_stitch.cpp collect_phase.o -lhts -lz -lpthread \
 *       -o test_phase_block_stitch && ./test_phase_block_stitch
 */

#include "collect_phase.hpp"
#include "collect_output.hpp"
#include "collect_phase_noisy.hpp"
#include "collect_types.hpp"
#include "gap_recovery.hpp"
#include "gap_evidence.hpp"

#include <cstdio>
#include <numeric>
#include <vector>
#include <sstream>
#include <unistd.h>

using namespace pgphase_collect;

static bool check(bool cond, const char* msg) {
    if (!cond) {
        std::printf("FAIL: %s\n", msg);
        return false;
    }
    return true;
}

static ReadRecord min_read() {
    ReadRecord r;
    r.is_skipped = false;
    return r;
}

static CandidateVariant dummy_cand(hts_pos_t ps) {
    CandidateVariant v;
    v.phase_set = ps;
    v.hap_to_cons_alle = {-1, 1, 0};
    return v;
}

static bool test_gap_edges_compose_before_relabelling() {
    PhasingChunk chunk;
    for (const hts_pos_t ps : {100, 200, 300}) {
        chunk.reads.push_back(min_read());
        chunk.haps.push_back(ps == 300 ? 2 : 1);
        chunk.phase_sets.push_back(ps);
        chunk.candidates.push_back(dummy_cand(ps));
    }
    std::vector<PhasingChunk> chunks;
    chunks.push_back(std::move(chunk));
    const std::vector<GapPhaseEdge> edges = {
        {100, 200, true},
        {200, 300, true},
        {100, 300, true},
    };
    const int conflicts = apply_gap_phase_edges(chunks, edges);
    bool ok = true;
    ok &= check(conflicts == 1, "inconsistent composed gap edge is rejected");
    ok &= check(chunks[0].phase_sets == std::vector<hts_pos_t>({100, 100, 100}),
                "gap edge chain receives one canonical phase set");
    ok &= check(chunks[0].haps == std::vector<int>({1, 2, 2}),
                "gap edge orientations compose before read relabelling");
    ok &= check(chunks[0].candidates[1].hap_to_cons_alle[1] == 0 &&
                    chunks[0].candidates[2].hap_to_cons_alle[1] == 1,
                "candidate orientations follow the composed parity");
    return ok;
}

// Two chunks, one overlapping pair: pre down index 0, cur up index 0.
static bool test_cross_chunk_flip_when_haps_disagree() {
    std::printf("--- test_cross_chunk_flip_when_haps_disagree ---\n");
    std::vector<PhasingChunk> chunks;
    chunks.resize(2);
    PhasingChunk& pre = chunks[0];
    PhasingChunk& cur = chunks[1];
    pre.region.tid = 0;
    cur.region.tid = 0;

    pre.reads.push_back(min_read());
    pre.haps.push_back(1);
    pre.phase_sets.push_back(500);
    pre.candidates.push_back(dummy_cand(500));
    pre.down_ovlp_read_i = {{0}};
    pre.up_ovlp_read_i = {{}};

    cur.reads.push_back(min_read());
    cur.haps.push_back(2);
    cur.phase_sets.push_back(900);
    cur.candidates.push_back(dummy_cand(900));
    cur.up_ovlp_read_i = {{0}};
    cur.down_ovlp_read_i = {{}};

    Options opts;
    stitch_chunk_haps(chunks, &opts, nullptr);

    bool ok = true;
    const CandidateVariant& v = chunks[1].candidates[0];
    ok &= check(v.phase_set == 500, "phase_set merged to pre anchor");
    ok &= check(v.hap_to_cons_alle[1] == 0 && v.hap_to_cons_alle[2] == 1, "hap consensus swapped");
    std::printf("%s\n", ok ? "PASS" : "FAIL");
    return ok;
}

// Equal votes: one agree (−1) and one disagree (+1) ⇒ flip_hap_score == 0 ⇒ no updates.
static bool test_flip_score_zero_no_merge() {
    std::printf("--- test_flip_score_zero_no_merge ---\n");
    PhasingChunk pre, cur;
    pre.region.tid = cur.region.tid = 0;

    for (int i = 0; i < 2; ++i) {
        pre.reads.push_back(min_read());
        pre.haps.push_back(1);
        pre.phase_sets.push_back(100);
    }
    pre.candidates.push_back(dummy_cand(100));
    pre.down_ovlp_read_i = {{0, 1}};

    for (int i = 0; i < 2; ++i) {
        cur.reads.push_back(min_read());
        cur.haps.push_back(0);
        cur.phase_sets.push_back(200);
    }
    cur.haps[0] = 1;  // agree with pre
    cur.haps[1] = 2;  // disagree
    cur.candidates.push_back(dummy_cand(200));
    cur.up_ovlp_read_i = {{0, 1}};

    std::vector<PhasingChunk> chunks;
    chunks.push_back(std::move(pre));
    chunks.push_back(std::move(cur));

    Options opts;
    stitch_chunk_haps(chunks, &opts, nullptr);

    bool ok = true;
    ok &= check(chunks[1].candidates[0].phase_set == 200, "phase_set unchanged when flip score ties");
    ok &= check(chunks[1].candidates[0].hap_to_cons_alle[1] == 1 &&
                    chunks[1].candidates[0].hap_to_cons_alle[2] == 0,
                "hap consensus unchanged when flip score ties");
    std::printf("%s\n", ok ? "PASS" : "FAIL");
    return ok;
}

static bool test_skipped_overlap_read_ignored() {
    std::printf("--- test_skipped_overlap_read_ignored ---\n");
    PhasingChunk pre, cur;
    pre.region.tid = cur.region.tid = 0;

    ReadRecord rskip = min_read();
    rskip.is_skipped = true;
    pre.reads.push_back(std::move(rskip));
    pre.haps.push_back(1);
    pre.phase_sets.push_back(10);
    pre.candidates.push_back(dummy_cand(10));
    pre.down_ovlp_read_i = {{0}};

    cur.reads.push_back(min_read());
    cur.haps.push_back(2);
    cur.phase_sets.push_back(20);
    cur.candidates.push_back(dummy_cand(20));
    cur.up_ovlp_read_i = {{0}};

    std::vector<PhasingChunk> chunks;
    chunks.push_back(std::move(pre));
    chunks.push_back(std::move(cur));

    Options opts;
    stitch_chunk_haps(chunks, &opts, nullptr);

    bool ok = check(chunks[1].candidates[0].phase_set == 20, "no informative overlap reads → no merge");
    std::printf("%s\n", ok ? "PASS" : "FAIL");
    return ok;
}

// Build two chunks joined by a single agreeing overlap read, so
// flip_hap_score == +1.  At margin 0 they merge; at margin 1 the boundary
// abstains (|1| <= 1) and the blocks stay separate.
static void make_single_agree_pair(std::vector<PhasingChunk>& chunks) {
    chunks.clear();
    chunks.resize(2);
    PhasingChunk& pre = chunks[0];
    PhasingChunk& cur = chunks[1];
    pre.region.tid = cur.region.tid = 0;

    pre.reads.push_back(min_read());
    pre.haps.push_back(1);
    pre.phase_sets.push_back(100);
    pre.candidates.push_back(dummy_cand(100));
    pre.down_ovlp_read_i = {{0}};

    cur.reads.push_back(min_read());
    cur.haps.push_back(1);  // agrees with pre → score +1 (no flip)
    cur.phase_sets.push_back(200);
    cur.candidates.push_back(dummy_cand(200));
    cur.up_ovlp_read_i = {{0}};
}

static bool test_margin_zero_merges_single_vote() {
    std::printf("--- test_margin_zero_merges_single_vote ---\n");
    std::vector<PhasingChunk> chunks;
    make_single_agree_pair(chunks);

    Options opts;
    opts.stitch_min_margin = 0;  // default behavior
    stitch_chunk_haps(chunks, &opts, nullptr);

    bool ok = check(chunks[1].candidates[0].phase_set == 100,
                    "single agreeing vote merges at margin 0");
    std::printf("%s\n", ok ? "PASS" : "FAIL");
    return ok;
}

static bool test_margin_one_abstains_single_vote() {
    std::printf("--- test_margin_one_abstains_single_vote ---\n");
    std::vector<PhasingChunk> chunks;
    make_single_agree_pair(chunks);

    Options opts;
    opts.stitch_min_margin = 1;  // |score|=1 <= 1 → abstain
    stitch_chunk_haps(chunks, &opts, nullptr);

    bool ok = check(chunks[1].candidates[0].phase_set == 200,
                    "single agreeing vote abstains at margin 1");
    std::printf("%s\n", ok ? "PASS" : "FAIL");
    return ok;
}

static CandidateVariant noisy_merge_cand(hts_pos_t pos, VariantCategory category,
                                         int alt_cov) {
    CandidateVariant cand;
    cand.key.tid = 0;
    cand.key.pos = pos;
    cand.key.type = VariantType::Deletion;
    cand.key.ref_len = 1;
    cand.counts.category = category;
    cand.counts.candvarcate_initial = category;
    cand.counts.alt_cov = alt_cov;
    cand.lcd_var_i_to_cate = category_to_flag(category);
    return cand;
}

static bool test_private_msa_merge_admits_only_whitelist() {
    std::printf("--- test_private_msa_merge_admits_only_whitelist ---\n");
    PhasingChunk chunk;
    chunk.region.tid = 0;
    chunk.reads.push_back(min_read());
    chunk.candidates.push_back(
        noisy_merge_cand(100, VariantCategory::RepeatHetIndel, 3));
    chunk.candidates[0].graph_site = true;
    ReadVariantProfile original;
    original.start_var_idx = original.end_var_idx = 0;
    original.alleles = {0};
    original.alt_qi = {kGraphConfirmedAltQi};
    original.bam_alleles = {-2};
    original.bam_qi = {12};
    original.graph_alleles = {0};
    chunk.read_var_profile.push_back(original);

    std::vector<CandidateVariant> msa_vars = {
        noisy_merge_cand(100, VariantCategory::NoisyCandHet, 8),
        noisy_merge_cand(200, VariantCategory::NoisyCandHet, 9),
        noisy_merge_cand(300, VariantCategory::NoisyCandHet, 10),
    };
    std::vector<VariantCategory> msa_cats(
        msa_vars.size(), VariantCategory::NoisyCandHet);
    std::vector<ReadVariantProfile> msa_profiles(1);
    msa_profiles[0].start_var_idx = 0;
    msa_profiles[0].end_var_idx = 2;
    msa_profiles[0].alleles = {1, 1, 0};
    msa_profiles[0].alt_qi = {-1, -1, -1};

    VariantKeySet whitelist;
    whitelist.insert(msa_vars[0].key);
    whitelist.insert(msa_vars[2].key);
    const int admitted = merge_var_profile(
        chunk, msa_vars, msa_cats, msa_profiles, &whitelist);

    bool ok = true;
    ok &= check(admitted == 2, "two whitelisted MSA sites admitted");
    ok &= check(chunk.candidates.size() == 2, "unlisted MSA site excluded");
    ok &= check(chunk.candidates[0].counts.category == VariantCategory::NoisyCandHet,
                "whitelisted repeat collision replaced by MSA call");
    ok &= check(chunk.candidates[0].counts.alt_cov == 8,
                "MSA counts replace repeat candidate counts");
    ok &= check(chunk.candidates[0].graph_site &&
                chunk.read_var_profile[0].bam_alleles == std::vector<int>({-2, -1}) &&
                chunk.read_var_profile[0].bam_qi == std::vector<int>({12, -1}) &&
                chunk.read_var_profile[0].graph_alleles == std::vector<int>({0, -1}),
                "MSA replacement retains catalog identity and conflicting source history");
    ok &= check(chunk.candidates[0].counts.candvarcate_initial ==
                    VariantCategory::RepeatHetIndel,
                "repeat candidate provenance preserved");
    ok &= check(chunk.read_var_profile[0].start_var_idx == 0 &&
                    chunk.read_var_profile[0].end_var_idx == 1 &&
                    chunk.read_var_profile[0].alleles == std::vector<int>({1, 0}),
                "MSA profile compacted onto admitted sites");
    std::printf("%s\n", ok ? "PASS" : "FAIL");
    return ok;
}

static bool test_allele_link_orientation_and_tie() {
    PhasingChunk chunk;
    chunk.region.tid = 0;
    for (int i = 0; i < 2; ++i) {
        auto cand = dummy_cand(100 + 100 * i);
        cand.key.pos = 100 + 100 * i;
        chunk.candidates.push_back(cand);
    }
    chunk.read_var_cr.reset(cr_init());
    for (int i = 0; i < 4; ++i) {
        chunk.reads.push_back(min_read());
        chunk.haps.push_back(0);
        ReadVariantProfile profile;
        profile.start_var_idx = 0;
        profile.end_var_idx = 1;
        profile.alleles = std::vector<int>{i % 2, 1 - i % 2};
        chunk.read_var_profile.push_back(profile);
        cr_add(chunk.read_var_cr.get(), "cr", 0, 2, i);
    }
    cr_index(chunk.read_var_cr.get());
    Options opts;
    opts.link_by_alleles = true;
    opts.min_block_link_reads = 2;
    iter_update_var_hap_cons_phase_set(chunk, {0, 1}, opts);
    bool ok = check(chunk.candidates[0].phase_set == chunk.candidates[1].phase_set,
                    "untagged opposite-allele reads join blocks");
    ok &= check(chunk.candidates[1].hap_to_cons_alle[1] == 0,
                "opposite link flips consensus exactly once");
    chunk.candidates[1].hap_to_cons_alle = chunk.candidates[0].hap_to_cons_alle;
    for (int i = 0; i < 2; ++i)
        chunk.read_var_profile[i].alleles[1] = chunk.read_var_profile[i].alleles[0];
    iter_update_var_hap_cons_phase_set(chunk, {0, 1}, opts);
    ok &= check(chunk.candidates[0].phase_set != chunk.candidates[1].phase_set,
                "two agree and two conflict abstain");
    return ok;
}

static bool test_region_msa_repeat_collision() {
    PhasingChunk chunk;
    chunk.reads.push_back(min_read());
    chunk.candidates.push_back(noisy_merge_cand(100, VariantCategory::RepeatHetIndel, 3));
    const auto msa = noisy_merge_cand(100, VariantCategory::NoisyCandHet, 8);
    ReadVariantProfile profile;
    profile.start_var_idx = profile.end_var_idx = 0;
    profile.alleles = std::vector<int>{1};
    profile.alt_qi = std::vector<int>{-1};
    VariantKeySet whitelist{msa.key};
    bool ok = check(merge_var_profile(chunk, {msa}, {VariantCategory::NoisyCandHet},
                                      {profile}, &whitelist, true, true) == 0,
                    "SNP tier excludes repeat indel");
    ok &= check(merge_var_profile(chunk, {msa}, {VariantCategory::NoisyCandHet},
                                  {profile}, &whitelist, true, false) == 1,
                "indel escalation replaces existing repeat in region mode");
    ok &= check(chunk.candidates[0].counts.category == VariantCategory::NoisyCandHet &&
                chunk.read_var_profile[0].alleles == std::vector<int>{1},
                "MSA category and observations enter second phasing pass together");
    return ok;
}

static bool test_recovery_links_join_earlier_components() {
    PhasingChunk chunk;
    for (int i = 0; i < 3; ++i) {
        auto candidate = dummy_cand(100 + 100 * i);
        candidate.key.pos = 100 + 100 * i;
        chunk.candidates.push_back(candidate);
    }
    chunk.read_var_cr.reset(cr_init());
    for (int i = 0; i < 4; ++i) {
        chunk.reads.push_back(min_read());
        chunk.haps.push_back(0);
        ReadVariantProfile profile;
        profile.start_var_idx = 0;
        profile.end_var_idx = 2;
        profile.alleles = i < 2 ? std::vector<int>({0, -1, 0}) : std::vector<int>({-1, 0, 1});
        chunk.read_var_profile.push_back(profile);
        cr_add(chunk.read_var_cr.get(), "cr", 0, 3, i);
    }
    cr_index(chunk.read_var_cr.get());
    Options opts;
    opts.link_by_alleles = true;
    opts.block_link_window = 8;
    opts.min_block_link_reads = 2;
    iter_update_var_hap_cons_phase_set(chunk, {0, 1, 2}, opts);
    bool ok = check(chunk.candidates[0].phase_set != chunk.candidates[2].phase_set,
                    "nearest-only linking reproduces the missed two-component bridge");
    opts.recover_gaps = true;
    opts.private_msa_admit_all_in_region = true;
    iter_update_var_hap_cons_phase_set(chunk, {0, 1, 2}, opts);
    ok &= check(chunk.candidates[0].phase_set == chunk.candidates[1].phase_set &&
                chunk.candidates[0].phase_set == chunk.candidates[2].phase_set,
                "a later verified site joins both earlier components");
    ok &= check(chunk.candidates[0].hap_to_cons_alle[1] == chunk.candidates[2].hap_to_cons_alle[1] &&
                chunk.candidates[0].hap_to_cons_alle[1] != chunk.candidates[1].hap_to_cons_alle[1],
                "component union composes both edge orientations");
    ok &= check(iter_update_var_hap_cons_phase_set(chunk, {0, 1, 2}, opts) == 0,
                "recovery component orientations converge");
    chunk.candidates[2].is_homopolymer_indel = true;
    chunk.candidates[2].gap_link_supported = true;
    chunk.candidates[2].lcd_var_i_to_cate = kCandNoisyCandHet;
    iter_update_var_hap_cons_phase_set(chunk, {0, 1, 2}, opts);
    ok &= check(chunk.candidates[0].phase_set != chunk.candidates[1].phase_set,
                "ordinary MSA round excludes homopolymer links");
    opts.gap_hp_link_beg = 100; opts.gap_hp_link_end = 250;
    iter_update_var_hap_cons_phase_set(chunk, {0, 1, 2}, opts);
    ok &= check(chunk.candidates[0].phase_set != chunk.candidates[1].phase_set,
                "homopolymer outside the unresolved gap stays excluded");
    opts.gap_hp_link_end = 300;
    iter_update_var_hap_cons_phase_set(chunk, {0, 1, 2}, opts);
    ok &= check(chunk.candidates[0].phase_set == chunk.candidates[1].phase_set,
                "verified in-gap homopolymer can bridge both components");
    chunk.candidates[2].lcd_var_i_to_cate = kLongcalldRepHetVar;
    iter_update_var_hap_cons_phase_set(chunk, {0, 1, 2}, opts);
    ok &= check(chunk.candidates[0].phase_set != chunk.candidates[1].phase_set,
                "unverified repeat cannot enter the fallback");
    chunk.candidates[2].is_homopolymer_indel = false;
    opts.gap_hp_link_beg = opts.gap_hp_link_end = -1;
    chunk.reads.clear();
    chunk.haps.clear();
    chunk.read_var_profile.clear();
    chunk.read_var_cr.reset(cr_init());
    for (auto& candidate : chunk.candidates)
        candidate.hap_to_cons_alle = dummy_cand(100).hap_to_cons_alle;
    for (int i = 0; i < 10; ++i) {
        chunk.reads.push_back(min_read());
        chunk.haps.push_back(0);
        ReadVariantProfile profile;
        profile.start_var_idx = 0;
        profile.end_var_idx = 2;
        const int allele = i % 2;
        profile.alleles = i < 4 ? std::vector<int>({allele, allele, -1}) :
                         i < 8 ? std::vector<int>({-1, allele, allele}) :
                                 std::vector<int>({allele, -1, 1 - allele});
        chunk.read_var_profile.push_back(profile);
        cr_add(chunk.read_var_cr.get(), "cr", 0, 3, i);
    }
    cr_index(chunk.read_var_cr.get());
    iter_update_var_hap_cons_phase_set(chunk, {0, 1, 2}, opts);
    ok &= check(chunk.candidates[0].hap_to_cons_alle[1] == chunk.candidates[1].hap_to_cons_alle[1] &&
                chunk.candidates[0].hap_to_cons_alle[1] == chunk.candidates[2].hap_to_cons_alle[1],
                "a weaker conflicting cycle does not overturn two stronger links");
    for (auto& profile : chunk.read_var_profile) profile.alleles = {-1, -1, -1};
    chunk.read_var_profile[0].alleles = {0, 0, -1};
    chunk.read_var_profile[1].alleles = {1, 1, -1};
    chunk.read_var_profile[2].alleles = {0, 1, -1};
    iter_update_var_hap_cons_phase_set(chunk, {0, 1, 2}, opts);
    ok &= check(chunk.candidates[0].phase_set == chunk.candidates[1].phase_set,
                "recovery preserves the existing majority-support threshold");
    auto& repeat = chunk.candidates[1];
    repeat.is_homopolymer_indel = true;
    repeat.gap_link_supported = true;
    repeat.lcd_var_i_to_cate = kCandNoisyCandHet;
    opts.gap_hp_link_beg = 100; opts.gap_hp_link_end = 300;
    iter_update_var_hap_cons_phase_set(chunk, {0, 1, 2}, opts);
    ok &= check(chunk.candidates[0].phase_set != repeat.phase_set,
                "a 2-to-1 repeat edge fails the existing net margin");
    chunk.read_var_profile[3].alleles = {0, 0, -1};
    chunk.read_var_profile[4].alleles = {1, 1, -1};
    repeat.hap_to_cons_alle[1] = repeat.hap_to_cons_alle[2] = 1;
    repeat.counts.ref_cov = repeat.counts.alt_cov = 5;
    repeat.counts.allele_fraction = 0.5;
    iter_update_var_hap_cons_phase_set(chunk, {0, 1, 2}, opts);
    ok &= check(chunk.candidates[0].phase_set == repeat.phase_set &&
                repeat.hap_to_cons_alle[1] != repeat.hap_to_cons_alle[2],
                "verified repeat genotype survives provisional labels and a 4-to-1 edge can join");
    repeat.hap_to_cons_alle[1] = repeat.hap_to_cons_alle[2] = 1;
    repeat.counts.ref_cov = 1;
    repeat.counts.allele_fraction = 0.9;
    iter_update_var_hap_cons_phase_set(chunk, {0, 1, 2}, opts);
    ok &= check(repeat.hap_to_cons_alle[1] == repeat.hap_to_cons_alle[2],
                "an unsupported repeat genotype is not forced heterozygous");
    return ok;
}

static bool test_read_phase_set_requires_observed_allele() {
    PhasingChunk chunk;
    for (int i = 0; i < 2; ++i) {
        auto v = dummy_cand(100 + 200 * i);
        v.key.pos = 100 + 200 * i;
        v.key.type = VariantType::Snp;
        v.lcd_var_i_to_cate = kCandCleanHetSnp;
        v.counts.category = VariantCategory::CleanHetSnp;
        v.counts.n_uniq_alles = 2;
        v.counts.alle_covs = {1, 1};
        v.counts.total_cov = 2;
        chunk.candidates.push_back(v);
    }
    chunk.read_var_cr.reset(cr_init());
    for (int i = 0; i < 4; ++i) {
        chunk.reads.push_back(min_read());
        ReadVariantProfile profile;
        profile.start_var_idx = 0; profile.end_var_idx = 1;
        profile.alleles = i < 2 ? std::vector<int>({i, -1}) : std::vector<int>({-1, i - 2});
        profile.alt_qi = {-1, -1};
        chunk.read_var_profile.push_back(profile);
        cr_add(chunk.read_var_cr.get(), "cr", 0, 2, i);
    }
    cr_index(chunk.read_var_cr.get());
    Options opts;
    opts.link_by_alleles = true;
    opts.min_block_link_reads = 2;
    assign_hap_based_on_germline_het_vars_kmeans(chunk, opts, kCandGermlineClean);
    return check(chunk.candidates[0].phase_set != chunk.candidates[1].phase_set &&
                 chunk.phase_sets[2] == chunk.candidates[1].phase_set &&
                 chunk.phase_sets[3] == chunk.candidates[1].phase_set,
                 "missing allele at a spanned variant cannot assign the read to its phase block");
}

static bool test_read_phase_set_ignores_non_scoring_repeat() {
    PhasingChunk chunk;
    for (int i = 0; i < 3; ++i) {
        auto v = dummy_cand(100 + 100 * i);
        v.key.pos = 100 + 100 * i;
        v.key.type = i == 1 ? VariantType::Deletion : VariantType::Snp;
        v.is_homopolymer_indel = i == 1;
        v.lcd_var_i_to_cate = i == 1 ? kCandNoisyCandHet : kCandCleanHetSnp;
        v.counts.n_uniq_alles = 2;
        v.counts.alle_covs = {1, 1};
        v.counts.total_cov = 2;
        chunk.candidates.push_back(v);
    }
    chunk.read_var_cr.reset(cr_init());
    for (int i = 0; i < 4; ++i) {
        chunk.reads.push_back(min_read());
        ReadVariantProfile profile;
        profile.start_var_idx = 0; profile.end_var_idx = 2;
        profile.alleles = i < 2 ? std::vector<int>({i, -1, -1}) :
                                 std::vector<int>({-1, i - 2, i - 2});
        profile.alt_qi = {-1, -1, -1};
        chunk.read_var_profile.push_back(profile);
        cr_add(chunk.read_var_cr.get(), "cr", 0, 3, i);
    }
    cr_index(chunk.read_var_cr.get());
    Options opts;
    opts.link_by_alleles = true;
    assign_hap_based_on_germline_het_vars_kmeans(chunk, opts, kCandGermlineVarCate);
    bool ok = check(chunk.candidates[0].phase_set != chunk.candidates[2].phase_set &&
                    chunk.phase_sets[2] == chunk.candidates[2].phase_set &&
                    chunk.phase_sets[3] == chunk.candidates[2].phase_set,
                    "a repeat excluded from read hap scoring cannot steal phase-set ownership");
    opts.recover_gaps = opts.private_msa_admit_all_in_region = true;
    opts.gap_hp_link_beg = opts.gap_hp_link_end = 199;
    opts.min_block_link_reads = 1;
    assign_hap_based_on_germline_het_vars_kmeans(chunk, opts, kCandGermlineVarCate);
    ok &= check(chunk.candidates[1].gap_link_supported,
                "a repeat separating both haplotypes within a local block is eligible");
    chunk.haps.assign(chunk.reads.size(), 0);
    chunk.phase_sets.assign(chunk.reads.size(), -1);
    assign_hap_based_on_germline_het_vars_kmeans(chunk, opts, kCandGermlineVarCate);
    ok &= check(chunk.candidates[1].gap_link_supported,
                "direct clean-site allele association does not require prior read HP tags");
    chunk.read_var_profile[3].alleles[1] = 0;
    assign_hap_based_on_germline_het_vars_kmeans(chunk, opts, kCandGermlineVarCate);
    ok &= check(!chunk.candidates[1].gap_link_supported,
                "a repeat without within-block haplotype association is not an anchor");
    return ok;
}

static bool test_msa_edge_requires_both_haplotype_votes() {
    bool ok = true;
    for (int mode = 0; mode < 8; ++mode) {
        const bool consistent = mode & 1;
        PhasingChunk chunk;
        for (int i = 0; i < 2; ++i) {
            auto var = dummy_cand(100 + 100 * i);
            var.key.pos = 100 + 100 * i;
            var.key.type = i ? VariantType::Insertion : VariantType::Snp;
            var.key.alt = i ? "AC" : "T";
            var.lcd_var_i_to_cate = i ? kCandNoisyCandHet : kCandCleanHetSnp;
            var.msa_verified = var.gap_link_supported = i;
            var.hap_to_cons_alle = {-1, 0, 1};
            if (i && (mode & 2)) std::swap(var.hap_to_cons_alle[1], var.hap_to_cons_alle[2]);
            chunk.candidates.push_back(var);
        }
        chunk.read_var_cr.reset(cr_init());
        // Pooled 9:3 support favors a join, but the second haplotype opposes
        // that orientation 3:1. Increasing reference depth cannot resolve it.
        for (int i = 0; i < 12; ++i) {
            auto read = min_read();
            read.qname = "allele" + std::to_string(i);
            chunk.reads.push_back(std::move(read));
            ReadVariantProfile profile;
            profile.start_var_idx = 0;
            profile.end_var_idx = 1;
            profile.alleles = {i >= 8 ? 1 : 0, (consistent ? i >= 8 : i == 11) ? 1 : 0};
            if (mode & 4) std::swap(profile.alleles[0], profile.alleles[1]);
            profile.alt_qi.assign(2, -1);
            chunk.read_var_profile.push_back(profile);
            cr_add(chunk.read_var_cr.get(), "cr", 0, 2, i);
        }
        cr_index(chunk.read_var_cr.get());
        Options opts;
        opts.recover_gaps = opts.private_msa_admit_all_in_region = opts.link_by_alleles = true;
        opts.gap_recovery_beg = 150;
        opts.gap_recovery_end = 250;
        opts.min_block_link_reads = 2;
        iter_update_var_hap_cons_phase_set(chunk, {0, 1}, opts);
        ok &= check((chunk.candidates[0].phase_set == chunk.candidates[1].phase_set) == consistent,
                    "pooled reference depth cannot override contradictory MSA allele linkage");
    }
    return ok;
}

static bool test_recovery_rejects_unsupported_biallelic_indel() {
    bool ok = true;
    for (int mode = 0; mode < 2; ++mode) {
        const bool supported = mode != 0;
        PhasingChunk chunk;
        for (int i = 0; i < 3; ++i) {
            auto var = dummy_cand(100 + 100 * i);
            var.key.pos = 100 + 100 * i;
            var.key.type = i == 1 ? VariantType::Insertion : VariantType::Snp;
            var.key.alt = i == 1 ? "AC" : "T";
            var.lcd_var_i_to_cate = i == 1 ? kCandNoisyCandHet : kCandCleanHetSnp;
            var.msa_verified = i == 1;
            var.gap_link_supported = supported;
            var.hap_to_cons_alle = {-1, 0, 1};
            chunk.candidates.push_back(var);
        }
        chunk.read_var_cr.reset(cr_init());
        // Separate read groups overlap the candidate bridge, never both flanks.
        // A rejected site must not connect them merely through reference reads.
        for (int i = 0; i < 8; ++i) {
            auto read = min_read();
            read.qname = "bridge" + std::to_string(i);
            chunk.reads.push_back(std::move(read));
            ReadVariantProfile profile;
            profile.start_var_idx = 0;
            profile.end_var_idx = 2;
            const int allele = supported ? i % 2 : 0;
            const int flank_allele = supported ? i % 2 : (i % 4 == 3 ? 1 : 0);
            profile.alleles = {i < 4 ? flank_allele : -1, allele, i >= 4 ? flank_allele : -1};
            profile.alt_qi.assign(3, -1);
            chunk.read_var_profile.push_back(profile);
            cr_add(chunk.read_var_cr.get(), "cr", 0, 3, i);
        }
        cr_index(chunk.read_var_cr.get());
        Options opts;
        opts.recover_gaps = opts.private_msa_admit_all_in_region = opts.link_by_alleles = true;
        opts.gap_recovery_beg = 150;
        opts.gap_recovery_end = 250;
        opts.min_block_link_reads = 2;
        iter_update_var_hap_cons_phase_set(chunk, {0, 1, 2}, opts);
        ok &= check((chunk.candidates[0].phase_set == chunk.candidates[2].phase_set) == (mode == 1),
                    "recovery honors MSA bridge validation");
    }
    return ok;
}

static bool test_repeat_anchor_prefers_complete_evidence() {
    PhasingChunk chunk;
    for (int i = 0; i < 3; ++i) {
        auto v = dummy_cand(100 + 100 * i);
        v.key.pos = 100 + 100 * i;
        v.key.type = i == 1 ? VariantType::Deletion : VariantType::Snp;
        v.is_homopolymer_indel = i == 1;
        v.lcd_var_i_to_cate = i == 1 ? kCandNoisyCandHet : kCandCleanHetSnp;
        v.hap_to_cons_alle[1] = 0; v.hap_to_cons_alle[2] = 1;
        chunk.candidates.push_back(v);
    }
    chunk.read_var_cr.reset(cr_init());
    for (int i = 0; i < 12; ++i) {
        chunk.reads.push_back(min_read());
        ReadVariantProfile p;
        p.start_var_idx = 0; p.end_var_idx = 2;
        // The full anchor has equal repeat-allele ratios on both haplotypes.
        // A sparse subset happens to separate the repeat perfectly.
        const int repeat_allele = i % 6 < 2 ? 1 : 0;
        p.alleles = {i < 6 ? 0 : 1, repeat_allele, i < 4 ? repeat_allele : -1};
        p.alt_qi = {-1, -1, -1};
        chunk.read_var_profile.push_back(p);
        cr_add(chunk.read_var_cr.get(), "cr", 0, 3, i);
    }
    cr_index(chunk.read_var_cr.get());
    chunk.haps.assign(12, 0); chunk.phase_sets.assign(12, -1);
    Options opts;
    opts.recover_gaps = opts.private_msa_admit_all_in_region = opts.link_by_alleles = true;
    opts.gap_hp_link_beg = opts.gap_hp_link_end = 199;
    opts.min_block_link_reads = 2;
    assign_hap_based_on_germline_het_vars_kmeans(chunk, opts, kCandGermlineVarCate);
    return check(!chunk.candidates[1].gap_link_supported,
                 "a favorable sparse subset cannot override the fuller clean-site evidence");
}

static bool test_msa_unsorted_profile_indices() {
    PhasingChunk chunk;
    chunk.reads.push_back(min_read());
    ReadVariantProfile profile;
    profile.start_var_idx = 0;
    profile.end_var_idx = 1;
    profile.alleles = std::vector<int>{1, 0};
    profile.alt_qi = std::vector<int>{17, -1};
    profile.graph_alleles = {0, 1};
    merge_var_profile(chunk,
                      {noisy_merge_cand(200, VariantCategory::NoisyCandHet, 8),
                       noisy_merge_cand(100, VariantCategory::NoisyCandHet, 8)},
                      {VariantCategory::NoisyCandHet, VariantCategory::NoisyCandHet},
                      {profile});
    return check(chunk.candidates[0].key.pos == 100 &&
                 chunk.read_var_profile[0].alleles == std::vector<int>({0, 1}) &&
                 chunk.read_var_profile[0].alt_qi == std::vector<int>({-1, 17}) &&
                 chunk.read_var_profile[0].graph_alleles == std::vector<int>({1, 0}),
                 "sorting MSA calls preserves allele and query-position ownership");
}

static void make_gap_fixture(std::vector<PhasingChunk>& chunks, PhasingChunk& proposal) {
    chunks.resize(2);
    for (int ci = 0; ci < 2; ++ci) {
        auto& chunk = chunks[ci];
        chunk.region.tid = 0;
        chunk.region.beg = ci == 0 ? 1 : 250;
        chunk.region.end = ci == 0 ? 249 : 600;
        const hts_pos_t ps = ci == 0 ? 100 : 300;
        auto candidate = dummy_cand(ps);
        candidate.key.tid = 0;
        candidate.key.pos = ps;
        candidate.key.type = VariantType::Snp;
        candidate.counts.category = VariantCategory::CleanHetSnp;
        candidate.lcd_var_i_to_cate = kCandCleanHetSnp;
        chunk.candidates.push_back(candidate);
        for (int j = 0; j < 2; ++j) {
            auto read = min_read();
            read.qname = "read" + std::to_string(ci * 2 + j);
            chunk.reads.push_back(std::move(read));
            chunk.haps.push_back(j + 1);
            chunk.phase_sets.push_back(ps);
        }
    }
    auto extra = min_read();
    extra.qname = "extra";
    chunks[0].reads.push_back(std::move(extra));
    chunks[0].haps.push_back(0);
    chunks[0].phase_sets.push_back(-1);
    auto unrelated = min_read();
    unrelated.qname = "unrelated";
    chunks[1].reads.push_back(std::move(unrelated));
    chunks[1].haps.push_back(1);
    chunks[1].phase_sets.push_back(500);
    auto other_site = chunks[1].candidates[0];
    other_site.key.pos = other_site.phase_set = 500;
    chunks[1].candidates.push_back(other_site);

    proposal.region.tid = 0;
    for (int i = 0; i < 5; ++i) {
        auto read = min_read();
        read.qname = i == 4 ? "extra" : "read" + std::to_string(i);
        read.n_clean_agree_snps = 2;
        proposal.reads.push_back(std::move(read));
        proposal.haps.push_back(i < 2 ? i + 1 : i < 4 ? 4 - i : 1);
        proposal.phase_sets.push_back(200);
        ReadVariantProfile profile;
        profile.start_var_idx = profile.end_var_idx = 0;
        profile.alleles.push_back(i % 2);
        profile.alt_qi.push_back(-1);
        proposal.read_var_profile.push_back(profile);
    }
    auto site = chunks[0].candidates[0];
    site.key.pos = site.phase_set = 200;
    proposal.candidates.push_back(site);
}

static bool test_gap_recovery_join_preserves_blocks() {
    std::vector<PhasingChunk> chunks;
    PhasingChunk proposal;
    make_gap_fixture(chunks, proposal);
    proposal.candidates[0].lcd_var_i_to_cate = kCandNonAnchorHet;
    const auto gaps = find_phase_gaps(chunks);
    bool ok = check(gaps.size() == 2 && gaps[0].left_end == 100 && gaps[0].right_beg == 300,
                    "find internal gaps across chunk boundaries");
    Options opts;
    opts.stitch_rule = kStitchRuleBothStrands;
    std::vector<GapLinkEvidence> evidence;
    const auto edge_only = stitch_gap_proposal(chunks, proposal, gaps[0], opts,
                                               nullptr, true, true, &evidence);
    ok &= check(evidence.size() == 1 && evidence[0].proposal_ps == 200 &&
                evidence[0].votes[0] == std::array<int, 4>{1, 0, 0, 1} &&
                evidence[0].votes[1] == std::array<int, 4>{0, 1, 1, 0},
                "read-only audit preserves both original flank vote matrices");
    ok &= check(edge_only.joined && edge_only.right_flip && edge_only.reads_added == 0 &&
                chunks[0].haps == std::vector<int>({1, 2, 0}) &&
                chunks[1].phase_sets == std::vector<hts_pos_t>({300, 300, 500}),
                "orientation-only recovery reports a flip without changing reads or blocks");
    Options hp_opts = opts;
    hp_opts.gap_hp_link_beg = gaps[0].left_end;
    hp_opts.gap_hp_link_end = gaps[0].right_beg;
    hp_opts.min_block_link_reads = 1;
    const auto hp_edge = stitch_gap_proposal(chunks, proposal, gaps[0], hp_opts,
                                            nullptr, true, true);
    ok &= check(hp_edge.joined && hp_edge.right_flip && hp_edge.reads_added == 0 &&
                chunks[0].haps == std::vector<int>({1, 2, 0}) &&
                chunks[1].phase_sets == std::vector<hts_pos_t>({300, 300, 500}) &&
                chunks[0].candidates.size() == 1,
                "homopolymer confirmation cannot mutate either flank before acceptance");
    const auto result = stitch_gap_proposal(chunks, proposal, gaps[0], opts);
    ok &= check(result.joined && result.reads_added == 1, "bridge joins both flanks and adds an unphased read");
    ok &= check(chunks[0].haps == std::vector<int>({1, 2, 1}), "left block orientation retained");
    ok &= check(chunks[1].haps == std::vector<int>({2, 1, 1}), "right block flips uniformly; unrelated read unchanged");
    ok &= check(chunks[1].phase_sets == std::vector<hts_pos_t>({100, 100, 500}), "only target phase set relabelled");
    ok &= check(chunks[0].candidates.size() == 2 && chunks[0].candidates[1].phase_set == 100,
                "recovery site enters output in the stitched orientation");
    ok &= check(chunks[0].candidates[1].lcd_var_i_to_cate == kCandNonAnchorHet,
                "recovery does not promote non-anchor candidates");
    ok &= check(chunks[1].candidates[0].hap_to_cons_alle[1] == 0,
                "candidate orientation follows the block flip");
    return ok;
}

static bool test_gap_recovery_partial_extension() {
    std::vector<PhasingChunk> chunks;
    PhasingChunk proposal;
    make_gap_fixture(chunks, proposal);
    proposal.haps[3] = 0;
    Options opts;
    opts.stitch_rule = kStitchRuleBothStrands;
    const auto result = stitch_gap_proposal(chunks, proposal, find_phase_gaps(chunks)[0], opts);
    return check(result.left_linked && !result.right_linked && !result.joined && result.reads_added == 1 &&
                 chunks[1].phase_sets[0] == 300 && chunks[1].haps[0] == 1,
                 "one-sided evidence preserves an extension but cannot close the gap");
}

static bool test_gap_independent_emits_every_qualifying_group() {
    // A wide gap can hold several independent internal blocks. Emitting only the
    // largest group discards the rest, and it also picks badly: the largest
    // gap-only group is usually the flank-adjacent component, whose reads the
    // same recovery round has already attached. Both groups here clear
    // min_reads, so both must be emitted.
    std::vector<PhasingChunk> chunks(1);
    PhasingChunk& chunk = chunks[0];
    chunk.region.tid = 0;
    chunk.region.beg = 1;
    chunk.region.end = 1000;
    for (int i = 0; i < 8; ++i) {
        auto read = min_read();
        read.qname = "gapread" + std::to_string(i);
        chunk.reads.push_back(std::move(read));
        chunk.haps.push_back(0);
        chunk.phase_sets.push_back(-1);
    }
    const GapReadIndex index(chunks);

    PhasingChunk proposal;
    proposal.region.tid = 0;
    for (int i = 0; i < 8; ++i) {
        auto read = min_read();
        read.qname = "gapread" + std::to_string(i);
        proposal.reads.push_back(std::move(read));
        proposal.haps.push_back(i % 2 + 1);
        proposal.phase_sets.push_back(i < 5 ? 400 : 700);
    }

    PhaseGap gap{};
    gap.tid = 0;
    gap.left_end = 100;
    gap.right_beg = 900;
    gap.left_ps = 100;
    gap.right_ps = 900;
    gap.region_beg = 1;
    gap.region_end = 1000;
    Options opts;
    hts_pos_t emitted = -1;
    std::set<hts_pos_t> emitted_round;
    const int applied = emit_independent_gap_block(chunks, proposal, gap, opts, index, 3,
                                                   &emitted, &emitted_round);
    bool ok = check(applied == 8 && emitted_round.size() == 2 &&
                    emitted_round.count(400) == 1 && emitted_round.count(700) == 1,
                    "every gap-only group clearing min_reads is emitted, not just the largest");
    ok &= check(emitted == 400, "emitted_ps still reports the largest group");
    ok &= check(chunk.phase_sets[0] == 400 && chunk.phase_sets[7] == 700 &&
                chunk.haps[0] == 1 && chunk.haps[7] == 2,
                "reads of both groups carry their own block and haplotype");

    // A group below min_reads is still screened out.
    std::vector<PhasingChunk> strict_chunks(1);
    PhasingChunk& strict = strict_chunks[0];
    strict.region = chunk.region;
    for (int i = 0; i < 8; ++i) {
        auto read = min_read();
        read.qname = "gapread" + std::to_string(i);
        strict.reads.push_back(std::move(read));
        strict.haps.push_back(0);
        strict.phase_sets.push_back(-1);
    }
    const GapReadIndex strict_index(strict_chunks);
    hts_pos_t strict_emitted = -1;
    std::set<hts_pos_t> strict_round;
    const int strict_applied = emit_independent_gap_block(
        strict_chunks, proposal, gap, opts, strict_index, 4, &strict_emitted, &strict_round);
    ok &= check(strict_applied == 5 && strict_round.size() == 1 && strict_round.count(400) == 1,
                "a group below min_reads is not emitted");
    return ok;
}

static bool test_gap_allele_attach_join_only_gates_one_sided() {
    // Isolate the allele-agreement path. `extra` is unphased in the chunk AND in
    // the proposal, so the proposal-based attachment cannot claim it and only
    // implied_assignment can; chunk 0 gains read profiles because the shared
    // fixture has none, which is why the other one-sided tests never reach this
    // path. Zeroing proposal read 3 leaves the right flank unlinked, so the join
    // is one-sided.
    auto build = [](std::vector<PhasingChunk>& chunks, PhasingChunk& proposal) {
        make_gap_fixture(chunks, proposal);
        proposal.haps[3] = 0;
        proposal.haps[4] = 0;
        chunks[0].read_var_profile.resize(chunks[0].reads.size());
        for (auto& prof : chunks[0].read_var_profile) {
            prof.start_var_idx = prof.end_var_idx = 0;
            prof.alleles.assign(1, -1);
            prof.alt_qi.assign(1, -1);
        }
        const size_t extra_i = chunks[0].reads.size() - 1;
        chunks[0].read_var_profile[extra_i].alleles[0] =
            chunks[0].candidates[0].hap_to_cons_alle[1];
        return extra_i;
    };

    std::vector<PhasingChunk> open_chunks;
    PhasingChunk open_proposal;
    const size_t extra_i = build(open_chunks, open_proposal);
    Options open_opts;
    open_opts.stitch_rule = kStitchRuleBothStrands;
    const auto open_result = stitch_gap_proposal(open_chunks, open_proposal,
                                                 find_phase_gaps(open_chunks)[0], open_opts);
    bool ok = check(!open_result.joined && open_result.left_linked &&
                    open_chunks[0].haps[extra_i] == 1 &&
                    open_chunks[0].phase_sets[extra_i] == 100,
                    "allele agreement attaches to a one-sided link by default");

    std::vector<PhasingChunk> gated_chunks;
    PhasingChunk gated_proposal;
    build(gated_chunks, gated_proposal);
    Options gated_opts = open_opts;
    gated_opts.gap_allele_attach_join_only = true;
    const auto gated_result = stitch_gap_proposal(gated_chunks, gated_proposal,
                                                  find_phase_gaps(gated_chunks)[0], gated_opts);
    ok &= check(!gated_result.joined && gated_result.left_linked &&
                gated_chunks[0].haps[extra_i] == 0 &&
                gated_chunks[0].phase_sets[extra_i] == -1,
                "allele attachment is withheld from a one-sided link when it must follow a join");
    return ok;
}

static bool test_gap_read_index_freezes_initial_assignments() {
    std::vector<PhasingChunk> chunks;
    PhasingChunk proposal;
    make_gap_fixture(chunks, proposal);
    const GapReadIndex index(chunks);
    Options opts;
    opts.stitch_rule = kStitchRuleBothStrands;
    proposal.haps[3] = 0;
    const auto first = stitch_gap_proposal(chunks, proposal, find_phase_gaps(chunks)[0], opts, &index);
    bool ok = check(first.reads_added == 1 && !first.joined,
                    "cached index permits an initial partial extension");
    // The next proposal disagrees on the newly tagged read. It must not turn
    // into anchor evidence merely because an earlier recovery touched it.
    proposal.haps[4] = 2;
    proposal.haps[3] = 1;
    const auto second = stitch_gap_proposal(chunks, proposal, find_phase_gaps(chunks)[0], opts, &index);
    ok &= check(second.joined && second.reads_added == 0 && chunks[0].haps[2] == 1 &&
                    index.assignments.size() == 5,
                "reused index keeps only pre-recovery anchor assignments");
    return ok;
}

static bool test_gap_msa_covers_clean_stretches() {
    PhasingChunk chunk;
    chunk.ref_beg = 1;
    chunk.ref_end = 2000;
    chunk.noisy_regions.push_back({700, 800, 3});
    prepare_gap_msa_regions(chunk, 100, 1300);
    bool ok = true;
    for (int pos = 100; pos <= 1300; ++pos) {
        bool covered = false;
        for (const auto& region : chunk.noisy_regions)
            covered |= region.beg <= pos && region.end >= pos;
        ok &= covered;
    }
    return check(ok && chunk.noisy_regions[0].beg == 700 && chunk.noisy_regions[0].end == 800,
                 "gap-driven MSA covers clean stretches and preserves noisy windows");
}

static bool test_gap_recovery_prefers_complete_bridge() {
    std::vector<PhasingChunk> chunks;
    PhasingChunk proposal;
    make_gap_fixture(chunks, proposal);
    for (int i = 0; i < 4; ++i) {
        auto core_read = min_read();
        core_read.qname = "left_only" + std::to_string(i);
        auto proposal_read = min_read();
        proposal_read.qname = core_read.qname;
        chunks[0].reads.push_back(std::move(core_read));
        chunks[0].haps.push_back(1 + i % 2);
        chunks[0].phase_sets.push_back(100);
        proposal.reads.push_back(std::move(proposal_read));
        proposal.haps.push_back(1 + i % 2);
        proposal.phase_sets.push_back(150);
        proposal.read_var_profile.emplace_back();
    }
    Options opts;
    opts.stitch_rule = kStitchRuleBothStrands;
    const auto result = stitch_gap_proposal(chunks, proposal, find_phase_gaps(chunks)[0], opts);
    return check(result.joined, "a complete bridge wins over stronger disconnected endpoint matches");
}

static bool test_gap_recovery_deduplicates_votes() {
    std::vector<PhasingChunk> chunks;
    PhasingChunk proposal;
    make_gap_fixture(chunks, proposal);
    auto duplicate = min_read();
    duplicate.qname = "read0";
    proposal.reads.push_back(std::move(duplicate));
    proposal.haps.push_back(1);
    proposal.phase_sets.push_back(200);
    proposal.read_var_profile.emplace_back();
    Options opts;
    opts.stitch_rule = kStitchRuleNetMargin;
    opts.stitch_min_margin = 2;
    const auto result = stitch_gap_proposal(chunks, proposal, find_phase_gaps(chunks)[0], opts);
    return check(!result.left_linked && !result.right_linked && result.reads_added == 0,
                 "duplicate alignments cannot manufacture sufficient endpoint votes");
}

static bool test_gap_tiers_are_cumulative() {
    PhasingChunk chunk;
    chunk.reads.push_back(min_read());
    auto snp = noisy_merge_cand(100, VariantCategory::NoisyCandHet, 8);
    snp.key.type = VariantType::Snp;
    snp.key.alt = "A";
    const auto indel = noisy_merge_cand(200, VariantCategory::NoisyCandHet, 8);
    ReadVariantProfile profile;
    profile.start_var_idx = 0;
    profile.end_var_idx = 1;
    profile.alleles = std::vector<int>{1, 0};
    profile.alt_qi = std::vector<int>{10, -1};
    const std::vector<CandidateVariant> sites{snp, indel};
    const std::vector<VariantCategory> categories(2, VariantCategory::NoisyCandHet);
    merge_var_profile(chunk, sites, categories, {profile}, nullptr, true, true);
    bool ok = check(chunk.candidates.size() == 1 && chunk.candidates[0].key.type == VariantType::Snp,
                    "tier two admits only SNPs");
    merge_var_profile(chunk, sites, categories, {profile}, nullptr, true, false);
    ok &= check(chunk.candidates.size() == 2 && chunk.read_var_profile[0].alleles == std::vector<int>({1, 0}),
                "tier three adds indels while retaining tier-two SNP observations");
    return ok;
}

static bool test_gap_recovery_respects_requested_regions() {
    std::vector<PhasingChunk> chunks;
    PhasingChunk proposal;
    make_gap_fixture(chunks, proposal);
    chunks[0].region.end = 150;
    const auto gaps = find_phase_gaps(chunks);
    return check(gaps.size() == 1 && gaps[0].left_end == 300 && gaps[0].region_beg == 250,
                 "do not recover excluded intervals between separate requested regions");
}

static bool test_gap_recovery_adopts_unresolved_site_and_read() {
    std::vector<PhasingChunk> chunks;
    PhasingChunk proposal;
    make_gap_fixture(chunks, proposal);
    chunks[0].haps[2] = 1; // An HP without a PS is still unresolved.
    auto unphased_site = proposal.candidates[0];
    unphased_site.phase_set = -1;
    unphased_site.counts.category = VariantCategory::LowCoverage;
    unphased_site.lcd_var_i_to_cate = kLongcalldLowCovVar;
    unphased_site.hap_to_cons_alle = {-1, -1, -1};
    chunks[0].candidates.push_back(unphased_site);
    Options opts;
    opts.stitch_rule = kStitchRuleBothStrands;
    const auto result = stitch_gap_proposal(chunks, proposal, find_phase_gaps(chunks)[0], opts);
    return check(result.joined && result.reads_added == 1 && chunks[0].phase_sets[2] == 100 &&
                 chunks[0].candidates[1].phase_set == 100 && chunks[0].candidates[1].hap_to_cons_alle[1] == 1 &&
                 chunks[0].candidates[1].counts.category == VariantCategory::CleanHetSnp,
                 "previously unphased exact sites and reads receive consistent recovered phase");
}

static bool test_gap_recovery_keeps_unphased_observations() {
    std::vector<PhasingChunk> chunks;
    PhasingChunk proposal;
    make_gap_fixture(chunks, proposal);
    proposal.haps[4] = 0;
    Options opts;
    opts.stitch_rule = kStitchRuleBothStrands;
    const auto result = stitch_gap_proposal(chunks, proposal, find_phase_gaps(chunks)[0], opts);
    const auto& profile = chunks[0].read_var_profile[2];
    return check(result.joined && result.reads_added == 0 && chunks[0].haps[2] == 0 &&
                 profile.start_var_idx == 1 && profile.alleles == std::vector<int>({0}),
                 "untagged reads retain raw observations without receiving an unsupported HP tag");
}

static bool test_gap_hp_trial_is_transactional() {
    std::vector<PhasingChunk> chunks;
    PhasingChunk proposal;
    make_gap_fixture(chunks, proposal);
    const auto gap = find_phase_gaps(chunks)[0];
    Options opts;
    opts.gap_hp_link_beg = gap.left_end;
    opts.gap_hp_link_end = gap.right_beg;
    opts.min_block_link_reads = 2;
    const auto rejected = stitch_gap_proposal(chunks, proposal, gap, opts);
    bool ok = check(!rejected.joined && rejected.reads_added == 0 &&
                    chunks[0].haps.back() == 0 && chunks[1].phase_sets[0] == gap.right_ps,
                    "weak homopolymer trial leaves original blocks and untagged reads untouched");
    opts.min_block_link_reads = 1;
    proposal.phase_sets[2] = proposal.phase_sets[3] = 400;
    const auto partial = stitch_gap_proposal(chunks, proposal, gap, opts);
    ok &= check(!partial.joined && partial.reads_added == 0 && chunks[0].haps.back() == 0,
                "two disconnected anchors cannot commit a homopolymer extension");
    proposal.phase_sets[2] = proposal.phase_sets[3] = 200;
    const auto joined = stitch_gap_proposal(chunks, proposal, gap, opts);
    ok &= check(joined.joined && joined.reads_added == 1 &&
                chunks[1].phase_sets[0] == gap.left_ps,
                "supported complete homopolymer proposal uses the existing stitch path");
    return ok;
}

static AlnStr site_alignment(const std::string& target, const std::string& query) {
    AlnStr aln;
    for (char c : target) aln.target_aln.push_back(c == '-' ? 5 : std::string("ACGTN").find(c));
    for (char c : query) aln.query_aln.push_back(c == '-' ? 5 : std::string("ACGTN").find(c));
    aln.aln_len = static_cast<int>(target.size());
    aln.target_end = aln.query_end = aln.aln_len - 1;
    return aln;
}

static bool test_msa_site_observations() {
    VariantKey key;
    key.pos = 103;
    key.type = VariantType::Snp;
    key.ref_len = 1;
    key.alt = "T";
    auto alt = site_alignment("ACGACGT", "ACGTCGT");
    auto ref = site_alignment("ACGACGT", "ACGACGT");
    bool ok = check(call_msa_site_allele({alt, alt}, key, 100) == 1,
                    "SNP observation does not require choosing a whole consensus");
    ok &= check(call_msa_site_allele({ref, ref}, key, 100) == 0, "direct reference observation");
    ok &= check(call_msa_site_allele({alt, ref}, key, 100) == -1, "alignment disagreement abstains");
    auto noisy = site_alignment("ACGACGT", "ATGTCGT");
    ok &= check(call_msa_site_allele({noisy, noisy}, key, 100) == -1, "nonmatching flank abstains");
    alt.query_end = 5;
    ok &= check(call_msa_site_allele({alt, alt}, key, 100) == -1, "partial flank abstains");
    key.type = VariantType::Insertion;
    key.ref_len = 0;
    key.alt = "TT";
    alt = site_alignment("ACG--CGT", "ACGTTCGT");
    ref = site_alignment("ACGCGT", "ACGCGT");
    ok &= check(call_msa_site_allele({alt, alt}, key, 100) == 1, "exact insertion observation");
    ok &= check(call_msa_site_allele({ref, ref}, key, 100) == 0, "insertion reference observation");
    auto other = site_alignment("ACG---CGT", "ACGTTTCGT");
    ok &= check(call_msa_site_allele({other, other}, key, 100) == -1, "third insertion allele abstains");
    key.type = VariantType::Deletion;
    key.ref_len = 2;
    key.alt.clear();
    alt = site_alignment("ACGTTCGT", "ACG--CGT");
    ok &= check(call_msa_site_allele({alt, alt}, key, 100) == 1, "exact deletion observation");
    return ok;
}

static bool test_composed_msa_repeat_placement() {
    VariantKey key;
    key.pos = 103; key.type = VariantType::Insertion; key.alt = "GT";
    auto left = site_alignment("ACG--GTGTGTCGA", "ACGGTGTGTGTCGA");
    auto right = site_alignment("ACGGTGT--GTCGA", "ACGGTGTGTGTCGA");
    bool ok = check(call_msa_site_allele({left, right}, key, 100) == -1,
                    "equivalent repeat placements reproduce the missing observation");
    left_normalize_msa_alignment(right);
    ok &= check(right.target_aln == left.target_aln && right.query_aln == left.query_aln &&
                call_msa_site_allele({left, right}, key, 100) == 1,
                "canonical indel placement restores the exact insertion observation");
    left = site_alignment("ACGATATATCGA", "ACG--ATATCGA");
    right = site_alignment("ACGATATATCGA", "ACGATAT--CGA");
    left_normalize_msa_alignment(right);
    ok &= check(right.query_aln == left.query_aln && right.target_aln == left.target_aln,
                "deletion normalization preserves both ungapped sequences");
    return ok;
}

static bool test_msa_cluster_membership_is_not_an_allele() {
    PhasingChunk chunk;
    chunk.region.tid = 0; chunk.ref_beg = 100; chunk.ref_seq = "ACGACGT";
    chunk.reads.resize(1);
    const auto ref = site_alignment("ACGACGT", "ACGACGT");
    const auto alt = site_alignment("ACGACGT", "ACGTCGT");
    std::array<std::vector<AlnStr>, 2> alignments = {{{alt}, {ref, alt, {}}}};
    std::vector<CandidateVariant> vars;
    std::vector<VariantCategory> categories;
    std::vector<ReadVariantProfile> profiles;
    Options opts;
    opts.recover_gaps = true;
    make_vars_from_msa_cons_aln(opts, chunk, 1, {}, 100, 2, {0, 1}, {{{}, {0}}},
                               alignments, vars, categories, profiles);
    bool ok = check(vars.size() == 1 && profiles[0].alleles == std::vector<int>({1}) &&
                    vars[0].counts.ref_cov == 0 && vars[0].counts.alt_cov == 1,
                    "a reference-cluster read carrying ALT is not imputed as reference");
    alignments[1][1] = site_alignment("ACGACGT", "ACTTCGT");
    make_vars_from_msa_cons_aln(opts, chunk, 1, {}, 100, 2, {0, 1}, {{{}, {0}}},
                               alignments, vars, categories, profiles);
    ok &= check(profiles[0].alleles == std::vector<int>({1}),
                "one flank error does not erase a strictly better local allele");
    alignments[1][1] = site_alignment("ACGACGT", "ACGGCGT");
    make_vars_from_msa_cons_aln(opts, chunk, 1, {}, 100, 2, {0, 1}, {{{}, {0}}},
                               alignments, vars, categories, profiles);
    ok &= check(profiles[0].alleles == std::vector<int>({-1}) && vars[0].counts.total_cov == 0,
                "a third allele is unknown rather than an invented reference vote");
    return ok;
}

static bool test_assigned_repeat_allows_one_local_error() {
    PhasingChunk chunk;
    chunk.region.tid = 0; chunk.ref_beg = 100; chunk.ref_seq = "ACGAAACGT";
    chunk.reads.resize(1);
    const auto short_del = site_alignment(chunk.ref_seq, "ACG-AACGT");
    const auto long_del = site_alignment(chunk.ref_seq, "ACG--ACGT");
    const auto read = site_alignment("ACGACGT", "ACG-CGT");
    const std::array<std::vector<AlnStr>, 2> alignments = {{{short_del}, {long_del, read, {}}}};
    std::vector<CandidateVariant> vars;
    std::vector<VariantCategory> categories;
    std::vector<ReadVariantProfile> profiles;
    Options opts;
    opts.recover_gaps = true;
    make_vars_from_msa_cons_aln(opts, chunk, 1, {}, 100, 2, {0, 1}, {{{}, {0}}},
                               alignments, vars, categories, profiles);
    return check(vars.size() == 2 && vars[1].key.pos == 104 &&
                 profiles[0].alleles[1] == 1 && vars[1].counts.alt_cov == 1,
                 "a one-base repeat error supports the strictly closer length allele");
}

static bool test_unassigned_msa_local_allele_evidence() {
    CandidateVariant var;
    var.key.pos = 103; var.key.type = VariantType::Snp;
    var.key.ref_len = 1; var.key.alt = "T";
    var.counts.category = VariantCategory::NoisyCandHet;
    var.counts.alle_covs = {2, 2}; var.counts.total_cov = 4;
    const auto ref = site_alignment("ACGACGT", "ACGACGT");
    const auto alt = site_alignment("ACGACGT", "ACGTCGT");
    const auto noisy_alt = site_alignment("ACGACGT", "ACTTCGT");
    const auto third = site_alignment("ACGACGT", "ACGGCGT");
    const std::array<AlnStr, 2> consensuses = {ref, alt};
    std::vector<CandidateVariant> vars{var};
    std::vector<ReadVariantProfile> profiles(3);
    Options opts;
    add_msa_site_observations(opts,
        {{0, {noisy_alt, noisy_alt}}, {1, {noisy_alt, ref}}, {2, {third, third}}},
        100, false, vars, profiles, &consensuses);
    return check(profiles[0].alleles == std::vector<int>({1}) &&
                 profiles[1].start_var_idx == -1 && profiles[2].start_var_idx == -1 &&
                 vars[0].counts.total_cov == 5,
                 "unassigned reads retain local evidence only when both paths uniquely agree");
}

static bool test_msa_counts_do_not_depend_on_recover_gaps() {
    // A candidate's allele counts describe the reads, so the same MSA input must
    // produce the same counts whether or not gap recovery is enabled. It does
    // not today: refresh_assigned_msa_observations runs only when recover_gaps
    // is set (collect_phase_noisy.cpp), so the hybrid arm re-scores MSA site
    // observations and the alignment-only channel never does. On
    // chr20:48,176,830-48,229,446 eight of the seventy-six candidates shared
    // between the two channels carry different DP/REF_COUNT/ALT_COUNT, which is
    // what this pins. The fixture deliberately has no nested deletion, so
    // split_nested_msa_deletions is a no-op and the refresh is the only
    // difference between the two runs.
    const auto build = [](bool recover_gaps,
                          std::vector<CandidateVariant>& vars,
                          std::vector<VariantCategory>& categories,
                          std::vector<ReadVariantProfile>& profiles) {
        PhasingChunk chunk;
        chunk.region.tid = 0;
        chunk.ref_beg = 100;
        chunk.ref_seq = "ACGTACGTACGT";
        chunk.reads.resize(2);
        const auto hap1 = site_alignment(chunk.ref_seq, "ACGTTCGTACGT");
        const auto hap2 = site_alignment(chunk.ref_seq, "ACGTACGTACGT");
        const std::array<std::vector<AlnStr>, 2> alignments = {{{hap1}, {hap2}}};
        Options opts;
        opts.recover_gaps = recover_gaps;
        // Assigned read clusters, because refresh_assigned_msa_observations only
        // re-scores reads that a cluster owns: with empty clusters it is a no-op
        // and the two arms agree trivially.
        chunk.reads.resize(4);
        const std::array<int, 2> clu_n_seqs = {2, 2};
        const std::array<std::vector<int>, 2> clu_read_ids = {{{0, 1}, {2, 3}}};
        const std::array<std::vector<AlnStr>, 2> read_alns = {{{hap1, hap1}, {hap2, hap2}}};
        make_vars_from_msa_cons_aln(opts, chunk, 4, {0, 1, 2, 3}, 100, 2,
                                   clu_n_seqs, clu_read_ids,
                                   read_alns, vars, categories, profiles);
        const std::array<AlnStr, 2> consensuses = {hap1, hap2};
        add_msa_site_observations(opts, {{0, {hap1, hap1}}, {1, {hap2, hap2}}},
                                  100, false, vars, profiles, &consensuses);
    };
    std::vector<CandidateVariant> off_vars, on_vars;
    std::vector<VariantCategory> off_cate, on_cate;
    std::vector<ReadVariantProfile> off_prof, on_prof;
    build(false, off_vars, off_cate, off_prof);
    build(true, on_vars, on_cate, on_prof);
    bool ok = check(off_vars.size() == on_vars.size(),
                    "the same MSA input yields the same candidates either way");
    if (!ok) return false;
    for (size_t vi = 0; vi < off_vars.size(); ++vi) {
        const VariantCounts& a = off_vars[vi].counts;
        const VariantCounts& b = on_vars[vi].counts;
        ok &= check(a.total_cov == b.total_cov && a.ref_cov == b.ref_cov &&
                    a.alt_cov == b.alt_cov && a.alle_covs == b.alle_covs &&
                    off_cate[vi] == on_cate[vi],
                    "candidate counts and category do not depend on recover_gaps");
    }
    return ok;
}

static bool test_nested_msa_deletions_share_common_event() {
    PhasingChunk chunk;
    chunk.region.tid = 0;
    chunk.ref_beg = 100;
    chunk.ref_seq = "ACGTATATACGT";
    chunk.reads.resize(2);
    const auto short_del = site_alignment(chunk.ref_seq, "ACG----TACGT");
    const auto long_del = site_alignment(chunk.ref_seq, "ACG------CGT");
    const std::array<std::vector<AlnStr>, 2> alignments = {{{short_del}, {long_del}}};
    std::vector<CandidateVariant> vars;
    std::vector<VariantCategory> categories;
    std::vector<ReadVariantProfile> profiles;
    Options opts;
    opts.recover_gaps = true;
    make_vars_from_msa_cons_aln(opts, chunk, 2, {}, 100, 2, {0, 0}, {},
                               alignments, vars, categories, profiles);
    bool ok = check(vars.size() == 2 && vars[0].key.pos == 103 &&
                    vars[0].key.ref_len == 4 && categories[0] == VariantCategory::NoisyCandHom &&
                    vars[1].key.pos == 107 && vars[1].key.ref_len == 2 &&
                    categories[1] == VariantCategory::NoisyCandHet,
                    "4/6 deletion genotype becomes common 4-base deletion plus 2-base difference");
    if (!ok) return false;
    const std::array<AlnStr, 2> consensuses = {short_del, long_del};
    add_msa_site_observations(opts, {{0, {short_del, short_del}}, {1, {long_del, long_del}}},
                              100, false, vars, profiles, &consensuses);
    ok &= check(profiles[0].alleles == std::vector<int>({0}) &&
                profiles[1].alleles == std::vector<int>({1}) &&
                vars[1].counts.ref_cov == 1 && vars[1].counts.alt_cov == 1,
                "the two repeat lengths give opposite observations at the heterozygous difference");
    return ok;
}

static bool test_msa_observation_heterozygosity_gate() {
    CandidateVariant var;
    var.key.pos = 103;
    var.key.type = VariantType::Snp;
    var.key.ref_len = 1;
    var.key.alt = "T";
    var.counts.category = VariantCategory::NoisyCandHet;
    var.counts.alle_covs = {1, 1};
    var.counts.total_cov = 2;
    std::vector<CandidateVariant> vars{var};
    std::vector<ReadVariantProfile> profiles(10);
    std::vector<UnassignedMsaRead> reads;
    const auto ref = site_alignment("ACGACGT", "ACGACGT");
    for (int i = 0; i < 10; ++i) reads.push_back({i, {ref, ref}});
    Options opts;
    add_msa_site_observations(opts, reads, 100, true, vars, profiles);
    bool ok = check(vars[0].counts.total_cov == 2 && profiles[0].start_var_idx == -1,
                    "reference-dominated evidence does not extend a false MSA het");
    reads.resize(1);
    add_msa_site_observations(opts, reads, 100, true, vars, profiles);
    ok &= check(vars[0].counts.total_cov == 3 && vars[0].counts.alle_covs[0] == 2 &&
                profiles[0].alleles == std::vector<int>({0}),
                "balanced site observation updates counts and profile exactly once");
    return ok;
}

static bool test_msa_two_alternate_insertions() {
    PhasingChunk chunk;
    chunk.region.tid = 0; chunk.ref_beg = 100; chunk.ref_seq = "ACGGTA";
    chunk.reads.resize(6);
    const auto six = site_alignment("ACG-------GTA", "ACGTTTTTT-GTA");
    const auto seven = site_alignment("ACG-------GTA", "ACGTTTTTTTGTA");
    const auto ref = site_alignment("ACG-------GTA", "ACG-------GTA");
    const auto eight = site_alignment("ACG--------GTA", "ACGTTTTTTTTGTA");
    const auto nine = site_alignment("ACG---------GTA", "ACGTTTTTTTTTGTA");
    const std::array<std::vector<AlnStr>, 2> alignments = {{{six}, {seven}}};
    std::vector<CandidateVariant> vars;
    std::vector<VariantCategory> categories;
    std::vector<ReadVariantProfile> profiles;
    Options opts; opts.recover_gaps = opts.private_msa_admit_all_in_region = true;
    make_vars_from_msa_cons_aln(opts, chunk, 6, {}, 100, 2, {0, 0}, {},
                               alignments, vars, categories, profiles);
    bool ok = check(vars.size() == 1 && vars[0].msa_insertion_alts ==
                    std::vector<std::string>({"TTTTTT", "TTTTTTT"}),
                    "two alternate insertions form one three-allele site");
    if (!ok) return false;
    const std::array<AlnStr, 2> consensuses = {six, seven};
    add_msa_site_observations(opts,
        {{0, {six, six}}, {1, {seven, seven}}, {2, {ref, ref}},
         {3, {nine, nine}}, {4, {six, seven}}, {5, {eight, eight}}},
        100, false, vars, profiles, &consensuses);
    ok &= check(vars[0].counts.alle_covs == std::vector<int>({1, 1, 2}) &&
                vars[0].counts.ref_cov == 1 && vars[0].counts.alt_cov == 3 &&
                profiles[0].alleles == std::vector<int>({1}) &&
                profiles[1].alleles == std::vector<int>({2}) &&
                profiles[2].alleles == std::vector<int>({0}) &&
                profiles[3].start_var_idx == -1 && profiles[4].start_var_idx == -1 &&
                profiles[5].alleles == std::vector<int>({2}),
                "insertion lengths retain distinct indices and unsupported lengths abstain");
    chunk.candidates = vars;
    auto& candidate = chunk.candidates[0];
    candidate.gap_link_supported = true;
    candidate.counts.alle_covs = {0, 3, 3};
    candidate.counts.total_cov = 6;
    candidate.hap_to_cons_alle = {-1, 2, 2};
    opts.link_by_alleles = opts.private_msa_admit_all_in_region = true;
    iter_update_var_hap_cons_phase_set(chunk, {0}, opts);
    ok &= check(candidate.hap_to_cons_alle[1] != candidate.hap_to_cons_alle[2] &&
                candidate.hap_to_cons_alle[1] > 0 && candidate.hap_to_cons_alle[2] > 0,
                "provisional labels cannot collapse a supported two-ALT genotype");
    PhasingChunk phased;
    phased.candidates = chunk.candidates;
    phased.reads.resize(6); phased.read_var_profile.resize(6);
    phased.read_var_cr.reset(cr_init());
    for (int ri = 0; ri < 6; ++ri) {
        auto& profile = phased.read_var_profile[ri];
        profile.start_var_idx = profile.end_var_idx = 0;
        profile.alleles = {ri < 3 ? 1 : 2}; profile.alt_qi = {-1};
        cr_add(phased.read_var_cr.get(), "cr", 0, 1, ri);
    }
    cr_index(phased.read_var_cr.get());
    assign_hap_based_on_germline_het_vars_kmeans(phased, opts, kCandGermlineVarCate);
    ok &= check(phased.candidates[0].hap_to_cons_alle[1] != phased.candidates[0].hap_to_cons_alle[2] &&
                phased.haps[0] != 0 && phased.haps[3] != 0 && phased.haps[0] != phased.haps[3],
                "the full k-means pass keeps the two insertion lengths on opposite haplotypes");
    char path[] = "/tmp/pgphase-multi-vcf-XXXXXX";
    const int fd = mkstemp(path);
    if (fd < 0) return check(false, "create multiallelic VCF reference fixture");
    FILE* fasta = fdopen(fd, "w");
    std::fputs(">chr1\nACGGTA\n", fasta);
    std::fclose(fasta);
    if (fai_build(path) != 0) return check(false, "index multiallelic VCF reference fixture");
    {
        std::unique_ptr<faidx_t, decltype(&fai_destroy)> fai(fai_load(path), fai_destroy);
        const std::string header_text = "@SQ\tSN:chr1\tLN:6\n";
        std::unique_ptr<bam_hdr_t, decltype(&bam_hdr_destroy)> header(
            sam_hdr_parse(header_text.size(), header_text.c_str()), bam_hdr_destroy);
        ReferenceCache reference(fai.get());
        candidate.key.pos = 4;
        candidate.phase_set = 3;
        candidate.counts.ref_cov = 0; candidate.counts.alt_cov = 6;
        opts.min_depth = opts.min_alt_depth = 1;
        std::ostringstream output;
        write_phased_variants_vcf_records(output, opts, header.get(), reference, chunk.candidates);
        ok &= check(output.str().find("\tG\tGTTTTTT,GTTTTTTT\t") != std::string::npos &&
                    output.str().find("\t1|2:6:0,3,3:0.5,0.5:0:3") != std::string::npos &&
                    output.str().find(";AF=0.5,0.5;") != std::string::npos,
                    "VCF preserves both ALT alleles, 1|2 genotype, allele depths and PS");
    }
    std::remove(path);
    std::remove((std::string(path) + ".fai").c_str());
    return ok;
}

// An MSA-verified private het SNP inside a gap earns a bridge verdict only when
// --gap-bridge-private-snps admits it; the default path leaves it unsupported.
// Reads between --min-mapq and --min-assign-mapq supply allele evidence but must
// not carry a haplotype tag; equal floors must behave exactly as one floor.
static bool test_assign_mapq_floor_gates_tags_only() {
    bool ok = true;
    Options opts;
    ok &= check(opts.min_assign_mapq == opts.min_mapq,
                "assignment floor defaults to the discovery floor");
    ok &= check(read_carries_phase_tags(opts.min_mapq, opts),
                "a read at the default floor carries tags");
    ok &= check(!read_carries_phase_tags(opts.min_mapq - 1, opts),
                "a read below the default floor carries no tags");
    opts.min_mapq = 1;
    ok &= check(!read_carries_phase_tags(3, opts),
                "an admitted low-MAPQ read carries no tags once the floors differ");
    ok &= check(read_carries_phase_tags(30, opts),
                "a confidently mapped read still carries tags");
    opts.min_assign_mapq = 5;
    ok &= check(read_carries_phase_tags(5, opts),
                "lowering the assignment floor admits reads at it");
    return ok;
}

static bool test_verified_msa_snp_bridges_without_flag() {
    bool ok = true;
    for (const bool enabled : {false, true}) {
        PhasingChunk chunk;
        CandidateVariant site;
        site.key.pos = 200; site.key.type = VariantType::Snp;
        site.key.alt = "A"; site.key.ref_len = 1;
        site.lcd_var_i_to_cate = kCandNoisyCandHet;
        site.msa_verified = true;
        CandidateVariant anchor;
        anchor.key.pos = 100; anchor.key.type = VariantType::Snp;
        anchor.key.alt = "T"; anchor.key.ref_len = 1;
        anchor.lcd_var_i_to_cate = kCandCleanHetSnp;
        anchor.hap_to_cons_alle = {-1, 0, 1};
        chunk.candidates = {anchor, site};
        chunk.read_var_cr.reset(cr_init());
        const std::array<int, 4> table{6, 0, 0, 6};
        for (int cell = 0; cell < 4; ++cell)
            for (int n = 0; n < table[cell]; ++n) {
                const int ri = static_cast<int>(chunk.reads.size());
                chunk.reads.emplace_back();
                ReadVariantProfile profile;
                profile.start_var_idx = 0; profile.end_var_idx = 1;
                profile.alleles = {cell / 2, cell % 2};
                profile.alt_qi = {-1, -1};
                chunk.read_var_profile.push_back(profile);
                cr_add(chunk.read_var_cr.get(), "cr", 0, 2, ri);
            }
        cr_index(chunk.read_var_cr.get());
        Options opts;
        opts.recover_gaps = opts.link_by_alleles = opts.private_msa_admit_all_in_region = true;
        opts.gap_recovery_beg = 150; opts.gap_recovery_end = 250;
        opts.gap_bridge_private_snps = enabled;
        assign_hap_based_on_germline_het_vars_kmeans(chunk, opts, kCandGermlineVarCate);
        ok &= check(chunk.candidates[1].gap_link_supported,
                    "MSA-verified het SNP earns a bridge verdict with the flag off too");
    }
    return ok;
}

static bool test_hp_gap_site_scores_reads_only_in_its_gap() {
    bool ok = true;
    // A homopolymer indel is excluded from read scores everywhere except the gap
    // the homopolymer tier admitted it for. Without that exception the tier could
    // earn link support and still leave the gap's reads unphased, because the
    // reads inside such a gap frequently have no other interior evidence.
    for (const bool hp_tier : {false, true}) {
        PhasingChunk chunk;
        CandidateVariant anchor;
        anchor.key.pos = 100; anchor.key.type = VariantType::Snp;
        anchor.key.alt = "T"; anchor.key.ref_len = 1;
        anchor.lcd_var_i_to_cate = kCandCleanHetSnp;
        anchor.hap_to_cons_alle = {-1, 0, 1};
        CandidateVariant hp_site;
        hp_site.key.pos = 200; hp_site.key.type = VariantType::Deletion;
        hp_site.key.alt = "."; hp_site.key.ref_len = 4;
        hp_site.lcd_var_i_to_cate = kCandNoisyCandHet;
        hp_site.msa_verified = true;
        hp_site.is_homopolymer_indel = true;
        chunk.candidates = {anchor, hp_site};
        chunk.read_var_cr.reset(cr_init());
        const std::array<int, 4> table{6, 0, 0, 6};
        for (int cell = 0; cell < 4; ++cell)
            for (int n = 0; n < table[cell]; ++n) {
                const int ri = static_cast<int>(chunk.reads.size());
                chunk.reads.emplace_back();
                ReadVariantProfile profile;
                profile.start_var_idx = 0; profile.end_var_idx = 1;
                profile.alleles = {cell / 2, cell % 2};
                profile.alt_qi = {-1, -1};
                chunk.read_var_profile.push_back(profile);
                cr_add(chunk.read_var_cr.get(), "cr", 0, 2, ri);
            }
        cr_index(chunk.read_var_cr.get());
        Options opts;
        opts.recover_gaps = opts.link_by_alleles = opts.private_msa_admit_all_in_region = true;
        opts.gap_recovery_beg = 150; opts.gap_recovery_end = 250;
        if (hp_tier) { opts.gap_hp_link_beg = 150; opts.gap_hp_link_end = 250; }
        assign_hap_based_on_germline_het_vars_kmeans(chunk, opts, kCandGermlineVarCate);
        ok &= check(chunk.candidates[1].hp_gap_scorable == hp_tier,
                    hp_tier ? "verified homopolymer indel is scorable inside the homopolymer tier's gap"
                            : "homopolymer indel stays unscorable outside the homopolymer tier");
    }
    return ok;
}

static bool test_msa_insertion_pair_clean_anchor() {
    bool ok = true;
    // The first pattern separates the two alleles despite one weak row.
    // The second is a homozygous insertion split by MSA length errors.
    for (const auto& votes : {std::array<int, 4>{9, 8, 7, 15},
                              std::array<int, 4>{6, 10, 8, 8}}) {
        PhasingChunk chunk;
        CandidateVariant anchor;
        anchor.key.pos = 100; anchor.key.type = VariantType::Snp;
        anchor.key.alt = "T"; anchor.key.ref_len = 1;
        anchor.lcd_var_i_to_cate = kCandCleanHetSnp;
        anchor.hap_to_cons_alle = {-1, 0, 1};
        anchor.counts.alle_covs = {votes[0] + votes[1], votes[2] + votes[3]};
        anchor.counts.total_cov = std::accumulate(votes.begin(), votes.end(), 0);
        CandidateVariant insertion;
        insertion.key.pos = 200; insertion.key.type = VariantType::Insertion;
        insertion.key.alt = "TTTTTT";
        insertion.msa_insertion_alts = {"TTTTTT", "TTTTTTT"};
        insertion.lcd_var_i_to_cate = kCandNoisyCandHet;
        insertion.counts.alle_covs = {0, votes[0] + votes[2], votes[1] + votes[3]};
        insertion.counts.total_cov = anchor.counts.total_cov;
        chunk.candidates = {anchor, insertion};
        chunk.read_var_cr.reset(cr_init());
        for (int cell = 0; cell < 4; ++cell) {
            for (int n = 0; n < votes[cell]; ++n) {
                const int ri = chunk.reads.size();
                chunk.reads.emplace_back();
                ReadVariantProfile profile;
                profile.start_var_idx = 0; profile.end_var_idx = 1;
                profile.alleles = {cell / 2, 1 + cell % 2};
                profile.alt_qi = {-1, -1};
                chunk.read_var_profile.push_back(profile);
                cr_add(chunk.read_var_cr.get(), "cr", 0, 2, ri);
            }
        }
        cr_index(chunk.read_var_cr.get());
        Options opts;
        opts.recover_gaps = opts.link_by_alleles = opts.private_msa_admit_all_in_region = true;
        assign_hap_based_on_germline_het_vars_kmeans(chunk, opts, kCandGermlineVarCate);
        ok &= check(chunk.candidates[1].gap_link_supported == (votes[0] == 9),
                    "clean anchor separates a real insertion pair from MSA length noise");
    }
    return ok;
}

static bool test_deletion_reference_with_overlapping_snp() {
    VariantKey key;
    key.pos = 103; key.type = VariantType::Deletion; key.ref_len = 4;
    const auto present = site_alignment("ACGTTCCGTA", "ACGTTCTGTA");
    const auto deleted = site_alignment("ACGTTCCGTA", "ACG----GTA");
    const auto third = site_alignment("ACGTTCCGTA", "ACGTTATGTA");
    const auto partial = site_alignment("ACGTTCCGTA", "ACG--CTGTA");
    const std::array<AlnStr, 2> consensuses = {present, deleted};
    bool ok = check(call_msa_site_allele({present, present}, key, 100, &consensuses) == 0,
                    "a verified SNP inside a deletion footprint does not erase the non-deletion allele");
    ok &= check(call_msa_site_allele({deleted, deleted}, key, 100, &consensuses) == 1 &&
                call_msa_site_allele({third, third}, key, 100, &consensuses) == -1 &&
                call_msa_site_allele({partial, partial}, key, 100, &consensuses) == -1 &&
                call_msa_site_allele({present, deleted}, key, 100, &consensuses) == -1,
                "unsupported substitutions, intermediate deletions and path conflicts remain unknown");
    return ok;
}

static bool test_msa_supported_flank_variant() {
    VariantKey key;
    key.pos = 103; key.type = VariantType::Snp; key.ref_len = 1; key.alt = "T";
    const auto read = site_alignment("ACGAC-GT", "ACGTCAGT");
    const std::array<AlnStr, 2> consensuses = {
        site_alignment("ACGAC-GT", "ACGACAGT"), read};
    bool ok = check(call_msa_site_allele({read, read}, key, 100) == -1,
                    "an unexplained flank insertion is rejected");
    ok &= check(call_msa_site_allele({read, read}, key, 100, &consensuses) == 1,
                "an insertion supported by fixed consensuses does not hide the nearby SNP");
    const auto error = site_alignment("ACGAC-GT", "ACGTCTGT");
    ok &= check(call_msa_site_allele({error, error}, key, 100, &consensuses) == -1,
                "an insertion sequence absent from both consensuses remains rejected");
    return ok;
}

static bool test_msa_snp_backfills_all_bam_reads() {
    PhasingChunk chunk;
    CandidateVariant clean;
    clean.key.pos = 100;
    clean.key.type = VariantType::Snp;
    clean.key.ref_len = 1;
    clean.key.alt = "T";
    clean.ref_base = 0;
    clean.lcd_var_i_to_cate = kCandCleanHetSnp;
    CandidateVariant recovered = clean;
    recovered.key.pos = 200;
    recovered.lcd_var_i_to_cate = kCandNoisyCandHet;
    recovered.msa_verified = true;
    CandidateVariant insertion;
    insertion.key.pos = 150;
    insertion.key.type = VariantType::Insertion;
    insertion.key.alt = "G";
    insertion.lcd_var_i_to_cate = kCandNoisyCandHet;
    insertion.msa_verified = true;
    chunk.candidates = {clean, insertion, recovered};

    auto read = min_read();
    read.mapq = 60;
    read.alignment.reset(bam_init1());
    std::string sequence(201, 'A');
    sequence[100] = 'T';
    const uint32_t cigar = bam_cigar_gen(sequence.size(), BAM_CMATCH);
    bam_set1(read.alignment.get(), 7, "backfill", 0, 0, 99, 60, 1, &cigar,
             -1, -1, 0, sequence.size(), sequence.c_str(), nullptr, 0);
    std::fill(bam_get_qual(read.alignment.get()),
              bam_get_qual(read.alignment.get()) + sequence.size(), 40);
    chunk.reads.push_back(std::move(read));
    ReadVariantProfile profile;
    profile.start_var_idx = 0;
    profile.end_var_idx = 0;
    profile.alleles = {0};
    profile.alt_qi = {0};
    chunk.read_var_profile.push_back(std::move(profile));

    Options opts;
    opts.min_bq = 10;
    const int added = backfill_msa_observations(chunk, opts, 150, 250);
    return check(added == 2 && chunk.read_var_profile[0].end_var_idx == 2 &&
                     chunk.read_var_profile[0].alleles[1] == 0 &&
                     chunk.read_var_profile[0].alleles[2] == 1,
                 "MSA SNP and exact indel observations are backfilled from BAM reads");
}

static bool test_gap_backfill_preserves_homozygous_calls() {
    bool ok = true;
    for (int mode = 0; mode < 2; ++mode) {
        PhasingChunk chunk;
        CandidateVariant var;
        var.key.pos = 200;
        var.key.type = VariantType::Deletion;
        var.key.ref_len = 1;
        var.msa_verified = true;
        var.is_homopolymer_indel = mode == 1;
        var.counts.category = VariantCategory::NoisyCandHom;
        var.lcd_var_i_to_cate = kCandNoisyCandHom;
        chunk.candidates.push_back(var);
        for (int ri = 0; ri < 4; ++ri) {
            auto read = min_read();
            read.qname = std::to_string(ri);
            read.mapq = 60;
            read.alignment.reset(bam_init1());
            const bool deletion = ri >= 2;
            const std::string sequence(deletion ? 200 : 201, 'A');
            const std::vector<uint32_t> cigar = deletion
                ? std::vector<uint32_t>{bam_cigar_gen(100, BAM_CMATCH),
                    bam_cigar_gen(1, BAM_CDEL), bam_cigar_gen(100, BAM_CMATCH)}
                : std::vector<uint32_t>{bam_cigar_gen(201, BAM_CMATCH)};
            bam_set1(read.alignment.get(), 7, "backfill", 0, 0, 99, 60,
                     cigar.size(), cigar.data(), -1, -1, 0,
                     sequence.size(), sequence.c_str(), nullptr, 0);
            std::fill(bam_get_qual(read.alignment.get()),
                      bam_get_qual(read.alignment.get()) + sequence.size(), 40);
            chunk.reads.push_back(std::move(read));
            chunk.read_var_profile.emplace_back();
        }
        Options opts;
        opts.min_alt_depth = 2;
        ok &= check(backfill_msa_observations(chunk, opts, 150, 250) == 0,
                    "ordinary backfill leaves homozygous candidates alone");
        ok &= check(chunk.candidates[0].lcd_var_i_to_cate == kCandNoisyCandHom &&
                    chunk.candidates[0].is_homopolymer_indel == (mode == 1),
                    "confirmation cannot promote a sparse genotype or remove its repeat flag");
    }
    return ok;
}

static bool test_gap_bridge_uses_equivalent_shifted_msa_insertion(int mode = 0) {
    PhasingChunk chunk;
    chunk.ref_beg = 1;
    chunk.ref_end = 400;
    chunk.ref_seq.assign(400, 'A');
    chunk.ref_seq.replace(199, 5, "TCTCT");
    for (int vi = 0; vi < 3; ++vi) {
        CandidateVariant var;
        var.key.pos = 100 + vi * 100;
        var.key.type = VariantType::Snp;
        var.key.ref_len = 1;
        var.key.alt = "T";
        var.ref_base = 0;
        var.lcd_var_i_to_cate = kCandCleanHetSnp;
        var.hap_to_cons_alle = {-1, 0, 1};
        chunk.candidates.push_back(var);
    }
    auto& insertion = chunk.candidates[1];
    insertion.key.type = VariantType::Insertion;
    insertion.key.ref_len = 0;
    insertion.key.alt = "TC";
    insertion.msa_insertion_alts = {"TC", "TCTC"};
    insertion.lcd_var_i_to_cate = kCandNoisyCandHet;
    insertion.msa_verified = true;
    insertion.gap_link_supported = true;
    insertion.hap_to_cons_alle = {-1, 1, 2};
    if (mode == 2) {
        insertion.key.alt = "T";
        insertion.msa_insertion_alts = {"T", "TT"};
    }

    std::vector<std::vector<int>> alleles = {
        {0, 1, -1}, {1, 2, -1}, {-1, -1, 0}, {-1, -1, 1},
        {-1, 1, 1}, {-1, 1, 1}};
    if (mode > 0) alleles.pop_back();
    chunk.read_var_cr.reset(cr_init());
    for (size_t ri = 0; ri < alleles.size(); ++ri) {
        auto read = min_read();
        read.mapq = 60;
        if (ri >= 4) {
            read.alignment.reset(bam_init1());
            std::string sequence(mode == 2 ? 103 : 104, 'A');
            if (mode == 2) sequence[1] = 'T';
            else { sequence[6] = 'C'; sequence[7] = 'T'; }
            sequence.back() = 'T';
            const std::array<uint32_t, 3> cigar = {
                bam_cigar_gen(mode == 2 ? 1 : 6, BAM_CMATCH),
                bam_cigar_gen(mode == 2 ? 1 : 2, BAM_CINS),
                bam_cigar_gen(mode == 2 ? 101 : 96, BAM_CMATCH)};
            bam_set1(read.alignment.get(), 7, "msaedge", 0, 0, 198, 60,
                     cigar.size(), cigar.data(),
                     -1, -1, 0, sequence.size(), sequence.c_str(), nullptr, 0);
            std::fill(bam_get_qual(read.alignment.get()),
                      bam_get_qual(read.alignment.get()) + sequence.size(), mode == 3 ? 10 : 40);
            if (mode >= 4) {
                auto* quality = bam_get_qual(read.alignment.get());
                quality[5] = 22;
                quality[7] = mode == 6 ? 19 : mode == 7 ? 255 : 22;
                if (mode == 5) quality[6] = 22;
            }
        }
        chunk.reads.push_back(std::move(read));
        ReadVariantProfile profile;
        profile.start_var_idx = 0;
        profile.end_var_idx = 2;
        profile.alleles = alleles[ri];
        profile.alt_qi.assign(3, -1);
        chunk.read_var_profile.push_back(std::move(profile));
        cr_add(chunk.read_var_cr.get(), "cr", 0, 3, ri);
    }
    cr_index(chunk.read_var_cr.get());
    Options opts;
    opts.link_by_alleles = opts.recover_gaps = opts.private_msa_admit_all_in_region = true;
    opts.min_block_link_reads = 2;
    opts.block_link_window = 8;
    iter_update_var_hap_cons_phase_set(chunk, {0, 1, 2}, opts);
    const bool joined = chunk.candidates[0].phase_set == chunk.candidates[2].phase_set;
    return check(joined == (mode < 2 || mode == 4),
                 "singleton MSA insertion bridges require separated alleles and confident BAM support");
}

static bool test_gap_bridge_validates_msa_deletion() {
    bool ok = true;
    for (int mode = 0; mode < 7; ++mode) {
        PhasingChunk chunk;
        for (int vi = 0; vi < 3; ++vi) {
            CandidateVariant var;
            var.key.pos = 100 + vi * 100;
            var.key.type = VariantType::Snp;
            var.key.ref_len = 1;
            var.key.alt = "T";
            var.ref_base = 0;
            var.lcd_var_i_to_cate = kCandCleanHetSnp;
            var.hap_to_cons_alle = {-1, 0, 1};
            chunk.candidates.push_back(var);
        }
        auto& deletion = chunk.candidates[1];
        deletion.key.type = VariantType::Deletion;
        deletion.key.ref_len = 2;
        deletion.key.alt.clear();
        deletion.msa_verified = deletion.gap_link_supported = true;
        deletion.is_homopolymer_indel = mode == 6;
        deletion.lcd_var_i_to_cate = kCandNoisyCandHet;
        const int bridge_allele = mode == 1 || mode == 3 ? 1 : 0;
        std::vector<std::vector<int>> alleles = {
            {0, 0, -1}, {1, 1, -1}, {-1, -1, 0}, {-1, -1, 1},
            {-1, bridge_allele, 1}};
        if (mode == 5) alleles.push_back({-1, 0, 0});
        chunk.read_var_cr.reset(cr_init());
        for (size_t ri = 0; ri < alleles.size(); ++ri) {
            auto read = min_read();
            read.mapq = 60;
            if (ri >= 4) {
                read.alignment.reset(bam_init1());
                const bool deleted = mode == 1 || mode == 3 || mode == 4;
                std::string sequence(deleted ? 100 : 102, 'A');
                sequence.back() = alleles[ri][2] ? 'T' : 'A';
                std::vector<uint32_t> cigar;
                if (deleted) cigar = {bam_cigar_gen(1, BAM_CMATCH),
                    bam_cigar_gen(2, mode == 3 ? BAM_CREF_SKIP : BAM_CDEL),
                    bam_cigar_gen(99, BAM_CMATCH)};
                else cigar = {bam_cigar_gen(102, BAM_CMATCH)};
                bam_set1(read.alignment.get(), 6, "bridge", 0, 0, 198, 60,
                         cigar.size(), cigar.data(), -1, -1, 0,
                         sequence.size(), sequence.c_str(), nullptr, 0);
                std::fill(bam_get_qual(read.alignment.get()),
                          bam_get_qual(read.alignment.get()) + sequence.size(), mode == 2 ? 10 : 40);
            }
            chunk.reads.push_back(std::move(read));
            ReadVariantProfile profile;
            profile.start_var_idx = 0;
            profile.end_var_idx = 2;
            profile.alleles = alleles[ri];
            chunk.read_var_profile.push_back(std::move(profile));
            cr_add(chunk.read_var_cr.get(), "cr", 0, 3, ri);
        }
        cr_index(chunk.read_var_cr.get());
        Options opts;
        opts.link_by_alleles = opts.recover_gaps = opts.private_msa_admit_all_in_region = true;
        opts.min_block_link_reads = 2;
        iter_update_var_hap_cons_phase_set(chunk, {0, 1, 2}, opts);
        const bool joined = chunk.candidates[0].phase_set == chunk.candidates[2].phase_set;
        ok &= check(joined == (mode < 2), "MSA deletion bridge requires the exact high-quality BAM allele");
        if (joined)
            ok &= check((chunk.candidates[0].hap_to_cons_alle[1] !=
                         chunk.candidates[2].hap_to_cons_alle[1]) == (bridge_allele == 0),
                        "MSA deletion bridge preserves REF/ALT orientation");
    }
    return ok;
}

static bool test_orphan_msa_site_uses_stitched_read_orientation() {
    bool ok = true;
    for (int mode = 0; mode < 4; ++mode) {
        std::vector<PhasingChunk> chunks(1);
        auto& chunk = chunks.front();
        chunk.region.tid = 0;
        CandidateVariant var;
        var.key.pos = 200;
        var.key.type = VariantType::Insertion;
        var.key.alt = "C";
        var.msa_insertion_alts = {"C", "CC"};
        var.msa_verified = mode != 3;
        var.hap_to_cons_alle = {-1, 1, 2};
        var.phase_set = mode == 2 ? 100 : 200;
        chunk.candidates.push_back(var);
        chunk.read_var_cr.reset(cr_init());
        for (int ri = 0; ri < 4; ++ri) {
            auto read = min_read();
            read.qname = std::to_string(ri);
            chunk.reads.push_back(std::move(read));
            chunk.haps.push_back(ri < 2 ? 1 : 2);
            chunk.phase_sets.push_back(100);
            ReadVariantProfile profile;
            profile.start_var_idx = profile.end_var_idx = 0;
            profile.alleles = {ri < 2 || mode == 1 ? 2 : 1};
            chunk.read_var_profile.push_back(std::move(profile));
            cr_add(chunk.read_var_cr.get(), "cr", 0, 1, ri);
        }
        cr_index(chunk.read_var_cr.get());
        Options opts;
        opts.min_block_link_reads = 2;
        const int attached = anchor_orphan_msa_sites(chunks, opts);
        ok &= check(attached == (mode == 0 ? 1 : 0),
                    "only verified orphan sites with support on both haplotypes are attached");
        if (mode == 0)
            ok &= check(chunk.candidates[0].phase_set == 100 &&
                        chunk.candidates[0].hap_to_cons_alle[1] == 2 &&
                        chunk.candidates[0].hap_to_cons_alle[2] == 1,
                        "alternate/alternate site orientation follows the stitched reads");
        ok &= check(chunk.haps == std::vector<int>({1, 1, 2, 2}) &&
                    chunk.phase_sets == std::vector<hts_pos_t>({100, 100, 100, 100}),
                    "orphan attachment does not change read assignments");
    }
    return ok;
}

static bool test_gap_bridge_separates_indel_boundary_signal() {
    bool ok = true;
    for (int mode = 0; mode < 4; ++mode) {
        PhasingChunk chunk;
        const std::array<int, 5> positions{100, 200, 300, mode == 1 ? 350 : 301, 400};
        for (const auto pos : positions) {
            CandidateVariant var;
            var.key.pos = pos;
            var.key.type = VariantType::Snp;
            var.key.ref_len = 1;
            var.key.alt = "T";
            var.ref_base = 0;
            var.lcd_var_i_to_cate = kCandCleanHetSnp;
            var.hap_to_cons_alle = {-1, 0, 1};
            chunk.candidates.push_back(var);
        }
        auto& insertion = chunk.candidates[1];
        insertion.key.type = VariantType::Insertion;
        insertion.key.ref_len = 0;
        insertion.msa_verified = insertion.gap_link_supported = true;
        insertion.lcd_var_i_to_cate = kCandNoisyCandHet;
        auto& deletion = chunk.candidates[2];
        deletion.key.type = VariantType::Deletion;
        deletion.key.alt.clear();
        deletion.msa_verified = mode != 2;
        deletion.gap_link_supported = true;
        deletion.lcd_var_i_to_cate = kCandNoisyCandHet;
        const std::vector<std::vector<int>> alleles = {
            {0, 0, -1, -1, -1}, {1, 1, -1, -1, -1},
            {-1, -1, 0, 0, 0}, {-1, -1, 1, 1, 1},
            {-1, 0, 0, 0, 1}};
        chunk.read_var_cr.reset(cr_init());
        for (size_t ri = 0; ri < alleles.size(); ++ri) {
            auto read = min_read();
            read.mapq = 60;
            if (ri == 4) {
                read.alignment.reset(bam_init1());
                std::string sequence(202, 'A');
                sequence.back() = 'T';
                const uint32_t cigar = bam_cigar_gen(sequence.size(), BAM_CMATCH);
                bam_set1(read.alignment.get(), 6, "bridge", 0, 0, 198, 60, 1, &cigar,
                         -1, -1, 0, sequence.size(), sequence.c_str(), nullptr, 0);
                std::fill(bam_get_qual(read.alignment.get()),
                          bam_get_qual(read.alignment.get()) + sequence.size(), 40);
                if (mode == 3) bam_get_qual(read.alignment.get())[sequence.size() - 1] = 10;
            }
            chunk.reads.push_back(std::move(read));
            ReadVariantProfile profile;
            profile.start_var_idx = 0;
            profile.end_var_idx = 4;
            profile.alleles = alleles[ri];
            chunk.read_var_profile.push_back(std::move(profile));
            cr_add(chunk.read_var_cr.get(), "cr", 0, 5, ri);
        }
        cr_index(chunk.read_var_cr.get());
        Options opts;
        opts.link_by_alleles = opts.recover_gaps = opts.private_msa_admit_all_in_region = true;
        opts.min_block_link_reads = 2;
        iter_update_var_hap_cons_phase_set(chunk, {0, 1, 2, 3, 4}, opts);
        const bool joined = chunk.candidates[0].phase_set == chunk.candidates[4].phase_set;
        ok &= check(joined == (mode == 0),
                    "indel-boundary signal cannot veto a separate confident clean SNP bridge");
        if (joined)
            ok &= check(chunk.candidates[0].hap_to_cons_alle[1] != chunk.candidates[4].hap_to_cons_alle[1],
                        "indel-boundary rescue uses the independent SNP orientation");
    }
    return ok;
}

static bool test_gap_clean_block_bridge() {
    bool ok = true;
    for (int mode = 0; mode < 11; ++mode) {
        PhasingChunk chunk;
        for (int vi = 0; vi < 4; ++vi) {
            CandidateVariant var;
            var.key.pos = 100 + vi * 100;
            var.key.type = VariantType::Snp;
            var.key.ref_len = 1; var.key.alt = "T"; var.ref_base = 0;
            var.lcd_var_i_to_cate = kCandCleanHetSnp;
            var.counts.alle_covs = {3, 3};
            var.hap_to_cons_alle = {-1, 0, 1};
            chunk.candidates.push_back(var);
        }
        std::vector<std::vector<int>> alleles = {
            {0, 0, -1, -1}, {1, 1, -1, -1},
            {-1, -1, 0, 0}, {-1, -1, 1, 1}, {0, 0, 1, 1}};
        if (mode == 1 || mode >= 4) alleles.back() = {0, -1, 1, -1};
        if (mode == 6 || mode == 7) alleles.back() = {0, -1, 1, 1};
        if (mode == 2) alleles.push_back({0, -1, 0, -1});
        chunk.read_var_cr.reset(cr_init());
        for (size_t ri = 0; ri < alleles.size(); ++ri) {
            auto read = min_read();
            read.mapq = mode == 3 && ri == 4 ? 0 : 60;
            if (mode >= 4 && ri == 4) {
                read.alignment.reset(bam_init1());
                std::string sequence(301, 'A');
                sequence[200] = sequence[300] = 'T';
                const uint32_t cigar = bam_cigar_gen(sequence.size(), BAM_CMATCH);
                bam_set1(read.alignment.get(), 6, "bridge", 0, 0, 99, 60, 1, &cigar,
                         -1, -1, 0, sequence.size(), sequence.c_str(), nullptr, 0);
                const int quality = mode == 4 ? 40 : mode == 6 ? 22 : 10;
                std::fill(bam_get_qual(read.alignment.get()),
                          bam_get_qual(read.alignment.get()) + sequence.size(), quality);
                if (mode == 10) bam_get_qual(read.alignment.get())[200] = 40;
            }
            chunk.reads.push_back(std::move(read));
            ReadVariantProfile profile;
            profile.start_var_idx = 0; profile.end_var_idx = 3;
            profile.alleles = alleles[ri];
            profile.alt_qi.assign(4, -1);
            if (ri == 4 && (mode == 8 || mode == 9 || mode == 10)) {
                profile.alt_qi[0] = kGraphConfirmedAltQi;
                if (mode == 8) profile.alt_qi[2] = kGraphConfirmedAltQi;
            }
            chunk.read_var_profile.push_back(profile);
            cr_add(chunk.read_var_cr.get(), "cr", 0, 4, ri);
        }
        cr_index(chunk.read_var_cr.get());
        Options opts;
        opts.link_by_alleles = opts.recover_gaps = opts.private_msa_admit_all_in_region = true;
        opts.min_block_link_reads = 2; opts.block_link_window = 8;
        iter_update_var_hap_cons_phase_set(chunk, {0, 1, 2, 3}, opts);
        const bool joined = chunk.candidates[0].phase_set == chunk.candidates[3].phase_set;
        const std::string bridge_message =
            "singleton block bridge requires strong clean SNPs or graph confirmation on both flanks (mode " +
            std::to_string(mode) + ")";
        ok &= check(joined == (mode == 0 || mode == 4 || mode == 8 || mode == 10),
                    bridge_message.c_str());
        if (joined) {
            ok &= check(chunk.candidates[0].hap_to_cons_alle[1] != chunk.candidates[3].hap_to_cons_alle[1],
                        "block bridge composes the supported opposite orientation");
            ok &= check(iter_update_var_hap_cons_phase_set(chunk, {0, 1, 2, 3}, opts) == 0,
                        "block bridge orientation converges");
        }
    }
    return ok;
}

static bool test_read_hp_matches_reported_phase_set() {
    PhasingChunk chunk;
    for (int vi = 0; vi < 4; ++vi) {
        CandidateVariant var;
        var.key.pos = 100 + vi * 100;
        var.key.type = VariantType::Snp; var.key.ref_len = 1; var.key.alt = "T";
        var.lcd_var_i_to_cate = kCandCleanHetSnp;
        var.counts.alle_covs = {4, 4}; var.counts.total_cov = 8;
        chunk.candidates.push_back(var);
    }
    const std::vector<std::vector<int>> alleles = {
        {0, 0, -1, -1}, {0, 0, -1, -1}, {1, 1, -1, -1}, {1, 1, -1, -1},
        {-1, -1, 0, 0}, {-1, -1, 0, 0}, {-1, -1, 1, 1}, {-1, -1, 1, 1},
        {0, 0, 0, 0}, {0, 0, 1, 1}};
    chunk.read_var_cr.reset(cr_init());
    for (size_t ri = 0; ri < alleles.size(); ++ri) {
        chunk.reads.push_back(min_read());
        ReadVariantProfile profile;
        profile.start_var_idx = 0; profile.end_var_idx = 3;
        profile.alleles = alleles[ri];
        chunk.read_var_profile.push_back(profile);
        cr_add(chunk.read_var_cr.get(), "cr", 0, 4, ri);
    }
    cr_index(chunk.read_var_cr.get());
    Options opts; opts.link_by_alleles = true; opts.min_block_link_reads = 2;
    assign_hap_based_on_germline_het_vars_kmeans(chunk, opts, kCandGermlineClean);
    return check(chunk.candidates[0].phase_set != chunk.candidates[3].phase_set &&
                 chunk.phase_sets[8] == chunk.candidates[0].phase_set &&
                 chunk.phase_sets[9] == chunk.candidates[0].phase_set &&
                 chunk.haps[8] != 0 && chunk.haps[8] == chunk.haps[9] &&
                 chunk.candidates[0].hap_to_cons_alle[chunk.haps[8]] == 0 &&
                 chunk.candidates[0].hap_to_alle_profile[chunk.haps[8]][0] == 4 &&
                 chunk.candidates[0].hap_to_alle_profile[3 - chunk.haps[8]][0] == 0,
                 "spanning reads update and report the left-block hap independently of right-block alleles");
}

static bool test_graph_gap_uses_bam_nucleotide() {
    PhasingChunk chunk;
    CandidateVariant site;
    site.key.pos = 100;
    site.key.type = VariantType::Snp;
    site.key.alt = "T";
    site.ref_base = 0;
    chunk.candidates.push_back(site);
    for (int i = 0; i < 2; ++i) {
        auto read = min_read();
        read.beg = 100; read.end = 110; read.mapq = 60;
        read.alignment.reset(bam_init1());
        const std::string sequence(11, 'A');
        const uint32_t cigar = bam_cigar_gen(11, BAM_CMATCH);
        bam_set1(read.alignment.get(), 4, "test", 0, 0, 99, 60, 1, &cigar,
                 -1, -1, 0, sequence.size(), sequence.c_str(), nullptr, 0);
        std::fill(bam_get_qual(read.alignment.get()), bam_get_qual(read.alignment.get()) + 11, 40);
        chunk.reads.push_back(std::move(read));
        ReadVariantProfile profile;
        profile.start_var_idx = profile.end_var_idx = 0;
        profile.alleles = {1};
        profile.alt_qi = {kGraphConfirmedAltQi};
        profile.graph_alleles = {i == 0 ? 1 : -1};
        chunk.read_var_profile.push_back(profile);
    }
    Options opts;
    const int selected = select_graph_gap_bam_reads(chunk, {0, 1, 2, 100, 105, 1, 200}, opts);
    return check(selected == 1 && chunk.read_var_profile[0].alleles[0] == 0 &&
                 chunk.read_var_profile[0].graph_alleles[0] == 1 && chunk.reads[1].is_skipped,
                 "graph selects the read, but BAM supplies its nucleotide");
}

static PhasingChunk gap_evidence_fixture(bool permute = false) {
    PhasingChunk chunk;
    chunk.region.tid = 0;
    chunk.region.beg = chunk.ref_beg = 50;
    chunk.region.end = chunk.ref_end = 350;
    chunk.ref_seq.assign(301, 'A');
    for (const int pos : {80, 100, 150, 200, 220, 290, 300}) {
        auto var = dummy_cand(pos);
        var.key = {0, pos, VariantType::Snp, 1, "T"};
        var.ref_base = 0;
        var.lcd_var_i_to_cate = kCandCleanHetSnp;
        var.counts.total_cov = 999;
        if (pos == 100 || pos == 300) {
            var.graph_site = true;
            var.phase_set = pos;
            var.hap_to_cons_alle = {-1, 0, 1};
        }
        if (pos == 200) {
            var.key.type = VariantType::Insertion;
            var.key.ref_len = 0;
            var.msa_verified = true;
            var.lcd_var_i_to_cate = kCandNoisyCandHet;
            var.msa_insertion_alts = permute ? std::vector<std::string>{"T", "GT"} : std::vector<std::string>{"GT", "T"};
        }
        if (pos == 220) {
            var.key.type = VariantType::Insertion;
            var.key.alt = "CG"; // a replacement, not a pure insertion
        }
        if (pos == 290) {
            var.key.type = VariantType::Deletion;
            var.key.ref_len = 20;
            var.key.alt.clear();
        }
        chunk.candidates.push_back(var);
    }
    for (int i = 0; i < 2; ++i) {
        auto read = min_read();
        read.qname = "molecule" + std::to_string(i);
        read.tid = 0;
        read.mapq = 60;
        read.beg = 50; read.end = 350;
        read.alignment.reset(bam_init1());
        const std::string sequence(301, 'A');
        const uint32_t cigar = bam_cigar_gen(301, BAM_CMATCH);
        bam_set1(read.alignment.get(), read.qname.size(), read.qname.c_str(), 0, 0, 49, 60,
                 1, &cigar, -1, -1, 0, sequence.size(), sequence.c_str(), nullptr, 0);
        std::fill(bam_get_qual(read.alignment.get()), bam_get_qual(read.alignment.get()) + 301, 40);
        chunk.reads.push_back(std::move(read));
        ReadVariantProfile profile;
        profile.start_var_idx = 0; profile.end_var_idx = 6;
        profile.alleles = {1, 0, 1, permute ? 2 - i : 1 + i, 1, 1, 1};
        profile.alt_qi.assign(7, kGraphConfirmedAltQi);
        profile.bam_alleles = {1, 1, 1, 1, -1, -1, 0};
        profile.bam_qi.assign(7, 50);
        profile.graph_alleles = {-1, 0, -1, -1, -1, -1, 1};
        chunk.read_var_profile.push_back(profile);
    }
    return chunk;
}

static bool test_gap_evidence_scope_and_sources() {
    const PhaseGap gap{0, 100, 300, 100, 300, 50, 350};
    const GapEvidence evidence(gap_evidence_fixture(), gap);
    const GapEvidence permuted(gap_evidence_fixture(true), gap);
    auto projected = evidence.project();
    bool ok = check(projected.candidates.size() == 4, "only gap-owned private events and graph anchors enter projection");
    ok &= check(evidence.events().size() == 6 && evidence.events()[3].role == GapEventRole::Unsupported &&
                evidence.events()[4].role == GapEventRole::Boundary,
                "complex and boundary events retain explicit non-projectable records");
    ok &= check(projected.candidates[0].key.pos == 100 && projected.candidates.back().key.pos == 300,
                "flank extraction context does not admit the outside-gap private SNP");
    ok &= check(evidence.events()[2].alleles == std::vector<std::string>({"", "T", "GT"}) &&
                evidence.events()[2].id == permuted.events()[2].id,
                "canonical event identity is invariant to MSA allele dictionary order");
    ok &= check(projected.read_var_profile[0].alleles == std::vector<int>({0, 1, 2, 1}) &&
                projected.read_var_profile[0].bam_alleles == std::vector<int>({1, 1, 1, 0}) &&
                projected.read_var_profile[0].alleles == permuted.project().read_var_profile[0].alleles,
                "source-local allele IDs map by sequence and disagreements survive projection");
    ok &= check(evidence.observations()[0].bam.query_index == 50 &&
                evidence.observations()[0].bam.base_quality == 40 &&
                evidence.observations()[0].graph.index == 0,
                "graph confirmation preserves original BAM coordinate and quality");
    ok &= check(projected.candidates[2].counts.total_cov == 2 &&
                projected.candidates[2].counts.alle_covs == std::vector<int>({0, 1, 1}),
                "counts derive from projected molecules, including genotype 1/2");
    projected.read_var_profile[0].alleles[0] = 99;
    projected.candidates[0].phase_set = 999;
    const auto replay = evidence.project();
    ok &= check(replay.read_var_profile[0].alleles[0] == 0 && replay.candidates[0].phase_set == 100 &&
                replay.haps == std::vector<int>({0, 0}), "rephasing cannot modify immutable gap evidence");
    const auto bam_view = evidence.project(Options{}, true);
    ok &= check(bam_view.read_var_profile[0].alleles[0] == 1 &&
                bam_view.read_var_profile[0].graph_alleles[0] == 0,
                "BAM-only views reuse stored observations without erasing graph conflicts");
    const auto second = evidence.project();
    for (size_t i = 0; i < second.candidates.size(); ++i)
        if (!second.candidates[i].graph_site)
            ok &= check(gap_owns_variant(gap, second.candidates[i].key), "private candidate footprint belongs to its gap");
    auto anchored = gap_evidence_fixture();
    anchored.candidates[3].graph_site = true;
    anchored.candidates[3].phase_set = gap.left_ps;
    anchored.candidates[3].hap_to_cons_alle = {-1, 1, 2};
    const auto anchor_view = GapEvidence(std::move(anchored), gap).project();
    ok &= check(anchor_view.candidates[2].hap_to_cons_alle == std::array<int, 3>({-1, 2, 1}),
                "anchor consensus follows the canonical MSA allele dictionary");
    return ok;
}

int main() {
    int failures = 0;
    failures += test_gap_evidence_scope_and_sources() ? 0 : 1;
    failures += test_graph_gap_uses_bam_nucleotide() ? 0 : 1;
    failures += test_gap_edges_compose_before_relabelling() ? 0 : 1;
    failures += test_msa_snp_backfills_all_bam_reads() ? 0 : 1;
    failures += test_gap_backfill_preserves_homozygous_calls() ? 0 : 1;
    failures += test_gap_bridge_uses_equivalent_shifted_msa_insertion() ? 0 : 1;
    for (int mode = 1; mode < 8; ++mode)
        failures += test_gap_bridge_uses_equivalent_shifted_msa_insertion(mode) ? 0 : 1;
    failures += test_read_hp_matches_reported_phase_set() ? 0 : 1;
    failures += test_gap_clean_block_bridge() ? 0 : 1;
    failures += test_gap_bridge_validates_msa_deletion() ? 0 : 1;
    failures += test_gap_bridge_separates_indel_boundary_signal() ? 0 : 1;
    failures += test_orphan_msa_site_uses_stitched_read_orientation() ? 0 : 1;
    failures += test_assign_mapq_floor_gates_tags_only() ? 0 : 1;
    failures += test_verified_msa_snp_bridges_without_flag() ? 0 : 1;
    failures += test_hp_gap_site_scores_reads_only_in_its_gap() ? 0 : 1;
    failures += test_msa_insertion_pair_clean_anchor() ? 0 : 1;
    failures += test_msa_two_alternate_insertions() ? 0 : 1;
    failures += test_deletion_reference_with_overlapping_snp() ? 0 : 1;
    failures += test_msa_edge_requires_both_haplotype_votes() ? 0 : 1;
    failures += test_recovery_rejects_unsupported_biallelic_indel() ? 0 : 1;
    failures += test_repeat_anchor_prefers_complete_evidence() ? 0 : 1;
    failures += test_unassigned_msa_local_allele_evidence() ? 0 : 1;
    failures += test_assigned_repeat_allows_one_local_error() ? 0 : 1;
    failures += test_msa_cluster_membership_is_not_an_allele() ? 0 : 1;
    failures += test_read_phase_set_ignores_non_scoring_repeat() ? 0 : 1;
    failures += test_composed_msa_repeat_placement() ? 0 : 1;
    failures += test_read_phase_set_requires_observed_allele() ? 0 : 1;
    failures += test_msa_counts_do_not_depend_on_recover_gaps() ? 0 : 1;
    failures += test_nested_msa_deletions_share_common_event() ? 0 : 1;
    failures += test_gap_hp_trial_is_transactional() ? 0 : 1;
    failures += test_recovery_links_join_earlier_components() ? 0 : 1;
    failures += test_gap_read_index_freezes_initial_assignments() ? 0 : 1;
    failures += test_msa_supported_flank_variant() ? 0 : 1;
    failures += test_msa_observation_heterozygosity_gate() ? 0 : 1;
    failures += test_msa_site_observations() ? 0 : 1;
    failures += test_gap_recovery_keeps_unphased_observations() ? 0 : 1;
    failures += test_gap_recovery_adopts_unresolved_site_and_read() ? 0 : 1;
    failures += test_gap_recovery_respects_requested_regions() ? 0 : 1;
    failures += test_gap_recovery_prefers_complete_bridge() ? 0 : 1;
    failures += test_gap_recovery_deduplicates_votes() ? 0 : 1;
    failures += test_gap_tiers_are_cumulative() ? 0 : 1;
    failures += test_gap_recovery_join_preserves_blocks() ? 0 : 1;
    failures += test_gap_recovery_partial_extension() ? 0 : 1;
    failures += test_gap_allele_attach_join_only_gates_one_sided() ? 0 : 1;
    failures += test_gap_independent_emits_every_qualifying_group() ? 0 : 1;
    failures += test_gap_msa_covers_clean_stretches() ? 0 : 1;
    failures += test_msa_unsorted_profile_indices() ? 0 : 1;
    failures += test_allele_link_orientation_and_tie() ? 0 : 1;
    failures += test_region_msa_repeat_collision() ? 0 : 1;
    failures += test_cross_chunk_flip_when_haps_disagree() ? 0 : 1;
    failures += test_flip_score_zero_no_merge() ? 0 : 1;
    failures += test_skipped_overlap_read_ignored() ? 0 : 1;
    failures += test_margin_zero_merges_single_vote() ? 0 : 1;
    failures += test_margin_one_abstains_single_vote() ? 0 : 1;
    failures += test_private_msa_merge_admits_only_whitelist() ? 0 : 1;
    if (failures == 0)
        std::printf("ALL PASS\n");
    else
        std::printf("FAILURES: %d\n", failures);
    return failures == 0 ? 0 : 1;
}
