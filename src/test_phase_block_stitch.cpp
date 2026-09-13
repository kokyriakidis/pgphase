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
#include "collect_phase_noisy.hpp"
#include "collect_types.hpp"

#include <cstdio>
#include <vector>

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

static bool test_msa_unsorted_profile_indices() {
    PhasingChunk chunk;
    chunk.reads.push_back(min_read());
    ReadVariantProfile profile;
    profile.start_var_idx = 0;
    profile.end_var_idx = 1;
    profile.alleles = std::vector<int>{1, 0};
    profile.alt_qi = std::vector<int>{17, -1};
    merge_var_profile(chunk,
                      {noisy_merge_cand(200, VariantCategory::NoisyCandHet, 8),
                       noisy_merge_cand(100, VariantCategory::NoisyCandHet, 8)},
                      {VariantCategory::NoisyCandHet, VariantCategory::NoisyCandHet},
                      {profile});
    return check(chunk.candidates[0].key.pos == 100 &&
                 chunk.read_var_profile[0].alleles == std::vector<int>({0, 1}) &&
                 chunk.read_var_profile[0].alt_qi == std::vector<int>({-1, 17}),
                 "sorting MSA calls preserves allele and query-position ownership");
}

int main() {
    int failures = 0;
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
