/**
 * Per-function parity tests: pin the behaviour of functions ported from
 * longcallD against the semantics of their upstream definition.
 *
 * Each TEST_CASE names the upstream function and line it is pinning. A shared
 * name is not evidence the bodies agree -- these cases are that evidence, and
 * they fail if a divergence is reintroduced.
 */
#define CATCH_CONFIG_MAIN
#include "../third_party/catch2/catch.hpp"

#include "collect_phase.hpp"
#include "phasing_types.hpp"

using namespace pgphase_collect;

namespace {

/// A clean het SNP with a known per-haplotype consensus and an empty profile.
CandidateVariant port_cand(hts_pos_t pos, hts_pos_t phase_set, int cons1, int cons2) {
    CandidateVariant v;
    v.key.pos = pos;
    v.key.type = VariantType::Snp;
    v.key.ref_len = 1;
    v.lcd_var_i_to_cate = kCandCleanHetSnp;
    v.phase_set = phase_set;
    v.hap_to_cons_alle[1] = cons1;
    v.hap_to_cons_alle[2] = cons2;
    v.hap_to_alle_profile[1].assign(2, 0);
    v.hap_to_alle_profile[2].assign(2, 0);
    return v;
}

/// One read whose profile covers every candidate, carrying `alleles`.
PhasingChunk two_block_chunk(const std::vector<int>& alleles,
                             hts_pos_t ps_a, hts_pos_t ps_b) {
    PhasingChunk c;
    c.candidates.push_back(port_cand(1000, ps_a, 0, 1));
    c.candidates.push_back(port_cand(2000, ps_b, 0, 1));
    c.reads.resize(1);
    c.haps.assign(1, 0);
    c.phase_sets.assign(1, -1);
    ReadVariantProfile p;
    p.read_id = 0;
    p.start_var_idx = 0;
    p.end_var_idx = 1;
    p.alleles = alleles;
    c.read_var_profile.push_back(p);
    return c;
}

} // namespace

// ── update_var_hap_profile_based_on_read_hap (longcallD assign_hap.c:292) ────

TEST_CASE("one hap per read reaches every variant it covers, across phase sets") {
    // Upstream takes a single hap for the read and applies it to every variant
    // the read covers; it has no notion of a phase set here. Two candidates in
    // DIFFERENT phase sets must therefore both be updated by one call.
    PhasingChunk c = two_block_chunk({1, 1}, 1000, 2000);
    update_var_hap_profile_based_on_read_hap(c, 0, 1, kCandGermlineClean);
    CHECK(c.candidates[0].hap_to_alle_profile[1][1] == 1);
    CHECK(c.candidates[1].hap_to_alle_profile[1][1] == 1);
    CHECK(c.candidates[0].hap_to_alle_profile[2][1] == 0);
    CHECK(c.candidates[1].hap_to_alle_profile[2][1] == 0);
}

TEST_CASE("an unassigned read counts toward both haplotypes") {
    // hap == 0 increments both profiles, matching upstream's handling of a read
    // the sweep has not labelled yet.
    PhasingChunk c = two_block_chunk({0, 1}, 1000, 1000);
    update_var_hap_profile_based_on_read_hap(c, 0, 0, kCandGermlineClean);
    CHECK(c.candidates[0].hap_to_alle_profile[1][0] == 1);
    CHECK(c.candidates[0].hap_to_alle_profile[2][0] == 1);
    CHECK(c.candidates[1].hap_to_alle_profile[1][1] == 1);
    CHECK(c.candidates[1].hap_to_alle_profile[2][1] == 1);
}

TEST_CASE("a category outside the mask is not counted") {
    PhasingChunk c = two_block_chunk({1, 1}, 1000, 1000);
    c.candidates[1].lcd_var_i_to_cate = kCandNoisyCandHet;
    update_var_hap_profile_based_on_read_hap(c, 0, 1, kCandGermlineClean);
    CHECK(c.candidates[0].hap_to_alle_profile[1][1] == 1);
    CHECK(c.candidates[1].hap_to_alle_profile[1][1] == 0);
}

TEST_CASE("an unobserved allele at a covered site is skipped") {
    PhasingChunk c = two_block_chunk({-1, 1}, 1000, 1000);
    update_var_hap_profile_based_on_read_hap(c, 0, 2, kCandGermlineClean);
    CHECK(c.candidates[0].hap_to_alle_profile[2][0] == 0);
    CHECK(c.candidates[0].hap_to_alle_profile[2][1] == 0);
    CHECK(c.candidates[1].hap_to_alle_profile[2][1] == 1);
}

TEST_CASE("the phase-set argument has no upstream counterpart and filters when given") {
    // Retained for the graph arm's recovery; the BAM path must never pass it,
    // which is what keeps that path faithful. Pinned so the parameter's effect
    // is documented rather than discovered.
    PhasingChunk c = two_block_chunk({1, 1}, 1000, 2000);
    update_var_hap_profile_based_on_read_hap(c, 0, 1, kCandGermlineClean,
                                             std::optional<hts_pos_t>(1000));
    CHECK(c.candidates[0].hap_to_alle_profile[1][1] == 1);
    CHECK(c.candidates[1].hap_to_alle_profile[1][1] == 0);
}

// ── check_agree_haps (longcallD assign_hap.c:307) ────────────────────────────

TEST_CASE("check_agree_haps returns agree, conflict and uninformative") {
    // cons is (hap1, hap2) = (0, 1) at both sites.
    SECTION("both alleles match the same haplotype -> agree") {
        PhasingChunk c = two_block_chunk({0, 0}, 1000, 1000);
        CHECK(check_agree_haps(c, 0, 1, 0, 1) == 1);
    }
    SECTION("the second allele matches the other haplotype -> conflict") {
        PhasingChunk c = two_block_chunk({0, 1}, 1000, 1000);
        CHECK(check_agree_haps(c, 0, 1, 0, 1) == 0);
    }
    SECTION("hap 0 carries no claim -> uninformative") {
        PhasingChunk c = two_block_chunk({0, 0}, 1000, 1000);
        CHECK(check_agree_haps(c, 0, 0, 0, 1) == -1);
    }
    SECTION("an unobserved allele -> uninformative") {
        PhasingChunk c = two_block_chunk({-1, 0}, 1000, 1000);
        CHECK(check_agree_haps(c, 0, 1, 0, 1) == -1);
    }
    SECTION("a variant outside the read's span -> uninformative") {
        PhasingChunk c = two_block_chunk({0, 0}, 1000, 1000);
        CHECK(check_agree_haps(c, 0, 1, 0, 5) == -1);
    }
}

// ── init_assign_read_hap_based_on_cons_alle (longcallD assign_hap.c:151) ─────

TEST_CASE("a read with no variant profile cannot be labelled") {
    PhasingChunk c = two_block_chunk({0, 0}, 1000, 1000);
    c.read_var_profile[0].start_var_idx = -1;
    CHECK(init_assign_read_hap_based_on_cons_alle(c, 0, kCandGermlineClean) == -1);
}

TEST_CASE("a read matching one haplotype's consensus is labelled with it") {
    PhasingChunk c = two_block_chunk({0, 0}, 1000, 1000);   // hap1 cons is 0 at both
    CHECK(init_assign_read_hap_based_on_cons_alle(c, 0, kCandGermlineClean) == 1);
    PhasingChunk d = two_block_chunk({1, 1}, 1000, 1000);   // hap2 cons is 1 at both
    CHECK(init_assign_read_hap_based_on_cons_alle(d, 0, kCandGermlineClean) == 2);
}
