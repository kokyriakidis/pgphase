// Unit tests for the phasing predicates: small, pure-ish functions whose
// verdicts gate what enters the solve.
//
// These exist because a real defect survived every integration test in this
// repo. `var_is_homopolymer_indel`'s insertion branch compared the reference as
// a raw FASTA byte against an nt4-coded alt base, so it returned false for
// every insertion ever passed to it (commit 510f865). The flag gates four
// consumers -- read scoring, the link list, phase-set eligibility and pivot
// choice -- and the leak cost 3.007 points of chr20 read hamming. It is
// expressible as three lines of synthetic reference here, and needs no BAM.
#define CATCH_CONFIG_MAIN
#include "../third_party/catch2/catch.hpp"

#include "collect_phase.hpp"
#include "collect_phase_noisy.hpp"
#include "collect_var.hpp"
#include "noise_filter.hpp"
#include "phasing_types.hpp"

#include "cgranges.h"

#include <unistd.h>

#include <cstdio>

using namespace pgphase_collect;

namespace {

/// A chunk carrying only the reference slice the predicate reads.
PhasingChunk ref_chunk(hts_pos_t ref_beg, const std::string& seq) {
    PhasingChunk chunk;
    chunk.ref_beg = ref_beg;
    chunk.ref_end = ref_beg + static_cast<hts_pos_t>(seq.size());
    chunk.ref_seq = seq;
    return chunk;
}

bool hp_ins(const PhasingChunk& c, hts_pos_t pos, const std::string& alt) {
    return var_is_homopolymer_indel(c, pos, VariantType::Insertion, 0, alt);
}
bool hp_del(const PhasingChunk& c, hts_pos_t pos, int ref_len) {
    return var_is_homopolymer_indel(c, pos, VariantType::Deletion, ref_len, {});
}


/// A minimal indel key. `alt` holds the inserted bases for an insertion.
VariantKey ins_key(hts_pos_t pos, const std::string& alt) {
    VariantKey k;
    k.pos = pos;
    k.type = VariantType::Insertion;
    k.ref_len = 0;
    k.alt = alt;
    return k;
}
VariantKey del_key(hts_pos_t pos, int ref_len) {
    VariantKey k;
    k.pos = pos;
    k.type = VariantType::Deletion;
    k.ref_len = ref_len;
    return k;
}

/// A candidate that passes every allele_depths_call_het test, so each SECTION
/// can break exactly one of them.
CandidateVariant het_candidate(hts_pos_t pos) {
    CandidateVariant v;
    v.key.pos = pos;
    v.key.type = VariantType::Snp;
    v.key.ref_len = 1;
    v.lcd_var_i_to_cate = kCandCleanHetSnp;
    v.alignment_verified = true;
    v.hap_to_alle_profile[1].assign(2, 10);
    v.hap_to_alle_profile[2].assign(2, 10);
    v.counts.ref_cov = 15;
    v.counts.alt_cov = 15;
    v.counts.allele_fraction = 0.5;
    return v;
}

} // namespace

TEST_CASE("base_to_nt4 is case-insensitive and rejects ambiguity") {
    CHECK(base_to_nt4('A') == 0);
    CHECK(base_to_nt4('a') == 0);
    CHECK(base_to_nt4('C') == 1);
    CHECK(base_to_nt4('c') == 1);
    CHECK(base_to_nt4('G') == 2);
    CHECK(base_to_nt4('g') == 2);
    CHECK(base_to_nt4('T') == 3);
    CHECK(base_to_nt4('t') == 3);
    CHECK(base_to_nt4('U') == 3);
    // Anything else must land outside 0..3 so the predicates can reject it.
    CHECK(base_to_nt4('N') > 3);
    CHECK(base_to_nt4('n') > 3);
    CHECK(base_to_nt4('-') > 3);
}

TEST_CASE("var_is_homopolymer_indel: insertions") {
    // 1000       1010
    // |          |
    // GGGGGAAAAAAAAAAGGGGG
    const PhasingChunk up = ref_chunk(1000, "GGGGGAAAAAAAAAAGGGGG");

    SECTION("single base into a run of >= 5 is a homopolymer indel") {
        CHECK(hp_ins(up, 1005, "A"));
    }
    SECTION("multiple copies of the same base still qualify") {
        CHECK(hp_ins(up, 1005, "AA"));
        CHECK(hp_ins(up, 1005, "AAA"));
    }
    SECTION("a mixed insertion does not, even in a run") {
        CHECK_FALSE(hp_ins(up, 1005, "AC"));
        CHECK_FALSE(hp_ins(up, 1005, "CA"));
    }
    SECTION("the inserted base must match the run") {
        CHECK_FALSE(hp_ins(up, 1005, "C"));
    }
    SECTION("a run shorter than 5 from the position does not qualify") {
        // Only 4 A's remain from 1011 onward before the G's resume.
        CHECK_FALSE(hp_ins(up, 1011, "A"));
    }
    SECTION("no run at all") {
        const PhasingChunk mixed = ref_chunk(1000, "ACGTACGTACGTACGTACGT");
        CHECK_FALSE(hp_ins(mixed, 1005, "A"));
    }
    SECTION("an empty alt is not an insertion we can judge") {
        CHECK_FALSE(hp_ins(up, 1005, ""));
    }

    SECTION("REGRESSION: a soft-masked (lowercase) run still qualifies") {
        // The reference is lowercase in exactly the repeat tracts this asks
        // about. Comparing raw bytes failed here twice over: 'a' (97) never
        // equals the nt4 code 0, and it would not have matched 'A' either.
        const PhasingChunk lower = ref_chunk(1000, "ggggg" "aaaaaaaaaa" "ggggg");
        CHECK(hp_ins(lower, 1005, "A"));
        CHECK(hp_ins(lower, 1005, "a"));
    }

    SECTION("REGRESSION: the chr20 locus that inverted a flank") {
        // chr20:55,871,832-55,871,845 reads `gtctcaaaaaaaaa`, and the 1 bp
        // insertion at 55,871,837 sits at the head of the 8 bp A-run. It
        // segregates at 0.525 against read truth -- noise -- yet was flagged 0
        // and therefore scored reads, setting the parity for a whole flank.
        const PhasingChunk chr20 = ref_chunk(55871832, "gtctcaaaaaaaaa");
        CHECK(hp_ins(chr20, 55871837, "A"));
    }
}

TEST_CASE("var_is_homopolymer_indel: deletions") {
    const PhasingChunk up = ref_chunk(1000, "GGGGGAAAAAAAAAAGGGGG");

    SECTION("a single-base deletion inside a run qualifies") {
        CHECK(hp_del(up, 1005, 1));
    }
    SECTION("a multi-base deletion of identical bases qualifies") {
        CHECK(hp_del(up, 1005, 3));
    }
    SECTION("a deletion spanning two different bases does not") {
        // From 1013 the span crosses the A -> G boundary.
        CHECK_FALSE(hp_del(up, 1013, 3));
    }
    SECTION("soft-masked runs qualify") {
        const PhasingChunk lower = ref_chunk(1000, "ggggg" "aaaaaaaaaa" "ggggg");
        CHECK(hp_del(lower, 1005, 1));
    }
}

TEST_CASE("var_is_homopolymer_indel: rejects what it cannot judge") {
    const PhasingChunk up = ref_chunk(1000, "GGGGGAAAAAAAAAAGGGGG");

    SECTION("a SNP is never a homopolymer indel") {
        CHECK_FALSE(var_is_homopolymer_indel(up, 1005, VariantType::Snp, 1, "A"));
    }
    SECTION("a position before the slice") {
        CHECK_FALSE(hp_ins(up, 999, "A"));
        CHECK_FALSE(hp_del(up, 999, 1));
    }
    SECTION("a position whose context runs off the end of the slice") {
        // The predicate reads 5 bases; the slice ends at 1020.
        CHECK_FALSE(hp_ins(up, 1018, "G"));
        CHECK_FALSE(hp_del(up, 1018, 1));
        CHECK_FALSE(hp_del(up, 1015, 8));
    }
    SECTION("an ambiguous base in the context") {
        const PhasingChunk n = ref_chunk(1000, "GGGGGAANAAAAAAAGGGGG");
        CHECK_FALSE(hp_ins(n, 1005, "A"));
        CHECK_FALSE(hp_del(n, 1005, 1));
    }
    SECTION("an empty reference slice") {
        const PhasingChunk empty = ref_chunk(1000, "");
        CHECK_FALSE(hp_ins(empty, 1000, "A"));
        CHECK_FALSE(hp_del(empty, 1000, 1));
    }
}

TEST_CASE("select_stitch_orientation: net-margin rule (the default)") {
    Options opts;
    opts.stitch_rule = kStitchRuleNetMargin;
    opts.stitch_min_margin = 0;
    bool flip = false;

    SECTION("a clear same-orientation vote merges without flipping") {
        // n11 + n22 dominate: flip_hap_score = n12 + n21 - n11 - n22 < 0.
        CHECK(select_stitch_orientation({30, 0, 0, 25}, &opts, flip));
        CHECK_FALSE(flip);
    }
    SECTION("a clear crossed vote merges with a flip") {
        CHECK(select_stitch_orientation({0, 30, 25, 0}, &opts, flip));
        CHECK(flip);
    }
    SECTION("a tie abstains") {
        CHECK_FALSE(select_stitch_orientation({10, 10, 10, 10}, &opts, flip));
        CHECK_FALSE(select_stitch_orientation({0, 0, 0, 0}, &opts, flip));
    }
    SECTION("margin 0 merges on a single net vote") {
        CHECK(select_stitch_orientation({1, 0, 0, 0}, &opts, flip));
    }
    SECTION("a margin abstains on a net vote that does not exceed it") {
        opts.stitch_min_margin = 2;
        CHECK_FALSE(select_stitch_orientation({2, 0, 0, 0}, &opts, flip));
        CHECK(select_stitch_orientation({3, 0, 0, 0}, &opts, flip));
    }
    SECTION("a null Options behaves as margin 0 with the default rule") {
        CHECK(select_stitch_orientation({1, 0, 0, 0}, nullptr, flip));
        CHECK_FALSE(select_stitch_orientation({0, 0, 0, 0}, nullptr, flip));
    }
}

TEST_CASE("select_stitch_orientation: both-strands rule") {
    Options opts;
    opts.stitch_rule = kStitchRuleBothStrands;
    bool flip = false;

    SECTION("requires evidence on BOTH links of the winning orientation") {
        // 30 no-flip votes, but only one of the two links carries any.
        CHECK_FALSE(select_stitch_orientation({30, 0, 0, 0}, &opts, flip));
        CHECK(select_stitch_orientation({30, 0, 0, 1}, &opts, flip));
        CHECK_FALSE(flip);
    }
    SECTION("when both orientations qualify, the larger net wins") {
        CHECK(select_stitch_orientation({5, 20, 20, 5}, &opts, flip));
        CHECK(flip);
        CHECK(select_stitch_orientation({20, 5, 5, 20}, &opts, flip));
        CHECK_FALSE(flip);
    }
    SECTION("when both qualify equally, it abstains") {
        CHECK_FALSE(select_stitch_orientation({10, 10, 10, 10}, &opts, flip));
    }
}

TEST_CASE("select_stitch_orientation: literal and both-strands-margin rules") {
    bool flip = false;

    SECTION("the literal rule stitches contested seams on purpose") {
        Options opts;
        opts.stitch_rule = kStitchRuleLiteral;
        CHECK(select_stitch_orientation({1, 1, 0, 0}, &opts, flip));
        // One vote each way: the net decides the orientation.
        CHECK(select_stitch_orientation({1, 5, 0, 0}, &opts, flip));
        CHECK(flip);
        // Nothing on one side: no merge.
        CHECK_FALSE(select_stitch_orientation({5, 0, 0, 5}, &opts, flip));
    }
    SECTION("both-strands-margin needs both links AND a net above the margin") {
        Options opts;
        opts.stitch_rule = kStitchRuleBothStrandsMargin;
        opts.stitch_min_margin = 3;
        // Net is 20 and both crossed links carry votes.
        CHECK(select_stitch_orientation({0, 10, 10, 0}, &opts, flip));
        CHECK(flip);
        // Net is only 2, below the margin.
        CHECK_FALSE(select_stitch_orientation({4, 5, 3, 4}, &opts, flip));
        // Net is large but one crossed link is empty.
        CHECK_FALSE(select_stitch_orientation({0, 20, 0, 0}, &opts, flip));
    }
}


TEST_CASE("var_is_homopolymer_pg: reference-only STR detection") {
    // 1000            1015
    const std::string up = "CCCCCAAAAAAAAAAGGGGGGGGGG";
    const hts_pos_t beg = 1000, end = beg + static_cast<hts_pos_t>(up.size());

    SECTION("a deletion inside a long run is detected") {
        CHECK(var_is_homopolymer_pg(del_key(1007, 1), up, beg, end, 5));
    }
    SECTION("soft-masked reference is detected too (nt4 comparison)") {
        std::string lower = up;
        for (char& c : lower) c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
        CHECK(var_is_homopolymer_pg(del_key(1007, 1), lower, beg, end, 5));
    }
    SECTION("a multi-base repeat unit is detected") {
        // AT x 6 is a unit-2 STR with more than three copies.
        const std::string str2 = "CCCCCATATATATATATCCCCC";
        CHECK(var_is_homopolymer_pg(del_key(1007, 2), str2, beg,
                                    beg + static_cast<hts_pos_t>(str2.size()), 5));
    }
    SECTION("a non-repeat context is not") {
        const std::string mixed = "ACGTACGGTTAACCGGTTACGT";
        CHECK_FALSE(var_is_homopolymer_pg(del_key(1010, 1), mixed, beg,
                                          beg + static_cast<hts_pos_t>(mixed.size()), 5));
    }
    SECTION("an indel longer than xid is not judged") {
        CHECK_FALSE(var_is_homopolymer_pg(del_key(1007, 6), up, beg, end, 5));
        CHECK_FALSE(var_is_homopolymer_pg(ins_key(1007, "AAAAAA"), up, beg, end, 5));
    }
    SECTION("an empty reference slice") {
        CHECK_FALSE(var_is_homopolymer_pg(del_key(1007, 1), "", beg, end, 5));
    }
}

TEST_CASE("var_is_repeat_region_pg: three tandem copies of the indel motif") {
    // The documented example: delete "AT" where the reference continues ATATAT.
    const std::string str2 = "GGGGG" "ATATATATATATATAT" "GGGGG";
    const hts_pos_t beg = 1000, end = beg + static_cast<hts_pos_t>(str2.size());

    SECTION("a deletion whose motif repeats three times downstream") {
        CHECK(var_is_repeat_region_pg(del_key(1005, 2), str2, beg, end, 5));
    }
    SECTION("an insertion consistent with the tandem repeat") {
        CHECK(var_is_repeat_region_pg(ins_key(1005, "AT"), str2, beg, end, 5));
    }
    SECTION("an insertion of a different motif is not") {
        CHECK_FALSE(var_is_repeat_region_pg(ins_key(1005, "GC"), str2, beg, end, 5));
    }
    SECTION("a non-repeat context is not") {
        const std::string mixed = "GGGGG" "ACGTTGCAACGTTGCA" "GGGGG";
        CHECK_FALSE(var_is_repeat_region_pg(del_key(1005, 2), mixed, beg,
                                            beg + static_cast<hts_pos_t>(mixed.size()), 5));
    }
    SECTION("an indel longer than xid is not judged") {
        CHECK_FALSE(var_is_repeat_region_pg(del_key(1005, 6), str2, beg, end, 5));
        CHECK_FALSE(var_is_repeat_region_pg(ins_key(1005, "ATATAT"), str2, beg, end, 5));
    }
    SECTION("a motif whose three copies run past the slice is not judged") {
        CHECK_FALSE(var_is_repeat_region_pg(del_key(end - 3, 2), str2, beg, end, 5));
        CHECK_FALSE(var_is_repeat_region_pg(ins_key(end - 3, "AT"), str2, beg, end, 5));
    }
    SECTION("a position before the slice") {
        CHECK_FALSE(var_is_repeat_region_pg(del_key(999, 2), str2, beg, end, 5));
    }

    SECTION("a deletion whose two windows straddle a soft-mask boundary") {
        // Half the tract masked: a byte comparison would call this no-repeat.
        std::string half = "GGGGG" "ATATAT" "atatatatat" "GGGGG";
        CHECK(var_is_repeat_region_pg(del_key(1005, 2), half, beg,
                                      beg + static_cast<hts_pos_t>(half.size()), 5));
    }
    SECTION("the call site's OR is insensitive to the insertion branch here") {
        // collect_var.cpp:1490-1492 classifies RepeatHetIndel on
        //   var_is_homopolymer_pg(...) || var_is_repeat_region_pg(...)
        // and the first is a unit-1..6 STR test read through nt4. For a tandem
        // repeat like this one it already returns true, case-insensitively, so
        // the insertion branch's case bug was LATENT at the call site: fixing it
        // left chr20 output byte-identical. This pins the redundancy, so a future
        // change that narrows var_is_homopolymer_pg cannot silently expose it.
        std::string lower = str2;
        for (char& c : lower) c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
        CHECK(var_is_homopolymer_pg(ins_key(1005, "AT"), lower, beg, end, 5));
        CHECK(var_is_homopolymer_pg(del_key(1005, 2), lower, beg, end, 5));
    }

    SECTION("a run of N is not a tandem repeat") {
        const std::string ns = "GGGGG" "NNNNNNNNNNNNNNNN" "GGGGG";
        CHECK_FALSE(var_is_repeat_region_pg(del_key(1005, 2), ns, beg,
                                            beg + static_cast<hts_pos_t>(ns.size()), 5));
        CHECK_FALSE(var_is_repeat_region_pg(ins_key(1005, "NN"), ns, beg,
                                            beg + static_cast<hts_pos_t>(ns.size()), 5));
    }

    SECTION("soft-masked reference must be judged the same as uppercase") {
        // The reference is lowercase in exactly the tandem repeats this asks
        // about, and candidate alt bases are uppercased upstream
        // (bam_digar.cpp:323). Its sibling var_is_homopolymer_pg compares
        // through nt4 and is case-insensitive, so the OR at
        // collect_var.cpp:1490-1492 must not depend on which of the two answers.
        std::string lower = str2;
        for (char& c : lower) c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
        CHECK(var_is_repeat_region_pg(del_key(1005, 2), lower, beg, end, 5));
        CHECK(var_is_repeat_region_pg(ins_key(1005, "AT"), lower, beg, end, 5));
    }
}

TEST_CASE("allele_depths_call_het: every exclusion in order") {
    Options opts;
    opts.min_alt_depth = 2;
    opts.min_af = 0.2;
    opts.max_af = 0.8;
    opts.retry_windows = {{1000, 2000}};

    SECTION("a clear het inside a retry window is admitted") {
        CHECK(allele_depths_call_het(het_candidate(1500), opts));
    }
    SECTION("off entirely with no window and no joint flag") {
        Options none = opts;
        none.retry_windows.clear();
        CHECK_FALSE(allele_depths_call_het(het_candidate(1500), none));
    }
    SECTION("a catalog-only row is not repaired") {
        CandidateVariant v = het_candidate(1500);
        v.alignment_verified = false;
        CHECK_FALSE(allele_depths_call_het(v, opts));
    }
    SECTION("the position must fall inside a window") {
        CHECK_FALSE(allele_depths_call_het(het_candidate(2500), opts));
        // The window is half-open: [beg, end).
        CHECK(allele_depths_call_het(het_candidate(1000), opts));
        CHECK_FALSE(allele_depths_call_het(het_candidate(2000), opts));
    }
    SECTION("joint orientation makes the verdict chunk-wide") {
        Options joint = opts;
        joint.retry_windows.clear();
        joint.joint_het_orientation = true;
        CHECK(allele_depths_call_het(het_candidate(999999), joint));
    }
    SECTION("category must be a het class") {
        CandidateVariant v = het_candidate(1500);
        v.lcd_var_i_to_cate = kCandCleanHom;
        CHECK_FALSE(allele_depths_call_het(v, opts));
        v.lcd_var_i_to_cate = kCandNoisyCandHom;
        CHECK_FALSE(allele_depths_call_het(v, opts));
        v.lcd_var_i_to_cate = kCandNoisyCandHet;
        CHECK(allele_depths_call_het(v, opts));
    }
    SECTION("a multiallelic record is excluded, being oriented jointly") {
        CandidateVariant v = het_candidate(1500);
        v.msa_insertion_alts = {"A", "AA"};
        CHECK_FALSE(allele_depths_call_het(v, opts));
    }
    SECTION("both haplotype profiles need at least two observations") {
        CandidateVariant v = het_candidate(1500);
        v.hap_to_alle_profile[2].assign(1, 10);
        CHECK_FALSE(allele_depths_call_het(v, opts));
        v.hap_to_alle_profile[2].clear();
        CHECK_FALSE(allele_depths_call_het(v, opts));
    }
    SECTION("reference and alternate depth must both reach min_alt_depth") {
        CandidateVariant v = het_candidate(1500);
        v.counts.ref_cov = 1;
        CHECK_FALSE(allele_depths_call_het(v, opts));
        v.counts.ref_cov = 15;
        v.counts.alt_cov = 1;
        CHECK_FALSE(allele_depths_call_het(v, opts));
    }
    SECTION("allele fraction must sit inside [min_af, max_af]") {
        CandidateVariant v = het_candidate(1500);
        v.counts.allele_fraction = 0.1;
        CHECK_FALSE(allele_depths_call_het(v, opts));
        v.counts.allele_fraction = 0.95;
        CHECK_FALSE(allele_depths_call_het(v, opts));
        // A homopolymer indel also needs the outside-in frontier to validate
        // its boundary link; depth alone cannot admit a repeat row.
        v.counts.allele_fraction = 0.5;
        v.is_homopolymer_indel = true;
        CHECK_FALSE(allele_depths_call_het(v, opts));
        v.gap_link_supported = true;
        CHECK(allele_depths_call_het(v, opts));
    }
}


TEST_CASE("MSA indel retry checks local and crossing dropout independently",
          "[recovery][msa][admission]") {
    PhasingChunk chunk;
    CandidateVariant left = het_candidate(1000);
    left.key = del_key(1000, 1);
    left.msa_verified = true;
    CandidateVariant contrast = het_candidate(1000);
    contrast.key.type = VariantType::Insertion;
    contrast.key.ref_len = 0;
    contrast.msa_verified = true;
    contrast.lcd_var_i_to_cate = kCandNoisyCandHet;
    chunk.candidates = {left, contrast, het_candidate(2000)};
    for (int ri = 0; ri < 65; ++ri) {
        ReadRecord read;
        read.beg = 900;
        read.end = ri < 8 ? 2100 : 1100;
        read.mapq = 60;
        chunk.reads.push_back(std::move(read));
        ReadVariantProfile profile;
        profile.start_var_idx = 0;
        profile.end_var_idx = 2;
        profile.alleles = {ri >= 8 && ri < 16 ? 0 : -1, 0, 0};
        chunk.read_var_profile.push_back(profile);
    }
    Options opts;
    SECTION("deep local dropout and eight missing crossing calls admit retry") {
        CHECK(msa_boundary_dropout_is_supported(chunk, {0, 2}, opts));
    }
    SECTION("normal local coverage cannot be hidden by crossing dropout") {
        for (size_t ri = 8; ri < chunk.reads.size(); ++ri)
            chunk.read_var_profile[ri].alleles[0] = 0;
        CHECK_FALSE(msa_boundary_dropout_is_supported(chunk, {0, 2}, opts));
    }
    SECTION("a crossing observation prevents a falsely decisive dropout") {
        for (int allele : {0, 1, 2}) {
            chunk.read_var_profile[0].alleles[0] = allele;
            CHECK_FALSE(msa_boundary_dropout_is_supported(chunk, {0, 2}, opts));
        }
    }
    SECTION("a missing profile counts as dropout") {
        chunk.read_var_profile[0] = ReadVariantProfile{};
        CHECK(msa_boundary_dropout_is_supported(chunk, {0, 2}, opts));
    }
    SECTION("unknown low and skipped mapping coverage is excluded") {
        for (int mapq : {255, 29}) {
            chunk.reads[0].mapq = mapq;
            CHECK_FALSE(msa_boundary_dropout_is_supported(chunk, {0, 2}, opts));
        }
        chunk.reads[0].mapq = 60;
        chunk.reads[0].is_skipped = true;
        CHECK_FALSE(msa_boundary_dropout_is_supported(chunk, {0, 2}, opts));
    }
    SECTION("low local coverage cannot admit retry") {
        chunk.reads.resize(19);
        chunk.read_var_profile.resize(19);
        CHECK_FALSE(msa_boundary_dropout_is_supported(chunk, {0, 2}, opts));
    }
    SECTION("configured depth also applies to crossing coverage") {
        opts.min_depth = 9;
        CHECK_FALSE(msa_boundary_dropout_is_supported(chunk, {0, 2}, opts));
    }
    SECTION("unverified sites and SNPs keep their existing admission path") {
        chunk.candidates[0].msa_verified = false;
        CHECK_FALSE(msa_boundary_dropout_is_supported(chunk, {0, 2}, opts));
        chunk.candidates[0].msa_verified = true;
        chunk.candidates[0].key.type = VariantType::Snp;
        CHECK_FALSE(msa_boundary_dropout_is_supported(chunk, {0, 2}, opts));
    }
    SECTION("a verified diploid contrast is required even when its genotype collapsed") {
        chunk.candidates[1].hap_to_cons_alle = {-1, 0, 0};
        CHECK(msa_boundary_dropout_is_supported(chunk, {0, 2}, opts));
        chunk.candidates[1].msa_verified = false;
        CHECK_FALSE(msa_boundary_dropout_is_supported(chunk, {0, 2}, opts));
        chunk.candidates[1].msa_verified = true;
        chunk.candidates[1].key.type = VariantType::Deletion;
        CHECK_FALSE(msa_boundary_dropout_is_supported(chunk, {0, 2}, opts));
        chunk.candidates[1].key.type = VariantType::Insertion;
        chunk.candidates[1].counts.n_uniq_alles = 3;
        CHECK_FALSE(msa_boundary_dropout_is_supported(chunk, {0, 2}, opts));
        chunk.candidates[1].counts.n_uniq_alles = 2;
        chunk.candidates[1].counts.alt_cov = opts.min_alt_depth - 1;
        CHECK_FALSE(msa_boundary_dropout_is_supported(chunk, {0, 2}, opts));
        chunk.candidates[1].counts.alt_cov = opts.min_alt_depth;
        chunk.candidates[1].counts.allele_fraction = 0;
        CHECK_FALSE(msa_boundary_dropout_is_supported(chunk, {0, 2}, opts));
    }
    SECTION("a complementary same-type MSA pair can request diploid retry") {
        chunk.candidates[0].phase_set = 1000;
        chunk.candidates[0].hap_to_cons_alle = {-1, 0, 1};
        chunk.candidates[1].key = del_key(1000, 2);
        chunk.candidates[1].phase_set = 1000;
        chunk.candidates[1].hap_to_cons_alle = {-1, 1, 0};
        CHECK_FALSE(msa_boundary_dropout_is_supported(chunk, {0, 2}, opts));
        CHECK(msa_boundary_dropout_is_supported(chunk, {0, 2}, opts, true));
        chunk.candidates[1].phase_set = 2000;
        CHECK_FALSE(msa_boundary_dropout_is_supported(chunk, {0, 2}, opts, true));
        chunk.candidates[1].phase_set = 1000;
        chunk.candidates[1].hap_to_cons_alle = {-1, 0, 1};
        CHECK_FALSE(msa_boundary_dropout_is_supported(chunk, {0, 2}, opts, true));
        chunk.candidates[1].hap_to_cons_alle = {-1, 0, 0};
        CHECK_FALSE(msa_boundary_dropout_is_supported(chunk, {0, 2}, opts, true));
        chunk.candidates[1].hap_to_cons_alle = {-1, 1, 0};
        chunk.candidates[1].key.ref_len = 1;
        CHECK_FALSE(msa_boundary_dropout_is_supported(chunk, {0, 2}, opts, true));
    }
    SECTION("either boundary can show the independent dropout") {
        chunk.candidates[0].key.type = VariantType::Snp;
        chunk.candidates[2].key = del_key(2000, 1);
        chunk.candidates[1].key.pos = 2000;
        chunk.candidates[2].msa_verified = true;
        for (size_t ri = 0; ri < chunk.reads.size(); ++ri) {
            chunk.reads[ri].beg = ri < 8 ? 900 : 1900;
            chunk.reads[ri].end = 2100;
            chunk.read_var_profile[ri].alleles[2] =
                chunk.read_var_profile[ri].alleles[0];
            chunk.read_var_profile[ri].alleles[0] = 0;
        }
        CHECK(msa_boundary_dropout_is_supported(chunk, {0, 2}, opts));
    }
}

TEST_CASE("joint genotype orientation retains a depth-supported BAM repeat bridge",
          "[recovery][genotype][phase-set]") {
    PhasingChunk chunk;
    CandidateVariant left = het_candidate(1000);
    CandidateVariant bridge = het_candidate(1500);
    CandidateVariant right = het_candidate(2000);
    bridge.key = del_key(1500, 1);
    bridge.lcd_var_i_to_cate = kCandNoisyCandHet;
    bridge.msa_verified = true;
    bridge.alignment_verified = false;
    bridge.is_homopolymer_indel = true;
    chunk.candidates = {left, bridge, right};
    for (auto& candidate : chunk.candidates)
        candidate.hap_to_cons_alle = {-1, 0, 1};
    chunk.read_var_cr.reset(cr_init());
    for (int side = 0; side < 2; ++side) {
        for (int allele = 0; allele < 2; ++allele) {
            for (int replicate = 0; replicate < 4; ++replicate) {
                const int ri = static_cast<int>(chunk.reads.size());
                chunk.reads.emplace_back();
                chunk.haps.push_back(allele + 1);
                ReadVariantProfile profile;
                profile.start_var_idx = side;
                profile.end_var_idx = side + 1;
                profile.alleles = {allele, allele};
                chunk.read_var_profile.push_back(profile);
                cr_add(chunk.read_var_cr.get(), "cr", side, side + 2, ri);
            }
        }
    }
    cr_index(chunk.read_var_cr.get());
    Options opts;
    opts.upstream_assign_hap = true;
    opts.joint_het_orientation = true;
    SECTION("explicit joint orientation lets the retained genotype participate") {
        REQUIRE(allele_depths_call_het(chunk.candidates[1], opts));
        iter_update_var_hap_cons_phase_set(chunk, {0, 1, 2}, opts);
        CHECK(chunk.candidates[2].phase_set == chunk.candidates[0].phase_set);
    }
    SECTION("allele linking preserves an opposite connection through the bridge") {
        opts.link_by_alleles = true;
        for (size_t ri = 8; ri < chunk.read_var_profile.size(); ++ri) {
            chunk.read_var_profile[ri].alleles[1] =
                1 - chunk.read_var_profile[ri].alleles[0];
            chunk.haps[ri] = 3 - chunk.haps[ri];
        }
        iter_update_var_hap_cons_phase_set(chunk, {0, 1, 2}, opts);
        CHECK(chunk.candidates[2].phase_set == chunk.candidates[0].phase_set);
        CHECK(chunk.candidates[0].hap_to_cons_alle[1] !=
              chunk.candidates[2].hap_to_cons_alle[1]);
    }
    SECTION("ordinary BAM mode keeps the upstream repeat exclusion") {
        opts.joint_het_orientation = false;
        iter_update_var_hap_cons_phase_set(chunk, {0, 1, 2}, opts);
        CHECK(chunk.candidates[2].phase_set != chunk.candidates[0].phase_set);
    }
    SECTION("an allele-depth failure still excludes the repeat bridge") {
        chunk.candidates[1].counts.alt_cov = 0;
        REQUIRE_FALSE(allele_depths_call_het(chunk.candidates[1], opts));
        iter_update_var_hap_cons_phase_set(chunk, {0, 1, 2}, opts);
        CHECK(chunk.candidates[2].phase_set != chunk.candidates[0].phase_set);
    }
}

TEST_CASE("is_repeat_indel: VCF-anchored tandem repeat test") {
    // Anchor at `pos`; indel content starts at pos+1.
    // 1000  1005
    const std::string str2 = "GGGGG" "ATATATATATATATAT" "GGGGG";
    const hts_pos_t beg = 1000, end = beg + static_cast<hts_pos_t>(str2.size());

    SECTION("a deletion of one repeat unit") {
        // anchor G at 1004, delete the AT at 1005-1006
        CHECK(is_repeat_indel(1004, "GAT", "G", str2, beg, end, 5));
    }
    SECTION("an insertion of one repeat unit") {
        CHECK(is_repeat_indel(1004, "G", "GAT", str2, beg, end, 5));
    }
    SECTION("a SNP is not an indel") {
        CHECK_FALSE(is_repeat_indel(1004, "G", "A", str2, beg, end, 5));
    }
    SECTION("a different motif is not a repeat here") {
        CHECK_FALSE(is_repeat_indel(1004, "G", "GCC", str2, beg, end, 5));
    }
    SECTION("longer than max_xgaps is not judged") {
        CHECK_FALSE(is_repeat_indel(1004, "G", "GATATATATATAT", str2, beg, end, 5));
    }
    SECTION("an empty reference slice") {
        CHECK_FALSE(is_repeat_indel(1004, "GAT", "G", "", beg, end, 5));
    }

    SECTION("REGRESSION: soft-masked reference, both branches") {
        std::string lower = str2;
        for (char& c : lower) c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
        CHECK(is_repeat_indel(1004, "gat", "g", lower, beg, end, 5));
        CHECK(is_repeat_indel(1004, "G", "GAT", lower, beg, end, 5));
    }
    SECTION("REGRESSION: a deletion straddling a soft-mask boundary") {
        // The deletion branch compared raw bytes, so a masked/unmasked join
        // read as no-repeat while the insertion branch (already nt4) did not.
        const std::string half = "GGGGG" "ATATAT" "atatatatat" "GGGGG";
        CHECK(is_repeat_indel(1004, "GAT", "G", half, beg,
                              beg + static_cast<hts_pos_t>(half.size()), 5));
    }
    SECTION("REGRESSION: a run of N is not a tandem repeat") {
        // memcmp compared N to N and returned equal.
        const std::string ns = "GGGGG" "NNNNNNNNNNNNNNNN" "GGGGG";
        CHECK_FALSE(is_repeat_indel(1004, "GNN", "G", ns, beg,
                                    beg + static_cast<hts_pos_t>(ns.size()), 5));
    }
}

TEST_CASE("find_low_complexity_intervals and pos_in_low_complexity") {
    const hts_pos_t beg = 1000;
    const std::string run = "ACGTCAGGTCAGT" + std::string(60, 'A') + "TGACTGCATGCAT";

    SECTION("a long homopolymer is reported as low complexity") {
        const auto iv = find_low_complexity_intervals(run, beg);
        REQUIRE_FALSE(iv.empty());
        CHECK(pos_in_low_complexity(beg + 40, iv));
        CHECK_FALSE(pos_in_low_complexity(beg + 1, iv));
    }
    SECTION("soft-masked input gives the same intervals") {
        // sdust's seq_nt4_table maps lowercase acgt to 0..3, so masking must not
        // change the verdict. This pins that assumption.
        std::string lower = run;
        for (char& c : lower) c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
        const auto a = find_low_complexity_intervals(run, beg);
        const auto b = find_low_complexity_intervals(lower, beg);
        REQUIRE(a.size() == b.size());
        for (size_t i = 0; i < a.size(); ++i) {
            CHECK(a[i].beg == b[i].beg);
            CHECK(a[i].end == b[i].end);
        }
    }
    SECTION("an empty slice yields nothing") {
        CHECK(find_low_complexity_intervals("", beg).empty());
        CHECK_FALSE(pos_in_low_complexity(beg, {}));
    }
}

TEST_CASE("trim_to_minimal_vcf") {
    SECTION("a right-trimmable deletion") {
        hts_pos_t pos = 100;
        std::string ref = "ATTTT", alt = "ATTT";
        trim_to_minimal_vcf(pos, ref, alt);
        CHECK(ref.size() > alt.size());
        CHECK(ref.size() - alt.size() == 1);
    }
    SECTION("a left-trimmable insertion keeps one anchor base") {
        hts_pos_t pos = 100;
        std::string ref = "GGA", alt = "GGAA";
        trim_to_minimal_vcf(pos, ref, alt);
        CHECK(alt.size() - ref.size() == 1);
        CHECK(pos >= 100);
    }
    SECTION("a SNP is left alone") {
        hts_pos_t pos = 100;
        std::string ref = "A", alt = "G";
        trim_to_minimal_vcf(pos, ref, alt);
        CHECK(ref == "A");
        CHECK(alt == "G");
        CHECK(pos == 100);
    }
}

// ---------------------------------------------------------------------------
// verify_chunk_invariants: the properties a chunk must still satisfy after the
// in-chunk recovery inserts candidates and reads into it. Each case here
// corresponds to a defect that shipped and was found by its symptom -- a phase
// block count, a missing record -- rather than by a check.
// ---------------------------------------------------------------------------

namespace {

PhasingChunk three_site_chunk() {
    PhasingChunk chunk;
    chunk.ref_beg = 1000;
    chunk.ref_end = 2000;
    for (int i = 0; i < 3; ++i) {
        CandidateVariant cand;
        cand.key.tid = 0;
        cand.key.pos = 1000 + 10 * i;
        cand.key.type = VariantType::Snp;
        cand.key.ref_len = 1;
        cand.key.alt = "A";
        chunk.candidates.push_back(cand);
        ReadRecord read;
        read.qname = std::string("read") + static_cast<char>('a' + i);
        read.beg = 1000;
        read.end = 1100;
        chunk.reads.push_back(std::move(read));
        ReadVariantProfile prof;
        prof.read_id = i;
        chunk.read_var_profile.push_back(prof);
    }
    chunk.haps.assign(3, 0);
    chunk.phase_sets.assign(3, -1);
    return chunk;
}

}  // namespace

TEST_CASE("verify_chunk_invariants accepts a consistent chunk", "[invariants]") {
    const PhasingChunk chunk = three_site_chunk();
    CHECK_NOTHROW(verify_chunk_invariants(chunk, 3, 3, 3, 0, 0));
    // A re-solved region inside the chunk is fine.
    CHECK_NOTHROW(verify_chunk_invariants(chunk, 3, 3, 3, 1200, 1400));
}

TEST_CASE("verify_chunk_invariants catches the defects that shipped", "[invariants]") {
    SECTION("a short per-site array shifts metadata onto the wrong candidate") {
        const PhasingChunk chunk = three_site_chunk();
        CHECK_THROWS_AS(verify_chunk_invariants(chunk, 3, 2, 3, 0, 0), std::runtime_error);
        CHECK_THROWS_AS(verify_chunk_invariants(chunk, 2, 3, 3, 0, 0), std::runtime_error);
        CHECK_THROWS_AS(verify_chunk_invariants(chunk, 3, 3, 2, 0, 0), std::runtime_error);
    }
    SECTION("reads out of qname order break the cross-chunk stitch's merge-join") {
        PhasingChunk chunk = three_site_chunk();
        std::swap(chunk.reads[0].qname, chunk.reads[2].qname);
        CHECK_THROWS_AS(verify_chunk_invariants(chunk, 3, 3, 3, 0, 0), std::runtime_error);
    }
    SECTION("candidates out of position order break the solve's sweep") {
        PhasingChunk chunk = three_site_chunk();
        std::swap(chunk.candidates[0].key.pos, chunk.candidates[2].key.pos);
        CHECK_THROWS_AS(verify_chunk_invariants(chunk, 3, 3, 3, 0, 0), std::runtime_error);
    }
    SECTION("a read-indexed vector of the wrong length") {
        PhasingChunk chunk = three_site_chunk();
        chunk.read_var_profile.pop_back();
        CHECK_THROWS_AS(verify_chunk_invariants(chunk, 3, 3, 3, 0, 0), std::runtime_error);
    }
    SECTION("a profile whose read_id no longer matches its slot") {
        PhasingChunk chunk = three_site_chunk();
        chunk.read_var_profile[1].read_id = 7;
        CHECK_THROWS_AS(verify_chunk_invariants(chunk, 3, 3, 3, 0, 0), std::runtime_error);
    }
    SECTION("BAM qualities must use the same candidate offsets as alleles") {
        PhasingChunk chunk = three_site_chunk();
        auto& profile = chunk.read_var_profile[0];
        profile.start_var_idx = 0;
        profile.end_var_idx = 1;
        profile.alleles = {0, 1};
        profile.alt_qi = {10, 11};
        profile.bam_base_qualities = {40};
        CHECK_THROWS_AS(verify_chunk_invariants(chunk, 3, 3, 3, 0, 0), std::runtime_error);
        profile.bam_base_qualities = {40, 5};
        CHECK_NOTHROW(verify_chunk_invariants(chunk, 3, 3, 3, 0, 0));
        profile.bam_base_qualities.clear();
        CHECK_NOTHROW(verify_chunk_invariants(chunk, 3, 3, 3, 0, 0));
    }
    SECTION("a re-solved region reaching outside the chunk") {
        const PhasingChunk chunk = three_site_chunk();
        CHECK_THROWS_AS(verify_chunk_invariants(chunk, 3, 3, 3, 900, 1400), std::runtime_error);
        CHECK_THROWS_AS(verify_chunk_invariants(chunk, 3, 3, 3, 1200, 2400), std::runtime_error);
    }
}


// ---------------------------------------------------------------------------
// drop_conflicting_haplotype_alleles
//
// One haplotype carries one allele. Chromosome 20 emitted 106 positions that
// broke that: 96 a SNP together with the longer allele containing it, 9 an
// insertion anchored on the reference base together with a SNP changing it.
// Read counts at those loci (288,018: 40 reference, 24 GAG, none carrying a
// bare G) say the alt-carrying reads carry the whole insertion, so the
// contained record is the one to drop.
//
// The records are built the way the pipeline holds them, confirmed by probing a
// real run at 288,018: both are VariantType::Snp at the SAME key.pos, differing
// in ref_len (1 against 3) and alt ('G' against 'GAG'). That matters, because
// grouping is by VariantKey::sort_pos(), which subtracts one for a non-SNP --
// building the longer allele as an Insertion at the same pos puts the two in
// DIFFERENT groups and nothing is compared.
// ---------------------------------------------------------------------------
namespace {

pgphase_collect::CandidateVariant conflict_record(hts_pos_t pos, int ref_len,
                                                  const std::string& alt, int hap,
                                                  int alt_depth,
                                                  pgphase_collect::VariantType type =
                                                      pgphase_collect::VariantType::Snp) {
    pgphase_collect::CandidateVariant c;
    c.key.tid = 0;
    c.key.pos = pos;
    c.key.ref_len = ref_len;
    c.key.alt = alt;
    c.key.type = type;
    c.hap_to_cons_alle = {0, 0, 0};
    c.hap_to_cons_alle[static_cast<size_t>(hap)] = 1;
    c.counts.alle_covs = {10, alt_depth};
    return c;
}

std::vector<std::string> alts_of(const pgphase_collect::CandidateTable& t) {
    std::vector<std::string> out;
    for (const auto& c : t) out.push_back(c.key.alt);
    return out;
}

}  // namespace

TEST_CASE("a contained allele loses to the complete one on the same haplotype",
          "[conflict]") {
    pgphase_collect::CandidateTable t{conflict_record(288018, 1, "G", 1, 21),
                                      conflict_record(288018, 3, "GAG", 1, 18)};
    pgphase_collect::drop_conflicting_haplotype_alleles(t);
    REQUIRE(t.size() == 1);
    CHECK(alts_of(t) == std::vector<std::string>{"GAG"});
}

TEST_CASE("the same two alleles on OPPOSITE haplotypes are both kept",
          "[conflict]") {
    // 1|2 is a real genotype: each haplotype has one allele and neither
    // contradicts the other, so this must survive untouched.
    pgphase_collect::CandidateTable t{conflict_record(288018, 1, "G", 1, 21),
                                      conflict_record(288018, 3, "GAG", 2, 18)};
    pgphase_collect::drop_conflicting_haplotype_alleles(t);
    CHECK(t.size() == 2);
}

TEST_CASE("when neither allele contains the other, allele depth decides",
          "[conflict]") {
    // 1,907,006: the insertion carries 10 against the SNP's 5.
    pgphase_collect::CandidateTable t{conflict_record(1907006, 1, "ATCCATCC", 2, 10),
                                      conflict_record(1907006, 1, "G", 2, 5)};
    pgphase_collect::drop_conflicting_haplotype_alleles(t);
    REQUIRE(t.size() == 1);
    CHECK(alts_of(t) == std::vector<std::string>{"ATCCATCC"});
}

TEST_CASE("equal depth breaks towards the shorter allele, deterministically",
          "[conflict]") {
    pgphase_collect::CandidateTable t{conflict_record(500, 1, "ATT", 1, 7),
                                      conflict_record(500, 1, "G", 1, 7)};
    pgphase_collect::drop_conflicting_haplotype_alleles(t);
    REQUIRE(t.size() == 1);
    CHECK(alts_of(t) == std::vector<std::string>{"G"});
}

TEST_CASE("records at different positions are never compared", "[conflict]") {
    pgphase_collect::CandidateTable t{conflict_record(100, 1, "G", 1, 9),
                                      conflict_record(200, 1, "GAG", 1, 9)};
    pgphase_collect::drop_conflicting_haplotype_alleles(t);
    CHECK(t.size() == 2);
}

TEST_CASE("an indel anchors one base earlier, and still groups with its SNP",
          "[conflict]") {
    // sort_pos() subtracts one for a non-SNP, so the insertion recorded at
    // pos+1 shares a group with the SNP at pos -- the real pairing in the data.
    pgphase_collect::CandidateTable t{
        conflict_record(700, 1, "G", 1, 12),
        conflict_record(701, 0, "GA", 1, 9, pgphase_collect::VariantType::Insertion)};
    pgphase_collect::drop_conflicting_haplotype_alleles(t);
    REQUIRE(t.size() == 1);
    CHECK(alts_of(t) == std::vector<std::string>{"GA"});
}

namespace {

PhasingChunk shifted_insertion_chunk(const std::string& cigar = "6M1I2M",
                                     const std::string& sequence = "CATTTTTGC") {
    PhasingChunk chunk = ref_chunk(100, "CATTTTGC");
    CandidateVariant candidate;
    candidate.key = ins_key(103, "T");
    candidate.msa_verified = true;
    candidate.lcd_var_i_to_cate = kCandNoisyCandHet;
    chunk.candidates.push_back(candidate);
    ReadRecord read;
    read.qname = "shifted";
    read.beg = 100;
    read.end = 107;
    read.mapq = 60;
    read.alignment.reset(bam_init1());
    const std::string header_text = "@SQ\tSN:chr\tLN:1000\n";
    std::unique_ptr<bam_hdr_t, decltype(&bam_hdr_destroy)> header(
        sam_hdr_parse(header_text.size(), header_text.c_str()), &bam_hdr_destroy);
    REQUIRE(header != nullptr);
    const std::string sam = "shifted\t0\tchr\t100\t60\t" + cigar +
        "\t*\t0\t0\t" + sequence + "\t" + std::string(sequence.size(), 'I');
    kstring_t line{0, 0, nullptr};
    kputs(sam.c_str(), &line);
    const int parsed = sam_parse1(&line, header.get(), read.alignment.get());
    std::free(line.s);
    REQUIRE(parsed >= 0);
    chunk.reads.push_back(std::move(read));
    ReadVariantProfile profile;
    profile.read_id = 0;
    chunk.read_var_profile.push_back(profile);
    return chunk;
}

} // namespace

TEST_CASE("equivalent deletion calls preserve separate allele rows",
          "[msa][recovery][representation][deletion]") {
    // Both placements delete the same A. The exact-position caller sees REF
    // at 103, but a shifted ALT must never vote REF for the two-base row.
    char path[] = "/tmp/pgphase-deletion-XXXXXX";
    const int fd = mkstemp(path);
    REQUIRE(fd >= 0);
    struct TemporaryReference {
        std::string path;
        ~TemporaryReference() {
            std::remove(path.c_str());
            std::remove((path + ".fai").c_str());
        }
    } temporary{path};
    const auto close_file = [](FILE* stream) { std::fclose(stream); };
    std::unique_ptr<FILE, decltype(close_file)> file(fdopen(fd, "w"), close_file);
    REQUIRE(file != nullptr);
    const std::string bases = std::string(99, 'N') + "CAAAAAGC" + std::string(93, 'N');
    REQUIRE(std::fprintf(file.get(), ">chr\n%s\n", bases.c_str()) > 0);
    file.reset();
    std::unique_ptr<faidx_t, decltype(&fai_destroy)> fai(fai_load(path), &fai_destroy);
    REQUIRE(fai != nullptr);
    ReferenceCache reference(fai.get());
    const std::string header_text = "@SQ\tSN:chr\tLN:200\n";
    std::unique_ptr<bam_hdr_t, decltype(&bam_hdr_destroy)> header(
        sam_hdr_parse(header_text.size(), header_text.c_str()), &bam_hdr_destroy);
    REQUIRE(header != nullptr);
    PhasingChunk chunk = shifted_insertion_chunk("5M1D2M", "CAAAAGC");
    CandidateVariant deletion;
    deletion.key.type = VariantType::Deletion;
    deletion.key.pos = 103;
    deletion.key.ref_len = 1;
    const auto call = [&](const CandidateVariant& candidate) {
        return bam_equivalent_deletion_allele(
            chunk.reads[0].alignment.get(), candidate, reference, 0, header.get(), 30);
    };

    SECTION("the shifted ALT reaches its original row") {
        int qi = -1;
        CHECK(bam_exact_indel_allele(chunk.reads[0].alignment.get(), deletion, 30, &qi) == 0);
        CHECK(call(deletion) == 1);
        CHECK(deletion.key.pos == 103);
        CHECK(deletion.key.ref_len == 1);
    }
    SECTION("a different deletion length is unknown, not REF") {
        deletion.key.ref_len = 2;
        CHECK(call(deletion) == -1);
    }
    SECTION("longer repeat deletions cannot impersonate either boundary ALT") {
        for (const auto& read : std::array<std::pair<std::string, std::string>, 2>{
                 std::make_pair("3M3D2M", "CAAGC"),
                 std::make_pair("2M4D2M", "CAGC")}) {
            chunk = shifted_insertion_chunk(read.first, read.second);
            deletion.key.ref_len = 1;
            CHECK(call(deletion) == -1);
            deletion.key.ref_len = 2;
            CHECK(call(deletion) == -1);
        }
    }
    SECTION("a two-base shifted ALT stays distinct from the one-base row") {
        chunk = shifted_insertion_chunk("4M2D2M", "CAAAGC");
        CHECK(call(deletion) == -1);
        deletion.key.ref_len = 2;
        CHECK(call(deletion) == 1);
    }
    SECTION("a left-shifted equivalent ALT is callable") {
        chunk = shifted_insertion_chunk("2M1D5M", "CAAAAGC");
        CHECK(call(deletion) == 1);
    }
    SECTION("an exact ALT is callable") {
        chunk = shifted_insertion_chunk("3M1D4M", "CAAAAGC");
        CHECK(call(deletion) == 1);
    }
    SECTION("clean reference is callable") {
        chunk = shifted_insertion_chunk("8M", "CAAAAAGC");
        CHECK(call(deletion) == 0);
    }
    SECTION("another edit inside the verified span is ambiguous") {
        chunk = shifted_insertion_chunk("3M1I2M1D2M", "CAATAAGC");
        CHECK(call(deletion) == -1);
    }
    SECTION("an unrelated edit outside the verified span does not hide ALT") {
        chunk = shifted_insertion_chunk("1M1I4M1D2M", "CTAAAAGC");
        CHECK(call(deletion) == 1);
    }
    SECTION("flank mismatches and missing or low qualities are not evidence") {
        chunk = shifted_insertion_chunk("5M1D2M", "CAACAGC");
        CHECK(call(deletion) == -1);
        chunk = shifted_insertion_chunk("5M1D2M", "CAAAAGC");
        bam_get_qual(chunk.reads[0].alignment.get())[4] = 29;
        CHECK(call(deletion) == -1);
        bam_get_qual(chunk.reads[0].alignment.get())[4] = 255;
        CHECK(call(deletion) == -1);
    }
    SECTION("the complete verification span must lie on the read") {
        chunk = shifted_insertion_chunk("6M", "CAAAAA");
        deletion.key.pos = 105;
        CHECK(call(deletion) == -1);
    }
}

TEST_CASE("recovery fills missing shifted deletion ALT observations",
          "[msa][recovery][representation][deletion-backfill]") {
    const auto make_chunk = [](const std::string& cigar,
                               const std::string& sequence) {
        PhasingChunk chunk = shifted_insertion_chunk(cigar + "100M", sequence + std::string(100, 'A'));
        chunk.ref_seq = "CAAAAAGC" + std::string(100, 'A');
        chunk.ref_end = chunk.ref_beg + static_cast<hts_pos_t>(chunk.ref_seq.size()) - 1;
        chunk.reads[0].end += 100;
        chunk.candidates[0].key = del_key(103, 1);
        chunk.candidates[0].phase_set = 100;
        chunk.candidates[0].hap_to_cons_alle[1] = 1;
        chunk.candidates[0].hap_to_cons_alle[2] = 0;
        CandidateVariant snp;
        snp.key.type = VariantType::Snp;
        snp.key.pos = 203;
        snp.key.ref_len = 1;
        snp.key.alt = "G";
        snp.ref_base = 0;
        snp.counts.category = VariantCategory::CleanHetSnp;
        snp.phase_set = 100;
        snp.hap_to_cons_alle[1] = 0;
        snp.hap_to_cons_alle[2] = 1;
        chunk.candidates.push_back(snp);
        return chunk;
    };
    Options opts;
    opts.min_bq = 30;
    const auto fill = [&opts](PhasingChunk& chunk) {
        return backfill_msa_observations(chunk, opts, 102, 102);
    };
    SECTION("an uncallable shifted ALT reaches its original row") {
        for (const std::string cigar : {"4M1D3M", "2M1D5M"}) {
            PhasingChunk chunk = make_chunk(cigar, "CAAAAGC");
            REQUIRE(fill(chunk) == 1);
            CHECK(chunk.read_var_profile[0].alleles == std::vector<int>{1});
            CHECK(chunk.candidates[0].key.pos == 103);
            CHECK(chunk.candidates[0].key.ref_len == 1);
            CHECK(chunk.read_var_cr != nullptr);
            REQUIRE(chunk.read_var_profile[0].alt_qi.size() == 1);
            CHECK(chunk.read_var_profile[0].alt_qi[0] >= 0);
            CHECK(fill(chunk) == 0);
        }
    }
    SECTION("a different deletion length stays unknown") {
        PhasingChunk chunk = make_chunk("4M1D3M", "CAAAAGC");
        chunk.candidates[0].key.ref_len = 2;
        CHECK(fill(chunk) == 0);
        CHECK(chunk.read_var_profile[0].start_var_idx == -1);
    }
    SECTION("callable exact contrasts remain unchanged") {
        PhasingChunk chunk = make_chunk("5M1D2M", "CAAAAGC");
        REQUIRE(fill(chunk) == 1);
        CHECK(chunk.read_var_profile[0].alleles == std::vector<int>{0});
        chunk = make_chunk("4M1I4M", "CAAAAAAGC");
        chunk.candidates[0].key.ref_len = 3;
        REQUIRE(fill(chunk) == 1);
        CHECK(chunk.read_var_profile[0].alleles == std::vector<int>{0});
    }
    SECTION("colocated MSA alternatives retain their binary source contrast") {
        PhasingChunk chunk = make_chunk("4M1I4M", "CAAAAAAGC");
        chunk.candidates[0].key.ref_len = 3;
        CandidateVariant alternative = chunk.candidates[0];
        alternative.key = ins_key(103, "A");
        chunk.candidates.push_back(alternative);
        REQUIRE(fill(chunk) >= 1);
        CHECK(chunk.read_var_profile[0].alleles[0] == 0);
        chunk = make_chunk("5M1D2M", "CAAAAGC");
        chunk.candidates[0].key.ref_len = 2;
        alternative = chunk.candidates[0];
        alternative.key.ref_len = 1;
        chunk.candidates.push_back(alternative);
        REQUIRE(fill(chunk) >= 1);
        CHECK(chunk.read_var_profile[0].alleles[0] == 0);
    }
    SECTION("a missing multi-row MSA contrast is not repaired in isolation") {
        PhasingChunk chunk = make_chunk("4M1D3M", "CAAAAGC");
        CandidateVariant alternative = chunk.candidates[0];
        alternative.key = ins_key(103, "A");
        chunk.candidates.push_back(alternative);
        fill(chunk);
        const auto& profile = chunk.read_var_profile[0];
        const int deletion_call = profile.start_var_idx == 0
            ? profile.alleles[0] : -1;
        CHECK(deletion_call == -1);
    }
    SECTION("clean REF and exact ALT remain callable") {
        for (const bool alt : {false, true}) {
            PhasingChunk chunk = make_chunk(alt ? "3M1D4M" : "8M",
                                             alt ? "CAAAAGC" : "CAAAAAGC");
            REQUIRE(fill(chunk) == 1);
            CHECK(chunk.read_var_profile[0].alleles ==
                  std::vector<int>{alt ? 1 : 0});
        }
    }
    SECTION("an independently callable SNP must preserve the source gauge") {
        PhasingChunk chunk = make_chunk("4M1D3M", "CAAAAGC");
        chunk.candidates[1].phase_set = 200;
        CHECK(fill(chunk) == 0);
        chunk = make_chunk("4M1D3M", "CAAAAGC");
        chunk.candidates[1].hap_to_cons_alle[1] = 1;
        chunk.candidates[1].hap_to_cons_alle[2] = 0;
        CHECK(fill(chunk) == 0);
        chunk = make_chunk("4M1D3M", "CAAAAGC");
        for (const int quality : {29, 255}) {
            bam_get_qual(chunk.reads[0].alignment.get())[102] = quality;
            CHECK(fill(chunk) == 0);
        }
    }
    SECTION("unknown and low mapping qualities do not admit a new ALT") {
        opts.min_mapq = 0;
        for (const int mapq : {29, 255}) {
            PhasingChunk chunk = make_chunk("4M1D3M", "CAAAAGC");
            chunk.reads[0].mapq = mapq;
            chunk.reads[0].alignment->core.qual = mapq;
            CHECK(fill(chunk) == 0);
            CHECK(chunk.read_var_profile[0].start_var_idx == -1);
        }
    }
    SECTION("missing and low-quality shifted flanks are not evidence") {
        opts.min_bq = 1;
        for (const int quality : {29, 255}) {
            PhasingChunk chunk = make_chunk("4M1D3M", "CAAAAGC");
            bam_get_qual(chunk.reads[0].alignment.get())[4] = quality;
            CHECK(fill(chunk) == 0);
            CHECK(chunk.read_var_profile[0].start_var_idx == -1);
        }
    }
    SECTION("existing MSA allele decisions are preserved") {
        for (const int allele : {0, 1, -2}) {
            PhasingChunk chunk = make_chunk("4M1D3M", "CAAAAGC");
            auto& profile = chunk.read_var_profile[0];
            profile.start_var_idx = profile.end_var_idx = 0;
            profile.alleles = {allele};
            profile.alt_qi = {-1};
            CHECK(fill(chunk) == 0);
            CHECK(profile.alleles == std::vector<int>{allele});
        }
    }
    SECTION("an independent flank substitution does not hide primitive REF") {
        PhasingChunk chunk = make_chunk("8M", "CACAAAGC");
        REQUIRE(fill(chunk) == 1);
        CHECK(chunk.read_var_profile[0].alleles == std::vector<int>{0});
    }
}

TEST_CASE("targeted recovery calls a shifted single-base insertion before phasing",
          "[msa][recovery][representation]") {
    PhasingChunk chunk = shifted_insertion_chunk();
    Options opts;
    opts.retry_windows = {{102, 102}};
    SECTION("the equivalent ALT creates an indexed sparse profile") {
        int qi = -1;
        // The exact-position caller calls this shifted event REF. Recovery
        // must recognize its reference-edit-equivalent ALT before solving.
        CHECK(bam_exact_indel_allele(chunk.reads[0].alignment.get(),
                                    chunk.candidates[0], 30, &qi) == 0);
        REQUIRE(backfill_shifted_msa_insertions(chunk, opts) == 1);
        const auto& profile = chunk.read_var_profile[0];
        CHECK(profile.start_var_idx == 0);
        CHECK(profile.alleles == std::vector<int>{1});
        CHECK(profile.alt_qi == std::vector<int>{6});
        REQUIRE(chunk.read_var_cr != nullptr);
        int64_t* hits = nullptr;
        int64_t capacity = 0;
        CHECK(cr_overlap(chunk.read_var_cr.get(), "cr", 0, 1, &hits, &capacity) == 1);
        std::free(hits);
        CHECK(backfill_shifted_msa_insertions(chunk, opts) == 0);
    }
    SECTION("a left-shifted event represents the same ALT") {
        chunk = shifted_insertion_chunk("2M1I6M");
        CHECK(backfill_shifted_msa_insertions(chunk, opts) == 1);
    }
    SECTION("existing MSA calls retain their own allele") {
        auto& profile = chunk.read_var_profile[0];
        profile.start_var_idx = profile.end_var_idx = 0;
        profile.alleles = {0};
        profile.alt_qi = {-1};
        CHECK(backfill_shifted_msa_insertions(chunk, opts) == 0);
        CHECK(profile.alleles == std::vector<int>{0});
    }
    SECTION("ordinary BAM and sites outside recovery remain unchanged") {
        opts.retry_windows.clear();
        CHECK(backfill_shifted_msa_insertions(chunk, opts) == 0);
        opts.retry_windows = {{103, 107}};
        CHECK(backfill_shifted_msa_insertions(chunk, opts) == 0);
        opts.retry_windows = {{100, 101}};
        CHECK(backfill_shifted_msa_insertions(chunk, opts) == 0);
    }
    SECTION("unknown and low mapping qualities cannot supply a bridge") {
        for (const int mapq : {29, 255}) {
            chunk.reads[0].mapq = mapq;
            CHECK(backfill_shifted_msa_insertions(chunk, opts) == 0);
        }
    }
    SECTION("inserted bases and flanks require known Q30 qualities") {
        for (const int qi : {3, 6}) {
            for (const int quality : {29, 255}) {
                bam_get_qual(chunk.reads[0].alignment.get())[qi] = quality;
                CHECK(backfill_shifted_msa_insertions(chunk, opts) == 0);
            }
            bam_get_qual(chunk.reads[0].alignment.get())[qi] = 40;
        }
    }
    SECTION("a different event is never converted to this ALT or REF") {
        chunk = shifted_insertion_chunk("6M1I2M", "CATTTTAGC");
        CHECK(backfill_shifted_msa_insertions(chunk, opts) == 0);
        chunk = shifted_insertion_chunk("6M2I2M", "CATTTTTTGC");
        CHECK(backfill_shifted_msa_insertions(chunk, opts) == 0);
        chunk = shifted_insertion_chunk("3M1I3M1I2M", "CATTTTTTGC");
        CHECK(backfill_shifted_msa_insertions(chunk, opts) == 0);
    }
    SECTION("a skipped reference base breaks equivalence") {
        chunk.ref_seq[4] = 'A';
        CHECK(backfill_shifted_msa_insertions(chunk, opts) == 0);
    }
    SECTION("an incomplete or ambiguous reference cannot certify equivalence") {
        chunk.ref_seq[4] = 'N';
        CHECK(backfill_shifted_msa_insertions(chunk, opts) == 0);
        chunk.ref_seq = "CATTTT";
        CHECK(backfill_shifted_msa_insertions(chunk, opts) == 0);
    }
    SECTION("an exact event and complex candidates keep existing handling") {
        chunk = shifted_insertion_chunk("3M1I5M");
        CHECK(backfill_shifted_msa_insertions(chunk, opts) == 0);
        chunk = shifted_insertion_chunk();
        chunk.candidates[0].key.alt = "TT";
        CHECK(backfill_shifted_msa_insertions(chunk, opts) == 0);
        chunk.candidates[0].key.alt = "T";
        chunk.candidates[0].msa_insertion_alts = {"T", "TT"};
        CHECK(backfill_shifted_msa_insertions(chunk, opts) == 0);
    }
}

TEST_CASE("MSA deletion rows retain a verified alternate-absent insertion allele",
          "[msa][recovery][complementary-indels]") {
    const auto alignment = [](const std::string& target, const std::string& query) {
        REQUIRE(target.size() == query.size());
        AlnStr result;
        for (char base : target) result.target_aln.push_back(base == '-' ? 5 : base_to_nt4(base));
        for (char base : query) result.query_aln.push_back(base == '-' ? 5 : base_to_nt4(base));
        result.aln_len = static_cast<int>(target.size());
        result.target_end = result.query_end = result.aln_len - 1;
        return result;
    };
    const AlnStr insertion = alignment("ACG-TGCA", "ACGTTGCA");
    const AlnStr deletion = alignment("ACGTGCA", "ACG-GCA");
    const VariantKey key = del_key(103, 1);
    std::array<AlnStr, 2> consensuses{insertion, deletion};
    SECTION("the other exact haplotype is zero for this deletion row") {
        CHECK(call_msa_site_allele({insertion, insertion}, key, 100, &consensuses) == 0);
        CHECK(call_msa_site_allele({deletion, deletion}, key, 100, &consensuses) == 1);
        CHECK(call_msa_site_allele({insertion, insertion}, key, 100) == -1);
    }
    SECTION("absent target or conflicting compositions do not invent a contrast") {
        consensuses = {insertion, insertion};
        CHECK(call_msa_site_allele({insertion, insertion}, key, 100, &consensuses) == -1);
        consensuses = {insertion, deletion};
        CHECK(call_msa_site_allele({insertion, deletion}, key, 100, &consensuses) == -1);
        const AlnStr different = alignment("ACG-TGCA", "ACGTAGCA");
        CHECK(call_msa_site_allele({different, different}, key, 100, &consensuses) == -1);
    }
    SECTION("partial or ambiguous context cannot supply the alternate-absent allele") {
        consensuses[1].target_beg = 1;
        CHECK(call_msa_site_allele({insertion, insertion}, key, 100, &consensuses) == -1);
        consensuses = {insertion, deletion};
        consensuses[0].query_aln[4] = 4;
        CHECK(call_msa_site_allele({insertion, insertion}, key, 100, &consensuses) == -1);
    }
    SECTION("a different deletion length retains its own representation") {
        const VariantKey longer = del_key(103, 2);
        const AlnStr partial = alignment("ACGTTGCA", "ACG-TGCA");
        const AlnStr complete = alignment("ACGTTGCA", "ACG--GCA");
        const std::array<AlnStr, 2> deletion_context{partial, complete};
        CHECK(call_msa_site_allele({partial, partial}, longer, 100,
                                   &deletion_context) == -1);
        CHECK(call_msa_site_allele({complete, complete}, longer, 100,
                                   &deletion_context) == 1);
    }
    SECTION("a mixed edit requires exact context instead of enabling fuzzy calls") {
        CandidateVariant candidate;
        candidate.key = key;
        candidate.counts.category = VariantCategory::NoisyCandHet;
        candidate.counts.alle_covs = {1, 1};
        candidate.counts.total_cov = 2;
        std::vector<CandidateVariant> candidates{candidate};
        std::vector<ReadVariantProfile> profiles(1);
        UnassignedMsaRead read;
        read.read_id = 0;
        const AlnStr mismatch = alignment("ACG-TGCA", "ACGTTGTA");
        read.ref_read = {mismatch, mismatch};
        Options opts;
        add_msa_site_observations(opts, {read}, 100, candidates, profiles, &consensuses);
        CHECK(profiles[0].start_var_idx == -1);
        CHECK(candidates[0].counts.alle_covs == std::vector<int>{1, 1});
    }
    SECTION("local observation recovery preserves the separate rows and counts") {
        CandidateVariant insertion_row;
        insertion_row.key = ins_key(103, "T");
        insertion_row.counts.category = VariantCategory::NoisyCandHet;
        insertion_row.counts.alle_covs = {1, 1};
        insertion_row.counts.total_cov = 2;
        CandidateVariant deletion_row = insertion_row;
        deletion_row.key = key;
        std::vector<CandidateVariant> candidates{insertion_row, deletion_row};
        std::vector<ReadVariantProfile> profiles(3);
        profiles[0].start_var_idx = profiles[1].start_var_idx = 0;
        profiles[0].end_var_idx = profiles[1].end_var_idx = 1;
        profiles[0].alleles = {1, 0};
        profiles[1].alleles = {0, 1};
        UnassignedMsaRead read;
        read.read_id = 2;
        read.ref_read = {insertion, insertion};
        Options opts;
        add_msa_site_observations(opts, {read}, 100, candidates, profiles, &consensuses);
        CHECK(profiles[2].alleles == std::vector<int>{1, 0});
        CHECK(candidates.size() == 2);
        CHECK(candidates[0].counts.alle_covs == std::vector<int>{1, 2});
        CHECK(candidates[1].counts.alle_covs == std::vector<int>{2, 1});
        CHECK(candidates[0].key.type == VariantType::Insertion);
        CHECK(candidates[0].key.alt == "T");
        CHECK(candidates[1].key.pos == key.pos);
        CHECK(candidates[1].key.type == VariantType::Deletion);
        CHECK(candidates[1].key.ref_len == 1);
    }
}

TEST_CASE("repeat insertion calls retain equivalent ALT beyond the old search bound",
          "[msa][recovery][representation][insertion]") {
    std::string reference = "G";
    for (int i = 0; i < 20; ++i) reference += "ATCT";
    reference += "C";
    const auto make_chunk = [&reference](int offset = 45, const std::string& inserted = "ATCT") {
        const std::string sequence = reference.substr(0, offset) + inserted + reference.substr(offset);
        PhasingChunk chunk = shifted_insertion_chunk(
            std::to_string(offset) + "M" + std::to_string(inserted.size()) + "I" +
            std::to_string(reference.size() - offset) + "M", sequence);
        chunk.ref_seq = reference;
        chunk.ref_end = chunk.ref_beg + reference.size() - 1;
        chunk.reads[0].end = bam_endpos(chunk.reads[0].alignment.get());
        chunk.candidates[0].key = ins_key(101, "ATCT");
        chunk.candidates[0].phase_set = 101;
        chunk.candidates[0].hap_to_cons_alle[1] = 1;
        chunk.candidates[0].hap_to_cons_alle[2] = 0;
        chunk.phase_sets = {101};
        chunk.haps = {2};
        return chunk;
    };
    PhasingChunk chunk = make_chunk();
    const auto call = [](const PhasingChunk& target) {
        return bam_shifted_repeat_insertion_query_index(
            target.reads[0].alignment.get(), target.candidates[0], target, 30);
    };
    SECTION("a 44-base shift calls ALT without changing the source gauge or row") {
        CHECK(insertion_equivalent_positions(101, "ATCT", chunk) ==
              std::make_pair(hts_pos_t{101}, hts_pos_t{181}));
        int qi = -1;
        CHECK(bam_exact_indel_allele(chunk.reads[0].alignment.get(),
                                    chunk.candidates[0], 30, &qi) == 0);
        CHECK(call(chunk) == 45);
        CHECK(chunk.candidates[0].key.pos == 101);
        CHECK(chunk.candidates[0].key.alt == "ATCT");
        CHECK(chunk.candidates[0].phase_set == 101);
        CHECK(chunk.candidates[0].hap_to_cons_alle[1] == 1);
        CHECK(chunk.candidates[0].hap_to_cons_alle[2] == 0);
        CHECK(chunk.phase_sets == std::vector<hts_pos_t>{101});
        CHECK(chunk.haps == std::vector<int>{2});
        CHECK(chunk.read_var_profile[0].alleles.empty());
    }
    SECTION("motif rotations and shifts in either direction preserve the edit") {
        chunk = make_chunk(44, "TATC");
        CHECK(call(chunk) == 44);
        chunk = make_chunk(1);
        chunk.candidates[0].key = ins_key(145, "ATCT");
        CHECK(call(chunk) == 1);
    }
    SECTION("all inserted and crossed reference bases require known Q30") {
        for (const int qi : {0, 22, 45, 48, 49}) {
            for (const int quality : {29, 255}) {
                chunk = make_chunk();
                bam_get_qual(chunk.reads[0].alignment.get())[qi] = quality;
                CHECK(call(chunk) == -1);
            }
        }
    }
    SECTION("different insertion sequences and lengths are not this ALT") {
        chunk = make_chunk(45, "ATCA");
        CHECK(call(chunk) == -1);
        chunk = make_chunk(45, "ATCTATCT");
        CHECK(call(chunk) == -1);
        chunk = make_chunk();
        chunk.candidates[0].key.alt = "ATCN";
        CHECK(call(chunk) == -1);
    }
    SECTION("a reference mismatch or ambiguity breaks the equivalence certificate") {
        for (const char base : {'C', 'N', 'M'}) {
            chunk = make_chunk();
            chunk.ref_seq[22] = base;
            CHECK(call(chunk) == -1);
        }
        chunk = make_chunk();
        chunk.ref_seq.resize(30);
        CHECK(call(chunk) == -1);
        chunk = make_chunk();
        for (char& base : chunk.ref_seq)
            base = static_cast<char>(std::tolower(static_cast<unsigned char>(base)));
        CHECK(call(chunk) == 45);
    }
    SECTION("reference boundaries and empty alleles terminate the interval") {
        CHECK(insertion_equivalent_positions(181, "ATCT", chunk) ==
              std::make_pair(hts_pos_t{101}, hts_pos_t{181}));
        CHECK(insertion_equivalent_positions(101, "", chunk) ==
              std::make_pair(hts_pos_t{101}, hts_pos_t{101}));
        chunk.ref_seq = "ATCTATCT";
        chunk.ref_beg = 1;
        CHECK(insertion_equivalent_positions(1, "ATCT", chunk) ==
              std::make_pair(hts_pos_t{1}, hts_pos_t{9}));
    }
    SECTION("compound insertion and deletion paths are not certified") {
        chunk = shifted_insertion_chunk("13M4I32M4I37M",
            reference.substr(0, 13) + "ATCT" + reference.substr(13, 32) +
            "ATCT" + reference.substr(45));
        chunk.ref_seq = reference;
        chunk.candidates[0].key = ins_key(101, "ATCT");
        CHECK(call(chunk) == -1);
        chunk = shifted_insertion_chunk("25M1D19M4I37M",
            reference.substr(0, 25) + reference.substr(26, 19) +
            "ATCT" + reference.substr(45));
        chunk.ref_seq = reference;
        chunk.candidates[0].key = ins_key(101, "ATCT");
        CHECK(call(chunk) == -1);
    }
    SECTION("the worker reference cache calls the same edit as the local fixture") {
        char path[] = "/tmp/pgphase-insertion-XXXXXX";
        const int fd = mkstemp(path);
        REQUIRE(fd >= 0);
        struct TemporaryReference {
            std::string path;
            ~TemporaryReference() {
                std::remove(path.c_str());
                std::remove((path + ".fai").c_str());
            }
        } temporary{path};
        const auto close_file = [](FILE* stream) { std::fclose(stream); };
        std::unique_ptr<FILE, decltype(close_file)> file(fdopen(fd, "w"), close_file);
        REQUIRE(file != nullptr);
        const std::string bases = std::string(99, 'N') + reference;
        REQUIRE(std::fprintf(file.get(), ">chr\n%s\n", bases.c_str()) > 0);
        file.reset();
        std::unique_ptr<faidx_t, decltype(&fai_destroy)> fai(fai_load(path), &fai_destroy);
        REQUIRE(fai != nullptr);
        ReferenceCache cache(fai.get());
        const std::string header_text = "@SQ\tSN:chr\tLN:181\n";
        std::unique_ptr<bam_hdr_t, decltype(&bam_hdr_destroy)> header(
            sam_hdr_parse(header_text.size(), header_text.c_str()), &bam_hdr_destroy);
        REQUIRE(header != nullptr);
        CHECK(insertion_equivalent_positions(101, "ATCT", cache, 0, header.get()) ==
              insertion_equivalent_positions(101, "ATCT", chunk));
        CHECK(bam_shifted_repeat_insertion_query_index(chunk.reads[0].alignment.get(),
            chunk.candidates[0], cache, 0, header.get(), 30) == call(chunk));
    }
    SECTION("exact placements retain their existing caller") {
        chunk = make_chunk(1);
        CHECK(call(chunk) == -1);
    }
}

TEST_CASE("postsolve MSA backfill includes indel VCF boundary anchors",
          "[msa][recovery][coordinates]") {
    Options opts;
    SECTION("an insertion at a singleton seam retains its exact ALT") {
        PhasingChunk chunk = shifted_insertion_chunk("3M1I5M");
        REQUIRE(chunk.candidates[0].key.sort_pos() == 102);
        REQUIRE(backfill_msa_observations(chunk, opts, 102, 102) == 1);
        CHECK(chunk.candidates[0].key.pos == 103);
        CHECK(chunk.read_var_profile[0].alleles == std::vector<int>{1});
        CHECK(chunk.read_var_profile[0].alt_qi == std::vector<int>{3});
        CHECK(chunk.read_var_cr != nullptr);
        CHECK(backfill_msa_observations(chunk, opts, 102, 102) == 0);
    }
    SECTION("a deletion at the right boundary retains its exact ALT") {
        PhasingChunk chunk = shifted_insertion_chunk("3M1D4M", "CATTTGC");
        chunk.candidates[0].key = del_key(103, 1);
        REQUIRE(chunk.candidates[0].key.sort_pos() == 102);
        REQUIRE(backfill_msa_observations(chunk, opts, 100, 102) == 1);
        CHECK(chunk.candidates[0].key.pos == 103);
        CHECK(chunk.read_var_profile[0].alleles == std::vector<int>{1});
        CHECK(chunk.read_var_cr != nullptr);
    }
    SECTION("reference observations also use the VCF anchor") {
        for (const bool deletion : {false, true}) {
            PhasingChunk chunk = shifted_insertion_chunk("8M", "CATTTTGC");
            if (deletion) chunk.candidates[0].key = del_key(103, 1);
            REQUIRE(backfill_msa_observations(chunk, opts, 102, 102) == 1);
            CHECK(chunk.read_var_profile[0].alleles == std::vector<int>{0});
        }
    }
    SECTION("the internal indel base cannot admit an out-of-window anchor") {
        PhasingChunk chunk = shifted_insertion_chunk("3M1I5M");
        CHECK(backfill_msa_observations(chunk, opts, 103, 107) == 0);
        CHECK(backfill_msa_observations(chunk, opts, 100, 101) == 0);
        CHECK(chunk.read_var_profile[0].start_var_idx == -1);
    }
    SECTION("SNP coordinates have no anchor offset") {
        PhasingChunk chunk = shifted_insertion_chunk("8M", "CATTATGC");
        chunk.candidates[0].key.type = VariantType::Snp;
        chunk.candidates[0].key.ref_len = 1;
        chunk.candidates[0].key.alt = "T";
        chunk.candidates[0].ref_base = base_to_nt4('A');
        CHECK(backfill_msa_observations(chunk, opts, 102, 102) == 0);
        REQUIRE(backfill_msa_observations(chunk, opts, 103, 103) == 1);
        CHECK(chunk.read_var_profile[0].alleles == std::vector<int>{1});
    }
    SECTION("existing MSA calls and ambiguity are preserved") {
        for (const int allele : {0, 1, -2}) {
            PhasingChunk chunk = shifted_insertion_chunk("3M1I5M");
            auto& profile = chunk.read_var_profile[0];
            profile.start_var_idx = profile.end_var_idx = 0;
            profile.alleles = {allele};
            profile.alt_qi = {-1};
            CHECK(backfill_msa_observations(chunk, opts, 102, 102) == 0);
            CHECK(profile.alleles == std::vector<int>{allele});
        }
    }
    SECTION("boundary admission does not bypass the allele quality check") {
        PhasingChunk chunk = shifted_insertion_chunk("3M1I5M");
        bam_get_qual(chunk.reads[0].alignment.get())[3] = 0;
        CHECK(backfill_msa_observations(chunk, opts, 102, 102) == 0);
        CHECK(chunk.read_var_profile[0].start_var_idx == -1);
    }
}

TEST_CASE("sparse recovery growth preserves the quality of each BAM site",
          "[recovery][profile][quality]") {
    ReadVariantProfile profile;
    profile.start_var_idx = 5;
    profile.end_var_idx = 6;
    profile.alleles = {0, 1};
    profile.alt_qi = {10, 11};
    profile.graph_alleles = {1, 0};
    profile.bam_alleles = {0, 1};
    profile.bam_qi = {10, 11};
    profile.bam_base_qualities = {40, 5};
    profile.bam_mapq = 60;
    SECTION("prepending recovery sites must not move Q40 onto another site") {
        update_read_var_profile_with_allele(3, 1, 8, profile);
        CHECK(profile.bam_base_qualities == std::vector<uint8_t>{0, 0, 40, 5});
        CHECK(profile.bam_alleles == std::vector<int>{-1, -1, 0, 1});
        CHECK(profile.graph_alleles == std::vector<int>{-1, -1, 1, 0});
        CHECK(profile.bam_qi == std::vector<int>{-1, -1, 10, 11});
        CHECK(profile.bam_mapq == 60);
    }
    SECTION("appending sites gives them no unmeasured base quality") {
        update_read_var_profile_with_allele(9, 0, 14, profile);
        CHECK(profile.bam_base_qualities == std::vector<uint8_t>{40, 5, 0, 0, 0});
        CHECK(profile.bam_alleles == std::vector<int>{0, 1, -1, -1, -1});
    }
    SECTION("growth on both sides preserves the original site offsets") {
        update_read_var_profile_with_allele(3, 1, 8, profile);
        update_read_var_profile_with_allele(9, 0, 14, profile);
        REQUIRE(profile.bam_base_qualities.size() == profile.alleles.size());
        CHECK(profile.bam_base_qualities[5 - profile.start_var_idx] == 40);
        CHECK(profile.bam_base_qualities[6 - profile.start_var_idx] == 5);
        CHECK(profile.bam_base_qualities.front() == 0);
        CHECK(profile.bam_base_qualities.back() == 0);
    }
    SECTION("profiles without BAM qualities keep that channel absent") {
        profile.bam_base_qualities.clear();
        update_read_var_profile_with_allele(3, 1, 8, profile);
        update_read_var_profile_with_allele(9, 0, 14, profile);
        CHECK(profile.bam_base_qualities.empty());
    }
    SECTION("an existing site retains its measured quality") {
        update_read_var_profile_with_allele(5, 1, 10, profile);
        CHECK(profile.bam_base_qualities == std::vector<uint8_t>{40, 5});
    }
}

TEST_CASE("BAM block corroboration uses any independent agreeing SNP",
          "[recovery][stitch][corroboration]") {
    PhasingChunk chunk;
    for (const hts_pos_t pos : {1000, 1950, 2000, 3000}) {
        CandidateVariant candidate;
        candidate.key.pos = pos;
        candidate.key.type = VariantType::Snp;
        candidate.counts.category = VariantCategory::CleanHetSnp;
        candidate.hap_to_cons_alle = {-1, 0, 1};
        chunk.candidates.push_back(candidate);
    }
    chunk.candidates[2].key.type = VariantType::Insertion;
    chunk.candidates[2].key.alt = "A";
    chunk.candidates[2].counts.category = VariantCategory::NoisyCandHet;
    chunk.candidates[2].msa_verified = true;
    ReadRecord read;
    read.mapq = 60;
    chunk.reads.push_back(std::move(read));
    ReadVariantProfile profile;
    profile.start_var_idx = 0;
    profile.end_var_idx = 3;
    profile.alleles = {0, 0, 0, 0};
    profile.bam_alleles = profile.alleles;
    profile.bam_base_qualities = {40, 40, 0, 40};
    profile.bam_mapq = 60;
    chunk.read_var_profile.push_back(profile);
    SECTION("a later nearby SNP cannot hide an earlier independent SNP") {
        const auto flip = corroborated_bam_block_flip(chunk, {0, 1, 2}, {3});
        REQUIRE(flip.has_value());
        CHECK_FALSE(*flip);
        CHECK(corroborated_bam_block_flip(chunk, {2, 1, 0}, {3}) == flip);
    }
    SECTION("a low-quality distant SNP does not provide corroboration") {
        chunk.read_var_profile[0].bam_base_qualities[0] = 5;
        CHECK_FALSE(corroborated_bam_block_flip(chunk, {0, 1, 2}, {3}));
    }
    SECTION("an opposing SNP remains a veto even beside an agreeing pair") {
        chunk.read_var_profile[0].bam_alleles[1] = 1;
        CHECK_FALSE(corroborated_bam_block_flip(chunk, {0, 1, 2}, {3}));
    }
    SECTION("the right block can supply the independent indel") {
        const auto flip = corroborated_bam_block_flip(chunk, {3}, {0, 1, 2});
        REQUIRE(flip.has_value());
        CHECK_FALSE(*flip);
    }
    SECTION("an earlier indel also retains its independent corroboration") {
        chunk.read_var_profile[0].bam_base_qualities[0] = 5;
        CandidateVariant candidate = chunk.candidates[2];
        candidate.key.pos = 1000;
        chunk.candidates.push_back(candidate);
        chunk.read_var_profile[0].end_var_idx = 4;
        chunk.read_var_profile[0].alleles.push_back(0);
        chunk.read_var_profile[0].bam_alleles.push_back(0);
        chunk.read_var_profile[0].bam_base_qualities.push_back(0);
        const auto flip = corroborated_bam_block_flip(chunk, {4, 1, 2}, {3});
        REQUIRE(flip.has_value());
        CHECK_FALSE(*flip);
    }
    SECTION("independent means at least the established hundred-base spacing") {
        chunk.read_var_profile[0].bam_base_qualities[0] = 5;
        // Insertion sort_pos is the VCF anchor one base before key.pos.
        chunk.candidates[2].key.pos = 2051;
        CHECK(corroborated_bam_block_flip(chunk, {0, 1, 2}, {3}).has_value());
        chunk.candidates[2].key.pos = 2050;
        CHECK_FALSE(corroborated_bam_block_flip(chunk, {0, 1, 2}, {3}));
    }
    SECTION("the qualifying molecule determines the relative block flip") {
        chunk.read_var_profile[0].bam_alleles[3] = 1;
        const auto flip = corroborated_bam_block_flip(chunk, {0, 1, 2}, {3});
        REQUIRE(flip.has_value());
        CHECK(*flip);
    }
}

TEST_CASE("recovery MSA calls stay inside physical read coverage", "[msa-coverage]") {
    std::vector<ReadRecord> reads(2);
    reads[0].beg = 100;
    reads[0].end = 110;
    reads[1].beg = 105;
    reads[1].end = 115;
    std::vector<CandidateVariant> vars;
    const auto add = [&](VariantType type, hts_pos_t pos, int ref_len) {
        CandidateVariant var;
        var.key.type = type;
        var.key.pos = pos;
        var.key.ref_len = ref_len;
        var.counts.alle_covs = {20, 20};
        var.counts.total_cov = 40;
        vars.push_back(std::move(var));
    };
    add(VariantType::Snp, 99, 1);
    add(VariantType::Snp, 100, 1);
    add(VariantType::Snp, 110, 1);
    add(VariantType::Insertion, 100, 0);
    add(VariantType::Insertion, 101, 0);
    add(VariantType::Insertion, 110, 0);
    add(VariantType::Insertion, 111, 0);
    add(VariantType::Deletion, 100, 1);
    add(VariantType::Deletion, 101, 9);
    add(VariantType::Deletion, 101, 10);
    std::vector<ReadVariantProfile> profiles(2);
    for (int ri = 0; ri < 2; ++ri) {
        ReadVariantProfile& profile = profiles[ri];
        profile.read_id = ri;
        profile.start_var_idx = ri == 0 ? 0 : 3;
        profile.end_var_idx = 9;
        profile.alleles.assign(10 - profile.start_var_idx, 1 - ri);
        profile.alt_qi.assign(profile.alleles.size(), 5);
    }
    restrict_msa_observations_to_read_coverage(reads, vars, profiles);
    CHECK(profiles[0].alleles == std::vector<int>{-1, 1, 1, -1, 1, 1, -1, -1, 1, -1});
    CHECK(profiles[1].alleles == std::vector<int>{-1, -1, 0, 0, -1, -1, -1});
    CHECK(profiles[0].alt_qi[0] == -1);
    CHECK(profiles[0].alt_qi[1] == 5);
    CHECK(vars[0].counts.total_cov == 0);
    CHECK(vars[2].counts.total_cov == 1);
    CHECK(vars[5].counts.alle_covs == std::vector<int>{1, 1});
    CHECK(vars[5].counts.total_cov == 2);
    CHECK(vars[5].counts.ref_cov == 1);
    CHECK(vars[5].counts.alt_cov == 1);
    CHECK(vars[5].counts.allele_fraction == Approx(0.5));
    CHECK(vars[8].counts.total_cov == 1);
    CHECK(vars[9].counts.total_cov == 0);
}
