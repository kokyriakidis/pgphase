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

TEST_CASE("padded graph indels use the emitted physical BAM key",
          "[recovery][representation]") {
    const auto check_key = [](hts_pos_t pos, const std::string& ref,
                              const std::string& alt, VariantType type,
                              hts_pos_t expected_pos, int ref_len,
                              const std::string& expected_alt) {
        const VariantKey raw = vcf_to_variant_key(11, pos, ref, alt);
        CHECK(raw.tid == 11);
        CHECK(raw.type == type);
        CHECK(raw.pos == expected_pos);
        CHECK(raw.ref_len == ref_len);
        CHECK(raw.alt == expected_alt);
        std::string emitted_ref = ref;
        std::string emitted_alt = alt;
        trim_to_minimal_vcf(pos, emitted_ref, emitted_alt);
        const VariantKey emitted =
            vcf_to_variant_key(11, pos, emitted_ref, emitted_alt);
        CHECK(exact_comp_var_site(&raw, &emitted) == 0);
    };
    SECTION("common suffix is context rather than inserted or deleted bases") {
        check_key(100, "ACG", "ATCG", VariantType::Insertion, 101, 0, "T");
        check_key(100, "ATCG", "ACG", VariantType::Deletion, 101, 1, "");
    }
    SECTION("complex replacements retain their nonmatching reference and ALT") {
        check_key(100, "ACG", "ATTG", VariantType::Insertion, 101, 1, "TT");
        check_key(100, "ATTG", "ACG", VariantType::Deletion, 101, 2, "C");
    }
    SECTION("minimal indels and padded SNPs keep their existing convention") {
        check_key(100, "A", "AT", VariantType::Insertion, 101, 0, "T");
        check_key(100, "AT", "A", VariantType::Deletion, 101, 1, "");
        check_key(100, "ACG", "ATG", VariantType::Snp, 101, 1, "T");
    }
}

TEST_CASE("graph indel context trimming preserves BAM repeat placement",
          "[recovery][representation]") {
    const VariantKey deletion = vcf_to_variant_key(11, 100, "TAA", "TA");
    CHECK(deletion.type == VariantType::Deletion);
    CHECK(deletion.pos == 102);
    CHECK(deletion.ref_len == 1);
    CHECK(deletion.alt.empty());
    const VariantKey insertion = vcf_to_variant_key(11, 100, "TAA", "TAAA");
    CHECK(insertion.type == VariantType::Insertion);
    CHECK(insertion.pos == 103);
    CHECK(insertion.ref_len == 0);
    CHECK(insertion.alt == "A");
}

TEST_CASE("recovery anchors require a complete heterozygous genotype",
          "[recovery-anchor]") {
    CandidateVariant site;
    site.phase_set = 100;
    site.counts.category = VariantCategory::CleanHetSnp;
    site.hap_to_cons_alle = {-1, 0, 1};
    CHECK(is_phase_set_anchor(site));
    site.hap_to_cons_alle = {-1, 1, 0};
    CHECK(is_phase_set_anchor(site));
    site.hap_to_cons_alle = {-1, 1, 2};
    CHECK(is_phase_set_anchor(site));
    site.hap_to_cons_alle = {-1, 2, 1};
    CHECK(is_phase_set_anchor(site));
    site.hap_to_cons_alle = {-1, 2, 3};
    CHECK_FALSE(is_phase_set_anchor(site));
    site.hap_to_cons_alle = {-1, 1, 1};
    CHECK_FALSE(is_phase_set_anchor(site));
    site.hap_to_cons_alle = {-1, 0, 0};
    CHECK_FALSE(is_phase_set_anchor(site));
    site.hap_to_cons_alle = {-1, -1, 1};
    CHECK_FALSE(is_phase_set_anchor(site));
    site.hap_to_cons_alle = {-1, 1, -1};
    CHECK_FALSE(is_phase_set_anchor(site));
    site.hap_to_cons_alle = {-1, 0, 1};
    site.counts.category = VariantCategory::CleanHom;
    CHECK_FALSE(is_phase_set_anchor(site));
    site.counts.category = VariantCategory::NoisyCandHom;
    CHECK_FALSE(is_phase_set_anchor(site));
    site.counts.category = VariantCategory::NoisyCandHet;
    CHECK(is_phase_set_anchor(site));
    site.phase_set = 0;
    CHECK_FALSE(is_phase_set_anchor(site));
    site.phase_set = -1;
    CHECK_FALSE(is_phase_set_anchor(site));
}

TEST_CASE("de novo MSA uses each clustered read's coverage flags",
          "[msa][coverage]") {
    Options opts;
    opts.min_af = 0.3;
    NoisyReadInfo info;
    info.n_reads = 6;
    const std::string reference = "ACGTACGTACGTACGTACGT";
    for (int i = 0; i < info.n_reads; ++i) {
        const std::string sequence = i == 5 ? reference.substr(1) : reference;
        std::vector<uint8_t> bases;
        for (const char base : sequence) bases.push_back(base_to_nt4(base));
        info.noisy_read_ids.push_back(100 + i);
        info.lens.push_back(static_cast<int>(bases.size()));
        info.seqs.push_back(std::move(bases));
        // The partial read is not in a cluster. Its flag must not be borrowed
        // by cluster zero's full-cover reads, particularly the leading DEL.
        info.fully_covers.push_back(i == 0 ? kNoisyRightCover : kNoisyBothCover);
    }
    std::vector<uint8_t> ref;
    for (const char base : reference) ref.push_back(base_to_nt4(base));
    std::array<int, 2> counts{};
    std::array<std::vector<int>, 2> ids;
    std::array<std::vector<AlnStr>, 2> alignments;
    const int n_cons = wfa_collect_noisy_aln_str_no_ps_hap(
        opts, info, ref.data(), static_cast<int>(ref.size()), true,
        counts, ids, alignments);
    REQUIRE(n_cons > 0);
    bool found = false;
    for (int ci = 0; ci < n_cons; ++ci) {
        for (int ri = 0; ri < counts[ci]; ++ri) {
            const AlnStr& aln = alignments[ci][2 * ri + 1];
            CHECK(ids[ci][ri] != 100);
            CHECK(aln.target_beg == 0);
            CHECK(aln.query_beg == 0);
            CHECK(aln.target_end == aln.aln_len - 1);
            CHECK(aln.query_end == aln.aln_len - 1);
            found |= ids[ci][ri] == 105;
        }
    }
    CHECK(found);
}

TEST_CASE("phased MSA retains excluded partial reads without inventing coverage",
          "[msa][coverage]") {
    Options opts;
    opts.noisy_reg_flank_len = 3;
    // Even a one-point assignment margin must not turn an uncovered end
    // into a haplotype vote.
    opts.msa_ambiguity_margin = 1;
    const std::string reference = "ACG" + std::string(20, 'T') + "GCA";
    NoisyReadInfo info;
    const auto add_read = [&](const std::string& sequence, int hap, int cover) {
        std::vector<uint8_t> bases;
        for (const char base : sequence) bases.push_back(base_to_nt4(base));
        info.noisy_read_ids.push_back(100 + info.n_reads++);
        info.lens.push_back(static_cast<int>(bases.size()));
        info.quals.emplace_back(bases.size(), 40);
        info.seqs.push_back(std::move(bases));
        info.fully_covers.push_back(cover);
        info.haps.push_back(hap);
        info.phase_sets.push_back(hap == 0 ? 0 : 100);
    };
    for (int i = 0; i < 3; ++i) {
        add_read(reference, 1, kNoisyBothCover);
        add_read("ACG" + std::string(24, 'T') + "GCA", 2, kNoisyBothCover);
    }
    int hap = 1;
    int cover = kNoisyRightCover;
    SECTION("a phased suffix was excluded from the MSA") {}
    SECTION("an unphased suffix has the same coverage limits") { hap = 0; }
    SECTION("a phased prefix was excluded from the MSA") { cover = kNoisyLeftCover; }
    SECTION("an unphased prefix has the same coverage limits") {
        hap = 0;
        cover = kNoisyLeftCover;
    }
    SECTION("a read with neither covered end contributes no alignment") {
        hap = 0;
        cover = kNoisyNoCover;
    }
    add_read(cover == kNoisyRightCover ? std::string(10, 'T') + "GCA"
                                      : "ACG" + std::string(10, 'T'), hap, cover);
    std::vector<uint8_t> ref;
    for (const char base : reference) ref.push_back(base_to_nt4(base));
    std::array<int, 2> counts{};
    std::array<std::vector<int>, 2> ids;
    std::array<std::vector<AlnStr>, 2> alignments;
    std::vector<UnassignedMsaRead> unassigned;
    REQUIRE(wfa_collect_noisy_aln_str_with_ps_hap(
        opts, false, info, 100, 3, 3, ref.data(), static_cast<int>(ref.size()),
        true, counts, ids, alignments, &unassigned) == 2);
    // The suffix fits both consensuses equally. Its missing repeat prefix
    // provides neither a haplotype vote nor deletion observations.
    REQUIRE(counts[0] == 3);
    REQUIRE(counts[1] == 3);
    if (cover == kNoisyNoCover) {
        CHECK(unassigned.empty());
        return;
    }
    REQUIRE(unassigned.size() == 1);
    CHECK(unassigned[0].read_id == 106);
    for (const AlnStr& aln : unassigned[0].ref_read) {
        CHECK(aln.target_beg == 0);
        CHECK(aln.target_end == aln.aln_len - 1);
        if (cover == kNoisyRightCover) {
            CHECK(aln.query_beg > 0);
            CHECK(aln.query_end == aln.aln_len - 1);
        } else {
            CHECK(aln.query_beg == 0);
            CHECK(aln.query_end < aln.aln_len - 1);
        }
    }
}

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

}  // namespace


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

} // namespace

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


