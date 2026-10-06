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

TEST_CASE("local MSA recall preserves fixed consensus membership",
          "[msa][recovery][observations]") {
    Options opts;
    opts.recall_unplaced_msa_insertions = true;
    opts.add_unplaced_msa_observations = false;
    opts.msa_ambiguity_margin = 1;
    const std::string reference = "ACG" + std::string(20, 'T') + "GCA";
    bool literal_ref = false, different_locus = false, tandem_insertion = false;
    SECTION("two alternate insertions share one reference anchor") {}
    SECTION("literal reference and one insertion do not establish the contrast") {
        literal_ref = true;
    }
    SECTION("insertions at different anchors do not establish the contrast") {
        different_locus = true;
    }
    SECTION("two tandem-repeat insertion lengths retain the same fixed consensuses") {
        tandem_insertion = true;
    }
    NoisyReadInfo info;
    for (int i = 0; i < 7; ++i) {
        const int hap = i == 6 ? 0 : i % 2 + 1;
        const std::string sequence = tandem_insertion
            ? "ACG" + std::string(20, 'T') + (hap == 1 ? "TATA" : "TATATATA") + "GCA"
            : different_locus && hap == 2
            ? "ACG" + std::string(20, 'T') + "G" + std::string(8, 'A') + "CA"
            : "ACG" + std::string(hap == 1 ? (literal_ref ? 20 : 24) : 28, 'T') + "GCA";
        std::vector<uint8_t> bases;
        for (const char base : sequence) bases.push_back(base_to_nt4(base));
        info.noisy_read_ids.push_back(100 + i);
        info.lens.push_back(static_cast<int>(bases.size()));
        info.quals.emplace_back(bases.size(), 40);
        info.seqs.push_back(std::move(bases));
        info.fully_covers.push_back(kNoisyBothCover);
        info.haps.push_back(hap);
        info.phase_sets.push_back(hap == 0 ? 0 : 100);
        ++info.n_reads;
    }
    std::vector<uint8_t> ref;
    for (const char base : reference) ref.push_back(base_to_nt4(base));
    std::array<int, 2> counts{};
    std::array<std::vector<int>, 2> ids;
    std::array<std::vector<AlnStr>, 2> alignments;
    std::vector<UnassignedMsaRead> recalled;
    REQUIRE(wfa_collect_noisy_aln_str_with_ps_hap(
        opts, false, info, 100, 3, 3, ref.data(), static_cast<int>(ref.size()),
        true, counts, ids, alignments, &recalled) == 2);
    CHECK(counts[0] == 3);
    CHECK(counts[1] == 3);
    if (literal_ref || different_locus) {
        CHECK(recalled.empty());
        return;
    }
    REQUIRE(recalled.size() == 1);
    CHECK(recalled[0].read_id == 106);
    const std::array<AlnStr, 2> consensuses{alignments[0][0], alignments[1][0]};
    if (tandem_insertion) return;
    VariantKey key;
    key.pos = 103;
    key.type = VariantType::Insertion;
    key.alt = "TTTT";
    CHECK(call_msa_site_allele(recalled[0].ref_read, key, 100, &consensuses) == 0);
    key.alt = "TTTTTTTT";
    CHECK(call_msa_site_allele(recalled[0].ref_read, key, 100, &consensuses) == 1);
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
        // Recovery seams include both VCF anchors: [beg, end].
        CHECK(allele_depths_call_het(het_candidate(1000), opts));
        CHECK(allele_depths_call_het(het_candidate(2000), opts));
        CHECK_FALSE(allele_depths_call_het(het_candidate(2001), opts));
    }
    SECTION("single-anchor seams admit their SNP and indel boundary rows") {
        constexpr hts_pos_t kAnchor = 2000;
        Options singleton = opts;
        singleton.retry_windows = {{kAnchor, kAnchor}};
        for (const VariantType type : {VariantType::Snp, VariantType::Insertion,
                                       VariantType::Deletion}) {
            CandidateVariant boundary = het_candidate(kAnchor);
            boundary.key.type = type;
            if (type != VariantType::Snp) {
                boundary.key.pos = kAnchor + 1;
                boundary.lcd_var_i_to_cate = kCandCleanHetIndel;
            }
            REQUIRE(boundary.key.sort_pos() == kAnchor);
            CHECK(allele_depths_call_het(boundary, singleton));
            --boundary.key.pos;
            CHECK_FALSE(allele_depths_call_het(boundary, singleton));
            boundary.key.pos += 2;
            CHECK_FALSE(allele_depths_call_het(boundary, singleton));
        }
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

TEST_CASE("MSA conflicts inside one source PS can request a focused retry",
          "[recovery][msa][admission][source-conflict]") {
    const auto make_chunk = [] {
        PhasingChunk chunk;
        CandidateVariant left = het_candidate(1000);
        CandidateVariant right = het_candidate(2000);
        left.key.type = VariantType::Insertion;
        left.lcd_var_i_to_cate = kCandNoisyCandHet;
        left.counts.category = VariantCategory::NoisyCandHet;
        left.msa_verified = true;
        right.key.type = VariantType::Deletion;
        right.counts.category = VariantCategory::NoisyCandHet;
        right.msa_verified = true;
        left.phase_set = right.phase_set = 1000;
        left.hap_to_cons_alle = right.hap_to_cons_alle = {-1, 0, 1};
        chunk.candidates = {left, right};
        for (int ri = 0; ri < 26; ++ri) {
            ReadRecord read;
            read.qname = "molecule" + std::to_string(ri);
            read.beg = 900;
            read.end = 2100;
            read.mapq = 60;
            chunk.reads.push_back(std::move(read));
            ReadVariantProfile profile;
            profile.read_id = ri;
            profile.start_var_idx = 0;
            profile.end_var_idx = 1;
            profile.alleles = {ri < 14 ? 1 : 0, 1};
            chunk.read_var_profile.push_back(std::move(profile));
            chunk.haps.push_back(ri < 14 ? 2 : 1);
            chunk.phase_sets.push_back(1000);
        }
        return chunk;
    };
    PhasingChunk chunk = make_chunk();
    Options opts;
    opts.min_depth = 6;
    const auto admits = [&opts](const PhasingChunk& input) {
        return msa_source_conflict_is_supported(input, {0, 1}, opts);
    };
    SECTION("a nominal PS with a conflicting weak internal edge is retried") {
        CHECK(admits(chunk));
    }
    SECTION("different source gauges do not establish an internal conflict") {
        chunk.candidates[1].phase_set = 2000;
        CHECK_FALSE(admits(chunk));
    }
    SECTION("homozygotes and unphased rows cannot name a path edge") {
        chunk.candidates[1].counts.category = VariantCategory::NoisyCandHom;
        CHECK_FALSE(admits(chunk));
        chunk = make_chunk();
        chunk.candidates[1].phase_set = 0;
        CHECK_FALSE(admits(chunk));
    }
    SECTION("clean SNP pairs do not request indel MSA") {
        for (auto& site : chunk.candidates) site.key.type = VariantType::Snp;
        CHECK_FALSE(admits(chunk));
    }
    SECTION("both consistent haplotypes already support the ordinary path") {
        chunk.read_var_profile[0].alleles = {0, 0};
        CHECK_FALSE(admits(chunk));
    }
    SECTION("small or nonsignificant conflicting cohorts are insufficient") {
        chunk.reads.resize(5);
        CHECK_FALSE(admits(chunk));
        chunk = make_chunk();
        for (size_t ri = 14; ri < 25; ++ri)
            chunk.read_var_profile[ri].alleles = {1, 1};
        CHECK_FALSE(admits(chunk));
    }
    SECTION("low and unknown MAPQ, skipped and other PS reads do not contribute") {
        for (int mode = 0; mode < 4; ++mode) {
            chunk = make_chunk();
            for (size_t ri = 14; ri < chunk.reads.size(); ++ri) {
                if (mode == 0) chunk.reads[ri].mapq = 29;
                if (mode == 1) chunk.reads[ri].mapq = 255;
                if (mode == 2) chunk.reads[ri].is_skipped = true;
                if (mode == 3) chunk.phase_sets[ri] = 2000;
            }
            CHECK_FALSE(admits(chunk));
        }
    }
    SECTION("a profile outside the physical molecule cannot supply a pair") {
        for (size_t ri = 14; ri < chunk.reads.size(); ++ri)
            chunk.reads[ri].end = 1500;
        CHECK_FALSE(admits(chunk));
    }
    SECTION("missing or third-allele observations are not conflicts") {
        for (const int allele : {-1, 2}) {
            chunk = make_chunk();
            for (size_t ri = 14; ri < chunk.reads.size(); ++ri)
                chunk.read_var_profile[ri].alleles[0] = allele;
            CHECK_FALSE(admits(chunk));
        }
    }
    SECTION("co-located alternatives are one locus, not an internal cut") {
        chunk.candidates[1].key.pos = chunk.candidates[0].key.pos;
        CHECK_FALSE(admits(chunk));
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
    SECTION("reference call error measures the whole verified footprint") {
        chunk = shifted_insertion_chunk("8M", "CAAAAAGC");
        bam1_t* read = chunk.reads[0].alignment.get();
        bam_get_qual(read)[3] = 17;
        double error = 1.0;
        CHECK(call(deletion) == -1);
        CHECK(bam_equivalent_deletion_allele(
            read, deletion, reference, 0, header.get(), 10, &error) == 0);
        CHECK(error > 0.019);
        CHECK(error < 0.05);
        bam_get_qual(read)[3] = 10;
        CHECK(bam_equivalent_deletion_allele(
            read, deletion, reference, 0, header.get(), 10, &error) == 0);
        CHECK(error > 0.05);
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
    SECTION("a shifted deletion and compensating mismatches preserve the exact allele") {
        // Deleting AG at 105 gives CAAAAC. The CIGAR instead deletes AA at
        // 104 and calls the surviving G as A: the final haplotype is identical.
        chunk = shifted_insertion_chunk("4M2D2M", "CAAAAC");
        deletion.key.pos = 105;
        deletion.key.ref_len = 2;
        CHECK(call(deletion) == -1);
        CHECK(bam_matches_deletion_sequence(chunk.reads[0].alignment.get(),
            deletion, reference, 0, header.get(), 30));
        CHECK(deletion.key.pos == 105);
        CHECK(deletion.key.ref_len == 2);
        deletion.key.ref_len = 1;
        CHECK(call(deletion) == -1);
        CHECK_FALSE(bam_matches_deletion_sequence(chunk.reads[0].alignment.get(),
            deletion, reference, 0, header.get(), 20));
        deletion.key.ref_len = 2;
        bam_get_qual(chunk.reads[0].alignment.get())[4] = 29;
        CHECK(call(deletion) == -1);
        CHECK_FALSE(bam_matches_deletion_sequence(chunk.reads[0].alignment.get(),
            deletion, reference, 0, header.get(), 30));
        CHECK(bam_matches_deletion_sequence(chunk.reads[0].alignment.get(),
            deletion, reference, 0, header.get(), 20));
        bam_get_qual(chunk.reads[0].alignment.get())[4] = 19;
        CHECK_FALSE(bam_matches_deletion_sequence(chunk.reads[0].alignment.get(),
            deletion, reference, 0, header.get(), 20));
        bam_get_qual(chunk.reads[0].alignment.get())[4] = 255;
        CHECK(call(deletion) == -1);
        CHECK_FALSE(bam_matches_deletion_sequence(chunk.reads[0].alignment.get(),
            deletion, reference, 0, header.get(), 20));
    }
    SECTION("equal deletion lengths with a different final sequence remain unknown") {
        chunk = shifted_insertion_chunk("4M2D2M", "CAAAGC");
        deletion.key.pos = 105;
        deletion.key.ref_len = 2;
        CHECK(call(deletion) == -1);
        CHECK_FALSE(bam_matches_deletion_sequence(chunk.reads[0].alignment.get(),
            deletion, reference, 0, header.get(), 20));
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

TEST_CASE("recovery jointly calls complementary MSA deletions",
          "[msa][recovery][representation][deletion-pair]") {
    const auto make_chunk = [](const std::string& cigar, const std::string& sequence) {
        PhasingChunk chunk = shifted_insertion_chunk(cigar, sequence);
        chunk.ref_seq = "CAAAAAGC";
        CandidateVariant& first = chunk.candidates.front();
        first.key = del_key(103, 1);
        first.msa_verified = true;
        first.is_homopolymer_indel = true;
        first.phase_set = 100;
        first.hap_to_cons_alle = {-1, 0, 1};
        first.counts.category = VariantCategory::NoisyCandHet;
        first.counts.n_uniq_alles = 2;
        first.counts.alle_covs = {20, 20};
        first.counts.total_cov = 40;
        CandidateVariant second = first;
        second.key.ref_len = 2;
        second.hap_to_cons_alle = {-1, 1, 0};
        chunk.candidates.push_back(second);
        return chunk;
    };
    Options opts;
    opts.min_bq = 30;
    PhasingChunk chunk = make_chunk("4M2D2M", "CAAAGC");
    const auto fill = [&]() { return backfill_msa_retry_deletions(chunk, opts, 102, 102, false); };
    SECTION("ordinary backfill does not promote physical pair calls to rescue evidence") {
        CHECK(backfill_msa_observations(chunk, opts, 102, 102) == 0);
        CHECK(chunk.read_var_profile[0].start_var_idx == -1);
        CHECK(fill() == 2);
    }
    SECTION("a verified two-base ALT supplies both missing binary contrasts") {
        CHECK(fill() == 2);
        CHECK(chunk.read_var_profile[0].alleles == std::vector<int>{0, 1});
        CHECK(chunk.candidates.size() == 2);
        CHECK(chunk.candidates[0].key.ref_len == 1);
        CHECK(chunk.candidates[1].key.ref_len == 2);
        CHECK(chunk.candidates[0].hap_to_cons_alle == std::array<int, 3>{-1, 0, 1});
        CHECK(chunk.candidates[1].hap_to_cons_alle == std::array<int, 3>{-1, 1, 0});
        for (const auto& candidate : chunk.candidates) {
            CHECK(candidate.counts.alle_covs == std::vector<int>{20, 20});
            CHECK(candidate.counts.total_cov == 40);
            CHECK(candidate.phase_set == 100);
        }
        CHECK(fill() == 0);
    }
    SECTION("exact and left-shifted one-base ALTs keep the other contrast zero") {
        for (const auto& cigar : {"3M1D4M", "2M1D5M"}) {
            chunk = make_chunk(cigar, "CAAAAGC");
            CHECK(fill() == 2);
            CHECK(chunk.read_var_profile[0].alleles == std::vector<int>{1, 0});
        }
    }
    SECTION("a third deletion length or literal REF does not identify a haplotype") {
        for (const auto& read : std::array<std::pair<std::string, std::string>, 2>{
                std::make_pair("3M3D2M", "CAAGC"),
                std::make_pair("8M", "CAAAAAGC")}) {
            chunk = make_chunk(read.first, read.second);
            CHECK(fill() == 0);
            CHECK(chunk.read_var_profile[0].start_var_idx == -1);
        }
    }
    SECTION("existing MSA contrasts are never overwritten or completed from another gauge") {
        update_read_var_profile_with_allele(0, 1, -1, chunk.read_var_profile[0]);
        CHECK(fill() == 0);
        CHECK(chunk.read_var_profile[0].alleles == std::vector<int>{1});
    }
    SECTION("independent blocks cannot manufacture a complementary locus") {
        chunk.candidates[1].phase_set = 200;
        CHECK(fill() == 0);
        chunk.candidates[1].phase_set = 100;
        chunk.candidates[1].hap_to_cons_alle = {-1, 0, 1};
        CHECK(fill() == 0);
    }
    SECTION("both rows must lie in the targeted seam") {
        CHECK(backfill_msa_retry_deletions(chunk, opts, 103, 106, false) == 0);
    }
    SECTION("compound edits and low or missing qualities remain unknown") {
        chunk = make_chunk("3M1I1M2D2M", "CAATAGC");
        CHECK(fill() == 0);
        chunk = make_chunk("4M2D2M", "CAAAGC");
        for (const int quality : {29, 255}) {
            bam_get_qual(chunk.reads[0].alignment.get())[3] = quality;
            CHECK(fill() == 0);
        }
        bam_get_qual(chunk.reads[0].alignment.get())[3] = 40;
        for (const int mapq : {5, 255}) {
            chunk.reads[0].mapq = mapq;
            CHECK(fill() == 0);
        }
    }
    SECTION("a third MSA alternative prevents a two-row projection") {
        CandidateVariant other = chunk.candidates[0];
        other.key.ref_len = 3;
        chunk.candidates.push_back(other);
        CHECK(fill() == 0);
    }
}

TEST_CASE("isolated MSA deletion calls diagnose a reversed source edge",
          "[msa][recovery][representation][isolated-deletion-retry]") {
    const auto make_chunk = [](const std::string& cigar, const std::string& sequence) {
        PhasingChunk chunk = shifted_insertion_chunk(cigar, sequence);
        chunk.ref_seq = "CAAAAAGC";
        CandidateVariant& deletion = chunk.candidates.front();
        deletion.key = del_key(103, 2);
        deletion.is_homopolymer_indel = true;
        deletion.phase_set = 100;
        deletion.hap_to_cons_alle = {-1, 0, 1};
        deletion.counts.category = VariantCategory::NoisyCandHet;
        deletion.counts.n_uniq_alles = 2;
        deletion.counts.alle_covs = {20, 20};
        deletion.counts.total_cov = 40;
        return chunk;
    };
    Options opts;
    PhasingChunk chunk = make_chunk("4M2D2M", "CAAAGC");
    const auto fill = [&]() { return backfill_msa_retry_deletions(chunk, opts, 102, 102, true); };
    SECTION("exact and shifted physical calls preserve the frozen source genotype") {
        CHECK(backfill_msa_observations(chunk, opts, 102, 102) == 0);
        for (const auto& cigar : {"3M2D3M", "4M2D2M"}) {
            chunk = make_chunk(cigar, "CAAAGC");
            CHECK(fill() == 1);
            CHECK(chunk.read_var_profile[0].alleles == std::vector<int>{1});
            CHECK(fill() == 0);
            CHECK(chunk.candidates[0].key.pos == 103);
            CHECK(chunk.candidates[0].key.ref_len == 2);
            CHECK(chunk.candidates[0].key.alt.empty());
            CHECK(chunk.candidates[0].hap_to_cons_alle == std::array<int, 3>{-1, 0, 1});
            CHECK(chunk.candidates[0].phase_set == 100);
            CHECK(chunk.candidates[0].counts.alle_covs == std::vector<int>{20, 20});
            CHECK(chunk.candidates[0].counts.total_cov == 40);
        }
    }
    SECTION("ordinary flanks retain the complementary-row admission rule") {
        CHECK(backfill_msa_retry_deletions(chunk, opts, 102, 102, false) == 0);
        CHECK(chunk.read_var_profile[0].start_var_idx == -1);
    }
    SECTION("verified literal reference can expose the opposite source conflict") {
        chunk = make_chunk("8M", "CAAAAAGC");
        CHECK(fill() == 1);
        CHECK(chunk.read_var_profile[0].alleles == std::vector<int>{0});
    }
    SECTION("third lengths and compound events do not identify this binary allele") {
        chunk = make_chunk("3M3D2M", "CAAGC");
        CHECK(fill() == 0);
        chunk = make_chunk("3M1I1M2D2M", "CAATAGC");
        CHECK(fill() == 0);
        chunk = make_chunk("4M2D2M", "CAAAGC");
        chunk.candidates[0].key.alt = "A";
        CHECK(fill() == 0);
    }
    SECTION("an additional MSA indel makes an isolated contrast ambiguous") {
        CandidateVariant alternative = chunk.candidates[0];
        alternative.key = ins_key(103, "A");
        chunk.candidates.push_back(alternative);
        CHECK(fill() == 0);
    }
    SECTION("existing calls and ineligible source rows remain untouched") {
        update_read_var_profile_with_allele(0, 0, -1, chunk.read_var_profile[0]);
        CHECK(fill() == 0);
        CHECK(chunk.read_var_profile[0].alleles == std::vector<int>{0});
        for (const hts_pos_t ps : {hts_pos_t{0}, hts_pos_t{-1}}) {
            chunk = make_chunk("4M2D2M", "CAAAGC");
            chunk.candidates[0].phase_set = ps;
            CHECK(fill() == 0);
        }
        chunk = make_chunk("4M2D2M", "CAAAGC");
        chunk.candidates[0].is_homopolymer_indel = false;
        CHECK(fill() == 0);
    }
    SECTION("qualities, mapping confidence and seam membership remain required") {
        for (const int quality : {29, 255}) {
            bam_get_qual(chunk.reads[0].alignment.get())[3] = quality;
            CHECK(fill() == 0);
        }
        bam_get_qual(chunk.reads[0].alignment.get())[3] = 40;
        for (const int mapq : {29, 255}) {
            chunk.reads[0].mapq = mapq;
            CHECK(fill() == 0);
        }
        chunk.reads[0].mapq = 60;
        CHECK(backfill_msa_retry_deletions(chunk, opts, 103, 107, true) == 0);
    }
    SECTION("restored pairs request a retry without changing source orientations") {
        CandidateVariant deletion = chunk.candidates.front();
        CandidateVariant snp = het_candidate(100);
        snp.key.alt = "T";
        snp.hap_to_cons_alle = {-1, 0, 1};
        snp.phase_set = 100;
        chunk.candidates = {snp, deletion};
        for (int ri = 1; ri < 6; ++ri) {
            ReadRecord read;
            read.qname = "spanner" + std::to_string(ri);
            read.beg = 100;
            read.end = 107;
            read.mapq = 60;
            read.alignment.reset(bam_dup1(chunk.reads.front().alignment.get()));
            chunk.reads.push_back(std::move(read));
            chunk.read_var_profile.push_back(ReadVariantProfile{});
        }
        for (auto& profile : chunk.read_var_profile) {
            profile.start_var_idx = profile.end_var_idx = 0;
            profile.alleles = {0};
        }
        chunk.phase_sets.assign(6, 100);
        chunk.haps.assign(6, 1);
        CHECK_FALSE(msa_source_conflict_is_supported(chunk, {0, 1}, opts));
        CHECK(fill() == 6);
        CHECK(msa_source_conflict_is_supported(chunk, {0, 1}, opts));
        CHECK(chunk.candidates[0].hap_to_cons_alle == chunk.candidates[1].hap_to_cons_alle);
        CHECK(chunk.phase_sets == std::vector<hts_pos_t>(6, 100));
        CHECK(chunk.haps == std::vector<int>(6, 1));
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
    SECTION("an exact-position REF must not hide a certified shifted ALT") {
        PhasingChunk chunk = make_chunk("5M1D2M", "CAAAAGC");
        int qi = -1;
        REQUIRE(bam_exact_indel_allele(chunk.reads[0].alignment.get(),
            chunk.candidates[0], 30, &qi) == 0);
        REQUIRE(fill(chunk) == 1);
        CHECK(chunk.read_var_profile[0].alleles == std::vector<int>{1});
        CHECK(chunk.candidates[0].phase_set == 100);
        CHECK(chunk.candidates[0].key.pos == 103);
    }
    SECTION("a different event keeps its exact binary contrast") {
        PhasingChunk chunk = make_chunk("4M1I4M", "CAAAAAAGC");
        chunk.candidates[0].key.ref_len = 3;
        REQUIRE(fill(chunk) == 1);
        CHECK(chunk.read_var_profile[0].alleles == std::vector<int>{0});
    }
    SECTION("a shifted ALT preserves an absent gauge and rejects a contradictory one") {
        for (const bool contradiction : {false, true}) {
            PhasingChunk chunk = make_chunk("5M1D2M", "CAAAAGC");
            if (contradiction) {
                chunk.candidates[1].hap_to_cons_alle[1] = 1;
                chunk.candidates[1].hap_to_cons_alle[2] = 0;
            } else {
                chunk.candidates[1].phase_set = 200;
            }
            CHECK(fill(chunk) == (contradiction ? 0 : 1));
            if (contradiction)
                CHECK(chunk.read_var_profile[0].start_var_idx == -1);
            else
                CHECK(chunk.read_var_profile[0].alleles == std::vector<int>{0});
        }
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

TEST_CASE("MSA insertion rows retain the other verified insertion allele",
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
    const AlnStr four = alignment("ACG----TGC", "ACGTTTTTGC");
    const AlnStr eight = alignment("ACG--------TGC", "ACGTTTTTTTTTGC");
    const AlnStr seven = alignment("ACG-------TGC", "ACGTTTTTTTTGC");
    const AlnStr six = alignment("ACG------TGC", "ACGTTTTTTTGC");
    const VariantKey four_key = ins_key(103, "TTTT");
    const VariantKey eight_key = ins_key(103, "TTTTTTTT");
    std::array<AlnStr, 2> consensuses{four, eight};
    SECTION("exact diploid alleles project to complementary rows") {
        CHECK(call_msa_site_allele({four, four}, four_key, 100, &consensuses) == 1);
        CHECK(call_msa_site_allele({eight, eight}, four_key, 100, &consensuses) == 0);
        CHECK(call_msa_site_allele({four, four}, eight_key, 100, &consensuses) == 0);
        CHECK(call_msa_site_allele({eight, eight}, eight_key, 100, &consensuses) == 1);
        CHECK(call_msa_site_allele({eight, eight}, four_key, 100) == -1);
        CHECK(call_msa_site_allele({seven, seven}, four_key, 100, &consensuses) == -1);
        CHECK(call_msa_site_allele({four, eight}, four_key, 100, &consensuses) == -1);
    }
    SECTION("missing selected allele or incomplete contrast stays unknown") {
        consensuses = {eight, eight};
        CHECK(call_msa_site_allele({eight, eight}, four_key, 100, &consensuses) == -1);
        consensuses = {four, eight};
        consensuses[0].target_beg = 1;
        CHECK(call_msa_site_allele({eight, eight}, four_key, 100, &consensuses) == -1);
        consensuses = {four, eight};
        consensuses[0].query_aln[0] = 4;
        CHECK(call_msa_site_allele({eight, eight}, four_key, 100, &consensuses) == -1);
    }
    SECTION("different consensus flanks do not certify an insertion contrast") {
        consensuses[0].query_aln[0] = base_to_nt4('T');
        CHECK(call_msa_site_allele({eight, eight}, four_key, 100, &consensuses) == -1);
    }
    SECTION("existing one-error recall uses the diploid contrast without merging") {
        CandidateVariant first;
        first.key = four_key;
        first.counts.category = VariantCategory::NoisyCandHet;
        first.counts.alle_covs = {1, 1};
        first.counts.total_cov = 2;
        CandidateVariant second = first;
        second.key = eight_key;
        std::vector<CandidateVariant> candidates{first, second};
        std::vector<ReadVariantProfile> profiles(2);
        UnassignedMsaRead supported;
        supported.read_id = 0;
        supported.ref_read = {seven, seven};
        UnassignedMsaRead ambiguous;
        ambiguous.read_id = 1;
        ambiguous.ref_read = {six, six};
        Options opts;
        add_msa_site_observations(opts, {supported, ambiguous}, 100,
                                  candidates, profiles, &consensuses);
        REQUIRE(candidates.size() == 2);
        CHECK(exact_comp_var_site(&candidates[0].key, &four_key) == 0);
        CHECK(exact_comp_var_site(&candidates[1].key, &eight_key) == 0);
        CHECK(profiles[0].alleles == std::vector<int>{0, 1});
        CHECK(profiles[1].start_var_idx == -1);
        CHECK(candidates[0].counts.alle_covs == std::vector<int>{2, 1});
        CHECK(candidates[1].counts.alle_covs == std::vector<int>{1, 2});
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
    reference += std::string(120, 'C');
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
    const auto add_snp = [](PhasingChunk& target) {
        CandidateVariant snp;
        snp.key.type = VariantType::Snp;
        snp.key.pos = 201;
        snp.key.ref_len = 1;
        snp.key.alt = "G";
        snp.ref_base = 1; // C in the physically aligned flank.
        snp.counts.category = VariantCategory::CleanHetSnp;
        snp.phase_set = 101;
        snp.hap_to_cons_alle = {-1, 0, 1};
        target.candidates.push_back(snp);
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
    SECTION("recovery backfill does not call an equivalent shifted ALT reference") {
        Options opts;
        add_snp(chunk);
        REQUIRE(backfill_msa_observations(chunk, opts, 100, 100) == 1);
        CHECK(chunk.read_var_profile[0].alleles == std::vector<int>{1});
        CHECK(chunk.candidates[0].key.pos == 101);
        CHECK(chunk.candidates[0].key.alt == "ATCT");
        CHECK(chunk.phase_sets == std::vector<hts_pos_t>{101});
    }
    SECTION("supplementary insertion calls preserve the discovery genotype census") {
        Options opts;
        add_snp(chunk);
        auto& counts = chunk.candidates[0].counts;
        counts.alle_covs = {3, 2};
        counts.total_cov = 5;
        counts.ref_cov = 3;
        counts.alt_cov = 2;
        counts.allele_fraction = 0.4;
        counts.forward_ref = 2;
        counts.reverse_ref = 1;
        counts.forward_alt = counts.reverse_alt = 1;
        REQUIRE(backfill_msa_observations(chunk, opts, 100, 100) == 1);
        CHECK(chunk.read_var_profile[0].alleles == std::vector<int>{1});
        CHECK(counts.alle_covs == std::vector<int>{3, 2});
        CHECK(counts.total_cov == 5);
        CHECK(counts.ref_cov == 3);
        CHECK(counts.alt_cov == 2);
        CHECK(counts.allele_fraction == Approx(0.4));
        CHECK(counts.forward_ref == 2);
        CHECK(counts.reverse_ref == 1);
        CHECK(counts.forward_alt == 1);
        CHECK(counts.reverse_alt == 1);
        CHECK(chunk.candidates[0].hap_to_cons_alle == std::array<int, 3>{-1, 1, 0});
        CHECK(chunk.haps == std::vector<int>{2});
        CHECK(chunk.phase_sets == std::vector<hts_pos_t>{101});
        CHECK(backfill_msa_observations(chunk, opts, 100, 100) == 0);
        CHECK(counts.alle_covs == std::vector<int>{3, 2});
        CHECK(counts.total_cov == 5);
    }
    SECTION("supplementary SNP calls do not reclassify the source genotype") {
        Options opts;
        add_snp(chunk);
        auto& snp = chunk.candidates.back();
        snp.msa_verified = true;
        snp.counts.category = VariantCategory::NoisyCandHet;
        snp.lcd_var_i_to_cate = kCandNoisyCandHet;
        snp.counts.alle_covs = {2, 3};
        snp.counts.ref_cov = 2;
        snp.counts.alt_cov = 3;
        snp.counts.total_cov = 5;
        snp.counts.allele_fraction = 0.6;
        REQUIRE(backfill_msa_observations(chunk, opts, 201, 201) == 1);
        CHECK(chunk.read_var_profile[0].alleles == std::vector<int>{0});
        CHECK(snp.counts.alle_covs == std::vector<int>{2, 3});
        CHECK(snp.counts.total_cov == 5);
        CHECK(snp.counts.ref_cov == 2);
        CHECK(snp.counts.alt_cov == 3);
        CHECK(snp.counts.allele_fraction == Approx(0.6));
        CHECK(snp.counts.category == VariantCategory::NoisyCandHet);
        CHECK(snp.phase_set == 101);
        CHECK(snp.hap_to_cons_alle == std::array<int, 3>{-1, 0, 1});
        CHECK(chunk.haps == std::vector<int>{2});
        CHECK(chunk.phase_sets == std::vector<hts_pos_t>{101});
    }
    SECTION("recovery preserves an existing MSA contrast and independent row gauges") {
        Options opts;
        update_read_var_profile_with_allele(0, 0, -1, chunk.read_var_profile[0]);
        CHECK(backfill_msa_observations(chunk, opts, 100, 100) == 0);
        CHECK(chunk.read_var_profile[0].alleles == std::vector<int>{0});
        chunk = make_chunk();
        CandidateVariant other = chunk.candidates[0];
        other.key.alt = "ATCTATCT";
        other.hap_to_cons_alle = {-1, 0, 1};
        chunk.candidates.push_back(other);
        add_snp(chunk);
        REQUIRE(backfill_msa_observations(chunk, opts, 100, 100) == 2);
        CHECK(chunk.read_var_profile[0].alleles == std::vector<int>{1, 0});
        CHECK(chunk.candidates.size() == 3);
        CHECK(chunk.candidates[0].hap_to_cons_alle == std::array<int, 3>{-1, 1, 0});
        CHECK(chunk.candidates[1].hap_to_cons_alle == std::array<int, 3>{-1, 0, 1});
        CHECK(chunk.haps == std::vector<int>{2});
    }
    SECTION("an absent SNP gauge does not reinterpret the existing source projection") {
        Options opts;
        REQUIRE(backfill_msa_observations(chunk, opts, 100, 100) == 1);
        CHECK(chunk.read_var_profile[0].alleles == std::vector<int>{0});
        CHECK(chunk.haps == std::vector<int>{2});
        CHECK(chunk.phase_sets == std::vector<hts_pos_t>{101});
        REQUIRE(chunk.pending_msa_observations.size() == 1);
        CHECK(chunk.pending_msa_observations[0].allele == 1);
        CHECK_FALSE(chunk.pending_msa_observations[0].update_counts);
    }
    SECTION("a certified third repeat length retracts both supplementary REF calls") {
        chunk = make_chunk(45, "ATCTATCT");
        chunk.haps = {0};
        chunk.phase_sets = {kUnphasedReadPhaseSet};
        chunk.candidates[0].counts.alle_covs = {3, 2};
        chunk.candidates[0].counts.total_cov = 5;
        CandidateVariant other = chunk.candidates[0];
        other.key.alt = "ATCTATCTATCT";
        other.hap_to_cons_alle = {-1, 0, 1};
        chunk.candidates.push_back(other);
        Options opts;
        REQUIRE(backfill_msa_observations(chunk, opts, 100, 100) == 2);
        CHECK(chunk.read_var_profile[0].alleles == std::vector<int>{0, 0});
        REQUIRE(chunk.pending_msa_observations.size() == 2);
        CHECK(apply_pending_msa_observations(chunk));
        CHECK(chunk.read_var_profile[0].alleles == std::vector<int>{-1, -1});
        CHECK(chunk.rejected_msa_observations.size() == 2);
        CHECK(chunk.candidates[0].counts.total_cov == 5);
        CHECK(chunk.candidates[1].counts.total_cov == 5);
    }
    SECTION("a contradictory independent SNP cannot certify the insertion source gauge") {
        Options opts;
        add_snp(chunk);
        chunk.candidates.back().hap_to_cons_alle = {-1, 1, 0};
        CHECK(backfill_msa_observations(chunk, opts, 100, 100) == 0);
        CHECK(chunk.read_var_profile[0].start_var_idx == -1);
        CHECK(chunk.haps == std::vector<int>{2});
        CHECK(chunk.pending_msa_observations.empty());
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
        chunk = shifted_insertion_chunk("13M4I32M4I" +
            std::to_string(reference.size() - 45) + "M",
            reference.substr(0, 13) + "ATCT" + reference.substr(13, 32) +
            "ATCT" + reference.substr(45));
        chunk.ref_seq = reference;
        chunk.candidates[0].key = ins_key(101, "ATCT");
        CHECK(call(chunk) == -1);
        chunk = shifted_insertion_chunk("25M1D19M4I" +
            std::to_string(reference.size() - 45) + "M",
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

TEST_CASE("pending fixed-consensus calls survive supplementary CIGAR projection",
          "[msa][recovery][observations]") {
    PhasingChunk chunk;
    for (const auto& alt : {std::string("TATA"), std::string("TATATATA")}) {
        CandidateVariant site = het_candidate(1500);
        site.key = ins_key(1500, alt);
        site.msa_verified = true;
        site.counts.category = VariantCategory::NoisyCandHet;
        site.lcd_var_i_to_cate = kCandNoisyCandHet;
        site.phase_set = 1000;
        site.hap_to_cons_alle = {-1, 1, 0};
        site.counts.alle_covs = {3, 3};
        site.counts.ref_cov = site.counts.alt_cov = 3;
        site.counts.total_cov = 6;
        chunk.candidates.push_back(site);
    }
    chunk.candidates[1].hap_to_cons_alle = {-1, 0, 1};
    chunk.reads.emplace_back();
    chunk.haps = {1};
    chunk.phase_sets = {1000};
    ReadVariantProfile profile;
    profile.read_id = 0;
    profile.start_var_idx = 0;
    profile.end_var_idx = 1;
    // The original MSA calls were missing when queued. Later CIGAR backfill
    // reports absence at both ALT anchors, losing the diploid contrast.
    profile.alleles = {0, 0};
    profile.alt_qi = {-1, -1};
    chunk.read_var_profile.push_back(profile);
    chunk.pending_msa_observations = {
        {chunk.candidates[0].key, 0, 1}, {chunk.candidates[1].key, 0, 0}};
    SECTION("verified complementary observations replace only the queued calls") {
        CHECK(apply_pending_msa_observations(chunk));
        CHECK(chunk.read_var_profile[0].alleles == std::vector<int>{1, 0});
        CHECK(chunk.candidates[0].counts.alle_covs == std::vector<int>{3, 4});
        CHECK(chunk.candidates[1].counts.alle_covs == std::vector<int>{4, 3});
        CHECK(chunk.candidates[0].counts.total_cov == 7);
        CHECK(chunk.candidates[1].counts.total_cov == 7);
        CHECK(chunk.read_var_cr != nullptr);
        CHECK(chunk.pending_msa_observations.empty());
        CHECK_FALSE(apply_pending_msa_observations(chunk));
        CHECK(chunk.candidates[0].counts.total_cov == 7);
        CHECK(chunk.haps == std::vector<int>{1});
        CHECK(chunk.phase_sets == std::vector<hts_pos_t>{1000});
        CHECK(chunk.candidates[0].hap_to_cons_alle == std::array<int, 3>{-1, 1, 0});
        CHECK(chunk.candidates[1].hap_to_cons_alle == std::array<int, 3>{-1, 0, 1});
        CHECK(chunk.candidates[0].phase_set == 1000);
        CHECK(chunk.candidates[1].phase_set == 1000);
    }
    SECTION("duplicate calls are counted once and conflicting calls abstain") {
        chunk.pending_msa_observations.push_back(chunk.pending_msa_observations.front());
        chunk.pending_msa_observations.push_back({chunk.candidates[1].key, 0, 1});
        CHECK(apply_pending_msa_observations(chunk));
        CHECK(chunk.candidates[0].counts.total_cov == 7);
        CHECK(chunk.candidates[1].counts.total_cov == 6);
        REQUIRE(chunk.conflicting_msa_observations.size() == 2);
        CHECK(chunk.conflicting_msa_observations[0].allele == 0);
        CHECK(chunk.conflicting_msa_observations[1].allele == 1);
        RecoveryPhaseGauge gauge;
        retain_recovery_bam_evidence(chunk, {{1000, 1000}}, gauge);
        REQUIRE(gauge.conflicting_recalls.size() == 2);
        CHECK(gauge.conflicting_recalls[0].allele == 0);
        CHECK(gauge.conflicting_recalls[1].allele == 1);
    }
    SECTION("supplementary physical calls preserve the discovery census") {
        for (auto& observation : chunk.pending_msa_observations)
            observation.update_counts = false;
        CHECK(apply_pending_msa_observations(chunk));
        CHECK(chunk.read_var_profile[0].alleles == std::vector<int>{1, 0});
        CHECK(chunk.candidates[0].counts.alle_covs == std::vector<int>{3, 3});
        CHECK(chunk.candidates[1].counts.alle_covs == std::vector<int>{3, 3});
        CHECK(chunk.candidates[0].counts.total_cov == 6);
        CHECK(chunk.candidates[1].counts.total_cov == 6);
        CHECK(chunk.haps == std::vector<int>{1});
        CHECK(chunk.phase_sets == std::vector<hts_pos_t>{1000});
        CHECK_FALSE(apply_pending_msa_observations(chunk));
    }
    SECTION("a third allele abstains from both complementary votes") {
        chunk.pending_msa_observations = {
            {chunk.candidates[0].key, 0, -1, false},
            {chunk.candidates[1].key, 0, -1, false}};
        CHECK(apply_pending_msa_observations(chunk));
        CHECK(chunk.read_var_profile[0].alleles == std::vector<int>{-1, -1});
        CHECK(chunk.candidates[0].counts.total_cov == 6);
        CHECK(chunk.candidates[1].counts.total_cov == 6);
        CHECK(chunk.haps == std::vector<int>{1});
        CHECK(chunk.phase_sets == std::vector<hts_pos_t>{1000});
    }
    SECTION("calls in another block preserve the independently oriented read label") {
        const int hap = GENERATE(1, 2);
        chunk.haps[0] = hap;
        chunk.phase_sets[0] = 2000;
        CHECK(apply_pending_msa_observations(chunk));
        CHECK(chunk.read_var_profile[0].alleles == std::vector<int>{1, 0});
        CHECK(chunk.candidates[0].counts.total_cov == 7);
        CHECK(chunk.candidates[1].counts.total_cov == 7);
        CHECK(chunk.haps == std::vector<int>{hap});
        CHECK(chunk.phase_sets == std::vector<hts_pos_t>{2000});
        CHECK(chunk.candidates[0].phase_set == 1000);
        CHECK(chunk.candidates[1].phase_set == 1000);
        CHECK(chunk.candidates[0].hap_to_cons_alle == std::array<int, 3>{-1, 1, 0});
        CHECK(chunk.candidates[1].hap_to_cons_alle == std::array<int, 3>{-1, 0, 1});
    }
    SECTION("verified calls preserve evidence independently of the prior read label") {
        chunk.haps[0] = 2;
        CHECK(apply_pending_msa_observations(chunk));
        CHECK(chunk.read_var_profile[0].alleles == std::vector<int>{1, 0});
        CHECK(chunk.candidates[0].counts.total_cov == 7);
        CHECK(chunk.candidates[1].counts.total_cov == 7);
        CHECK(chunk.haps == std::vector<int>{2});
        CHECK(chunk.phase_sets == std::vector<hts_pos_t>{1000});
    }
    SECTION("a verified MSA call takes precedence over supplementary projection") {
        chunk.pending_msa_observations.push_back({chunk.candidates[0].key, 0, 0, false});
        CHECK(apply_pending_msa_observations(chunk));
        CHECK(chunk.read_var_profile[0].alleles == std::vector<int>{1, 0});
        CHECK(chunk.candidates[0].counts.total_cov == 7);
        CHECK(chunk.candidates[1].counts.total_cov == 7);
        REQUIRE(chunk.conflicting_msa_observations.size() == 2);
        CHECK_FALSE(chunk.conflicting_msa_observations[0].update_counts);
        CHECK(chunk.conflicting_msa_observations[1].update_counts);
        CHECK(chunk.haps == std::vector<int>{1});
        CHECK(chunk.phase_sets == std::vector<hts_pos_t>{1000});
    }
    SECTION("lost or demoted source candidates cannot gain an observation") {
        chunk.candidates[0].phase_set = 0;
        chunk.candidates[1].msa_verified = false;
        CHECK_FALSE(apply_pending_msa_observations(chunk));
        CHECK(chunk.read_var_profile[0].alleles == std::vector<int>{0, 0});
        CHECK(chunk.candidates[0].counts.total_cov == 6);
        CHECK(chunk.candidates[1].counts.total_cov == 6);
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
    SECTION("missing SNP quality cannot provide corroboration") {
        chunk.read_var_profile[0].bam_base_qualities[0] = 255;
        CHECK_FALSE(corroborated_bam_block_flip(chunk, {0, 1, 2}, {3}));
    }
    SECTION("unknown BAM mapping quality cannot provide corroboration") {
        chunk.read_var_profile[0].bam_mapq = 255;
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

TEST_CASE("insertion edit identity includes the common haplotype background",
          "[msa][recovery][representation][insertion-equivalence]") {
    CHECK(insertion_edits_are_equivalent(100, "ACGT", 102, "GTAC", "AC"));
    CHECK(insertion_edits_are_equivalent(100, "acgt", 102, "gtac", "ac"));
    // The same insertions disagree on the genomic background and agree only
    // after the independently called common C->T substitution is included.
    CHECK_FALSE(insertion_edits_are_equivalent(100, "AT", 102, "AT", "AC"));
    CHECK(insertion_edits_are_equivalent(100, "AT", 102, "AT", "AT"));
    CHECK(insertion_edits_are_equivalent(100, "AC", 100, "ac", ""));
    CHECK_FALSE(insertion_edits_are_equivalent(100, "AC", 100, "AT", ""));
    CHECK_FALSE(insertion_edits_are_equivalent(102, "AC", 100, "AC", "AC"));
    CHECK_FALSE(insertion_edits_are_equivalent(0, "AC", 2, "AC", "AC"));
    CHECK_FALSE(insertion_edits_are_equivalent(100, "", 100, "", ""));
    CHECK_FALSE(insertion_edits_are_equivalent(100, "AC", 102, "ACA", "AC"));
    CHECK_FALSE(insertion_edits_are_equivalent(100, "AC", 102, "AC", "A"));
    CHECK_FALSE(insertion_edits_are_equivalent(100, "AN", 102, "AN", "AN"));
    CHECK_FALSE(insertion_edits_are_equivalent(100, "AC", 102, "AC", "AN"));
}

TEST_CASE("core phase-set joins preserve independent output-only read gauges",
          "[recovery][stitch][read-only-gauge]") {
    for (const bool flip : {false, true}) {
        PhasingChunk chunk;
        chunk.reads.resize(4);
        chunk.haps = {1, 2, 0, 0};
        chunk.phase_sets = {100, 200, 0, 0};
        chunk.gap_haps = {0, 0, 1, 2};
        chunk.gap_phase_sets = {0, 0, 100 + kGapFillPsOffset, 200 + kGapFillPsOffset};
        chunk.bam_fallback_haps = {0, 0, 2, 1};
        chunk.bam_fallback_phase_sets = {0, 0, 100 + kBamFallbackPsOffset,
                                        200 + kBamFallbackPsOffset};
        CandidateVariant left, right;
        left.key = ins_key(100, "AC");
        right.key = ins_key(200, "AC");
        left.phase_set = 100;
        right.phase_set = 200;
        left.hap_to_cons_alle = {-1, 0, 1};
        right.hap_to_cons_alle = {-1, flip ? 1 : 0, flip ? 0 : 1};
        chunk.candidates = {left, right};
        REQUIRE(merge_phase_sets_in_place(chunk, 100, 200, flip));
        CHECK(chunk.candidates[0].hap_to_cons_alle == chunk.candidates[1].hap_to_cons_alle);
        CHECK(chunk.candidates[1].phase_set == 100);
        CHECK(chunk.phase_sets == std::vector<hts_pos_t>{100, 100, 0, 0});
        CHECK(chunk.haps[1] == (flip ? 1 : 2));
        CHECK(chunk.gap_haps == std::vector<int>{0, 0, 1, 2});
        CHECK(chunk.gap_phase_sets == std::vector<hts_pos_t>{
            0, 0, 100 + kGapFillPsOffset, 200 + kGapFillPsOffset});
        CHECK(chunk.bam_fallback_haps == std::vector<int>{0, 0, 2, 1});
        CHECK(chunk.bam_fallback_phase_sets == std::vector<hts_pos_t>{
            0, 0, 100 + kBamFallbackPsOffset, 200 + kBamFallbackPsOffset});
    }
}


TEST_CASE("complete BAM evidence retains context independently of graph ownership",
          "[recovery][complete-block]") {
    PhasingChunk source;
    for (int ci = 0; ci < 4; ++ci) {
        CandidateVariant site;
        site.key.pos = 100 + ci * 100;
        site.key.type = VariantType::Snp;
        site.phase_set = ci < 2 ? 100 : 300;
        site.hap_to_cons_alle = {-1, 0, 1};
        site.counts.category = VariantCategory::CleanHetSnp;
        source.candidates.push_back(site);
    }
    for (int ri = 0; ri < 20; ++ri) {
        ReadRecord read;
        read.qname = "molecule-" + std::to_string(ri);
        read.mapq = 60;
        source.reads.push_back(std::move(read));
        ReadVariantProfile profile;
        profile.start_var_idx = 0;
        profile.end_var_idx = 3;
        // Only the context sites call both blocks. The live gap rows (1,2)
        // cannot establish this connection without the independent snapshot.
        profile.alleles = {ri % 2, -1, -1, ri % 2};
        profile.bam_base_qualities = {40, 0, 0, 40};
        source.read_var_profile.push_back(std::move(profile));
    }
    RecoveryPhaseGauge gauge;
    retain_recovery_bam_evidence(source, {{100, 101}, {300, 301}}, gauge);
    PhasingChunk graph;
    graph.candidates = {source.candidates[1], source.candidates[2]};
    graph.candidates[0].phase_set = 101;
    graph.candidates[1].phase_set = 301;
    const auto flip = [&]() {
        return complete_recovery_block_flip(graph, gauge, 101, 301, {0}, {1}, 5);
    };
    REQUIRE(gauge.bam_sites.size() == 4);
    REQUIRE(gauge.bam_reads.size() == 20);
    REQUIRE(flip().has_value());
    CHECK_FALSE(*flip());
    CHECK(graph.candidates.size() == 2);
    CHECK(graph.candidates[0].phase_set == 101);
    CHECK(source.candidates[0].phase_set == 100);
    SECTION("clean labels cannot replace physical allele certificates") {
        for (auto& read : gauge.bam_reads) read.base_qualities.clear();
        CHECK_FALSE(flip());
    }
    SECTION("low-quality SNP bases cannot certify a whole-block join") {
        for (auto& read : gauge.bam_reads)
            std::fill(read.base_qualities.begin(), read.base_qualities.end(), 5);
        CHECK_FALSE(flip());
    }
    SECTION("low MAPQ calls remain observations but cannot certify a join") {
        for (auto& read : gauge.bam_reads) read.mapq = 5;
        CHECK_FALSE(flip());
    }
    SECTION("missing base qualities cannot masquerade as high quality") {
        for (auto& read : gauge.bam_reads)
            std::fill(read.base_qualities.begin(), read.base_qualities.end(), 255);
        CHECK_FALSE(flip());
    }
    SECTION("physical calls cannot contradict the block's existing read HP gauge") {
        for (const auto& saved : gauge.bam_reads) {
            ReadRecord read;
            read.qname = saved.qname;
            read.mapq = saved.mapq;
            graph.reads.push_back(std::move(read));
            graph.haps.push_back(2 - saved.observations.front().second);
            graph.phase_sets.push_back(graph.phase_sets.size() < 10 ? 101 : 301);
        }
        CHECK_FALSE(flip());
        for (int& hap : graph.haps) hap = 3 - hap;
        REQUIRE(flip().has_value());
        CHECK_FALSE(*flip());
    }
    SECTION("one agreeing block cannot hide its neighbor's contradictory HP tags") {
        for (size_t ri = 0; ri < gauge.bam_reads.size(); ++ri) {
            const auto& saved = gauge.bam_reads[ri];
            ReadRecord read;
            read.qname = saved.qname;
            graph.reads.push_back(std::move(read));
            const int allele = saved.observations.front().second;
            graph.haps.push_back(ri < 10 ? 2 - allele : 1 + allele);
            graph.phase_sets.push_back(ri < 10 ? 101 : 301);
        }
        CHECK_FALSE(flip());
    }
    SECTION("independent gauges can need a flip") {
        for (auto& site : gauge.bam_sites)
            if (site.phase_set == 301) std::swap(site.hap1_allele, site.hap2_allele);
        std::swap(graph.candidates[1].hap_to_cons_alle[1],
                  graph.candidates[1].hap_to_cons_alle[2]);
        REQUIRE(flip().has_value());
        CHECK(*flip());
    }
    SECTION("one haplotype cannot certify a whole-block connection") {
        for (auto& read : gauge.bam_reads)
            for (auto& observation : read.observations) observation.second = 0;
        CHECK_FALSE(flip());
    }
    SECTION("balanced evidence abstains") {
        for (size_t ri = 0; ri < gauge.bam_reads.size(); ++ri)
            if (ri % 2 == 0) gauge.bam_reads[ri].observations.back().second ^= 1;
        CHECK_FALSE(flip());
    }
    SECTION("unknown and low mapping qualities do not vote") {
        for (size_t ri = 0; ri < gauge.bam_reads.size(); ++ri)
            gauge.bam_reads[ri].mapq = ri % 2 == 0 ? 255 : 4;
        CHECK_FALSE(flip());
    }
    SECTION("alternative alignments do not multiply support") {
        for (auto& read : source.reads) read.qname = "one-molecule";
        retain_recovery_bam_evidence(source, {{100, 101}, {300, 301}}, gauge);
        CHECK(gauge.bam_reads.empty());
        CHECK_FALSE(flip());
    }
    SECTION("flank-only blocks preserve reads without aliasing graph PS labels") {
        source.candidates[0].hap_to_cons_alle = {-1, 0, 0};
        retain_recovery_bam_evidence(source, {{100, 101}}, gauge);
        REQUIRE(gauge.bam_sites.size() == 3);
        CHECK(gauge.bam_sites[0].key.pos == 200);
        CHECK(gauge.bam_sites[0].phase_set == 101);
        CHECK(gauge.bam_sites[0].source_phase_set == 100);
        CHECK(gauge.bam_sites[1].phase_set == 0);
        CHECK(gauge.bam_sites[1].source_phase_set == 300);
        CHECK(gauge.bam_sites[2].phase_set == 0);
        CHECK(gauge.bam_sites[2].source_phase_set == 300);
        REQUIRE(gauge.bam_reads.size() == 20);
        for (const auto& read : gauge.bam_reads) {
            REQUIRE(read.observations.size() == 1);
            CHECK(read.observations.front().first == 2);
        }
        // A raw source coordinate equal to a graph PS does not grant ownership.
        CHECK_FALSE(complete_recovery_block_flip(
            graph, gauge, 300, 301, {0}, {1}, 5));
    }
    SECTION("clean SNPs cannot be outvoted by repetitive indels") {
        RecoveryBamSite indel = gauge.bam_sites.front();
        indel.clean_snp = false;
        indel.key.type = VariantType::Deletion;
        for (int i = 0; i < 10; ++i) {
            indel.phase_set = i < 5 ? 101 : 301;
            gauge.bam_sites.push_back(indel);
        }
        for (auto& read : gauge.bam_reads) {
            const int hap = read.observations.front().second;
            for (size_t ci = 4; ci < gauge.bam_sites.size(); ++ci)
                read.observations.emplace_back(ci, ci < 9 ? hap : 1 - hap);
        }
        REQUIRE(flip().has_value());
        CHECK_FALSE(*flip());
    }
    SECTION("snapshot survives source destruction and numeric PS collisions") {
        source.candidates.clear();
        source.reads.clear();
        source.read_var_profile.clear();
        graph.candidates[0].phase_set = 100;
        REQUIRE(flip().has_value());
        CHECK_FALSE(*flip());
        CHECK_FALSE(complete_recovery_block_flip(graph, gauge, 101, 101, {}, {}, 5));
    }
    SECTION("graph flanks are scored by molecule name in their current gauge") {
        for (const auto& source_read : source.reads) {
            ReadRecord read;
            read.qname = source_read.qname;
            read.mapq = source_read.mapq;
            graph.reads.push_back(std::move(read));
        }
        std::sort(graph.reads.begin(), graph.reads.end(),
                  [](const ReadRecord& a, const ReadRecord& b) { return a.qname < b.qname; });
        graph.read_var_profile.clear();
        for (const auto& read : graph.reads) {
            const int hap = std::stoi(read.qname.substr(9)) % 2;
            ReadVariantProfile profile;
            profile.start_var_idx = 0;
            profile.end_var_idx = 1;
            profile.alleles = {hap, -1};
            profile.bam_alleles = {hap, -1};
            profile.bam_base_qualities = {40, 0};
            graph.read_var_profile.push_back(std::move(profile));
        }
        const auto graph_flip = complete_recovery_block_flip(
            graph, gauge, 900, 301, {0}, {1}, 5);
        REQUIRE(graph_flip.has_value());
        CHECK_FALSE(*graph_flip);
        std::swap(graph.candidates[0].hap_to_cons_alle[1],
                  graph.candidates[0].hap_to_cons_alle[2]);
        const auto reversed = complete_recovery_block_flip(
            graph, gauge, 900, 301, {0}, {1}, 5);
        REQUIRE(reversed.has_value());
        CHECK(*reversed);
        for (auto& profile : graph.read_var_profile)
            profile.bam_alleles[0] = 1 - profile.alleles[0];
        CHECK_FALSE(complete_recovery_block_flip(
            graph, gauge, 900, 301, {0}, {1}, 5));
        for (auto& profile : graph.read_var_profile)
            profile.bam_alleles.clear();
        CHECK_FALSE(complete_recovery_block_flip(
            graph, gauge, 900, 301, {0}, {1}, 5));
    }
}

TEST_CASE("complete BAM blocks stitch atomically without injecting flank rows",
          "[recovery][complete-block][stitch]") {
    PhasingChunk chunk;
    for (int ci = 0; ci < 4; ++ci) {
        CandidateVariant site;
        site.key.pos = 500 + ci * 1000;
        site.key.type = VariantType::Snp;
        site.counts.category = VariantCategory::CleanHetSnp;
        site.phase_set = 100 + ci * 100;
        site.hap_to_cons_alle = {-1, ci == 3 ? 1 : 0, ci == 3 ? 0 : 1};
        site.bam_injected = ci == 1 || ci == 2;
        chunk.candidates.push_back(site);
    }
    RecoveryPhaseGauge gauge;
    gauge.beg = 0;
    gauge.end = 4000;
    gauge.imported_phase_sets = {200, 300};
    RecoveryBlockGaugeVote left;
    left.graph_phase_set = 100;
    left.bam_phase_set = 200;
    left.counts = {{{10, 0}, {0, 10}}};
    left.shared_candidate_same = 8;
    RecoveryBlockGaugeVote right;
    right.graph_phase_set = 400;
    right.bam_phase_set = 300;
    right.counts = {{{0, 10}, {10, 0}}};
    right.shared_candidate_cross = 8;
    gauge.block_votes = {left, right};
    PhasingChunk source;
    source.candidates = {chunk.candidates[1], chunk.candidates[2]};
    for (int ri = 0; ri < 40; ++ri) {
        ReadRecord read;
        read.qname = "bridge-" + std::to_string(ri);
        read.mapq = 60;
        source.reads.push_back(std::move(read));
        ReadVariantProfile profile;
        profile.start_var_idx = 0;
        profile.end_var_idx = 1;
        profile.alleles = {ri % 2, ri % 2};
        profile.bam_base_qualities = {40, 40};
        source.read_var_profile.push_back(std::move(profile));
    }
    for (size_t ri = 0; ri < source.reads.size(); ++ri) {
        ReadRecord read;
        read.qname = source.reads[ri].qname;
        read.mapq = 60;
        chunk.reads.push_back(std::move(read));
    }
    std::sort(chunk.reads.begin(), chunk.reads.end(),
              [](const ReadRecord& a, const ReadRecord& b) { return a.qname < b.qname; });
    for (const auto& read : chunk.reads) {
        const int allele = std::stoi(read.qname.substr(7)) % 2;
        ReadVariantProfile profile;
        profile.start_var_idx = 0;
        profile.end_var_idx = 3;
        profile.alleles = {allele, -1, -1, allele};
        profile.bam_alleles = profile.alleles;
        profile.bam_base_qualities = {40, 0, 0, 40};
        chunk.read_var_profile.push_back(std::move(profile));
        const int original_i = std::stoi(read.qname.substr(7));
        chunk.phase_sets.push_back(100 + (original_i / 10) * 100);
        chunk.haps.push_back(original_i < 30 ? 1 + allele : 2 - allele);
    }
    retain_recovery_bam_evidence(source, {{200, 200}, {300, 300}}, gauge);
    Options opts;
    const std::unordered_map<hts_pos_t, bool> complete{{200, true}, {300, true}};
    const RecoverySeam seam{500, 3500, 100, 400};
    SECTION("complete source observations close the full chain") {
        CHECK(stitch_complete_recovery_phase_blocks(
            chunk, {seam}, {gauge}, opts, complete, {}) == 3);
        REQUIRE(chunk.candidates.size() == 4);
        for (const auto& site : chunk.candidates) {
            CHECK(site.phase_set == 100);
            CHECK(site.hap_to_cons_alle == std::array<int, 3>{-1, 0, 1});
        }
    }
    SECTION("an earlier attachment flips live blocks but keeps source evidence immutable") {
        chunk.candidates.front().hap_to_cons_alle = {-1, 1, 0};
        for (size_t ri = 0; ri < chunk.haps.size(); ++ri)
            if (chunk.phase_sets[ri] == 100) chunk.haps[ri] = 3 - chunk.haps[ri];
        gauge.block_votes.front().counts = {{{0, 10}, {10, 0}}};
        gauge.block_votes.front().shared_candidate_same = 0;
        gauge.block_votes.front().shared_candidate_cross = 8;
        CHECK(stitch_complete_recovery_phase_blocks(
            chunk, {seam}, {gauge}, opts, complete, {}) == 3);
        for (const auto& site : chunk.candidates) {
            CHECK(site.phase_set == 100);
            CHECK(site.hap_to_cons_alle == std::array<int, 3>{-1, 1, 0});
        }
        CHECK(gauge.bam_sites.front().hap1_allele == 0);
    }
    SECTION("an internal source cut cannot be hidden by a large block vote") {
        const std::unordered_map<hts_pos_t, bool> broken{{200, false}, {300, true}};
        stitch_complete_recovery_phase_blocks(chunk, {seam}, {gauge}, opts, broken, {});
        CHECK(chunk.candidates.front().phase_set != chunk.candidates.back().phase_set);
    }
    SECTION("local flank votes cannot certify an internally unsupported graph block") {
        CandidateVariant context = chunk.candidates.front();
        context.key.pos = 100;
        chunk.candidates.insert(chunk.candidates.begin(), context);
        for (auto& profile : chunk.read_var_profile) {
            profile.end_var_idx = 4;
            profile.alleles.insert(profile.alleles.begin(), -1);
            profile.bam_alleles.insert(profile.bam_alleles.begin(), -1);
            profile.bam_base_qualities.insert(profile.bam_base_qualities.begin(), 0);
        }
        stitch_complete_recovery_phase_blocks(chunk, {seam}, {gauge}, opts, complete, {});
        CHECK(chunk.candidates.front().phase_set != chunk.candidates.back().phase_set);
    }
    SECTION("overlapping reads certify a long graph flank without one end-to-end read") {
        CandidateVariant context = chunk.candidates.front();
        context.key.pos = 100;
        chunk.candidates.insert(chunk.candidates.begin(), context);
        context.key.pos = 300;
        chunk.candidates.insert(chunk.candidates.begin() + 1, context);
        for (auto& profile : chunk.read_var_profile) {
            profile.end_var_idx += 2;
            profile.alleles.insert(profile.alleles.begin(), 2, -1);
            profile.bam_alleles.insert(profile.bam_alleles.begin(), 2, -1);
            profile.bam_base_qualities.insert(profile.bam_base_qualities.begin(), 2, 0);
        }
        // A--B and B--C are supported by different molecules on both
        // haplotypes. No molecule calls both ends of the graph block.
        for (int ri = 15; ri >= 0; --ri) {
            ReadRecord read;
            read.qname = "aaa-context-" + std::to_string(ri);
            read.mapq = 60;
            chunk.reads.insert(chunk.reads.begin(), std::move(read));
            ReadVariantProfile profile;
            profile.start_var_idx = 0;
            profile.end_var_idx = 5;
            profile.alleles.assign(6, -1);
            const int first = ri < 8 ? 0 : 1;
            profile.alleles[first] = ri % 2;
            profile.alleles[first + 1] = ri % 2;
            profile.bam_alleles = profile.alleles;
            chunk.read_var_profile.insert(chunk.read_var_profile.begin(), profile);
            chunk.phase_sets.insert(chunk.phase_sets.begin(), 100);
            chunk.haps.insert(chunk.haps.begin(), 1 + ri % 2);
        }
        SECTION("connected chains preserve the full flank orientation") {
            CHECK(stitch_complete_recovery_phase_blocks(
                chunk, {seam}, {gauge}, opts, complete, {}) == 3);
            for (const auto& site : chunk.candidates) {
                CHECK(site.phase_set == 100);
                CHECK(site.hap_to_cons_alle == std::array<int, 3>{-1, 0, 1});
            }
        }
        SECTION("interleaved disconnected components cannot masquerade as one path") {
            context.key.pos = 700;
            chunk.candidates.insert(chunk.candidates.begin() + 3, context);
            for (auto& profile : chunk.read_var_profile) {
                profile.end_var_idx += 1;
                profile.alleles.insert(profile.alleles.begin() + 3, -1);
                profile.bam_alleles.insert(profile.bam_alleles.begin() + 3, -1);
                if (!profile.bam_base_qualities.empty())
                    profile.bam_base_qualities.insert(profile.bam_base_qualities.begin() + 3, 0);
            }
            // A--C and B--D cross every coordinate cut, but no molecule
            // establishes the orientation between those two components.
            for (size_t ri = 0; ri < 16; ++ri) {
                auto& alleles = chunk.read_var_profile[ri].alleles;
                alleles.assign(7, -1);
                const int first = ri < 8 ? 0 : 1;
                alleles[first] = static_cast<int>(ri % 2);
                alleles[first + 2] = static_cast<int>(ri % 2);
            }
            CHECK(stitch_complete_recovery_phase_blocks(
                chunk, {seam}, {gauge}, opts, complete, {}) == 0);
        }
        SECTION("a disconnected cut cannot inherit a certificate from its neighbors") {
            for (size_t ri = 8; ri < 16; ++ri)
                chunk.read_var_profile[ri].alleles[2] = -1;
            CHECK(stitch_complete_recovery_phase_blocks(
                chunk, {seam}, {gauge}, opts, complete, {}) == 0);
        }
        SECTION("an internal switch cannot be hidden by correct boundary votes") {
            std::swap(chunk.candidates[1].hap_to_cons_alle[1],
                      chunk.candidates[1].hap_to_cons_alle[2]);
            CHECK(stitch_complete_recovery_phase_blocks(
                chunk, {seam}, {gauge}, opts, complete, {}) == 0);
        }
        SECTION("unknown mapping quality cannot certify internal connectivity") {
            for (size_t ri = 0; ri < 16; ++ri) chunk.reads[ri].mapq = 255;
            CHECK(stitch_complete_recovery_phase_blocks(
                chunk, {seam}, {gauge}, opts, complete, {}) == 0);
        }
    }
    SECTION("an imported source certificate does not cover unsupported absorbed rows") {
        CandidateVariant context = chunk.candidates[1];
        context.key.pos = 1700;
        chunk.candidates.insert(chunk.candidates.begin() + 2, context);
        for (auto& profile : chunk.read_var_profile) {
            profile.end_var_idx += 1;
            profile.alleles.insert(profile.alleles.begin() + 2, -1);
            profile.bam_alleles.insert(profile.bam_alleles.begin() + 2, -1);
            profile.bam_base_qualities.insert(profile.bam_base_qualities.begin() + 2, 0);
        }
        CHECK(stitch_complete_recovery_phase_blocks(
            chunk, {seam}, {gauge}, opts, complete, {}) == 0);
        CHECK(chunk.candidates.front().phase_set != chunk.candidates.back().phase_set);
    }
    SECTION("a homozygous row cannot fabricate an imported block anchor") {
        CandidateVariant hom = chunk.candidates[1];
        hom.key.pos = 2000;
        hom.phase_set = 500;
        hom.hap_to_cons_alle = {-1, 0, 0};
        chunk.candidates.insert(chunk.candidates.begin() + 2, hom);
        for (auto& profile : chunk.read_var_profile) {
            profile.end_var_idx = 4;
            profile.alleles.insert(profile.alleles.begin() + 2, -1);
            profile.bam_alleles.insert(profile.bam_alleles.begin() + 2, -1);
            profile.bam_base_qualities.insert(profile.bam_base_qualities.begin() + 2, 0);
        }
        gauge.imported_phase_sets.push_back(500);
        CHECK(stitch_complete_recovery_phase_blocks(
            chunk, {seam}, {gauge}, opts, complete, {}) == 3);
        CHECK(chunk.candidates.front().phase_set == chunk.candidates.back().phase_set);
        CHECK(chunk.candidates[2].hap_to_cons_alle == std::array<int, 3>{-1, 0, 0});
    }
    SECTION("a retry without the old block cannot hide earlier complete evidence") {
        RecoveryPhaseGauge retry;
        retry.beg = gauge.beg;
        retry.end = gauge.end;
        CHECK(stitch_complete_recovery_phase_blocks(
            chunk, {seam}, {retry, gauge}, opts, complete, {}) == 3);
        CHECK(chunk.candidates.front().phase_set == chunk.candidates.back().phase_set);
    }
    SECTION("missing current read gauge cannot certify a new whole-block join") {
        std::fill(chunk.phase_sets.begin(), chunk.phase_sets.end(), 0);
        CHECK(stitch_complete_recovery_phase_blocks(
            chunk, {seam}, {gauge}, opts, complete, {}) == 0);
        CHECK(chunk.candidates.front().phase_set != chunk.candidates.back().phase_set);
    }
    SECTION("without the source matrix the live gap sites lack a connection") {
        gauge.bam_sites.clear();
        gauge.bam_reads.clear();
        stitch_complete_recovery_phase_blocks(chunk, {seam}, {gauge}, opts, complete, {});
        CHECK(chunk.candidates.front().phase_set != chunk.candidates.back().phase_set);
    }
}
