/*
 * Compile this against longcallD's original assign_hap.c, then compare the
 * BAM-mode C++ ports at their function boundaries. No reference algorithm is
 * transcribed into this test.
 */
#define CATCH_CONFIG_MAIN
#include "../third_party/catch2/catch.hpp"

#include "collect_phase.hpp"
#include "collect_phase_noisy.hpp"
#include "test_upstream_phase_bridge.h"

#include <array>
#include <vector>

using namespace pgphase_collect;

static CandidateVariant candidate(int category, int n_alleles, int homopolymer) {
    CandidateVariant var;
    var.lcd_var_i_to_cate = static_cast<uint32_t>(category);
    var.counts.n_uniq_alles = n_alleles;
    var.counts.alle_covs.assign(static_cast<size_t>(n_alleles), 0);
    var.is_homopolymer_indel = homopolymer != 0;
    return var;
}

TEST_CASE("MSA homopolymer detector matches original C in BAM mode") {
    for (const std::string ref : {"AAAAAAAA", "aaaaaaaa", "CCCCCCCC",
                                  "cccccccc", "TTTTTTTT", "ACGTACGT"}) {
        PhasingChunk chunk;
        chunk.ref_beg = 100;
        chunk.ref_seq = ref;
        for (int offset : {0, 1, 2}) {
            for (const std::string alt : {"A", "C", "G", "T", "AA", "CC"}) {
                CAPTURE(ref, offset, alt);
                const int expected = upstream_msa_homopolymer(
                    ref.c_str(), offset, 1, 0, alt.c_str());
                CHECK(var_is_homopolymer_indel(
                    chunk, 100 + offset, VariantType::Insertion, 0, alt, true) ==
                    (expected != 0));
            }
            for (int len : {1, 3, 6}) {
                CAPTURE(ref, offset, len);
                const int expected = upstream_msa_homopolymer(
                    ref.c_str(), offset, 0, len, "");
                CHECK(var_is_homopolymer_indel(
                    chunk, 100 + offset, VariantType::Deletion, len, "", true) ==
                    (expected != 0));
            }
        }
    }
}

TEST_CASE("var_init_hap_profile_cons_allele matches original C") {
    const int categories[] = {kCandCleanHetSnp, kCandCleanHetIndel,
                              kCandCleanHom, kCandNoisyCandHet, kCandNoisyCandHom};
    for (int category : categories)
    for (int n_alleles : {2, 3})
    for (int is_ont : {0, 1})
    for (int homopolymer : {0, 1})
    for (int had_profile : {0, 1})
    for (int cov0 : {0, 1, 3})
    for (int cov1 : {0, 1, 3}) {
        int coverage[3] = {cov0, cov1, n_alleles == 3 ? 2 : 0};
        int expected_cons[3] = {-1, 1, 0};
        int expected_profile[3][3] = {{4, 3, 2}, {3, 2, 1}, {2, 1, 0}};
        upstream_var_init(is_ont, homopolymer, category, n_alleles, coverage,
                          had_profile, expected_cons, expected_profile);

        CandidateVariant var = candidate(category, n_alleles, homopolymer);
        var.counts.alle_covs.assign(coverage, coverage + n_alleles);
        var.hap_to_cons_alle = {-1, 1, 0};
        if (had_profile) {
            const int initial_profile[3][3] = {{4, 3, 2}, {3, 2, 1}, {2, 1, 0}};
            for (int h = 0; h < 3; ++h)
                var.hap_to_alle_profile[h].assign(
                    initial_profile[h], initial_profile[h] + n_alleles);
        }
        CandidateTable variants = {var};
        var_init_hap_profile_cons_allele(is_ont != 0, variants, {0}, false);
        const CandidateVariant& got = variants[0];
        for (int h = 0; h < 3; ++h) {
            CHECK(got.hap_to_cons_alle[h] == expected_cons[h]);
            for (int a = 0; a < n_alleles; ++a)
                CHECK(got.hap_to_alle_profile[h][a] == expected_profile[h][a]);
        }
    }
}

TEST_CASE("update_var_hap_to_cons_alle matches original C") {
    for (int n_alleles : {2, 3})
    for (int is_ont : {0, 1})
    for (int homopolymer : {0, 1})
    for (int hap : {0, 1, 2})
    for (int p0 : {0, 1, 3, 5})
    for (int p1 : {0, 1, 3, 5})
    for (int p2 : {0, 1, 3, 5}) {
        int expected_profile[3][3] = {{0, 0, 0}, {p0, p1, p2}, {p2, p1, p0}};
        int expected_cons[3] = {0, 1, 0};
        upstream_update_cons(is_ont, homopolymer, kCandCleanHetSnp, n_alleles,
                             hap, expected_profile, expected_cons);
        CandidateVariant var = candidate(kCandCleanHetSnp, n_alleles, homopolymer);
        var.hap_to_cons_alle = {0, 1, 0};
        var.hap_to_alle_profile[1] = {p0, p1};
        var.hap_to_alle_profile[2] = {p2, p1};
        if (n_alleles == 3) {
            var.hap_to_alle_profile[1].push_back(p2);
            var.hap_to_alle_profile[2].push_back(p0);
        }
        update_var_hap_to_cons_alle(is_ont != 0, var, hap);
        for (int h = 0; h < 3; ++h)
            CHECK(var.hap_to_cons_alle[h] == expected_cons[h]);
    }
}

TEST_CASE("read_to_cons_allele_score matches original C in BAM mode") {
    const int categories[] = {kCandCleanHetSnp, kCandCleanHetIndel,
                              kCandCleanHom, kCandNoisyCandHet, kCandNoisyCandHom};
    for (int category : categories)
    for (int n_alleles : {2, 3})
    for (int cons1 : {-1, 0, 1})
    for (int cons2 : {-1, 0, 1})
    for (int hap : {1, 2})
    for (int allele = 0; allele < n_alleles; ++allele) {
        int expected_cons[3] = {-1, cons1, cons2};
        const int expected = upstream_read_score(hap, category, allele, expected_cons);
        CandidateVariant var = candidate(category, n_alleles, 0);
        var.hap_to_cons_alle = {-1, cons1, cons2};
        const int got = read_to_cons_allele_score(var, hap, allele, true, true, true);
        CHECK(got == expected);
        for (int h = 0; h < 3; ++h)
            CHECK(var.hap_to_cons_alle[h] == expected_cons[h]);
    }
}

TEST_CASE("init_assign_read_hap_based_on_cons_alle matches original C in BAM mode") {
    const std::array<std::array<int, 2>, 9> category_pairs = {{
        {{kCandCleanHetSnp, kCandCleanHetSnp}},
        {{kCandCleanHetSnp, kCandCleanHetIndel}},
        {{kCandCleanHetSnp, kCandCleanHom}},
        {{kCandCleanHetSnp, kCandNoisyCandHet}},
        {{kCandCleanHetSnp, kCandNoisyCandHom}},
        {{kCandNoisyCandHet, kCandNoisyCandHet}},
        {{kCandCleanHom, kCandCleanHom}},
        {{kCandCleanHetIndel, kCandNoisyCandHet}},
        {{kCandCleanHetSnp, kLongcalldLowCovVar}}
    }};
    const std::array<std::array<int, 2>, 6> cons_pairs = {{
        {{-1, -1}}, {{-1, 0}}, {{0, -1}},
        {{0, 1}}, {{1, 0}}, {{1, 1}}
    }};
    const std::array<std::array<int, 2>, 5> allele_pairs = {{
        {{0, 0}}, {{0, 1}}, {{1, 0}}, {{1, 1}}, {{-1, 1}}
    }};
    const std::array<std::array<int, 2>, 3> snp_pairs = {{
        {{1, 1}}, {{1, 0}}, {{0, 1}}
    }};
    const std::array<std::array<int, 2>, 3> hp_pairs = {{
        {{0, 0}}, {{1, 0}}, {{0, 1}}
    }};
    for (const auto& categories : category_pairs)
    for (const auto& cons_a : cons_pairs)
    for (const auto& cons_b : cons_pairs)
    for (const auto& alleles : allele_pairs)
    for (const auto& snps : snp_pairs)
    for (const auto& hp : hp_pairs) {
        int upstream_cons[3][3] = {{-1, cons_a[0], cons_a[1]},
                                   {-1, cons_b[0], cons_b[1]},
                                   {-1, -1, -1}};
        int expected_agree = 99, expected_conflict = 99;
        const int upstream_hap = upstream_init_read_hap(
            2, categories.data(), snps.data(), hp.data(), alleles.data(),
            upstream_cons, kCandGermlineVarCate,
            &expected_agree, &expected_conflict);

        PhasingChunk chunk;
        for (int i = 0; i < 2; ++i) {
            CandidateVariant var = candidate(categories[i], 2, hp[i]);
            var.key.type = snps[i] ? VariantType::Snp : VariantType::Insertion;
            var.hap_to_cons_alle = i == 0
                ? std::array<int, 3>{-1, cons_a[0], cons_a[1]}
                : std::array<int, 3>{-1, cons_b[0], cons_b[1]};
            chunk.candidates.push_back(var);
        }
        chunk.reads.resize(1);
        ReadVariantProfile profile;
        profile.read_id = 0;
        profile.start_var_idx = 0;
        profile.end_var_idx = 1;
        profile.alleles = {alleles[0], alleles[1]};
        chunk.read_var_profile.push_back(profile);
        const int got_hap = init_assign_read_hap_based_on_cons_alle(
            chunk, 0, kCandGermlineVarCate, std::nullopt, true, true, true);
        CAPTURE(categories[0], categories[1], cons_a[0], cons_a[1],
                cons_b[0], cons_b[1], alleles[0], alleles[1],
                snps[0], snps[1], hp[0], hp[1]);
        CHECK(got_hap == upstream_hap);
        CHECK(chunk.reads[0].n_clean_agree_snps == expected_agree);
        CHECK(chunk.reads[0].n_clean_conflict_snps == expected_conflict);
        for (int i = 0; i < 2; ++i)
        for (int h = 0; h < 3; ++h)
            CHECK(chunk.candidates[i].hap_to_cons_alle[h] == upstream_cons[i][h]);
    }
}

TEST_CASE("iter_update_var_hap_cons_phase_set matches original C in BAM mode") {
    for (int n_vars : {2, 3})
    for (int hap_mode : {0, 1, 2})
    for (int n_reads = 0; n_reads <= 4; ++n_reads)
    for (int conflict_mask_1 = 0; conflict_mask_1 < (1 << n_reads); ++conflict_mask_1)
    for (int conflict_mask_2 = 0;
         conflict_mask_2 < (n_vars == 3 ? (1 << n_reads) : 1);
         ++conflict_mask_2) {
        int haps[4] = {1, 1, 1, 1};
        for (int i = 0; i < n_reads; ++i)
            haps[i] = hap_mode == 0 ? 1 : (hap_mode == 1 ? 2 : 1 + (i & 1));
        int alleles[4][3] = {{0, 0, 0}, {0, 0, 0},
                             {0, 0, 0}, {0, 0, 0}};
        for (int i = 0; i < n_reads; ++i) {
            alleles[i][1] = (conflict_mask_1 >> i) & 1;
            alleles[i][2] = (conflict_mask_2 >> i) & 1;
        }
        int expected_cons[3][3] = {{-1, 0, 1}, {-1, 0, 1}, {-1, 0, 1}};
        int expected_ps[3] = {0, 0, 0};
        const int expected_changed = upstream_phase_link(
            n_vars, n_reads, haps, alleles, expected_cons, expected_ps);

        PhasingChunk chunk;
        for (int vi = 0; vi < n_vars; ++vi) {
            CandidateVariant var = candidate(kCandCleanHetSnp, 2, 0);
            var.key.pos = 100 * (vi + 1);
            var.key.type = VariantType::Snp;
            var.hap_to_cons_alle = {-1, 0, 1};
            chunk.candidates.push_back(var);
        }
        chunk.reads.resize(static_cast<size_t>(n_reads));
        chunk.haps.assign(haps, haps + n_reads);
        chunk.read_var_cr.reset(cr_init());
        for (int i = 0; i < n_reads; ++i) {
            ReadVariantProfile profile;
            profile.read_id = i;
            profile.start_var_idx = 0;
            profile.end_var_idx = n_vars - 1;
            profile.alleles.assign(alleles[i], alleles[i] + n_vars);
            chunk.read_var_profile.push_back(profile);
            cr_add(chunk.read_var_cr.get(), "cr", 0, n_vars, i);
        }
        cr_index(chunk.read_var_cr.get());
        Options opts;
        opts.upstream_assign_hap = true;
        std::vector<int> valid_idx;
        for (int i = 0; i < n_vars; ++i) valid_idx.push_back(i);
        const int got_changed = iter_update_var_hap_cons_phase_set(chunk, valid_idx, opts);
        CAPTURE(n_vars, hap_mode, n_reads, conflict_mask_1, conflict_mask_2);
        CHECK(got_changed == expected_changed);
        for (int vi = 0; vi < n_vars; ++vi) {
            CHECK(chunk.candidates[vi].phase_set == expected_ps[vi]);
            for (int h = 0; h < 3; ++h)
                CHECK(chunk.candidates[vi].hap_to_cons_alle[h] == expected_cons[vi][h]);
        }
    }
}

TEST_CASE("complete k-means pass matches original C on small BAM matrices") {
    for (int n_vars : {2, 3})
    for (int n_reads = 0; n_reads <= 4; ++n_reads)
    for (int mask1 = 0; mask1 < (1 << n_reads); ++mask1)
    for (int mask2 = 0; mask2 < (n_vars == 3 ? (1 << n_reads) : 1); ++mask2) {
        int alleles[4][3] = {{0, 0, 0}, {0, 0, 0},
                             {0, 0, 0}, {0, 0, 0}};
        for (int ri = 0; ri < n_reads; ++ri) {
            alleles[ri][1] = (mask1 >> ri) & 1;
            alleles[ri][2] = (mask2 >> ri) & 1;
        }
        int expected_cons[3][3] = {{-1, -1, -1},
                                   {-1, -1, -1},
                                   {-1, -1, -1}};
        int64_t expected_candidate_ps[3] = {0, 0, 0};
        int expected_haps[4] = {0, 0, 0, 0};
        int64_t expected_read_ps[4] = {0, 0, 0, 0};
        upstream_full_kmeans(n_vars, n_reads, alleles, 0, expected_cons,
                             expected_candidate_ps, expected_haps,
                             expected_read_ps);

        PhasingChunk chunk;
        for (int vi = 0; vi < n_vars; ++vi) {
            CandidateVariant var = candidate(kCandCleanHetSnp, 2, 0);
            var.key.pos = 100 * (vi + 1);
            var.key.type = VariantType::Snp;
            var.counts.total_cov = n_reads;
            for (int ri = 0; ri < n_reads; ++ri)
                ++var.counts.alle_covs[alleles[ri][vi]];
            chunk.candidates.push_back(var);
        }
        chunk.reads.resize(static_cast<size_t>(n_reads));
        chunk.read_var_cr.reset(cr_init());
        for (int ri = 0; ri < n_reads; ++ri) {
            ReadVariantProfile profile;
            profile.read_id = ri;
            profile.start_var_idx = 0;
            profile.end_var_idx = n_vars - 1;
            profile.alleles.assign(alleles[ri], alleles[ri] + n_vars);
            chunk.read_var_profile.push_back(profile);
            cr_add(chunk.read_var_cr.get(), "cr", 0, n_vars, ri);
        }
        cr_index(chunk.read_var_cr.get());
        Options opts;
        opts.msa_sites_vote_without_gap_link = true;
        opts.infer_complement_at_multiallelic = true;
        opts.upstream_read_scoring = true;
        opts.upstream_assign_hap = true;
        assign_hap_based_on_germline_het_vars_kmeans(
            chunk, opts, kCandGermlineClean);
        CAPTURE(n_vars, n_reads, mask1, mask2);
        for (int vi = 0; vi < n_vars; ++vi) {
            CHECK(chunk.candidates[vi].phase_set == expected_candidate_ps[vi]);
            for (int h = 0; h < 3; ++h) {
                CAPTURE(vi, h);
                CHECK(chunk.candidates[vi].hap_to_cons_alle[h] == expected_cons[vi][h]);
            }
            const int c1 = expected_cons[vi][1] < 0 ? 0 : expected_cons[vi][1];
            const int c2 = expected_cons[vi][2] < 0 ? 0 : expected_cons[vi][2];
            const int expected_hap_alt = c1 != 0 && c2 != 0 ? 3 : c1 != 0 ? 1 : c2 != 0 ? 2 : 0;
            const int expected_hap_ref = c1 != 0 && c2 == 0 ? 2 : c1 == 0 && c2 != 0 ? 1 : 0;
            CHECK(chunk.candidates[vi].hap_alt == expected_hap_alt);
            CHECK(chunk.candidates[vi].hap_ref == expected_hap_ref);
        }
        for (int ri = 0; ri < n_reads; ++ri) {
            CHECK(chunk.haps[ri] == expected_haps[ri]);
            CHECK(chunk.phase_sets[ri] == expected_read_ps[ri]);
        }
    }
}
