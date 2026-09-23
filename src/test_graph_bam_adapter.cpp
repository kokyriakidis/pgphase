#include "graph_bam_adapter.hpp"

#include "collect_phase.hpp"

#include <algorithm>
#include <iostream>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

using namespace pgphase_collect;

static bool check(bool cond, const std::string& msg) {
    if (!cond) std::cerr << "FAIL: " << msg << "\n";
    return cond;
}

static const ReadVariantProfile* profile_for_read(const PhasingChunk& chunk,
                                                  const std::string& qname) {
    for (const ReadRecord& read : chunk.reads) {
        if (read.qname != qname) continue;
        const int read_i = static_cast<int>(&read - chunk.reads.data());
        for (const ReadVariantProfile& profile : chunk.read_var_profile) {
            if (profile.read_id == read_i) {
                return &profile;
            }
        }
    }
    return nullptr;
}

static GraphSite make_site(const std::string& id,
                           hts_pos_t pos,
                           const std::string& allele0,
                           const std::string& allele1) {
    GraphSite site;
    site.chrom = "chr1";
    site.ref_contig = "chr1";
    site.pos = pos;
    site.ref_beg = pos;
    site.ref_end = pos + 1;
    site.id = id;
    site.allele_traversals = {allele0, allele1};
    site.allele_walks = {parse_graph_walk(allele0), parse_graph_walk(allele1)};
    site.skip_reason = graph_site_validation_skip_reason(site);
    site.eligible = site.skip_reason.empty();
    return site;
}

int main() {
    bool ok = true;

    GraphSiteCatalog catalog;
    catalog.sites.push_back(make_site("s1", 100, ">1>2>3", ">1>4>3"));
    catalog.sites.push_back(make_site("s2", 200, ">5>6>7", ">5>8>7"));

    std::vector<GraphReadAllele> rows = {
        {"s1", "chr1", 100, "read_a", 0},
        {"s2", "chr1", 200, "read_a", 0},
        {"s1", "chr1", 100, "read_b", 0},
        {"s2", "chr1", 200, "read_b", 0},
        {"s1", "chr1", 100, "read_c", 1},
        {"s2", "chr1", 200, "read_c", 1},
        {"s1", "chr1", 100, "read_d", 1},
        {"s2", "chr1", 200, "read_d", 1},
    };

    Options build_opts;
    build_opts.min_depth = 1;
    build_opts.min_alt_depth = 1;

    std::vector<GraphChunkBuildResult> chunks;
    chunks.push_back(build_graph_chunk(catalog.view_all(), rows, "chr1", 0, 300, 0, build_opts));
    ok &= check(chunks[0].chunk.candidates.size() == 2, "adapter builds two candidates");
    ok &= check(chunks[0].chunk.reads.size() == 4, "adapter builds four reads");
    ok &= check(chunks[0].chunk.read_var_profile.size() == 4, "adapter builds read profiles");
    ok &= check(chunks[0].chunk.read_var_cr != nullptr, "adapter builds read-var cgranges");
    ok &= check(chunks[0].chunk.candidates[0].phase_set == kUnsetCandidatePhaseSet,
                "graph candidate uses longcallD unset phase-set sentinel");
    ok &= check(chunks[0].chunk.phase_sets[0] == kUnphasedReadPhaseSet,
                "graph read uses longcallD unphased phase-set sentinel");

    Options opts;
    opts.read_technology = ReadTechnology::Hifi;
    phase_graph_chunks(chunks, opts);
    ok &= check(chunks[0].chunk.candidates[0].phase_set > 0, "BAM k-means phases graph candidate");
    ok &= check(chunks[0].chunk.haps.size() == 4, "BAM k-means assigns graph read hap vector");
    for (int hap : chunks[0].chunk.haps) {
        ok &= check(hap == 1 || hap == 2, "graph read hap assigned by BAM k-means");
    }

    GraphSiteCatalog conditional_catalog;
    conditional_catalog.sites.push_back(make_site("parent", 100, ">10>11>12", ">10>13>12"));
    GraphSite child = make_site("child", 110, ">20>21>22", ">20>23>22");
    child.parent = "parent";
    child.conditional_parent_alleles = {1};
    conditional_catalog.sites.push_back(child);
    std::vector<GraphReadAllele> conditional_rows = {
        {"child", "chr1", 110, "child_only", 1},
        {"parent", "chr1", 100, "with_parent", 1},
        {"child", "chr1", 110, "with_parent", 1},
    };
    Options no_af_opts;
    no_af_opts.min_depth = 1;
    no_af_opts.min_alt_depth = 1;
    no_af_opts.min_af = 0.0;
    no_af_opts.max_af = 1.0;
    GraphChunkBuildResult conditional_chunk =
        build_graph_chunk(conditional_catalog.view_all(), conditional_rows, "chr1", 0, 200, 0, no_af_opts);
    ok &= check(conditional_chunk.chunk.reads.size() == 1,
                "conditional child-only graph read is missing");
    ok &= check(!conditional_chunk.chunk.reads.empty() &&
                conditional_chunk.chunk.reads[0].qname == "with_parent",
                "conditional child with parent is retained");

    // --- AF / depth filter ---
    // Four sites: one het (passes), one hom-alt (AF=1.0 > max_af), one low-depth
    // (total=4 < min_depth=5), one low-alt (alt=1 < min_alt_depth=2).
    // A cross-read observing both the het site and the hom-alt site verifies that
    // its profile is truncated to the surviving site only.
    {
        GraphSiteCatalog af_catalog;
        af_catalog.sites.push_back(make_site("af_het",       100, ">1>2>3",    ">1>4>3"));
        af_catalog.sites.push_back(make_site("af_hom_alt",   200, ">5>6>7",    ">5>8>7"));
        af_catalog.sites.push_back(make_site("af_low_depth", 300, ">9>10>11",  ">9>12>11"));
        af_catalog.sites.push_back(make_site("af_low_alt",   400, ">13>14>15", ">13>16>15"));

        std::vector<GraphReadAllele> af_rows;
        // af_het: 4 ref + 4 alt → total=8, AF=0.5 → passes
        for (int i = 0; i < 4; ++i) {
            af_rows.push_back({"af_het", "chr1", 100, "het_ref_" + std::to_string(i), 0});
            af_rows.push_back({"af_het", "chr1", 100, "het_alt_" + std::to_string(i), 1});
        }
        // af_hom_alt: 6 alt reads → AF=1.0 > max_af → filtered
        for (int i = 0; i < 6; ++i)
            af_rows.push_back({"af_hom_alt", "chr1", 200, "hom_" + std::to_string(i), 1});
        // af_low_depth: 2 ref + 2 alt → total=4 < min_depth=5 → filtered
        for (int i = 0; i < 2; ++i) {
            af_rows.push_back({"af_low_depth", "chr1", 300, "ld_ref_" + std::to_string(i), 0});
            af_rows.push_back({"af_low_depth", "chr1", 300, "ld_alt_" + std::to_string(i), 1});
        }
        // af_low_alt: 6 ref + 1 alt → alt=1 < min_alt_depth=2 → filtered
        for (int i = 0; i < 6; ++i)
            af_rows.push_back({"af_low_alt", "chr1", 400, "la_ref_" + std::to_string(i), 0});
        af_rows.push_back({"af_low_alt", "chr1", 400, "la_alt_0", 1});
        // cross_read: observes af_het (allele 0) and af_hom_alt (allele 1);
        // after filter only the af_het observation survives in its profile
        af_rows.push_back({"af_het",     "chr1", 100, "cross_read", 0});
        af_rows.push_back({"af_hom_alt", "chr1", 200, "cross_read", 1});

        Options default_opts;  // min_depth=5, min_alt_depth=2, min_af=0.20, max_af=0.80
        GraphChunkBuildResult af_chunk =
            build_graph_chunk(af_catalog.view_all(), af_rows, "chr1", 0, 500, 0, default_opts);

        ok &= check(af_chunk.chunk.candidates.size() == 1,
                    "af filter: only het site survives");
        ok &= check(!af_chunk.site_ids.empty() &&
                    af_chunk.site_ids[0] == "af_het",
                    "af filter: surviving site id is af_het");
        // het_ref_0..3 + het_alt_0..3 + cross_read = 9 reads with surviving observations
        ok &= check(af_chunk.chunk.reads.size() == 9,
                    "af filter: reads with only filtered-site observations are dropped");
        ok &= check(af_chunk.chunk.candidates[0].counts.total_cov == 9,
                    "af filter: surviving site depth includes cross_read");
    }

    // --- Per-allele min_alt_depth drop ---
    // A triallelic snarl where allele 2 has only 1 supporting read (below min_alt_depth=2).
    // After dropping allele 2, the site becomes biallelic (ref=5, alt1=5) with AF=0.5 → passes.
    // Reads that observed allele 2 lose that observation and are dropped if it was their only site.
    {
        GraphSiteCatalog tri_catalog;
        // Three-allele snarl: walks ">1>2>3" (ref), ">1>4>3" (alt1), ">1>5>3" (alt2)
        GraphSite tri_site;
        tri_site.chrom = "chr1";
        tri_site.ref_contig = "chr1";
        tri_site.pos = 100;
        tri_site.ref_beg = 100;
        tri_site.ref_end = 101;
        tri_site.id = "tri";
        tri_site.allele_traversals = {">1>2>3", ">1>4>3", ">1>5>3"};
        tri_site.allele_walks = {parse_graph_walk(">1>2>3"),
                                 parse_graph_walk(">1>4>3"),
                                 parse_graph_walk(">1>5>3")};
        tri_site.skip_reason = graph_site_validation_skip_reason(tri_site);
        tri_site.eligible = tri_site.skip_reason.empty();
        tri_catalog.sites.push_back(tri_site);

        std::vector<GraphReadAllele> tri_rows;
        // allele 0 (ref): 5 reads
        for (int i = 0; i < 5; ++i)
            tri_rows.push_back({"tri", "chr1", 100, "ref_" + std::to_string(i), 0});
        // allele 1 (alt1): 5 reads — passes min_alt_depth=2
        for (int i = 0; i < 5; ++i)
            tri_rows.push_back({"tri", "chr1", 100, "alt1_" + std::to_string(i), 1});
        // allele 2 (alt2): 1 read — below min_alt_depth=2, must be dropped
        tri_rows.push_back({"tri", "chr1", 100, "noise_read", 2});

        Options default_opts;
        GraphChunkBuildResult tri_chunk =
            build_graph_chunk(tri_catalog.view_all(), tri_rows, "chr1", 0, 200, 0, default_opts);

        // Site survives: after dropping allele 2, alle_covs=[5,5], AF=0.5, total=10 >= min_depth=5
        ok &= check(tri_chunk.chunk.candidates.size() == 1,
                    "per-allele drop: triallelic site survives as biallelic after noise drop");
        ok &= check(!tri_chunk.chunk.candidates.empty() &&
                    tri_chunk.chunk.candidates[0].counts.n_uniq_alles == 2,
                    "per-allele drop: surviving site is biallelic");
        ok &= check(!tri_chunk.chunk.candidates.empty() &&
                    tri_chunk.chunk.candidates[0].counts.total_cov == 10,
                    "per-allele drop: noise read excluded from total_cov");
        // noise_read had only the dropped allele observation → no surviving site → not in reads
        ok &= check(tri_chunk.chunk.reads.size() == 10,
                    "per-allele drop: noise read with only dropped-allele observation is excluded");
    }

    // --- Multiallelic biallelic decomposition ---
    // A triallelic site (ref + alt1 + alt2) is split into two biallelic (ref vs alt_i)
    // pairs. Each pair gets its own AF/depth filter.
    // Sub-test A: both pairs pass — ref reads fan out to both pairs (allele 0 each);
    //             alt reads contribute allele 1 to their own pair and allele 0 to
    //             the other pair when snarl allele phasing is enabled.
    {
        GraphSite tri_both;
        tri_both.chrom = tri_both.ref_contig = "chr1";
        tri_both.pos = tri_both.ref_beg = 100; tri_both.ref_end = 101;
        tri_both.id = "tri_both";
        tri_both.allele_traversals = {">1>2>3", ">1>4>3", ">1>5>3"};
        tri_both.allele_walks = {parse_graph_walk(">1>2>3"),
                                 parse_graph_walk(">1>4>3"),
                                 parse_graph_walk(">1>5>3")};
        tri_both.skip_reason = graph_site_validation_skip_reason(tri_both);
        tri_both.eligible = tri_both.skip_reason.empty();
        GraphSiteCatalog tri_both_cat; tri_both_cat.sites.push_back(tri_both);

        std::vector<GraphReadAllele> tb_rows;
        for (int i = 0; i < 5; ++i)
            tb_rows.push_back({"tri_both", "chr1", 100, "ref_" + std::to_string(i), 0});
        tb_rows.push_back({"tri_both", "chr1", 100, "ref_cross", 0});
        for (int i = 0; i < 5; ++i)
            tb_rows.push_back({"tri_both", "chr1", 100, "a1_" + std::to_string(i), 1});
        for (int i = 0; i < 5; ++i)
            tb_rows.push_back({"tri_both", "chr1", 100, "a2_" + std::to_string(i), 2});

        Options default_opts;
        default_opts.snarl_allele_phasing = true;
        auto tb = build_graph_chunk(tri_both_cat.view_all(), tb_rows, "chr1", 0, 200, 0, default_opts);

        ok &= check(tb.chunk.candidates.size() == 2,
                    "multiallelic decomp: triallelic → 2 biallelic pairs");
        ok &= check(tb.site_ids.size() == 2 &&
                    tb.site_ids[0] == "tri_both:1" && tb.site_ids[1] == "tri_both:2",
                    "multiallelic decomp: pair IDs carry original alt index");
        // Under snarl allele phasing, "ref" for each pair means all reads not
        // carrying that alt: 6 literal-ref reads + 5 reads on the other alt.
        ok &= check(tb.chunk.candidates[0].counts.ref_cov == 11 &&
                    tb.chunk.candidates[0].counts.alt_cov == 5,
                    "multiallelic decomp: pair 0 counts correct");
        ok &= check(tb.chunk.candidates[1].counts.ref_cov == 11 &&
                    tb.chunk.candidates[1].counts.alt_cov == 5,
                    "multiallelic decomp: pair 1 counts correct");
        // All 16 reads (6 ref + 5 alt1 + 5 alt2) survive
        ok &= check(tb.chunk.reads.size() == 16,
                    "multiallelic decomp: all reads retained when both pairs pass");
        const ReadVariantProfile* a1_profile = profile_for_read(tb.chunk, "a1_0");
        ok &= check(a1_profile != nullptr &&
                    a1_profile->start_var_idx == 0 &&
                    a1_profile->end_var_idx == 1 &&
                    a1_profile->alleles.size() == 2 &&
                    a1_profile->alleles[0] == 1 &&
                    a1_profile->alleles[1] == 0,
                    "multiallelic decomp: alt1 read votes ref for alt2 pair");
        const ReadVariantProfile* a2_profile = profile_for_read(tb.chunk, "a2_0");
        ok &= check(a2_profile != nullptr &&
                    a2_profile->start_var_idx == 0 &&
                    a2_profile->end_var_idx == 1 &&
                    a2_profile->alleles.size() == 2 &&
                    a2_profile->alleles[0] == 0 &&
                    a2_profile->alleles[1] == 1,
                    "multiallelic decomp: alt2 read votes ref for alt1 pair");
    }

    // Sub-test A2: keeping a no-reference alt1/alt2 snarl whole must remain a
    // heterozygous n-allelic anchor, not collapse to homozygous non-reference.
    {
        GraphSite tri_alt_alt;
        tri_alt_alt.chrom = tri_alt_alt.ref_contig = "chr1";
        tri_alt_alt.pos = tri_alt_alt.ref_beg = 100; tri_alt_alt.ref_end = 101;
        tri_alt_alt.id = "tri_alt_alt";
        tri_alt_alt.allele_traversals = {">1>2>3", ">1>4>3", ">1>5>3"};
        tri_alt_alt.allele_walks = {parse_graph_walk(">1>2>3"),
                                    parse_graph_walk(">1>4>3"),
                                    parse_graph_walk(">1>5>3")};
        tri_alt_alt.skip_reason = graph_site_validation_skip_reason(tri_alt_alt);
        tri_alt_alt.eligible = tri_alt_alt.skip_reason.empty();
        GraphSiteCatalog tri_alt_alt_cat; tri_alt_alt_cat.sites.push_back(tri_alt_alt);

        std::vector<GraphReadAllele> aa_rows;
        for (int i = 0; i < 5; ++i)
            aa_rows.push_back({"tri_alt_alt", "chr1", 100, "a1_" + std::to_string(i), 1});
        for (int i = 0; i < 5; ++i)
            aa_rows.push_back({"tri_alt_alt", "chr1", 100, "a2_" + std::to_string(i), 2});

        Options keep_whole_opts;
        keep_whole_opts.snarl_keep_whole = true;
        keep_whole_opts.snarl_allele_phasing = true;
        auto aa = build_graph_chunk(tri_alt_alt_cat.view_all(), aa_rows, "chr1", 0, 200, 0,
                                    keep_whole_opts);

        ok &= check(aa.chunk.candidates.size() == 1,
                    "whole snarl: alt1/alt2 site retained as one candidate");
        ok &= check(!aa.chunk.candidates.empty() &&
                    aa.chunk.candidates[0].counts.n_uniq_alles == 3,
                    "whole snarl: candidate keeps all allele slots");
        ok &= check(!aa.chunk.candidates.empty() &&
                    aa.chunk.candidates[0].counts.category == VariantCategory::CleanHetIndel,
                    "whole snarl: alt1/alt2 site classified as heterozygous");
    }

    // Sub-test B: alt2 has AF > max_af (hom-alt), its pair is filtered.
    // Reads observing alt2 are dropped; ref reads appear only at the surviving pair.
    {
        GraphSite tri_pf;
        tri_pf.chrom = tri_pf.ref_contig = "chr1";
        tri_pf.pos = tri_pf.ref_beg = 100; tri_pf.ref_end = 101;
        tri_pf.id = "tri_pf";
        tri_pf.allele_traversals = {">1>2>3", ">1>4>3", ">1>5>3"};
        tri_pf.allele_walks = {parse_graph_walk(">1>2>3"),
                               parse_graph_walk(">1>4>3"),
                               parse_graph_walk(">1>5>3")};
        tri_pf.skip_reason = graph_site_validation_skip_reason(tri_pf);
        tri_pf.eligible = tri_pf.skip_reason.empty();
        GraphSiteCatalog tri_pf_cat; tri_pf_cat.sites.push_back(tri_pf);

        std::vector<GraphReadAllele> pf_rows;
        for (int i = 0; i < 5; ++i)
            pf_rows.push_back({"tri_pf", "chr1", 100, "ref_" + std::to_string(i), 0});
        for (int i = 0; i < 5; ++i)
            pf_rows.push_back({"tri_pf", "chr1", 100, "a1_" + std::to_string(i), 1});
        // alt2: 50 reads → AF = 50/55 > max_af=0.80 → pair filtered
        for (int i = 0; i < 50; ++i)
            pf_rows.push_back({"tri_pf", "chr1", 100, "a2_" + std::to_string(i), 2});

        Options default_opts;
        auto pf = build_graph_chunk(tri_pf_cat.view_all(), pf_rows, "chr1", 0, 200, 0, default_opts);

        ok &= check(pf.chunk.candidates.size() == 1,
                    "per-pair AF filter: only alt1 pair survives (alt2 high_af filtered)");
        ok &= check(!pf.site_ids.empty() && pf.site_ids[0] == "tri_pf:1",
                    "per-pair AF filter: surviving pair id is tri_pf:1");
        ok &= check(!pf.filtered_sites.empty() &&
                    pf.filtered_sites[0].filter_reason == "high_af",
                    "per-pair AF filter: alt2 pair recorded as high_af");
        // Only ref (5) and alt1 (5) reads survive; alt2 reads dropped
        ok &= check(pf.chunk.reads.size() == 10,
                    "per-pair AF filter: alt2 reads dropped, ref+alt1 retained");
    }

    // --- Chunk boundary: site assigned to exactly one chunk ---
    // A site whose ref_beg is at 0-based position 199 (1-based pos=200) spans the
    // boundary between chunk0=[0,200) and chunk1=[200,400).  The start-position
    // check assigns it only to chunk0.  The test pre-filters the catalog per chunk
    // (mirroring the production path) so build_graph_chunk receives only the
    // sites that belong to each chunk.
    {
        GraphSiteCatalog full_catalog;
        full_catalog.sites.push_back(make_site("left",     100, ">1>2>3",  ">1>4>3"));
        full_catalog.sites.push_back(make_site("boundary", 200, ">5>6>7",  ">5>8>7"));
        full_catalog.sites.push_back(make_site("right",    300, ">9>10>11",">9>12>11"));

        // Pre-filter catalog per chunk interval (same as production path).
        auto filter_catalog = [](const GraphSiteCatalog& cat,
                                 hts_pos_t beg, hts_pos_t end) {
            GraphSiteCatalog out;
            for (const GraphSite& site : cat.sites) {
                const hts_pos_t site_beg0 =
                    (site.ref_beg > 0 ? site.ref_beg : site.pos) - 1;
                if (site_beg0 >= beg && site_beg0 < end)
                    out.sites.push_back(site);
            }
            return out;
        };
        GraphSiteCatalog cat0 = filter_catalog(full_catalog, 0, 200);
        GraphSiteCatalog cat1 = filter_catalog(full_catalog, 200, 400);

        // Give every site 6 ref + 6 alt reads so all pass default depth/AF filters.
        std::vector<GraphReadAllele> boundary_rows;
        for (const auto& [sid, pos] : std::vector<std::pair<std::string,int>>{
                 {"left",100}, {"boundary",200}, {"right",300}}) {
            for (int i = 0; i < 6; ++i) {
                boundary_rows.push_back({sid, "chr1", pos,
                                         sid + "_ref_" + std::to_string(i), 0});
                boundary_rows.push_back({sid, "chr1", pos,
                                         sid + "_alt_" + std::to_string(i), 1});
            }
        }

        Options default_opts;
        GraphChunkBuildResult chunk0 =
            build_graph_chunk(cat0.view_all(), boundary_rows, "chr1", 0,   200, 0, default_opts);
        GraphChunkBuildResult chunk1 =
            build_graph_chunk(cat1.view_all(), boundary_rows, "chr1", 200, 400, 1, default_opts);

        // "left" (beg0=99) and "boundary" (beg0=199) both start in [0,200)
        ok &= check(chunk0.chunk.candidates.size() == 2,
                    "boundary: chunk0 contains left and boundary sites");
        // "right" (beg0=299) starts in [200,400); "boundary" must NOT appear again
        ok &= check(chunk1.chunk.candidates.size() == 1,
                    "boundary: chunk1 contains only right site, boundary not duplicated");
        ok &= check(!chunk1.site_ids.empty() && chunk1.site_ids[0] == "right",
                    "boundary: chunk1 site is right");
    }

    // --- Compact stitch replay from the chr20:4.76-4.79 Mb seam ---
    // This synthetic matrix uses candidate coordinates and support patterns from
    // the seam containing the 4,778,793-4,785,720 short-gap target. Recovery
    // first phases BAM-only gap reads into an independent local block.
    // Stitching then propagates the left gauge through that block and into the
    // right block. Parental labels are assertion data and never enter either
    // decision.
    {
        const hts_pos_t left_phase_set = 4642430;
        const hts_pos_t right_phase_set = 4785719;
        const auto make_candidate = [](hts_pos_t pos, VariantType type,
                                       std::array<int, 3> consensus,
                                       hts_pos_t phase_set,
                                       bool homopolymer, bool injected) {
            CandidateVariant candidate;
            candidate.key.pos = pos;
            candidate.key.type = type;
            candidate.key.ref_len = type == VariantType::Deletion ? 20 : 0;
            candidate.counts.n_uniq_alles = 2;
            candidate.counts.alle_covs = {32, 32};
            candidate.counts.total_cov = 64;
            candidate.counts.ref_cov = 32;
            candidate.counts.alt_cov = 32;
            candidate.counts.allele_fraction = 0.5;
            candidate.counts.category = type == VariantType::Snp
                                            ? VariantCategory::CleanHetSnp
                                            : VariantCategory::CleanHetIndel;
            candidate.lcd_var_i_to_cate = type == VariantType::Snp
                                              ? kCandCleanHetSnp
                                              : kCandCleanHetIndel;
            candidate.hap_to_cons_alle = consensus;
            candidate.phase_set = phase_set;
            candidate.is_homopolymer_indel = homopolymer;
            candidate.bam_injected = injected;
            candidate.alignment_verified = injected;
            return candidate;
        };

        // Phase the gap reads with the production BAM solver. Truth labels are
        // retained only for the assertion below and never enter the solve.
        PhasingChunk local;
        local.candidates = {
            make_candidate(4766929, VariantType::Insertion, {-1, 0, 1},
                           kUnsetCandidatePhaseSet, false, false),
            make_candidate(4778793, VariantType::Insertion, {-1, 0, 1},
                           kUnsetCandidatePhaseSet, false, false),
        };
        for (CandidateVariant& candidate : local.candidates) {
            candidate.counts.alle_covs = {4, 5};
            candidate.counts.total_cov = 9;
            candidate.counts.ref_cov = 4;
            candidate.counts.alt_cov = 5;
        }
        std::vector<char> local_truth;
        const auto add_local_read = [&](const std::string& name, char truth,
                                        int allele) {
            const int read_i = static_cast<int>(local.reads.size());
            ReadRecord read;
            read.qname = name;
            read.mapq = 60;
            local.reads.push_back(std::move(read));
            ReadVariantProfile profile;
            profile.read_id = read_i;
            profile.start_var_idx = 0;
            profile.end_var_idx = 1;
            profile.alleles = {allele, allele};
            profile.alt_qi = {60, 60};
            local.read_var_profile.push_back(std::move(profile));
            local_truth.push_back(truth);
        };
        for (int i = 0; i < 4; ++i)
            add_local_read("gap_m_" + std::to_string(i), 'M', 1);
        add_local_read("gap_to_right_m", 'M', 1);
        for (int i = 0; i < 4; ++i)
            add_local_read("gap_p_" + std::to_string(i), 'P', 0);
        rebuild_read_var_cr(local);
        Options local_opts;
        local_opts.msa_sites_vote_without_gap_link = true;
        local_opts.infer_complement_at_multiallelic = true;
        local_opts.upstream_read_scoring = true;
        local_opts.upstream_assign_hap = true;
        assign_hap_based_on_germline_het_vars_kmeans(
            local, local_opts, kCandGermlineClean);

        int local_straight = 0;
        int local_flipped = 0;
        for (size_t read_i = 0; read_i < local.reads.size(); ++read_i) {
            const bool maternal = local_truth[read_i] == 'M';
            local_straight += local.haps[read_i] == (maternal ? 2 : 1);
            local_flipped += local.haps[read_i] == (maternal ? 1 : 2);
        }
        ok &= check(std::max(local_straight, local_flipped) == 9 &&
                    local.candidates[0].phase_set > 0 &&
                    local.candidates[0].phase_set == local.candidates[1].phase_set,
                    "4.78 Mb fixture: production solver phases all local reads");
        const hts_pos_t gap_phase_set = local.candidates[0].phase_set;
        for (CandidateVariant& candidate : local.candidates) {
            candidate.bam_injected = true;
            candidate.alignment_verified = true;
        }

        PhasingChunk boundary;
        boundary.candidates = {
            make_candidate(4761808, VariantType::Snp, {-1, 0, 1},
                           left_phase_set, false, false),
            local.candidates[0],
            local.candidates[1],
            // The right block has the opposite numeric HP gauge.
            make_candidate(4785720, VariantType::Deletion, {-1, 1, 0},
                           right_phase_set, false, false),
            // No read connects this far site to the seam. It must still flip
            // because PS orientation is an atomic block property.
            make_candidate(4790000, VariantType::Snp, {-1, 0, 1},
                           right_phase_set, false, false),
        };

        std::vector<char> parental_truth;
        const auto add_read = [&](const std::string& name, char truth,
                                  int incoming_hap, hts_pos_t incoming_phase_set,
                                  int first, std::vector<int> alleles) {
            const int read_i = static_cast<int>(boundary.reads.size());
            ReadRecord read;
            read.qname = name;
            read.mapq = 60;
            boundary.reads.push_back(std::move(read));
            ReadVariantProfile profile;
            profile.read_id = read_i;
            profile.start_var_idx = first;
            profile.end_var_idx = first + static_cast<int>(alleles.size()) - 1;
            profile.alleles = std::move(alleles);
            profile.alt_qi.assign(profile.alleles.size(), 60);
            boundary.read_var_profile.push_back(std::move(profile));
            boundary.haps.push_back(incoming_hap);
            boundary.phase_sets.push_back(incoming_phase_set);
            parental_truth.push_back(truth);
        };
        for (int i = 0; i < 3; ++i)
            add_read("left_m_" + std::to_string(i), 'M', 2, left_phase_set,
                     0, {1, 1});
        for (int i = 0; i < 3; ++i)
            add_read("left_p_" + std::to_string(i), 'P', 1, left_phase_set,
                     0, {0, 0});
        for (size_t i = 0; i < local.reads.size(); ++i) {
            if (local.reads[i].qname == "gap_to_right_m") continue;
            const int allele = local_truth[i] == 'M' ? 1 : 0;
            add_read(local.reads[i].qname, local_truth[i], local.haps[i],
                     gap_phase_set, 1, {allele, allele});
        }
        // One local-block molecule reaches the right-block boundary and
        // establishes a crossed edge between their independent gauges.
        add_read("gap_to_right_m", 'M', local.haps[4], gap_phase_set,
                 2, {1, 1});
        for (int i = 0; i < 13; ++i)
            add_read("right_m_" + std::to_string(i), 'M', 1, right_phase_set,
                     3, {1});
        for (int i = 0; i < 31; ++i)
            add_read("right_p_" + std::to_string(i), 'P', 2, right_phase_set,
                     3, {0});
        for (int i = 0; i < 4; ++i)
            add_read("right_tail_m_" + std::to_string(i), 'M', 1,
                     right_phase_set, 4, {0});
        for (int i = 0; i < 4; ++i)
            add_read("right_tail_p_" + std::to_string(i), 'P', 2,
                     right_phase_set, 4, {1});
        rebuild_read_var_cr(boundary);

        // Validate local gap phasing independently of its arbitrary HP gauge.
        int gap_straight = 0;
        int gap_flipped = 0;
        for (size_t read_i = 0; read_i < boundary.reads.size(); ++read_i) {
            if (boundary.phase_sets[read_i] != gap_phase_set) continue;
            const bool maternal = parental_truth[read_i] == 'M';
            gap_straight += boundary.haps[read_i] == (maternal ? 2 : 1);
            gap_flipped += boundary.haps[read_i] == (maternal ? 1 : 2);
        }
        ok &= check(std::max(gap_straight, gap_flipped) == 9,
                    "4.78 Mb fixture: local gap PS phases all gap reads correctly");

        Options stitch_opts;
        stitch_opts.block_link_window = 128;
        stitch_opts.min_block_link_reads = 1;
        // A balanced same/cross vote must leave the downstream block intact.
        for (size_t read_i = 0; read_i < 3; ++read_i)
            boundary.read_var_profile[read_i].alleles[1] = 0;
        rebuild_read_var_cr(boundary);
        ok &= check(!stitch_phase_sets_by_alleles(
                        boundary, left_phase_set, gap_phase_set, stitch_opts) &&
                    boundary.candidates[1].phase_set == gap_phase_set &&
                    boundary.candidates[1].hap_to_cons_alle ==
                        local.candidates[0].hap_to_cons_alle,
                    "4.78 Mb fixture: tied edge leaves the local block unchanged");
        for (size_t read_i = 0; read_i < 3; ++read_i)
            boundary.read_var_profile[read_i].alleles[1] = 1;
        rebuild_read_var_cr(boundary);

        const std::vector<RecoverySeam> recovery_windows = {
            {4761808, 4785720, left_phase_set, right_phase_set},
        };
        RecoveryPhaseGauge recovery_gauge;
        recovery_gauge.beg = 4700000;
        recovery_gauge.end = 4850000;
        recovery_gauge.imported_phase_sets = {gap_phase_set};
        const bool bam_matches_left_gauge = local_straight == 9;
        recovery_gauge.graph_votes = {
            PhaseSetGaugeVote{left_phase_set,
                              bam_matches_left_gauge ? 45 : 0,
                              bam_matches_left_gauge ? 0 : 45},
            PhaseSetGaugeVote{right_phase_set,
                              bam_matches_left_gauge ? 0 : 63,
                              bam_matches_left_gauge ? 63 : 0},
        };
        if (bam_matches_left_gauge) {
            recovery_gauge.block_votes = {
                RecoveryBlockGaugeVote{
                    left_phase_set, gap_phase_set, {{{8, 0}, {0, 8}}}},
                RecoveryBlockGaugeVote{
                    right_phase_set, gap_phase_set, {{{0, 8}, {8, 0}}}},
            };
        } else {
            recovery_gauge.block_votes = {
                RecoveryBlockGaugeVote{
                    left_phase_set, gap_phase_set, {{{0, 8}, {8, 0}}}},
                RecoveryBlockGaugeVote{
                    right_phase_set, gap_phase_set, {{{8, 0}, {0, 8}}}},
            };
        }
        ok &= check(stitch_recovery_phase_sets_left_to_right(
                        boundary, recovery_windows, {recovery_gauge},
                        stitch_opts) == 0,
                    "4.78 Mb fixture: unvalidated one-block seam remains independent");

        ok &= check(boundary.candidates[1].phase_set == gap_phase_set &&
                    boundary.candidates[2].phase_set == gap_phase_set &&
                    boundary.candidates[3].phase_set == right_phase_set &&
                    boundary.candidates[4].phase_set == right_phase_set,
                    "4.78 Mb fixture: failed whole-seam validation is atomic");

        int assigned = 0;
        for (size_t read_i = 0; read_i < boundary.reads.size(); ++read_i) {
            assigned += (boundary.haps[read_i] == 1 ||
                         boundary.haps[read_i] == 2) &&
                        boundary.phase_sets[read_i] > 0;
        }
        ok &= check(assigned == static_cast<int>(boundary.reads.size()),
                    "4.78 Mb fixture: independent blocks preserve read assignments");
    }

    // The ordinary adjacent-block stitch cannot cross two disjoint read
    // cohorts. The DP fallback must find the supported injected site between
    // them, solve its two orientations, and join the right block atomically.
    {
        constexpr hts_pos_t kLeft = 100;
        constexpr hts_pos_t kRight = 300;
        PhasingChunk replay;
        const auto add_candidate = [&](hts_pos_t pos, hts_pos_t phase_set,
                                       std::array<int, 3> consensus,
                                       bool injected) {
            CandidateVariant candidate;
            candidate.key.pos = pos;
            candidate.key.type = VariantType::Snp;
            candidate.counts.n_uniq_alles = 2;
            candidate.counts.alle_covs = {12, 12};
            candidate.lcd_var_i_to_cate = kCandCleanHetSnp;
            candidate.hap_to_cons_alle = consensus;
            candidate.phase_set = phase_set;
            candidate.bam_injected = injected;
            replay.candidates.push_back(std::move(candidate));
        };
        add_candidate(1000, kLeft, {-1, 0, 1}, false);
        add_candidate(1100, kLeft, {-1, 0, 1}, false);
        add_candidate(10000, 0,
                      {-1, -1, -1}, true);
        add_candidate(20000, kRight, {-1, 0, 1}, false);
        add_candidate(20100, kRight, {-1, 0, 1}, false);

        const auto add_read = [&](std::string name, int first,
                                  std::vector<int> alleles, int hap,
                                  hts_pos_t phase_set) {
            ReadRecord read;
            read.qname = std::move(name);
            read.mapq = 60;
            replay.reads.push_back(std::move(read));
            ReadVariantProfile profile;
            profile.read_id = static_cast<int>(replay.read_var_profile.size());
            profile.start_var_idx = first;
            profile.end_var_idx = first + static_cast<int>(alleles.size()) - 1;
            profile.alleles = std::move(alleles);
            profile.alt_qi.assign(profile.alleles.size(), 60);
            replay.read_var_profile.push_back(std::move(profile));
            replay.haps.push_back(hap);
            replay.phase_sets.push_back(phase_set);
        };
        for (int allele = 0; allele <= 1; ++allele) {
            for (int i = 0; i < 6; ++i) {
                add_read("left_path_" + std::to_string(allele) + "_" +
                             std::to_string(i),
                         1, {allele, allele}, allele + 1, kLeft);
                add_read("right_path_" + std::to_string(allele) + "_" +
                             std::to_string(i),
                         2, {allele, allele}, allele + 1, kRight);
            }
        }
        add_read("unassigned_internal", 2, {0}, 0,
                 kUnphasedReadPhaseSet);
        rebuild_read_var_cr(replay);

        Options stitch_opts;
        stitch_opts.min_block_link_reads = 1;
        stitch_opts.block_link_window = 8;
        const size_t joined = stitch_recovery_phase_sets_left_to_right(
            replay, {{1100, 20000, kLeft, kRight}}, {}, stitch_opts);
        ok &= check(joined == 1 &&
                    std::all_of(replay.candidates.begin(),
                                replay.candidates.end(),
                                [](const CandidateVariant& candidate) {
                                    return candidate.phase_set == kLeft;
                                }) &&
                    replay.candidates[2].hap_to_cons_alle[1] == 0 &&
                    replay.candidates[2].hap_to_cons_alle[2] == 1,
                    "recovery DP: exact supported path joins disjoint read cohorts");
        ok &= check(replay.haps.back() == 1 &&
                    replay.phase_sets.back() == kLeft,
                    "recovery DP: exact path allele assigns an unphased read");
    }

    // Every local edge has six agreeing reads, below the exact-path DP's
    // ten-read edge threshold. Taken together they define one unique diploid
    // minimum-error chain. The final MEC pass must solve the whole seam without
    // weakening the established direct-edge thresholds.
    {
        constexpr hts_pos_t kLeft = 100;
        constexpr hts_pos_t kRight = 300;
        PhasingChunk replay;
        const auto add_candidate = [&](hts_pos_t pos, hts_pos_t phase_set,
                                       bool injected) {
            CandidateVariant candidate;
            candidate.key.pos = pos;
            candidate.key.type = VariantType::Snp;
            candidate.counts.n_uniq_alles = 2;
            candidate.counts.alle_covs = {6, 6};
            candidate.lcd_var_i_to_cate = kCandCleanHetSnp;
            candidate.hap_to_cons_alle = injected
                ? std::array<int, 3>{-1, -1, -1}
                : std::array<int, 3>{-1, 0, 1};
            candidate.phase_set = phase_set;
            candidate.bam_injected = injected;
            candidate.alignment_verified = injected;
            replay.candidates.push_back(std::move(candidate));
        };
        add_candidate(1000, kLeft, false);
        add_candidate(1100, kLeft, false);
        add_candidate(8000, kUnsetCandidatePhaseSet, true);
        add_candidate(12000, kUnsetCandidatePhaseSet, true);
        add_candidate(20000, kRight, false);
        add_candidate(20100, kRight, false);

        const auto add_edge_reads = [&](const std::string& prefix, int first) {
            for (int allele = 0; allele <= 1; ++allele) {
                for (int copy = 0; copy < 3; ++copy) {
                    ReadRecord read;
                    read.qname = prefix + "_" + std::to_string(allele) +
                                 "_" + std::to_string(copy);
                    read.mapq = 60;
                    replay.reads.push_back(std::move(read));
                    ReadVariantProfile profile;
                    profile.read_id = static_cast<int>(
                        replay.read_var_profile.size());
                    profile.start_var_idx = first;
                    profile.end_var_idx = first + 1;
                    profile.alleles = {allele, allele};
                    profile.alt_qi = {60, 60};
                    replay.read_var_profile.push_back(std::move(profile));
                    replay.haps.push_back(0);
                    replay.phase_sets.push_back(kUnphasedReadPhaseSet);
                }
            }
        };
        add_edge_reads("mec_left", 1);
        add_edge_reads("mec_middle", 2);
        add_edge_reads("mec_right", 3);
        rebuild_read_var_cr(replay);

        Options stitch_opts;
        stitch_opts.min_block_link_reads = 8;
        stitch_opts.block_link_window = 8;
        const size_t joined = stitch_recovery_phase_sets_left_to_right(
            replay, {{1100, 20000, kLeft, kRight}}, {}, stitch_opts);
        ok &= check(joined == 1 &&
                    std::all_of(replay.candidates.begin(),
                                replay.candidates.end(),
                                [](const CandidateVariant& candidate) {
                                    return candidate.phase_set == kLeft;
                                }) &&
                    replay.candidates[2].hap_to_cons_alle[1] == 0 &&
                    replay.candidates[3].hap_to_cons_alle[1] == 0,
                    "recovery MEC: globally supported weak-edge chain joins exactly");
    }

    // No single locus pair reaches the ordinary stitch threshold here:
    // four read cohorts cover four different boundary-site pairs. Aggregating
    // each molecule's consensus across the two phase sets recovers the one
    // well-supported parity without counting a read more than once.
    {
        constexpr hts_pos_t kLeft = 100;
        constexpr hts_pos_t kRight = 300;
        PhasingChunk replay;
        const auto add_candidate = [&](hts_pos_t pos, hts_pos_t phase_set) {
            CandidateVariant candidate;
            candidate.key.pos = pos;
            candidate.key.type = VariantType::Snp;
            candidate.counts.n_uniq_alles = 2;
            candidate.counts.alle_covs = {8, 8};
            candidate.lcd_var_i_to_cate = kCandCleanHetSnp;
            candidate.hap_to_cons_alle = {-1, 0, 1};
            candidate.phase_set = phase_set;
            replay.candidates.push_back(std::move(candidate));
        };
        add_candidate(1000, kLeft);
        add_candidate(1100, kLeft);
        add_candidate(2000, kRight);
        add_candidate(2100, kRight);

        for (int left_i = 0; left_i < 2; ++left_i) {
            for (int right_i = 2; right_i < 4; ++right_i) {
                for (int allele = 0; allele < 2; ++allele) {
                    for (int copy = 0; copy < 2; ++copy) {
                        ReadRecord read;
                        read.qname = "aggregate_" + std::to_string(left_i) +
                                     "_" + std::to_string(right_i) + "_" +
                                     std::to_string(allele) + "_" +
                                     std::to_string(copy);
                        read.mapq = 60;
                        replay.reads.push_back(std::move(read));

                        ReadVariantProfile profile;
                        profile.read_id = static_cast<int>(
                            replay.read_var_profile.size());
                        profile.start_var_idx = 0;
                        profile.end_var_idx = 3;
                        profile.alleles.assign(4, -1);
                        profile.alleles[static_cast<size_t>(left_i)] = allele;
                        profile.alleles[static_cast<size_t>(right_i)] = allele;
                        profile.alt_qi.assign(4, 60);
                        replay.read_var_profile.push_back(std::move(profile));
                        replay.haps.push_back(allele + 1);
                        replay.phase_sets.push_back(kLeft);
                    }
                }
            }
        }
        rebuild_read_var_cr(replay);

        Options stitch_opts;
        stitch_opts.min_block_link_reads = 8;
        stitch_opts.block_link_window = 8;
        ok &= check(!stitch_phase_sets_by_alleles(
                        replay, kLeft, kRight, stitch_opts),
                    "recovery aggregate: every individual site pair abstains");
        const size_t joined = stitch_recovery_phase_sets_left_to_right(
            replay, {{1100, 2000, kLeft, kRight}}, {}, stitch_opts);
        ok &= check(joined == 1 &&
                    replay.candidates[2].phase_set == kLeft &&
                    replay.candidates[3].phase_set == kLeft,
                    "recovery aggregate: distributed read evidence joins the blocks");
    }

    // A targeted solve may put several MSA candidates in one local phase set.
    // Stitching preserves the BAM subsolve's read assignment and does not
    // assign an unrelated, previously unphased read from the completed chain.
    {
        constexpr hts_pos_t kLeft = 100;
        constexpr hts_pos_t kLocal = 200;
        constexpr hts_pos_t kRight = 300;
        PhasingChunk replay;
        const auto add_candidate = [&](hts_pos_t pos, hts_pos_t phase_set,
                                       std::array<int, 3> consensus,
                                       bool injected) {
            CandidateVariant candidate;
            candidate.key.pos = pos;
            candidate.key.type = VariantType::Deletion;
            candidate.key.ref_len = 1;
            candidate.counts.n_uniq_alles = 2;
            candidate.counts.alle_covs = {5, 5};
            candidate.lcd_var_i_to_cate = kCandNoisyCandHet;
            candidate.hap_to_cons_alle = consensus;
            candidate.phase_set = phase_set;
            candidate.bam_injected = injected;
            replay.candidates.push_back(std::move(candidate));
        };
        add_candidate(1000, kLeft, {-1, 0, 1}, false);
        add_candidate(2000, kLocal, {-1, 0, 1}, true);   // unrelated MSA row
        add_candidate(3000, kLocal, {-1, 0, 1}, true);   // selected bridge
        add_candidate(4000, kRight, {-1, 0, 1}, false);

        const auto add_read = [&](std::string name, int first,
                                  std::vector<int> alleles, int hap,
                                  hts_pos_t phase_set) {
            ReadRecord read;
            read.qname = std::move(name);
            read.mapq = 60;
            replay.reads.push_back(std::move(read));
            ReadVariantProfile profile;
            profile.read_id = static_cast<int>(replay.read_var_profile.size());
            profile.start_var_idx = first;
            profile.end_var_idx = first + static_cast<int>(alleles.size()) - 1;
            profile.alleles = std::move(alleles);
            profile.alt_qi.assign(profile.alleles.size(), 60);
            replay.read_var_profile.push_back(std::move(profile));
            replay.haps.push_back(hap);
            replay.phase_sets.push_back(phase_set);
        };
        for (int i = 0; i < 4; ++i)
            add_read("edge_" + std::to_string(i), 2, {0, 1}, 0,
                     kUnphasedReadPhaseSet);
        // The unsupported row says hap1 while the selected bridge says hap2.
        add_read("read_to_rescore", 1, {0, 1}, 0,
                 kUnphasedReadPhaseSet);
        rebuild_read_var_cr(replay);

        RecoveryPhaseGauge gauge;
        gauge.beg = 900;
        gauge.end = 4100;
        gauge.imported_phase_sets = {kLocal};
        gauge.graph_votes = {
            PhaseSetGaugeVote{kLeft, 4, 0},
            PhaseSetGaugeVote{kRight, 0, 4},
        };
        gauge.block_votes = {
            RecoveryBlockGaugeVote{kLeft, kLocal, {{{8, 0}, {0, 8}}}},
            RecoveryBlockGaugeVote{kRight, kLocal, {{{0, 8}, {8, 0}}}},
        };
        Options stitch_opts;
        stitch_opts.min_block_link_reads = 1;
        stitch_opts.block_link_window = 8;
        const size_t joined = stitch_recovery_phase_sets_left_to_right(
            replay, {{1000, 4000, kLeft, kRight}}, {gauge}, stitch_opts);
        ok &= check(joined == 0 &&
                    !replay.candidates[1].gap_link_supported &&
                    !replay.candidates[2].gap_link_supported &&
                    replay.haps.back() == 0 &&
                    replay.phase_sets.back() == kUnphasedReadPhaseSet,
                    "recovery chain: stitch preserves BAM read assignments");
    }

    // An imported block without direct block-to-flank evidence remains
    // independent even when aggregate allele evidence would join it. The same
    // post-break graph shortcut also abstains.
    {
        const auto run_post_break = [](bool injected) {
            constexpr hts_pos_t kLeft = 100;
            constexpr hts_pos_t kLocal = 200;
            constexpr hts_pos_t kRight = 300;
            PhasingChunk replay;
            for (const auto [pos, phase_set, is_injected] : {
                     std::tuple<hts_pos_t, hts_pos_t, bool>{1000, kLeft, false},
                     {2000, kLocal, injected},
                     {3000, kRight, false}}) {
                CandidateVariant candidate;
                candidate.key.pos = pos;
                candidate.key.type = VariantType::Snp;
                candidate.counts.n_uniq_alles = 2;
                candidate.counts.alle_covs = {8, 8};
                candidate.lcd_var_i_to_cate = kCandCleanHetSnp;
                candidate.hap_to_cons_alle = {-1, 0, 1};
                candidate.phase_set = phase_set;
                candidate.bam_injected = is_injected;
                candidate.alignment_verified = is_injected;
                replay.candidates.push_back(std::move(candidate));
            }
            for (int allele = 0; allele <= 1; ++allele) {
                for (int copy = 0; copy < 4; ++copy) {
                    ReadRecord read;
                    read.qname = "post_break_" + std::to_string(allele) +
                                 "_" + std::to_string(copy);
                    read.mapq = 60;
                    replay.reads.push_back(std::move(read));
                    ReadVariantProfile profile;
                    profile.read_id = static_cast<int>(
                        replay.read_var_profile.size());
                    profile.start_var_idx = 1;
                    profile.end_var_idx = 2;
                    profile.alleles = {allele, allele};
                    profile.alt_qi = {60, 60};
                    replay.read_var_profile.push_back(std::move(profile));
                    replay.haps.push_back(0);
                    replay.phase_sets.push_back(kUnphasedReadPhaseSet);
                }
            }
            rebuild_read_var_cr(replay);
            Options stitch_opts;
            stitch_opts.min_block_link_reads = 1;
            stitch_opts.block_link_window = 8;
            const size_t joined = stitch_recovery_phase_sets_left_to_right(
                replay, {{1000, 3000, kLeft, kRight}}, {}, stitch_opts);
            return std::make_pair(joined, std::move(replay));
        };

        auto injected = run_post_break(true);
        ok &= check(injected.first == 0 &&
                    injected.second.candidates[0].phase_set == 100 &&
                    injected.second.candidates[1].phase_set == 200 &&
                    injected.second.candidates[2].phase_set == 300,
                    "recovery chain: BAM block without direct gauge stays independent");
        auto graph_only = run_post_break(false);
        ok &= check(graph_only.first == 0 &&
                    graph_only.second.candidates[0].phase_set == 100 &&
                    graph_only.second.candidates[1].phase_set == 200 &&
                    graph_only.second.candidates[2].phase_set == 300,
                    "recovery chain: graph shortcut abstains after break");
    }

    // Distinct centered SNPs can carry a graph/BAM relation even when
    // their VariantKeys differ. Require significant support independently in
    // both deterministic read halves; the exact-key candidate count stays zero.
    {
        constexpr hts_pos_t kLeft = 101;
        constexpr hts_pos_t kLocal = 202;
        constexpr hts_pos_t kRight = 303;
        PhasingChunk replay;
        for (const auto [pos, phase_set, injected] : {
                 std::tuple<hts_pos_t, hts_pos_t, bool>{1000, kLeft, false},
                 {2000, kLocal, true},
                 {3000, kRight, false}}) {
            CandidateVariant candidate;
            candidate.key.pos = pos;
            candidate.key.type = VariantType::Snp;
            candidate.counts.ref_cov = 16;
            candidate.counts.alt_cov = 16;
            candidate.counts.allele_fraction = 0.5;
            candidate.hap_to_cons_alle = {-1, 0, 1};
            candidate.phase_set = phase_set;
            candidate.bam_injected = injected;
            candidate.alignment_verified = injected;
            replay.candidates.push_back(std::move(candidate));
        }
        const auto stable_fold = [](const std::string& name) {
            uint64_t hash = 14695981039346656037ULL;
            for (const unsigned char byte : name) {
                hash ^= byte;
                hash *= 1099511628211ULL;
            }
            return static_cast<int>(hash & 1ULL);
        };
        std::array<int, 2> fold_counts{0, 0};
        for (int copy = 0; fold_counts[0] < 8 || fold_counts[1] < 8; ++copy) {
            const std::string name = "distinct_snp_" + std::to_string(copy);
            const int fold = stable_fold(name);
            if (fold_counts[fold] >= 8) continue;
            const int allele = fold_counts[fold]++ & 1;
            ReadRecord read;
            read.qname = name;
            read.mapq = 60;
            replay.reads.push_back(std::move(read));
            ReadVariantProfile profile;
            profile.read_id = static_cast<int>(
                replay.read_var_profile.size());
            profile.start_var_idx = 0;
            profile.end_var_idx = 1;
            profile.alleles = {allele, allele};
            profile.alt_qi = {60, 60};
            replay.read_var_profile.push_back(std::move(profile));
            replay.haps.push_back(0);
            replay.phase_sets.push_back(kUnphasedReadPhaseSet);
        }
        rebuild_read_var_cr(replay);
        RecoveryPhaseGauge gauge;
        gauge.beg = 900;
        gauge.end = 3100;
        gauge.imported_phase_sets = {kLocal};
        RecoveryBlockGaugeVote distinct_snp_vote{
            kLeft, kLocal, {{{8, 0}, {0, 8}}}};
        distinct_snp_vote.shared_candidate_same = 2;
        gauge.block_votes = {distinct_snp_vote};
        Options stitch_opts;
        stitch_opts.min_block_link_reads = 1;
        stitch_opts.block_link_window = 8;
        const size_t joined = stitch_recovery_phase_sets_left_to_right(
            replay, {{1000, 3000, kLeft, kRight}}, {gauge}, stitch_opts);
        ok &= check(joined == 1 &&
                    replay.candidates[0].phase_set == kLeft &&
                    replay.candidates[1].phase_set == kLeft &&
                    replay.candidates[2].phase_set == kRight,
                    "recovery representation: split-stable distinct SNPs replace exact key");
    }

    // Keep complementary multi-allelic BAM rows separate. A selected verified
    // boundary indel may be off-center because row-wise AF measures that one
    // alternate against all other alleles; the source gauge remains mandatory.
    {
        constexpr hts_pos_t kLeft = 101;
        constexpr hts_pos_t kLocal = 202;
        constexpr hts_pos_t kRight = 303;
        PhasingChunk replay;
        CandidateVariant left;
        left.key.pos = 1000;
        left.key.type = VariantType::Snp;
        left.counts.ref_cov = 16;
        left.counts.alt_cov = 16;
        left.counts.allele_fraction = 0.5;
        left.hap_to_cons_alle = {-1, 0, 1};
        left.phase_set = kLeft;
        CandidateVariant local;
        local.key.pos = 2000;
        local.key.type = VariantType::Deletion;
        local.key.ref_len = 5;
        local.counts.ref_cov = 8;
        local.counts.alt_cov = 24;
        local.counts.allele_fraction = 0.75;
        local.hap_to_cons_alle = {-1, 0, 1};
        local.phase_set = kLocal;
        local.bam_injected = true;
        local.alignment_verified = true;
        CandidateVariant right = left;
        right.key.pos = 3000;
        right.phase_set = kRight;
        replay.candidates = {left, local, right};
        for (int allele = 0; allele <= 1; ++allele) {
            for (int copy = 0; copy < 12; ++copy) {
                ReadRecord read;
                read.qname = "off_center_indel_" +
                             std::to_string(allele) + "_" +
                             std::to_string(copy);
                read.mapq = 60;
                replay.reads.push_back(std::move(read));
                ReadVariantProfile profile;
                profile.read_id = static_cast<int>(
                    replay.read_var_profile.size());
                profile.start_var_idx = 1;
                profile.end_var_idx = 2;
                profile.alleles = {allele, allele};
                profile.alt_qi = {60, 60};
                replay.read_var_profile.push_back(std::move(profile));
                replay.haps.push_back(0);
                replay.phase_sets.push_back(kUnphasedReadPhaseSet);
            }
        }
        rebuild_read_var_cr(replay);
        RecoveryPhaseGauge gauge;
        gauge.beg = 900;
        gauge.end = 3100;
        gauge.imported_phase_sets = {kLocal};
        RecoveryBlockGaugeVote vote{
            kRight, kLocal, {{{12, 0}, {0, 12}}}};
        vote.shared_candidate_same = 1;
        gauge.block_votes = {vote};
        Options stitch_opts;
        stitch_opts.min_block_link_reads = 1;
        stitch_opts.block_link_window = 8;
        const size_t joined = stitch_recovery_phase_sets_left_to_right(
            replay, {{1000, 3000, kLeft, kRight}}, {gauge}, stitch_opts);
        ok &= check(joined == 1 &&
                    replay.candidates[0].phase_set == kLeft &&
                    replay.candidates[1].phase_set == kLocal &&
                    replay.candidates[2].phase_set == kLocal,
                    "recovery representation: verified off-center boundary indel remains separate");
    }

    // The exact solve validates the complete neighboring blocks. Boundary
    // scope is available only when that problem exceeds the variable budget;
    // a statistically supported candidate gauge must authorize it, while one
    // shared anchor must still abstain.
    {
        constexpr hts_pos_t kLeft = 101;
        constexpr hts_pos_t kLocalA = 202;
        constexpr hts_pos_t kLocal = 203;
        constexpr hts_pos_t kRight = 303;
        PhasingChunk replay;
        const auto add_candidate = [&](hts_pos_t pos, hts_pos_t phase_set,
                                       bool injected, bool is_oriented) {
            CandidateVariant candidate;
            candidate.key.pos = pos;
            candidate.key.type = VariantType::Snp;
            candidate.counts.ref_cov = 16;
            candidate.counts.alt_cov = 16;
            candidate.counts.allele_fraction = 0.5;
            candidate.hap_to_cons_alle =
                is_oriented ? std::array<int, 3>{-1, 0, 1}
                            : std::array<int, 3>{-1, -1, -1};
            candidate.phase_set = phase_set;
            candidate.bam_injected = injected;
            candidate.alignment_verified = injected;
            replay.candidates.push_back(std::move(candidate));
        };
        add_candidate(1000, kLeft, false, true);
        add_candidate(1500, kLocalA, true, true);
        add_candidate(2000, kLocal, true, true);
        add_candidate(3000, kRight, false, true);
        constexpr int kOutsideVariables = 21;
        for (int i = 0; i < kOutsideVariables; ++i)
            add_candidate(3100 + i, kUnsetCandidatePhaseSet, false, false);
        add_candidate(6000, kRight, false, true);

        const auto add_read = [&](std::string name, int first,
                                  std::vector<int> alleles) {
            ReadRecord read;
            read.qname = std::move(name);
            read.mapq = 60;
            replay.reads.push_back(std::move(read));
            ReadVariantProfile profile;
            profile.read_id = static_cast<int>(replay.read_var_profile.size());
            profile.start_var_idx = first;
            profile.end_var_idx = first + static_cast<int>(alleles.size()) - 1;
            profile.alleles = std::move(alleles);
            profile.alt_qi.assign(profile.alleles.size(), 60);
            replay.read_var_profile.push_back(std::move(profile));
            replay.haps.push_back(0);
            replay.phase_sets.push_back(kUnphasedReadPhaseSet);
        };
        for (int allele = 0; allele <= 1; ++allele)
            for (int copy = 0; copy < 12; ++copy)
                add_read("local_edge_" + std::to_string(allele) + "_" +
                             std::to_string(copy),
                         2, {allele, allele});
        for (int variable = 0; variable < kOutsideVariables; ++variable) {
            std::vector<int> ref(static_cast<size_t>(variable + 2), -1);
            std::vector<int> alt(static_cast<size_t>(variable + 2), -1);
            ref[0] = 0;
            ref[static_cast<size_t>(variable + 1)] = 0;
            alt[0] = 1;
            alt[static_cast<size_t>(variable + 1)] = 1;
            add_read("outside_ref_" + std::to_string(variable), 3,
                     std::move(ref));
            add_read("outside_alt_" + std::to_string(variable), 3,
                     std::move(alt));
        }
        std::vector<int> far_ref(
            static_cast<size_t>(kOutsideVariables + 2), -1);
        std::vector<int> far_alt(
            static_cast<size_t>(kOutsideVariables + 2), -1);
        far_ref.front() = 0;
        far_ref.back() = 0;
        far_alt.front() = 1;
        far_alt.back() = 1;
        // The distant end of the graph block has the opposite local polarity.
        // Equal molecule support makes whole-block orientation ambiguous, so a
        // lone shared candidate at the boundary must not hide this older break.
        replay.candidates.back().hap_to_cons_alle = {-1, 1, 0};
        for (int copy = 0; copy < 12; ++copy) {
            add_read("outside_far_ref_" + std::to_string(copy), 3, far_ref);
            add_read("outside_far_alt_" + std::to_string(copy), 3, far_alt);
        }
        rebuild_read_var_cr(replay);

        RecoveryPhaseGauge gauge;
        gauge.beg = 900;
        gauge.end = 6100;
        gauge.imported_phase_sets = {kLocalA, kLocal};
        RecoveryBlockGaugeVote vote{
            kRight, kLocal, {{{12, 0}, {0, 12}}}};
        vote.shared_candidate_same = 1;
        gauge.block_votes = {vote};
        Options stitch_opts;
        stitch_opts.min_block_link_reads = 1;
        stitch_opts.block_link_window = 8;
        const std::vector<RecoverySeam> windows = {
            {1000, 3000, kLeft, kRight},
        };
        const size_t weak_joined = stitch_recovery_phase_sets_left_to_right(
            replay, windows, {gauge}, stitch_opts);
        ok &= check(weak_joined == 0 &&
                    replay.candidates[2].phase_set == kLocal &&
                    replay.candidates[3].phase_set == kRight,
                    "recovery edge: one anchor cannot authorize boundary scope");

        gauge.block_votes[0].shared_candidate_same = 6;
        const size_t joined = stitch_recovery_phase_sets_left_to_right(
            replay, windows, {gauge}, stitch_opts);
        ok &= check(joined == 1 &&
                    replay.candidates[2].phase_set == kLocal &&
                    replay.candidates[3].phase_set == kLocal &&
                    replay.candidates.back().phase_set == kLocal,
                    "recovery edge: significant anchors permit bounded boundary solve");
    }

    // A strong local candidate gauge cannot override a complete-block tie.
    // The former any-failure retry narrowed this matrix and joined it; the
    // resource-aware flow retries only when the variable limit was exceeded.
    {
        constexpr hts_pos_t kLeft = 101;
        constexpr hts_pos_t kLocal = 202;
        constexpr hts_pos_t kRight = 303;
        PhasingChunk replay;
        const auto add_candidate = [&](hts_pos_t pos, hts_pos_t phase_set,
                                       bool injected,
                                       std::array<int, 3> consensus) {
            CandidateVariant candidate;
            candidate.key.pos = pos;
            candidate.key.type = VariantType::Snp;
            candidate.counts.ref_cov = 16;
            candidate.counts.alt_cov = 16;
            candidate.counts.allele_fraction = 0.5;
            candidate.hap_to_cons_alle = consensus;
            candidate.phase_set = phase_set;
            candidate.bam_injected = injected;
            replay.candidates.push_back(std::move(candidate));
        };
        add_candidate(1000, kLeft, false, {-1, 0, 1});
        add_candidate(2000, kLocal, true, {-1, 0, 1});
        add_candidate(3000, kRight, false, {-1, 0, 1});
        add_candidate(6000, kRight, false, {-1, 1, 0});

        const auto stable_fold = [](const std::string& value) {
            uint64_t hash = 14695981039346656037ULL;
            for (const unsigned char byte : value) {
                hash ^= byte;
                hash *= 1099511628211ULL;
            }
            return static_cast<int>(hash & 1ULL);
        };
        const auto add_read = [&](std::string name, int first,
                                  std::vector<int> alleles) {
            ReadRecord read;
            read.qname = std::move(name);
            read.mapq = 60;
            replay.reads.push_back(std::move(read));
            ReadVariantProfile profile;
            profile.read_id = static_cast<int>(replay.read_var_profile.size());
            profile.start_var_idx = first;
            profile.end_var_idx = first + static_cast<int>(alleles.size()) - 1;
            profile.alleles = std::move(alleles);
            profile.alt_qi.assign(profile.alleles.size(), 60);
            replay.read_var_profile.push_back(std::move(profile));
            replay.haps.push_back(0);
            replay.phase_sets.push_back(kUnphasedReadPhaseSet);
        };
        for (const std::string prefix : {"near", "far"}) {
            std::array<int, 2> fold_counts{};
            for (int copy = 0; fold_counts[0] < 16 || fold_counts[1] < 16;
                 ++copy) {
                const std::string name =
                    "full_tie_" + prefix + "_" + std::to_string(copy);
                const int fold = stable_fold(name);
                if (fold_counts[fold] >= 16) continue;
                const int allele = fold_counts[fold]++ & 1;
                add_read(name, 1,
                         prefix == "near"
                             ? std::vector<int>{allele, allele}
                             : std::vector<int>{allele, -1, allele});
            }
        }
        rebuild_read_var_cr(replay);

        RecoveryPhaseGauge gauge;
        gauge.beg = 900;
        gauge.end = 6100;
        gauge.imported_phase_sets = {kLocal};
        RecoveryBlockGaugeVote vote{
            kRight, kLocal, {{{16, 0}, {0, 16}}}};
        vote.shared_candidate_same = 6;
        gauge.block_votes = {vote};
        Options stitch_opts;
        stitch_opts.min_block_link_reads = 1;
        stitch_opts.block_link_window = 8;
        const size_t joined = stitch_recovery_phase_sets_left_to_right(
            replay, {{1000, 3000, kLeft, kRight}}, {gauge}, stitch_opts);
        ok &= check(joined == 0 &&
                    replay.candidates[1].phase_set == kLocal &&
                    replay.candidates[2].phase_set == kRight &&
                    replay.candidates[3].phase_set == kRight,
                    "recovery edge: full-block tie cannot trigger boundary scope");
    }

    // --- Compact explicit-seam replay for every selected chr20 target gap ---
    // The expensive discovery and allele matrices are reduced to the invariant
    // that failed in production: stitching must use the phase-set identities
    // captured by the canonical seam detector, even when candidate keys use a
    // shifted indel representation. Keep the previously fixed short gaps here
    // as regressions, then append every noncentromeric gap the current full
    // chr20 run leaves open while HiPhase crosses it at >=98% local truth
    // purity. The table keeps all selected coordinates in a sub-millisecond
    // test; integration measurements and support counts remain in evaluations/.
    {
        const std::vector<std::pair<hts_pos_t, hts_pos_t>> target_gaps = {
            // Previously fixed short-gap regressions.
            {1903021, 1903046}, {1907007, 1907045},
            {4778793, 4785720}, {32168699, 32169151},
            {34094605, 34102868}, {35613763, 35613765},
            {36059715, 36069697}, {55882617, 55883020},
            {58416190, 58418868}, {61763858, 61770862},
            {65509649, 65509898}, {65511683, 65513879},
            {65953554, 65953582}, {65994906, 65995588},

            // Current full-chr20 misses from remaining_hiphase_correct_gaps.tsv.
            {865572, 882277}, {882277, 890261},
            {4778792, 4785719}, {11495992, 11497008},
            {11796979, 11813446}, {14264549, 14272741},
            {14679241, 14679247}, {14719378, 14729749},
            {15023123, 15039543}, {15056025, 15071132},
            {18194808, 18218259}, {21669098, 21688374},
            {21823066, 21844359}, {25985709, 25986123},
            {30813841, 30814005}, {31886648, 31901501},
            {32035459, 32050364}, {32234664, 32246127},
            {32910395, 32910411}, {33792638, 33797964},
            {35281130, 35305645}, {35328965, 35347817},
            {38283561, 38303247}, {45230668, 45252114},
            {46727050, 46747598}, {46895038, 46896191},
            {55381221, 55381472}, {55381472, 55382729},
            {56064697, 56083708}, {57854341, 57866713},
            {58834248, 58836037}, {60971452, 60980031},
            {61373859, 61375344}, {62623253, 62642316},
        };
        ok &= check(target_gaps.size() == 48,
                    "target-gap replay covers 14 historical and 34 current coordinate cases");
        Options stitch_opts;
        stitch_opts.min_block_link_reads = 1;
        for (const auto& [beg, end] : target_gaps) {
            const hts_pos_t left_phase_set = beg + 100000000;
            const hts_pos_t right_phase_set = end + 200000000;
            PhasingChunk replay;
            CandidateVariant left;
            left.key.pos = beg - 7;  // deliberately differs from canonical beg
            left.hap_to_cons_alle = {-1, 0, 1};
            left.phase_set = left_phase_set;
            CandidateVariant right;
            right.key.pos = end + 9;  // deliberately differs from canonical end
            right.hap_to_cons_alle = {-1, 1, 0};
            right.phase_set = right_phase_set;
            replay.candidates = {left, right};

            const RecoverySeam seam{
                beg, end, left_phase_set, right_phase_set};
            RecoveryPhaseGauge gauge;
            gauge.beg = beg - 50000;
            gauge.end = end + 50000;
            gauge.graph_votes = {
                PhaseSetGaugeVote{left_phase_set, 3, 0},
                PhaseSetGaugeVote{right_phase_set, 3, 0},
            };
            const size_t joined = stitch_recovery_phase_sets_left_to_right(
                replay, {seam}, {gauge}, stitch_opts);
            ok &= check(joined == 1 &&
                        replay.candidates[0].phase_set == left_phase_set &&
                        replay.candidates[1].phase_set == left_phase_set,
                        "target-gap replay: explicit seam identities join " +
                            std::to_string(beg) + "-" + std::to_string(end));
        }
    }

    // A left-to-right merge can absorb the phase set named by the following
    // seam. The following targeted BAM solve votes for the surviving label, so
    // the stitcher must resolve the old detector ID before reading that vote.
    {
        constexpr hts_pos_t kLeft = 101;
        constexpr hts_pos_t kMiddle = 202;
        constexpr hts_pos_t kRight = 303;
        Options stitch_opts;
        stitch_opts.min_block_link_reads = 1;
        PhasingChunk replay;
        for (const auto [pos, phase_set] :
             {std::pair<hts_pos_t, hts_pos_t>{1000, kLeft},
              std::pair<hts_pos_t, hts_pos_t>{2000, kMiddle},
              std::pair<hts_pos_t, hts_pos_t>{3000, kRight}}) {
            CandidateVariant candidate;
            candidate.key.pos = pos;
            candidate.hap_to_cons_alle = {-1, 0, 1};
            candidate.phase_set = phase_set;
            replay.candidates.push_back(std::move(candidate));
        }
        const std::vector<RecoverySeam> seams = {
            {1000, 2000, kLeft, kMiddle},
            {2000, 3000, kMiddle, kRight},
        };
        RecoveryPhaseGauge first;
        first.beg = 900;
        first.end = 2100;
        first.graph_votes = {
            PhaseSetGaugeVote{kLeft, 4, 0},
            PhaseSetGaugeVote{kMiddle, 4, 0},
        };
        RecoveryPhaseGauge second;
        second.beg = 1900;
        second.end = 3100;
        // kMiddle has already been absorbed. This is the exact vote shape from
        // the consecutive 61.76 and 65.51 Mb recovery seams.
        second.graph_votes = {
            PhaseSetGaugeVote{kLeft, 4, 0},
            PhaseSetGaugeVote{kRight, 4, 0},
        };
        ok &= check(stitch_recovery_phase_sets_left_to_right(
                        replay, seams, {first, second}, stitch_opts) == 2 &&
                    std::all_of(replay.candidates.begin(),
                                replay.candidates.end(),
                                [](const CandidateVariant& candidate) {
                                    return candidate.phase_set == kLeft;
                                }),
                    "consecutive seam replay: surviving phase-set gauge is reused");
    }

    // Separate BAM phase sets have independent HP gauges even when they came
    // from one targeted solve. Direct per-block votes must take precedence over
    // a pooled solve-wide vote that happens to point in the opposite direction.
    {
        constexpr hts_pos_t kLeft = 101;
        constexpr hts_pos_t kLocal = 202;
        constexpr hts_pos_t kRight = 303;
        PhasingChunk replay;
        for (const auto [pos, phase_set, injected] : {
                 std::tuple<hts_pos_t, hts_pos_t, bool>{1000, kLeft, false},
                 std::tuple<hts_pos_t, hts_pos_t, bool>{2000, kLocal, true},
                 std::tuple<hts_pos_t, hts_pos_t, bool>{3000, kRight, false}}) {
            CandidateVariant candidate;
            candidate.key.pos = pos;
            candidate.hap_to_cons_alle = {-1, 0, 1};
            candidate.phase_set = phase_set;
            candidate.bam_injected = injected;
            replay.candidates.push_back(std::move(candidate));
        }
        RecoveryPhaseGauge gauge;
        gauge.beg = 900;
        gauge.end = 3100;
        gauge.imported_phase_sets = {kLocal};
        gauge.graph_votes = {
            PhaseSetGaugeVote{kLeft, 0, 10},
            PhaseSetGaugeVote{kRight, 10, 0},
        };
        gauge.block_votes = {
            RecoveryBlockGaugeVote{kLeft, kLocal, {{{8, 0}, {0, 8}}}},
            RecoveryBlockGaugeVote{kRight, kLocal, {{{8, 0}, {0, 8}}}},
        };
        Options stitch_opts;
        stitch_opts.min_block_link_reads = 1;
        PhasingChunk uncertain;
        uncertain.candidates = replay.candidates;
        RecoveryPhaseGauge uncertain_gauge = gauge;
        uncertain_gauge.graph_votes.clear();
        uncertain_gauge.block_votes = {
            RecoveryBlockGaugeVote{kLeft, kLocal, {{{52, 48}, {48, 52}}}},
            RecoveryBlockGaugeVote{kRight, kLocal, {{{52, 48}, {48, 52}}}},
        };
        ok &= check(stitch_recovery_phase_sets_left_to_right(
                        uncertain, {{1000, 3000, kLeft, kRight}},
                        {uncertain_gauge}, stitch_opts) == 0 &&
                    uncertain.candidates[1].phase_set == kLocal &&
                    uncertain.candidates[2].phase_set == kRight,
                    "recovery block gauge: high-depth near-tie abstains");

        // Shared clean candidates validate a decisive read orientation but
        // never create a join on their own. This prevents an internal switch in
        // a long BAM block from being promoted through two distant anchors.
        PhasingChunk candidate_anchored;
        candidate_anchored.candidates = replay.candidates;
        RecoveryPhaseGauge candidate_gauge;
        candidate_gauge.beg = 900;
        candidate_gauge.end = 3100;
        candidate_gauge.imported_phase_sets = {kLocal};
        RecoveryBlockGaugeVote left_anchor{
            kLeft, kLocal, {{{4, 0}, {0, 4}}}};
        left_anchor.shared_candidate_same = 1;
        RecoveryBlockGaugeVote right_anchor{
            kRight, kLocal, {{{4, 0}, {0, 4}}}};
        right_anchor.shared_candidate_same = 1;
        candidate_gauge.block_votes = {left_anchor, right_anchor};
        ok &= check(stitch_recovery_phase_sets_left_to_right(
                        candidate_anchored,
                        {{1000, 3000, kLeft, kRight}},
                        {candidate_gauge}, stitch_opts) == 0 &&
                    candidate_anchored.candidates[0].phase_set == kLeft &&
                    candidate_anchored.candidates[1].phase_set == kLocal &&
                    candidate_anchored.candidates[2].phase_set == kRight,
                    "recovery block gauge: one-block validation is atomic");

        PhasingChunk conflicting_candidates;
        conflicting_candidates.candidates = replay.candidates;
        left_anchor.shared_candidate_cross = 1;
        right_anchor.shared_candidate_cross = 1;
        candidate_gauge.block_votes = {left_anchor, right_anchor};
        ok &= check(stitch_recovery_phase_sets_left_to_right(
                        conflicting_candidates,
                        {{1000, 3000, kLeft, kRight}},
                        {candidate_gauge}, stitch_opts) == 0,
                    "recovery block gauge: conflicting shared candidates abstain");

        const size_t joined = stitch_recovery_phase_sets_left_to_right(
            replay, {{1000, 3000, kLeft, kRight}}, {gauge}, stitch_opts);
        ok &= check(joined == 0 &&
                    replay.candidates[0].phase_set == kLeft &&
                    replay.candidates[1].phase_set == kLocal &&
                    replay.candidates[2].phase_set == kRight,
                    "recovery block gauge: incomplete one-block bridge changes nothing");
    }

    // Alleles shared by a graph block and an imported block can be strongly
    // correlated even when their phase gauges are incompatible. Without the
    // source-specific 2x2 vote, keep the BAM block local.
    {
        constexpr hts_pos_t kLeft = 101;
        constexpr hts_pos_t kLocal = 202;
        PhasingChunk replay;
        for (const auto [pos, phase_set, injected] : {
                 std::tuple<hts_pos_t, hts_pos_t, bool>{1000, kLeft, false},
                 {2000, kLocal, true}}) {
            CandidateVariant candidate;
            candidate.key.pos = pos;
            candidate.key.type = VariantType::Snp;
            candidate.counts.n_uniq_alles = 2;
            candidate.counts.alle_covs = {8, 8};
            candidate.lcd_var_i_to_cate = kCandCleanHetSnp;
            candidate.hap_to_cons_alle = {-1, 0, 1};
            candidate.phase_set = phase_set;
            candidate.bam_injected = injected;
            replay.candidates.push_back(std::move(candidate));
        }
        for (int allele = 0; allele <= 1; ++allele) {
            for (int copy = 0; copy < 8; ++copy) {
                ReadRecord read;
                read.qname = "aggregate_only_" + std::to_string(allele) +
                             "_" + std::to_string(copy);
                read.mapq = 60;
                replay.reads.push_back(std::move(read));
                ReadVariantProfile profile;
                profile.read_id = static_cast<int>(
                    replay.read_var_profile.size());
                profile.start_var_idx = 0;
                profile.end_var_idx = 1;
                profile.alleles = {allele, allele};
                profile.alt_qi = {60, 60};
                replay.read_var_profile.push_back(std::move(profile));
                replay.haps.push_back(0);
                replay.phase_sets.push_back(kUnphasedReadPhaseSet);
            }
        }
        rebuild_read_var_cr(replay);
        RecoveryPhaseGauge gauge;
        gauge.beg = 900;
        gauge.end = 2100;
        gauge.imported_phase_sets = {kLocal};
        Options stitch_opts;
        stitch_opts.min_block_link_reads = 1;
        stitch_opts.block_link_window = 8;
        ok &= check(stitch_recovery_phase_sets_left_to_right(
                        replay, {{1000, 2000, kLeft, kLocal}},
                        {gauge}, stitch_opts) == 0 &&
                    replay.candidates[0].phase_set == kLeft &&
                    replay.candidates[1].phase_set == kLocal,
                    "recovery blocks: aggregate-only graph attachment abstains");
    }

    // Multiple BAM phase blocks from one subsolve begin with unrelated HP
    // gauges. A final stitch may join them only when reads spanning that exact
    // adjacent pair pass the parity test; each outer block independently uses
    // its own graph-flank vote.
    {
        constexpr hts_pos_t kLeft = 101;
        constexpr hts_pos_t kLocalA = 202;
        constexpr hts_pos_t kLocalB = 303;
        constexpr hts_pos_t kRight = 404;
        PhasingChunk replay;
        for (const auto [pos, phase_set, injected] : {
                 std::tuple<hts_pos_t, hts_pos_t, bool>{1000, kLeft, false},
                 {2000, kLocalA, true},
                 {3000, kLocalB, true},
                 {4000, kRight, false}}) {
            CandidateVariant candidate;
            candidate.key.pos = pos;
            candidate.key.type = VariantType::Snp;
            candidate.counts.n_uniq_alles = 2;
            candidate.counts.alle_covs = {16, 16};
            candidate.lcd_var_i_to_cate = kCandCleanHetSnp;
            candidate.hap_to_cons_alle = {-1, 0, 1};
            candidate.phase_set = phase_set;
            candidate.bam_injected = injected;
            replay.candidates.push_back(std::move(candidate));
        }
        const auto add_read = [&](std::string name, int first,
                                  std::vector<int> alleles, int hap,
                                  hts_pos_t phase_set) {
            ReadRecord read;
            read.qname = std::move(name);
            read.mapq = 60;
            replay.reads.push_back(std::move(read));
            ReadVariantProfile profile;
            profile.read_id = static_cast<int>(
                replay.read_var_profile.size());
            profile.start_var_idx = first;
            profile.end_var_idx = first +
                static_cast<int>(alleles.size()) - 1;
            profile.alleles = std::move(alleles);
            profile.alt_qi.assign(profile.alleles.size(), 60);
            replay.read_var_profile.push_back(std::move(profile));
            replay.haps.push_back(hap);
            replay.phase_sets.push_back(phase_set);
        };
        add_read("local_a", 1, {0}, 1, kLocalA);
        add_read("local_b", 2, {0}, 1, kLocalB);
        add_read("right", 3, {0}, 1, kRight);
        for (int allele = 0; allele <= 1; ++allele) {
            for (int copy = 0; copy < 8; ++copy) {
                add_read("between_" + std::to_string(allele) + "_" +
                             std::to_string(copy),
                         1, {allele, allele}, 0,
                         kUnphasedReadPhaseSet);
                // The complete imported chain is allowed to commit only when
                // an independent read relation between the graph flanks agrees
                // with its accumulated parity.
                add_read("outer_" + std::to_string(allele) + "_" +
                             std::to_string(copy),
                         0, {allele, -1, -1, 1 - allele}, 0,
                         kUnphasedReadPhaseSet);
            }
        }
        rebuild_read_var_cr(replay);

        RecoveryPhaseGauge gauge;
        gauge.beg = 900;
        gauge.end = 4100;
        gauge.imported_phase_sets = {kLocalA, kLocalB};
        RecoveryBlockGaugeVote left_vote{
            kLeft, kLocalA, {{{8, 0}, {0, 8}}}};
        left_vote.shared_candidate_same = 6;
        RecoveryBlockGaugeVote right_vote{
            kRight, kLocalB, {{{0, 8}, {8, 0}}}};
        right_vote.shared_candidate_cross = 6;
        gauge.block_votes = {left_vote, right_vote};
        Options stitch_opts;
        stitch_opts.min_block_link_reads = 1;
        stitch_opts.block_link_window = 8;
        const size_t joined = stitch_recovery_phase_sets_left_to_right(
            replay, {{1000, 4000, kLeft, kRight}}, {gauge}, stitch_opts);
        ok &= check(joined == 3 &&
                    replay.candidates[0].phase_set == kLeft &&
                    replay.candidates[1].phase_set == kLeft &&
                    replay.candidates[2].phase_set == kLeft &&
                    replay.candidates[3].phase_set == kLeft,
                    "recovery blocks: direct evidence stitches independent blocks");
        ok &= check(replay.phase_sets[0] == kLeft &&
                    replay.haps[0] == 1 &&
                    replay.phase_sets[1] == kLeft &&
                    replay.haps[1] == 1 &&
                    replay.phase_sets[2] == kLeft &&
                    replay.haps[2] == 2,
                    "recovery blocks: flank attachment preserves and flips reads atomically");
    }

    // Reads with only excluded-site observations are invisible to the clean
    // solve. Seven already phased molecules orient two sites statistically;
    // an eighth molecule can then be assigned without changing the candidates or
    // joining phase sets.
    {
        PhasingChunk rescue;
        CandidateVariant excluded;
        excluded.key.pos = 1000;
        excluded.key.type = VariantType::Insertion;
        excluded.counts.n_uniq_alles = 2;
        excluded.counts.category = VariantCategory::RepeatHetIndel;
        excluded.lcd_var_i_to_cate = kLongcalldRepHetVar;
        rescue.candidates.push_back(excluded);
        excluded.key.pos = 1100;
        rescue.candidates.push_back(excluded);

        constexpr hts_pos_t kPhaseSet = 900;
        for (int read_i = 0; read_i < 8; ++read_i) {
            ReadRecord read;
            read.qname = "rescue_" + std::to_string(read_i);
            rescue.reads.push_back(std::move(read));

            ReadVariantProfile profile;
            profile.read_id = read_i;
            profile.start_var_idx = 0;
            profile.end_var_idx = 1;
            profile.alleles = {read_i < 3 ? 0 : 1,
                               read_i < 3 ? 0 : 1};
            rescue.read_var_profile.push_back(std::move(profile));

            rescue.haps.push_back(read_i < 3 ? 1 :
                                  read_i < 7 ? 2 : 0);
            rescue.phase_sets.push_back(read_i < 7 ? kPhaseSet :
                                        kUnphasedReadPhaseSet);
        }

        const size_t rescued = rescue_unphased_graph_reads(rescue);
        ok &= check(rescued == 1 && rescue.haps[7] == 0 &&
                    rescue.gap_haps[7] == 2 &&
                    rescue.gap_phase_sets[7] ==
                        kPhaseSet + kGapFillPsOffset,
                    "graph read rescue: oriented excluded site tags read");
        ok &= check(rescue.candidates[0].phase_set ==
                        kUnsetCandidatePhaseSet &&
                    rescue.candidates[1].phase_set ==
                        kUnsetCandidatePhaseSet,
                    "graph read rescue: candidate phasing is unchanged");
        ok &= check(rescue_unphased_graph_reads(rescue) == 0 &&
                    rescue.gap_phase_sets[7] ==
                        kPhaseSet + kGapFillPsOffset,
                    "graph read rescue: fixed point preserves one PS offset");
    }

    // A read can appear in two adjacent chunks while only the upstream chunk
    // has informative alleles. The downstream visit must not erase that valid
    // HP/PS assignment merely because it is unphased. A later phased visit
    // still owns the read, preserving ordinary downstream ownership.
    {
        std::unordered_map<std::string, PhaseReadOutputRow> output_rows;
        merge_graph_chunk_into_read_rows(output_rows, chunks[0], 0);
        const PhaseReadOutputRow initial = output_rows.at("read_a");
        ok &= check(initial.has_phased_assignment,
                    "graph output merge: fixture starts phased");

        const std::vector<int> saved_haps = chunks[0].chunk.haps;
        const std::vector<hts_pos_t> saved_phase_sets =
            chunks[0].chunk.phase_sets;
        const int saved_chunk_id = chunks[0].chunk.region.chunk_id;
        std::fill(chunks[0].chunk.haps.begin(),
                  chunks[0].chunk.haps.end(), 0);
        std::fill(chunks[0].chunk.phase_sets.begin(),
                  chunks[0].chunk.phase_sets.end(),
                  kUnphasedReadPhaseSet);
        chunks[0].chunk.region.chunk_id = saved_chunk_id + 1;
        merge_graph_chunk_into_read_rows(output_rows, chunks[0], 0);
        ok &= check(output_rows.at("read_a").hap == initial.hap &&
                    output_rows.at("read_a").phase_set == initial.phase_set &&
                    output_rows.at("read_a").has_phased_assignment,
                    "graph output merge: unphased overlap preserves assignment");

        const int replacement_hap = initial.hap == 1 ? 2 : 1;
        chunks[0].chunk.gap_haps.assign(chunks[0].chunk.reads.size(), 0);
        chunks[0].chunk.gap_phase_sets.assign(
            chunks[0].chunk.reads.size(), kUnphasedReadPhaseSet);
        chunks[0].chunk.gap_haps[0] = replacement_hap;
        chunks[0].chunk.gap_phase_sets[0] = initial.phase_set + 500;
        merge_graph_chunk_into_read_rows(output_rows, chunks[0], 0);
        ok &= check(output_rows.at("read_a").hap == initial.hap &&
                    output_rows.at("read_a").phase_set == initial.phase_set,
                    "graph output merge: rescue cannot replace primary assignment");
        chunks[0].chunk.gap_haps.clear();
        chunks[0].chunk.gap_phase_sets.clear();

        chunks[0].chunk.haps[0] = replacement_hap;
        chunks[0].chunk.phase_sets[0] = initial.phase_set + 1000;
        chunks[0].chunk.region.chunk_id = saved_chunk_id + 2;
        merge_graph_chunk_into_read_rows(output_rows, chunks[0], 0);
        ok &= check(output_rows.at("read_a").hap == replacement_hap &&
                    output_rows.at("read_a").phase_set ==
                        initial.phase_set + 1000,
                    "graph output merge: later phased overlap owns assignment");

        chunks[0].chunk.haps = saved_haps;
        chunks[0].chunk.phase_sets = saved_phase_sets;
        chunks[0].chunk.region.chunk_id = saved_chunk_id;
    }

    std::ostringstream sites;
    write_graph_phase_sites_tsv(sites, chunks);
    ok &= check(sites.str().find("HAP1_ALLELE") != std::string::npos,
                "BAM-derived graph phase site output");
    if (ok) {
        std::cout << "ALL PASS\n";
        return 0;
    }
    return 1;
}
