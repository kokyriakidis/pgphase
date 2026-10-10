#include "graph_bam_adapter.hpp"

#include "collect_phase.hpp"
#include "collect_output.hpp"
#include "collect_var.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <memory>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include <htslib/faidx.h>
#include <htslib/sam.h>
#include <unistd.h>

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

static bool test_tandem_insertion_motif_length() {
    bool ok = true;
    ok &= check(tandem_insertion_motif_length("", "GAAGGAAG") == 4, "complete source deletion uses a four-base fundamental motif");
    ok &= check(tandem_insertion_motif_length("TTCC", "TTCCTTCC") == 4, "four-base repeat motif is callable");
    ok &= check(tandem_insertion_motif_length("", "TTTAT") == 5, "reference versus five-base repeat is callable");
    ok &= check(tandem_insertion_motif_length("AT", "ATATAT") == 2, "short motif uses its fundamental period");
    ok &= check(tandem_insertion_motif_length("ATATATAGAG", "ATATATATAGAGAG") == 0, "compound MSA contrast is not a pure repeat");
    ok &= check(tandem_insertion_motif_length("TTCC", "TTCCAAAA") == 0, "incompatible motifs abstain");
    ok &= check(tandem_insertion_motif_length("TTCC", "TTCC") == 0, "identical insertion alleles abstain");
    return ok;
}

static bool test_joint_candidate_loci() {
    bool ok = true;
    const std::string reference = "CAAATG";
    const auto base = [&](hts_pos_t pos) {
        return pos >= 1 && pos <= static_cast<hts_pos_t>(reference.size())
            ? reference[static_cast<size_t>(pos - 1)] : 'N';
    };
    const auto fixture = [&](size_t copies) {
        GraphChunkBuildResult gc;
        gc.chunk.region.tid = 0;
        gc.chunk.ref_beg = 1;
        gc.chunk.ref_end = reference.size();
        for (size_t ci = 0; ci < copies; ++ci) {
            CandidateVariant row;
            row.key = vcf_to_variant_key(0, std::min<hts_pos_t>(2 + ci, 4), "A", "AA");
            row.bam_injected = ci % 2 != 0;
            row.graph_site = !row.bam_injected;
            row.msa_verified = row.bam_injected;
            if (row.graph_site) row.key.alt = ">graph" + std::to_string(ci);
            row.counts.category = VariantCategory::CleanHetIndel;
            row.counts.n_uniq_alles = 2;
            row.counts.ref_cov = 4;
            row.counts.alt_cov = 4;
            row.counts.total_cov = 8;
            row.counts.alle_covs = {4, 4};
            row.lcd_var_i_to_cate = kCandCleanHetIndel;
            gc.chunk.candidates.push_back(row);
            gc.site_ids.push_back("source" + std::to_string(ci));
            GraphSiteMeta meta;
            meta.chrom = "chr1";
            meta.pos = std::min<hts_pos_t>(2 + ci, 4);
            meta.ref = "A";
            meta.alts = {"AA"};
            gc.site_meta.push_back(meta);
            gc.site_allele_orig_idx.push_back({0, 1});
        }
        for (int ri = 0; ri < 8; ++ri) {
            ReadRecord read;
            read.qname = "read" + std::to_string(ri);
            read.beg = 1;
            read.end = reference.size();
            gc.chunk.reads.push_back(std::move(read));
            ReadVariantProfile profile;
            profile.read_id = ri;
            profile.start_var_idx = 0;
            profile.end_var_idx = copies - 1;
            profile.alleles.assign(copies, ri % 2);
            profile.alt_qi.assign(copies, 42);
            profile.graph_alleles = profile.alleles;
            profile.bam_alleles = profile.alleles;
            profile.bam_qi.assign(copies, 43);
            profile.bam_base_qualities.assign(copies, 31);
            gc.chunk.read_var_profile.push_back(profile);
        }
        rebuild_read_var_cr(gc.chunk);
        return gc;
    };
    Options opts;
    auto single = fixture(1);
    phase_joint_graph_candidates(single, opts, base);
    for (const size_t copies : {2, 4}) {
        auto gc = fixture(copies);
        const auto profiles = gc.chunk.read_var_profile;
        const auto candidates = gc.chunk.candidates;
        const auto ids = gc.site_ids;
        phase_joint_graph_candidates(gc, opts, base);
        ok &= check(gc.joint_candidate_loci.size() == 1 &&
                    gc.joint_candidate_loci[0].second.size() == copies,
                    "joint locus: graph and physical BAM aliases share one complete contrast");
        ok &= check(gc.chunk.haps == single.chunk.haps && gc.chunk.phase_sets == single.chunk.phase_sets,
                    "joint locus: additional descriptions cannot change the single-locus read phase");
        for (size_t ri = 0; ri < profiles.size(); ++ri) {
            const auto& p = gc.chunk.read_var_profile[ri];
            const auto& old = profiles[ri];
            ok &= check(gc.chunk.reads[ri].n_vars_scored == single.chunk.reads[ri].n_vars_scored &&
                        gc.chunk.reads[ri].hap_score_margin == single.chunk.reads[ri].hap_score_margin,
                        "joint locus: aliases cannot increase vote count or haplotype margin");
            ok &= check(p.alleles == old.alleles && p.alt_qi == old.alt_qi &&
                        p.graph_alleles == old.graph_alleles && p.bam_alleles == old.bam_alleles &&
                        p.bam_qi == old.bam_qi && p.bam_base_qualities == old.bam_base_qualities,
                        "joint locus: all original observation channels and quality certificates survive");
        }
        for (size_t ci = 0; ci < copies; ++ci) {
            const auto& row = gc.chunk.candidates[ci];
            const auto& old = candidates[ci];
            ok &= check(row.key.alt == old.key.alt && row.key.pos == old.key.pos &&
                        row.bam_injected == old.bam_injected && row.graph_site == old.graph_site &&
                        row.msa_verified == old.msa_verified && row.counts.alle_covs == old.counts.alle_covs &&
                        row.lcd_var_i_to_cate == old.lcd_var_i_to_cate && gc.site_ids == ids,
                        "joint locus: physical/topology identities, depths and provenance remain distinct");
            ok &= check(row.hap_to_cons_alle == gc.chunk.candidates[0].hap_to_cons_alle &&
                        row.phase_set == gc.chunk.candidates[0].phase_set,
                        "joint locus: aliases receive one consistently oriented result");
        }
    }
    auto conflict = fixture(2);
    conflict.chunk.read_var_profile[0].alleles = {0, 1};
    rebuild_read_var_cr(conflict.chunk);
    phase_joint_graph_candidates(conflict, opts, base);
    ok &= check(conflict.chunk.haps[0] == 0 && conflict.chunk.reads[0].n_vars_scored == 0,
                "joint locus: opposing descriptions abstain instead of choosing the first call");
    auto exact = fixture(2);
    exact.site_meta[1] = exact.site_meta[0];
    exact.chunk.candidates[1].key = vcf_to_variant_key(0, 2, "A", "AA");
    auto& sparse = exact.chunk.read_var_profile[1];
    sparse.start_var_idx = sparse.end_var_idx = 1;
    for (auto* values : {&sparse.alleles, &sparse.alt_qi, &sparse.graph_alleles,
                        &sparse.bam_alleles, &sparse.bam_qi}) values->erase(values->begin());
    sparse.bam_base_qualities.erase(sparse.bam_base_qualities.begin());
    rebuild_read_var_cr(exact.chunk);
    phase_joint_graph_candidates(exact, opts, base);
    ok &= check(exact.chunk.haps[1] != 0 && exact.chunk.read_var_profile[1].start_var_idx == 1 &&
                exact.chunk.read_var_profile[1].alleles == std::vector<int>{1},
                "joint locus: alias-only sparse coverage votes once and restores its original extent");
    auto filtered = fixture(2);
    filtered.chunk.candidates[1].lcd_var_i_to_cate = kCandNoisyCandHet;
    filtered.chunk.candidates[1].counts.category = VariantCategory::NoisyCandHet;
    phase_joint_graph_candidates(filtered, opts, base);
    ok &= check(filtered.chunk.candidates[1].phase_set <= 0 &&
                filtered.chunk.candidates[1].lcd_var_i_to_cate == kCandNoisyCandHet,
                "joint locus: identity cannot promote an excluded noisy row into the clean solve");
    auto multiallelic = fixture(2);
    multiallelic.site_meta[0].alts.push_back("AAA");
    multiallelic.site_meta[1].alts.push_back("AAA");
    opts.snarl_allele_phasing = true;
    phase_joint_graph_candidates(multiallelic, opts, base);
    ok &= check(multiallelic.chunk.reads[0].n_vars_scored == 2,
                "joint locus: non-selected ALT classes cannot claim an actual reference contrast");
    auto partial = fixture(2);
    partial.chunk.reads[1].beg = 4;
    opts.snarl_allele_phasing = false;
    phase_joint_graph_candidates(partial, opts, base);
    ok &= check(partial.chunk.reads[0].n_vars_scored == 2 &&
                partial.chunk.read_var_profile[1].alleles == std::vector<int>({1, 1}),
                "joint locus: incomplete molecule coverage cannot project a shifted alias leftward");
    auto unsupported = fixture(2);
    unsupported.chunk.read_var_profile[1].alleles[0] = -1;
    unsupported.chunk.read_var_profile[1].graph_alleles[0] = -1;
    phase_joint_graph_candidates(unsupported, opts, base);
    ok &= check(unsupported.chunk.reads[0].n_vars_scored == 2 &&
                unsupported.chunk.read_var_profile[1].alleles == std::vector<int>({-1, 1}),
                "joint locus: a different allele context cannot certify a missing parent REF class");
    conflict.site_allele_orig_idx[1] = {1, 2};
    rebuild_joint_candidate_loci(conflict, base);
    ok &= check(conflict.joint_candidate_loci.empty(), "joint locus: ALT/ALT is not a REF/ALT alias");
    conflict.site_allele_orig_idx[1] = {0, 1};
    conflict.site_meta[1].alts = {"AAA"};
    rebuild_joint_candidate_loci(conflict, base);
    ok &= check(conflict.joint_candidate_loci.empty(), "joint locus: distinct repeat lengths remain separate");
    conflict.site_meta[1].alts = {"AA"};
    conflict.site_meta[1].ref = "C";
    rebuild_joint_candidate_loci(conflict, base);
    ok &= check(conflict.joint_candidate_loci.empty(), "joint locus: an invalid reference cannot join aliases");
    return ok;
}

int main(int argc, char** argv) {
    if (argc == 2 && std::strcmp(argv[1], "--joint-loci") == 0) {
        const bool ok = test_joint_candidate_loci();
        if (ok) std::cout << "ALL PASS\n";
        return ok ? 0 : 1;
    }
    bool ok = true;
    ok &= test_tandem_insertion_motif_length();
    ok &= test_joint_candidate_loci();

    {
        GraphSiteCatalog catalog;
        GraphSite site = make_site("snp_branch", 100, ">1>2>3>4>5", ">1>8>3>4>5");
        site.ref = "A";
        site.alts = {"C", "CT", "AT", "G", "CC", "T"};
        for (const std::string walk : {">1>8>3>9>5", ">1>2>3>9>5",
                                        ">1>6>3>9>5", ">1>8>3>1>8>3",
                                        ">1<8>3>9>5"})
            site.allele_walks.push_back(parse_graph_walk(walk));
        catalog.sites.push_back(site);
        const auto make_chunk = [&] {
            GraphChunkBuildResult graph;
            CandidateVariant snp;
            snp.counts.n_uniq_alles = 2;
            snp.counts.category = VariantCategory::CleanHetSnp;
            snp.counts.alle_covs = {10, 10};
            snp.counts.ref_cov = 10;
            snp.counts.alt_cov = 10;
            snp.phase_set = 100;
            snp.hap_to_cons_alle = {-1, 0, 1};
            graph.chunk.candidates = {snp};
            graph.site_ids = {"snp_branch:1"};
            graph.site_meta = {{"chr1", 100, "A", site.alts}};
            graph.site_allele_orig_idx = {{0, 1}};
            for (const std::string name : {"alt", "ref", "unknown", "repeat", "reverse_node",
                                            "duplicate", "conflict", "same_branch_conflict",
                                            "low_mapq", "known"}) {
                ReadRecord read;
                read.qname = name;
                read.mapq = 60;
                graph.chunk.reads.push_back(std::move(read));
                ReadVariantProfile profile;
                profile.read_id = static_cast<int>(graph.chunk.reads.size() - 1);
                profile.start_var_idx = 0;
                profile.end_var_idx = 0;
                profile.alleles = {name == "known" ? 0 : -1};
                profile.alt_qi = {-1};
                graph.chunk.read_var_profile.push_back(profile);
            }
            graph.chunk.haps.assign(graph.chunk.reads.size(), 1);
            graph.chunk.phase_sets.assign(graph.chunk.reads.size(), 100);
            return graph;
        };
        const auto row = [](const std::string& name, int allele, int mapq = 60) {
            return GraphReadAllele{"snp_branch", "chr1", 100, name, allele, mapq, false};
        };
        const std::vector<GraphReadAllele> rows{
            row("alt", 2), row("ref", 3), row("unknown", 4), row("repeat", 5),
            row("reverse_node", 6), row("duplicate", 2), row("duplicate", 2),
            row("conflict", 2), row("conflict", 3),
            row("same_branch_conflict", 1), row("same_branch_conflict", 2),
            row("low_mapq", 2, 5), row("known", 2)};
        Options opts;
        opts.min_mapq = 30;
        GraphChunkBuildResult graph = make_chunk();
        ok &= check(supplement_phased_snp_branches(catalog.view_all(), rows, graph, opts) == 3,
                    "SNP branch: only unique eligible missing calls are added once");
        const auto allele = [&](const std::string& name) {
            return profile_for_read(graph.chunk, name)->alleles[0];
        };
        ok &= check(allele("alt") == 1 && allele("ref") == 0 && allele("duplicate") == 1,
                    "SNP branch: exact local branch survives another catalog indel");
        for (const std::string name : {"unknown", "repeat", "reverse_node", "conflict",
                                        "same_branch_conflict", "low_mapq"})
            ok &= check(allele(name) == -1, "SNP branch: ambiguous calls abstain: " + name);
        ok &= check(allele("known") == 0 && graph.chunk.candidates[0].counts.alle_covs ==
                    std::vector<int>{10, 10} && graph.chunk.candidates[0].phase_set == 100 &&
                    graph.chunk.candidates[0].hap_to_cons_alle == std::array<int, 3>{-1, 0, 1} &&
                    graph.chunk.haps == std::vector<int>(10, 1) &&
                    graph.chunk.phase_sets == std::vector<hts_pos_t>(10, 100),
                    "SNP branch: preserve existing calls, genotype counts and gauges");
        ok &= check(supplement_phased_snp_branches(catalog.view_all(), rows, graph, opts) == 0,
                    "SNP branch: applying supplementation twice is idempotent");
        graph = make_chunk();
        auto reversed = rows;
        std::reverse(reversed.begin(), reversed.end());
        ok &= check(supplement_phased_snp_branches(catalog.view_all(), reversed, graph, opts) == 3 &&
                    allele("conflict") == -1 && allele("same_branch_conflict") == -1,
                    "SNP branch: source order cannot choose a conflicted allele");
        graph = make_chunk();
        Options narrowed = opts;
        narrowed.max_af = 0.52;
        ok &= check(supplement_phased_snp_branches(catalog.view_all(), rows, graph, narrowed) == 0 &&
                    allele("alt") == -1 && allele("ref") == -1,
                    "SNP branch: evidence contradicting the frozen genotype cannot strengthen it");
        graph = make_chunk();
        graph.chunk.candidates[0].phase_set = 0;
        ok &= check(supplement_phased_snp_branches(catalog.view_all(), rows, graph, opts) == 0,
                    "SNP branch: unphased candidates cannot gain stitching evidence");
        graph = make_chunk();
        catalog.sites[0].conditional_parent_alleles = {1};
        ok &= check(supplement_phased_snp_branches(catalog.view_all(), rows, graph, opts) == 0,
                    "SNP branch: conditional child cannot bypass parent gating");
    }


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
    ok &= check(std::all_of(chunks[0].chunk.candidates.begin(),
                            chunks[0].chunk.candidates.end(),
                            [](const CandidateVariant& candidate) {
                                return candidate.graph_site;
                            }),
                "adapter marks graph-origin candidates");
    ok &= check(chunks[0].chunk.reads.size() == 4, "adapter builds four reads");
    ok &= check(chunks[0].chunk.read_var_profile.size() == 4, "adapter builds read profiles");
    ok &= check(chunks[0].chunk.read_var_cr != nullptr, "adapter builds read-var cgranges");
    ok &= check(chunks[0].chunk.candidates[0].phase_set == kUnsetCandidatePhaseSet,
                "graph candidate uses longcallD unset phase-set sentinel");
    ok &= check(chunks[0].chunk.phase_sets[0] == kUnphasedReadPhaseSet,
                "graph read uses longcallD unphased phase-set sentinel");

    {
        GraphSiteCatalog physical_catalog = catalog;
        physical_catalog.sites[0].bam_homozygous_alt = true;
        auto physical = build_graph_chunk(physical_catalog.view_all(), rows,
            "chr1", 0, 300, 0, build_opts);
        const auto original_categories = physical.chunk.candidates;
        assign_hap_based_on_germline_het_vars_kmeans(physical.chunk, build_opts, kCandGermlineClean);
        const auto original_haps = physical.chunk.haps;
        const auto original_phase_sets = physical.chunk.phase_sets;
        const auto original_neighbor = physical.chunk.candidates[1];
        reclassify_physically_validated_graph_snps(physical_catalog.view_all(), physical);
        ok &= check(physical.chunk.candidates.size() == 2 &&
                    physical.chunk.candidates[0].counts.ref_cov == 2 &&
                    physical.chunk.candidates[0].counts.alt_cov == 2,
                    "physical homozygote: catalog counts and neighbor survive rebuilding");
        ok &= check(physical.chunk.candidates[0].counts.category == VariantCategory::CleanHom &&
                    physical.chunk.candidates[0].lcd_var_i_to_cate == kCandCleanHom &&
                    !is_phase_set_anchor(physical.chunk.candidates[0]) &&
                    physical.chunk.candidates[0].hap_to_cons_alle[1] == 1 &&
                    physical.chunk.candidates[0].hap_to_cons_alle[2] == 1 &&
                    physical.chunk.candidates[1].counts.category == VariantCategory::CleanHetSnp,
                    "physical homozygote: graph walk segregation cannot recreate heterozygosity");
        ok &= check(physical.chunk.haps == original_haps &&
                    physical.chunk.phase_sets == original_phase_sets &&
                    physical.chunk.candidates[1].hap_to_cons_alle == original_neighbor.hap_to_cons_alle,
                    "physical homozygote: surviving anchors and reads retain their initial gauge");
        GraphSiteCatalog deletion_catalog = catalog;
        deletion_catalog.sites[0].bam_alt_deletion_no_ref = true;
        auto deletion = build_graph_chunk(deletion_catalog.view_all(), rows,
            "chr1", 0, 300, 0, build_opts);
        assign_hap_based_on_germline_het_vars_kmeans(deletion.chunk, build_opts, kCandGermlineClean);
        const auto deletion_neighbor = deletion.chunk.candidates[1];
        const auto deletion_haps = deletion.chunk.haps;
        reclassify_physically_validated_graph_snps(deletion_catalog.view_all(), deletion);
        ok &= check(deletion.site_meta[0].bam_alt_deletion_no_ref &&
                    !is_phase_set_anchor(deletion.chunk.candidates[0]) &&
                    deletion.chunk.candidates[0].counts.ref_cov == 2 &&
                    deletion.chunk.candidates[1].hap_to_cons_alle == deletion_neighbor.hap_to_cons_alle &&
                    deletion.chunk.haps == deletion_haps,
                    "REF-absent deletion: demote after solving without changing surviving gauges or evidence");
        GraphSiteCatalog minor_catalog = catalog;
        minor_catalog.sites[0].bam_low_fraction_snp = true;
        auto minor = build_graph_chunk(minor_catalog.view_all(), rows,
            "chr1", 0, 300, 0, build_opts);
        assign_hap_based_on_germline_het_vars_kmeans(minor.chunk, build_opts, kCandGermlineClean);
        const auto minor_haps = minor.chunk.haps;
        const auto minor_neighbor = minor.chunk.candidates[1];
        minor.site_ids[0] += ":4";
        CandidateVariant retained_deletion = minor.chunk.candidates[0];
        retained_deletion.key.type = VariantType::Deletion;
        retained_deletion.counts.category = VariantCategory::CleanHetIndel;
        minor.chunk.candidates.push_back(retained_deletion);
        minor.site_ids.push_back(graph_site_key_str(minor_catalog.sites[0]) + ":3");
        minor.site_meta.push_back(minor.site_meta[0]);
        reclassify_physically_validated_graph_snps(minor_catalog.view_all(), minor);
        ok &= check(minor.site_meta[0].bam_low_fraction_snp &&
                    minor.chunk.candidates[0].counts.category == VariantCategory::LowAlleleFraction &&
                    !is_phase_set_anchor(minor.chunk.candidates[0]),
                    "physical minor SNP: decomposed ID retains catalog provenance");
        ok &= check(minor.chunk.candidates[2].phase_set == retained_deletion.phase_set &&
                    minor.chunk.candidates[2].hap_to_cons_alle == retained_deletion.hap_to_cons_alle &&
                    !minor.site_meta[2].bam_low_fraction_snp &&
                    minor.chunk.candidates[1].hap_to_cons_alle == minor_neighbor.hap_to_cons_alle &&
                    minor.chunk.haps == minor_haps,
                    "physical minor SNP: other ALT rows and surviving read gauges stay intact");
        GraphChunkBuildResult singleton;
        singleton.site_ids.push_back(graph_site_key_str(physical_catalog.sites[0]));
        singleton.chunk.candidates.push_back(original_categories[0]);
        singleton.chunk.candidates[0].phase_set = 100;
        singleton.chunk.reads.resize(1);
        singleton.chunk.haps = {1};
        singleton.chunk.phase_sets = {100};
        reclassify_physically_validated_graph_snps(physical_catalog.view_all(), singleton);
        ok &= check(singleton.chunk.haps[0] == 0 &&
                    singleton.chunk.phase_sets[0] == kUnphasedReadPhaseSet,
                    "physical homozygote: singleton read labels require reassignment");
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

        ok &= check(tb.site_meta.size() == 2 && tb.site_meta[0].non_selected_alt_class &&
                    tb.site_meta[1].non_selected_alt_class,
                    "multiallelic decomp: collapsed other-ALT classes remain explicit");
        ok &= check(tb.chunk.candidates.size() == 2,
                    "multiallelic decomp: triallelic → 2 biallelic pairs");
        ok &= check(tb.site_ids.size() == 2 &&
                    tb.site_ids[0] == "tri_both:1" && tb.site_ids[1] == "tri_both:2",
                    "multiallelic decomp: pair IDs carry original alt index");
        // Both split rows retain the full VCF ALT list. Sequence matching must
        // select the allele carried by that row, including when only ALT2 survives.
        tb.site_meta[0].ref = tb.site_meta[1].ref = "A";
        tb.site_meta[0].alts = tb.site_meta[1].alts = {"C", "G"};
        const std::string* selected_first = selected_graph_candidate_alt(tb, 0);
        const std::string* selected_second = selected_graph_candidate_alt(tb, 1);
        ok &= check(selected_first != nullptr && *selected_first == "C" &&
                    selected_second != nullptr && *selected_second == "G",
                    "multiallelic decomp: each pair matches only its selected ALT");
        GraphChunkBuildResult alt2_only;
        alt2_only.chunk.candidates.push_back(tb.chunk.candidates[1]);
        alt2_only.site_meta.push_back(tb.site_meta[1]);
        alt2_only.site_allele_orig_idx.push_back(tb.site_allele_orig_idx[1]);
        const std::string* selected_only = selected_graph_candidate_alt(alt2_only, 0);
        ok &= check(selected_only != nullptr && *selected_only == "G",
                    "multiallelic decomp: ALT1 cannot claim an ALT2-only candidate");
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

    // In the default REF-vs-ALT projection, a read on another ALT cannot
    // testify that it carries the literal reference allele. This occurs at
    // chr20:21,742,440: CT->CTT reads are unknown on the CT->C deletion row.
    {
        GraphSite multi;
        multi.chrom = multi.ref_contig = "chr1";
        multi.pos = multi.ref_beg = 100;
        multi.ref_end = 101;
        multi.id = "multi_default";
        multi.allele_traversals = {">1>2>3", ">1>4>3", ">1>5>3"};
        multi.allele_walks = {parse_graph_walk(">1>2>3"),
                              parse_graph_walk(">1>4>3"),
                              parse_graph_walk(">1>5>3")};
        multi.skip_reason = graph_site_validation_skip_reason(multi);
        multi.eligible = multi.skip_reason.empty();
        GraphSiteCatalog catalog;
        catalog.sites.push_back(multi);
        std::vector<GraphReadAllele> observations;
        for (int i = 0; i < 6; ++i)
            observations.push_back({"multi_default", "chr1", 100,
                                    "ref_" + std::to_string(i), 0});
        for (int i = 0; i < 5; ++i) {
            observations.push_back({"multi_default", "chr1", 100,
                                    "alt1_" + std::to_string(i), 1});
            observations.push_back({"multi_default", "chr1", 100,
                                    "alt2_" + std::to_string(i), 2});
        }
        Options opts;
        auto result = build_graph_chunk(catalog.view_all(), observations,
                                        "chr1", 0, 200, 0, opts);
        const ReadVariantProfile* alt2 =
            profile_for_read(result.chunk, "alt2_0");
        ok &= check(result.chunk.candidates.size() == 2 && alt2 != nullptr &&
                    alt2->start_var_idx == 1 && alt2->alleles.size() == 1 &&
                    alt2->alleles[0] == 1,
                    "default multiallelic projection: other ALT is unknown, not REF");
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
        ok &= check(selected_graph_candidate_alt(aa, 0) == nullptr,
                    "whole snarl: binary BAM allele cannot claim a multiallelic row");
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

    // Equal support in independent blocks prefers SNPs over indels, while
    // equal support from the same evidence tier remains ambiguous.
    {
        const auto tied_blocks = [](bool second_snp, bool certify = true) {
            PhasingChunk chunk;
            for (int i = 0; i < 2; ++i) {
                CandidateVariant site;
                const bool is_snp = i == 0 || second_snp;
                site.key.pos = 1000 + 100 * i;
                site.key.type = is_snp ? VariantType::Snp : VariantType::Insertion;
                site.key.ref_len = is_snp ? 1 : 0;
                site.key.alt = "A";
                site.counts.n_uniq_alles = 2;
                site.counts.category = is_snp ? VariantCategory::CleanHetSnp
                                              : VariantCategory::CleanHetIndel;
                site.lcd_var_i_to_cate = is_snp ? kCandCleanHetSnp : kCandCleanHetIndel;
                site.phase_set = 900 + 100 * i;
                site.hap_to_cons_alle = {-1, 0, 1};
                chunk.candidates.push_back(site);
            }
            chunk.reads.emplace_back();
            ReadVariantProfile profile;
            profile.read_id = 0;
            profile.start_var_idx = 0;
            profile.end_var_idx = 1;
            profile.alleles = {0, 1};
            profile.alt_qi = {-1, -1};
            chunk.read_var_profile.push_back(profile);
            chunk.haps = {0};
            chunk.phase_sets = {kUnphasedReadPhaseSet};
            if (certify) {
                for (int site_i = 0; site_i < 2; ++site_i) {
                    constexpr int kPrimarySupport = 32;
                    for (int i = 0; i < kPrimarySupport; ++i) {
                        const int allele = i % 2;
                        chunk.reads.emplace_back();
                        ReadVariantProfile primary;
                        primary.read_id = static_cast<int>(chunk.reads.size()) - 1;
                        primary.start_var_idx = primary.end_var_idx = site_i;
                        primary.alleles = {allele};
                        primary.alt_qi = {-1};
                        chunk.read_var_profile.push_back(primary);
                        chunk.haps.push_back(allele + 1);
                        chunk.phase_sets.push_back(chunk.candidates[site_i].phase_set);
                    }
                }
            }
            return chunk;
        };
        auto mixed = tied_blocks(false);
        ok &= check(rescue_unphased_graph_reads(mixed) == 1 &&
                    mixed.gap_haps[0] == 1 &&
                    mixed.gap_phase_sets[0] == 900 + kGapFillPsOffset,
                    "read rescue: equal SNP and indel support prefers the SNP block");
        ok &= check(mixed.candidates[0].phase_set == 900 &&
                    mixed.candidates[1].phase_set == 1000,
                    "read rescue: evidence tier preference never joins blocks");
        auto snps = tied_blocks(true);
        ok &= check(rescue_unphased_graph_reads(snps) == 0,
                    "read rescue: equal SNP blocks remain ambiguous");
        auto weak = tied_blocks(false, false);
        ok &= check(rescue_unphased_graph_reads(weak) == 0,
                    "read rescue: an uncertified SNP cannot break a block tie");
        auto inferred = tied_blocks(false);
        inferred.candidates[0].phase_set = kUnsetCandidatePhaseSet;
        inferred.candidates[0].hap_to_cons_alle = {-1, -1, -1};
        ok &= check(rescue_unphased_graph_reads(inferred) == 0,
                    "read rescue: an inferred SNP association cannot break a block tie");
        auto coincident = tied_blocks(false);
        CandidateVariant inferred_snp = coincident.candidates[0];
        inferred_snp.phase_set = kUnsetCandidatePhaseSet;
        inferred_snp.hap_to_cons_alle = {-1, -1, -1};
        coincident.candidates.push_back(inferred_snp);
        auto& direct_indel = coincident.candidates[0];
        direct_indel.key.type = VariantType::Insertion;
        direct_indel.key.ref_len = 0;
        direct_indel.counts.category = VariantCategory::CleanHetIndel;
        direct_indel.lcd_var_i_to_cate = kCandCleanHetIndel;
        auto& target = coincident.read_var_profile[0];
        target.end_var_idx = 2;
        target.alleles.push_back(0);
        target.alt_qi.push_back(-1);
        for (size_t ri = 1; ri < coincident.read_var_profile.size(); ++ri) {
            auto& primary = coincident.read_var_profile[ri];
            if (primary.start_var_idx != 0) continue;
            const int allele = primary.alleles[0];
            primary.end_var_idx = 2;
            primary.alleles = {allele, -1, allele};
            primary.alt_qi = {-1, -1, -1};
        }
        ok &= check(rescue_unphased_graph_reads(coincident) == 0,
                    "read rescue: co-located indel and inferred SNP cannot combine certificates");
        auto indels = tied_blocks(false);
        indels.candidates[0] = indels.candidates[1];
        indels.candidates[0].key.pos = 1000;
        indels.candidates[0].phase_set = 900;
        ok &= check(rescue_unphased_graph_reads(indels) == 0,
                    "read rescue: equal indel blocks remain ambiguous");
    }

    // One inferred site may tag a read only when primary graph assignments
    // establish a low-discordance allele/haplotype relation. Thirty-one
    // perfectly separated primary reads pass the Wilson bound; seven do not.
    {
        const auto singleton_rescue = [](int primary_reads) {
            PhasingChunk rescue;
            CandidateVariant excluded;
            excluded.key.pos = 1200;
            excluded.key.type = VariantType::Insertion;
            excluded.counts.n_uniq_alles = 2;
            excluded.counts.category = VariantCategory::RepeatHetIndel;
            excluded.lcd_var_i_to_cate = kLongcalldRepHetVar;
            rescue.candidates.push_back(std::move(excluded));

            constexpr hts_pos_t kPhaseSet = 950;
            for (int read_i = 0; read_i <= primary_reads; ++read_i) {
                const bool unphased = read_i == primary_reads;
                const bool hap1 = read_i < primary_reads / 2;
                ReadRecord read;
                read.qname = "singleton_" + std::to_string(read_i);
                rescue.reads.push_back(std::move(read));

                ReadVariantProfile profile;
                profile.read_id = read_i;
                profile.start_var_idx = 0;
                profile.end_var_idx = 0;
                profile.alleles = {hap1 ? 0 : 1};
                rescue.read_var_profile.push_back(std::move(profile));
                rescue.haps.push_back(unphased ? 0 : (hap1 ? 1 : 2));
                rescue.phase_sets.push_back(
                    unphased ? kUnphasedReadPhaseSet : kPhaseSet);
            }
            return std::make_pair(
                rescue_unphased_graph_reads(rescue), std::move(rescue));
        };

        auto [strong_count, strong] = singleton_rescue(31);
        ok &= check(
            strong_count == 1 && strong.gap_haps.back() == 2 &&
                strong.gap_phase_sets.back() == 950 + kGapFillPsOffset,
            "graph read rescue: strong inferred singleton tags read");

        auto [weak_count, weak] = singleton_rescue(7);
        ok &= check(
            weak_count == 0 && weak.gap_haps.back() == 0,
            "graph read rescue: weak inferred singleton abstains");
    }

    // A read can appear in two adjacent chunks while only the upstream chunk
    // has informative alleles. The downstream visit must not erase that valid
    // HP/PS assignment merely because it is unphased. A later phased visit
    // still owns the read, preserving ordinary downstream ownership.
    {
        chunks[0].chunk.haps.assign(chunks[0].chunk.reads.size(), 1);
        chunks[0].chunk.phase_sets.assign(chunks[0].chunk.reads.size(), 100);
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

    // A stale heterozygous internal consensus does not turn a homozygous
    // category into a direct read marker. Independent primary association
    // remains available through the ordinary inferred-marker path.
    {
        const auto singleton = [](VariantCategory category) {
            PhasingChunk chunk;
            CandidateVariant site;
            site.key.type = VariantType::Insertion;
            site.key.pos = 100;
            site.counts.category = category;
            site.phase_set = 100;
            site.hap_to_cons_alle = {-1, 0, 1};
            chunk.candidates.push_back(site);
            chunk.reads.emplace_back();
            chunk.haps = {0};
            chunk.phase_sets = {0};
            ReadVariantProfile profile;
            profile.read_id = 0;
            profile.start_var_idx = profile.end_var_idx = 0;
            profile.alleles = {1};
            chunk.read_var_profile.push_back(profile);
            return chunk;
        };
        PhasingChunk noisy_hom = singleton(VariantCategory::NoisyCandHom);
        ok &= check(rescue_unphased_graph_reads(noisy_hom) == 0,
                    "hom singleton: stale noisy het orientation is not direct evidence");
        PhasingChunk clean_hom = singleton(VariantCategory::CleanHom);
        ok &= check(rescue_unphased_graph_reads(clean_hom) == 0,
                    "hom singleton: stale clean het orientation is not direct evidence");
        PhasingChunk het = singleton(VariantCategory::NoisyCandHet);
        ok &= check(rescue_unphased_graph_reads(het) == 1 && het.gap_haps[0] == 2,
                    "hom singleton: ordinary phased heterozygote remains usable");
    }

    // A BAM-observation rescue is additive across overlapping chunks. Once an
    // earlier graph-supported rescue assigned the read, a later fill-only
    // assignment cannot replace its HP/PS gauge.
    {
        GraphChunkBuildResult fill_chunk;
        fill_chunk.chunk.region.chunk_id = 10;
        ReadRecord read;
        read.qname = "fill_only";
        fill_chunk.chunk.reads.push_back(std::move(read));
        ReadVariantProfile profile;
        profile.read_id = 0;
        fill_chunk.chunk.read_var_profile.push_back(std::move(profile));
        fill_chunk.chunk.haps = {0};
        fill_chunk.chunk.phase_sets = {kUnphasedReadPhaseSet};
        fill_chunk.chunk.gap_haps = {1};
        fill_chunk.chunk.gap_phase_sets = {1100 + kGapFillPsOffset};
        fill_chunk.chunk.gap_from_bam_observation = {false};

        std::unordered_map<std::string, PhaseReadOutputRow> rows;
        merge_graph_chunk_into_read_rows(rows, fill_chunk, 0);
        fill_chunk.chunk.region.chunk_id = 11;
        fill_chunk.chunk.gap_haps[0] = 2;
        fill_chunk.chunk.gap_phase_sets[0] = 1200 + kGapFillPsOffset;
        fill_chunk.chunk.gap_from_bam_observation[0] = true;
        merge_graph_chunk_into_read_rows(rows, fill_chunk, 0);

        ok &= check(
            rows.at("fill_only").hap == 1 &&
                rows.at("fill_only").phase_set ==
                    1100 + kGapFillPsOffset,
            "graph output merge: BAM-observation rescue is fill-only");
    }

    // A graph MNP has ref_len > 1 even though its internal type is SNP.
    // The VCF must retain the whole replaced reference sequence.
    {
        char fasta_path[] = "/tmp/pgphase_mnp_writer_XXXXXX";
        const int fd = mkstemp(fasta_path);
        ok &= check(fd >= 0, "MNP writer: create temporary reference");
        if (fd >= 0) {
            close(fd);
            {
                std::ofstream fasta(fasta_path);
                fasta << ">chr1\nATGCCATGCCAT\n";
            }
            ok &= check(fai_build(fasta_path) == 0,
                        "MNP writer: index temporary reference");
            std::unique_ptr<faidx_t, decltype(&fai_destroy)> fai(
                fai_load(fasta_path), fai_destroy);
            const char header_text[] =
                "@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:12\n";
            std::unique_ptr<bam_hdr_t, decltype(&sam_hdr_destroy)> header(
                sam_hdr_parse(std::strlen(header_text), header_text),
                sam_hdr_destroy);
            ok &= check(fai != nullptr && header != nullptr,
                        "MNP writer: load reference and header");
            if (fai != nullptr && header != nullptr) {
                CandidateVariant candidate;
                candidate.key.tid = 0;
                candidate.key.type = VariantType::Snp;
                candidate.key.pos = 2;
                candidate.key.ref_len = 2;
                candidate.key.alt = "CA";
                candidate.counts.category = VariantCategory::CleanHetSnp;
                candidate.counts.ref_cov = 5;
                candidate.counts.alt_cov = 5;
                candidate.counts.total_cov = 10;
                candidate.counts.allele_fraction = 0.5;
                candidate.hap_to_cons_alle = {-1, 1, 0};
                candidate.phase_set = 2;
                ReferenceCache ref(fai.get());
                std::ostringstream vcf;
                write_phased_variants_vcf_records(
                    vcf, Options{}, header.get(), ref, {candidate});
                ok &= check(vcf.str().find("\t2\t.\tTG\tCA\t") !=
                                std::string::npos,
                            "MNP writer: complete REF describes the graph allele");
                Options replacement_opts;
                replacement_opts.min_sv_len = 1;
                candidate.key.type = VariantType::Deletion;
                candidate.key.pos = 3;
                candidate.key.ref_len = 2;
                candidate.key.alt = "A";
                candidate.counts.category = VariantCategory::CleanHetIndel;
                std::ostringstream deletion_vcf;
                write_phased_variants_vcf_records(
                    deletion_vcf, replacement_opts, header.get(), ref, {candidate});
                ok &= check(deletion_vcf.str().find("\t2\t.\tTGC\tTA\t") !=
                                std::string::npos &&
                            deletion_vcf.str().find("SVLEN=-1") !=
                                std::string::npos,
                            "deletion writer: replacement ALT and net length survive");

                candidate.key.type = VariantType::Insertion;
                candidate.key.ref_len = 1;
                candidate.key.alt = "AA";
                std::ostringstream insertion_vcf;
                write_phased_variants_vcf_records(
                    insertion_vcf, replacement_opts, header.get(), ref, {candidate});
                ok &= check(insertion_vcf.str().find("\t2\t.\tTG\tTAA\t") !=
                                std::string::npos &&
                            insertion_vcf.str().find("SVLEN=1") !=
                                std::string::npos,
                            "insertion writer: replacement net length is reported");
            }
            std::remove(fasta_path);
            std::remove((std::string(fasta_path) + ".fai").c_str());
        }
    }

    if (ok) {
        std::cout << "ALL PASS\n";
        return 0;
    }
    return 1;
}
