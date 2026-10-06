#include "graph_bam_adapter.hpp"

#include "collect_phase.hpp"
#include "collect_output.hpp"

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

int main() {
    bool ok = true;

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

    {
        CandidateVariant graph;
        graph.graph_site = true;
        graph.key.pos = 90;
        graph.key.alt = ">10>20>30";
        graph.counts.n_uniq_alles = 2;
        graph.counts.category = VariantCategory::RepeatHetIndel;
        graph.lcd_var_i_to_cate = kLongcalldRepHetVar;
        CandidateVariant source;
        source.key.pos = 101;
        source.key.type = VariantType::Deletion;
        source.key.ref_len = 1;
        source.counts.category = VariantCategory::NoisyCandHet;
        source.counts.n_uniq_alles = 2;
        source.counts.ref_cov = 11;
        source.counts.alt_cov = 13;
        source.lcd_var_i_to_cate = kCandNoisyCandHet;
        source.msa_verified = true;
        source.phase_set = 95;
        source.hap_to_cons_alle = {-1, 1, 0};
        CandidateVariant left = source;
        left.key.pos = 91;
        CandidateVariant right = source;
        right.key.pos = 111;
        std::vector<CandidateVariant> sites{right, source, left};
        ok &= check(bam_site_has_only_weak_links(sites, 1, {90, 100}),
                    "isolated BAM allele: both weak edges permit independent ownership");
        ok &= check(!bam_site_has_only_weak_links(sites, 1, {90}) &&
                    !bam_site_has_only_weak_links(sites, 1, {100}),
                    "isolated BAM allele: either supported edge retains source component");
        sites[0].phase_set = 200;
        ok &= check(bam_site_has_only_weak_links(sites, 1, {90}),
                    "isolated BAM allele: another PS does not define a source edge");
        sites[0] = source;
        sites[0].key.ref_len = 2;
        ok &= check(!bam_site_has_only_weak_links(sites, 1, {90, 100}),
                    "isolated BAM allele: co-located contrasts remain independent rows");
        ok &= check(!bam_site_has_only_weak_links({source}, 0, {90}),
                    "isolated BAM allele: a singleton source has no internal weak edge");
        sites[0] = right;
        sites[1].phase_set = kUnsetCandidatePhaseSet;
        ok &= check(!bam_site_has_only_weak_links(sites, 1, {90, 100}),
                    "isolated BAM allele: an unphased row supplies no local genotype");

        const CandidateVariant original_graph = graph;
        ok &= check(adopt_unphased_graph_allele_from_bam(graph, source, 96),
                    "shared MSA: unphased catalog allele retains verified BAM genotype");
        ok &= check(graph.key.pos == 101 && graph.key.alt.empty() &&
                    graph.key.type == VariantType::Deletion && graph.key.ref_len == 1 &&
                    graph.counts.category == VariantCategory::NoisyCandHet &&
                    graph.counts.alle_covs == std::vector<int>{11, 13} &&
                    graph.hap_to_cons_alle == std::array<int, 3>{-1, 1, 0} &&
                    graph.phase_set == 96 && graph.msa_verified &&
                    graph.alignment_verified && graph.bam_injected && graph.graph_site,
                    "shared MSA: key, counts, phase gauge and validation stay together");
        ok &= check(!adopt_unphased_graph_allele_from_bam(graph, source, 97) &&
                    graph.phase_set == 96,
                    "shared MSA: cannot overwrite an already phased catalog row");
        graph = original_graph;
        source.msa_verified = false;
        ok &= check(!adopt_unphased_graph_allele_from_bam(graph, source, 96),
                    "shared MSA: an unverified noisy candidate cannot own the allele");
        source.msa_verified = true;
        source.hap_to_cons_alle = {-1, 1, 1};
        ok &= check(!adopt_unphased_graph_allele_from_bam(graph, source, 96),
                    "shared MSA: homozygotes cannot supply an oriented genotype");
        source.hap_to_cons_alle = {-1, 1, 2};
        source.counts.n_uniq_alles = 3;
        ok &= check(!adopt_unphased_graph_allele_from_bam(graph, source, 96),
                    "shared MSA: a multiallelic source cannot replace a binary row");
        source.hap_to_cons_alle = {-1, 1, 0};
        source.counts.n_uniq_alles = 2;
        graph.counts.n_uniq_alles = 3;
        ok &= check(!adopt_unphased_graph_allele_from_bam(graph, source, 96),
                    "shared MSA: binary source cannot claim a whole snarl");
        graph = original_graph;
        ok &= check(!adopt_unphased_graph_allele_from_bam(graph, source, 0),
                    "shared MSA: mapped source label must be phased");
        graph.counts.category = VariantCategory::CleanHetIndel;
        ok &= check(!adopt_unphased_graph_allele_from_bam(graph, source, 96),
                    "shared MSA: only a demoted repeat can change ownership");
    }


    {
        // Exercise the production BAM query, graph projection and diploid
        // count/quality checks together, without comparator or truth inputs.
        char fasta_path[] = "/tmp/pgphase_retry_snp_XXXXXX";
        const int fd = mkstemp(fasta_path);
        ok &= check(fd >= 0, "SNP retry: create synthetic fixture");
        if (fd >= 0) {
            close(fd);
            const std::string bam_path = std::string(fasta_path) + ".bam";
            {
                std::ofstream fasta(fasta_path);
                // Soft masking must not remove physical REF calls.
                fasta << ">chr1\n" << std::string(300, 'a') << '\n';
            }
            ok &= check(fai_build(fasta_path) == 0,
                        "SNP retry: index synthetic reference");
            const auto make_graph = [] {
                GraphChunkBuildResult graph;
                for (const hts_pos_t pos : {100, 200}) {
                    CandidateVariant site;
                    site.key.tid = 0;
                    site.key.type = VariantType::Snp;
                    site.key.pos = pos - 10;
                    site.key.ref_len = 2;
                    site.key.alt = ">10>20>30";
                    site.counts.n_uniq_alles = 2;
                    site.counts.category = VariantCategory::CleanHetSnp;
                    site.phase_set = pos - 10;
                    site.hap_to_cons_alle = {-1, 0, 1};
                    graph.chunk.candidates.push_back(site);
                    graph.site_meta.push_back(
                        {"chr1", pos - 1, "AA", {pos == 100 ? "AC" : "AG"}});
                    graph.site_allele_orig_idx.push_back({0, 1});
                }
                return graph;
            };
            const auto graph = make_graph();
            const RecoverySeam seam{100, 200, 90, 190};
            // counts: REF/REF, REF/ALT, ALT/REF, ALT/ALT.
            const auto query = [&](const GraphChunkBuildResult& input,
                                   const std::array<int, 4>& counts,
                                   int mapq = 60, int flags = 0,
                                   bool repeat_names = false,
                                   int left_baseq = 40, int right_baseq = 40,
                                   std::optional<bool>* graph_parity = nullptr) {
                const char header_text[] =
                    "@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:300\n";
                std::unique_ptr<bam_hdr_t, HeaderDeleter> header(
                    sam_hdr_parse(std::strlen(header_text), header_text));
                std::unique_ptr<samFile, decltype(&hts_close)> output(
                    sam_open(bam_path.c_str(), "wb"), hts_close);
                std::unique_ptr<bam1_t, AlignmentDeleter> alignment(bam_init1());
                if (!header || !output || !alignment ||
                    sam_hdr_write(output.get(), header.get()) < 0)
                    throw std::runtime_error("SNP retry fixture: write BAM header");
                for (size_t cell = 0; cell < counts.size(); ++cell) {
                    for (int i = 0; i < counts[cell]; ++i) {
                        std::string bases(140, 'A');
                        bases[19] = cell >= 2 ? 'C' : 'A';
                        bases[119] = cell % 2 != 0 ? 'G' : 'A';
                        std::string qualities(140, 'I');
                        // 255 is the BAM missing-quality sentinel, not a high Q.
                        qualities[19] = static_cast<char>(33 +
                            (left_baseq == 255 ? 40 : left_baseq));
                        qualities[119] = static_cast<char>(33 +
                            (right_baseq == 255 ? 40 : right_baseq));
                        std::string sam = "read" + std::to_string(cell) + "_" +
                            std::to_string(repeat_names ? 0 : i) + "\t" +
                            std::to_string(flags) + "\tchr1\t81\t" +
                            std::to_string(mapq) + "\t140M\t*\t0\t0\t" +
                            bases + "\t" + qualities;
                        kstring_t line{sam.size(), sam.size() + 1, sam.data()};
                        if (sam_parse1(&line, header.get(), alignment.get()) < 0)
                            throw std::runtime_error("SNP retry fixture: parse read");
                        if (left_baseq == 255) bam_get_qual(alignment.get())[19] = 255;
                        if (right_baseq == 255) bam_get_qual(alignment.get())[119] = 255;
                        if (sam_write1(output.get(), header.get(), alignment.get()) < 0)
                            throw std::runtime_error("SNP retry fixture: write read");
                    }
                }
                output.reset();
                if (sam_index_build(bam_path.c_str(), 0) < 0)
                    throw std::runtime_error("SNP retry fixture: index BAM");
                Options opts;
                opts.ref_fasta = fasta_path;
                opts.bam_files = {bam_path};
                WorkerContext context(opts);
                return has_direct_snp_parity_for_retry(input, seam, context, 0, graph_parity);
            };
            ok &= check(query(graph, {8, 0, 0, 8}),
                        "SNP retry: padded graph walks use selected physical bases");
            ok &= check(query(graph, {0, 8, 8, 0}),
                        "SNP retry: alternate connection remains eligible");
            ok &= check(!query(graph, {2, 0, 0, 2}),
                        "SNP retry: Q40 does not certify four unanimous molecules");
            ok &= check(!query(graph, {0, 0, 0, 28}),
                        "SNP retry: one haplotype does not certify diploid parity");
            ok &= check(!query(graph, {8, 8, 8, 8}),
                        "SNP retry: tied parity stays ineligible");
            ok &= check(query(graph, {15, 1, 1, 15}),
                        "SNP retry: decisive nonunanimous calls remain eligible");
            ok &= check(!query(graph, {3, 1, 1, 3}),
                        "SNP retry: weak 3:1 count evidence stays ineligible");
            ok &= check(!query(graph, {8, 0, 0, 8}, 29),
                        "SNP retry: low MAPQ is not parity evidence");
            ok &= check(!query(graph, {8, 0, 0, 8}, 255),
                        "SNP retry: missing MAPQ is not parity evidence");
            ok &= check(!query(graph, {8, 0, 0, 8}, 60, 0, false, 29),
                        "SNP retry: low left base quality is excluded");
            ok &= check(!query(graph, {8, 0, 0, 8}, 60, 0, false, 40, 29),
                        "SNP retry: low right base quality is excluded");
            ok &= check(!query(graph, {8, 0, 0, 8}, 60, 0, false, 255),
                        "SNP retry: missing left base quality is excluded");
            ok &= check(!query(graph, {8, 0, 0, 8}, 60, 0, false, 40, 255),
                        "SNP retry: missing right base quality is excluded");
            for (const int flag : {BAM_FSECONDARY, BAM_FSUPPLEMENTARY,
                                   BAM_FDUP, BAM_FQCFAIL})
                ok &= check(!query(graph, {8, 0, 0, 8}, 60, flag),
                            "SNP retry: nonprimary/failed molecules are excluded");
            ok &= check(!query(graph, {8, 0, 0, 8}, 60, 0, true),
                        "SNP retry: repeated read names cannot multiply support");
            std::optional<bool> graph_parity = true;
            ok &= check(query(graph, {8, 0, 0, 8}, 60, 0, false, 40, 40,
                              &graph_parity) && graph_parity && !*graph_parity,
                        "SNP retry: same graph relation reaches the stitch constraint");
            ok &= check(query(graph, {0, 8, 8, 0}, 60, 0, false, 40, 40,
                              &graph_parity) && graph_parity && *graph_parity,
                        "SNP retry: swapped graph relation reaches the stitch constraint");
            ok &= check(!query(graph, {2, 0, 0, 2}, 60, 0, false, 40, 40,
                               &graph_parity) && !graph_parity,
                        "SNP retry: failed admission clears an earlier relation");
            auto modified = make_graph();
            modified.site_meta[0].alts = {"AT", "AC"};
            modified.site_allele_orig_idx[0] = {0, 2};
            ok &= check(query(modified, {8, 0, 0, 8}),
                        "SNP retry: decomposed graph row selects original ALT2");
            modified.site_allele_orig_idx[0] = {0, 1};
            ok &= check(!query(modified, {8, 0, 0, 8}),
                        "SNP retry: a different ALT is not usable evidence");
            modified = make_graph();
            modified.site_meta[0].ref = "AAA";
            modified.site_meta[0].alts = {"ACC"};
            ok &= check(!query(modified, {8, 0, 0, 8}),
                        "SNP retry: graph MNP is not a single-base SNP");
            modified = make_graph();
            modified.site_meta[0].alts = {"AAC"};
            ok &= check(!query(modified, {8, 0, 0, 8}),
                        "SNP retry: graph insertion is not a single-base SNP");
            modified = make_graph();
            modified.site_allele_orig_idx[0].clear();
            ok &= check(!query(modified, {8, 0, 0, 8}),
                        "SNP retry: missing graph allele mapping cannot certify parity");
            modified = make_graph();
            modified.chunk.candidates[0].counts.n_uniq_alles = 3;
            ok &= check(!query(modified, {8, 0, 0, 8}),
                        "SNP retry: whole multiallelic graph row is excluded");
            modified = make_graph();
            modified.chunk.candidates[0].hap_to_cons_alle = {-1, 0, 0};
            ok &= check(!query(modified, {8, 0, 0, 8}),
                        "SNP retry: homozygote cannot orient a boundary");
            modified = make_graph();
            modified.chunk.candidates[0].counts.category = VariantCategory::NoisyCandHet;
            ok &= check(!query(modified, {8, 0, 0, 8}),
                        "SNP retry: unverified noisy SNP cannot orient a boundary");
            modified.chunk.candidates[0].msa_verified = true;
            modified.chunk.candidates[0].alignment_verified = true;
            ok &= check(query(modified, {8, 0, 0, 8}),
                        "SNP retry: alignment/MSA-verified SNP remains eligible");
            modified = make_graph();
            for (size_t i = 0; i < 2; ++i) {
                auto& site = modified.chunk.candidates[i];
                site.bam_injected = true;
                site.key.pos = i == 0 ? 100 : 200;
                site.key.ref_len = 1;
                site.key.alt = i == 0 ? "C" : "G";
            }
            modified.site_meta.clear();
            modified.site_allele_orig_idx.clear();
            ok &= check(query(modified, {8, 0, 0, 8}),
                        "SNP retry: injected BAM SNPs retain exact sequence keys");
            ok &= check(query(modified, {2, 0, 0, 2}),
                        "SNP retry: existing BAM-only four-read admission is preserved");
            ok &= check(query(modified, {0, 0, 0, 28}),
                        "SNP retry: existing BAM-only source evidence is preserved");
            ok &= check(query(modified, {3, 1, 1, 3}),
                        "SNP retry: existing BAM-only dominant relation is preserved");
            graph_parity = true;
            ok &= check(query(modified, {8, 0, 0, 8}, 60, 0, false, 40, 40,
                              &graph_parity) && !graph_parity,
                        "SNP retry: BAM-only admission does not create a graph constraint");
            modified.chunk.candidates[0] = graph.chunk.candidates[0];
            modified.site_meta = graph.site_meta;
            modified.site_allele_orig_idx = graph.site_allele_orig_idx;
            ok &= check(!query(modified, {2, 0, 0, 2}),
                        "SNP retry: mixed graph/BAM four-read vote needs count confidence");
            ok &= check(!query(modified, {0, 0, 0, 28}),
                        "SNP retry: mixed graph/BAM one-haplotype cohort stays ineligible");
            ok &= check(!query(modified, {3, 1, 1, 3}),
                        "SNP retry: mixed graph/BAM 3:1 split needs count confidence");
            ok &= check(!query(modified, {20, 0, 2, 0}),
                        "SNP retry: one haplotype's reversal cannot supply diploid support");
            ok &= check(!query(modified, {20, 0, 5, 2}),
                        "SNP retry: aggregate majority cannot override a haplotype reversal");
            std::remove(bam_path.c_str());
            std::remove((bam_path + ".bai").c_str());
            std::remove((std::string(fasta_path) + ".fai").c_str());
            std::remove(fasta_path);
        }
    }

    {
        const auto make_boundary = [] {
            PhasingChunk chunk;
            for (const hts_pos_t pos : {100, 200}) {
                CandidateVariant site;
                site.key.type = VariantType::Snp;
                site.key.pos = pos;
                site.key.ref_len = 1;
                site.phase_set = pos;
                site.hap_to_cons_alle = {-1, 0, 1};
                chunk.candidates.push_back(site);
            }
            return chunk;
        };
        const std::vector<RecoverySeam> seams = {{100, 200, 100, 200}};
        const auto make_gauge = [](bool flip) {
            RecoveryPhaseGauge gauge;
            gauge.beg = 90;
            gauge.end = 210;
            gauge.physical_snp_bridges.push_back({100, 200, flip});
            return gauge;
        };
        Options opts;
        auto boundary = make_boundary();
        auto gauge = make_gauge(false);
        ok &= check(stitch_recovery_phase_sets_left_to_right(
                        boundary, seams, {gauge}, opts) == 1 &&
                    boundary.candidates[1].phase_set == 100 &&
                    boundary.candidates[1].hap_to_cons_alle[1] == 0,
                    "physical constraint: same parity joins without reversing alleles");
        boundary = make_boundary();
        gauge = make_gauge(true);
        ok &= check(stitch_recovery_phase_sets_left_to_right(
                        boundary, seams, {gauge}, opts) == 1 &&
                    boundary.candidates[1].phase_set == 100 &&
                    boundary.candidates[1].hap_to_cons_alle[1] == 1,
                    "physical constraint: alternate parity applies exactly one flip");
        for (const bool first_flip : {false, true}) {
            boundary = make_boundary();
            gauge = make_gauge(first_flip);
            gauge.physical_snp_bridges.push_back({100, 200, !first_flip});
            ok &= check(stitch_recovery_phase_sets_left_to_right(
                            boundary, seams, {gauge}, opts) == 0 &&
                        boundary.candidates[1].phase_set == 200 &&
                        boundary.candidates[1].hap_to_cons_alle[1] == 0,
                        "physical constraint: contradictory relations veto either row order");
        }
        boundary = make_boundary();
        gauge = make_gauge(false);
        gauge.physical_snp_bridges.push_back({100, 200, false});
        ok &= check(stitch_recovery_phase_sets_left_to_right(
                        boundary, seams, {gauge}, opts) == 1,
                    "physical constraint: agreeing duplicate evidence remains usable");
    }

    {
        PhasingChunk pair;
        pair.candidates.resize(2);
        for (CandidateVariant& site : pair.candidates) {
            site.bam_injected = true;
            site.hap_to_cons_alle = {-1, 0, 1};
        }
        constexpr int kVoteCount = 32;
        constexpr int kMinMapq = 30;
        constexpr double kMaxP = 0.001;
        for (int ri = 0; ri < kVoteCount; ++ri) {
            ReadRecord read;
            read.qname = "paired-molecule-" + std::to_string(ri);
            read.mapq = 10;
            pair.reads.push_back(std::move(read));
            ReadVariantProfile profile;
            profile.start_var_idx = 0;
            profile.end_var_idx = 1;
            profile.alleles = {ri % 2, ri % 2};
            profile.bam_alleles = profile.alleles;
            profile.bam_mapq = 60;
            pair.read_var_profile.push_back(std::move(profile));
        }
        const auto supported = [&] {
            return local_run_boundary_flip(pair, 0, 1, kMinMapq, kMaxP);
        };
        const auto same = supported();
        ok &= check(same && !*same,
                    "BAM pair: low GAF MAPQ cannot discard mapped BAM votes");
        for (ReadVariantProfile& profile : pair.read_var_profile)
            profile.bam_alleles[1] = 1 - profile.bam_alleles[0];
        const auto cross = supported();
        ok &= check(cross && *cross,
                    "BAM pair: parity uses BAM calls rather than working calls");
        for (ReadRecord& read : pair.reads) read.mapq = 60;
        for (const int mapq : {10, 255}) {
            for (ReadVariantProfile& profile : pair.read_var_profile)
                profile.bam_mapq = mapq;
            ok &= check(!supported(),
                        "BAM pair: GAF MAPQ cannot certify low or unknown BAM MAPQ");
        }
        for (ReadVariantProfile& profile : pair.read_var_profile) {
            profile.bam_mapq = 60;
            profile.bam_alleles = {-1, -1};
        }
        ok &= check(!supported(),
                    "BAM pair: missing source calls cannot borrow working calls");
        pair.candidates[1].bam_injected = false;
        const auto mixed = supported();
        ok &= check(mixed && !*mixed,
                    "mixed pair: retain the established working-matrix parity");
        for (ReadRecord& read : pair.reads) read.mapq = 10;
        ok &= check(!supported(),
                    "mixed pair: BAM MAPQ cannot certify the graph observations");
        pair.candidates[1].bam_injected = true;
        for (ReadVariantProfile& profile : pair.read_var_profile)
            profile.bam_alleles = {0, 0};
        ok &= check(!supported(),
                    "BAM pair: one allele class cannot orient both haplotypes");
    }

    {
        // The complementary pair's downstream certificate must not strand
        // its cut-free BAM prefix. This proof cannot borrow a graph gauge,
        // ambiguous deletion REF calls, or a favorable earlier SNP.
        const auto make_prefix = [] {
            GraphChunkBuildResult result;
            for (size_t ci = 0; ci < 4; ++ci) {
                CandidateVariant site;
                site.key.pos = ci == 0 ? 121 : ci == 1 ? 150 : 201;
                site.key.type = ci == 0 ? VariantType::Insertion :
                    ci == 1 ? VariantType::Snp : VariantType::Deletion;
                site.key.ref_len = ci == 0 ? 0 : ci == 3 ? 2 : 1;
                site.key.alt = ci < 2 ? "C" : "";
                site.counts.category = ci == 1 ? VariantCategory::CleanHetSnp :
                    VariantCategory::CleanHetIndel;
                site.phase_set = 80;
                site.hap_to_cons_alle = {-1, ci == 2 ? 0 : 1, ci == 2 ? 1 : 0};
                site.bam_injected = site.msa_verified = site.alignment_verified = true;
                result.chunk.candidates.push_back(site);
                RecoverySourceSite origin;
                origin.candidate_index = ci;
                origin.phase_set = 40;
                origin.hap1_allele = site.hap_to_cons_alle[1];
                origin.hap2_allele = site.hap_to_cons_alle[2];
                origin.can_adopt = true;
                result.recovery_source_sites.push_back(origin);
            }
            result.recovery_source_weak_cuts[40] = {250};
            for (int ri = 0; ri < 18; ++ri) {
                ReadRecord read;
                read.qname = "prefix-molecule-" + std::to_string(ri);
                read.mapq = 5;
                result.chunk.reads.push_back(std::move(read));
                ReadVariantProfile profile;
                profile.start_var_idx = 0;
                profile.end_var_idx = 3;
                const int snp = ri < 12 ? 0 : 1;
                profile.bam_alleles = {-1, snp, 1 - snp, snp};
                profile.alleles = {-1, -1, -1, -1};
                profile.bam_mapq = 60;
                result.chunk.read_var_profile.push_back(std::move(profile));
            }
            return result;
        };
        auto prefix = make_prefix();
        const auto supported = [&] {
            return bam_prefix_before_deletion_pair(prefix, 2, 3, 110);
        };
        ok &= check(supported() == std::optional<hts_pos_t>(120),
                    "BAM prefix: retain all audited rows before the supported pair");
        for (CandidateVariant& row : prefix.chunk.candidates)
            std::swap(row.hap_to_cons_alle[1], row.hap_to_cons_alle[2]);
        ok &= check(supported() == std::optional<hts_pos_t>(120),
                    "BAM prefix: a uniform source gauge reversal remains valid");
        prefix = make_prefix();
        prefix.chunk.read_var_profile[0].bam_alleles[1] = 1;
        ok &= check(!supported(), "BAM prefix: one opposing molecule vetoes transfer");
        for (const int mapq : {10, 255}) {
            prefix = make_prefix();
            for (ReadVariantProfile& read : prefix.chunk.read_var_profile)
                read.bam_mapq = mapq;
            ok &= check(!supported(), "BAM prefix: reject low or unknown BAM MAPQ");
        }
        prefix = make_prefix();
        for (ReadVariantProfile& read : prefix.chunk.read_var_profile)
            read.bam_alleles = {-1, 0, 1, 0};
        ok &= check(!supported(), "BAM prefix: one haplotype does not certify diploidy");
        for (const int deletion_call : {0, 1}) {
            prefix = make_prefix();
            for (ReadVariantProfile& read : prefix.chunk.read_var_profile)
                read.bam_alleles[2] = read.bam_alleles[3] = deletion_call;
            ok &= check(!supported(), "BAM prefix: double REF or ALT abstains");
        }
        prefix = make_prefix();
        for (size_t ri = 0; ri < prefix.chunk.reads.size(); ++ri)
            prefix.chunk.reads[ri].qname = "duplicate-" + std::to_string(ri % 2);
        ok &= check(!supported(), "BAM prefix: duplicate names cannot inflate support");
        prefix = make_prefix();
        prefix.recovery_source_weak_cuts[40] = {150};
        ok &= check(!supported(), "BAM prefix: an internal weak cut vetoes transfer");
        prefix = make_prefix();
        prefix.recovery_source_quality_cuts[40] = {150};
        ok &= check(!supported(), "BAM prefix: an internal quality cut vetoes transfer");
        prefix = make_prefix();
        prefix.recovery_source_weak_cuts.clear();
        ok &= check(!supported(), "BAM prefix: an unaudited source is not certified");
        prefix = make_prefix();
        prefix.recovery_source_sites.erase(prefix.recovery_source_sites.begin());
        ok &= check(!supported(), "BAM prefix: missing row provenance vetoes transfer");
        prefix = make_prefix();
        prefix.recovery_source_sites[0].phase_set = 41;
        ok &= check(!supported(), "BAM prefix: another source cannot inherit the certificate");
        prefix = make_prefix();
        std::swap(prefix.chunk.candidates[0].hap_to_cons_alle[1],
                  prefix.chunk.candidates[0].hap_to_cons_alle[2]);
        ok &= check(!supported(), "BAM prefix: inconsistent original gauges veto transfer");
        prefix = make_prefix();
        prefix.chunk.candidates[0].bam_injected = false;
        ok &= check(!supported(), "BAM prefix: catalog rows require their own certificate");
        prefix = make_prefix();
        prefix.chunk.candidates[0].hap_to_cons_alle = {-1, 0, 0};
        ok &= check(supported() == std::optional<hts_pos_t>(150),
                    "BAM prefix: homozygotes do not extend the certified interval");
        prefix = make_prefix();
        prefix.chunk.candidates.push_back(prefix.chunk.candidates[2]);
        ok &= check(!supported(), "BAM prefix: a third oriented deletion row vetoes transfer");
        prefix = make_prefix();
        auto nearer = prefix.chunk.candidates[1];
        nearer.key.pos = 180;
        prefix.chunk.candidates.push_back(nearer);
        auto origin = prefix.recovery_source_sites[1];
        origin.candidate_index = 4;
        prefix.recovery_source_sites.push_back(origin);
        ok &= check(!supported(), "BAM prefix: do not bypass an uncalled nearest clean SNP");
    }

    {
        // Homozygotes inherit the BAM source's PS, but cannot orient it or
        // extend its certified heterozygous interval across a weak cut.
        constexpr hts_pos_t kRunPhaseSet = 80;
        constexpr hts_pos_t kSourcePhaseSet = 40;
        const auto make_run = [kSourcePhaseSet] {
            GraphChunkBuildResult result;
            for (const hts_pos_t pos : {100, 200}) {
                CandidateVariant site;
                site.key.type = VariantType::Snp;
                site.key.pos = pos;
                site.key.ref_len = 1;
                site.phase_set = kRunPhaseSet;
                site.hap_to_cons_alle = {-1, 0, 1};
                site.bam_injected = true;
                result.chunk.candidates.push_back(site);
                RecoverySourceSite source;
                source.candidate_index = result.chunk.candidates.size() - 1;
                source.phase_set = kSourcePhaseSet;
                source.hap1_allele = 0;
                source.hap2_allele = 1;
                source.can_adopt = true;
                result.recovery_source_sites.push_back(source);
            }
            result.recovery_source_weak_cuts[kSourcePhaseSet] = {};
            return result;
        };
        auto run = make_run();
        ok &= check(bam_source_run_supported(run, kRunPhaseSet),
                    "BAM run: exact source gauge is supported");
        CandidateVariant homozygote;
        homozygote.key.type = VariantType::Snp;
        homozygote.key.pos = 150;
        homozygote.phase_set = kRunPhaseSet;
        homozygote.hap_to_cons_alle = {-1, 0, 0};
        homozygote.bam_injected = true;
        run.chunk.candidates.push_back(homozygote);
        ok &= check(bam_source_run_supported(run, kRunPhaseSet),
                    "BAM run: source-labelled homozygote is not a gauge veto");
        run.chunk.candidates.back().key.pos = 50;
        run.chunk.candidates.back().hap_to_cons_alle = {-1, 1, 1};
        run.recovery_source_weak_cuts[kSourcePhaseSet] = {75, 200};
        ok &= check(bam_source_run_supported(run, kRunPhaseSet),
                    "BAM run: homozygote does not extend the run across a cut");
        RecoverySourceSite homozygous_source;
        homozygous_source.candidate_index = run.chunk.candidates.size() - 1;
        homozygous_source.phase_set = kSourcePhaseSet;
        homozygous_source.hap1_allele = 1;
        homozygous_source.hap2_allele = 1;
        run.recovery_source_sites.push_back(homozygous_source);
        run.recovery_source_sites.push_back(homozygous_source);
        ok &= check(bam_source_run_supported(run, kRunPhaseSet),
                    "BAM run: duplicate homozygous provenance is not a gauge veto");
        run.chunk.candidates.back().counts.category = VariantCategory::NoisyCandHom;
        run.chunk.candidates.back().hap_to_cons_alle = {-1, 0, 1};
        ok &= check(bam_source_run_supported(run, kRunPhaseSet),
                    "BAM run: stale homozygous orientation cannot extend a path cut");
        run.chunk.candidates.back().counts.category = VariantCategory::CleanHom;
        ok &= check(bam_source_run_supported(run, kRunPhaseSet),
                    "BAM run: clean homozygous orientation is not provenance");
        run.recovery_source_sites.resize(2);
        run.chunk.candidates.back().counts.category = VariantCategory::CleanHetSnp;
        run.chunk.candidates.back().hap_to_cons_alle = {-1, -1, -1};
        ok &= check(!bam_source_run_supported(run, kRunPhaseSet),
                    "BAM run: malformed source-labelled row remains a veto");
        run = make_run();
        run.chunk.candidates[0].hap_to_cons_alle = {-1, 1, 0};
        run.chunk.candidates[1].hap_to_cons_alle = {-1, 1, 0};
        ok &= check(bam_source_run_supported(run, kRunPhaseSet),
                    "BAM run: a consistent source flip is supported");
        run.chunk.candidates[1].hap_to_cons_alle = {-1, 0, 1};
        ok &= check(!bam_source_run_supported(run, kRunPhaseSet),
                    "BAM run: conflicting heterozygous gauges remain a veto");
        run = make_run();
        run.recovery_source_weak_cuts[kSourcePhaseSet] = {150};
        ok &= check(!bam_source_run_supported(run, kRunPhaseSet),
                    "BAM run: an interior weak cut remains a veto");
        run.recovery_source_weak_cuts.clear();
        ok &= check(!bam_source_run_supported(run, kRunPhaseSet),
                    "BAM run: missing path metadata is not a certificate");
        run = make_run();
        run.recovery_source_sites[1].phase_set = kSourcePhaseSet + 1;
        ok &= check(!bam_source_run_supported(run, kRunPhaseSet),
                    "BAM run: two independent sources remain separate");
        run = make_run();
        run.recovery_source_sites.push_back(run.recovery_source_sites.front());
        ok &= check(!bam_source_run_supported(run, kRunPhaseSet),
                    "BAM run: duplicate heterozygous provenance remains ambiguous");
        run = make_run();
        run.chunk.candidates[0].bam_injected = false;
        ok &= check(!bam_source_run_supported(run, kRunPhaseSet),
                    "BAM run: a graph heterozygote needs a graph path");
        // A physical insertion bridge can independently orient an attached
        // component. Its shared graph rows still need exact source lineage,
        // not the component's new numeric phase-set label.
        run.recovery_source_weak_cuts[kSourcePhaseSet] = {250};
        ok &= check(bam_source_run_supported(run, kRunPhaseSet, true),
                    "shared BAM run: downstream source cut lies outside the component");
        std::swap(run.chunk.candidates[0].hap_to_cons_alle[1],
                  run.chunk.candidates[0].hap_to_cons_alle[2]);
        ok &= check(!bam_source_run_supported(run, kRunPhaseSet, true),
                    "shared BAM run: a reversed graph row cannot borrow the source gauge");
        std::swap(run.chunk.candidates[0].hap_to_cons_alle[1],
                  run.chunk.candidates[0].hap_to_cons_alle[2]);
        run.recovery_source_sites[0].can_adopt = false;
        ok &= check(!bam_source_run_supported(run, kRunPhaseSet, true),
                    "shared BAM run: an untranslatable graph allele remains a veto");
        run.recovery_source_sites[0].can_adopt = true;
        run.recovery_source_sites.erase(run.recovery_source_sites.begin());
        ok &= check(!bam_source_run_supported(run, kRunPhaseSet, true),
                    "shared BAM run: missing graph provenance remains a veto");
        run = make_run();
        run.recovery_source_quality_cuts[kSourcePhaseSet] = {150};
        ok &= check(!bam_source_run_supported(run, kRunPhaseSet),
                    "BAM run: an internal quality cut cannot certify a component");
        run.chunk.candidates[0].bam_injected = false;
        ok &= check(!bam_source_run_supported(run, kRunPhaseSet, true),
                    "shared BAM run: an internal quality cut remains a veto");
        run.recovery_source_quality_cuts[kSourcePhaseSet] = {250};
        ok &= check(bam_source_run_supported(run, kRunPhaseSet, true),
                    "shared BAM run: a later quality cut does not disconnect earlier rows");
        run.recovery_source_weak_cuts[kSourcePhaseSet] = {150};
        ok &= check(!bam_source_run_supported(run, kRunPhaseSet, true),
                    "shared BAM run: an internal weak cut remains a veto");
        run = make_run();
        run.recovery_source_sites[0].can_adopt = false;
        ok &= check(!bam_source_run_supported(run, kRunPhaseSet),
                    "BAM run: an untranslatable allele gauge is rejected");
        run = make_run();
        run.chunk.candidates.resize(1);
        ok &= check(!bam_source_run_supported(run, kRunPhaseSet),
                    "BAM run: a singleton cannot certify a path");
        run.chunk.candidates.front() = homozygote;
        ok &= check(!bam_source_run_supported(run, kRunPhaseSet),
                    "BAM run: homozygotes alone cannot certify a path");
    }

    {
        const auto make_island = [] {
            GraphChunkBuildResult gc;
            for (size_t ci = 0; ci < 4; ++ci) {
                CandidateVariant row;
                row.key.pos = ci == 3 ? 200 : 100 + 10 * ci;
                row.key.type = VariantType::Snp;
                row.counts.category = VariantCategory::CleanHetSnp;
                row.phase_set = ci == 2 ? 90 : 80;
                row.hap_to_cons_alle = ci == 3 ? std::array<int, 3>{-1, 0, 1} :
                    std::array<int, 3>{-1, 1, 0};
                row.bam_injected = ci != 2;
                gc.chunk.candidates.push_back(row);
                RecoverySourceSite source;
                source.candidate_index = ci;
                source.phase_set = 40;
                source.hap1_allele = 1;
                source.hap2_allele = 0;
                source.can_adopt = true;
                gc.recovery_source_sites.push_back(source);
            }
            CandidateVariant graph_anchor = gc.chunk.candidates[0];
            graph_anchor.key.pos = 105;
            graph_anchor.bam_injected = false;
            gc.chunk.candidates.push_back(graph_anchor);
            RecoverySourceSite shared = gc.recovery_source_sites[0];
            shared.candidate_index = 4;
            shared.clean_shared_snp = true;
            gc.recovery_source_sites.push_back(shared);
            gc.recovery_source_weak_cuts[40] = {130};
            for (int ri = 0; ri < 2; ++ri) {
                ReadRecord read;
                read.qname = "island-read-" + std::to_string(ri);
                gc.chunk.reads.push_back(std::move(read));
                ReadVariantProfile profile;
                profile.start_var_idx = 0;
                profile.end_var_idx = 3;
                profile.alleles = {-1, -1, ri == 1 ? 0 : -1, 1};
                gc.chunk.read_var_profile.push_back(std::move(profile));
                RecoverySourceRead source;
                source.read_index = ri;
                source.phase_set = 40;
                source.hap = 1;
                gc.recovery_source_reads.push_back(source);
            }
            gc.chunk.haps = {2, 2};
            gc.chunk.phase_sets = {80, 80};
            return gc;
        };
        auto gc = make_island();
        detach_bam_sites_across_weak_cuts(gc);
        const hts_pos_t detached_ps = gc.chunk.candidates[3].phase_set;
        ok &= check(detached_ps > 0 && detached_ps != 80 && detached_ps != 90,
                    "source island: the closest owner's label cannot hide an earlier owner's tail");
        ok &= check(gc.chunk.candidates[0].phase_set == 80 &&
                    gc.chunk.candidates[1].phase_set == 80 &&
                    gc.chunk.candidates[2].phase_set == 90,
                    "source island: established near components remain separate");
        ok &= check(gc.chunk.candidates[3].hap_to_cons_alle[1] == 1 &&
                    gc.chunk.candidates[3].hap_to_cons_alle[2] == 0,
                    "source island: detached allele gauge matches its original source");
        ok &= check(gc.chunk.phase_sets[0] == detached_ps && gc.chunk.haps[0] == 1,
                    "source island: an exclusive source read follows the restored gauge");
        ok &= check(gc.chunk.phase_sets[1] == 80 && gc.chunk.haps[1] == 2,
                    "source island: a read supporting another component keeps its near assignment");
        gc = make_island();
        gc.chunk.read_var_profile[1].alleles[0] = 0;
        detach_bam_sites_across_weak_cuts(gc);
        ok &= check(gc.chunk.candidates[3].phase_set == 80,
                    "source island: a physical bridge preserves an additional owner's join");
        gc = make_island();
        gc.chunk.candidates[3].bam_injected = false;
        detach_bam_sites_across_weak_cuts(gc);
        ok &= check(gc.chunk.candidates[3].phase_set == 80,
                    "source island: a certified far catalog anchor preserves its existing join");
        gc = make_island();
        gc.recovery_source_sites.back().clean_shared_snp = false;
        detach_bam_sites_across_weak_cuts(gc);
        ok &= check(gc.chunk.candidates[3].phase_set == 80,
                    "source island: an additional owner requires a shared clean SNP in this component");
        gc = make_island();
        gc.chunk.candidates[2].phase_set = 40;
        gc.recovery_source_sites.back().clean_shared_snp = false;
        detach_bam_sites_across_weak_cuts(gc);
        ok &= check(gc.chunk.candidates[3].phase_set == 80,
                    "source island: the original source label still owns the nearest row");
        gc = make_island();
        gc.recovery_source_weak_cuts[40] = {115, 130};
        detach_bam_sites_across_weak_cuts(gc);
        ok &= check(gc.chunk.candidates[3].phase_set == 80,
                    "source island: an older component's owner is not reopened at a later cut");
        gc = make_island();
        gc.recovery_source_weak_cuts[40].clear();
        detach_bam_sites_across_weak_cuts(gc);
        ok &= check(gc.chunk.candidates[3].phase_set == 80,
                    "source island: an uncut source is unchanged");
    }

    {
        constexpr hts_pos_t kCurrentPhaseSet = 80;
        constexpr hts_pos_t kOriginalPhaseSet = 40;
        const auto make_source = [kOriginalPhaseSet] {
            GraphChunkBuildResult result;
            for (const hts_pos_t pos : {100, 200}) {
                CandidateVariant row;
                row.key.type = pos == 100 ? VariantType::Deletion : VariantType::Snp;
                row.key.pos = pos;
                row.key.ref_len = 1;
                row.phase_set = kCurrentPhaseSet;
                row.hap_to_cons_alle = {-1, 0, 1};
                row.bam_injected = pos == 100;
                result.chunk.candidates.push_back(row);
                RecoverySourceSite source;
                source.candidate_index = result.chunk.candidates.size() - 1;
                source.phase_set = kOriginalPhaseSet;
                source.hap1_allele = 0;
                source.hap2_allele = 1;
                source.can_adopt = true;
                result.recovery_source_sites.push_back(source);
            }
            result.recovery_source_path_supported[kOriginalPhaseSet] = true;
            result.recovery_source_weak_cuts[kOriginalPhaseSet] = {};
            result.recovery_source_quality_cuts[kOriginalPhaseSet] = {};
            return result;
        };
        auto source = make_source();
        ok &= check(bam_source_site_path_supported(source, 0),
                    "BAM site path: original source anchors a relabeled deletion");
        source.chunk.candidates[0].hap_to_cons_alle = {-1, 1, 0};
        source.chunk.candidates[1].hap_to_cons_alle = {-1, 1, 0};
        ok &= check(bam_source_site_path_supported(source, 0),
                    "BAM site path: consistent source reversal remains valid");
        source.chunk.candidates[1].hap_to_cons_alle = {-1, 0, 1};
        ok &= check(!bam_source_site_path_supported(source, 0),
                    "BAM site path: deletion cannot reverse alone in its block");
        source = make_source();
        source.recovery_source_path_supported[kOriginalPhaseSet] = false;
        source.recovery_source_path_supported[kCurrentPhaseSet] = true;
        ok &= check(!bam_source_site_path_supported(source, 0),
                    "BAM site path: destination certificate cannot hide an incomplete source");
        source.recovery_source_path_supported.erase(kOriginalPhaseSet);
        ok &= check(!bam_source_site_path_supported(source, 0),
                    "BAM site path: destination certificate cannot replace missing provenance");
        source = make_source();
        source.recovery_source_weak_cuts[kOriginalPhaseSet] = {150};
        source.recovery_source_path_supported[kCurrentPhaseSet] = true;
        ok &= check(!bam_source_site_path_supported(source, 0),
                    "BAM site path: original weak cut vetoes a complete destination");
        source = make_source();
        source.recovery_source_quality_cuts[kOriginalPhaseSet] = {150};
        ok &= check(!bam_source_site_path_supported(source, 0),
                    "BAM site path: original quality cut remains a veto");
        source = make_source();
        source.recovery_source_weak_cuts[kCurrentPhaseSet] = {150};
        ok &= check(bam_source_site_path_supported(source, 0),
                    "BAM site path: unrelated destination cuts do not describe this source");
        source = make_source();
        source.recovery_source_sites.push_back(source.recovery_source_sites[0]);
        ok &= check(!bam_source_site_path_supported(source, 0),
                    "BAM site path: duplicate deletion provenance is ambiguous");
        source = make_source();
        source.recovery_source_sites[0].can_adopt = false;
        ok &= check(!bam_source_site_path_supported(source, 0),
                    "BAM site path: untranslatable deletion cannot certify a join");
        source = make_source();
        source.recovery_source_sites[1].can_adopt = false;
        ok &= check(!bam_source_site_path_supported(source, 0),
                    "BAM site path: untranslatable anchor cannot certify a gauge");
        source = make_source();
        source.recovery_source_sites[1].phase_set = kOriginalPhaseSet + 1;
        ok &= check(!bam_source_site_path_supported(source, 0),
                    "BAM site path: a different source cannot anchor the deletion");
        source = make_source();
        source.chunk.candidates[1].phase_set = kCurrentPhaseSet + 1;
        ok &= check(!bam_source_site_path_supported(source, 0),
                    "BAM site path: original source anchor must be in the current block");
        source = make_source();
        source.chunk.candidates[1].hap_to_cons_alle = {-1, 1, 1};
        ok &= check(!bam_source_site_path_supported(source, 0),
                    "BAM site path: a homozygote cannot certify orientation");
        source = make_source();
        source.chunk.candidates[1].counts.category = VariantCategory::NoisyCandHom;
        ok &= check(!bam_source_site_path_supported(source, 0),
                    "BAM site path: stale homozygous labels cannot anchor a gauge");
        source = make_source();
        source.chunk.candidates[0].counts.category = VariantCategory::CleanHom;
        ok &= check(!bam_source_site_path_supported(source, 0),
                    "BAM site path: stale homozygous deletion cannot certify parity");
        source = make_source();
        source.chunk.candidates[1].key.type = VariantType::Deletion;
        source.chunk.candidates[1].key.pos = 100;
        ok &= check(!bam_source_site_path_supported(source, 0),
                    "BAM site path: another allele at the same locus is not an anchor");
        source = make_source();
        source.chunk.candidates[0].hap_to_cons_alle = {-1, 0, 0};
        ok &= check(!bam_source_site_path_supported(source, 0),
                    "BAM site path: a homozygous deletion cannot certify parity");
        source = make_source();
        source.recovery_source_sites[0].hap2_allele = 2;
        ok &= check(!bam_source_site_path_supported(source, 0),
                    "BAM site path: a nonbiallelic source gauge is rejected");
        source = make_source();
        source.chunk.candidates[0].bam_injected = false;
        ok &= check(!bam_source_site_path_supported(source, 0),
                    "BAM site path: a graph-owned row needs graph certification");
        source = make_source();
        source.chunk.candidates[0].phase_set = kUnsetCandidatePhaseSet;
        ok &= check(!bam_source_site_path_supported(source, 0),
                    "BAM site path: an unphased row cannot orient a block");
        ok &= check(!bam_source_site_path_supported(source, 2),
                    "BAM site path: missing candidate is rejected");
        constexpr size_t kLongInsertionLength = 65;
        source = make_source();
        source.chunk.candidates[0].key.type = VariantType::Insertion;
        source.chunk.candidates[0].key.alt.assign(kLongInsertionLength, 'A');
        ok &= check(bam_source_site_path_supported(source, 0),
                    "BAM site path: a relabeled long insertion uses its original source");
        source.recovery_source_path_supported[kOriginalPhaseSet] = false;
        source.recovery_source_path_supported[kCurrentPhaseSet] = true;
        ok &= check(!bam_source_site_path_supported(source, 0),
                    "BAM site path: long insertion cannot borrow its destination certificate");
        source.recovery_source_path_supported[kOriginalPhaseSet] = true;
        source.recovery_source_quality_cuts[kOriginalPhaseSet] = {150};
        ok &= check(!bam_source_site_path_supported(source, 0),
                    "BAM site path: long insertion retains its original quality-cut veto");
        source.recovery_source_quality_cuts[kOriginalPhaseSet].clear();
        source.chunk.candidates[1].hap_to_cons_alle = {-1, 1, 0};
        ok &= check(!bam_source_site_path_supported(source, 0),
                    "BAM site path: long insertion must share its source anchor gauge");

        const auto make_prefix = [&make_source, kOriginalPhaseSet] {
            auto result = make_source();
            result.chunk.candidates[0].key.type = VariantType::Insertion;
            result.chunk.candidates[0].key.alt = "G";
            result.chunk.candidates[1].counts.category = VariantCategory::CleanHetSnp;
            result.recovery_source_sites[1].clean_shared_snp = true;
            result.recovery_source_path_supported[kOriginalPhaseSet] = false;
            result.recovery_source_weak_cuts[kOriginalPhaseSet] = {250};
            return result;
        };
        source = make_prefix();
        ok &= check(bam_source_prefix_to_graph_supported(source, 0),
                    "BAM prefix: a later source cut does not break the shared-SNP prefix");
        source.recovery_source_weak_cuts[kOriginalPhaseSet] = {150};
        ok &= check(!bam_source_prefix_to_graph_supported(source, 0),
                    "BAM prefix: an internal weak cut vetoes attachment");
        source = make_prefix();
        source.recovery_source_quality_cuts[kOriginalPhaseSet] = {150};
        ok &= check(!bam_source_prefix_to_graph_supported(source, 0),
                    "BAM prefix: an internal quality cut vetoes attachment");
        source = make_prefix();
        source.recovery_source_sites[1].clean_shared_snp = false;
        ok &= check(!bam_source_prefix_to_graph_supported(source, 0),
                    "BAM prefix: a nonshared graph row cannot anchor the source");
        source = make_prefix();
        source.chunk.candidates[1].hap_to_cons_alle = {-1, 1, 0};
        ok &= check(!bam_source_prefix_to_graph_supported(source, 0),
                    "BAM prefix: the shared SNP cannot reverse independently");
        source.chunk.candidates[0].hap_to_cons_alle = {-1, 1, 0};
        ok &= check(bam_source_prefix_to_graph_supported(source, 0),
                    "BAM prefix: consistent source reversal preserves attachment");
        source = make_prefix();
        source.recovery_source_sites.push_back(source.recovery_source_sites[1]);
        ok &= check(!bam_source_prefix_to_graph_supported(source, 0),
                    "BAM prefix: duplicate claims from the same source are ambiguous");
        source = make_prefix();
        auto foreign = source.recovery_source_sites[1];
        foreign.phase_set = kCurrentPhaseSet;
        foreign.hap1_allele = 1;
        foreign.hap2_allele = 0;
        source.recovery_source_sites.push_back(foreign);
        ok &= check(bam_source_prefix_to_graph_supported(source, 0),
                    "BAM prefix: a second solve does not invalidate the exact original gauge");
        source = make_prefix();
        auto foreign_row = source.chunk.candidates[1];
        foreign_row.key.pos = 150;
        source.chunk.candidates.push_back(foreign_row);
        ok &= check(!bam_source_prefix_to_graph_supported(source, 0),
                    "BAM prefix: every intervening anchor needs original provenance");
        source = make_prefix();
        source.recovery_source_weak_cuts.erase(kOriginalPhaseSet);
        ok &= check(!bam_source_prefix_to_graph_supported(source, 0),
                    "BAM prefix: missing path evidence is not a cut-free certificate");
        ok &= check(!bam_source_prefix_to_graph_supported(source, 9),
                    "BAM prefix: missing marker is rejected");
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

    // A broad graph/BAM gauge orients blocks but cannot bridge a local gap
    // with no read observing both sides. The next, directly supported seam
    // must still be eligible after that abstention.
    {
        constexpr hts_pos_t kLeft = 100;
        constexpr hts_pos_t kMiddle = 200;
        constexpr hts_pos_t kRight = 300;
        PhasingChunk replay;
        for (const auto [pos, phase_set] :
             {std::pair<hts_pos_t, hts_pos_t>{1000, kLeft},
              {2000, kMiddle}, {3000, kRight}}) {
            CandidateVariant candidate;
            candidate.key.pos = pos;
            candidate.key.type = VariantType::Snp;
            candidate.counts.n_uniq_alles = 2;
            candidate.counts.alle_covs = {12, 12};
            candidate.lcd_var_i_to_cate = kCandCleanHetSnp;
            candidate.hap_to_cons_alle = {-1, 0, 1};
            candidate.phase_set = phase_set;
            replay.candidates.push_back(std::move(candidate));
        }
        const auto add_read = [&](const std::string& name, int first,
                                  std::vector<int> alleles) {
            const int read_i = static_cast<int>(replay.reads.size());
            ReadRecord read;
            read.qname = name;
            read.mapq = 60;
            replay.reads.push_back(std::move(read));
            ReadVariantProfile profile;
            profile.read_id = read_i;
            profile.start_var_idx = first;
            profile.end_var_idx = first + static_cast<int>(alleles.size()) - 1;
            profile.alleles = std::move(alleles);
            profile.alt_qi.assign(profile.alleles.size(), 60);
            replay.read_var_profile.push_back(std::move(profile));
            replay.haps.push_back(0);
            replay.phase_sets.push_back(kUnphasedReadPhaseSet);
        };
        for (int allele = 0; allele <= 1; ++allele) {
            for (int copy = 0; copy < 6; ++copy) {
                const std::string suffix = std::to_string(allele) + "_" +
                                           std::to_string(copy);
                add_read("left_" + suffix, 0, {allele});
                add_read("middle_right_" + suffix, 1,
                         {allele, allele});
            }
        }
        rebuild_read_var_cr(replay);

        RecoveryPhaseGauge gauge;
        gauge.beg = 1000;
        gauge.end = 3000;
        gauge.graph_votes = {
            PhaseSetGaugeVote{kLeft, 12, 0},
            PhaseSetGaugeVote{kMiddle, 12, 0},
            PhaseSetGaugeVote{kRight, 12, 0},
        };
        Options stitch_opts;
        stitch_opts.min_block_link_reads = 1;
        stitch_opts.block_link_window = 8;
        const size_t joined = stitch_recovery_phase_sets_left_to_right(
            replay, {{1000, 2000, kLeft, kMiddle},
                     {2000, 3000, kMiddle, kRight}},
            {gauge}, stitch_opts);
        ok &= check(joined == 1 &&
                    replay.candidates[0].phase_set == kLeft &&
                    replay.candidates[1].phase_set == kMiddle &&
                    replay.candidates[2].phase_set == kMiddle,
                    "recovery gauge: zero-read seam stays split, later clean edge joins");
    }

    // A broad graph/BAM read gauge can disagree with the direct allele edge
    // across one seam. The molecule calls must determine the joined parity.
    {
        constexpr hts_pos_t kLeft = 100;
        constexpr hts_pos_t kRight = 200;
        PhasingChunk replay;
        for (const auto [pos, phase_set] :
             {std::pair<hts_pos_t, hts_pos_t>{1000, kLeft},
              {2000, kRight}}) {
            CandidateVariant candidate;
            candidate.key.pos = pos;
            candidate.key.type = VariantType::Snp;
            candidate.counts.n_uniq_alles = 2;
            candidate.counts.alle_covs = {6, 6};
            candidate.lcd_var_i_to_cate = kCandCleanHetSnp;
            candidate.hap_to_cons_alle = {-1, 0, 1};
            candidate.phase_set = phase_set;
            replay.candidates.push_back(std::move(candidate));
        }
        for (int allele = 0; allele <= 1; ++allele) {
            for (int copy = 0; copy < 6; ++copy) {
                const int read_i = static_cast<int>(replay.reads.size());
                ReadRecord read;
                read.qname = "cross_" + std::to_string(allele) + "_" +
                             std::to_string(copy);
                read.mapq = 60;
                replay.reads.push_back(std::move(read));
                ReadVariantProfile profile;
                profile.read_id = read_i;
                profile.start_var_idx = 0;
                profile.end_var_idx = 1;
                profile.alleles = {allele, 1 - allele};
                profile.alt_qi = {60, 60};
                replay.read_var_profile.push_back(std::move(profile));
                replay.haps.push_back(0);
                replay.phase_sets.push_back(kUnphasedReadPhaseSet);
            }
        }
        rebuild_read_var_cr(replay);
        RecoveryPhaseGauge gauge;
        gauge.beg = 1000;
        gauge.end = 2000;
        gauge.graph_votes = {
            PhaseSetGaugeVote{kLeft, 12, 0},
            PhaseSetGaugeVote{kRight, 12, 0},
        };
        Options stitch_opts;
        stitch_opts.min_block_link_reads = 1;
        stitch_opts.block_link_window = 8;
        const size_t joined = stitch_recovery_phase_sets_left_to_right(
            replay, {{1000, 2000, kLeft, kRight}}, {gauge}, stitch_opts);
        ok &= check(joined == 1 &&
                    replay.candidates[1].phase_set == kLeft &&
                    replay.candidates[1].hap_to_cons_alle[1] == 1,
                    "recovery gauge: contradictory direct SNP edge sets parity");
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
            left.key.type = VariantType::Snp;
            left.lcd_var_i_to_cate = kCandCleanHetSnp;
            left.hap_to_cons_alle = {-1, 0, 1};
            left.phase_set = left_phase_set;
            CandidateVariant right;
            right.key.pos = end + 9;  // deliberately differs from canonical end
            right.key.type = VariantType::Snp;
            right.lcd_var_i_to_cate = kCandCleanHetSnp;
            right.hap_to_cons_alle = {-1, 1, 0};
            right.phase_set = right_phase_set;
            replay.candidates = {left, right};
            for (int allele = 0; allele <= 1; ++allele) {
                for (int copy = 0; copy < 4; ++copy) {
                    const int read_i = static_cast<int>(replay.reads.size());
                    ReadRecord read;
                    read.qname = "target_" + std::to_string(allele) +
                                 "_" + std::to_string(copy);
                    read.mapq = 60;
                    replay.reads.push_back(std::move(read));
                    ReadVariantProfile profile;
                    profile.read_id = read_i;
                    profile.start_var_idx = 0;
                    profile.end_var_idx = 1;
                    profile.alleles = {allele, 1 - allele};
                    profile.alt_qi = {60, 60};
                    replay.read_var_profile.push_back(std::move(profile));
                    replay.haps.push_back(0);
                    replay.phase_sets.push_back(kUnphasedReadPhaseSet);
                }
            }
            rebuild_read_var_cr(replay);

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
            candidate.key.type = VariantType::Snp;
            candidate.lcd_var_i_to_cate = kCandCleanHetSnp;
            candidate.hap_to_cons_alle = {-1, 0, 1};
            candidate.phase_set = phase_set;
            replay.candidates.push_back(std::move(candidate));
        }
        for (int first = 0; first < 2; ++first) {
            for (int allele = 0; allele <= 1; ++allele) {
                for (int copy = 0; copy < 4; ++copy) {
                    const int read_i = static_cast<int>(replay.reads.size());
                    ReadRecord read;
                    read.qname = "consecutive_" + std::to_string(first) +
                                 "_" + std::to_string(allele) + "_" +
                                 std::to_string(copy);
                    read.mapq = 60;
                    replay.reads.push_back(std::move(read));
                    ReadVariantProfile profile;
                    profile.read_id = read_i;
                    profile.start_var_idx = first;
                    profile.end_var_idx = first + 1;
                    profile.alleles = {allele, allele};
                    profile.alt_qi = {60, 60};
                    replay.read_var_profile.push_back(std::move(profile));
                    replay.haps.push_back(0);
                    replay.phase_sets.push_back(kUnphasedReadPhaseSet);
                }
            }
        }
        rebuild_read_var_cr(replay);
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
        // The graph votes look decisive, but spanning reads already assigned
        // to the first BAM block contradict its boundary allele consensus.
        PhasingChunk inconsistent_source;
        inconsistent_source.candidates = replay.candidates;
        for (const ReadRecord& source : replay.reads) {
            ReadRecord read;
            read.qname = source.qname;
            read.mapq = source.mapq;
            read.is_skipped = source.qname.rfind("outer_", 0) == 0;
            inconsistent_source.reads.push_back(std::move(read));
        }
        inconsistent_source.read_var_profile = replay.read_var_profile;
        inconsistent_source.haps = replay.haps;
        inconsistent_source.phase_sets = replay.phase_sets;
        for (size_t read_i = 0; read_i < inconsistent_source.reads.size(); ++read_i) {
            const std::string& name = inconsistent_source.reads[read_i].qname;
            if (name.rfind("between_", 0) != 0) continue;
            const int allele =
                inconsistent_source.read_var_profile[read_i].alleles.front();
            inconsistent_source.haps[read_i] = allele == 0 ? 2 : 1;
            inconsistent_source.phase_sets[read_i] = kLocalA;
        }
        rebuild_read_var_cr(inconsistent_source);
        RecoveryPhaseGauge source_vote = gauge;
        source_vote.block_votes[0].shared_candidate_same = 20;
        source_vote.block_votes[1].counts = {{{8, 0}, {0, 8}}};
        source_vote.block_votes[1].shared_candidate_same = 20;
        source_vote.block_votes[1].shared_candidate_cross = 0;
        source_vote.graph_votes = {
            PhaseSetGaugeVote{kLeft, 16, 0},
            PhaseSetGaugeVote{kRight, 16, 0},
        };
        stitch_recovery_phase_sets_left_to_right(
            inconsistent_source, {{1000, 4000, kLeft, kRight}},
            {source_vote}, stitch_opts);
        ok &= check(
            inconsistent_source.candidates[0].phase_set !=
                inconsistent_source.candidates[3].phase_set,
            "recovery blocks: source HP conflict vetoes a false graph join");
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

    // Independent BAM blocks are accepted only when every informative graph
    // link is diploid, non-random, and has a narrow block-wide discordance
    // interval. This prevents high-depth weak links from passing on p-value
    // alone.
    {
        IndependentBamBlockLink clean;
        clean.counts = {{{15, 0}, {0, 15}}};
        ok &= check(independent_bam_block_is_supported({clean}),
                    "independent BAM block: clean diploid link passes");

        IndependentBamBlockLink random;
        random.counts = {{{8, 7}, {7, 8}}};
        ok &= check(!independent_bam_block_is_supported({clean, random}),
                    "independent BAM block: one mixed graph link rejects block");

        IndependentBamBlockLink high_depth_but_noisy;
        high_depth_but_noisy.counts = {{{45, 5}, {5, 45}}};
        ok &= check(!independent_bam_block_is_supported(
                        {high_depth_but_noisy}),
                    "independent BAM block: confidence bound rejects noise");

        PhasingChunk fallback;
        fallback.reads.resize(3);
        fallback.haps = {1, 0, 0};
        fallback.phase_sets = {100, 0, 0};
        fallback.gap_haps = {0, 2, 0};
        fallback.gap_phase_sets = {0, 200 + kGapFillPsOffset, 0};
        fallback.bam_fallback_haps = {2, 1, 2};
        fallback.bam_fallback_phase_sets = {
            300 + kBamFallbackPsOffset,
            300 + kBamFallbackPsOffset,
            300 + kBamFallbackPsOffset};
        ok &= check(apply_independent_bam_read_blocks(fallback) == 1 &&
                    fallback.haps[0] == 1 && fallback.phase_sets[0] == 100 &&
                    fallback.gap_haps[1] == 2 &&
                    fallback.gap_phase_sets[1] ==
                        200 + kGapFillPsOffset &&
                    fallback.gap_haps[2] == 2 &&
                    fallback.gap_phase_sets[2] ==
                        300 + kBamFallbackPsOffset,
                    "independent BAM block: primary and rescue assignments win");
    }

    // Keeping an MSA homopolymer genotype for linking must not promote one
    // noisy locus into unconditional read evidence. Use the same singleton
    // confidence rule as inferred markers, without changing its genotype.
    {
        const auto msa_singleton = [](int concordant, int conflicting) {
            PhasingChunk chunk;
            CandidateVariant site;
            site.key.pos = 1000;
            site.key.type = VariantType::Deletion;
            site.counts.n_uniq_alles = 2;
            site.lcd_var_i_to_cate = kCandNoisyCandHet;
            site.phase_set = 900;
            site.hap_to_cons_alle = {-1, 0, 1};
            site.bam_injected = true;
            site.msa_verified = true;
            site.is_homopolymer_indel = true;
            site.read_rescue_requires_validation = true;
            chunk.candidates.push_back(site);
            const int primary_reads = concordant + conflicting;
            for (int ri = 0; ri <= primary_reads; ++ri) {
                chunk.reads.emplace_back();
                const int hap = ri % 2 + 1;
                chunk.haps.push_back(ri < primary_reads ? hap : 0);
                chunk.phase_sets.push_back(ri < primary_reads ? 900 : 0);
                ReadVariantProfile profile;
                profile.read_id = ri;
                profile.start_var_idx = 0;
                profile.end_var_idx = 0;
                profile.alleles = {ri < concordant ? hap - 1 :
                                   ri < primary_reads ? 2 - hap : 0};
                chunk.read_var_profile.push_back(profile);
            }
            return chunk;
        };
        PhasingChunk noisy = msa_singleton(36, 12);
        ok &= check(rescue_unphased_graph_reads(noisy) == 0 &&
                    noisy.gap_haps.back() == 0,
                    "MSA singleton: significant but noisy locus stays unassigned");
        ok &= check(noisy.candidates[0].phase_set == 900 &&
                    noisy.candidates[0].hap_to_cons_alle[1] == 0 &&
                    noisy.candidates[0].hap_to_cons_alle[2] == 1,
                    "MSA singleton: rejection preserves linked genotype");
        PhasingChunk noisy_nonrepeat = msa_singleton(36, 12);
        noisy_nonrepeat.candidates[0].is_homopolymer_indel = false;
        ok &= check(rescue_unphased_graph_reads(noisy_nonrepeat) == 0,
                    "MSA singleton: noisy nonrepeat row needs the same confidence");
        PhasingChunk duplicate = msa_singleton(36, 12);
        duplicate.candidates.push_back(duplicate.candidates.front());
        for (auto& profile : duplicate.read_var_profile) {
            profile.end_var_idx = 1;
            profile.alleles.push_back(profile.alleles.front());
        }
        ok &= check(rescue_unphased_graph_reads(duplicate) == 0,
                    "MSA singleton: duplicate rows at one locus cannot certify rescue");
        PhasingChunk paired = msa_singleton(36, 12);
        paired.candidates.push_back(paired.candidates.front());
        paired.candidates.back().key.pos = 1100;
        for (auto& profile : paired.read_var_profile) {
            profile.end_var_idx = 1;
            profile.alleles.push_back(profile.alleles.front());
        }
        ok &= check(rescue_unphased_graph_reads(paired) == 1,
                    "MSA singleton: two independently located markers retain rescue");
        PhasingChunk clean = msa_singleton(32, 0);
        ok &= check(rescue_unphased_graph_reads(clean) == 1 &&
                    clean.gap_haps.back() == 1,
                    "MSA singleton: clean diploid primary support permits rescue");

        // The BAM genotype owns an independent block. Graph observations can
        // still rescue a read into an established block, in that block's gauge.
        const auto independent_msa = [&] {
            PhasingChunk independent = msa_singleton(64, 0);
            independent.candidates[0].graph_site = true;
            independent.candidates[0].bam_independent_genotype = true;
            independent.candidates[0].hap_to_cons_alle = {-1, 1, 0};
            for (size_t ri = 0; ri < independent.reads.size(); ++ri) {
                if (ri + 1 < independent.reads.size()) independent.phase_sets[ri] = 700;
                auto& profile = independent.read_var_profile[ri];
                profile.graph_alleles = profile.alleles;
                profile.alleles = {0};
                profile.bam_alleles = {0};
            }
            return independent;
        };
        PhasingChunk independent = independent_msa();
        PhasingChunk ambiguous = independent_msa();
        for (size_t ri = 32; ri + 1 < ambiguous.reads.size(); ++ri)
            ambiguous.phase_sets[ri] = 800;
        ok &= check(rescue_unphased_graph_reads(ambiguous) == 0 &&
                    ambiguous.gap_haps.back() == 0,
                    "independent MSA: two supported graph gauges cannot be pooled");
        PhasingChunk missing_graph = independent_msa();
        for (auto& profile : missing_graph.read_var_profile) profile.graph_alleles = {-1};
        ok &= check(rescue_unphased_graph_reads(missing_graph) == 0,
                    "independent MSA: graph rescue cannot borrow absent graph calls");
        ok &= check(rescue_unphased_graph_reads(independent) == 1 &&
                    independent.gap_haps.back() == 1 &&
                    independent.gap_phase_sets.back() == 700 + kGapFillPsOffset,
                    "independent MSA: graph association preserves established read gauge");
        ok &= check(independent.candidates[0].phase_set == 900 &&
                    independent.candidates[0].hap_to_cons_alle ==
                        std::array<int, 3>{-1, 1, 0},
                    "independent MSA: read rescue neither reorients nor stitches source row");
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

    // A repeat indel may supply one trustworthy read-only singleton when a
    // nearby phased clean SNP confirms its allele orientation in both read
    // halves. An unstable half or a site outside recovery cannot use the rule.
    {
        const auto make_direct_snp_fixture = [](bool unstable) {
            PhasingChunk chunk;
            CandidateVariant anchor;
            anchor.key.pos = 1000;
            anchor.key.type = VariantType::Snp;
            anchor.lcd_var_i_to_cate = kCandCleanHetSnp;
            anchor.phase_set = 1000;
            anchor.hap_to_cons_alle = {-1, 0, 1};
            chunk.candidates.push_back(anchor);
            CandidateVariant excluded;
            excluded.key.pos = 1100;
            excluded.key.type = VariantType::Deletion;
            excluded.counts.n_uniq_alles = 2;
            excluded.counts.category = VariantCategory::RepeatHetIndel;
            excluded.lcd_var_i_to_cate = kLongcalldRepHetVar;
            chunk.candidates.push_back(excluded);

            const auto fold_of = [](const std::string& qname) {
                uint64_t hash = 14695981039346656037ULL;
                for (const unsigned char byte : qname) {
                    hash ^= byte;
                    hash *= 1099511628211ULL;
                }
                return static_cast<size_t>(hash & 1ULL);
            };
            std::array<int, 2> fold_counts{};
            for (int name_i = 0;
                 fold_counts[0] < 18 || fold_counts[1] < 18;
                 ++name_i) {
                const std::string qname =
                    "direct_snp_" + std::to_string(name_i);
                const size_t fold = fold_of(qname);
                if (fold_counts[fold] == 18) continue;
                const int copy = fold_counts[fold]++;
                const int hap = copy % 2 == 0 ? 1 : 2;
                const int anchor_allele = hap - 1;
                const bool error = unstable
                    ? fold == 1 && copy >= 9
                    : copy >= 16;
                ReadRecord read;
                read.qname = qname;
                chunk.reads.push_back(std::move(read));
                ReadVariantProfile profile;
                profile.read_id = static_cast<int>(chunk.read_var_profile.size());
                profile.start_var_idx = 0;
                profile.end_var_idx = 1;
                profile.alleles = {
                    anchor_allele,
                    error ? 1 - anchor_allele : anchor_allele};
                chunk.read_var_profile.push_back(std::move(profile));
                chunk.haps.push_back(hap);
                chunk.phase_sets.push_back(1000);
            }
            for (int allele = 0; allele < 2; ++allele) {
                ReadRecord read;
                read.qname = "direct_unphased_" + std::to_string(allele);
                chunk.reads.push_back(std::move(read));
                ReadVariantProfile profile;
                profile.read_id = static_cast<int>(chunk.read_var_profile.size());
                profile.start_var_idx = 1;
                profile.end_var_idx = 1;
                profile.alleles = {allele};
                chunk.read_var_profile.push_back(std::move(profile));
                chunk.haps.push_back(0);
                chunk.phase_sets.push_back(kUnphasedReadPhaseSet);
            }
            return chunk;
        };
        const std::vector<RecoverySeam> windows = {{1050, 1150, 1000, 2000}};
        PhasingChunk outside = make_direct_snp_fixture(false);
        ok &= check(rescue_unphased_graph_reads(outside) == 0,
                    "graph direct SNP rescue: outside recovery abstains");
        PhasingChunk supported = make_direct_snp_fixture(false);
        ok &= check(rescue_unphased_graph_reads(supported, windows) == 2 &&
                    supported.gap_haps[36] == 1 &&
                    supported.gap_haps[37] == 2 &&
                    supported.gap_phase_sets[36] ==
                        1000 + kGapFillPsOffset &&
                    supported.candidates[1].phase_set ==
                        kUnsetCandidatePhaseSet,
                    "graph direct SNP rescue: split-fold link tags singleton reads");
        PhasingChunk unstable = make_direct_snp_fixture(true);
        ok &= check(rescue_unphased_graph_reads(unstable, windows) == 0,
                    "graph direct SNP rescue: unstable read half abstains");
        PhasingChunk discordant = make_direct_snp_fixture(false);
        discordant.read_var_profile[0].alleles[1] ^= 1;
        discordant.read_var_profile[1].alleles[1] ^= 1;
        ok &= check(rescue_unphased_graph_reads(discordant, windows) == 0,
                    "graph direct SNP rescue: uncertain link abstains");
        PhasingChunk ambiguous = make_direct_snp_fixture(false);
        CandidateVariant alternative = ambiguous.candidates[1];
        alternative.key.type = VariantType::Insertion;
        ambiguous.candidates.push_back(std::move(alternative));
        ok &= check(rescue_unphased_graph_reads(ambiguous, windows) == 0,
                    "graph direct SNP rescue: co-located alleles abstain");
        PhasingChunk sparse = make_direct_snp_fixture(false);
        for (int read_i = 0; read_i < 50; ++read_i) {
            ReadRecord read;
            read.qname = "sparse_site_" + std::to_string(read_i);
            sparse.reads.push_back(std::move(read));
            ReadVariantProfile profile;
            profile.read_id = static_cast<int>(sparse.read_var_profile.size());
            profile.start_var_idx = 1;
            profile.end_var_idx = 1;
            profile.alleles = {read_i % 2};
            sparse.read_var_profile.push_back(std::move(profile));
            sparse.haps.push_back(0);
            sparse.phase_sets.push_back(kUnphasedReadPhaseSet);
        }
        ok &= check(rescue_unphased_graph_reads(sparse, windows) == 0,
                    "graph direct SNP rescue: sparse anchor coverage abstains");
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

        // A BAM-only row must reach output even though it has no graph profile.
        // The same secondary channel must not replace a graph-primary row.
        constexpr hts_pos_t kBamOnlyPhaseSet =
            300 + kBamFallbackPsOffset;
        chunks[0].chunk.bam_output_fallback_reads = {
            {"bam_only", 2, kBamOnlyPhaseSet},
            {"read_a", initial.hap, kBamOnlyPhaseSet}};
        merge_graph_chunk_into_read_rows(output_rows, chunks[0], 0);
        ok &= check(output_rows.at("bam_only").hap == 2 &&
                    output_rows.at("bam_only").phase_set ==
                        kBamOnlyPhaseSet &&
                    output_rows.at("bam_only").has_phased_assignment &&
                    !output_rows.at("bam_only").has_primary_assignment,
                    "graph output merge: BAM-only assignment creates row");
        ok &= check(output_rows.at("read_a").hap == replacement_hap &&
                    output_rows.at("read_a").phase_set ==
                        initial.phase_set + 1000,
                    "graph output merge: primary assignment beats BAM-only fallback");
        chunks[0].chunk.bam_output_fallback_reads.clear();

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

    {
        const std::string header_text = "@SQ\tSN:chr\tLN:1000\n";
        const std::unique_ptr<bam_hdr_t, decltype(&bam_hdr_destroy)> header(
            sam_hdr_parse(header_text.size(), header_text.c_str()), bam_hdr_destroy);
        const std::unique_ptr<bam1_t, decltype(&bam_destroy1)> alignment(
            bam_init1(), bam_destroy1);
        const std::string sam = "quality\t0\tchr\t100\t60\t5M1D4M\t*\t0\t0\tACGTCGTTA\tIIIIIIIII";
        kstring_t line{0, 0, nullptr};
        kputs(sam.c_str(), &line);
        const int parsed = sam_parse1(&line, header.get(), alignment.get());
        std::free(line.s);
        ok &= check(parsed >= 0, "physical SNP quality: parse fixture");
        if (parsed >= 0) {
            const auto quality = [&](hts_pos_t pos, char ref, char alt, int allele) {
                return bam_snp_observation_quality(alignment.get(), pos, ref, alt, allele);
            };
            ok &= check(quality(100, 'a', 'g', 0) == 40 &&
                        quality(102, 'a', 'g', 1) == 40,
                        "physical SNP quality: matching REF/ALT keeps measured quality");
            ok &= check(quality(100, 'A', 'G', 1) == 0 &&
                        quality(102, 'A', 'G', 0) == 0,
                        "physical SNP quality: MSA cannot borrow opposite-allele quality");
            ok &= check(quality(101, 'A', 'G', 0) == 0 &&
                        quality(101, 'A', 'G', 1) == 0,
                        "physical SNP quality: a third BAM base is not a certificate");
            ok &= check(quality(105, 'A', 'G', 0) == 0 &&
                        quality(105, 'A', 'G', 1) == 0 &&
                        quality(99, 'A', 'G', 0) == 0 &&
                        quality(110, 'A', 'G', 0) == 0,
                        "physical SNP quality: deleted and uncovered positions abstain");
            bam_get_qual(alignment.get())[0] = 255;
            ok &= check(quality(100, 'A', 'G', 0) == 0,
                        "physical SNP quality: absent quality cannot certify an allele");
            bam_get_qual(alignment.get())[0] = 40;
            ok &= check(quality(100, 'A', 'G', -1) == 0 &&
                        quality(100, 'A', 'G', 2) == 0 &&
                        quality(100, 'A', 'A', 0) == 0 &&
                        quality(100, 'N', 'G', 0) == 0 &&
                        quality(100, 'A', 'N', 0) == 0 &&
                        bam_snp_observation_quality(nullptr, 100, 'A', 'G', 0) == 0,
                        "physical SNP quality: only a known binary nucleotide contrast certifies");
        }
    }

    {
        const auto make_seed = [] {
            PhasingChunk seed;
            seed.reads.resize(1);
            seed.haps = {2};
            seed.phase_sets = {100};
            seed.read_var_profile.resize(1);
            auto& profile = seed.read_var_profile[0];
            profile.read_id = 0;
            profile.start_var_idx = 0;
            profile.end_var_idx = 1;
            profile.alleles = {0, 0};
            profile.bam_alleles = {0, 0};
            profile.bam_base_qualities = {40, 40};
            profile.bam_mapq = 60;
            for (const hts_pos_t pos : {100, 250}) {
                CandidateVariant candidate;
                candidate.key.pos = pos;
                candidate.key.type = VariantType::Snp;
                candidate.key.ref_len = 1;
                candidate.key.alt = "T";
                candidate.counts.category = VariantCategory::CleanHetSnp;
                candidate.hap_to_cons_alle = {-1, 0, 1};
                candidate.phase_set = 100;
                seed.candidates.push_back(candidate);
            }
            return seed;
        };
        const auto seed = make_seed();
        auto corrected = make_seed();
        ok &= check(refresh_recovered_read_haps_from_bam_snps(corrected) == 1 &&
                    corrected.haps == std::vector<int>{1} &&
                    corrected.phase_sets == seed.phase_sets &&
                    corrected.reads[0].n_clean_agree_snps == 2 &&
                    corrected.reads[0].n_clean_conflict_snps == 0 &&
                    corrected.candidates[0].hap_to_cons_alle == seed.candidates[0].hap_to_cons_alle &&
                    corrected.candidates[1].phase_set == seed.candidates[1].phase_set,
                    "physical SNP refresh: repair inherited HP without changing the gauge");
        const auto make_certified_seed = [&] {
            auto chunk = make_seed();
            chunk.candidates[0].counts.category = VariantCategory::NoisyCandHet;
            chunk.candidates[0].msa_verified = true;
            chunk.candidates[0].bam_injected = true;
            chunk.read_var_profile[0].bam_base_qualities[1] = 0;
            // Sixteen conflict-free molecules pass the existing singleton
            // association/Wilson gates; both alleles need physical support.
            for (int i = 0; i < 16; ++i) {
                chunk.reads.emplace_back();
                chunk.reads.back().qname = "support_" + std::to_string(i);
                chunk.haps.push_back(0);
                chunk.phase_sets.push_back(100);
                auto profile = chunk.read_var_profile[0];
                profile.read_id = i + 1;
                profile.alleles = {i % 2, i % 2};
                profile.bam_alleles = profile.alleles;
                profile.bam_base_qualities = {40, 40};
                chunk.read_var_profile.push_back(std::move(profile));
            }
            return chunk;
        };
        auto msa = make_certified_seed();
        msa.candidates[0].counts.category = VariantCategory::NoisyCandHet;
        msa.candidates[0].msa_verified = true;
        msa.candidates[0].bam_injected = true;
        msa.read_var_profile[0].bam_base_qualities[1] = 0;
        ok &= check(refresh_recovered_read_haps_from_bam_snps(msa) == 1 &&
                    msa.haps[0] == 1 && msa.phase_sets[0] == seed.phase_sets[0] &&
                    std::all_of(msa.haps.begin() + 1, msa.haps.end(), [](int hap) { return hap == 0; }) &&
                    msa.reads[0].n_clean_agree_snps == 0,
                    "physical SNP refresh: independently verified MSA SNP corrects indel HP");
        msa = make_certified_seed();
        msa.candidates[0].counts.category = VariantCategory::NoisyCandHet;
        msa.candidates[0].msa_verified = true;
        msa.candidates[0].bam_injected = true;
        msa.read_var_profile[0].bam_alleles[1] = 1;
        msa.read_var_profile[0].alleles[1] = 1;
        msa.read_var_profile[0].bam_base_qualities[1] = 40;
        ok &= check(refresh_recovered_read_haps_from_bam_snps(msa) == 0 && msa.haps[0] == 2,
                    "physical SNP refresh: contradictory clean SNP vetoes MSA witness");
        msa = make_certified_seed();
        msa.read_var_profile[0].bam_alleles[1] = 1;
        msa.read_var_profile[0].alleles[1] = 1;
        msa.read_var_profile[0].bam_base_qualities[1] = 10;
        ok &= check(refresh_recovered_read_haps_from_bam_snps(msa) == 0 && msa.haps[0] == 2,
                    "physical SNP refresh: weak clean contradiction vetoes a noisy singleton");
        msa = make_certified_seed();
        msa.read_var_profile[0].bam_alleles[1] = 1;
        msa.read_var_profile[0].alleles[1] = 1;
        msa.read_var_profile[0].bam_base_qualities[1] = 0;
        ok &= check(refresh_recovered_read_haps_from_bam_snps(msa) == 1 && msa.haps[0] == 1,
                    "physical SNP refresh: an uncertified clean call cannot veto correction");
        msa = make_certified_seed();
        auto extra_clean = msa.candidates[1];
        extra_clean.key.pos = 500;
        msa.candidates.push_back(extra_clean);
        extra_clean.key.pos = 750;
        msa.candidates.push_back(extra_clean);
        for (auto& profile : msa.read_var_profile) {
            profile.end_var_idx = 3;
            profile.bam_base_qualities[1] = 40;
            profile.alleles.push_back(profile.alleles[0]);
            profile.bam_alleles.push_back(profile.bam_alleles[0]);
            profile.bam_base_qualities.push_back(40);
            profile.alleles.push_back(-1);
            profile.bam_alleles.push_back(-1);
            profile.bam_base_qualities.push_back(0);
        }
        msa.read_var_profile[0].alleles[3] = 1;
        msa.read_var_profile[0].bam_alleles[3] = 1;
        msa.read_var_profile[0].bam_base_qualities[3] = 10;
        ok &= check(refresh_recovered_read_haps_from_bam_snps(msa) == 1 &&
                    msa.haps[0] == 1 && msa.reads[0].n_clean_agree_snps == 2,
                    "physical SNP refresh: spaced Q30 clean witnesses retain their certificate");
        msa = make_certified_seed();
        msa.candidates[0].counts.category = VariantCategory::NoisyCandHet;
        msa.candidates[0].msa_verified = true;
        msa.candidates[0].bam_injected = false;
        msa.read_var_profile[0].bam_base_qualities[1] = 0;
        ok &= check(refresh_recovered_read_haps_from_bam_snps(msa) == 0,
                    "physical SNP refresh: MSA flag without imported source cannot certify HP");
        msa = make_certified_seed();
        auto farther = msa.candidates[1];
        farther.key.pos = 500;
        msa.candidates.push_back(std::move(farther));
        for (size_t ri = 0; ri < msa.read_var_profile.size(); ++ri) {
            auto& profile = msa.read_var_profile[ri];
            profile.end_var_idx = 2;
            profile.alleles.push_back(profile.alleles[0]);
            profile.bam_alleles.push_back(profile.bam_alleles[0]);
            profile.bam_base_qualities[1] = 0;
            profile.bam_base_qualities.push_back(ri == 0 ? 0 : 40);
        }
        ok &= check(refresh_recovered_read_haps_from_bam_snps(msa) == 1 && msa.haps[0] == 1,
                    "physical SNP refresh: farther clean anchors retain missing nearest-site evidence");
        msa = make_certified_seed();
        msa.candidates[1].hap_to_cons_alle = {-1, 1, 0};
        ok &= check(refresh_recovered_read_haps_from_bam_snps(msa) == 0 && msa.haps[0] == 2,
                    "physical SNP refresh: an MSA gauge reversed against clean SNPs abstains");
        msa = make_certified_seed();
        for (auto& profile : msa.read_var_profile) {
            profile.alleles = {0, 0};
            profile.bam_alleles = {0, 0};
        }
        ok &= check(refresh_recovered_read_haps_from_bam_snps(msa) == 0,
                    "physical SNP refresh: one-sided site support cannot certify a singleton");
        msa = make_certified_seed();
        msa.reads.resize(2);
        msa.haps.resize(2);
        msa.phase_sets.resize(2);
        msa.read_var_profile.resize(2);
        ok &= check(refresh_recovered_read_haps_from_bam_snps(msa) == 0,
                    "physical SNP refresh: one agreeing molecule cannot certify an MSA gauge");
        auto rejected = make_seed();
        rejected.read_var_profile[0].bam_base_qualities[1] = 29;
        ok &= check(refresh_recovered_read_haps_from_bam_snps(rejected) == 0 &&
                    rejected.haps == seed.haps,
                    "physical SNP refresh: one quality-bearing locus is insufficient");
        rejected = make_seed();
        rejected.read_var_profile[0].bam_mapq = 255;
        ok &= check(refresh_recovered_read_haps_from_bam_snps(rejected) == 0,
                    "physical SNP refresh: unknown MAPQ abstains");
        rejected = make_seed();
        rejected.read_var_profile[0].bam_mapq = 29;
        ok &= check(refresh_recovered_read_haps_from_bam_snps(rejected) == 0,
                    "physical SNP refresh: low MAPQ abstains");
        rejected = make_seed();
        rejected.read_var_profile[0].bam_alleles[1] = 1;
        rejected.read_var_profile[0].alleles[1] = 1;
        ok &= check(refresh_recovered_read_haps_from_bam_snps(rejected) == 0 &&
                    rejected.haps == seed.haps,
                    "physical SNP refresh: conflicting physical calls abstain");
        rejected = make_seed();
        rejected.read_var_profile[0].alleles[1] = 1;
        ok &= check(refresh_recovered_read_haps_from_bam_snps(rejected) == 0,
                    "physical SNP refresh: provenance disagreements cannot add support");
        rejected = make_seed();
        rejected.candidates[1].key.pos = 100;
        ok &= check(refresh_recovered_read_haps_from_bam_snps(rejected) == 0,
                    "physical SNP refresh: duplicate descriptions are one locus");
        rejected = make_seed();
        rejected.candidates[1].key.pos = 199;
        ok &= check(refresh_recovered_read_haps_from_bam_snps(rejected) == 0,
                    "physical SNP refresh: adjacent loci are not independent witnesses");
        rejected = make_seed();
        rejected.candidates[1].counts.category = VariantCategory::NoisyCandHet;
        ok &= check(refresh_recovered_read_haps_from_bam_snps(rejected) == 0,
                    "physical SNP refresh: noisy SNPs cannot replace clean witnesses");
        rejected = make_seed();
        rejected.candidates[1].phase_set = 200;
        ok &= check(refresh_recovered_read_haps_from_bam_snps(rejected) == 0,
                    "physical SNP refresh: independent phase sets cannot supply the gauge");
        rejected = make_seed();
        rejected.haps[0] = 0;
        ok &= check(refresh_recovered_read_haps_from_bam_snps(rejected) == 0,
                    "physical SNP refresh: unassigned reads remain for ordinary rescue");
    }

    {
        PhasingChunk coverage;
        coverage.reads.resize(3);
        for (ReadRecord& read : coverage.reads) {
            read.beg = 100;
            read.end = 200;
        }
        coverage.haps = {1, 2, 0};
        coverage.phase_sets = {7, 7, 0};
        coverage.gap_haps = {0, 0, 1};
        coverage.gap_phase_sets = {0, 0, 7 + kGapFillPsOffset};
        ok &= check(core_dominates_rescue_coverage(coverage, 7),
                    "partial rescue transfer: established core covers the cohort");
        coverage.reads[0].end = 150;
        coverage.reads[1].end = 150;
        ok &= check(!core_dominates_rescue_coverage(coverage, 7),
                    "partial rescue transfer: local rescue dominance vetoes splitting");
        coverage.reads[2].is_skipped = true;
        ok &= check(core_dominates_rescue_coverage(coverage, 7),
                    "partial rescue transfer: skipped reads do not define a cohort");
        coverage.reads[2].is_skipped = false;
        coverage.gap_phase_sets[2] = 8 + kGapFillPsOffset;
        ok &= check(core_dominates_rescue_coverage(coverage, 7),
                    "partial rescue transfer: independent cohorts are not combined");
    }

    {
        PhasingChunk masked;
        masked.candidates.resize(3);
        CandidateVariant& deletion = masked.candidates[0];
        deletion.key.type = VariantType::Deletion;
        deletion.key.pos = 100;
        deletion.key.ref_len = 6;
        deletion.bam_injected = true;
        deletion.msa_verified = true;
        deletion.counts.category = VariantCategory::NoisyCandHet;
        deletion.phase_set = 7;
        deletion.hap_to_cons_alle = {-1, 1, 0};
        for (size_t ci = 1; ci < masked.candidates.size(); ++ci) {
            CandidateVariant& snp = masked.candidates[ci];
            snp.key.type = VariantType::Snp;
            snp.key.pos = 100 + 3 * (ci - 1);
            snp.key.ref_len = 1;
            snp.key.alt = "T";
            snp.msa_verified = true;
            snp.counts.category = VariantCategory::NoisyCandHet;
            snp.phase_set = 7;
            snp.hap_to_cons_alle = {-1, 0, 1};
        }
        ReadVariantProfile profile;
        profile.start_var_idx = 0;
        profile.end_var_idx = 2;
        profile.alleles = {1, -1, -1};
        profile.bam_alleles = profile.alleles;
        ok &= check(masked_snp_deletion_haplotype(masked, profile, 0) == 1,
                    "masked SNP deletion: ALT identifies the missing SNP haplotype");
        for (CandidateVariant& candidate : masked.candidates)
            std::swap(candidate.hap_to_cons_alle[1], candidate.hap_to_cons_alle[2]);
        ok &= check(masked_snp_deletion_haplotype(masked, profile, 0) == 2,
                    "masked SNP deletion: flipped core gauge is preserved");
        for (CandidateVariant& candidate : masked.candidates)
            std::swap(candidate.hap_to_cons_alle[1], candidate.hap_to_cons_alle[2]);
        profile.alleles[1] = 1;
        ok &= check(masked_snp_deletion_haplotype(masked, profile, 0) == 0,
                    "masked SNP deletion: contrary phased allele vetoes filling");
        profile.alleles[1] = 0;
        masked.candidates[1].phase_set = 8;
        ok &= check(masked_snp_deletion_haplotype(masked, profile, 0) == 0,
                    "masked SNP deletion: another observed phase block vetoes filling");
        profile.alleles[1] = -1;
        masked.candidates[1].phase_set = 7;
        profile.bam_alleles[1] = 1;
        ok &= check(masked_snp_deletion_haplotype(masked, profile, 0) == 0,
                    "masked SNP deletion: a contrary BAM-only call vetoes filling");
        profile.bam_alleles[1] = -1;
        profile.bam_alleles[0] = -1;
        ok &= check(masked_snp_deletion_haplotype(masked, profile, 0) == 0,
                    "masked SNP deletion: both channels must retain deletion ALT");
        profile.bam_alleles[0] = 1;
        masked.candidates[1].key.pos = 110;
        masked.candidates[2].key.pos = 120;
        ok &= check(masked_snp_deletion_haplotype(masked, profile, 0) == 0,
                    "masked SNP deletion: nearby unobserved SNPs are not masked");
    }

    {
        GraphChunkBuildResult gc;
        gc.chunk.candidates.resize(1);
        CandidateVariant& marker = gc.chunk.candidates[0];
        marker.key.tid = 11;
        marker.key.pos = 100;
        marker.key.type = VariantType::Snp;
        marker.key.ref_len = 1;
        marker.key.alt = "T";
        marker.phase_set = 7;
        marker.counts.category = VariantCategory::CleanHetSnp;
        marker.hap_to_cons_alle = {-1, 0, 1};
        marker.counts.n_uniq_alles = 2;
        gc.site_meta.push_back({"chr20", 100, "A", {"T"}});
        gc.site_allele_orig_idx.push_back({0, 1});
        RecoverySourceSite source;
        source.candidate_index = 0;
        source.phase_set = 9;
        source.hap1_allele = 0;
        source.hap2_allele = 1;
        source.can_adopt = true;
        source.msa_key = marker.key;
        gc.recovery_source_sites.push_back(source);
        gc.recovery_source_path_supported[9] = true;
        gc.recovery_source_weak_cuts[9] = {};
        RecoveryPhaseGauge gauge;
        for (int si = 0; si < 3; ++si) {
            VariantKey key = marker.key;
            key.pos += si * 10;
            gauge.bam_sites.push_back({key, 9, 0, 1, si != 0, 9});
        }
        for (int ri = 0; ri < 12; ++ri) {
            const int allele = ri % 2;
            gauge.bam_reads.push_back({std::to_string(ri), 60,
                {{0, allele}, {1, allele}, {2, allele}}, {40, 40, 40}});
        }
        gc.recovery_phase_gauges.push_back(gauge);
        ok &= check(retained_source_snp_path_anchor_supported(gc, 0),
                    "source SNP anchor: recalled allele uses independent physical clean loci");
        gc.recovery_source_weak_cuts[9] = {105};
        ok &= check(!retained_source_snp_path_anchor_supported(gc, 0),
                    "source SNP anchor: an unsupported source cut cannot certify the path");
        gc.recovery_source_weak_cuts[9].clear();
        gc.recovery_source_quality_cuts[9] = {105};
        ok &= check(!retained_source_snp_path_anchor_supported(gc, 0),
                    "source SNP anchor: source quality cuts remain a veto");
        gc.recovery_source_quality_cuts[9].clear();
        gc.recovery_source_path_supported[9] = false;
        ok &= check(!retained_source_snp_path_anchor_supported(gc, 0),
                    "source SNP anchor: complete independent source is required");
        gc.recovery_source_path_supported[9] = true;
        gc.recovery_source_sites.push_back(source);
        ok &= check(!retained_source_snp_path_anchor_supported(gc, 0),
                    "source SNP anchor: competing source claims are rejected");
        gc.recovery_source_sites.pop_back();
        gc.recovery_source_sites[0].msa_key.reset();
        ok &= check(!retained_source_snp_path_anchor_supported(gc, 0),
                    "source SNP anchor: recalled genotype provenance is required");
        gc.recovery_source_sites[0].msa_key = marker.key;
        gc.site_meta[0].alts[0] = "G";
        ok &= check(!retained_source_snp_path_anchor_supported(gc, 0),
                    "source SNP anchor: graph and source must describe the same physical SNP");
        gc.site_meta[0].alts[0] = "T";
        gc.recovery_phase_gauges[0].bam_reads[0].observations[1].second = 1;
        ok &= check(!retained_source_snp_path_anchor_supported(gc, 0),
                    "source SNP anchor: a contrary clean locus vetoes the recalled gauge");
        gc.recovery_phase_gauges[0] = gauge;
        gc.recovery_phase_gauges[0].bam_sites[2].key.pos = 110;
        ok &= check(!retained_source_snp_path_anchor_supported(gc, 0),
                    "source SNP anchor: duplicate descriptions are not independent loci");
        gc.recovery_phase_gauges[0] = gauge;
        for (RecoveryBamRead& read : gc.recovery_phase_gauges[0].bam_reads)
            read.base_qualities[0] = 0;
        ok &= check(!retained_source_snp_path_anchor_supported(gc, 0),
                    "source SNP anchor: deletion-masked or missing physical bases abstain");
        gc.recovery_phase_gauges[0] = gauge;
        for (RecoveryBamRead& read : gc.recovery_phase_gauges[0].bam_reads)
            if (read.observations[0].second == 1) read.mapq = 10;
        ok &= check(!retained_source_snp_path_anchor_supported(gc, 0),
                    "source SNP anchor: both physical allele classes need mapping support");
        gc.recovery_phase_gauges[0] = gauge;
        gc.recovery_phase_gauges[0].bam_reads.resize(8);
        ok &= check(!retained_source_snp_path_anchor_supported(gc, 0),
                    "source SNP anchor: weak physical counts cannot certify parity");
    }

    {
        IndependentBamBlockLink gauge;
        gauge.counts = {{{30, 0}, {0, 30}}};
        const double q17_error = std::pow(10.0, -1.7) + 0.000302;
        const double q17_odds = std::log((1.0 - q17_error) / q17_error);
        const auto crossed = calibrated_indel_bridge_flip(gauge, {0, 1}, q17_odds);
        ok &= check(crossed && *crossed,
                    "shared deletion parity: calibrated Q17 primary pair meets the 5% error bound");
        const auto same = calibrated_indel_bridge_flip(gauge, {1, 0}, -q17_odds);
        ok &= check(same && !*same,
                    "shared deletion parity: the same-gauge orientation is retained");
        ok &= check(!calibrated_indel_bridge_flip(gauge, {1, 1}, 10.0),
                    "shared deletion parity: a contrary molecule vetoes the bridge");
        ok &= check(!calibrated_indel_bridge_flip(gauge, {0, 0}, 10.0),
                    "shared deletion parity: no physical pair means no bridge");
        ok &= check(!calibrated_indel_bridge_flip(gauge, {0, 1}, 2.0),
                    "shared deletion parity: weak physical quality cannot certify the gauge");
        ok &= check(!calibrated_indel_bridge_flip(gauge, {0, 1}, -q17_odds),
                    "shared deletion parity: likelihood and observed orientation must agree");
        ok &= check(calibrated_indel_bridge_flip(gauge, {0, 2}, 2.0 * q17_odds).has_value(),
                    "shared deletion parity: added agreeing evidence preserves acceptance");
        gauge.counts = {{{7, 0}, {0, 14}}};
        ok &= check(calibrated_indel_bridge_flip(gauge, {0, 1}, q17_odds).has_value(),
                    "shared deletion parity: sparse precise evidence meets the joint 80% bound");
        gauge.counts = {{{11, 0}, {1, 13}}};
        const auto insertion = calibrated_indel_bridge_flip(gauge, {1, 0}, -5.9507);
        ok &= check(insertion && !*insertion,
                    "insertion parity: diploid calibration tolerates one discordant donor within the 80% bound");
        ok &= check(!calibrated_indel_bridge_flip(gauge, {1, 1}, -5.9507),
                    "insertion parity: a contrary bridge cannot be hidden by calibration depth");
        gauge.counts = {{{10, 0}, {1, 12}}};
        const auto tandem = calibrated_indel_bridge_flip(gauge, {0, 1}, 7.23801);
        ok &= check(tandem && *tandem,
                    "tandem insertion: separate diploid calibration and a precise bridge meet the joint bound");
        ok &= check(!calibrated_indel_bridge_flip(gauge, {1, 1}, 7.23801),
                    "tandem insertion: a contrary physical bridge vetoes the join");
        gauge.counts = {{{10, 0}, {5, 12}}};
        ok &= check(!calibrated_indel_bridge_flip(gauge, {0, 1}, 7.23801),
                    "tandem insertion: precise parity cannot replace adequate allele calibration");
        gauge.counts = {{{10, 0}, {0, 2}}};
        const auto mixed = calibrated_indel_bridge_flip(gauge, {1, 0}, -7.564577);
        ok &= check(mixed && !*mixed,
                    "mixed repeat: complementary insertion/deletion calibration meets the joint 80% bound");
        ok &= check(!calibrated_indel_bridge_flip(gauge, {1, 1}, -7.564577),
                    "mixed repeat: a contrary physical bridge vetoes the join");
        gauge.counts = {{{10, 0}, {0, 1}}};
        ok &= check(!calibrated_indel_bridge_flip(gauge, {1, 0}, -7.564577),
                    "mixed repeat: nearest-SNP-only calibration exceeds the joint error budget");
        gauge.counts = {{{2, 0}, {0, 8}}};
        ok &= check(!calibrated_indel_bridge_flip(gauge, {0, 1}, q17_odds),
                    "shared deletion parity: insufficient calibration fails the joint error bound");
        gauge.counts = {{{30, 10}, {10, 30}}};
        ok &= check(!calibrated_indel_bridge_flip(gauge, {0, 1}, 20.0),
                    "shared deletion parity: high depth cannot conceal excessive calibration discordance");
        gauge.counts = {{{30, 0}, {0, 0}}};
        ok &= check(!calibrated_indel_bridge_flip(gauge, {0, 1}, q17_odds),
                    "shared deletion parity: calibration requires both haplotypes");
        gauge.counts = {{{0, 30}, {30, 0}}};
        ok &= check(!calibrated_indel_bridge_flip(gauge, {0, 1}, q17_odds),
                    "shared deletion parity: a reversed marker gauge cannot supply the bridge");
    }

    {
        IndependentBamBlockLink gauge;
        gauge.counts = {{{7, 0}, {0, 4}}};
        const auto same = calibrated_source_deletion_bridge_flip(gauge, {1, 0}, 0.135193);
        ok &= check(same && !*same, "source deletion: diploid calibration permits a precise sparse bridge");
        const auto reverse = calibrated_source_deletion_bridge_flip(gauge, {0, 1}, 0.135193);
        ok &= check(reverse && *reverse, "source deletion: reverse orientation is preserved");
        ok &= check(!calibrated_source_deletion_bridge_flip(gauge, {1, 1}, 0.01),
                    "source deletion: any contrary bridge vetoes the union");
        ok &= check(!calibrated_source_deletion_bridge_flip(gauge, {0, 0}, 0.01),
                    "source deletion: calibration alone cannot join blocks");
        ok &= check(!calibrated_source_deletion_bridge_flip(gauge, {1, 0}, 0.200001),
                    "source deletion: the actual joint quality bound cannot exceed 20%");
        ok &= check(!calibrated_source_deletion_bridge_flip(gauge, {1, 0}, -0.1),
                    "source deletion: an invalid quality bound abstains");
        gauge.counts[0][1] = 1;
        ok &= check(!calibrated_source_deletion_bridge_flip(gauge, {1, 0}, 0.01),
                    "source deletion: contrary calibration cannot be outvoted");
        gauge.counts = {{{7, 0}, {0, 1}}};
        ok &= check(!calibrated_source_deletion_bridge_flip(gauge, {1, 0}, 0.01),
                    "source deletion: both allele classes need independent calibration");
        gauge.counts = {{{2, 0}, {0, 2}}};
        ok &= check(!calibrated_source_deletion_bridge_flip(gauge, {1, 0}, 0.01),
                    "source deletion: weak association cannot certify the gauge");
    }

    {
        GraphChunkBuildResult shared;
        shared.chunk.candidates.resize(1);
        CandidateVariant& marker = shared.chunk.candidates[0];
        marker.key.type = VariantType::Deletion;
        marker.key.pos = 100;
        marker.alignment_verified = true;
        marker.counts.n_uniq_alles = 2;
        marker.counts.category = VariantCategory::NoisyCandHet;
        marker.phase_set = 7;
        marker.hap_to_cons_alle = {-1, 0, 1};
        shared.site_meta.push_back({"chr20", 100, "AT", {"A"}});
        shared.site_allele_orig_idx.push_back({0, 1});
        RecoverySourceSite source;
        source.candidate_index = 0;
        source.can_adopt = true;
        VariantKey key;
        key.pos = 100;
        key.type = VariantType::Deletion;
        key.ref_len = 1;
        source.msa_key = key;
        shared.recovery_source_sites.push_back(source);
        const auto physical = retained_shared_deletion_key(shared, 0);
        ok &= check(physical && physical->pos == 101 && physical->ref_len == 1 &&
                    physical->alt.empty() && !marker.msa_verified,
                    "shared MSA deletion: catalog keeps its physical edit and source proof");
        marker.alignment_verified = false;
        ok &= check(!retained_shared_deletion_key(shared, 0),
                    "shared MSA deletion: unverified graph allele is excluded");
        marker.alignment_verified = true;
        shared.recovery_source_sites[0].msa_key.reset();
        ok &= check(!retained_shared_deletion_key(shared, 0),
                    "shared MSA deletion: retained MSA provenance is required");
        shared.recovery_source_sites[0].msa_key = key;
        shared.recovery_source_sites.push_back(source);
        ok &= check(!retained_shared_deletion_key(shared, 0),
                    "shared MSA deletion: competing source claims are excluded");
        shared.recovery_source_sites.pop_back();
        shared.site_meta[0].alts[0] = "GC";
        ok &= check(!retained_shared_deletion_key(shared, 0),
                    "shared MSA deletion: a replacement is not a pure deletion");
    }

    {
        PhasingChunk chunk;
        chunk.candidates.resize(3);
        for (CandidateVariant& site : chunk.candidates) {
            site.key.type = VariantType::Snp;
            site.key.ref_len = 1;
            site.counts.n_uniq_alles = 2;
            site.counts.category = VariantCategory::CleanHetSnp;
            site.phase_set = 7;
            site.hap_to_cons_alle = {-1, 0, 1};
        }
        ReadVariantProfile profile;
        profile.start_var_idx = 0;
        profile.alleles = {0, 1, -1};
        profile.graph_alleles = profile.alleles;
        profile.bam_alleles = profile.alleles;
        profile.bam_base_qualities = {0, 40, 0};
        ok &= check(masked_bam_snp_haplotype(chunk, profile, 1, {0}) == 2,
                    "masked BAM SNP: a deleted REF can abstain behind a physical witness");
        ok &= check(masked_bam_snp_haplotype(chunk, profile, 1, {}) == 0,
                    "masked BAM SNP: an ordinary contrary REF remains a veto");
        profile.bam_base_qualities[1] = 17;
        ok &= check(masked_bam_snp_haplotype(chunk, profile, 1, {0}) == 0,
                    "masked BAM SNP: a low quality witness cannot connect the read");
        profile.bam_base_qualities[1] = 40;
        profile.graph_alleles[1] = 0;
        ok &= check(masked_bam_snp_haplotype(chunk, profile, 1, {0}) == 0,
                    "masked BAM SNP: graph and BAM witness alleles must agree");
        profile.graph_alleles[1] = 1;
        profile.alleles[2] = 0;
        ok &= check(masked_bam_snp_haplotype(chunk, profile, 1, {0}) == 0,
                    "masked BAM SNP: another unmasked contrary call vetoes connection");
        profile.alleles[2] = -1;
        chunk.candidates[0].phase_set = 8;
        ok &= check(masked_bam_snp_haplotype(chunk, profile, 1, {0}) == 0,
                    "masked BAM SNP: deleted calls cannot connect a different phase set");
        chunk.candidates[0].phase_set = 7;
        profile.bam_base_qualities[0] = 40;
        ok &= check(masked_bam_snp_haplotype(chunk, profile, 1, {0}) == 0,
                    "masked BAM SNP: a quality-bearing REF cannot be masked");
    }

    {
        std::array<IndependentBamBlockLink, 2> cohorts;
        for (auto& cohort : cohorts) cohort.counts = {{{20, 0}, {0, 20}}};
        ok &= check(calibrated_terminal_insertion_hap1(cohorts) == 0,
                    "terminal insertion: both independent cohorts orient REF on hap1");
        cohorts[0].counts = {{{14, 0}, {0, 7}}};
        cohorts[1].counts = {{{10, 0}, {0, 7}}};
        ok &= check(calibrated_terminal_insertion_hap1(cohorts) == 0,
                    "terminal insertion: independent physical calls certify the endpoint");
        cohorts[0].counts = {{{7, 0}, {0, 10}}};
        cohorts[1].counts = {{{9, 0}, {0, 8}}};
        ok &= check(calibrated_terminal_insertion_hap1(cohorts) == 0,
                    "terminal mixed repeat: separate joint ALT classes certify the source gauge");
        cohorts[1].counts = {{{0, 9}, {8, 0}}};
        ok &= check(!calibrated_terminal_insertion_hap1(cohorts),
                    "terminal mixed repeat: graph/source gauge disagreement vetoes the endpoint");
        for (auto& cohort : cohorts) cohort.counts = {{{7, 0}, {0, 6}}};
        ok &= check(!calibrated_terminal_insertion_hap1(cohorts),
                    "terminal insertion: adequate separate cohorts still need the combined bound");
        for (auto& cohort : cohorts) cohort.counts = {{{0, 20}, {20, 0}}};
        ok &= check(calibrated_terminal_insertion_hap1(cohorts) == 1,
                    "terminal insertion: both independent cohorts orient ALT on hap1");
        cohorts[1].counts = {{{20, 0}, {0, 20}}};
        ok &= check(!calibrated_terminal_insertion_hap1(cohorts),
                    "terminal insertion: opposite cohort gauges cannot phase the endpoint");
        cohorts[1].counts = {{{0, 40}, {0, 0}}};
        ok &= check(!calibrated_terminal_insertion_hap1(cohorts),
                    "terminal insertion: one observed allele cannot establish a heterozygote");
        cohorts[1].counts = {{{0, 4}, {4, 0}}};
        ok &= check(!calibrated_terminal_insertion_hap1(cohorts),
                    "terminal insertion: sparse evidence cannot certify an endpoint");
        cohorts[1].counts = {{{10, 30}, {30, 10}}};
        ok &= check(!calibrated_terminal_insertion_hap1(cohorts),
                    "terminal insertion: discordance above the calibrated bound abstains");
    }

    {
        std::array<IndependentBamBlockLink, 2> cohorts;
        cohorts[0].counts = {{{1, 8}, {8, 0}}};
        cohorts[1].counts = {{{0, 8}, {8, 0}}};
        ok &= check(calibrated_verified_insertion_hap1(cohorts) == 1,
                    "verified insertion: independent shifted calls retain the source gauge");
        ok &= check(!calibrated_terminal_insertion_hap1(cohorts),
                    "verified insertion: retaining a source does not weaken the terminal gate");
        cohorts[1].counts = {{{8, 0}, {0, 8}}};
        ok &= check(!calibrated_verified_insertion_hap1(cohorts),
                    "verified insertion: opposite independent gauges abstain");
        cohorts[1].counts = {{{0, 16}, {0, 0}}};
        ok &= check(!calibrated_verified_insertion_hap1(cohorts),
                    "verified insertion: each cohort must contain both alleles");
        for (auto& cohort : cohorts) cohort.counts = {{{10, 30}, {30, 10}}};
        ok &= check(!calibrated_verified_insertion_hap1(cohorts),
                    "verified insertion: confidence bound must retain eighty percent");
        for (auto& cohort : cohorts) cohort.counts = {{{0, 3}, {3, 0}}};
        ok &= check(!calibrated_verified_insertion_hap1(cohorts),
                    "verified insertion: sparse associations cannot certify a gauge");
    }

    {
        std::array<IndependentBamBlockLink, 2> gauges;
        gauges[0].counts = {{{12, 3}, {0, 8}}};
        gauges[1].counts = {{{9, 2}, {1, 9}}};
        ok &= check(calibrated_deletion_chain_flip(gauges, {4, 1}) == false,
                    "deletion chain: calibrated length errors permit an 80 percent bridge");
        ok &= check(calibrated_deletion_chain_flip(gauges, {1, 4}) == true,
                    "deletion chain: opposite parity preserves the calibrated gauge");
        ok &= check(!calibrated_deletion_chain_flip(gauges, {3, 2}),
                    "deletion chain: a bridge below 80 percent remains open");
        ok &= check(!calibrated_deletion_chain_flip(gauges, {2, 0}),
                    "deletion chain: two molecules cannot establish the connection");
        ok &= check(!calibrated_deletion_chain_flip(gauges, {0, 0}),
                    "deletion chain: calibration alone cannot join noncrossing components");
        gauges[1].counts = {{{9, 4}, {2, 9}}};
        ok &= check(!calibrated_deletion_chain_flip(gauges, {4, 1}),
                    "deletion chain: a poorly calibrated repeat remains excluded");
        gauges[1].counts = {{{20, 0}, {0, 0}}};
        ok &= check(!calibrated_deletion_chain_flip(gauges, {4, 1}),
                    "deletion chain: both diploid classes need independent calibration");
        gauges[1].counts = {{{2, 0}, {0, 2}}};
        ok &= check(!calibrated_deletion_chain_flip(gauges, {4, 1}),
                    "deletion chain: sparse calibration cannot certify the gauge");
    }

    {
        IndependentBamBlockLink gauge;
        gauge.counts = {{{28, 1}, {5, 28}}};
        const std::array<std::array<int, 2>, 2> bridge{{{5, 0}, {0, 2}}};
        ok &= check(calibrated_repeat_chain_flip(gauge, bridge, {5, 1}, -30.0) == false,
                    "repeat chain: eighty-percent SNP agreement tolerates one noisy length");
        ok &= check(calibrated_repeat_chain_flip(gauge, bridge, {1, 5}, 30.0) == true,
                    "repeat chain: independent SNPs can reverse the downstream gauge");
        ok &= check(!calibrated_repeat_chain_flip(gauge, bridge, {4, 2}, -30.0),
                    "repeat chain: below eighty-percent agreement cannot join");
        ok &= check(!calibrated_repeat_chain_flip(gauge, bridge, {2, 0}, -30.0),
                    "repeat chain: two right molecules cannot replace calibration");
        ok &= check(!calibrated_repeat_chain_flip(gauge, bridge, {5, 1}, -2.0),
                    "repeat chain: weak physical quality cannot certify orientation");
        ok &= check(!calibrated_repeat_chain_flip(gauge, {{{5, 1}, {0, 2}}}, {5, 1}, -30.0),
                    "repeat chain: a contrary molecule across the interior vetoes joining");
        ok &= check(!calibrated_repeat_chain_flip(gauge, {{{5, 0}, {0, 0}}}, {5, 1}, -30.0),
                    "repeat chain: both physical repeat classes must cross");
        gauge.counts = {{{28, 20}, {20, 28}}};
        ok &= check(!calibrated_repeat_chain_flip(gauge, bridge, {5, 1}, -30.0),
                    "repeat chain: an unsupported left gauge cannot seed recovery");
    }

    {
        std::array<IndependentBamBlockLink, 2> gauges;
        for (auto& gauge : gauges) gauge.counts = {{{50, 0}, {0, 50}}};
        ok &= check(calibrated_repeat_snp_bridge_flip(gauges, {1, 0}, 0.0004) == false,
                    "repeat SNP bridge: calibrated sparse same-gauge molecule joins");
        ok &= check(calibrated_repeat_snp_bridge_flip(gauges, {0, 1}, 0.0004) == true,
                    "repeat SNP bridge: calibrated reverse-gauge molecule joins");
        ok &= check(!calibrated_repeat_snp_bridge_flip(gauges, {1, 1}, 0.0004),
                    "repeat SNP bridge: contrary independent molecule vetoes union");
        ok &= check(!calibrated_repeat_snp_bridge_flip(gauges, {0, 0}, 0.0004),
                    "repeat SNP bridge: calibration without a physical bridge abstains");
        ok &= check(!calibrated_repeat_snp_bridge_flip(gauges, {1, 0}, 0.002),
                    "repeat SNP bridge: weak base or mapping quality abstains");
        auto unsupported = gauges;
        unsupported[1].counts = {{{100, 0}, {0, 0}}};
        ok &= check(!calibrated_repeat_snp_bridge_flip(unsupported, {1, 0}, 0.0004),
                    "repeat SNP bridge: both diploid classes must calibrate each flank");
        unsupported[1].counts = {{{0, 50}, {50, 0}}};
        ok &= check(!calibrated_repeat_snp_bridge_flip(unsupported, {1, 0}, 0.0004),
                    "repeat SNP bridge: physical calls cannot reverse an established gauge");
        unsupported[1].counts = {{{2, 0}, {0, 2}}};
        ok &= check(!calibrated_repeat_snp_bridge_flip(unsupported, {1, 0}, 0.0004),
                    "repeat SNP bridge: tiny calibration cannot replace a full path proof");
        unsupported[1].counts = {{{40, 10}, {10, 40}}};
        ok &= check(!calibrated_repeat_snp_bridge_flip(unsupported, {1, 0}, 0.0004),
                    "repeat SNP bridge: joint error above twenty percent abstains");
    }

    {
        ok &= check(complementary_insertion_length_class(7, 4, 8) == 1,
                    "repeat insertion: seven bases select the eight-base class");
        ok &= check(complementary_insertion_length_class(1, 4, 8) == 0,
                    "repeat insertion: slippage remains a non-reference short class");
        ok &= check(complementary_insertion_length_class(6, 4, 8) == -1 &&
                    complementary_insertion_length_class(0, 4, 8) == -1,
                    "repeat insertion: ties and reference observations abstain");
        IndependentBamBlockLink gauge;
        gauge.counts = {{{4, 0}, {0, 5}}};
        const std::array<std::array<int, 2>, 2> reverse{{{0, 1}, {0, 1}}};
        const std::array<double, 2> errors{0.0031, 0.1004};
        ok &= check(calibrated_repeat_insertion_bridge_flip(gauge, reverse, errors) == true,
                    "repeat insertion: diploid physical witnesses calibrate a reverse union");
        ok &= check(!calibrated_repeat_insertion_bridge_flip(gauge, {{{0, 1}, {0, 0}}}, errors),
                    "repeat insertion: one class cannot replace the exact-ALT certificate");
        ok &= check(!calibrated_repeat_insertion_bridge_flip(gauge, {{{0, 1}, {1, 0}}}, errors),
                    "repeat insertion: opposing class gauges veto the union");
        ok &= check(!calibrated_repeat_insertion_bridge_flip(gauge, {{{1, 1}, {0, 1}}}, errors),
                    "repeat insertion: any contrary independent molecule vetoes the union");
        gauge.counts[1][1] = 3;
        ok &= check(!calibrated_repeat_insertion_bridge_flip(gauge, reverse, errors),
                    "repeat insertion: insufficient calibration abstains");
        gauge.counts = {{{4, 0}, {0, 5}}};
        ok &= check(!calibrated_repeat_insertion_bridge_flip(gauge, reverse, {0.2, 0.2}),
                    "repeat insertion: joint error above twenty percent abstains");
        gauge.counts = {{{0, 4}, {5, 0}}};
        ok &= check(!calibrated_repeat_insertion_bridge_flip(gauge, reverse, errors),
                    "repeat insertion: a reversed upstream gauge remains a veto");
    }

    {
        IndependentBamBlockLink left, right;
        left.counts = {{{27, 0}, {3, 14}}};
        right.counts = {{{15, 1}, {1, 4}}};
        const std::array<int, 2> parity{0, 4};
        const std::array<std::array<int, 2>, 2> classes{{{0, 1}, {0, 3}}};
        const std::array<double, 2> errors{0.000502, 1e-9};
        ok &= check(calibrated_indel_bridge_flip(left, parity, 28.6298970081) == true &&
                    calibrated_repeat_insertion_bridge_flip(right, classes, errors) == true,
                    "compound prefix: separate diploid gauges certify the physical bridge");
        ok &= check(!calibrated_repeat_insertion_bridge_flip(right, {{{0, 0}, {0, 3}}}, errors),
                    "compound prefix: shared-prefix quality cannot replace a missing allele class");
        ok &= check(!calibrated_repeat_insertion_bridge_flip(right, {{{0, 1}, {1, 3}}}, errors),
                    "compound prefix: contrary physical molecules veto the union");
    }

    {
        ok &= check(graph_snp_low_alt_fraction_supported(59, 8, 0, 0, 0.2, 1000),
                    "physical minor repeat SNP fails the diploid balance test");
        ok &= check(!graph_snp_low_alt_fraction_supported(59, 20, 0, 0, 0.2, 1000) &&
                    !graph_snp_low_alt_fraction_supported(6, 1, 0, 0, 0.2, 1000),
                    "valid ALT fraction and weak evidence preserve graph SNPs");
        ok &= check(!graph_snp_low_alt_fraction_supported(59, 8, 1, 0, 0.2, 1000) &&
                    !graph_snp_low_alt_fraction_supported(59, 8, 0, 1, 0.2, 1000) &&
                    !graph_snp_low_alt_fraction_supported(59, 8, 0, 0, 0.2, 0),
                    "uncallable alternate alleles and missing family size veto rejection");
        ok &= check(graph_snp_padded_deletion_supported(0, 28, 9, 0, 1000),
                    "padded SNP: significant physical ALT/deletion contrast survives normalization");
        ok &= check(!graph_snp_padded_deletion_supported(1, 28, 9, 0, 1000) &&
                    !graph_snp_padded_deletion_supported(0, 28, 9, 1, 1000) &&
                    !graph_snp_padded_deletion_supported(0, 28, 9, 0, 0),
                    "padded SNP: REF, third base or absent family vetoes contradiction");
        ok &= check(!graph_snp_padded_deletion_supported(0, 28, 1, 0, 1000) &&
                    !graph_snp_padded_deletion_supported(0, 40, 9, 0, 1000) &&
                    !graph_snp_padded_deletion_supported(0, 10, 9, 0, 1000),
                    "padded SNP: singleton, incidental deletion and weak significance abstain");
        ok &= check(physical_deletion_gauge_haplotype({2, 0}, 0.0001) == 1 &&
                    physical_deletion_gauge_haplotype({0, 3}, 0.001) == 2,
                    "deletion gauge: independent pairs orient either haplotype");
        ok &= check(!physical_deletion_gauge_haplotype({2, 1}, 0.0001) &&
                    !physical_deletion_gauge_haplotype({1, 0}, 0.0001) &&
                    !physical_deletion_gauge_haplotype({0, 0}, 0.0) &&
                    !physical_deletion_gauge_haplotype({2, 0}, 0.0011) &&
                    !physical_deletion_gauge_haplotype({2, 0}, -1.0),
                    "deletion gauge: contrary, weak, missing and invalid evidence abstains");
        ok &= check(graph_snp_ref_absence_supported(0, 54, 0, 0, 1000),
                    "REF absence: homozygous ALT does not require a deletion");
        ok &= check(!graph_snp_ref_absence_supported(1, 54, 0, 0, 1000) &&
                    !graph_snp_ref_absence_supported(0, 54, 0, 1, 1000),
                    "REF absence: any REF or third base vetoes exclusion");
        ok &= check(!graph_snp_ref_absence_supported(0, 10, 0, 0, 1000) &&
                    !graph_snp_ref_absence_supported(0, 0, 0, 0, 1000) &&
                    !graph_snp_ref_absence_supported(0, 54, 0, 0, 0),
                    "REF absence: weak, missing and untested evidence abstains");
        ok &= check(graph_snp_ref_absence_supported(0, 40, 10, 0, 1000) &&
                    !graph_snp_ref_absence_supported(0, 40, 9, 0, 1000) &&
                    !graph_snp_ref_absence_supported(0, 41, 10, 0, 1000),
                    "REF absence: ALT/deletion keeps its count and fraction gates");
    }

    {
        const GraphSnpReferenceEvidence terminal{1, 0, 4, 1.0};
        GraphSnpReferenceEvidence partner{1, 4, 0, 0.000032};
        ok &= check(graph_snp_cohort_is_physically_contradicted(terminal, partner),
                    "graph SNP cohort: deletion and reference-only evidence excludes false anchors");
        ok &= check(!graph_snp_cohort_is_physically_contradicted(std::nullopt, partner) &&
                    !graph_snp_cohort_is_physically_contradicted(terminal, std::nullopt),
                    "graph SNP cohort: credible physical ALT or missing evidence vetoes exclusion");
        auto weak = terminal;
        weak.alternate_class_deletions = 1;
        ok &= check(!graph_snp_cohort_is_physically_contradicted(weak, partner),
                    "graph SNP cohort: a single deleted molecule cannot exclude anchors");
        weak = terminal;
        weak.reference_class_reads = 0;
        ok &= check(!graph_snp_cohort_is_physically_contradicted(weak, partner),
                    "graph SNP cohort: an unverified reference class cannot exclude anchors");
        partner.alternate_class_reference_reads = 1;
        ok &= check(!graph_snp_cohort_is_physically_contradicted(terminal, partner),
                    "graph SNP cohort: both sites need independent contrary molecules");
        partner.alternate_class_reference_reads = 4;
        partner.wrong_alternate_bound = 0.0011;
        ok &= check(!graph_snp_cohort_is_physically_contradicted(terminal, partner),
                    "graph SNP cohort: weak mapping or base evidence cannot exclude anchors");
        partner.wrong_alternate_bound = 0.000032;
        partner.reference_class_reads = 0;
        ok &= check(!graph_snp_cohort_is_physically_contradicted(terminal, partner),
                    "graph SNP cohort: partner reference class also needs a physical witness");
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
