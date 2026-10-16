#include "graph_sites.hpp"

#include <cstdio>
#include <fstream>
#include <iostream>
#include <set>
#include <tuple>
#include <sstream>
#include <string>

using namespace pgphase_collect;

static bool check(bool cond, const std::string& msg) {
    if (!cond) std::cerr << "FAIL: " << msg << "\n";
    return cond;
}

static const GraphSite* find_site(const GraphSiteCatalog& catalog, const std::string& id) {
    for (const GraphSite& site : catalog.sites) {
        if (site.id == id) return &site;
    }
    return nullptr;
}

static bool test_loader_rejects_bad_records_without_throwing() {
    // A catalog is external data: one corrupt record must cost that record, not
    // the run, and every discard must be counted rather than silent.
    const char* lines[] = {
        "chr20\t100\ts1\tA\tG\t60\t.\tAT=>1>2,>1>3\n",          // fine
        "chr20\tNOTANUMBER\ts2\tA\tG\t60\t.\tAT=>1>2,>1>3\n",   // POS not a number
        "chr20\t0\ts3\tA\tG\t60\t.\tAT=>1>2,>1>3\n",            // POS not positive
        "chr20\t300\ts4\tA\t<DEL>\t60\t.\tAT=>1>2,>1>3\n",      // symbolic allele
        "chr20\t400\ts5\ta\tg\t60\t.\tAT=>1>2,>1>3\n",          // lower case
        "chr20\t500\ts6\tA\n",                                     // short line
        "chr20\t600\ts7\tA\tG\t60\t.\tAT=>1>2,>1>3;END=NOPE\n", // END not a number
    };
    const std::string path = "/tmp/pgphase_loader_test.vcf";
    {
        std::ofstream out(path);
        out << "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n";
        for (const char* l : lines) out << l;
    }
    const GraphSiteCatalog cat = load_graph_site_catalog_from_vcf(path, {}, true);
    bool ok = true;
    ok &= check(cat.stats.data_lines == 7, "every data line is counted");
    ok &= check(cat.stats.sites_parsed == 2, "the two well-formed records survive parsing");
    ok &= check(cat.stats.bad_position == 3, "two bad POS and one bad END are counted");
    ok &= check(cat.stats.unsupported_allele == 1, "the symbolic allele is rejected");
    ok &= check(cat.stats.short_line == 1, "the short line is counted");
    // Case is normalised rather than carried through, so a lower-case allele
    // still matches a candidate derived from an upper-case reference.
    bool saw_upper = false;
    for (const GraphSite& s : cat.sites)
        if (s.pos == 400) saw_upper = (s.ref == "A" && s.alts.size() == 1 && s.alts[0] == "G");
    ok &= check(saw_upper, "a lower-case allele is upper-cased at load");
    return ok;
}

using SiteKey = std::tuple<hts_pos_t, std::string, std::string, std::set<std::pair<std::string, std::string>>, int,
                           std::string>;

// Sites as a set: ALT order is the builder's choice (vg's, or pgcat's by
// haplotype count), and must not change what is loaded.
static std::set<SiteKey> site_set(const GraphSiteCatalog& cat) {
    std::set<SiteKey> out;
    for (const GraphSite& s : cat.sites) {
        std::set<std::pair<std::string, std::string>> alleles;
        for (size_t i = 0; i < s.alts.size(); ++i)
            alleles.emplace(s.alts[i], i + 1 < s.allele_traversals.size() ? s.allele_traversals[i + 1] : "");
        out.emplace(s.pos, s.id, s.ref, alleles, s.level, s.parent);
    }
    return out;
}

// The micb catalog as pgcat writes it (a folder of <contig>.pgsites, or one
// file) and as vg deconstruct does (a tabix VCF): the same sites for whole
// contigs and for regions, records overlapping a region included as tabix does.
static bool test_pgsites_catalog_matches_vcf() {
    const std::string dir = "test_data/graph/micb-kir3dl1.catalog";
    const std::string vcf = "test_data/graph/micb-kir3dl1.sites.vcf.gz";
    bool ok = true;
    require_site_catalog(dir);
    require_site_catalog(dir + "/chr6.pgsites");
    SiteCatalogHandle folder(dir), file(dir + "/chr6.pgsites"), tabix(vcf);
    struct Query { const char* contig; hts_pos_t beg, end; size_t expect; };
    const Query queries[] = {
        {"CHM13#0#chr6", 1, 400000000, 202},
        {"CHM13#0#chr19", 1, 400000000, 715},
        {"CHM13#0#chr6", 31351180, 31351181, 0},  // inside the >16>32 record
        {"CHM13#0#chr6", 31355000, 31357000, 0},
        {"chr19", 18930000, 18935000, 0},         // bare contig name
    };
    for (const Query& q : queries) {
        const GraphSiteCatalog a = load_sites_for_region(folder, q.contig, q.beg, q.end);
        const GraphSiteCatalog b = load_sites_for_region(tabix, q.contig, q.beg, q.end);
        const std::string what = std::string(q.contig) + ":" + std::to_string(q.beg) + "-" + std::to_string(q.end);
        ok &= check(!a.sites.empty(), what + ": sites loaded from the folder");
        ok &= check(q.expect == 0 || a.sites.size() == q.expect, what + ": every site of the contig");
        ok &= check(site_set(a) == site_set(b), what + ": the same sites as the VCF");
        if (std::string(q.contig).find("chr6") != std::string::npos)
            ok &= check(site_set(load_sites_for_region(file, q.contig, q.beg, q.end)) == site_set(a),
                        what + ": the single file loads as the folder does");
    }
    ok &= check(load_sites_for_region(file, "CHM13#0#chr19", 1, 400000000).sites.empty(),
                "a single contig's file has no sites on another contig");
    ok &= check(load_sites_for_region(folder, "CHM13#0#chrM", 1, 400000000).sites.empty(),
                "a contig without a file has no sites");
    ok &= check(load_sites_for_region(folder, "GRCh38#0#chr6", 1, 400000000).sites.empty() &&
                    load_sites_for_region(tabix, "GRCh38#0#chr6", 1, 400000000).sites.empty(),
                "another sample's path has no sites in a CHM13 catalog");
    return ok;
}

int main() {
    const std::string path = "/tmp/pgphase_graph_sites_test.vcf";
    {
        std::ofstream out(path);
        out << "##fileformat=VCFv4.2\n";
        out << "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n";
        out << "chr20\t100\tparent\tA\tC,*\t.\tPASS\tLV=0;AT=>1>2>4,>1<3>4;X=Y\n";
        out << "chr20\t110\tchild\tG\tT\t.\tPASS\tLV=1;PS=parent;AT=>2>4>6,>2>5>6\n";
        out << "chrNested\t120\tnested\tG\tA,*\t.\tPASS\tLV=1;PS=parent;PA=1;RC=chr20;END=130;AT=>6>7>9,>6>8>9,*\n";
        out << "chr20\t140\tbad\tA\tC\t.\tPASS\tLV=0;AT=>10>11,bad_walk\n";
    }

    bool ok = true;
    GraphSiteCatalog catalog = load_graph_site_catalog_from_vcf(path);
    const GraphSite* parent = find_site(catalog, "parent");
    const GraphSite* child = find_site(catalog, "child");
    const GraphSite* nested = find_site(catalog, "nested");
    const GraphSite* bad = find_site(catalog, "bad");
    ok &= check(catalog.sites.size() == 4, "loaded four graph sites");
    ok &= check(parent != nullptr, "parent site found");
    ok &= check(child != nullptr, "child site found");
    ok &= check(nested != nullptr, "nested site found");
    ok &= check(bad != nullptr, "malformed site found");
    ok &= check(parent != nullptr && parent->level == 0, "parent LV parsed");
    ok &= check(parent != nullptr && parent->parent.empty(), "parent has no PS");
    ok &= check(parent != nullptr && parent->has_spanning_deletion, "star allele detected");
    ok &= check(parent != nullptr && parent->allele_traversals.size() == 2, "parent AT traversals parsed");
    ok &= check(parent != nullptr && parent->allele_walks.size() == 2, "parent AT walks parsed");
    ok &= check(parent != nullptr && parent->allele_walks[1].size() == 3, "parent second walk length");
    ok &= check(parent != nullptr && parent->allele_walks[1][1].node == 3, "walk node id parsed");
    ok &= check(parent != nullptr && parent->allele_walks[1][1].reverse, "walk orientation parsed");
    ok &= check(parent != nullptr && graph_site_between_query(*parent) == "1+:4+",
                "between query from boundary nodes");
    ok &= check(child != nullptr && child->level == 1, "child LV parsed");
    ok &= check(child != nullptr && child->parent == "parent", "child PS parsed");
    ok &= check(child != nullptr && child->allele_traversals[1] == ">2>5>6", "child AT traversal parsed");
    ok &= check(nested != nullptr && nested->ref_contig == "chr20", "nested RC parsed");
    ok &= check(nested != nullptr && nested->ref_end == 130, "nested END parsed");
    ok &= check(nested != nullptr && nested->conditional_parent_alleles.size() == 1 &&
                nested->conditional_parent_alleles[0] == 1,
                "nested PA parsed");
    ok &= check(nested != nullptr && nested->allele_walks.size() == 3 &&
                nested->allele_walks[2].empty(),
                "star traversal kept as missing graph walk");
    ok &= check(bad != nullptr && !bad->eligible, "malformed traversal makes site ineligible");
    ok &= check(bad != nullptr && bad->skip_reason == "malformed_allele_traversal",
                "malformed traversal skip reason");

    std::ostringstream tsv;
    write_graph_site_catalog_tsv(tsv, catalog);
    ok &= check(tsv.str().find("HAS_STAR") != std::string::npos, "TSV header written");
    ok &= check(tsv.str().find("child") != std::string::npos, "TSV row written");

    const GraphWalk walk = parse_graph_walk(">11<12>13");
    ok &= check(walk.size() == 3, "standalone walk length");
    ok &= check(graph_walk_to_string(walk) == ">11<12>13", "walk roundtrip");
    ok &= check(graph_walk_to_string(reverse_graph_walk(walk)) == "<13>12<11",
                "walk reverse orientation");
    ok &= check(child != nullptr &&
                match_graph_allele_exact(parse_graph_walk(">2>5>6"),
                                         child->allele_walks) == 1,
                "exact allele match");
    ok &= check(child != nullptr &&
                match_graph_allele_exact(parse_graph_walk("<6<5<2"),
                                         child->allele_walks) == 1,
                "reverse exact allele match");
    ok &= check(child != nullptr &&
                match_graph_allele_exact(parse_graph_walk(">2>6>6"),
                                         child->allele_walks) == kGraphAlleleMissing,
                "missing allele match");

    const std::string gfa =
        "S\t1\tA\n"
        "S\t2\tC\n"
        "S\t3\tG\n"
        "S\t4\tT\n"
        "W\tREF\t0\tchr1\t0\t3\t>1>2>4\n"
        "W\tALT\t0\tchr1\t0\t3\t>1>3>4\n";
    GraphSiteCatalog auto_catalog = load_graph_site_catalog_from_gfa_text(gfa, "REF", "chr1");
    ok &= check(auto_catalog.sites.size() == 1, "automatic GFA catalog finds one anchored event");
    ok &= check(auto_catalog.sites[0].allele_walks.size() == 2, "automatic event has ref and alt walks");
    ok &= check(auto_catalog.sites[0].allele_traversals[0] == ">1>2>4",
                "automatic event ref traversal");
    ok &= check(auto_catalog.sites[0].allele_traversals[1] == ">1>3>4",
                "automatic event alt traversal");
    ok &= check(graph_site_between_query(auto_catalog.sites[0]) == "1+:4+",
                "automatic event between boundary");

    GraphSite incompatible;
    incompatible.id = "incompatible";
    incompatible.allele_traversals = {">1>2", ">3>4"};
    incompatible.allele_walks = {parse_graph_walk(">1>2"), parse_graph_walk(">3>4")};
    ok &= check(graph_site_validation_skip_reason(incompatible) == "incompatible_allele_boundaries",
                "incompatible boundaries are invalid");
    ok &= check(graph_site_between_query(incompatible).empty(),
                "incompatible boundaries have no between query");

    ok &= test_loader_rejects_bad_records_without_throwing();
    ok &= test_pgsites_catalog_matches_vcf();

    std::remove(path.c_str());
    if (ok) {
        std::cout << "ALL PASS\n";
        return 0;
    }
    return 1;
}
