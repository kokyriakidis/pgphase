#include "graph_sites.hpp"
#include "site_catalog.hpp"

#include <cstdio>
#include <cstdlib>
#include <random>
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

// Sites as a set: ALT order is the builder's choice (vg's, or pgcatalog's by
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

static std::string tmp_path(const std::string& name) {
    const char* dir = std::getenv("TMPDIR");
    return std::string(dir != nullptr && *dir != '\0' ? dir : "/tmp") + "/" + name;
}

static std::string read_file(const std::string& path) {
    std::ifstream in(path, std::ios::binary);
    std::ostringstream s;
    s << in.rdbuf();
    return s.str();
}

static void write_file(const std::string& path, const std::string& data) {
    std::ofstream out(path, std::ios::binary);
    out << data;
}

static const char* const kCatalog = "test_data/graph/micb-kir3dl1.CHM13.sites";

// The micb catalog as dromos catalog writes it (one file for both contigs) and
// as vg deconstruct does (a tabix VCF): the same sites for whole contigs and
// for regions, records overlapping a region included as tabix does.
static bool test_site_catalog_matches_vcf() {
    const std::string vcf = "test_data/graph/micb-kir3dl1.sites.vcf.gz";
    bool ok = true;
    require_site_catalog(kCatalog);
    SiteCatalogHandle catalog(kCatalog), tabix(vcf);
    struct Query { const char* contig; hts_pos_t beg, end; size_t expect; };
    const Query queries[] = {
        {"CHM13#0#chr6", 1, 400000000, 202},
        {"CHM13#0#chr19", 1, 400000000, 715},
        {"CHM13#0#chr6", 31351180, 31351181, 0},  // inside the >16>32 record
        {"CHM13#0#chr6", 31355000, 31357000, 0},
        {"chr19", 18930000, 18935000, 0},         // bare contig name
    };
    for (const Query& q : queries) {
        const GraphSiteCatalog a = load_sites_for_region(catalog, q.contig, q.beg, q.end);
        const GraphSiteCatalog b = load_sites_for_region(tabix, q.contig, q.beg, q.end);
        const std::string what = std::string(q.contig) + ":" + std::to_string(q.beg) + "-" + std::to_string(q.end);
        ok &= check(!a.sites.empty(), what + ": sites loaded from the catalog");
        ok &= check(q.expect == 0 || a.sites.size() == q.expect, what + ": every site of the contig");
        ok &= check(site_set(a) == site_set(b), what + ": the same sites as the VCF");
    }
    ok &= check(load_sites_for_region(catalog, "CHM13#0#chrM", 1, 400000000).sites.empty(),
                "a contig the catalog lacks has no sites");
    ok &= check(load_sites_for_region(catalog, "GRCh38#0#chr6", 1, 400000000).sites.empty() &&
                    load_sites_for_region(tabix, "GRCh38#0#chr6", 1, 400000000).sites.empty(),
                "another sample's path has no sites in a CHM13 catalog");
    // Every thread's handle reads the one file opened for the process.
    SiteCatalogHandle second(kCatalog);
    ok &= check(&catalog.catalog->file() == &second.catalog->file(), "handles share the opened catalog");
    return ok;
}

// The header: the graph's identity and each contig's.
static bool test_site_catalog_header() {
    bool ok = true;
    const std::shared_ptr<const SiteCatalogFile> f = SiteCatalogFile::open(kCatalog);
    const std::vector<SiteCatalogContig>& c = f->contigs();
    ok &= check(c.size() == 2 && c[0].name == "chr6" && c[1].name == "chr19", "contigs in natural order");
    ok &= check(c.size() == 2 && c[0].chrom == "CHM13#0#chr6" && c[0].length == 31363898 &&
                    c[0].records == 202 && c[0].md5_hex() == "a5b49a5f3e75bedbb6ce1a71871a525c",
                "chr6's CHROM, length, records and MD5");
    ok &= check(c.size() == 2 && c[1].length == 18941251 && c[1].records == 715 &&
                    c[1].md5_hex() == "987a61a2c458a9ed426d5257295cb293",
                "chr19's length, records and MD5");
    ok &= check(f->find_contig("chr19") == 1 && f->find_contig("CHM13#0#chr19") == 1 &&
                    f->find_contig("GRCh38#0#chr19") == -1 && f->find_contig("chr20") == -1,
                "contigs found by name or by the reference sample's path name");
    ok &= check(f->property("graph.fingerprint") == "size=73920 fnv1a=9736178bc136e029" &&
                    f->property("reference.sample") == "CHM13" && f->property("no.such.key").empty(),
                "properties");
    return ok;
}

// A file that is not a version 1 catalog, or is damaged, is refused, saying why.
static bool test_site_catalog_refuses_damage() {
    bool ok = true;
    const std::string good = read_file(kCatalog);
    if (!check(good.size() > 200, "the test catalog is readable")) return false;
    const auto refused = [&](const std::string& name, const std::string& data, const std::string& expect,
                             bool on_query) {
        const std::string p = tmp_path("pgphase_test_" + name + ".sites");
        write_file(p, data);
        std::string error;
        try {
            SiteCatalogHandle h(p);
            if (on_query) {
                load_sites_for_region(h, "chr6", 1, 400000000);
                load_sites_for_region(h, "chr19", 1, 400000000);
            }
        } catch (const std::exception& e) {
            error = e.what();
        }
        std::remove(p.c_str());
        return check(error.find(expect) != std::string::npos,
                     name + ": refused with '" + expect + "' (got '" + error + "')");
    };
    const uint64_t header_at = [&]() {
        uint64_t v = 0;
        for (int i = 7; i >= 0; --i) v = (v << 8) | static_cast<unsigned char>(good[24 + i]);
        return v;
    }();
    std::string d = good;
    d[8] = 2;
    ok &= refused("version", d, "format version 2", false);
    d = good;
    d[17] ^= 1;
    ok &= refused("preamble", d, "preamble's hash", false);
    ok &= refused("truncated", good.substr(0, good.size() - 10), "truncated", false);
    d = good;
    d[header_at + 20] ^= 0x40;
    ok &= refused("header", d, "its header", false);
    d = good;
    d[64] ^= 1;  // the first block's first byte
    ok &= refused("block", d, "CRC differs", true);
    d = good;
    d.replace(0, 8, "PGSITES2");
    ok &= refused("pgsites", d, "pgcatalog .pgsites", false);
    std::string error;
    try {
        require_site_catalog("test_data/graph");
    } catch (const std::exception& e) {
        error = e.what();
    }
    ok &= check(error.find("a folder") != std::string::npos, "a folder is refused (got '" + error + "')");
    return ok;
}

// The index's binary searches pick the blocks a linear scan does.
static bool test_site_catalog_index_search() {
    std::mt19937_64 rng(7);
    bool ok = true;
    for (int trial = 0; trial < 50; ++trial) {
        SiteCatalogIndex ix;
        uint64_t pos = 1 + rng() % 50, reach = 0;
        const size_t n = 1 + rng() % 40;
        for (size_t b = 0; b < n; ++b) {
            ix.first_pos.push_back(pos);
            // Mostly short records; now and then one reaching far ahead.
            const uint64_t end = pos + (rng() % 8 == 0 ? rng() % 2000 : rng() % 60);
            ix.max_end.push_back(end);
            reach = std::max(reach, end);
            ix.reach.push_back(reach);
            pos += rng() % 3 == 0 ? 0 : rng() % 100;
        }
        for (int q = 0; q < 200; ++q) {
            const uint64_t beg = 1 + rng() % (pos + 100), end = beg + rng() % 300;
            std::vector<size_t> linear, searched;
            for (size_t b = 0; b < n; ++b)
                if (ix.first_pos[b] <= end && ix.max_end[b] >= beg) linear.push_back(b);
            const auto [first, last] = ix.overlapping(beg, end);
            for (size_t b = first; b < last; ++b)
                if (ix.max_end[b] >= beg) searched.push_back(b);
            ok &= check(linear == searched, "binary search on trial " + std::to_string(trial));
            if (!ok) return false;
        }
    }
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
    ok &= test_site_catalog_matches_vcf();
    ok &= test_site_catalog_header();
    ok &= test_site_catalog_refuses_damage();
    ok &= test_site_catalog_index_search();

    std::remove(path.c_str());
    if (ok) {
        std::cout << "ALL PASS\n";
        return 0;
    }
    return 1;
}
