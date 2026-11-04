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

// A catalog's sites as a set (ALT order is the builder's choice: vg's, or
// dromos's by haplotype count), hashed (FNV-1a): POS, ID, REF and its walk,
// and each ALT with its walk.
static uint64_t site_digest(const GraphSiteCatalog& cat) {
    std::set<std::tuple<hts_pos_t, std::string, std::string, std::set<std::pair<std::string, std::string>>>> sites;
    for (const GraphSite& s : cat.sites) {
        std::set<std::pair<std::string, std::string>> alleles;
        for (size_t i = 0; i < s.alts.size(); ++i)
            alleles.emplace(s.alts[i],
                            i + 1 < s.allele_walks.size() ? graph_walk_to_string(s.allele_walks[i + 1]) : "");
        sites.emplace(s.pos, s.id,
                      s.ref + "|" + graph_walk_to_string(s.allele_walks.empty() ? GraphWalk{} : s.allele_walks[0]),
                      alleles);
    }
    uint64_t h = 1469598103934665603ULL;
    const auto put = [&](const std::string& t) {
        for (const unsigned char c : t) {
            h ^= c;
            h *= 1099511628211ULL;
        }
        h ^= 0xff;
        h *= 1099511628211ULL;
    };
    for (const auto& [pos, id, ref, alleles] : sites) {
        put(std::to_string(pos));
        put(id);
        put(ref);
        for (const auto& [a, w] : alleles) {
            put(a);
            put(w);
        }
    }
    return h;
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
static const char* const kRegionCatalog = "test_data/graph_chr20/chr20_25M.CHM13.sites";

// The micb catalog as dromos catalog writes it (one file for both contigs):
// the sites vg deconstruct -a reported for whole contigs and for regions,
// records overlapping a region included as tabix did (the digests are of
// vg's sites VCF, read by the sites VCF reader pgphase had until dromos's
// catalog replaced it).
static bool test_site_catalog_sites() {
    bool ok = true;
    require_site_catalog(kCatalog);
    SiteCatalogHandle catalog(kCatalog);
    struct Query { const char* contig; hts_pos_t beg, end; size_t sites; uint64_t digest; };
    const Query queries[] = {
        {"CHM13#0#chr6", 1, 400000000, 202, 0xeb4c80d5cc53f3b8ULL},
        {"CHM13#0#chr19", 1, 400000000, 715, 0xd13c1413c839370dULL},
        {"CHM13#0#chr6", 31351180, 31351181, 1, 0x03c93e9f26c8ae80ULL},  // inside the >16>32 record
        {"CHM13#0#chr6", 31355000, 31357000, 27, 0xa84f33dd1118909fULL},
        {"chr19", 18930000, 18935000, 94, 0x2efe4d9cb837305fULL},         // bare contig name
    };
    for (const Query& q : queries) {
        const GraphSiteCatalog a = load_sites_for_region(catalog, q.contig, q.beg, q.end);
        const std::string what = std::string(q.contig) + ":" + std::to_string(q.beg) + "-" + std::to_string(q.end);
        ok &= check(a.sites.size() == q.sites, what + ": " + std::to_string(q.sites) + " sites (got " +
                                                   std::to_string(a.sites.size()) + ")");
        ok &= check(site_digest(a) == q.digest, what + ": the sites vg deconstruct reported");
    }
    ok &= check(load_sites_for_region(catalog, "CHM13#0#chrM", 1, 400000000).sites.empty(),
                "a contig the catalog lacks has no sites");
    ok &= check(load_sites_for_region(catalog, "GRCh38#0#chr6", 1, 400000000).sites.empty(),
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
    uint64_t lo = 0, hi = 0;
    ok &= check(f->node_range(lo, hi) && lo == 1 && hi == 2891, "the graph's node ids");
    ok &= check(f->property("graph.fingerprint") == "size=73920 fnv1a=9736178bc136e029" &&
                    f->property("reference.sample") == "CHM13" && f->property("no.such.key").empty(),
                "properties");
    const std::string line = f->vcf_header_line();
    ok &= check(line.rfind("##pgphase_site_catalog=<File=", 0) == 0 &&
                    line.find("Fingerprint=\"size=73920 fnv1a=9736178bc136e029\"") != std::string::npos &&
                    line.back() == '>',
                "the VCF header line names the graph");
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
    error.clear();
    try {
        SiteCatalogHandle h("test_data/graph/micb-kir3dl1_HG003.coord.gaf.gz");
    } catch (const std::exception& e) {
        error = e.what();
    }
    ok &= check(error.find("not dromos's site catalog") != std::string::npos,
                "another file is refused (got '" + error + "')");
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

// The catalog's contigs against a FASTA's: the same reference, one that runs
// past the graph's path, one with bases where the graph has none (N), and
// ones that are not the catalog's reference.
static bool test_site_catalog_reference_check() {
    bool ok = true;
    const std::shared_ptr<const SiteCatalogFile> f = SiteCatalogFile::open(kCatalog);
    const std::string ref = "test_data/graph/ref.fa.gz";
    for (const SiteCatalogReferenceCheck& c :
         check_site_catalog_reference(f, ref, {"CHM13#0#chr6", "CHM13#0#chr19"}, 2))
        ok &= check(c.md5_matches && !c.refuse && c.catalog_contig >= 0, c.fasta_contig + ": " + c.message);

    std::unique_ptr<faidx_t, void (*)(faidx_t*)> fai(fai_load(ref.c_str()), fai_destroy);
    hts_pos_t len = 0;
    char* raw = faidx_fetch_seq64(fai.get(), "CHM13#0#chr19", 0, 18941250, &len);
    const std::string seq = raw != nullptr ? std::string(raw, static_cast<size_t>(len)) : std::string();
    std::free(raw);
    if (!check(seq.size() == 18941251, "chr19 read from the test FASTA")) return false;
    SiteCatalogHandle h(kCatalog);
    const GraphSiteCatalog sites = load_sites_for_region(h, "chr19", 1, 18941251);
    if (!check(!sites.sites.empty(), "chr19's sites")) return false;
    const GraphSite& first = sites.sites.front();

    const auto checked = [&](const std::string& name, const std::string& contig_seq) {
        const std::string p = tmp_path("pgphase_test_" + name + ".fa");
        write_file(p, ">CHM13#0#chr19\n" + contig_seq + "\n");
        std::remove((p + ".fai").c_str());
        if (fai_build(p.c_str()) != 0) std::cerr << "cannot index " << p << "\n";
        const SiteCatalogReferenceCheck c = check_site_catalog_reference(f, p, {"CHM13#0#chr19"}, 1).at(0);
        std::remove(p.c_str());
        std::remove((p + ".fai").c_str());
        return c;
    };
    SiteCatalogReferenceCheck c = checked("longer", seq + "ACGT");
    ok &= check(c.md5_matches && !c.refuse && c.message.find("ends at") != std::string::npos,
                "a FASTA contig running past the graph's path: a warning (" + c.message + ")");
    std::string s = seq;
    ok &= check(s[100] == 'N', "the test FASTA's chr19 starts with N");
    s[100] = 'A';
    c = checked("gap", s);
    ok &= check(!c.md5_matches && !c.refuse && c.sampled_records > 0 && c.sampled_agree == c.sampled_records,
                "bases where the graph has none: a warning (" + c.message + ")");
    s = seq;
    const size_t at = static_cast<size_t>(first.pos - 1);
    s[at] = s[at] == 'A' ? 'C' : 'A';
    c = checked("ref", s);
    ok &= check(!c.md5_matches && c.refuse && c.sampled_agree < c.sampled_records,
                "a REF the FASTA does not spell: refused (" + c.message + ")");
    c = checked("short", seq.substr(0, 1000));
    ok &= check(c.refuse, "a FASTA contig shorter than the graph's path: refused (" + c.message + ")");

    // A catalog of a region (dromos catalog --region): chr20 of the HPRC graph
    // (66,210,255 bp) against a FASTA that ends with the region, at 25.5 Mb.
    const std::shared_ptr<const SiteCatalogFile> slice = SiteCatalogFile::open(kRegionCatalog);
    ok &= check(slice->contigs().size() == 1 && slice->contigs()[0].regions.size() == 1 &&
                    slice->contigs()[0].regions[0] == std::make_pair<uint64_t, uint64_t>(25000001, 25500000) &&
                    slice->contigs()[0].records == 6963,
                "the region catalog's region and records");
    c = check_site_catalog_reference(slice, "test_data/graph_chr20/ref.fa.gz", {"chr20"}, 1).at(0);
    ok &= check(!c.refuse && c.sampled_records > 0 && c.sampled_agree == c.sampled_records &&
                    c.message.find("holding the catalog's regions") != std::string::npos,
                "a FASTA contig reaching the catalog's regions: a warning (" + c.message + ")");
    return ok;
}

int main() {
    bool ok = true;
    GraphSite parent;
    parent.id = "parent";
    parent.alts = {"C", "*"};
    parent.allele_walks = {parse_graph_walk(">1>2>4"), parse_graph_walk(">1<3>4")};
    parent.skip_reason = graph_site_validation_skip_reason(parent);
    parent.eligible = parent.skip_reason.empty();
    ok &= check(parent.eligible, "a site with two walks between the same boundaries is eligible");
    ok &= check(parent.allele_walks[1].size() == 3, "walk length");
    ok &= check(parent.allele_walks[1][1].node == 3, "walk node id parsed");
    ok &= check(parent.allele_walks[1][1].reverse, "walk orientation parsed");

    const GraphWalk walk = parse_graph_walk(">11<12>13");
    ok &= check(walk.size() == 3, "standalone walk length");
    ok &= check(graph_walk_to_string(walk) == ">11<12>13", "walk roundtrip");

    GraphSite incompatible;
    incompatible.id = "incompatible";
    incompatible.allele_walks = {parse_graph_walk(">1>2"), parse_graph_walk(">3>4")};
    ok &= check(graph_site_validation_skip_reason(incompatible) == "incompatible_allele_boundaries",
                "incompatible boundaries are invalid");

    ok &= test_site_catalog_sites();
    ok &= test_site_catalog_header();
    ok &= test_site_catalog_refuses_damage();
    ok &= test_site_catalog_index_search();
    ok &= test_site_catalog_reference_check();

    if (ok) {
        std::cout << "ALL PASS\n";
        return 0;
    }
    return 1;
}
