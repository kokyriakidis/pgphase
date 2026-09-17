// Do the alignment channel's sites reach the hybrid, and arrive intact?
//
// Three separate questions, three separate test cases, because they fail for
// different reasons and a single "injection works" assertion would hide which:
//
//   1. COMPLETENESS   -- is every site the alignment channel finds present in
//                        the hybrid's candidate table? A site the hybrid never
//                        receives cannot be phased however good the solver is.
//   2. REPRESENTATION -- is each shared site described by the same alleles?
//                        A locus whose two haplotypes are both non-reference
//                        needs one record carrying both alternates; replacing
//                        it with a single-allele record from the graph catalog
//                        collapses its depth and drives its allele fraction to
//                        1, which reads as homozygous. Measured at
//                        chr20:55,896,396, where the merged REF=T(16) record at
//                        DP 69 was displaced by a claim's version at DP 37,
//                        2/35, AF 0.946, classified CleanHom -- a real
//                        heterozygote at purity 1.000 over 28 reads.
//   3. COUNTS         -- are the per-site read counts right, both against the
//                        alignment channel's own numbers and internally? The
//                        strand tallies are load-bearing rather than
//                        informational: the ONT strand-bias screen in
//                        classify_graph_only_candidates computes an expected
//                        value from forward_alt + reverse_alt and declines to
//                        test when it is not positive, so tallies left at zero
//                        silently exempt every injected candidate from a screen
//                        every alignment-derived candidate faces.
//
// Cheap by construction: both channels are run once per window and cached, so
// adding a test case costs assertions rather than pipeline runs.

#define CATCH_CONFIG_MAIN
#include "../third_party/catch2/catch.hpp"

#include <htslib/faidx.h>
#include <htslib/tbx.h>
#include <htslib/sam.h>

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <map>
#include <set>
#include <sstream>
#include <string>
#include <tuple>
#include <vector>

namespace {

std::string env_or(const char* key, const std::string& fallback) {
    const char* v = std::getenv(key);
    return (v != nullptr && *v != '\0') ? std::string(v) : fallback;
}

/// Threads for every pipeline run these tests make, never fewer than four.
///
/// Four is a floor rather than a default because below it these tests get
/// slower for no reason, and above it they gain nothing. Measured on
/// chr20:30,000,000-35,000,000, varying only -t: 161.6 s at 1 thread, 46.7 s at
/// 4 (3.46x), 43.6 s at 12 and 44.0 s at 20 -- the ceiling is ~4.2x because
/// collect_pipeline.cpp caps workers at the chunks in a batch and joins between
/// batches, and chunk costs are uneven. Total CPU work is flat across all four
/// and the VCF is byte-identical, so this affects only wall time.
///
/// PGPHASE_TEST_THREADS can raise it; a lower value is clamped up rather than
/// honoured, so a run of these tests cannot be made accidentally serial.
int test_threads() {
    constexpr int kFloor = 4;
    const char* v = std::getenv("PGPHASE_TEST_THREADS");
    if (v == nullptr || *v == '\0') return kFloor;
    int requested = 0;
    try {
        requested = std::stoi(v);
    } catch (const std::exception&) {
        return kFloor;
    }
    return std::max(kFloor, requested);
}

bool file_exists(const std::string& path) {
    std::ifstream in(path);
    return in.good();
}

std::vector<std::string> split_char(const std::string& s, char sep) {
    std::vector<std::string> out;
    std::string field;
    std::istringstream in(s);
    while (std::getline(in, field, sep)) out.push_back(field);
    return out;
}

std::vector<std::string> split_tabs(const std::string& line) {
    std::vector<std::string> out;
    std::string field;
    std::istringstream in(line);
    while (std::getline(in, field, '\t')) {
        while (!field.empty() && (field.back() == '\r' || field.back() == '\n')) field.pop_back();
        out.push_back(field);
    }
    return out;
}

/// One row of a candidate table, by name rather than by column index: the
/// column order is not part of any contract and has changed before.
struct Candidate {
    std::map<std::string, std::string> f;
    const std::string& get(const std::string& k) const {
        static const std::string empty;
        const auto it = f.find(k);
        return it == f.end() ? empty : it->second;
    }
    long long num(const std::string& k) const {
        const std::string& v = get(k);
        return v.empty() || v == "." ? 0 : std::stoll(v);
    }
    double real(const std::string& k) const {
        const std::string& v = get(k);
        return v.empty() || v == "." ? 0.0 : std::stod(v);
    }
};

/// A site's identity: position, event type, and the exact alleles. The alleles
/// are part of the key on purpose -- two records at one position describing
/// different events are different sites, and collapsing them by position is how
/// a multiallelic locus gets silently reduced to one of its alternates.
using Key = std::tuple<long long, std::string, std::string, std::string>;

Key key_of(const Candidate& c) {
    return {c.num("POS"), c.get("TYPE"), c.get("REF"), c.get("ALT")};
}

std::string show(const Key& k) {
    std::ostringstream o;
    o << std::get<0>(k) << " " << std::get<1>(k) << " " << std::get<2>(k).substr(0, 24)
      << ">" << std::get<3>(k).substr(0, 24);
    return o.str();
}

std::map<Key, Candidate> load_candidates(const std::string& path) {
    std::map<Key, Candidate> out;
    std::ifstream in(path);
    std::string line;
    if (!std::getline(in, line)) return out;
    const auto header = split_tabs(line);
    while (std::getline(in, line)) {
        if (line.empty()) continue;
        const auto vals = split_tabs(line);
        Candidate c;
        for (size_t i = 0; i < header.size() && i < vals.size(); ++i) c.f[header[i]] = vals[i];
        out[key_of(c)] = c;
    }
    return out;
}

/// One emitted VCF record, reduced to what these tests judge: where it is, what
/// alleles it names, the genotype, and the per-allele depths.
struct Record {
    long long pos = 0;
    std::string ref, alt, gt;
    std::vector<int> ad;
    bool multiallelic() const { return alt.find(',') != std::string::npos; }
    /// Both genotype calls are the same allele.
    bool homozygous_call() const {
        const size_t bar = gt.find_first_of("|/");
        return bar != std::string::npos && gt.substr(0, bar) == gt.substr(bar + 1);
    }
    /// Reads behind two or more ALTs, which no homozygous call can describe.
    int alts_with_reads() const {
        int n = 0;
        for (size_t i = 1; i < ad.size(); ++i) if (ad[i] > 0) ++n;
        return n;
    }
};

std::vector<Record> load_records(const std::string& path) {
    std::vector<Record> out;
    std::ifstream in(path);
    std::string line;
    while (std::getline(in, line)) {
        if (line.empty() || line[0] == '#') continue;
        const auto f = split_tabs(line);
        if (f.size() < 10) continue;
        Record r;
        r.pos = std::stoll(f[1]);
        r.ref = f[3];
        r.alt = f[4];
        const auto sample = split_char(f[9], ':');
        if (sample.empty()) continue;
        r.gt = sample[0];
        // AD is the third sub-field of the sample column in this writer's layout
        // (GT:DP:AD:...); parsed positionally because the FORMAT string is fixed.
        if (sample.size() > 2)
            for (const std::string& v : split_char(sample[2], ','))
                r.ad.push_back(v.empty() || !isdigit(static_cast<unsigned char>(v[0])) ? 0
                                                                                       : std::stoi(v));
        out.push_back(std::move(r));
    }
    return out;
}

struct Window { long long gap_left = 0, gap_right = 0, gap_bp = 0; };

/// Known, documented defects these tests find and that are not yet fixed.
///
/// A test that is red for a known reason gates nothing, and a test that simply
/// does not look would not have found these. So the known cases are listed and
/// the assertion is "no NEW ones": a duplicate at a position not in this file
/// fails, and a strandless count above the recorded ceiling fails.
/// evaluations/2026-09-17-injection-tests/README.md carries the mechanism for
/// each. Removing a row here is how a fix gets locked in.
struct Allowance {
    std::set<std::pair<long long, long long>> duplicates;   // (window, position)
    std::set<std::pair<long long, long long>> overdepth;    // (window, position)
    std::set<std::pair<long long, long long>> genotype;     // (window, position)
    std::map<long long, int> emitted_multiallelic;          // window -> floor
    std::set<std::pair<long long, long long>> claims;       // (window, position)
    std::set<std::pair<long long, long long>> ref_convention;  // (window, position)
    std::map<long long, int> strandless_ceiling;            // window -> count
};

Allowance load_allowance(const std::string& path) {
    Allowance a;
    std::ifstream in(path);
    std::string line;
    while (std::getline(in, line)) {
        if (line.empty() || line[0] == '#') continue;
        const auto f = split_tabs(line);
        if (f.size() < 3 || f[0] == "kind") continue;
        if (f[0] == "duplicate") a.duplicates.insert({std::stoll(f[1]), std::stoll(f[2])});
        else if (f[0] == "depth") a.overdepth.insert({std::stoll(f[1]), std::stoll(f[2])});
        else if (f[0] == "genotype") a.genotype.insert({std::stoll(f[1]), std::stoll(f[2])});
        else if (f[0] == "emitted_multi") a.emitted_multiallelic[std::stoll(f[1])] = std::stoi(f[2]);
        else if (f[0] == "claim") a.claims.insert({std::stoll(f[1]), std::stoll(f[2])});
        else if (f[0] == "ref_convention")
            a.ref_convention.insert({std::stoll(f[1]), std::stoll(f[2])});
        else if (f[0] == "strandless") a.strandless_ceiling[std::stoll(f[1])] = std::stoi(f[2]);
    }
    return a;
}

std::vector<Window> load_panel(const std::string& path) {
    std::vector<Window> out;
    std::ifstream in(path);
    std::string line;
    bool header = true;
    while (std::getline(in, line)) {
        if (header) { header = false; continue; }
        if (line.empty()) continue;
        const auto f = split_tabs(line);
        if (f.size() < 3) continue;
        out.push_back({std::stoll(f[0]), std::stoll(f[1]), std::stoll(f[2])});
    }
    return out;
}

struct Paths {
    std::string binary = env_or("PGPHASE_BIN", "./pgphase");
    std::string test_data = env_or("PGPHASE_TEST_DATA", "test_data");
    std::string panel = env_or("PGPHASE_PANEL", "evaluations/2026-09-16-test-panel/panel.tsv");
    std::string workdir = env_or("PGPHASE_TEST_WORKDIR", "/tmp/pgphase-injection-tests");
    std::string allowance = env_or("PGPHASE_INJECTION_ALLOW",
                                   "src/test_bam_site_injection_allow.tsv");
    bool complete() const {
        return file_exists(binary) && file_exists(panel) &&
               file_exists(test_data + "/chm13v2.0.chr20.renamed.fa");
    }
    std::string missing() const {
        std::string m;
        if (!file_exists(binary)) m += "pgphase binary (" + binary + ", run make)";
        if (!file_exists(test_data + "/chm13v2.0.chr20.renamed.fa"))
            m += (m.empty() ? "" : ", ") + std::string("test_data/");
        if (!file_exists(panel)) m += (m.empty() ? "" : ", ") + std::string("panel (" + panel + ")");
        return m;
    }
};

/// Both channels over one window. `hybrid` is the alignment channel plus the
/// graph catalog injected into it; `alignment` is the same machinery with no
/// graph channel at all, which is the reference the injection has to preserve.
struct Pair {
    std::map<Key, Candidate> alignment, hybrid;
    std::vector<Record> alignment_records, hybrid_records;
};

bool run(const Paths& p, const std::string& subcommand, const Window& w,
         const std::string& extra, std::string& outdir) {
    std::ostringstream dir;
    dir << p.workdir << "/" << subcommand << "/w" << w.gap_left;
    outdir = dir.str();
    std::system(("mkdir -p '" + outdir + "'").c_str());
    std::ostringstream cmd;
    cmd << "mkdir -p '" << outdir << "' && '" << p.binary << "' " << subcommand
        << " --ref '" << p.test_data << "/chm13v2.0.chr20.renamed.fa'"
        << " --bam '" << p.test_data
        << "/HG002_chr20_hifi_mapped_to_CHM13_chr20_annotated.bam'" << extra
        << " -r 'CHM13#0#chr20:" << (w.gap_left - 50000) << "-" << (w.gap_right + 50000) << "'"
        << " -t " << test_threads()
        << " -o '" << outdir << "/candidates.tsv'"
        << " --phased-vcf-out '" << outdir << "/native.vcf'"
        << " -b '" << outdir << "/phased.bam'"
        << " > '" << outdir << "/stdout.log' 2> '" << outdir << "/stderr.log'";
    // The exact invocation is written beside its logs: a failure message names
    // the directory, and the first thing worth seeing there is the command that
    // produced it -- including the thread count actually used.
    {
        std::ofstream rec(outdir + "/cmd.txt");
        rec << cmd.str() << "\n";
    }
    return std::system(cmd.str().c_str()) == 0;
}

const Pair& channels(const Paths& p, const Window& w) {
    static std::map<long long, Pair> cache;
    const auto hit = cache.find(w.gap_left);
    if (hit != cache.end()) return hit->second;

    std::string dir_a, dir_h;
    const bool ok_a = run(p, "collect-bam-variation", w, "", dir_a);
    const std::string graph = " --graph-sites '" + p.test_data + "/chr20.sites.striped.vcf.gz'" +
                              " --gaf '" + p.test_data + "/HG002.chr20.annotated.coord.gaf.gz'";
    const bool ok_h = run(p, "collect-hybrid-variation", w, graph, dir_h);
    INFO("logs under " << dir_a << " and " << dir_h);
    REQUIRE(ok_a);
    REQUIRE(ok_h);
    Pair pr;
    pr.alignment = load_candidates(dir_a + "/candidates.tsv");
    pr.hybrid = load_candidates(dir_h + "/candidates.tsv");
    pr.alignment_records = load_records(dir_a + "/native.vcf");
    pr.hybrid_records = load_records(dir_h + "/native.vcf");
    REQUIRE(!pr.alignment.empty());
    REQUIRE(!pr.hybrid.empty());
    return cache.emplace(w.gap_left, std::move(pr)).first->second;
}

/// The read intervals in the window, from the BAM itself.
///
/// A candidate's depth cannot exceed the number of reads that OVERLAP its
/// reference span -- a read has to touch the site to be observed at it. Counting
/// coverage at the anchor base only is not that bound: a read overlapping part
/// of a 24 bp deletion can contribute an observation while not covering the
/// anchor, which is why the anchor count was exceeded by one read at
/// chr20:55,905,752 (DP+LOW_QUAL 57 against 56) by a record that is not in fact
/// over-counted.
///
/// This is the check that catches double counting, which is worth a bound that
/// cannot be argued with: before the backfill was rewritten to derive counts
/// rather than accumulate them, a read could be counted once as low-quality
/// depth and again as an allele vote, and depth_with_low_quality then exceeded
/// the site's own coverage -- the quantity the minimum-depth gate reads.
const std::vector<std::pair<long long, long long>>& read_spans(const Paths& p,
                                                               const Window& w) {
    static std::map<long long, std::vector<std::pair<long long, long long>>> cache;
    const auto hit = cache.find(w.gap_left);
    if (hit != cache.end()) return hit->second;
    std::vector<std::pair<long long, long long>> spans;
    const std::string bam =
        p.test_data + "/HG002_chr20_hifi_mapped_to_CHM13_chr20_annotated.bam";
    samFile* fp = sam_open(bam.c_str(), "r");
    REQUIRE(fp != nullptr);
    bam_hdr_t* hdr = sam_hdr_read(fp);
    REQUIRE(hdr != nullptr);
    hts_idx_t* idx = sam_index_load(fp, bam.c_str());
    REQUIRE(idx != nullptr);
    std::ostringstream region;
    region << "CHM13#0#chr20:" << (w.gap_left - 51000) << "-" << (w.gap_right + 51000);
    hts_itr_t* itr = sam_itr_querys(idx, hdr, region.str().c_str());
    REQUIRE(itr != nullptr);
    bam1_t* rec = bam_init1();
    while (sam_itr_next(fp, itr, rec) >= 0) {
        if (rec->core.flag & (BAM_FUNMAP | BAM_FSECONDARY | BAM_FSUPPLEMENTARY)) continue;
        if (rec->core.qual < 1) continue;
        spans.emplace_back(rec->core.pos + 1, bam_endpos(rec));  // 1-based inclusive
    }
    bam_destroy1(rec);
    hts_itr_destroy(itr);
    hts_idx_destroy(idx);
    bam_hdr_destroy(hdr);
    sam_close(fp);
    REQUIRE(!spans.empty());
    return cache.emplace(w.gap_left, std::move(spans)).first->second;
}

int overlapping(const std::vector<std::pair<long long, long long>>& spans,
                long long beg, long long end) {
    int n = 0;
    for (const auto& [b, e] : spans)
        if (b <= end && e >= beg) ++n;
    return n;
}

/// The reference base(s) at a 1-based position, from the FASTA the run used.
/// A candidate's REF is a claim about the reference; checking it against the
/// reference itself is the one way to catch a coordinate convention error that
/// every internally-consistent check would pass.
std::string reference_at(const Paths& p, long long pos, int len) {
    static faidx_t* fai = nullptr;
    if (fai == nullptr) {
        fai = fai_load((p.test_data + "/chm13v2.0.chr20.renamed.fa").c_str());
        REQUIRE(fai != nullptr);
    }
    std::ostringstream reg;
    reg << "CHM13#0#chr20:" << pos << "-" << (pos + len - 1);
    int got = 0;
    char* seq = fai_fetch(fai, reg.str().c_str(), &got);
    if (seq == nullptr || got <= 0) return {};
    std::string out(seq, static_cast<size_t>(got));
    free(seq);
    for (char& c : out) c = static_cast<char>(toupper(static_cast<unsigned char>(c)));
    return out;
}

/// The same window run again with a wider flank. Retrieval of a site in the
/// interior is a statement about the reads there, so it must not depend on how
/// much sequence was asked for around it.
const std::map<Key, Candidate>& widened(const Paths& p, const Window& w) {
    static std::map<long long, std::map<Key, Candidate>> cache;
    const auto hit = cache.find(w.gap_left);
    if (hit != cache.end()) return hit->second;
    Window wide = w;
    wide.gap_left = w.gap_left - 20000;
    wide.gap_right = w.gap_right + 20000;
    std::string dir;
    const std::string graph = " --graph-sites '" + p.test_data + "/chr20.sites.striped.vcf.gz'" +
                              " --gaf '" + p.test_data + "/HG002.chr20.annotated.coord.gaf.gz'";
    // A distinct output directory, or it would collide with the primary run.
    const std::string keep = p.workdir;
    Paths wp = p;
    wp.workdir = p.workdir + "/widened";
    REQUIRE(run(wp, "collect-hybrid-variation", wide, graph, dir));
    auto loaded = load_candidates(dir + "/candidates.tsv");
    REQUIRE(!loaded.empty());
    return cache.emplace(w.gap_left, std::move(loaded)).first->second;
}

/// The catalog's claims over a region, as (pos, ref, alts).
struct Claim { long long pos = 0; std::string ref; std::vector<std::string> alts; };

std::vector<Claim> load_claims(const Paths& p, long long lo, long long hi) {
    std::vector<Claim> out;
    const std::string path = p.test_data + "/chr20.sites.striped.vcf.gz";
    htsFile* fp = hts_open(path.c_str(), "r");
    if (fp == nullptr) return out;
    tbx_t* tbx = tbx_index_load(path.c_str());
    if (tbx == nullptr) { hts_close(fp); return out; }
    std::ostringstream reg;
    reg << "CHM13#0#chr20:" << lo << "-" << hi;
    hts_itr_t* itr = tbx_itr_querys(tbx, reg.str().c_str());
    kstring_t line = {0, 0, nullptr};
    while (itr != nullptr && tbx_itr_next(fp, tbx, itr, &line) >= 0) {
        const auto f = split_tabs(std::string(line.s, line.l));
        if (f.size() < 5) continue;
        Claim c;
        c.pos = std::stoll(f[1]);
        c.ref = f[3];
        for (char& ch : c.ref) ch = static_cast<char>(toupper(static_cast<unsigned char>(ch)));
        for (std::string a : split_char(f[4], ',')) {
            if (a == "*") continue;
            for (char& ch : a) ch = static_cast<char>(toupper(static_cast<unsigned char>(ch)));
            c.alts.push_back(a);
        }
        out.push_back(std::move(c));
    }
    free(line.s);
    if (itr != nullptr) tbx_itr_destroy(itr);
    tbx_destroy(tbx);
    hts_close(fp);
    return out;
}

bool inputs_ready(const Paths& p) {
    if (p.complete()) return true;
    WARN("injection tests need inputs that are absent: " << p.missing());
    return false;
}

}  // namespace

// 1. Every site the alignment channel finds must survive injection.
TEST_CASE("injection: the alignment channel's sites all reach the hybrid",
          "[injection][completeness]") {
    const Paths p;
    if (!inputs_ready(p)) { SUCCEED("skipped: inputs absent"); return; }
    for (const auto& w : load_panel(p.panel)) {
        DYNAMIC_SECTION("window " << w.gap_left) {
            const Pair& ch = channels(p, w);
            std::vector<std::string> lost;
            for (const auto& [k, c] : ch.alignment)
                if (ch.hybrid.find(k) == ch.hybrid.end()) lost.push_back(show(k));
            INFO("alignment channel candidates: " << ch.alignment.size()
                 << ", hybrid: " << ch.hybrid.size());
            if (!lost.empty()) {
                std::ostringstream o;
                for (const auto& s : lost) o << "\n    " << s;
                FAIL("injection dropped " << lost.size() << " candidate(s) the alignment "
                     "channel found:" << o.str());
            }
            // Injection may only ADD. A hybrid table smaller than the alignment
            // table means something was displaced rather than merged.
            CHECK(ch.hybrid.size() >= ch.alignment.size());
        }
    }
}

// 2. A shared site must keep its allele set.
TEST_CASE("injection: shared sites keep their alleles", "[injection][representation]") {
    const Paths p;
    if (!inputs_ready(p)) { SUCCEED("skipped: inputs absent"); return; }
    for (const auto& w : load_panel(p.panel)) {
        DYNAMIC_SECTION("window " << w.gap_left) {
            const Pair& ch = channels(p, w);
            // Positions where the alignment channel holds a record carrying two
            // alternates. The hybrid must still hold one there: reducing it to a
            // single alternate is the collapse that reads as homozygous.
            int multiallelic = 0, preserved = 0;
            std::vector<std::string> collapsed;
            for (const auto& [k, c] : ch.alignment) {
                const std::string& alt = std::get<3>(k);
                const bool two_alleles = alt.find(',') != std::string::npos;
                if (!two_alleles) continue;
                ++multiallelic;
                bool still_there = false;
                for (const auto& [hk, hc] : ch.hybrid) {
                    if (std::get<0>(hk) != std::get<0>(k)) continue;
                    if (std::get<3>(hk).find(',') != std::string::npos) { still_there = true; break; }
                }
                if (still_there) ++preserved; else collapsed.push_back(show(k));
            }
            INFO("multiallelic records in the alignment channel: " << multiallelic
                 << ", still multiallelic in the hybrid: " << preserved);
            if (!collapsed.empty()) {
                std::ostringstream o;
                for (const auto& s : collapsed) o << "\n    " << s;
                FAIL("injection reduced " << collapsed.size()
                     << " multiallelic locus/loci to a single alternate:" << o.str());
            }
            // And for every shared key the alleles are identical by construction
            // (they are part of the key), so what remains to check is that the
            // hybrid did not also add a SECOND record at the same position with
            // a different allele -- two descriptions of one locus.
            std::map<long long, int> het_in_hybrid, het_in_alignment;
            for (const auto& [k, c] : ch.hybrid)
                if (c.get("CATEGORY").find("HET") != std::string::npos)
                    ++het_in_hybrid[std::get<0>(k)];
            for (const auto& [k, c] : ch.alignment)
                if (c.get("CATEGORY").find("HET") != std::string::npos)
                    ++het_in_alignment[std::get<0>(k)];
            std::vector<std::string> doubled;
            for (const auto& [pos, n] : het_in_hybrid) {
                const auto a = het_in_alignment.find(pos);
                const int in_alignment = a == het_in_alignment.end() ? 0 : a->second;
                if (n > 1 && n > in_alignment) {
                    std::ostringstream o;
                    o << pos << " (" << in_alignment << " het record(s) in the alignment channel, "
                      << n << " in the hybrid)";
                    doubled.push_back(o.str());
                }
            }
            const Allowance allow = load_allowance(p.allowance);
            std::vector<std::string> unexpected;
            for (const auto& [pos, n] : het_in_hybrid) {
                const auto a = het_in_alignment.find(pos);
                const int in_alignment = a == het_in_alignment.end() ? 0 : a->second;
                if (n > 1 && n > in_alignment &&
                    allow.duplicates.count({w.gap_left, pos}) == 0) {
                    std::ostringstream o;
                    o << pos << " (" << in_alignment << " het record(s) in the alignment channel, "
                      << n << " in the hybrid)";
                    unexpected.push_back(o.str());
                }
            }
            INFO("known duplicate positions in this window: "
                 << std::count_if(allow.duplicates.begin(), allow.duplicates.end(),
                                  [&](const std::pair<long long, long long>& d) {
                                      return d.first == w.gap_left; }));
            if (!unexpected.empty()) {
                std::ostringstream o;
                for (const auto& s : unexpected) o << "\n    " << s;
                FAIL("injection added a second het description at " << unexpected.size()
                     << " position(s) not in the allowance file:" << o.str());
            }
            // The doubled list is kept for the message above; every entry in it
            // that is allowed is a known defect, not a passing case.
            (void)doubled;

            // A record naming two alternates with reads behind both cannot be
            // genotyped homozygous. Keeping the alleles is only half of correct
            // representation: 55,883,019 carries 29 reads on one alternate and
            // 33 on the other with ZERO reference, and a homozygous call there
            // describes neither haplotype. This is the half the merge exists for
            // -- split into two biallelic records, each allele would be measured
            // against a reference no read carries, its allele fraction would run
            // to 1, and both would classify homozygous.
            //
            // Checked in BOTH channels. The alignment channel is where these
            // records are built, so a wrong genotype there is the origin; the
            // hybrid is where it would be carried. Looking only at the hybrid
            // hides the defect entirely whenever the hybrid does not emit the
            // record at all, which is currently every one of them.
            std::vector<std::string> hom_at_multiallelic;
            for (const auto& [channel, records] :
                 {std::make_pair("alignment", &ch.alignment_records),
                  std::make_pair("hybrid", &ch.hybrid_records)}) {
                for (const Record& r : *records) {
                    if (!r.multiallelic() || !r.homozygous_call()) continue;
                    if (r.alts_with_reads() < 2) continue;
                    if (allow.genotype.count({w.gap_left, r.pos}) > 0) continue;
                    std::ostringstream o;
                    o << channel << " " << r.pos << " " << r.ref.substr(0, 18) << ">"
                      << r.alt.substr(0, 22) << " GT=" << r.gt << " AD=";
                    for (size_t i = 0; i < r.ad.size(); ++i) o << (i ? "," : "") << r.ad[i];
                    hom_at_multiallelic.push_back(o.str());
                }
            }

            // How many of the alignment channel's multiallelic loci the hybrid
            // also EMITS. Surviving into the candidate table is not the same as
            // reaching the output: the hybrid's solve excludes the noisy class
            // (skip_noisy_kmeans, a hybrid-only override), so a merged
            // NOISY_CAND_HET record gets no phase set and is never written. The
            // count is recorded rather than asserted to be equal, because
            // admitting that class chunk-wide is a measured bad trade -- but it
            // is recorded so that a change in either direction is visible.
            int align_multi = 0, hybrid_multi = 0;
            for (const Record& r : ch.alignment_records) if (r.multiallelic()) ++align_multi;
            for (const Record& r : ch.hybrid_records) if (r.multiallelic()) ++hybrid_multi;
            const auto emitted = allow.emitted_multiallelic.find(w.gap_left);
            INFO("multiallelic records emitted: alignment " << align_multi << ", hybrid "
                 << hybrid_multi << " (recorded " 
                 << (emitted == allow.emitted_multiallelic.end() ? -1 : emitted->second) << ")");
            CHECK(emitted != allow.emitted_multiallelic.end());
            if (emitted != allow.emitted_multiallelic.end())
                CHECK(hybrid_multi >= emitted->second);
            if (!hom_at_multiallelic.empty()) {
                std::ostringstream o;
                for (const auto& x : hom_at_multiallelic) o << "\n    " << x;
                FAIL("multiallelic record(s) genotyped homozygous with reads on both "
                     "alternates: " << hom_at_multiallelic.size() << o.str());
            }
        }
    }
}

// 3. A shared site's read counts must be the alignment channel's, and internally
//    consistent; an injected site's counts must at least be internally consistent.
TEST_CASE("injection: read counts and alleles per site are correct",
          "[injection][counts]") {
    const Paths p;
    if (!inputs_ready(p)) { SUCCEED("skipped: inputs absent"); return; }
    for (const auto& w : load_panel(p.panel)) {
        DYNAMIC_SECTION("window " << w.gap_left) {
            const Pair& ch = channels(p, w);
            const auto& spans = read_spans(p, w);
            int shared = 0, injected = 0, drifted = 0;
            std::vector<std::string> bad_depth, bad_internal;

            const Allowance allow = load_allowance(p.allowance);
            int strandless_alignment = 0;
            auto internally_consistent = [&](const Candidate& c, const Key& k, bool is_injected,
                                             std::vector<std::string>& out) {
                const long long dp = c.num("DP"), ref = c.num("REF_COUNT"),
                                alt = c.num("ALT_COUNT"), lq = c.num("LOW_QUAL_COUNT");
                const long long fr = c.num("FORWARD_REF"), rr = c.num("REVERSE_REF"),
                                fa = c.num("FORWARD_ALT"), ra = c.num("REVERSE_ALT");
                std::ostringstream o;
                if (dp != ref + alt)
                    o << " DP " << dp << " != ref+alt " << (ref + alt) << ";";
                if (ref + alt > 0) {
                    const double af = static_cast<double>(alt) / (ref + alt);
                    if (std::abs(af - c.real("AF")) > 1e-4)
                        o << " AF " << c.real("AF") << " != alt/(ref+alt) " << af << ";";
                }
                // The strand tallies must sum to the counts they accompany: they
                // are mirrored at each accumulation site rather than re-derived,
                // so a mismatch means a site that adds coverage without strand.
                // Asserted for INJECTED candidates, whose accumulation sites were
                // fixed to carry the strand. Alignment-derived candidates built
                // by the MSA path still leave these at zero -- the same defect in
                // a different path, counted below against a ceiling rather than
                // asserted away.
                const bool strand_missing =
                    (ref > 0 && fr + rr == 0) || (alt > 0 && fa + ra == 0);
                if (strand_missing && !is_injected) ++strandless_alignment;
                else if (is_injected || !strand_missing) {
                    if (ref > 0 && fr + rr != ref)
                        o << " ref strand " << fr << "+" << rr << " != REF_COUNT " << ref << ";";
                    if (alt > 0 && fa + ra != alt)
                        o << " alt strand " << fa << "+" << ra << " != ALT_COUNT " << alt << ";";
                }
                if (lq < 0) o << " negative LOW_QUAL_COUNT;";
                if (!o.str().empty()) out.push_back(show(k) + " --" + o.str());
            };

            for (const auto& [k, hc] : ch.hybrid) {
                const auto it = ch.alignment.find(k);
                const bool is_injected = it == ch.alignment.end();
                if (is_injected) ++injected; else ++shared;
                internally_consistent(hc, k, is_injected, bad_internal);

                // Depth cannot exceed the reads covering the site. This is the
                // check that catches double counting: before the backfill was
                // rewritten to derive counts instead of accumulating them, a
                // read could be counted once as low-quality depth and again as
                // an allele vote, and depth_with_low_quality then exceeded the
                // site's own coverage -- which is what the minimum-depth gate
                // reads.
                const long long pos = std::get<0>(k);
                const long long ref_len =
                    std::max<long long>(1, static_cast<long long>(std::get<2>(k).size()));
                const int over = overlapping(spans, pos, pos + ref_len - 1);
                const long long dp = hc.num("DP") + hc.num("LOW_QUAL_COUNT");
                if (over > 0 && dp > over && allow.overdepth.count({w.gap_left, pos}) == 0) {
                    std::ostringstream o;
                    o << " DP+LOW_QUAL " << dp << " exceeds " << over << " overlapping reads;";
                    bad_depth.push_back(show(k) + " --" + o.str());
                }

                // Cross-channel count drift is REPORTED, not asserted. The two
                // channels legitimately differ by a read or two: measured at
                // chr20:12,754,263 the hybrid's 32/16 matches the reads while
                // the alignment channel reports 32/18, and at 12,797,914 it is
                // the other way round. Neither is the reference for the other.
                if (it != ch.alignment.end() &&
                    (hc.num("DP") != it->second.num("DP") ||
                     hc.num("REF_COUNT") != it->second.num("REF_COUNT") ||
                     hc.num("ALT_COUNT") != it->second.num("ALT_COUNT"))) ++drifted;
            }

            INFO("shared with the alignment channel: " << shared << ", injected: " << injected
                 << ", shared sites whose counts differ between channels: " << drifted);
            auto report = [](const char* what, const std::vector<std::string>& v) {
                if (v.empty()) return;
                std::ostringstream o;
                for (size_t i = 0; i < v.size() && i < 12; ++i) o << "\n    " << v[i];
                if (v.size() > 12) o << "\n    ... " << (v.size() - 12) << " more";
                FAIL(what << ": " << v.size() << o.str());
            };
            report("candidates whose own fields are inconsistent", bad_internal);
            report("candidates whose depth exceeds their coverage", bad_depth);
            const auto ceiling = allow.strandless_ceiling.find(w.gap_left);
            INFO("alignment-derived candidates with no strand tally: " << strandless_alignment
                 << " (documented ceiling "
                 << (ceiling == allow.strandless_ceiling.end() ? -1 : ceiling->second) << ")");
            CHECK(ceiling != allow.strandless_ceiling.end());
            if (ceiling != allow.strandless_ceiling.end())
                CHECK(strandless_alignment <= ceiling->second);
            CHECK(shared > 0);
        }
    }
}


// 4. A candidate's REF must be what the reference actually says -- under the
//    convention that candidate type actually uses.
//
//    The candidate table does NOT use one convention for all types, and the
//    difference is a trap worth pinning down in a test rather than rediscovering:
//
//      SNP, DEL : POS is the site, REF is the reference string starting AT POS.
//      INS      : POS is one past the anchor, REF holds the ANCHOR base, which
//                 is the reference base at POS-1, and ALT holds only the
//                 INSERTED bases (not anchor+inserted as VCF would write them).
//
//    Measured on the panel: 186 SNP and 21 DEL candidates match the reference at
//    POS, and 30 of 36 INS candidates match at POS-1 with the remaining 6
//    ambiguous because the anchor sits in a homopolymer and both positions carry
//    the same base. Asserting the VCF convention for insertions -- which an
//    earlier version of this test did -- fails all 30 and says nothing true.
//
//    The emitted VCF is checked separately against the VCF convention, because
//    that is the output a consumer reads.
TEST_CASE("retrieval: REF matches the reference sequence", "[injection][alleles]") {
    const Paths p;
    if (!inputs_ready(p)) { SUCCEED("skipped: inputs absent"); return; }
    for (const auto& w : load_panel(p.panel)) {
        DYNAMIC_SECTION("window " << w.gap_left) {
            const Pair& ch = channels(p, w);
            const Allowance allow = load_allowance(p.allowance);
            std::vector<std::string> mismatched;
            int checked = 0;
            for (const auto& [k, c] : ch.hybrid) {
                const std::string& ref = std::get<2>(k);
                const std::string& type = std::get<1>(k);
                if (ref.empty() || ref == ".") continue;
                // An insertion has TWO legitimate encodings and the table carries
                // no ref_len column to tell them apart, so either is accepted:
                //   ref_len == 0 (a clean insertion, replacing nothing): REF is
                //       the ANCHOR base, the reference at POS-1.
                //   ref_len  > 0 (a claim whose REF ran past the anchor, so those
                //       bases are replaced): REF is the CONSUMED reference, at POS.
                const long long at = (type == "INS") ? std::get<0>(k) - 1 : std::get<0>(k);
                const std::string actual = reference_at(p, at, static_cast<int>(ref.size()));
                if (actual.empty()) continue;
                ++checked;
                bool ok = (actual == ref);
                if (!ok && type == "INS")
                    ok = (reference_at(p, std::get<0>(k), static_cast<int>(ref.size())) == ref);
                if (!ok) {
                    std::ostringstream o;
                    o << show(k) << " -- reference at " << at << " says " << actual;
                    mismatched.push_back(o.str());
                }
            }
            INFO("candidates whose REF was checked against the FASTA: " << checked);
            REQUIRE(checked > 0);
            if (!mismatched.empty()) {
                std::ostringstream o;
                for (size_t i = 0; i < mismatched.size() && i < 10; ++i) o << "\n    " << mismatched[i];
                if (mismatched.size() > 10) o << "\n    ... " << (mismatched.size() - 10) << " more";
                FAIL("REF disagrees with the reference at " << mismatched.size()
                     << " candidate(s):" << o.str());
            }

            // The emitted records, against the VCF convention a consumer relies
            // on: REF is the reference string at POS, and for an indel every ALT
            // begins with REF.
            std::vector<std::string> bad_vcf;
            int vcf_checked = 0;
            for (const Record& r : ch.hybrid_records) {
                if (r.ref.empty()) continue;
                const std::string actual =
                    reference_at(p, r.pos, static_cast<int>(r.ref.size()));
                if (actual.empty()) continue;
                ++vcf_checked;
                std::ostringstream why;
                if (actual != r.ref) why << " REF!=" << actual << ";";
                for (const std::string& alt : split_char(r.alt, ',')) {
                    if (alt.size() == r.ref.size()) continue;  // substitution
                    // A length-changing record is anchored: REF and ALT share
                    // the base at POS. Requiring ALT to begin with the WHOLE REF
                    // is right only for a simple indel -- a complex event, where
                    // the claim replaces the bases it consumes, shares just the
                    // anchor. chr20:55,905,389 CT>CCG is valid VCF and an earlier
                    // version of this check called it broken.
                    if (alt.empty() || r.ref.empty() || alt[0] != r.ref[0])
                        why << " ALT " << alt.substr(0, 14) << " shares no anchor with REF;";
                }
                if (!why.str().empty())
                    bad_vcf.push_back(std::to_string(r.pos) + " " + r.ref.substr(0, 14) + ">" +
                                      r.alt.substr(0, 14) + " --" + why.str());
            }
            INFO("emitted records checked against the VCF convention: " << vcf_checked);
            if (!bad_vcf.empty()) {
                std::ostringstream o;
                for (size_t i = 0; i < bad_vcf.size() && i < 10; ++i) o << "\n    " << bad_vcf[i];
                if (bad_vcf.size() > 10) o << "\n    ... " << (bad_vcf.size() - 10) << " more";
                FAIL("emitted record(s) breaking the VCF REF/ALT convention: "
                     << bad_vcf.size() << o.str());
            }
        }
    }
}

// 5. The candidate table and the emitted records must describe one thing.
TEST_CASE("retrieval: emitted records agree with the candidate table",
          "[injection][consistency]") {
    const Paths p;
    if (!inputs_ready(p)) { SUCCEED("skipped: inputs absent"); return; }
    for (const auto& w : load_panel(p.panel)) {
        DYNAMIC_SECTION("window " << w.gap_left) {
            const Pair& ch = channels(p, w);
            const long long lo = w.gap_left - 50000, hi = w.gap_right + 50000;

            // Every candidate lies inside the region that was asked for.
            std::vector<std::string> outside;
            for (const auto& [k, c] : ch.hybrid) {
                const long long pos = std::get<0>(k);
                if (pos < lo || pos > hi) outside.push_back(show(k));
            }
            if (!outside.empty()) {
                std::ostringstream o;
                for (size_t i = 0; i < outside.size() && i < 8; ++i) o << "\n    " << outside[i];
                FAIL("candidate(s) outside the requested region " << lo << "-" << hi << ": "
                     << outside.size() << o.str());
            }

            // An emitted record's depths must be the candidate's depths. The VCF
            // is written from the candidate, so a disagreement means one of the
            // two was rewritten after the other was produced.
            int matched = 0;
            std::vector<std::string> disagree, orphan;
            for (const Record& r : ch.hybrid_records) {
                if (r.multiallelic() || r.ad.size() < 2) continue;  // 1-ALT records only
                const Key k{r.pos, r.ref.size() == r.alt.size() && r.ref.size() == 1 ? "SNP"
                            : (r.alt.size() > r.ref.size() ? "INS" : "DEL"), r.ref, r.alt};
                const auto it = ch.hybrid.find(k);
                if (it == ch.hybrid.end()) continue;  // representation differs; not this test
                ++matched;
                const long long cref = it->second.num("REF_COUNT");
                const long long calt = it->second.num("ALT_COUNT");
                if (r.ad[0] != cref || r.ad[1] != calt) {
                    std::ostringstream o;
                    o << r.pos << " VCF AD=" << r.ad[0] << "," << r.ad[1]
                      << " against candidate " << cref << "/" << calt;
                    disagree.push_back(o.str());
                }
            }
            INFO("emitted records matched to their candidate: " << matched);
            if (!disagree.empty()) {
                std::ostringstream o;
                for (size_t i = 0; i < disagree.size() && i < 10; ++i) o << "\n    " << disagree[i];
                if (disagree.size() > 10) o << "\n    ... " << (disagree.size() - 10) << " more";
                FAIL("emitted depths disagree with the candidate table at " << disagree.size()
                     << " record(s):" << o.str());
            }
            CHECK(matched > 0);
        }
    }
}

// 6. Retrieval must not depend on how much sequence was asked for around a site.
TEST_CASE("retrieval: a wider window finds the same interior sites",
          "[injection][stability]") {
    const Paths p;
    if (!inputs_ready(p)) { SUCCEED("skipped: inputs absent"); return; }
    for (const auto& w : load_panel(p.panel)) {
        DYNAMIC_SECTION("window " << w.gap_left) {
            const Pair& ch = channels(p, w);
            const std::map<Key, Candidate>& wide = widened(p, w);
            // Compare only the interior, well inside BOTH runs' flanks: a site
            // near an edge legitimately sees fewer reads in the narrower run.
            const long long lo = w.gap_left - 25000, hi = w.gap_right + 25000;
            auto interior = [&](const Key& k) {
                const long long pos = std::get<0>(k);
                return pos >= lo && pos <= hi;
            };
            std::vector<std::string> only_narrow, only_wide;
            for (const auto& [k, c] : ch.hybrid)
                if (interior(k) && wide.find(k) == wide.end()) only_narrow.push_back(show(k));
            for (const auto& [k, c] : wide)
                if (interior(k) && ch.hybrid.find(k) == ch.hybrid.end()) only_wide.push_back(show(k));
            INFO("interior " << lo << "-" << hi << "; narrow run " << ch.hybrid.size()
                 << " candidates, wide run " << wide.size());
            auto report = [](const char* what, const std::vector<std::string>& v) {
                if (v.empty()) return;
                std::ostringstream o;
                for (size_t i = 0; i < v.size() && i < 10; ++i) o << "\n    " << v[i];
                if (v.size() > 10) o << "\n    ... " << (v.size() - 10) << " more";
                FAIL(what << ": " << v.size() << o.str());
            };
            report("interior sites the narrow window found and the wider one did not", only_narrow);
            report("interior sites the wider window found and the narrow one did not", only_wide);
        }
    }
}


// 7. A claim's alleles have to be carried whole, reference bases included.
//
// A claim whose REF runs past the anchor describes a locus where those extra
// reference bases are REPLACED. Emitting the claim's ALT while consuming only
// the anchor leaves the rest of the claim's REF in place, so the record asserts
// a haplotype one or more bases longer than the claim did -- the alleles look
// carried, and the sequence is wrong.
//
// Found at chr20:12,680,256, where the catalog claims AT > AAATAAAATAAAATA and
// we emit A > AAATAAAATAAAATA, leaving the T: our record means
// AAATAAAATAAAATA-T where the claim means AAATAAAATAAAATA. The candidate table
// shows the same locus as INS REF=T at 12,680,257 -- the only candidate on the
// panel whose REF is the base at its own POS rather than the anchor, which is
// how the discrepancy surfaces.
//
// Deliberately narrow: only a record carrying the claim's ALT VERBATIM while
// consuming fewer reference bases is flagged. A record describing a different
// event at the same position is the alignment channel's own call and is not
// this defect.
TEST_CASE("retrieval: an injected claim's reference bases are consumed whole",
          "[injection][claims]") {
    const Paths p;
    if (!inputs_ready(p)) { SUCCEED("skipped: inputs absent"); return; }
    const Allowance allow = load_allowance(p.allowance);
    for (const auto& w : load_panel(p.panel)) {
        DYNAMIC_SECTION("window " << w.gap_left) {
            const Pair& ch = channels(p, w);
            const long long lo = w.gap_left - 50000, hi = w.gap_right + 50000;
            std::map<long long, std::vector<const Record*>> by_pos;
            for (const Record& r : ch.hybrid_records) by_pos[r.pos].push_back(&r);

            const auto claims = load_claims(p, lo, hi);
            REQUIRE(!claims.empty());
            int examined = 0;
            std::vector<std::string> partial;
            for (const Claim& c : claims) {
                if (c.ref.size() < 2) continue;
                bool has_insertion = false;
                for (const std::string& a : c.alts) has_insertion |= a.size() > c.ref.size();
                if (!has_insertion) continue;
                const auto at = by_pos.find(c.pos);
                if (at == by_pos.end()) continue;
                ++examined;
                for (const Record* r : at->second) {
                    if (r->ref.size() >= c.ref.size()) continue;
                    for (const std::string& ours : split_char(r->alt, ',')) {
                        for (const std::string& a : c.alts) {
                            if (a.size() <= c.ref.size() || ours != a) continue;
                            if (allow.claims.count({w.gap_left, c.pos}) > 0) continue;
                            std::ostringstream o;
                            o << c.pos << " claim " << c.ref << ">" << a.substr(0, 22)
                              << " emitted as " << r->ref << ">" << ours.substr(0, 22)
                              << " (consumed " << r->ref.size() << " of " << c.ref.size()
                              << " reference bases)";
                            partial.push_back(o.str());
                        }
                    }
                }
            }
            INFO("claims with REF past the anchor and an emitted record at the same position: "
                 << examined);
            if (!partial.empty()) {
                std::ostringstream o;
                for (size_t i = 0; i < partial.size() && i < 10; ++i) o << "\n    " << partial[i];
                if (partial.size() > 10) o << "\n    ... " << (partial.size() - 10) << " more";
                FAIL("claim(s) emitted with their reference bases only partly consumed: "
                     << partial.size() << o.str());
            }
        }
    }
}


// 8. A site the alignment channel already had must arrive as it was.
//
// "As they are" is not one rule for every site, because injection legitimately
// re-measures some of them. Measured across the panel on 1,517 shared keys:
//
//   clean classes (CLEAN_HET_SNP, CLEAN_HET_INDEL, CLEAN_HOM): 1,201 sites,
//       ZERO count changes and ZERO demotions. Injection does not touch them.
//   noisy classes: 80 of 316 sites have their counts re-measured, and 105 are
//       promoted to a clean class by the graph claim.
//
// The noisy re-measurement is a repair, not damage, which is why this test does
// not demand identity there. At the four largest changes the alignment channel
// reports DP 7-11 where 71-78 reads overlap -- the depth-starvation signature --
// and the hybrid reports 69-77. Requiring identity would pin the starved value.
// The bound that does apply to them, DP within the overlapping reads, is
// asserted in [counts].
//
// So: clean sites must be byte-identical and must never be demoted; noisy drift
// is reported, and its direction is held to promotions only.
TEST_CASE("retrieval: a site the alignment channel had arrives unchanged",
          "[injection][fidelity]") {
    const Paths p;
    if (!inputs_ready(p)) { SUCCEED("skipped: inputs absent"); return; }
    static const char* kCountFields[] = {"DP",          "REF_COUNT",   "ALT_COUNT",
                                         "LOW_QUAL_COUNT", "FORWARD_REF", "REVERSE_REF",
                                         "FORWARD_ALT", "REVERSE_ALT", "AF"};
    auto is_clean = [](const std::string& c) { return c.rfind("CLEAN", 0) == 0; };
    for (const auto& w : load_panel(p.panel)) {
        DYNAMIC_SECTION("window " << w.gap_left) {
            const Pair& ch = channels(p, w);
            int clean = 0, noisy = 0, promoted = 0, noisy_recounted = 0;
            std::vector<std::string> altered, demoted;
            for (const auto& [k, a] : ch.alignment) {
                const auto it = ch.hybrid.find(k);
                if (it == ch.hybrid.end()) continue;  // absence is [completeness]'s job
                const Candidate& h = it->second;
                const std::string acat = a.get("CATEGORY"), hcat = h.get("CATEGORY");
                std::string changed;
                for (const char* f : kCountFields)
                    if (a.get(f) != h.get(f))
                        changed += std::string(" ") + f + " " + a.get(f) + "->" + h.get(f) + ";";
                if (is_clean(acat)) {
                    ++clean;
                    if (!changed.empty()) altered.push_back(show(k) + " [" + acat + "] --" + changed);
                    if (!is_clean(hcat))
                        demoted.push_back(show(k) + " " + acat + " -> " + hcat);
                } else {
                    ++noisy;
                    if (!changed.empty()) ++noisy_recounted;
                    if (is_clean(hcat)) ++promoted;
                    else if (hcat != acat)
                        demoted.push_back(show(k) + " " + acat + " -> " + hcat);
                }
            }
            INFO("clean-class shared sites " << clean << "; noisy-class " << noisy
                 << ", of which " << noisy_recounted << " re-measured and " << promoted
                 << " promoted");
            REQUIRE(clean > 0);
            auto report = [](const char* what, const std::vector<std::string>& v) {
                if (v.empty()) return;
                std::ostringstream o;
                for (size_t i = 0; i < v.size() && i < 10; ++i) o << "\n    " << v[i];
                if (v.size() > 10) o << "\n    ... " << (v.size() - 10) << " more";
                FAIL(what << ": " << v.size() << o.str());
            };
            report("clean-class site(s) whose counts injection changed", altered);
            report("site(s) injection demoted", demoted);
        }
    }
}
