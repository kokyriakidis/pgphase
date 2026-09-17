// Regression tests over the chr20 gap windows.
//
// These are integration tests, not unit tests: each one runs the pipeline binary
// on one window of the real test data, parses the phased VCF and phased BAM it
// produced, and scores the result against the parental truth. A window costs
// about 1.2 s, so the whole panel is a few seconds and belongs in a normal test
// run rather than in a benchmark script.
//
// WHY A WINDOW TEST AND NOT A UNIT TEST. Everything this session found in the
// gap windows was invisible to unit tests and to read-level accuracy alike: a
// phantom site bridging the k-means, a verdict blocking its own correction, a
// duplicate record beside a correctly merged locus, an inverted flank behind a
// join no read spans. Those are properties of the pipeline's output on real
// data. Tests that assert them have to run the pipeline.
//
// WHAT IS ASSERTED, per window and per arm:
//
//   spans        -- does one phase set bracket the whole gap? A window that
//                   spans when the expectation says it should not is as much a
//                   failure as the reverse: a join across an interval no read
//                   crosses is a coin flip, and chr20:48,176,830-48,229,446
//                   was once reported CLOSED at 100% read accuracy while its
//                   two halves sat on opposite haplotypes.
//   in-gap hets  -- phased heterozygotes strictly inside the gap. This is the
//                   quantity the retry moves; the default leaves the noisy
//                   class unoriented, so it phases only the boundary sites.
//   tagged       -- reads carrying HP, i.e. coverage gained.
//   concordance  -- of the reads that both carry HP and appear in truth, the
//                   fraction agreeing with the majority orientation of their
//                   own phase set. Scored per phase set, because a phase set's
//                   labels are only meaningful relative to itself.
//   discordant   -- the absolute count, floored as well as the rate: a rate can
//                   be held up by coverage while reads go wrong.
//
// Expectations are committed in src/test_gap_windows_expect.tsv and are FLOORS
// and CEILINGS, not equalities: a change that phases more sites correctly must
// not have to edit this file, while one that loses coverage or flips a read
// must fail. Regenerate with tools/refresh_gap_window_expectations.sh only when
// a measured improvement is intended, and say so in the commit.

#define CATCH_CONFIG_MAIN
#include "../third_party/catch2/catch.hpp"

#include <htslib/sam.h>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <map>
#include <sstream>
#include <string>
#include <unordered_map>
#include <vector>

namespace {

std::string env_or(const char* key, const std::string& fallback) {
    const char* v = std::getenv(key);
    return (v != nullptr && *v != '\0') ? std::string(v) : fallback;
}

/// Threads for every pipeline run these tests make, never fewer than four.
/// PGPHASE_TEST_THREADS may raise it; a lower value is clamped up, so a test run
/// cannot be made accidentally serial. Four is where the scaling stops paying:
/// on chr20:30-35 Mb, 161.6 s at 1 thread against 46.7 s at 4, and 43.6 s at 12.
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

/// One row of the committed panel: the window, and what a competitor achieves on
/// it. The competitor columns are context for a failure message, not assertions
/// -- the tests assert against our own committed expectations.
struct Window {
    long long gap_left = 0;
    long long gap_right = 0;
    long long gap_bp = 0;
    std::string best_competitor;
    double competitor_acc = 0.0;
};

/// A window's expected outcome under one arm. Floors and ceilings, never
/// equalities, so an improvement does not have to edit the file.
struct Expectation {
    bool spans = false;        // per-window rows: asserted as an EQUALITY
    int min_spanned = 0;       // TOTAL rows only: the spans column read as a count
    int min_in_gap_hets = 0;
    int min_tagged = 0;
    double min_concordance = 0.0;
    int max_discordant = 0;
};

/// What one run of the pipeline produced on one window.
struct Outcome {
    bool spans = false;
    int in_gap_hets = 0;
    int tagged = 0;
    int scored = 0;
    int correct = 0;
    int blocks = 0;
    double concordance() const {
        return scored > 0 ? static_cast<double>(correct) / scored : 0.0;
    }
    int discordant() const { return scored - correct; }
};

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
        Window w;
        w.gap_left = std::stoll(f[0]);
        w.gap_right = std::stoll(f[1]);
        w.gap_bp = std::stoll(f[2]);
        if (f.size() > 4) w.best_competitor = f[4];
        if (f.size() > 5) w.competitor_acc = std::stod(f[5]);
        out.push_back(w);
    }
    return out;
}

/// Keyed by "<arm>\t<gap_left>" so one file holds every arm.
std::map<std::string, Expectation> load_expectations(const std::string& path) {
    std::map<std::string, Expectation> out;
    std::ifstream in(path);
    std::string line;
    while (std::getline(in, line)) {
        if (line.empty() || line[0] == '#') continue;
        const auto f = split_tabs(line);
        if (f.size() < 7 || f[0] == "arm") continue;
        Expectation e;
        // On a TOTAL row the spans column is a COUNT of windows spanned; on a
        // per-window row it is a boolean. One column, two readings, because the
        // alternative is a second file that can fall out of step with this one.
        if (f[1] == "TOTAL") e.min_spanned = std::stoi(f[2]);
        else e.spans = (f[2] == "1" || f[2] == "yes" || f[2] == "true");
        e.min_in_gap_hets = std::stoi(f[3]);
        e.min_tagged = std::stoi(f[4]);
        e.min_concordance = std::stod(f[5]);
        e.max_discordant = std::stoi(f[6]);
        out[f[0] + "\t" + f[1]] = e;
    }
    return out;
}

std::unordered_map<std::string, char> load_truth(const std::string& path) {
    std::unordered_map<std::string, char> out;
    std::ifstream in(path);
    std::string line;
    while (std::getline(in, line)) {
        const auto tab = line.find('\t');
        if (tab == std::string::npos) continue;
        const std::string hap = line.substr(tab + 1);
        if (hap.rfind("MATERNAL", 0) == 0) out[line.substr(0, tab)] = 'M';
        else if (hap.rfind("PATERNAL", 0) == 0) out[line.substr(0, tab)] = 'P';
    }
    return out;
}

/// Whether a VCF genotype string describes a phased heterozygote. 1|2 counts:
/// a locus whose two haplotypes are both non-reference is exactly the case the
/// multiallelic work exists to represent, and excluding it would make the test
/// blind to losing it.
bool is_phased_het(const std::string& gt) {
    const auto bar = gt.find('|');
    if (bar == std::string::npos || bar == 0 || bar + 1 >= gt.size()) return false;
    return gt.substr(0, bar) != gt.substr(bar + 1);
}

/// Per phase set, the first and last phased heterozygote, and how many fall
/// strictly inside the gap.
void parse_vcf(const std::string& path, const Window& w, Outcome& out) {
    std::ifstream in(path);
    std::string line;
    std::map<std::string, std::pair<long long, long long>> extent;
    while (std::getline(in, line)) {
        if (line.empty() || line[0] == '#') continue;
        const auto f = split_tabs(line);
        if (f.size() < 10) continue;
        const std::string sample = f[9];
        const std::string gt = sample.substr(0, sample.find(':'));
        if (!is_phased_het(gt)) continue;
        const long long pos = std::stoll(f[1]);
        // PS is located through FORMAT rather than assumed to be a fixed field.
        std::string ps;
        {
            std::istringstream ks(f[8]), vs(sample);
            std::string k, v;
            while (std::getline(ks, k, ':') && std::getline(vs, v, ':'))
                if (k == "PS") { ps = v; break; }
        }
        if (ps.empty() || ps == "." || ps == "0") continue;
        auto it = extent.find(ps);
        if (it == extent.end()) extent[ps] = {pos, pos};
        else {
            it->second.first = std::min(it->second.first, pos);
            it->second.second = std::max(it->second.second, pos);
        }
        if (pos > w.gap_left && pos < w.gap_right) ++out.in_gap_hets;
    }
    out.blocks = static_cast<int>(extent.size());
    for (const auto& [ps, span] : extent)
        if (span.first <= w.gap_left && span.second >= w.gap_right) out.spans = true;
}

/// Score the phased BAM per phase set: a phase set's HP labels are arbitrary up
/// to a global flip, so the orientation truth prefers is chosen per phase set
/// and the reads that disagree with it are the discordant ones.
void score_bam(const std::string& path,
               const std::unordered_map<std::string, char>& truth,
               Outcome& out) {
    samFile* fp = sam_open(path.c_str(), "r");
    REQUIRE(fp != nullptr);
    bam_hdr_t* hdr = sam_hdr_read(fp);
    REQUIRE(hdr != nullptr);
    bam1_t* rec = bam_init1();

    // per phase set: [hap1-with-MAT, hap1-with-PAT] as the two orientations
    std::unordered_map<long long, std::pair<int, int>> votes;
    while (sam_read1(fp, hdr, rec) >= 0) {
        if (rec->core.flag & BAM_FUNMAP) continue;
        const uint8_t* hp = bam_aux_get(rec, "HP");
        if (hp == nullptr) continue;
        ++out.tagged;
        const uint8_t* ps = bam_aux_get(rec, "PS");
        if (ps == nullptr) continue;
        const auto found = truth.find(bam_get_qname(rec));
        if (found == truth.end()) continue;
        const long long set_id = bam_aux2i(ps);
        const long long hap = bam_aux2i(hp);
        auto& v = votes[set_id];
        const bool hap1 = (hap == 1);
        if ((hap1 && found->second == 'M') || (!hap1 && found->second == 'P')) ++v.first;
        else ++v.second;
    }
    for (const auto& [set_id, v] : votes) {
        (void)set_id;
        out.scored += v.first + v.second;
        out.correct += std::max(v.first, v.second);
    }
    bam_destroy1(rec);
    bam_hdr_destroy(hdr);
    sam_close(fp);
}

struct Paths {
    std::string binary, test_data, panel, expectations, truth_map, workdir;
    bool complete() const {
        const bool emitting = std::getenv("PGPHASE_EMIT_EXPECTATIONS") != nullptr;
        return file_exists(binary) && file_exists(panel) && file_exists(truth_map) &&
               (emitting || file_exists(expectations)) &&
               file_exists(test_data + "/chm13v2.0.chr20.renamed.fa");
    }
    std::string missing() const {
        std::string m;
        auto note = [&m](bool ok, const std::string& what) {
            if (!ok) m += (m.empty() ? "" : ", ") + what;
        };
        note(file_exists(binary), "pgphase binary (" + binary + ", run make)");
        note(file_exists(test_data + "/chm13v2.0.chr20.renamed.fa"), "test_data/");
        note(file_exists(panel), "panel (" + panel + ")");
        note(file_exists(expectations), "expectations (" + expectations + ")");
        note(file_exists(truth_map),
             "truth map (" + truth_map + ", run scripts/make_truth_hap_map.sh)");
        return m;
    }
};

Paths paths() {
    Paths p;
    p.binary = env_or("PGPHASE_BIN", "./pgphase");
    p.test_data = env_or("PGPHASE_TEST_DATA", "test_data");
    p.panel = env_or("PGPHASE_PANEL", "evaluations/2026-09-16-test-panel/panel.tsv");
    p.expectations = env_or("PGPHASE_EXPECT", "src/test_gap_windows_expect.tsv");
    p.truth_map = env_or("PGPHASE_TRUTH_MAP", "test_data/derived/chr20_truth_hap.tsv");
    p.workdir = env_or("PGPHASE_TEST_WORKDIR", "/tmp/pgphase-window-tests");
    return p;
}

/// Run one arm over one window. Returns false when the binary failed, leaving
/// its stderr on disk for the failure message.
bool run_arm(const Paths& p, const Window& w, const std::string& arm,
             const std::string& flags, std::string& outdir) {
    std::ostringstream dir;
    dir << p.workdir << "/" << arm << "/w" << w.gap_left;
    outdir = dir.str();
    std::ostringstream cmd;
    cmd << "mkdir -p '" << outdir << "' && '" << p.binary << "' collect-hybrid-variation"
        << " --ref '" << p.test_data << "/chm13v2.0.chr20.renamed.fa'"
        << " --bam '" << p.test_data
        << "/HG002_chr20_hifi_mapped_to_CHM13_chr20_annotated.bam'"
        << " --graph-sites '" << p.test_data << "/chr20.sites.striped.vcf.gz'"
        << " --gaf '" << p.test_data << "/HG002.chr20.annotated.coord.gaf.gz'"
        << " -r 'CHM13#0#chr20:" << (w.gap_left - 50000) << "-" << (w.gap_right + 50000) << "'"
        << " -t " << test_threads() << " " << flags
        << " -o '" << outdir << "/candidates.tsv'"
        << " --phased-vcf-out '" << outdir << "/native.vcf'"
        << " -b '" << outdir << "/phased.bam'"
        << " > '" << outdir << "/stdout.log' 2> '" << outdir << "/stderr.log'";
    return std::system(cmd.str().c_str()) == 0;
}

/// Memoised across test cases: Catch2 runs them in one process, and the panel
/// totals are the per-window outcomes summed, so re-running the pipeline for
/// them would double the suite's runtime for no extra coverage.
std::map<std::string, Outcome>& outcome_cache() {
    static std::map<std::string, Outcome> cache;
    return cache;
}

Outcome measure_uncached(const Paths& p, const Window& w, const std::string& arm,
                         const std::string& flags,
                         const std::unordered_map<std::string, char>& truth) {
    std::string dir;
    const bool ok = run_arm(p, w, arm, flags, dir);
    INFO("arm '" << arm << "' on window " << w.gap_left << "-" << w.gap_right
         << "; outputs and logs under " << dir);
    REQUIRE(ok);
    Outcome out;
    parse_vcf(dir + "/native.vcf", w, out);
    score_bam(dir + "/phased.bam", truth, out);
    return out;
}

Outcome measure(const Paths& p, const Window& w, const std::string& arm,
                const std::string& flags,
                const std::unordered_map<std::string, char>& truth) {
    const std::string key = arm + "\t" + std::to_string(w.gap_left);
    auto& cache = outcome_cache();
    const auto hit = cache.find(key);
    if (hit != cache.end()) return hit->second;
    const Outcome out = measure_uncached(p, w, arm, flags, truth);
    cache.emplace(key, out);
    return out;
}

void check_against(const Window& w, const std::string& arm, const Outcome& got,
                   const Expectation& want) {
    INFO("window " << w.gap_left << "-" << w.gap_right << " (" << w.gap_bp
         << " bp), arm '" << arm << "'; best competitor " << w.best_competitor
         << " at " << (100.0 * w.competitor_acc) << "%");
    INFO("measured: spans=" << (got.spans ? "yes" : "no")
         << " in_gap_hets=" << got.in_gap_hets << " blocks=" << got.blocks
         << " tagged=" << got.tagged << " scored=" << got.scored
         << " concordance=" << got.concordance() << " discordant=" << got.discordant());

    // Spanning is an equality, in both directions: a span that appears where the
    // expectation says there is none is the coin-flip join, not an improvement.
    CHECK(got.spans == want.spans);
    CHECK(got.in_gap_hets >= want.min_in_gap_hets);
    CHECK(got.tagged >= want.min_tagged);
    CHECK(got.scored > 0);
    CHECK(got.concordance() >= want.min_concordance);
    CHECK(got.discordant() <= want.max_discordant);
}

}  // namespace

/// Emit the expectations file instead of asserting, so a refresh cannot drift
/// from the measurement. Floors are set to what was measured and the discordant
/// ceiling likewise, which makes any later loss of coverage or accuracy a
/// failure while leaving room for improvement.
///
/// Written to a file rather than stdout on purpose: Catch2 captures stdout
/// during a test case and reports it indented, so a refresh that piped stdout
/// would silently produce an unparseable file.
/// A floor printed to four decimals has to be rounded DOWN. Rounding to nearest
/// produced 0.9976 for a measured 0.99757869, i.e. a floor above the
/// measurement, so a freshly generated file failed against the very run that
/// produced it.
double floor4(double v) { return std::floor(v * 10000.0) / 10000.0; }

void emit_expectations(const std::string& out_path, const Paths& p,
                       const std::vector<Window>& panel,
                       const std::unordered_map<std::string, char>& truth,
                       const std::vector<std::pair<std::string, std::string>>& arms) {
    std::FILE* out = std::fopen(out_path.c_str(), "w");
    INFO("cannot write expectations to " << out_path);
    REQUIRE(out != nullptr);
    std::fprintf(out, "# Expected outcome per arm and window for src/test_gap_windows.cpp.\n");
    std::fprintf(out, "# Floors and ceilings, not equalities, EXCEPT the spans column, which is\n");
    std::fprintf(out, "# asserted exactly per window (a span appearing where none is expected is\n");
    std::fprintf(out, "# a join across an interval no read crosses, not an improvement) and read\n");
    std::fprintf(out, "# as a count of spanned windows on the TOTAL rows.\n");
    std::fprintf(out, "# Regenerate with scripts/refresh_gap_window_expectations.sh.\n");
    std::fprintf(out, "arm\twindow\tspans\tmin_in_gap_hets\tmin_tagged\tmin_concordance\tmax_discordant\n");
    for (const auto& [arm, flags] : arms) {
        int spanned = 0, hets = 0, tagged = 0, scored = 0, correct = 0;
        for (const auto& w : panel) {
            const Outcome got = measure(p, w, arm, flags, truth);
            std::fprintf(out, "%s\t%lld\t%d\t%d\t%d\t%.4f\t%d\n", arm.c_str(), w.gap_left,
                        got.spans ? 1 : 0, got.in_gap_hets, got.tagged,
                        floor4(got.concordance()), got.discordant());
            spanned += got.spans ? 1 : 0;
            hets += got.in_gap_hets; tagged += got.tagged;
            scored += got.scored; correct += got.correct;
        }
        std::fprintf(out, "%s\tTOTAL\t%d\t%d\t%d\t%.4f\t%d\n", arm.c_str(), spanned, hets,
                    tagged, floor4(scored ? double(correct) / scored : 0.0),
                    scored - correct);
    }
    std::fclose(out);
    WARN("wrote expectations to " << out_path);
}

TEST_CASE("chr20 gap windows", "[gap][windows]") {
    const Paths p = paths();
    if (!p.complete()) {
        WARN("gap-window tests need inputs that are absent: " << p.missing());
        SUCCEED("skipped: inputs absent");
        return;
    }
    const auto panel = load_panel(p.panel);
    const auto truth = load_truth(p.truth_map);
    REQUIRE(!panel.empty());
    REQUIRE(truth.size() > 1000);

    // The two arms that exist. Adding a third means adding its rows to the
    // expectations file; a window with no row for an arm fails loudly rather
    // than being skipped, so the file cannot silently fall behind the panel.
    const std::vector<std::pair<std::string, std::string>> arms = {
        {"default", ""},
        {"retry", "--retry-unphased-with-bam"},
    };

    // PGPHASE_EMIT_EXPECTATIONS holds the path to write; "1" means the default.
    const std::string emit = env_or("PGPHASE_EMIT_EXPECTATIONS", "");
    if (!emit.empty()) {
        emit_expectations(emit == "1" ? p.expectations : emit, p, panel, truth, arms);
        SUCCEED("emitted expectations");
        return;
    }
    const auto expect = load_expectations(p.expectations);
    REQUIRE(!expect.empty());

    for (const auto& [arm, flags] : arms) {
        for (const auto& w : panel) {
            DYNAMIC_SECTION(arm << " / window " << w.gap_left) {
                const auto key = arm + "\t" + std::to_string(w.gap_left);
                const auto it = expect.find(key);
                if (it == expect.end())
                    FAIL("no expectation row for arm '" << arm << "' window "
                         << w.gap_left << " in " << p.expectations
                         << " -- run scripts/refresh_gap_window_expectations.sh");
                const Outcome got = measure(p, w, arm, flags, truth);
                check_against(w, arm, got, it->second);
            }
        }
    }
}

TEST_CASE("chr20 gap windows: panel totals", "[gap][windows][totals]") {
    const Paths p = paths();
    if (!p.complete()) {
        WARN("gap-window tests need inputs that are absent: " << p.missing());
        SUCCEED("skipped: inputs absent");
        return;
    }
    // The first test case does the emitting; this one has nothing to add to the
    // file and must not assert against expectations that do not exist yet.
    if (std::getenv("PGPHASE_EMIT_EXPECTATIONS") != nullptr) {
        SUCCEED("skipped while emitting expectations");
        return;
    }
    const auto panel = load_panel(p.panel);
    const auto expect = load_expectations(p.expectations);
    const auto truth = load_truth(p.truth_map);

    // A per-window test can pass everywhere while the panel as a whole moves,
    // because each window's floor is generous on its own. These totals are the
    // numbers the work is actually steered by.
    struct Totals { int spanned = 0, in_gap_hets = 0, tagged = 0, scored = 0, correct = 0; };
    std::map<std::string, Totals> totals;
    const std::vector<std::pair<std::string, std::string>> arms = {
        {"default", ""},
        {"retry", "--retry-unphased-with-bam"},
    };
    for (const auto& [arm, flags] : arms) {
        for (const auto& w : panel) {
            const Outcome got = measure(p, w, arm, flags, truth);
            auto& t = totals[arm];
            t.spanned += got.spans ? 1 : 0;
            t.in_gap_hets += got.in_gap_hets;
            t.tagged += got.tagged;
            t.scored += got.scored;
            t.correct += got.correct;
        }
    }
    for (const auto& [arm, t] : totals) {
        const auto it = expect.find(arm + "\tTOTAL");
        if (it == expect.end())
            FAIL("no TOTAL row for arm '" << arm << "' in " << p.expectations);
        INFO("arm '" << arm << "' totals: spanned=" << t.spanned
             << " in_gap_hets=" << t.in_gap_hets << " tagged=" << t.tagged
             << " scored=" << t.scored
             << " concordance=" << (t.scored ? double(t.correct) / t.scored : 0.0)
             << " discordant=" << (t.scored - t.correct));
        const Expectation& want = it->second;
        CHECK(t.spanned >= want.min_spanned);
        CHECK(t.in_gap_hets >= want.min_in_gap_hets);
        CHECK(t.tagged >= want.min_tagged);
        CHECK((t.scored ? double(t.correct) / t.scored : 0.0) >= want.min_concordance);
        CHECK((t.scored - t.correct) <= want.max_discordant);
    }
}
