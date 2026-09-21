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
#include <set>
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
/// What a window is expected to do.
///
/// Deliberately NOT a read count. How many reads end up tagged is a consequence
/// of decisions the solver is free to make differently -- a margin, a link
/// threshold, which of two equally good sites anchors a block -- so a floor on
/// it fails on changes that are improvements, and pinning it invites a refresh
/// that launders a real regression. The two properties worth asserting are
/// whether the gap is closed and whether the sites we hold and need are used.
/// Rounded DOWN to two decimals. A floor printed with round-to-nearest can sit
/// ABOVE the measurement that produced it, so a freshly emitted file fails
/// against its own run.
inline double floor2(double v) { return std::floor(v * 100.0) / 100.0; }

struct Expectation {
    bool spans = false;        // per-window rows: asserted as an EQUALITY
    int min_spanned = 0;       // TOTAL rows only: the spans column read as a count
    int min_in_gap_hets = 0;
    /// Per-arm read-concordance floor. Not a read COUNT -- a rate, and the only
    /// numeric quality bound here. It is per arm because the arms are not
    /// equally mature: the graph-first configuration pays a read-PLACEMENT cost
    /// (a third of the clean het anchors outside a gap are alignment-only,
    /// including every clean het indel in the regions examined), and recording
    /// its measured floor states that deficit in the baseline instead of either
    /// hiding it behind a loose global bound or letting it block the suite.
    double min_concordance = 0.0;
    /// Floor on Outcome::separated() -- the fraction of the window's scorable
    /// reads that ONE block places correctly. The headline quantity: a change
    /// that splits a window keeps concordance at 1.00 and drops this.
    double min_separated = 0.0;
};

/// What one run of the pipeline produced on one window.
struct Outcome {
    std::string diagnosis;  // why a red span check is red; see explain_gap
    /// Positions whose records put two DIFFERENT alleles on the SAME haplotype.
    /// Asserted at zero, not floored: one haplotype carries one allele.
    int hap_allele_conflicts = 0;
    bool spans = false;
    int in_gap_hets = 0;
    int tagged = 0;
    int scored = 0;
    int correct = 0;
    int blocks = 0;
    /// Reads correctly separated into ONE block, over every read in the window
    /// that truth can score. This is the product the pipeline exists to make:
    /// purity alone hides fragmentation (two immaculate half-blocks separate
    /// nobody), and block count alone hides switches. Measured on the largest
    /// phase set touching the window.
    int dominant_correct = 0;
    int window_scorable = 0;
    /// Candidates inside the gap that are in an admitted class -- a CLEAN het,
    /// which the stage-1 mask accepts -- and carry no phase set. A site we hold,
    /// and that the solve is allowed to use, left unused.
    int unused_clean_hets = 0;
    /// Required sites -- the ones the closing arm phases inside this gap -- that
    /// the run failed to retrieve at all, and ones it retrieved but left without
    /// a phase set. The first is a discovery or injection regression, the second
    /// an admission regression, and a `spans` check alone can pass through both
    /// by finding some other way across.
    int required_missing = 0;
    int required_unused = 0;
    std::string required_detail;
    /// Phased heterozygotes strictly inside the gap, as (pos, ref, alt), and the
    /// category of every in-gap candidate. Together these let the binary WRITE
    /// the required-sites file it later asserts against, so the file cannot
    /// drift from a second implementation of "which sites close this gap".
    std::vector<std::array<std::string, 3>> in_gap_sites;
    std::map<long long, std::string> gap_cats;
    /// A block that reaches both sides of the gap but puts a different parent on
    /// haplotype 1 at each end. Read-level accuracy cannot see this when no read
    /// crosses the gap, so it is asserted directly.
    bool switched = false;
    double concordance() const {
        return scored > 0 ? static_cast<double>(correct) / scored : 0.0;
    }
    int discordant() const { return scored - correct; }
    /// Fraction of the window's truth-scorable reads that one block places
    /// correctly. Contiguity and correctness in one number.
    double separated() const {
        return window_scorable > 0
                   ? static_cast<double>(dominant_correct) / window_scorable
                   : 0.0;
    }
};

std::vector<Window> load_panel(const std::string& path) {
    std::vector<Window> out;
    std::ifstream in(path);
    std::string line;
    bool header = true;
    while (std::getline(in, line)) {
        // Comments precede the header, so skip them BEFORE consuming it --
        // otherwise the first '#' line is eaten as the header and the real
        // header line reaches std::stoll. The panel documents its own two
        // selection bases in that block, so it has to survive the parse.
        if (!line.empty() && line[0] == '#') continue;
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
        if (f.size() < 4 || f[0] == "arm") continue;
        Expectation e;
        // On a TOTAL row the spans column is a COUNT of windows spanned; on a
        // per-window row it is a boolean. One column, two readings, because the
        // alternative is a second file that can fall out of step with this one.
        if (f[1] == "TOTAL") e.min_spanned = std::stoi(f[2]);
        else e.spans = (f[2] == "1" || f[2] == "yes" || f[2] == "true");
        e.min_in_gap_hets = std::stoi(f[3]);
        e.min_concordance = f.size() > 4 ? std::stod(f[4]) : 0.0;
        e.min_separated = f.size() > 5 ? std::stod(f[5]) : 0.0;
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
/// Count sites inside the gap that the solve was allowed to use and did not.
void parse_candidates(const std::string& path, const Window& w, Outcome& out) {
    std::ifstream in(path);
    std::string line;
    int pos_i = -1, cat_i = -1, ps_i = -1;
    bool header = true;
    while (std::getline(in, line)) {
        if (line.empty()) continue;
        const auto f = split_tabs(line);
        if (header) {
            header = false;
            for (size_t i = 0; i < f.size(); ++i) {
                if (f[i] == "POS") pos_i = static_cast<int>(i);
                else if (f[i] == "CATEGORY") cat_i = static_cast<int>(i);
                else if (f[i] == "PHASE_SET") ps_i = static_cast<int>(i);
            }
            continue;
        }
        if (pos_i < 0 || cat_i < 0 || ps_i < 0) return;
        if (f.size() <= static_cast<size_t>(std::max(pos_i, std::max(cat_i, ps_i)))) continue;
        const long long pos = std::stoll(f[static_cast<size_t>(pos_i)]);
        if (pos <= w.gap_left || pos >= w.gap_right) continue;
        out.gap_cats[pos] = f[static_cast<size_t>(cat_i)];
        if (f[static_cast<size_t>(cat_i)].rfind("CLEAN_HET", 0) != 0) continue;
        const std::string& ps = f[static_cast<size_t>(ps_i)];
        if (ps == "0" || ps.empty() || ps == ".") ++out.unused_clean_hets;
    }
}

/// The sites a closing arm's chain rests on, keyed by arm and window.
struct RequiredSite { long long pos = 0; std::string ref, alt, cat; };

std::map<std::string, std::vector<RequiredSite>> load_required(const std::string& path) {
    std::map<std::string, std::vector<RequiredSite>> out;
    std::ifstream in(path);
    std::string line;
    while (std::getline(in, line)) {
        if (line.empty() || line[0] == '#') continue;
        const auto f = split_tabs(line);
        if (f.size() < 6 || f[0] == "arm") continue;
        RequiredSite r;
        r.pos = std::stoll(f[2]);
        r.ref = f[3]; r.alt = f[4]; r.cat = f[5];
        out[f[0] + "\t" + f[1]].push_back(r);
    }
    return out;
}

/// Assert-side of the required list: is each site retrieved, and is it used?
///
/// Matching allows +/-2 bp because an insertion's candidate anchors one base
/// past the position the VCF reports, so an exact match would report a site as
/// missing that is present under the other convention.
void check_required_sites(const std::string& candidates_path, const std::string& arm,
                          const Window& w, Outcome& out) {
    const char* env = std::getenv("PGPHASE_REQUIRED");
    const std::string path = env != nullptr ? env : "src/test_gap_windows_required.tsv";
    if (!file_exists(path)) return;
    static const auto required = load_required(path);
    const auto it = required.find(arm + "\t" + std::to_string(w.gap_left));
    if (it == required.end()) return;

    std::map<long long, std::string> found;  // pos -> phase set
    std::ifstream in(candidates_path);
    std::string line;
    int pos_i = -1, ps_i = -1;
    bool header = true;
    while (std::getline(in, line)) {
        if (line.empty()) continue;
        const auto f = split_tabs(line);
        if (header) {
            header = false;
            for (size_t i = 0; i < f.size(); ++i) {
                if (f[i] == "POS") pos_i = static_cast<int>(i);
                else if (f[i] == "PHASE_SET") ps_i = static_cast<int>(i);
            }
            continue;
        }
        if (pos_i < 0 || ps_i < 0) return;
        if (f.size() <= static_cast<size_t>(std::max(pos_i, ps_i))) continue;
        found[std::stoll(f[static_cast<size_t>(pos_i)])] = f[static_cast<size_t>(ps_i)];
    }
    std::ostringstream detail;
    for (const RequiredSite& r : it->second) {
        const std::string* ps = nullptr;
        for (long long d = -2; d <= 2 && ps == nullptr; ++d) {
            const auto hit = found.find(r.pos + d);
            if (hit != found.end()) ps = &hit->second;
        }
        if (ps == nullptr) {
            ++out.required_missing;
            detail << "\n    NOT RETRIEVED " << r.pos << " " << r.ref << ">" << r.alt
                   << " (recorded as " << r.cat << ")";
        } else if (*ps == "0" || ps->empty() || *ps == ".") {
            ++out.required_unused;
            detail << "\n    RETRIEVED BUT UNUSED " << r.pos << " " << r.ref << ">"
                   << r.alt << " (recorded as " << r.cat << ")";
        }
    }
    out.required_detail = detail.str();
}

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
        if (pos > w.gap_left && pos < w.gap_right) {
            ++out.in_gap_hets;
            out.in_gap_sites.push_back({std::to_string(pos), f[3], f[4]});
        }
    }
    out.blocks = static_cast<int>(extent.size());
    for (const auto& [ps, span] : extent)
        if (span.first <= w.gap_left && span.second >= w.gap_right) out.spans = true;
}

/// Score the phased BAM per phase set: a phase set's HP labels are arbitrary up
/// to a global flip, so the orientation truth prefers is chosen per phase set
/// and the reads that disagree with it are the discordant ones.
/// Read spans by qname, taken from the INPUT alignment.
///
/// The graph arm writes its phased BAM UNALIGNED -- every record has flag 4, no
/// reference and no position -- because the haplotype call belongs to the read,
/// not to a placement. Scoring still needs coordinates: concordance needs to
/// know a read is in the window, and the switch check needs to know which side
/// of the gap it sits on. Those come from the input BAM, matched by name, which
/// is also the only source that cannot disagree with what the pipeline read.
using ReadSpans = std::unordered_map<std::string, std::pair<long long, long long>>;


/// Count positions whose records claim the same haplotype with two different
/// alleles. No truth needed: one haplotype carries one allele, so a SNP and the
/// insertion containing it (T>G with T>GAG) contradict each other on their face.
/// Asserted at zero rather than floored -- it is never acceptable.
int count_hap_allele_conflicts(const std::string& vcf_path) {
    std::ifstream in(vcf_path);
    if (!in) return 0;
    std::string line;
    std::map<long long, std::array<std::set<std::string>, 3>> claimed;
    while (std::getline(in, line)) {
        if (line.empty() || line[0] == '#') continue;
        const auto f = split_tabs(line);
        if (f.size() < 10) continue;
        const std::string gt = f[9].substr(0, f[9].find(':'));
        const long long pos = std::stoll(f[1]);
        if (gt == "1|0") claimed[pos][1].insert(f[4]);
        else if (gt == "0|1") claimed[pos][2].insert(f[4]);
    }
    int bad = 0;
    for (const auto& kv : claimed)
        if (kv.second[1].size() > 1 || kv.second[2].size() > 1) ++bad;
    return bad;
}

void score_bam(const std::string& path, const Window& w,
               const std::unordered_map<std::string, char>& truth,
               const ReadSpans& spans, Outcome& out) {
    samFile* fp = sam_open(path.c_str(), "r");
    REQUIRE(fp != nullptr);
    bam_hdr_t* hdr = sam_hdr_read(fp);
    REQUIRE(hdr != nullptr);
    bam1_t* rec = bam_init1();

    // per phase set: [hap1-with-MAT, hap1-with-PAT] as the two orientations
    std::unordered_map<long long, std::pair<int, int>> votes;
    // the same tally restricted to reads wholly on one side of the gap, so a
    // block reaching both sides can be checked for an internal switch
    std::unordered_map<long long, std::pair<int, int>> left_votes, right_votes;
    // Restricted to reads overlapping the WINDOW itself. `votes` spans the whole
    // flanked run, so using it as the numerator against a window-sized
    // denominator produced separated() above 1.0.
    std::unordered_map<long long, std::pair<int, int>> win_votes;
    while (sam_read1(fp, hdr, rec) >= 0) {
        if ((rec->core.flag & BAM_FUNMAP) && spans.empty()) continue;
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
        const bool mat_on_hap1 =
            (hap1 && found->second == 'M') || (!hap1 && found->second == 'P');
        if (mat_on_hap1) ++v.first; else ++v.second;
        long long beg = rec->core.pos;
        long long end = bam_endpos(rec);
        if (rec->core.flag & BAM_FUNMAP) {
            const auto sp = spans.find(bam_get_qname(rec));
            if (sp == spans.end()) continue;
            beg = sp->second.first;
            end = sp->second.second;
        }
        if (end >= w.gap_left && beg <= w.gap_right) {
            auto& wv = win_votes[set_id];
            if (mat_on_hap1) ++wv.first; else ++wv.second;
        }
        std::pair<int, int>* side = nullptr;
        if (end <= w.gap_left) side = &left_votes[set_id];
        else if (beg >= w.gap_right) side = &right_votes[set_id];
        if (side != nullptr) {
            if (mat_on_hap1) ++side->first; else ++side->second;
        }
    }
    // A side counts as placed only when it is confident: at least five scored
    // reads and at least 90% of them agreeing. Two confidently placed ends that
    // disagree are a switch inside one block.
    // The dominant block: the phase set scoring the most reads in this window,
    // and how many of those it places on the right parent.
    {
        int best_n = 0;
        for (const auto& kv : win_votes) {
            const int n = kv.second.first + kv.second.second;
            if (n > best_n) {
                best_n = n;
                out.dominant_correct = std::max(kv.second.first, kv.second.second);
            }
        }
        out.window_scorable = 0;
        for (const auto& [name, se] : spans) {
            if (se.second < w.gap_left || se.first > w.gap_right) continue;
            if (truth.count(name)) ++out.window_scorable;
        }
    }

    const auto placed = [](const std::pair<int, int>& v) -> int {
        const int n = v.first + v.second;
        if (n < 5) return -1;
        const int top = std::max(v.first, v.second);
        if (static_cast<double>(top) < 0.90 * n) return -1;
        return v.first >= v.second ? 1 : 0;
    };
    for (const auto& [set_id, lv] : left_votes) {
        const auto rit = right_votes.find(set_id);
        if (rit == right_votes.end()) continue;
        const int l = placed(lv), r = placed(rit->second);
        if (l >= 0 && r >= 0 && l != r) out.switched = true;
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
    std::string binary, test_data, panel, expectations, required, truth_map, workdir;
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
    p.required = env_or("PGPHASE_REQUIRED", "src/test_gap_windows_required.tsv");
    p.truth_map = env_or("PGPHASE_TRUTH_MAP", "test_data/derived/chr20_truth_hap.tsv");
    p.workdir = env_or("PGPHASE_TEST_WORKDIR", "/tmp/pgphase-window-tests");
    return p;
}

ReadSpans& input_read_spans(const Paths& p, const Window& w) {
    static std::map<hts_pos_t, ReadSpans> cache;
    auto it = cache.find(w.gap_left);
    if (it != cache.end()) return it->second;
    ReadSpans spans;
    const std::string bam = p.test_data +
        "/HG002_chr20_hifi_mapped_to_CHM13_chr20_annotated.bam";
    samFile* fp = sam_open(bam.c_str(), "r");
    REQUIRE(fp != nullptr);
    bam_hdr_t* hdr = sam_hdr_read(fp);
    REQUIRE(hdr != nullptr);
    hts_idx_t* idx = sam_index_load(fp, bam.c_str());
    REQUIRE(idx != nullptr);
    std::ostringstream reg;
    reg << "CHM13#0#chr20:" << (w.gap_left - 50000) << "-" << (w.gap_right + 50000);
    hts_itr_t* itr = sam_itr_querys(idx, hdr, reg.str().c_str());
    REQUIRE(itr != nullptr);
    bam1_t* rec = bam_init1();
    while (sam_itr_next(fp, itr, rec) >= 0) {
        if (rec->core.flag & (BAM_FUNMAP | BAM_FSECONDARY | BAM_FSUPPLEMENTARY)) continue;
        spans[bam_get_qname(rec)] = {rec->core.pos, bam_endpos(rec)};
    }
    bam_destroy1(rec);
    hts_itr_destroy(itr);
    hts_idx_destroy(idx);
    bam_hdr_destroy(hdr);
    sam_close(fp);
    return cache.emplace(w.gap_left, std::move(spans)).first->second;
}

/// Run one arm over one window. Returns false when the binary failed, leaving
/// its stderr on disk for the failure message.
bool run_arm(const Paths& p, const Window& w, const std::string& arm,
             const std::string& flags, std::string& outdir) {
    std::ostringstream dir;
    dir << p.workdir << "/" << arm << "/w" << w.gap_left;
    outdir = dir.str();
    std::ostringstream cmd;
    cmd << "mkdir -p '" << outdir << "' && '" << p.binary << "' collect-graph-variation"
        << " --ref '" << p.test_data << "/chm13v2.0.chr20.renamed.fa'"
        << " --bam '" << p.test_data
        << "/HG002_chr20_hifi_mapped_to_CHM13_chr20_annotated.bam'"
        << " --sites '" << p.test_data << "/chr20.sites.striped.vcf.gz'"
        << " --gaf '" << p.test_data << "/HG002.chr20.annotated.coord.gaf.gz'"
        << " -r 'CHM13#0#chr20:" << (w.gap_left - 50000) << "-" << (w.gap_right + 50000) << "'"
        << " -t " << test_threads() << " " << flags
        << " -o '" << outdir << "/candidates.tsv'"
        << " --phased-vcf-out '" << outdir << "/native.vcf'"
        << " --phased-bam-out '" << outdir << "/phased.bam'"
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

/// Why a window failed, in the three terms the cause can take.
///
/// A red span check says the gap did not close; it does not say whether the
/// gap CAN close, which is the first thing anyone iterating on recovery needs
/// and the thing that used to cost a manual drilldown each time:
///
///   UNCLOSEABLE -- some interior position is crossed by no read. No admission
///                  or linking change can help; a competitor spanning it is
///                  making a join its own reads do not support.
///   NO SITES    -- reads cross, but the gap interior holds no phased het, so
///                  the chain has nothing to step through. Admission.
///   NOT LINKED  -- reads cross AND interior sites are phased, but they did not
///                  end up in one phase set. Linking or orientation.
///
/// Computed from what the run already wrote, so it costs no extra pipeline run.
std::string explain_gap(const Paths& p, const Window& w, const Outcome& got) {
    const ReadSpans& spans = input_read_spans(p, w);
    // thinnest read coverage strictly inside the gap, sampled on a grid
    const long long step = std::max<long long>(1, (w.gap_right - w.gap_left) / 200);
    int thinnest = INT_MAX;
    long long thinnest_at = w.gap_left;
    for (long long x = w.gap_left + step; x < w.gap_right; x += step) {
        int n = 0;
        for (const auto& [name, se] : spans)
            if (se.first <= x && se.second >= x) ++n;
        if (n < thinnest) { thinnest = n; thinnest_at = x; }
    }
    std::ostringstream o;
    o << "gap diagnosis: thinnest interior read coverage " << thinnest
      << " at " << thinnest_at << "; phased hets inside the gap " << got.in_gap_hets
      << "; distinct phase sets touching the window " << got.blocks << " -- ";
    if (thinnest == 0)
        o << "UNCLOSEABLE (a position inside the gap is crossed by no read)";
    else if (got.in_gap_hets == 0)
        o << "NO SITES (reads cross it, but nothing inside is phased: admission)";
    else
        o << "NOT LINKED (reads cross and interior sites are phased, but they did"
             " not join: linking or orientation)";
    return o.str();
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
    parse_candidates(dir + "/candidates.tsv", w, out);
    check_required_sites(dir + "/candidates.tsv", arm, w, out);
    score_bam(dir + "/phased.bam", w, truth, input_read_spans(p, w), out);
    out.hap_allele_conflicts = count_hap_allele_conflicts(dir + "/native.vcf");
    out.diagnosis = explain_gap(p, w, out);
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
         << " concordance=" << got.concordance() << " discordant=" << got.discordant()
         << " separated=" << got.separated()
         << " (" << got.dominant_correct << "/" << got.window_scorable << ")");

    // Spanning is an equality, in both directions: a span that appears where the
    // expectation says there is none is the coin-flip join, not an improvement.
    // On a red span check, say WHY in the same breath: a failing window that
    // cannot be closed by anyone is a different fact from one whose sites were
    // refused, and iterating on recovery needs that distinction immediately.
    INFO(got.diagnosis);
    // Never acceptable, and cheap to notice: one haplotype with two alleles at
    // one position. Not floored against a recorded baseline -- asserted at zero.
    INFO("arm '" << arm << "' window " << w.gap_left << ": positions putting two "
         "alleles on one haplotype = " << got.hap_allele_conflicts);
    CHECK(got.hap_allele_conflicts == 0);

    CHECK(got.spans == want.spans);
    // Sites phased INSIDE the gap -- a count of sites, not of reads.
    CHECK(got.in_gap_hets >= want.min_in_gap_hets);
    // Every site in the gap that the solve was allowed to use must be used.
    INFO("clean hets inside the gap left without a phase set: " << got.unused_clean_hets);
    CHECK(got.unused_clean_hets == 0);
    // Every site the closing arm's chain rests on must still be retrieved and
    // still be used. A gap that closes some other way while one of these has
    // gone missing is not the same closure.
    INFO("required sites: " << got.required_missing << " not retrieved, "
         << got.required_unused << " retrieved but unused" << got.required_detail);
    CHECK(got.required_missing == 0);
    CHECK(got.required_unused == 0);
    // A block reaching both sides of the gap must not switch across it. This is
    // the hazard a read-level number cannot see: when no read crosses the gap,
    // an inverted join scores 100% and only the two ends disagree.
    CHECK_FALSE(got.switched);
    CHECK(got.scored > 0);
    // A floor, not a pinned measurement: it fires on a collapse, not on a
    // decision that shifts a handful of reads. Per arm, because the arms are not
    // equally mature -- see Expectation::min_concordance.
    CHECK(got.concordance() >= want.min_concordance);
    // Reads correctly separated into one block: contiguity and correctness at
    // once, and the number to move.
    CHECK(got.separated() >= want.min_separated);
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

void emit_expectations(const std::string& out_path, const Paths& p,
                       const std::vector<Window>& panel,
                       const std::unordered_map<std::string, char>& truth,
                       const std::vector<std::pair<std::string, std::string>>& arms) {
    std::FILE* out = std::fopen(out_path.c_str(), "w");
    INFO("cannot write expectations to " << out_path);
    REQUIRE(out != nullptr);
    std::fprintf(out, "# Expected outcome per arm and window for src/test_gap_windows.cpp.\n");
    std::fprintf(out, "# spans is an EQUALITY per window (a span appearing where none is expected is\n");
    std::fprintf(out, "# a join across an interval no read crosses, not an improvement) and a count\n");
    std::fprintf(out, "# of spanned windows on the TOTAL rows. min_in_gap_hets is a floor on SITES\n");
    std::fprintf(out, "# phased inside the gap.\n");
    std::fprintf(out, "#\n");
    std::fprintf(out, "# There is deliberately no read-count column. How many reads end up tagged\n");
    std::fprintf(out, "# follows from decisions the solver may make differently, so a floor on it\n");
    std::fprintf(out, "# fails on improvements and pinning it invites a refresh that launders a\n");
    std::fprintf(out, "# regression. The assertions that do not live in this file and cannot drift:\n");
    std::fprintf(out, "# no clean het inside a gap is left without a phase set, no block switches\n");
    std::fprintf(out, "# across a gap. min_concordance is a per-arm read-concordance floor -- a rate,\n");
    std::fprintf(out, "# not a count -- rounded DOWN, per arm because the arms are not equally\n");
    std::fprintf(out, "# mature: graph-first pays a read-placement cost and recording its measured\n");
    std::fprintf(out, "# floor states that in the baseline rather than hiding it behind a loose\n");
    std::fprintf(out, "# global bound or letting it block the suite.\n");
    std::fprintf(out, "# Regenerate with scripts/refresh_gap_window_expectations.sh.\n");
    std::fprintf(out, "# min_separated is the headline: the fraction of the window's truth-scorable\n");
    std::fprintf(out, "# reads that ONE block places correctly. Purity alone hides fragmentation --\n");
    std::fprintf(out, "# two immaculate half-blocks separate nobody -- and block count alone hides\n");
    std::fprintf(out, "# switches. Rounded DOWN, like the concordance floor.\n");
    std::fprintf(out, "arm\twindow\tspans\tmin_in_gap_hets\tmin_concordance\tmin_separated\n");
    for (const auto& [arm, flags] : arms) {
        int spanned = 0, hets = 0;
        double worst = 1.0, worst_sep = 1.0;
        for (const auto& w : panel) {
            const Outcome got = measure(p, w, arm, flags, truth);
            std::fprintf(out, "%s\t%lld\t%d\t%d\t%.2f\t%.2f\n", arm.c_str(), w.gap_left,
                        got.spans ? 1 : 0, got.in_gap_hets, floor2(got.concordance()),
                        floor2(got.separated()));
            spanned += got.spans ? 1 : 0;
            hets += got.in_gap_hets;
            worst = std::min(worst, got.concordance());
            worst_sep = std::min(worst_sep, got.separated());
        }
        std::fprintf(out, "%s\tTOTAL\t%d\t%d\t%.2f\t%.2f\n", arm.c_str(), spanned, hets,
                     floor2(worst), floor2(worst_sep));
    }
    std::fclose(out);
    WARN("wrote expectations to " << out_path);

    // The sites a closing arm's chain rests on, written from the same run that
    // produced the expectations so the two cannot disagree.
    std::FILE* req = std::fopen(p.required.c_str(), "w");
    INFO("cannot write required sites to " << p.required);
    REQUIRE(req != nullptr);
    std::fprintf(req, "# Sites whose retrieval and use close a gap.\n#\n");
    std::fprintf(req, "# Derived, not asserted from memory: for every window the closing arm\n");
    std::fprintf(req, "# spans, every heterozygote it phases STRICTLY INSIDE the gap is listed.\n");
    std::fprintf(req, "# Those are the sites the closure rests on -- remove one and the chain\n");
    std::fprintf(req, "# across the gap loses a step.\n#\n");
    std::fprintf(req, "# The test asserts two things per row, in the named arm:\n");
    std::fprintf(req, "#   RETRIEVED -- a candidate exists at the position (+/-2 bp, since an\n");
    std::fprintf(req, "#                insertion anchors its candidate one base past the VCF\n");
    std::fprintf(req, "#                position)\n");
    std::fprintf(req, "#   USED      -- that candidate carries a non-zero PHASE_SET\n#\n");
    std::fprintf(req, "# A site that stops being retrieved is a discovery or injection\n");
    std::fprintf(req, "# regression; one retrieved and no longer used is an admission\n");
    std::fprintf(req, "# regression. Either way the gap has lost the evidence that closes it,\n");
    std::fprintf(req, "# which a spans check alone can still pass by finding another way across.\n");
    std::fprintf(req, "arm\twindow\tpos\tref\talt\tcategory_when_recorded\n");
    for (const auto& [arm, flags] : arms) {
        for (const auto& w : panel) {
            const Outcome got = measure(p, w, arm, flags, truth);
            if (!got.spans) continue;
            for (const auto& site : got.in_gap_sites) {
                const long long pos = std::stoll(site[0]);
                std::string cat = "UNKNOWN";
                for (long long d = -2; d <= 2; ++d) {
                    const auto hit = got.gap_cats.find(pos + d);
                    if (hit != got.gap_cats.end()) { cat = hit->second; break; }
                }
                std::fprintf(req, "%s\t%lld\t%lld\t%s\t%s\t%s\n", arm.c_str(), w.gap_left,
                             pos, site[1].c_str(), site[2].c_str(), cat.c_str());
            }
        }
    }
    std::fclose(req);
    WARN("wrote required sites to " << p.required);
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
    // One arm. Recovery is not a mode any more: it runs inside every chunk's
    // own first solve, so there is nothing to switch on and nothing to compare
    // against but the competitor and the committed expectations.
    const std::vector<std::pair<std::string, std::string>> arms = {
        {"graph", ""},
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

TEST_CASE("recovery preserves complementary BAM deletion rows", "[gap][representation]") {
    const Paths p = paths();
    if (!p.complete()) {
        WARN("gap-window tests need inputs that are absent: " << p.missing());
        SUCCEED("skipped: inputs absent");
        return;
    }
    Window w;
    w.gap_left = 11235279;
    w.gap_right = 11262361;
    std::string dir;
    REQUIRE(run_arm(p, w, "graph", "", dir));

    std::ifstream in(dir + "/native.vcf");
    REQUIRE(in.good());
    std::set<std::pair<std::string, std::string>> alleles;
    int rows = 0;
    std::string line;
    while (std::getline(in, line)) {
        if (line.empty() || line[0] == '#') continue;
        const auto fields = split_tabs(line);
        if (fields.size() < 5 || fields[1] != "11255369") continue;
        ++rows;
        alleles.emplace(fields[3], fields[4]);
    }
    CHECK(rows == 2);
    CHECK(alleles == std::set<std::pair<std::string, std::string>>{
                         {"GA", "G"}, {"GAA", "G"}});
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
    struct Totals {
        int spanned = 0, in_gap_hets = 0, tagged = 0, scored = 0, correct = 0;
        int unused_clean_hets = 0;
        bool switched = false;
    };
    std::map<std::string, Totals> totals;
    // One arm. Recovery is not a mode any more: it runs inside every chunk's
    // own first solve, so there is nothing to switch on and nothing to compare
    // against but the competitor and the committed expectations.
    const std::vector<std::pair<std::string, std::string>> arms = {
        {"graph", ""},
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
            t.unused_clean_hets += got.unused_clean_hets;
            t.switched = t.switched || got.switched;
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
        CHECK(t.unused_clean_hets == 0);
        CHECK_FALSE(t.switched);
        CHECK((t.scored ? double(t.correct) / t.scored : 0.0) >= 0.95);
    }
}
