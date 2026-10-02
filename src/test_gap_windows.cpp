// Regression tests over the chr20 gap windows.
//
// These are integration tests, not unit tests: each one runs the pipeline binary
// on one window of the real test data, parses the phased VCF and phased BAM it
// produced, and scores the result against the parental truth. Owning-chunk
// replays are shared when multiple windows request identical pipeline inputs.
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
//   spans        -- does one phase set bracket the whole gap? An unexpected
//                   span needs physical read coverage and parental-orientation
//                   review before its expectation changes: chr20:48,176,830-
//                   48,229,446 was once reported CLOSED at 100% read accuracy
//                   while its halves sat on opposite haplotypes.
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
// Expectations are committed in src/test_gap_windows_expect.tsv. Coverage and
// accuracy use floors; spans are exact until a new join's orientation is
// reviewed. Regenerate with scripts/refresh_gap_window_expectations.sh only
// for an intended, measured improvement.

#define CATCH_CONFIG_MAIN
#include "../third_party/catch2/catch.hpp"

#include <htslib/sam.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <map>
#include <set>
#include <sstream>
#include <string>
#include <unordered_map>
#include <unordered_set>
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

/// Both boundaries identify a case: distinct gaps can start at one variant.
std::string window_key(const Window& w) {
    return std::to_string(w.gap_left) + "-" + std::to_string(w.gap_right);
}

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
    std::string diagnosis;  // span status and missing-link evidence; see explain_gap
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

/// Keyed by "<arm>\t<gap_left>-<gap_right>" so shared starts remain distinct.
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

/// Return the first changed reference base, matching vcf_to_variant_key().
/// VCF indels retain a shared left anchor, while the gap manifest and candidate
/// table use the unanchored BAM coordinate. Comparing raw VCF POS against those
/// boundaries made a correctly spanning deletion appear one base too short.
long long normalized_vcf_position(long long pos, const std::string& ref,
                                  const std::string& alt) {
    size_t shared = 0;
    const size_t limit = std::min(ref.size(), alt.size());
    while (shared < limit && ref[shared] == alt[shared]) ++shared;
    return pos + static_cast<long long>(shared);
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
    const auto it = required.find(arm + "\t" + window_key(w));
    if (it == required.end()) return;

    // Several alleles can share a normalized position. Preserve whether any
    // row at that position is phased instead of letting the last row win.
    std::map<long long, bool> found;
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
        const std::string& ps = f[static_cast<size_t>(ps_i)];
        found[std::stoll(f[static_cast<size_t>(pos_i)])] |=
            !ps.empty() && ps != "0" && ps != ".";
    }
    std::ostringstream detail;
    for (const RequiredSite& r : it->second) {
        bool retrieved = false;
        bool used = false;
        for (long long d = -2; d <= 2; ++d) {
            const auto hit = found.find(r.pos + d);
            if (hit == found.end()) continue;
            retrieved = true;
            used |= hit->second;
        }
        if (!retrieved) {
            ++out.required_missing;
            detail << "\n    NOT RETRIEVED " << r.pos << " " << r.ref << ">" << r.alt
                   << " (recorded as " << r.cat << ")";
        } else if (!used) {
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
        const long long raw_pos = std::stoll(f[1]);
        const long long pos = normalized_vcf_position(
            raw_pos, f[3], f[4]);
        // A VCF indel occupies both its mandatory anchor coordinate and its
        // canonical first-changed coordinate. Use that closed interval for
        // phase-block extent, while candidate membership below continues to
        // use the canonical coordinate. Reducing an insertion to only pos + 1
        // creates a false one-base gap at its left breakpoint.
        const long long record_beg = std::min(raw_pos, pos);
        const long long record_end = std::max(raw_pos, pos);
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
        if (it == extent.end()) extent[ps] = {record_beg, record_end};
        else {
            it->second.first = std::min(it->second.first, record_beg);
            it->second.second = std::max(it->second.second, record_end);
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
    static std::map<std::pair<long long, long long>, ReadSpans> cache;
    const auto key = std::make_pair(w.gap_left, w.gap_right);
    auto it = cache.find(key);
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
    return cache.emplace(key, std::move(spans)).first->second;
}

/// Run one arm over one window. Returns false when the binary failed, leaving
/// its stderr on disk for the failure message.
bool run_arm(const Paths& p, const Window& w, const std::string& arm,
             const std::string& flags, std::string& outdir) {
    std::ostringstream dir;
    dir << p.workdir << "/" << arm << "/w" << window_key(w);
    outdir = dir.str();
    std::ostringstream cmd;
    cmd << "mkdir -p '" << outdir << "' && '" << p.binary << "' collect-graph-variation"
        << " --ref '" << p.test_data << "/chm13v2.0.chr20.renamed.fa'"
        << " --bam '" << p.test_data
        << "/HG002_chr20_hifi_mapped_to_CHM13_chr20_annotated.bam'"
        << " --sites '" << p.test_data << "/chr20.sites.striped.vcf.gz'"
        << " --gaf '" << p.test_data << "/HG002.chr20.annotated.coord.gaf.gz'"
        // These graph source blocks begin beyond the 50-kb test padding.
        // Replay their owning chunk so the tested flank gauges match production.
        << " -r 'CHM13#0#chr20:"
        << ((w.gap_left == 48971192 || w.gap_left == 48929511)
                ? 48000001 :
            w.gap_left == 7047080 ? 7000001 :
            w.gap_left == 8166027 ? 8000001 :
            (w.gap_left == 514902 || w.gap_left == 528850) ? 1 :
            w.gap_left == 1180618 ? 1000001 :
            w.gap_left == 10727690 ? 10000001 :
            w.gap_left == 11573074 ? 11000001 :
            w.gap_left == 14446295 ? 14000001 :
            w.gap_left == 15351845 ? 15000001 :
            (w.gap_left == 53945086 || w.gap_left == 53947946)
                ? 53000001 :
            w.gap_left == 20707556 ? 20000001 :
            w.gap_left == 33211072 ? 33000001 :
            w.gap_left == 30673476 ? 30000001 :
            w.gap_left == 35328965 ? 35000001 :
            w.gap_left == 6513891 ? 6000001 :
            w.gap_left == 32431751 ? 32300001 :
            w.gap_left == 32234664 ? 32000001 :
            w.gap_left == 23421003 ? 23000001 :
            w.gap_left == 13784702 ? 13000001 :
            w.gap_left == 12256072 ? 12000001 :
            (w.gap_left == 19373922 || w.gap_left == 19395544 ||
             w.gap_left == 19403172) ? 19000001 :
            w.gap_left == 21594343 ? 21000001 :
            (w.gap_left == 36332599 || w.gap_left == 36611593)
                ? 36000001 :
            w.gap_left == 38331110 ? 38000001 :
            (w.gap_left == 39147805 || w.gap_left == 39848887)
                ? 39000001 :
            (w.gap_left == 64138752 || w.gap_left == 64144256 ||
             w.gap_left == 64140314) ? 64000001 :
            (w.gap_left == 56662188 || w.gap_left == 56064697)
                ? 56000001 :
            w.gap_left == 47003897 ? 47000001 :
            w.gap_left == 61757551 ? 61000001 :
            w.gap_left == 37458817 ? 37000001 :
            w.gap_left == 46727050 ? 46000001 :
            (w.gap_left == 54000001 || w.gap_left == 54483506 ||
             w.gap_left == 54684758)
                ? 54000001 :
            w.gap_left == 3000001 ? 3000001 :
            (w.gap_left == 21159070 || w.gap_left == 21179807 ||
             w.gap_left == 21377985 || w.gap_left == 21435750 ||
             w.gap_left == 21514518)
                ? 21000001 : w.gap_left - 50000)
        << "-" << ((w.gap_left == 48971192 || w.gap_left == 48929511)
                ? 49000000 :
            w.gap_left == 7047080 ? 8000000 :
            w.gap_left == 8166027 ? 9000000 :
            (w.gap_left == 514902 || w.gap_left == 528850) ? 1000000 :
            w.gap_left == 1180618 ? 2000000 :
            w.gap_left == 10727690 ? 11000000 :
            w.gap_left == 11573074 ? 12000000 :
            w.gap_left == 14446295 ? 15000000 :
            w.gap_left == 15351845 ? 16000000 :
            (w.gap_left == 53945086 || w.gap_left == 53947946)
                ? 54000000 :
            w.gap_left == 20707556 ? 21000000 :
            w.gap_left == 33211072 ? 34000000 :
            w.gap_left == 30673476 ? 31000000 :
            w.gap_left == 35328965 ? 36000000 :
            w.gap_left == 6513891 ? 7000000 :
            w.gap_left == 32431751 ? 32500000 :
            w.gap_left == 32234664 ? 33000000 :
            w.gap_left == 23421003 ? 24000000 :
            w.gap_left == 13784702 ? 14000000 :
            w.gap_left == 12256072 ? 13000000 :
            (w.gap_left == 19373922 || w.gap_left == 19395544 ||
             w.gap_left == 19403172) ? 20000000 :
            w.gap_left == 21594343 ? 22000000 :
            (w.gap_left == 36332599 || w.gap_left == 36611593)
                ? 37000000 :
            w.gap_left == 38331110 ? 39000000 :
            (w.gap_left == 39147805 || w.gap_left == 39848887)
                ? 40000000 :
            (w.gap_left == 64138752 || w.gap_left == 64144256 ||
             w.gap_left == 64140314) ? 65000000 :
            (w.gap_left == 56662188 || w.gap_left == 56064697)
                ? 57000000 :
            w.gap_left == 47003897 ? 48000000 :
            w.gap_left == 61757551 ? 62000000 :
            w.gap_left == 37458817 ? 38000000 :
            w.gap_left == 46727050 ? 47000000 :
            (w.gap_left == 54000001 || w.gap_left == 54483506 ||
             w.gap_left == 54684758)
                ? 55000000 :
            w.gap_left == 3000001 ? 4000000 :
            (w.gap_left == 21159070 || w.gap_left == 21179807 ||
             w.gap_left == 21377985 || w.gap_left == 21435750 ||
             w.gap_left == 21514518)
                ? 22000000 : w.gap_right + 50000) << "'"
        << " -t " << test_threads() << " " << flags
        << " -o '" << outdir << "/candidates.tsv'"
        << " --phased-vcf-out '" << outdir << "/native.vcf'"
        << " --phased-bam-out '" << outdir << "/phased.bam'"
        << " > '" << outdir << "/stdout.log' 2> '" << outdir << "/stderr.log'";
    // Several gap cases replay the same owning chunk with identical inputs.
    // Keep each case's expected output path while sharing the completed run.
    const std::string command = cmd.str();
    const std::string output_placeholder = "<output>";
    std::string cache_key = p.workdir + "\n" + command;
    for (size_t pos = 0; (pos = cache_key.find(outdir, pos)) !=
                          std::string::npos; pos += output_placeholder.size())
        cache_key.replace(pos, outdir.size(), output_placeholder);
    const auto link_outputs = [](const std::filesystem::path& source_dir,
                                 const std::filesystem::path& target_dir) {
        std::error_code error;
        std::filesystem::create_directories(target_dir, error);
        if (error) return false;
        for (const char* name : {"candidates.tsv", "native.vcf",
                                 "phased.bam", "stdout.log", "stderr.log"}) {
            const std::filesystem::path source = source_dir / name;
            const std::filesystem::path target = target_dir / name;
            if (source == target) continue;
            std::filesystem::remove(target, error);
            if (error) return false;
            std::filesystem::create_hard_link(source, target, error);
            if (error) return false;
        }
        return true;
    };
    static std::map<std::string, std::string> completed_runs;
    const auto cached = completed_runs.find(cache_key);
    if (cached != completed_runs.end())
        return link_outputs(cached->second, outdir);
    if (std::system(command.c_str()) != 0) return false;
    const std::filesystem::path cache_dir =
        std::filesystem::path(p.workdir) / ".replay-cache" /
        ("run" + std::to_string(completed_runs.size()));
    if (!link_outputs(outdir, cache_dir)) return false;
    completed_runs.emplace(std::move(cache_key), cache_dir.string());
    return true;
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
/// Diagnose the span result using read coverage and phased interior sites.
/// A failed span alone does not explain whether recovery has a molecular path:
///
///   UNCLOSEABLE -- some interior position is crossed by no read. No admission
///                  or linking change can help; a competitor spanning it is
///                  making a join its own reads do not support.
///   ADJACENT   -- no interior reference base exists; only a direct allele
///                  link can join the neighboring boundary sites.
///   NO SITES    -- reads cross, but the gap interior holds no phased het.
///                  A direct flank link or another site would be needed.
///   NOT LINKED  -- reads cross AND interior sites are phased, but they did not
///                  end up in one phase set. Linking or orientation.
///
/// Connected windows are reported first: an already joined pair may need no
/// interior heterozygote at all. Computed from existing output, without a
/// second pipeline run.
std::string explain_gap(const Paths& p, const Window& w, const Outcome& got) {
    if (w.gap_right <= w.gap_left + 1) {
        std::ostringstream o;
        o << "gap diagnosis: adjacent boundary sites; distinct phase sets"
             " touching the window " << got.blocks << " -- "
          << (got.spans
                  ? "CONNECTED (one phase set reaches both boundary positions)"
                  : "ADJACENT SITES (a direct allele link is needed)");
        return o.str();
    }
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
    if (got.spans)
        o << "CONNECTED (one phase set reaches both boundary positions)";
    else if (thinnest == 0)
        o << "UNCLOSEABLE (a position inside the gap is crossed by no read)";
    else if (got.in_gap_hets == 0)
        o << "NO SITES (reads cross it, but a direct flank link or another"
             " phased site is needed)";
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
    const std::string key = arm + "\t" + window_key(w);
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

    // Spanning is exact until molecule support and parental orientation of a
    // new join have been reviewed. Several longer gaps do have crossing reads.
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
    std::fprintf(out, "# spans is an EQUALITY per window until new molecule support and\n");
    std::fprintf(out, "# parental orientation are reviewed; TOTAL rows count spans.\n");
    std::fprintf(out, "# min_in_gap_hets is a floor on SITES\n");
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
            std::fprintf(out, "%s\t%s\t%d\t%d\t%.2f\t%.2f\n", arm.c_str(), window_key(w).c_str(),
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
                std::fprintf(req, "%s\t%s\t%lld\t%s\t%s\t%s\n", arm.c_str(), window_key(w).c_str(),
                             pos, site[1].c_str(), site[2].c_str(), cat.c_str());
            }
        }
    }
    std::fclose(req);
    WARN("wrote required sites to " << p.required);
}

TEST_CASE("anchored indel spans its VCF breakpoint",
          "[gap][representation][unit]") {
    const std::string path = "/tmp/pgphase-anchored-span-regression.vcf";
    {
        std::ofstream out(path);
        REQUIRE(out.good());
        out << "##fileformat=VCFv4.2\n"
            << "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n"
            << "chr20\t100\t.\tA\tATC\t.\tPASS\t.\tGT:PS\t0|1:100\n"
            << "chr20\t200\t.\tT\tC\t.\tPASS\t.\tGT:PS\t1|0:100\n";
    }

    Window window;
    window.gap_left = 100;
    window.gap_right = 200;
    Outcome outcome;
    parse_vcf(path, window, outcome);
    std::remove(path.c_str());

    CHECK(outcome.spans);
    CHECK(outcome.in_gap_hets == 1);
}

TEST_CASE("adjacent gap diagnosis has no interior coverage", "[gap][unit]") {
    Window window;
    window.gap_left = 528827;
    window.gap_right = 528828;
    Outcome outcome;
    outcome.blocks = 2;
    CHECK(explain_gap(Paths{}, window, outcome).find("ADJACENT SITES") !=
          std::string::npos);
    outcome.spans = true;
    outcome.blocks = 1;
    CHECK(explain_gap(Paths{}, window, outcome).find("CONNECTED") !=
          std::string::npos);
}

TEST_CASE("original MSA SNP dropout is retried before CIGAR backfill",
          "[gap][stitch-connectivity][msa]") {
    const Paths p = paths();
    if (!p.complete()) {
        WARN("gap-window tests need inputs that are absent: " << p.missing());
        SUCCEED("skipped: inputs absent");
        return;
    }
    // A short solve already connects this pair. Its owning chunk preserves the
    // two source gauges that leave 35/36 crossing reads without the MSA SNP.
    // CIGAR backfill restores paired calls after solving but cannot repair the
    // source HP/PS solution; admission must inspect the original SNP matrix.
    Window region;
    region.gap_left = 39050001;
    region.gap_right = 39950000;
    std::string dir;
    REQUIRE(run_arm(p, region, "msa_snp_dropout_39m", "", dir));
    const std::set<std::string> boundaries{
        "39838293", "39848887", "39856144", "39868811"};
    std::map<std::string, std::string> genotypes;
    std::map<std::string, std::string> phase_sets;
    long long snp_depth = 0;
    std::ifstream vcf(dir + "/native.vcf");
    REQUIRE(vcf.good());
    std::string line;
    while (std::getline(vcf, line)) {
        if (line.empty() || line[0] == '#') continue;
        const auto fields = split_tabs(line);
        if (fields.size() < 10 || boundaries.count(fields[1]) == 0) continue;
        const size_t separator = fields[9].rfind(':');
        REQUIRE(separator != std::string::npos);
        genotypes[fields[1]] = fields[9].substr(0, fields[9].find(':'));
        phase_sets[fields[1]] = fields[9].substr(separator + 1);
        if (fields[1] == "39856144") {
            const size_t depth = fields[7].find("DP=");
            REQUIRE(depth != std::string::npos);
            snp_depth = std::stoll(fields[7].substr(depth + 3));
        }
    }
    REQUIRE(phase_sets.size() == boundaries.size());
    for (const auto& [pos, ps] : phase_sets) {
        CHECK(ps == phase_sets.at("39838293"));
        CHECK((genotypes.at(pos) == "0|1" || genotypes.at(pos) == "1|0"));
    }
    // Independent observed alleles support the same DEL/SNP polarity. A PS
    // equality without this allele relation could conceal a wrong join.
    CHECK(genotypes.at("39848887") == genotypes.at("39856144"));
    CHECK(snp_depth >= 60);
    Window gap;
    gap.gap_left = 39848887;
    gap.gap_right = 39856144;
    Outcome reads;
    score_bam(dir + "/phased.bam", gap, load_truth(p.truth_map),
              input_read_spans(p, gap), reads);
    CHECK_FALSE(reads.switched);
    CHECK(reads.concordance() >= 0.98);
    CHECK(reads.separated() == 1.0);
}

TEST_CASE("conflicting sparse pairs retry with their observed binomial tail",
          "[gap][stitch-connectivity][msa][sparse-pair]") {
    const Paths p = paths();
    if (!p.complete()) {
        WARN("gap-window tests need inputs that are absent: " << p.missing());
        SUCCEED("skipped: inputs absent");
        return;
    }
    // Seven paired calls contain one opposite vote. Their 50:50 tail is
    // 0.0625, not the 0.0078125 probability of seven unanimous votes. The
    // owning chunk is needed to exercise this grouped solve's focused retry.
    Window gap;
    gap.gap_left = 54684758;
    gap.gap_right = 54702072;
    std::string dir;
    REQUIRE(run_arm(p, gap, "sparse_pair_54m", "", dir));
    using Key = std::array<std::string, 3>;
    const std::set<Key> expected{
        {"54585084", "G", "GCACACA"},
        {"54585084", "GCA", "G"},
        {"54590762", "C", "T"},
        {"54596124", "C", "T"},
        {"54607447", "G", "GGAAGGAAAGAAAGAAAGAAA"},
        {"54684758", "CCT", "C"},
        {"54702072", "G", "A"},
        {"54894127", "G", "C"},
        {"54912022", "G", "A"}};
    std::map<Key, std::pair<std::string, std::string>> sites;
    std::ifstream vcf(dir + "/native.vcf");
    REQUIRE(vcf.good());
    std::string line;
    while (std::getline(vcf, line)) {
        if (line.empty() || line[0] == '#') continue;
        const auto fields = split_tabs(line);
        if (fields.size() < 10) continue;
        const Key key{fields[1], fields[3], fields[4]};
        if (expected.count(key) == 0) continue;
        const size_t separator = fields[9].rfind(':');
        REQUIRE(separator != std::string::npos);
        sites[key] = {fields[9].substr(0, fields[9].find(':')),
                      fields[9].substr(separator + 1)};
    }
    REQUIRE(sites.size() == expected.size());
    const Key deletion{"54684758", "CCT", "C"};
    const Key snp{"54702072", "G", "A"};
    const auto& boundary = sites.at(deletion);
    for (const auto& [key, site] : sites) {
        CHECK((site.first == "0|1" || site.first == "1|0"));
        // The private context rows belong to the selected BAM block, even
        // outside the seam. Keeping only the boundaries truncates its evidence.
        CHECK(site.second == boundary.second);
        if (key == Key{"54585084", "G", "GCACACA"})
            CHECK(site.first != sites.at(snp).first);
        else
            CHECK(site.first == sites.at(snp).first);
    }
    Outcome reads;
    score_bam(dir + "/phased.bam", gap, load_truth(p.truth_map),
              input_read_spans(p, gap), reads);
    CHECK_FALSE(reads.switched);
    CHECK(reads.concordance() >= 0.98);
    CHECK(reads.separated() >= 0.83);
    CHECK(count_hap_allele_conflicts(dir + "/native.vcf") == 0);
}

TEST_CASE("a moved deletion does not certify a whole right-block join",
          "[gap][stitch-connectivity][msa][moved-deletion]") {
    const Paths p = paths();
    if (!p.complete()) {
        WARN("gap-window tests need inputs that are absent: " << p.missing());
        SUCCEED("skipped: inputs absent");
        return;
    }
    // Padding this replay gives exactly the owning 58--59 Mb chunk. In this
    // context the deletion has moved into the right PS while its original
    // source still has a weak cut. The right PS's own source certificate is
    // insufficient to join the left SNP through that moved row.
    Window replay;
    replay.gap_left = 58050001;
    replay.gap_right = 58950000;
    std::string dir;
    REQUIRE(run_arm(p, replay, "moved_deletion_58m", "", dir));
    using Key = std::array<std::string, 3>;
    const Key left{"58366458", "T", "C"};
    const Key deletion{"58385702", "CA", "C"};
    const Key right{"58391091", "C", "G"};
    std::map<Key, std::string> phase_sets;
    std::ifstream vcf(dir + "/native.vcf");
    REQUIRE(vcf.good());
    std::string line;
    while (std::getline(vcf, line)) {
        if (line.empty() || line[0] == '#') continue;
        const auto fields = split_tabs(line);
        if (fields.size() < 10) continue;
        const Key key{fields[1], fields[3], fields[4]};
        if (key != left && key != deletion && key != right) continue;
        const size_t separator = fields[9].rfind(':');
        REQUIRE(separator != std::string::npos);
        phase_sets.emplace(key, fields[9].substr(separator + 1));
    }
    REQUIRE(phase_sets.size() == 3);
    CHECK(phase_sets.at(left) != phase_sets.at(deletion));
    CHECK(phase_sets.at(deletion) == phase_sets.at(right));
    Window gap;
    gap.gap_left = 58366458;
    gap.gap_right = 58385702;
    Outcome reads;
    score_bam(dir + "/phased.bam", gap, load_truth(p.truth_map),
              input_read_spans(p, gap), reads);
    CHECK_FALSE(reads.switched);
    CHECK(reads.concordance() >= 0.98);
    CHECK(reads.separated() >= 0.39);
}

TEST_CASE("focused retry certifies its path after CIGAR backfill",
          "[gap][stitch-connectivity][msa][backfill-certificate]") {
    const Paths p = paths();
    if (!p.complete()) {
        WARN("gap-window tests need inputs that are absent: " << p.missing());
        SUCCEED("skipped: inputs absent");
        return;
    }
    // The raw focused MSA path has no cut, but its transferred backfilled
    // matrix has a cut at 4.778 Mb. Accepting the earlier certificate removes
    // the existing 4.767 Mb bridge in the owning chromosome chunk.
    Window region;
    region.gap_left = 4050001;
    region.gap_right = 4950000;
    std::string dir;
    REQUIRE(run_arm(p, region, "backfill_certificate_4m", "", dir));
    std::map<std::string, std::pair<std::string, std::string>> sites;
    std::ifstream vcf(dir + "/native.vcf");
    REQUIRE(vcf.good());
    std::string line;
    while (std::getline(vcf, line)) {
        if (line.empty() || line[0] == '#') continue;
        const auto fields = split_tabs(line);
        if (fields.size() < 10 ||
            (fields[1] != "4766928" && fields[1] != "4792960")) continue;
        const size_t separator = fields[9].rfind(':');
        REQUIRE(separator != std::string::npos);
        sites[fields[1]] = {fields[9].substr(0, fields[9].find(':')),
                           fields[9].substr(separator + 1)};
    }
    REQUIRE(sites.size() == 2);
    CHECK(sites.at("4766928").second == sites.at("4792960").second);
    CHECK(sites.at("4766928").first == sites.at("4792960").first);
    CHECK((sites.at("4766928").first == "0|1" ||
           sites.at("4766928").first == "1|0"));
    Window gap;
    gap.gap_left = 4766928;
    gap.gap_right = 4792960;
    Outcome reads;
    score_bam(dir + "/phased.bam", gap, load_truth(p.truth_map),
              input_read_spans(p, gap), reads);
    CHECK_FALSE(reads.switched);
    CHECK(reads.concordance() >= 0.98);
    // The owning-chunk replay retains the established 64% separation floor;
    // whole-chromosome stitching has a different read denominator here.
    CHECK(reads.separated() >= 0.64);
}

TEST_CASE("an earlier focused retry does not hide a later MSA dropout",
          "[gap][msa][retry-scheduling]") {
    const Paths p = paths();
    if (!p.complete()) {
        WARN("gap-window tests need inputs that are absent: " << p.missing());
        SUCCEED("skipped: inputs absent");
        return;
    }
    // The owning chunk groups several seams. Its first focused retry fails;
    // the later complementary deletion pair must still receive its own solve.
    Window replay;
    replay.gap_left = 57050001;
    replay.gap_right = 57950000;
    const std::string arm = "retry_scheduling_57m";
    const std::string output = p.workdir + "/" + arm + "/w" + window_key(replay);
    std::filesystem::create_directories(output);
    // Remove old debug files: the assertion must describe this invocation.
    for (const auto& entry : std::filesystem::directory_iterator(output)) {
        if (entry.path().filename().string().find("matrix.") == 0)
            std::filesystem::remove(entry.path());
    }
    std::string dir;
    REQUIRE(run_arm(p, replay, arm,
                    "--phase-matrix-dump '" + output + "/matrix'", dir));
    size_t focused_trials = 0;
    bool found_later_pair = false;
    for (const auto& entry : std::filesystem::directory_iterator(dir)) {
        const std::string name = entry.path().filename().string();
        if (name.find(".focused.") == std::string::npos ||
            name.find(".trial-source.tsv") == std::string::npos)
            continue;
        ++focused_trials;
        std::ifstream matrix(entry.path());
        REQUIRE(matrix.good());
        // Matrix keys retain the internal deletion coordinate and length.
        std::map<std::pair<std::string, std::string>,
                 std::array<std::string, 3>> deletions;
        std::string right_ps;
        std::string line;
        while (std::getline(matrix, line)) {
            const auto fields = split_tabs(line);
            if (fields.size() < 11 || fields[0] != "VAR") continue;
            if (fields[2] == "57854342" && fields[3] == "D")
                deletions[{fields[2], fields[6]}] =
                    {fields[8], fields[9], fields[10]};
            if (fields[2] == "57866713" && fields[3] == "X")
                right_ps = fields[8];
        }
        const auto one = deletions.find({"57854342", "1"});
        const auto two = deletions.find({"57854342", "2"});
        if (one == deletions.end() || two == deletions.end() || right_ps.empty())
            continue;
        found_later_pair = true;
        CHECK(one->second[0] == two->second[0]);
        CHECK(one->second[1] != one->second[2]);
        CHECK(one->second[1] == two->second[2]);
        CHECK(one->second[2] == two->second[1]);
        // Retrying restores calls but does not authorize a conflicting join.
        CHECK(one->second[0] != right_ps);
    }
    CHECK(focused_trials >= 2);
    CHECK(found_later_pair);
    Window gap;
    gap.gap_left = 57854341;
    gap.gap_right = 57866713;
    Outcome outcome;
    parse_vcf(dir + "/native.vcf", gap, outcome);
    CHECK_FALSE(outcome.spans);
    CHECK(count_hap_allele_conflicts(dir + "/native.vcf") == 0);
}

TEST_CASE("complete recovery MSA blocks retain the complex left flank",
          "[gap][msa-complete-transfer]") {
    const Paths p = paths();
    if (!p.complete()) {
        WARN("gap-window tests need inputs that are absent: " << p.missing());
        SUCCEED("skipped: inputs absent");
        return;
    }
    Window gap;
    gap.gap_left = 528850;
    gap.gap_right = 542052;
    std::string dir;
    REQUIRE(run_arm(p, gap, "msa_complete_transfer", "", dir));
    std::ifstream vcf(dir + "/native.vcf");
    REQUIRE(vcf.good());
    std::map<std::string, std::pair<std::string, std::string>> sites;
    std::string line;
    while (std::getline(vcf, line)) {
        if (line.empty() || line[0] == '#') continue;
        const auto fields = split_tabs(line);
        if (fields.size() < 10) continue;
        const std::string key = fields[1] + ":" + fields[3] + ">" + fields[4];
        if (key != "528728:C>A" && key != "528825:AAAAATATAT>A" &&
            key != "528827:A>T" && key != "528850:TGA>T" &&
            key != "542052:CT>C" && key != "545002:G>A")
            continue;
        const size_t colon = fields[9].rfind(':');
        REQUIRE(colon != std::string::npos);
        sites.emplace(key, std::make_pair(fields[9].substr(0, 3),
                                         fields[9].substr(colon + 1)));
    }
    REQUIRE(sites.size() == 6);
    const auto& left = sites.at("528850:TGA>T");
    CHECK(left.second != ".");
    for (const auto& entry : sites) {
        CHECK(is_phased_het(entry.second.first));
        CHECK(entry.second.second == left.second);
    }
    CHECK(left.first == sites.at("542052:CT>C").first);
    CHECK(left.first != sites.at("545002:G>A").first);
    CHECK(left.first != sites.at("528827:A>T").first);
    CHECK(sites.at("528825:AAAAATATAT>A").first == sites.at("528827:A>T").first);
    CHECK(sites.at("528728:C>A").first == sites.at("528827:A>T").first);
    Outcome outcome;
    parse_vcf(dir + "/native.vcf", gap, outcome);
    CHECK(outcome.spans);
    CHECK(count_hap_allele_conflicts(dir + "/native.vcf") == 0);
    const auto truth = load_truth(p.truth_map);
    const auto& spans = input_read_spans(p, gap);
    Outcome flank_reads;
    score_bam(dir + "/phased.bam", gap, truth, spans, flank_reads);
    CHECK_FALSE(flank_reads.switched);
    CHECK(flank_reads.concordance() >= 0.99);
    std::unordered_map<std::string, char> local_truth;
    for (const auto& entry : spans) {
        if (entry.second.second < gap.gap_left || entry.second.first > gap.gap_right)
            continue;
        const auto found = truth.find(entry.first);
        if (found != truth.end()) local_truth.emplace(*found);
    }
    Outcome local_reads;
    score_bam(dir + "/phased.bam", gap, local_truth, spans, local_reads);
    CHECK(local_reads.scored >= 123);
    CHECK(local_reads.correct >= 123);
    CHECK(local_reads.concordance() >= 1.0);
    CHECK(local_reads.discordant() == 0);
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
            DYNAMIC_SECTION(arm << " / window " << window_key(w)) {
                const auto key = arm + "\t" + window_key(w);
                const auto it = expect.find(key);
                if (it == expect.end())
                    FAIL("no expectation row for arm '" << arm << "' window "
                         << window_key(w) << " in " << p.expectations
                         << " -- run scripts/refresh_gap_window_expectations.sh");
                const Outcome got = measure(p, w, arm, flags, truth);
                check_against(w, arm, got, it->second);
                if (arm == "graph" && w.gap_left == 21159070) {
                    // The 21-22 Mb recovery group includes several seams. A
                    // broad MSA retry can add deletion observations while
                    // silently removing this established boundary allele.
                    std::ifstream vcf(p.workdir + "/graph/w" + window_key(w) +
                                      "/native.vcf");
                    REQUIRE(vcf.good());
                    std::set<std::string> sites;
                    std::string line;
                    while (std::getline(vcf, line)) {
                        if (line.empty() || line[0] == '#') continue;
                        const auto fields = split_tabs(line);
                        if (fields.size() < 10) continue;
                        const std::string key = fields[1] + ":" + fields[3] +
                            ">" + fields[4];
                        if (key == "21159070:C>T" ||
                            key == "21172487:CA>C" ||
                            key == "21172501:A>C") {
                            CHECK(fields[9].find('|') != std::string::npos);
                            const size_t colon = fields[9].rfind(':');
                            REQUIRE(colon != std::string::npos);
                            CHECK(fields[9].substr(colon + 1) != ".");
                            sites.insert(key);
                        }
                    }
                    CHECK(sites.size() == 3);
                }
                if (arm == "graph" && w.gap_left == 46727050) {
                    // A callable clean SNP pair already joins this block.
                    // A nearer noisy SNP cannot preempt that evidence.
                    std::ifstream vcf(p.workdir + "/graph/w" + window_key(w) +
                                      "/native.vcf");
                    REQUIRE(vcf.good());
                    std::map<std::string, std::pair<std::string, std::string>> sites;
                    std::string line;
                    while (std::getline(vcf, line)) {
                        if (line.empty() || line[0] == '#') continue;
                        const auto fields = split_tabs(line);
                        if (fields.size() < 10) continue;
                        const std::string key = fields[1] + ":" + fields[3] +
                            ">" + fields[4];
                        if (key != "46727050:A>G" &&
                            key != "46747598:C>CATATATAT" &&
                            key != "46748667:T>C")
                            continue;
                        const size_t colon = fields[9].rfind(':');
                        REQUIRE(colon != std::string::npos);
                        sites.emplace(key, std::make_pair(
                            fields[9].substr(0, 3), fields[9].substr(colon + 1)));
                    }
                    REQUIRE(sites.size() == 3);
                    const auto& left = sites.at("46727050:A>G");
                    const auto& right = sites.at("46748667:T>C");
                    CHECK(left.second != ".");
                    CHECK(left.second == right.second);
                    CHECK(left.first != right.first);
                    CHECK(sites.at("46747598:C>CATATATAT") == right);
                }
                if (arm == "graph" && w.gap_left == 37458817) {
                    // A noisy right SNP has an inconsistent link to its
                    // next clean SNP. The indel fallback must not use that
                    // noisy SNP to flip and absorb the whole right block.
                    std::ifstream vcf(p.workdir + "/graph/w" + window_key(w) +
                                      "/native.vcf");
                    REQUIRE(vcf.good());
                    std::map<std::string, std::string> phase_sets;
                    std::string line;
                    while (std::getline(vcf, line)) {
                        if (line.empty() || line[0] == '#') continue;
                        const auto fields = split_tabs(line);
                        if (fields.size() < 10) continue;
                        const std::string key = fields[1] + ":" + fields[3] +
                            ">" + fields[4];
                        if (key != "37458817:CA>C" &&
                            key != "37461999:CTT>C" &&
                            key != "37466820:T>A" &&
                            key != "37482721:G>A")
                            continue;
                        const size_t colon = fields[9].rfind(':');
                        REQUIRE(colon != std::string::npos);
                        phase_sets.emplace(key, fields[9].substr(colon + 1));
                    }
                    REQUIRE(phase_sets.size() == 4);
                    CHECK(phase_sets.at("37458817:CA>C") !=
                          phase_sets.at("37461999:CTT>C"));
                    CHECK(phase_sets.at("37461999:CTT>C") ==
                          phase_sets.at("37466820:T>A"));
                    CHECK(phase_sets.at("37466820:T>A") ==
                          phase_sets.at("37482721:G>A"));
                }
                if (arm == "graph" && w.gap_left == 61757551) {
                    // The demoted right SNP must agree with its downstream
                    // clean SNP before its physical bridge joins whole blocks.
                    // Also retain the previously closed left-hand gap.
                    std::ifstream vcf(p.workdir + "/graph/w" + window_key(w) +
                                      "/native.vcf");
                    REQUIRE(vcf.good());
                    std::map<std::string, std::pair<std::string, std::string>> sites;
                    std::string line;
                    while (std::getline(vcf, line)) {
                        if (line.empty() || line[0] == '#') continue;
                        const auto fields = split_tabs(line);
                        if (fields.size() < 10) continue;
                        const std::string key = fields[1] + ":" + fields[3] +
                            ">" + fields[4];
                        if (key != "61747506:ATT>A" &&
                            key != "61757551:G>A" &&
                            key != "61773799:A>G" &&
                            key != "61782778:C>T")
                            continue;
                        const size_t colon = fields[9].rfind(':');
                        REQUIRE(colon != std::string::npos);
                        sites.emplace(key, std::make_pair(
                            fields[9].substr(0, 3), fields[9].substr(colon + 1)));
                    }
                    REQUIRE(sites.size() == 4);
                    const auto& left = sites.at("61757551:G>A");
                    const auto& right = sites.at("61773799:A>G");
                    CHECK(left.second != ".");
                    CHECK(sites.at("61747506:ATT>A") == left);
                    CHECK(left == right);
                    CHECK(sites.at("61782778:C>T").second == right.second);
                    CHECK(sites.at("61782778:C>T").first != right.first);
                }
                if (arm == "graph" && w.gap_left == 54483506) {
                    // The clean SNP pair crosses the recovered insertion. Its
                    // relative phase and the earlier downstream join must
                    // both survive the owning-chunk stitch order.
                    std::ifstream vcf(p.workdir + "/graph/w" + window_key(w) +
                                      "/native.vcf");
                    REQUIRE(vcf.good());
                    std::map<std::string, std::pair<std::string, std::string>> sites;
                    std::string line;
                    while (std::getline(vcf, line)) {
                        if (line.empty() || line[0] == '#') continue;
                        const auto fields = split_tabs(line);
                        if (fields.size() < 10) continue;
                        const std::string key = fields[1] + ":" + fields[3] +
                            ">" + fields[4];
                        if (key != "54483506:T>G" && key != "54490229:G>C" &&
                            key != "54547514:T>A" && key != "54569072:C>T")
                            continue;
                        const size_t colon = fields[9].rfind(':');
                        REQUIRE(colon != std::string::npos);
                        sites.emplace(key, std::make_pair(
                            fields[9].substr(0, 3), fields[9].substr(colon + 1)));
                    }
                    REQUIRE(sites.size() == 4);
                    const auto& left = sites.at("54483506:T>G");
                    const auto& right = sites.at("54490229:G>C");
                    const auto& downstream = sites.at("54547514:T>A");
                    const auto& terminal = sites.at("54569072:C>T");
                    CHECK(left.second != ".");
                    CHECK(left.second == right.second);
                    CHECK(right.second == downstream.second);
                    CHECK(downstream.second == terminal.second);
                    CHECK(left.first != right.first);
                    CHECK(right.first == downstream.first);
                    CHECK(downstream.first != terminal.first);
                }
                if (arm == "graph" && w.gap_left == 61747506) {
                    // Ten physical SNP pairs agree with the stored alleles.
                    // A retry must carry their relation into the stitch rather
                    // than letting a stronger indel edge reverse the right PS.
                    std::ifstream vcf(p.workdir + "/graph/w" + window_key(w) +
                                      "/native.vcf");
                    REQUIRE(vcf.good());
                    std::map<std::string, std::pair<std::string, std::string>> sites;
                    std::string line;
                    while (std::getline(vcf, line)) {
                        if (line.empty() || line[0] == '#') continue;
                        const auto fields = split_tabs(line);
                        if (fields.size() < 10 || fields[3].size() != 1 ||
                            fields[4].size() != 1 ||
                            (fields[1] != "61757551" && fields[1] != "61773799"))
                            continue;
                        const auto colon = fields[9].rfind(':');
                        REQUIRE(colon != std::string::npos);
                        sites.emplace(fields[1], std::make_pair(
                            fields[9].substr(0, 3), fields[9].substr(colon + 1)));
                    }
                    REQUIRE(sites.size() == 2);
                    CHECK(sites.at("61757551").second == sites.at("61773799").second);
                    CHECK(sites.at("61757551").first == sites.at("61773799").first);
                }
                if (arm == "graph" && w.gap_left == 30673476) {
                    // Parental SNP calls put these adjacent blocks on opposite
                    // ALT haplotypes. The right orphan block must retain its
                    // internal orientation after joining the left read block.
                    std::ifstream vcf(p.workdir + "/graph/w" + window_key(w) +
                                      "/native.vcf");
                    REQUIRE(vcf.good());
                    std::map<std::string, std::pair<std::string, std::string>> sites;
                    std::string line;
                    while (std::getline(vcf, line)) {
                        if (line.empty() || line[0] == '#') continue;
                        const auto fields = split_tabs(line);
                        if (fields.size() < 10) continue;
                        const std::string key = fields[1] + ":" + fields[3] +
                            ">" + fields[4];
                        if (key != "30673476:G>C" &&
                            key != "30673709:T>A" &&
                            key != "30676069:T>C") continue;
                        const size_t colon = fields[9].rfind(':');
                        REQUIRE(colon != std::string::npos);
                        sites.emplace(key, std::make_pair(
                            fields[9].substr(0, 3), fields[9].substr(colon + 1)));
                    }
                    REQUIRE(sites.size() == 3);
                    const auto& left = sites.at("30673476:G>C");
                    const auto& right = sites.at("30673709:T>A");
                    const auto& terminal = sites.at("30676069:T>C");
                    CHECK(left.second == right.second);
                    CHECK(right.second == terminal.second);
                    CHECK(left.first != right.first);
                    CHECK(left.first == terminal.first);
                }
                if (arm == "graph" && w.gap_left == 64144256) {
                    // The BAM source has a weak cut before the left deletion.
                    // Close the two boundary alleles without importing the
                    // established left graph block into that local component.
                    std::ifstream vcf(p.workdir + "/graph/w" + window_key(w) +
                                      "/native.vcf");
                    REQUIRE(vcf.good());
                    std::string left_graph_ps, left_boundary_ps, right_boundary_ps;
                    std::string line;
                    while (std::getline(vcf, line)) {
                        if (line.empty() || line[0] == '#') continue;
                        const auto fields = split_tabs(line);
                        if (fields.size() < 10) continue;
                        const auto colon = fields[9].rfind(':');
                        if (colon == std::string::npos) continue;
                        const std::string ps = fields[9].substr(colon + 1);
                        if (fields[1] == "64118182" && fields[3] == "T" &&
                            fields[4] == "C") left_graph_ps = ps;
                        if (fields[1] == "64144256" && fields[3] == "CA" &&
                            fields[4] == "C") left_boundary_ps = ps;
                        if (fields[1] == "64144722" &&
                            fields[3] == "ATGGTGGGGG" && fields[4] == "A")
                            right_boundary_ps = ps;
                    }
                    REQUIRE(!left_graph_ps.empty());
                    REQUIRE(!left_boundary_ps.empty());
                    REQUIRE(!right_boundary_ps.empty());
                    CHECK(left_graph_ps != left_boundary_ps);
                    CHECK(left_boundary_ps == right_boundary_ps);
                }
                if (arm == "graph" &&
                    (w.gap_left == 19395544 || w.gap_left == 19403172)) {
                    // The two SNP bridges must preserve one ALT haplotype
                    // through the BAM-only source run and graph singleton.
                    std::ifstream vcf(p.workdir + "/graph/w" + window_key(w) +
                                      "/native.vcf");
                    REQUIRE(vcf.good());
                    std::map<std::string, std::pair<std::string, std::string>> sites;
                    std::string line;
                    while (std::getline(vcf, line)) {
                        if (line.empty() || line[0] == '#') continue;
                        const auto fields = split_tabs(line);
                        if (fields.size() < 10) continue;
                        const std::string key = fields[1] + ":" + fields[3] +
                            ">" + fields[4];
                        if (key != "19395544:T>C" &&
                            key != "19403172:T>C" &&
                            key != "19414720:G>A") continue;
                        const auto colon = fields[9].rfind(':');
                        REQUIRE(colon != std::string::npos);
                        sites.emplace(key, std::make_pair(
                            fields[9].substr(0, 3), fields[9].substr(colon + 1)));
                    }
                    REQUIRE(sites.count("19395544:T>C") == 1);
                    REQUIRE(sites.count("19403172:T>C") == 1);
                    REQUIRE(sites.count("19414720:G>A") == 1);
                    const auto& first = sites.at("19395544:T>C");
                    const auto& middle = sites.at("19403172:T>C");
                    const auto& last = sites.at("19414720:G>A");
                    CHECK(first.second == middle.second);
                    CHECK(middle.second == last.second);
                    CHECK(first.first == middle.first);
                    CHECK(middle.first == last.first);
                }
                if (arm == "graph" && w.gap_left == 1180618) {
                    // Read HP calls orient the short BAM insertion run, while
                    // the later graph insertion stays behind its weak edge.
                    std::ifstream vcf(p.workdir + "/graph/w" + window_key(w) +
                                      "/native.vcf");
                    REQUIRE(vcf.good());
                    std::map<std::string, std::pair<std::string, std::string>> sites;
                    std::string line;
                    while (std::getline(vcf, line)) {
                        if (line.empty() || line[0] == '#') continue;
                        const auto fields = split_tabs(line);
                        if (fields.size() < 10) continue;
                        const std::string key = fields[1] + ":" + fields[3] +
                            ">" + fields[4];
                        if (key != "1180618:T>A" &&
                            key != "1194189:T>TATC" &&
                            key != "1194227:TCACCAC>T" &&
                            key != "1194230:C>T" &&
                            key != "1194233:C>T" &&
                            fields[1] != "1196967")
                            continue;
                        const auto colon = fields[9].rfind(':');
                        REQUIRE(colon != std::string::npos);
                        sites.emplace(key, std::make_pair(
                            fields[9].substr(0, 3), fields[9].substr(colon + 1)));
                    }
                    REQUIRE(sites.size() == 6);
                    const auto& left = sites.at("1180618:T>A");
                    CHECK(sites.at("1194189:T>TATC") == left);
                    CHECK(sites.at("1194230:C>T") == left);
                    CHECK(sites.at("1194233:C>T") == left);
                    CHECK(sites.at("1194227:TCACCAC>T").second == left.second);
                    CHECK(sites.at("1194227:TCACCAC>T").first != left.first);
                    const auto graph = std::find_if(
                        sites.begin(), sites.end(),
                        [](const auto& item) {
                            return item.first.substr(0, 8) == "1196967:";
                        });
                    REQUIRE(graph != sites.end());
                    CHECK(graph->second.second != left.second);
                }
                if (arm == "graph" && w.gap_left == 12256072) {
                    // The complete BAM source resolves tied graph edges, but
                    // the physical SNP pair fixes the allele orientation.
                    // Pin the exact left SNP, first recovered deletion, and
                    // right SNP; a PS extent alone can hide a wrong join.
                    std::ifstream vcf(p.workdir + "/graph/w" + window_key(w) +
                                      "/native.vcf");
                    REQUIRE(vcf.good());
                    std::map<std::string, std::pair<std::string, std::string>> sites;
                    std::string line;
                    while (std::getline(vcf, line)) {
                        if (line.empty() || line[0] == '#') continue;
                        const auto fields = split_tabs(line);
                        if (fields.size() < 10) continue;
                        const std::string key = fields[1] + ":" + fields[3] +
                            ">" + fields[4];
                        if (key != "12243047:C>G" &&
                            key != "12256072:A>G" &&
                            key != "12269535:TA>T" &&
                            key != "12277078:G>T")
                            continue;
                        const auto colon = fields[9].rfind(':');
                        REQUIRE(colon != std::string::npos);
                        sites.emplace(key, std::make_pair(
                            fields[9].substr(0, 3), fields[9].substr(colon + 1)));
                    }
                    REQUIRE(sites.size() == 4);
                    const auto& left = sites.at("12256072:A>G");
                    CHECK(left.second != ".");
                    CHECK(sites.at("12243047:C>G").second == left.second);
                    CHECK(sites.at("12243047:C>G").first != left.first);
                    CHECK(sites.at("12269535:TA>T") == left);
                    CHECK(sites.at("12277078:G>T") == left);
                }
                if (arm == "graph" &&
                    (w.gap_left == 48929511 || w.gap_left == 56064697)) {
                    // Both joins depend on their owning 1-Mb chunk. Pin the
                    // boundary alleles so a distant PS extent cannot mask a
                    // switch or a missing physical variant at either end.
                    const std::string left_key = w.gap_left == 48929511
                        ? "48929511:GTGGA>G" : "56064697:G>A";
                    const std::string right_key = w.gap_left == 48929511
                        ? "48950388:C>T" : "56083708:TCA>T";
                    std::ifstream vcf(p.workdir + "/graph/w" + window_key(w) +
                                      "/native.vcf");
                    REQUIRE(vcf.good());
                    std::map<std::string, std::pair<std::string, std::string>> sites;
                    std::string line;
                    while (std::getline(vcf, line)) {
                        if (line.empty() || line[0] == '#') continue;
                        const auto fields = split_tabs(line);
                        if (fields.size() < 10) continue;
                        const std::string key = fields[1] + ":" + fields[3] +
                            ">" + fields[4];
                        if (key != left_key && key != right_key) continue;
                        const auto colon = fields[9].rfind(':');
                        REQUIRE(colon != std::string::npos);
                        sites.emplace(key, std::make_pair(
                            fields[9].substr(0, 3), fields[9].substr(colon + 1)));
                    }
                    REQUIRE(sites.size() == 2);
                    CHECK(sites.at(left_key).second != ".");
                    CHECK(sites.at(right_key) == sites.at(left_key));
                }
                if (arm == "graph" && w.gap_left == 58366458) {
                    // The two clean SNP ALTs belong to opposite parents. The
                    // intervening repeat deletion has nonsegregating CIGAR
                    // calls, so it cannot certify the joined block's gauge.
                    std::ifstream vcf(p.workdir + "/graph/w" + window_key(w) +
                                      "/native.vcf");
                    REQUIRE(vcf.good());
                    std::map<std::string, std::pair<std::string, std::string>> sites;
                    std::string line;
                    while (std::getline(vcf, line)) {
                        if (line.empty() || line[0] == '#') continue;
                        const auto fields = split_tabs(line);
                        if (fields.size() < 10) continue;
                        const std::string key = fields[1] + ":" + fields[3] +
                            ">" + fields[4];
                        if (key != "58366458:T>C" &&
                            key != "58391091:C>G") continue;
                        const auto colon = fields[9].rfind(':');
                        REQUIRE(colon != std::string::npos);
                        sites.emplace(key, std::make_pair(
                            fields[9].substr(0, 3), fields[9].substr(colon + 1)));
                    }
                    REQUIRE(sites.count("58366458:T>C") == 1);
                    REQUIRE(sites.count("58391091:C>G") == 1);
                    const auto& left = sites.at("58366458:T>C");
                    const auto& right = sites.at("58391091:C>G");
                    if (left.second == right.second)
                        CHECK(left.first != right.first);
                }
                if (arm == "graph" &&
                    (w.gap_left == 7047080 || w.gap_left == 48971192 ||
                     w.gap_left == 56662188)) {
                    // A block extent can span the interval without containing
                    // both intended boundary alleles. Assert their exact VCF
                    // keys share one phased PS in the owning production chunk.
                    const bool first = w.gap_left == 7047080;
                    const bool normalized_snp = w.gap_left == 56662188;
                    const std::string left_ref = first ? "C" : "G";
                    const std::string left_alt = first ? "CA" :
                        normalized_snp ? "T" : "GT";
                    const std::string right_ref = first || normalized_snp ? "C" : "T";
                    const std::string right_alt = first || normalized_snp ? "T" : "C";
                    std::ifstream vcf(p.workdir + "/graph/w" + window_key(w) +
                                      "/native.vcf");
                    REQUIRE(vcf.good());
                    std::string left_ps, right_ps, line;
                    while (std::getline(vcf, line)) {
                        if (line.empty() || line[0] == '#') continue;
                        const auto fields = split_tabs(line);
                        if (fields.size() < 10) continue;
                        const bool left = fields[1] ==
                                std::to_string(w.gap_left) &&
                            fields[3] == left_ref && fields[4] == left_alt;
                        const bool right = fields[1] ==
                                std::to_string(w.gap_right) &&
                            fields[3] == right_ref && fields[4] == right_alt;
                        if (!left && !right) continue;
                        CHECK(is_phased_het(fields[9].substr(0,
                            fields[9].find(':'))));
                        std::istringstream formats(fields[8]), sample(fields[9]);
                        std::string format, value;
                        while (std::getline(formats, format, ':') &&
                               std::getline(sample, value, ':'))
                            if (format == "PS") {
                                if (left) left_ps = value;
                                if (right) right_ps = value;
                                break;
                            }
                    }
                    REQUIRE((!left_ps.empty() && left_ps != "." && left_ps != "0"));
                    REQUIRE((!right_ps.empty() && right_ps != "." && right_ps != "0"));
                    CHECK(left_ps == right_ps);
                }
                if (arm == "graph" && w.gap_left == 60033052) {
                    // This graph REF/ALT SNP has no callable REF base in BAM:
                    // physical reads carry ALT or a deletion. It must not
                    // orient either neighboring phase block.
                    std::ifstream vcf(p.workdir + "/graph/w" + window_key(w) + "/native.vcf");
                    REQUIRE(vcf.good());
                    std::string line;
                    while (std::getline(vcf, line)) {
                        if (line.empty() || line[0] == '#') continue;
                        const auto fields = split_tabs(line);
                        REQUIRE(fields.size() >= 2);
                        CHECK(fields[1] != "60033350");
                    }
                }
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

TEST_CASE("verified MSA insertion survives graph output classification",
          "[gap][representation]") {
    const Paths p = paths();
    if (!p.complete()) {
        WARN("gap-window tests need inputs that are absent: " << p.missing());
        SUCCEED("skipped: inputs absent");
        return;
    }
    // The 50-kb padding makes this exactly the 35–36 Mb owning chunk.
    // The narrow panel replay keeps its existing accuracy floor and span gate.
    Window w;
    w.gap_left = 35050001;
    w.gap_right = 35950000;
    std::string dir;
    REQUIRE(run_arm(p, w, "graph", "", dir));

    std::ifstream vcf(dir + "/native.vcf");
    REQUIRE(vcf.good());
    std::map<std::string, std::string> phase_sets;
    std::string line;
    while (std::getline(vcf, line)) {
        if (line.empty() || line[0] == '#') continue;
        const auto fields = split_tabs(line);
        if (fields.size() < 10) continue;
        const std::string key = fields[1] + ":" + fields[3] +
            ">" + fields[4];
        if (key != "35328965:A>AT" &&
            key != "35342608:G>GAGATAGAT" &&
            key != "35347817:CT>C")
            continue;
        const size_t colon = fields[9].rfind(':');
        REQUIRE(colon != std::string::npos);
        CHECK(fields[9].find('|') != std::string::npos);
        phase_sets.emplace(key, fields[9].substr(colon + 1));
    }
    REQUIRE(phase_sets.size() == 3);
    CHECK(phase_sets.at("35342608:G>GAGATAGAT") ==
          phase_sets.at("35347817:CT>C"));
    // The graph repeat allele independently links both imported BAM blocks.
    CHECK(phase_sets.at("35328965:A>AT") ==
          phase_sets.at("35342608:G>GAGATAGAT"));
}

TEST_CASE("a phased BAM deletion survives an unphased graph duplicate",
          "[gap][representation]") {
    const Paths p = paths();
    if (!p.complete()) {
        WARN("gap-window tests need inputs that are absent: " << p.missing());
        SUCCEED("skipped: inputs absent");
        return;
    }
    Window w;
    w.gap_left = 20707556;
    w.gap_right = 20731688;
    std::string dir;
    REQUIRE(run_arm(p, w, "graph", "", dir));

    std::ifstream in(dir + "/native.vcf");
    REQUIRE(in.good());
    int recovered_rows = 0;
    std::string line;
    while (std::getline(in, line)) {
        if (line.empty() || line[0] == '#') continue;
        const auto fields = split_tabs(line);
        if (fields.size() < 10 || fields[1] != "20711883" ||
            fields[3] != "CA" || fields[4] != "C")
            continue;
        ++recovered_rows;
        const std::string& sample = fields[9];
        CHECK((sample.substr(0, 3) == "0|1" ||
               sample.substr(0, 3) == "1|0"));
        const auto colon = sample.rfind(':');
        REQUIRE(colon != std::string::npos);
        CHECK(sample.substr(colon + 1) != "0");
        CHECK(sample.substr(colon + 1) != ".");
    }
    CHECK(recovered_rows == 1);
}

TEST_CASE("graph SNP retry preserves a verified BAM deletion connection",
          "[gap][stitch-confidence]") {
    const Paths p = paths();
    if (!p.complete()) {
        WARN("gap-window tests need inputs that are absent: " << p.missing());
        SUCCEED("skipped: inputs absent");
        return;
    }
    // This existing window replays the complete 15-Mb owning chunk. New
    // graph-allele admission checks must retain the BAM-only source connection.
    Window w;
    w.gap_left = 15351845;
    w.gap_right = 15367755;
    std::string dir;
    REQUIRE(run_arm(p, w, "graph", "", dir));
    std::ifstream vcf(dir + "/native.vcf");
    REQUIRE(vcf.good());
    std::map<std::string, std::pair<std::string, std::string>> sites;
    std::map<std::string, std::pair<std::string, std::string>> bridge;
    std::string line;
    while (std::getline(vcf, line)) {
        if (line.empty() || line[0] == '#') continue;
        const auto fields = split_tabs(line);
        if (fields.size() < 10) continue;
        const std::string key = fields[1] + ":" + fields[3] + ">" + fields[4];
        if (key == "15329499:A>T" || key == "15351845:CT>C" ||
            key == "15367755:A>AATCT" ||
            key == "15380509:C>T") {
            const auto colon = fields[9].rfind(':');
            REQUIRE(colon != std::string::npos);
            bridge.emplace(fields[1], std::make_pair(
                fields[9].substr(0, fields[9].find(':')),
                fields[9].substr(colon + 1)));
        }
        if (key != "15039543:A>G" && key != "15056025:G>A" &&
            key != "15071132:ATGTATACACACACGTGTGTGTATACGCACACCTGTG>A" &&
            key != "15100456:CT>C")
            continue;
        const auto colon = fields[9].rfind(':');
        REQUIRE(colon != std::string::npos);
        sites.emplace(fields[1], std::make_pair(
            fields[9].substr(0, 3), fields[9].substr(colon + 1)));
    }
    REQUIRE(sites.size() == 4);
    CHECK(sites.at("15039543").second == sites.at("15056025").second);
    CHECK(sites.at("15056025").second == sites.at("15071132").second);
    CHECK(sites.at("15071132").second == sites.at("15100456").second);
    CHECK(sites.at("15039543").first == sites.at("15056025").first);
    CHECK(sites.at("15071132").first == sites.at("15100456").first);
    // The MSA deletion/ATCT insertion must be on opposite haplotypes. A shared
    // label alone would also pass after the upstream double-flip defect.
    REQUIRE(bridge.size() == 4);
    CHECK(is_phased_het(bridge.at("15329499").first));
    CHECK(is_phased_het(bridge.at("15351845").first));
    CHECK(is_phased_het(bridge.at("15367755").first));
    CHECK(is_phased_het(bridge.at("15380509").first));
    CHECK(bridge.at("15351845").second == bridge.at("15329499").second);
    CHECK(bridge.at("15367755").second == bridge.at("15351845").second);
    CHECK(bridge.at("15380509").second == bridge.at("15351845").second);
    CHECK(bridge.at("15351845").first != bridge.at("15367755").first);
    CHECK(bridge.at("15367755").first == bridge.at("15380509").first);
    CHECK(bridge.at("15329499").first == bridge.at("15367755").first);

    const auto truth = load_truth(p.truth_map);
    const auto& spans = input_read_spans(p, w);
    Outcome flank_reads;
    score_bam(dir + "/phased.bam", w, truth, spans, flank_reads);
    CHECK_FALSE(flank_reads.switched);
    std::unordered_map<std::string, char> gap_truth;
    for (const auto& entry : spans) {
        if (entry.second.second < w.gap_left || entry.second.first > w.gap_right)
            continue;
        const auto found = truth.find(entry.first);
        if (found != truth.end()) gap_truth.emplace(*found);
    }
    Outcome local_reads;
    score_bam(dir + "/phased.bam", w, gap_truth, spans, local_reads);
    CHECK(local_reads.scored >= 92);
    CHECK(local_reads.correct >= 89);
    CHECK(local_reads.concordance() >= 0.96);
    CHECK(local_reads.scored - local_reads.correct <= 3);
}

TEST_CASE("an unsupported BAM source cut keeps the far deletion independent",
          "[gap][stitch-connectivity]") {
    const Paths p = paths();
    if (!p.complete()) {
        WARN("gap-window tests need inputs that are absent: " << p.missing());
        SUCCEED("skipped: inputs absent");
        return;
    }
    // The BAM sub-solve assigns both sites one source PS, but no read calls
    // both alleles. The downstream deletion is supported by the right block.
    Window w;
    w.gap_left = 33211072;
    w.gap_right = 33227050;
    std::string dir;
    REQUIRE(run_arm(p, w, "weak_source_cut", "", dir));

    std::ifstream vcf(dir + "/native.vcf");
    REQUIRE(vcf.good());
    std::map<std::string, std::pair<std::string, std::string>> sites;
    std::string line;
    while (std::getline(vcf, line)) {
        if (line.empty() || line[0] == '#') continue;
        const auto fields = split_tabs(line);
        if (fields.size() < 10 ||
            (fields[1] != "33211072" && fields[1] != "33227050" &&
             fields[1] != "33233843"))
            continue;
        const std::string& sample = fields[9];
        const size_t colon = sample.rfind(':');
        REQUIRE(colon != std::string::npos);
        sites.emplace(fields[1], std::make_pair(
            sample.substr(0, sample.find(':')), sample.substr(colon + 1)));
    }
    REQUIRE(sites.size() == 3);
    CHECK(sites.at("33211072").second != sites.at("33227050").second);
    CHECK(sites.at("33227050").second == sites.at("33233843").second);
    CHECK(sites.at("33227050").first == sites.at("33233843").first);
}

TEST_CASE("graph SNPs bridge supported recovery blocks", "[gap][graph-bridge]") {
    const Paths p = paths();
    if (!p.complete()) {
        WARN("gap-window tests need inputs that are absent: " << p.missing());
        SUCCEED("skipped: inputs absent");
        return;
    }
    // run_arm adds 50 kb on each side; this selects the audited
    // chr20:38,233,000-38,343,000 region.
    Window w;
    w.gap_left = 38283000;
    w.gap_right = 38293000;
    std::string dir;
    REQUIRE(run_arm(p, w, "graph_bridge", "", dir));

    std::ifstream in(dir + "/native.vcf");
    REQUIRE(in.good());
    std::set<std::string> phase_sets;
    std::map<std::string, std::string> snp_genotypes;
    int boundary_deletions = 0;
    std::string line;
    while (std::getline(in, line)) {
        if (line.empty() || line[0] == '#') continue;
        const auto fields = split_tabs(line);
        if (fields.size() < 10) continue;
        const bool left_snp = fields[1] == "38259286";
        const bool right_snp = fields[1] == "38279170";
        const bool boundary_deletion = fields[1] == "38283561" &&
                                       fields[3].size() > fields[4].size();
        if (!left_snp && !right_snp && !boundary_deletion) continue;
        const std::string& sample = fields[9];
        const size_t separator = sample.rfind(':');
        REQUIRE(separator != std::string::npos);
        phase_sets.insert(sample.substr(separator + 1));
        if (left_snp || right_snp)
            snp_genotypes[fields[1]] = sample.substr(0, 3);
        if (boundary_deletion) ++boundary_deletions;
    }
    REQUIRE(snp_genotypes.size() == 2);
    CHECK(snp_genotypes.at("38259286") !=
          snp_genotypes.at("38279170"));
    CHECK(boundary_deletions == 2);
    REQUIRE(phase_sets.size() == 1);
    CHECK(*phase_sets.begin() != ".");
}

TEST_CASE("recovery stitch preserves the next flank across an unlinked seam",
          "[gap][stitch-connectivity]") {
    const Paths p = paths();
    if (!p.complete()) {
        WARN("gap-window tests need inputs that are absent: " << p.missing());
        SUCCEED("skipped: inputs absent");
        return;
    }
    // At 51.235-51.262 Mb, high-MAPQ reads cover both clean SNPs but no
    // molecule observes both. The next 8.7 kb has 28 clean-SNP allele links.
    Window w;
    w.gap_left = 51235063;
    w.gap_right = 51262081;
    std::string dir;
    REQUIRE(run_arm(p, w, "stitch_connectivity", "", dir));

    std::ifstream in(dir + "/native.vcf");
    REQUIRE(in.good());
    std::map<std::string, std::string> phase_sets;
    std::string line;
    while (std::getline(in, line)) {
        if (line.empty() || line[0] == '#') continue;
        const auto fields = split_tabs(line);
        if (fields.size() < 10) continue;
        if (fields[1] != "51235063" && fields[1] != "51262081" &&
            fields[1] != "51270774") continue;
        const size_t separator = fields[9].rfind(':');
        REQUIRE(separator != std::string::npos);
        phase_sets[fields[1]] = fields[9].substr(separator + 1);
    }
    REQUIRE(phase_sets.size() == 3);
    CHECK(phase_sets.at("51235063") != ".");
    CHECK(phase_sets.at("51262081") != ".");
    CHECK(phase_sets.at("51270774") != ".");
    CHECK(phase_sets.at("51235063") != phase_sets.at("51262081"));
    CHECK(phase_sets.at("51262081") == phase_sets.at("51270774"));

    // The production 1 Mb chunk has another recovery seam to the right;
    // its larger transaction must not restore the unsupported left join.
    Window full_chunk;
    full_chunk.gap_left = 51050001;
    full_chunk.gap_right = 51950000;
    REQUIRE(run_arm(p, full_chunk, "stitch_connectivity_full", "", dir));
    std::ifstream full_vcf(dir + "/native.vcf");
    REQUIRE(full_vcf.good());
    std::map<std::string, std::string> full_phase_sets;
    while (std::getline(full_vcf, line)) {
        if (line.empty() || line[0] == '#') continue;
        const auto fields = split_tabs(line);
        if (fields.size() < 10 ||
            (fields[1] != "51235063" && fields[1] != "51262081" &&
             fields[1] != "51270774" && fields[1] != "51286463" &&
             fields[1] != "51287372"))
            continue;
        const size_t separator = fields[9].rfind(':');
        REQUIRE(separator != std::string::npos);
        full_phase_sets[fields[1]] = fields[9].substr(separator + 1);
    }
    REQUIRE(full_phase_sets.size() == 5);
    CHECK(full_phase_sets.at("51235063") !=
          full_phase_sets.at("51262081"));
    CHECK(full_phase_sets.at("51262081") ==
          full_phase_sets.at("51270774"));
    CHECK(full_phase_sets.at("51270774") ==
          full_phase_sets.at("51286463"));
    CHECK(full_phase_sets.at("51286463") ==
          full_phase_sets.at("51287372"));

    // The BAM source has a weak cut before its 19.397 Mb local run.
    // Keep the earlier prefix separate while the two clean SNP bridges
    // connect that run to the graph singleton in the correct orientation.
    Window weak_source;
    weak_source.gap_left = 19050000;
    weak_source.gap_right = 19950000;
    REQUIRE(run_arm(p, weak_source, "stitch_weak_source", "", dir));
    std::ifstream weak_vcf(dir + "/native.vcf");
    REQUIRE(weak_vcf.good());
    std::map<std::string, std::pair<std::string, std::string>> weak_sites;
    while (std::getline(weak_vcf, line)) {
        if (line.empty() || line[0] == '#') continue;
        const auto fields = split_tabs(line);
        if (fields.size() < 10) continue;
        const std::string key = fields[1] + ":" + fields[3] +
            ">" + fields[4];
        if (key != "19377345:T>TATATATAGAG" &&
            key != "19395544:T>C" && key != "19403172:T>C" &&
            key != "19414720:G>A") continue;
        const size_t separator = fields[9].rfind(':');
        REQUIRE(separator != std::string::npos);
        weak_sites[key] = {fields[9].substr(0, 3),
                           fields[9].substr(separator + 1)};
    }
    REQUIRE(weak_sites.size() == 4);
    const auto& prefix = weak_sites.at("19377345:T>TATATATAGAG");
    const auto& first = weak_sites.at("19395544:T>C");
    const auto& middle = weak_sites.at("19403172:T>C");
    const auto& last = weak_sites.at("19414720:G>A");
    CHECK(prefix.second != first.second);
    CHECK(first.second == middle.second);
    CHECK(middle.second == last.second);
    CHECK(first.first == middle.first);
    CHECK(middle.first == last.first);
}

TEST_CASE("a clean-SNP gap bridge keeps the prior graph gauge",
          "[gap][stitch-connectivity]") {
    const Paths p = paths();
    if (!p.complete()) {
        WARN("gap-window tests need inputs that are absent: " << p.missing());
        SUCCEED("skipped: inputs absent");
        return;
    }
    // Match the 47–48 Mb chromosome chunk. A 47.636-Mb source cut has two
    // quality-40 SNP molecules, but its left graph block was already flipped
    // at the preceding seam. Reusing the saved gauge without translating that
    // flip joins 2,200-plus reads to the opposite parent.
    Window w;
    w.gap_left = 47050001;
    w.gap_right = 47950000;
    std::string dir;
    REQUIRE(run_arm(p, w, "quality_snp_bridge_chunk", "", dir));

    std::ifstream vcf(dir + "/native.vcf");
    REQUIRE(vcf.good());
    const std::set<std::string> positions = {
        "47018892", "47118148", "47636774", "47659910"};
    std::map<std::string, std::pair<std::string, std::string>> sites;
    std::string line;
    while (std::getline(vcf, line)) {
        if (line.empty() || line[0] == '#') continue;
        const auto fields = split_tabs(line);
        if (fields.size() < 10 || positions.count(fields[1]) == 0) continue;
        const std::string gt = fields[9].substr(0, fields[9].find(':'));
        const size_t separator = fields[9].rfind(':');
        REQUIRE(separator != std::string::npos);
        sites[fields[1]] = {gt, fields[9].substr(separator + 1)};
    }
    REQUIRE(sites.size() == positions.size());
    const std::string& phase_set = sites.at("47018892").second;
    for (const std::string& pos : positions)
        CHECK(sites.at(pos).second == phase_set);
    CHECK(sites.at("47018892").first != sites.at("47118148").first);
    CHECK(sites.at("47636774").first != sites.at("47659910").first);

    samFile* bam = sam_open((dir + "/phased.bam").c_str(), "r");
    REQUIRE(bam != nullptr);
    bam_hdr_t* header = sam_hdr_read(bam);
    REQUIRE(header != nullptr);
    bam1_t* record = bam_init1();
    const auto truth = load_truth(p.truth_map);
    std::unordered_set<std::string> seen;
    int maternal_on_hap1 = 0;
    int paternal_on_hap1 = 0;
    while (sam_read1(bam, header, record) >= 0) {
        const std::string qname = bam_get_qname(record);
        const auto known = truth.find(qname);
        if (known == truth.end()) continue;
        const uint8_t* ps = bam_aux_get(record, "PS");
        const uint8_t* hp = bam_aux_get(record, "HP");
        if (ps == nullptr || hp == nullptr ||
            std::to_string(bam_aux2i(ps)) != phase_set)
            continue;
        const int hap = bam_aux2i(hp);
        if ((hap != 1 && hap != 2) || !seen.insert(qname).second) continue;
        const bool maternal = (hap == 1) == (known->second == 'M');
        if (maternal) ++maternal_on_hap1;
        else ++paternal_on_hap1;
    }
    bam_destroy1(record);
    bam_hdr_destroy(header);
    sam_close(bam);
    const int total = maternal_on_hap1 + paternal_on_hap1;
    CHECK(total >= 2000);
    // The physically certified 47 Mb join changes one local truth assignment;
    // 2,967 of 2,998 reads still keep the established graph gauge.
    CHECK(std::max(maternal_on_hap1, paternal_on_hap1) >=
          0.989 * static_cast<double>(total));
}

TEST_CASE("a local source path does not lose its graph bridge",
          "[gap][stitch-connectivity]") {
    const Paths p = paths();
    if (!p.complete()) {
        WARN("gap-window tests need inputs that are absent: " << p.missing());
        SUCCEED("skipped: inputs absent");
        return;
    }
    Window w;
    w.gap_left = 56323427;
    w.gap_right = 56343002;
    std::string dir;
    REQUIRE(run_arm(p, w, "stitch_local_source_span", "", dir));

    std::ifstream vcf(dir + "/native.vcf");
    REQUIRE(vcf.good());
    std::map<std::string, std::pair<std::string, std::string>> sites;
    std::string line;
    while (std::getline(vcf, line)) {
        if (line.empty() || line[0] == '#') continue;
        const auto fields = split_tabs(line);
        if (fields.size() < 10 ||
            (fields[1] != "56323427" && fields[1] != "56343002"))
            continue;
        const size_t separator = fields[9].rfind(':');
        REQUIRE(separator != std::string::npos);
        sites[fields[1]] = {fields[9].substr(0, 3),
                            fields[9].substr(separator + 1)};
    }
    REQUIRE(sites.size() == 2);
    CHECK(sites.at("56323427").second == sites.at("56343002").second);
    CHECK(sites.at("56323427").first == sites.at("56343002").first);

    Outcome reads;
    score_bam(dir + "/phased.bam", w, load_truth(p.truth_map),
              input_read_spans(p, w), reads);
    CHECK_FALSE(reads.switched);
    CHECK(reads.concordance() >= 0.99);
}

TEST_CASE("corroborated SNP molecules join sparse graph seams",
          "[gap][stitch-connectivity]") {
    const Paths p = paths();
    if (!p.complete()) {
        WARN("gap-window tests need inputs that are absent: " << p.missing());
        SUCCEED("skipped: inputs absent");
        return;
    }
    struct Bridge { long long left; long long right; bool same_allele; };
    const std::array<Bridge, 3> bridges{{
        {20875468, 20895876, true},
        {24357615, 24379633, false},
        {48022132, 48043584, true},
    }};
    for (const Bridge& bridge : bridges) {
        INFO("boundary " << bridge.left << " -> " << bridge.right);
        Window w;
        w.gap_left = bridge.left;
        w.gap_right = bridge.right;
        std::string dir;
        REQUIRE(run_arm(p, w, "stitch_corroborated_snp", "", dir));
        std::ifstream vcf(dir + "/native.vcf");
        REQUIRE(vcf.good());
        std::map<long long, std::pair<std::string, std::string>> sites;
        std::string line;
        while (std::getline(vcf, line)) {
            if (line.empty() || line[0] == '#') continue;
            const auto fields = split_tabs(line);
            if (fields.size() < 10) continue;
            const long long pos = std::stoll(fields[1]);
            if (pos != bridge.left && pos != bridge.right) continue;
            const size_t separator = fields[9].rfind(':');
            REQUIRE(separator != std::string::npos);
            sites[pos] = {fields[9].substr(0, 3),
                          fields[9].substr(separator + 1)};
        }
        REQUIRE(sites.size() == 2);
        CHECK(sites.at(bridge.left).second == sites.at(bridge.right).second);
        CHECK((sites.at(bridge.left).first == sites.at(bridge.right).first) ==
              bridge.same_allele);
        Outcome reads;
        score_bam(dir + "/phased.bam", w, load_truth(p.truth_map),
                  input_read_spans(p, w), reads);
        CHECK_FALSE(reads.switched);
        CHECK(reads.concordance() >= 0.99);
    }
}

TEST_CASE("an indel boundary with allele dropout gets an MSA retry",
          "[gap][stitch-connectivity]") {
    const Paths p = paths();
    if (!p.complete()) {
        WARN("gap-window tests need inputs that are absent: " << p.missing());
        SUCCEED("skipped: inputs absent");
        return;
    }
    Window w;
    w.gap_left = 24581764;
    w.gap_right = 24601281;
    std::string dir;
    REQUIRE(run_arm(p, w, "stitch_indel_dropout", "", dir));
    std::ifstream vcf(dir + "/native.vcf");
    REQUIRE(vcf.good());
    std::map<std::string, std::string> phase_sets;
    std::string line;
    while (std::getline(vcf, line)) {
        if (line.empty() || line[0] == '#') continue;
        const auto fields = split_tabs(line);
        if (fields.size() < 10 ||
            (fields[1] != "24581764" && fields[1] != "24601281"))
            continue;
        const size_t separator = fields[9].rfind(':');
        REQUIRE(separator != std::string::npos);
        phase_sets[fields[1]] = fields[9].substr(separator + 1);
    }
    REQUIRE(phase_sets.size() == 2);
    CHECK(phase_sets.at("24581764") == phase_sets.at("24601281"));
    Outcome reads;
    score_bam(dir + "/phased.bam", w, load_truth(p.truth_map),
              input_read_spans(p, w), reads);
    CHECK_FALSE(reads.switched);
    CHECK(reads.concordance() >= 0.99);
}

TEST_CASE("a lone clean SNP pair cannot join whole blocks",
          "[gap][stitch-connectivity]") {
    const Paths p = paths();
    if (!p.complete()) {
        WARN("gap-window tests need inputs that are absent: " << p.missing());
        SUCCEED("skipped: inputs absent");
        return;
    }
    // One MAPQ-60, Q40 molecule crosses these two SNPs. Each is the only
    // clean SNP it observes in its block, and joining them reverses 731 reads.
    Window w;
    w.gap_left = 62050001;
    w.gap_right = 62950000;
    std::string dir;
    REQUIRE(run_arm(p, w, "stitch_lone_snp_pair", "", dir));
    std::ifstream vcf(dir + "/native.vcf");
    REQUIRE(vcf.good());
    std::map<std::string, std::string> phase_sets;
    std::string line;
    while (std::getline(vcf, line)) {
        if (line.empty() || line[0] == '#') continue;
        const auto fields = split_tabs(line);
        if (fields.size() < 10 ||
            (fields[1] != "62623253" && fields[1] != "62645168"))
            continue;
        const size_t separator = fields[9].rfind(':');
        REQUIRE(separator != std::string::npos);
        phase_sets[fields[1]] = fields[9].substr(separator + 1);
    }
    REQUIRE(phase_sets.size() == 2);
    CHECK(phase_sets.at("62623253") != phase_sets.at("62645168"));
}

TEST_CASE("physical bridges survive their owning graph chunks",
          "[gap][stitch-connectivity]") {
    const Paths p = paths();
    if (!p.complete()) {
        WARN("gap-window tests need inputs that are absent: " << p.missing());
        SUCCEED("skipped: inputs absent");
        return;
    }
    struct BridgeCase {
        long long chunk_left;
        long long chunk_right;
        long long gap_left;
        long long gap_right;
        const char* left_site;
        const char* right_site;
        double min_chunk_concordance;
        double min_separated = 0.75;
    };
    const std::array<BridgeCase, 13> bridges{{
        // The two four-base CIGAR deletion placements have one haplotype.
        {1050001, 1950000, 1508171, 1523721,
         "1508171", "1523721", 0.99, 0.70},
        {3050001, 3950000, 3573979, 3596663,
         "3573979", "3596663", 0.98},
        // Physical deletion calls join only the suffix after a weak graph
        // SNP edge, preserving the earlier block as a separate phase set.
        {13050001, 13950000, 13784702, 13800824,
         "13784702", "13800824", 0.93, 0.87},
        {21050001, 21950000, 21736539, 21757517,
         "21736539", "21757517", 0.93},
        // Two repeat-length observations join the boundary insertions.
        // Two clean molecules independently certify the earlier weak edge.
        {21050001, 21950000, 21377985, 21395286,
         "21377985", "21395286", 0.93, 0.84},
        // Two independent high-quality molecules carry the same SNP pair.
        {32000001, 33000000, 32215055, 32233534,
         "32215055", "32233534", 0.85},
        // Two insertion-ALT molecules join the next complete BAM source
        // without reopening the already joined 32.215 Mb graph boundary.
        {32050001, 32950000, 32234664, 32246127,
         "32234664", "32246127", 0.85},
        // The right BAM deletion has a complete source and graph SNP path;
        // its single physical bridge meets the wrong-parity probability bound.
        {23050001, 23950000, 23421003, 23445252,
         "23421003", "23445252", 0.96},
        {54000001, 55000000, 54547514, 54569072,
         "54547514", "54569072", 0.99},
        {54000001, 55000000, 54894127, 54912022,
         "54894127", "54912022", 0.99},
        // A shifted A insertion orients the next BAM block from a lone graph SNP.
        {60050001, 60950000, 60033052, 60048237,
         "60033052", "60048237", 0.98, 0.65},
        // The 50 kb test padding replays the exact 60–61 Mb owning chunk.
        {60050001, 60950000, 60453499, 60467115,
         "60453499", "60467115", 0.98},
        // The right graph block has one SNP before the 19 Mb chunk edge.
        {18050001, 18950000, 18983414, 18999993,
         "18983414", "18999993", 0.98, 0.50},
    }};
    for (const BridgeCase& bridge : bridges) {
        DYNAMIC_SECTION("chr20:" << bridge.gap_left) {
            // A short replay can join a boundary that the full graph block
            // later splits. Keep the complete owning chunk in this test.
            Window chunk;
            chunk.gap_left = bridge.chunk_left;
            chunk.gap_right = bridge.chunk_right;
            std::string dir;
            REQUIRE(run_arm(p, chunk,
                            "stitch_single_molecule_" +
                                std::to_string(bridge.gap_left),
                            "", dir));
            std::ifstream vcf(dir + "/native.vcf");
            REQUIRE(vcf.good());
            std::map<std::string, std::string> phase_sets;
            std::map<std::string, std::string> genotypes;
            std::map<std::string, std::string> earlier_phase_sets;
            std::map<std::string, std::string> earlier_genotypes;
            std::string line;
            while (std::getline(vcf, line)) {
                if (line.empty() || line[0] == '#') continue;
                const auto fields = split_tabs(line);
                if (fields.size() < 10) continue;
                const size_t separator = fields[9].rfind(':');
                if (bridge.gap_left == 13784702 &&
                    (fields[1] == "13752640" ||
                     fields[1] == "13773452")) {
                    REQUIRE(separator != std::string::npos);
                    earlier_phase_sets[fields[1]] =
                        fields[9].substr(separator + 1);
                }
                if (bridge.gap_left == 21377985 &&
                    fields[1] == "21343641") {
                    REQUIRE(separator != std::string::npos);
                    earlier_phase_sets[fields[1]] =
                        fields[9].substr(separator + 1);
                }
                if (((bridge.gap_left == 3573979) &&
                     (fields[1] == "3529324" || fields[1] == "3542977")) ||
                    ((bridge.gap_left == 32234664) &&
                     (fields[1] == "32215055" || fields[1] == "32233534"))) {
                    REQUIRE(separator != std::string::npos);
                    earlier_phase_sets[fields[1]] =
                        fields[9].substr(separator + 1);
                    earlier_genotypes[fields[1]] = fields[9].substr(0, 3);
                }
                if (fields[1] != bridge.left_site &&
                    fields[1] != bridge.right_site)
                    continue;
                if (bridge.gap_left == 13784702) {
                    CHECK(fields[3] == (fields[1] == bridge.left_site ?
                        "AT" : "GATA"));
                    CHECK(fields[4] == (fields[1] == bridge.left_site ?
                        "A" : "G"));
                }
                if (bridge.gap_left == 21377985) {
                    CHECK(fields[3] == (fields[1] == bridge.left_site ?
                        "C" : "A"));
                    CHECK(fields[4] == (fields[1] == bridge.left_site ?
                        "CTGTGTG" : "AT"));
                }
                if (bridge.gap_left == 23421003) {
                    CHECK(fields[3] == (fields[1] == bridge.left_site ?
                        "C" : "TGAAAGAAGA"));
                    CHECK(fields[4] == "T");
                }
                if (bridge.gap_left == 32234664) {
                    CHECK(fields[3] == (fields[1] == bridge.left_site ? "A" : "T"));
                    CHECK(fields[4] == (fields[1] == bridge.left_site ? "T" :
                        "TGGAATGGAATGGAATGGAATGGAATGTAGTCAACGCGAGT"
                        "GGAATGGATTGGAATGGAATGGAAG"));
                }
                REQUIRE(separator != std::string::npos);
                phase_sets[fields[1]] = fields[9].substr(separator + 1);
                genotypes[fields[1]] = fields[9].substr(0, 3);
            }
            REQUIRE(phase_sets.size() == 2);
            CHECK(phase_sets.at(bridge.left_site) != ".");
            CHECK(phase_sets.at(bridge.left_site) ==
                  phase_sets.at(bridge.right_site));
            if (bridge.gap_left == 13784702) {
                REQUIRE(earlier_phase_sets.size() == 2);
                CHECK(earlier_phase_sets.at("13752640") !=
                      phase_sets.at(bridge.left_site));
                CHECK(earlier_phase_sets.at("13773452") ==
                      phase_sets.at(bridge.left_site));
                CHECK(genotypes.at(bridge.left_site) ==
                      genotypes.at(bridge.right_site));
            }
            if (bridge.gap_left == 21377985) {
                REQUIRE(earlier_phase_sets.size() == 1);
                CHECK(earlier_phase_sets.at("21343641") ==
                      phase_sets.at(bridge.left_site));
                CHECK(genotypes.at(bridge.left_site) ==
                      genotypes.at(bridge.right_site));
            }
            if (bridge.gap_left == 32234664) {
                REQUIRE(earlier_phase_sets.size() == 2);
                CHECK(earlier_phase_sets.at("32215055") ==
                      phase_sets.at(bridge.left_site));
                CHECK(earlier_phase_sets.at("32233534") ==
                      phase_sets.at(bridge.right_site));
                CHECK(earlier_genotypes.at("32215055") !=
                      earlier_genotypes.at("32233534"));
            }
            if (bridge.gap_left == 3573979) {
                REQUIRE(earlier_phase_sets.size() == 2);
                // The BAM-certified path keeps the intermediate deletion
                // bridge intact when the earlier SNP seam joins this block.
                CHECK(earlier_phase_sets.at("3529324") ==
                      earlier_phase_sets.at("3542977"));
                CHECK(earlier_genotypes.at("3529324") !=
                      earlier_genotypes.at("3542977"));
                Window earlier_gap;
                earlier_gap.gap_left = 3529324;
                earlier_gap.gap_right = 3542977;
                Outcome earlier_reads;
                score_bam(dir + "/phased.bam", earlier_gap,
                          load_truth(p.truth_map),
                          input_read_spans(p, earlier_gap), earlier_reads);
                CHECK_FALSE(earlier_reads.switched);
                CHECK(earlier_reads.concordance() >= 0.98);
            }
            if (bridge.gap_left == 54894127 ||
                bridge.gap_left == 18983414 ||
                bridge.gap_left == 1508171)
                CHECK(genotypes.at(bridge.left_site) ==
                      genotypes.at(bridge.right_site));
            if (bridge.gap_left == 54547514 ||
                bridge.gap_left == 32215055 ||
                bridge.gap_left == 32234664 ||
                bridge.gap_left == 23421003 ||
                bridge.gap_left == 60453499 ||
                bridge.gap_left == 60033052)
                CHECK(genotypes.at(bridge.left_site) !=
                      genotypes.at(bridge.right_site));

            Window gap;
            gap.gap_left = bridge.gap_left;
            gap.gap_right = bridge.gap_right;
            Outcome reads;
            score_bam(dir + "/phased.bam", gap, load_truth(p.truth_map),
                      input_read_spans(p, gap), reads);
            CHECK_FALSE(reads.switched);
            CHECK(reads.separated() >= bridge.min_separated);
            // Concordance covers the whole million-base chunk. The 21 Mb
            // chunk has unrelated discordant reads outside this seam.
            CHECK(reads.concordance() >= bridge.min_chunk_concordance);
        }
    }
}

TEST_CASE("phased right reads bridge only the 64.140 Mb deletion row",
          "[gap][representation]") {
    const Paths p = paths();
    if (!p.complete()) {
        WARN("gap-window tests need inputs that are absent: " << p.missing());
        SUCCEED("skipped: inputs absent");
        return;
    }
    Window chunk;
    chunk.gap_left = 64050001;
    chunk.gap_right = 64950000;
    std::string dir;
    REQUIRE(run_arm(p, chunk, "stitch_right_read_gauge_64m", "", dir));

    std::ifstream vcf(dir + "/native.vcf");
    REQUIRE(vcf.good());
    std::map<std::string, std::pair<std::string, std::string>> rows;
    std::string line;
    while (std::getline(vcf, line)) {
        if (line.empty() || line[0] == '#') continue;
        const auto fields = split_tabs(line);
        if (fields.size() < 10) continue;
        const std::string key =
            fields[1] + ":" + fields[3] + ">" + fields[4];
        if (key != "64118182:T>C" &&
            key != "64128828:G>C" &&
            key != "64134226:CT>C" &&
            key != "64134226:CTT>C" &&
            key != "64138752:CA>C" &&
            key != "64140314:AT>A" &&
            key != "64144256:CA>C" &&
            key != "64144256:CAA>C")
            continue;
        const size_t separator = fields[9].rfind(':');
        REQUIRE(separator != std::string::npos);
        rows[key] = {fields[9].substr(0, 3),
                     fields[9].substr(separator + 1)};
    }
    REQUIRE(rows.size() == 8);
    const auto& left = rows.at("64140314:AT>A");
    const auto& right = rows.at("64144256:CA>C");
    const auto& complement = rows.at("64144256:CAA>C");
    CHECK(left.second == right.second);
    CHECK(right.second == complement.second);
    CHECK(rows.at("64118182:T>C").second != left.second);
    CHECK(rows.at("64128828:G>C").second != left.second);
    // The complementary BAM deletion pair certifies the intervening run;
    // the graph SNP before it still lacks a paired allele edge.
    CHECK(rows.at("64134226:CT>C").second == left.second);
    CHECK(rows.at("64134226:CTT>C").second == left.second);
    CHECK(rows.at("64138752:CA>C").second == left.second);
    CHECK(left.first != right.first);
    CHECK(right.first != complement.first);

    Window gap;
    gap.gap_left = 64140314;
    gap.gap_right = 64144256;
    Outcome reads;
    score_bam(dir + "/phased.bam", gap, load_truth(p.truth_map),
              input_read_spans(p, gap), reads);
    CHECK_FALSE(reads.switched);
    CHECK(reads.concordance() >= 0.99);
    CHECK(reads.separated() >= 0.72);

    gap.gap_left = 64138752;
    gap.gap_right = 64140314;
    Outcome joined_reads;
    score_bam(dir + "/phased.bam", gap, load_truth(p.truth_map),
              input_read_spans(p, gap), joined_reads);
    CHECK_FALSE(joined_reads.switched);
    CHECK(joined_reads.concordance() >= 0.99);
    CHECK(joined_reads.separated() >= 0.71);
}

TEST_CASE("source allele repair preserves complementary repeat rows",
          "[gap][representation]") {
    const Paths p = paths();
    if (!p.complete()) {
        WARN("gap-window tests need inputs that are absent: " << p.missing());
        SUCCEED("skipped: inputs absent");
        return;
    }
    // The source repair can connect this locus using its restored allele
    // observations. Keep both complementary rows; insertion length alone
    // still cannot reinterpret the deletion as insertion REF.
    Window chunk;
    chunk.gap_left = 23421003;
    chunk.gap_right = 23445252;
    std::string dir;
    REQUIRE(run_arm(p, chunk, "repeat_multiallelic_guard", "", dir));
    std::ifstream vcf(dir + "/native.vcf");
    REQUIRE(vcf.good());
    std::map<std::string, std::pair<std::string, std::string>> rows;
    std::string line;
    while (std::getline(vcf, line)) {
        if (line.empty() || line[0] == '#') continue;
        const auto fields = split_tabs(line);
        if (fields.size() < 10) continue;
        if (fields[1] != "23792418" &&
            fields[1] != "23806565")
            continue;
        const std::string key =
            fields[1] + ":" + fields[3] + ">" + fields[4];
        if (key != "23792418:T>TT" &&
            key != "23806565:T>TAC" &&
            key != "23806565:TACACACAC>T")
            continue;
        const size_t separator = fields[9].rfind(':');
        REQUIRE(separator != std::string::npos);
        rows[key] = {fields[9].substr(0, 3),
                     fields[9].substr(separator + 1)};
    }
    REQUIRE(rows.size() == 3);
    CHECK(rows.at("23792418:T>TT").second ==
          rows.at("23806565:T>TAC").second);
    CHECK(rows.at("23792418:T>TT").first ==
          rows.at("23806565:T>TAC").first);
    CHECK(rows.at("23806565:T>TAC").second ==
          rows.at("23806565:TACACACAC>T").second);
    CHECK(rows.at("23806565:T>TAC").first !=
          rows.at("23806565:TACACACAC>T").first);
}

TEST_CASE("graph haplotypes bridge equivalent complementary BAM deletions",
          "[gap][stitch-connectivity]") {
    const Paths p = paths();
    if (!p.complete()) {
        WARN("gap-window tests need inputs that are absent: " << p.missing());
        SUCCEED("skipped: inputs absent");
        return;
    }
    // Use the whole owning chunk: the right deletion pair belongs to a BAM
    // phase block, while the left SNP's read haplotype comes from the graph.
    Window chunk;
    chunk.gap_left = 59050001;
    chunk.gap_right = 59950000;
    std::string dir;
    REQUIRE(run_arm(p, chunk, "stitch_graph_snp_deletions", "", dir));
    std::ifstream vcf(dir + "/native.vcf");
    REQUIRE(vcf.good());
    std::string left_ps, short_ps, long_ps;
    std::string left_gt, short_gt, long_gt;
    std::string line;
    while (std::getline(vcf, line)) {
        if (line.empty() || line[0] == '#') continue;
        const auto fields = split_tabs(line);
        if (fields.size() < 10) continue;
        const size_t separator = fields[9].rfind(':');
        REQUIRE(separator != std::string::npos);
        const std::string gt = fields[9].substr(0, 3);
        const std::string ps = fields[9].substr(separator + 1);
        if (fields[1] == "59825454" && fields[3] == "C" &&
            fields[4] == "T") {
            left_gt = gt;
            left_ps = ps;
        } else if (fields[1] == "59842960" && fields[3] == "GTGTT" &&
                   fields[4] == "G") {
            short_gt = gt;
            short_ps = ps;
        } else if (fields[1] == "59842960" &&
                   fields[3] == "GTGTTTGTT" && fields[4] == "G") {
            long_gt = gt;
            long_ps = ps;
        }
    }
    REQUIRE_FALSE(left_ps.empty());
    REQUIRE_FALSE(short_ps.empty());
    REQUIRE_FALSE(long_ps.empty());
    CHECK(left_ps != ".");
    CHECK(left_ps == short_ps);
    CHECK(left_ps == long_ps);
    CHECK(left_gt == short_gt);
    CHECK(left_gt != long_gt);

    Window gap;
    gap.gap_left = 59825454;
    gap.gap_right = 59842960;
    Outcome reads;
    score_bam(dir + "/phased.bam", gap, load_truth(p.truth_map),
              input_read_spans(p, gap), reads);
    CHECK_FALSE(reads.switched);
    CHECK(reads.concordance() >= 0.99);
}

TEST_CASE("physical links close the long 47 Mb gap without a switch",
          "[gap][stitch-connectivity]") {
    const Paths p = paths();
    if (!p.complete()) {
        WARN("gap-window tests need inputs that are absent: " << p.missing());
        SUCCEED("skipped: inputs absent");
        return;
    }
    Window gap;
    gap.gap_left = 47003897;
    gap.gap_right = 47713869;
    std::string dir;
    REQUIRE(run_arm(p, gap, "stitch_long_47m", "", dir));
    std::ifstream vcf(dir + "/native.vcf");
    REQUIRE(vcf.good());
    std::map<std::string, std::pair<std::string, std::string>> sites;
    std::string line;
    while (std::getline(vcf, line)) {
        if (line.empty() || line[0] == '#') continue;
        const auto fields = split_tabs(line);
        if (fields.size() < 10 ||
            (fields[1] != "47003897" && fields[1] != "47689418" &&
             fields[1] != "47694119" && fields[1] != "47713869"))
            continue;
        const size_t separator = fields[9].rfind(':');
        REQUIRE(separator != std::string::npos);
        sites[fields[1]] = {fields[9].substr(0, 3),
                            fields[9].substr(separator + 1)};
    }
    REQUIRE(sites.size() == 4);
    const std::string& phase_set = sites.at("47003897").second;
    REQUIRE(phase_set != ".");
    CHECK(phase_set == sites.at("47689418").second);
    CHECK(phase_set == sites.at("47694119").second);
    CHECK(phase_set == sites.at("47713869").second);
    CHECK(sites.at("47003897").first != sites.at("47713869").first);

    Outcome reads;
    score_bam(dir + "/phased.bam", gap, load_truth(p.truth_map),
              input_read_spans(p, gap), reads);
    CHECK_FALSE(reads.switched);
    CHECK(reads.concordance() >= 0.98);
    CHECK(reads.separated() >= 0.90);
}

TEST_CASE("validated graph bridge survives a BAM phase-set ID change",
          "[gap][stitch-connectivity]") {
    const Paths p = paths();
    if (!p.complete()) {
        WARN("gap-window tests need inputs that are absent: " << p.missing());
        SUCCEED("skipped: inputs absent");
        return;
    }
    // The complete left graph block starts before the local BAM source block.
    // Its validation solve and the targeted injection solve therefore assign
    // different numeric PS IDs to the physical SNP/deletion bridge.
    Window region;
    region.gap_left = 22950001;
    region.gap_right = 23050000;
    std::string dir;
    REQUIRE(run_arm(p, region, "stitch_validation_ps_remap", "", dir));
    std::ifstream vcf(dir + "/native.vcf");
    REQUIRE(vcf.good());
    std::map<std::string, std::string> phase_sets;
    std::string line;
    while (std::getline(vcf, line)) {
        if (line.empty() || line[0] == '#') continue;
        const auto fields = split_tabs(line);
        if (fields.size() < 10 ||
            (fields[1] != "22980600" && fields[1] != "23008891"))
            continue;
        const size_t separator = fields[9].rfind(':');
        REQUIRE(separator != std::string::npos);
        phase_sets[fields[1]] = fields[9].substr(separator + 1);
    }
    REQUIRE(phase_sets.size() == 2);
    CHECK(phase_sets.at("22980600") != ".");
    CHECK(phase_sets.at("22980600") == phase_sets.at("23008891"));

    Window gap;
    gap.gap_left = 22980600;
    gap.gap_right = 23008891;
    Outcome reads;
    score_bam(dir + "/phased.bam", gap, load_truth(p.truth_map),
              input_read_spans(p, gap), reads);
    CHECK_FALSE(reads.switched);
    CHECK(reads.concordance() >= 0.99);
}

TEST_CASE("one shared BAM anchor cannot absorb the earlier 19.4 Mb block",
          "[gap][stitch-connectivity]") {
    const Paths p = paths();
    if (!p.complete()) {
        WARN("gap-window tests need inputs that are absent: " << p.missing());
        SUCCEED("skipped: inputs absent");
        return;
    }
    // The two direct clean-SNP bridges join the local BAM run to the right
    // block. One shared BAM anchor still cannot absorb the earlier graph
    // prefix across the source's weak cut.
    Window region;
    region.gap_left = 19050001;
    region.gap_right = 19950000;
    std::string dir;
    REQUIRE(run_arm(p, region, "stitch_one_anchor_19m", "", dir));
    std::ifstream vcf(dir + "/native.vcf");
    REQUIRE(vcf.good());
    std::map<std::string, std::pair<std::string, std::string>> sites;
    std::string line;
    while (std::getline(vcf, line)) {
        if (line.empty() || line[0] == '#') continue;
        const auto fields = split_tabs(line);
        if (fields.size() < 10) continue;
        const std::string key = fields[1] + ":" + fields[3] +
            ">" + fields[4];
        if (key != "19377345:T>TATATATAGAG" &&
            key != "19403172:T>C" && key != "19414720:G>A")
            continue;
        const size_t separator = fields[9].rfind(':');
        REQUIRE(separator != std::string::npos);
        sites[key] = {fields[9].substr(0, 3),
                      fields[9].substr(separator + 1)};
    }
    REQUIRE(sites.size() == 3);
    const auto& prefix = sites.at("19377345:T>TATATATAGAG");
    const auto& middle = sites.at("19403172:T>C");
    const auto& last = sites.at("19414720:G>A");
    CHECK(prefix.second != middle.second);
    CHECK(middle.second == last.second);
    CHECK(middle.first == last.first);
}

TEST_CASE("direct SNP proof reuses a graph block at 56.15 Mb",
          "[gap][stitch-connectivity]") {
    const Paths p = paths();
    if (!p.complete()) {
        WARN("gap-window tests need inputs that are absent: " << p.missing());
        SUCCEED("skipped: inputs absent");
        return;
    }
    // The left graph block is already joined upstream. Its next BAM edge has
    // an exact MEC SNP certificate in both read halves, so reuse is valid.
    Window region;
    region.gap_left = 56050001;
    region.gap_right = 56950000;
    std::string dir;
    REQUIRE(run_arm(p, region, "stitch_reused_direct_snp_56m", "", dir));
    std::ifstream vcf(dir + "/native.vcf");
    REQUIRE(vcf.good());
    std::map<std::string, std::string> phase_sets;
    std::string line;
    while (std::getline(vcf, line)) {
        if (line.empty() || line[0] == '#') continue;
        const auto fields = split_tabs(line);
        if (fields.size() < 10 ||
            (fields[1] != "56150368" && fields[1] != "56156525"))
            continue;
        const size_t separator = fields[9].rfind(':');
        REQUIRE(separator != std::string::npos);
        phase_sets[fields[1]] = fields[9].substr(separator + 1);
    }
    REQUIRE(phase_sets.size() == 2);
    CHECK(phase_sets.at("56150368") != ".");
    CHECK(phase_sets.at("56150368") == phase_sets.at("56156525"));

    Window gap;
    gap.gap_left = 56150368;
    gap.gap_right = 56156525;
    Outcome reads;
    score_bam(dir + "/phased.bam", gap, load_truth(p.truth_map),
              input_read_spans(p, gap), reads);
    CHECK_FALSE(reads.switched);
    CHECK(reads.concordance() >= 0.95);
}

TEST_CASE("physical SNP bridge crosses the 23 Mb graph chunk boundary",
          "[gap][stitch-connectivity]") {
    const Paths p = paths();
    if (!p.complete()) {
        WARN("gap-window tests need inputs that are absent: " << p.missing());
        SUCCEED("skipped: inputs absent");
        return;
    }
    // Run the two owning 1-Mb chunks: a single-window replay already joined
    // these alleles, but the full chromosome used to split them at 23 Mb.
    Window region;
    region.gap_left = 22050001;
    region.gap_right = 23950000;
    std::string dir;
    REQUIRE(run_arm(p, region, "stitch_cross_chunk_23m",
                    "--chunk-size 1000000", dir));
    std::ifstream vcf(dir + "/native.vcf");
    REQUIRE(vcf.good());
    std::map<std::string, std::string> phase_sets;
    std::string line;
    while (std::getline(vcf, line)) {
        if (line.empty() || line[0] == '#') continue;
        const auto fields = split_tabs(line);
        if (fields.size() < 10 ||
            (fields[1] != "22980600" && fields[1] != "23008891"))
            continue;
        const size_t separator = fields[9].rfind(':');
        REQUIRE(separator != std::string::npos);
        phase_sets[fields[1]] = fields[9].substr(separator + 1);
    }
    REQUIRE(phase_sets.size() == 2);
    CHECK(phase_sets.at("22980600") != ".");
    CHECK(phase_sets.at("22980600") == phase_sets.at("23008891"));

    Window gap;
    gap.gap_left = 22980600;
    gap.gap_right = 23008891;
    Outcome reads;
    score_bam(dir + "/phased.bam", gap, load_truth(p.truth_map),
              input_read_spans(p, gap), reads);
    CHECK_FALSE(reads.switched);
    CHECK(reads.concordance() >= 0.97);
}

TEST_CASE("multiple BAM source labels preserve a 22.98 Mb graph bridge",
          "[gap][stitch-connectivity]") {
    const Paths p = paths();
    if (!p.complete()) {
        WARN("gap-window tests need inputs that are absent: " << p.missing());
        SUCCEED("skipped: inputs absent");
        return;
    }
    // The complete left graph block has exact matches to two BAM phase sets.
    // Its boundary source remains coherent with the graph and the physical
    // SNP/deletion bridge, even though the earlier BAM source has another ID.
    Window region;
    region.gap_left = 22550001;
    region.gap_right = 23450000;
    std::string dir;
    REQUIRE(run_arm(p, region, "stitch_multiple_source_labels_23m", "", dir));
    std::ifstream vcf(dir + "/native.vcf");
    REQUIRE(vcf.good());
    std::map<std::string, std::string> phase_sets;
    std::string line;
    while (std::getline(vcf, line)) {
        if (line.empty() || line[0] == '#') continue;
        const auto fields = split_tabs(line);
        if (fields.size() < 10 ||
            (fields[1] != "22980600" && fields[1] != "23008891"))
            continue;
        const size_t separator = fields[9].rfind(':');
        REQUIRE(separator != std::string::npos);
        phase_sets[fields[1]] = fields[9].substr(separator + 1);
    }
    REQUIRE(phase_sets.size() == 2);
    CHECK(phase_sets.at("22980600") == phase_sets.at("23008891"));

    Window gap;
    gap.gap_left = 22980600;
    gap.gap_right = 23008891;
    Outcome reads;
    score_bam(dir + "/phased.bam", gap, load_truth(p.truth_map),
              input_read_spans(p, gap), reads);
    CHECK_FALSE(reads.switched);
    CHECK(reads.concordance() >= 0.98);
}

TEST_CASE("a different deletion locus cannot validate a remapped bridge",
          "[gap][stitch-connectivity]") {
    const Paths p = paths();
    if (!p.complete()) {
        WARN("gap-window tests need inputs that are absent: " << p.missing());
        SUCCEED("skipped: inputs absent");
        return;
    }
    // The 37.56-Mb boundary has one physical deletion bridge, but the other
    // deletion row is 173 bp away. Accepting that single allele after its BAM
    // source PS changes joins two graph blocks in the wrong orientation.
    Window region;
    region.gap_left = 37050001;
    region.gap_right = 37950000;
    std::string dir;
    REQUIRE(run_arm(p, region, "stitch_distinct_deletion_loci", "", dir));
    std::ifstream vcf(dir + "/native.vcf");
    REQUIRE(vcf.good());
    std::map<std::string, std::string> phase_sets;
    std::string line;
    while (std::getline(vcf, line)) {
        if (line.empty() || line[0] == '#') continue;
        const auto fields = split_tabs(line);
        if (fields.size() < 10 ||
            (fields[1] != "37556826" && fields[1] != "37606443"))
            continue;
        const size_t separator = fields[9].rfind(':');
        REQUIRE(separator != std::string::npos);
        phase_sets[fields[1]] = fields[9].substr(separator + 1);
    }
    REQUIRE(phase_sets.size() == 2);
    CHECK(phase_sets.at("37556826") != ".");
    CHECK(phase_sets.at("37606443") != ".");
    CHECK(phase_sets.at("37556826") != phase_sets.at("37606443"));
}

TEST_CASE("BAM transfer preserves an already connected graph phase set",
          "[gap][stitch-connectivity]") {
    const Paths p = paths();
    if (!p.complete()) {
        WARN("gap-window tests need inputs that are absent: " << p.missing());
        SUCCEED("skipped: inputs absent");
        return;
    }
    Window w;
    w.gap_left = 11050001;
    w.gap_right = 11950000;
    std::string dir;
    REQUIRE(run_arm(p, w, "stitch_preserve_graph", "", dir));

    std::ifstream in(dir + "/native.vcf");
    REQUIRE(in.good());
    std::map<std::string, std::string> phase_sets;
    std::string line;
    while (std::getline(in, line)) {
        if (line.empty() || line[0] == '#') continue;
        const auto fields = split_tabs(line);
        if (fields.size() < 10 ||
            (fields[1] != "11357244" && fields[1] != "11360353"))
            continue;
        const size_t separator = fields[9].rfind(':');
        REQUIRE(separator != std::string::npos);
        phase_sets[fields[1]] = fields[9].substr(separator + 1);
    }
    REQUIRE(phase_sets.size() == 2);
    CHECK(phase_sets.at("11357244") != ".");
    CHECK(phase_sets.at("11357244") == phase_sets.at("11360353"));
}

TEST_CASE("complementary BAM boundary rows close the 41.900 Mb seam",
          "[gap][stitch-connectivity]") {
    const Paths p = paths();
    if (!p.complete()) {
        WARN("gap-window tests need inputs that are absent: " << p.missing());
        SUCCEED("skipped: inputs absent");
        return;
    }
    // The preceding seam shares a BAM solve region with this deletion pair.
    // Its block must stay independent while the focused MSA closes this one.
    Window full_chunk;
    full_chunk.gap_left = 41050001;
    full_chunk.gap_right = 41950000;
    std::string dir;
    REQUIRE(run_arm(p, full_chunk, "stitch_complementary_41m", "", dir));

    std::ifstream in(dir + "/native.vcf");
    REQUIRE(in.good());
    std::map<std::string, std::string> phase_sets;
    std::string line;
    while (std::getline(in, line)) {
        if (line.empty() || line[0] == '#') continue;
        const auto fields = split_tabs(line);
        if (fields.size() < 10 ||
            (fields[1] != "41866917" && fields[1] != "41898323" &&
             fields[1] != "41900800" && fields[1] != "41919471"))
            continue;
        const size_t separator = fields[9].rfind(':');
        REQUIRE(separator != std::string::npos);
        phase_sets[fields[1]] = fields[9].substr(separator + 1);
    }
    REQUIRE(phase_sets.size() == 4);
    CHECK(phase_sets.at("41866917") != phase_sets.at("41898323"));
    CHECK(phase_sets.at("41900800") == phase_sets.at("41919471"));

    Window gap;
    gap.gap_left = 41900800;
    gap.gap_right = 41919471;
    Outcome reads;
    score_bam(dir + "/phased.bam", gap, load_truth(p.truth_map),
              input_read_spans(p, gap), reads);
    CHECK_FALSE(reads.switched);
    CHECK(reads.concordance() >= 0.98);
}

TEST_CASE("certified deletion bridge preserves the 11.599 Mb source alleles",
          "[gap][stitch-connectivity][msa][source-path]") {
    const Paths p = paths();
    if (!p.complete()) {
        WARN("gap-window tests need inputs that are absent: " << p.missing());
        SUCCEED("skipped: inputs absent");
        return;
    }
    // The incomplete focused retry must still leave the original BAM rows
    // intact. A separately certified physical deletion bridge can now join
    // those blocks without replacing their solve or reversing their alleles.
    Window full_chunk;
    full_chunk.gap_left = 11050001;
    full_chunk.gap_right = 11950000;
    std::string dir;
    REQUIRE(run_arm(p, full_chunk, "stitch_incomplete_11m", "", dir));

    const std::map<std::string, std::string> expected = {
        {"11517849:C>T", "0|1"},
        {"11573074:G>C", "1|0"},
        {"11586531:T>TACACACACACACAC", "0|1"},
        {"11586531:T>TACACACACACACACACAC", "1|0"},
        {"11591585:CT>C", "1|0"},
        {"11599138:G>A", "1|0"}};
    std::ifstream in(dir + "/native.vcf");
    REQUIRE(in.good());
    std::map<std::string, std::pair<std::string, std::string>> sites;
    std::string line;
    while (std::getline(in, line)) {
        if (line.empty() || line[0] == '#') continue;
        const auto fields = split_tabs(line);
        if (fields.size() < 10) continue;
        const std::string key = fields[1] + ":" + fields[3] + ">" + fields[4];
        if (expected.count(key) == 0) continue;
        const size_t separator = fields[9].rfind(':');
        REQUIRE(separator != std::string::npos);
        sites.emplace(key, std::make_pair(fields[9].substr(0, 3),
                                          fields[9].substr(separator + 1)));
    }
    REQUIRE(sites.size() == expected.size());
    const std::string phase_set = sites.at("11573074:G>C").second;
    CHECK(phase_set != ".");
    for (const auto& [key, genotype] : expected) {
        CHECK(sites.at(key).first == genotype);
        CHECK(sites.at(key).second == phase_set);
    }
    CHECK(count_hap_allele_conflicts(dir + "/native.vcf") == 0);

    const auto truth = load_truth(p.truth_map);
    for (const auto& interval : {std::pair<long long, long long>{11586531, 11599138},
                                {11573074, 11586531}}) {
        Window gap;
        gap.gap_left = interval.first;
        gap.gap_right = interval.second;
        Outcome reads;
        score_bam(dir + "/phased.bam", gap, truth, input_read_spans(p, gap), reads);
        CHECK_FALSE(reads.switched);
        CHECK(reads.concordance() >= 0.97);
        if (gap.gap_left == 11573074) {
            CHECK(reads.separated() >= 0.64);
            // score_bam counts the whole owning chunk; its accepted baseline
            // has 90 discordant reads, including sites outside this gap.
            CHECK(reads.discordant() <= 90);
        }
    }
}

TEST_CASE("complete BAM source path closes the 48.929 Mb seam in its owning chunk",
          "[gap][stitch-connectivity]") {
    const Paths p = paths();
    if (!p.complete()) {
        WARN("gap-window tests need inputs that are absent: " << p.missing());
        SUCCEED("skipped: inputs absent");
        return;
    }
    // The focused BAM solve spans both graph flanks, but its only private
    // in-gap row is an indel. Keep the full adjacent phase-set context so
    // the source path and both graph/BAM gauges can validate the join.
    Window full_chunk;
    full_chunk.gap_left = 48050001;
    full_chunk.gap_right = 48950000;
    std::string dir;
    REQUIRE(run_arm(p, full_chunk, "stitch_full_ps_48m", "", dir));

    std::ifstream in(dir + "/native.vcf");
    REQUIRE(in.good());
    std::map<std::string, std::string> phase_sets;
    std::string line;
    while (std::getline(in, line)) {
        if (line.empty() || line[0] == '#') continue;
        const auto fields = split_tabs(line);
        if (fields.size() < 10 ||
            (fields[1] != "48929511" && fields[1] != "48950388"))
            continue;
        const size_t separator = fields[9].rfind(':');
        REQUIRE(separator != std::string::npos);
        phase_sets[fields[1]] = fields[9].substr(separator + 1);
    }
    REQUIRE(phase_sets.size() == 2);
    CHECK(phase_sets.at("48929511") != ".");
    CHECK(phase_sets.at("48929511") == phase_sets.at("48950388"));

    Window gap;
    gap.gap_left = 48929511;
    gap.gap_right = 48950388;
    Outcome reads;
    score_bam(dir + "/phased.bam", gap, load_truth(p.truth_map),
              input_read_spans(p, gap), reads);
    CHECK_FALSE(reads.switched);
    CHECK(reads.concordance() >= 0.98);
}

TEST_CASE("complete adjacent graph phase sets close the 56 Mb seam in its owning chunk",
          "[gap][stitch-connectivity]") {
    const Paths p = paths();
    if (!p.complete()) {
        WARN("gap-window tests need inputs that are absent: " << p.missing());
        SUCCEED("skipped: inputs absent");
        return;
    }
    // The 120 kb replay has one BAM source block, but the ordinary grouped
    // full-chunk solve contains an unrelated later seam. The focused retry
    // must keep the established BAM rows and connect both complete graph PSs.
    Window full_chunk;
    full_chunk.gap_left = 56050001;
    full_chunk.gap_right = 56950000;
    std::string dir;
    REQUIRE(run_arm(p, full_chunk, "stitch_full_ps_56m", "", dir));

    std::ifstream in(dir + "/native.vcf");
    REQUIRE(in.good());
    std::map<std::string, std::string> phase_sets;
    std::string line;
    while (std::getline(in, line)) {
        if (line.empty() || line[0] == '#') continue;
        const auto fields = split_tabs(line);
        if (fields.size() < 10 ||
            (fields[1] != "56064697" && fields[1] != "56083708" &&
             fields[1] != "56107742"))
            continue;
        const size_t separator = fields[9].rfind(':');
        REQUIRE(separator != std::string::npos);
        phase_sets[fields[1]] = fields[9].substr(separator + 1);
    }
    REQUIRE(phase_sets.size() == 3);
    CHECK(phase_sets.at("56064697") != ".");
    CHECK(phase_sets.at("56064697") == phase_sets.at("56083708"));
    CHECK(phase_sets.at("56083708") == phase_sets.at("56107742"));

    Window gap;
    gap.gap_left = 56064697;
    gap.gap_right = 56083708;
    Outcome reads;
    score_bam(dir + "/phased.bam", gap, load_truth(p.truth_map),
              input_read_spans(p, gap), reads);
    CHECK_FALSE(reads.switched);
    CHECK(reads.concordance() >= 0.98);
}

TEST_CASE("complete BAM path and verified deletion close the 5.31 Mb seam",
          "[gap][stitch-connectivity]") {
    const Paths p = paths();
    if (!p.complete()) {
        WARN("gap-window tests need inputs that are absent: " << p.missing());
        SUCCEED("skipped: inputs absent");
        return;
    }
    // Use the entire owning chunk: an earlier seam has already reused the
    // left graph block, and the right BAM source has a weak cut farther away.
    Window chunk;
    chunk.gap_left = 5050001;
    chunk.gap_right = 5950000;
    std::string dir;
    REQUIRE(run_arm(p, chunk, "stitch_complete_graph_path_5m", "", dir));

    std::ifstream vcf(dir + "/native.vcf");
    REQUIRE(vcf.good());
    std::map<std::string, std::string> phase_sets;
    std::string line;
    while (std::getline(vcf, line)) {
        if (line.empty() || line[0] == '#') continue;
        const auto fields = split_tabs(line);
        if (fields.size() < 10 ||
            (fields[1] != "5309406" && fields[1] != "5345085" &&
             fields[1] != "5350509" && fields[1] != "5393615"))
            continue;
        const size_t separator = fields[9].rfind(':');
        REQUIRE(separator != std::string::npos);
        phase_sets[fields[1]] = fields[9].substr(separator + 1);
    }
    REQUIRE(phase_sets.size() == 4);
    CHECK(phase_sets.at("5309406") != ".");
    CHECK(phase_sets.at("5309406") == phase_sets.at("5345085"));
    // Direct physical SNP calls join the local edge. The known source cut
    // farther right still leaves the distant graph block independent.
    CHECK(phase_sets.at("5345085") == phase_sets.at("5350509"));
    CHECK(phase_sets.at("5350509") != phase_sets.at("5393615"));

    Window gap;
    gap.gap_left = 5309406;
    gap.gap_right = 5345085;
    Outcome reads;
    score_bam(dir + "/phased.bam", gap, load_truth(p.truth_map),
              input_read_spans(p, gap), reads);
    CHECK_FALSE(reads.switched);
    CHECK(reads.concordance() >= 0.95);

    Window local_edge;
    local_edge.gap_left = 5345085;
    local_edge.gap_right = 5350509;
    Outcome local_reads;
    score_bam(dir + "/phased.bam", local_edge, load_truth(p.truth_map),
              input_read_spans(p, local_edge), local_reads);
    CHECK_FALSE(local_reads.switched);
    CHECK(local_reads.concordance() >= 0.99);
}

TEST_CASE("BAM-left MEC cannot reverse the 34.1 Mb graph block",
          "[gap][stitch-connectivity]") {
    const Paths p = paths();
    if (!p.complete()) {
        WARN("gap-window tests need inputs that are absent: " << p.missing());
        SUCCEED("skipped: inputs absent");
        return;
    }
    Window chunk;
    chunk.gap_left = 34050001;
    chunk.gap_right = 34950000;
    std::string dir;
    REQUIRE(run_arm(p, chunk, "stitch_bam_left_34m", "", dir));

    std::ifstream vcf(dir + "/native.vcf");
    REQUIRE(vcf.good());
    std::map<std::string, std::string> phase_sets;
    bool phased_off_edge_insertion = false;
    std::string line;
    while (std::getline(vcf, line)) {
        if (line.empty() || line[0] == '#') continue;
        const auto fields = split_tabs(line);
        if (fields.size() >= 2 && fields[1] == "34146891")
            phased_off_edge_insertion = true;
        if (fields.size() < 10 ||
            (fields[1] != "34102867" && fields[1] != "34122094" &&
             fields[1] != "34134418"))
            continue;
        const size_t separator = fields[9].rfind(':');
        REQUIRE(separator != std::string::npos);
        phase_sets[fields[1]] = fields[9].substr(separator + 1);
    }
    REQUIRE(phase_sets.size() == 3);
    CHECK(phase_sets.at("34102867") != phase_sets.at("34134418"));
    CHECK(phase_sets.at("34122094") == phase_sets.at("34134418"));
    CHECK_FALSE(phased_off_edge_insertion);
}

TEST_CASE("BAM recovery admits missing MSA pairs at 17.62 Mb",
          "[gap][stitch-connectivity]") {
    const Paths p = paths();
    if (!p.complete()) {
        WARN("gap-window tests need inputs that are absent: " << p.missing());
        SUCCEED("skipped: inputs absent");
        return;
    }
    Window w;
    w.gap_left = 17616778;
    w.gap_right = 17625527;
    std::string dir;
    REQUIRE(run_arm(p, w, "stitch_msa_pair", "", dir));

    std::ifstream in(dir + "/native.vcf");
    REQUIRE(in.good());
    std::map<std::string, std::string> phase_sets;
    std::string line;
    while (std::getline(in, line)) {
        if (line.empty() || line[0] == '#') continue;
        const auto fields = split_tabs(line);
        if (fields.size() < 10 ||
            (fields[1] != "17616778" && fields[1] != "17625527"))
            continue;
        const size_t separator = fields[9].rfind(':');
        REQUIRE(separator != std::string::npos);
        phase_sets[fields[1]] = fields[9].substr(separator + 1);
    }
    REQUIRE(phase_sets.size() == 2);
    CHECK(phase_sets.at("17616778") != ".");
    CHECK(phase_sets.at("17616778") == phase_sets.at("17625527"));

    Outcome reads;
    score_bam(dir + "/phased.bam", w, load_truth(p.truth_map),
              input_read_spans(p, w), reads);
    CHECK_FALSE(reads.switched);
    CHECK(reads.concordance() >= 0.80);
}

TEST_CASE("BAM block attaches to one supported graph flank",
          "[gap][stitch-connectivity]") {
    const Paths p = paths();
    if (!p.complete()) {
        WARN("gap-window tests need inputs that are absent: " << p.missing());
        SUCCEED("skipped: inputs absent");
        return;
    }
    Window w;
    w.gap_left = 14050001;
    w.gap_right = 14950000;
    std::string dir;
    REQUIRE(run_arm(p, w, "stitch_one_flank", "", dir));

    std::ifstream in(dir + "/native.vcf");
    REQUIRE(in.good());
    std::map<std::string, std::string> phase_sets;
    std::string line;
    while (std::getline(in, line)) {
        if (line.empty() || line[0] == '#') continue;
        const auto fields = split_tabs(line);
        if (fields.size() < 10 ||
            (fields[1] != "14577646" && fields[1] != "14584522" &&
             fields[1] != "14612632"))
            continue;
        const size_t separator = fields[9].rfind(':');
        REQUIRE(separator != std::string::npos);
        phase_sets[fields[1]] = fields[9].substr(separator + 1);
    }
    REQUIRE(phase_sets.size() == 3);
    CHECK(phase_sets.at("14577646") != ".");
    CHECK(phase_sets.at("14577646") == phase_sets.at("14584522"));
    CHECK(phase_sets.at("14584522") != phase_sets.at("14612632"));
}

TEST_CASE("BAM recovery uses a supported run before a weak source cut",
          "[gap][stitch-connectivity]") {
    const Paths p = paths();
    if (!p.complete()) {
        WARN("gap-window tests need inputs that are absent: " << p.missing());
        SUCCEED("skipped: inputs absent");
        return;
    }
    Window w;
    w.gap_left = 8050001;
    w.gap_right = 8950000;
    std::string dir;
    REQUIRE(run_arm(p, w, "stitch_source_run", "", dir));

    std::ifstream in(dir + "/native.vcf");
    REQUIRE(in.good());
    std::map<std::string, std::string> phase_sets;
    std::string line;
    while (std::getline(in, line)) {
        if (line.empty() || line[0] == '#') continue;
        const auto fields = split_tabs(line);
        if (fields.size() < 10 ||
            (fields[1] != "8632381" && fields[1] != "8638940" &&
             fields[1] != "8662670"))
            continue;
        const size_t separator = fields[9].rfind(':');
        REQUIRE(separator != std::string::npos);
        phase_sets[fields[1]] = fields[9].substr(separator + 1);
    }
    REQUIRE(phase_sets.size() == 3);
    CHECK(phase_sets.at("8632381") != ".");
    CHECK(phase_sets.at("8632381") == phase_sets.at("8638940"));
    CHECK(phase_sets.at("8638940") == phase_sets.at("8662670"));
    Window gap;
    gap.gap_left = 8638940;
    gap.gap_right = 8662670;
    Outcome reads;
    score_bam(dir + "/phased.bam", gap, load_truth(p.truth_map),
              input_read_spans(p, gap), reads);
    CHECK_FALSE(reads.switched);
    CHECK(reads.concordance() >= 0.99);
}

TEST_CASE("BAM source path with one-haplotype molecule support closes 34.844 Mb",
          "[gap][stitch-connectivity]") {
    const Paths p = paths();
    if (!p.complete()) {
        WARN("gap-window tests need inputs that are absent: " << p.missing());
        SUCCEED("skipped: inputs absent");
        return;
    }
    Window w;
    w.gap_left = 34050001;
    w.gap_right = 34950000;
    std::string dir;
    REQUIRE(run_arm(p, w, "stitch_one_hap_source", "", dir));

    std::ifstream in(dir + "/native.vcf");
    REQUIRE(in.good());
    std::map<std::string, std::string> phase_sets;
    std::string line;
    while (std::getline(in, line)) {
        if (line.empty() || line[0] == '#') continue;
        const auto fields = split_tabs(line);
        if (fields.size() < 10 ||
            (fields[1] != "34844194" && fields[1] != "34844579"))
            continue;
        const size_t separator = fields[9].rfind(':');
        REQUIRE(separator != std::string::npos);
        phase_sets[fields[1]] = fields[9].substr(separator + 1);
    }
    REQUIRE(phase_sets.size() == 2);
    CHECK(phase_sets.at("34844194") != ".");
    CHECK(phase_sets.at("34844194") == phase_sets.at("34844579"));

    Window bridge;
    bridge.gap_left = 34844194;
    bridge.gap_right = 34844579;
    Outcome reads;
    score_bam(dir + "/phased.bam", bridge, load_truth(p.truth_map),
              input_read_spans(p, bridge), reads);
    CHECK_FALSE(reads.switched);
    CHECK(reads.concordance() >= 0.98);
}

TEST_CASE("clean indel boundaries use the phased read path",
          "[gap][stitch-connectivity]") {
    const Paths p = paths();
    if (!p.complete()) {
        WARN("gap-window tests need inputs that are absent: " << p.missing());
        SUCCEED("skipped: inputs absent");
        return;
    }
    struct Bridge {
        long long region_left;
        long long region_right;
        long long left;
        long long right;
        bool same_allele;
    };
    const std::array<Bridge, 2> bridges{{
        {14050001, 14950000, 14235594, 14239930, true},
        {36050001, 36950000, 36016138, 36018966, false},
    }};
    for (const Bridge& bridge : bridges) {
        INFO("boundary " << bridge.left << " -> " << bridge.right);
        Window region;
        region.gap_left = bridge.region_left;
        region.gap_right = bridge.region_right;
        std::string dir;
        REQUIRE(run_arm(p, region, "stitch_clean_indel", "", dir));

        std::ifstream in(dir + "/native.vcf");
        REQUIRE(in.good());
        std::map<long long, std::pair<std::string, std::string>> sites;
        std::string line;
        while (std::getline(in, line)) {
            if (line.empty() || line[0] == '#') continue;
            const auto fields = split_tabs(line);
            if (fields.size() < 10) continue;
            const long long pos = std::stoll(fields[1]);
            if (pos != bridge.left && pos != bridge.right) continue;
            const size_t separator = fields[9].rfind(':');
            REQUIRE(separator != std::string::npos);
            sites[pos] = {fields[9].substr(0, 3),
                          fields[9].substr(separator + 1)};
        }
        REQUIRE(sites.size() == 2);
        CHECK(sites.at(bridge.left).second != ".");
        CHECK(sites.at(bridge.left).second == sites.at(bridge.right).second);
        CHECK((sites.at(bridge.left).first == sites.at(bridge.right).first) ==
              bridge.same_allele);

        Window gap;
        gap.gap_left = bridge.left;
        gap.gap_right = bridge.right;
        Outcome reads;
        score_bam(dir + "/phased.bam", gap, load_truth(p.truth_map),
                  input_read_spans(p, gap), reads);
        CHECK_FALSE(reads.switched);
        CHECK(reads.concordance() >= 0.95);
    }
}

TEST_CASE("one weak-cut BAM run attaches to one supported neighbor",
          "[gap][stitch-connectivity]") {
    const Paths p = paths();
    if (!p.complete()) {
        WARN("gap-window tests need inputs that are absent: " << p.missing());
        SUCCEED("skipped: inputs absent");
        return;
    }
    struct Bridge {
        long long region_left;
        long long region_right;
        long long left;
        long long right;
        size_t left_rows;
    };
    const std::array<Bridge, 2> bridges{{
        {50001, 950000, 542052, 545002, 1},
        {37050001, 37950000, 37461999, 37466820, 2},
    }};
    for (const Bridge& bridge : bridges) {
        INFO("weak-cut boundary " << bridge.left << " -> " << bridge.right);
        Window region;
        region.gap_left = bridge.region_left;
        region.gap_right = bridge.region_right;
        std::string dir;
        REQUIRE(run_arm(p, region, "stitch_weak_run", "", dir));
        std::ifstream in(dir + "/native.vcf");
        REQUIRE(in.good());
        std::vector<std::pair<std::string, std::string>> left_rows;
        std::vector<std::pair<std::string, std::string>> right_rows;
        std::map<long long, std::pair<std::string, std::string>> adjacent_pair;
        std::string line;
        while (std::getline(in, line)) {
            if (line.empty() || line[0] == '#') continue;
            const auto fields = split_tabs(line);
            if (fields.size() < 10) continue;
            const long long pos = std::stoll(fields[1]);
            if (bridge.left == 542052 && (pos == 528827 || pos == 528828)) {
                const size_t separator = fields[9].rfind(':');
                REQUIRE(separator != std::string::npos);
                adjacent_pair[pos] = {fields[9].substr(0, 3),
                                      fields[9].substr(separator + 1)};
            }
            if (pos != bridge.left && pos != bridge.right) continue;
            const size_t separator = fields[9].rfind(':');
            REQUIRE(separator != std::string::npos);
            auto& rows = pos == bridge.left ? left_rows : right_rows;
            rows.emplace_back(fields[9].substr(0, 3),
                              fields[9].substr(separator + 1));
        }
        REQUIRE(left_rows.size() == bridge.left_rows);
        REQUIRE(right_rows.size() == 1);
        for (const auto& left : left_rows) {
            CHECK(left.second != ".");
            CHECK(left.second == right_rows.front().second);
        }
        CHECK(left_rows.front().first != right_rows.front().first);
        if (bridge.left_rows == 2)
            CHECK(left_rows[0].first != left_rows[1].first);
        if (bridge.left == 542052) {
            REQUIRE(adjacent_pair.size() == 2);
            // The child-snarl SNP is now projected onto the BAM deletion
            // block. Its ALT and the insertion ALT must occupy opposite
            // haplotypes in the same phase set.
            const auto& snp = adjacent_pair.at(528827);
            const auto& insertion = adjacent_pair.at(528828);
            CHECK(snp.second == insertion.second);
            CHECK(snp.first != insertion.first);
            Window gap;
            gap.gap_left = 528827;
            gap.gap_right = 528828;
            Outcome reads;
            score_bam(dir + "/phased.bam", gap, load_truth(p.truth_map),
                      input_read_spans(p, gap), reads);
            CHECK_FALSE(reads.switched);
            CHECK(reads.concordance() >= 0.95);
        }
    }
}

TEST_CASE("complete BAM blocks use split-stable allele votes without a read gauge",
          "[gap][stitch-connectivity]") {
    const Paths p = paths();
    if (!p.complete()) {
        WARN("gap-window tests need inputs that are absent: " << p.missing());
        SUCCEED("skipped: inputs absent");
        return;
    }
    Window region;
    region.gap_left = 32050001;
    region.gap_right = 32950000;
    std::string dir;
    REQUIRE(run_arm(p, region, "stitch_aggregate", "", dir));

    std::ifstream in(dir + "/native.vcf");
    REQUIRE(in.good());
    std::map<long long, std::pair<std::string, std::string>> sites;
    std::string line;
    while (std::getline(in, line)) {
        if (line.empty() || line[0] == '#') continue;
        const auto fields = split_tabs(line);
        if (fields.size() < 10) continue;
        const long long pos = std::stoll(fields[1]);
        if (pos != 32490058 && pos != 32490150) continue;
        const size_t separator = fields[9].rfind(':');
        REQUIRE(separator != std::string::npos);
        sites[pos] = {fields[9].substr(0, 3),
                      fields[9].substr(separator + 1)};
    }
    REQUIRE(sites.size() == 2);
    CHECK(sites.at(32490058).second != ".");
    CHECK(sites.at(32490058).second == sites.at(32490150).second);
    CHECK(sites.at(32490058).first == sites.at(32490150).first);
}

TEST_CASE("BAM-supported inner block connects after an invalid graph SNP is excluded",
          "[gap][stitch-connectivity]") {
    const Paths p = paths();
    if (!p.complete()) {
        WARN("gap-window tests need inputs that are absent: " << p.missing());
        SUCCEED("skipped: inputs absent");
        return;
    }
    Window w;
    w.gap_left = 55300000;
    w.gap_right = 55500000;
    std::string dir;
    REQUIRE(run_arm(p, w, "stitch_bam_inner", "", dir));

    std::ifstream in(dir + "/native.vcf");
    REQUIRE(in.good());
    const std::set<std::string> sites{
        "55331033", "55360776", "55373606", "55381221", "55381472",
        "55382729"};
    std::map<std::string, std::string> phase_sets;
    std::string line;
    while (std::getline(in, line)) {
        if (line.empty() || line[0] == '#') continue;
        const auto fields = split_tabs(line);
        if (fields.size() < 10 || sites.count(fields[1]) == 0) continue;
        if (fields[8].find("PS") == std::string::npos) {
            phase_sets[fields[1]] = ".";
        } else {
            const size_t separator = fields[9].rfind(':');
            REQUIRE(separator != std::string::npos);
            phase_sets[fields[1]] = fields[9].substr(separator + 1);
        }
    }
    REQUIRE(phase_sets.count("55331033") == 1);
    REQUIRE(phase_sets.count("55360776") == 1);
    REQUIRE(phase_sets.count("55382729") == 1);
    CHECK(phase_sets.at("55331033") != phase_sets.at("55360776"));
    // Physical reads carry T or a deletion at this graph G>T site, with no
    // callable G. Exclude the false anchor while retaining the BAM bridge.
    CHECK(phase_sets.count("55373606") == 0);
    REQUIRE(phase_sets.count("55381221") == 1);
    REQUIRE(phase_sets.count("55381472") == 1);
    CHECK(phase_sets.at("55360776") == phase_sets.at("55381221"));
    CHECK(phase_sets.at("55381221") == phase_sets.at("55381472"));
    CHECK(phase_sets.at("55381472") == phase_sets.at("55382729"));

    Window bridge;
    bridge.gap_left = 55360776;
    bridge.gap_right = 55382729;
    Outcome reads;
    score_bam(dir + "/phased.bam", bridge, load_truth(p.truth_map),
              input_read_spans(p, bridge), reads);
    CHECK_FALSE(reads.switched);
    CHECK(reads.concordance() >= 0.95);
}

TEST_CASE("low-MAPQ BAM evidence cannot veto a graph heterozygote",
          "[gap][anchor-quality]") {
    const Paths p = paths();
    if (!p.complete()) {
        WARN("gap-window tests need inputs that are absent: " << p.missing());
        SUCCEED("skipped: inputs absent");
        return;
    }
    Window w;
    w.gap_left = 30750000;
    w.gap_right = 30810000;
    std::string dir;
    REQUIRE(run_arm(p, w, "anchor_quality", "", dir));

    std::ifstream in(dir + "/native.vcf");
    REQUIRE(in.good());
    std::string line;
    bool found = false;
    while (std::getline(in, line)) {
        if (line.empty() || line[0] == '#') continue;
        const auto fields = split_tabs(line);
        if (fields.size() < 10 || fields[1] != "30794399") continue;
        found = true;
        CHECK(fields[9].find('|') != std::string::npos);
        CHECK(fields[8].find("PS") != std::string::npos);
    }
    CHECK(found);
}

TEST_CASE("shifted BAM deletion joins a certified graph suffix",
          "[gap][stitch-connectivity]") {
    const Paths p = paths();
    if (!p.complete()) {
        WARN("gap-window tests need inputs that are absent: " << p.missing());
        SUCCEED("skipped: inputs absent");
        return;
    }
    // The shifted BAM deletion joins the right graph SNP. Independent
    // primary calls also certify the next graph deletion and the local SNP
    // suffix; the earlier weak graph edge remains split.
    Window region;
    region.gap_left = 13050001;
    region.gap_right = 13950000;
    std::string dir;
    REQUIRE(run_arm(p, region, "stitch_shifted_suffix_13m", "", dir));
    const std::set<std::string> boundaries{
        "13752640", "13773452", "13804846", "13830800", "13844727"};
    std::map<std::string, std::string> phase_sets;
    std::ifstream vcf(dir + "/native.vcf");
    REQUIRE(vcf.good());
    std::string line;
    while (std::getline(vcf, line)) {
        if (line.empty() || line[0] == '#') continue;
        const auto fields = split_tabs(line);
        if (fields.size() < 10 || boundaries.count(fields[1]) == 0)
            continue;
        const size_t separator = fields[9].rfind(':');
        REQUIRE(separator != std::string::npos);
        phase_sets[fields[1]] = fields[9].substr(separator + 1);
    }
    REQUIRE(phase_sets.size() == boundaries.size());
    CHECK(phase_sets.at("13752640") != phase_sets.at("13773452"));
    CHECK(phase_sets.at("13773452") == phase_sets.at("13804846"));
    CHECK(phase_sets.at("13804846") == phase_sets.at("13830800"));
    CHECK(phase_sets.at("13830800") == phase_sets.at("13844727"));

    Window gap;
    gap.gap_left = 13830800;
    gap.gap_right = 13844727;
    Outcome reads;
    score_bam(dir + "/phased.bam", gap, load_truth(p.truth_map),
              input_read_spans(p, gap), reads);
    CHECK_FALSE(reads.switched);
    CHECK(reads.concordance() >= 0.93);
}

TEST_CASE("a directly supported BAM suffix closes the 17.865 Mb gap",
          "[gap][stitch-connectivity]") {
    const Paths p = paths();
    if (!p.complete()) {
        WARN("gap-window tests need inputs that are absent: " << p.missing());
        SUCCEED("skipped: inputs absent");
        return;
    }
    // The BAM source's final weak cut lies between this clean SNP and MSA
    // insertion. Molecules calling both sites certify their common phase;
    // the earlier source prefix must retain its independent gauge.
    Window region;
    region.gap_left = 17050001;
    region.gap_right = 17950000;
    std::string dir;
    REQUIRE(run_arm(p, region, "stitch_bam_suffix_17m", "", dir));
    const std::set<std::string> positions{
        "17839399", "17852024", "17865146", "17883198"};
    std::map<std::string, std::string> phase_sets;
    std::ifstream vcf(dir + "/native.vcf");
    REQUIRE(vcf.good());
    std::string line;
    while (std::getline(vcf, line)) {
        if (line.empty() || line[0] == '#') continue;
        const auto fields = split_tabs(line);
        if (fields.size() < 10 || positions.count(fields[1]) == 0)
            continue;
        const size_t separator = fields[9].rfind(':');
        REQUIRE(separator != std::string::npos);
        phase_sets[fields[1]] = fields[9].substr(separator + 1);
    }
    REQUIRE(phase_sets.size() == positions.size());
    CHECK(phase_sets.at("17839399") != phase_sets.at("17852024"));
    CHECK(phase_sets.at("17852024") == phase_sets.at("17865146"));
    CHECK(phase_sets.at("17865146") == phase_sets.at("17883198"));

    Window gap;
    gap.gap_left = 17865146;
    gap.gap_right = 17883198;
    Outcome reads;
    score_bam(dir + "/phased.bam", gap, load_truth(p.truth_map),
              input_read_spans(p, gap), reads);
    CHECK_FALSE(reads.switched);
    CHECK(reads.concordance() >= 0.99);
}

TEST_CASE("a newly imported BAM seam closes the 15.056 Mb gap",
          "[gap][stitch-connectivity]") {
    const Paths p = paths();
    if (!p.complete()) {
        WARN("gap-window tests need inputs that are absent: " << p.missing());
        SUCCEED("skipped: inputs absent");
        return;
    }
    // Reproduce the owning chunk. The first solve imports two independent
    // BAM blocks inside a wider graph seam; the physical SNP pair permits a
    // second solve across their newly exposed 15 kb boundary.
    Window region;
    region.gap_left = 15000001;
    region.gap_right = 15950000;
    std::string dir;
    REQUIRE(run_arm(p, region, "stitch_new_bam_seam_15m", "", dir));

    std::ifstream in(dir + "/native.vcf");
    REQUIRE(in.good());
    const std::set<std::string> boundaries{
        "15056025", "15071132", "15071145"};
    std::map<std::string, std::string> phase_sets;
    std::string line;
    while (std::getline(in, line)) {
        if (line.empty() || line[0] == '#') continue;
        const auto fields = split_tabs(line);
        if (fields.size() < 10 || boundaries.count(fields[1]) == 0)
            continue;
        const size_t separator = fields[9].rfind(':');
        REQUIRE(separator != std::string::npos);
        phase_sets[fields[1]] = fields[9].substr(separator + 1);
    }
    REQUIRE(phase_sets.size() == boundaries.size());
    CHECK(phase_sets.at("15056025") != ".");
    CHECK(phase_sets.at("15056025") == phase_sets.at("15071132"));
    CHECK(phase_sets.at("15071132") == phase_sets.at("15071145"));

    Window gap;
    gap.gap_left = 15056025;
    gap.gap_right = 15071132;
    Outcome reads;
    score_bam(dir + "/phased.bam", gap, load_truth(p.truth_map),
              input_read_spans(p, gap), reads);
    CHECK_FALSE(reads.switched);
    CHECK(reads.concordance() >= 0.98);
}

TEST_CASE("recovered right deletion joins a source-backed graph block",
          "[gap][stitch-connectivity]") {
    const Paths p = paths();
    if (!p.complete()) {
        WARN("gap-window tests need inputs that are absent: " << p.missing());
        SUCCEED("skipped: inputs absent");
        return;
    }
    Window gap;
    gap.gap_left = 8166027;
    gap.gap_right = 8172072;
    std::string dir;
    REQUIRE(run_arm(p, gap, "stitch_source_graph_deletion_8m", "", dir));

    std::ifstream in(dir + "/native.vcf");
    REQUIRE(in.good());
    std::map<std::string, std::pair<std::string, std::string>> calls;
    std::string line;
    while (std::getline(in, line)) {
        if (line.empty() || line[0] == '#') continue;
        const auto fields = split_tabs(line);
        if (fields.size() < 10) continue;
        const std::string key = fields[1] + ":" + fields[3] + ":" + fields[4];
        if (key != "8165210:T:A" && key != "8166027:C:CTTTT" &&
            key != "8166027:C:CTTTTTTT" && key != "8172072:GAGGAA:G" &&
            key != "8176552:G:A")
            continue;
        const size_t separator = fields[9].rfind(':');
        REQUIRE(separator != std::string::npos);
        calls[key] = {fields[9].substr(0, 3), fields[9].substr(separator + 1)};
    }
    REQUIRE(calls.size() == 5);
    const std::string phase_set = calls.at("8165210:T:A").second;
    CHECK(phase_set != ".");
    for (const auto& [key, call] : calls) {
        INFO(key);
        CHECK(is_phased_het(call.first));
        CHECK(call.second == phase_set);
    }
    CHECK(calls.at("8165210:T:A").first !=
          calls.at("8176552:G:A").first);
    CHECK(calls.at("8166027:C:CTTTT").first !=
          calls.at("8166027:C:CTTTTTTT").first);

    Outcome reads;
    score_bam(dir + "/phased.bam", gap, load_truth(p.truth_map),
              input_read_spans(p, gap), reads);
    CHECK_FALSE(reads.switched);
    CHECK(reads.concordance() >= 0.95);
}

TEST_CASE("adjacent anchors do not trigger a second BAM recovery solve",
          "[gap][stitch-connectivity]") {
    const Paths p = paths();
    if (!p.complete()) {
        WARN("gap-window tests need inputs that are absent: " << p.missing());
        SUCCEED("skipped: inputs absent");
        return;
    }
    Window region;
    region.gap_left = 17000001;
    region.gap_right = 17950000;
    std::string dir;
    REQUIRE(run_arm(p, region, "stitch_no_adjacent_retry_17m", "", dir));

    std::ifstream in(dir + "/native.vcf");
    REQUIRE(in.good());
    std::map<std::string, std::string> phase_sets;
    std::string line;
    while (std::getline(in, line)) {
        if (line.empty() || line[0] == '#') continue;
        const auto fields = split_tabs(line);
        if (fields.size() < 10 ||
            (fields[1] != "17839398" && fields[1] != "17839399"))
            continue;
        const size_t separator = fields[9].rfind(':');
        REQUIRE(separator != std::string::npos);
        phase_sets[fields[1]] = fields[9].substr(separator + 1);
    }
    REQUIRE(phase_sets.size() == 2);
    CHECK(phase_sets.at("17839398") != phase_sets.at("17839399"));

    Window gap;
    gap.gap_left = 17839398;
    gap.gap_right = 17839399;
    Outcome reads;
    score_bam(dir + "/phased.bam", gap, load_truth(p.truth_map),
              input_read_spans(p, gap), reads);
    CHECK_FALSE(reads.switched);
    CHECK(reads.concordance() >= 0.99);
}

TEST_CASE("shifted single-base BAM insertion closes its owning 23 Mb gap",
          "[gap][stitch-connectivity][representation]") {
    const Paths p = paths();
    if (!p.complete()) {
        WARN("gap-window tests need inputs that are absent: " << p.missing());
        SUCCEED("skipped: inputs absent");
        return;
    }
    Window chunk;
    chunk.gap_left = 23050001;
    chunk.gap_right = 23950000;
    std::string dir;
    REQUIRE(run_arm(p, chunk, "stitch_shifted_insertion_23m", "", dir));
    std::ifstream vcf(dir + "/native.vcf");
    REQUIRE(vcf.good());
    std::map<std::string, std::pair<std::string, std::string>> alleles;
    std::string line;
    while (std::getline(vcf, line)) {
        if (line.empty() || line[0] == '#') continue;
        const auto fields = split_tabs(line);
        if (fields.size() < 10) continue;
        const std::string key = fields[1] + ":" + fields[3] + ">" + fields[4];
        if (key != "23792418:T>TT" && key != "23806565:T>TAC" &&
            key != "23806565:TACACACAC>T") continue;
        const size_t colon = fields[9].rfind(':');
        REQUIRE(colon != std::string::npos);
        alleles[key] = {fields[9].substr(0, fields[9].find(':')),
                        fields[9].substr(colon + 1)};
    }
    REQUIRE(alleles.size() == 3);
    const auto& left = alleles.at("23792418:T>TT");
    const auto& insertion = alleles.at("23806565:T>TAC");
    const auto& deletion = alleles.at("23806565:TACACACAC>T");
    REQUIRE(is_phased_het(left.first));
    CHECK(left.second != ".");
    CHECK(left.second == insertion.second);
    CHECK(left.second == deletion.second);
    CHECK(left.first == insertion.first);
    CHECK(insertion.first != deletion.first);

    Window gap;
    gap.gap_left = 23792418;
    gap.gap_right = 23806565;
    Outcome reads;
    score_bam(dir + "/phased.bam", gap, load_truth(p.truth_map),
              input_read_spans(p, gap), reads);
    CHECK_FALSE(reads.switched);
    CHECK(reads.concordance() >= 0.97);
    CHECK(reads.separated() >= 0.50);
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
