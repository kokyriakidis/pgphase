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
#include <sys/stat.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <map>
#include <memory>
#include <set>
#include <sstream>
#include <string>
#include <tuple>
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
    int primary_scorable = 0;
    int primary_correct = 0;
    int core_correct = 0;
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

const std::unordered_map<std::string, char>& load_truth(const std::string& path) {
    struct stat metadata;
    REQUIRE(stat(path.c_str(), &metadata) == 0);
    const auto key = std::make_tuple(path, metadata.st_dev, metadata.st_ino,
        metadata.st_size, metadata.st_mtim.tv_sec, metadata.st_mtim.tv_nsec,
        metadata.st_ctim.tv_sec, metadata.st_ctim.tv_nsec);
    static std::map<decltype(key), std::unordered_map<std::string, char>> cache;
    const auto hit = cache.find(key);
    if (hit != cache.end()) return hit->second;
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
    return cache.emplace(key, std::move(out)).first->second;
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
    std::unordered_map<long long, std::pair<int, int>> primary_votes;
    std::unordered_map<std::string, std::pair<long long, bool>> primary_tags;
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
        if (!(rec->core.flag & (BAM_FSECONDARY | BAM_FSUPPLEMENTARY)) &&
            (hap == 1 || hap == 2) && set_id > 0) {
            const std::string name = bam_get_qname(rec);
            if (!primary_tags.emplace(name, std::make_pair(set_id, mat_on_hap1)).second)
                FAIL("duplicate scored primary read: " << name);
            auto& pv = primary_votes[set_id];
            if (mat_on_hap1) ++pv.first; else ++pv.second;
        }
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

    // Use each whole replay block's parental orientation, then score exactly
    // the primary input molecules overlapping the gap. Abstentions stay in the
    // denominator; rescue phase sets cannot count as a connected core.
    constexpr long long kReadRescuePhaseSetOffset = 1000000000;
    std::unordered_map<long long, int> core_counts;
    for (const auto& [name, span] : spans) {
        if (span.second < w.gap_left || span.first >= w.gap_right || !truth.count(name))
            continue;
        ++out.primary_scorable;
        const auto tag = primary_tags.find(name);
        if (tag == primary_tags.end()) continue;
        const auto& vote = primary_votes.at(tag->second.first);
        if (tag->second.second != (vote.first >= vote.second)) continue;
        ++out.primary_correct;
        if (tag->second.first < kReadRescuePhaseSetOffset)
            ++core_counts[tag->second.first];
    }
    for (const auto& [ps, count] : core_counts) {
        (void)ps;
        out.core_correct = std::max(out.core_correct, count);
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

/// Keep production-context overrides reviewable beside the gap panel.
std::pair<long long, long long> replay_region(const Window& w) {
    static const auto regions = [] {
        std::map<long long, std::pair<long long, long long>> out;
        std::ifstream in("src/test_gap_replays.tsv");
        REQUIRE(in.good());
        std::string line;
        while (std::getline(in, line)) {
            if (line.empty() || line[0] == '#' || line.rfind("gap_left", 0) == 0)
                continue;
            const auto fields = split_tabs(line);
            REQUIRE(fields.size() == 3);
            const long long left = std::stoll(fields[0]);
            const auto region = std::make_pair(std::stoll(fields[1]), std::stoll(fields[2]));
            REQUIRE(region.first <= region.second);
            REQUIRE(out.emplace(left, region).second);
        }
        return out;
    }();
    const auto it = regions.find(w.gap_left);
    return it == regions.end()
        ? std::make_pair(w.gap_left - 50000, w.gap_right + 50000) : it->second;
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
        << " -r 'CHM13#0#chr20:" << replay_region(w).first
        << "-" << replay_region(w).second << "'"
        << " -t " << test_threads() << " " << flags
        << " -o '" << outdir << "/candidates.tsv'"
        << " --phased-vcf-out '" << outdir << "/native.vcf'"
        << " --phased-bam-out '" << outdir << "/phased.bam'"
        << " > '" << outdir << "/stdout.log' 2> '" << outdir << "/stderr.log'";
    // Persist completed pipeline outputs, but always execute the assertions.
    // The helper fingerprints the binary and inputs and locks across shards.
    const auto shell_quote = [](const std::string& value) {
        std::string quoted = "'";
        for (const char c : value)
            quoted += c == '\'' ? "'\\''" : std::string(1, c);
        return quoted + "'";
    };
    const std::string cached_command =
        "env -u PGPHASE_GAP_FILTER -u PGPHASE_GAP_BENCHMARK -u PGPHASE_GAP_CERTIFIED "
        "python3 scripts/cache_gap_replay.py --outdir " + shell_quote(outdir) +
        " --command " + shell_quote(cmd.str());
    return std::system(cached_command.c_str()) == 0;
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

/// Read bounds shared by owning-context regressions; '-' means no such bound.
void check_read_floors(const std::string& name, const Outcome& got) {
    static const auto floors = [] {
        std::map<std::string, std::vector<std::string>> out;
        std::ifstream in("src/test_gap_read_floors.tsv");
        REQUIRE(in.good());
        std::string line;
        while (std::getline(in, line)) {
            if (line.empty() || line[0] == '#' || line.rfind("check", 0) == 0) continue;
            auto fields = split_tabs(line);
            REQUIRE(fields.size() == 7);
            const std::string name = fields.front();
            fields.erase(fields.begin());
            REQUIRE(out.emplace(name, std::move(fields)).second);
        }
        return out;
    }();
    const auto it = floors.find(name);
    REQUIRE(it != floors.end());
    const auto& row = it->second;
    const auto value = [](const std::string& text) {
        const auto slash = text.find('/');
        return slash == std::string::npos ? std::stod(text) :
            std::stod(text.substr(0, slash)) / std::stod(text.substr(slash + 1));
    };
    INFO("read regression floors for " << name);
    if (row[0] != "-") CHECK(got.scored >= std::stoi(row[0]));
    if (row[1] != "-") CHECK(got.correct >= std::stoi(row[1]));
    if (row[2] != "-") CHECK(got.discordant() <= std::stoi(row[2]));
    if (row[3] != "-") CHECK(got.concordance() >= value(row[3]));
    if (row[4] != "-") CHECK(got.separated() >= value(row[4]));
    if (row[5] != "-") CHECK(got.window_scorable == std::stoi(row[5]));
}

struct GapBenchmark {
    int scorable = 0;
    int correct = 0;
    int core_correct = 0;
};

const std::map<std::string, GapBenchmark>& gap_benchmarks() {
    static const auto rows = [] {
        std::map<std::string, GapBenchmark> out;
        std::ifstream in(env_or("PGPHASE_GAP_BENCHMARK", "src/test_gap_hiphase.tsv"));
        REQUIRE(in.good());
        std::string line;
        while (std::getline(in, line)) {
            if (line.empty() || line[0] == '#' || line.rfind("window", 0) == 0) continue;
            const auto fields = split_tabs(line);
            REQUIRE(fields.size() == 4);
            REQUIRE(out.emplace(fields[0], GapBenchmark{std::stoi(fields[1]),
                std::stoi(fields[2]), std::stoi(fields[3])}).second);
        }
        return out;
    }();
    return rows;
}

const std::set<std::string>& certified_gaps() {
    static const auto gaps = [] {
        std::set<std::string> out;
        std::ifstream in(env_or("PGPHASE_GAP_CERTIFIED", "src/test_gap_certified.tsv"));
        REQUIRE(in.good());
        std::string line;
        while (std::getline(in, line)) {
            if (line.empty() || line[0] == '#' || line == "window") continue;
            REQUIRE(out.insert(line).second);
        }
        return out;
    }();
    return gaps;
}

void check_gap_contract(const Paths& p, const Window& w, const Outcome& got) {
    const auto key = window_key(w);
    const auto benchmark = gap_benchmarks().find(key);
    REQUIRE(benchmark != gap_benchmarks().end());
    REQUIRE(got.primary_scorable == benchmark->second.scorable);
    constexpr double kMinCorrectPrimaryFraction = 0.80;
    const bool quality = got.primary_scorable > 0 &&
        static_cast<double>(got.primary_correct) / got.primary_scorable >=
        kMinCorrectPrimaryFraction;
    const bool parity = got.primary_correct >= benchmark->second.correct &&
        got.core_correct >= benchmark->second.core_correct;
    INFO("primary correctness " << got.primary_correct << "/" << got.primary_scorable
         << "; core correct " << got.core_correct << "; HiPhase correct/core "
         << benchmark->second.correct << "/" << benchmark->second.core_correct);
    // Keep historical regression floors; the strict certification manifest
    // grows when a closure has been reviewed under the current contract.
    if (certified_gaps().count(key)) {
        CHECK(quality);
        CHECK(parity);
    }
    std::filesystem::create_directories(p.workdir);
    const std::string report_path = p.workdir + "/gap-contract.tsv";
    static std::set<std::string> initialized_reports;
    const bool first = initialized_reports.insert(report_path).second;
    std::ofstream report(report_path, first ? std::ios::trunc : std::ios::app);
    REQUIRE(report.good());
    if (first)
        report << "window\tspans\tscorable\tcorrect\tcore_correct\thiphase_correct"
                  "\thiphase_core_correct\tpasses_80pct\tpasses_parity\n";
    report << key << '\t' << got.spans << '\t' << got.primary_scorable << '\t'
           << got.primary_correct << '\t' << got.core_correct << '\t'
           << benchmark->second.correct << '\t' << benchmark->second.core_correct
           << '\t' << quality << '\t' << parity << '\n';
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

static void gap_original_msa_snp_dropout_is_retried_before_cigar_backfill(const Paths& p) {
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

static void gap_conflicting_sparse_pairs_retry_with_their_observed_binomial_tail(const Paths& p) {
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

static void gap_a_moved_deletion_does_not_certify_a_whole_right_block_join(const Paths& p) {
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

static void gap_focused_retry_certifies_its_path_after_cigar_backfill(const Paths& p) {
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

static void gap_an_earlier_focused_retry_does_not_hide_a_later_msa_dropout(const Paths& p) {
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

static void gap_complete_recovery_msa_blocks_retain_the_complex_left_flank(const Paths& p) {
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
    const auto& truth = load_truth(p.truth_map);
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

static void gap_isolated_shared_bam_genotype_survives_graph_repeat_demotion(const Paths& p) {
    // Discovery needs the complete owning chunk, including both source gauges.
    Window owner;
    owner.gap_left = 21050001;
    owner.gap_right = 21950000;
    std::string dir;
    REQUIRE(run_arm(p, owner, "isolated_shared_msa_21m", "", dir));

    std::ifstream vcf(dir + "/native.vcf");
    REQUIRE(vcf.good());
    std::map<std::string, std::map<std::string, std::string>> rows;
    std::string line;
    while (std::getline(vcf, line)) {
        if (line.empty() || line[0] == '#') continue;
        const auto fields = split_tabs(line);
        if (fields.size() < 10) continue;
        const std::string key = fields[1] + ":" + fields[3] + ">" + fields[4];
        if (key != "21823066:GT>G" && key != "21831480:CT>C" &&
            key != "21844359:T>TTATATATA") continue;
        std::vector<std::string> names, values;
        std::string value;
        std::istringstream format(fields[8]), sample(fields[9]);
        while (std::getline(format, value, ':')) names.push_back(value);
        while (std::getline(sample, value, ':')) values.push_back(value);
        REQUIRE(names.size() == values.size());
        REQUIRE(rows.count(key) == 0);
        auto& row = rows[key];
        for (size_t i = 0; i < names.size(); ++i) row[names[i]] = values[i];
    }
    REQUIRE(rows.size() == 3);
    const auto& left = rows.at("21823066:GT>G");
    const auto& middle = rows.at("21831480:CT>C");
    const auto& right = rows.at("21844359:T>TTATATATA");
    CHECK(is_phased_het(middle.at("GT")));
    CHECK(middle.at("AD") == "11,4");
    CHECK(middle.at("DP") == "15");
    CHECK(std::stoll(middle.at("PS")) > 0);
    // A shared source label does not certify the zero-pair left edge.
    CHECK(middle.at("PS") != left.at("PS"));
    CHECK(middle.at("PS") != right.at("PS"));
    CHECK(left.at("PS") != right.at("PS"));

    Window gap;
    gap.gap_left = 21823066;
    gap.gap_right = 21844359;
    Outcome reads;
    parse_vcf(dir + "/native.vcf", gap, reads);
    score_bam(dir + "/phased.bam", gap, load_truth(p.truth_map),
              input_read_spans(p, gap), reads);
    CHECK_FALSE(reads.spans);
    CHECK(reads.scored >= 3321);
    CHECK(reads.correct >= 3110);
    CHECK(reads.discordant() <= 211);
}

static void gap_chr20_gap_windows(const Paths& p) {
    const auto panel = load_panel(p.panel);
    const auto& truth = load_truth(p.truth_map);
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
                check_gap_contract(p, w, got);
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
                    // The audited BAM prefix reaches the complementary pair.
                    // Its local transfer enables the ordinary stitch to join
                    // the catalog block with the opposite deletion orientation.
                    std::ifstream vcf(p.workdir + "/graph/w" + window_key(w) +
                                      "/native.vcf");
                    REQUIRE(vcf.good());
                    std::string left_graph_ps, left_boundary_ps, right_boundary_ps;
                    std::string left_graph_gt, left_boundary_gt;
                    std::string line;
                    while (std::getline(vcf, line)) {
                        if (line.empty() || line[0] == '#') continue;
                        const auto fields = split_tabs(line);
                        if (fields.size() < 10) continue;
                        const auto colon = fields[9].rfind(':');
                        if (colon == std::string::npos) continue;
                        const std::string ps = fields[9].substr(colon + 1);
                        if (fields[1] == "64118182" && fields[3] == "T" &&
                            fields[4] == "C") {
                            left_graph_ps = ps;
                            left_graph_gt = fields[9].substr(0, 3);
                        }
                        if (fields[1] == "64144256" && fields[3] == "CA" &&
                            fields[4] == "C") {
                            left_boundary_ps = ps;
                            left_boundary_gt = fields[9].substr(0, 3);
                        }
                        if (fields[1] == "64144722" &&
                            fields[3] == "ATGGTGGGGG" && fields[4] == "A")
                            right_boundary_ps = ps;
                    }
                    REQUIRE(!left_graph_ps.empty());
                    REQUIRE(!left_boundary_ps.empty());
                    REQUIRE(!right_boundary_ps.empty());
                    CHECK(left_graph_ps == left_boundary_ps);
                    CHECK(left_graph_gt != left_boundary_gt);
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
                    // Read HP calls orient the short BAM run. The repeat-shifted
                    // graph insertion now joins through exact edit identity and
                    // independently supported graph SNP continuity.
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
                    CHECK(graph->second == left);
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

static void gap_recovery_preserves_complementary_bam_deletion_rows(const Paths& p) {
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

static void gap_verified_msa_insertion_survives_graph_output_classification(const Paths& p) {
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

static void gap_a_phased_bam_deletion_survives_an_unphased_graph_duplicate(const Paths& p) {
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

static void gap_graph_snp_retry_preserves_a_verified_bam_deletion_connection(const Paths& p) {
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

    const auto& truth = load_truth(p.truth_map);
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

static void gap_an_unsupported_bam_source_cut_keeps_the_far_deletion_independent(const Paths& p) {
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

static void gap_graph_snps_bridge_supported_recovery_blocks(const Paths& p) {
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

static void gap_recovery_stitch_preserves_the_next_flank_across_an_unlinked_seam(const Paths& p) {
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

static void gap_a_clean_snp_gap_bridge_keeps_the_prior_graph_gauge(const Paths& p) {
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
    const auto& truth = load_truth(p.truth_map);
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

static void gap_a_local_source_path_does_not_lose_its_graph_bridge(const Paths& p) {
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

static void gap_corroborated_snp_molecules_join_sparse_graph_seams(const Paths& p) {
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

static void gap_an_indel_boundary_with_allele_dropout_gets_an_msa_retry(const Paths& p) {
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

static void gap_a_lone_boundary_snp_pair_needs_a_repaired_complete_graph_flank(const Paths& p) {
    // One MAPQ-60, Q40 molecule crosses these two SNPs. Each is the only
    // clean SNP it observes in its block. Joining without repairing the
    // internal switch reverses hundreds of reads; the boundary pair alone
    // cannot supply that missing internal certificate.
    Window w;
    w.gap_left = 62050001;
    w.gap_right = 62950000;
    std::string dir;
    REQUIRE(run_arm(p, w, "stitch_lone_snp_pair", "", dir));
    std::ifstream vcf(dir + "/native.vcf");
    REQUIRE(vcf.good());
    std::map<std::string, std::pair<std::string, std::string>> sites;
    std::string line;
    while (std::getline(vcf, line)) {
        if (line.empty() || line[0] == '#') continue;
        const auto fields = split_tabs(line);
        if (fields.size() < 10 ||
            (fields[1] != "62623253" && fields[1] != "62645168" &&
             fields[1] != "62718395" && fields[1] != "62722021"))
            continue;
        const size_t separator = fields[9].rfind(':');
        REQUIRE(separator != std::string::npos);
        sites[fields[1]] = {fields[9].substr(0, 3),
                            fields[9].substr(separator + 1)};
    }
    REQUIRE(sites.size() == 4);
    if (sites.at("62623253").second == sites.at("62645168").second) {
        CHECK(is_phased_het(sites.at("62718395").first));
        CHECK(is_phased_het(sites.at("62722021").first));
        CHECK(sites.at("62718395").second == sites.at("62623253").second);
        CHECK(sites.at("62722021").second == sites.at("62623253").second);
        CHECK(sites.at("62718395").first != sites.at("62722021").first);
        Outcome reads;
        score_bam(dir + "/phased.bam", w, load_truth(p.truth_map),
                  input_read_spans(p, w), reads);
        CHECK_FALSE(reads.switched);
        CHECK(reads.scored >= 3750);
        CHECK(reads.correct >= 3722);
        CHECK(reads.discordant() <= 28);
    }
}

static void gap_physical_bridges_survive_their_owning_graph_chunks(const Paths& p) {
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
        // Physical deletion calls retain the suffix bridge. Independent
        // quality-backed source evidence now also certifies its prefix.
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
                CHECK(earlier_phase_sets.at("13752640") ==
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

static void gap_phased_right_reads_retain_the_certified_64_mb_bam_prefix(const Paths& p) {
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
    // Source transfer first preserves the cut-free BAM prefix; the ordinary
    // stitch then independently connects the preceding catalog block.
    CHECK(rows.at("64118182:T>C").second == left.second);
    CHECK(rows.at("64128828:G>C").second == left.second);
    CHECK(rows.at("64128828:G>C").first != rows.at("64134226:CT>C").first);
    CHECK(rows.at("64128828:G>C").first == rows.at("64134226:CTT>C").first);
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

    gap.gap_left = 64128828;
    gap.gap_right = 64134226;
    Outcome prefix_reads;
    parse_vcf(dir + "/native.vcf", gap, prefix_reads);
    score_bam(dir + "/phased.bam", gap, load_truth(p.truth_map),
              input_read_spans(p, gap), prefix_reads);
    CHECK(prefix_reads.spans);
    CHECK_FALSE(prefix_reads.switched);
    CHECK(prefix_reads.scored >= 3013);
    CHECK(prefix_reads.discordant() <= 3);
    CHECK(prefix_reads.concordance() >= 0.999);
    CHECK(prefix_reads.separated() >= 0.72);

    gap.gap_left = 64138752;
    gap.gap_right = 64140314;
    Outcome joined_reads;
    score_bam(dir + "/phased.bam", gap, load_truth(p.truth_map),
              input_read_spans(p, gap), joined_reads);
    CHECK_FALSE(joined_reads.switched);
    CHECK(joined_reads.concordance() >= 0.99);
    CHECK(joined_reads.separated() >= 0.71);
}

static void gap_source_allele_repair_preserves_complementary_repeat_rows(const Paths& p) {
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

static void gap_graph_haplotypes_bridge_equivalent_complementary_bam_deletions(const Paths& p) {
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

static void gap_physical_links_close_the_long_47_mb_gap_without_a_switch(const Paths& p) {
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

static void gap_validated_graph_bridge_survives_a_bam_phase_set_id_change(const Paths& p) {
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

static void gap_one_shared_bam_anchor_cannot_absorb_the_earlier_19_4_mb_block(const Paths& p) {
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

static void gap_direct_snp_proof_reuses_a_graph_block_at_56_15_mb(const Paths& p) {
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

static void gap_physical_snp_bridge_crosses_the_23_mb_graph_chunk_boundary(const Paths& p) {
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

static void gap_multiple_bam_source_labels_preserve_a_22_98_mb_graph_bridge(const Paths& p) {
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

static void gap_a_different_deletion_locus_cannot_validate_a_remapped_bridge(const Paths& p) {
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

static void gap_bam_transfer_preserves_an_already_connected_graph_phase_set(const Paths& p) {
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

static void gap_complementary_bam_boundary_rows_close_the_41_900_mb_seam(const Paths& p) {
    // The preceding seam shares a BAM solve region with this deletion pair.
    // The verified noisy-SNP bridge now also joins its outer graph blocks;
    // both connections must keep their allele gauge and parental orientation.
    Window full_chunk;
    full_chunk.gap_left = 41050001;
    full_chunk.gap_right = 41950000;
    std::string dir;
    REQUIRE(run_arm(p, full_chunk, "stitch_complementary_41m", "", dir));

    std::ifstream in(dir + "/native.vcf");
    REQUIRE(in.good());
    std::map<std::string, std::string> phase_sets;
    std::map<std::string, std::pair<std::string, std::string>> repeat_sites;
    std::string line;
    while (std::getline(in, line)) {
        if (line.empty() || line[0] == '#') continue;
        const auto fields = split_tabs(line);
        if (fields.size() < 10 ||
            (fields[1] != "41866917" && fields[1] != "41898323" &&
             fields[1] != "41900800" && fields[1] != "41919471" &&
             fields[1] != "41879449" && fields[1] != "41880908")) continue;
        const size_t separator = fields[9].rfind(':');
        REQUIRE(separator != std::string::npos);
        if (fields[1] == "41879449" || fields[1] == "41880908")
            repeat_sites.emplace(fields[1] + ":" + fields[3] + ">" + fields[4],
                std::make_pair(fields[9].substr(0, 3), fields[9].substr(separator + 1)));
        if (fields[1] != "41866917" && fields[1] != "41898323" &&
            fields[1] != "41900800" && fields[1] != "41919471") continue;
        phase_sets[fields[1]] = fields[9].substr(separator + 1);
    }
    REQUIRE(phase_sets.size() == 4);
    CHECK(phase_sets.at("41866917") == phase_sets.at("41898323"));
    CHECK(phase_sets.at("41900800") == phase_sets.at("41919471"));

    Window gap;
    gap.gap_left = 41900800;
    gap.gap_right = 41919471;
    Outcome reads;
    score_bam(dir + "/phased.bam", gap, load_truth(p.truth_map),
              input_read_spans(p, gap), reads);
    CHECK_FALSE(reads.switched);
    CHECK(reads.concordance() >= 0.98);

    // Local tandem-repeat calls must close the preceding private-site gap
    // without displacing the independently certified downstream source solve.
    REQUIRE(repeat_sites.size() == 3);
    const auto& left = repeat_sites.at("41879449:A>AGA");
    const auto& short_repeat = repeat_sites.at("41880908:T>TCACACACACACACA");
    const auto& long_repeat = repeat_sites.at("41880908:T>TCACACACACACACACACA");
    CHECK(is_phased_het(left.first));
    CHECK(left.second != ".");
    CHECK(left.second == short_repeat.second);
    CHECK(left.second == long_repeat.second);
    CHECK(left.first != short_repeat.first);
    CHECK(left.first == long_repeat.first);
    Window recovered;
    recovered.gap_left = 41879449;
    recovered.gap_right = 41880908;
    Outcome connectivity;
    parse_vcf(dir + "/native.vcf", recovered, connectivity);
    CHECK(connectivity.spans);
    Outcome whole_chunk;
    score_bam(dir + "/phased.bam", full_chunk, load_truth(p.truth_map),
              input_read_spans(p, full_chunk), whole_chunk);
    CHECK(whole_chunk.scored >= 4154);
    CHECK(whole_chunk.correct >= 4126);
    CHECK(whole_chunk.discordant() <= 28);
    CHECK_FALSE(whole_chunk.switched);
    CHECK(whole_chunk.concordance() >= 0.993);
}

static void gap_certified_deletion_bridge_preserves_the_11_599_mb_source_alleles(const Paths& p) {
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

    const auto& truth = load_truth(p.truth_map);
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

static void gap_complete_bam_source_path_closes_the_48_929_mb_seam_in_its_owning_chunk(const Paths& p) {
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

static void gap_complete_adjacent_graph_phase_sets_close_the_56_mb_seam_in_its_owning_chunk(const Paths& p) {
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

static void gap_complete_bam_path_and_verified_deletion_close_the_5_31_mb_seam(const Paths& p) {
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
    std::map<std::string, std::string> genotypes;
    std::string line;
    while (std::getline(vcf, line)) {
        if (line.empty() || line[0] == '#') continue;
        const auto fields = split_tabs(line);
        if (fields.size() < 10 ||
            (fields[1] != "5309406" && fields[1] != "5345085" &&
             fields[1] != "5350509" && fields[1] != "5393615" &&
             fields[1] != "5511231" && fields[1] != "5531924"))
            continue;
        const size_t separator = fields[9].rfind(':');
        REQUIRE(separator != std::string::npos);
        phase_sets[fields[1]] = fields[9].substr(separator + 1);
        genotypes[fields[1]] = fields[9].substr(0, 3);
    }
    REQUIRE(phase_sets.size() == 6);
    CHECK(phase_sets.at("5309406") != ".");
    CHECK(phase_sets.at("5309406") == phase_sets.at("5345085"));
    // Direct physical SNP calls join the local edge. The known source cut
    // farther right still leaves the distant graph block independent.
    CHECK(phase_sets.at("5345085") == phase_sets.at("5350509"));
    CHECK(phase_sets.at("5350509") != phase_sets.at("5393615"));

    // The right MSA deletion loses two partial-read calls in the larger
    // validation solve. Both original CIGARs agree with the source contrast;
    // backfill must recover them without crossing the upstream source cut.
    CHECK(phase_sets.at("5511231") == phase_sets.at("5531924"));
    CHECK(genotypes.at("5511231") != genotypes.at("5531924"));
    Window ng50_gap;
    ng50_gap.gap_left = 5511231;
    ng50_gap.gap_right = 5531924;
    Outcome ng50_reads;
    parse_vcf(dir + "/native.vcf", ng50_gap, ng50_reads);
    score_bam(dir + "/phased.bam", ng50_gap, load_truth(p.truth_map),
              input_read_spans(p, ng50_gap), ng50_reads);
    CHECK(ng50_reads.spans);
    CHECK_FALSE(ng50_reads.switched);
    CHECK(ng50_reads.scored >= 4024);
    CHECK(ng50_reads.correct >= 4001);
    CHECK(ng50_reads.discordant() <= 23);
    CHECK(ng50_reads.dominant_correct >= 131);

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

static void gap_repeat_deletion_recovery_preserves_the_61_738_mb_flank_gauges(const Paths& p) {
    // The raw source's deletion orientation conflicts with physical SNPs on
    // the left. Whole-chunk context catches the false joins from broad recall.
    Window chunk;
    chunk.gap_left = 61050001;
    chunk.gap_right = 61950000;
    std::string dir;
    REQUIRE(run_arm(p, chunk, "repeat_deletion_61m", "", dir));
    Window gap;
    gap.gap_left = 61738239;
    gap.gap_right = 61747506;
    Outcome reads;
    parse_vcf(dir + "/native.vcf", gap, reads);
    score_bam(dir + "/phased.bam", gap, load_truth(p.truth_map),
              input_read_spans(p, gap), reads);
    CHECK(reads.spans);
    CHECK(reads.scored >= 3468);
    CHECK(reads.correct >= 3447);
    CHECK(reads.discordant() <= 21);
    CHECK_FALSE(reads.switched);
    // Retain the independently verified Q17 reference call in the connected
    // core; a per-base Q30 cutoff loses parity despite a <5% call error.
    CHECK(reads.primary_scorable == 96);
    CHECK(reads.primary_correct >= 94);
    CHECK(reads.core_correct >= 91);

    std::ifstream vcf(dir + "/native.vcf");
    REQUIRE(vcf.good());
    std::map<std::string, std::pair<std::string, std::string>> calls;
    std::string line;
    while (std::getline(vcf, line)) {
        if (line.empty() || line[0] == '#') continue;
        const auto fields = split_tabs(line);
        if (fields.size() < 10) continue;
        const bool boundary = fields[1] == "61732321" || fields[1] == "61738239" ||
            fields[1] == "61747506" || fields[1] == "61757551";
        if (!boundary) continue;
        const size_t separator = fields[9].rfind(':');
        REQUIRE(separator != std::string::npos);
        const auto call = std::make_pair(fields[9].substr(0, 3),
                                        fields[9].substr(separator + 1));
        calls[fields[1]] = call;
    }
    REQUIRE(calls.size() == 4);
    for (const auto& call : calls) {
        REQUIRE((call.second.first == "0|1" || call.second.first == "1|0"));
        REQUIRE(call.second.second != ".");
        CHECK(std::stoll(call.second.second) > 0);
    }
    CHECK(calls.at("61732321").second == calls.at("61757551").second);
    CHECK(calls.at("61732321").first == calls.at("61757551").first);
    // The source matrix has two ALT/ALT pairs here despite assigning opposite
    // haplotypes. A nominal BAM PS must not import that internal reversal.
    CHECK(calls.at("61738239").second == calls.at("61747506").second);
    CHECK(calls.at("61738239").first == calls.at("61747506").first);
}

static void gap_bam_left_mec_cannot_reverse_the_34_1_mb_graph_block(const Paths& p) {
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

static void gap_bam_recovery_admits_missing_msa_pairs_at_17_62_mb(const Paths& p) {
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

static void gap_bam_block_attaches_to_one_supported_graph_flank(const Paths& p) {
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

static void gap_bam_recovery_uses_a_supported_run_before_a_weak_source_cut(const Paths& p) {
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

static void gap_bam_source_path_with_one_haplotype_molecule_support_closes_34_844_mb(const Paths& p) {
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

static void gap_clean_indel_boundaries_use_the_phased_read_path(const Paths& p) {
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

static void gap_one_weak_cut_bam_run_attaches_to_one_supported_neighbor(const Paths& p) {
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

static void gap_complete_bam_blocks_use_split_stable_allele_votes_without_a_read_gauge(const Paths& p) {
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

static void gap_bam_supported_inner_block_connects_after_an_invalid_graph_snp_is_excluded(const Paths& p) {
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

static void gap_low_mapq_bam_evidence_cannot_veto_a_graph_heterozygote(const Paths& p) {
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

static void gap_shifted_bam_deletion_joins_a_certified_graph_suffix(const Paths& p) {
    // The shifted BAM deletion joins the right graph SNP. Independent
    // primary calls also certify the next graph deletion and the local SNP
    // suffix. Independent source quality and chain certificates now also
    // retain the earlier shared-SNP edge.
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
    CHECK(phase_sets.at("13752640") == phase_sets.at("13773452"));
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

static void gap_a_directly_supported_bam_suffix_closes_the_17_865_mb_gap(const Paths& p) {
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

static void gap_a_newly_imported_bam_seam_closes_the_15_056_mb_gap(const Paths& p) {
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
    // A source MSA SNP can have a wrong local gauge. Its singleton read
    // certificate must preserve the established clean-SNP read orientation.
    CHECK(reads.concordance() >= 0.99);
}

static void gap_recovered_right_deletion_joins_a_source_backed_graph_block(const Paths& p) {
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

static void gap_adjacent_anchors_do_not_trigger_a_second_bam_recovery_solve(const Paths& p) {
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

static void gap_shifted_single_base_bam_insertion_closes_its_owning_23_mb_gap(const Paths& p) {
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

static void gap_a_single_seam_admits_focused_recovery_of_complementary_bam_rows(const Paths& p) {
    Window gap;
    gap.gap_left = 10727690;
    gap.gap_right = 10746628;
    std::string dir;
    REQUIRE(run_arm(p, gap, "stitch_single_seam_10m", "", dir));
    const std::set<std::string> keys{
        "10706319:G>T", "10727690:T>C", "10746473:CCTTT>C",
        "10746473:CCTTTCTTTCTTTCTTT>C", "10746628:C>CTCTT",
        "10760419:G>A"};
    std::map<std::string, std::pair<std::string, std::string>> calls;
    std::ifstream vcf(dir + "/native.vcf");
    REQUIRE(vcf.good());
    std::string line;
    while (std::getline(vcf, line)) {
        if (line.empty() || line[0] == '#') continue;
        const auto fields = split_tabs(line);
        if (fields.size() < 10) continue;
        const std::string key = fields[1] + ":" + fields[3] + ">" + fields[4];
        if (keys.count(key) == 0) continue;
        REQUIRE(fields[8].find("PS") != std::string::npos);
        calls[key] = {fields[9].substr(0, 3),
                      fields[9].substr(fields[9].rfind(':') + 1)};
    }
    REQUIRE(calls.size() == keys.size());
    const std::string ps = calls.at("10727690:T>C").second;
    CHECK(ps != ".");
    for (const auto& [key, call] : calls) {
        INFO(key);
        CHECK(is_phased_het(call.first));
        CHECK(call.second == ps);
    }
    // Both deletion alleles survive separately. The 16-base deletion and
    // insertion share the left SNP's ALT haplotype; the 4-base deletion does
    // not. Pin parity rather than arbitrary numeric HP labels.
    CHECK(calls.at("10727690:T>C").first == calls.at("10746628:C>CTCTT").first);
    CHECK(calls.at("10746473:CCTTTCTTTCTTTCTTT>C").first ==
          calls.at("10746628:C>CTCTT").first);
    CHECK(calls.at("10746473:CCTTT>C").first !=
          calls.at("10746473:CCTTTCTTTCTTTCTTT>C").first);
    CHECK(calls.at("10706319:G>T").first == calls.at("10760419:G>A").first);

    const auto& truth = load_truth(p.truth_map);
    Window graph_seam;
    graph_seam.gap_left = 10706319;
    graph_seam.gap_right = 10760419;
    Outcome flanks;
    score_bam(dir + "/phased.bam", graph_seam, truth,
              input_read_spans(p, graph_seam), flanks);
    CHECK_FALSE(flanks.switched);
    CHECK(flanks.concordance() >= 0.98);

    const auto& spans = input_read_spans(p, gap);
    std::unordered_map<std::string, char> local_truth;
    for (const auto& [name, span] : spans) {
        if (span.second < gap.gap_left || span.first > gap.gap_right) continue;
        const auto found = truth.find(name);
        if (found != truth.end()) local_truth.emplace(*found);
    }
    Outcome local_reads;
    score_bam(dir + "/phased.bam", gap, local_truth, spans, local_reads);
    CHECK(local_reads.scored >= 143);
    CHECK(local_reads.correct >= 136);
    CHECK(local_reads.discordant() <= 7);
    CHECK(local_reads.concordance() >= 0.95);
}

static void gap_an_internal_msa_conflict_does_not_authorize_an_unsupported_block_join(const Paths& p) {
    Window gap;
    gap.gap_left = 19373922;
    gap.gap_right = 19395544;
    std::string dir;
    REQUIRE(run_arm(p, gap, "stitch_internal_conflict_guard_19m", "", dir));
    const std::set<std::string> keys{
        "19373922:T>TTTCC", "19373922:T>TTTCCTTCC", "19395544:T>C"};
    std::map<std::string, std::pair<std::string, std::string>> calls;
    std::ifstream vcf(dir + "/native.vcf");
    REQUIRE(vcf.good());
    std::string line;
    while (std::getline(vcf, line)) {
        if (line.empty() || line[0] == '#') continue;
        const auto fields = split_tabs(line);
        if (fields.size() < 10) continue;
        const std::string key = fields[1] + ":" + fields[3] + ">" + fields[4];
        if (keys.count(key) == 0) continue;
        REQUIRE(fields[8].find("PS") != std::string::npos);
        calls[key] = {fields[9].substr(0, 3),
                      fields[9].substr(fields[9].rfind(':') + 1)};
    }
    REQUIRE(calls.size() == keys.size());
    for (const auto& [key, call] : calls) {
        INFO(key);
        CHECK(is_phased_het(call.first));
        CHECK(call.second != ".");
    }
    const auto& short_insertion = calls.at("19373922:T>TTTCC");
    const auto& long_insertion = calls.at("19373922:T>TTTCCTTCC");
    CHECK(short_insertion.second == long_insertion.second);
    CHECK(short_insertion.first != long_insertion.first);
    // A nominal PS is insufficient: the final source path has a weak cut.
    CHECK(short_insertion.second != calls.at("19395544:T>C").second);
}

static void gap_focused_recovery_retains_a_supported_partial_path_inside_a_graph_seam(const Paths& p) {
    Window gap;
    gap.gap_left = 36332599;
    gap.gap_right = 36354890;
    std::string dir;
    REQUIRE(run_arm(p, gap, "stitch_partial_source_36m", "", dir));
    const std::set<std::string> keys{
        "36332599:G>A", "36343992:C>CA", "36354890:CATAT>C",
        "36620864:G>A", "36623545:A>AT"};
    std::map<std::string, std::pair<std::string, std::string>> calls;
    std::ifstream vcf(dir + "/native.vcf");
    REQUIRE(vcf.good());
    std::string line;
    while (std::getline(vcf, line)) {
        if (line.empty() || line[0] == '#') continue;
        const auto fields = split_tabs(line);
        if (fields.size() < 10) continue;
        const std::string key = fields[1] + ":" + fields[3] + ">" + fields[4];
        if (keys.count(key) == 0) continue;
        if (key != "36620864:G>A") REQUIRE(fields[8].find("PS") != std::string::npos);
        calls[key] = {fields[9].substr(0, 3),
                      fields[8].find("PS") == std::string::npos ? "." :
                      fields[9].substr(fields[9].rfind(':') + 1)};
    }
    REQUIRE(calls.size() == keys.size());
    const auto& left = calls.at("36332599:G>A");
    const auto& insertion = calls.at("36343992:C>CA");
    const auto& deletion = calls.at("36354890:CATAT>C");
    CHECK(left.second != ".");
    for (const auto& [key, call] : calls) {
        INFO(key);
        if (key != "36620864:G>A") CHECK(is_phased_het(call.first));
        if (key.compare(0, 3, "363") == 0)
            CHECK(call.second == left.second);
    }
    // All original BAM bases are ALT; the catalog call survives as a
    // homozygote and cannot split the independently supported indel path.
    CHECK(calls.at("36620864:G>A").first == "1/1");
    CHECK(calls.at("36620864:G>A").second == ".");
    CHECK(calls.at("36623545:A>AT").second == left.second);

    // Pin the chain's relative allele orientation, not arbitrary HP labels.
    CHECK(left.first == deletion.first);
    CHECK(left.first != insertion.first);

    const auto& truth = load_truth(p.truth_map);
    const auto& spans = input_read_spans(p, gap);
    Outcome chunk_reads;
    score_bam(dir + "/phased.bam", gap, truth, spans, chunk_reads);
    CHECK_FALSE(chunk_reads.switched);
    CHECK(chunk_reads.correct >= 3637);
    CHECK(chunk_reads.discordant() <= 153);

    std::unordered_map<std::string, char> local_truth;
    for (const auto& [name, span] : spans) {
        if (span.second < gap.gap_left || span.first > gap.gap_right) continue;
        const auto found = truth.find(name);
        if (found != truth.end()) local_truth.emplace(*found);
    }
    Outcome local_reads;
    score_bam(dir + "/phased.bam", gap, local_truth, spans, local_reads);
    CHECK(local_reads.scored >= 172);
    CHECK(local_reads.correct >= 166);
    CHECK(local_reads.discordant() <= 6);
    CHECK(local_reads.concordance() >= 0.96);
}

static void gap_focused_bam_recovery_retains_insertion_representations_of_graph_child_snps(const Paths& p) {
    Window chunk;
    chunk.gap_left = 59050001;
    chunk.gap_right = 59950000;
    std::string dir;
    REQUIRE(run_arm(p, chunk, "stitch_child_insertions_59m", "", dir));
    std::ifstream vcf(dir + "/native.vcf");
    REQUIRE(vcf.good());
    std::map<std::string, std::pair<std::string, std::string>> calls;
    std::string line;
    while (std::getline(vcf, line)) {
        if (line.empty() || line[0] == '#') continue;
        const auto fields = split_tabs(line);
        if (fields.size() < 10 ||
            (fields[1] != "59679069" && fields[1] != "59757500")) continue;
        const std::string key = fields[1] + ":" + fields[3] + ">" + fields[4];
        REQUIRE(fields[8].find("PS") != std::string::npos);
        calls[key] = {fields[9].substr(0, 3),
                      fields[9].substr(fields[9].rfind(':') + 1)};
    }
    REQUIRE(calls.size() == 2);
    REQUIRE(calls.count("59679069:A>AT") == 1);
    REQUIRE(calls.count("59757500:T>TC") == 1);
    const auto& left = calls.at("59679069:A>AT");
    const auto& right = calls.at("59757500:T>TC");
    CHECK(is_phased_het(left.first));
    CHECK(is_phased_het(right.first));
    CHECK(left.first != right.first);
    CHECK(left.second != ".");
    CHECK(left.second == right.second);

    Window gap;
    gap.gap_left = 59679069;
    gap.gap_right = 59757500;
    Outcome reads;
    score_bam(dir + "/phased.bam", gap, load_truth(p.truth_map),
              input_read_spans(p, gap), reads);
    CHECK_FALSE(reads.switched);
    CHECK(reads.correct >= 3681);
    CHECK(reads.discordant() <= 2);
    CHECK(reads.concordance() >= 0.999);
}

static void gap_chr20_gap_windows_panel_totals(const Paths& p) {
    // The first test case does the emitting; this one has nothing to add to the
    // file and must not assert against expectations that do not exist yet.
    if (std::getenv("PGPHASE_EMIT_EXPECTATIONS") != nullptr) {
        SUCCEED("skipped while emitting expectations");
        return;
    }
    const auto panel = load_panel(p.panel);
    const auto expect = load_expectations(p.expectations);
    const auto& truth = load_truth(p.truth_map);

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

static void gap_observed_bam_insertion_runs_retain_long_alleles_before_the_graph_boundary(const Paths& p) {
    Window gap;
    gap.gap_left = 1194233;
    gap.gap_right = 1196894;
    std::string dir;
    REQUIRE(run_arm(p, gap, "stitch_long_insertion_run_1m", "", dir));
    std::map<long long, std::pair<std::string, std::string>> calls;
    std::ifstream vcf(dir + "/native.vcf");
    REQUIRE(vcf.good());
    std::string line;
    while (std::getline(vcf, line)) {
        if (line.empty() || line[0] == '#') continue;
        const auto fields = split_tabs(line);
        if (fields.size() < 10) continue;
        const long long pos = std::stoll(fields[1]);
        if (pos != gap.gap_left && pos != gap.gap_right &&
            pos != 1196967 && pos != 1197105)
            continue;
        REQUIRE(fields[8].find("PS") != std::string::npos);
        if (pos == gap.gap_left) {
            CHECK(fields[3] == "C");
            CHECK(fields[4] == "T");
        } else if (pos == 1197105) {
            CHECK(fields[3] == "T");
            CHECK(fields[4] == "C");
        } else {
            CHECK(fields[4].size() - fields[3].size() == 155);
            if (pos == gap.gap_right) {
                CHECK(fields[3] == "A");
                CHECK(fields[4] == "ACCAGCCTGGGCAACATGGTGAAACTCTGTCTCTATAAATTAGCTGGATGTGGTGGTGTGACCCTGGAGTCCCAACTACTTGAGAGGCTGAGGTGGGAGGATTGCTTGAGCCCAGGAGGTGGAGGTTGTGGTGAGCCATGATCACACCACTGTACT");
            }
        }
        REQUIRE(calls.emplace(pos, std::make_pair(fields[9].substr(0, 3),
            fields[9].substr(fields[9].rfind(':') + 1))).second);
    }
    REQUIRE(calls.size() == 4);
    const auto& left = calls.at(gap.gap_left);
    const auto& insertion = calls.at(gap.gap_right);
    CHECK(is_phased_het(left.first));
    CHECK(is_phased_het(insertion.first));
    CHECK(left.second != ".");
    // Gauge labels may flip globally; the SNP and insertion ALT must agree.
    CHECK(insertion == left);
    // Both insertion descriptions and the first physical graph SNP must
    // retain the verified allele relation, not merely overlapping PS extents.
    CHECK(calls.at(1196967) == left);
    CHECK(calls.at(1197105) == left);

    const auto& truth = load_truth(p.truth_map);
    const auto& spans = input_read_spans(p, gap);
    Outcome chunk_reads;
    score_bam(dir + "/phased.bam", gap, truth, spans, chunk_reads);
    CHECK_FALSE(chunk_reads.switched);
    CHECK(chunk_reads.discordant() <= 10);
    std::unordered_map<std::string, char> local_truth;
    for (const auto& [name, span] : spans) {
        if (span.second < gap.gap_left || span.first > gap.gap_right) continue;
        const auto found = truth.find(name);
        if (found != truth.end()) local_truth.emplace(*found);
    }
    Outcome local_reads;
    score_bam(dir + "/phased.bam", gap, local_truth, spans, local_reads);
    CHECK(local_reads.scored >= 90);
    CHECK(local_reads.correct >= 90);
    CHECK(local_reads.discordant() == 0);
    CHECK(local_reads.concordance() == 1.0);
}

static void gap_equivalent_bam_and_graph_insertions_close_their_repeat_seam(const Paths& p) {
    Window gap;
    gap.gap_left = 1196894;
    gap.gap_right = 1196967;
    std::string dir;
    REQUIRE(run_arm(p, gap, "stitch_equivalent_insertion_1m", "", dir));
    Outcome variants;
    parse_vcf(dir + "/native.vcf", gap, variants);
    CHECK(variants.spans);
    const auto& truth = load_truth(p.truth_map);
    const auto& spans = input_read_spans(p, gap);
    std::unordered_map<std::string, char> local_truth;
    for (const auto& [name, span] : spans) {
        if (span.second < gap.gap_left || span.first > gap.gap_right) continue;
        const auto found = truth.find(name);
        if (found != truth.end()) local_truth.emplace(*found);
    }
    Outcome local_reads;
    score_bam(dir + "/phased.bam", gap, local_truth, spans, local_reads);
    CHECK(local_reads.scored >= 78);
    CHECK(local_reads.correct >= 78);
    CHECK(local_reads.discordant() == 0);
    CHECK(local_reads.separated() == 1.0);
    CHECK_FALSE(local_reads.switched);
}

static void gap_an_attached_bam_component_keeps_its_source_identity_at_an_insertion_bridge(const Paths& p) {
    Window gap;
    gap.gap_left = 34046350;
    gap.gap_right = 34055920;
    std::string dir;
    REQUIRE(run_arm(p, gap, "stitch_shared_source_34m", "", dir));
    std::map<long long, std::pair<std::string, std::string>> rows;
    std::ifstream vcf(dir + "/native.vcf");
    REQUIRE(vcf.good());
    std::string line;
    while (std::getline(vcf, line)) {
        if (line.empty() || line[0] == '#') continue;
        const auto fields = split_tabs(line);
        if (fields.size() < 10) continue;
        const long long pos = std::stoll(fields[1]);
        if (pos != 34026369 && pos != gap.gap_left &&
            pos != gap.gap_right && pos != 34072055 && pos != 34079290)
            continue;
        const size_t separator = fields[9].rfind(':');
        REQUIRE(separator != std::string::npos);
        if (pos == gap.gap_left) {
            CHECK(fields[3] == "T");
            CHECK(fields[4] == "TTATATATTTATATATATTTATATATATATTTA");
        }
        rows[pos] = {fields[9].substr(0, 3), fields[9].substr(separator + 1)};
    }
    REQUIRE(rows.size() == 5);
    const auto& insertion = rows.at(gap.gap_left);
    CHECK(insertion == rows.at(34026369));
    for (const long long pos : {gap.gap_right, 34072055LL, 34079290LL}) {
        CHECK(rows.at(pos).second == insertion.second);
        CHECK(rows.at(pos).first != insertion.first);
    }
    Outcome reads;
    parse_vcf(dir + "/native.vcf", gap, reads);
    score_bam(dir + "/phased.bam", gap, load_truth(p.truth_map),
              input_read_spans(p, gap), reads);
    CHECK(reads.spans);
    CHECK_FALSE(reads.switched);
    CHECK(reads.scored >= 3567);
    CHECK(reads.discordant() <= 28);
    CHECK(reads.concordance() >= 0.992);
    CHECK(reads.separated() >= 118.0 / 123.0);
}


static void gap_independent_bam_snp_pairs_certify_both_sides_of_a_deletion_seam(const Paths& p) {
    // Both graph blocks contain an agreeing one-haplotype edge. The indel
    // bridge cannot certify those edges; independent BAM SNP pairs must.
    Window gap;
    gap.gap_left = 47751480;
    gap.gap_right = 47762233;
    std::string dir;
    REQUIRE(run_arm(p, gap, "stitch_deletion_physical_paths_47m", "", dir));
    std::ifstream vcf(dir + "/native.vcf");
    REQUIRE(vcf.good());
    std::map<long long, std::pair<std::string, std::string>> rows;
    std::string line;
    while (std::getline(vcf, line)) {
        if (line.empty() || line[0] == '#') continue;
        const auto fields = split_tabs(line);
        if (fields.size() < 10) continue;
        const long long pos = std::stoll(fields[1]);
        if (pos != 47738673 && pos != gap.gap_left &&
            pos != gap.gap_right && pos != 47764156 && pos != 47781656)
            continue;
        const size_t separator = fields[9].rfind(':');
        REQUIRE(separator != std::string::npos);
        if (pos == gap.gap_left) {
            CHECK(fields[3] == "CT");
            CHECK(fields[4] == "C");
        }
        rows[pos] = {fields[9].substr(0, 3), fields[9].substr(separator + 1)};
    }
    REQUIRE(rows.size() == 5);
    const auto& deletion = rows.at(gap.gap_left);
    CHECK(deletion == rows.at(47738673));
    CHECK(deletion == rows.at(gap.gap_right));
    CHECK(deletion == rows.at(47781656));
    CHECK(rows.at(47764156).second == deletion.second);
    CHECK(rows.at(47764156).first != deletion.first);
    Outcome reads;
    parse_vcf(dir + "/native.vcf", gap, reads);
    score_bam(dir + "/phased.bam", gap, load_truth(p.truth_map),
              input_read_spans(p, gap), reads);
    CHECK(reads.spans);
    CHECK_FALSE(reads.switched);
    CHECK(reads.scored >= 3993);
    CHECK(reads.discordant() <= 36);
    CHECK(reads.concordance() >= 0.9909);
    // Continuity must preserve the existing local reads; their repeat-only
    // assignment errors remain a separate issue from this connection.
    CHECK(reads.separated() >= 70.0 / 84.0);
}

static void gap_padded_graph_repeats_retain_the_verified_bam_genotype(const Paths& p) {
    struct Replay {
        long long chunk_beg;
        int scored;
        int correct;
        int errors;
        std::vector<std::string> alleles;
    };
    const std::vector<Replay> replays{
        // An extra repeat call must not join the owning flank blocks with a
        // reversed gauge; a small replay can hide that whole-block switch.
        {2000001, 4088, 4026, 62,
            {"2311483:A>ATTTCTTTCTTTC", "2311483:A>ATTTCTTTCTTTCTTTCTTTC"}},
        {10000001, 4404, 4359, 45, {"10488935:GCTTTTTT>G", "10891999:CA>C"}},
        {35000001, 3362, 3192, 170, {"35490917:TA>T"}},
        {41000001, 4154, 4109, 45, {"41885034:TA>T"}}
    };
    const auto& truth = load_truth(p.truth_map);
    for (const Replay& replay : replays) {
        DYNAMIC_SECTION("owning chunk beginning " << replay.chunk_beg) {
            // Standard padding reproduces the full source gauge. A narrow
            // window can retain the call through a different recovery seam.
            Window w;
            w.gap_left = replay.chunk_beg + 50000;
            w.gap_right = replay.chunk_beg + 949999;
            std::string dir;
            REQUIRE(run_arm(p, w, "padded_repeat_transfer", "", dir));
            std::ifstream vcf(dir + "/native.vcf");
            REQUIRE(vcf.good());
            std::set<std::string> expected(replay.alleles.begin(), replay.alleles.end());
            std::set<std::string> found;
            std::string line;
            while (std::getline(vcf, line)) {
                if (line.empty() || line[0] == '#') continue;
                const auto fields = split_tabs(line);
                if (fields.size() < 10) continue;
                const std::string key = fields[1] + ":" + fields[3] + ">" + fields[4];
                if (expected.count(key) == 0) continue;
                CHECK(found.insert(key).second);
                CHECK(fields[7].find("CAT=NOISY_CAND_HET") != std::string::npos);
                CHECK(is_phased_het(fields[9].substr(0, fields[9].find(':'))));
                std::istringstream format(fields[8]), sample(fields[9]);
                std::string field, value, ps;
                while (std::getline(format, field, ':') && std::getline(sample, value, ':'))
                    if (field == "PS") ps = value;
                REQUIRE_FALSE(ps.empty());
                CHECK(std::stoll(ps) > 0);
            }
            CHECK(found == expected);
            Outcome owning;
            score_bam(dir + "/phased.bam", w, truth, input_read_spans(p, w), owning);
            CHECK(owning.scored >= replay.scored);
            CHECK(owning.correct >= replay.correct);
            CHECK(owning.discordant() <= replay.errors);
        }
    }
}

static void gap_one_haplotype_repeat_deletion_votes_preserve_parental_block_orientation(const Paths& p) {
    // Two Q30 physical ALT/ALT pairs suggest a flip here, but their repeat
    // deletion calls disagree with the independently phased flanks. A future
    // connection is allowed only in the parental orientation, with the owning
    // block context retained so a short replay cannot hide the inversion.
    Window gap;
    gap.gap_left = 4866153;
    gap.gap_right = 4874129;
    std::string dir;
    REQUIRE(run_arm(p, gap, "repeat_deletion_parental_orientation_4m", "", dir));
    std::ifstream vcf(dir + "/native.vcf");
    REQUIRE(vcf.good());
    std::map<long long, std::pair<std::string, std::string>> boundaries;
    std::string line;
    while (std::getline(vcf, line)) {
        if (line.empty() || line[0] == '#') continue;
        const auto fields = split_tabs(line);
        if (fields.size() < 10) continue;
        const long long pos = std::stoll(fields[1]);
        if (pos != gap.gap_left && pos != gap.gap_right) continue;
        CHECK(fields[3] == (pos == gap.gap_left ? "CA" : "ATTT"));
        CHECK(fields[4] == (pos == gap.gap_left ? "C" : "A"));
        const size_t separator = fields[9].rfind(':');
        REQUIRE(separator != std::string::npos);
        boundaries[pos] = {fields[9].substr(0, 3), fields[9].substr(separator + 1)};
    }
    REQUIRE(boundaries.size() == 2);
    const auto& left = boundaries.at(gap.gap_left);
    const auto& right = boundaries.at(gap.gap_right);
    if (left.second == right.second)
        CHECK(left.first != right.first);
    const auto& truth = load_truth(p.truth_map);
    const auto& spans = input_read_spans(p, gap);
    Outcome owning;
    parse_vcf(dir + "/native.vcf", gap, owning);
    score_bam(dir + "/phased.bam", gap, truth, spans, owning);
    CHECK_FALSE(owning.switched);
    // Recovered insertion observations add 26 phased and 16 correct reads;
    // retain that coverage gain and its measured ten-error cost explicitly.


    check_read_floors("gap_one_haplotype_repeat_deletion_votes_preserve_parental_block_orientation", owning);
    std::unordered_map<std::string, char> local_truth;
    for (const auto& entry : spans) {
        if (entry.second.second < gap.gap_left || entry.second.first > gap.gap_right)
            continue;
        const auto found = truth.find(entry.first);
        if (found != truth.end()) local_truth.emplace(*found);
    }
    Outcome local;
    score_bam(dir + "/phased.bam", gap, local_truth, spans, local);
    CHECK(local.scored >= 52);
    CHECK(local.correct >= 51);
    CHECK(local.discordant() <= 1);
}

static void gap_exact_snp_branches_preserve_the_path_through_a_complex_catalog_allele(const Paths& p) {
    // The right graph SNP also occurs in a catalog allele carrying another
    // substitution. Losing that exact branch leaves the internal path one-sided.
    // Replay the owning chunk to retain the complete left and right gauges.
    Window gap;
    gap.gap_left = 20506159;
    gap.gap_right = 20517197;
    std::string dir;
    REQUIRE(run_arm(p, gap, "stitch_catalog_snp_branch_20m", "", dir));
    std::ifstream vcf(dir + "/native.vcf");
    REQUIRE(vcf.good());
    std::map<std::string, std::pair<std::string, std::string>> rows;
    const std::set<std::string> keys{"20445143:A>C", "20506159:C>CT", "20517197:A>AAT",
                                     "20517551:C>T", "20523889:A>T", "20540928:C>T"};
    std::string line;
    while (std::getline(vcf, line)) {
        if (line.empty() || line[0] == '#') continue;
        const auto fields = split_tabs(line);
        REQUIRE(fields.size() >= 10);
        const std::string key = fields[1] + ":" + fields[3] + ">" + fields[4];
        if (keys.count(key) == 0) continue;
        const size_t separator = fields[9].rfind(':');
        REQUIRE(separator != std::string::npos);
        CHECK(rows.emplace(key, std::make_pair(fields[9].substr(0, 3),
                                              fields[9].substr(separator + 1))).second);
    }
    REQUIRE(rows.size() == keys.size());
    const auto& left = rows.at("20506159:C>CT");
    CHECK(left.second != ".");
    for (const auto& [key, row] : rows) {
        INFO(key);
        CHECK(is_phased_het(row.first));
        CHECK(row.second == left.second);
    }
    CHECK(left.first != rows.at("20517197:A>AAT").first);
    CHECK(left.first == rows.at("20517551:C>T").first);
    CHECK(left.first == rows.at("20523889:A>T").first);
    CHECK(left.first != rows.at("20445143:A>C").first);
    CHECK(left.first != rows.at("20540928:C>T").first);
    const auto& truth = load_truth(p.truth_map);
    const auto& spans = input_read_spans(p, gap);
    Outcome owning;
    parse_vcf(dir + "/native.vcf", gap, owning);
    score_bam(dir + "/phased.bam", gap, truth, spans, owning);
    CHECK(owning.spans);
    CHECK_FALSE(owning.switched);



    check_read_floors("gap_exact_snp_branches_preserve_the_path_through_a_complex_catalog_allele", owning);
    std::unordered_map<std::string, char> local_truth;
    for (const auto& [name, span] : spans) {
        if (span.second < gap.gap_left || span.first > gap.gap_right) continue;
        const auto found = truth.find(name);
        if (found != truth.end()) local_truth.emplace(*found);
    }
    Outcome local;
    score_bam(dir + "/phased.bam", gap, local_truth, spans, local);
    CHECK(local.scored >= 110);
    CHECK(local.correct >= 109);
    CHECK(local.discordant() <= 1);
}

static void gap_complementary_insertion_recall_cannot_invert_its_snp_flanks(const Paths& p) {
    // Missing insertion calls can hide a weak cut inside one BAM source PS.
    // A retry must not attach that source to both graph blocks in the wrong
    // parental orientation. Retain the complete owning-chunk flank context.
    Window gap;
    gap.gap_left = 24121713;
    gap.gap_right = 24131707;
    std::string dir;
    REQUIRE(run_arm(p, gap, "msa_insertion_parental_orientation_24m", "", dir));
    const std::set<std::string> keys{"24103779:C>T", "24105188:A>G",
        "24121713:C>CTTTT", "24121713:C>CTTTTTTTT", "24131707:CT>C",
        "24142287:A>G", "24142446:G>A"};
    std::map<std::string, std::pair<std::string, std::string>> rows;
    std::ifstream vcf(dir + "/native.vcf");
    REQUIRE(vcf.good());
    std::string line;
    while (std::getline(vcf, line)) {
        if (line.empty() || line[0] == '#') continue;
        const auto fields = split_tabs(line);
        REQUIRE(fields.size() >= 10);
        const std::string key = fields[1] + ":" + fields[3] + ">" + fields[4];
        if (keys.count(key) == 0) continue;
        const size_t separator = fields[9].rfind(':');
        REQUIRE(separator != std::string::npos);
        CHECK(rows.emplace(key, std::make_pair(fields[9].substr(0, 3),
                                              fields[9].substr(separator + 1))).second);
    }
    REQUIRE(rows.size() == keys.size());
    const auto& left = rows.at("24103779:C>T");
    const auto& right = rows.at("24142287:A>G");
    CHECK(is_phased_het(left.first));
    CHECK(is_phased_het(right.first));
    CHECK(left.second == right.second);
    CHECK(left.first != right.first);
    const auto& four = rows.at("24121713:C>CTTTT");
    const auto& eight = rows.at("24121713:C>CTTTTTTTT");
    if (four.second == eight.second) CHECK(four.first != eight.first);
    const auto& truth = load_truth(p.truth_map);
    const auto& spans = input_read_spans(p, gap);
    Outcome owning;
    parse_vcf(dir + "/native.vcf", gap, owning);
    score_bam(dir + "/phased.bam", gap, truth, spans, owning);

    CHECK(owning.spans);
    CHECK_FALSE(owning.switched);
    CHECK(owning.primary_scorable == 115);
    CHECK(owning.primary_correct >= 93);
    CHECK(owning.core_correct >= 93);
    check_gap_contract(p, gap, owning);


    check_read_floors("gap_complementary_insertion_recall_cannot_invert_its_snp_flanks", owning);
    std::unordered_map<std::string, char> local_truth;
    for (const auto& [name, span] : spans) {
        if (span.second < gap.gap_left || span.first > gap.gap_right) continue;
        const auto found = truth.find(name);
        if (found != truth.end()) local_truth.emplace(*found);
    }
    Outcome local;
    score_bam(dir + "/phased.bam", gap, local_truth, spans, local);
    CHECK(local.scored >= 89);
    CHECK(local.correct >= 93);
    CHECK(local.discordant() == 0);

    std::array<std::array<int, 2>, 2> parents{};
    const std::unique_ptr<samFile, decltype(&hts_close)> bam(
        sam_open((dir + "/phased.bam").c_str(), "r"), hts_close);
    REQUIRE(bam != nullptr);
    const std::unique_ptr<bam_hdr_t, decltype(&bam_hdr_destroy)> header(
        sam_hdr_read(bam.get()), bam_hdr_destroy);
    REQUIRE(header != nullptr);
    const std::unique_ptr<bam1_t, decltype(&bam_destroy1)> record(bam_init1(), bam_destroy1);
    REQUIRE(record != nullptr);
    while (sam_read1(bam.get(), header.get(), record.get()) >= 0) {
        if (record->core.flag & (BAM_FSECONDARY | BAM_FSUPPLEMENTARY)) continue;
        const std::string name = bam_get_qname(record.get());
        const auto parent = truth.find(name);
        const auto span = spans.find(name);
        const uint8_t* hp = bam_aux_get(record.get(), "HP");
        const uint8_t* ps = bam_aux_get(record.get(), "PS");
        if (parent == truth.end() || span == spans.end() || hp == nullptr || ps == nullptr ||
            bam_aux2i(ps) != std::stoll(left.second) || (bam_aux2i(hp) != 1 && bam_aux2i(hp) != 2)) continue;
        const bool mat_on_hap1 = (bam_aux2i(hp) == 1) == (parent->second == 'M');
        if (span->second.second < gap.gap_left && span->second.second > gap.gap_left - 50000)
            ++parents[0][mat_on_hap1];
        if (span->second.first >= gap.gap_right && span->second.first < gap.gap_right + 50000)
            ++parents[1][mat_on_hap1];
    }
    const int orientation = parents[0][1] >= parents[0][0] ? 1 : 0;
    for (const auto& flank : parents) {
        const int total = flank[0] + flank[1];
        REQUIRE(total >= 5);
        CHECK(static_cast<double>(flank[orientation]) / total >= 0.90);
    }
}

static void gap_overlapping_recovery_solves_cannot_exchange_snp_quality_certificates(const Paths& p) {
    // Two source solves disagree. The working BAM slot abstains, the graph's
    // independent REF survives, and each source keeps its own SNP certificate.
    Window region;
    region.gap_left = 3000001;
    region.gap_right = 4000000;
    const std::string prefix = p.workdir + "/retained_snp_quality";
    std::string dir;
    REQUIRE(run_arm(p, region, "retained_snp_quality",
                    "--phase-matrix-dump '" + prefix + "'", dir));
    std::ifstream matrix(prefix + ".chunk0.recovery-final.tsv");
    REQUIRE(matrix.good());
    const std::string qname = "m84031_231217_034919_s2/16389926/ccs";
    int site_index = -1;
    bool observation_seen = false, quality_seen = false;
    std::string line;
    while (std::getline(matrix, line)) {
        const auto fields = split_tabs(line);
        if (fields.size() >= 4 && fields[0] == "VAR" &&
            fields[2] == "3597791" && fields[3] == "X")
            site_index = std::stoi(fields[1]);
        if (fields.size() < 4 || fields[1] != qname || site_index < 0 ||
            std::stoi(fields[2]) != site_index) continue;
        if (fields[0] == "OBS") {
            observation_seen = true;
            CHECK(fields[3] == "0");
            CHECK(fields[4] == "0");
            CHECK(fields[5] == "-2");
        } else if (fields[0] == "BAMQ") {
            quality_seen = true;
            CHECK(fields[3] == "0");
        }
    }
    REQUIRE(site_index >= 0);
    CHECK(observation_seen);
    CHECK(quality_seen);
    std::ifstream evidence(prefix + ".chunk0.bam-source-evidence.tsv");
    REQUIRE(evidence.good());
    std::set<std::pair<int, int>> source_calls;
    while (std::getline(evidence, line)) {
        const auto fields = split_tabs(line);
        if (fields.size() < 9 || fields[1] != "3597791" || fields[5] != qname ||
            fields[8] != "call") continue;
        source_calls.emplace(std::stoi(fields[6]), std::stoi(fields[7]));
    }
    CHECK(source_calls == std::set<std::pair<int, int>>{{0, 40}, {1, 0}});
    // Neither the recovery retry nor the late whole-chunk solve may refill
    // the disputed BAM slot or borrow a quality certificate from the graph.
    for (const char* stage : {"bam-overlay-input", "bam-overlay-output"}) {
        std::ifstream overlay(prefix + ".chunk0." + stage + ".tsv");
        REQUIRE(overlay.good());
        int index = -1;
        bool blocked = false;
        while (std::getline(overlay, line)) {
            const auto fields = split_tabs(line);
            if (fields.size() >= 4 && fields[0] == "VAR" &&
                fields[2] == "3597791" && fields[3] == "X")
                index = std::stoi(fields[1]);
            if (fields.size() >= 6 && fields[0] == "OBS" && fields[1] == qname &&
                index >= 0 && std::stoi(fields[2]) == index) {
                blocked = true;
                CHECK(fields[3] == "0");
                CHECK(fields[5] == "-2");
            }
        }
        CHECK(blocked);
    }
}

static void gap_complementary_insertion_evidence_preserves_the_17_50_mb_connection(const Paths& p) {
    // The left block starts in the preceding chunk. A short replay cannot
    // certify its complete gauge, so keep both owning chunks in this test.
    Window gap;
    gap.gap_left = 17502614;
    gap.gap_right = 17521341;
    const auto& truth = load_truth(p.truth_map);
    const auto expect = load_expectations(p.expectations);
    const Outcome panel = measure(p, gap, "graph", "", truth);
    const auto expected = expect.find("graph\t" + window_key(gap));
    REQUIRE(expected != expect.end());
    check_against(gap, "graph", panel, expected->second);
    std::string dir;
    REQUIRE(run_arm(p, gap, "complementary_insertion_connection_17m", "", dir));
    std::map<long long, std::pair<std::string, std::string>> rows;
    std::ifstream vcf(dir + "/native.vcf");
    REQUIRE(vcf.good());
    std::string line;
    while (std::getline(vcf, line)) {
        if (line.empty() || line[0] == '#') continue;
        const auto fields = split_tabs(line);
        if (fields.size() < 10) continue;
        const long long pos = std::stoll(fields[1]);
        if (pos != gap.gap_left && pos != gap.gap_right) continue;
        CHECK(fields[3] == (pos == gap.gap_left ? "CT" : "A"));
        CHECK(fields[4] == (pos == gap.gap_left ? "C" : "G"));
        const size_t separator = fields[9].rfind(':');
        REQUIRE(separator != std::string::npos);
        rows[pos] = {fields[9].substr(0, 3), fields[9].substr(separator + 1)};
    }
    REQUIRE(rows.size() == 2);
    const auto& left = rows.at(gap.gap_left);
    const auto& right = rows.at(gap.gap_right);
    CHECK(is_phased_het(left.first));
    CHECK(left.second == right.second);
    CHECK(left.first == right.first);

    const auto& spans = input_read_spans(p, gap);
    std::unordered_map<std::string, char> local_truth;
    for (const auto& [name, span] : spans) {
        if (span.second < gap.gap_left || span.first > gap.gap_right) continue;
        const auto found = truth.find(name);
        if (found != truth.end()) local_truth.emplace(*found);
    }
    Outcome local;
    score_bam(dir + "/phased.bam", gap, local_truth, spans, local);
    CHECK(local.scored >= 136);
    CHECK(local.correct >= 131);
    CHECK(local.discordant() <= 5);
    // Individual noisy reads may disagree, but the two flanks must still
    // choose the same parent. Test the relative gauge independently of GT.
    std::array<std::array<int, 2>, 2> votes{};
    const std::unique_ptr<samFile, decltype(&hts_close)> bam(
        sam_open((dir + "/phased.bam").c_str(), "r"), hts_close);
    REQUIRE(bam != nullptr);
    const std::unique_ptr<bam_hdr_t, decltype(&bam_hdr_destroy)> hdr(
        sam_hdr_read(bam.get()), bam_hdr_destroy);
    REQUIRE(hdr != nullptr);
    const std::unique_ptr<bam1_t, decltype(&bam_destroy1)> record(bam_init1(), bam_destroy1);
    while (sam_read1(bam.get(), hdr.get(), record.get()) >= 0) {
        if ((record->core.flag & (BAM_FSECONDARY | BAM_FSUPPLEMENTARY)) != 0)
            continue;
        const std::string name = bam_get_qname(record.get());
        const auto span = spans.find(name);
        const auto parent = truth.find(name);
        const uint8_t* hp = bam_aux_get(record.get(), "HP");
        const uint8_t* ps = bam_aux_get(record.get(), "PS");
        if (span == spans.end() || parent == truth.end() || hp == nullptr || ps == nullptr ||
            std::to_string(bam_aux2i(ps)) != left.second) continue;
        if (bam_aux2i(hp) != 1 && bam_aux2i(hp) != 2) continue;
        const int side = span->second.second < gap.gap_left ? 0 :
                         span->second.first > gap.gap_right ? 1 : -1;
        if (side < 0) continue;
        const bool maternal_on_hap1 = (bam_aux2i(hp) == 1) == (parent->second == 'M');
        ++votes[static_cast<size_t>(side)][maternal_on_hap1];
    }
    REQUIRE(votes[0][0] + votes[0][1] >= 5);
    REQUIRE(votes[1][0] + votes[1][1] >= 5);
    REQUIRE(votes[0][0] != votes[0][1]);
    REQUIRE(votes[1][0] != votes[1][1]);
    const bool left_parent = votes[0][1] > votes[0][0];
    CHECK(left_parent == (votes[1][1] > votes[1][0]));
}

static void gap_complete_bam_evidence_preserves_finalized_37_mb_block_orientation(const Paths& p) {
    // A correct inner-edge vote must not run before source attachment, which
    // can replay a different orientation across the rest of the owning block.
    Window chunk;
    chunk.gap_left = 37050001;
    chunk.gap_right = 37950000;
    std::string dir;
    REQUIRE(run_arm(p, chunk, "complete_source_37m", "", dir));
    Window gap;
    gap.gap_left = 37598255;
    gap.gap_right = 37606443;
    Outcome reads;
    parse_vcf(dir + "/native.vcf", gap, reads);
    score_bam(dir + "/phased.bam", gap, load_truth(p.truth_map),
              input_read_spans(p, gap), reads);
    // A future correct closure is allowed. The full owning chunk must retain
    // its read coverage and parental orientation when these blocks are joined.
    CHECK(reads.scored >= 3607);
    CHECK(reads.correct >= 3561);
    CHECK(reads.discordant() <= 46);
    CHECK_FALSE(reads.switched);
}

static void gap_complementary_deletion_recovery_closes_the_owning_50_mb_gap(const Paths& p) {
    // The short replay lacks the full adjacent source blocks. Its missing
    // allele calls must not veto a connection proven in the owning chunk.
    Window gap;
    gap.gap_left = 50548245;
    gap.gap_right = 50562066;
    std::string dir;
    REQUIRE(run_arm(p, gap, "complementary_deletion_50m", "", dir));
    Outcome reads;
    parse_vcf(dir + "/native.vcf", gap, reads);
    score_bam(dir + "/phased.bam", gap, load_truth(p.truth_map),
              input_read_spans(p, gap), reads);
    CHECK(reads.spans);
    CHECK_FALSE(reads.switched);
    CHECK(reads.scored >= 3887);
    CHECK(reads.correct >= 3876);
    CHECK(reads.discordant() <= 11);
    CHECK(count_hap_allele_conflicts(dir + "/native.vcf") == 0);

    std::map<std::string, std::pair<std::string, std::string>> rows;
    std::ifstream vcf(dir + "/native.vcf");
    REQUIRE(vcf.good());
    std::string line;
    while (std::getline(vcf, line)) {
        if (line.empty() || line[0] == '#') continue;
        const auto fields = split_tabs(line);
        if (fields.size() < 10) continue;
        const long long pos = std::stoll(fields[1]);
        if (pos != gap.gap_left && pos != gap.gap_right) continue;
        const size_t separator = fields[9].rfind(':');
        REQUIRE(separator != std::string::npos);
        rows[fields[3] + ">" + fields[4]] = {
            fields[9].substr(0, 3), fields[9].substr(separator + 1)};
    }
    REQUIRE(rows.count("CA>C") == 1);
    REQUIRE(rows.count("CAA>C") == 1);
    REQUIRE(rows.count("A>G") == 1);
    const auto& shorter = rows.at("CA>C");
    const auto& longer = rows.at("CAA>C");
    const auto& snp = rows.at("A>G");
    CHECK(is_phased_het(shorter.first));
    CHECK(is_phased_het(longer.first));
    CHECK(is_phased_het(snp.first));
    CHECK(shorter.second == longer.second);
    CHECK(longer.second == snp.second);
    // The two deletion ALTs are complementary. The two-base deletion and
    // right SNP ALT must land on the same parent, independently of HP labels.
    CHECK(shorter.first != longer.first);
    CHECK(longer.first == snp.first);
}

static void gap_shifted_compound_cigar_deletions_close_the_owning_7_9_mb_gap(const Paths& p) {
    Window gap;
    gap.gap_left = 7901413;
    gap.gap_right = 7918883;
    std::string dir;
    REQUIRE(run_arm(p, gap, "compound_deletion_parental_orientation_7m", "", dir));
    // Read the exact separate deletion rows without constructing a merged
    // allele. The graph SNP gauges provide the parental-orientation check.
    std::map<int64_t, std::vector<std::pair<std::string, std::string>>> rows;
    std::set<std::string> deletion_keys;
    std::ifstream vcf(dir + "/native.vcf");
    REQUIRE(vcf.good());
    std::string line;
    while (std::getline(vcf, line)) {
        if (line.empty() || line[0] == '#') continue;
        const auto fields = split_tabs(line);
        REQUIRE(fields.size() >= 10);
        const auto pos = std::stoll(fields[1]);
        if (pos != 7901325 && pos != 7901413 && pos != 7918883 &&
            pos != 7924485) continue;
        if (pos == 7918883) deletion_keys.insert(fields[3] + ">" + fields[4]);
        const size_t separator = fields[9].rfind(':');
        REQUIRE(separator != std::string::npos);
        rows[pos].emplace_back(fields[9].substr(0, 3),
                              fields[9].substr(separator + 1));
    }
    REQUIRE(rows[7901325].size() == 1);
    REQUIRE(rows[7901413].size() == 1);
    REQUIRE(rows[7918883].size() == 2);
    CHECK(deletion_keys == std::set<std::string>{
        "ATATTTTTTTTTTT>A", "ATATTTTTTTTTTTTT>A"});
    REQUIRE(rows[7924485].size() == 1);
    const auto& left = rows[7901325].front();
    const auto& right = rows[7924485].front();
    CHECK(is_phased_het(left.first));
    CHECK(left.second == right.second);
    CHECK(left.first == right.first);
    CHECK(rows[7901413].front() == left);
    CHECK(rows[7918883][0].second == left.second);
    CHECK(rows[7918883][1].second == left.second);
    CHECK(rows[7918883][0].first == left.first);
    CHECK(rows[7918883][0].first != rows[7918883][1].first);
    const auto& truth = load_truth(p.truth_map);
    const auto& spans = input_read_spans(p, gap);
    Outcome owning;
    score_bam(dir + "/phased.bam", gap, truth, spans, owning);


    check_read_floors("gap_shifted_compound_cigar_deletions_close_the_owning_7_9_mb_gap", owning);
    std::unordered_map<std::string, char> local_truth;
    for (const auto& [name, span] : spans) {
        if (span.second < gap.gap_left || span.first > gap.gap_right) continue;
        const auto found = truth.find(name);
        if (found != truth.end()) local_truth.emplace(*found);
    }
    Outcome local;
    score_bam(dir + "/phased.bam", gap, local_truth, spans, local);
    CHECK(local.scored >= 116);
    CHECK(local.correct >= 114);
    CHECK(local.discordant() <= 2);
    CHECK(local.dominant_correct >= 104);
}

static void gap_whole_chunk_bam_fallback_preserves_verified_recovery_observations(const Paths& p) {
    // The owning chunk contains separate 19/20-base MSA insertion rows whose
    // calls differ from the later whole-chunk MSA. Existing calls must survive
    // without disabling fallback for other, previously missing observations.
    Window gap;
    gap.gap_left = 61738239;
    gap.gap_right = 61757551;
    const std::string prefix = p.workdir + "/verified_recovery_transfer";
    std::string dir;
    REQUIRE(run_arm(p, gap, "verified_recovery_transfer_61m",
                    "--phase-matrix-dump '" + prefix + "'", dir));
    using CallKey = std::pair<std::string, int>;
    std::map<CallKey, int> before;
    std::set<int> msa_sites;
    std::ifstream input(prefix + ".chunk0.bam-overlay-input.tsv");
    REQUIRE(input.good());
    std::string line;
    while (std::getline(input, line)) {
        const auto fields = split_tabs(line);
        if (fields.size() >= 3 && fields[0] == "#META" && fields[2] == "msa_verified=1")
            msa_sites.insert(std::stoi(fields[1]));
        if (fields.size() >= 6 && fields[0] == "OBS" &&
            (std::stoi(fields[5]) >= 0 || fields[5] == "-2"))
            before.emplace(CallKey{fields[1], std::stoi(fields[2])}, std::stoi(fields[5]));
    }
    std::ifstream output(prefix + ".chunk0.bam-overlay-output.tsv");
    REQUIRE(output.good());
    size_t checked_msa_calls = 0;
    size_t missing_calls_filled = 0;
    size_t changed_calls = 0;
    size_t retained_calls = 0;
    while (std::getline(output, line)) {
        const auto fields = split_tabs(line);
        if (fields.size() < 6 || fields[0] != "OBS") continue;
        const CallKey key{fields[1], std::stoi(fields[2])};
        const int allele = std::stoi(fields[5]);
        const auto old = before.find(key);
        if (old == before.end()) {
            missing_calls_filled += allele >= 0;
            continue;
        }
        ++retained_calls;
        checked_msa_calls += msa_sites.count(key.second);
        changed_calls += allele != old->second;
    }
    CHECK(checked_msa_calls > 0);
    CHECK(missing_calls_filled > 0);
    CHECK(retained_calls == before.size());
    CHECK(changed_calls == 0);

    // The transfer trace resolves BAM sequence keys to graph candidate indices;
    // raw graph-walk strings cannot be compared to the BAM alleles directly.
    std::ifstream transfer(prefix + ".chunk0.transfer.tsv");
    REQUIRE(transfer.good());
    using SourceCallKey = std::tuple<std::string, std::string, std::string,
                                    std::string, std::string, std::string, int>;
    std::set<SourceCallKey> independent_calls;
    std::ifstream evidence(prefix + ".chunk0.bam-source-evidence.tsv");
    REQUIRE(evidence.good());
    while (std::getline(evidence, line)) {
        const auto fields = split_tabs(line);
        if (fields.size() < 9 || fields[8] != "call") continue;
        independent_calls.emplace(fields[0], fields[1], fields[2], fields[3],
                                  fields[4], fields[5], std::stoi(fields[6]));
    }
    REQUIRE(!independent_calls.empty());
    size_t checked_independent_calls = 0, lost_independent_calls = 0;
    std::map<CallKey, std::set<int>> source_calls;
    std::map<CallKey, int> transferred_calls;
    while (std::getline(transfer, line)) {
        const auto fields = split_tabs(line);
        if (fields.size() < 15) continue;
        if (fields[14] == "1") {
            ++checked_independent_calls;
            lost_independent_calls += independent_calls.count(SourceCallKey{
                fields[0], fields[1], fields[2], fields[3], fields[4],
                fields[7], std::stoi(fields[8])}) == 0;
        }
        if (fields[12] != "mapped") continue;
        const CallKey key{fields[7], std::stoi(fields[9])};
        source_calls[key].insert(std::stoi(fields[8]));
        transferred_calls[key] = std::stoi(fields[10]);
    }
    REQUIRE(!source_calls.empty());
    CHECK(checked_independent_calls > 0);
    CHECK(lost_independent_calls == 0);
    size_t lost_source_calls = 0;
    for (const auto& [key, alleles] : source_calls) {
        // Shared conflicts abstain in the working slot. The check above proves
        // that every source allele survives independently, including context
        // sites that do not own an injected working candidate.
        if (alleles.size() == 1 && transferred_calls.at(key) != -2)
            lost_source_calls += transferred_calls.at(key) != *alleles.begin();
    }
    CHECK(lost_source_calls == 0);
}

static void gap_selected_bam_blocks_retain_verified_msa_calls_across_graph_gauge_disagreement(const Paths& p) {
    struct Replay {
        long long gap_left, gap_right, pos;
        const char* ref;
        std::array<std::string, 2> alts;
        int depth, correct, discordant;
    };
    // Replay complete owning chunks. The old admission veto lost six and five
    // verified call pairs respectively even though the selected BAM genotype
    // remains complementary. Added evidence must not weaken read accuracy.
    const std::array<Replay, 2> replays{{
        {11573074, 11586531, 11862429, "T", {"TGATA", "TGATAGATA"}, 68, 4004, 90},
        {34046350, 34055920, 34835166, "C",
         {"CTATATATATATATATATA", "CTATATATATATATATATATATATA"}, 56, 3539, 28}
    }};
    for (const Replay& replay : replays) {
        DYNAMIC_SECTION("verified source at " << replay.pos) {
            Window gap;
            gap.gap_left = replay.gap_left;
            gap.gap_right = replay.gap_right;
            std::string dir;
            REQUIRE(run_arm(p, gap, "verified_graph_gauge_" +
                            std::to_string(replay.pos), "", dir));
            std::ifstream vcf(dir + "/native.vcf");
            REQUIRE(vcf.good());
            std::map<std::string, std::map<std::string, std::string>> rows;
            std::string line;
            while (std::getline(vcf, line)) {
                if (line.empty() || line[0] == '#') continue;
                const auto fields = split_tabs(line);
                if (fields.size() < 10 || std::stoll(fields[1]) != replay.pos ||
                    fields[3] != replay.ref ||
                    std::find(replay.alts.begin(), replay.alts.end(), fields[4]) ==
                        replay.alts.end()) continue;
                REQUIRE(rows.count(fields[4]) == 0);
                std::istringstream format(fields[8]), sample(fields[9]);
                std::string name, value;
                while (std::getline(format, name, ':') && std::getline(sample, value, ':'))
                    rows[fields[4]][name] = value;
            }
            REQUIRE(rows.size() == 2);
            for (const auto& [alt, row] : rows) {
                CHECK(is_phased_het(row.at("GT")));
                CHECK(std::stoi(row.at("DP")) >= replay.depth);
                CHECK(std::stoll(row.at("PS")) > 0);
            }
            CHECK(rows.at(replay.alts[0]).at("GT") != rows.at(replay.alts[1]).at("GT"));
            CHECK(rows.at(replay.alts[0]).at("PS") == rows.at(replay.alts[1]).at("PS"));
            Outcome reads;
            score_bam(dir + "/phased.bam", gap, load_truth(p.truth_map),
                      input_read_spans(p, gap), reads);
            CHECK(reads.correct >= replay.correct);
            CHECK(reads.discordant() <= replay.discordant);
        }
    }
}

static void gap_bam_private_blocks_retain_verified_msa_calls_without_graph_anchors(const Paths& p) {
    Window gap;
    gap.gap_left = 21594343;
    gap.gap_right = 21612458;
    std::string dir;
    REQUIRE(run_arm(p, gap, "verified_private_block_21m", "", dir));
    std::ifstream vcf(dir + "/native.vcf");
    REQUIRE(vcf.good());
    std::map<std::string, std::map<std::string, std::string>> rows;
    std::string line;
    while (std::getline(vcf, line)) {
        if (line.empty() || line[0] == '#') continue;
        const auto fields = split_tabs(line);
        if (fields.size() < 10 || fields[1] != "21594343" || fields[3] != "C" ||
            (fields[4] != "CTCTT" && fields[4] != "CTCTTTCTT")) continue;
        REQUIRE(rows.count(fields[4]) == 0);
        std::istringstream format(fields[8]), sample(fields[9]);
        std::string name, value;
        while (std::getline(format, name, ':') && std::getline(sample, value, ':'))
            rows[fields[4]][name] = value;
    }
    REQUIRE(rows.size() == 2);
    const auto& shorter = rows.at("CTCTT");
    const auto& longer = rows.at("CTCTTTCTT");
    // The independent BAM block has no shared graph SNP. Its two fixed MSA
    // consensuses still classify 62 reads, including 51 deferred call pairs.
    for (const auto& [alt, row] : rows) {
        CHECK(is_phased_het(row.at("GT")));
        CHECK(std::stoll(row.at("DP")) >= 62);
        CHECK(std::stoll(row.at("PS")) > 0);
    }
    CHECK(shorter.at("GT") != longer.at("GT"));
    CHECK(shorter.at("PS") == longer.at("PS"));
    Outcome reads;
    parse_vcf(dir + "/native.vcf", gap, reads);
    score_bam(dir + "/phased.bam", gap, load_truth(p.truth_map),
              input_read_spans(p, gap), reads);
    // Observation admission must not turn an independent gauge into a join.
    CHECK_FALSE(reads.spans);
    CHECK(reads.scored >= 3321);
    CHECK(reads.correct >= 3111);
    CHECK(reads.discordant() <= 210);
}

static void gap_verified_noisy_snp_and_symmetric_graph_path_close_the_41_881_mb_gap(const Paths& p) {
    // The right block includes a repeat SNP whose incoming edge has both
    // haplotypes but whose outgoing edge does not. Keep the complete block:
    // the direct edge around that SNP, and a previously certified BAM join,
    // are both needed to validate the boundary's inherited allele gauge.
    Window gap;
    gap.gap_left = 41880908;
    gap.gap_right = 41885033;
    std::string dir;
    REQUIRE(run_arm(p, gap, "verified_noisy_snp_graph_path_41m", "", dir));
    Outcome reads;
    parse_vcf(dir + "/native.vcf", gap, reads);
    score_bam(dir + "/phased.bam", gap, load_truth(p.truth_map),
              input_read_spans(p, gap), reads);
    CHECK(reads.spans);
    CHECK_FALSE(reads.switched);
    CHECK(reads.scored >= 4156);
    CHECK(reads.correct >= 4129);
    CHECK(reads.discordant() <= 27);
    CHECK(reads.separated() >= 59.0 / 91.0);
    CHECK(count_hap_allele_conflicts(dir + "/native.vcf") == 0);

    std::map<std::string, std::map<std::string, std::string>> rows;
    std::ifstream vcf(dir + "/native.vcf");
    REQUIRE(vcf.good());
    std::string line;
    while (std::getline(vcf, line)) {
        if (line.empty() || line[0] == '#') continue;
        const auto fields = split_tabs(line);
        if (fields.size() < 10) continue;
        const long long pos = std::stoll(fields[1]);
        if (pos != gap.gap_left && pos != gap.gap_right && pos != 41898323)
            continue;
        const std::string key = fields[1] + ":" + fields[3] + ">" + fields[4];
        REQUIRE(rows.count(key) == 0);
        std::istringstream format(fields[8]), sample(fields[9]);
        std::string name, value;
        while (std::getline(format, name, ':') && std::getline(sample, value, ':'))
            rows[key][name] = value;
    }
    REQUIRE(rows.size() == 4);
    const auto& shorter = rows.at("41880908:T>TCACACACACACACA");
    const auto& longer = rows.at("41880908:T>TCACACACACACACACACA");
    const auto& noisy = rows.at("41885033:A>T");
    const auto& clean = rows.at("41898323:C>A");
    for (const auto& [key, row] : rows) {
        CHECK(is_phased_het(row.at("GT")));
        CHECK(row.at("PS") == shorter.at("PS"));
    }
    // Complementary insertion rows keep their exact alleles and depths. The
    // shorter ALT belongs to the same parent as both right SNP ALTs.
    CHECK(shorter.at("GT") != longer.at("GT"));
    CHECK(shorter.at("GT") == noisy.at("GT"));
    CHECK(noisy.at("GT") == clean.at("GT"));
    CHECK(shorter.at("AD") == "22,31");
    CHECK(longer.at("AD") == "32,21");
    CHECK(noisy.at("AD") == "10,10");
}


static void gap_composed_physical_bridges_close_23_461_mb_and_preserve_the_preceding_join(const Paths& p) {
    Window gap;
    gap.gap_left = 23460963;
    gap.gap_right = 23480815;
    std::string dir;
    REQUIRE(run_arm(p, gap, "composed_physical_bridge", "", dir));
    std::ifstream vcf(dir + "/native.vcf");
    REQUIRE(vcf.good());
    std::map<std::string, std::pair<std::string, std::string>> rows;
    std::string line;
    const std::set<std::string> positions{
        "23421003", "23445252", "23449834", "23460963", "23480815", "23503225"};
    while (std::getline(vcf, line)) {
        if (line.empty() || line[0] == '#') continue;
        const auto fields = split_tabs(line);
        if (fields.size() < 10 || positions.count(fields[1]) == 0) continue;
        const size_t separator = fields[9].rfind(':');
        REQUIRE(separator != std::string::npos);
        rows.emplace(fields[1] + ":" + fields[3] + ">" + fields[4],
                     std::make_pair(fields[9].substr(0, 3), fields[9].substr(separator + 1)));
    }
    REQUIRE(rows.size() == 7);
    const std::string phase_set = rows.at("23460963:G>A").second;
    REQUIRE(phase_set != ".");
    for (const auto& [key, row] : rows) {
        INFO(key);
        CHECK(is_phased_het(row.first));
        CHECK(row.second == phase_set);
    }
    // The two insertion descriptions retain their complementary source alleles.
    CHECK(rows.at("23480815:A>AC").first != rows.at("23480815:A>ACC").first);
    CHECK(rows.at("23460963:G>A").first == rows.at("23503225:G>A").first);
    Outcome connection;
    parse_vcf(dir + "/native.vcf", gap, connection);
    CHECK(connection.spans);
    Window preceding;
    preceding.gap_left = 23421003;
    preceding.gap_right = 23445252;
    Outcome old_connection;
    parse_vcf(dir + "/native.vcf", preceding, old_connection);
    CHECK(old_connection.spans);
    const auto& truth = load_truth(p.truth_map);
    Outcome orientation;
    score_bam(dir + "/phased.bam", gap, truth, input_read_spans(p, gap), orientation);
    CHECK_FALSE(orientation.switched);
    CHECK(orientation.concordance() >= 0.98);
    Window owning;
    owning.gap_left = 23000001;
    owning.gap_right = 24000000;
    Outcome reads;
    score_bam(dir + "/phased.bam", owning, truth, input_read_spans(p, owning), reads);
    CHECK(reads.scored >= 3602);
    CHECK(reads.correct >= 3542);
    CHECK(reads.discordant() <= 60);
}

static void gap_deferred_physical_bridge_preserves_read_rescue_across_22_24_mb(const Paths& p) {
    Window gap;
    gap.gap_left = 23460963;
    gap.gap_right = 23480815;
    std::string dir;
    REQUIRE(run_arm(p, gap, "deferred_physical_bridge_paired_chunks",
                    "-r 'CHM13#0#chr20:22000001-24000000'", dir));
    Outcome connection;
    parse_vcf(dir + "/native.vcf", gap, connection);
    CHECK(connection.spans);
    Window owning;
    owning.gap_left = 22000001;
    owning.gap_right = 24000000;
    Outcome reads;
    score_bam(dir + "/phased.bam", owning, load_truth(p.truth_map),
              input_read_spans(p, owning), reads);
    // Changing the marker cohorts before rescue moves five read-only labels
    // between unrelated gauges, turning three correct calls into discordance.
    CHECK(reads.scored >= 7428);
    CHECK(reads.correct >= 7307);
    CHECK(reads.discordant() <= 121);
    // The core bridge does not certify these independent rescue labels.
    const std::map<std::string, int> retained_haps{
        {"m84031_231217_034919_s2/56692951/ccs", 1},
        {"m84031_231217_034919_s2/99026879/ccs", 1},
        {"m84031_231217_062403_s3/74782355/ccs", 2}};
    constexpr long long kRetainedRescuePhaseSet = 1023480815;
    const std::unique_ptr<samFile, decltype(&hts_close)> bam(
        sam_open((dir + "/phased.bam").c_str(), "r"), hts_close);
    REQUIRE(bam != nullptr);
    const std::unique_ptr<bam_hdr_t, decltype(&bam_hdr_destroy)> hdr(
        sam_hdr_read(bam.get()), bam_hdr_destroy);
    REQUIRE(hdr != nullptr);
    const std::unique_ptr<bam1_t, decltype(&bam_destroy1)> record(
        bam_init1(), bam_destroy1);
    std::set<std::string> seen;
    while (sam_read1(bam.get(), hdr.get(), record.get()) >= 0) {
        const auto expected = retained_haps.find(bam_get_qname(record.get()));
        if (expected == retained_haps.end()) continue;
        const uint8_t* hp = bam_aux_get(record.get(), "HP");
        const uint8_t* ps = bam_aux_get(record.get(), "PS");
        REQUIRE(hp != nullptr);
        REQUIRE(ps != nullptr);
        CHECK(bam_aux2i(hp) == expected->second);
        CHECK(bam_aux2i(ps) == kRetainedRescuePhaseSet);
        seen.insert(expected->first);
    }
    CHECK(seen.size() == retained_haps.size());
}

static void gap_complementary_insertion_alts_bridge_the_9_mb_chunk_boundary(const Paths& p) {
    Window gap;
    gap.gap_left = 8977829;
    gap.gap_right = 9014032;
    std::string dir;
    REQUIRE(run_arm(p, gap, "complementary_insertion_boundary", "", dir));
    Outcome connection;
    parse_vcf(dir + "/native.vcf", gap, connection);
    CHECK(connection.spans);
    const auto& truth = load_truth(p.truth_map);
    score_bam(dir + "/phased.bam", gap, truth, input_read_spans(p, gap), connection);
    CHECK_FALSE(connection.switched);
    CHECK(connection.concordance() >= 0.98);
    CHECK(connection.separated() >= 0.56);
    CHECK(connection.window_scorable == 223);
    CHECK(connection.dominant_correct >= 126);
    CHECK(connection.scored >= 8525);
    CHECK(connection.correct >= 8517);
    CHECK(connection.discordant() <= 8);
    std::map<long long, std::pair<std::string, std::string>> anchors;
    std::ifstream vcf(dir + "/native.vcf");
    REQUIRE(vcf.good());
    std::string line;
    while (std::getline(vcf, line)) {
        if (line.empty() || line[0] == '#') continue;
        const auto fields = split_tabs(line);
        if (fields.size() < 10) continue;
        const long long pos = std::stoll(fields[1]);
        if (pos != gap.gap_left && pos != gap.gap_right) continue;
        const size_t separator = fields[9].rfind(':');
        REQUIRE(separator != std::string::npos);
        anchors.emplace(pos, std::make_pair(fields[9].substr(0, 3),
                                            fields[9].substr(separator + 1)));
    }
    REQUIRE(anchors.size() == 2);
    CHECK(is_phased_het(anchors.at(gap.gap_left).first));
    CHECK(is_phased_het(anchors.at(gap.gap_right).first));
    CHECK(anchors.at(gap.gap_left).second == anchors.at(gap.gap_right).second);
    CHECK(anchors.at(gap.gap_left).first != anchors.at(gap.gap_right).first);
}

static void gap_physical_snp_switch_repair_closes_the_owning_62_mb_gap(const Paths& p) {
    Window gap;
    gap.gap_left = 62623253;
    gap.gap_right = 62642316;
    std::string dir;
    REQUIRE(run_arm(p, gap, "physical_switch_62m", "", dir));
    const auto& truth = load_truth(p.truth_map);
    const auto& spans = input_read_spans(p, gap);
    Outcome owning;
    parse_vcf(dir + "/native.vcf", gap, owning);
    score_bam(dir + "/phased.bam", gap, truth, spans, owning);
    CHECK(owning.spans);
    CHECK_FALSE(owning.switched);


    check_read_floors("gap_physical_snp_switch_repair_closes_the_owning_62_mb_gap", owning);
    std::map<long long, std::pair<std::string, std::string>> rows;
    std::ifstream vcf(dir + "/native.vcf");
    REQUIRE(vcf.good());
    std::string line;
    while (std::getline(vcf, line)) {
        if (line.empty() || line[0] == '#') continue;
        const auto fields = split_tabs(line);
        REQUIRE(fields.size() >= 10);
        const long long pos = std::stoll(fields[1]);
        if (pos != 62716326 && pos != 62718395 && pos != 62722021 &&
            pos != 62722416) continue;
        const size_t separator = fields[9].rfind(':');
        REQUIRE(separator != std::string::npos);
        REQUIRE(rows.emplace(pos, std::make_pair(fields[9].substr(0, 3),
            fields[9].substr(separator + 1))).second);
    }
    REQUIRE(rows.size() == 4);
    for (const auto& [pos, row] : rows) {
        INFO(pos);
        CHECK(is_phased_het(row.first));
        CHECK(row.second == rows.at(62718395).second);
    }
    CHECK(rows.at(62716326).first != rows.at(62718395).first);
    CHECK(rows.at(62718395).first != rows.at(62722021).first);
    CHECK(rows.at(62722021).first == rows.at(62722416).first);
}

static void gap_quality_backed_source_path_closes_the_owning_13_mb_gap(const Paths& p) {
    Window gap;
    gap.gap_left = 13752640;
    gap.gap_right = 13773452;
    std::string dir;
    REQUIRE(run_arm(p, gap, "source_quality_path_13m", "", dir));
    Outcome owning;
    parse_vcf(dir + "/native.vcf", gap, owning);
    score_bam(dir + "/phased.bam", gap, load_truth(p.truth_map),
              input_read_spans(p, gap), owning);
    CHECK(owning.spans);
    CHECK_FALSE(owning.switched);





    check_read_floors("gap_quality_backed_source_path_closes_the_owning_13_mb_gap", owning);
    std::map<long long, std::pair<std::string, std::string>> rows;
    std::ifstream vcf(dir + "/native.vcf");
    REQUIRE(vcf.good());
    std::string line;
    while (std::getline(vcf, line)) {
        if (line.empty() || line[0] == '#') continue;
        const auto fields = split_tabs(line);
        REQUIRE(fields.size() >= 10);
        const long long pos = std::stoll(fields[1]);
        if (pos != 13751618 && pos != gap.gap_left && pos != gap.gap_right) continue;
        const size_t separator = fields[9].rfind(':');
        REQUIRE(separator != std::string::npos);
        REQUIRE(rows.emplace(pos, std::make_pair(fields[9].substr(0, 3),
            fields[9].substr(separator + 1))).second);
    }
    REQUIRE(rows.size() == 3);
    for (const auto& [pos, row] : rows) {
        INFO(pos);
        CHECK(is_phased_het(row.first));
        CHECK(row.first == rows.at(gap.gap_left).first);
        CHECK(row.second == rows.at(gap.gap_left).second);
    }
}

static void gap_a_recovered_short_insertion_retains_its_independent_source_path(const Paths& p) {
    Window gap;
    gap.gap_left = 33747591;
    gap.gap_right = 33749688;
    std::string dir;
    REQUIRE(run_arm(p, gap, "short_insertion_source_33m", "", dir));
    Outcome owning;
    parse_vcf(dir + "/native.vcf", gap, owning);
    score_bam(dir + "/phased.bam", gap, load_truth(p.truth_map),
              input_read_spans(p, gap), owning);
    CHECK(owning.spans);
    CHECK_FALSE(owning.switched);




    // HiPhase places 64 of these 67 molecules into one connected block.

    check_read_floors("gap_a_recovered_short_insertion_retains_its_independent_source_path", owning);
    std::map<long long, std::pair<std::string, std::string>> rows;
    std::ifstream vcf(dir + "/native.vcf");
    REQUIRE(vcf.good());
    std::string line;
    while (std::getline(vcf, line)) {
        if (line.empty() || line[0] == '#') continue;
        const auto fields = split_tabs(line);
        REQUIRE(fields.size() >= 10);
        const long long pos = std::stoll(fields[1]);
        if (pos != gap.gap_left && pos != gap.gap_right) continue;
        const size_t separator = fields[9].rfind(':');
        REQUIRE(separator != std::string::npos);
        REQUIRE(rows.emplace(pos, std::make_pair(fields[9].substr(0, 3),
            fields[9].substr(separator + 1))).second);
    }
    REQUIRE(rows.size() == 2);
    for (const auto& [pos, row] : rows) {
        INFO(pos);
        CHECK(is_phased_het(row.first));
        CHECK(row.second == rows.at(gap.gap_left).second);
    }
    CHECK(rows.at(gap.gap_left).first != rows.at(gap.gap_right).first);
    const std::map<std::string, int> core_haps{
        {"m84031_231217_034919_s2/104792409/ccs", 1},
        {"m84031_231217_034919_s2/199365264/ccs", 2},
        {"m84031_231217_034919_s2/58000534/ccs", 2},
        {"m84031_231217_062403_s3/205719379/ccs", 2}};
    const std::set<std::string> conflicting_reads{
        "m84031_231217_034919_s2/181539002/ccs",
        "m84031_231217_034919_s2/261293766/ccs"};
    const std::unique_ptr<samFile, decltype(&hts_close)> bam(
        sam_open((dir + "/phased.bam").c_str(), "r"), hts_close);
    REQUIRE(bam != nullptr);
    const std::unique_ptr<bam_hdr_t, decltype(&bam_hdr_destroy)> header(
        sam_hdr_read(bam.get()), bam_hdr_destroy);
    REQUIRE(header != nullptr);
    const std::unique_ptr<bam1_t, decltype(&bam_destroy1)> record(
        bam_init1(), bam_destroy1);
    REQUIRE(record != nullptr);
    std::set<std::string> seen;
    while (sam_read1(bam.get(), header.get(), record.get()) >= 0) {
        if (record->core.flag & (BAM_FSECONDARY | BAM_FSUPPLEMENTARY)) continue;
        const std::string name = bam_get_qname(record.get());
        const auto expected = core_haps.find(name);
        if (expected == core_haps.end() && conflicting_reads.count(name) == 0)
            continue;
        const uint8_t* hp = bam_aux_get(record.get(), "HP");
        const uint8_t* ps = bam_aux_get(record.get(), "PS");
        if (expected != core_haps.end()) {
            REQUIRE(hp != nullptr);
            REQUIRE(ps != nullptr);
            CHECK(bam_aux2i(hp) == expected->second);
            CHECK(bam_aux2i(ps) == std::stoll(rows.at(gap.gap_left).second));
        } else {
            CHECK((hp == nullptr || bam_aux2i(hp) == 0));
        }
        seen.insert(name);
    }
    CHECK(seen.size() == core_haps.size() + conflicting_reads.size());
}


static void gap_shared_graph_rows_retain_their_independent_bam_run_certificate(const Paths& p) {
    Window gap;
    gap.gap_left = 49720311;
    gap.gap_right = 49742973;
    std::string dir;
    REQUIRE(run_arm(p, gap, "shared_source_snp_49m", "", dir));
    Outcome owning;
    parse_vcf(dir + "/native.vcf", gap, owning);
    score_bam(dir + "/phased.bam", gap, load_truth(p.truth_map),
              input_read_spans(p, gap), owning);
    CHECK(owning.spans);
    CHECK_FALSE(owning.switched);




    // Native HiPhase puts the same 137 correct molecules in one core block.
    check_read_floors("gap_shared_graph_rows_retain_their_independent_bam_run_certificate", owning);
    std::map<long long, std::pair<std::string, std::string>> rows;
    std::ifstream vcf(dir + "/native.vcf");
    REQUIRE(vcf.good());
    std::string line;
    while (std::getline(vcf, line)) {
        if (line.empty() || line[0] == '#') continue;
        const auto fields = split_tabs(line);
        REQUIRE(fields.size() >= 10);
        const long long pos = std::stoll(fields[1]);
        if (pos != gap.gap_left && pos != gap.gap_right) continue;
        const size_t separator = fields[9].rfind(':');
        REQUIRE(separator != std::string::npos);
        REQUIRE(rows.emplace(pos, std::make_pair(fields[9].substr(0, 3),
            fields[9].substr(separator + 1))).second);
    }
    REQUIRE(rows.size() == 2);
    CHECK(is_phased_het(rows.at(gap.gap_left).first));
    CHECK(rows.at(gap.gap_left) == rows.at(gap.gap_right));
}

static void gap_a_clean_snp_retry_certifies_a_complete_independent_bam_suffix(const Paths& p) {
    Window gap;
    gap.gap_left = 57764235;
    gap.gap_right = 57785224;
    std::string dir;
    REQUIRE(run_arm(p, gap, "clean_snp_source_retry_57m", "", dir));
    Outcome owning;
    parse_vcf(dir + "/native.vcf", gap, owning);
    score_bam(dir + "/phased.bam", gap, load_truth(p.truth_map),
              input_read_spans(p, gap), owning);
    CHECK(owning.spans);
    CHECK_FALSE(owning.switched);




    // Native HiPhase puts the same 128 correct molecules in one core block.
    check_read_floors("gap_a_clean_snp_retry_certifies_a_complete_independent_bam_suffix", owning);
    std::map<long long, std::pair<std::string, std::string>> rows;
    std::ifstream vcf(dir + "/native.vcf");
    REQUIRE(vcf.good());
    std::string line;
    while (std::getline(vcf, line)) {
        if (line.empty() || line[0] == '#') continue;
        const auto fields = split_tabs(line);
        REQUIRE(fields.size() >= 10);
        const long long pos = std::stoll(fields[1]);
        if (pos != gap.gap_left && pos != 57785772) continue;
        const size_t separator = fields[9].rfind(':');
        REQUIRE(separator != std::string::npos);
        REQUIRE(rows.emplace(pos, std::make_pair(fields[9].substr(0, 3),
            fields[9].substr(separator + 1))).second);
    }
    REQUIRE(rows.size() == 2);
    CHECK(is_phased_het(rows.at(gap.gap_left).first));
    CHECK(is_phased_het(rows.at(57785772).first));
    CHECK(rows.at(gap.gap_left).first != rows.at(57785772).first);
    CHECK(rows.at(gap.gap_left).second == rows.at(57785772).second);
}

TEST_CASE("gap certification counts abstentions and excludes rescue cores", "[gap][unit]") {
    const auto directory = std::filesystem::temp_directory_path() /
        "pgphase-gap-contract-unit";
    std::filesystem::create_directories(directory);
    const auto path = directory / "reads.sam";
    {
        std::ofstream sam(path);
        REQUIRE(sam.good());
        sam << "@HD\tVN:1.6\n@SQ\tSN:chr20\tLN:1000\n"
            << "maternal\t0\tchr20\t100\t60\t10M\t*\t0\t0\tAAAAAAAAAA\t*\tHP:i:1\tPS:i:100\n"
            << "paternal\t0\tchr20\t100\t60\t10M\t*\t0\t0\tAAAAAAAAAA\t*\tHP:i:2\tPS:i:100\n"
            << "rescued\t0\tchr20\t100\t60\t10M\t*\t0\t0\tAAAAAAAAAA\t*\tHP:i:1\tPS:i:1000000001\n"
            << "unphased\t0\tchr20\t100\t60\t10M\t*\t0\t0\tAAAAAAAAAA\t*\n";
    }
    const std::unordered_map<std::string, char> truth{
        {"maternal", 'M'}, {"paternal", 'P'}, {"rescued", 'M'}, {"unphased", 'M'}};
    const ReadSpans spans{{"maternal", {99, 109}}, {"paternal", {99, 109}},
        {"rescued", {99, 109}}, {"unphased", {99, 109}}};
    Window gap;
    gap.gap_left = 100;
    gap.gap_right = 108;
    Outcome got;
    score_bam(path.string(), gap, truth, spans, got);
    CHECK(got.primary_scorable == 4);
    CHECK(got.primary_correct == 3);
    CHECK(got.core_correct == 2);
    std::filesystem::remove(path);
}

static void gap_retained_shared_deletion_closes_the_62_408_mb_gap(const Paths& p) {
    Window gap;
    gap.gap_left = 62408056;
    gap.gap_right = 62432427;
    const Outcome got = measure(p, gap, "graph", "", load_truth(p.truth_map));
    CHECK(got.spans);
    CHECK_FALSE(got.switched);
    CHECK(got.primary_scorable == 177);
    CHECK(got.primary_correct >= 157);
    CHECK(got.core_correct >= 157);
    check_gap_contract(p, gap, got);
    std::string dir;
    REQUIRE(run_arm(p, gap, "graph", "", dir));
    std::map<long long, std::pair<std::string, std::string>> calls;
    std::ifstream vcf(dir + "/native.vcf");
    REQUIRE(vcf.good());
    std::string line;
    while (std::getline(vcf, line)) {
        if (line.empty() || line[0] == '#') continue;
        const auto fields = split_tabs(line);
        if (fields.size() < 10) continue;
        const long long pos = std::stoll(fields[1]);
        if (pos != gap.gap_left && pos != gap.gap_right) continue;
        calls[pos] = {fields[9].substr(0, 3), fields[9].substr(fields[9].rfind(':') + 1)};
    }
    REQUIRE(calls.size() == 2);
    REQUIRE(is_phased_het(calls.at(gap.gap_left).first));
    CHECK(calls.at(gap.gap_left) == calls.at(gap.gap_right));
    const long long core = std::stoll(calls.at(gap.gap_left).second);
    const int alt_hap = calls.at(gap.gap_left).first == "1|0" ? 1 : 2;
    const std::set<std::string> retained{
        "m84031_231217_034919_s2/69075309/ccs",
        "m84031_231217_062403_s3/58004984/ccs",
        "m84031_231217_034919_s2/36898742/ccs",
        "m84031_231217_062403_s3/80809387/ccs"};
    const std::unique_ptr<samFile, decltype(&hts_close)> bam(
        sam_open((dir + "/phased.bam").c_str(), "r"), hts_close);
    REQUIRE(bam != nullptr);
    const std::unique_ptr<bam_hdr_t, decltype(&bam_hdr_destroy)> header(
        sam_hdr_read(bam.get()), bam_hdr_destroy);
    const std::unique_ptr<bam1_t, decltype(&bam_destroy1)> record(bam_init1(), bam_destroy1);
    REQUIRE(header != nullptr);
    REQUIRE(record != nullptr);
    std::set<std::string> seen;
    while (sam_read1(bam.get(), header.get(), record.get()) >= 0) {
        const std::string name = bam_get_qname(record.get());
        if (!retained.count(name)) continue;
        const uint8_t* hp = bam_aux_get(record.get(), "HP");
        const uint8_t* ps = bam_aux_get(record.get(), "PS");
        REQUIRE(hp != nullptr);
        REQUIRE(ps != nullptr);
        CHECK(bam_aux2i(hp) == alt_hap);
        CHECK(bam_aux2i(ps) == core);
        seen.insert(name);
    }
    CHECK(seen == retained);
}

static void gap_cut_free_source_paths_close_the_17_634_mb_gap(const Paths& p) {
    Window gap;
    gap.gap_left = 17634393;
    gap.gap_right = 17667022;
    const auto& truth = load_truth(p.truth_map);
    const Outcome got = measure(p, gap, "graph", "", truth);
    check_read_floors("gap_cut_free_source_paths_close_the_17_634_mb_gap", got);
    CHECK(got.spans);
    CHECK_FALSE(got.switched);
    CHECK(got.primary_scorable == 202);
    CHECK(got.primary_correct >= 164);
    CHECK(got.core_correct >= 164);
    check_gap_contract(p, gap, got);
    std::string dir;
    REQUIRE(run_arm(p, gap, "graph", "", dir));
    std::map<std::string, std::pair<std::string, std::string>> calls;
    std::ifstream vcf(dir + "/native.vcf");
    REQUIRE(vcf.good());
    const std::set<long long> positions{
        17487837, 17521341, 17633256, 17634393, 17667022, 17723058, 17744061};
    std::string line;
    while (std::getline(vcf, line)) {
        if (line.empty() || line[0] == '#') continue;
        const auto fields = split_tabs(line);
        if (fields.size() < 10 || !positions.count(std::stoll(fields[1]))) continue;
        const std::string key = fields[1] + ":" + fields[3] + ">" + fields[4];
        REQUIRE(calls.emplace(key, std::make_pair(fields[9].substr(0, 3),
            fields[9].substr(fields[9].rfind(':') + 1))).second);
    }
    REQUIRE(calls.size() == 8);
    const auto& left_snp = calls.at("17633256:G>C");
    const long long core = std::stoll(left_snp.second);
    for (const auto& [key, call] : calls) {
        INFO(key);
        CHECK(is_phased_het(call.first));
        CHECK(call.second == left_snp.second);
    }
    CHECK(calls.at("17487837:C>G").first == left_snp.first);
    CHECK(calls.at("17521341:A>G").first != left_snp.first);
    CHECK(calls.at("17667022:G>T").first == left_snp.first);
    CHECK(calls.at("17634393:CT>C").first == left_snp.first);
    CHECK(calls.at("17634393:CTTTTTTTTTTTT>C").first != left_snp.first);
    CHECK(calls.at("17723058:A>T").first == calls.at("17744061:G>A").first);

    // Disjoint primary flanks must independently establish the same parent;
    // a pooled major orientation alone could conceal a reversed join.
    const auto& spans = input_read_spans(p, gap);
    std::array<std::array<int, 2>, 2> parents{};
    const std::unique_ptr<samFile, decltype(&hts_close)> bam(
        sam_open((dir + "/phased.bam").c_str(), "r"), hts_close);
    REQUIRE(bam != nullptr);
    const std::unique_ptr<bam_hdr_t, decltype(&bam_hdr_destroy)> header(
        sam_hdr_read(bam.get()), bam_hdr_destroy);
    const std::unique_ptr<bam1_t, decltype(&bam_destroy1)> record(bam_init1(), bam_destroy1);
    REQUIRE(header != nullptr);
    REQUIRE(record != nullptr);
    while (sam_read1(bam.get(), header.get(), record.get()) >= 0) {
        if (record->core.flag & (BAM_FSECONDARY | BAM_FSUPPLEMENTARY)) continue;
        const std::string name = bam_get_qname(record.get());
        const auto parent = truth.find(name);
        const auto span = spans.find(name);
        const uint8_t* hp = bam_aux_get(record.get(), "HP");
        const uint8_t* ps = bam_aux_get(record.get(), "PS");
        if (parent == truth.end() || span == spans.end() || hp == nullptr || ps == nullptr ||
            bam_aux2i(ps) != core || (bam_aux2i(hp) != 1 && bam_aux2i(hp) != 2)) continue;
        const bool mat_on_hap1 = (bam_aux2i(hp) == 1) == (parent->second == 'M');
        if (span->second.second < gap.gap_left) ++parents[0][mat_on_hap1];
        if (span->second.first >= gap.gap_right) ++parents[1][mat_on_hap1];
    }
    const int orientation = parents[0][1] >= parents[0][0] ? 1 : 0;
    for (const auto& flank : parents) {
        const int total = flank[0] + flank[1];
        REQUIRE(total >= 5);
        CHECK(static_cast<double>(flank[orientation]) / total >= 0.90);
    }
}

static void gap_calibrated_repeat_snps_join_the_largest_hiphase_block(const Paths& p) {
    Window gap;
    gap.gap_left = 65509355;
    gap.gap_right = 65509406;
    const auto& truth = load_truth(p.truth_map);
    const Outcome got = measure(p, gap, "graph", "", truth);
    check_read_floors("gap_calibrated_repeat_snps_join_the_largest_hiphase_block", got);
    CHECK(got.spans);
    CHECK_FALSE(got.switched);
    CHECK(got.primary_scorable == 6);
    CHECK(got.primary_correct >= 6);
    CHECK(got.core_correct >= 6);
    check_gap_contract(p, gap, got);
    std::string dir;
    REQUIRE(run_arm(p, gap, "graph", "", dir));
    std::map<std::string, std::pair<std::string, std::string>> calls;
    std::ifstream vcf(dir + "/native.vcf");
    REQUIRE(vcf.good());
    const std::set<long long> positions{
        65000495, 65505704, 65509355, 65509406, 65515063, 65996492};
    std::string line;
    while (std::getline(vcf, line)) {
        if (line.empty() || line[0] == '#') continue;
        const auto fields = split_tabs(line);
        if (fields.size() < 10 || !positions.count(std::stoll(fields[1]))) continue;
        const std::string key = fields[1] + ":" + fields[3] + ">" + fields[4];
        REQUIRE(calls.emplace(key, std::make_pair(fields[9].substr(0, 3),
            fields[9].substr(fields[9].rfind(':') + 1))).second);
    }
    REQUIRE(calls.size() == positions.size());
    const auto& left_snp = calls.at("65505704:C>G");
    const long long core = std::stoll(left_snp.second);
    for (const auto& [key, call] : calls) {
        INFO(key);
        CHECK(is_phased_het(call.first));
        CHECK(call.second == left_snp.second);
    }
    CHECK(calls.at("65515063:C>A").first == left_snp.first);

    // Disjoint primary flanks must independently establish the same parent;
    // a pooled major orientation alone could conceal a reversed join.
    const auto& spans = input_read_spans(p, gap);
    std::array<std::array<int, 2>, 2> parents{};
    const std::unique_ptr<samFile, decltype(&hts_close)> bam(
        sam_open((dir + "/phased.bam").c_str(), "r"), hts_close);
    REQUIRE(bam != nullptr);
    const std::unique_ptr<bam_hdr_t, decltype(&bam_hdr_destroy)> header(
        sam_hdr_read(bam.get()), bam_hdr_destroy);
    const std::unique_ptr<bam1_t, decltype(&bam_destroy1)> record(bam_init1(), bam_destroy1);
    REQUIRE(header != nullptr);
    REQUIRE(record != nullptr);
    while (sam_read1(bam.get(), header.get(), record.get()) >= 0) {
        if (record->core.flag & (BAM_FSECONDARY | BAM_FSUPPLEMENTARY)) continue;
        const std::string name = bam_get_qname(record.get());
        const auto parent = truth.find(name);
        const auto span = spans.find(name);
        const uint8_t* hp = bam_aux_get(record.get(), "HP");
        const uint8_t* ps = bam_aux_get(record.get(), "PS");
        if (parent == truth.end() || span == spans.end() || hp == nullptr || ps == nullptr ||
            bam_aux2i(ps) != core || (bam_aux2i(hp) != 1 && bam_aux2i(hp) != 2)) continue;
        const bool mat_on_hap1 = (bam_aux2i(hp) == 1) == (parent->second == 'M');
        if (span->second.second < gap.gap_left) ++parents[0][mat_on_hap1];
        if (span->second.first >= gap.gap_right) ++parents[1][mat_on_hap1];
    }
    const int orientation = parents[0][1] >= parents[0][0] ? 1 : 0;
    for (const auto& flank : parents) {
        const int total = flank[0] + flank[1];
        REQUIRE(total >= 5);
        CHECK(static_cast<double>(flank[orientation]) / total >= 0.90);
    }
}

static void gap_calibrated_source_deletion_joins_the_fourth_largest_hiphase_block(const Paths& p) {
    Window gap;
    gap.gap_left = 52696940;
    gap.gap_right = 52711825;
    const auto& truth = load_truth(p.truth_map);
    const Outcome got = measure(p, gap, "graph", "", truth);
    CHECK(got.spans);
    CHECK_FALSE(got.switched);
    CHECK(got.primary_scorable == 136);
    CHECK(got.primary_correct >= 119);
    CHECK(got.core_correct >= 119);
    check_gap_contract(p, gap, got);
    check_read_floors("gap_calibrated_source_deletion_joins_the_fourth_largest_hiphase_block", got);
    std::string dir;
    REQUIRE(run_arm(p, gap, "graph", "", dir));
    const std::set<std::string> keys{"52696410:GGA>G", "52696940:G>A",
        "52711825:G>GA", "52715881:TTGTG>T", "52727684:A>C"};
    std::map<std::string, std::pair<std::string, std::string>> rows;
    std::ifstream vcf(dir + "/native.vcf");
    REQUIRE(vcf.good());
    std::string line;
    while (std::getline(vcf, line)) {
        if (line.empty() || line[0] == '#') continue;
        const auto fields = split_tabs(line);
        REQUIRE(fields.size() >= 10);
        const std::string key = fields[1] + ":" + fields[3] + ">" + fields[4];
        if (!keys.count(key)) continue;
        REQUIRE(rows.emplace(key, std::make_pair(fields[9].substr(0, 3),
            fields[9].substr(fields[9].rfind(':') + 1))).second);
    }
    REQUIRE(rows.size() == keys.size());
    const auto& left = rows.at("52696940:G>A");
    const long long core = std::stoll(left.second);
    for (const auto& [key, row] : rows) {
        INFO(key);
        CHECK(is_phased_het(row.first));
        CHECK(row.second == left.second);
        CHECK(row.first == (key == "52711825:G>GA" ?
            std::string(1, left.first[2]) + "|" + left.first[0] : left.first));
    }
    // Only two reads end in the disjoint left flank. Instead use disjoint
    // marker-bearing cohorts: left-SNP reads ending before the deletion and
    // downstream SNP reads starting beyond the left SNP. The bridge molecule
    // belongs to neither cohort, so a pooled majority cannot hide inversion.
    const auto& spans = input_read_spans(p, gap);
    std::array<std::array<int, 2>, 2> parents{};
    const std::unique_ptr<samFile, decltype(&hts_close)> bam(
        sam_open((dir + "/phased.bam").c_str(), "r"), hts_close);
    REQUIRE(bam != nullptr);
    const std::unique_ptr<bam_hdr_t, decltype(&bam_hdr_destroy)> header(
        sam_hdr_read(bam.get()), bam_hdr_destroy);
    const std::unique_ptr<bam1_t, decltype(&bam_destroy1)> record(bam_init1(), bam_destroy1);
    REQUIRE(header != nullptr);
    REQUIRE(record != nullptr);
    bool bridge_seen = false;
    while (sam_read1(bam.get(), header.get(), record.get()) >= 0) {
        const std::string name = bam_get_qname(record.get());
        const auto parent = truth.find(name);
        const auto span = spans.find(name);
        const uint8_t* hp = bam_aux_get(record.get(), "HP");
        const uint8_t* ps = bam_aux_get(record.get(), "PS");
        if (parent == truth.end() || span == spans.end() || !hp || !ps ||
            bam_aux2i(ps) != core || (bam_aux2i(hp) != 1 && bam_aux2i(hp) != 2)) continue;
        const bool mat_on_hap1 = (bam_aux2i(hp) == 1) == (parent->second == 'M');
        if (span->second.first < 52696940 && span->second.second >= 52696940 &&
            span->second.second < 52715882) ++parents[0][mat_on_hap1];
        if (span->second.first >= 52696940 && span->second.second >= 52727684)
            ++parents[1][mat_on_hap1];
        if (name == "m84031_231217_062403_s3/163778931/ccs") {
            bridge_seen = true;
            CHECK(bam_aux2i(hp) == (left.first[0] == '1' ? 1 : 2));
            CHECK(parent->second == 'P');
        }
    }
    REQUIRE(bridge_seen);
    const int orientation = parents[0][1] >= parents[0][0] ? 1 : 0;
    for (const auto& cohort : parents) {
        const int total = cohort[0] + cohort[1];
        REQUIRE(total >= 20);
        CHECK(static_cast<double>(cohort[orientation]) / total >= 0.90);
    }
}

static void gap_physically_contradicted_graph_snps_close_the_25_855_mb_seam(const Paths& p) {
    // Copy-specific graph SNPs on the same molecules are contradicted by the
    // original CIGAR and sequence. Keep the native owner: a short replay loses
    // the left phase set and cannot reproduce the false terminal anchor.
    Window gap;
    gap.gap_left = 25855631;
    gap.gap_right = 25855633;
    const auto& truth = load_truth(p.truth_map);
    const Outcome got = measure(p, gap, "graph", "", truth);
    CHECK(got.spans);
    CHECK_FALSE(got.switched);
    CHECK(got.primary_scorable == 91);
    CHECK(got.primary_correct >= 88);
    CHECK(got.core_correct >= 88);
    check_gap_contract(p, gap, got);
    check_read_floors("gap_physically_contradicted_graph_snps_close_the_25_855_mb_seam", got);
    std::string dir;
    REQUIRE(run_arm(p, gap, "graph", "", dir));
    const std::set<std::string> keys{"25842903:G>C", "25856023:T>G",
        "25859488:T>C", "25891621:G>A"};
    std::map<std::string, std::pair<std::string, std::string>> rows;
    std::ifstream vcf(dir + "/native.vcf");
    REQUIRE(vcf.good());
    std::string line;
    while (std::getline(vcf, line)) {
        if (line.empty() || line[0] == '#') continue;
        const auto fields = split_tabs(line);
        REQUIRE(fields.size() >= 10);
        const std::string key = fields[1] + ":" + fields[3] + ">" + fields[4];
        if (key == "25853217:G>A" || key == "25855631:T>C")
            CHECK_FALSE(is_phased_het(fields[9].substr(0, 3)));
        if (keys.count(key) == 0) continue;
        CHECK(rows.emplace(key, std::make_pair(fields[9].substr(0, 3),
            fields[9].substr(fields[9].rfind(':') + 1))).second);
    }
    REQUIRE(rows.size() == keys.size());
    const auto& left = rows.at("25842903:G>C");
    REQUIRE(is_phased_het(left.first));
    const long long core = std::stoll(left.second);
    for (const auto& [key, row] : rows) {
        INFO(key);
        CHECK(is_phased_het(row.first));
        CHECK(row.second == left.second);
        CHECK(row.first == left.first);
    }
    const auto& spans = input_read_spans(p, gap);
    std::array<std::array<int, 2>, 2> parents{};
    const std::unique_ptr<samFile, decltype(&hts_close)> bam(
        sam_open((dir + "/phased.bam").c_str(), "r"), hts_close);
    REQUIRE(bam != nullptr);
    const std::unique_ptr<bam_hdr_t, decltype(&bam_hdr_destroy)> header(
        sam_hdr_read(bam.get()), bam_hdr_destroy);
    REQUIRE(header != nullptr);
    const std::unique_ptr<bam1_t, decltype(&bam_destroy1)> record(bam_init1(), bam_destroy1);
    REQUIRE(record != nullptr);
    while (sam_read1(bam.get(), header.get(), record.get()) >= 0) {
        if (record->core.flag & (BAM_FSECONDARY | BAM_FSUPPLEMENTARY)) continue;
        const std::string name = bam_get_qname(record.get());
        const auto parent = truth.find(name);
        const auto span = spans.find(name);
        const uint8_t* hp = bam_aux_get(record.get(), "HP");
        const uint8_t* ps = bam_aux_get(record.get(), "PS");
        if (parent == truth.end() || span == spans.end() || hp == nullptr || ps == nullptr ||
            bam_aux2i(ps) != core || (bam_aux2i(hp) != 1 && bam_aux2i(hp) != 2)) continue;
        const bool mat_on_hap1 = (bam_aux2i(hp) == 1) == (parent->second == 'M');
        if (span->second.second < gap.gap_left && span->second.second > gap.gap_left - 50000)
            ++parents[0][mat_on_hap1];
        if (span->second.first >= gap.gap_right && span->second.first < gap.gap_right + 50000)
            ++parents[1][mat_on_hap1];
    }
    const int orientation = parents[0][1] >= parents[0][0] ? 1 : 0;
    for (const auto& flank : parents) {
        const int total = flank[0] + flank[1];
        REQUIRE(total >= 5);
        CHECK(static_cast<double>(flank[orientation]) / total >= 0.90);
    }
}

static void gap_complementary_insertions_join_the_second_largest_hiphase_block(const Paths& p) {
    Window gap;
    gap.gap_left = 10325039;
    gap.gap_right = 10337661;
    const auto& truth = load_truth(p.truth_map);
    const Outcome got = measure(p, gap, "graph", "", truth);
    check_read_floors("gap_complementary_insertions_join_the_second_largest_hiphase_block", got);
    CHECK(got.spans);
    CHECK_FALSE(got.switched);
    CHECK(got.primary_scorable == 127);
    CHECK(got.primary_correct >= 108);
    CHECK(got.core_correct >= 98);
    check_gap_contract(p, gap, got);
    std::string dir;
    REQUIRE(run_arm(p, gap, "graph", "", dir));
    std::map<std::string, std::pair<std::string, std::string>> calls;
    std::ifstream vcf(dir + "/native.vcf");
    REQUIRE(vcf.good());
    const std::set<long long> positions{
        10002997, 10318335, 10325039, 10337661, 10343595, 10999362};
    std::string line;
    while (std::getline(vcf, line)) {
        if (line.empty() || line[0] == '#') continue;
        const auto fields = split_tabs(line);
        if (fields.size() < 10 || !positions.count(std::stoll(fields[1]))) continue;
        const std::string key = fields[1] + ":" + fields[3] + ">" + fields[4];
        REQUIRE(calls.emplace(key, std::make_pair(fields[9].substr(0, 3),
            fields[9].substr(fields[9].rfind(':') + 1))).second);
    }
    REQUIRE(calls.size() == positions.size() + 1);
    const auto& left_snp = calls.at("10318335:C>A");
    const long long core = std::stoll(left_snp.second);
    for (const auto& [key, call] : calls) {
        INFO(key);
        CHECK(is_phased_het(call.first));
        CHECK(call.second == left_snp.second);
    }
    CHECK(calls.at("10325039:T>TGGAAGGAA").first == left_snp.first);
    CHECK(calls.at("10337661:G>GA").first == left_snp.first);
    CHECK(calls.at("10325039:T>TGGAA").first != left_snp.first);
    CHECK(calls.at("10343595:A>G").first != left_snp.first);

    // Disjoint primary flanks must independently establish the same parent;
    // a pooled major orientation alone could conceal a reversed join.
    const auto& spans = input_read_spans(p, gap);
    std::array<std::array<int, 2>, 2> parents{};
    const std::unique_ptr<samFile, decltype(&hts_close)> bam(
        sam_open((dir + "/phased.bam").c_str(), "r"), hts_close);
    REQUIRE(bam != nullptr);
    const std::unique_ptr<bam_hdr_t, decltype(&bam_hdr_destroy)> header(
        sam_hdr_read(bam.get()), bam_hdr_destroy);
    const std::unique_ptr<bam1_t, decltype(&bam_destroy1)> record(bam_init1(), bam_destroy1);
    REQUIRE(header != nullptr);
    REQUIRE(record != nullptr);
    while (sam_read1(bam.get(), header.get(), record.get()) >= 0) {
        if (record->core.flag & (BAM_FSECONDARY | BAM_FSUPPLEMENTARY)) continue;
        const std::string name = bam_get_qname(record.get());
        const auto parent = truth.find(name);
        const auto span = spans.find(name);
        const uint8_t* hp = bam_aux_get(record.get(), "HP");
        const uint8_t* ps = bam_aux_get(record.get(), "PS");
        if (parent == truth.end() || span == spans.end() || hp == nullptr || ps == nullptr ||
            bam_aux2i(ps) != core || (bam_aux2i(hp) != 1 && bam_aux2i(hp) != 2)) continue;
        const bool mat_on_hap1 = (bam_aux2i(hp) == 1) == (parent->second == 'M');
        if (span->second.second < gap.gap_left) ++parents[0][mat_on_hap1];
        if (span->second.first >= gap.gap_right) ++parents[1][mat_on_hap1];
    }
    const int orientation = parents[0][1] >= parents[0][0] ? 1 : 0;
    for (const auto& flank : parents) {
        const int total = flank[0] + flank[1];
        REQUIRE(total >= 5);
        CHECK(static_cast<double>(flank[orientation]) / total >= 0.90);
    }
}

static void gap_homozygous_graph_snp_does_not_split_a_supported_block(const Paths& p) {
    Window gap;
    gap.gap_left = 36614185;
    gap.gap_right = 36623545;
    const auto& truth = load_truth(p.truth_map);
    const Outcome got = measure(p, gap, "graph", "", truth);
    CHECK(got.spans);
    CHECK_FALSE(got.switched);
    CHECK(got.primary_scorable == 102);
    CHECK(got.primary_correct >= 93);
    CHECK(got.core_correct >= 92);
    check_gap_contract(p, gap, got);
    std::string dir;
    REQUIRE(run_arm(p, gap, "graph", "", dir));
    const std::set<std::string> keys{"36597548:A>G", "36614185:A>AT",
        "36620864:G>A", "36623545:A>AT", "36629321:C>T"};
    std::map<std::string, std::pair<std::string, std::string>> calls;
    std::ifstream vcf(dir + "/native.vcf");
    REQUIRE(vcf.good());
    std::string line;
    while (std::getline(vcf, line)) {
        if (line.empty() || line[0] == '#') continue;
        const auto fields = split_tabs(line);
        const std::string key = fields[1] + ":" + fields[3] + ">" + fields[4];
        if (!keys.count(key)) continue;
        calls.emplace(key, std::make_pair(fields[9].substr(0, 3),
            fields[8].find("PS") == std::string::npos ? "." :
            fields[9].substr(fields[9].rfind(':') + 1)));
        if (key == "36620864:G>A") {
            CHECK(fields[9].find(":55:11,44:") != std::string::npos);
            CHECK(fields[7].find("CAT=CLEAN_HOM") != std::string::npos);
        }
    }
    REQUIRE(calls.size() == keys.size());
    const auto& left = calls.at("36597548:A>G");
    CHECK(left.second != ".");
    CHECK(calls.at("36620864:G>A").first == "1/1");
    CHECK(calls.at("36620864:G>A").second == ".");
    for (const std::string key : {"36614185:A>AT", "36623545:A>AT", "36629321:C>T"}) {
        const auto& call = calls.at(key);
        CHECK(is_phased_het(call.first));
        CHECK(call.second == left.second);
        CHECK(call.first == (key == "36629321:C>T" ? left.first :
            std::string(1, left.first[2]) + "|" + left.first[0]));
    }
}

static void gap_cut_free_insertion_prefix_closes_the_36_268_mb_gap(const Paths& p) {
    Window gap;
    gap.gap_left = 36268558;
    gap.gap_right = 36286778;
    const auto& truth = load_truth(p.truth_map);
    const Outcome got = measure(p, gap, "graph", "", truth);
    CHECK(got.spans);
    CHECK_FALSE(got.switched);
    CHECK(got.primary_scorable == 109);
    CHECK(got.primary_correct >= 100);
    CHECK(got.core_correct >= 100);
    check_gap_contract(p, gap, got);
    check_read_floors("gap_cut_free_insertion_prefix_closes_the_36_268_mb_gap", got);
    std::string dir;
    REQUIRE(run_arm(p, gap, "graph", "", dir));
    const std::set<long long> positions{36268558, 36286778, 36299817, 36317511};
    std::map<long long, std::pair<std::string, long long>> calls;
    std::ifstream vcf(dir + "/native.vcf");
    REQUIRE(vcf.good());
    std::string line;
    while (std::getline(vcf, line)) {
        if (line.empty() || line[0] == '#') continue;
        const auto fields = split_tabs(line);
        const long long pos = std::stoll(fields[1]);
        if (!positions.count(pos)) continue;
        REQUIRE(fields[8].find("PS") != std::string::npos);
        REQUIRE(calls.emplace(pos, std::make_pair(fields[9].substr(0, 3),
            std::stoll(fields[9].substr(fields[9].rfind(':') + 1)))).second);
    }
    REQUIRE(calls.size() == positions.size());
    const auto& left = calls.at(gap.gap_left);
    for (const auto& [pos, call] : calls) {
        INFO(pos);
        CHECK(is_phased_het(call.first));
        CHECK(call.second == left.second);
        CHECK(call.first == (pos == 36317511 ?
            std::string(1, left.first[2]) + "|" + left.first[0] : left.first));
    }
    const auto& spans = input_read_spans(p, gap);
    std::array<std::array<int, 2>, 2> parents{};
    const std::unique_ptr<samFile, decltype(&hts_close)> bam(
        sam_open((dir + "/phased.bam").c_str(), "r"), hts_close);
    REQUIRE(bam != nullptr);
    const std::unique_ptr<bam_hdr_t, decltype(&bam_hdr_destroy)> header(
        sam_hdr_read(bam.get()), bam_hdr_destroy);
    const std::unique_ptr<bam1_t, decltype(&bam_destroy1)> record(bam_init1(), bam_destroy1);
    REQUIRE(header != nullptr);
    REQUIRE(record != nullptr);
    while (sam_read1(bam.get(), header.get(), record.get()) >= 0) {
        if (record->core.flag & (BAM_FSECONDARY | BAM_FSUPPLEMENTARY)) continue;
        const std::string name = bam_get_qname(record.get());
        const auto parent = truth.find(name);
        const auto span = spans.find(name);
        const uint8_t* hp = bam_aux_get(record.get(), "HP");
        const uint8_t* ps = bam_aux_get(record.get(), "PS");
        if (parent == truth.end() || span == spans.end() || hp == nullptr || ps == nullptr ||
            bam_aux2i(ps) != left.second || (bam_aux2i(hp) != 1 && bam_aux2i(hp) != 2)) continue;
        const bool mat_on_hap1 = (bam_aux2i(hp) == 1) == (parent->second == 'M');
        if (span->second.second < gap.gap_left) ++parents[0][mat_on_hap1];
        if (span->second.first >= gap.gap_right) ++parents[1][mat_on_hap1];
    }
    const int orientation = parents[0][1] >= parents[0][0] ? 1 : 0;
    for (const auto& flank : parents) {
        const int total = flank[0] + flank[1];
        REQUIRE(total >= 5);
        CHECK(static_cast<double>(flank[orientation]) / total >= 0.90);
    }
}

static void gap_ref_absent_graph_snp_recovers_the_40_633_mb_terminal_deletion(const Paths& p) {
    Window gap;
    gap.gap_left = 40633644;
    gap.gap_right = 40636353;
    const auto& truth = load_truth(p.truth_map);
    const Outcome got = measure(p, gap, "graph", "", truth);
    CHECK(got.spans);
    CHECK_FALSE(got.switched);
    CHECK(got.primary_scorable == 45);
    CHECK(got.primary_correct >= 45);
    CHECK(got.core_correct >= 41);
    check_gap_contract(p, gap, got);
    check_read_floors("gap_ref_absent_graph_snp_recovers_the_40_633_mb_terminal_deletion", got);
    std::string dir;
    REQUIRE(run_arm(p, gap, "graph", "", dir));
    std::map<std::string, std::pair<std::string, long long>> calls;
    std::ifstream vcf(dir + "/native.vcf");
    REQUIRE(vcf.good());
    std::string line;
    while (std::getline(vcf, line)) {
        if (line.empty() || line[0] == '#') continue;
        const auto fields = split_tabs(line);
        const std::string key = fields[1] + ":" + fields[3] + ">" + fields[4];
        if (key != "40633644:G>A" && key != "40636353:TG>T") continue;
        REQUIRE(calls.emplace(key, std::make_pair(fields[9].substr(0, 3),
            std::stoll(fields[9].substr(fields[9].rfind(':') + 1)))).second);
    }
    REQUIRE(calls.size() == 2);
    const auto& left = calls.at("40633644:G>A");
    const auto& deletion = calls.at("40636353:TG>T");
    CHECK(is_phased_het(left.first));
    CHECK(is_phased_het(deletion.first));
    CHECK(deletion.second == left.second);
    CHECK(deletion.first == std::string(1, left.first[2]) + "|" + left.first[0]);
    const auto& spans = input_read_spans(p, gap);
    std::array<std::array<int, 2>, 2> parents{};
    const std::unique_ptr<samFile, decltype(&hts_close)> bam(
        sam_open((dir + "/phased.bam").c_str(), "r"), hts_close);
    REQUIRE(bam != nullptr);
    const std::unique_ptr<bam_hdr_t, decltype(&bam_hdr_destroy)> header(
        sam_hdr_read(bam.get()), bam_hdr_destroy);
    const std::unique_ptr<bam1_t, decltype(&bam_destroy1)> record(bam_init1(), bam_destroy1);
    REQUIRE(header != nullptr);
    REQUIRE(record != nullptr);
    while (sam_read1(bam.get(), header.get(), record.get()) >= 0) {
        if (record->core.flag & (BAM_FSECONDARY | BAM_FSUPPLEMENTARY)) continue;
        const std::string name = bam_get_qname(record.get());
        const auto parent = truth.find(name);
        const auto span = spans.find(name);
        const uint8_t* hp = bam_aux_get(record.get(), "HP");
        const uint8_t* ps = bam_aux_get(record.get(), "PS");
        if (parent == truth.end() || span == spans.end() || hp == nullptr || ps == nullptr ||
            bam_aux2i(ps) != left.second || (bam_aux2i(hp) != 1 && bam_aux2i(hp) != 2)) continue;
        const bool mat_on_hap1 = (bam_aux2i(hp) == 1) == (parent->second == 'M');
        if (span->second.second < gap.gap_left) ++parents[0][mat_on_hap1];
        // These reads start beyond the upstream marker and observe only the
        // terminal allele; they are disjoint from the upstream cohort.
        if (span->second.first > gap.gap_left && span->second.second >= gap.gap_right)
            ++parents[1][mat_on_hap1];
    }
    const int orientation = parents[0][1] >= parents[0][0] ? 1 : 0;
    for (const auto& flank : parents) {
        const int total = flank[0] + flank[1];
        REQUIRE(total >= 3);
        CHECK(static_cast<double>(flank[orientation]) / total >= 0.90);
    }
}

static void gap_physically_validated_repeat_snp_closes_the_11_796_mb_gap(const Paths& p) {
    Window gap;
    gap.gap_left = 11796979;
    gap.gap_right = 11813446;
    const auto& truth = load_truth(p.truth_map);
    const Outcome got = measure(p, gap, "graph", "", truth);
    CHECK(got.spans);
    CHECK_FALSE(got.switched);
    CHECK(got.primary_scorable == 136);
    CHECK(got.primary_correct >= 126);
    CHECK(got.core_correct >= 126);
    check_gap_contract(p, gap, got);
    check_read_floors("gap_physically_validated_repeat_snp_closes_the_11_796_mb_gap", got);
    std::string dir;
    REQUIRE(run_arm(p, gap, "graph", "", dir));
    std::map<std::string, std::pair<std::string, long long>> calls;
    std::ifstream vcf(dir + "/native.vcf");
    REQUIRE(vcf.good());
    std::string line;
    while (std::getline(vcf, line)) {
        if (line.empty() || line[0] == '#') continue;
        const auto fields = split_tabs(line);
        const std::string key = fields[1] + ":" + fields[3] + ">" + fields[4];
        CHECK(key != "11813622:T>C");
        if (key != "11796969:GTGTGTGTGTGTGTA>G" && key != "11796979:GTGTA>G" &&
            key != "11813446:AT>A" && key != "11813622:T>TAC" && key != "11813668:T>C") continue;
        REQUIRE(calls.emplace(key, std::make_pair(fields[9].substr(0, 3),
            std::stoll(fields[9].substr(fields[9].rfind(':') + 1)))).second);
    }
    REQUIRE(calls.size() == 5);
    const auto& left = calls.at("11796979:GTGTA>G");
    for (const auto& [key, call] : calls) {
        CHECK(is_phased_het(call.first));
        CHECK(call.second == left.second);
        CHECK(call.first == (key == "11796979:GTGTA>G" || key == "11813668:T>C" ? left.first :
            std::string(1, left.first[2]) + "|" + left.first[0]));
    }
    const auto& spans = input_read_spans(p, gap);
    std::array<std::array<int, 2>, 2> parents{};
    const std::unique_ptr<samFile, decltype(&hts_close)> bam(
        sam_open((dir + "/phased.bam").c_str(), "r"), hts_close);
    REQUIRE(bam != nullptr);
    const std::unique_ptr<bam_hdr_t, decltype(&bam_hdr_destroy)> header(
        sam_hdr_read(bam.get()), bam_hdr_destroy);
    const std::unique_ptr<bam1_t, decltype(&bam_destroy1)> record(bam_init1(), bam_destroy1);
    REQUIRE(header != nullptr);
    REQUIRE(record != nullptr);
    while (sam_read1(bam.get(), header.get(), record.get()) >= 0) {
        if (record->core.flag & (BAM_FSECONDARY | BAM_FSUPPLEMENTARY)) continue;
        const std::string name = bam_get_qname(record.get());
        const auto parent = truth.find(name);
        const auto span = spans.find(name);
        const uint8_t* hp = bam_aux_get(record.get(), "HP");
        const uint8_t* ps = bam_aux_get(record.get(), "PS");
        if (parent == truth.end() || span == spans.end() || hp == nullptr || ps == nullptr ||
            bam_aux2i(ps) != left.second || (bam_aux2i(hp) != 1 && bam_aux2i(hp) != 2)) continue;
        const bool mat_on_hap1 = (bam_aux2i(hp) == 1) == (parent->second == 'M');
        if (span->second.first <= gap.gap_left && span->second.second >= gap.gap_left &&
            span->second.second < gap.gap_right) ++parents[0][mat_on_hap1];
        if (span->second.first > gap.gap_left && span->second.second >= gap.gap_right)
            ++parents[1][mat_on_hap1];
    }
    const int orientation = parents[0][1] >= parents[0][0] ? 1 : 0;
    for (const auto& flank : parents) {
        const int total = flank[0] + flank[1];
        REQUIRE(total >= 3);
        CHECK(static_cast<double>(flank[orientation]) / total >= 0.90);
    }
    // A late local union must also retain the already stitched continuation
    // in the next owning chunk.
    Window continuation;
    continuation.gap_left = 11050001;
    continuation.gap_right = 12950000;
    std::string continuation_dir;
    REQUIRE(run_arm(p, continuation, "graph", "", continuation_dir));
    std::ifstream continuation_vcf(continuation_dir + "/native.vcf");
    REQUIRE(continuation_vcf.good());
    std::map<std::string, long long> continuation_ps;
    while (std::getline(continuation_vcf, line)) {
        if (line.empty() || line[0] == '#') continue;
        const auto fields = split_tabs(line);
        const std::string key = fields[1] + ":" + fields[3] + ">" + fields[4];
        if (key != "11796979:GTGTA>G" && key != "11813446:AT>A" && key != "12717796:C>T") continue;
        CHECK(is_phased_het(fields[9].substr(0, 3)));
        continuation_ps.emplace(key, std::stoll(fields[9].substr(fields[9].rfind(':') + 1)));
    }
    REQUIRE(continuation_ps.size() == 3);
    CHECK(continuation_ps.at("11796979:GTGTA>G") == continuation_ps.at("11813446:AT>A"));
    CHECK(continuation_ps.at("11796979:GTGTA>G") == continuation_ps.at("12717796:C>T"));
}

static void gap_graph_repeat_indel_chain_closes_the_12_717_mb_gap(const Paths& p) {
    Window gap;
    gap.gap_left = 12717796;
    gap.gap_right = 12740002;
    const auto& truth = load_truth(p.truth_map);
    const Outcome got = measure(p, gap, "graph", "", truth);
    CHECK(got.spans);
    CHECK_FALSE(got.switched);
    CHECK(got.primary_scorable == 178);
    CHECK(got.primary_correct >= 164);
    CHECK(got.core_correct >= 164);
    check_gap_contract(p, gap, got);
    check_read_floors("gap_graph_repeat_indel_chain_closes_the_12_717_mb_gap", got);
    std::string dir;
    REQUIRE(run_arm(p, gap, "graph", "", dir));
    std::map<std::string, std::pair<std::string, long long>> calls;
    std::ifstream vcf(dir + "/native.vcf");
    REQUIRE(vcf.good());
    std::string line;
    while (std::getline(vcf, line)) {
        if (line.empty() || line[0] == '#') continue;
        const auto fields = split_tabs(line);
        const std::string key = fields[1] + ":" + fields[3] + ">" + fields[4];
        if (key != "12717796:C>T" && key != "12735894:TA>T" && key != "12752291:C>T") continue;
        REQUIRE(calls.emplace(key, std::make_pair(fields[9].substr(0, 3),
            std::stoll(fields[9].substr(fields[9].rfind(':') + 1)))).second);
    }
    REQUIRE(calls.size() == 3);
    const auto& left = calls.at("12717796:C>T");
    for (const auto& [key, call] : calls) {
        CHECK(is_phased_het(call.first));
        CHECK(call.second == left.second);
        CHECK(call.first == (key == "12752291:C>T" ?
            std::string(1, left.first[2]) + "|" + left.first[0] : left.first));
    }
    const auto& spans = input_read_spans(p, gap);
    std::array<std::array<int, 2>, 2> parents{};
    const std::unique_ptr<samFile, decltype(&hts_close)> bam(
        sam_open((dir + "/phased.bam").c_str(), "r"), hts_close);
    REQUIRE(bam != nullptr);
    const std::unique_ptr<bam_hdr_t, decltype(&bam_hdr_destroy)> header(
        sam_hdr_read(bam.get()), bam_hdr_destroy);
    const std::unique_ptr<bam1_t, decltype(&bam_destroy1)> record(bam_init1(), bam_destroy1);
    REQUIRE(header != nullptr);
    REQUIRE(record != nullptr);
    while (sam_read1(bam.get(), header.get(), record.get()) >= 0) {
        if (record->core.flag & (BAM_FSECONDARY | BAM_FSUPPLEMENTARY)) continue;
        const std::string name = bam_get_qname(record.get());
        const auto parent = truth.find(name);
        const auto span = spans.find(name);
        const uint8_t* hp = bam_aux_get(record.get(), "HP");
        const uint8_t* ps = bam_aux_get(record.get(), "PS");
        if (parent == truth.end() || span == spans.end() || hp == nullptr || ps == nullptr ||
            bam_aux2i(ps) != left.second || (bam_aux2i(hp) != 1 && bam_aux2i(hp) != 2)) continue;
        const bool mat_on_hap1 = (bam_aux2i(hp) == 1) == (parent->second == 'M');
        if (span->second.first <= gap.gap_left && span->second.second >= gap.gap_left &&
            span->second.second < gap.gap_right) ++parents[0][mat_on_hap1];
        if (span->second.first > gap.gap_left && span->second.second >= gap.gap_right)
            ++parents[1][mat_on_hap1];
    }
    const int orientation = parents[0][1] >= parents[0][0] ? 1 : 0;
    for (const auto& flank : parents) {
        const int total = flank[0] + flank[1];
        REQUIRE(total >= 3);
        CHECK(static_cast<double>(flank[orientation]) / total >= 0.90);
    }
    // A late local union must also retain the already stitched continuation
    // in the next owning chunk.
    Window continuation;
    continuation.gap_left = 11050001;
    continuation.gap_right = 12950000;
    std::string continuation_dir;
    REQUIRE(run_arm(p, continuation, "graph", "", continuation_dir));
    std::ifstream continuation_vcf(continuation_dir + "/native.vcf");
    REQUIRE(continuation_vcf.good());
    std::map<std::string, long long> continuation_ps;
    while (std::getline(continuation_vcf, line)) {
        if (line.empty() || line[0] == '#') continue;
        const auto fields = split_tabs(line);
        const std::string key = fields[1] + ":" + fields[3] + ">" + fields[4];
        if (key != "11796979:GTGTA>G" && key != "12717796:C>T" && key != "12752291:C>T") continue;
        CHECK(is_phased_het(fields[9].substr(0, 3)));
        continuation_ps.emplace(key, std::stoll(fields[9].substr(fields[9].rfind(':') + 1)));
    }
    REQUIRE(continuation_ps.size() == 3);
    CHECK(continuation_ps.at("11796979:GTGTA>G") == continuation_ps.at("12717796:C>T"));
    CHECK(continuation_ps.at("12717796:C>T") == continuation_ps.at("12752291:C>T"));
}

static void gap_physical_terminal_insertion_closes_the_12_954_mb_gap(const Paths& p) {
    Window gap;
    gap.gap_left = 12954878;
    gap.gap_right = 12955838;
    const auto& truth = load_truth(p.truth_map);
    const Outcome got = measure(p, gap, "graph", "", truth);
    CHECK(got.spans);
    CHECK_FALSE(got.switched);
    CHECK(got.primary_scorable == 84);
    CHECK(got.primary_correct >= 84);
    CHECK(got.core_correct >= 84);
    check_gap_contract(p, gap, got);
    check_read_floors("gap_physical_terminal_insertion_closes_the_12_954_mb_gap", got);
    std::string dir;
    REQUIRE(run_arm(p, gap, "graph", "", dir));
    std::map<std::string, std::pair<std::string, long long>> calls;
    std::ifstream vcf(dir + "/native.vcf");
    REQUIRE(vcf.good());
    std::string line;
    while (std::getline(vcf, line)) {
        if (line.empty() || line[0] == '#') continue;
        const auto fields = split_tabs(line);
        const std::string key = fields[1] + ":" + fields[3] + ">" + fields[4];
        if (key != "12954878:A>G" && key != "12955838:T>TC" && key != "12954322:G>A") continue;
        REQUIRE(calls.emplace(key, std::make_pair(fields[9].substr(0, 3),
            std::stoll(fields[9].substr(fields[9].rfind(':') + 1)))).second);
    }
    REQUIRE(calls.size() == 3);
    const auto& left = calls.at("12954878:A>G");
    for (const auto& [key, call] : calls) {
        CHECK(is_phased_het(call.first));
        CHECK(call.second == left.second);
        CHECK(call.first == left.first);
    }
    const auto& spans = input_read_spans(p, gap);
    std::array<std::array<int, 2>, 2> parents{};
    const std::unique_ptr<samFile, decltype(&hts_close)> bam(
        sam_open((dir + "/phased.bam").c_str(), "r"), hts_close);
    REQUIRE(bam != nullptr);
    const std::unique_ptr<bam_hdr_t, decltype(&bam_hdr_destroy)> header(
        sam_hdr_read(bam.get()), bam_hdr_destroy);
    const std::unique_ptr<bam1_t, decltype(&bam_destroy1)> record(bam_init1(), bam_destroy1);
    REQUIRE(header != nullptr);
    REQUIRE(record != nullptr);
    while (sam_read1(bam.get(), header.get(), record.get()) >= 0) {
        if (record->core.flag & (BAM_FSECONDARY | BAM_FSUPPLEMENTARY)) continue;
        const std::string name = bam_get_qname(record.get());
        const auto parent = truth.find(name);
        const auto span = spans.find(name);
        const uint8_t* hp = bam_aux_get(record.get(), "HP");
        const uint8_t* ps = bam_aux_get(record.get(), "PS");
        if (parent == truth.end() || span == spans.end() || hp == nullptr || ps == nullptr ||
            bam_aux2i(ps) != left.second || (bam_aux2i(hp) != 1 && bam_aux2i(hp) != 2)) continue;
        const bool mat_on_hap1 = (bam_aux2i(hp) == 1) == (parent->second == 'M');
        if (span->second.first <= gap.gap_left && span->second.second >= gap.gap_left &&
            span->second.second < gap.gap_right) ++parents[0][mat_on_hap1];
        if (span->second.first > gap.gap_left && span->second.second >= gap.gap_right)
            ++parents[1][mat_on_hap1];
    }
    const int orientation = parents[0][1] >= parents[0][0] ? 1 : 0;
    for (const auto& flank : parents) {
        const int total = flank[0] + flank[1];
        REQUIRE(total >= 3);
        CHECK(static_cast<double>(flank[orientation]) / total >= 0.90);
    }
    // The terminal insertion must retain the final stitched SNP gauge
    // across the preceding owning chunk.
    Window continuation;
    continuation.gap_left = 11050001;
    continuation.gap_right = 12955838;
    std::string continuation_dir;
    REQUIRE(run_arm(p, continuation, "graph", "", continuation_dir));
    std::ifstream continuation_vcf(continuation_dir + "/native.vcf");
    REQUIRE(continuation_vcf.good());
    std::map<std::string, long long> continuation_ps;
    while (std::getline(continuation_vcf, line)) {
        if (line.empty() || line[0] == '#') continue;
        const auto fields = split_tabs(line);
        const std::string key = fields[1] + ":" + fields[3] + ">" + fields[4];
        if (key != "11796979:GTGTA>G" && key != "12735894:TA>T" && key != "12955838:T>TC") continue;
        CHECK(is_phased_het(fields[9].substr(0, 3)));
        continuation_ps.emplace(key, std::stoll(fields[9].substr(fields[9].rfind(':') + 1)));
    }
    REQUIRE(continuation_ps.size() == 3);
    CHECK(continuation_ps.at("11796979:GTGTA>G") == continuation_ps.at("12735894:TA>T"));
    CHECK(continuation_ps.at("12735894:TA>T") == continuation_ps.at("12955838:T>TC"));
}

static void gap_calibrated_tandem_insertions_close_the_0_865_mb_gap(const Paths& p) {
    Window gap;
    gap.gap_left = 865572;
    gap.gap_right = 882277;
    const auto& truth = load_truth(p.truth_map);
    const Outcome got = measure(p, gap, "graph", "", truth);
    CHECK(got.spans);
    CHECK_FALSE(got.switched);
    CHECK(got.primary_scorable == 71);
    CHECK(got.primary_correct >= 59);
    CHECK(got.core_correct >= 58);
    check_gap_contract(p, gap, got);
    check_read_floors("gap_calibrated_tandem_insertions_close_the_0_865_mb_gap", got);
    std::string dir;
    REQUIRE(run_arm(p, gap, "graph", "", dir));
    std::map<std::string, std::pair<std::string, long long>> calls;
    std::ifstream vcf(dir + "/native.vcf");
    REQUIRE(vcf.good());
    std::string line;
    while (std::getline(vcf, line)) {
        if (line.empty() || line[0] == '#') continue;
        const auto fields = split_tabs(line);
        const std::string key = fields[1] + ":" + fields[3] + ">" + fields[4];
        if (key != "863406:C>A" && key != "882277:A>ATC" && key != "882277:A>ATCTC" && key != "890261:TT>AC") continue;
        if (key == "882277:A>ATC" || key == "882277:A>ATCTC") {
            CHECK(fields[9].find(key == "882277:A>ATC" ? ":28:17,11:" : ":28:13,15:") != std::string::npos);
        }
        REQUIRE(calls.emplace(key, std::make_pair(fields[9].substr(0, 3),
            std::stoll(fields[9].substr(fields[9].rfind(':') + 1)))).second);
    }
    REQUIRE(calls.size() == 4);
    const auto& left = calls.at("863406:C>A");
    for (const auto& [key, call] : calls) {
        CHECK(is_phased_het(call.first));
        CHECK(call.second == left.second);
        CHECK(call.first == (key == "882277:A>ATC" ? (left.first == "0|1" ? "1|0" : "0|1") : left.first));
    }
    const auto& spans = input_read_spans(p, gap);
    std::array<std::array<int, 2>, 2> parents{};
    const std::unique_ptr<samFile, decltype(&hts_close)> bam(
        sam_open((dir + "/phased.bam").c_str(), "r"), hts_close);
    REQUIRE(bam != nullptr);
    const std::unique_ptr<bam_hdr_t, decltype(&bam_hdr_destroy)> header(
        sam_hdr_read(bam.get()), bam_hdr_destroy);
    const std::unique_ptr<bam1_t, decltype(&bam_destroy1)> record(bam_init1(), bam_destroy1);
    REQUIRE(header != nullptr);
    REQUIRE(record != nullptr);
    while (sam_read1(bam.get(), header.get(), record.get()) >= 0) {
        if (record->core.flag & (BAM_FSECONDARY | BAM_FSUPPLEMENTARY)) continue;
        const std::string name = bam_get_qname(record.get());
        const auto parent = truth.find(name);
        const auto span = spans.find(name);
        const uint8_t* hp = bam_aux_get(record.get(), "HP");
        const uint8_t* ps = bam_aux_get(record.get(), "PS");
        if (parent == truth.end() || span == spans.end() || hp == nullptr || ps == nullptr ||
            bam_aux2i(ps) != left.second || (bam_aux2i(hp) != 1 && bam_aux2i(hp) != 2)) continue;
        const bool mat_on_hap1 = (bam_aux2i(hp) == 1) == (parent->second == 'M');
        if (span->second.first <= gap.gap_left && span->second.second >= gap.gap_left &&
            span->second.second < gap.gap_right) ++parents[0][mat_on_hap1];
        if (span->second.first > gap.gap_left && span->second.second >= gap.gap_right)
            ++parents[1][mat_on_hap1];
    }
    const int orientation = parents[0][1] >= parents[0][0] ? 1 : 0;
    for (const auto& flank : parents) {
        const int total = flank[0] + flank[1];
        REQUIRE(total >= 3);
        CHECK(static_cast<double>(flank[orientation]) / total >= 0.90);
    }
    // A serial union must retain the stitched continuation into the next chunk.
    Window continuation;
    continuation.gap_left = 100001;
    continuation.gap_right = 1117887;
    std::string continuation_dir;
    REQUIRE(run_arm(p, continuation, "graph", "", continuation_dir));
    std::ifstream continuation_vcf(continuation_dir + "/native.vcf");
    REQUIRE(continuation_vcf.good());
    std::map<std::string, long long> continuation_ps;
    std::map<std::string, std::string> continuation_gt;
    while (std::getline(continuation_vcf, line)) {
        if (line.empty() || line[0] == '#') continue;
        const auto fields = split_tabs(line);
        const std::string key = fields[1] + ":" + fields[3] + ">" + fields[4];
        if (key != "130540:T>C" && key != "882277:A>ATC" && key != "1117887:T>A") continue;
        CHECK(is_phased_het(fields[9].substr(0, 3)));
        continuation_gt.emplace(key, fields[9].substr(0, 3));
        continuation_ps.emplace(key, std::stoll(fields[9].substr(fields[9].rfind(':') + 1)));
    }
    REQUIRE(continuation_ps.size() == 3);
    CHECK(continuation_gt.at("130540:T>C") == continuation_gt.at("882277:A>ATC"));
    CHECK(continuation_gt.at("882277:A>ATC") == continuation_gt.at("1117887:T>A"));
    CHECK(continuation_ps.at("130540:T>C") == continuation_ps.at("882277:A>ATC"));
    CHECK(continuation_ps.at("882277:A>ATC") == continuation_ps.at("1117887:T>A"));
}

static void gap_calibrated_mixed_repeat_closes_the_57_085_mb_gap(const Paths& p) {
    Window gap;
    gap.gap_left = 57085410;
    gap.gap_right = 57104654;
    const auto& truth = load_truth(p.truth_map);
    const Outcome got = measure(p, gap, "graph", "", truth);
    CHECK(got.spans);
    CHECK_FALSE(got.switched);
    CHECK(got.primary_scorable == 155);
    CHECK(got.primary_correct >= 143);
    CHECK(got.core_correct >= 140);
    check_gap_contract(p, gap, got);
    check_read_floors("gap_calibrated_mixed_repeat_closes_the_57_085_mb_gap", got);
    std::string dir;
    REQUIRE(run_arm(p, gap, "graph", "", dir));
    std::map<std::string, std::pair<std::string, long long>> calls;
    std::ifstream vcf(dir + "/native.vcf");
    REQUIRE(vcf.good());
    std::string line;
    while (std::getline(vcf, line)) {
        if (line.empty() || line[0] == '#') continue;
        const auto fields = split_tabs(line);
        const std::string key = fields[1] + ":" + fields[3] + ">" + fields[4];
        if (key != "57085410:A>C" && key != "57104654:T>TAAA" && key != "57104654:TAA>T" && key != "57123860:A>C") continue;
        if (key == "57104654:T>TAAA" || key == "57104654:TAA>T") {
            CHECK(fields[9].find(key == "57104654:T>TAAA" ? ":67:44,23:" : ":67:45,22:") != std::string::npos);
        }
        REQUIRE(calls.emplace(key, std::make_pair(fields[9].substr(0, 3),
            std::stoll(fields[9].substr(fields[9].rfind(':') + 1)))).second);
    }
    REQUIRE(calls.size() == 4);
    const auto& left = calls.at("57085410:A>C");
    for (const auto& [key, call] : calls) {
        CHECK(is_phased_het(call.first));
        CHECK(call.second == left.second);
        CHECK(call.first == (key == "57104654:TAA>T" ? (left.first == "0|1" ? "1|0" : "0|1") : left.first));
    }
    const auto& spans = input_read_spans(p, gap);
    std::array<std::array<int, 2>, 2> parents{};
    const std::unique_ptr<samFile, decltype(&hts_close)> bam(
        sam_open((dir + "/phased.bam").c_str(), "r"), hts_close);
    REQUIRE(bam != nullptr);
    const std::unique_ptr<bam_hdr_t, decltype(&bam_hdr_destroy)> header(
        sam_hdr_read(bam.get()), bam_hdr_destroy);
    const std::unique_ptr<bam1_t, decltype(&bam_destroy1)> record(bam_init1(), bam_destroy1);
    REQUIRE(header != nullptr);
    REQUIRE(record != nullptr);
    while (sam_read1(bam.get(), header.get(), record.get()) >= 0) {
        if (record->core.flag & (BAM_FSECONDARY | BAM_FSUPPLEMENTARY)) continue;
        const std::string name = bam_get_qname(record.get());
        const auto parent = truth.find(name);
        const auto span = spans.find(name);
        const uint8_t* hp = bam_aux_get(record.get(), "HP");
        const uint8_t* ps = bam_aux_get(record.get(), "PS");
        if (parent == truth.end() || span == spans.end() || hp == nullptr || ps == nullptr ||
            bam_aux2i(ps) != left.second || (bam_aux2i(hp) != 1 && bam_aux2i(hp) != 2)) continue;
        const bool mat_on_hap1 = (bam_aux2i(hp) == 1) == (parent->second == 'M');
        if (span->second.first <= gap.gap_left && span->second.second >= gap.gap_left &&
            span->second.second < gap.gap_right) ++parents[0][mat_on_hap1];
        if (span->second.first > gap.gap_left && span->second.second >= gap.gap_right)
            ++parents[1][mat_on_hap1];
    }
    const int orientation = parents[0][1] >= parents[0][0] ? 1 : 0;
    for (const auto& flank : parents) {
        const int total = flank[0] + flank[1];
        REQUIRE(total >= 3);
        CHECK(static_cast<double>(flank[orientation]) / total >= 0.90);
    }
    // A serial union must retain the stitched continuation into the next chunk.
    Window continuation;
    continuation.gap_left = 56542357;
    continuation.gap_right = 57465653;
    std::string continuation_dir;
    REQUIRE(run_arm(p, continuation, "graph", "", continuation_dir));
    std::ifstream continuation_vcf(continuation_dir + "/native.vcf");
    REQUIRE(continuation_vcf.good());
    std::map<std::string, long long> continuation_ps;
    std::map<std::string, std::string> continuation_gt;
    while (std::getline(continuation_vcf, line)) {
        if (line.empty() || line[0] == '#') continue;
        const auto fields = split_tabs(line);
        const std::string key = fields[1] + ":" + fields[3] + ">" + fields[4];
        if (key != "56542357:C>T" && key != "57104654:T>TAAA" && key != "57465653:C>A") continue;
        CHECK(is_phased_het(fields[9].substr(0, 3)));
        continuation_gt.emplace(key, fields[9].substr(0, 3));
        continuation_ps.emplace(key, std::stoll(fields[9].substr(fields[9].rfind(':') + 1)));
    }
    REQUIRE(continuation_ps.size() == 3);
    CHECK(continuation_gt.at("56542357:C>T") != continuation_gt.at("57104654:T>TAAA"));
    CHECK(continuation_gt.at("57104654:T>TAAA") == continuation_gt.at("57465653:C>A"));
    CHECK(continuation_ps.at("56542357:C>T") == continuation_ps.at("57104654:T>TAAA"));
    CHECK(continuation_ps.at("57104654:T>TAAA") == continuation_ps.at("57465653:C>A"));
}

static void gap_equivalent_mixed_repeat_closes_the_55_309_mb_endpoint(const Paths& p) {
    Window gap;
    gap.gap_left = 55309789;
    gap.gap_right = 55309794;
    const auto& truth = load_truth(p.truth_map);
    const Outcome got = measure(p, gap, "graph", "", truth);
    CHECK(got.spans);
    CHECK_FALSE(got.switched);
    CHECK(got.primary_scorable == 67);
    CHECK(got.primary_correct >= 66);
    CHECK(got.core_correct >= 66);
    check_gap_contract(p, gap, got);
    check_read_floors("gap_equivalent_mixed_repeat_closes_the_55_309_mb_endpoint", got);
    std::string dir;
    REQUIRE(run_arm(p, gap, "graph", "", dir));
    std::map<std::string, std::pair<std::string, long long>> calls;
    std::ifstream vcf(dir + "/native.vcf");
    REQUIRE(vcf.good());
    std::string line;
    while (std::getline(vcf, line)) {
        if (line.empty() || line[0] == '#') continue;
        const auto fields = split_tabs(line);
        const std::string key = fields[1] + ":" + fields[3] + ">" + fields[4];
        if (key != "55309475:G>A" && key != "55309789:C>CT" && key != "55309789:CTT>C" && key != "55309794:T>TT" && key != "55309794:TTT>T") continue;
        if (key == "55309794:T>TT" || key == "55309794:TTT>T") {
            CHECK(fields[9].find(key == "55309794:T>TT" ? ":41:16,25:" : ":42:16,26:") != std::string::npos);
        }
        REQUIRE(calls.emplace(key, std::make_pair(fields[9].substr(0, 3),
            std::stoll(fields[9].substr(fields[9].rfind(':') + 1)))).second);
    }
    REQUIRE(calls.size() == 5);
    const auto& left = calls.at("55309475:G>A");
    for (const auto& [key, call] : calls) {
        CHECK(is_phased_het(call.first));
        CHECK(call.second == left.second);
        CHECK(call.first == ((key == "55309794:T>TT" || key == "55309789:C>CT") ? (left.first == "0|1" ? "1|0" : "0|1") : left.first));
    }
    const auto& spans = input_read_spans(p, gap);
    std::array<std::array<int, 2>, 2> parents{};
    const std::unique_ptr<samFile, decltype(&hts_close)> bam(
        sam_open((dir + "/phased.bam").c_str(), "r"), hts_close);
    REQUIRE(bam != nullptr);
    const std::unique_ptr<bam_hdr_t, decltype(&bam_hdr_destroy)> header(
        sam_hdr_read(bam.get()), bam_hdr_destroy);
    const std::unique_ptr<bam1_t, decltype(&bam_destroy1)> record(bam_init1(), bam_destroy1);
    REQUIRE(header != nullptr);
    REQUIRE(record != nullptr);
    while (sam_read1(bam.get(), header.get(), record.get()) >= 0) {
        if (record->core.flag & (BAM_FSECONDARY | BAM_FSUPPLEMENTARY)) continue;
        const std::string name = bam_get_qname(record.get());
        const auto parent = truth.find(name);
        const auto span = spans.find(name);
        const uint8_t* hp = bam_aux_get(record.get(), "HP");
        const uint8_t* ps = bam_aux_get(record.get(), "PS");
        if (parent == truth.end() || span == spans.end() || hp == nullptr || ps == nullptr ||
            bam_aux2i(ps) != left.second || (bam_aux2i(hp) != 1 && bam_aux2i(hp) != 2)) continue;
        const bool mat_on_hap1 = (bam_aux2i(hp) == 1) == (parent->second == 'M');
        if (span->second.first <= 55300140 && span->second.second >= 55300140 &&
            span->second.second < gap.gap_left) ++parents[0][mat_on_hap1];
        if (span->second.first > 55309475 && span->second.second >= gap.gap_right)
            ++parents[1][mat_on_hap1];
    }
    const int orientation = parents[0][1] >= parents[0][0] ? 1 : 0;
    for (const auto& flank : parents) {
        const int total = flank[0] + flank[1];
        // Two held-out reads start beyond the last SNP and carry only the endpoint.
        REQUIRE(total >= 2);
        CHECK(static_cast<double>(flank[orientation]) / total >= 0.90);
    }
    // A serial union must retain the stitched continuation into the next chunk.
    Window continuation;
    continuation.gap_left = 54453020;
    continuation.gap_right = 55309794;
    std::string continuation_dir;
    REQUIRE(run_arm(p, continuation, "graph", "", continuation_dir));
    std::ifstream continuation_vcf(continuation_dir + "/native.vcf");
    REQUIRE(continuation_vcf.good());
    std::map<std::string, long long> continuation_ps;
    std::map<std::string, std::string> continuation_gt;
    while (std::getline(continuation_vcf, line)) {
        if (line.empty() || line[0] == '#') continue;
        const auto fields = split_tabs(line);
        const std::string key = fields[1] + ":" + fields[3] + ">" + fields[4];
        if (key != "54453020:G>A" && key != "55309794:T>TT" && key != "55309789:C>CT") continue;
        CHECK(is_phased_het(fields[9].substr(0, 3)));
        continuation_gt.emplace(key, fields[9].substr(0, 3));
        continuation_ps.emplace(key, std::stoll(fields[9].substr(fields[9].rfind(':') + 1)));
    }
    REQUIRE(continuation_ps.size() == 3);
    CHECK(continuation_ps.at("54453020:G>A") == continuation_ps.at("55309794:T>TT"));
    CHECK(continuation_ps.at("55309794:T>TT") == continuation_ps.at("55309789:C>CT"));
}

/// One fixture and selector for every panel window and mechanism regression.
static void gap_compound_insertion_prefix_closes_the_45_876_mb_gap(const Paths& p) {
    Window gap;
    gap.gap_left = 45876157;
    gap.gap_right = 45896820;
    const auto& truth = load_truth(p.truth_map);
    const Outcome got = measure(p, gap, "graph", "", truth);
    CHECK(got.spans);
    CHECK_FALSE(got.switched);
    CHECK(got.primary_scorable == 141);
    CHECK(got.primary_correct >= 129);
    CHECK(got.core_correct >= 129);
    check_gap_contract(p, gap, got);
    check_read_floors("gap_compound_insertion_prefix_closes_the_45_876_mb_gap", got);
    const auto check_markers = [](const std::string& dir, bool require_end) {
        std::map<std::string, std::pair<std::string, long long>> calls;
        std::ifstream vcf(dir + "/native.vcf");
        REQUIRE(vcf.good());
        std::string line;
        while (std::getline(vcf, line)) {
            if (line.empty() || line[0] == '#') continue;
            const auto fields = split_tabs(line);
            const std::string key = fields[1] + ":" + fields[3] + ">" + fields[4];
            if (key != "45866904:G>GACAGACAGACACACACAC" && key != "45866904:G>GACAGACAGACACACACACAC" &&
                key != "45876157:G>GT" && key != "45896820:A>G" && key != "46636707:C>T") continue;
            REQUIRE(calls.emplace(key, std::make_pair(fields[9].substr(0, 3),
                std::stoll(fields[9].substr(fields[9].rfind(':') + 1)))).second);
        }
        REQUIRE(calls.count("45866904:G>GACAGACAGACACACACAC") == 1);
        REQUIRE(calls.count("45866904:G>GACAGACAGACACACACACAC") == 1);
        REQUIRE(calls.count("45876157:G>GT") == 1);
        REQUIRE(calls.count("45896820:A>G") == 1);
        if (require_end) REQUIRE(calls.count("46636707:C>T") == 1);
        const auto& left = calls.at("45866904:G>GACAGACAGACACACACACAC");
        for (const auto& [key, call] : calls) {
            CHECK(is_phased_het(call.first));
            CHECK(call.second == left.second);
        }
        CHECK(calls.at("45896820:A>G").first == left.first);
        CHECK(calls.at("45866904:G>GACAGACAGACACACACAC").first != left.first);
        return left.second;
    };
    std::string dir;
    REQUIRE(run_arm(p, gap, "graph", "", dir));
    const long long core = check_markers(dir, false);
    const auto& spans = input_read_spans(p, gap);
    std::array<std::array<int, 2>, 2> parents{};
    const std::unique_ptr<samFile, decltype(&hts_close)> bam(sam_open((dir + "/phased.bam").c_str(), "r"), hts_close);
    REQUIRE(bam != nullptr);
    const std::unique_ptr<bam_hdr_t, decltype(&bam_hdr_destroy)> header(sam_hdr_read(bam.get()), bam_hdr_destroy);
    const std::unique_ptr<bam1_t, decltype(&bam_destroy1)> record(bam_init1(), bam_destroy1);
    REQUIRE(header != nullptr);
    REQUIRE(record != nullptr);
    while (sam_read1(bam.get(), header.get(), record.get()) >= 0) {
        if (record->core.flag & (BAM_FSECONDARY | BAM_FSUPPLEMENTARY)) continue;
        const std::string name = bam_get_qname(record.get());
        const auto parent = truth.find(name);
        const auto span = spans.find(name);
        const uint8_t* hp = bam_aux_get(record.get(), "HP");
        const uint8_t* ps = bam_aux_get(record.get(), "PS");
        if (parent == truth.end() || span == spans.end() || hp == nullptr || ps == nullptr ||
            bam_aux2i(ps) != core || (bam_aux2i(hp) != 1 && bam_aux2i(hp) != 2)) continue;
        const bool mat_on_hap1 = (bam_aux2i(hp) == 1) == (parent->second == 'M');
        // The compound-only cohort has noisy repeat lengths. The original
        // bridge molecules and an independent SNP cohort verify the gauge.
        if (span->second.first <= 45866905 && span->second.second >= 45883707) ++parents[0][mat_on_hap1];
        if (span->second.first > 45866905 && span->second.second >= gap.gap_right) ++parents[1][mat_on_hap1];
    }
    const int orientation = parents[0][1] >= parents[0][0] ? 1 : 0;
    for (const auto& flank : parents) {
        const int total = flank[0] + flank[1];
        REQUIRE(total >= 3);
        CHECK(static_cast<double>(flank[orientation]) / total >= 0.90);
    }
    Window continuation;
    continuation.gap_left = 45866904;
    continuation.gap_right = 46636707;
    std::string continuation_dir;
    REQUIRE(run_arm(p, continuation, "graph", "", continuation_dir));
    check_markers(continuation_dir, true);
}


static void gap_long_repeat_anchors_close_the_15_mb_block(const Paths& p) {
    const auto& truth = load_truth(p.truth_map);
    for (const auto& bounds : {std::array<long long, 5>{15023123, 15039543, 126, 101, 101},
                              std::array<long long, 5>{15100456, 15101262, 73, 70, 70},
                              std::array<long long, 5>{15095642, 15101261, 101, 95, 95}}) {
        Window gap;
        gap.gap_left = bounds[0]; gap.gap_right = bounds[1];
        const Outcome got = measure(p, gap, "graph", "", truth);
        CHECK(got.spans);
        CHECK_FALSE(got.switched);
        CHECK(got.primary_scorable == bounds[2]);
        CHECK(got.primary_correct >= bounds[3]);
        CHECK(got.core_correct >= bounds[4]);
        check_gap_contract(p, gap, got);
        const bool owning_chunk = gap.gap_left != 15095642;
        if (owning_chunk) check_read_floors("gap_long_repeat_anchors_close_the_15_mb_block", got);
        const long long joined_ps = owning_chunk ? 15003748 : 15071132;
        std::string dir;
        REQUIRE(run_arm(p, gap, "graph", "", dir));
        const auto& spans = input_read_spans(p, gap);
        std::array<std::array<int, 2>, 2> parents{};
        std::map<std::string, bool> affected_parents;
        const std::unique_ptr<samFile, decltype(&hts_close)> bam(sam_open((dir + "/phased.bam").c_str(), "r"), hts_close);
        REQUIRE(bam != nullptr);
        const std::unique_ptr<bam_hdr_t, decltype(&bam_hdr_destroy)> header(sam_hdr_read(bam.get()), bam_hdr_destroy);
        const std::unique_ptr<bam1_t, decltype(&bam_destroy1)> record(bam_init1(), bam_destroy1);
        REQUIRE(header != nullptr); REQUIRE(record != nullptr);
        while (sam_read1(bam.get(), header.get(), record.get()) >= 0) {
            if (record->core.flag & (BAM_FSECONDARY | BAM_FSUPPLEMENTARY)) continue;
            const std::string name = bam_get_qname(record.get());
            const auto parent = truth.find(name);
            const auto span = spans.find(name);
            const uint8_t* hp = bam_aux_get(record.get(), "HP");
            const uint8_t* ps = bam_aux_get(record.get(), "PS");
            if (parent == truth.end() || span == spans.end() || hp == nullptr || ps == nullptr ||
                bam_aux2i(ps) != joined_ps || (bam_aux2i(hp) != 1 && bam_aux2i(hp) != 2)) continue;
            const bool mat_on_hap1 = (bam_aux2i(hp) == 1) == (parent->second == 'M');
            if (name == "m84031_231217_034919_s2/237175194/ccs" ||
                name == "m84031_231217_034919_s2/82248351/ccs" ||
                name == "m84031_231217_062403_s3/235673333/ccs")
                affected_parents.emplace(name, mat_on_hap1);
            // Independent cohorts on either side of each physical repeat
            // check that the joined core retains one parental orientation.
            const long long marker = gap.gap_left == 15023123 ? 15019256 : 15101263;
            if (span->second.first <= marker && span->second.second < gap.gap_right) ++parents[0][mat_on_hap1];
            if (span->second.first > marker && span->second.second >= gap.gap_right) ++parents[1][mat_on_hap1];
        }
        const int orientation = parents[0][1] >= parents[0][0] ? 1 : 0;
        if (gap.gap_left != 15023123) {
            REQUIRE(affected_parents.size() == 3);
            for (const auto& [name, mat_on_hap1] : affected_parents) {
                INFO("retained source insertion and tandem veto: " << name);
                CHECK(mat_on_hap1 == static_cast<bool>(orientation));
            }
        }
        for (const auto& cohort : parents) {
            const int total = cohort[0] + cohort[1];
            REQUIRE(total >= 3);
            CHECK(static_cast<double>(cohort[orientation]) / total >= 0.90);
        }
    }
    Window continuation;
    continuation.gap_left = 14719378; continuation.gap_right = 15335938;
    std::string dir;
    REQUIRE(run_arm(p, continuation, "graph", "", dir));
    std::map<std::string, std::pair<std::string, long long>> calls;
    const std::set<std::string> required{
        "14719378:C>T", "15019255:TACACAC>T", "15023123:CT>C", "15039543:A>G",
        "15100456:CT>C", "15100456:CTT>C", "15101262:A>ATG", "15101262:A>ATGTGTGTG", "15109300:C>CA", "15335938:G>A"};
    std::ifstream vcf(dir + "/native.vcf");
    REQUIRE(vcf.good());
    std::string line;
    while (std::getline(vcf, line)) {
        if (line.empty() || line[0] == '#') continue;
        const auto fields = split_tabs(line);
        const std::string key = fields[1] + ":" + fields[3] + ">" + fields[4];
        if (required.count(key)) REQUIRE(calls.emplace(key, std::make_pair(fields[9].substr(0, 3),
            std::stoll(fields[9].substr(fields[9].rfind(':') + 1)))).second);
    }
    REQUIRE(calls.size() == required.size());
    for (const auto& [key, call] : calls) {
        CHECK(is_phased_het(call.first));
        CHECK(call.second == 14719378);
    }
    CHECK(calls.at("15019255:TACACAC>T").first == calls.at("15039543:A>G").first);
    CHECK(calls.at("15100456:CT>C").first != calls.at("15100456:CTT>C").first);
    CHECK(calls.at("15101262:A>ATG").first != calls.at("15101262:A>ATGTGTGTG").first);
    CHECK(calls.at("15101262:A>ATGTGTGTG").first != calls.at("15039543:A>G").first);
    CHECK(calls.at("15109300:C>CA").first == calls.at("15101262:A>ATGTGTGTG").first);
}


static void gap_calibrated_deletion_chain_closes_the_56_mb_boundary(const Paths& p) {
    Window gap;
    gap.gap_left = 55999194; gap.gap_right = 56040612;
    const auto& truth = load_truth(p.truth_map);
    const Outcome got = measure(p, gap, "graph", "", truth);
    CHECK(got.spans);
    CHECK_FALSE(got.switched);
    CHECK(got.primary_scorable == 246);
    CHECK(got.primary_correct >= 214);
    CHECK(got.core_correct >= 214);
    check_gap_contract(p, gap, got);
    check_read_floors("gap_calibrated_deletion_chain_closes_the_56_mb_boundary", got);
    std::string dir;
    REQUIRE(run_arm(p, gap, "graph", "", dir));
    const auto& spans = input_read_spans(p, gap);
    std::array<std::array<int, 2>, 2> parents{};
    const std::set<std::string> restored{
        "m84031_231217_034919_s2/208800100/ccs", "m84031_231217_062403_s3/152046853/ccs",
        "m84031_231217_062403_s3/254677189/ccs", "m84031_231217_034919_s2/80151528/ccs"};
    std::set<std::string> checked;
    const std::unique_ptr<samFile, decltype(&hts_close)> bam(sam_open((dir + "/phased.bam").c_str(), "r"), hts_close);
    REQUIRE(bam != nullptr);
    const std::unique_ptr<bam_hdr_t, decltype(&bam_hdr_destroy)> header(sam_hdr_read(bam.get()), bam_hdr_destroy);
    const std::unique_ptr<bam1_t, decltype(&bam_destroy1)> record(bam_init1(), bam_destroy1);
    REQUIRE(header != nullptr); REQUIRE(record != nullptr);
    while (sam_read1(bam.get(), header.get(), record.get()) >= 0) {
        const std::string name = bam_get_qname(record.get());
        const uint8_t* hp = bam_aux_get(record.get(), "HP");
        const uint8_t* ps = bam_aux_get(record.get(), "PS");
        if (name == "m84031_231217_062403_s3/151328165/ccs") CHECK((hp == nullptr || bam_aux2i(hp) == 0));
        const auto parent = truth.find(name);
        const auto span = spans.find(name);
        if (parent == truth.end() || span == spans.end() || hp == nullptr || ps == nullptr ||
            bam_aux2i(ps) != 55883019 || (bam_aux2i(hp) != 1 && bam_aux2i(hp) != 2)) continue;
        const bool mat_on_hap1 = (bam_aux2i(hp) == 1) == (parent->second == 'M');
        if (restored.count(name)) { CHECK(mat_on_hap1); checked.insert(name); }
        if (span->second.first <= gap.gap_left && span->second.second < 56007501) ++parents[0][mat_on_hap1];
        if (span->second.first > 56027379 && span->second.second >= gap.gap_right) ++parents[1][mat_on_hap1];
    }
    CHECK(checked == restored);
    for (const auto& cohort : parents) {
        REQUIRE(cohort[0] + cohort[1] >= 3);
        CHECK(static_cast<double>(cohort[1]) / (cohort[0] + cohort[1]) >= 0.90);
    }
}

TEST_CASE("all gaps", "[gap][windows][integration]") {
    const Paths p = paths();
    if (!p.complete()) {
        WARN("gap-window tests need inputs that are absent: " << p.missing());
        SUCCEED("skipped: inputs absent");
        return;
    }
    struct GapCheck {
        const char* name;
        const char* tags;
        void (*run)(const Paths&);
    };
    static const GapCheck checks[] = {
        {"calibrated deletion chain closes the 56 Mb boundary", "[gap][deletion-chain][orientation]", gap_calibrated_deletion_chain_closes_the_56_mb_boundary},
        {"Calibrated tandem insertions close the 0.865 Mb gap", "[gap][tandem-insertion][orientation]", gap_calibrated_tandem_insertions_close_the_0_865_mb_gap},
        {"calibrated mixed repeat closes the 57.085 Mb gap", "[gap][stitch-connectivity][msa][orientation]", gap_calibrated_mixed_repeat_closes_the_57_085_mb_gap},
        {"equivalent mixed repeat closes the 55.309 Mb endpoint", "[gap][terminal-repeat][orientation]", gap_equivalent_mixed_repeat_closes_the_55_309_mb_endpoint},
        {"long repeat anchors close the 15.023 and 15.100 Mb gaps", "[gap][long-repeat][orientation]", gap_long_repeat_anchors_close_the_15_mb_block},
        {"compound insertion prefix closes the 45.876 Mb gap", "[gap][compound-prefix][orientation]", gap_compound_insertion_prefix_closes_the_45_876_mb_gap},
        {"Physical terminal insertion completes the 12.954 Mb block", "[gap][terminal-insertion][orientation]", gap_physical_terminal_insertion_closes_the_12_954_mb_gap},
        {"Physical graph repeat chain closes the 12.717 Mb gap", "[gap][repeat-indel-chain][orientation]", gap_graph_repeat_indel_chain_closes_the_12_717_mb_gap},
        {"Physical repeat SNP validation closes the 11.796 Mb gap", "[gap][repeat-snp-validation][orientation]", gap_physically_validated_repeat_snp_closes_the_11_796_mb_gap},
        {"REF-absent graph SNP recovers the 40.633 Mb terminal deletion", "[gap][terminal-deletion][orientation]", gap_ref_absent_graph_snp_recovers_the_40_633_mb_terminal_deletion},
        {"cut-free insertion prefix closes the 36.268 Mb gap", "[gap][insertion-prefix][orientation]", gap_cut_free_insertion_prefix_closes_the_36_268_mb_gap},
        {"homozygous graph SNP does not split a supported block", "[gap][homozygous-graph-snp][orientation]", gap_homozygous_graph_snp_does_not_split_a_supported_block},
        {"calibrated source deletion joins the fourth largest HiPhase block", "[gap][fourth-largest-block][source-deletion][orientation]", gap_calibrated_source_deletion_joins_the_fourth_largest_hiphase_block},
        {"physically contradicted graph SNPs close the 25.855 Mb seam", "[gap][graph-snp-contradiction][orientation]", gap_physically_contradicted_graph_snps_close_the_25_855_mb_seam},
        {"complementary insertions join the second largest HiPhase block", "[gap][second-largest-block][orientation]", gap_complementary_insertions_join_the_second_largest_hiphase_block},
        {"calibrated repeat SNPs join the largest HiPhase block", "[gap][repeat][largest-block][orientation]", gap_calibrated_repeat_snps_join_the_largest_hiphase_block},
        {"cut-free source paths close the 17.634 Mb gap", "[gap][shared-deletion][source-path][orientation]", gap_cut_free_source_paths_close_the_17_634_mb_gap},
        {"retained shared deletion closes the 62.408 Mb gap", "[gap][shared-deletion][masked-bam-snp][orientation]", gap_retained_shared_deletion_closes_the_62_408_mb_gap},
        {"original MSA SNP dropout is retried before CIGAR backfill", "[gap][stitch-connectivity][msa]", gap_original_msa_snp_dropout_is_retried_before_cigar_backfill},
        {"conflicting sparse pairs retry with their observed binomial tail", "[gap][stitch-connectivity][msa][sparse-pair]", gap_conflicting_sparse_pairs_retry_with_their_observed_binomial_tail},
        {"a moved deletion does not certify a whole right-block join", "[gap][stitch-connectivity][msa][moved-deletion]", gap_a_moved_deletion_does_not_certify_a_whole_right_block_join},
        {"focused retry certifies its path after CIGAR backfill", "[gap][stitch-connectivity][msa][backfill-certificate]", gap_focused_retry_certifies_its_path_after_cigar_backfill},
        {"an earlier focused retry does not hide a later MSA dropout", "[gap][msa][retry-scheduling]", gap_an_earlier_focused_retry_does_not_hide_a_later_msa_dropout},
        {"complete recovery MSA blocks retain the complex left flank", "[gap][msa-complete-transfer]", gap_complete_recovery_msa_blocks_retain_the_complex_left_flank},
        {"isolated shared BAM genotype survives graph repeat demotion", "[gap][representation]", gap_isolated_shared_bam_genotype_survives_graph_repeat_demotion},
        {"chr20 gap windows", "[gap][windows]", gap_chr20_gap_windows},
        {"recovery preserves complementary BAM deletion rows", "[gap][representation]", gap_recovery_preserves_complementary_bam_deletion_rows},
        {"verified MSA insertion survives graph output classification", "[gap][representation]", gap_verified_msa_insertion_survives_graph_output_classification},
        {"a phased BAM deletion survives an unphased graph duplicate", "[gap][representation]", gap_a_phased_bam_deletion_survives_an_unphased_graph_duplicate},
        {"graph SNP retry preserves a verified BAM deletion connection", "[gap][stitch-confidence]", gap_graph_snp_retry_preserves_a_verified_bam_deletion_connection},
        {"an unsupported BAM source cut keeps the far deletion independent", "[gap][stitch-connectivity]", gap_an_unsupported_bam_source_cut_keeps_the_far_deletion_independent},
        {"graph SNPs bridge supported recovery blocks", "[gap][graph-bridge]", gap_graph_snps_bridge_supported_recovery_blocks},
        {"recovery stitch preserves the next flank across an unlinked seam", "[gap][stitch-connectivity]", gap_recovery_stitch_preserves_the_next_flank_across_an_unlinked_seam},
        {"a clean-SNP gap bridge keeps the prior graph gauge", "[gap][stitch-connectivity]", gap_a_clean_snp_gap_bridge_keeps_the_prior_graph_gauge},
        {"a local source path does not lose its graph bridge", "[gap][stitch-connectivity]", gap_a_local_source_path_does_not_lose_its_graph_bridge},
        {"corroborated SNP molecules join sparse graph seams", "[gap][stitch-connectivity]", gap_corroborated_snp_molecules_join_sparse_graph_seams},
        {"an indel boundary with allele dropout gets an MSA retry", "[gap][stitch-connectivity]", gap_an_indel_boundary_with_allele_dropout_gets_an_msa_retry},
        {"a lone boundary SNP pair needs a repaired complete graph flank", "[gap][stitch-connectivity]", gap_a_lone_boundary_snp_pair_needs_a_repaired_complete_graph_flank},
        {"physical bridges survive their owning graph chunks", "[gap][stitch-connectivity]", gap_physical_bridges_survive_their_owning_graph_chunks},
        {"phased right reads retain the certified 64 Mb BAM prefix", "[gap][representation]", gap_phased_right_reads_retain_the_certified_64_mb_bam_prefix},
        {"source allele repair preserves complementary repeat rows", "[gap][representation]", gap_source_allele_repair_preserves_complementary_repeat_rows},
        {"graph haplotypes bridge equivalent complementary BAM deletions", "[gap][stitch-connectivity]", gap_graph_haplotypes_bridge_equivalent_complementary_bam_deletions},
        {"physical links close the long 47 Mb gap without a switch", "[gap][stitch-connectivity]", gap_physical_links_close_the_long_47_mb_gap_without_a_switch},
        {"validated graph bridge survives a BAM phase-set ID change", "[gap][stitch-connectivity]", gap_validated_graph_bridge_survives_a_bam_phase_set_id_change},
        {"one shared BAM anchor cannot absorb the earlier 19.4 Mb block", "[gap][stitch-connectivity]", gap_one_shared_bam_anchor_cannot_absorb_the_earlier_19_4_mb_block},
        {"direct SNP proof reuses a graph block at 56.15 Mb", "[gap][stitch-connectivity]", gap_direct_snp_proof_reuses_a_graph_block_at_56_15_mb},
        {"physical SNP bridge crosses the 23 Mb graph chunk boundary", "[gap][stitch-connectivity]", gap_physical_snp_bridge_crosses_the_23_mb_graph_chunk_boundary},
        {"multiple BAM source labels preserve a 22.98 Mb graph bridge", "[gap][stitch-connectivity]", gap_multiple_bam_source_labels_preserve_a_22_98_mb_graph_bridge},
        {"a different deletion locus cannot validate a remapped bridge", "[gap][stitch-connectivity]", gap_a_different_deletion_locus_cannot_validate_a_remapped_bridge},
        {"BAM transfer preserves an already connected graph phase set", "[gap][stitch-connectivity]", gap_bam_transfer_preserves_an_already_connected_graph_phase_set},
        {"complementary BAM boundary rows close the 41.900 Mb seam", "[gap][stitch-connectivity]", gap_complementary_bam_boundary_rows_close_the_41_900_mb_seam},
        {"certified deletion bridge preserves the 11.599 Mb source alleles", "[gap][stitch-connectivity][msa][source-path]", gap_certified_deletion_bridge_preserves_the_11_599_mb_source_alleles},
        {"complete BAM source path closes the 48.929 Mb seam in its owning chunk", "[gap][stitch-connectivity]", gap_complete_bam_source_path_closes_the_48_929_mb_seam_in_its_owning_chunk},
        {"complete adjacent graph phase sets close the 56 Mb seam in its owning chunk", "[gap][stitch-connectivity]", gap_complete_adjacent_graph_phase_sets_close_the_56_mb_seam_in_its_owning_chunk},
        {"complete BAM path and verified deletion close the 5.31 Mb seam", "[gap][stitch-connectivity]", gap_complete_bam_path_and_verified_deletion_close_the_5_31_mb_seam},
        {"repeat-deletion recovery preserves the 61.738 Mb flank gauges", "[gap][recovery][stitch-connectivity]", gap_repeat_deletion_recovery_preserves_the_61_738_mb_flank_gauges},
        {"BAM-left MEC cannot reverse the 34.1 Mb graph block", "[gap][stitch-connectivity]", gap_bam_left_mec_cannot_reverse_the_34_1_mb_graph_block},
        {"BAM recovery admits missing MSA pairs at 17.62 Mb", "[gap][stitch-connectivity]", gap_bam_recovery_admits_missing_msa_pairs_at_17_62_mb},
        {"BAM block attaches to one supported graph flank", "[gap][stitch-connectivity]", gap_bam_block_attaches_to_one_supported_graph_flank},
        {"BAM recovery uses a supported run before a weak source cut", "[gap][stitch-connectivity]", gap_bam_recovery_uses_a_supported_run_before_a_weak_source_cut},
        {"BAM source path with one-haplotype molecule support closes 34.844 Mb", "[gap][stitch-connectivity]", gap_bam_source_path_with_one_haplotype_molecule_support_closes_34_844_mb},
        {"clean indel boundaries use the phased read path", "[gap][stitch-connectivity]", gap_clean_indel_boundaries_use_the_phased_read_path},
        {"one weak-cut BAM run attaches to one supported neighbor", "[gap][stitch-connectivity]", gap_one_weak_cut_bam_run_attaches_to_one_supported_neighbor},
        {"complete BAM blocks use split-stable allele votes without a read gauge", "[gap][stitch-connectivity]", gap_complete_bam_blocks_use_split_stable_allele_votes_without_a_read_gauge},
        {"BAM-supported inner block connects after an invalid graph SNP is excluded", "[gap][stitch-connectivity]", gap_bam_supported_inner_block_connects_after_an_invalid_graph_snp_is_excluded},
        {"low-MAPQ BAM evidence cannot veto a graph heterozygote", "[gap][anchor-quality]", gap_low_mapq_bam_evidence_cannot_veto_a_graph_heterozygote},
        {"shifted BAM deletion joins a certified graph suffix", "[gap][stitch-connectivity]", gap_shifted_bam_deletion_joins_a_certified_graph_suffix},
        {"a directly supported BAM suffix closes the 17.865 Mb gap", "[gap][stitch-connectivity]", gap_a_directly_supported_bam_suffix_closes_the_17_865_mb_gap},
        {"a newly imported BAM seam closes the 15.056 Mb gap", "[gap][stitch-connectivity]", gap_a_newly_imported_bam_seam_closes_the_15_056_mb_gap},
        {"recovered right deletion joins a source-backed graph block", "[gap][stitch-connectivity]", gap_recovered_right_deletion_joins_a_source_backed_graph_block},
        {"adjacent anchors do not trigger a second BAM recovery solve", "[gap][stitch-connectivity]", gap_adjacent_anchors_do_not_trigger_a_second_bam_recovery_solve},
        {"shifted single-base BAM insertion closes its owning 23 Mb gap", "[gap][stitch-connectivity][representation]", gap_shifted_single_base_bam_insertion_closes_its_owning_23_mb_gap},
        {"a single seam admits focused recovery of complementary BAM rows", "[gap][stitch-connectivity][representation]", gap_a_single_seam_admits_focused_recovery_of_complementary_bam_rows},
        {"an internal MSA conflict does not authorize an unsupported block join", "[gap][stitch-connectivity][source-conflict-guard]", gap_an_internal_msa_conflict_does_not_authorize_an_unsupported_block_join},
        {"focused recovery retains a supported partial path inside a graph seam", "[gap][stitch-connectivity][representation][partial-source]", gap_focused_recovery_retains_a_supported_partial_path_inside_a_graph_seam},
        {"focused BAM recovery retains insertion representations of graph child SNPs", "[gap][representation]", gap_focused_bam_recovery_retains_insertion_representations_of_graph_child_snps},
        {"chr20 gap windows: panel totals", "[gap][windows][totals]", gap_chr20_gap_windows_panel_totals},
        {"observed BAM insertion runs retain long alleles before the graph boundary", "[gap][stitch-connectivity][representation][long-insertion-run]", gap_observed_bam_insertion_runs_retain_long_alleles_before_the_graph_boundary},
        {"equivalent BAM and graph insertions close their repeat seam", "[gap][stitch-connectivity][representation][insertion-equivalence]", gap_equivalent_bam_and_graph_insertions_close_their_repeat_seam},
        {"an attached BAM component keeps its source identity at an insertion bridge", "[gap][representation][source-component]", gap_an_attached_bam_component_keeps_its_source_identity_at_an_insertion_bridge},
        {"independent BAM SNP pairs certify both sides of a deletion seam", "[gap][deletion][physical-path]", gap_independent_bam_snp_pairs_certify_both_sides_of_a_deletion_seam},
        {"padded graph repeats retain the verified BAM genotype", "[gap][representation][transfer]", gap_padded_graph_repeats_retain_the_verified_bam_genotype},
        {"one-haplotype repeat deletion votes preserve parental block orientation", "[gap][deletion][orientation]", gap_one_haplotype_repeat_deletion_votes_preserve_parental_block_orientation},
        {"exact SNP branches preserve the path through a complex catalog allele", "[gap][representation][snp-branch]", gap_exact_snp_branches_preserve_the_path_through_a_complex_catalog_allele},
        {"complementary insertion recall cannot invert its SNP flanks", "[gap][msa][source-gap][orientation]", gap_complementary_insertion_recall_cannot_invert_its_snp_flanks},
        {"overlapping recovery solves cannot exchange SNP quality certificates", "[gap][msa][quality-provenance]", gap_overlapping_recovery_solves_cannot_exchange_snp_quality_certificates},
        {"complementary insertion evidence preserves the 17.50 Mb connection", "[gap][msa][insertion-connection][orientation]", gap_complementary_insertion_evidence_preserves_the_17_50_mb_connection},
        {"complete BAM evidence preserves finalized 37 Mb block orientation", "[gap][recovery][complete-block][stitch-connectivity]", gap_complete_bam_evidence_preserves_finalized_37_mb_block_orientation},
        {"complementary deletion recovery closes the owning 50 Mb gap", "[gap][recovery][deletion-pair][stitch-connectivity]", gap_complementary_deletion_recovery_closes_the_owning_50_mb_gap},
        {"shifted compound CIGAR deletions close the owning 7.9 Mb gap", "[gap][msa][representation][deletion-pair]", gap_shifted_compound_cigar_deletions_close_the_owning_7_9_mb_gap},
        {"whole-chunk BAM fallback preserves verified recovery observations", "[gap][msa][transfer][observations]", gap_whole_chunk_bam_fallback_preserves_verified_recovery_observations},
        {"selected BAM blocks retain verified MSA calls across graph gauge disagreement", "[gap][msa][transfer][source-admission]", gap_selected_bam_blocks_retain_verified_msa_calls_across_graph_gauge_disagreement},
        {"BAM-private blocks retain verified MSA calls without graph anchors", "[gap][msa][transfer][private-block]", gap_bam_private_blocks_retain_verified_msa_calls_without_graph_anchors},
        {"verified noisy SNP and symmetric graph path close the 41.881 Mb gap", "[gap][msa][stitch-connectivity][graph-path]", gap_verified_noisy_snp_and_symmetric_graph_path_close_the_41_881_mb_gap},
        {"composed physical bridges close 23.461 Mb and preserve the preceding join", "[gap][stitch-connectivity][msa][source-path][composed-bridge]", gap_composed_physical_bridges_close_23_461_mb_and_preserve_the_preceding_join},
        {"deferred physical bridge preserves read rescue across 22-24 Mb", "[gap][stitch-connectivity][composed-bridge][read-rescue]", gap_deferred_physical_bridge_preserves_read_rescue_across_22_24_mb},
        {"complementary insertion ALTs bridge the 9 Mb chunk boundary", "[gap][stitch-connectivity][msa][complementary-insertion][chunk-boundary]", gap_complementary_insertion_alts_bridge_the_9_mb_chunk_boundary},
        {"physical SNP switch repair closes the owning 62 Mb gap", "[gap][physical-switch][orientation]", gap_physical_snp_switch_repair_closes_the_owning_62_mb_gap},
        {"quality-backed source path closes the owning 13 Mb gap", "[gap][source-quality-path][orientation]", gap_quality_backed_source_path_closes_the_owning_13_mb_gap},
        {"a recovered short insertion retains its independent source path", "[gap][short-insertion-source][orientation]", gap_a_recovered_short_insertion_retains_its_independent_source_path},
        {"shared graph rows retain their independent BAM run certificate", "[gap][shared-source-snp][orientation]", gap_shared_graph_rows_retain_their_independent_bam_run_certificate},
        {"a clean SNP retry certifies a complete independent BAM suffix", "[gap][clean-snp-source-retry][orientation]", gap_a_clean_snp_retry_certifies_a_complete_independent_bam_suffix},
    };
    const std::string filter = env_or("PGPHASE_GAP_FILTER", "");
    const bool emitting = std::getenv("PGPHASE_EMIT_EXPECTATIONS") != nullptr;
    int selected = 0;
    for (const auto& check : checks) {
        if (emitting && std::string(check.name) != "chr20 gap windows") continue;
        if (!filter.empty() && std::string(check.name).find(filter) == std::string::npos &&
            std::string(check.tags).find(filter) == std::string::npos) continue;
        ++selected;
        DYNAMIC_SECTION(check.name) {
            INFO("regression tags: " << check.tags);
            check.run(p);
        }
    }
    REQUIRE(selected > 0);
}
