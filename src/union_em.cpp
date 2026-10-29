// Union gap phasing: the global EM over the read x site matrix, with block
// cuts and the weak-bridge guard.

#include "trace.hpp"
#include "union_internal.hpp"

#include <array>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <numeric>

namespace pgphase_collect {

// Error-rate bounds; starting error for clean and other sites; the learned
// error at which a site is neither reported phased nor used to label reads;
// iteration limits; the read posterior needed for a label; the boundary
// log-likelihood below which blocks are cut; and the MAPQ from which a read
// teaches the model (lower reads are labelled from it but do not shape it).
constexpr double kEmMinError = 0.01;
// A site called by more than this many times the chunk's median callers sits
// in a collapsed duplication: its "alleles" are often the copies' differences,
// which read in lockstep and carried joins at margins in the hundreds
// (chr18:21.05 Mb, 140-214 callers against a median near 60). It still
// labels reads but holds no join. Over chr12/18/20: +72 correct / -69
// discordant arm reads (all reads +333 / -421); 1.6 and 2.5 do less.
constexpr double kEmJoinMaxDepthRatio = 2.0;
// A haplotype duplicated at a site: callers above kEmDupDepthRatio x the
// chunk's median and the minority allele under kEmDupMinorFraction of them
// (one parent's two copies differ, the other parent matches one copy). The
// site separates copies, not parents, and holds no join. chr19:8.274-8.284 Mb,
// 1.5x depth, sites of error 0.01 splitting maternal + half the paternal reads
// from the other paternal half, carried a join that flipped 29 kb of reads.
// chr19 +83 correct / -83 discordant, one arm switch fewer; chr12 -1 / +1;
// chr20, chr18, chr17 unchanged. At 1.2 and 0.42: +1 switch on chr17.
constexpr double kEmDupDepthRatio = 1.3;
constexpr double kEmDupMinorFraction = 0.4;
// A site whose callers' mean MAPQ is below this sits where reads cannot be
// placed confidently -- again mostly copies' differences: on chr12, 433 of
// the 450 noise sites (parents not split) that labelled discordant reads,
// and on chr18 all 26, against 0.3% of real hets. It labels reads but holds
// no join, like a deep site. Over chr12/18/20: +58 correct / -61 discordant
// arm reads, switches -2, flips -2, all reads +886 / -634 (30-50 alike;
// muting its labels instead cost 1,832 correct reads).
constexpr double kEmJoinMinMeanMapq = 40.0;
constexpr double kEmMaxError = 0.45;
constexpr double kEmCleanStartError = 0.05;
constexpr double kEmOtherStartError = 0.2;
constexpr double kEmUnreliableError = 0.3;
// The learned error up to which an unverified call (unresolved_noisy) is
// written phased. Such calls label no reads, so this decides only the VCF:
// 0.1 against 0.3 cut arm switches 12 / 24 / 10 -> 10 / 20 / 6 on chr12,
// chr18 and chr20 and kept 60-75% of the truth heterozygotes they add.
constexpr double kEmUnverifiedReportError = 0.1;
// A phased site written unphased in the VCF above these learned errors (it
// still labels reads and marks them for the rescue). A real SNP is rarely
// read wrong: over chr12/18/20 the 324 truth-het SNPs phased at error
// 0.1-0.3 held 10 of the 61 arm flips (3%, against 0.007% below 0.1), and
// indels at 0.25-0.3 another 7 flips and 2 switches in 1,117 (0.8%).
// Written unphased: switches 24 -> 17, flips 61 -> 36, phased truth hets
// -1,660 (0.6%), reads unchanged.
constexpr double kEmVcfMaxSnpError = 0.1;
// How far an injected call the EM phased may lie from an unphased one it phases through the alignment solve's
// haplotypes (phase_by_alignment_solve).
constexpr hts_pos_t kEmAlignmentAnchorSpan = 5000;
constexpr double kEmVcfMaxIndelError = 0.25;
constexpr int kEmIterations = 15;
constexpr int kEmSwitchRounds = 30;
constexpr double kEmReadLabelPosterior = 0.8;
// A read below kEmReadLabelPosterior whose calls in a block all agree takes
// that block's label from this posterior, applied last and fill-only (with the
// island labels). Typically one call at an injected repeat site of learned
// error 0.2-0.3, which can never reach 0.8 alone: such calls on unphased reads
// were right 301 times in 315 (chr12/18/20). Applied inside the EM, the same
// labels steered windows, stitching and rescue: +261 correct / +101
// discordant, chr20 blocks flipping. Deferred: +369 / +4. 0.6 adds 8 / 5.
constexpr double kEmDeferredLabelPosterior = 0.70;
// A read with no deferred label either takes one from its unverified calls
// (noisy regions too long for MSA; sites phased for the VCF, learned error
// under kEmUnverifiedReportError) at this posterior. They never label a read
// that has any other way (see join_excluded); as the last resort, measured
// on chr12/18/20 on top of the echo sites: +22 correct / -3 discordant / -19
// unphased (chr18 alone +19 / -31 unphased).
constexpr double kEmUnverifiedLabelPosterior = 0.8;
constexpr double kEmBlockCutLogLikelihood = 4.0;
// The cut for joins carried only by weak sites (see bridge_weak_sites).
constexpr double kEmWeakBridgeCutLogLikelihood = 10.0;
constexpr int kEmLearnMinMapq = 20;
// Pseudo-reads with which each site's learned error is shrunk toward the
// pooled error of its class in the chunk (empirical Bayes; see m_step).
constexpr double kEmClassPriorReads = 5.0;
// The class pool is the median measured error of the class's sites that have
// at least kEmPoolMinWeight reads' worth of teaching weight, when there are
// at least kEmPoolMinSites of them (otherwise the weighted mean).
constexpr double kEmPoolMinWeight = 10.0;
constexpr size_t kEmPoolMinSites = 5;
// Reads that must span a block boundary, with evidence on both sides, for the
// blocks to join: one molecule can be chimeric, mis-placed or mis-called, and
// its likelihood alone can clear the cut (measured on chr18 at 38.34 Mb: one
// 26 kb read, one call per side, joined 428 hets in the wrong orientation).
constexpr int kEmMinJoinReads = 2;
// Reads that must carry a join's strong margin (its evidence at reliable,
// non-weak sites on both sides) for the join to count as strong; with fewer
// it is a weak bridge and needs kEmWeakBridgeCutLogLikelihood. Measured at
// chr12:118.48-118.50 Mb: a margin of 5.2 without weak sites, carried by two
// reads with a wrong call each, joined 1,114 reads in the wrong orientation
// through a homopolymer call of learned error 0.27. Three are still too few:
// at chr17:11.54-11.56 Mb three reads, their only calls left of a 16 kb gap
// in a noisy MSA region and partly wrong, carried a margin of 7.1 and joined
// 121 confidently labelled reads in the wrong orientation (HiPhase, WhatsHap
// start a block there). Four: chr17 +144 correct / -144 discordant, one
// switch fewer; chr20 -4 / +4, chr19 -1 / +1; chr12, chr18 unchanged.
constexpr int kEmMinStrongJoinReads = 4;
// A site whose less observed allele carries under this fraction of its reads
// does not take part: a skewed split is the signature of a paralog or error
// call, and in a homozygous stretch such a site labels reads at random.
constexpr double kEmMinMinorAlleleFraction = 0.25;
// An MSA-genotyped het such a split keeps out is written unphased when the MSA's reads are at least this fraction
// of the chunk's median site depth (write_skewed_msa_hets). HG002 HiFi chr12/18/20, plegma: +152 true het regions,
// -156 missed, -76 wrong, +28 false; reads unchanged. Without the depth rule +154 true for +83 false; a floor on
// the minor allele's share (5-20%) on top of it changes nothing.
constexpr double kEmMinMsaReadsFraction = 0.5;
// Echo sites: a site within kEmEchoMaxGap bp after another whose shared reads
// (at least kEmEchoMinShare of the larger site's callers) split the same way
// (kEmEchoMinConcordance of them) repeats that site's evidence. One event
// written as several records -- an MSA splitting reads by how their aligner
// placed a shift in a homopolymer (chr18:56,758,778-781, four injected sites,
// a 1|1 change in the truth), catalog SNPs inside a tandem repeat that follow
// where an insertion was placed (chr12:55,857,714 and 732) -- otherwise votes
// once per record, and two or four such votes outweigh a read's real hets
// whatever the parents say. An echo site is learned and phased but does not
// vote on read labels; block joins and switch moves still count it (muting it
// there too flipped blocks: chr12 -19, chr20 -24 correct). Measured on
// chr12/18/20 against 20293c8: +81 correct, -76 discordant, -5 unphased.
constexpr hts_pos_t kEmEchoMaxGap = 30;
constexpr double kEmEchoMinConcordance = 0.95;
constexpr double kEmEchoMinShare = 0.8;
constexpr int kEmEchoMinReads = 4;
// Island blocks (see GlobalEm::island_blocks). Every site the EM
// rejects offers islands, whatever its learned error (it is capped at
// kEmMaxError): a cap at 0.40 left +155 correct / +16 discordant / -171
// unphased arms reads over chr12/18/20 (all reads +2,249 / +1,476 / -3,725),
// 0.35 lost -464 / -104 / +568.
// A candidate outside the solve offers islands below this raw disagreement
// with the labelled reads.
constexpr double kIslandOutsideMaxError = 0.30;
// Unlabelled reads that must call a candidate for it to hold an island.
constexpr int kIslandMinReads = 2;

static double log_sum_exp(double a, double b) {
    const double m = std::max(a, b);
    return m + std::log(std::exp(a - m) + std::exp(b - m));
}

static void set_site_phase(CandidateVariant& c, hts_pos_t ps, int hap1_allele, int hap2_allele) {
    c.phase_set = ps;
    c.hap_to_cons_alle[1] = hap1_allele;
    c.hap_to_cons_alle[2] = hap2_allele;
    const bool h1 = hap1_allele != 0, h2 = hap2_allele != 0;
    c.hap_alt = h1 && h2 ? 3 : h1 ? 1 : h2 ? 2 : 0;
    c.hap_ref = h1 && h2 ? 0 : h1 ? 2 : h2 ? 1 : 0;
}

// Probability that a read is on haplotype 1 given its log-likelihoods on the
// two haplotypes (E-step form) and given their difference (label form).
static double read_posterior(double l1, double l2) {
    const double d = l2 - l1;
    return d > 700 ? 0.0 : d < -700 ? 1.0 : 1.0 / (1.0 + std::exp(d));
}

static double hap1_posterior(double margin) {
    return margin > 700 ? 1.0 : margin < -700 ? 0.0 : 1.0 / (1.0 + std::exp(-margin));
}

namespace {

// A voting site: a chunk candidate (locus < 0) contrasting alleles x and y, or
// a locus window contrasting its two haplotype sequences (0 and 1).
struct EmSite { long locus; size_t ci; int x; int y; hts_pos_t pos; };
// A read's call at site s: side b (0: x, 1: y) and its own error e.
struct EmObs { int s; int b; double e; };

// The global EM over one chunk's read x site matrix: site selection, the EM
// with switch moves, block cuts, read labels and the last-resort labels
// (deferred, unverified, one-edit windows, islands), applied in that order by
// run().
// A het whose orientation the MSA or alignment-solve passes below can read and set: REF and ALT (alleles 0 and 1), or
// a 1/2 site's two ALTs. Its "first" allele is allele 1 either way; its orientation, whether haplotype 1 carries it.
static bool phaseable_het(const std::array<int, 3>& h) {
    const int lo = std::min(h[1], h[2]), hi = std::max(h[1], h[2]);
    return (lo == 0 && hi == 1) || (lo == 1 && hi == 2);
}

static bool first_allele_on_hap1(const CandidateVariant& c) { return c.hap_to_cons_alle[1] == 1; }

static void set_first_allele_side(CandidateVariant& c, int on_hap1) {
    const int other = std::max(c.hap_to_cons_alle[1], c.hap_to_cons_alle[2]) == 2 ? 2 : 0;
    c.hap_to_cons_alle[1] = on_hap1 ? 1 : other;
    c.hap_to_cons_alle[2] = on_hap1 ? other : 1;
}

class GlobalEm {
public:
    GlobalEm(PhasingChunk& chunk, const std::vector<LocusWindowSite>* loci, const std::vector<char>* bridge_weak,
             const ObservationErrors* errors)
        : chunk_(chunk), loci_(loci), bridge_weak_(bridge_weak), errors_(errors),
          n_locus_(loci != nullptr ? loci->size() : 0) {}

    size_t run() {
        select_sites();
        n_ = sites_.size();
        if (n_ < 2) return 0;
        collect_observations();
        find_echo_sites();
        setup_joins();
        solve();
        cut_blocks();
        write_site_phases();
        label_reads();
        if (loci_ != nullptr) island_blocks();
        write_diagnostics();
        return blocks_;
    }

private:
    PhasingChunk& chunk_;
    const std::vector<LocusWindowSite>* loci_;
    const std::vector<char>* bridge_weak_;
    const ObservationErrors* errors_;
    const size_t n_locus_;

    std::vector<std::map<int, int>> allele_counts_;  // per candidate: allele -> reads
    std::vector<int> undecided_;                      // reads in range without an allele (diagnostic)
    std::vector<EmSite> sites_;                       // in position order
    size_t n_ = 0;
    std::vector<int> site_of_;     // candidate -> site, or -1
    std::vector<int> locus_site_;  // locus window -> site, or -1 (island-only)
    std::vector<int> phase_;       // 0: hap1 carries x
    std::vector<double> err_;      // learned error
    std::vector<std::array<double, 2>> site_log_;  // log(1 - e), log(e) of a call without its own error
    std::vector<std::vector<EmObs>> robs_;  // per read, in site order
    std::vector<double> post_;              // P(read on hap1)
    std::vector<int> start_label_;          // label the read started with (diagnostic)
    std::vector<char> learns_;              // read teaches the model (MAPQ)
    std::vector<char> echo_;
    // Holds no join: called by > kEmJoinMaxDepthRatio x the chunk's median
    // callers, or by reads of mean MAPQ < kEmJoinMinMeanMapq.
    std::vector<char> join_deep_;
    std::vector<char> site_weak_, join_weak_, join_excluded_;
    std::vector<std::vector<char>> join_skip_;
    std::vector<int> site_class_;
    bool m_started_ = false;
    std::vector<double> delta_, delta_strong_;
    std::vector<int> spanning_, spanning_strong_;
    std::vector<hts_pos_t> block_;
    std::vector<char> weak_bridge_;
    size_t blocks_ = 1;

    // Sites in position order with their two-allele contrast; a locus window
    // contrasts its two local haplotype sequences (0 and 1).
    void select_sites() {
        allele_counts_.assign(chunk_.candidates.size(), {});
        undecided_.assign(chunk_.candidates.size(), 0);
        for (size_t ri = 0; ri < chunk_.read_var_profile.size(); ++ri) {
            const ReadVariantProfile& prof = chunk_.read_var_profile[ri];
            if (prof.start_var_idx < 0 || chunk_.reads[ri].is_skipped) continue;
            for (size_t k = 0; k < prof.alleles.size(); ++k) {
                if (prof.alleles[k] >= 0)
                    ++allele_counts_[static_cast<size_t>(prof.start_var_idx) + k][prof.alleles[k]];
                else
                    ++undecided_[static_cast<size_t>(prof.start_var_idx) + k];
            }
        }
        for (size_t ci = 0; ci < chunk_.candidates.size(); ++ci) {
            const CandidateVariant& c = chunk_.candidates[ci];
            if ((c.lcd_var_i_to_cate & kCandGermlineVarCate) == 0) continue;
            int x = c.hap_to_cons_alle[1], y = c.hap_to_cons_alle[2];
            // An MSA-genotyped site of two insertion alleles (a 1/2): its alleles are the MSA's, not the two its reads
            // call most -- reads call one letter more in a long run as the reference allele as often as not
            // (chr12:20,069,335, T+1 and T+15: written as the T+15 het alone). HG002 HiFi chr12/18/20, plegma: true
            // het regions +39, wrong -35, missed -4; arm reads +38 correct, -37 discordant.
            if (c.msa_genotyped && c.msa_insertion_alts.size() == 2) x = 1, y = 2;
            if (x < 0 || y < 0 || x == y) {
                std::vector<std::pair<int, int>> by_count;
                for (const auto& [a, n] : allele_counts_[ci]) by_count.emplace_back(n, a);
                std::sort(by_count.rbegin(), by_count.rend());
                if (by_count.size() < 2) {
                    pgphase_trace::candidate("em", "not_site", c, "one_allele_called",
                                             pgphase_trace::genotype_text(c));
                    continue;
                }
                x = by_count[0].second;
                y = by_count[1].second;
            }
            const auto count_of = [&](int a) {
                const auto it = allele_counts_[ci].find(a);
                return it == allele_counts_[ci].end() ? 0 : it->second;
            };
            const int nx = count_of(x), ny = count_of(y);
            if (nx + ny == 0 || std::min(nx, ny) < kEmMinMinorAlleleFraction * (nx + ny)) {
                if (pgphase_trace::on())
                    pgphase_trace::candidate("em", "not_site", c, "minor_allele_below_fraction",
                                             "alleles " + std::to_string(x) + ":" + std::to_string(nx) + " " +
                                             std::to_string(y) + ":" + std::to_string(ny));
                continue;
            }
            sites_.push_back(EmSite{-1, ci, x, y, c.key.sort_pos()});
        }
        for (size_t li = 0; li < n_locus_; ++li)
            if (!(*loci_)[li].island_only) sites_.push_back(EmSite{static_cast<long>(li), 0, 0, 1, (*loci_)[li].pos});
        std::stable_sort(sites_.begin(), sites_.end(),
                         [](const EmSite& a, const EmSite& b) { return a.pos < b.pos; });
    }

    // Read observations as (site, side, own error) in site order; each read's
    // starting posterior from its incoming label.
    void collect_observations() {
        site_of_.assign(chunk_.candidates.size(), -1);
        locus_site_.assign(n_locus_, -1);
        for (size_t k = 0; k < n_; ++k) {
            if (sites_[k].locus >= 0) locus_site_[static_cast<size_t>(sites_[k].locus)] = static_cast<int>(k);
            else site_of_[sites_[k].ci] = static_cast<int>(k);
        }
        phase_.assign(n_, 0);
        err_.assign(n_, 0.0);
        for (size_t k = 0; k < n_; ++k) {
            const bool clean = sites_[k].locus < 0 &&
                (chunk_.candidates[sites_[k].ci].lcd_var_i_to_cate & (kCandCleanHetSnp | kCandCleanHetIndel)) != 0;
            err_[k] = clean ? kEmCleanStartError : kEmOtherStartError;
        }
        refresh_site_logs();
        robs_.assign(chunk_.reads.size(), {});
        for (size_t ri = 0; ri < chunk_.read_var_profile.size(); ++ri) {
            const ReadVariantProfile& prof = chunk_.read_var_profile[ri];
            if (prof.start_var_idx < 0 || chunk_.reads[ri].is_skipped) continue;
            for (size_t k = 0; k < prof.alleles.size(); ++k) {
                const int s = site_of_[static_cast<size_t>(prof.start_var_idx) + k];
                if (s < 0 || prof.alleles[k] < 0) continue;
                const EmSite& site = sites_[static_cast<size_t>(s)];
                double e = 0.0;
                if (errors_ != nullptr) {
                    const auto it = errors_->find(observation_key(ri, static_cast<size_t>(prof.start_var_idx) + k));
                    if (it != errors_->end()) e = it->second;
                }
                if (prof.alleles[k] == site.x) robs_[ri].push_back(EmObs{s, 0, e});
                else if (prof.alleles[k] == site.y) robs_[ri].push_back(EmObs{s, 1, e});
            }
        }
        for (size_t li = 0; li < n_locus_; ++li) {
            if (locus_site_[li] < 0) continue;  // island-only window
            for (const auto& [ri, side, e] : (*loci_)[li].observations)
                if (ri < robs_.size()) robs_[ri].push_back(EmObs{locus_site_[li], side, static_cast<double>(e)});
        }
        post_.assign(chunk_.reads.size(), 0.5);
        start_label_.assign(chunk_.reads.size(), 0);
        learns_.assign(chunk_.reads.size(), 0);
        for (size_t ri = 0; ri < robs_.size(); ++ri) {
            std::sort(robs_[ri].begin(), robs_[ri].end(),
                      [](const EmObs& a, const EmObs& b) { return a.s < b.s; });
            const int hap = ri < chunk_.haps.size() ? chunk_.haps[ri] : 0;
            post_[ri] = hap == 1 ? 0.95 : hap == 2 ? 0.05 : 0.5;
            if (ri < start_label_.size()) start_label_[ri] = hap;
            learns_[ri] = chunk_.reads[ri].mapq >= kEmLearnMinMapq;
        }
    }

    // Echo sites (see kEmEchoMaxGap): each candidate pair within the gap,
    // counted over the reads that call both.
    void find_echo_sites() {
        echo_.assign(n_, 0);
        std::vector<int> callers(n_, 0);
        std::map<std::pair<int, int>, std::array<int, 2>> pair_calls;  // (k, j): same / opposite sides
        for (const std::vector<EmObs>& obs : robs_) {
            for (size_t i = 0; i < obs.size(); ++i) {
                ++callers[static_cast<size_t>(obs[i].s)];
                const EmSite& si = sites_[static_cast<size_t>(obs[i].s)];
                if (si.locus >= 0) continue;
                for (size_t j = i + 1; j < obs.size(); ++j) {
                    const EmSite& sj = sites_[static_cast<size_t>(obs[j].s)];
                    if (sj.pos - si.pos > kEmEchoMaxGap) break;
                    if (obs[j].s == obs[i].s || sj.locus >= 0) continue;
                    ++pair_calls[{obs[i].s, obs[j].s}][obs[i].b == obs[j].b ? 0 : 1];
                }
            }
        }
        for (const auto& [kj, c] : pair_calls) {
            const int shared = c[0] + c[1];
            const int larger = std::max(callers[static_cast<size_t>(kj.first)], callers[static_cast<size_t>(kj.second)]);
            if (shared >= kEmEchoMinReads && shared >= kEmEchoMinShare * larger &&
                std::max(c[0], c[1]) >= kEmEchoMinConcordance * shared)
                echo_[static_cast<size_t>(kj.second)] = 1;
        }
    }

    // Which sites and calls the block-cut test counts, and each site's class
    // for the error prior.
    void setup_joins() {
        // Weak sites (see bridge_weak_sites) are phased but cannot hold blocks
        // together: the final cuts are judged on the other sites' evidence.
        site_weak_.assign(n_, 0);
        if (bridge_weak_ != nullptr)
            for (size_t k = 0; k < n_; ++k)
                if (sites_[k].locus < 0 && sites_[k].ci < bridge_weak_->size())
                    site_weak_[k] = (*bridge_weak_)[sites_[k].ci];
        // The join test counts a locus once. A locus window and a candidate
        // inside it are read from the same bases, so a read's two calls there
        // agree by construction; counted twice they outvote a real heterozygote
        // read once. And a window over a homopolymer length call is that call
        // re-read, so it is weak too. At chr12:122,284,562 a window and the +1
        // indel inside it (one locus, purity 0.53 against the parents) carried a
        // 30 kb join past a true het (122,290,121, purity 0.95) with a margin of
        // 5.65 from "strong" sites, and 165 reads went discordant. Labels and
        // learning keep every observation: only the decision to join uses this.
        join_weak_ = site_weak_;
        join_skip_.assign(robs_.size(), {});
        // Calls kept from noisy regions too long for MSA (unresolved_noisy) are
        // left out of the join test altogether: a divergent copy's mismatches
        // (695 such calls in 24 kb at chr20:49.65 Mb) clear even the weak-bridge
        // margin by sheer number, and flipped a block of 110 reads. Clean hets no
        // confidently mapped read carries (low_mapq_only) are treated the same:
        // as ordinary sites they added two arm switches over chr12/18/20 (chr12
        // 34.64 Mb, chr20 32.25 Mb) for three correct reads.
        join_excluded_.assign(n_, 0);
        for (size_t k = 0; k < n_; ++k)
            if (sites_[k].locus < 0 && (chunk_.candidates[sites_[k].ci].unresolved_noisy ||
                                        chunk_.candidates[sites_[k].ci].low_mapq_only)) join_excluded_[k] = 1;
        // Sites in collapsed duplications or poorly mapped reads hold no join
        // (kEmJoinMaxDepthRatio, kEmJoinMinMeanMapq).
        join_deep_.assign(n_, 0);
        {
            std::vector<int> callers(n_, 0);
            std::vector<double> mapq_sum(n_, 0.0);
            for (size_t ri = 0; ri < robs_.size(); ++ri) {
                const std::vector<EmObs>& obs = robs_[ri];
                for (size_t i = 0; i < obs.size(); ++i) {
                    if (i > 0 && obs[i].s == obs[i - 1].s) continue;
                    ++callers[static_cast<size_t>(obs[i].s)];
                    mapq_sum[static_cast<size_t>(obs[i].s)] += chunk_.reads[ri].mapq;
                }
            }
            for (size_t k = 0; k < n_; ++k)
                if (callers[k] > 0 && mapq_sum[k] / callers[k] < kEmJoinMinMeanMapq) join_deep_[k] = 1;
            std::vector<int> called;
            for (const int c : callers) if (c > 0) called.push_back(c);
            if (!called.empty()) {
                std::nth_element(called.begin(), called.begin() + static_cast<long>(called.size() / 2), called.end());
                const double median = called[called.size() / 2];
                for (size_t k = 0; k < n_; ++k)
                    if (callers[k] > kEmJoinMaxDepthRatio * median) join_deep_[k] = 1;
                std::vector<std::array<int, 2>> sides(n_, {0, 0});
                for (const std::vector<EmObs>& obs : robs_)
                    for (size_t i = 0; i < obs.size(); ++i) {
                        if (i > 0 && obs[i].s == obs[i - 1].s) continue;
                        ++sides[static_cast<size_t>(obs[i].s)][obs[i].b ? 1 : 0];
                    }
                for (size_t k = 0; k < n_; ++k) {
                    const int tot = sides[k][0] + sides[k][1];
                    if (tot > 0 && callers[k] > kEmDupDepthRatio * median &&
                        std::min(sides[k][0], sides[k][1]) < kEmDupMinorFraction * tot) join_deep_[k] = 1;
                }
            }
        }
        std::vector<hts_pos_t> site_pos(n_);
        for (size_t k = 0; k < n_; ++k) site_pos[k] = sites_[k].pos;
        for (size_t k = 0; k < n_; ++k) {
            if (sites_[k].locus < 0) continue;
            const LocusWindowSite& w = (*loci_)[static_cast<size_t>(sites_[k].locus)];
            const auto lo = std::lower_bound(site_pos.begin(), site_pos.end(), w.pos);
            const auto hi = std::upper_bound(lo, site_pos.end(), std::max(w.pos, w.end));
            for (auto it = lo; it != hi; ++it) {
                const size_t j = static_cast<size_t>(it - site_pos.begin());
                if (sites_[j].locus < 0 && site_weak_[j]) join_weak_[k] = 1;
            }
        }
        std::vector<std::pair<hts_pos_t, hts_pos_t>> spans;
        for (size_t ri = 0; ri < robs_.size(); ++ri) {
            spans.clear();
            for (const EmObs& o : robs_[ri]) {
                const EmSite& site = sites_[static_cast<size_t>(o.s)];
                if (site.locus < 0) continue;
                const LocusWindowSite& w = (*loci_)[static_cast<size_t>(site.locus)];
                spans.emplace_back(w.pos, std::max(w.pos, w.end));
            }
            if (spans.empty()) continue;
            join_skip_[ri].assign(robs_[ri].size(), 0);
            for (size_t i = 0; i < robs_[ri].size(); ++i) {
                const EmSite& site = sites_[static_cast<size_t>(robs_[ri][i].s)];
                if (site.locus >= 0) continue;
                for (const auto& [beg, end] : spans)
                    if (beg <= site.pos && site.pos <= end) join_skip_[ri][i] = 1;
            }
        }
        // Site classes for the error prior: clean catalog het, other catalog
        // site, alignment-injected site, locus window -- each split by weak --
        // and calls from noisy regions too long for MSA, pooled apart: they are
        // flagged weak, and sharing the homopolymer calls' pool raised those
        // calls' errors with their own (chr12 lost 232 labels at 37 and 133 Mb).
        site_class_.assign(n_, 0);
        for (size_t k = 0; k < n_; ++k) {
            int kind = 3;
            if (sites_[k].locus < 0) {
                const CandidateVariant& c = chunk_.candidates[sites_[k].ci];
                kind = c.unresolved_noisy ? 4
                     : c.bam_injected ? 2
                     : (c.lcd_var_i_to_cate & (kCandCleanHetSnp | kCandCleanHetIndel)) != 0 ? 0 : 1;
            }
            site_class_[k] = 2 * kind + (site_weak_[k] ? 1 : 0);
        }
    }

    double obs_error(const EmObs& o) const {
        return std::min(kEmMaxError, 1.0 - (1.0 - err_[static_cast<size_t>(o.s)]) * (1.0 - o.e));
    }

    // Most calls carry no error of their own: their two log terms depend on
    // the site's error alone and are computed once per error update (the
    // E-steps, switch moves and cut tests took ~4% of a run in std::log).
    void refresh_site_logs() {
        site_log_.resize(n_);
        for (size_t k = 0; k < n_; ++k) {
            const double e = std::min(kEmMaxError, 1.0 - (1.0 - err_[k]) * (1.0 - 0.0));
            site_log_[k] = {std::log(1.0 - e), std::log(e)};
        }
    }

    std::array<double, 2> read_log_likelihoods(const EmObs& o) const {
        const bool match1 = o.b == phase_[static_cast<size_t>(o.s)];
        if (o.e == 0.0) {
            const auto& l = site_log_[static_cast<size_t>(o.s)];
            return {match1 ? l[0] : l[1], match1 ? l[1] : l[0]};
        }
        const double e = obs_error(o);
        return {std::log(match1 ? 1.0 - e : e), std::log(match1 ? e : 1.0 - e)};
    }

    // A call's vote on its read's label: none at an echo site.
    std::array<double, 2> label_log_likelihoods(const EmObs& o) const {
        if (echo_[static_cast<size_t>(o.s)]) return {0.0, 0.0};
        return read_log_likelihoods(o);
    }

    // A site's error is its reads' disagreement rate, shrunk toward the pooled
    // rate of its class by kEmClassPriorReads pseudo-reads. Two sites that
    // share their reads and disagree are otherwise symmetric: whichever the
    // starting labels follow reads error 0.01 and the other ~0.3, so a
    // homopolymer length call could reject a clean SNP (chr12:6,134,662, 35
    // paternal and 24 maternal reads all split by the SNP, rejected at 0.316
    // against a +T/+TT call that does not separate the parents). The class
    // pools say how often each kind of site is wrong in this chunk.
    void m_step() {
        std::vector<std::array<double, 2>> w(n_, {0.0, 0.0});  // hap1 weight per side
        for (size_t ri = 0; ri < robs_.size(); ++ri) {
            if (!learns_[ri] || robs_[ri].empty()) continue;
            // A read whose observations all sit at one site has a posterior
            // made from that site alone: counted, it returns the site's own
            // current error and pins it wherever it started. Weakly connected
            // injected hets stayed rejected so -- chr12:105,144,892, 82 reads
            // at purity 0.98 against the parents, error 0.38 -- and their
            // reads unphased. Such a read is still labelled and still teaches
            // a clean catalog het, whose start (kEmCleanStartError) is the
            // better prior: a lone clean SNP taught only by reads crossing
            // noisy neighbours was rejected (chr18:71,515,222, 74 reads at
            // purity 0.99, 21 reads went discordant).
            const bool single_site = robs_[ri].front().s == robs_[ri].back().s;
            if (single_site && site_class_[static_cast<size_t>(robs_[ri].front().s)] / 2 != 0) continue;
            for (const EmObs& o : robs_[ri]) {
                w[static_cast<size_t>(o.s)][static_cast<size_t>(o.b)] += post_[ri];
                w[static_cast<size_t>(o.s)][static_cast<size_t>(1 - o.b)] += 1.0 - post_[ri];
            }
        }
        std::array<double, 10> class_minor{}, class_total{};
        // A class's pool is the median error of its sites whose error is
        // measured (>= kEmPoolMinWeight reads' worth of teaching weight), not
        // the weighted mean: noise sites at ~0.4 dragged the mean up for the
        // class's real hets, typically at 0.02-0.05, and a real het whose own
        // error is not pinned by its reads was shrunk there (chr12:79,657,521,
        // purity 0.99 against the parents, error 0.20, 41 reads just under
        // the label posterior). A site with no teaching weight at all takes
        // the pool instead of keeping its start error.
        std::array<std::vector<double>, 10> class_rates;
        for (size_t k = 0; k < n_; ++k) {
            class_minor[static_cast<size_t>(site_class_[k])] += std::min(w[k][0], w[k][1]);
            class_total[static_cast<size_t>(site_class_[k])] += w[k][0] + w[k][1];
            if (w[k][0] + w[k][1] >= kEmPoolMinWeight)
                class_rates[static_cast<size_t>(site_class_[k])].push_back(std::min(w[k][0], w[k][1]) / (w[k][0] + w[k][1]));
        }
        std::array<double, 10> class_pool{};
        for (size_t c = 0; c < 10; ++c) {
            class_pool[c] = class_total[c] > 0.0 ? class_minor[c] / class_total[c] : kEmOtherStartError;
            if (class_rates[c].size() >= kEmPoolMinSites) {
                std::vector<double>& r = class_rates[c];
                std::nth_element(r.begin(), r.begin() + static_cast<long>(r.size() / 2), r.end());
                class_pool[c] = r[r.size() / 2];
            }
        }
        const double shrink = m_started_ ? kEmClassPriorReads : 0.0;
        for (size_t k = 0; k < n_; ++k) {
            const double total = w[k][0] + w[k][1];
            if (total < 1e-6) {
                if (m_started_)
                    err_[k] = std::min(kEmMaxError, std::max(kEmMinError, class_pool[static_cast<size_t>(site_class_[k])]));
                continue;
            }
            phase_[k] = w[k][0] >= w[k][1] ? 0 : 1;
            const size_t c = static_cast<size_t>(site_class_[k]);
            const double pooled = class_pool[c];
            const double e = (std::min(w[k][0], w[k][1]) + shrink * pooled) / (total + shrink);
            err_[k] = std::min(kEmMaxError, std::max(kEmMinError, e));
        }
        refresh_site_logs();
        m_started_ = true;
    }

    // Posteriors of reads; with only_unlabelled, of the reads that started
    // without a label.
    void e_step(bool only_unlabelled = false) {
        for (size_t ri = 0; ri < robs_.size(); ++ri) {
            if (robs_[ri].empty() || (only_unlabelled && start_label_[ri] != 0)) continue;
            double l1 = 0.0, l2 = 0.0;
            for (const EmObs& o : robs_[ri]) {
                const auto ll = label_log_likelihoods(o);
                l1 += ll[0];
                l2 += ll[1];
            }
            post_[ri] = read_posterior(l1, l2);
        }
    }

    // delta[k] = log L(current) - log L(every site right of k flipped), over
    // the reads that teach the model or, for block cuts, over every read: a
    // read too ambiguously mapped to shape the model still shows the molecule
    // continues across the boundary.
    // joining: the block-cut test, with each locus counted once (join_skip,
    // join_weak) rather than the switch moves' every observation.
    std::vector<double> boundary_delta(bool every_read, bool strong_only = false, bool joining = false) const {
        static const std::vector<char> no_skip;
        std::vector<double> diff(n_ + 1, 0.0);
        std::vector<EmObs> strong;
        for (size_t ri = 0; ri < robs_.size(); ++ri) {
            const std::vector<EmObs>* items_ptr = &robs_[ri];
            const std::vector<char>& skip = joining ? join_skip_[ri] : no_skip;
            const std::vector<char>& weak_of = joining ? join_weak_ : site_weak_;
            if (strong_only || joining) {
                strong.clear();
                for (size_t i = 0; i < robs_[ri].size(); ++i) {
                    const EmObs& o = robs_[ri][i];
                    if (joining && ((!skip.empty() && skip[i]) || join_excluded_[static_cast<size_t>(o.s)] ||
                                    join_deep_[static_cast<size_t>(o.s)])) continue;
                    if (strong_only && weak_of[static_cast<size_t>(o.s)]) continue;
                    strong.push_back(o);
                }
                items_ptr = &strong;
            }
            const auto& items = *items_ptr;
            if ((!every_read && !learns_[ri]) || items.size() < 2 || items.front().s == items.back().s)
                continue;
            std::vector<std::array<double, 2>> ll(items.size());
            double tot1 = 0.0, tot2 = 0.0;
            for (size_t i = 0; i < items.size(); ++i) {
                ll[i] = read_log_likelihoods(items[i]);
                tot1 += ll[i][0];
                tot2 += ll[i][1];
            }
            double c1 = 0.0, c2 = 0.0;
            for (size_t i = 0; i + 1 < items.size(); ++i) {
                c1 += ll[i][0];
                c2 += ll[i][1];
                const int k = items[i].s, k_next = items[i + 1].s;
                if (k_next == k) continue;
                const double gain = log_sum_exp(tot1, tot2) -
                                    log_sum_exp(c1 + (tot2 - c2), c2 + (tot1 - c1));
                diff[static_cast<size_t>(k)] += gain;
                diff[static_cast<size_t>(k_next)] -= gain;
            }
        }
        std::vector<double> delta(n_, 0.0);
        double run = 0.0;
        for (size_t k = 0; k < n_; ++k) {
            run += diff[k];
            delta[k] = run;
        }
        return delta;
    }

    // EM, then switch moves: flip everything right of the worst boundary while
    // that raises the likelihood. Reads the clean stage left unlabelled start
    // from an E-step under the sites' initial orientation and starting errors,
    // not from 0.5. At 0.5 every site's first M-step error is 0.5 -- a saddle --
    // and the class prior holds an island there: two sites shared by 75 reads
    // at chr18:60,372,121, 97% pure against the parents, stayed at error 0.45
    // with no read labelled.
    void solve() {
        e_step(/*only_unlabelled=*/true);
        for (int round = 0; round < kEmSwitchRounds; ++round) {
            for (int it = 0; it < kEmIterations; ++it) {
                m_step();
                e_step();
            }
            const std::vector<double> delta = boundary_delta(false);
            size_t worst = n_;
            for (size_t k = 0; k + 1 < n_; ++k)
                if (delta[k] < -1e-6 && (worst == n_ || delta[k] < delta[worst])) worst = k;
            if (worst == n_) break;
            for (size_t k = worst + 1; k < n_; ++k) phase_[k] ^= 1;
            for (size_t ri = 0; ri < robs_.size(); ++ri)
                if (!robs_[ri].empty() && static_cast<size_t>(robs_[ri].front().s) > worst)
                    post_[ri] = 1.0 - post_[ri];
        }
    }

    // Blocks: cut where flipping the rest would cost less than the threshold,
    // or where too few reads cross.
    void cut_blocks() {
        delta_ = boundary_delta(true, false, true);
        // A join that the non-weak sites alone do not carry (its margin comes
        // from weak sites: homopolymer length calls) needs a larger margin.
        // Measured on chr20: wrong joins of this kind sit just above the
        // ordinary cut (median 6), right ones well above it (median 19).
        delta_strong_ = bridge_weak_ != nullptr ? boundary_delta(true, true, true) : delta_;
        // Reads whose calls the join test counts (join_skip, join_excluded,
        // join_deep as in delta_) lie on both sides of each boundary: a read whose calls on
        // one side are all left out of the margin cannot stand for a molecule
        // crossing it.
        spanning_.assign(n_, 0);
        {
            std::vector<int> diff(n_ + 1, 0);
            for (size_t ri = 0; ri < robs_.size(); ++ri) {
                const std::vector<EmObs>& obs = robs_[ri];
                int first = -1, last = -1;
                for (size_t i = 0; i < obs.size(); ++i) {
                    if ((!join_skip_[ri].empty() && join_skip_[ri][i]) || join_excluded_[static_cast<size_t>(obs[i].s)] ||
                        join_deep_[static_cast<size_t>(obs[i].s)])
                        continue;
                    if (first < 0) first = obs[i].s;
                    last = obs[i].s;
                }
                if (first < 0 || first == last) continue;
                ++diff[static_cast<size_t>(first)];
                --diff[static_cast<size_t>(last)];
            }
            int run = 0;
            for (size_t k = 0; k < n_; ++k) spanning_[k] = run += diff[k];
        }
        // Reads whose strong evidence -- reliable, join-eligible, non-weak
        // sites -- lies on both sides of each boundary. The strong margin is
        // only as independent as the reads behind it (see kEmMinStrongJoinReads).
        spanning_strong_.assign(n_, 0);
        if (bridge_weak_ != nullptr) {
            std::vector<int> diff(n_ + 1, 0);
            for (const std::vector<EmObs>& obs : robs_) {
                int first = -1, last = -1;
                for (const EmObs& o : obs) {
                    const size_t k = static_cast<size_t>(o.s);
                    if (err_[k] >= kEmUnreliableError || join_excluded_[k] || join_weak_[k] || join_deep_[k]) continue;
                    if (first < 0) first = o.s;
                    last = o.s;
                }
                if (first < 0 || first == last) continue;
                ++diff[static_cast<size_t>(first)];
                --diff[static_cast<size_t>(last)];
            }
            int run = 0;
            for (size_t k = 0; k < n_; ++k) spanning_strong_[k] = run += diff[k];
        }
        block_.assign(n_, 0);
        weak_bridge_.assign(n_, 0);
        for (size_t k = 0; k < n_; ++k)
            weak_bridge_[k] = delta_strong_[k] < kEmBlockCutLogLikelihood ||
                              (bridge_weak_ != nullptr && spanning_strong_[k] < kEmMinStrongJoinReads);
        // A site that alone connects its neighbours -- no read observes
        // join-eligible sites on both sides of it -- with weak joins on both
        // sides: the block would hinge on that one site's phase. Cut the
        // weaker of its two joins. chr6:152,636,704, an injected call whose
        // reads split the parents 0.58 (learned error 0.11), was the only
        // connection between 152,625,051 and 152,649,680 and flipped the block.
        // chr6 +930 correct / -927 discordant, one arm switch fewer; chr16
        // -1 / +1; chr12, chr17, chr18, chr19, chr20, chr21 unchanged.
        std::vector<char> hinge_cut(n_, 0);
        if (n_ >= 3) {
            std::vector<int> through(n_, 0);   // reads with join-eligible sites on both sides of k
            std::vector<int> diff(n_ + 1, 0);
            for (size_t ri = 0; ri < robs_.size(); ++ri) {
                const std::vector<EmObs>& obs = robs_[ri];
                int first = -1, last = -1;
                for (size_t i = 0; i < obs.size(); ++i) {
                    const size_t k = static_cast<size_t>(obs[i].s);
                    if ((!join_skip_[ri].empty() && join_skip_[ri][i]) || join_excluded_[k] || join_deep_[k]) continue;
                    if (first < 0) first = obs[i].s;
                    last = obs[i].s;
                }
                if (first < 0 || last - first < 2) continue;
                ++diff[static_cast<size_t>(first) + 1];
                --diff[static_cast<size_t>(last)];
            }
            int run = 0;
            for (size_t k = 0; k < n_; ++k) through[k] = run += diff[k];
            for (size_t k = 1; k + 1 < n_; ++k) {
                if (through[k] > 0 || !weak_bridge_[k - 1] || !weak_bridge_[k]) continue;
                hinge_cut[delta_[k - 1] <= delta_[k] ? k - 1 : k] = 1;
            }
        }
        hts_pos_t current = sites_[0].pos;
        blocks_ = 1;
        for (size_t k = 0; k < n_; ++k) {
            block_[k] = current;
            if (k + 1 < n_ && (delta_[k] < kEmBlockCutLogLikelihood || spanning_[k] < kEmMinJoinReads ||
                               (weak_bridge_[k] && delta_[k] < kEmWeakBridgeCutLogLikelihood) || hinge_cut[k])) {
                current = sites_[k + 1].pos;
                ++blocks_;
            }
        }
    }

    // Each site's phase and block onto its candidate (or locus window); an
    // unreliable site is left unphased.
    void write_site_phases() {
        for (CandidateVariant& c : chunk_.candidates)
            if ((c.lcd_var_i_to_cate & kCandGermlineVarCate) != 0) c.phase_set = kUnsetCandidatePhaseSet, c.vcf_unphased = false;
        for (size_t k = 0; k < n_; ++k) {
            if (sites_[k].locus >= 0) {
                const LocusWindowSite& w = (*loci_)[static_cast<size_t>(sites_[k].locus)];
                w.phase = phase_[k];
                w.block = block_[k];
                w.error = err_[k];
                continue;
            }
            CandidateVariant& c = chunk_.candidates[sites_[k].ci];
            if (err_[k] >= (join_excluded_[k] ? kEmUnverifiedReportError : kEmUnreliableError)) {
                c.hap_to_cons_alle[1] = c.hap_to_cons_alle[2] = -1;
                c.hap_alt = c.hap_ref = 0;
                // An MSA-genotyped het the reads cannot phase is still a het: its two haplotype consensus
                // sequences say so, the reads only fail to tell which is which (a long homopolymer that
                // HiFi reads one short as often as not: chr12:953,298, A24/A25, error 0.33). Written
                // unphased, with no phase set, so it labels, joins and rephases nothing. HG002 HiFi
                // chr12/18/20, plegma: +267 true het regions, -282 missed, +11 false; reads unchanged.
                if (!join_excluded_[k] && c.msa_genotyped && err_[k] < kEmMaxError) {
                    set_site_phase(c, kUnsetCandidatePhaseSet, sites_[k].x, sites_[k].y);
                    c.vcf_unphased = true;
                }
                if (pgphase_trace::on())
                    pgphase_trace::candidate("em", "site", c, c.vcf_unphased ? "msa_het_unphasable" : "unreliable",
                                             "err=" + std::to_string(err_[k]).substr(0, 5) +
                                             (join_excluded_[k] ? " join_excluded" : ""));
                continue;
            }
            const int hap1 = phase_[k] == 0 ? sites_[k].x : sites_[k].y;
            const int hap2 = phase_[k] == 0 ? sites_[k].y : sites_[k].x;
            set_site_phase(c, block_[k], hap1, hap2);
            // Unverified calls from noisy regions too long for MSA: inside VNTRs
            // they and the truth describe different copies of the repeat
            // (chr18:79.877 Mb, six switches over 32 correctly phased reads).
            c.vcf_unphased = c.unresolved_noisy ||
                             err_[k] >= (c.key.type == VariantType::Snp ? kEmVcfMaxSnpError : kEmVcfMaxIndelError);
            if (pgphase_trace::on())
                pgphase_trace::candidate("em", "site", c, c.vcf_unphased ? "phased_written_unphased" : "phased",
                                         "block=" + std::to_string(block_[k]) + " err=" +
                                         std::to_string(err_[k]).substr(0, 5) + " " + pgphase_trace::genotype_text(c) +
                                         (c.unresolved_noisy ? " unresolved_noisy" : "") + (join_deep_[k] ? " join_deep" : ""));
        }
        confirm_phase_by_msa();
        write_skewed_msa_hets();
        phase_unphasable_by_msa();
        phase_by_alignment_solve();
    }

    // What phase_unphasable_by_msa and confirm_phase_by_msa do within an MSA region, by the alignment solve's phase
    // instead, for an injected het (or a catalog row its complement call describes) still written unphased (the region's lone call: chr20:782,916, a homopolymer
    // insertion at error 0.32 with no other call in its region): the injected calls the EM phased within
    // kEmAlignmentAnchorSpan of it in the same alignment phase set relate to it as that solve's haplotypes do. Phased
    // when every such call implies one orientation in one block (and, if the EM phased it, the EM's own).
    void phase_by_alignment_solve() {
        const auto alt_on_hap1 = [](const CandidateVariant& c) { return first_allele_on_hap1(c); };
        const auto usable = [](const CandidateVariant& c) {
            const auto& h = c.hap_to_cons_alle;
            return c.aln_phase_set > 0 && c.aln_hap1_allele >= 0 && h.size() > 2 && phaseable_het(h);
        };
        std::vector<size_t> anchors;
        for (size_t k = 0; k < n_; ++k) {
            if (sites_[k].locus >= 0) continue;
            const CandidateVariant& d = chunk_.candidates[sites_[k].ci];
            if (usable(d) && d.phase_set > 0 && !d.vcf_unphased) anchors.push_back(sites_[k].ci);
        }
        if (anchors.empty()) return;
        for (CandidateVariant& c : chunk_.candidates) {
            if (!usable(c) || !c.vcf_unphased || c.unresolved_noisy) continue;
            hts_pos_t ps = 0;
            int implied = -1;
            bool consistent = true;
            for (const size_t di : anchors) {
                const CandidateVariant& d = chunk_.candidates[di];
                if (d.aln_phase_set != c.aln_phase_set || std::abs(d.key.pos - c.key.pos) > kEmAlignmentAnchorSpan)
                    continue;
                const bool cis = (c.aln_hap1_allele == 1) == (d.aln_hap1_allele == 1);
                const int here = cis == alt_on_hap1(d) ? 1 : 0;
                consistent = consistent && (ps == 0 || ps == d.phase_set) && (implied < 0 || implied == here);
                ps = d.phase_set;
                implied = here;
            }
            if (implied < 0 || !consistent) continue;
            if (c.phase_set > 0) {
                if (c.phase_set != ps || (alt_on_hap1(c) ? 1 : 0) != implied) continue;
            } else {
                set_first_allele_side(c, implied);
                c.phase_set = ps;
            }
            c.vcf_unphased = false;
            pgphase_trace::candidate("em", "site", c, "phased_by_alignment_solve", pgphase_trace::genotype_text(c));
        }
    }

    // An MSA het written unphased with no phase set -- the EM could not phase it (msa_het_unphasable), or its reads
    // split too unevenly to be a site (write_skewed_msa_hets) -- takes its phase from the MSA where its het region holds
    // sites the EM phased: the MSA's consensus sequences say whether it is cis or trans to each. Only when every such
    // site is in one block and implies one orientation.
    void phase_unphasable_by_msa() {
        const auto alt_on_hap1 = [](const CandidateVariant& c) { return first_allele_on_hap1(c); };
        const auto usable = [](const CandidateVariant& c) {
            const auto& h = c.hap_to_cons_alle;
            return c.msa_region >= 0 && c.msa_hap1_allele >= 0 && h.size() > 2 && phaseable_het(h);
        };
        std::vector<size_t> anchors;
        for (size_t k = 0; k < n_; ++k) {
            if (sites_[k].locus >= 0) continue;
            const CandidateVariant& d = chunk_.candidates[sites_[k].ci];
            if (usable(d) && d.phase_set > 0 && !d.vcf_unphased) anchors.push_back(sites_[k].ci);
        }
        if (anchors.empty()) return;
        for (CandidateVariant& c : chunk_.candidates) {
            if (!usable(c) || !c.vcf_unphased || c.phase_set > 0 || c.unresolved_noisy) continue;
            hts_pos_t ps = 0;
            int implied = -1;
            bool consistent = true;
            for (const size_t di : anchors) {
                const CandidateVariant& d = chunk_.candidates[di];
                if (d.msa_region != c.msa_region) continue;
                const bool msa_cis = (c.msa_hap1_allele == 1) == (d.msa_hap1_allele == 1);
                const int here = msa_cis == alt_on_hap1(d) ? 1 : 0;
                consistent = consistent && (ps == 0 || ps == d.phase_set) && (implied < 0 || implied == here);
                ps = d.phase_set;
                implied = here;
            }
            if (implied < 0 || !consistent) continue;
            set_first_allele_side(c, implied);
            c.phase_set = ps;
            c.vcf_unphased = false;
            pgphase_trace::candidate("em", "site", c, "phased_by_msa", pgphase_trace::genotype_text(c));
        }
    }

    // A site written unphased for its learned error (kEmVcfMaxSnpError, kEmVcfMaxIndelError) whose orientation
    // against every site of its het MSA region phased in its block is the MSA's own (cis when one consensus carries
    // both ALTs, trans otherwise) is written phased: the MSA's consensus sequences, built from the region's reads,
    // confirm what the EM learned less surely (chr20:15,508,905 T>A at error 0.14 beside 15,508,906 T>A at 0.01, one
    // consensus carrying both, the truth 1|0 1|0).
    void confirm_phase_by_msa() {
        std::vector<size_t> by_region;
        for (size_t k = 0; k < n_; ++k) {
            if (sites_[k].locus >= 0) continue;
            const CandidateVariant& c = chunk_.candidates[sites_[k].ci];
            if (c.msa_region >= 0 && c.msa_hap1_allele >= 0 && c.phase_set > 0 && phaseable_het(c.hap_to_cons_alle))
                by_region.push_back(sites_[k].ci);
        }
        const auto alt_on_hap1 = [](const CandidateVariant& c) { return first_allele_on_hap1(c); };
        if (pgphase_trace::on())
            for (const size_t ci : by_region) {
                const CandidateVariant& c = chunk_.candidates[ci];
                if (c.vcf_unphased || c.unresolved_noisy) continue;
                int agree = 0, disagree = 0;
                for (const size_t di : by_region) {
                    const CandidateVariant& d = chunk_.candidates[di];
                    if (di == ci || d.vcf_unphased || d.msa_region != c.msa_region || d.phase_set != c.phase_set) continue;
                    const bool em_cis = alt_on_hap1(c) == alt_on_hap1(d);
                    const bool msa_cis = (c.msa_hap1_allele == 1) == (d.msa_hap1_allele == 1);
                    ++(em_cis == msa_cis ? agree : disagree);
                }
                if (agree + disagree > 0)
                    pgphase_trace::candidate("em", "msa_check", c, disagree > 0 ? "disagrees" : "agrees",
                                             "agree=" + std::to_string(agree) + " disagree=" + std::to_string(disagree) +
                                             " err=" + std::to_string(err_[static_cast<size_t>(site_of_[ci])]).substr(0, 5) +
                                             " " + pgphase_trace::genotype_text(c));
            }
        for (const size_t ci : by_region) {
            CandidateVariant& c = chunk_.candidates[ci];
            if (!c.vcf_unphased || c.unresolved_noisy) continue;
            int agree = 0, disagree = 0;
            for (const size_t di : by_region) {
                const CandidateVariant& d = chunk_.candidates[di];
                if (di == ci || d.vcf_unphased || d.msa_region != c.msa_region || d.phase_set != c.phase_set) continue;
                const bool em_cis = alt_on_hap1(c) == alt_on_hap1(d);
                const bool msa_cis = (c.msa_hap1_allele == 1) == (d.msa_hap1_allele == 1);
                ++(em_cis == msa_cis ? agree : disagree);
            }
            if (agree == 0 || disagree > 0) continue;
            c.vcf_unphased = false;
            pgphase_trace::candidate("em", "site", c, "phased_msa_confirmed",
                                     "agree=" + std::to_string(agree) + " " + pgphase_trace::genotype_text(c));
        }
    }

    // An MSA-genotyped het whose reads split too unevenly to be a site (kEmMinMinorAlleleFraction) is still a
    // het: the reads call a long homopolymer one short as often as not (chr20:2,060,967, T22 against T24, the MSA's
    // clusters of 49 and 28 reads, the reads' calls 48 REF to 13 ALT). Written unphased, as an unphasable one above,
    // when the reads the MSA genotyped it from are at least kEmMinMsaReadsFraction of the chunk's median site
    // depth; fewer is an MSA of a handful of phased reads (chr20:15,363,779, 8 of 79), mostly a false het.
    void write_skewed_msa_hets() {
        std::vector<int> depths;
        for (size_t k = 0; k < n_; ++k) {
            if (sites_[k].locus >= 0) continue;
            int n = 0;
            for (const auto& [a, m] : allele_counts_[sites_[k].ci]) n += m;
            depths.push_back(n);
        }
        if (depths.empty()) return;
        std::nth_element(depths.begin(), depths.begin() + depths.size() / 2, depths.end());
        const double min_reads = kEmMinMsaReadsFraction * depths[depths.size() / 2];
        for (size_t ci = 0; ci < chunk_.candidates.size(); ++ci) {
            CandidateVariant& c = chunk_.candidates[ci];
            if (site_of_[ci] >= 0 || !c.msa_genotyped || c.island_only || c.low_mapq_only || c.unresolved_noisy ||
                (c.lcd_var_i_to_cate & kCandGermlineVarCate) == 0)
                continue;
            const std::vector<int>& covs = c.counts.alle_covs;
            const int msa_reads = covs.empty() ? c.counts.ref_cov + c.counts.alt_cov
                                               : std::accumulate(covs.begin(), covs.end(), 0);
            if (msa_reads < min_reads) continue;
            std::vector<std::pair<int, int>> by_count;
            for (const auto& [a, n] : allele_counts_[ci]) by_count.emplace_back(n, a);
            std::sort(by_count.rbegin(), by_count.rend());
            const bool two_alts = c.msa_insertion_alts.size() == 2;  // its alleles the MSA's, as in select_sites
            const int x = two_alts ? 1 : by_count.size() >= 2 ? by_count[0].second : 0;
            const int y = two_alts ? 2 : by_count.size() >= 2 ? by_count[1].second : 1;
            set_site_phase(c, kUnsetCandidatePhaseSet, x, y);
            c.vcf_unphased = true;
            pgphase_trace::candidate("em", "site", c, "skewed_msa_het_unphased", pgphase_trace::genotype_text(c));
        }
    }

    // A last-resort label: applied after stitching and the rescue, fill-only
    // (apply_island_labels), so it steers no window, stitch or rescue.
    void defer_label(size_t ri, int hap, hts_pos_t ps) {
        if (chunk_.island_haps.size() < chunk_.reads.size()) {
            chunk_.island_haps.resize(chunk_.reads.size(), 0);
            chunk_.island_phase_sets.resize(chunk_.reads.size(), 0);
        }
        chunk_.island_haps[ri] = hap;
        chunk_.island_phase_sets[ri] = ps;
    }

    // Reads: label in the block where their own evidence is strongest. A read
    // below the label posterior takes, in turn, a deferred label (its calls in
    // a block all agree, kEmDeferredLabelPosterior) or one from its unverified
    // calls (kEmUnverifiedLabelPosterior). Labels from one-edit window calls
    // came after these until fill labels were written under their block's PS:
    // scored so they were +6 correct / +16 discordant, and were removed.
    void label_reads() {
        for (size_t ri = 0; ri < robs_.size() && ri < chunk_.haps.size(); ++ri) {
            chunk_.haps[ri] = 0;
            chunk_.phase_sets[ri] = kUnphasedReadPhaseSet;
            std::map<hts_pos_t, std::array<int, 2>> directions;  // calls favouring hap 1 / hap 2
            std::map<hts_pos_t, std::array<double, 2>> per_block;
            for (const EmObs& o : robs_[ri]) {
                if (err_[static_cast<size_t>(o.s)] >= kEmUnreliableError) continue;
                // Unverified calls (join_excluded) are phased for the VCF but do
                // not label reads: on chr12, chr18 and chr20 leaving them out took
                // 90, 125 and 12 reads off the discordant count and added correct
                // ones on chr12 and chr18.
                if (join_excluded_[static_cast<size_t>(o.s)]) continue;
                const auto ll = label_log_likelihoods(o);
                auto& l = per_block[block_[static_cast<size_t>(o.s)]];
                l[0] += ll[0];
                l[1] += ll[1];
                if (ll[0] != ll[1]) ++directions[block_[static_cast<size_t>(o.s)]][ll[0] > ll[1] ? 0 : 1];
            }
            double best = 0.0;
            for (const auto& [ps, l] : per_block) {
                const double margin = l[0] - l[1];
                const double p1 = hap1_posterior(margin);
                if (std::max(p1, 1.0 - p1) < kEmReadLabelPosterior || std::abs(margin) <= best) continue;
                best = std::abs(margin);
                chunk_.haps[ri] = p1 >= 0.5 ? 1 : 2;
                chunk_.phase_sets[ri] = ps;
            }
            if (chunk_.haps[ri] != 0 || loci_ == nullptr) continue;
            // Deferred: all of the read's calls in a block point one way.
            double weak_best = 0.0;
            int weak_hap = 0;
            hts_pos_t weak_ps = 0;
            for (const auto& [ps, l] : per_block) {
                const double margin = l[0] - l[1];
                const double p1 = hap1_posterior(margin);
                const auto& dir = directions[ps];
                if (dir[0] != 0 && dir[1] != 0) continue;
                if (std::max(p1, 1.0 - p1) < kEmDeferredLabelPosterior || std::abs(margin) <= weak_best) continue;
                weak_best = std::abs(margin);
                weak_hap = p1 >= 0.5 ? 1 : 2;
                weak_ps = ps;
            }
            // Failing that, its unverified calls.
            if (weak_hap == 0) {
                std::map<hts_pos_t, std::array<double, 2>> unverified;
                for (const EmObs& o : robs_[ri]) {
                    const size_t k = static_cast<size_t>(o.s);
                    if (!join_excluded_[k] || err_[k] >= kEmUnverifiedReportError) continue;
                    const auto ll = label_log_likelihoods(o);
                    auto& l = unverified[block_[k]];
                    l[0] += ll[0];
                    l[1] += ll[1];
                }
                for (const auto& [ps, l] : unverified) {
                    const double margin = l[0] - l[1];
                    const double p1 = hap1_posterior(margin);
                    if (std::max(p1, 1.0 - p1) < kEmUnverifiedLabelPosterior || std::abs(margin) <= weak_best) continue;
                    weak_best = std::abs(margin);
                    weak_hap = p1 >= 0.5 ? 1 : 2;
                    weak_ps = ps;
                }
            }
            if (weak_hap != 0) defer_label(ri, weak_hap, weak_ps);
        }
    }

    // Island blocks. A read the solve leaves unlabelled is phased locally by
    // its own call at one site that cannot vote: a block of that site alone,
    // never joined to anything, so a wrong island costs its own reads and
    // never flips a neighbour. Candidates: sites the EM rejects (learned error
    // kEmUnreliableError or more), repeat het indels it never sees
    // (REP_HET_INDEL) and MSA-verified hets injection kept out (island_only),
    // the last two judged by their raw disagreement with the labelled reads,
    // and island-only length windows (repeat loci genotyped from read lengths;
    // ungated). Each read takes its lowest-error candidate. Recorded here,
    // applied after the rescue to reads still unphased (apply_island_labels).
    // Measured on chr12/18/20 against c8391f1: +2,349 correct, +457
    // discordant, -2,806 unphased reads, one arm switch; islands at truth hets
    // run 18-20 : 1, at repeat-length noise 1.5 : 1.
    void island_blocks() {
        struct IslandCand { size_t ci; int x; int y; hts_pos_t pos; int k; double err; int phase; char kind; };
        std::vector<IslandCand> cands;
        std::vector<int> cand_of_site(n_, -1);
        for (size_t k = 0; k < n_; ++k) {
            if (sites_[k].locus >= 0 || join_excluded_[k] || err_[k] < kEmUnreliableError) continue;
            cand_of_site[k] = static_cast<int>(cands.size());
            cands.push_back(IslandCand{sites_[k].ci, sites_[k].x, sites_[k].y, sites_[k].pos, static_cast<int>(k),
                                       err_[k], phase_[k], 'E'});
        }
        std::vector<int> cand_of_ci(chunk_.candidates.size(), -1);
        for (size_t ci = 0; ci < chunk_.candidates.size(); ++ci) {
            const CandidateVariant& c = chunk_.candidates[ci];
            if ((c.counts.category != VariantCategory::RepeatHetIndel && !c.island_only) || site_of_[ci] >= 0) continue;
            std::vector<std::pair<int, int>> by_count;
            for (const auto& [al, cnt] : allele_counts_[ci]) by_count.emplace_back(cnt, al);
            std::sort(by_count.rbegin(), by_count.rend());
            if (by_count.size() < 2) continue;
            const int nx = by_count[0].first, ny = by_count[1].first;
            if (ny < kEmMinMinorAlleleFraction * (nx + ny)) continue;
            cand_of_ci[ci] = static_cast<int>(cands.size());
            cands.push_back(IslandCand{ci, by_count[0].second, by_count[1].second, c.key.sort_pos(), -1, 1.0, 0,
                                       c.island_only ? 'M' : 'R'});
        }
        // Island-only length windows (kind W; build_length_window_sites): sides 0
        // and 1 are the shorter and longer length.
        std::vector<int> cand_of_locus(n_locus_, -1);
        for (size_t li = 0; li < n_locus_; ++li) {
            if (!(*loci_)[li].island_only) continue;
            cand_of_locus[li] = static_cast<int>(cands.size());
            cands.push_back(IslandCand{SIZE_MAX, 0, 1, (*loci_)[li].pos, -1, 1.0, 0, 'W'});
        }
        // A site outside the solve: its raw disagreement with the labelled
        // reads, block by block; too few labelled reads to judge counts as the
        // unreliable cut.
        {
            std::vector<std::map<hts_pos_t, std::array<int, 2>>> agree(cands.size());
            for (size_t ri = 0; ri < chunk_.read_var_profile.size() && ri < chunk_.haps.size(); ++ri) {
                const int hap = chunk_.haps[ri];
                if (hap != 1 && hap != 2) continue;
                const ReadVariantProfile& prof = chunk_.read_var_profile[ri];
                if (prof.start_var_idx < 0 || chunk_.reads[ri].is_skipped) continue;
                for (size_t j = 0; j < prof.alleles.size(); ++j) {
                    const int ic = cand_of_ci[static_cast<size_t>(prof.start_var_idx) + j];
                    if (ic < 0) continue;
                    const int al = prof.alleles[j];
                    const IslandCand& cand = cands[static_cast<size_t>(ic)];
                    if (al != cand.x && al != cand.y) continue;
                    ++agree[static_cast<size_t>(ic)][chunk_.phase_sets[ri]][((al == cand.x) == (hap == 1)) ? 0 : 1];
                }
            }
            for (size_t li = 0; li < n_locus_; ++li) {
                const int ic = cand_of_locus[li];
                if (ic < 0) continue;
                for (const auto& [ri, side, e] : (*loci_)[li].observations) {
                    if (ri >= chunk_.haps.size()) continue;
                    const int hap = chunk_.haps[ri];
                    if (hap != 1 && hap != 2) continue;
                    ++agree[static_cast<size_t>(ic)][chunk_.phase_sets[ri]][((side == 0) == (hap == 1)) ? 0 : 1];
                }
            }
            for (size_t ic = 0; ic < cands.size(); ++ic) {
                if (cands[ic].kind == 'E') continue;
                int best = 0, total = 0;
                for (const auto& [ps, v] : agree[ic]) { best += std::max(v[0], v[1]); total += v[0] + v[1]; }
                cands[ic].err = total >= 3 ? 1.0 - static_cast<double>(best) / total : kEmUnreliableError;
            }
        }
        // A length window is not gated by agreement with the labelled reads: its
        // two lengths are measured from the reads themselves. Gated, 5 of ~270
        // unphased reads it calls were labelled; ungated +232 correct / +36
        // discordant over chr12/18/20.
        const auto usable = [&](size_t ic) {
            return cands[ic].kind == 'E' || cands[ic].kind == 'W' || cands[ic].err < kIslandOutsideMaxError;
        };
        // Each unlabelled read's calls at the candidates.
        std::vector<std::vector<std::pair<int, int>>> iobs(robs_.size());
        for (size_t ri = 0; ri < robs_.size() && ri < chunk_.haps.size(); ++ri) {
            if (chunk_.haps[ri] != 0) continue;
            if (ri < chunk_.island_haps.size() && chunk_.island_haps[ri] != 0) continue;  // deferred block label
            for (const EmObs& o : robs_[ri])
                if (cand_of_site[static_cast<size_t>(o.s)] >= 0)
                    iobs[ri].emplace_back(cand_of_site[static_cast<size_t>(o.s)], o.b);
            if (ri >= chunk_.read_var_profile.size()) continue;
            const ReadVariantProfile& prof = chunk_.read_var_profile[ri];
            if (prof.start_var_idx < 0 || chunk_.reads[ri].is_skipped) continue;
            for (size_t j = 0; j < prof.alleles.size(); ++j) {
                const int ic = cand_of_ci[static_cast<size_t>(prof.start_var_idx) + j];
                if (ic < 0 || !usable(static_cast<size_t>(ic))) continue;
                const int al = prof.alleles[j];
                if (al == cands[static_cast<size_t>(ic)].x) iobs[ri].emplace_back(ic, 0);
                else if (al == cands[static_cast<size_t>(ic)].y) iobs[ri].emplace_back(ic, 1);
            }
        }
        for (size_t li = 0; li < n_locus_; ++li) {
            const int ic = cand_of_locus[li];
            if (ic < 0 || !usable(static_cast<size_t>(ic))) continue;
            for (const auto& [ri, side, e] : (*loci_)[li].observations) {
                if (ri >= iobs.size() || ri >= chunk_.haps.size() || chunk_.haps[ri] != 0) continue;
                if (ri < chunk_.island_haps.size() && chunk_.island_haps[ri] != 0) continue;
                iobs[ri].emplace_back(ic, side);
            }
        }
        std::vector<int> unlabelled_callers(cands.size(), 0);
        for (const auto& v : iobs)
            for (const auto& [ic, side] : v) ++unlabelled_callers[static_cast<size_t>(ic)];
        std::vector<int> pick_cand(robs_.size(), -1), pick_side(robs_.size(), 0);
        for (size_t ri = 0; ri < iobs.size(); ++ri) {
            int best = -1, best_side = 0;
            for (const auto& [ic, side] : iobs[ri]) {
                const size_t u = static_cast<size_t>(ic);
                if (unlabelled_callers[u] < kIslandMinReads) continue;
                if (best >= 0) {
                    const size_t b = static_cast<size_t>(best);
                    if (cands[u].err > cands[b].err) continue;
                    if (cands[u].err == cands[b].err && unlabelled_callers[u] <= unlabelled_callers[b]) continue;
                }
                best = ic;
                best_side = side;
            }
            pick_cand[ri] = best;
            pick_side[ri] = best_side;
        }
        std::set<hts_pos_t> used_ps(block_.begin(), block_.end());
        std::map<size_t, hts_pos_t> island_ps;  // candidate -> its own phase set
        const auto ps_of = [&](size_t ic) {
            auto it = island_ps.find(ic);
            if (it != island_ps.end()) return it->second;
            hts_pos_t ps = cands[ic].pos;
            while (used_ps.count(ps)) ++ps;
            used_ps.insert(ps);
            island_ps.emplace(ic, ps);
            return ps;
        };
        for (size_t ri = 0; ri < robs_.size(); ++ri) {
            if (pick_cand[ri] < 0) continue;
            const size_t ic = static_cast<size_t>(pick_cand[ri]);
            defer_label(ri, pick_side[ri] == cands[ic].phase ? 1 : 2, ps_of(ic));
        }
        chunk_.island_sites.clear();
        for (const auto& [ic, ps] : island_ps) {
            const IslandCand& c = cands[ic];
            if (c.kind == 'W') continue;  // a window writes no VCF record
            const int hap1 = c.phase == 0 ? c.x : c.y;
            const int hap2 = c.phase == 0 ? c.y : c.x;
            chunk_.island_sites.push_back(PhasingChunk::IslandSite{c.ci, ps, hap1, hap2});
        }
        // Diagnostic, off by default: PGPHASE_DIAG=DIR writes every read an
        // island claims (islands.TID.BEG.tsv): read, site position, island PS,
        // the read's side, its island haplotype, its read error at the site,
        // the site's error (learned, or raw outside the solve), unlabelled
        // callers, the most member reads sharing one other verified site, 1,
        // kind (E rejected EM site, R repeat het indel, M island-only MSA het,
        // W length window).
        static const char* const island_diag = std::getenv("PGPHASE_DIAG");
        if (island_diag == nullptr) return;
        std::vector<int> linked_reads(cands.size(), 0);
        std::map<std::pair<int, int>, int> shared;  // (island, other EM site) -> member reads calling it
        for (size_t ri = 0; ri < robs_.size(); ++ri) {
            if (pick_cand[ri] < 0) continue;
            const int own = cands[static_cast<size_t>(pick_cand[ri])].k;
            for (const EmObs& o : robs_[ri]) {
                const size_t k = static_cast<size_t>(o.s);
                if (o.s == own || sites_[k].locus >= 0 || chunk_.candidates[sites_[k].ci].unresolved_noisy) continue;
                int& c = ++shared[{pick_cand[ri], o.s}];
                int& m = linked_reads[static_cast<size_t>(pick_cand[ri])];
                m = std::max(m, c);
            }
        }
        const std::string path = std::string(island_diag) + "/islands." + std::to_string(chunk_.region.tid) + "." +
                                 std::to_string(chunk_.region.beg) + ".tsv";
        if (FILE* fh = std::fopen(path.c_str(), "w")) {
            for (size_t ri = 0; ri < robs_.size(); ++ri) {
                if (pick_cand[ri] < 0) continue;
                const size_t ic = static_cast<size_t>(pick_cand[ri]);
                double read_err = 0.0;
                for (const EmObs& o : robs_[ri])
                    if (o.s == cands[ic].k) read_err = o.e;
                std::fprintf(fh, "%s\t%lld\t%lld\t%d\t%d\t%.4f\t%.3f\t%d\t%d\t1\t%c\n", chunk_.reads[ri].qname.c_str(),
                             static_cast<long long>(cands[ic].pos), static_cast<long long>(island_ps.at(ic)),
                             pick_side[ri], pick_side[ri] == cands[ic].phase ? 1 : 2, read_err, cands[ic].err,
                             unlabelled_callers[ic], linked_reads[ic], cands[ic].kind);
            }
            std::fclose(fh);
        }
    }

    // Diagnostics, off by default (only for the union solve, loci given).
    void write_diagnostics() const {
        if (loci_ == nullptr) return;
        // PGPHASE_DIAG=DIR writes the final solve's voting candidates
        // (emsites.TID.BEG.tsv): VCF-convention position, type, REF length, ALT,
        // origin (C catalog clean, O catalog other, I injected, U unverified),
        // allele read counts (allele:reads;...), reads in range without an
        // allele, EM site or "-" (left out: skewed alleles), the two alleles
        // contrasted, learned error, phase (which allele haplotype 1 carries),
        // block.
        static const char* const diag_dir = std::getenv("PGPHASE_DIAG");
        if (diag_dir != nullptr) {
            const std::string path = std::string(diag_dir) + "/emsites." + std::to_string(chunk_.region.tid) + "." +
                                     std::to_string(chunk_.region.beg) + ".tsv";
            if (FILE* fh = std::fopen(path.c_str(), "w")) {
                for (size_t ci = 0; ci < chunk_.candidates.size(); ++ci) {
                    const CandidateVariant& c = chunk_.candidates[ci];
                    if ((c.lcd_var_i_to_cate & kCandGermlineVarCate) == 0) continue;
                    const char type = c.key.type == VariantType::Snp ? 'S' : c.key.type == VariantType::Insertion ? 'I' : 'D';
                    const char origin = c.unresolved_noisy ? 'U' : c.bam_injected ? 'I'
                                      : (c.lcd_var_i_to_cate & (kCandCleanHetSnp | kCandCleanHetIndel)) != 0 ? 'C' : 'O';
                    std::string counts;
                    for (const auto& [a, n] : allele_counts_[ci])
                        counts += (counts.empty() ? "" : ";") + std::to_string(a) + ":" + std::to_string(n);
                    std::fprintf(fh, "%lld\t%c\t%d\t%s\t%c\t%s\t%d", static_cast<long long>(c.key.sort_pos()), type,
                                 c.key.ref_len, c.key.alt.substr(0, 40).c_str(), origin,
                                 counts.empty() ? "." : counts.c_str(), undecided_[ci]);
                    const int k = site_of_[ci];
                    if (k < 0) {
                        std::fprintf(fh, "\t-\t.\t.\t.\t.\t.\n");
                    } else {
                        const EmSite& site = sites_[static_cast<size_t>(k)];
                        std::fprintf(fh, "\tem\t%d\t%d\t%.3f\t%d\t%lld\n", site.x, site.y, err_[static_cast<size_t>(k)],
                                     phase_[static_cast<size_t>(k)] == 0 ? site.x : site.y,
                                     static_cast<long long>(block_[static_cast<size_t>(k)]));
                    }
                }
                std::fclose(fh);
            }
            // Every read's profile as the final solve saw it (emreads.TID.BEG.tsv):
            // name, MAPQ, skipped flag, then each call as VCF-convention position
            // of the candidate : allele, including candidates the EM left out.
            const std::string rpath = std::string(diag_dir) + "/emreads." + std::to_string(chunk_.region.tid) + "." +
                                      std::to_string(chunk_.region.beg) + ".tsv";
            if (FILE* fh = std::fopen(rpath.c_str(), "w")) {
                for (size_t ri = 0; ri < chunk_.reads.size(); ++ri) {
                    const ReadRecord& r = chunk_.reads[ri];
                    std::fprintf(fh, "%s\t%d\t%d", r.qname.c_str(), r.mapq, r.is_skipped ? 1 : 0);
                    const ReadVariantProfile* prof = nullptr;
                    for (const ReadVariantProfile& p2 : chunk_.read_var_profile)
                        if (p2.read_id == static_cast<int>(ri)) { prof = &p2; break; }
                    if (prof == nullptr) {
                        std::fputs("\tno_profile\n", fh);
                        continue;
                    }
                    if (prof->start_var_idx < 0) {
                        std::fputs("\tempty\n", fh);
                        continue;
                    }
                    for (size_t k = 0; k < prof->alleles.size(); ++k) {
                        const size_t ci = static_cast<size_t>(prof->start_var_idx) + k;
                        if (ci >= chunk_.candidates.size()) break;
                        std::fprintf(fh, "\t%lld:%d", static_cast<long long>(chunk_.candidates[ci].key.sort_pos()),
                                     prof->alleles[k]);
                    }
                    std::fputc('\n', fh);
                }
                std::fclose(fh);
            }
            // Every block-cut decision (joins.TID.BEG.tsv): between each EM site
            // and the next, their VCF-convention positions, the join margin (log L
            // of the solve minus log L with every site to the right flipped, each
            // locus counted once), the margin without weak sites, reads with
            // observations on both sides, the decision -- J joined, C cut (margin
            // below the cut), S cut (too few reads cross), W cut (a weak bridge
            // below its larger cut) -- and reads with strong evidence on both
            // sides.
            const std::string jpath = std::string(diag_dir) + "/joins." + std::to_string(chunk_.region.tid) + "." +
                                      std::to_string(chunk_.region.beg) + ".tsv";
            if (FILE* fh = std::fopen(jpath.c_str(), "w")) {
                for (size_t k = 0; k + 1 < n_; ++k) {
                    const char decision = delta_[k] < kEmBlockCutLogLikelihood ? 'C'
                                        : spanning_[k] < kEmMinJoinReads ? 'S'
                                        : weak_bridge_[k] && delta_[k] < kEmWeakBridgeCutLogLikelihood ? 'W' : 'J';
                    std::fprintf(fh, "%lld\t%lld\t%.2f\t%.2f\t%d\t%c\t%d\n", static_cast<long long>(sites_[k].pos),
                                 static_cast<long long>(sites_[k + 1].pos), delta_[k], delta_strong_[k], spanning_[k],
                                 decision, spanning_strong_[k]);
                }
                std::fclose(fh);
            }
        }
        // PGPHASE_EM_DUMP=DIR writes the final solve's evidence per read, one
        // file per chunk. A line is the read, MAPQ, label (HP, PS), per-block
        // log-likelihood margins (PS=margin), the label the read started the
        // solve with, then each observation as
        // POS:KIND:SIDE:SITE_ERROR:BLOCK:WEAK:READ_ERROR, with KIND C clean het,
        // O other, I injected from the alignment, U injected from a noisy
        // region too long for MSA, L window.
        static const char* const dump_dir = std::getenv("PGPHASE_EM_DUMP");
        if (dump_dir == nullptr) return;
        const std::string path = std::string(dump_dir) + "/em." + std::to_string(chunk_.region.tid) + "." +
                                 std::to_string(chunk_.region.beg) + ".tsv";
        FILE* fh = std::fopen(path.c_str(), "w");
        if (fh == nullptr) return;
        for (size_t ri = 0; ri < chunk_.reads.size() && ri < robs_.size(); ++ri) {
            const ReadRecord& r = chunk_.reads[ri];
            std::fprintf(fh, "%s\t%d\t%d\t%lld\t", r.qname.c_str(), r.mapq, chunk_.haps[ri],
                         static_cast<long long>(chunk_.phase_sets[ri]));
            std::map<hts_pos_t, double> margin;
            for (const EmObs& o : robs_[ri]) {
                if (err_[static_cast<size_t>(o.s)] >= kEmUnreliableError) continue;
                if (join_excluded_[static_cast<size_t>(o.s)]) continue;
                const auto ll = label_log_likelihoods(o);
                margin[block_[static_cast<size_t>(o.s)]] += ll[0] - ll[1];
            }
            const char* sep = "";
            for (const auto& [ps, m] : margin) {
                std::fprintf(fh, "%s%lld=%.2f", sep, static_cast<long long>(ps), m);
                sep = ",";
            }
            std::fprintf(fh, "\t%d", start_label_[ri]);
            for (const EmObs& o : robs_[ri]) {
                const EmSite& site = sites_[static_cast<size_t>(o.s)];
                char kind = 'L';
                if (site.locus < 0) {
                    const CandidateVariant& c = chunk_.candidates[site.ci];
                    kind = c.unresolved_noisy ? 'U'
                         : c.bam_injected ? 'I'
                         : (c.lcd_var_i_to_cate & (kCandCleanHetSnp | kCandCleanHetIndel)) != 0 ? 'C' : 'O';
                }
                std::fprintf(fh, "\t%lld:%c:%d:%.3f:%lld:%d:%.4f", static_cast<long long>(site.pos), kind, o.b,
                             err_[static_cast<size_t>(o.s)], static_cast<long long>(block_[static_cast<size_t>(o.s)]),
                             static_cast<int>(site_weak_[static_cast<size_t>(o.s)]), o.e);
            }
            std::fputc('\n', fh);
        }
        std::fclose(fh);
    }
};

}  // namespace

size_t phase_chunk_by_global_em(PhasingChunk& chunk, const std::vector<LocusWindowSite>* loci,
                                const std::vector<char>* bridge_weak,
                                const ObservationErrors* errors) {
    return GlobalEm(chunk, loci, bridge_weak, errors).run();
}

// ── Island labels ───────────────────────────────────────────────────────────

size_t apply_island_labels(PhasingChunk& chunk) {
    if (chunk.island_haps.empty()) return 0;
    if (chunk.gap_haps.size() < chunk.reads.size()) chunk.gap_haps.resize(chunk.reads.size(), 0);
    if (chunk.gap_phase_sets.size() < chunk.reads.size()) chunk.gap_phase_sets.resize(chunk.reads.size(), 0);
    if (chunk.gap_from_bam_observation.size() < chunk.reads.size())
        chunk.gap_from_bam_observation.resize(chunk.reads.size(), false);
    size_t labelled = 0;
    for (size_t ri = 0; ri < chunk.island_haps.size() && ri < chunk.reads.size(); ++ri) {
        if (chunk.island_haps[ri] == 0) continue;
        if (ri < chunk.haps.size() && chunk.haps[ri] != 0) continue;
        if (chunk.gap_haps[ri] != 0) continue;
        chunk.gap_haps[ri] = chunk.island_haps[ri];
        chunk.gap_phase_sets[ri] = chunk.island_phase_sets[ri] + kGapFillPsOffset;
        chunk.gap_from_bam_observation[ri] = true;  // fill only: never replaces another chunk's label
        ++labelled;
    }
    for (const PhasingChunk::IslandSite& s : chunk.island_sites) {
        if (s.candidate >= chunk.candidates.size()) continue;
        CandidateVariant& c = chunk.candidates[s.candidate];
        if (c.hap_to_cons_alle[1] >= 0 && c.hap_to_cons_alle[2] >= 0) continue;  // phased since
        set_site_phase(c, s.phase_set, s.hap1, s.hap2);
    }
    return labelled;
}

} // namespace pgphase_collect
