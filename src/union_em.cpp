// Union gap phasing: the global EM over the read x site matrix, with block
// cuts and the weak-bridge guard.

#include "union_internal.hpp"

#include <array>
#include <cstdint>
#include <cstdio>
#include <cstdlib>

namespace pgphase_collect {

// Error-rate bounds; starting error for clean and other sites; the learned
// error at which a site is neither reported phased nor used to label reads;
// iteration limits; the read posterior needed for a label; the boundary
// log-likelihood below which blocks are cut; and the MAPQ from which a read
// teaches the model (lower reads are labelled from it but do not shape it).
constexpr double kEmMinError = 0.01;
constexpr double kEmMaxError = 0.45;
constexpr double kEmCleanStartError = 0.05;
constexpr double kEmOtherStartError = 0.2;
constexpr double kEmUnreliableError = 0.3;
// The learned error up to which an unverified call (unresolved_noisy) is
// written phased. Such calls label no reads, so this decides only the VCF:
// 0.1 against 0.3 cut arm switches 12 / 24 / 10 -> 10 / 20 / 6 on chr12,
// chr18 and chr20 and kept 60-75% of the truth heterozygotes they add.
constexpr double kEmUnverifiedReportError = 0.1;
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
// through a homopolymer call of learned error 0.27.
constexpr int kEmMinStrongJoinReads = 3;
// A site whose less observed allele carries under this fraction of its reads
// does not take part: a skewed split is the signature of a paralog or error
// call, and in a homozygous stretch such a site labels reads at random.
constexpr double kEmMinMinorAlleleFraction = 0.25;
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
// Island blocks (see the end of phase_chunk_by_global_em). A rejected EM site
// offers islands below this learned error: 0.35 gave +975/+191 reads over
// chr12/18/20, 0.40 +1,458/+415, before the other candidates and
// lowest-error choice.
constexpr double kIslandMaxError = 0.40;
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

size_t phase_chunk_by_global_em(PhasingChunk& chunk, const std::vector<LocusWindowSite>* loci,
                                const std::vector<char>* bridge_weak,
                                const ObservationErrors* errors) {
    struct EmSite { long locus; size_t ci; int x; int y; hts_pos_t pos; };  // locus < 0: candidate
    struct EmObs { int s; int b; double e; };
    // Sites in position order with their two-allele contrast; a locus window
    // contrasts its two local haplotype sequences (0 and 1).
    std::vector<std::map<int, int>> allele_counts(chunk.candidates.size());
    std::vector<int> undecided(chunk.candidates.size(), 0);  // reads in range without an allele (diagnostic)
    for (size_t ri = 0; ri < chunk.read_var_profile.size(); ++ri) {
        const ReadVariantProfile& prof = chunk.read_var_profile[ri];
        if (prof.start_var_idx < 0 || chunk.reads[ri].is_skipped) continue;
        for (size_t k = 0; k < prof.alleles.size(); ++k) {
            if (prof.alleles[k] >= 0)
                ++allele_counts[static_cast<size_t>(prof.start_var_idx) + k][prof.alleles[k]];
            else
                ++undecided[static_cast<size_t>(prof.start_var_idx) + k];
        }
    }
    std::vector<EmSite> sites;
    for (size_t ci = 0; ci < chunk.candidates.size(); ++ci) {
        const CandidateVariant& c = chunk.candidates[ci];
        if ((c.lcd_var_i_to_cate & kCandGermlineVarCate) == 0) continue;
        int x = c.hap_to_cons_alle[1], y = c.hap_to_cons_alle[2];
        if (x < 0 || y < 0 || x == y) {
            std::vector<std::pair<int, int>> by_count;
            for (const auto& [a, n] : allele_counts[ci]) by_count.emplace_back(n, a);
            std::sort(by_count.rbegin(), by_count.rend());
            if (by_count.size() < 2) continue;
            x = by_count[0].second;
            y = by_count[1].second;
        }
        const auto count_of = [&](int a) {
            const auto it = allele_counts[ci].find(a);
            return it == allele_counts[ci].end() ? 0 : it->second;
        };
        const int nx = count_of(x), ny = count_of(y);
        if (nx + ny == 0 || std::min(nx, ny) < kEmMinMinorAlleleFraction * (nx + ny)) continue;
        sites.push_back(EmSite{-1, ci, x, y, c.key.sort_pos()});
    }
    const size_t n_locus = loci != nullptr ? loci->size() : 0;
    for (size_t li = 0; li < n_locus; ++li)
        if (!(*loci)[li].island_only) sites.push_back(EmSite{static_cast<long>(li), 0, 0, 1, (*loci)[li].pos});
    std::stable_sort(sites.begin(), sites.end(),
                     [](const EmSite& a, const EmSite& b) { return a.pos < b.pos; });
    const size_t n = sites.size();
    if (n < 2) return 0;
    std::vector<int> site_of(chunk.candidates.size(), -1);
    std::vector<int> locus_site(n_locus, -1);
    for (size_t k = 0; k < n; ++k) {
        if (sites[k].locus >= 0) locus_site[static_cast<size_t>(sites[k].locus)] = static_cast<int>(k);
        else site_of[sites[k].ci] = static_cast<int>(k);
    }
    std::vector<int> phase(n, 0);  // 0: hap1 carries x
    std::vector<double> err(n);
    for (size_t k = 0; k < n; ++k) {
        const bool clean = sites[k].locus < 0 &&
            (chunk.candidates[sites[k].ci].lcd_var_i_to_cate & (kCandCleanHetSnp | kCandCleanHetIndel)) != 0;
        err[k] = clean ? kEmCleanStartError : kEmOtherStartError;
    }
    // Read observations as (site, side, own error) in site order.
    std::vector<std::vector<EmObs>> robs(chunk.reads.size());
    for (size_t ri = 0; ri < chunk.read_var_profile.size(); ++ri) {
        const ReadVariantProfile& prof = chunk.read_var_profile[ri];
        if (prof.start_var_idx < 0 || chunk.reads[ri].is_skipped) continue;
        for (size_t k = 0; k < prof.alleles.size(); ++k) {
            const int s = site_of[static_cast<size_t>(prof.start_var_idx) + k];
            if (s < 0 || prof.alleles[k] < 0) continue;
            const EmSite& site = sites[static_cast<size_t>(s)];
            double e = 0.0;
            if (errors != nullptr) {
                const auto it = errors->find(observation_key(ri, static_cast<size_t>(prof.start_var_idx) + k));
                if (it != errors->end()) e = it->second;
            }
            if (prof.alleles[k] == site.x) robs[ri].push_back(EmObs{s, 0, e});
            else if (prof.alleles[k] == site.y) robs[ri].push_back(EmObs{s, 1, e});
        }
    }
    for (size_t li = 0; li < n_locus; ++li) {
        if (locus_site[li] < 0) continue;  // island-only window
        for (const auto& [ri, side, e] : (*loci)[li].observations)
            if (ri < robs.size()) robs[ri].push_back(EmObs{locus_site[li], side, static_cast<double>(e)});
    }
    std::vector<double> post(chunk.reads.size(), 0.5);
    std::vector<int> start_label(chunk.reads.size(), 0);  // for the diagnostic dump
    std::vector<char> learns(chunk.reads.size(), 0);
    for (size_t ri = 0; ri < robs.size(); ++ri) {
        std::sort(robs[ri].begin(), robs[ri].end(),
                  [](const EmObs& a, const EmObs& b) { return a.s < b.s; });
        const int hap = ri < chunk.haps.size() ? chunk.haps[ri] : 0;
        post[ri] = hap == 1 ? 0.95 : hap == 2 ? 0.05 : 0.5;
        if (ri < start_label.size()) start_label[ri] = hap;
        learns[ri] = chunk.reads[ri].mapq >= kEmLearnMinMapq;
    }
    // Echo sites (see kEmEchoMaxGap): each candidate pair within the gap,
    // counted over the reads that call both.
    std::vector<char> echo(n, 0);
    {
        std::vector<int> callers(n, 0);
        std::map<std::pair<int, int>, std::array<int, 2>> pair_calls;  // (k, j): same / opposite sides
        for (const std::vector<EmObs>& obs : robs) {
            for (size_t i = 0; i < obs.size(); ++i) {
                ++callers[static_cast<size_t>(obs[i].s)];
                const EmSite& si = sites[static_cast<size_t>(obs[i].s)];
                if (si.locus >= 0) continue;
                for (size_t j = i + 1; j < obs.size(); ++j) {
                    const EmSite& sj = sites[static_cast<size_t>(obs[j].s)];
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
                echo[static_cast<size_t>(kj.second)] = 1;
        }
    }
    const auto obs_error = [&](const EmObs& o) {
        return std::min(kEmMaxError, 1.0 - (1.0 - err[static_cast<size_t>(o.s)]) * (1.0 - o.e));
    };
    // Weak sites (see bridge_weak_sites) are phased but cannot hold blocks
    // together: the final cuts are judged on the other sites' evidence.
    std::vector<char> site_weak(n, 0);
    if (bridge_weak != nullptr)
        for (size_t k = 0; k < n; ++k)
            if (sites[k].locus < 0 && sites[k].ci < bridge_weak->size()) site_weak[k] = (*bridge_weak)[sites[k].ci];
    // The join test counts a locus once. A locus window and a candidate
    // inside it are read from the same bases, so a read's two calls there
    // agree by construction; counted twice they outvote a real heterozygote
    // read once. And a window over a homopolymer length call is that call
    // re-read, so it is weak too. At chr12:122,284,562 a window and the +1
    // indel inside it (one locus, purity 0.53 against the parents) carried a
    // 30 kb join past a true het (122,290,121, purity 0.95) with a margin of
    // 5.65 from "strong" sites, and 165 reads went discordant. Labels and
    // learning keep every observation: only the decision to join uses this.
    std::vector<char> join_weak = site_weak;
    std::vector<std::vector<char>> join_skip(robs.size());
    // Calls kept from noisy regions too long for MSA (unresolved_noisy) are
    // left out of the join test altogether: a divergent copy's mismatches
    // (695 such calls in 24 kb at chr20:49.65 Mb) clear even the weak-bridge
    // margin by sheer number, and flipped a block of 110 reads. Clean hets no
    // confidently mapped read carries (low_mapq_only) are treated the same:
    // as ordinary sites they added two arm switches over chr12/18/20 (chr12
    // 34.64 Mb, chr20 32.25 Mb) for three correct reads.
    std::vector<char> join_excluded(n, 0);
    for (size_t k = 0; k < n; ++k)
        if (sites[k].locus < 0 && (chunk.candidates[sites[k].ci].unresolved_noisy ||
                                   chunk.candidates[sites[k].ci].low_mapq_only)) join_excluded[k] = 1;
    {
        std::vector<hts_pos_t> site_pos(n);
        for (size_t k = 0; k < n; ++k) site_pos[k] = sites[k].pos;
        for (size_t k = 0; k < n; ++k) {
            if (sites[k].locus < 0) continue;
            const LocusWindowSite& w = (*loci)[static_cast<size_t>(sites[k].locus)];
            const auto lo = std::lower_bound(site_pos.begin(), site_pos.end(), w.pos);
            const auto hi = std::upper_bound(lo, site_pos.end(), std::max(w.pos, w.end));
            for (auto it = lo; it != hi; ++it) {
                const size_t j = static_cast<size_t>(it - site_pos.begin());
                if (sites[j].locus < 0 && site_weak[j]) join_weak[k] = 1;
            }
        }
        std::vector<std::pair<hts_pos_t, hts_pos_t>> spans;
        for (size_t ri = 0; ri < robs.size(); ++ri) {
            spans.clear();
            for (const EmObs& o : robs[ri]) {
                const EmSite& site = sites[static_cast<size_t>(o.s)];
                if (site.locus < 0) continue;
                const LocusWindowSite& w = (*loci)[static_cast<size_t>(site.locus)];
                spans.emplace_back(w.pos, std::max(w.pos, w.end));
            }
            if (spans.empty()) continue;
            join_skip[ri].assign(robs[ri].size(), 0);
            for (size_t i = 0; i < robs[ri].size(); ++i) {
                const EmSite& site = sites[static_cast<size_t>(robs[ri][i].s)];
                if (site.locus >= 0) continue;
                for (const auto& [beg, end] : spans)
                    if (beg <= site.pos && site.pos <= end) join_skip[ri][i] = 1;
            }
        }
    }
    // Site classes for the error prior: clean catalog het, other catalog
    // site, alignment-injected site, locus window -- each split by weak --
    // and calls from noisy regions too long for MSA, pooled apart: they are
    // flagged weak, and sharing the homopolymer calls' pool raised those
    // calls' errors with their own (chr12 lost 232 labels at 37 and 133 Mb).
    std::vector<int> site_class(n);
    for (size_t k = 0; k < n; ++k) {
        int kind = 3;
        if (sites[k].locus < 0) {
            const CandidateVariant& c = chunk.candidates[sites[k].ci];
            kind = c.unresolved_noisy ? 4
                 : c.bam_injected ? 2
                 : (c.lcd_var_i_to_cate & (kCandCleanHetSnp | kCandCleanHetIndel)) != 0 ? 0 : 1;
        }
        site_class[k] = 2 * kind + (site_weak[k] ? 1 : 0);
    }
    bool m_started = false;
    // A site's error is its reads' disagreement rate, shrunk toward the pooled
    // rate of its class by kEmClassPriorReads pseudo-reads. Two sites that
    // share their reads and disagree are otherwise symmetric: whichever the
    // starting labels follow reads error 0.01 and the other ~0.3, so a
    // homopolymer length call could reject a clean SNP (chr12:6,134,662, 35
    // paternal and 24 maternal reads all split by the SNP, rejected at 0.316
    // against a +T/+TT call that does not separate the parents). The class
    // pools say how often each kind of site is wrong in this chunk.
    const auto m_step = [&]() {
        std::vector<std::array<double, 2>> w(n, {0.0, 0.0});  // hap1 weight per side
        for (size_t ri = 0; ri < robs.size(); ++ri) {
            if (!learns[ri]) continue;
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
            const bool single_site = robs[ri].empty() || robs[ri].front().s == robs[ri].back().s;
            if (single_site && site_class[static_cast<size_t>(robs[ri].front().s)] / 2 != 0) continue;
            for (const EmObs& o : robs[ri]) {
                w[static_cast<size_t>(o.s)][static_cast<size_t>(o.b)] += post[ri];
                w[static_cast<size_t>(o.s)][static_cast<size_t>(1 - o.b)] += 1.0 - post[ri];
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
        for (size_t k = 0; k < n; ++k) {
            class_minor[static_cast<size_t>(site_class[k])] += std::min(w[k][0], w[k][1]);
            class_total[static_cast<size_t>(site_class[k])] += w[k][0] + w[k][1];
            if (w[k][0] + w[k][1] >= kEmPoolMinWeight)
                class_rates[static_cast<size_t>(site_class[k])].push_back(std::min(w[k][0], w[k][1]) / (w[k][0] + w[k][1]));
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
        const double shrink = m_started ? kEmClassPriorReads : 0.0;
        for (size_t k = 0; k < n; ++k) {
            const double total = w[k][0] + w[k][1];
            if (total < 1e-6) {
                if (m_started)
                    err[k] = std::min(kEmMaxError, std::max(kEmMinError, class_pool[static_cast<size_t>(site_class[k])]));
                continue;
            }
            phase[k] = w[k][0] >= w[k][1] ? 0 : 1;
            const size_t c = static_cast<size_t>(site_class[k]);
            const double pooled = class_pool[c];
            const double e = (std::min(w[k][0], w[k][1]) + shrink * pooled) / (total + shrink);
            err[k] = std::min(kEmMaxError, std::max(kEmMinError, e));
        }
        m_started = true;
    };
    const auto read_log_likelihoods = [&](const EmObs& o) -> std::array<double, 2> {
        const double e = obs_error(o);
        const bool match1 = o.b == phase[static_cast<size_t>(o.s)];
        return {std::log(match1 ? 1.0 - e : e), std::log(match1 ? e : 1.0 - e)};
    };
    // A call's vote on its read's label: none at an echo site.
    const auto label_log_likelihoods = [&](const EmObs& o) -> std::array<double, 2> {
        if (echo[static_cast<size_t>(o.s)]) return {0.0, 0.0};
        return read_log_likelihoods(o);
    };
    const auto e_step = [&]() {
        for (size_t ri = 0; ri < robs.size(); ++ri) {
            if (robs[ri].empty()) continue;
            double l1 = 0.0, l2 = 0.0;
            for (const EmObs& o : robs[ri]) {
                const auto ll = label_log_likelihoods(o);
                l1 += ll[0];
                l2 += ll[1];
            }
            const double d = l2 - l1;
            post[ri] = d > 700 ? 0.0 : d < -700 ? 1.0 : 1.0 / (1.0 + std::exp(d));
        }
    };
    // delta[k] = log L(current) - log L(every site right of k flipped), over
    // the reads that teach the model or, for block cuts, over every read: a
    // read too ambiguously mapped to shape the model still shows the molecule
    // continues across the boundary.
    // joining: the block-cut test, with each locus counted once (join_skip,
    // join_weak) rather than the switch moves' every observation.
    const auto boundary_delta = [&](bool every_read, bool strong_only = false, bool joining = false) {
        static const std::vector<char> no_skip;
        std::vector<double> diff(n + 1, 0.0);
        std::vector<EmObs> strong;
        for (size_t ri = 0; ri < robs.size(); ++ri) {
            const std::vector<EmObs>* items_ptr = &robs[ri];
            const std::vector<char>& skip = joining ? join_skip[ri] : no_skip;
            const std::vector<char>& weak_of = joining ? join_weak : site_weak;
            if (strong_only || joining) {
                strong.clear();
                for (size_t i = 0; i < robs[ri].size(); ++i) {
                    const EmObs& o = robs[ri][i];
                    if (joining && ((!skip.empty() && skip[i]) || join_excluded[static_cast<size_t>(o.s)])) continue;
                    if (strong_only && weak_of[static_cast<size_t>(o.s)]) continue;
                    strong.push_back(o);
                }
                items_ptr = &strong;
            }
            const auto& items = *items_ptr;
            if ((!every_read && !learns[ri]) || items.size() < 2 || items.front().s == items.back().s)
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
        std::vector<double> delta(n, 0.0);
        double run = 0.0;
        for (size_t k = 0; k < n; ++k) {
            run += diff[k];
            delta[k] = run;
        }
        return delta;
    };
    // EM, then switch moves: flip everything right of the worst boundary while
    // that raises the likelihood.
    std::vector<double> delta;
    // Reads the clean stage left unlabelled start from an E-step under the
    // sites' initial orientation and starting errors, not from 0.5. At 0.5
    // every site's first M-step error is 0.5 -- a saddle -- and the class
    // prior holds an island there: two sites shared by 75 reads at
    // chr18:60,372,121, 97% pure against the parents, stayed at error 0.45
    // with no read labelled.
    {
        for (size_t ri = 0; ri < robs.size(); ++ri) {
            if (robs[ri].empty() || start_label[ri] != 0) continue;
            double l1 = 0.0, l2 = 0.0;
            for (const EmObs& o : robs[ri]) {
                const auto ll = label_log_likelihoods(o);
                l1 += ll[0];
                l2 += ll[1];
            }
            const double d = l2 - l1;
            post[ri] = d > 700 ? 0.0 : d < -700 ? 1.0 : 1.0 / (1.0 + std::exp(d));
        }
    }
    for (int round = 0; round < kEmSwitchRounds; ++round) {
        for (int it = 0; it < kEmIterations; ++it) {
            m_step();
            e_step();
        }
        delta = boundary_delta(false);
        size_t worst = n;
        for (size_t k = 0; k + 1 < n; ++k)
            if (delta[k] < -1e-6 && (worst == n || delta[k] < delta[worst])) worst = k;
        if (worst == n) break;
        for (size_t k = worst + 1; k < n; ++k) phase[k] ^= 1;
        for (size_t ri = 0; ri < robs.size(); ++ri)
            if (!robs[ri].empty() && static_cast<size_t>(robs[ri].front().s) > worst)
                post[ri] = 1.0 - post[ri];
    }
    delta = boundary_delta(true, false, true);
    // A join that the non-weak sites alone do not carry (its margin comes from
    // weak sites: homopolymer length calls) needs a larger margin. Measured on
    // chr20: wrong joins of this kind sit just above the ordinary cut
    // (median 6), right ones well above it (median 19).
    const std::vector<double> delta_strong =
        bridge_weak != nullptr ? boundary_delta(true, true, true) : delta;
    // Reads whose observations lie on both sides of each boundary.
    std::vector<int> spanning(n, 0);
    {
        std::vector<int> diff(n + 1, 0);
        for (const std::vector<EmObs>& obs : robs) {
            if (obs.size() < 2 || obs.front().s == obs.back().s) continue;
            ++diff[static_cast<size_t>(obs.front().s)];
            --diff[static_cast<size_t>(obs.back().s)];
        }
        int run = 0;
        for (size_t k = 0; k < n; ++k) spanning[k] = run += diff[k];
    }
    // Reads whose strong evidence -- reliable, join-eligible, non-weak sites --
    // lies on both sides of each boundary. The strong margin is only as
    // independent as the reads behind it (see kEmMinStrongJoinReads).
    std::vector<int> spanning_strong(n, 0);
    if (bridge_weak != nullptr) {
        std::vector<int> diff(n + 1, 0);
        for (const std::vector<EmObs>& obs : robs) {
            int first = -1, last = -1;
            for (const EmObs& o : obs) {
                const size_t k = static_cast<size_t>(o.s);
                if (err[k] >= kEmUnreliableError || join_excluded[k] || join_weak[k]) continue;
                if (first < 0) first = o.s;
                last = o.s;
            }
            if (first < 0 || first == last) continue;
            ++diff[static_cast<size_t>(first)];
            --diff[static_cast<size_t>(last)];
        }
        int run = 0;
        for (size_t k = 0; k < n; ++k) spanning_strong[k] = run += diff[k];
    }
    // Blocks: cut where flipping the rest would cost less than the threshold,
    // or where too few reads cross.
    std::vector<hts_pos_t> block(n);
    std::vector<char> weak_bridge(n, 0);
    hts_pos_t current = sites[0].pos;
    size_t blocks = 1;
    for (size_t k = 0; k < n; ++k) {
        block[k] = current;
        weak_bridge[k] = delta_strong[k] < kEmBlockCutLogLikelihood ||
                         (bridge_weak != nullptr && spanning_strong[k] < kEmMinStrongJoinReads);
        if (k + 1 < n && (delta[k] < kEmBlockCutLogLikelihood || spanning[k] < kEmMinJoinReads ||
                          (weak_bridge[k] && delta[k] < kEmWeakBridgeCutLogLikelihood))) {
            current = sites[k + 1].pos;
            ++blocks;
        }
    }
    for (CandidateVariant& c : chunk.candidates)
        if ((c.lcd_var_i_to_cate & kCandGermlineVarCate) != 0) c.phase_set = kUnsetCandidatePhaseSet;
    for (size_t k = 0; k < n; ++k) {
        if (sites[k].locus >= 0) {
            const LocusWindowSite& w = (*loci)[static_cast<size_t>(sites[k].locus)];
            w.phase = phase[k];
            w.block = block[k];
            w.error = err[k];
            continue;
        }
        CandidateVariant& c = chunk.candidates[sites[k].ci];
        if (err[k] >= (join_excluded[k] ? kEmUnverifiedReportError : kEmUnreliableError)) {
            c.hap_to_cons_alle[1] = c.hap_to_cons_alle[2] = -1;
            c.hap_alt = c.hap_ref = 0;
            continue;
        }
        const int hap1 = phase[k] == 0 ? sites[k].x : sites[k].y;
        const int hap2 = phase[k] == 0 ? sites[k].y : sites[k].x;
        set_site_phase(c, block[k], hap1, hap2);
    }
    // Weak window calls per read (site, side), for the last-resort labels.
    std::vector<std::vector<std::pair<size_t, int>>> weak_window_obs(robs.size());
    for (size_t li = 0; li < n_locus; ++li) {
        if (locus_site[li] < 0) continue;
        for (const auto& [ri, side] : (*loci)[li].weak_observations)
            if (ri < weak_window_obs.size()) weak_window_obs[ri].emplace_back(static_cast<size_t>(locus_site[li]), side);
    }
    // Reads: label in the block where their own evidence is strongest.
    for (size_t ri = 0; ri < robs.size() && ri < chunk.haps.size(); ++ri) {
        chunk.haps[ri] = 0;
        chunk.phase_sets[ri] = kUnphasedReadPhaseSet;
        std::map<hts_pos_t, std::array<int, 2>> directions;  // calls favouring hap 1 / hap 2
        std::map<hts_pos_t, std::array<double, 2>> per_block;
        for (const EmObs& o : robs[ri]) {
            if (err[static_cast<size_t>(o.s)] >= kEmUnreliableError) continue;
            // Unverified calls (join_excluded) are phased for the VCF but do
            // not label reads: on chr12, chr18 and chr20 leaving them out took
            // 90, 125 and 12 reads off the discordant count and added correct
            // ones on chr12 and chr18.
            if (join_excluded[static_cast<size_t>(o.s)]) continue;
            const auto ll = label_log_likelihoods(o);
            auto& l = per_block[block[static_cast<size_t>(o.s)]];
            l[0] += ll[0];
            l[1] += ll[1];
            if (ll[0] != ll[1]) ++directions[block[static_cast<size_t>(o.s)]][ll[0] > ll[1] ? 0 : 1];
        }
        double best = 0.0;
        for (const auto& [ps, l] : per_block) {
            const double margin = l[0] - l[1];
            const double p1 = margin > 700 ? 1.0 : margin < -700 ? 0.0 : 1.0 / (1.0 + std::exp(-margin));
            if (std::max(p1, 1.0 - p1) < kEmReadLabelPosterior || std::abs(margin) <= best) continue;
            best = std::abs(margin);
            chunk.haps[ri] = p1 >= 0.5 ? 1 : 2;
            chunk.phase_sets[ri] = ps;
        }
        // Deferred: a read below the label posterior whose calls in a block all
        // point one way takes that block's label last (apply_island_labels,
        // fill-only), so it steers no window, stitch or rescue; failing that,
        // its unverified calls may label it (kEmUnverifiedLabelPosterior).
        if (chunk.haps[ri] != 0 || loci == nullptr) continue;
        double weak_best = 0.0;
        int weak_hap = 0;
        hts_pos_t weak_ps = 0;
        for (const auto& [ps, l] : per_block) {
            const double margin = l[0] - l[1];
            const double p1 = margin > 700 ? 1.0 : margin < -700 ? 0.0 : 1.0 / (1.0 + std::exp(-margin));
            const auto& dir = directions[ps];
            if (dir[0] != 0 && dir[1] != 0) continue;
            if (std::max(p1, 1.0 - p1) < kEmDeferredLabelPosterior || std::abs(margin) <= weak_best) continue;
            weak_best = std::abs(margin);
            weak_hap = p1 >= 0.5 ? 1 : 2;
            weak_ps = ps;
        }
        if (weak_hap == 0) {
            std::map<hts_pos_t, std::array<double, 2>> unverified;
            for (const EmObs& o : robs[ri]) {
                const size_t k = static_cast<size_t>(o.s);
                if (!join_excluded[k] || err[k] >= kEmUnverifiedReportError) continue;
                const auto ll = label_log_likelihoods(o);
                auto& l = unverified[block[k]];
                l[0] += ll[0];
                l[1] += ll[1];
            }
            for (const auto& [ps, l] : unverified) {
                const double margin = l[0] - l[1];
                const double p1 = margin > 700 ? 1.0 : margin < -700 ? 0.0 : 1.0 / (1.0 + std::exp(-margin));
                if (std::max(p1, 1.0 - p1) < kEmUnverifiedLabelPosterior || std::abs(margin) <= weak_best) continue;
                weak_best = std::abs(margin);
                weak_hap = p1 >= 0.5 ? 1 : 2;
                weak_ps = ps;
            }
        }
        // Failing those, its one-edit window calls (weak_observations) at
        // reliable windows, all agreeing within a block: +25 correct / +4
        // discordant / -29 unphased over chr12/18/20, switches unchanged.
        // Refining the windows' haplotype medoids from their assigned reads
        // first cost chr20 67 correct reads and 15% of the run time.
        if (weak_hap == 0 && ri < weak_window_obs.size() && !weak_window_obs[ri].empty()) {
            std::map<hts_pos_t, std::array<int, 2>> wv;
            for (const auto& [k, side] : weak_window_obs[ri]) {
                if (err[k] >= kEmUnreliableError) continue;
                ++wv[block[k]][side == phase[k] ? 0 : 1];
            }
            int best = 0;
            for (const auto& [ps, v] : wv) {
                if (v[0] != 0 && v[1] != 0) continue;
                if (v[0] + v[1] <= best) continue;
                best = v[0] + v[1];
                weak_hap = v[0] > 0 ? 1 : 2;
                weak_ps = ps;
            }
        }
        if (weak_hap != 0) {
            if (chunk.island_haps.size() < chunk.reads.size()) {
                chunk.island_haps.resize(chunk.reads.size(), 0);
                chunk.island_phase_sets.resize(chunk.reads.size(), 0);
            }
            chunk.island_haps[ri] = weak_hap;
            chunk.island_phase_sets[ri] = weak_ps;
        }
    }
    // Island blocks. A read the solve leaves unlabelled is phased locally by
    // its own call at one site that cannot vote: a block of that site alone,
    // never joined to anything, so a wrong island costs its own reads and
    // never flips a neighbour. Candidates: sites the EM rejects (learned error
    // in [kEmUnreliableError, kIslandMaxError)), repeat het indels it never
    // sees (REP_HET_INDEL) and MSA-verified hets injection kept out
    // (island_only), the last two judged by their raw disagreement with the
    // labelled reads, and island-only length windows (repeat loci genotyped
    // from read lengths; ungated). Each read takes its lowest-error candidate. Recorded
    // here, applied after the rescue to reads still unphased
    // (apply_island_labels). Measured on chr12/18/20 against c8391f1:
    // +2,349 correct, +457 discordant, -2,806 unphased reads, one arm switch;
    // islands at truth hets run 18-20 : 1, at repeat-length noise 1.5 : 1.
    if (loci != nullptr) {
        struct IslandCand { size_t ci; int x; int y; hts_pos_t pos; int k; double err; int phase; char kind; };
        std::vector<IslandCand> cands;
        std::vector<int> cand_of_site(n, -1);
        for (size_t k = 0; k < n; ++k) {
            if (sites[k].locus >= 0 || join_excluded[k] || err[k] < kEmUnreliableError || err[k] >= kIslandMaxError) continue;
            cand_of_site[k] = static_cast<int>(cands.size());
            cands.push_back(IslandCand{sites[k].ci, sites[k].x, sites[k].y, sites[k].pos, static_cast<int>(k), err[k], phase[k], 'E'});
        }
        std::vector<int> cand_of_ci(chunk.candidates.size(), -1);
        for (size_t ci = 0; ci < chunk.candidates.size(); ++ci) {
            const CandidateVariant& c = chunk.candidates[ci];
            if ((c.counts.category != VariantCategory::RepeatHetIndel && !c.island_only) || site_of[ci] >= 0) continue;
            std::vector<std::pair<int, int>> by_count;
            for (const auto& [al, cnt] : allele_counts[ci]) by_count.emplace_back(cnt, al);
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
        std::vector<int> cand_of_locus(n_locus, -1);
        for (size_t li = 0; li < n_locus; ++li) {
            if (!(*loci)[li].island_only) continue;
            cand_of_locus[li] = static_cast<int>(cands.size());
            cands.push_back(IslandCand{SIZE_MAX, 0, 1, (*loci)[li].pos, -1, 1.0, 0, 'W'});
        }
        // A site outside the solve: its raw disagreement with the labelled
        // reads, block by block; too few labelled reads to judge counts as the
        // unreliable cut.
        {
            std::vector<std::map<hts_pos_t, std::array<int, 2>>> agree(cands.size());
            for (size_t ri = 0; ri < chunk.read_var_profile.size() && ri < chunk.haps.size(); ++ri) {
                const int hap = chunk.haps[ri];
                if (hap != 1 && hap != 2) continue;
                const ReadVariantProfile& prof = chunk.read_var_profile[ri];
                if (prof.start_var_idx < 0 || chunk.reads[ri].is_skipped) continue;
                for (size_t j = 0; j < prof.alleles.size(); ++j) {
                    const int ic = cand_of_ci[static_cast<size_t>(prof.start_var_idx) + j];
                    if (ic < 0) continue;
                    const int al = prof.alleles[j];
                    const IslandCand& cand = cands[static_cast<size_t>(ic)];
                    if (al != cand.x && al != cand.y) continue;
                    ++agree[static_cast<size_t>(ic)][chunk.phase_sets[ri]][((al == cand.x) == (hap == 1)) ? 0 : 1];
                }
            }
            for (size_t li = 0; li < n_locus; ++li) {
                const int ic = cand_of_locus[li];
                if (ic < 0) continue;
                for (const auto& [ri, side, e] : (*loci)[li].observations) {
                    if (ri >= chunk.haps.size()) continue;
                    const int hap = chunk.haps[ri];
                    if (hap != 1 && hap != 2) continue;
                    ++agree[static_cast<size_t>(ic)][chunk.phase_sets[ri]][((side == 0) == (hap == 1)) ? 0 : 1];
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
        std::vector<std::vector<std::pair<int, int>>> iobs(robs.size());
        for (size_t ri = 0; ri < robs.size() && ri < chunk.haps.size(); ++ri) {
            if (chunk.haps[ri] != 0) continue;
            if (ri < chunk.island_haps.size() && chunk.island_haps[ri] != 0) continue;  // deferred block label
            for (const EmObs& o : robs[ri])
                if (cand_of_site[static_cast<size_t>(o.s)] >= 0) iobs[ri].emplace_back(cand_of_site[static_cast<size_t>(o.s)], o.b);
            if (ri >= chunk.read_var_profile.size()) continue;
            const ReadVariantProfile& prof = chunk.read_var_profile[ri];
            if (prof.start_var_idx < 0 || chunk.reads[ri].is_skipped) continue;
            for (size_t j = 0; j < prof.alleles.size(); ++j) {
                const int ic = cand_of_ci[static_cast<size_t>(prof.start_var_idx) + j];
                if (ic < 0 || !usable(static_cast<size_t>(ic))) continue;
                const int al = prof.alleles[j];
                if (al == cands[static_cast<size_t>(ic)].x) iobs[ri].emplace_back(ic, 0);
                else if (al == cands[static_cast<size_t>(ic)].y) iobs[ri].emplace_back(ic, 1);
            }
        }
        for (size_t li = 0; li < n_locus; ++li) {
            const int ic = cand_of_locus[li];
            if (ic < 0 || !usable(static_cast<size_t>(ic))) continue;
            for (const auto& [ri, side, e] : (*loci)[li].observations) {
                if (ri >= iobs.size() || ri >= chunk.haps.size() || chunk.haps[ri] != 0) continue;
                if (ri < chunk.island_haps.size() && chunk.island_haps[ri] != 0) continue;
                iobs[ri].emplace_back(ic, side);
            }
        }
        std::vector<int> unlabelled_callers(cands.size(), 0);
        for (const auto& v : iobs)
            for (const auto& [ic, side] : v) ++unlabelled_callers[static_cast<size_t>(ic)];
        std::vector<int> pick_cand(robs.size(), -1), pick_side(robs.size(), 0);
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
        std::set<hts_pos_t> used_ps(block.begin(), block.end());
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
        for (size_t ri = 0; ri < robs.size(); ++ri) {
            if (pick_cand[ri] < 0) continue;
            const size_t ic = static_cast<size_t>(pick_cand[ri]);
            if (chunk.island_haps.size() < chunk.reads.size()) {
                chunk.island_haps.resize(chunk.reads.size(), 0);
                chunk.island_phase_sets.resize(chunk.reads.size(), 0);
            }
            chunk.island_haps[ri] = pick_side[ri] == cands[ic].phase ? 1 : 2;
            chunk.island_phase_sets[ri] = ps_of(ic);
        }
        chunk.island_sites.clear();
        for (const auto& [ic, ps] : island_ps) {
            const IslandCand& c = cands[ic];
            if (c.kind == 'W') continue;  // a window writes no VCF record
            const int hap1 = c.phase == 0 ? c.x : c.y;
            const int hap2 = c.phase == 0 ? c.y : c.x;
            chunk.island_sites.push_back(PhasingChunk::IslandSite{c.ci, ps, hap1, hap2});
        }
        // Diagnostic, off by default: PGPHASE_DIAG=DIR writes every read an
        // island claims (islands.TID.BEG.tsv): read, site position, island PS,
        // the read's side, its island haplotype, its read error at the site,
        // the site's error (learned, or raw outside the solve), unlabelled
        // callers, the most member reads sharing one other verified site, 1,
        // kind (E rejected EM site, R repeat het indel, M island-only MSA het,
        // W length window).
        static const char* const island_diag = std::getenv("PGPHASE_DIAG");
        if (island_diag != nullptr) {
            std::vector<int> linked_reads(cands.size(), 0);
            std::map<std::pair<int, int>, int> shared;  // (island, other EM site) -> member reads calling it
            for (size_t ri = 0; ri < robs.size(); ++ri) {
                if (pick_cand[ri] < 0) continue;
                const int own = cands[static_cast<size_t>(pick_cand[ri])].k;
                for (const EmObs& o : robs[ri]) {
                    const size_t k = static_cast<size_t>(o.s);
                    if (o.s == own || sites[k].locus >= 0 || chunk.candidates[sites[k].ci].unresolved_noisy) continue;
                    int& c = ++shared[{pick_cand[ri], o.s}];
                    int& m = linked_reads[static_cast<size_t>(pick_cand[ri])];
                    m = std::max(m, c);
                }
            }
            const std::string path = std::string(island_diag) + "/islands." + std::to_string(chunk.region.tid) + "." +
                                     std::to_string(chunk.region.beg) + ".tsv";
            if (FILE* fh = std::fopen(path.c_str(), "w")) {
                for (size_t ri = 0; ri < robs.size(); ++ri) {
                    if (pick_cand[ri] < 0) continue;
                    const size_t ic = static_cast<size_t>(pick_cand[ri]);
                    double read_err = 0.0;
                    for (const EmObs& o : robs[ri])
                        if (o.s == cands[ic].k) read_err = o.e;
                    std::fprintf(fh, "%s\t%lld\t%lld\t%d\t%d\t%.4f\t%.3f\t%d\t%d\t1\t%c\n", chunk.reads[ri].qname.c_str(),
                                 static_cast<long long>(cands[ic].pos), static_cast<long long>(island_ps.at(ic)),
                                 pick_side[ri], pick_side[ri] == cands[ic].phase ? 1 : 2, read_err, cands[ic].err,
                                 unlabelled_callers[ic], linked_reads[ic], cands[ic].kind);
                }
                std::fclose(fh);
            }
        }
    }
    // Diagnostic, off by default: PGPHASE_EM_DUMP=DIR writes the final solve's
    // evidence per read, one file per chunk. A line is the read, MAPQ, label
    // (HP, PS), per-block log-likelihood margins (PS=margin), the label the
    // read started the solve with, then each observation as
    // POS:KIND:SIDE:SITE_ERROR:BLOCK:WEAK:READ_ERROR, with KIND C clean het,
    // O other, I injected from the alignment, U injected from a noisy region
    // too long for MSA, L window.
    // Diagnostic, off by default: PGPHASE_DIAG=DIR writes the final solve's
    // voting candidates (emsites.TID.BEG.tsv): VCF-convention position, type,
    // REF length, ALT, origin (C catalog clean, O catalog other, I injected, U
    // unverified), allele read counts (allele:reads;...), reads in range
    // without an allele, EM site or "-" (left out: skewed alleles), the two
    // alleles contrasted, learned error, phase (which allele haplotype 1
    // carries), block.
    static const char* const diag_dir = std::getenv("PGPHASE_DIAG");
    if (diag_dir != nullptr && loci != nullptr) {
        const std::string path = std::string(diag_dir) + "/emsites." + std::to_string(chunk.region.tid) + "." +
                                 std::to_string(chunk.region.beg) + ".tsv";
        if (FILE* fh = std::fopen(path.c_str(), "w")) {
            for (size_t ci = 0; ci < chunk.candidates.size(); ++ci) {
                const CandidateVariant& c = chunk.candidates[ci];
                if ((c.lcd_var_i_to_cate & kCandGermlineVarCate) == 0) continue;
                const char type = c.key.type == VariantType::Snp ? 'S' : c.key.type == VariantType::Insertion ? 'I' : 'D';
                const char origin = c.unresolved_noisy ? 'U' : c.bam_injected ? 'I'
                                  : (c.lcd_var_i_to_cate & (kCandCleanHetSnp | kCandCleanHetIndel)) != 0 ? 'C' : 'O';
                std::string counts;
                for (const auto& [a, n] : allele_counts[ci])
                    counts += (counts.empty() ? "" : ";") + std::to_string(a) + ":" + std::to_string(n);
                std::fprintf(fh, "%lld\t%c\t%d\t%s\t%c\t%s\t%d", static_cast<long long>(c.key.sort_pos()), type,
                             c.key.ref_len, c.key.alt.substr(0, 40).c_str(), origin,
                             counts.empty() ? "." : counts.c_str(), undecided[ci]);
                const int k = site_of[ci];
                if (k < 0) {
                    std::fprintf(fh, "\t-\t.\t.\t.\t.\t.\n");
                } else {
                    const EmSite& site = sites[static_cast<size_t>(k)];
                    std::fprintf(fh, "\tem\t%d\t%d\t%.3f\t%d\t%lld\n", site.x, site.y, err[static_cast<size_t>(k)],
                                 phase[static_cast<size_t>(k)] == 0 ? site.x : site.y,
                                 static_cast<long long>(block[static_cast<size_t>(k)]));
                }
            }
            std::fclose(fh);
        }
        // Every read's profile as the final solve saw it (emreads.TID.BEG.tsv):
        // name, MAPQ, skipped flag, then each call as VCF-convention position
        // of the candidate : allele, including candidates the EM left out.
        const std::string rpath = std::string(diag_dir) + "/emreads." + std::to_string(chunk.region.tid) + "." +
                                  std::to_string(chunk.region.beg) + ".tsv";
        if (FILE* fh = std::fopen(rpath.c_str(), "w")) {
            for (size_t ri = 0; ri < chunk.reads.size(); ++ri) {
                const ReadRecord& r = chunk.reads[ri];
                std::fprintf(fh, "%s\t%d\t%d", r.qname.c_str(), r.mapq, r.is_skipped ? 1 : 0);
                const ReadVariantProfile* prof = nullptr;
                for (const ReadVariantProfile& p2 : chunk.read_var_profile)
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
                    if (ci >= chunk.candidates.size()) break;
                    std::fprintf(fh, "\t%lld:%d", static_cast<long long>(chunk.candidates[ci].key.sort_pos()),
                                 prof->alleles[k]);
                }
                std::fputc('\n', fh);
            }
            std::fclose(fh);
        }
    }
    // Every block-cut decision (joins.TID.BEG.tsv): between each EM site and
    // the next, their VCF-convention positions, the join margin (log L of the
    // solve minus log L with every site to the right flipped, each locus
    // counted once), the margin without weak sites, reads with observations on
    // both sides, the decision -- J joined, C cut (margin below the cut), S cut
    // (too few reads cross), W cut (a weak bridge below its larger cut) -- and
    // reads with strong evidence on both sides.
    if (diag_dir != nullptr && loci != nullptr) {
        const std::string path = std::string(diag_dir) + "/joins." + std::to_string(chunk.region.tid) + "." +
                                 std::to_string(chunk.region.beg) + ".tsv";
        if (FILE* fh = std::fopen(path.c_str(), "w")) {
            for (size_t k = 0; k + 1 < n; ++k) {
                const char decision = delta[k] < kEmBlockCutLogLikelihood ? 'C'
                                    : spanning[k] < kEmMinJoinReads ? 'S'
                                    : weak_bridge[k] && delta[k] < kEmWeakBridgeCutLogLikelihood ? 'W' : 'J';
                std::fprintf(fh, "%lld\t%lld\t%.2f\t%.2f\t%d\t%c\t%d\n", static_cast<long long>(sites[k].pos),
                             static_cast<long long>(sites[k + 1].pos), delta[k], delta_strong[k], spanning[k], decision,
                             spanning_strong[k]);
            }
            std::fclose(fh);
        }
    }
    static const char* const dump_dir = std::getenv("PGPHASE_EM_DUMP");
    if (dump_dir != nullptr && loci != nullptr) {
        const std::string path = std::string(dump_dir) + "/em." + std::to_string(chunk.region.tid) + "." +
                                 std::to_string(chunk.region.beg) + ".tsv";
        if (FILE* fh = std::fopen(path.c_str(), "w")) {
            for (size_t ri = 0; ri < chunk.reads.size() && ri < robs.size(); ++ri) {
                const ReadRecord& r = chunk.reads[ri];
                std::fprintf(fh, "%s\t%d\t%d\t%lld\t", r.qname.c_str(), r.mapq, chunk.haps[ri],
                             static_cast<long long>(chunk.phase_sets[ri]));
                std::map<hts_pos_t, double> margin;
                for (const EmObs& o : robs[ri]) {
                    if (err[static_cast<size_t>(o.s)] >= kEmUnreliableError) continue;
                    if (join_excluded[static_cast<size_t>(o.s)]) continue;
                    const auto ll = label_log_likelihoods(o);
                    margin[block[static_cast<size_t>(o.s)]] += ll[0] - ll[1];
                }
                const char* sep = "";
                for (const auto& [ps, m] : margin) {
                    std::fprintf(fh, "%s%lld=%.2f", sep, static_cast<long long>(ps), m);
                    sep = ",";
                }
                std::fprintf(fh, "\t%d", start_label[ri]);
                for (const EmObs& o : robs[ri]) {
                    const EmSite& site = sites[static_cast<size_t>(o.s)];
                    char kind = 'L';
                    if (site.locus < 0) {
                        const CandidateVariant& c = chunk.candidates[site.ci];
                        kind = c.unresolved_noisy ? 'U'
                             : c.bam_injected ? 'I'
                             : (c.lcd_var_i_to_cate & (kCandCleanHetSnp | kCandCleanHetIndel)) != 0 ? 'C' : 'O';
                    }
                    std::fprintf(fh, "\t%lld:%c:%d:%.3f:%lld:%d:%.4f", static_cast<long long>(site.pos), kind, o.b,
                                 err[static_cast<size_t>(o.s)], static_cast<long long>(block[static_cast<size_t>(o.s)]),
                                 static_cast<int>(site_weak[static_cast<size_t>(o.s)]), o.e);
                }
                std::fputc('\n', fh);
            }
            std::fclose(fh);
        }
    }
    return blocks;
}

// ── Local haplotype windows ─────────────────────────────────────────────────

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
