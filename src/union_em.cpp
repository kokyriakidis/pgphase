// Union gap phasing: the global EM over the read x site matrix, with block
// cuts and the weak-bridge guard.

#include "union_internal.hpp"

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
constexpr int kEmIterations = 15;
constexpr int kEmSwitchRounds = 30;
constexpr double kEmReadLabelPosterior = 0.8;
constexpr double kEmBlockCutLogLikelihood = 4.0;
// The cut for joins carried only by weak sites (see bridge_weak_sites).
constexpr double kEmWeakBridgeCutLogLikelihood = 10.0;
constexpr int kEmLearnMinMapq = 20;
// A site whose less observed allele carries under this fraction of its reads
// does not take part: a skewed split is the signature of a paralog or error
// call, and in a homozygous stretch such a site labels reads at random.
constexpr double kEmMinMinorAlleleFraction = 0.25;

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
                                const std::vector<char>* bridge_weak) {
    struct EmSite { long locus; size_t ci; int x; int y; hts_pos_t pos; };  // locus < 0: candidate
    struct EmObs { int s; int b; double e; };
    // Sites in position order with their two-allele contrast; a locus window
    // contrasts its two local haplotype sequences (0 and 1).
    std::vector<std::map<int, int>> allele_counts(chunk.candidates.size());
    for (size_t ri = 0; ri < chunk.read_var_profile.size(); ++ri) {
        const ReadVariantProfile& prof = chunk.read_var_profile[ri];
        if (prof.start_var_idx < 0 || chunk.reads[ri].is_skipped) continue;
        for (size_t k = 0; k < prof.alleles.size(); ++k)
            if (prof.alleles[k] >= 0)
                ++allele_counts[static_cast<size_t>(prof.start_var_idx) + k][prof.alleles[k]];
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
        sites.push_back(EmSite{static_cast<long>(li), 0, 0, 1, (*loci)[li].pos});
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
            if (prof.alleles[k] == site.x) robs[ri].push_back(EmObs{s, 0, 0.0});
            else if (prof.alleles[k] == site.y) robs[ri].push_back(EmObs{s, 1, 0.0});
        }
    }
    for (size_t li = 0; li < n_locus; ++li)
        for (const auto& [ri, side, e] : (*loci)[li].observations)
            if (ri < robs.size()) robs[ri].push_back(EmObs{locus_site[li], side, static_cast<double>(e)});
    std::vector<double> post(chunk.reads.size(), 0.5);
    std::vector<char> learns(chunk.reads.size(), 0);
    for (size_t ri = 0; ri < robs.size(); ++ri) {
        std::sort(robs[ri].begin(), robs[ri].end(),
                  [](const EmObs& a, const EmObs& b) { return a.s < b.s; });
        const int hap = ri < chunk.haps.size() ? chunk.haps[ri] : 0;
        post[ri] = hap == 1 ? 0.95 : hap == 2 ? 0.05 : 0.5;
        learns[ri] = chunk.reads[ri].mapq >= kEmLearnMinMapq;
    }
    const auto obs_error = [&](const EmObs& o) {
        return std::min(kEmMaxError, 1.0 - (1.0 - err[static_cast<size_t>(o.s)]) * (1.0 - o.e));
    };
    const auto m_step = [&]() {
        std::vector<std::array<double, 2>> w(n, {0.0, 0.0});  // hap1 weight per side
        for (size_t ri = 0; ri < robs.size(); ++ri) {
            if (!learns[ri]) continue;
            for (const EmObs& o : robs[ri]) {
                w[static_cast<size_t>(o.s)][static_cast<size_t>(o.b)] += post[ri];
                w[static_cast<size_t>(o.s)][static_cast<size_t>(1 - o.b)] += 1.0 - post[ri];
            }
        }
        for (size_t k = 0; k < n; ++k) {
            const double total = w[k][0] + w[k][1];
            if (total < 1e-6) continue;
            phase[k] = w[k][0] >= w[k][1] ? 0 : 1;
            err[k] = std::min(kEmMaxError, std::max(kEmMinError, std::min(w[k][0], w[k][1]) / total));
        }
    };
    const auto read_log_likelihoods = [&](const EmObs& o) -> std::array<double, 2> {
        const double e = obs_error(o);
        const bool match1 = o.b == phase[static_cast<size_t>(o.s)];
        return {std::log(match1 ? 1.0 - e : e), std::log(match1 ? e : 1.0 - e)};
    };
    const auto e_step = [&]() {
        for (size_t ri = 0; ri < robs.size(); ++ri) {
            if (robs[ri].empty()) continue;
            double l1 = 0.0, l2 = 0.0;
            for (const EmObs& o : robs[ri]) {
                const auto ll = read_log_likelihoods(o);
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
    // Weak sites (see bridge_weak_sites) are phased but cannot hold blocks
    // together: the final cuts are judged on the other sites' evidence.
    std::vector<char> site_weak(n, 0);
    if (bridge_weak != nullptr)
        for (size_t k = 0; k < n; ++k)
            if (sites[k].locus < 0 && sites[k].ci < bridge_weak->size()) site_weak[k] = (*bridge_weak)[sites[k].ci];
    const auto boundary_delta = [&](bool every_read, bool strong_only = false) {
        std::vector<double> diff(n + 1, 0.0);
        std::vector<EmObs> strong;
        for (size_t ri = 0; ri < robs.size(); ++ri) {
            const std::vector<EmObs>* items_ptr = &robs[ri];
            if (strong_only) {
                strong.clear();
                for (const EmObs& o : robs[ri])
                    if (!site_weak[static_cast<size_t>(o.s)]) strong.push_back(o);
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
    delta = boundary_delta(true);
    // A join that the non-weak sites alone do not carry (its margin comes from
    // weak sites: homopolymer length calls) needs a larger margin. Measured on
    // chr20: wrong joins of this kind sit just above the ordinary cut
    // (median 6), right ones well above it (median 19).
    const std::vector<double> delta_strong =
        bridge_weak != nullptr ? boundary_delta(true, true) : delta;
    // Blocks: cut where flipping the rest would cost less than the threshold.
    std::vector<hts_pos_t> block(n);
    hts_pos_t current = sites[0].pos;
    size_t blocks = 1;
    for (size_t k = 0; k < n; ++k) {
        block[k] = current;
        const bool weak_bridge = delta_strong[k] < kEmBlockCutLogLikelihood;
        if (k + 1 < n && (delta[k] < kEmBlockCutLogLikelihood ||
                          (weak_bridge && delta[k] < kEmWeakBridgeCutLogLikelihood))) {
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
        if (err[k] >= kEmUnreliableError) {
            c.hap_to_cons_alle[1] = c.hap_to_cons_alle[2] = -1;
            c.hap_alt = c.hap_ref = 0;
            continue;
        }
        const int hap1 = phase[k] == 0 ? sites[k].x : sites[k].y;
        const int hap2 = phase[k] == 0 ? sites[k].y : sites[k].x;
        set_site_phase(c, block[k], hap1, hap2);
    }
    // Reads: label in the block where their own evidence is strongest.
    for (size_t ri = 0; ri < robs.size() && ri < chunk.haps.size(); ++ri) {
        chunk.haps[ri] = 0;
        chunk.phase_sets[ri] = kUnphasedReadPhaseSet;
        std::map<hts_pos_t, std::array<double, 2>> per_block;
        for (const EmObs& o : robs[ri]) {
            if (err[static_cast<size_t>(o.s)] >= kEmUnreliableError) continue;
            const auto ll = read_log_likelihoods(o);
            auto& l = per_block[block[static_cast<size_t>(o.s)]];
            l[0] += ll[0];
            l[1] += ll[1];
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
    }
    // Diagnostic, off by default: PGPHASE_EM_DUMP=DIR writes the final solve's
    // evidence per read, one file per chunk. A line is the read, MAPQ, label
    // (HP, PS), per-block log-likelihood margins (PS=margin), then each
    // observation as POS:KIND:SIDE:SITE_ERROR:BLOCK, with KIND C clean het,
    // O other, I injected from the alignment, L window.
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
                    const auto ll = read_log_likelihoods(o);
                    margin[block[static_cast<size_t>(o.s)]] += ll[0] - ll[1];
                }
                const char* sep = "";
                for (const auto& [ps, m] : margin) {
                    std::fprintf(fh, "%s%lld=%.2f", sep, static_cast<long long>(ps), m);
                    sep = ",";
                }
                for (const EmObs& o : robs[ri]) {
                    const EmSite& site = sites[static_cast<size_t>(o.s)];
                    char kind = 'L';
                    if (site.locus < 0) {
                        const CandidateVariant& c = chunk.candidates[site.ci];
                        kind = c.bam_injected ? 'I'
                             : (c.lcd_var_i_to_cate & (kCandCleanHetSnp | kCandCleanHetIndel)) != 0 ? 'C' : 'O';
                    }
                    std::fprintf(fh, "\t%lld:%c:%d:%.3f:%lld", static_cast<long long>(site.pos), kind, o.b,
                                 err[static_cast<size_t>(o.s)], static_cast<long long>(block[static_cast<size_t>(o.s)]));
                }
                std::fputc('\n', fh);
            }
            std::fclose(fh);
        }
    }
    return blocks;
}

// ── Local haplotype windows ─────────────────────────────────────────────────

} // namespace pgphase_collect
