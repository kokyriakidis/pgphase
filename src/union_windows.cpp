// Union gap phasing: haplotype (medoid) windows and allele windows over noisy
// loci, as extra EM sites.

#include <unordered_map>
#include <string_view>
#include "union_internal.hpp"
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <map>

#include "edlib.h"

namespace pgphase_collect {

// Query index aligned at (or, inside a deletion, after) a 0-based reference
// position; -1 outside the alignment (the walk itself ends at the read's end).
hts_pos_t locus_query_index(const bam1_t* b, hts_pos_t ref_pos) {
    if (ref_pos < b->core.pos) return -1;
    hts_pos_t r = b->core.pos, q = 0;
    const uint32_t* cigar = bam_get_cigar(b);
    for (uint32_t i = 0; i < b->core.n_cigar; ++i) {
        const int op = bam_cigar_op(cigar[i]);
        const hts_pos_t len = bam_cigar_oplen(cigar[i]);
        if (op == BAM_CMATCH || op == BAM_CEQUAL || op == BAM_CDIFF) {
            if (r + len > ref_pos) return q + (ref_pos - r);
            r += len;
            q += len;
        } else if (op == BAM_CDEL || op == BAM_CREF_SKIP) {
            if (r + len > ref_pos) return q;
            r += len;
        } else if (op == BAM_CINS || op == BAM_CSOFT_CLIP) {
            q += len;
        }
    }
    return -1;
}

int locus_edit_distance(const std::string& a, const std::string& b) {
    EdlibAlignResult r = edlibAlign(a.data(), static_cast<int>(a.size()), b.data(),
                                    static_cast<int>(b.size()),
                                    edlibNewAlignConfig(-1, EDLIB_MODE_NW, EDLIB_TASK_DISTANCE, nullptr, 0));
    const int d = r.editDistance;
    edlibFreeAlignResult(r);
    return d;
}

// The sequence with the least total edit distance to the others; the first
// such in input order. Identical sequences are aligned once, each pair once.
static const std::string* locus_medoid(const std::vector<const std::string*>& seqs) {
    std::vector<const std::string*> uniq;
    std::vector<long> count;
    for (const std::string* s : seqs) {
        size_t u = 0;
        while (u < uniq.size() && *uniq[u] != *s) ++u;
        if (u == uniq.size()) {
            uniq.push_back(s);
            count.push_back(0);
        }
        ++count[u];
    }
    std::vector<long> cost(uniq.size(), 0);
    for (size_t i = 0; i < uniq.size(); ++i)
        for (size_t j = i + 1; j < uniq.size(); ++j) {
            const long d = locus_edit_distance(*uniq[i], *uniq[j]);
            cost[i] += d * count[j];
            cost[j] += d * count[i];
        }
    size_t best = 0;
    for (size_t i = 1; i < uniq.size(); ++i)
        if (cost[i] < cost[best]) best = i;
    return uniq.empty() ? nullptr : uniq[best];
}

// (haplotype, phase set) per read from the chunk's phased sites.
static std::vector<std::pair<int, hts_pos_t>> site_derived_read_labels(const PhasingChunk& chunk) {
    std::vector<std::pair<int, hts_pos_t>> labels(chunk.reads.size(), {0, kUnphasedReadPhaseSet});
    for (size_t ri = 0; ri < chunk.read_var_profile.size() && ri < chunk.reads.size(); ++ri) {
        const ReadVariantProfile& prof = chunk.read_var_profile[ri];
        if (prof.start_var_idx < 0 || chunk.reads[ri].is_skipped) continue;
        std::map<hts_pos_t, std::array<int, 2>> votes;
        for (size_t k = 0; k < prof.alleles.size(); ++k) {
            const CandidateVariant& c = chunk.candidates[static_cast<size_t>(prof.start_var_idx) + k];
            const int a = prof.alleles[k];
            if (a < 0 || c.phase_set <= 0 || c.hap_to_cons_alle[1] < 0 ||
                c.hap_to_cons_alle[1] == c.hap_to_cons_alle[2]) continue;
            if (a == c.hap_to_cons_alle[1]) ++votes[c.phase_set][0];
            else if (a == c.hap_to_cons_alle[2]) ++votes[c.phase_set][1];
        }
        int best_margin = 0;
        for (const auto& [ps, v] : votes) {
            const int margin = std::abs(v[0] - v[1]);
            if (margin < kSeedLabelMinMargin ||
                std::max(v[0], v[1]) < kSeedLabelMinAgreement * (v[0] + v[1])) continue;
            if (margin > best_margin) {
                best_margin = margin;
                labels[ri] = {v[0] > v[1] ? 1 : 2, ps};
            }
        }
    }
    return labels;
}

std::vector<LocusWindowSite> build_locus_window_sites(const PhasingChunk& bam,
                                                      const PhasingChunk& graph) {
    std::vector<LocusWindowSite> out;
    if (bam.ref_seq.empty()) return out;
    const std::vector<std::pair<int, hts_pos_t>> seed_labels = site_derived_read_labels(graph);
    const ReferenceView ref{bam};
    // Loci: the alignment's noisy-region calls, merged.
    std::vector<std::pair<hts_pos_t, hts_pos_t>> calls;
    for (const CandidateVariant& c : bam.candidates) {
        if (c.counts.category != VariantCategory::NoisyCandHet &&
            c.counts.category != VariantCategory::NoisyCandHom) continue;
        const hts_pos_t p0 = c.key.pos - 1;
        const hts_pos_t len = std::max<hts_pos_t>(
            1, std::max<hts_pos_t>(c.key.ref_len, static_cast<hts_pos_t>(c.key.alt.size())));
        calls.emplace_back(p0, p0 + len);
    }
    std::sort(calls.begin(), calls.end());
    std::vector<std::pair<hts_pos_t, hts_pos_t>> merged;
    for (const auto& [a, b] : calls) {
        if (!merged.empty() && a - merged.back().second <= kLocusMergeGap)
            merged.back().second = std::max(merged.back().second, b);
        else
            merged.emplace_back(a, b);
    }
    const SpanningReadIndex index(bam, graph);
    // Diagnostic, off by default: PGPHASE_DIAG=DIR writes each locus and what
    // became of it (windows.TID.BEG.tsv): window start, end, spanning reads,
    // seeded reads per side, haplotype edit distance, reads assigned per side,
    // decision.
    static const char* const diag_dir = std::getenv("PGPHASE_DIAG");
    FILE* diag = nullptr;
    if (diag_dir != nullptr)
        diag = std::fopen((std::string(diag_dir) + "/windows." + std::to_string(graph.region.tid) + "." +
                           std::to_string(graph.region.beg) + ".tsv").c_str(), "w");
    const auto note = [&](hts_pos_t l0, hts_pos_t l1, size_t nseg, size_t s1, size_t s2, int dh, int n0, int n1,
                          const char* why) {
        if (diag != nullptr)
            std::fprintf(diag, "%lld\t%lld\t%zu\t%zu\t%zu\t%d\t%d\t%d\t%s\n", static_cast<long long>(l0 + 1),
                         static_cast<long long>(l1), nseg, s1, s2, dh, n0, n1, why);
    };
    for (const auto& [a, b] : merged) {
        const hts_pos_t l0 = std::max<hts_pos_t>(0, a - kLocusPad);
        const hts_pos_t l1 = b + ref.repeat_right(b) + kLocusPad;
        if (l1 - l0 > kLocusMaxLength) { note(l0, l1, 0, 0, 0, -1, 0, 0, "too_long"); continue; }
        const std::vector<std::pair<size_t, std::string>> segs = index.segments(l0, l1, kLocusMaxReads);
        if (segs.size() < static_cast<size_t>(2 * kLocusMinSide)) { note(l0, l1, segs.size(), 0, 0, -1, 0, 0, "few_reads"); continue; }
        // Seeds: reads already labelled in the dominant phase set here.
        std::map<hts_pos_t, int> ps_count;
        for (const auto& [ri, seq] : segs)
            if (seed_labels[ri].first != 0) ++ps_count[seed_labels[ri].second];
        std::vector<const std::string*> side1, side2;
        if (!ps_count.empty()) {
            const hts_pos_t ps = std::max_element(ps_count.begin(), ps_count.end(),
                [](const auto& x, const auto& y) { return x.second < y.second; })->first;
            for (const auto& [ri, seq] : segs) {
                if (seed_labels[ri].second != ps) continue;
                if (seed_labels[ri].first == 1 && side1.size() < kLocusSeedSample) side1.push_back(&seq);
                if (seed_labels[ri].first == 2 && side2.size() < kLocusSeedSample) side2.push_back(&seq);
            }
        }
        const std::string* hap1 = nullptr;
        const std::string* hap2 = nullptr;
        if (side1.size() >= 2 && side2.size() >= 2) {
            hap1 = locus_medoid(side1);
            hap2 = locus_medoid(side2);
        }
        // Seeds can all come from one haplotype (the other's reads unlabelled
        // here); identical seeded medoids then say nothing, and the reads
        // themselves are split instead.
        if (hap1 == nullptr || *hap1 == *hap2) {
            hap1 = hap2 = nullptr;
            // No usable seeds: split the spanning reads themselves. The EM
            // orients the locus through the reads it shares with other sites.
            std::vector<const std::string*> pool;
            for (const auto& [ri, seq] : segs)
                if (pool.size() < 2 * kLocusSeedSample) pool.push_back(&seq);
            hap1 = locus_medoid(pool);
            int far = -1;
            for (const std::string* seq : pool) {
                const int d = locus_edit_distance(*seq, *hap1);
                if (d > far) { far = d; hap2 = seq; }
            }
            for (int round = 0; round < kLocusSplitRounds && hap2 != nullptr; ++round) {
                std::vector<const std::string*> a, b;
                for (const std::string* seq : pool)
                    (locus_edit_distance(*seq, *hap1) <= locus_edit_distance(*seq, *hap2) ? a : b).push_back(seq);
                if (a.size() < 2 || b.size() < 2) break;
                if (a.size() > kLocusSeedSample) a.resize(kLocusSeedSample);
                if (b.size() > kLocusSeedSample) b.resize(kLocusSeedSample);
                hap1 = locus_medoid(a);
                hap2 = locus_medoid(b);
            }
            if (hap2 == nullptr) { note(l0, l1, segs.size(), side1.size(), side2.size(), -1, 0, 0, "no_split"); continue; }
        }
        if (*hap1 == *hap2) { note(l0, l1, segs.size(), side1.size(), side2.size(), 0, 0, 0, "homozygous"); continue; }  // homozygous here
        LocusWindowSite site;
        site.pos = l0 + 1;
        site.end = l1;
        int n0 = 0, n1 = 0;
        // A read is assigned at kLocusMinMargin edits, or at one edit when the
        // two haplotypes differ by a single substitution: they can give no read
        // a two-edit margin, and the window was dropped with every read
        // unassigned. A one-base length difference (a homopolymer call) keeps
        // the two-edit margin: allowing it too gained 622 correct reads for 260
        // discordant on chr12, chr18 and chr20.
        const bool one_substitution =
            hap1->size() == hap2->size() && locus_edit_distance(*hap1, *hap2) == 1;
        const int min_margin = one_substitution ? 1 : kLocusMinMargin;
        // Spanning reads often share a sequence: each is aligned once.
        std::unordered_map<std::string_view, std::pair<int, int>> seen;
        for (const auto& [ri, seq] : segs) {
            auto hit = seen.find(seq);
            if (hit == seen.end())
                hit = seen.emplace(seq, std::make_pair(locus_edit_distance(seq, *hap1), locus_edit_distance(seq, *hap2))).first;
            const int d1 = hit->second.first, d2 = hit->second.second;
            const int margin = std::abs(d1 - d2);
            if (margin < min_margin) {
                // One edit where two are required: kept for the last-resort
                // labels only (weak_observations).
                if (margin == 1) site.weak_observations.emplace_back(ri, d1 < d2 ? 0 : 1);
                continue;
            }
            const int side = d1 < d2 ? 0 : 1;
            (side == 0 ? n0 : n1)++;
            site.observations.emplace_back(
                ri, side, static_cast<float>(std::min(kLocusMaxObsError, 0.5 * std::exp(-margin))));
        }
        const int called = n0 + n1;
        const int dh = locus_edit_distance(*hap1, *hap2);
        if (std::min(n0, n1) < kLocusMinSide || std::min(n0, n1) < kLocusMinMinorFraction * called) {
            note(l0, l1, segs.size(), side1.size(), side2.size(), dh, n0, n1, "one_sided");
            continue;
        }
        note(l0, l1, segs.size(), side1.size(), side2.size(), dh, n0, n1, "window");
        out.push_back(std::move(site));
    }
    if (diag != nullptr) std::fclose(diag);
    return out;
}

std::vector<LocusWindowSite> build_allele_window_sites(const PhasingChunk& bam,
                                                       const GraphChunkBuildResult& graph_chunk,
                                                       const std::vector<LocusWindowSite>& taken) {
    const PhasingChunk& graph = graph_chunk.chunk;
    std::vector<LocusWindowSite> out;
    if (bam.ref_seq.empty()) return out;
    const ReferenceView ref{bam};
    const auto is_acgt = [](const std::string& s) {
        return s.find_first_not_of("ACGTacgt") == std::string::npos;
    };
    // Windows already observed by the EM: haplotype windows, and voting indel
    // rows (their pileup calls are the site's evidence already).
    std::vector<std::pair<hts_pos_t, hts_pos_t>> busy;
    for (const LocusWindowSite& w : taken) busy.emplace_back(w.pos - 1, w.end);
    for (const CandidateVariant& c : graph.candidates) {
        if ((c.lcd_var_i_to_cate & kCandGermlineVarCate) == 0 || c.key.type == VariantType::Snp) continue;
        busy.emplace_back(c.key.pos - 1, c.key.pos - 1 + std::max(1, c.key.ref_len));
    }
    std::sort(busy.begin(), busy.end());
    hts_pos_t widest = 0;
    for (const auto& b : busy) widest = std::max(widest, b.second - b.first);
    const auto overlaps_busy = [&](hts_pos_t a, hts_pos_t b) {
        auto it = std::lower_bound(busy.begin(), busy.end(), std::make_pair(a - widest, hts_pos_t{0}));
        for (; it != busy.end() && it->first < b; ++it)
            if (it->second > a) return true;
        return false;
    };
    // Heterozygous indels the alignment called (noisy-region MSA or repeat
    // context) that the EM does not otherwise observe: injected ones are now
    // voting rows and fall under `busy`.
    // Each candidate is a reference span [p0, p0 + ref_len) and the two allele
    // sequences that replace it on the two haplotypes.
    struct Allele { hts_pos_t p0; int ref_len; std::string seq0; std::string seq1; };
    std::vector<Allele> alleles;
    for (const CandidateVariant& c : bam.candidates) {
        if (c.key.type == VariantType::Snp || !c.msa_insertion_alts.empty()) continue;
        const VariantCategory cat = c.counts.category;
        if (cat != VariantCategory::RepeatHetIndel && cat != VariantCategory::NoisyCandHet) continue;
        if (!is_acgt(c.key.alt) || (c.key.ref_len == 0 && c.key.alt.empty())) continue;
        alleles.push_back(Allele{c.key.pos - 1, c.key.ref_len, ref.slice(c.key.pos - 1, c.key.pos - 1 + c.key.ref_len),
                                 c.key.alt});
    }
    // Two-allele heterozygotes (1/2): at a repeat, the two haplotypes can carry
    // two different indels and no reference allele. Called against REF one at a
    // time, each looks homozygous; together they are one heterozygous site.
    // Pair the two best-supported REF-less indels within one repeat window and
    // contrast the reference carrying one with the reference carrying the other.
    {
        struct Lone { hts_pos_t p0; int ref_len; std::string alt; int reads; hts_pos_t w0; hts_pos_t w1; };
        std::vector<Lone> lone;
        for (const CandidateVariant& c : bam.candidates) {
            if (c.key.type == VariantType::Snp || !c.msa_insertion_alts.empty()) continue;
            if (!is_acgt(c.key.alt) || (c.key.ref_len == 0 && c.key.alt.empty())) continue;
            const int alt_reads = c.counts.alt_cov, ref_reads = c.counts.ref_cov;
            if (alt_reads < kLocusMinSide || ref_reads > kPairMaxRefFraction * (alt_reads + ref_reads)) continue;
            const hts_pos_t p0 = c.key.pos - 1, end0 = p0 + c.key.ref_len;
            lone.push_back(Lone{p0, c.key.ref_len, c.key.alt, alt_reads,
                                p0 - ref.repeat_left(p0), end0 + ref.repeat_right(end0)});
        }
        std::sort(lone.begin(), lone.end(), [](const Lone& a, const Lone& b) { return a.w0 < b.w0; });
        for (size_t i = 0; i < lone.size();) {
            // One group: indels whose repeat spans overlap.
            size_t j = i + 1;
            hts_pos_t group_end = lone[i].w1;
            while (j < lone.size() && lone[j].w0 <= group_end) group_end = std::max(group_end, lone[j++].w1);
            std::vector<const Lone*> group;
            for (size_t k = i; k < j; ++k) group.push_back(&lone[k]);
            i = j;
            if (group.size() < 2) continue;
            std::sort(group.begin(), group.end(), [](const Lone* a, const Lone* b) { return a->reads > b->reads; });
            const Lone& a = *group[0];
            const Lone& b = *group[1];
            if (b.reads < kLocusMinMinorFraction * (a.reads + b.reads)) continue;
            const hts_pos_t lo = std::min(a.p0, b.p0);
            const hts_pos_t hi = std::max(a.p0 + a.ref_len, b.p0 + b.ref_len);
            const auto apply = [&](const Lone& e) {
                return ref.slice(lo, e.p0) + e.alt + ref.slice(e.p0 + e.ref_len, hi);
            };
            std::string seq_a = apply(a), seq_b = apply(b);
            for (char& ch : seq_a) ch = static_cast<char>(std::toupper(static_cast<unsigned char>(ch)));
            for (char& ch : seq_b) ch = static_cast<char>(std::toupper(static_cast<unsigned char>(ch)));
            if (seq_a == seq_b) continue;
            alleles.push_back(Allele{lo, static_cast<int>(hi - lo), std::move(seq_a), std::move(seq_b)});
        }
    }
    // The catalog's non-voting repeat indels, from their VCF-form alleles.
    if (graph_chunk.site_meta.size() == graph.candidates.size() &&
        graph_chunk.site_allele_orig_idx.size() == graph.candidates.size()) {
        for (size_t ci = 0; ci < graph.candidates.size(); ++ci) {
            const CandidateVariant& c = graph.candidates[ci];
            if ((c.lcd_var_i_to_cate & kCandGermlineVarCate) != 0) continue;
            if (c.counts.category != VariantCategory::RepeatHetIndel &&
                c.counts.category != VariantCategory::NoisyCandHet) continue;
            const GraphSiteMeta& meta = graph_chunk.site_meta[ci];
            const std::vector<int>& orig = graph_chunk.site_allele_orig_idx[ci];
            if (meta.pos <= 0 || meta.ref.empty() || meta.non_selected_alt_class || orig.size() < 2) continue;
            const auto allele_seq = [&meta](int k) -> const std::string* {
                if (k == 0) return &meta.ref;
                return k >= 1 && static_cast<size_t>(k) <= meta.alts.size() ? &meta.alts[static_cast<size_t>(k) - 1]
                                                                            : nullptr;
            };
            const std::string* a0 = allele_seq(orig[0]);
            const std::string* a1 = allele_seq(orig[1]);
            if (a0 == nullptr || a1 == nullptr || !is_acgt(*a0) || !is_acgt(*a1) || *a0 == *a1) continue;
            const hts_pos_t p0 = meta.pos - 1;
            const int ref_len = static_cast<int>(meta.ref.size());
            std::string expected = meta.ref;
            for (char& ch : expected) ch = static_cast<char>(std::toupper(static_cast<unsigned char>(ch)));
            if (ref.slice(p0, p0 + ref_len) != expected) continue;
            std::string s0 = *a0, s1 = *a1;
            for (char& ch : s0) ch = static_cast<char>(std::toupper(static_cast<unsigned char>(ch)));
            for (char& ch : s1) ch = static_cast<char>(std::toupper(static_cast<unsigned char>(ch)));
            alleles.push_back(Allele{p0, ref_len, std::move(s0), std::move(s1)});
        }
    }
    std::sort(alleles.begin(), alleles.end(), [](const Allele& a, const Allele& b) { return a.p0 < b.p0; });
    const SpanningReadIndex index(bam, graph);
    hts_pos_t last_end = -1;
    for (const Allele& al : alleles) {
        const hts_pos_t end0 = al.p0 + al.ref_len;
        const hts_pos_t w0 = std::max<hts_pos_t>(0, al.p0 - ref.repeat_left(al.p0) - kAlleleWindowFlank);
        const hts_pos_t w1 = end0 + ref.repeat_right(end0) + kAlleleWindowFlank;
        if (w1 - w0 > kLocusMaxLength || w0 < last_end || overlaps_busy(w0, w1)) continue;
        const std::string left = ref.slice(w0, al.p0), right = ref.slice(end0, w1);
        const std::string hap_ref = left + al.seq0 + right;
        const std::string hap_alt = left + al.seq1 + right;
        if (hap_ref == hap_alt || hap_ref.find('N') != std::string::npos) continue;
        LocusWindowSite site;
        site.pos = w0 + 1;
        site.end = w1;
        int n0 = 0, n1 = 0;
        // The two-edit margin is kept although the alleles are exact: for a
        // one-base indel it assigns no read, so no window is built there, and
        // that is the better outcome measured -- a one-edit margin (HiFi length
        // errors put many reads one edit from the wrong allele) cost 445
        // correct reads and left 453 more unphased on chr12, chr18 and chr20.
        for (const auto& [ri, seq] : index.segments(w0, w1, kAlleleWindowMaxReads)) {
            const int d0 = locus_edit_distance(seq, hap_ref), d1 = locus_edit_distance(seq, hap_alt);
            if (std::abs(d0 - d1) < kLocusMinMargin) continue;
            const int side = d0 < d1 ? 0 : 1;
            (side == 0 ? n0 : n1)++;
            site.observations.emplace_back(
                ri, side, static_cast<float>(std::min(kLocusMaxObsError, 0.5 * std::exp(-std::abs(d0 - d1)))));
        }
        const int called = n0 + n1;
        if (std::min(n0, n1) < kLocusMinSide || std::min(n0, n1) < kLocusMinMinorFraction * called) continue;
        last_end = w1;
        out.push_back(std::move(site));
    }
    return out;
}

// Length windows (see build_length_window_sites).
namespace {
// Run-length buckets of the stutter model (bp of the reference tract) and the
// largest length change modelled; beyond it the floor probability applies.
constexpr int kLengthBuckets[] = {6, 10, 15, 20, 30};
constexpr int kLengthMaxDev = 8;
constexpr int kLengthMinTract = 6;          // homopolymer / STR tract, bp
constexpr int kLengthMinReads = 8;          // spanning reads to genotype a locus
constexpr double kLengthMinHetLLR = 6.0;    // two lengths over one, log units
constexpr double kLengthMinReadLLR = 0.5;   // a read nearer one length than this calls nothing

int length_bucket(hts_pos_t tract) {
    int b = 0;
    for (int i = 0; i < static_cast<int>(sizeof(kLengthBuckets) / sizeof(int)); ++i)
        if (tract >= kLengthBuckets[i]) b = i;
    return b;
}
}  // namespace

std::vector<LocusWindowSite> build_length_window_sites(const PhasingChunk& bam,
                                                       const GraphChunkBuildResult& graph_chunk,
                                                       const std::vector<LocusWindowSite>& taken) {
    const PhasingChunk& graph = graph_chunk.chunk;
    std::vector<LocusWindowSite> out;
    if (bam.ref_seq.empty()) return out;
    const ReferenceView ref{bam};
    // Loci the EM already observes: windows, and voting indel rows.
    std::vector<std::pair<hts_pos_t, hts_pos_t>> busy;
    for (const LocusWindowSite& w : taken) busy.emplace_back(w.pos - 1, w.end);
    for (const CandidateVariant& c : graph.candidates) {
        if ((c.lcd_var_i_to_cate & kCandGermlineVarCate) == 0 || c.key.type == VariantType::Snp) continue;
        busy.emplace_back(c.key.pos - 1, c.key.pos - 1 + std::max(1, c.key.ref_len));
    }
    std::sort(busy.begin(), busy.end());
    hts_pos_t widest = 0;
    for (const auto& b : busy) widest = std::max(widest, b.second - b.first);
    const auto overlaps_busy = [&](hts_pos_t a, hts_pos_t b) {
        auto it = std::lower_bound(busy.begin(), busy.end(), std::make_pair(a - widest, hts_pos_t{0}));
        for (; it != busy.end() && it->first < b; ++it)
            if (it->second > a) return true;
        return false;
    };
    // Tracts: the repeat around every indel candidate of the alignment solve,
    // merged where they overlap.
    struct Tract { hts_pos_t t0; hts_pos_t t1; int period; };
    std::vector<Tract> tracts;
    for (const CandidateVariant& c : bam.candidates) {
        if (c.key.type == VariantType::Snp) continue;
        const hts_pos_t p0 = c.key.pos - 1;
        const hts_pos_t end0 = p0 + (c.key.type == VariantType::Deletion ? c.key.ref_len : 0);
        const hts_pos_t t0 = p0 - ref.repeat_left(p0), t1 = end0 + ref.repeat_right(end0);
        if (t1 - t0 < kLengthMinTract) continue;
        const std::string edit = c.key.type == VariantType::Deletion ? ref.slice(p0, end0) : c.key.alt;
        int period = static_cast<int>(edit.size());
        for (int q = 1; q <= static_cast<int>(edit.size()); ++q) {
            if (edit.size() % static_cast<size_t>(q) != 0) continue;
            bool rep = true;
            for (size_t i = static_cast<size_t>(q); i < edit.size() && rep; ++i) rep = edit[i] == edit[i - static_cast<size_t>(q)];
            if (rep) { period = q; break; }
        }
        tracts.push_back(Tract{t0, t1, std::max(1, period)});
    }
    std::sort(tracts.begin(), tracts.end(), [](const Tract& a, const Tract& b) { return a.t0 < b.t0; });
    std::vector<Tract> merged;
    for (const Tract& t : tracts) {
        if (!merged.empty() && t.t0 <= merged.back().t1) {
            merged.back().t1 = std::max(merged.back().t1, t.t1);
            merged.back().period = std::min(merged.back().period, t.period);
        } else {
            merged.push_back(t);
        }
    }
    // Each tract's reads: their length change across the tract.
    const SpanningReadIndex index(bam, graph);
    struct Locus { hts_pos_t w0; hts_pos_t w1; int period; hts_pos_t tract; std::vector<std::pair<size_t, int>> reads; };
    std::vector<Locus> loci;
    for (const Tract& t : merged) {
        const hts_pos_t w0 = std::max<hts_pos_t>(0, t.t0 - kAlleleWindowFlank), w1 = t.t1 + kAlleleWindowFlank;
        if (w1 - w0 > kLocusMaxLength || overlaps_busy(w0, w1)) continue;
        Locus l{w0, w1, t.period, t.t1 - t.t0, {}};
        for (const auto& [ri, seq] : index.segments(w0, w1, kAlleleWindowMaxReads))
            l.reads.emplace_back(ri, static_cast<int>(seq.size()) - static_cast<int>(w1 - w0));
        if (static_cast<int>(l.reads.size()) >= kLengthMinReads) loci.push_back(std::move(l));
    }
    if (loci.empty()) return out;
    // Stutter: deviation from each locus's most common length, by tract bucket.
    constexpr int nb = static_cast<int>(sizeof(kLengthBuckets) / sizeof(int));
    std::vector<std::vector<double>> dev_counts(nb, std::vector<double>(2 * kLengthMaxDev + 1, 1.0));
    for (const Locus& l : loci) {
        std::map<int, int> hist;
        for (const auto& [ri, d] : l.reads) ++hist[d];
        int mode = hist.begin()->first, best = 0;
        for (const auto& [d, n] : hist) if (n > best) { best = n; mode = d; }
        for (const auto& [ri, d] : l.reads) {
            const int dev = d - mode;
            if (std::abs(dev) <= kLengthMaxDev) dev_counts[length_bucket(l.tract)][dev + kLengthMaxDev] += 1.0;
        }
    }
    std::vector<std::vector<double>> logp(nb, std::vector<double>(2 * kLengthMaxDev + 1));
    double floor = 0.0;
    for (int b = 0; b < nb; ++b) {
        double tot = 0.0;
        for (double v : dev_counts[b]) tot += v;
        for (int i = 0; i <= 2 * kLengthMaxDev; ++i) {
            logp[b][i] = std::log(dev_counts[b][i] / tot);
            floor = std::min(floor, logp[b][i]);
        }
    }
    floor -= 2.0;
    const auto lp = [&](int b, int d, int len) {
        const int dev = d - len;
        return std::abs(dev) <= kLengthMaxDev ? logp[b][dev + kLengthMaxDev] : floor;
    };
    // Diagnostic, off by default: PGPHASE_DIAG=DIR writes every length window
    // (lengthwin.TID.BEG.tsv): 1-based window start and end, repeat period,
    // tract bp, the two lengths (change from the reference), het score,
    // reads on each side, one repeat unit apart (1/0), then read:side per
    // called read.
    static const char* const diag_dir = std::getenv("PGPHASE_DIAG");
    FILE* dump = diag_dir == nullptr ? nullptr
        : std::fopen((std::string(diag_dir) + "/lengthwin." + std::to_string(graph.region.tid) + "." +
                      std::to_string(graph.region.beg) + ".tsv").c_str(), "w");
    for (const Locus& l : loci) {
        const int b = length_bucket(l.tract);
        int lo = l.reads.front().second, hi = lo;
        for (const auto& [ri, d] : l.reads) { lo = std::min(lo, d); hi = std::max(hi, d); }
        double ll1 = -1e300;
        for (int c = lo - 1; c <= hi + 1; ++c) {
            double s = 0.0;
            for (const auto& [ri, d] : l.reads) s += lp(b, d, c);
            ll1 = std::max(ll1, s);
        }
        double ll2 = -1e300;
        int c1 = 0, c2 = 0;
        for (int x = lo - 1; x <= hi + 1; ++x) {
            for (int y = x + 1; y <= hi + 1; ++y) {
                double s = 0.0;
                for (const auto& [ri, d] : l.reads) {
                    const double a = lp(b, d, x), bb = lp(b, d, y);
                    const double m = std::max(a, bb);
                    s += m + std::log(0.5 * std::exp(a - m) + 0.5 * std::exp(bb - m));
                }
                if (s > ll2) { ll2 = s; c1 = x; c2 = y; }
            }
        }
        if (ll2 - ll1 < kLengthMinHetLLR) continue;
        LocusWindowSite site;
        site.pos = l.w0 + 1;
        site.end = l.w1;
        int n0 = 0, n1 = 0;
        for (const auto& [ri, d] : l.reads) {
            const double r = lp(b, d, c1) - lp(b, d, c2);
            if (std::abs(r) < kLengthMinReadLLR) continue;
            const int side = r > 0 ? 0 : 1;
            (side == 0 ? n0 : n1)++;
            site.observations.emplace_back(
                ri, side, static_cast<float>(std::min(kLocusMaxObsError, 1.0 / (1.0 + std::exp(std::abs(r))))));
        }
        if (std::min(n0, n1) < kLocusMinSide || std::min(n0, n1) < kLocusMinMinorFraction * (n0 + n1)) continue;
        if (dump != nullptr) {
            std::fprintf(dump, "%lld\t%lld\t%d\t%lld\t%d\t%d\t%.2f\t%d\t%d\t%d", static_cast<long long>(site.pos),
                         static_cast<long long>(site.end), l.period, static_cast<long long>(l.tract), c1, c2,
                         ll2 - ll1, n0, n1, std::abs(c2 - c1) <= l.period ? 1 : 0);
            for (const auto& [ri, side, e] : site.observations)
                std::fprintf(dump, "\t%s:%d", graph.reads[ri].qname.c_str(), side);
            std::fputc('\n', dump);
        }
        out.push_back(std::move(site));
    }
    if (dump != nullptr) std::fclose(dump);
    return out;
}

} // namespace pgphase_collect
