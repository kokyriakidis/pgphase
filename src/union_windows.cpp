// Union gap phasing: haplotype (medoid) windows and allele windows over noisy
// loci, as extra EM sites.

#include "union_internal.hpp"

#include "edlib.h"

namespace pgphase_collect {

// Query index aligned at (or, inside a deletion, after) a 0-based reference position.
hts_pos_t locus_query_index(const bam1_t* b, hts_pos_t ref_pos) {
    if (ref_pos < b->core.pos || ref_pos >= bam_endpos(b)) return -1;
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

static const std::string* locus_medoid(const std::vector<const std::string*>& seqs) {
    const std::string* best = nullptr;
    long best_cost = -1;
    for (size_t i = 0; i < seqs.size(); ++i) {
        long cost = 0;
        for (size_t j = 0; j < seqs.size(); ++j)
            if (j != i) cost += locus_edit_distance(*seqs[i], *seqs[j]);
        if (best_cost < 0 || cost < best_cost) { best = seqs[i]; best_cost = cost; }
    }
    return best;
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
    for (const auto& [a, b] : merged) {
        const hts_pos_t l0 = std::max<hts_pos_t>(0, a - kLocusPad);
        const hts_pos_t l1 = b + ref.repeat_right(b) + kLocusPad;
        if (l1 - l0 > kLocusMaxLength) continue;
        const std::vector<std::pair<size_t, std::string>> segs = index.segments(l0, l1, kLocusMaxReads);
        if (segs.size() < static_cast<size_t>(2 * kLocusMinSide)) continue;
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
            if (hap2 == nullptr) continue;
        }
        if (*hap1 == *hap2) continue;  // homozygous here
        LocusWindowSite site;
        site.pos = l0 + 1;
        site.end = l1;
        int n0 = 0, n1 = 0;
        for (const auto& [ri, seq] : segs) {
            const int d1 = locus_edit_distance(seq, *hap1), d2 = locus_edit_distance(seq, *hap2);
            const int margin = std::abs(d1 - d2);
            if (margin < kLocusMinMargin) continue;
            const int side = d1 < d2 ? 0 : 1;
            (side == 0 ? n0 : n1)++;
            site.observations.emplace_back(
                ri, side, static_cast<float>(std::min(kLocusMaxObsError, 0.5 * std::exp(-margin))));
        }
        const int called = n0 + n1;
        if (std::min(n0, n1) < kLocusMinSide || std::min(n0, n1) < kLocusMinMinorFraction * called) continue;
        out.push_back(std::move(site));
    }
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

} // namespace pgphase_collect
