// Union gap phasing, last resort: labels from MSA-verified indels and
// haplotype markers, applied through anchor reads after stitching.

#include "union_internal.hpp"

namespace pgphase_collect {

// One alignment event of a read against the reference: a mismatching base, an
// insertion or a deletion, with indels left-normalized so placements of one
// homopolymer or repeat length change compare equal.
struct ReadEvent {
    hts_pos_t pos;    // 0-based: the base (SNP), the first deleted base, or the base an insertion precedes
    std::string what; // "S" + base, "I" + inserted sequence, "D" + deleted length
    hts_pos_t end;    // 0-based exclusive reference end the event touches
    bool operator<(const ReadEvent& o) const { return pos != o.pos ? pos < o.pos : what < o.what; }
    bool operator==(const ReadEvent& o) const { return pos == o.pos && what == o.what; }
};

static std::vector<ReadEvent> read_events(const bam1_t* b, const ReferenceView& ref, hts_pos_t lo, hts_pos_t hi) {
    std::vector<ReadEvent> out;
    const uint32_t* cigar = bam_get_cigar(b);
    const uint8_t* seq = bam_get_seq(b);
    const uint8_t* qual = bam_get_qual(b);
    hts_pos_t r = b->core.pos, q = 0;
    for (uint32_t i = 0; i < b->core.n_cigar && r < hi; ++i) {
        const int op = bam_cigar_op(cigar[i]);
        const hts_pos_t len = bam_cigar_oplen(cigar[i]);
        if (op == BAM_CMATCH || op == BAM_CEQUAL || op == BAM_CDIFF) {
            for (hts_pos_t k = 0; k < len; ++k) {
                const hts_pos_t p = r + k;
                if (p < lo || p >= hi) continue;
                const char base = seq_nt16_str[bam_seqi(seq, q + k)];
                if (base == 'N' || qual[q + k] < kEventMinBaseQuality) continue;
                const char rb = ref.base(p);
                if (rb != 'N' && base != rb) out.push_back(ReadEvent{p, std::string("S") + base, p + 1});
            }
            r += len;
            q += len;
        } else if (op == BAM_CINS) {
            if (r >= lo && r < hi) {
                std::string ins(static_cast<size_t>(len), 'N');
                for (hts_pos_t k = 0; k < len; ++k) ins[static_cast<size_t>(k)] = seq_nt16_str[bam_seqi(seq, q + k)];
                hts_pos_t p = r;
                while (p > 0 && ref.base(p - 1) != 'N' && ref.base(p - 1) == ins.back()) {
                    ins = std::string(1, ref.base(p - 1)) + ins.substr(0, ins.size() - 1);
                    --p;
                }
                out.push_back(ReadEvent{p, "I" + ins, p});
            }
            q += len;
        } else if (op == BAM_CDEL) {
            if (r >= lo && r < hi) {
                hts_pos_t p = r;
                while (p > 0 && ref.base(p - 1) != 'N' && ref.base(p - 1) == ref.base(p + len - 1)) --p;
                out.push_back(ReadEvent{p, "D" + std::to_string(len), p + len});
            }
            r += len;
        } else if (op == BAM_CREF_SKIP) {
            r += len;
        } else if (op == BAM_CSOFT_CLIP) {
            q += len;
        }
    }
    std::sort(out.begin(), out.end());
    out.erase(std::unique(out.begin(), out.end()), out.end());
    return out;
}

size_t label_reads_from_haplotype_consensus(GraphChunkBuildResult& graph_chunk, const PhasingChunk& bam) {
    PhasingChunk& chunk = graph_chunk.chunk;
    if (bam.ref_seq.empty()) return 0;
    const ReferenceView ref{bam};
    const SpanningReadIndex index(bam, chunk);
    struct Span { size_t i; hts_pos_t beg; hts_pos_t end; };  // index entry, aligned span
    std::unordered_set<size_t> decided_before;
    for (const auto& d : graph_chunk.deferred_read_labels) decided_before.insert(std::get<0>(d));
    std::vector<Span> unlabelled, labelled;
    for (size_t i = 0; i < index.reads.size(); ++i) {
        const size_t ri = index.graph_index[i];
        if (ri >= chunk.haps.size() || chunk.reads[ri].is_skipped) continue;
        if (chunk.haps[ri] == 0 && decided_before.count(ri) != 0) continue;
        const bam1_t* aln = index.reads[i]->alignment.get();
        const Span sp{i, aln->core.pos, bam_endpos(aln)};
        if (chunk.haps[ri] == 0) unlabelled.push_back(sp);
        else if (chunk.phase_sets[ri] > 0) labelled.push_back(sp);
    }
    if (unlabelled.empty() || labelled.empty()) return 0;
    // Only the stretches unlabelled reads cover matter.
    std::vector<std::pair<hts_pos_t, hts_pos_t>> regions;
    for (const Span& u : unlabelled) regions.emplace_back(u.beg, u.end);
    std::sort(regions.begin(), regions.end());
    std::vector<std::pair<hts_pos_t, hts_pos_t>> merged;
    for (const auto& r : regions) {
        if (!merged.empty() && r.first <= merged.back().second) merged.back().second = std::max(merged.back().second, r.second);
        else merged.push_back(r);
    }
    // Events of the labelled reads inside those stretches, by (phase set, haplotype).
    std::map<std::pair<hts_pos_t, ReadEvent>, std::array<int, 2>> carriers;  // (ps, event) -> reads per haplotype
    std::map<hts_pos_t, std::vector<std::pair<Span, int>>> by_set;           // ps -> (span, hap index)
    for (const Span& l : labelled) {
        const size_t ri = index.graph_index[l.i];
        const hts_pos_t ps = chunk.phase_sets[ri];
        const int h = chunk.haps[ri] == 1 ? 0 : 1;
        auto it = std::upper_bound(merged.begin(), merged.end(), std::make_pair(l.end, hts_pos_t{0}));
        bool touches = false;
        for (auto m = merged.begin(); m != it; ++m)
            if (m->second > l.beg) { touches = true; break; }
        if (!touches) continue;
        by_set[ps].emplace_back(l, h);
        for (const auto& [a, b] : merged) {
            if (b <= l.beg || a >= l.end) continue;
            for (const ReadEvent& e : read_events(index.reads[l.i]->alignment.get(), ref, std::max(a, l.beg), std::min(b, l.end)))
                ++carriers[{ps, e}][h];
        }
    }
    // Discriminating events: carried by one haplotype's covering reads and
    // absent from the other's.
    struct Marker { hts_pos_t ps; ReadEvent e; int carrier; };
    std::vector<Marker> markers;
    for (const auto& [key, cnt] : carriers) {
        if (std::max(cnt[0], cnt[1]) < kConsensusMinReads) continue;
        const auto& [ps, e] = key;
        std::array<int, 2> cover{0, 0};
        for (const auto& [sp, h] : by_set[ps])
            if (sp.beg <= e.pos - kConsensusFlank && sp.end >= e.end + kConsensusFlank) ++cover[h];
        if (cover[0] < kConsensusMinReads || cover[1] < kConsensusMinReads) continue;
        const double f0 = static_cast<double>(cnt[0]) / cover[0], f1 = static_cast<double>(cnt[1]) / cover[1];
        if (f0 >= kConsensusCarrierFraction && f1 <= 1.0 - kConsensusCarrierFraction) markers.push_back(Marker{ps, e, 0});
        else if (f1 >= kConsensusCarrierFraction && f0 <= 1.0 - kConsensusCarrierFraction) markers.push_back(Marker{ps, e, 1});
    }
    if (markers.empty()) return 0;
    std::sort(markers.begin(), markers.end(), [](const Marker& a, const Marker& b) { return a.e.pos < b.e.pos; });
    std::map<hts_pos_t, size_t> anchor_of;
    for (size_t ri = 0; ri < chunk.haps.size(); ++ri)
        if (chunk.haps[ri] != 0 && chunk.phase_sets[ri] > 0) anchor_of.try_emplace(chunk.phase_sets[ri], ri);
    size_t decided = 0;
    for (const Span& u : unlabelled) {
        const auto first = std::lower_bound(markers.begin(), markers.end(), u.beg + kConsensusFlank,
            [](const Marker& m, hts_pos_t p) { return m.e.pos < p; });
        std::map<hts_pos_t, std::array<int, 2>> votes;  // ps -> votes for haplotype index
        std::vector<ReadEvent> own;
        bool have_own = false;
        for (auto m = first; m != markers.end() && m->e.pos < u.end - kConsensusFlank; ++m) {
            if (m->e.end + kConsensusFlank > u.end) continue;
            if (!have_own) {
                own = read_events(index.reads[u.i]->alignment.get(), ref, u.beg, u.end);
                have_own = true;
            }
            const bool carries = std::binary_search(own.begin(), own.end(), m->e);
            ++votes[m->ps][carries ? m->carrier : 1 - m->carrier];
        }
        if (votes.empty()) continue;
        const auto best = std::max_element(votes.begin(), votes.end(), [](const auto& a, const auto& b) {
            return a.second[0] + a.second[1] < b.second[0] + b.second[1];
        });
        const int n = best->second[0] + best->second[1];
        const int lead = std::max(best->second[0], best->second[1]);
        if (lead < kConsensusMinAgreement * n) continue;
        const int hap = best->second[0] >= best->second[1] ? 1 : 2;
        const auto anchor = anchor_of.find(best->first);
        if (anchor == anchor_of.end()) continue;
        const size_t ri = index.graph_index[u.i];
        graph_chunk.deferred_read_labels.emplace_back(ri, anchor->second, chunk.haps[anchor->second] == hap);
        ++decided;
    }
    return decided;
}

size_t label_reads_from_verified_indels(GraphChunkBuildResult& graph_chunk, const PhasingChunk& bam) {
    PhasingChunk& chunk = graph_chunk.chunk;
    if (bam.ref_seq.empty()) return 0;
    const ReferenceView ref{bam};
    const SpanningReadIndex index(bam, chunk);
    const auto is_acgt = [](const std::string& s) {
        return s.find_first_not_of("ACGTacgt") == std::string::npos;
    };
    // (phase set, haplotype) votes per unlabelled read; a read with conflicting
    // votes stays unlabelled.
    std::map<size_t, std::map<std::pair<hts_pos_t, int>, int>> votes;
    std::vector<LocusAllele> locus_alleles;
    for (const CandidateVariant& c : bam.candidates) {
        if (c.key.type == VariantType::Snp || !c.msa_insertion_alts.empty()) continue;
        if (!is_acgt(c.key.alt) || (c.key.ref_len == 0 && c.key.alt.empty())) continue;
        if (c.counts.alt_cov < kLocusMinSide) continue;
        locus_alleles.push_back(LocusAllele{c.key.pos - 1, c.key.pos - 1 + c.key.ref_len, c.key.alt});
    }
    std::sort(locus_alleles.begin(), locus_alleles.end(), [](const LocusAllele& a, const LocusAllele& b) { return a.p0 < b.p0; });
    hts_pos_t last_end = -1;
    for (const CandidateVariant& c : bam.candidates) {
        if (c.counts.category != VariantCategory::NoisyCandHet || !c.msa_verified) continue;
        if (c.key.type == VariantType::Snp || !c.msa_insertion_alts.empty()) continue;
        if (!is_acgt(c.key.alt) || (c.key.ref_len == 0 && c.key.alt.empty())) continue;
        const hts_pos_t p0 = c.key.pos - 1;
        const hts_pos_t end0 = p0 + c.key.ref_len;
        const hts_pos_t w0 = std::max<hts_pos_t>(0, p0 - ref.repeat_left(p0) - kAlleleWindowFlank);
        const hts_pos_t w1 = end0 + ref.repeat_right(end0) + kAlleleWindowFlank;
        if (w1 - w0 > kLocusMaxLength || w0 < last_end) continue;  // one call per repeat
        const std::string hap_ref = ref.slice(w0, w1);
        const std::string hap_alt = ref.slice(w0, p0) + c.key.alt + ref.slice(end0, w1);
        if (hap_ref == hap_alt || hap_ref.find('N') != std::string::npos) continue;
        std::vector<std::pair<size_t, int>> calls;  // read, allele
        const std::vector<std::string> competitors =
            competing_haplotypes(ref, w0, w1, LocusAllele{p0, end0, c.key.alt}, locus_alleles);
        for (const auto& [ri, seq] : index.segments(w0, w1, kAlleleWindowMaxReads)) {
            const int a = closest_allele_call(seq, hap_alt, competitors, 1);
            if (a >= 0) calls.emplace_back(ri, a);
        }
        // The site's orientation from the reads already labelled: which allele
        // haplotype 1 carries, in the phase set most of them share.
        std::map<hts_pos_t, std::array<int, 2>> by_set;  // (hap1 carries ALT, hap1 carries REF)
        for (const auto& [ri, a] : calls) {
            if (ri >= chunk.haps.size() || chunk.haps[ri] == 0 || chunk.phase_sets[ri] <= 0) continue;
            ++by_set[chunk.phase_sets[ri]][(chunk.haps[ri] == 1) == (a == 1) ? 0 : 1];
        }
        if (by_set.empty()) continue;
        const auto best = std::max_element(by_set.begin(), by_set.end(), [](const auto& x, const auto& y) {
            return x.second[0] + x.second[1] < y.second[0] + y.second[1];
        });
        const int n = best->second[0] + best->second[1];
        const int agree = std::max(best->second[0], best->second[1]);
        if (n < kLastResortMinLabelled || agree < kLastResortMinAgreement * n) continue;
        last_end = w1;
        const bool hap1_alt = best->second[0] >= best->second[1];
        for (const auto& [ri, a] : calls) {
            if (ri >= chunk.haps.size() || chunk.haps[ri] != 0) continue;
            const int hap = (a == 1) == hap1_alt ? 1 : 2;
            ++votes[ri][{best->first, hap}];
        }
    }
    // Decisions are tied to an anchor read of the same phase set and applied
    // after chunk stitching, so a last-resort label never votes in it and
    // follows any flip the stitch applies to its block.
    std::map<hts_pos_t, size_t> anchor_of;
    for (size_t ri = 0; ri < chunk.haps.size(); ++ri)
        if (chunk.haps[ri] != 0 && chunk.phase_sets[ri] > 0) anchor_of.try_emplace(chunk.phase_sets[ri], ri);
    size_t labelled = 0;
    for (const auto& [ri, v] : votes) {
        if (v.size() != 1) continue;  // the read's verified indels disagree
        const auto& [ps, hap] = v.begin()->first;
        const auto anchor = anchor_of.find(ps);
        if (anchor == anchor_of.end()) continue;
        graph_chunk.deferred_read_labels.emplace_back(ri, anchor->second, chunk.haps[anchor->second] == hap);
        ++labelled;
    }
    return labelled;
}

size_t apply_deferred_read_labels(GraphChunkBuildResult& graph_chunk) {
    PhasingChunk& chunk = graph_chunk.chunk;
    size_t applied = 0;
    for (const auto& [ri, anchor, same] : graph_chunk.deferred_read_labels) {
        if (ri >= chunk.haps.size() || anchor >= chunk.haps.size() || chunk.haps[ri] != 0) continue;
        const int anchor_hap = chunk.haps[anchor];
        if ((anchor_hap != 1 && anchor_hap != 2) || chunk.phase_sets[anchor] <= 0) continue;
        chunk.haps[ri] = same ? anchor_hap : 3 - anchor_hap;
        chunk.phase_sets[ri] = chunk.phase_sets[anchor];
        ++applied;
    }
    graph_chunk.deferred_read_labels.clear();
    return applied;
}

} // namespace pgphase_collect
