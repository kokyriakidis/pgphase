// Union gap phasing: alignment-only reads, their fills, and realignment of
// uncalled reads at injected indel sites.

#include "union_internal.hpp"

namespace pgphase_collect {

std::vector<std::string> competing_haplotypes(const ReferenceView& ref, hts_pos_t w0, hts_pos_t w1,
                                                     const LocusAllele& self,
                                                     const std::vector<LocusAllele>& sorted_alleles) {
    std::vector<std::string> haps{ref.slice(w0, w1)};
    // Only alleles of the same locus: edits inside this allele's own repeat.
    // A different variant elsewhere in the window would otherwise stand in
    // for the reference, and reads carrying both would be miscalled.
    const hts_pos_t lo = self.p0 - ref.repeat_left(self.p0);
    const hts_pos_t hi = self.end0 + ref.repeat_right(self.end0);
    auto it = std::lower_bound(sorted_alleles.begin(), sorted_alleles.end(), lo,
                               [](const LocusAllele& a, hts_pos_t p) { return a.p0 < p; });
    for (; it != sorted_alleles.end() && it->p0 <= hi; ++it) {
        if (it->end0 > hi || it->p0 < w0 || it->end0 > w1) continue;
        if (it->p0 == self.p0 && it->end0 == self.end0 && it->alt == self.alt) continue;
        std::string h = ref.slice(w0, it->p0) + it->alt + ref.slice(it->end0, w1);
        if (std::find(haps.begin(), haps.end(), h) == haps.end()) haps.push_back(std::move(h));
    }
    return haps;
}

int closest_allele_call(const std::string& seq, const std::string& hap_alt,
                               const std::vector<std::string>& competitors, int min_margin) {
    const int d1 = locus_edit_distance(seq, hap_alt);
    int d0 = -1;
    for (const std::string& h : competitors) {
        if (h == hap_alt) continue;
        const int d = locus_edit_distance(seq, h);
        if (d0 < 0 || d < d0) d0 = d;
    }
    if (d0 < 0 || std::abs(d0 - d1) < min_margin) return -1;
    return d1 < d0 ? 1 : 0;
}

// Index of the closest of several haplotype sequences, or -1 when the best
// does not beat the runner-up by min_margin edits.
static int closest_haplotype(const std::string& seq, const std::vector<std::string>& haps, int min_margin) {
    int best = -1, best_d = -1, second_d = -1;
    for (size_t i = 0; i < haps.size(); ++i) {
        const int d = locus_edit_distance(seq, haps[i]);
        if (best_d < 0 || d < best_d) {
            second_d = best_d;
            best_d = d;
            best = static_cast<int>(i);
        } else if (second_d < 0 || d < second_d) {
            second_d = d;
        }
    }
    if (second_d >= 0 && second_d - best_d < min_margin) return -1;
    return best;
}

// Haplotype sequences of a multi-allelic MSA insertion row over [w0, w1):
// index 0 the reference, index i the reference with msa_insertion_alts[i-1].
static std::vector<std::string> msa_allele_haplotypes(const ReferenceView& ref, const CandidateVariant& c,
                                                      hts_pos_t w0, hts_pos_t w1) {
    const hts_pos_t p0 = c.key.pos - 1;
    std::vector<std::string> haps{ref.slice(w0, w1)};
    for (const std::string& ins : c.msa_insertion_alts) {
        if (ins.find_first_not_of("ACGTacgt") != std::string::npos) return {};
        haps.push_back(ref.slice(w0, p0) + ins + ref.slice(p0, w1));
    }
    return haps;
}

// Widen a read profile's candidate range to [lo, hi], keeping every channel aligned.
static void widen_profile(ReadVariantProfile& p, int lo, int hi) {
    if (lo == p.start_var_idx && hi == p.end_var_idx) return;
    const size_t span = static_cast<size_t>(hi - lo + 1);
    if (p.start_var_idx < 0) {  // no calls yet: every channel starts empty
        p.alleles.clear();
        p.alt_qi.clear();
        p.graph_alleles.clear();
        p.bam_alleles.clear();
        p.bam_qi.clear();
        p.bam_base_qualities.clear();
    }
    const size_t shift = p.start_var_idx >= 0 ? static_cast<size_t>(p.start_var_idx - lo) : 0;
    const auto widen = [&](auto& v, auto fill, bool required) {
        if (v.empty() && !required) return;
        std::remove_reference_t<decltype(v)> w(span, fill);
        for (size_t k = 0; k < v.size(); ++k) w[shift + k] = v[k];
        v = std::move(w);
    };
    widen(p.alleles, -1, true);
    widen(p.alt_qi, 0, true);
    widen(p.graph_alleles, -1, false);
    widen(p.bam_alleles, -1, false);
    widen(p.bam_qi, 0, false);
    widen(p.bam_base_qualities, static_cast<uint8_t>(0), false);
    p.start_var_idx = lo;
    p.end_var_idx = hi;
}

size_t realign_indel_observations(PhasingChunk& chunk, const PhasingChunk& bam) {
    if (bam.ref_seq.empty()) return 0;
    const ReferenceView ref{bam};
    const SpanningReadIndex index(bam, chunk);
    const auto is_acgt = [](const std::string& s) {
        return s.find_first_not_of("ACGTacgt") == std::string::npos;
    };
    // Read -> its profile, for writing calls at a candidate index.
    std::vector<ReadVariantProfile*> profile_of(chunk.reads.size(), nullptr);
    for (ReadVariantProfile& p : chunk.read_var_profile)
        if (p.read_id >= 0 && static_cast<size_t>(p.read_id) < profile_of.size()) profile_of[p.read_id] = &p;
    std::map<size_t, std::vector<std::pair<size_t, int>>> calls;  // read -> (candidate, allele or -1)
    std::vector<LocusAllele> locus_alleles;
    for (const CandidateVariant& c : chunk.candidates) {
        if (!c.bam_injected || c.key.type == VariantType::Snp || !c.msa_insertion_alts.empty()) continue;
        if (!is_acgt(c.key.alt) || (c.key.ref_len == 0 && c.key.alt.empty())) continue;
        locus_alleles.push_back(LocusAllele{c.key.pos - 1, c.key.pos - 1 + c.key.ref_len, c.key.alt});
    }
    std::sort(locus_alleles.begin(), locus_alleles.end(), [](const LocusAllele& a, const LocusAllele& b) { return a.p0 < b.p0; });
    for (size_t ci = 0; ci < chunk.candidates.size(); ++ci) {
        const CandidateVariant& c = chunk.candidates[ci];
        if (!c.bam_injected || c.key.type == VariantType::Snp) continue;
        if (!c.msa_insertion_alts.empty()) {
            // One row, several insertion alleles: the closest allele by index.
            if (c.key.type != VariantType::Insertion || (c.lcd_var_i_to_cate & kCandGermlineVarCate) == 0) continue;
            const hts_pos_t p0 = c.key.pos - 1;
            const hts_pos_t w0 = std::max<hts_pos_t>(0, p0 - ref.repeat_left(p0) - kAlleleWindowFlank);
            const hts_pos_t w1 = p0 + ref.repeat_right(p0) + kAlleleWindowFlank;
            if (w1 - w0 > kLocusMaxLength) continue;
            const std::vector<std::string> haps = msa_allele_haplotypes(ref, c, w0, w1);
            if (haps.size() < 3) continue;
            for (const auto& [ri, seq] : index.segments(w0, w1, kAlleleWindowMaxReads))
                if (profile_of[ri] != nullptr) calls[ri].emplace_back(ci, closest_haplotype(seq, haps, 1));
            continue;
        }
        if ((c.lcd_var_i_to_cate & kCandGermlineVarCate) == 0) continue;
        if (!is_acgt(c.key.alt) || (c.key.ref_len == 0 && c.key.alt.empty())) continue;
        const hts_pos_t p0 = c.key.pos - 1;
        const hts_pos_t end0 = p0 + c.key.ref_len;
        const hts_pos_t w0 = std::max<hts_pos_t>(0, p0 - ref.repeat_left(p0) - kAlleleWindowFlank);
        const hts_pos_t w1 = end0 + ref.repeat_right(end0) + kAlleleWindowFlank;
        if (w1 - w0 > kLocusMaxLength) continue;
        const std::string hap_ref = ref.slice(w0, w1);
        const std::string hap_alt = ref.slice(w0, p0) + c.key.alt + ref.slice(end0, w1);
        if (hap_ref == hap_alt || hap_ref.find('N') != std::string::npos) continue;
        // Allele 0 is the closest of REF and every other allele of the locus.
        const std::vector<std::string> competitors =
            competing_haplotypes(ref, w0, w1, LocusAllele{p0, end0, c.key.alt}, locus_alleles);
        for (const auto& [ri, seq] : index.segments(w0, w1, kAlleleWindowMaxReads)) {
            if (profile_of[ri] == nullptr) continue;
            calls[ri].emplace_back(ci, closest_allele_call(seq, hap_alt, competitors, 1));
        }
    }
    // Every read spanning the site's window is re-called there, its range
    // widened if its other calls stop short of the site.
    size_t changed = 0;
    for (const auto& [ri, list] : calls) {
        ReadVariantProfile* prof = profile_of[ri];
        int lo = prof->start_var_idx, hi = prof->end_var_idx;
        for (const auto& [ci, a] : list) {
            if (a < 0 && (lo < 0 || static_cast<int>(ci) < lo || static_cast<int>(ci) > hi)) continue;
            lo = lo < 0 ? static_cast<int>(ci) : std::min(lo, static_cast<int>(ci));
            hi = std::max(hi, static_cast<int>(ci));
        }
        if (lo < 0) continue;
        widen_profile(*prof, lo, hi);
        for (const auto& [ci, a] : list) {
            if (static_cast<int>(ci) < lo || static_cast<int>(ci) > hi) continue;
            int& slot = prof->alleles[ci - static_cast<size_t>(lo)];
            // The MSA already called this read here; realignment only fills reads
            // it did not call.
            if (slot >= 0) continue;
            if (slot != a) { slot = a; ++changed; }
        }
    }
    if (changed > 0) rebuild_read_var_cr(chunk);
    return changed;
}

size_t fill_missing_observations(GraphChunkBuildResult& graph_chunk, const PhasingChunk& bam,
                                 const std::vector<char>& only_reads) {
    PhasingChunk& chunk = graph_chunk.chunk;
    if (bam.ref_seq.empty()) return 0;
    const ReferenceView ref{bam};
    const SpanningReadIndex index(bam, chunk);
    const auto upper = [](std::string s) {
        for (char& ch : s) ch = static_cast<char>(std::toupper(static_cast<unsigned char>(ch)));
        return s;
    };
    const auto is_acgt = [](const std::string& s) {
        return s.find_first_not_of("ACGTacgt") == std::string::npos;
    };
    const bool have_meta = graph_chunk.site_meta.size() == chunk.candidates.size() &&
                           graph_chunk.site_allele_orig_idx.size() == chunk.candidates.size();
    std::vector<ReadVariantProfile*> profile_of(chunk.reads.size(), nullptr);
    for (ReadVariantProfile& p : chunk.read_var_profile)
        if (p.read_id >= 0 && static_cast<size_t>(p.read_id) < profile_of.size()) profile_of[p.read_id] = &p;
    const auto has_call = [&](size_t ri, size_t ci) {
        const ReadVariantProfile* p = profile_of[ri];
        if (p == nullptr || p->start_var_idx < 0 || static_cast<int>(ci) < p->start_var_idx ||
            static_cast<int>(ci) > p->end_var_idx) return false;
        return p->alleles[ci - static_cast<size_t>(p->start_var_idx)] >= 0;
    };
    std::map<size_t, std::vector<std::pair<size_t, int>>> fills;  // read -> (candidate, allele)
    std::vector<LocusAllele> locus_alleles;
    for (const CandidateVariant& c : chunk.candidates) {
        if (!c.bam_injected || c.key.type == VariantType::Snp || !c.msa_insertion_alts.empty()) continue;
        if (!is_acgt(c.key.alt) || (c.key.ref_len == 0 && c.key.alt.empty())) continue;
        locus_alleles.push_back(LocusAllele{c.key.pos - 1, c.key.pos - 1 + c.key.ref_len, c.key.alt});
    }
    std::sort(locus_alleles.begin(), locus_alleles.end(), [](const LocusAllele& a, const LocusAllele& b) { return a.p0 < b.p0; });
    for (size_t ci = 0; ci < chunk.candidates.size(); ++ci) {
        const CandidateVariant& c = chunk.candidates[ci];
        if ((c.lcd_var_i_to_cate & kCandGermlineVarCate) == 0) continue;
        if (!c.msa_insertion_alts.empty()) {
            if (!c.bam_injected || c.key.type != VariantType::Insertion) continue;
            const hts_pos_t p0 = c.key.pos - 1;
            const hts_pos_t w0 = std::max<hts_pos_t>(0, p0 - ref.repeat_left(p0) - kAlleleWindowFlank);
            const hts_pos_t w1 = p0 + ref.repeat_right(p0) + kAlleleWindowFlank;
            if (w1 - w0 > kLocusMaxLength) continue;
            const std::vector<std::string> haps = msa_allele_haplotypes(ref, c, w0, w1);
            if (haps.size() < 3) continue;
            for (const auto& [ri, seq] : index.segments(w0, w1, kAlleleWindowMaxReads)) {
                if (ri >= only_reads.size() || !only_reads[ri] || chunk.reads[ri].is_skipped || has_call(ri, ci)) continue;
                const int a = closest_haplotype(seq, haps, kLocusMinMargin);
                if (a >= 0) fills[ri].emplace_back(ci, a);
            }
            continue;
        }
        // The site's reference span and its two candidate alleles as sequence.
        hts_pos_t p0 = 0;
        int ref_len = 0;
        std::string seq0, seq1;
        std::vector<std::string> decomposed_others;
        if (c.bam_injected) {
            if (!is_acgt(c.key.alt) || (c.key.ref_len == 0 && c.key.alt.empty())) continue;
            p0 = c.key.pos - 1;
            ref_len = c.key.ref_len;
            seq0 = ref.slice(p0, p0 + ref_len);
            seq1 = upper(c.key.alt);
        } else {
            if (!have_meta) continue;
            const GraphSiteMeta& meta = graph_chunk.site_meta[ci];
            const std::vector<int>& orig = graph_chunk.site_allele_orig_idx[ci];
            if (meta.pos <= 0 || meta.ref.empty() || orig.size() < 2) continue;
            const auto allele_seq = [&meta](int k) -> const std::string* {
                if (k == 0) return &meta.ref;
                return k >= 1 && static_cast<size_t>(k) <= meta.alts.size() ? &meta.alts[static_cast<size_t>(k) - 1]
                                                                            : nullptr;
            };
            const std::string* a0 = allele_seq(orig[0]);
            const std::string* a1 = allele_seq(orig[1]);
            if (a0 == nullptr || a1 == nullptr || !is_acgt(*a0) || !is_acgt(*a1)) continue;
            p0 = meta.pos - 1;
            ref_len = static_cast<int>(meta.ref.size());
            if (ref.slice(p0, p0 + ref_len) != upper(meta.ref)) continue;
            seq0 = upper(*a0);
            seq1 = upper(*a1);
            // Allele 0 of a decomposed multi-allelic row is "any other allele of
            // the site": every one of them competes with the selected allele.
            if (meta.non_selected_alt_class) {
                decomposed_others.push_back(upper(meta.ref));
                bool usable = is_acgt(meta.ref);
                for (size_t a = 0; a < meta.alts.size(); ++a) {
                    if (static_cast<int>(a) + 1 == orig[1]) continue;
                    if (!is_acgt(meta.alts[a])) { usable = false; break; }
                    decomposed_others.push_back(upper(meta.alts[a]));
                }
                if (!usable) continue;
            }
        }
        if (seq0 == seq1) continue;
        const hts_pos_t end0 = p0 + ref_len;
        const bool snp = seq0.size() == 1 && seq1.size() == 1 && ref_len == 1;
        const hts_pos_t w0 = std::max<hts_pos_t>(0, p0 - (snp ? 0 : ref.repeat_left(p0)) - kAlleleWindowFlank);
        const hts_pos_t w1 = end0 + (snp ? 0 : ref.repeat_right(end0)) + kAlleleWindowFlank;
        if (w1 - w0 > kLocusMaxLength) continue;
        const std::string left = ref.slice(w0, p0), right = ref.slice(end0, w1);
        const std::string hap0 = left + seq0 + right, hap1 = left + seq1 + right;
        if (hap0.find('N') != std::string::npos) continue;
        std::vector<std::string> competitors{hap0};
        if (c.bam_injected && !snp)
            competitors = competing_haplotypes(ref, w0, w1, LocusAllele{p0, end0, seq1}, locus_alleles);
        for (const std::string& o : decomposed_others) competitors.push_back(left + o + right);
        for (const auto& [ri, seq] : index.segments(w0, w1, kAlleleWindowMaxReads)) {
            if (ri >= only_reads.size() || !only_reads[ri] || chunk.reads[ri].is_skipped || has_call(ri, ci)) continue;
            // A one-edit difference decides a SNP; at an indel it is as often a
            // homopolymer length error as the allele.
            const int a = closest_allele_call(seq, hap1, competitors, snp ? 1 : kLocusMinMargin);
            if (a >= 0) fills[ri].emplace_back(ci, a);
        }
    }
    // Write the calls, widening a profile's index range where needed.
    size_t added = 0;
    for (auto& [ri, calls] : fills) {
        ReadVariantProfile* p = profile_of[ri];
        if (p == nullptr) continue;
        int lo = p->start_var_idx, hi = p->end_var_idx;
        for (const auto& [ci, a] : calls) {
            lo = lo < 0 ? static_cast<int>(ci) : std::min(lo, static_cast<int>(ci));
            hi = std::max(hi, static_cast<int>(ci));
        }
        widen_profile(*p, lo, hi);
        for (const auto& [ci, a] : calls) {
            p->alleles[ci - static_cast<size_t>(lo)] = a;
            ++added;
        }
    }
    if (added > 0) rebuild_read_var_cr(chunk);
    return added;
}

std::vector<char> add_alignment_only_reads(PhasingChunk& chunk, const PhasingChunk& bam) {
    std::unordered_set<std::string_view> present;
    present.reserve(chunk.reads.size());
    for (const ReadRecord& r : chunk.reads) present.insert(r.qname);
    size_t added = 0;
    std::vector<char> is_added(chunk.reads.size(), 0);
    for (const ReadRecord& source : bam.reads) {
        if (!source.alignment || source.is_skipped) continue;
        if ((source.alignment->core.flag & (BAM_FSECONDARY | BAM_FSUPPLEMENTARY | BAM_FDUP | BAM_FUNMAP)) != 0)
            continue;
        if (source.alignment->core.qual < kAddedReadMinMapq || present.count(source.qname) != 0) continue;
        ReadRecord read;
        read.tid = chunk.region.tid;
        read.qname = source.qname;
        read.mapq = source.mapq;
        read.beg = source.beg;
        read.end = source.end;
        read.reverse = source.reverse;
        ReadVariantProfile prof;
        prof.read_id = static_cast<int>(chunk.reads.size());
        chunk.reads.push_back(std::move(read));
        chunk.read_var_profile.push_back(std::move(prof));
        present.insert(chunk.reads.back().qname);
        is_added.push_back(1);
        ++added;
    }
    chunk.haps.resize(chunk.reads.size(), 0);
    chunk.phase_sets.resize(chunk.reads.size(), kUnphasedReadPhaseSet);
    if (added > 0) rebuild_read_var_cr(chunk);
    return is_added;
}

} // namespace pgphase_collect
