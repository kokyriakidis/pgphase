// Union gap phasing: alignment-only reads, their fills, and realignment of
// uncalled reads at injected indel sites.

#include <unordered_map>
#include <string_view>
#include <array>
#include <map>
#include <set>
#include <tuple>
#include "union_internal.hpp"
#include "trace.hpp"

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

// A multi-allelic MSA row: co-located insertions (merge_msa_insertion_alleles,
// empty reference span) or co-located deletions (merge_msa_colocated_deletions:
// the longer deletion's span, each allele the bases it keeps). Its window is the
// span extended through tandem repeats plus the flank; false when too long.
static bool msa_row_window(const ReferenceView& ref, const CandidateVariant& c, hts_pos_t* w0, hts_pos_t* w1) {
    if (c.key.type != VariantType::Insertion && c.key.type != VariantType::Deletion) return false;
    const hts_pos_t p0 = c.key.pos - 1;
    const hts_pos_t end0 = c.key.type == VariantType::Deletion ? p0 + c.key.ref_len : p0;
    *w0 = std::max<hts_pos_t>(0, p0 - ref.repeat_left(p0) - kAlleleWindowFlank);
    *w1 = end0 + ref.repeat_right(end0) + kAlleleWindowFlank;
    return *w1 - *w0 <= kLocusMaxLength;
}

// Haplotype sequences of a multi-allelic MSA row over [w0, w1): index 0 the
// reference, index i the reference with the row's span replaced by
// msa_insertion_alts[i-1].
static std::vector<std::string> msa_allele_haplotypes(const ReferenceView& ref, const CandidateVariant& c,
                                                      hts_pos_t w0, hts_pos_t w1) {
    const hts_pos_t p0 = c.key.pos - 1;
    const hts_pos_t end0 = c.key.type == VariantType::Deletion ? p0 + c.key.ref_len : p0;
    std::vector<std::string> haps{ref.slice(w0, w1)};
    for (const std::string& allele : c.msa_insertion_alts) {
        if (allele.find_first_not_of("ACGTacgt") != std::string::npos) return {};
        haps.push_back(ref.slice(w0, p0) + allele + ref.slice(end0, w1));
    }
    return haps;
}

// An MSA call's two haplotypes over its region: each consensus sequence (MsaRegion::spelled) between the reference
// flanks, for a call that is one allele of each consensus (msa_cons_allele) in a region the two differ in, short enough
// to realign reads over. -1 otherwise, else the region's index.
static long msa_region_haplotypes(const ReferenceView& ref, const PhasingChunk& bam, const CandidateVariant& c,
                                  hts_pos_t* w0, hts_pos_t* w1, std::array<std::string, 2>* haps) {
    const std::array<int8_t, 2>& allele = c.msa_cons_allele;
    if (!c.msa_verified || allele[0] < 0 || allele[1] < 0 || allele[0] + allele[1] != 1) return -1;
    for (size_t i = 0; i < bam.msa_regions.size(); ++i) {
        const PhasingChunk::MsaRegion& region = bam.msa_regions[i];
        if (!region.het() || region.beg >= c.key.pos || c.key.pos >= region.end) continue;
        const hts_pos_t r0 = region.beg - 1, r1 = region.end;  // 0-based, half open
        *w0 = std::max<hts_pos_t>(0, r0 - kAlleleWindowFlank);
        *w1 = r1 + kAlleleWindowFlank;
        if (*w1 - *w0 > kLocusMaxLength) return -1;
        for (size_t k = 0; k < 2; ++k) {
            std::string seq = ref.slice(*w0, r0);
            for (const std::string& at : region.spelled[k]) seq += at;
            seq += ref.slice(r1, *w1);
            (*haps)[k] = std::move(seq);
        }
        return (*haps)[0] != (*haps)[1] ? static_cast<long>(i) : -1;
    }
    return -1;
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

// A read's realignment margin m (edits between the called allele and the
// closest competitor) sets its error 1 / (1 + exp(lambda * m)); where the MSA
// called the read and realignment cannot decide, or decides the other way,
// the MSA call is kept at these errors.
constexpr double kMarginErrorLambda = 3.0;
constexpr double kMsaUndecidedError = 0.25;
constexpr double kMsaContradictedError = 0.35;

static double margin_error(int margin) {
    return 1.0 / (1.0 + std::exp(kMarginErrorLambda * margin));
}

size_t realign_indel_observations(PhasingChunk& chunk, const PhasingChunk& bam, ObservationErrors* errors) {
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
    std::map<std::pair<size_t, size_t>, int> margins;               // (read, candidate) -> edit margin
    std::vector<LocusAllele> locus_alleles;
    for (const CandidateVariant& c : chunk.candidates) {
        if (!c.bam_injected || c.key.type == VariantType::Snp || !c.msa_insertion_alts.empty()) continue;
        if (!is_acgt(c.key.alt) || (c.key.ref_len == 0 && c.key.alt.empty())) continue;
        locus_alleles.push_back(LocusAllele{c.key.pos - 1, c.key.pos - 1 + c.key.ref_len, c.key.alt});
    }
    std::sort(locus_alleles.begin(), locus_alleles.end(), [](const LocusAllele& a, const LocusAllele& b) { return a.p0 < b.p0; });
    // Per resolved MSA region: each spanning read's closer consensus (or -1) and its edit margin.
    std::map<long, std::vector<std::tuple<size_t, int, int>>> region_reads;
    for (size_t ci = 0; ci < chunk.candidates.size(); ++ci) {
        const CandidateVariant& c = chunk.candidates[ci];
        // An MSA call inside a region whose two consensus sequences differ: a read the MSA did not cluster is called
        // by the consensus it is closer to over the whole region, not by REF against this call's own allele over a
        // window. chr20:4,022,917, A>G on one haplotype inside the other's 6 bp deletion of a CA repeat (the truth
        // A>G,*): realigned site by site, the deletion haplotype's reads read the G of 4,022,923 (three CA units on)
        // there, 49 ALT to 9, and the EM's minor-allele floor dropped a het the MSA had right (clusters 12 / 9).
        if (c.bam_injected && (c.lcd_var_i_to_cate & kCandGermlineVarCate) != 0) {
            hts_pos_t w0 = 0, w1 = 0;
            std::array<std::string, 2> haps;
            const long region = msa_region_haplotypes(ref, bam, c, &w0, &w1, &haps);
            if (region >= 0) {
                auto [it, fresh] = region_reads.try_emplace(region);
                if (fresh)
                    for (const auto& [ri, seq] : index.segments(w0, w1, kAlleleWindowMaxReads)) {
                        if (profile_of[ri] == nullptr) continue;
                        const int d0 = locus_edit_distance(seq, haps[0]);
                        const int d1 = locus_edit_distance(seq, haps[1]);
                        it->second.emplace_back(ri, d0 == d1 ? -1 : d0 < d1 ? 0 : 1, std::abs(d0 - d1));
                    }
                for (const auto& [ri, k, margin] : it->second) {
                    calls[ri].emplace_back(ci, k < 0 ? -1 : c.msa_cons_allele[static_cast<size_t>(k)]);
                    if (errors != nullptr && k >= 0) margins[{ri, ci}] = margin;
                }
                continue;
            }
        }
        // MSA SNP calls too: they carry only the reads inside the MSA's
        // clusters (chr20:35,211,758: 4 of 55 crossing reads), and realigning
        // a read to REF and to the SNP's haplotype fills the rest the same way.
        // A pileup SNP is not filled: its reads without a call are mostly the
        // ones whose base there failed the base-quality filter.
        if (!c.bam_injected || (c.key.type == VariantType::Snp && !c.msa_verified)) continue;
        if (!c.msa_insertion_alts.empty()) {
            // One row, several insertion or deletion alleles: the closest
            // allele by index. Deletion rows were skipped before, leaving a
            // merged 1/2 deletion with only the MSA's clustered reads
            // (chr18:77,670,828: 11 of ~60 crossing reads).
            if ((c.lcd_var_i_to_cate & kCandGermlineVarCate) == 0) continue;
            hts_pos_t w0 = 0, w1 = 0;
            if (!msa_row_window(ref, c, &w0, &w1)) continue;
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
        // Distances to the ALT haplotype and the closest competitor, once per
        // distinct read sequence: the call (closest_allele_call, margin 1) and
        // the margin for the read's error come from the same pair.
        std::unordered_map<std::string_view, std::pair<int, int>> seen;
        for (const auto& [ri, seq] : index.segments(w0, w1, kAlleleWindowMaxReads)) {
            if (profile_of[ri] == nullptr) continue;
            auto hit = seen.find(seq);
            if (hit == seen.end()) {
                const int d1 = locus_edit_distance(seq, hap_alt);
                int d0 = -1;
                for (const std::string& h : competitors) {
                    if (h == hap_alt) continue;
                    const int d = locus_edit_distance(seq, h);
                    if (d0 < 0 || d < d0) d0 = d;
                }
                hit = seen.emplace(seq, std::make_pair(d1, d0)).first;
            }
            const auto [d1, d0] = hit->second;
            calls[ri].emplace_back(ci, d0 < 0 || std::abs(d0 - d1) < 1 ? -1 : d1 < d0 ? 1 : 0);
            if (errors != nullptr && d0 >= 0) margins[{ri, ci}] = std::abs(d0 - d1);
        }
    }
    // Every read spanning the site's window is re-called there, its range
    // widened if its other calls stop short of the site.
    size_t changed = 0;
    // Trace: per site, calls set (the MSA left the read uncalled) and MSA calls graded agreeing, undecided or
    // contradicted.
    std::map<size_t, std::array<int, 4>> site_tally;
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
            const auto mg = margins.find({ri, ci});
            // The MSA already called this read here; realignment only fills reads
            // it did not call, and grades the calls it did.
            if (slot >= 0) {
                // A SNP's existing call is a base call; realignment over a
                // window that also spans repeat indels is often undecided
                // there, and grading the call by it down-weighted good SNPs.
                if (chunk.candidates[ci].key.type == VariantType::Snp) continue;
                if (errors != nullptr && mg != margins.end())
                    (*errors)[observation_key(ri, ci)] =
                        a < 0 ? kMsaUndecidedError : a != slot ? kMsaContradictedError : margin_error(mg->second);
                if (pgphase_trace::on() && errors != nullptr && mg != margins.end()) {
                    ++site_tally[ci][a < 0 ? 2 : a != slot ? 3 : 1];
                    ++pgphase_trace::read_trace(chunk, ri)->graded;
                }
                continue;
            }
            if (slot != a) {
                slot = a;
                ++changed;
                if (pgphase_trace::on()) {
                    ++site_tally[ci][0];
                    ++pgphase_trace::read_trace(chunk, ri)->realigned;
                }
            }
            if (errors != nullptr && a >= 0 && mg != margins.end())
                (*errors)[observation_key(ri, ci)] = margin_error(mg->second);
        }
    }
    for (const auto& [ci, n] : site_tally) {
        const auto [beg, end] = pgphase_trace::span(chunk.candidates[ci]);
        pgphase_trace::emit("reads", "realign", beg, end, pgphase_trace::key_text(chunk.candidates[ci].key),
                            "calls_rewritten",
                            "set=" + std::to_string(n[0]) + " graded_agree=" + std::to_string(n[1]) +
                                " undecided=" + std::to_string(n[2]) + " contradicted=" + std::to_string(n[3]));
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
            if (!c.bam_injected) continue;
            hts_pos_t w0 = 0, w1 = 0;
            if (!msa_row_window(ref, c, &w0, &w1)) continue;
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
        // One call per distinct read sequence.
        std::unordered_map<std::string_view, int> seen;
        for (const auto& [ri, seq] : index.segments(w0, w1, kAlleleWindowMaxReads)) {
            if (ri >= only_reads.size() || !only_reads[ri] || chunk.reads[ri].is_skipped || has_call(ri, ci)) continue;
            // A one-edit difference decides a SNP; at an indel it is as often a
            // homopolymer length error as the allele.
            auto hit = seen.find(seq);
            if (hit == seen.end())
                hit = seen.emplace(seq, closest_allele_call(seq, hap1, competitors, snp ? 1 : kLocusMinMargin)).first;
            if (hit->second >= 0) fills[ri].emplace_back(ci, hit->second);
        }
    }
    // Write the calls, widening a profile's index range where needed.
    size_t added = 0;
    std::map<size_t, int> site_fills;  // trace: calls given per site
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
            if (pgphase_trace::on()) {
                ++site_fills[ci];
                ++pgphase_trace::read_trace(chunk, ri)->filled;
            }
        }
    }
    for (const auto& [ci, n] : site_fills) {
        const auto [beg, end] = pgphase_trace::span(chunk.candidates[ci]);
        pgphase_trace::emit("reads", "fill", beg, end, pgphase_trace::key_text(chunk.candidates[ci].key),
                            "alignment_only_reads_called", "n=" + std::to_string(n));
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
        if (pgphase_trace::on()) {
            PhasingChunk::ReadTrace* t = pgphase_trace::read_trace(chunk, chunk.reads.size() - 1);
            t->alignment_only = true;
            const size_t si = static_cast<size_t>(&source - bam.reads.data());
            t->aln_hap = si < bam.haps.size() ? bam.haps[si] : -1;
            t->aln_ps = si < bam.phase_sets.size() ? bam.phase_sets[si] : 0;
        }
    }
    if (pgphase_trace::on() && added > 0)
        pgphase_trace::emit("reads", "chunk", chunk.region.beg, chunk.region.end, "-", "alignment_only_reads_added",
                            "n=" + std::to_string(added));
    chunk.haps.resize(chunk.reads.size(), 0);
    chunk.phase_sets.resize(chunk.reads.size(), kUnphasedReadPhaseSet);
    if (added > 0) rebuild_read_var_cr(chunk);
    return is_added;
}

// The alignment's call at a catalog SNP outranks the graph's. Where both
// called a read and disagreed, the graph's base was right 7-12% of the time
// (chr12 / chr18 / chr20, reads we phase correctly); where the alignment
// calls the site but not that read (base quality, read end), the graph's
// call was right 77-79% of the time, against 99.9% where the two agree.
constexpr double kGraphOnlyCallError = 0.25;

// Below it a read's base is not its call.
constexpr int kAlignedBaseMinQuality = 10;

// The read's base at 1-based reference position pos by its alignment record, upper case; 0 if the alignment does not
// cover it (or the record is gone), deletes it, or reads it below kAlignedBaseMinQuality.
static char aligned_base(const ReadRecord& read, hts_pos_t pos) {
    const bam1_t* b = read.alignment.get();
    if (b == nullptr) return 0;
    const hts_pos_t ref_pos = pos - 1;
    if (ref_pos < b->core.pos) return 0;
    hts_pos_t r = b->core.pos, q = 0;
    const uint32_t* cigar = bam_get_cigar(b);
    for (uint32_t i = 0; i < b->core.n_cigar; ++i) {
        const int op = bam_cigar_op(cigar[i]);
        const hts_pos_t len = bam_cigar_oplen(cigar[i]);
        if (op == BAM_CMATCH || op == BAM_CEQUAL || op == BAM_CDIFF) {
            if (r + len > ref_pos) {
                const hts_pos_t qi = q + (ref_pos - r);
                if (bam_get_qual(b)[qi] < kAlignedBaseMinQuality) return 0;
                return seq_nt16_str[bam_seqi(bam_get_seq(b), qi)];
            }
            r += len;
            q += len;
        } else if (op == BAM_CDEL || op == BAM_CREF_SKIP) {
            if (r + len > ref_pos) return 0;
            r += len;
        } else if (op == BAM_CINS || op == BAM_CSOFT_CLIP) {
            q += len;
        }
    }
    return 0;
}

size_t confirm_catalog_snp_calls(GraphChunkBuildResult& graph_chunk, const PhasingChunk& bam,
                                 ObservationErrors* errors) {
    PhasingChunk& chunk = graph_chunk.chunk;
    if (graph_chunk.site_meta.size() != chunk.candidates.size() ||
        graph_chunk.site_allele_orig_idx.size() != chunk.candidates.size()) return 0;
    std::unordered_map<hts_pos_t, std::vector<size_t>> bam_snps;  // VCF position -> alignment SNP rows
    for (size_t bi = 0; bi < bam.candidates.size(); ++bi) {
        const CandidateVariant& c = bam.candidates[bi];
        if (c.key.type == VariantType::Snp && c.key.alt.size() == 1) bam_snps[c.key.pos].push_back(bi);
    }
    std::unordered_map<std::string_view, const ReadVariantProfile*> bam_profile;
    for (const ReadVariantProfile& p : bam.read_var_profile)
        if (p.read_id >= 0 && p.start_var_idx >= 0)
            bam_profile.emplace(bam.reads[static_cast<size_t>(p.read_id)].qname, &p);
    // A read's alignment records (a supplementary one beside the primary): its base is the primary's where that covers
    // the position, else the covering record's that starts first -- whichever record came first in the chunk differed
    // with the thread count.
    std::unordered_map<std::string_view, std::vector<const ReadRecord*>> bam_read;
    for (const ReadRecord& r : bam.reads) bam_read[r.qname].push_back(&r);
    const auto read_base = [&bam_read](std::string_view qname, hts_pos_t pos) -> char {
        const auto it = bam_read.find(qname);
        if (it == bam_read.end()) return 0;
        const ReadRecord* best = nullptr;
        for (const ReadRecord* r : it->second) {
            const bam1_t* b = r->alignment.get();
            if (b == nullptr || b->core.pos > pos - 1 || bam_endpos(b) < pos) continue;
            const bool primary = (b->core.flag & (BAM_FSECONDARY | BAM_FSUPPLEMENTARY)) == 0;
            const bool best_primary = best != nullptr &&
                                      (best->alignment->core.flag & (BAM_FSECONDARY | BAM_FSUPPLEMENTARY)) == 0;
            if (best == nullptr || (primary && !best_primary) ||
                (primary == best_primary && b->core.pos < best->alignment->core.pos))
                best = r;
        }
        return best != nullptr ? aligned_base(*best, pos) : 0;
    };
    const auto upper = [](char b) { return static_cast<char>(std::toupper(static_cast<unsigned char>(b))); };
    size_t changed = 0;
    std::map<size_t, std::array<int, 2>> site_tally;  // trace: calls replaced, graph-only calls down-weighted
    std::set<size_t> recount;  // sites whose calls the alignment's replaced (below)
    // trace, sites the alignment has no SNP at: reads absent from it, not read there, read
    std::map<size_t, std::array<int, 3>> aligned_tally;
    for (size_t ri = 0; ri < chunk.read_var_profile.size(); ++ri) {
        ReadVariantProfile& prof = chunk.read_var_profile[ri];
        if (prof.start_var_idx < 0 || prof.read_id < 0) continue;
        const std::string& qname = chunk.reads[static_cast<size_t>(prof.read_id)].qname;
        const auto found = bam_profile.find(qname);
        const ReadVariantProfile* bp = found == bam_profile.end() ? nullptr : found->second;
        for (size_t k = 0; k < prof.alleles.size(); ++k) {
            const int a = prof.alleles[k];
            if (a < 0) continue;
            const size_t ci = static_cast<size_t>(prof.start_var_idx) + k;
            if (chunk.candidates[ci].bam_injected) continue;
            const GraphSiteMeta& meta = graph_chunk.site_meta[ci];
            const std::vector<int>& orig = graph_chunk.site_allele_orig_idx[ci];
            if (meta.ref.size() != 1 || static_cast<size_t>(a) >= orig.size()) continue;
            bool snp_row = !meta.alts.empty();
            for (const std::string& alt : meta.alts) snp_row = snp_row && alt.size() == 1;
            if (!snp_row) continue;
            const auto site = bam_snps.find(meta.pos);
            // The read's base by the alignment, when it called the read here. Where the alignment has no SNP here at
            // all, its reads carry one base, the reference's as a rule: the read's own aligned base then
            // (chr20:1,012,629 C>G, 20 of 68 reads on the graph's ALT path, all 68 C in their alignments; written a
            // false 0/1).
            char base = 0;
            if (site == bam_snps.end()) {
                base = read_base(qname, meta.pos);
                if (pgphase_trace::on())
                    ++aligned_tally[ci][bam_read.count(qname) == 0 ? 0 : base == 0 ? 1 : 2];
                if (base == 0) continue;
            }
            for (const size_t bi : site == bam_snps.end() ? std::vector<size_t>{} : site->second) {
                const int b = bp != nullptr && static_cast<int>(bi) >= bp->start_var_idx &&
                              static_cast<int>(bi) <= bp->end_var_idx
                                  ? bp->alleles[bi - static_cast<size_t>(bp->start_var_idx)] : -1;
                if (b == 1) { base = upper(bam.candidates[bi].key.alt[0]); break; }
                if (b == 0 && site->second.size() == 1) base = upper(meta.ref[0]);
            }
            if (base == 0) {
                if (errors != nullptr) {
                    double& e = (*errors)[observation_key(ri, ci)];
                    e = std::max(e, kGraphOnlyCallError);
                    if (pgphase_trace::on()) ++site_tally[ci][1];
                }
                continue;
            }
            // Does graph allele j carry `base`?
            const auto carries = [&](size_t j) {
                const int o = orig[j];
                if (o == 0) return upper(meta.ref[0]) == base;
                return o >= 1 && static_cast<size_t>(o) <= meta.alts.size() &&
                       upper(meta.alts[static_cast<size_t>(o) - 1][0]) == base;
            };
            if (carries(static_cast<size_t>(a))) continue;
            int replacement = -1;
            for (size_t j = 0; j < orig.size(); ++j)
                if (carries(j)) { replacement = static_cast<int>(j); break; }
            prof.alleles[k] = replacement;
            ++changed;
            recount.insert(ci);
            if (pgphase_trace::on()) {
                ++site_tally[ci][0];
                ++pgphase_trace::read_trace(chunk, static_cast<size_t>(prof.read_id))->snp_replaced;
            }
        }
    }
    // A site whose calls the alignment's replaced is genotyped at output from its reads' calls as they now stand, not
    // the graph's: written from the graph's counts, chr20:319,724 T>C (10 of 35 reads on the graph's REF path, all 35
    // T>C in their alignments) was a false 0/1 beside a truth 1/1, and 1,012,629 C>G (20 of 68 on its ALT path, none
    // in their alignments) a false 0/1.
    if (!recount.empty()) {
        std::map<size_t, std::vector<int>> covs;
        for (const size_t ci : recount)
            covs[ci].assign(std::max<size_t>(2, chunk.candidates[ci].counts.alle_covs.size()), 0);
        for (const ReadVariantProfile& prof : chunk.read_var_profile) {
            if (prof.start_var_idx < 0 || prof.read_id < 0) continue;
            for (size_t k = 0; k < prof.alleles.size(); ++k) {
                const auto it = covs.find(static_cast<size_t>(prof.start_var_idx) + k);
                const int a = prof.alleles[k];
                if (it != covs.end() && a >= 0 && static_cast<size_t>(a) < it->second.size()) ++it->second[a];
            }
        }
        for (auto& [ci, c] : covs) {
            VariantCounts& counts = chunk.candidates[ci].counts;
            if (pgphase_trace::on()) {
                std::string before, after;
                for (const int n : counts.alle_covs) before += (before.empty() ? "" : ",") + std::to_string(n);
                for (const int n : c) after += (after.empty() ? "" : ",") + std::to_string(n);
                const auto [beg, end] = pgphase_trace::span(chunk.candidates[ci]);
                pgphase_trace::emit("reads", "snp_call", beg, end, pgphase_trace::key_text(chunk.candidates[ci].key),
                                    "recounted", "before=" + before + " after=" + after);
            }
            counts.alle_covs = std::move(c);
        }
    }
    for (const auto& [ci, n] : site_tally) {
        const auto [beg, end] = pgphase_trace::span(chunk.candidates[ci]);
        pgphase_trace::emit("reads", "snp_call", beg, end, pgphase_trace::key_text(chunk.candidates[ci].key),
                            "graph_call_vs_alignment",
                            "replaced=" + std::to_string(n[0]) + " graph_only_downweighted=" + std::to_string(n[1]));
    }
    for (const auto& [ci, n] : aligned_tally) {
        const auto [beg, end] = pgphase_trace::span(chunk.candidates[ci]);
        pgphase_trace::emit("reads", "snp_call", beg, end, pgphase_trace::key_text(chunk.candidates[ci].key),
                            "graph_call_vs_read_alignment",
                            "no_read=" + std::to_string(n[0]) + " not_read=" + std::to_string(n[1]) +
                                " read=" + std::to_string(n[2]) +
                                " replaced=" + std::to_string(site_tally.count(ci) ? site_tally[ci][0] : 0));
    }
    return changed;
}

} // namespace pgphase_collect
