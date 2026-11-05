/**
 * @file noisy_regions.cpp
 * Interval trees (cgranges) and noisy-region construction, extension, splitting, and containment filtering.
 */

#include "collect_var.hpp"
#include "collect_phase.hpp"
#include "noise_filter.hpp"
#include "plegma_msa.hpp"
#include "trace.hpp"

#include <algorithm>
#include <cctype>
#include <cinttypes>
#include <cstdio>
#include <cstdlib>
#include <limits>
#include <string>
#include <utility>
#include <vector>

namespace pgphase_collect {

// ════════════════════════════════════════════════════════════════════════════
// Interval utilities
// ════════════════════════════════════════════════════════════════════════════

/**
 * @brief Builds an unindexed `cgranges_t` from 1-based inclusive intervals (stored 0-based half-open).
 *
 * Returns nullptr if empty or all intervals have `beg > end`. Caller must `cr_index` before overlap queries.
 * @param intervals 1-based inclusive source intervals.
 * @return New `cgranges_t*` or nullptr.
 */
cgranges_t* intervals_to_cr(const std::vector<Interval>& intervals) {
    if (intervals.empty()) return nullptr;
    cgranges_t* cr = cr_init();
    int added = 0;
    for (const Interval& iv : intervals) {
        if (iv.beg > iv.end) continue;
        const int32_t st = static_cast<int32_t>(iv.beg - 1);
        const int32_t en = static_cast<int32_t>(iv.end);
        const int32_t label = iv.label;
        cr_add(cr, "cr", st, en, label);
        ++added;
    }
    if (added == 0) {
        cr_destroy(cr);
        return nullptr;
    }
    // Leave unindexed: callers run cr_index once before overlap/merge.
    return cr;
}

/**
 * @brief Reads intervals from the synthetic `"cr"` contig into \a out (inverse of `intervals_to_cr`).
 * @param cr Indexed or unindexed tree using contig `"cr"`.
 * @param out Destination vector (cleared first).
 */
void intervals_from_cr(const cgranges_t* cr, std::vector<Interval>& out) {
    out.clear();
    if (cr == nullptr || cr->n_r == 0) return;
    const int32_t ctg = cr_get_ctg(cr, "cr");
    if (ctg < 0) return;
    const cr_ctg_t* c = &cr->ctg[ctg];
    for (int64_t i = c->off; i < c->off + c->n; ++i) {
        const hts_pos_t beg = static_cast<hts_pos_t>(cr_start(cr, i) + 1);
        const hts_pos_t end = static_cast<hts_pos_t>(cr_end(cr, i));
        out.push_back(Interval{beg, end, cr_label(cr, i)});
    }
}

// Post-process noisy regions: read intervals from tree, convert to 1-based inclusive.
// noisy_reg_end[reg_i], ...)` stores `st == noisy_reg_start`, not `start - 1`.
static void intervals_from_cr_noisy_post_merge(const cgranges_t* cr, std::vector<Interval>& out) {
    out.clear();
    if (cr == nullptr || cr->n_r == 0) return;
    const int32_t ctg = cr_get_ctg(cr, "cr");
    if (ctg < 0) return;
    const cr_ctg_t* c = &cr->ctg[ctg];
    for (int64_t i = c->off; i < c->off + c->n; ++i) {
        const hts_pos_t beg = static_cast<hts_pos_t>(cr_start(cr, i));
        const hts_pos_t end = static_cast<hts_pos_t>(cr_end(cr, i));
        out.push_back(Interval{beg, end, cr_label(cr, i)});
    }
}

static cgranges_t* intervals_to_cr_noisy_post_merge(const std::vector<Interval>& intervals) {
    if (intervals.empty()) return nullptr;
    cgranges_t* cr = cr_init();
    int added = 0;
    for (const Interval& iv : intervals) {
        if (iv.beg > iv.end) continue;
        const int32_t st = static_cast<int32_t>(iv.beg);
        const int32_t en = static_cast<int32_t>(iv.end);
        const int32_t label = iv.label;
        cr_add(cr, "cr", st, en, label);
        ++added;
    }
    if (added == 0) {
        cr_destroy(cr);
        return nullptr;
    }
    return cr;
}

/**
 * @brief Indexed `cgranges_t` for one read's `noisy_regions` (overlap tests in noisy preprocessing).
 * @param read Read whose `noisy_regions` may be empty.
 * @return Indexed tree or nullptr.
 */
cgranges_t* build_read_noisy_cr(const ReadRecord& read) {
    if (read.noisy_regions.empty()) return nullptr;
    cgranges_t* cr = intervals_to_cr(read.noisy_regions);
    if (cr != nullptr) cr_index(cr);
    return cr;
}

/**
 * @brief Reference span for overlap tests (1-based inclusive).
 *
 * Insertions yield var_end == key.pos - 1 (zero-length span).
 * @param key Variant locus.
 * @param var_start Output first reference coordinate.
 * @param var_end Output last reference coordinate.
 */
void variant_genomic_span(const VariantKey& key, hts_pos_t& var_start, hts_pos_t& var_end) {
    if (key.type == VariantType::Snp) {
        var_start = var_end = key.pos;
        return;
    }
    if (key.type == VariantType::Deletion) {
        var_start = key.pos;
        var_end = key.pos + key.ref_len - 1;
        return;
    }
    // Insertions: use [pos, pos] for overlap/noisy-ratio checks.
    var_start = key.pos;
    var_end = key.pos;
}

// ════════════════════════════════════════════════════════════════════════════
// Noisy region processing
// ════════════════════════════════════════════════════════════════════════════

/**
 * @brief Expands \a start/\a end to the union with any overlapping low-complexity intervals.
 * @param low_comp_cr Indexed low-complexity tree or nullptr.
 * @param start Reference start (1-based, converted internally).
 * @param end Reference end.
 * @param new_start Output expanded start.
 * @param new_end Output expanded end.
 */
static void low_comp_cr_start_end(cgranges_t* low_comp_cr,
                                  int32_t start, int32_t end,
                                  int32_t* new_start, int32_t* new_end) {
    *new_start = start;
    *new_end = end;
    if (low_comp_cr == nullptr || low_comp_cr->n_r == 0) return;
    int64_t* low_comp_b = nullptr;
    int64_t max_low_comp_n = 0;
    const int64_t low_comp_n =
        cr_overlap(low_comp_cr, "cr", start - 1, end, &low_comp_b, &max_low_comp_n);
    for (int64_t j = 0; j < low_comp_n; ++j) {
        const int32_t s = cr_start(low_comp_cr, low_comp_b[j]) + 1;
        const int32_t e = cr_end(low_comp_cr, low_comp_b[j]);
        if (s < *new_start) *new_start = s;
        if (e > *new_end) *new_end = e;
    }
    std::free(low_comp_b);
}

/**
 * Extend noisy regions into adjacent low-complexity intervals, then merge.
 * @param chunk_noisy_regs Unindexed or indexed noisy tree (consumed/replaced when low-comp extends).
 * @param low_comp_cr Low-complexity intervals or nullptr.
 * @param merge_dis Merge distance for `cr_merge`.
 * @param min_sv_len `min_sv_len` argument to `cr_merge`.
 * @return Merged `cgranges_t*` (caller owns prior tree if replaced).
 */
static cgranges_t* cr_extend_noisy_regs_with_low_comp(cgranges_t* chunk_noisy_regs,
                                                      cgranges_t* low_comp_cr,
                                                      int merge_dis,
                                                      int min_sv_len) {
    if (chunk_noisy_regs == nullptr) return nullptr;
    cgranges_t* cur = chunk_noisy_regs;
    if (low_comp_cr != nullptr && low_comp_cr->n_r > 0) {
        cgranges_t* new_noisy_regs = cr_init();
        for (int i = 0; i < cur->n_r; ++i) {
            const int32_t start = cr_start(cur, i) + 1;
            const int32_t end = cr_end(cur, i);
            int32_t new_start = start;
            int32_t new_end = end;
            low_comp_cr_start_end(low_comp_cr, start, end, &new_start, &new_end);
            cr_add(new_noisy_regs, "cr", new_start - 1, new_end, cr_label(cur, i));
        }
        cr_index(new_noisy_regs);
        cr_destroy(cur);
        cur = new_noisy_regs;
    }
    return cr_merge(cur, -1, merge_dis, min_sv_len);
}

/**
 * @brief Categories excluded when flanking noisy regions with nearby clean variants.
 * @param c Final or intermediate category.
 */
static bool category_skipped_for_noisy_flank(VariantCategory c) {
    return c == VariantCategory::LowCoverage || c == VariantCategory::StrandBias ||
           c == VariantCategory::NonVariant;
}

/**
 * Per-noisy-region flanking bounds including nearest clean candidates.
 * @param chunk_noisy_regs Indexed noisy regions.
 * @param cand Candidates (sorted as in chunk).
 * @param cand_cate Parallel category per candidate (same length as `cand`).
 * @param start_out Output left bounds per noisy region.
 * @param end_out Output right bounds per noisy region.
 */
static void collect_noisy_reg_start_end_pgphase(const cgranges_t* chunk_noisy_regs,
                                                const CandidateTable& cand,
                                                const std::vector<VariantCategory>& cand_cate,
                                                std::vector<int32_t>& start_out,
                                                std::vector<int32_t>& end_out) {
    const auto flank_span = [](const VariantKey& key, hts_pos_t& var_start, hts_pos_t& var_end) {
        // Extend flanks to nearest clean candidate outside the noisy span:
        // var_start = var->pos; var_end = var->pos + var->ref_len - 1.
        var_start = key.pos;
        var_end = key.pos + key.ref_len - 1;
    };

    const int n_noisy = chunk_noisy_regs->n_r;
    start_out.assign(n_noisy, 0);
    end_out.assign(n_noisy, 0);
    if (cand.empty()) {
        for (int reg_i = 0; reg_i < n_noisy; ++reg_i) {
            const int32_t ori_reg_start = cr_start(chunk_noisy_regs, reg_i) + 1;
            const int32_t ori_reg_end = cr_end(chunk_noisy_regs, reg_i);
            start_out[static_cast<size_t>(reg_i)] = ori_reg_start - kNoisyRegFlankLen;
            end_out[static_cast<size_t>(reg_i)] = ori_reg_end + kNoisyRegFlankLen;
        }
        return;
    }
    std::vector<int> max_left_var_i(n_noisy, -1);
    std::vector<int> min_right_var_i(n_noisy, -1);

    int reg_i = 0;
    int var_i = 0;
    while (reg_i < n_noisy && var_i < static_cast<int>(cand.size())) {
        if (category_skipped_for_noisy_flank(cand_cate[static_cast<size_t>(var_i)])) {
            ++var_i;
            continue;
        }
        hts_pos_t var_start = 0;
        hts_pos_t var_end = 0;
        flank_span(cand[static_cast<size_t>(var_i)].key, var_start, var_end);
        const int32_t reg_start = cr_start(chunk_noisy_regs, reg_i) + 1;
        const int32_t reg_end = cr_end(chunk_noisy_regs, reg_i);
        if (static_cast<int32_t>(var_start) > reg_end) {
            if (min_right_var_i[static_cast<size_t>(reg_i)] == -1) {
                min_right_var_i[static_cast<size_t>(reg_i)] = var_i;
            }
            ++reg_i;
        } else if (static_cast<int32_t>(var_end) < reg_start) {
            max_left_var_i[static_cast<size_t>(reg_i)] = var_i;
            ++var_i;
        } else {
            ++var_i;
        }
    }

    for (reg_i = 0; reg_i < n_noisy; ++reg_i) {
        if (max_left_var_i[static_cast<size_t>(reg_i)] == -1) {
            max_left_var_i[static_cast<size_t>(reg_i)] =
                std::min(std::max(0, static_cast<int>(cand.size()) - 1), 0);
        }
        if (min_right_var_i[static_cast<size_t>(reg_i)] == -1) {
            min_right_var_i[static_cast<size_t>(reg_i)] =
                std::max(0, static_cast<int>(cand.size()) - 1);
        }
        const int32_t ori_reg_start = cr_start(chunk_noisy_regs, reg_i) + 1;
        const int32_t ori_reg_end = cr_end(chunk_noisy_regs, reg_i);
        int32_t cur_start = ori_reg_start - kNoisyRegFlankLen;
        for (int v = max_left_var_i[static_cast<size_t>(reg_i)]; v >= 0; --v) {
            if (category_skipped_for_noisy_flank(cand_cate[static_cast<size_t>(v)])) continue;
            hts_pos_t vs = 0;
            hts_pos_t ve = 0;
            flank_span(cand[static_cast<size_t>(v)].key, vs, ve);
            if (static_cast<int32_t>(ve) < cur_start - 1) break;
            if (static_cast<int32_t>(vs) - kNoisyRegFlankLen < cur_start) {
                cur_start = static_cast<int32_t>(vs) - kNoisyRegFlankLen;
            }
        }
        start_out[static_cast<size_t>(reg_i)] = cur_start;
        int32_t cur_end = ori_reg_end + kNoisyRegFlankLen;
        for (int v = min_right_var_i[static_cast<size_t>(reg_i)]; v < static_cast<int>(cand.size()); ++v) {
            if (category_skipped_for_noisy_flank(cand_cate[static_cast<size_t>(v)])) continue;
            hts_pos_t vs = 0;
            hts_pos_t ve = 0;
            flank_span(cand[static_cast<size_t>(v)].key, vs, ve);
            if (static_cast<int32_t>(vs) > cur_end + 1) break;
            if (static_cast<int32_t>(ve) + kNoisyRegFlankLen > cur_end) {
                cur_end = static_cast<int32_t>(ve) + kNoisyRegFlankLen;
            }
        }
        end_out[static_cast<size_t>(reg_i)] = cur_end;
    }
}

// Perfect tandem-repeat tracts of the chunk reference: period 1-6, three copies and at least
// kStrTractMinLen bp, at most kStrTractMaxLen bp (longer arrays are satellites, left to sdust).
// sdust's 20 bp window misses dinucleotide and longer repeats -- a pure TG tract scores about 4.2
// against its threshold of 5 -- yet the aligner may put an indel anywhere along such a tract. A
// noisy region ending inside one leaves the reads that put the indel at the other end with none
// in the MSA: at chr20:11,837,334 13 of 28 reads carrying +TG put it outside the 21 bp region
// and the cluster consensus was decided by read order.
static std::vector<Interval> tandem_repeat_tracts(const PhasingChunk& chunk) {
    std::vector<Interval> out;
    const std::string& s = chunk.ref_seq;
    const int n = static_cast<int>(s.size());
    const auto at = [&s](int i) { return std::toupper(static_cast<unsigned char>(s[static_cast<size_t>(i)])); };
    for (int period = 1; period <= kStrTractMaxPeriod; ++period) {
        int i = period;
        while (i < n) {
            if (at(i) == 'N' || at(i) != at(i - period)) {
                ++i;
                continue;
            }
            int j = i;
            while (j < n && at(j) != 'N' && at(j) == at(j - period)) ++j;
            const int len = j - i + period;  // the tract is [i - period, j)
            // One base of padding each side: an insertion at either end of the tract (keyed at the base
            // after it) is in the tract.
            if (len >= kStrTractMinLen && len >= 3 * period && len <= kStrTractMaxLen)
                out.push_back(Interval{chunk.ref_beg + i - period - 1, chunk.ref_beg + j, len});
            i = j;
        }
    }
    return out;
}

// Low-complexity sequence noisy regions extend into: sdust's intervals and, with plegma, the tandem-repeat tracts
// (abPOA keeps longcallD's regions: its outputs are the port's parity check).
std::vector<Interval> noisy_extension_intervals(const PhasingChunk& chunk) {
    std::vector<Interval> out = chunk.low_complexity_regions;
    if (!use_plegma_msa()) return out;
    const std::vector<Interval> tracts = tandem_repeat_tracts(chunk);
    out.insert(out.end(), tracts.begin(), tracts.end());
    return out;
}

/**
 * @brief Merges per-read noisy intervals into unified chunk-level regions and extends low-comp windows.
 *
 * Pre-process noisy regions: merge read-level noisy intervals, extend
 * It sweeps across individually flagged noisy segments (dense mismatches),
 * merges overlapping intervals from multiple reads via `cgranges_t` merges,
 * extends them to encompass neighboring low-complexity sequence (like homopolymers),
 * and finally drops regions that lack sufficient read support (min_alt_depth or min_af)
 * to be deemed true structural variations rather than isolated noise bursts.
 *
 * into low-complexity spans, and filter by depth/ratio thresholds.
 * collection, immediately before `classify_cand_vars`.
 *
 * @param chunk Collection of read records and current candidates in this block.
 * @param opts User-defined filtering thresholds.
 */
void pre_process_noisy_regs_pgphase(PhasingChunk& chunk, const Options& opts) {
    if (chunk.noisy_regions.empty()) return;

    CrangesOwner noisy_own;
    noisy_own.adopt(intervals_to_cr(chunk.noisy_regions));
    if (noisy_own.cr == nullptr) {
        chunk.noisy_regions.clear();
        return;
    }
    CrangesOwner low_own;
    low_own.adopt(intervals_to_cr(noisy_extension_intervals(chunk)));

    cr_index(noisy_own.cr);
    if (low_own.cr != nullptr) cr_index(low_own.cr);

    const int merge_dis = opts.noisy_reg_merge_dis;
    const int msv = opts.min_sv_len;

    // Extend noisy regions into adjacent low-complexity intervals, then re-merge.
    cgranges_t* noisy = noisy_own.release();
    noisy = cr_extend_noisy_regs_with_low_comp(noisy, low_own.cr, merge_dis, msv);
    noisy = cr_merge(noisy, -1, merge_dis, msv);
    noisy_own.adopt(noisy);

    cgranges_t* noisy_regs = noisy_own.cr;
    int64_t* ovlp_b = nullptr;
    int64_t max_b = 0;
    std::vector<uint8_t> skip_noisy_reg(static_cast<size_t>(noisy_regs->n_r), 0);
    std::vector<int> noisy_reg_to_total(static_cast<size_t>(noisy_regs->n_r), 0);
    std::vector<int> noisy_reg_to_noisy(static_cast<size_t>(noisy_regs->n_r), 0);

    const std::vector<int>& read_order = chunk.ordered_read_ids;
    const auto visit_read = [&](const ReadRecord& read) {
        if (read.is_skipped) return;
        const int64_t beg = read.beg;
        const int64_t end = read.end;
        const int64_t ovlp_n = cr_overlap(
            noisy_regs, "cr", static_cast<int32_t>(beg - 1), static_cast<int32_t>(end), &ovlp_b, &max_b);
        CrangesOwner read_noisy_own;
        read_noisy_own.adopt(build_read_noisy_cr(read));
        int64_t* read_ovlp_b = nullptr;
        int64_t read_max_b = 0;
        for (int64_t ovlp_i = 0; ovlp_i < ovlp_n; ++ovlp_i) {
            const int r_i = static_cast<int>(ovlp_b[ovlp_i]);
            noisy_reg_to_total[static_cast<size_t>(r_i)]++;
            const int32_t noisy_reg_start = cr_start(noisy_regs, r_i) + 1;
            const int32_t noisy_reg_end = cr_end(noisy_regs, r_i);
            int64_t noisy_digar_ovlp_n = 0;
            if (read_noisy_own.cr != nullptr) {
                noisy_digar_ovlp_n =
                    cr_overlap(read_noisy_own.cr, "cr", noisy_reg_start - 1, noisy_reg_end,
                               &read_ovlp_b, &read_max_b);
            }
            if (noisy_digar_ovlp_n > 0) {
                noisy_reg_to_noisy[static_cast<size_t>(r_i)]++;
            }
        }
        std::free(read_ovlp_b);
    };
    if (read_order.empty() && !chunk.reads.empty()) {
        for (const ReadRecord& r : chunk.reads) visit_read(r);
    } else {
        for (int ord : read_order) {
            if (ord < 0 || static_cast<size_t>(ord) >= chunk.reads.size()) continue;
            visit_read(chunk.reads[static_cast<size_t>(ord)]);
        }
    }
    std::free(ovlp_b);

    const int min_noisy_reg_reads = opts.min_alt_depth;
    const float min_noisy_reg_ratio = static_cast<float>(opts.min_af);
    int n_skipped = 0;
    for (int i = 0; i < noisy_regs->n_r; ++i) {
        const int n_noisy = noisy_reg_to_noisy[static_cast<size_t>(i)];
        const int n_total = noisy_reg_to_total[static_cast<size_t>(i)];
        // Skip if too few noisy reads
        // or if noisy/total ratio is below threshold.
        if (n_noisy < min_noisy_reg_reads ||
            (n_total > 0 &&
             static_cast<float>(n_noisy) / static_cast<float>(n_total) < min_noisy_reg_ratio)) {
            skip_noisy_reg[static_cast<size_t>(i)] = 1;
            ++n_skipped;
        }
    }

    if (n_skipped > 0) {
        cgranges_t* new_noisy_regs = cr_init();
        for (int i = 0; i < noisy_regs->n_r; ++i) {
            if (skip_noisy_reg[static_cast<size_t>(i)]) continue;
            cr_add(new_noisy_regs, "cr", cr_start(noisy_regs, i), cr_end(noisy_regs, i),
                   cr_label(noisy_regs, i));
        }
        cr_index(new_noisy_regs);
        noisy_own.adopt(new_noisy_regs);
    }

    intervals_from_cr(noisy_own.cr, chunk.noisy_regions);
}

/**
 * Widen noisy intervals using classified candidates, then re-merge.
 * @param chunk Chunk whose `noisy_regions` are rewritten.
 * @param cand Classified candidates (categories drive flanking).
 */
void post_process_noisy_regs_pgphase(PhasingChunk& chunk, const CandidateTable& cand) {
    if (chunk.noisy_regions.empty()) return;
    CrangesOwner noisy_own;
    noisy_own.adopt(intervals_to_cr(chunk.noisy_regions));
    if (noisy_own.cr == nullptr) {
        chunk.noisy_regions.clear();
        return;
    }
    cr_index(noisy_own.cr);
    cgranges_t* chunk_noisy_regs = noisy_own.release();
    const int n_noisy_regs = chunk_noisy_regs->n_r;
    std::vector<int32_t> noisy_reg_start(static_cast<size_t>(n_noisy_regs));
    std::vector<int32_t> noisy_reg_end(static_cast<size_t>(n_noisy_regs));
    std::vector<VariantCategory> cand_cate;
    cand_cate.reserve(cand.size());
    for (const auto& cv : cand) cand_cate.push_back(cv.counts.category);
    collect_noisy_reg_start_end_pgphase(
        chunk_noisy_regs, cand, cand_cate, noisy_reg_start, noisy_reg_end);

    cgranges_t* noisy_regs = cr_init();
    for (int reg_i = 0; reg_i < n_noisy_regs; ++reg_i) {
        cr_add(noisy_regs, "cr",
               noisy_reg_start[static_cast<size_t>(reg_i)],
               noisy_reg_end[static_cast<size_t>(reg_i)],
               cr_label(chunk_noisy_regs, reg_i));
    }
    cr_index(noisy_regs);
    cr_destroy(chunk_noisy_regs);
    cgranges_t* merged = cr_merge(noisy_regs, 0, -1, -1);
    noisy_own.adopt(merged);
    intervals_from_cr_noisy_post_merge(noisy_own.cr, chunk.noisy_regions);
}

/**
 * @brief Traverses surviving variants and drops those entirely swallowed by a noisy block.
 *
 * Final masking pass: set NonVariant for candidates fully inside noisy spans.
 * Like classify compaction, skips sites whose bitmask intersects
 * Skips candidates already in kLongcalldNotCandVarCate.
 * Once noisy intervals are merged and widened via `post_process_noisy_regs_pgphase()`, this
 * marks surviving germline candidates fully contained in noisy spans as `NON_VAR`.
 *
 * This ensures that a single structural variance window gets reported cleanly
 * instead of hundreds of raw SNVs generated through misalignment.
 *
 * @param chunk Chunk tracking the current region's evaluated candidates.
 * @note No-op when `chunk.noisy_regions` is empty.
 */
// Noisy regions longer than the MSA takes (max_noisy_reg_len), split where the reads agree with the reference rather
// than skipped whole (plegma; abPOA keeps longcallD's regions). A cut point lies kNoisySplitFlank or more from every
// candidate the flanks count (category_skipped_for_noisy_flank) and from low-complexity sequence and tandem tracts;
// each piece is cut at the last such point within the cap. A region a few bases over the cap lost every call in it:
// chr20:50,273,460-50,283,511, 10,051 bp, its 686 bp deletion and the SNPs beside it, which a region of 9.7 kb from
// another alignment of the same reads resolved. A stretch with no cut point within the cap stays whole (and skipped).
// Shorter pieces are spanned by more reads, too. A piece of a collapsed repeat (a satellite) runs on its confidently
// placed reads or not at all (collect_noisy_vars1): with all their reads, chr20's centromere (59 regions up to 165 kb,
// pieces of up to 478 reads) took 30 minutes.
void split_long_noisy_regions(PhasingChunk& chunk, const Options& opts) {
    if (!use_plegma_msa() || opts.max_noisy_reg_len <= 0) return;
    const hts_pos_t cap = opts.max_noisy_reg_len;
    bool any = false;
    for (const Interval& iv : chunk.noisy_regions) any = any || iv.end - iv.beg + 1 > cap;
    if (!any) return;
    // Blocked: candidate spans and low-complexity intervals, padded by kNoisySplitFlank; merged.
    std::vector<std::pair<hts_pos_t, hts_pos_t>> blocked;
    for (const CandidateVariant& c : chunk.candidates) {
        if (category_skipped_for_noisy_flank(c.counts.category)) continue;
        const hts_pos_t e = c.key.pos + std::max(1, c.key.ref_len) - 1;
        blocked.emplace_back(c.key.pos - kNoisySplitFlank, e + kNoisySplitFlank);
    }
    for (const Interval& iv : noisy_extension_intervals(chunk))
        blocked.emplace_back(iv.beg - kNoisySplitFlank, iv.end + kNoisySplitFlank);
    std::sort(blocked.begin(), blocked.end());
    std::vector<std::pair<hts_pos_t, hts_pos_t>> merged;
    for (const auto& b : blocked) {
        if (!merged.empty() && b.first <= merged.back().second + 1)
            merged.back().second = std::max(merged.back().second, b.second);
        else
            merged.push_back(b);
    }
    // The last free position at or before p (and after lo), or -1.
    const auto last_free = [&merged](hts_pos_t lo, hts_pos_t p) -> hts_pos_t {
        auto it = std::upper_bound(merged.begin(), merged.end(), std::make_pair(p, std::numeric_limits<hts_pos_t>::max()));
        while (p > lo) {
            if (it == merged.begin() || std::prev(it)->second < p) return p;  // p lies in no blocked interval
            --it;
            p = it->first - 1;  // left of this blocked interval
        }
        return -1;
    };
    const auto count_in = [&chunk](hts_pos_t a, hts_pos_t b) {
        int n = 0;
        for (const CandidateVariant& c : chunk.candidates)
            n += !category_skipped_for_noisy_flank(c.counts.category) && c.key.pos >= a && c.key.pos <= b;
        return n;
    };
    std::vector<Interval> out;
    for (const Interval& iv : chunk.noisy_regions) {
        if (iv.end - iv.beg + 1 <= cap) {
            out.push_back(iv);
            continue;
        }
        hts_pos_t seg = iv.beg;
        int pieces = 0;
        while (iv.end - seg + 1 > cap) {
            const hts_pos_t cut = last_free(seg + kNoisySplitFlank, seg + cap - 1);
            if (cut < 0) break;
            out.push_back(Interval{seg, cut, count_in(seg, cut), true});
            ++pieces;
            seg = cut + 1;
        }
        out.push_back(Interval{seg, iv.end, count_in(seg, iv.end), pieces > 0});
        if (pgphase_trace::on())
            pgphase_trace::emit("msa", "region", iv.beg, iv.end, "-", "split",
                                "pieces=" + std::to_string(pieces + 1) + " last=" + std::to_string(seg) + "-" +
                                    std::to_string(iv.end));
    }
    std::sort(out.begin(), out.end(), [](const Interval& a, const Interval& b) { return a.beg < b.beg; });
    chunk.noisy_regions = std::move(out);
}

void apply_noisy_containment_filter(PhasingChunk& chunk, int max_suppressing_len) {
    if (chunk.noisy_regions.empty()) return;
    // A span longer than the MSA's limit is never resolved, so suppressing its
    // candidates would leave the locus with no calls at all. Union phasing
    // keeps them: windows widened to nearby candidates and merged can chain
    // across a whole SV-rich stretch -- chr18:20,779,364-20,877,543, 98 kb,
    // held dozens of clean heterozygous SNPs that separate the parents and
    // that LongPhase phases, and the alignment solve called none of them.
    // Classification already masked candidates overlapping the per-read
    // noisy windows (classify_cand_vars_pgphase, pass 2), before those were
    // merged into the spans measured here; a clean call masked only for that
    // keeps its classification in candvarcate_initial and is restored.
    std::vector<Interval> suppressing, unresolved;
    for (const Interval& r : chunk.noisy_regions) {
        if (max_suppressing_len <= 0 || r.end - r.beg + 1 <= max_suppressing_len) suppressing.push_back(r);
        else unresolved.push_back(r);
    }
    if (!unresolved.empty()) {
        CrangesOwner keep;
        keep.adopt(intervals_to_cr_noisy_post_merge(unresolved));
        if (keep.cr != nullptr && keep.cr->n_r > 0) {
            cr_index(keep.cr);
            int64_t* kb = nullptr;
            int64_t km = 0;
            for (auto& cv : chunk.candidates) {
                const VariantCategory init = cv.counts.candvarcate_initial;
                if (cv.counts.category != VariantCategory::NonVariant) continue;
                if (init != VariantCategory::CleanHetSnp) continue;
                const VariantKey& k = cv.key;
                const int32_t beg = static_cast<int32_t>(k.pos - 1);
                const int32_t end = static_cast<int32_t>(k.type == VariantType::Insertion ? k.pos : k.pos + k.ref_len);
                if (cr_is_contained(keep.cr, "cr", beg, end, &kb, &km) > 0) {
                    cv.counts.category = init;
                    cv.lcd_var_i_to_cate = category_to_flag(init);
                    cv.unresolved_noisy = true;
                }
            }
            std::free(kb);
        }
    }
    if (max_suppressing_len > 0) {
        // A repeat het indel is kept in any noisy region, short ones included.
        // It is never injected, but union phasing realigns reads to its two
        // exact alleles (build_allele_window_sites) where nothing else covers
        // the locus. A short region's MSA separates haplotypes by reads already
        // labelled, so inside a gap it calls nothing and the indel vanished:
        // chr12:108,955,590, 52 REF / 23 ALT, the het HiPhase bridges two of
        // our blocks with.
        CrangesOwner all_spans;
        all_spans.adopt(intervals_to_cr_noisy_post_merge(chunk.noisy_regions));
        if (all_spans.cr != nullptr && all_spans.cr->n_r > 0) {
            cr_index(all_spans.cr);
            int64_t* kb = nullptr;
            int64_t km = 0;
            for (auto& cv : chunk.candidates) {
                if (cv.counts.category != VariantCategory::NonVariant ||
                    cv.counts.candvarcate_initial != VariantCategory::RepeatHetIndel) continue;
                const VariantKey& k = cv.key;
                const int32_t beg = static_cast<int32_t>(k.pos - 1);
                const int32_t end = static_cast<int32_t>(k.type == VariantType::Insertion ? k.pos : k.pos + k.ref_len);
                if (cr_is_contained(all_spans.cr, "cr", beg, end, &kb, &km) > 0) {
                    cv.counts.category = VariantCategory::RepeatHetIndel;
                    cv.lcd_var_i_to_cate = category_to_flag(VariantCategory::RepeatHetIndel);
                }
            }
            std::free(kb);
        }
    }
    if (suppressing.empty()) return;
    CrangesOwner own;
    own.adopt(intervals_to_cr_noisy_post_merge(suppressing));
    if (own.cr == nullptr || own.cr->n_r == 0) return;
    cr_index(own.cr);
    int64_t* b = nullptr;
    int64_t m = 0;
    for (auto& cv : chunk.candidates) {
        if ((cv.lcd_var_i_to_cate & kLongcalldNotCandVarCate) != 0) continue;
        if (max_suppressing_len > 0 && cv.counts.category == VariantCategory::RepeatHetIndel) continue;
        const VariantKey& k = cv.key;
        int64_t n = 0;
        if (k.type == VariantType::Insertion) {
            n = cr_is_contained(
                own.cr, "cr", static_cast<int32_t>(k.pos - 1), static_cast<int32_t>(k.pos), &b, &m);
        } else {
            n = cr_is_contained(own.cr, "cr",
                                static_cast<int32_t>(k.pos - 1),
                                static_cast<int32_t>(k.pos + k.ref_len),
                                &b, &m);
        }
        if (n > 0) {
            cv.counts.category = VariantCategory::NonVariant;
            cv.lcd_var_i_to_cate = category_to_flag(VariantCategory::NonVariant);
        }
    }
    std::free(b);
}

/**
 * @brief Fills `chunk.low_complexity_regions` via sdust on the loaded reference slice.
 * @param chunk Must have `ref_seq` and region/ref bounds aligned with the slice.
 */
void populate_low_complexity_intervals(PhasingChunk& chunk) {
    chunk.low_complexity_regions.clear();
    if (chunk.ref_seq.empty()) return;

    const hts_pos_t beg = std::max(chunk.region.beg, chunk.ref_beg);
    const hts_pos_t end = std::min(chunk.region.end, chunk.ref_end);
    if (beg > end) return;

    const size_t offset = static_cast<size_t>(beg - chunk.ref_beg);
    const int len = static_cast<int>(end - beg + 1);
    const std::string slice = chunk.ref_seq.substr(offset, static_cast<size_t>(len));
    chunk.low_complexity_regions = find_low_complexity_intervals(slice, beg);
}

/**
 * @brief Sets `ordered_read_ids` and unions per-read `noisy_regions` into chunk-level intervals.
 * @param chunk Chunk with `reads` populated.
 */
void populate_chunk_read_indexes(PhasingChunk& chunk) {
    chunk.ordered_read_ids.clear();
    chunk.noisy_regions.clear();
    chunk.ordered_read_ids.reserve(chunk.reads.size());
    for (size_t read_i = 0; read_i < chunk.reads.size(); ++read_i) {
        const ReadRecord& read = chunk.reads[read_i];
        chunk.ordered_read_ids.push_back(static_cast<int>(read_i));
        if (read.is_skipped) continue;
        for (const Interval& iv : read.noisy_regions) {
            if (iv.end < chunk.region.beg || iv.beg > chunk.region.end) continue;
            chunk.noisy_regions.push_back(iv);
        }
    }
    // Sort reads: pos asc, end desc, NM asc, qname asc.
    // Use std::sort (non-stable) to mirror qsort tie behavior as closely as possible.
    std::sort(chunk.ordered_read_ids.begin(), chunk.ordered_read_ids.end(),
              [&](int lhs_i, int rhs_i) {
                  const ReadRecord& lhs = chunk.reads[static_cast<size_t>(lhs_i)];
                  const ReadRecord& rhs = chunk.reads[static_cast<size_t>(rhs_i)];
                  if (lhs.beg != rhs.beg) return lhs.beg < rhs.beg;
                  if (lhs.end != rhs.end) return lhs.end > rhs.end;
                  if (lhs.nm != rhs.nm) return lhs.nm < rhs.nm;
                  return lhs.qname < rhs.qname;
              });
}

/**
 * Print noisy-region summary at --verbose 2.
 */
void dump_all_noisy_regions(const PhasingChunk& chunk, const Options& opts, const bam_hdr_t* header) {
    if (opts.verbose < 2) return;
    const char* chrom = "unknown";
    if (header != nullptr && chunk.region.tid >= 0 && chunk.region.tid < header->n_targets) {
        chrom = header->target_name[chunk.region.tid];
    }
    std::fprintf(stderr, "AllNoisyRegions: %s:%" PRId64 "-%" PRId64 " (%zu)\n",
                 chrom,
                 static_cast<int64_t>(chunk.region.beg),
                 static_cast<int64_t>(chunk.region.end),
                 chunk.noisy_regions.size());
    for (const Interval& region : chunk.noisy_regions) {
        std::fprintf(stderr, "NoisyRegion: %s:%" PRId64 "-%" PRId64 " (%d)\n",
                     chrom,
                     static_cast<int64_t>(region.beg),
                     static_cast<int64_t>(region.end),
                     region.label);
    }
}

} // namespace pgphase_collect
