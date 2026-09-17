/**
 * @file collect_pipeline.cpp
 * @brief Region chunking, parallel candidate collection, streaming writers, and CLI for collect-bam-variation.
 *
 * @details Coordinates are 1-based inclusive (BED is converted in `load_bed_regions`). Chunks are batched by
 * `reg_chunk_i` in `run_collect_bam_variation` so TSV/VCF can be streamed without holding
 * all candidates in memory.
 */

#include "collect_pipeline.hpp"

#include "arg_parse.hpp"
#include "bam_digar.hpp"
#include "collect_bam_output.hpp"
#include "collect_output.hpp"
#include "collect_phase.hpp"
#include "collect_phase_pgbam.hpp"
#include "collect_var.hpp"
#include "gap_recovery.hpp"
#include "gap_evidence.hpp"

#include <algorithm>
#include <atomic>
#include <chrono>
#include <cctype>
#include <cstdint>
#include <cstring>
#include <fstream>
#include <filesystem>
#include <getopt.h>
#include <iostream>
#include <iterator>
#include <map>
#include <memory>
#include <mutex>
#include <numeric>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <thread>
#include <unordered_map>
#include <vector>

#include <htslib/sam.h>

namespace pgphase_collect {

using BamAuthorityIntervals = std::unordered_map<int, std::vector<Interval>>;

static std::string authority_bed_contig(std::string contig) {
    constexpr const char* kMaternalSuffix = "_MATERNAL";
    constexpr const char* kPaternalSuffix = "_PATERNAL";
    for (const char* suffix : {kMaternalSuffix, kPaternalSuffix}) {
        const size_t n = std::char_traits<char>::length(suffix);
        if (contig.size() >= n &&
            contig.compare(contig.size() - n, n, suffix) == 0) {
            contig.resize(contig.size() - n);
            break;
        }
    }
    return contig;
}

static int authority_bed_tid(const std::string& bed_contig,
                             const bam_hdr_t* header) {
    const std::string normalized = authority_bed_contig(bed_contig);
    int match = -1;
    for (int tid = 0; tid < header->n_targets; ++tid) {
        const std::string target = header->target_name[tid];
        const size_t hash = target.rfind('#');
        const std::string suffix = hash == std::string::npos
                                       ? target
                                       : target.substr(hash + 1);
        if (target != normalized && suffix != normalized) continue;
        if (match >= 0) return -1;
        match = tid;
    }
    return match;
}

static BamAuthorityIntervals load_bam_authority_intervals(
        const std::string& path,
        const bam_hdr_t* header) {
    std::ifstream in(path);
    if (!in) throw std::runtime_error("failed to open BAM-authoritative BED: " + path);

    BamAuthorityIntervals intervals;
    std::string line;
    while (std::getline(in, line)) {
        if (line.empty() || line[0] == '#' || line.rfind("track", 0) == 0)
            continue;
        std::istringstream fields(line);
        std::string contig;
        hts_pos_t bed_beg = 0;
        hts_pos_t bed_end = 0;
        if (!(fields >> contig >> bed_beg >> bed_end) ||
            bed_beg < 0 || bed_end <= bed_beg) {
            throw std::runtime_error("invalid interval in BAM-authoritative BED: " + line);
        }
        const int tid = authority_bed_tid(contig, header);
        if (tid < 0) continue;
        intervals[tid].push_back(Interval{bed_beg + 1, bed_end});
    }

    for (auto& [tid, rows] : intervals) {
        (void)tid;
        std::sort(rows.begin(), rows.end(), [](const Interval& a, const Interval& b) {
            return a.beg < b.beg || (a.beg == b.beg && a.end < b.end);
        });
        std::vector<Interval> merged;
        for (const Interval& row : rows) {
            if (!merged.empty() && row.beg <= merged.back().end + 1) {
                merged.back().end = std::max(merged.back().end, row.end);
            } else {
                merged.push_back(row);
            }
        }
        rows = std::move(merged);
    }
    return intervals;
}

static bool is_bam_authoritative_position(
        const BamAuthorityIntervals* intervals,
        int tid,
        hts_pos_t pos) {
    if (intervals == nullptr) return false;
    const auto found = intervals->find(tid);
    if (found == intervals->end()) return false;
    const std::vector<Interval>& rows = found->second;
    const auto it = std::upper_bound(
        rows.begin(), rows.end(), pos,
        [](hts_pos_t value, const Interval& row) { return value < row.beg; });
    if (it == rows.begin()) return false;
    const Interval& row = *std::prev(it);
    return pos >= row.beg && pos <= row.end;
}

// ════════════════════════════════════════════════════════════════════════════
// Region chunking
// ════════════════════════════════════════════════════════════════════════════

/**
 * @brief Parses a genomic region string into a `RegionFilter`.
 *
 * Accepts `chrom`, `chrom:pos`, or `chrom:start-end` with optional comma thousands separators.
 * An empty string yields `enabled == false`.
 *
 * @param region Region literal from `-r` or positional arguments.
 * @return Parsed filter; throws if coordinates are invalid.
 */
RegionFilter parse_region(const std::string& region) {
    RegionFilter filter;
    if (region.empty()) return filter;

    auto strip_commas = [](std::string value) {
        value.erase(std::remove(value.begin(), value.end(), ','), value.end());
        return value;
    };

    filter.enabled = true;
    const size_t colon = region.find(':');
    if (colon == std::string::npos) {
        filter.chrom = region;
        return filter;
    }

    filter.chrom = region.substr(0, colon);
    const std::string range = region.substr(colon + 1);
    const size_t dash = range.find('-');
    const std::string beg = strip_commas(dash == std::string::npos ? range : range.substr(0, dash));
    const std::string end = strip_commas(dash == std::string::npos ? "" : range.substr(dash + 1));
    if (!beg.empty()) filter.beg = std::stoll(beg);
    if (!end.empty()) filter.end = std::stoll(end);
    if (filter.chrom.empty() || filter.beg < 1 || (filter.end >= 0 && filter.end < filter.beg)) {
        throw std::runtime_error("invalid region: " + region);
    }
    return filter;
}

/**
 * @brief Loads 3-column BED intervals as `RegionFilter` entries.
 *
 * Skips blank and `#` lines. Converts BED 0-based half-open `[bed_beg, bed_end)` to
 * 1-based inclusive `[bed_beg+1, bed_end]` on the reference.
 *
 * @param path Path to `--region-file`.
 * @return List of enabled filters; throws on malformed lines.
 */
std::vector<RegionFilter> load_bed_regions(const std::string& path) {
    std::ifstream in(path);
    if (!in) throw std::runtime_error("failed to open region file: " + path);
    std::vector<RegionFilter> regions;
    std::string line;
    while (std::getline(in, line)) {
        if (line.empty() || line[0] == '#') continue;
        std::istringstream fields(line);
        std::string chrom;
        hts_pos_t bed_beg = 0;
        hts_pos_t bed_end = 0;
        if (!(fields >> chrom >> bed_beg >> bed_end)) {
            throw std::runtime_error("invalid BED line in region file: " + line);
        }
        if (bed_beg < 0 || bed_end <= bed_beg) {
            throw std::runtime_error("invalid BED interval in region file: " + line);
        }
        regions.push_back(RegionFilter{true, chrom, bed_beg + 1, bed_end});
    }
    return regions;
}

/**
 * @brief Resolves a sequence name to BAM target id (`@SQ`).
 *
 * @param header BAM/CRAM header.
 * @param chrom Contig name (e.g. `chr1`).
 * @return Target index, or -1 if not found.
 */
static int tid_for_name(const bam_hdr_t* header, const std::string& chrom) {
    for (int tid = 0; tid < header->n_targets; ++tid) {
        if (chrom == header->target_name[tid]) return tid;
    }
    return -1;
}

/**
 * @brief Partitions a 1-based inclusive interval into fixed-width chunks.
 *
 * @param tid Contig index for every emitted `RegionChunk`.
 * @param beg First reference position (1-based inclusive).
 * @param end Last reference position (1-based inclusive).
 * @param chunk_size Maximum chunk width in bp.
 * @return Chunks with `beg`/`end` set; `chunk_id` and neighbors are filled later.
 */
static std::vector<RegionChunk> split_region(int tid,
                                             hts_pos_t beg,
                                             hts_pos_t end,
                                             hts_pos_t chunk_size) {
    std::vector<RegionChunk> chunks;
    for (hts_pos_t chunk_beg = beg; chunk_beg <= end; chunk_beg += chunk_size) {
        const hts_pos_t chunk_end = std::min(end, chunk_beg + chunk_size - 1);
        RegionChunk chunk;
        chunk.tid = tid;
        chunk.beg = chunk_beg;
        chunk.end = chunk_end;
        chunks.push_back(chunk);
    }
    return chunks;
}

/**
 * @brief Sorts chunks and fills ids used for batching and boundary overlap logic.
 *
 * Assigns `chunk_id`, per-contig `reg_chunk_i` and `reg_i`, and `prev_*` / `next_*` neighbor
 * fields when the adjacent chunk is on the same contig.
 *
 * @param chunks In/out list (typically from `split_region` / `add_filter_chunks`).
 */
static void annotate_chunk_neighbors(std::vector<RegionChunk>& chunks) {
    std::sort(chunks.begin(), chunks.end(), [](const RegionChunk& lhs, const RegionChunk& rhs) {
        if (lhs.tid != rhs.tid) return lhs.tid < rhs.tid;
        if (lhs.beg != rhs.beg) return lhs.beg < rhs.beg;
        return lhs.end < rhs.end;
    });

    int reg_chunk_i = -1;
    int reg_i = 0;
    int last_tid = -1;
    for (size_t i = 0; i < chunks.size(); ++i) {
        RegionChunk& chunk = chunks[i];
        chunk.chunk_id = static_cast<int>(i);
        if (i == 0 || chunk.tid != last_tid) {
            ++reg_chunk_i;
            reg_i = 0;
        }
        chunk.reg_chunk_i = reg_chunk_i;
        chunk.reg_i = reg_i++;
        last_tid = chunk.tid;
    }
    for (size_t i = 0; i < chunks.size(); ++i) {
        RegionChunk& chunk = chunks[i];
        if (i > 0 && chunks[i - 1].tid == chunk.tid) {
            chunk.prev_chunk_id = chunks[i - 1].chunk_id;
            chunk.prev_tid = chunks[i - 1].tid;
            chunk.prev_beg = chunks[i - 1].beg;
            chunk.prev_end = chunks[i - 1].end;
        }
        if (i + 1 < chunks.size() && chunks[i + 1].tid == chunk.tid) {
            chunk.next_chunk_id = chunks[i + 1].chunk_id;
            chunk.next_tid = chunks[i + 1].tid;
            chunk.next_beg = chunks[i + 1].beg;
            chunk.next_end = chunks[i + 1].end;
        }
    }
}

/**
 * @brief Appends `RegionChunk` tiles for one filter to `chunks`.
 *
 * Clips `filter.end` to the contig length when `end == -1` (full contig). Validates that the
 * contig exists in both BAM header and FASTA index.
 *
 * @param region Enabled filter with `chrom` and 1-based bounds.
 * @param header BAM header for contig lengths and name resolution.
 * @param fai Reference index for `faidx_has_seq`.
 * @param chunk_size Passed to `split_region`.
 * @param chunks Destination vector.
 */
static void add_filter_chunks(const RegionFilter& region,
                              const bam_hdr_t* header,
                              const faidx_t* fai,
                              hts_pos_t chunk_size,
                              std::vector<RegionChunk>& chunks) {
    if (!region.enabled) return;
    const int tid = tid_for_name(header, region.chrom);
    if (tid < 0)
        throw std::runtime_error("region contig is not present in BAM header: " + region.chrom);
    if (!faidx_has_seq(fai, region.chrom.c_str()))
        throw std::runtime_error("region contig is not present in FASTA: " + region.chrom);

    const hts_pos_t contig_end = static_cast<hts_pos_t>(header->target_len[tid]);
    const hts_pos_t end = region.end < 0 ? contig_end : std::min(region.end, contig_end);
    if (region.beg > end) return;
    auto region_chunks = split_region(tid, region.beg, end, chunk_size);
    chunks.insert(chunks.end(), region_chunks.begin(), region_chunks.end());
}

/**
 * @brief Builds the full list of region filters from CLI options.
 *
 * Concatenates `-r` entries, BED from `--region-file`, and autosome `chrN`/`N` whole-chromosome
 * filters when `--autosome` is set (requires contig in both BAM and FASTA).
 *
 * @param opts User options.
 * @param header Primary BAM header.
 * @param fai Reference FASTA index.
 * @return Ordered list of filters (may be empty if no `-r`/BED/autosome).
 */
static std::vector<RegionFilter> collect_region_filters(const Options& opts,
                                                        const bam_hdr_t* header,
                                                        const faidx_t* fai) {
    std::vector<RegionFilter> filters;
    for (const std::string& region : opts.regions) {
        filters.push_back(parse_region(region));
    }
    if (!opts.region_file.empty()) {
        auto bed_regions = load_bed_regions(opts.region_file);
        filters.insert(filters.end(), bed_regions.begin(), bed_regions.end());
    }
    if (opts.autosome) {
        for (int i = 1; i <= 22; ++i) {
            const std::string no_prefix = std::to_string(i);
            const std::string with_prefix = "chr" + no_prefix;
            if (tid_for_name(header, with_prefix) >= 0 && faidx_has_seq(fai, with_prefix.c_str())) {
                filters.push_back(RegionFilter{true, with_prefix, 1, -1});
            } else if (tid_for_name(header, no_prefix) >= 0 &&
                       faidx_has_seq(fai, no_prefix.c_str())) {
                filters.push_back(RegionFilter{true, no_prefix, 1, -1});
            }
        }
    }
    return filters;
}

/**
 * @brief Produces all `RegionChunk` tiles for the run.
 *
 * If any filter is present, only those intervals are tiled; otherwise every BAM contig that
 * exists in the FASTA index is split. Always ends with `annotate_chunk_neighbors`.
 *
 * @param opts Chunk size and region inputs.
 * @param header BAM header.
 * @param fai Reference index.
 * @return Sorted, annotated chunk list (may be empty if no valid contigs).
 */
std::vector<RegionChunk> build_region_chunks(const Options& opts,
                                             const bam_hdr_t* header,
                                             const faidx_t* fai) {
    const std::vector<RegionFilter> filters = collect_region_filters(opts, header, fai);
    std::vector<RegionChunk> chunks;

    if (!filters.empty()) {
        for (const RegionFilter& filter : filters) {
            add_filter_chunks(filter, header, fai, opts.chunk_size, chunks);
        }
        annotate_chunk_neighbors(chunks);
        return chunks;
    }

    for (int tid = 0; tid < header->n_targets; ++tid) {
        if (!faidx_has_seq(fai, header->target_name[tid])) continue;
        const hts_pos_t contig_end = static_cast<hts_pos_t>(header->target_len[tid]);
        if (contig_end <= 0) continue;
        auto contig_chunks = split_region(tid, 1, contig_end, opts.chunk_size);
        chunks.insert(chunks.end(), contig_chunks.begin(), contig_chunks.end());
    }

    annotate_chunk_neighbors(chunks);
    return chunks;
}

std::vector<RegionChunk> build_region_chunks(const Options& opts,
                                             const bam_hdr_t* header,
                                             const faidx_t* fai,
                                             const std::vector<RegionFilter>& filters) {
    std::vector<RegionChunk> chunks;
    if (!filters.empty()) {
        for (const RegionFilter& filter : filters) {
            add_filter_chunks(filter, header, fai, opts.chunk_size, chunks);
        }
        annotate_chunk_neighbors(chunks);
        return chunks;
    }
    for (int tid = 0; tid < header->n_targets; ++tid) {
        if (!faidx_has_seq(fai, header->target_name[tid])) continue;
        const hts_pos_t contig_end = static_cast<hts_pos_t>(header->target_len[tid]);
        if (contig_end <= 0) continue;
        auto contig_chunks = split_region(tid, 1, contig_end, opts.chunk_size);
        chunks.insert(chunks.end(), contig_chunks.begin(), contig_chunks.end());
    }
    annotate_chunk_neighbors(chunks);
    return chunks;
}

/**
 * @brief Opens primary alignment + reference and returns chunks from `build_region_chunks`.
 *
 * @param opts Must include `primary_bam_file()`, `ref_fasta`, and region-related fields.
 * @return Chunk list; throws if BAM header, index, or FAI cannot be opened.
 */
std::vector<RegionChunk> load_region_chunks(const Options& opts) {
    SamFile bam(opts.primary_bam_file(), 1, opts.ref_fasta);
    std::unique_ptr<bam_hdr_t, HeaderDeleter> header(sam_hdr_read(bam.get()));
    if (!header) throw std::runtime_error("failed to read BAM header");
    std::unique_ptr<hts_idx_t, IndexDeleter> index(
        sam_index_load(bam.get(), opts.primary_bam_file().c_str()));
    if (!index) throw std::runtime_error("region chunking requires an indexed BAM/CRAM");
    std::unique_ptr<faidx_t, FaiDeleter> fai(load_reference_index(opts.ref_fasta));
    return build_region_chunks(opts, header.get(), fai.get());
}

// ════════════════════════════════════════════════════════════════════════════
// Pipeline
// ════════════════════════════════════════════════════════════════════════════

// Initialize per-chunk bookkeeping after reads are loaded: identify which
// reads overlap adjacent chunks (for stitching), and zero-fill hap/PS arrays
// before k-means populates them.
static void initialize_chunk_overlap_state(PhasingChunk& chunk, size_t n_bams) {
    chunk.up_ovlp_read_i.assign(n_bams, {});
    chunk.down_ovlp_read_i.assign(n_bams, {});
    chunk.n_up_ovlp_skip_reads.assign(n_bams, 0);
    chunk.n_down_ovlp_skip_reads.assign(n_bams, 0);
    chunk.haps.assign(chunk.reads.size(), 0);
    chunk.phase_sets.assign(chunk.reads.size(), -1);

    for (size_t read_i = 0; read_i < chunk.reads.size(); ++read_i) {
        const ReadRecord& read = chunk.reads[read_i];
        if (read.input_index < 0 || static_cast<size_t>(read.input_index) >= n_bams) continue;
        if (read_overlaps_prev_region(chunk.region, read.tid, read.beg, read.end)) {
            chunk.up_ovlp_read_i[read.input_index].push_back(static_cast<int>(read_i));
        }
        if (read_overlaps_next_region(chunk.region, read.tid, read.beg, read.end)) {
            chunk.down_ovlp_read_i[read.input_index].push_back(static_cast<int>(read_i));
        }
    }
}

// Load reads from all input BAMs for one region chunk, identify overlap reads
// for stitching, and finalize the chunk (populate reference slice, compute
// quality stats, build digar-derived data structures).
static void load_and_prepare_chunk(PhasingChunk& chunk, const Options& opts, WorkerContext& context) {
    chunk.reads.clear();
    std::vector<OverlapSkipCounts> overlap_skip_counts(context.bams.size());
    for (size_t input_i = 0; input_i < context.bams.size(); ++input_i) {
        std::vector<ReadRecord> reads = load_read_records_for_chunk(
            opts,
            chunk.region,
            static_cast<int>(input_i),
            *context.bams[input_i],
            context.headers[input_i].get(),
            context.indexes[input_i].get(),
            context.ref,
            &overlap_skip_counts[input_i]);
        chunk.reads.reserve(chunk.reads.size() + reads.size());
        chunk.reads.insert(chunk.reads.end(),
                           std::make_move_iterator(reads.begin()),
                           std::make_move_iterator(reads.end()));
    }
    initialize_chunk_overlap_state(chunk, context.bams.size());
    for (size_t input_i = 0; input_i < overlap_skip_counts.size(); ++input_i) {
        chunk.n_up_ovlp_skip_reads[input_i] = overlap_skip_counts[input_i].upstream;
        chunk.n_down_ovlp_skip_reads[input_i] = overlap_skip_counts[input_i].downstream;
    }
    finalize_bam_chunk(chunk, context.ref, context.primary_header());
}

// Release bulky per-chunk intermediates after variant calling + k-means
// phasing but before stitching/output.  Frees digars, quality arrays,
// noisy-region interval trees, and low-complexity regions.  Keeps reads
// (for phased-BAM output), candidates, read profiles, haps/phase_sets,
// and overlap indices (all needed by stitching and VCF emission).
static void mid_free_chunk(PhasingChunk& chunk, const Options& opts) {
    // Per-read digar ops, base quals, and noisy sub-intervals.
    const bool keep_digars_for_refine = !opts.output_aln.empty() && opts.refine_aln;
    if (!keep_digars_for_refine) {
        for (ReadRecord& r : chunk.reads) {
            r.digars.clear();
            r.digars.shrink_to_fit();
            r.qual.clear();
            r.qual.shrink_to_fit();
            r.noisy_regions.clear();
            r.noisy_regions.shrink_to_fit();
        }
    }
    // Low-complexity interval list (used only during classification).
    chunk.low_complexity_regions.clear();
    chunk.low_complexity_regions.shrink_to_fit();
    // Noisy-read coverage/error interval trees and dedup marks.
    chunk.var_noisy_read_cov_cr.reset();
    chunk.var_noisy_read_err_cr.reset();
    chunk.var_noisy_read_marks.clear();
    chunk.var_noisy_read_marks.shrink_to_fit();
    chunk.var_noisy_read_mark_id = 0;
    // Chunk-level noisy region list.
    chunk.noisy_regions.clear();
    chunk.noisy_regions.shrink_to_fit();
    // Read-variant overlap interval tree (rebuilt if needed later).
    chunk.read_var_cr.reset();
}

// Process one genomic region: load reads from BAM, run variant calling
// (candidate collection, classification, k-means phasing), then free
// heavy intermediates.  Returns a PhasingChunk with candidates, read
// profiles, and hap assignments ready for stitching.
static PhasingChunk process_chunk(const RegionChunk& region,
                              const Options& opts,
                              WorkerContext& context) {
    PhasingChunk chunk;
    chunk.region = region;
    load_and_prepare_chunk(chunk, opts, context);
    collect_var_main(chunk, opts, context.primary_header());
    mid_free_chunk(chunk, opts);
    return chunk;
}

// Merge per-chunk candidate tables into a single sorted table, deduplicating
// by exact variant key (tid, pos, type, ref_len, alt).  Only exact-key dedup
// is applied here — fuzzy insertion collapse was already done within each
// chunk during candidate collection.
//
// When tiling overlap produces duplicate keys, the copy whose position falls
// inside its chunk's active region is preferred; ties keep stream order.
static CandidateTable merge_chunk_candidates(std::vector<PhasingChunk>& chunks) {
    struct TaggedRow {
        CandidateVariant v;
        int chunk_i = 0;
        int ord_in_chunk = 0;
    };
    std::vector<TaggedRow> rows;
    size_t reserve_n = 0;
    for (const PhasingChunk& chunk : chunks) reserve_n += chunk.candidates.size();
    rows.reserve(reserve_n);

    for (size_t ci = 0; ci < chunks.size(); ++ci) {
        CandidateTable& table = chunks[ci].candidates;
        for (size_t j = 0; j < table.size(); ++j) {
            rows.push_back(TaggedRow{std::move(table[j]), static_cast<int>(ci), static_cast<int>(j)});
        }
    }
    if (rows.empty()) return {};

    std::sort(rows.begin(), rows.end(), [](const TaggedRow& a, const TaggedRow& b) {
        const int c = exact_comp_var_site(&a.v.key, &b.v.key);
        if (c != 0) return c < 0;
        if (a.chunk_i != b.chunk_i) return a.chunk_i < b.chunk_i;
        return a.ord_in_chunk < b.ord_in_chunk;
    });

    CandidateTable merged;
    merged.reserve(rows.size());
    for (size_t i = 0; i < rows.size(); ++i) {
        if (!merged.empty() &&
            exact_comp_var_site(&merged.back().key, &rows[i].v.key) == 0) {
            const bool prev_pass = merged.back().lcd_make_variants_region_pass;
            const bool cur_pass = rows[i].v.lcd_make_variants_region_pass;
            if (!prev_pass && cur_pass) {
                merged.back() = std::move(rows[i].v);
            } else {
                merged.back().lcd_make_variants_region_pass = prev_pass || cur_pass;
            }
            continue;
        }
        merged.push_back(std::move(rows[i].v));
    }
    return merged;
}

// Parallel batch result: one PhasingChunk per chunk offset.
struct ChunkBatchResult {
    std::vector<PhasingChunk> chunks;
};

// Run process_chunk on chunks[batch_begin..batch_end) using a thread pool.
// Each worker opens its own BAM/FAI handles.  First exception is rethrown
// after all workers join.
static ChunkBatchResult collect_chunk_batch_parallel(const Options& opts,
                                                     const std::vector<RegionChunk>& chunks,
                                                     size_t batch_begin,
                                                     size_t batch_end) {
    if (batch_begin > batch_end || batch_end > chunks.size()) {
        throw std::runtime_error("invalid chunk batch range");
    }
    const size_t batch_size = batch_end - batch_begin;
    ChunkBatchResult result;
    result.chunks.resize(batch_size);
    if (batch_size == 0) return result;

    const size_t worker_count = std::min<size_t>(static_cast<size_t>(opts.threads), batch_size);
    std::atomic<size_t> next_offset{0};
    std::exception_ptr first_error;
    std::mutex error_mutex;
    std::vector<std::thread> workers;
    workers.reserve(worker_count);

    for (size_t worker_i = 0; worker_i < worker_count; ++worker_i) {
        workers.emplace_back([&, worker_i]() {
            (void)worker_i;
            try {
                WorkerContext context(opts);
                while (true) {
                    const size_t offset = next_offset.fetch_add(1);
                    if (offset >= batch_size) break;
                    result.chunks[offset] =
                        process_chunk(chunks[batch_begin + offset], opts, context);
                }
            } catch (...) {
                std::lock_guard<std::mutex> lock(error_mutex);
                if (!first_error) first_error = std::current_exception();
            }
        });
    }

    for (std::thread& worker : workers) worker.join();
    if (first_error) std::rethrow_exception(first_error);

    return result;
}

/**
 * @brief End-to-end collect-bam-variation driver with streaming output.
 *
 * Groups chunks by `reg_chunk_i`, processes each batch in parallel, merges candidates in memory
 * only within the batch, then appends TSV rows and optional VCF lines. Does not
 * hold the full genome candidate set in RAM.
 *
 * @param opts Output paths, reference, BAM list, and threading configuration.
 */
void run_collect_bam_variation(const Options& opts) {
    std::unique_ptr<PgbamSidecarData> pgbam_sidecar;
    if (!opts.pgbam_file.empty()) {
        pgbam_sidecar = std::make_unique<PgbamSidecarData>(load_pgbam_sidecar(opts.pgbam_file));
    }

    std::unique_ptr<faidx_t, FaiDeleter> fai(load_reference_index(opts.ref_fasta));

    const std::vector<RegionChunk> chunks = load_region_chunks(opts);
    SamFile bam(opts.primary_bam_file(), 1, opts.ref_fasta);
    std::unique_ptr<bam_hdr_t, HeaderDeleter> header(sam_hdr_read(bam.get()));
    if (!header) throw std::runtime_error("failed to read BAM header");
    ReferenceCache ref(fai.get());

    std::ofstream variant_out(opts.output_tsv);
    if (!variant_out) throw std::runtime_error("failed to open output: " + opts.output_tsv);
    write_variants_tsv_header(variant_out);

    std::ofstream vcf_out;
    if (!opts.output_vcf.empty()) {
        vcf_out.open(opts.output_vcf);
        if (!vcf_out) throw std::runtime_error("failed to open VCF output: " + opts.output_vcf);
        write_variants_vcf_header(vcf_out, opts, header.get());
    }
    std::ofstream phased_vcf_out;
    if (!opts.output_phased_vcf.empty()) {
        phased_vcf_out.open(opts.output_phased_vcf);
        if (!phased_vcf_out) {
            throw std::runtime_error("failed to open phased VCF output: " + opts.output_phased_vcf);
        }
        write_phased_variants_vcf_header(phased_vcf_out, opts, header.get());
    }

    std::unique_ptr<PhasedAlignmentWriter> phased_aln_writer;
    if (!opts.output_aln.empty()) {
        phased_aln_writer = std::make_unique<PhasedAlignmentWriter>(opts, header.get());
    }

    size_t n_variants = 0;
    size_t n_out_aln_reads = 0;
    size_t batch_begin = 0;
    while (batch_begin < chunks.size()) {
        size_t batch_end = batch_begin + 1;
        while (batch_end < chunks.size() &&
               chunks[batch_end].reg_chunk_i == chunks[batch_begin].reg_chunk_i) {
            ++batch_end;
        }

        ChunkBatchResult batch = collect_chunk_batch_parallel(
            opts, chunks, batch_begin, batch_end);
        stitch_chunk_haps(batch.chunks, &opts, pgbam_sidecar.get());
        CandidateTable variants = merge_chunk_candidates(batch.chunks);
        n_variants += variants.size();
        write_variants_tsv_records(variant_out, header.get(), ref, variants);
        if (!opts.output_vcf.empty()) {
            write_variants_vcf_records(vcf_out, opts, header.get(), ref, variants);
        }
        if (!opts.output_phased_vcf.empty()) {
            write_phased_variants_vcf_records(phased_vcf_out, opts, header.get(), ref, variants);
        }
        if (phased_aln_writer) {
            n_out_aln_reads += static_cast<size_t>(phased_aln_writer->write_chunks(batch.chunks));
        }

        batch_begin = batch_end;
    }

    std::cerr << "Processed " << chunks.size() << " region chunks with " << opts.threads
              << " worker thread(s)\n";
    std::cerr << "Collected " << n_variants << " candidate variant sites into "
              << opts.output_tsv << "\n";
    if (!opts.output_vcf.empty()) {
        std::cerr << "Wrote candidate VCF to " << opts.output_vcf << "\n";
    }
    if (!opts.output_phased_vcf.empty()) {
        std::cerr << "Wrote phased candidate VCF to " << opts.output_phased_vcf << "\n";
    }
    if (!opts.output_aln.empty()) {
        // Ensure output handle is closed/flushed before indexing.
        phased_aln_writer.reset();
        if (opts.output_aln_format == OutputAlignmentFormat::Cram) {
            std::cerr << "Output " << n_out_aln_reads << " reads to CRAM\n";
        } else if (opts.output_aln_format == OutputAlignmentFormat::Sam) {
            std::cerr << "Output " << n_out_aln_reads << " reads to SAM\n";
        } else {
            std::cerr << "Output " << n_out_aln_reads << " reads to BAM\n";
        }
        if (opts.refine_aln) {
            std::cerr << "Coordinate-sorting refined alignment (samtools sort -@"
                      << std::max(1, opts.threads) << ")...\n";
            coordinate_sort_refined_alignment_file_or_throw(opts);
        }
        if (opts.output_aln_format == OutputAlignmentFormat::Bam ||
            opts.output_aln_format == OutputAlignmentFormat::Cram) {
            if (sam_index_build(opts.output_aln.c_str(), 0) != 0) {
                throw std::runtime_error("failed to index output alignment: " + opts.output_aln);
            }
            if (opts.output_aln_format == OutputAlignmentFormat::Cram) {
                std::cerr << "Indexed output CRAM: " << opts.output_aln << ".crai\n";
            } else {
                std::cerr << "Indexed output BAM: " << opts.output_aln << ".bai\n";
            }
        }
    }
}

// ════════════════════════════════════════════════════════════════════════════
// Hybrid BAM+graph pipeline
// ════════════════════════════════════════════════════════════════════════════

} // namespace pgphase_collect

#include "hybrid_inject.hpp"
#include "graph_sites.hpp"
#include "graph_query.hpp"

namespace pgphase_collect {

/// Reference windows where this chunk's solve failed, either way it can fail.
///
/// Two failures, both observed rather than predicted from site spacing:
///
///   NOT PHASED    -- reads cover the window and no eligible het candidate there
///                    carries a phase set, so nothing in it was phased at all.
///   NOT CONNECTED -- it was phased, but into a different phase set from the
///                    block before it, so the solve could not link the two.
///
/// The second matters as much as the first: a window phased into its own block
/// leaves the same unusable result as one not phased, and it is the shape most of
/// these gaps actually have. Scattered unphased reads inside an otherwise phased
/// stretch are not a failure and are left alone.
static std::vector<std::pair<hts_pos_t, hts_pos_t>> collect_unphased_windows(
        const PhasingChunk& chunk, int min_reads, hts_pos_t min_bp) {
    constexpr hts_pos_t kBin = 1000;
    // A read's phase set comes only from an eligible heterozygous candidate
    // (update_read_phase_set, collect_phase.cpp): homopolymer indels,
    // kCandNoisyCandHom and unsupported MSA insertions are skipped, and a read
    // with no such candidate gets ps = -1. So a window the solve could not phase
    // is one holding no such candidate -- and it is measured that way rather than
    // by asking whether reads carry a phase set, because a read reaching in from
    // a flanking block carries one earned outside the window and would mask it.
    // The eligibility test mirrors update_read_phase_set exactly, including the
    // hap_to_cons_alle[1] / [2] indices; [0] is not a haplotype allele.
    std::set<hts_pos_t> phasing_bins;
    for (const auto& cand : chunk.candidates) {
        if (cand.phase_set < 0) continue;
        if (cand.is_homopolymer_indel ||
            cand.lcd_var_i_to_cate == kCandNoisyCandHom ||
            (!cand.msa_insertion_alts.empty() && !cand.gap_link_supported)) continue;
        if (cand.hap_to_cons_alle[1] == -1 || cand.hap_to_cons_alle[2] == -1 ||
            cand.hap_to_cons_alle[1] == cand.hap_to_cons_alle[2]) continue;
        phasing_bins.insert(cand.key.pos / kBin);
    }
    std::map<hts_pos_t, std::pair<int, int>> bins;  // bin -> {phased, unphased}
    const size_t n = std::min(chunk.reads.size(), chunk.haps.size());
    for (size_t i = 0; i < n; ++i) {
        const ReadRecord& read = chunk.reads[i];
        if (read.is_skipped || read.end <= read.beg) continue;
        for (hts_pos_t b = read.beg / kBin; b <= read.end / kBin; ++b) {
            auto& e = bins[b];
            if (phasing_bins.count(b) != 0) ++e.first; else ++e.second;
        }
    }
    std::vector<std::pair<hts_pos_t, hts_pos_t>> out;
    // NOT CONNECTED: consecutive phasing candidates in different phase sets. The
    // span between them is where the link failed, whatever its width -- a break
    // is a break -- so it is admitted without the width floor that the
    // not-phased runs carry.
    std::vector<std::pair<hts_pos_t, hts_pos_t>> breaks;
    {
        std::vector<std::pair<hts_pos_t, hts_pos_t>> ordered;  // pos -> phase set
        for (const auto& cand : chunk.candidates) {
            if (cand.phase_set < 0) continue;
            if (cand.is_homopolymer_indel ||
                cand.lcd_var_i_to_cate == kCandNoisyCandHom ||
                (!cand.msa_insertion_alts.empty() && !cand.gap_link_supported)) continue;
            if (cand.hap_to_cons_alle[1] == -1 || cand.hap_to_cons_alle[2] == -1 ||
                cand.hap_to_cons_alle[1] == cand.hap_to_cons_alle[2]) continue;
            ordered.emplace_back(cand.key.pos, cand.phase_set);
        }
        std::sort(ordered.begin(), ordered.end());
        for (size_t i = 1; i < ordered.size(); ++i) {
            if (ordered[i].second == ordered[i - 1].second) continue;
            breaks.emplace_back(ordered[i - 1].first, ordered[i].first + 1);
        }
    }
    hts_pos_t run_beg = -1;
    hts_pos_t prev = -2;
    for (const auto& [b, counts] : bins) {
        const bool dead = counts.first == 0 && counts.second >= min_reads;
        if (dead) {
            if (run_beg < 0 || b != prev + 1) {
                if (run_beg >= 0 && (prev + 1) * kBin - run_beg * kBin >= min_bp)
                    out.emplace_back(run_beg * kBin, (prev + 1) * kBin);
                run_beg = b;
            }
            prev = b;
        } else if (run_beg >= 0) {
            if ((prev + 1) * kBin - run_beg * kBin >= min_bp)
                out.emplace_back(run_beg * kBin, (prev + 1) * kBin);
            run_beg = -1;
        }
    }
    if (run_beg >= 0 && (prev + 1) * kBin - run_beg * kBin >= min_bp)
        out.emplace_back(run_beg * kBin, (prev + 1) * kBin);
    out.insert(out.end(), breaks.begin(), breaks.end());
    std::sort(out.begin(), out.end());
    // Merge overlaps so a candidate is not admitted twice and the report reads
    // as one failure per region.
    std::vector<std::pair<hts_pos_t, hts_pos_t>> merged;
    for (const auto& w : out) {
        if (!merged.empty() && w.first <= merged.back().second)
            merged.back().second = std::max(merged.back().second, w.second);
        else
            merged.push_back(w);
    }
    return merged;
}

/// Hybrid per-chunk processing: BAM classification → graph site injection →
/// BAM profile build → graph read injection → unified k-means.
static PhasingChunk process_chunk_hybrid(
        const RegionChunk& region,
        const Options& opts,
        WorkerContext& context,
        SitesVcfHandle& sites_handle,
        IndexedGafHandle& gaf_handle,
        const std::string& graph_query_contig,
        const std::unordered_map<std::string, std::string>& chrom_remap,
        const VariantKeySet* private_keys,
        const BamAuthorityIntervals* bam_authority,
        bool keep_intermediates = false) {
    PhasingChunk chunk;
    chunk.region = region;
    load_and_prepare_chunk(chunk, opts, context);

    // Steps 1-2: BAM candidate discovery + classification.
    collect_var_classify(chunk, opts, context.primary_header());
    // In additive mode the whitelist exists only to scope WHERE noisy-region
    // MSA runs (collect_noisy_vars_step4's containment check) and which MSA
    // calls escalate through the SNP/indel tiers; it is not an ownership
    // boundary, so every real BAM candidate is kept exactly as when no
    // --private-sites is given at all. Restricting retention to whitelist
    // keys here was silently dropping every non-whitelisted clean BAM
    // candidate inside the region -- including ones already trustworthy on
    // their own -- which is the opposite of "add extra sites".
    if (private_keys != nullptr && !opts.private_msa_admit_all_in_region) {
        CandidateTable retained;
        retained.reserve(chunk.candidates.size());
        for (CandidateVariant& candidate : chunk.candidates) {
            if (private_keys->find(candidate.key) != private_keys->end() ||
                is_bam_authoritative_position(
                    bam_authority, candidate.key.tid, candidate.key.sort_pos())) {
                retained.push_back(std::move(candidate));
            }
        }
        chunk.candidates = std::move(retained);
    }

    // Load graph sites and GAF reads for this region.
    GraphSiteCatalog chunk_catalog = load_sites_for_region(
        sites_handle, graph_query_contig, region.beg, region.end);
    for (GraphSite& s : chunk_catalog.sites) {
        auto it = chrom_remap.find(s.chrom);
        if (it != chrom_remap.end()) s.chrom = it->second;
        if (!s.ref_contig.empty()) {
            auto it2 = chrom_remap.find(s.ref_contig);
            if (it2 != chrom_remap.end()) s.ref_contig = it2->second;
        }
    }
    GraphSiteCatalogView chunk_view;
    chunk_view.source = &chunk_catalog.sites;
    chunk_view.indices.reserve(chunk_catalog.sites.size());
    for (size_t site_i = 0; site_i < chunk_catalog.sites.size(); ++site_i) {
        if (!is_bam_authoritative_position(
                bam_authority, region.tid, chunk_catalog.sites[site_i].pos)) {
            chunk_view.indices.push_back(site_i);
        }
    }

    std::vector<GraphReadAllele> chunk_rows;
    if (!chunk_view.empty()) {
        chunk_rows = scan_indexed_gaf_chunk(
            gaf_handle, graph_query_contig,
            region.beg - 1, region.end,
            chunk_view, opts.min_mapq);
    }

    // Phase A: add graph-only sites to candidate table.  Sites are added
    // unclassified (flag 0) and only gated into k-means later, once their
    // allele counts are final (see classify_graph_only_candidates below).
    SiteToCandidateMap site_map;
    std::unordered_set<int> graph_only_cands;
    std::unordered_set<int> all_graph_cands;
    GraphOnlyVcfAlleles graph_only_vcf_alleles;
    int bridged = 0, added = 0;
    if (!chunk_view.empty()) {
        site_map = inject_graph_sites(
            chunk, chunk_view, chrom_remap, opts, &bridged, &added,
            &graph_only_cands, &graph_only_vcf_alleles, &all_graph_cands);
    }

    // In authoritative mode, graph observations own every catalog-matched
    // candidate, including exact BAM/graph matches, and BAM's OWN read
    // evidence at those matches is discarded in favor of the graph's
    // (clear_bam_evidence_at_graph_candidates below). That is a real,
    // chromosome-wide change to how every graph/BAM match is resolved, not
    // something implied by merely supplying a whitelist -- in additive mode
    // (private_msa_admit_all_in_region) the whitelist exists only to scope
    // MSA and must not silently switch every other candidate's evidence
    // source too. Confirmed by measurement: enabling it unconditionally here
    // was clearing BAM evidence chromosome-wide and cost +720 discordant
    // reads on chr20 for a change that was supposed to be purely additive.
    const bool graph_authoritative =
        opts.graph_authoritative ||
        (private_keys != nullptr && !opts.private_msa_admit_all_in_region);
    const std::unordered_set<int>& graph_owned_cands =
        graph_authoritative ? all_graph_cands : graph_only_cands;

    // Step 3.1: build BAM read profiles against augmented candidate table.
    collect_var_build_profiles(chunk, opts);

    // When the graph owns these sites, the BAM's observations at them are
    // discarded rather than counted, so the graph's own evidence is not diluted
    // by an alignment call at the same position. The counts for the
    // non-authoritative case are derived after Phase B, below.
    if (graph_authoritative)
        clear_bam_evidence_at_graph_candidates(chunk, graph_owned_cands);

    // Phase B: inject graph-only reads and extend doubly-mapped profiles.
    const size_t n_bam_reads = chunk.reads.size();
    int reads_injected = 0;
    int reads_extended = 0;
    if (!chunk_rows.empty() && !site_map.empty()) {
        reads_injected = inject_graph_reads(
            chunk, chunk_rows, site_map, graph_owned_cands, opts,
            &reads_extended);
    }

    // Append graph-only reads to ordered_read_ids so k-means Phase 2
    // and update_read_phase_set visit them (they use ordered_read_ids).
    if (reads_injected > 0 && !chunk.ordered_read_ids.empty()) {
        for (size_t i = n_bam_reads; i < chunk.reads.size(); ++i)
            chunk.ordered_read_ids.push_back(static_cast<int>(i));
    }

    // Derive the graph-only candidates' counts from the final read profiles.
    // This has to follow Phase B: the injected and extended profiles carry
    // observations that did not exist when the BAM profiles were built, and
    // deriving the counts in one sweep afterwards is what lets the injection
    // sites write alleles only. Reading the profiles before Phase B is what
    // made each of those sites keep counts of its own.
    if (!graph_authoritative)
        backfill_graph_candidate_counts(chunk, graph_only_cands);

    // Gate graph-only candidates with the BAM pipeline's depth/AF/het
    // thresholds now that counts are final, then run the indel noise filter
    // on the promoted CleanHetIndel candidates.  Order matters: the noise
    // filter only acts on CleanHetIndel, so it must follow classification.
    int promoted = 0;
    if (!graph_owned_cands.empty()) {
        promoted = classify_graph_only_candidates(chunk, graph_owned_cands, opts);
        if (!chunk.ref_seq.empty()) {
            apply_hybrid_noise_filter(
                chunk, chunk.ref_seq, chunk.ref_beg, chunk.ref_end,
                graph_owned_cands, opts.noisy_reg_max_xgaps,
                &graph_only_vcf_alleles, opts.exp_hybrid_trim);
        }
    }

    if (opts.verbose >= 1 && (added > 0 || reads_injected > 0 || reads_extended > 0)) {
        std::fprintf(stderr,
            "hybrid chunk %d: bridged=%d added=%d graph_owned=%zu promoted=%d "
            "reads_injected=%d reads_extended=%d\n",
            region.chunk_id, bridged, added, graph_owned_cands.size(), promoted,
            reads_injected, reads_extended);
    }

    std::vector<uint32_t> discovery_flags;
    if (opts.recover_gaps) {
        discovery_flags.reserve(chunk.candidates.size());
        for (auto& candidate : chunk.candidates) {
            discovery_flags.push_back(candidate.lcd_var_i_to_cate);
            if (!candidate.graph_site) candidate.lcd_var_i_to_cate = 0;
        }
    }
    // Steps 3.2-4: k-means + noisy-region MSA. Private-site MSA is experimental:
    // when enabled, admit only exact-whitelist calls and phase them with the
    // clean graph core. The default preserves the private ownership boundary
    // without running MSA recall.
    if (private_keys != nullptr) {
        Options private_phase_opts = opts;
        if (opts.private_msa) {
            private_phase_opts.skip_noisy_kmeans = false;
            collect_var_run_phasing(chunk, private_phase_opts, private_keys);
        } else {
            private_phase_opts.max_noisy_reg_len = 0;
            collect_var_run_phasing(chunk, private_phase_opts);
        }
    } else if (opts.private_msa) {
        // No whitelist: private_msa's consensus-rescoring bridge-read admission
        // (align.cpp) still applies to every noisy MSA region, and there is no
        // ownership boundary to enforce, so run noisy MSA + k-means unrestricted.
        Options global_msa_opts = opts;
        global_msa_opts.skip_noisy_kmeans = false;
        collect_var_run_phasing(chunk, global_msa_opts);
    } else {
        collect_var_run_phasing(chunk, opts);
    }

    // The solve above saw only catalog sites. Where it left reads unphased, admit
    // the BAM's own sites IN THOSE WINDOWS and solve again, rather than leaving
    // the region to a later recovery pass. The admission is confined to the
    // failed windows because admitting these sites chunk-wide is a bad trade:
    // chromosome-wide on chr20 it doubled the read Hamming error (0.878% ->
    // 1.837%, 1,831 -> 3,415 discordant reads) while spanning 14 of 196 gaps.
    // `discovery_flags` is only populated under recover_gaps, where every
    // non-graph category is zeroed before the solve. Gating the retry on it made
    // --retry-unphased-with-bam a no-op on its own: without recovery there is
    // nothing to restore, `readmitted` stayed 0, and the re-solve below -- the
    // part that actually admits the noisy class -- never ran. The two concerns
    // are separate. Restoring zeroed categories is recovery-specific; asking for
    // the noisy-region pass inside a window the first solve could not phase is
    // not, and is exactly what the default path needs: on
    // chr20:55,843,827-55,889,113 the chunk already holds all three of the
    // interior sites a competitor crosses on (55,862,240 AATGGC>. at 31/29,
    // 55,862,270 T>CAGTAAATTAATTATC at 31/29, 55,883,020 ATAT>. at 17/18), each
    // classified NoisyCandHet with correct het depths and each left at
    // phase_set 0, because the hybrid subcommand sets skip_noisy_kmeans = true
    // (hybrid_collect.cpp:128) and that class never enters phasing.
    if (opts.retry_unphased_with_bam) {
        const auto windows = collect_unphased_windows(
            chunk, opts.retry_min_unphased_reads, opts.retry_min_window_bp);
        int readmitted = 0;
        if (!windows.empty() && !discovery_flags.empty()) {
            for (size_t vi = 0; vi < chunk.candidates.size(); ++vi) {
                const hts_pos_t pos = chunk.candidates[vi].key.pos;
                bool inside = false;
                for (const auto& [beg, end] : windows) {
                    if (pos >= beg && pos < end) { inside = true; break; }
                }
                if (!inside) continue;
                if (chunk.candidates[vi].graph_site) {
                    // The graph's sites are what left this window unphased, so
                    // the arm asks what the window looks like without them:
                    // discard them inside it and re-solve on verified alignment
                    // evidence alone. Zeroing the category removes the candidate
                    // from every het and link mask, as the second pass does.
                    if (opts.gap_bam_only) chunk.candidates[vi].lcd_var_i_to_cate = 0;
                    continue;
                }
                if (opts.gap_bam_only && !chunk.candidates[vi].msa_verified) {
                    chunk.candidates[vi].lcd_var_i_to_cate = 0;
                    continue;
                }
                chunk.candidates[vi].lcd_var_i_to_cate = discovery_flags[vi];
                ++readmitted;
            }
        }
        if (opts.verbose >= 1 && !windows.empty()) {
            for (const auto& [wbeg, wend] : windows) {
                int in_win = 0, in_win_graph = 0;
                for (const auto& cand : chunk.candidates) {
                    if (cand.key.pos < wbeg || cand.key.pos >= wend) continue;
                    ++in_win;
                    in_win_graph += cand.graph_site ? 1 : 0;
                }
                std::fprintf(stderr,
                    "hybrid chunk %d retry window %ld-%ld (%.1f kb): "
                    "%d candidate(s) inside, %d of them graph sites\n",
                    region.chunk_id, (long)wbeg, (long)wend,
                    (double)(wend - wbeg) / 1000.0, in_win, in_win_graph);
                std::map<uint32_t, std::pair<int, int>> by_cate;  // saved cate -> {total, phased}
                for (size_t vi = 0; vi < chunk.candidates.size(); ++vi) {
                    const auto& cand = chunk.candidates[vi];
                    if (cand.key.pos < wbeg || cand.key.pos >= wend) continue;
                    const uint32_t saved =
                        vi < discovery_flags.size() ? discovery_flags[vi] : cand.lcd_var_i_to_cate;
                    if (saved == 0) continue;
                    auto& e = by_cate[saved];
                    ++e.first;
                    e.second += cand.phase_set >= 0 ? 1 : 0;
                }
                for (const auto& [cate, counts] : by_cate) {
                    std::fprintf(stderr,
                        "      cate=0x%03x total=%d phased=%d\n",
                        cate, counts.first, counts.second);
                }
            }
            std::fprintf(stderr,
                "hybrid chunk %d retry: %zu unphased window(s), %d site(s) admitted\n",
                region.chunk_id, windows.size(), readmitted);
        }
        // Re-solve whenever a window failed, not only when a category had to be
        // restored: under stock defaults nothing was zeroed, so `readmitted` is
        // 0 while the sites are sitting there unphased.
        if (!windows.empty()) {
            Options retry_opts = opts;
            // collect_var_run_phasing skips the noisy-region MSA outright while
            // recover_gaps is set (collect_var.cpp), deferring it to the recovery
            // pass, so the noisy het class never exists in the first solve: inside
            // chr20:48,176,831-48,229,447 the chunk holds 1 clean het SNP, 2 clean
            // het indels and no NoisyCandHet at all, while the BAM channel run on
            // the same interval calls 8 of them and phases the window into one
            // block. The retry is the second try the deferral assumes, so it runs
            // that step here rather than leaving the region to recovery -- and the
            // noisy k-means with it, since orienting those candidates is the point
            // of admitting them.
            // Ask for the MSA step by name. Clearing recover_gaps would also
            // turn off every other guard gated on it -- it silently disabled
            // split_nested_msa_deletions, which then emitted both nested forms
            // of one tandem-repeat deletion as independent hets.
            retry_opts.force_noisy_msa = true;
            retry_opts.skip_noisy_kmeans = false;
            retry_opts.retry_windows = windows;
            collect_var_run_phasing(chunk, retry_opts);
        }
    }

    if (opts.recover_gaps) {
        for (size_t vi = 0; vi < discovery_flags.size(); ++vi)
            chunk.candidates[vi].lcd_var_i_to_cate = discovery_flags[vi];
    }

    // Recovery still consumes the candidate-indexed profiles. Pruning here
    // would shift candidates without remapping those profiles.
    if (!opts.recover_gaps) prune_not_candidate_variants(chunk);
    if (!keep_intermediates) mid_free_chunk(chunk, opts);
    return chunk;
}


/// Parallel batch for hybrid pipeline.
static ChunkBatchResult collect_hybrid_chunk_batch_parallel(
        const Options& opts,
        const std::vector<RegionChunk>& chunks,
        size_t batch_begin,
        size_t batch_end,
        const std::string& graph_query_contig,
        const std::unordered_map<std::string, std::string>& chrom_remap,
        const VariantKeySet* private_keys,
        const BamAuthorityIntervals* bam_authority) {
    const size_t batch_size = batch_end - batch_begin;
    ChunkBatchResult result;
    result.chunks.resize(batch_size);

    const size_t worker_count = std::min<size_t>(
        static_cast<size_t>(opts.threads), batch_size);
    std::atomic<size_t> next_offset{0};
    std::exception_ptr first_error;
    std::mutex error_mutex;

    std::vector<std::thread> workers;
    workers.reserve(worker_count);

    for (size_t worker_i = 0; worker_i < worker_count; ++worker_i) {
        workers.emplace_back([&, worker_i]() {
            (void)worker_i;
            try {
                WorkerContext context(opts);
                SitesVcfHandle sites_handle(opts.graph_sites_vcf);
                IndexedGafHandle gaf_handle(opts.gaf_file);

                while (true) {
                    const size_t offset = next_offset.fetch_add(1);
                    if (offset >= batch_size) break;
                    result.chunks[offset] = process_chunk_hybrid(
                        chunks[batch_begin + offset], opts, context,
                        sites_handle, gaf_handle,
                        graph_query_contig, chrom_remap, private_keys,
                        bam_authority, opts.recover_gaps);
                }
            } catch (...) {
                std::lock_guard<std::mutex> lock(error_mutex);
                if (!first_error) first_error = std::current_exception();
            }
        });
    }
    for (std::thread& w : workers) w.join();
    if (first_error) std::rethrow_exception(first_error);
    return result;
}

static void filter_hybrid_reads_by_margin(std::vector<PhasingChunk>& chunks,
                                          int min_margin,
                                          bool credit_bridge_snps) {
    if (min_margin <= 0) return;
    for (PhasingChunk& chunk : chunks) {
        for (size_t i = 0; i < chunk.reads.size(); ++i) {
            const ReadRecord& read = chunk.reads[i];
            if (read.n_clean_agree_snps - read.n_clean_conflict_snps >= min_margin)
                continue;
            // A read whose only informative sites in a stretch are
            // MSA-admitted bridge SNPs (see --private-msa-admit-all-in-region)
            // has zero clean-SNP margin by construction -- CleanHetSnp is a
            // different category -- and would be silently stripped here even
            // when its haplotype assignment already used that evidence via
            // hap_scores.  Credit it only when the admission mode that
            // produced it is active, and only using the bridge SNP counters
            // (never indels: an MSA indel can be placed at several equivalent
            // positions inside a repeat run and is not trustworthy enough to
            // rescue an otherwise-thin read).
            if (credit_bridge_snps &&
                (read.n_bridge_agree_snps > 0 || read.n_hp_gap_agree > 0)) {
                // hp_gap_scorable observations are credited here for the same
                // reason as bridge SNPs: the read's haplotype was assigned using
                // that evidence, so judging it on clean-SNP counts alone strips a
                // read the solve had already phased. This is the narrow case the
                // indel caveat above does not cover -- the site is MSA-verified
                // and inside the homopolymer tier's own gap window, and it is
                // frequently the only interior evidence such a gap has.
                const int bridge_margin =
                    (read.n_clean_agree_snps + read.n_bridge_agree_snps +
                     read.n_hp_gap_agree) -
                    (read.n_clean_conflict_snps + read.n_bridge_conflict_snps +
                     read.n_hp_gap_conflict);
                if (bridge_margin >= min_margin) continue;
            }
            if (i < chunk.haps.size()) chunk.haps[i] = 0;
            if (i < chunk.phase_sets.size()) chunk.phase_sets[i] = -1;
        }
    }
}

static void filter_hybrid_small_phase_sets(std::vector<PhasingChunk>& chunks,
                                           int min_reads) {
    if (min_reads <= 0) return;
    std::unordered_map<hts_pos_t, std::unordered_set<std::string>> phase_set_reads;
    for (const PhasingChunk& chunk : chunks) {
        for (size_t i = 0; i < chunk.haps.size() &&
                           i < chunk.phase_sets.size() &&
                           i < chunk.reads.size(); ++i) {
            if ((chunk.haps[i] == 1 || chunk.haps[i] == 2) &&
                chunk.phase_sets[i] >= 0) {
                std::string read_key = std::to_string(chunk.reads[i].input_index);
                read_key.push_back('\0');
                read_key += chunk.reads[i].qname;
                phase_set_reads[chunk.phase_sets[i]].insert(std::move(read_key));
            }
        }
    }
    for (PhasingChunk& chunk : chunks) {
        for (size_t i = 0; i < chunk.haps.size() && i < chunk.phase_sets.size(); ++i) {
            const hts_pos_t phase_set = chunk.phase_sets[i];
            const auto found = phase_set_reads.find(phase_set);
            if (phase_set < 0 ||
                (found != phase_set_reads.end() &&
                 static_cast<int>(found->second.size()) >= min_reads)) {
                continue;
            }
            chunk.haps[i] = 0;
            chunk.phase_sets[i] = -1;
        }
    }
}

static constexpr hts_pos_t kGapRecoveryFlank = 50000;
static constexpr hts_pos_t kGapRecoveryMsaFlank = 5000;
static constexpr hts_pos_t kGapRecoveryMaxMsaSpan = 250000;
static constexpr uint32_t kGapEvidenceCacheVersion = 6;
static constexpr char kGapEvidenceCacheMagic[] = "PGGAPEV";

static void gap_cache_hash_bytes(uint64_t& hash, const void* data, size_t size) {
    constexpr uint64_t kFnvPrime = 1099511628211ULL;
    const auto* bytes = static_cast<const unsigned char*>(data);
    for (size_t i = 0; i < size; ++i) {
        hash ^= bytes[i];
        hash *= kFnvPrime;
    }
}

template <typename T>
static void gap_cache_hash_value(uint64_t& hash, const T& value) {
    gap_cache_hash_bytes(hash, &value, sizeof(value));
}

static void gap_cache_hash_string(uint64_t& hash, const std::string& value) {
    gap_cache_hash_value(hash, value.size());
    gap_cache_hash_bytes(hash, value.data(), value.size());
}

static uint64_t gap_evidence_input_signature(
        const std::vector<PhasingChunk>& chunks,
        const std::vector<PhaseGap>& gaps, const Options& opts) {
    uint64_t hash = 14695981039346656037ULL;
    for (const auto value : {
             opts.min_mapq, opts.min_bq, opts.min_depth, opts.min_alt_depth,
             opts.max_noisy_reg_len, opts.max_noisy_reg_cov,
             opts.min_hap_full_reads, opts.min_hap_reads,
             opts.min_noisy_reg_size_to_sample_reads, opts.noisy_reg_flank_len,
             opts.match, opts.mismatch, opts.gap_open1, opts.gap_ext1,
             opts.gap_open2, opts.gap_ext2, opts.gap_aln,
             opts.private_msa_margin}) {
        gap_cache_hash_value(hash, value);
    }
    gap_cache_hash_value(hash, opts.partial_aln_ratio);
    gap_cache_hash_value(hash, opts.read_technology);
    for (const PhaseGap& gap : gaps) {
        gap_cache_hash_value(hash, gap.tid);
        gap_cache_hash_value(hash, gap.left_ps);
        gap_cache_hash_value(hash, gap.right_ps);
        gap_cache_hash_value(hash, gap.left_end);
        gap_cache_hash_value(hash, gap.right_beg);
    }
    for (const PhasingChunk& chunk : chunks) {
        gap_cache_hash_value(hash, chunk.region.tid);
        gap_cache_hash_value(hash, chunk.region.beg);
        gap_cache_hash_value(hash, chunk.region.end);
        for (const ReadRecord& read : chunk.reads) {
            gap_cache_hash_value(hash, read.input_index);
            gap_cache_hash_value(hash, read.beg);
            gap_cache_hash_value(hash, read.end);
            gap_cache_hash_string(hash, read.qname);
        }
        for (const int hap : chunk.haps) gap_cache_hash_value(hash, hap);
        for (const hts_pos_t phase_set : chunk.phase_sets)
            gap_cache_hash_value(hash, phase_set);
        for (const CandidateVariant& candidate : chunk.candidates) {
            gap_cache_hash_value(hash, candidate.key.tid);
            gap_cache_hash_value(hash, candidate.key.pos);
            gap_cache_hash_value(hash, candidate.key.type);
            gap_cache_hash_value(hash, candidate.key.ref_len);
            gap_cache_hash_string(hash, candidate.key.alt);
            gap_cache_hash_value(hash, candidate.counts.total_cov);
            gap_cache_hash_value(hash, candidate.counts.ref_cov);
            gap_cache_hash_value(hash, candidate.counts.alt_cov);
            gap_cache_hash_value(hash, candidate.counts.category);
            gap_cache_hash_value(hash, candidate.lcd_var_i_to_cate);
            gap_cache_hash_value(hash, candidate.phase_set);
            gap_cache_hash_value(hash, candidate.hap_alt);
            gap_cache_hash_value(hash, candidate.hap_ref);
            for (const int allele : candidate.hap_to_cons_alle)
                gap_cache_hash_value(hash, allele);
        }
    }
    return hash;
}

template <typename T>
static void write_gap_cache_value(std::ostream& out, const T& value) {
    out.write(reinterpret_cast<const char*>(&value), sizeof(value));
}

template <typename T>
static void read_gap_cache_value(std::istream& in, T& value) {
    in.read(reinterpret_cast<char*>(&value), sizeof(value));
}

static void write_gap_cache_string(std::ostream& out, const std::string& value) {
    const uint64_t size = value.size();
    write_gap_cache_value(out, size);
    out.write(value.data(), static_cast<std::streamsize>(size));
}

static void read_gap_cache_string(std::istream& in, std::string& value) {
    uint64_t size = 0;
    read_gap_cache_value(in, size);
    constexpr uint64_t kMaxCacheStringSize = 1ULL << 30;
    if (size > kMaxCacheStringSize)
        throw std::runtime_error("invalid gap evidence cache string length");
    value.resize(static_cast<size_t>(size));
    in.read(value.data(), static_cast<std::streamsize>(size));
}

template <typename T>
static void write_gap_cache_vector(std::ostream& out,
                                   const std::vector<T>& values) {
    const uint64_t size = values.size();
    write_gap_cache_value(out, size);
    if (size > 0) {
        out.write(reinterpret_cast<const char*>(values.data()),
                  static_cast<std::streamsize>(size * sizeof(T)));
    }
}

template <typename T>
static void read_gap_cache_vector(std::istream& in, std::vector<T>& values) {
    uint64_t size = 0;
    read_gap_cache_value(in, size);
    constexpr uint64_t kMaxCacheVectorElements = 1ULL << 32;
    if (size > kMaxCacheVectorElements)
        throw std::runtime_error("invalid gap evidence cache vector length");
    values.resize(static_cast<size_t>(size));
    if (size > 0) {
        in.read(reinterpret_cast<char*>(values.data()),
                static_cast<std::streamsize>(size * sizeof(T)));
    }
}

static void write_gap_cache_candidate(std::ostream& out,
                                      const CandidateVariant& candidate) {
    write_gap_cache_value(out, candidate.key.tid);
    write_gap_cache_value(out, candidate.key.pos);
    write_gap_cache_value(out, candidate.key.type);
    write_gap_cache_value(out, candidate.key.ref_len);
    write_gap_cache_string(out, candidate.key.alt);
    write_gap_cache_value(out, candidate.counts.total_cov);
    write_gap_cache_value(out, candidate.counts.ref_cov);
    write_gap_cache_value(out, candidate.counts.alt_cov);
    write_gap_cache_value(out, candidate.counts.low_qual_cov);
    write_gap_cache_value(out, candidate.counts.forward_ref);
    write_gap_cache_value(out, candidate.counts.reverse_ref);
    write_gap_cache_value(out, candidate.counts.forward_alt);
    write_gap_cache_value(out, candidate.counts.reverse_alt);
    write_gap_cache_value(out, candidate.counts.n_uniq_alles);
    write_gap_cache_vector(out, candidate.counts.alle_covs);
    write_gap_cache_value(out, candidate.counts.category);
    write_gap_cache_value(out, candidate.counts.candvarcate_initial);
    write_gap_cache_value(out, candidate.counts.allele_fraction);
    const uint64_t n_alts = candidate.msa_insertion_alts.size();
    write_gap_cache_value(out, n_alts);
    for (const std::string& alt : candidate.msa_insertion_alts)
        write_gap_cache_string(out, alt);
    write_gap_cache_value(out, candidate.ref_base);
    write_gap_cache_value(out, candidate.alt_ref_base);
    write_gap_cache_value(out, candidate.phase_set);
    write_gap_cache_value(out, candidate.hap_alt);
    write_gap_cache_value(out, candidate.hap_ref);
    write_gap_cache_value(out, candidate.is_homopolymer_indel);
    write_gap_cache_value(out, candidate.gap_link_supported);
    write_gap_cache_value(out, candidate.msa_verified);
    write_gap_cache_value(out, candidate.graph_site);
    write_gap_cache_value(out, candidate.lcd_make_variants_region_pass);
    write_gap_cache_value(out, candidate.lcd_var_i_to_cate);
    for (const auto& profile : candidate.hap_to_alle_profile)
        write_gap_cache_vector(out, profile);
    for (const int allele : candidate.hap_to_cons_alle)
        write_gap_cache_value(out, allele);
}

static CandidateVariant read_gap_cache_candidate(std::istream& in) {
    CandidateVariant candidate;
    read_gap_cache_value(in, candidate.key.tid);
    read_gap_cache_value(in, candidate.key.pos);
    read_gap_cache_value(in, candidate.key.type);
    read_gap_cache_value(in, candidate.key.ref_len);
    read_gap_cache_string(in, candidate.key.alt);
    read_gap_cache_value(in, candidate.counts.total_cov);
    read_gap_cache_value(in, candidate.counts.ref_cov);
    read_gap_cache_value(in, candidate.counts.alt_cov);
    read_gap_cache_value(in, candidate.counts.low_qual_cov);
    read_gap_cache_value(in, candidate.counts.forward_ref);
    read_gap_cache_value(in, candidate.counts.reverse_ref);
    read_gap_cache_value(in, candidate.counts.forward_alt);
    read_gap_cache_value(in, candidate.counts.reverse_alt);
    read_gap_cache_value(in, candidate.counts.n_uniq_alles);
    read_gap_cache_vector(in, candidate.counts.alle_covs);
    read_gap_cache_value(in, candidate.counts.category);
    read_gap_cache_value(in, candidate.counts.candvarcate_initial);
    read_gap_cache_value(in, candidate.counts.allele_fraction);
    uint64_t n_alts = 0;
    read_gap_cache_value(in, n_alts);
    if (n_alts > 1000000)
        throw std::runtime_error("invalid gap evidence cache alternate count");
    candidate.msa_insertion_alts.resize(static_cast<size_t>(n_alts));
    for (std::string& alt : candidate.msa_insertion_alts)
        read_gap_cache_string(in, alt);
    read_gap_cache_value(in, candidate.ref_base);
    read_gap_cache_value(in, candidate.alt_ref_base);
    read_gap_cache_value(in, candidate.phase_set);
    read_gap_cache_value(in, candidate.hap_alt);
    read_gap_cache_value(in, candidate.hap_ref);
    read_gap_cache_value(in, candidate.is_homopolymer_indel);
    read_gap_cache_value(in, candidate.gap_link_supported);
    read_gap_cache_value(in, candidate.msa_verified);
    read_gap_cache_value(in, candidate.graph_site);
    read_gap_cache_value(in, candidate.lcd_make_variants_region_pass);
    read_gap_cache_value(in, candidate.lcd_var_i_to_cate);
    for (auto& profile : candidate.hap_to_alle_profile)
        read_gap_cache_vector(in, profile);
    for (int& allele : candidate.hap_to_cons_alle)
        read_gap_cache_value(in, allele);
    return candidate;
}

static void write_gap_evidence_cache(const std::string& path,
                                     const std::vector<PhasingChunk>& chunks,
                                     uint64_t input_signature) {
    const std::string temporary = path + ".tmp";
    std::ofstream out(temporary, std::ios::binary | std::ios::trunc);
    if (!out) throw std::runtime_error("failed to create gap evidence cache: " + path);
    out.write(kGapEvidenceCacheMagic, sizeof(kGapEvidenceCacheMagic));
    write_gap_cache_value(out, kGapEvidenceCacheVersion);
    write_gap_cache_value(out, input_signature);
    const uint64_t n_chunks = chunks.size();
    write_gap_cache_value(out, n_chunks);
    for (const PhasingChunk& chunk : chunks) {
        write_gap_cache_value(out, chunk.region.tid);
        write_gap_cache_value(out, chunk.region.beg);
        write_gap_cache_value(out, chunk.region.end);
        const uint64_t n_reads = chunk.reads.size();
        write_gap_cache_value(out, n_reads);
        for (const ReadRecord& read : chunk.reads) {
            write_gap_cache_value(out, read.input_index);
            write_gap_cache_value(out, read.beg);
            write_gap_cache_value(out, read.end);
            write_gap_cache_string(out, read.qname);
        }
        const uint64_t n_candidates = chunk.candidates.size();
        write_gap_cache_value(out, n_candidates);
        for (const CandidateVariant& candidate : chunk.candidates)
            write_gap_cache_candidate(out, candidate);
        const uint64_t n_profiles = chunk.read_var_profile.size();
        write_gap_cache_value(out, n_profiles);
        for (const ReadVariantProfile& profile : chunk.read_var_profile) {
            write_gap_cache_value(out, profile.read_id);
            write_gap_cache_value(out, profile.start_var_idx);
            write_gap_cache_value(out, profile.end_var_idx);
            write_gap_cache_vector(out, profile.alleles);
            write_gap_cache_vector(out, profile.alt_qi);
            write_gap_cache_vector(out, profile.graph_alleles);
            write_gap_cache_vector(out, profile.bam_alleles);
            write_gap_cache_vector(out, profile.bam_qi);
        }
    }
    out.close();
    if (!out) throw std::runtime_error("failed to write gap evidence cache: " + path);
    std::error_code error;
    std::filesystem::rename(temporary, path, error);
    if (error)
        throw std::runtime_error("failed to finalize gap evidence cache " + path +
                                 ": " + error.message());
}

static void read_gap_evidence_cache(const std::string& path,
                                    std::vector<PhasingChunk>& chunks,
                                    uint64_t expected_signature) {
    std::ifstream in(path, std::ios::binary);
    if (!in) throw std::runtime_error("failed to open gap evidence cache: " + path);
    char magic[sizeof(kGapEvidenceCacheMagic)]{};
    in.read(magic, sizeof(magic));
    uint32_t version = 0;
    read_gap_cache_value(in, version);
    if (std::memcmp(magic, kGapEvidenceCacheMagic, sizeof(magic)) != 0 ||
        version != kGapEvidenceCacheVersion) {
        throw std::runtime_error("incompatible gap evidence cache: " + path);
    }
    uint64_t input_signature = 0;
    read_gap_cache_value(in, input_signature);
    if (input_signature != expected_signature) {
        throw std::runtime_error(
            "gap evidence cache inputs or MSA settings do not match this run; "
            "remove and rebuild: " + path);
    }
    uint64_t n_chunks = 0;
    read_gap_cache_value(in, n_chunks);
    if (n_chunks != chunks.size())
        throw std::runtime_error("gap evidence cache chunk count does not match this run: " + path);
    for (PhasingChunk& chunk : chunks) {
        int tid = -1;
        hts_pos_t beg = 0, end = 0;
        read_gap_cache_value(in, tid);
        read_gap_cache_value(in, beg);
        read_gap_cache_value(in, end);
        if (tid != chunk.region.tid || beg != chunk.region.beg || end != chunk.region.end)
            throw std::runtime_error("gap evidence cache regions do not match this run: " + path);
        uint64_t n_reads = 0;
        read_gap_cache_value(in, n_reads);
        if (n_reads != chunk.reads.size())
            throw std::runtime_error("gap evidence cache reads do not match this run: " + path);
        for (const ReadRecord& read : chunk.reads) {
            int input_index = -1;
            hts_pos_t read_beg = 0, read_end = 0;
            std::string qname;
            read_gap_cache_value(in, input_index);
            read_gap_cache_value(in, read_beg);
            read_gap_cache_value(in, read_end);
            read_gap_cache_string(in, qname);
            if (input_index != read.input_index || read_beg != read.beg ||
                read_end != read.end || qname != read.qname) {
                throw std::runtime_error("gap evidence cache read order does not match this run: " + path);
            }
        }
        uint64_t n_candidates = 0;
        read_gap_cache_value(in, n_candidates);
        if (n_candidates > (1ULL << 32))
            throw std::runtime_error("invalid gap evidence cache candidate count: " + path);
        CandidateTable candidates;
        candidates.reserve(static_cast<size_t>(n_candidates));
        for (uint64_t i = 0; i < n_candidates; ++i)
            candidates.push_back(read_gap_cache_candidate(in));
        uint64_t n_profiles = 0;
        read_gap_cache_value(in, n_profiles);
        if (n_profiles != chunk.reads.size())
            throw std::runtime_error("gap evidence cache profiles do not match this run: " + path);
        std::vector<ReadVariantProfile> profiles(static_cast<size_t>(n_profiles));
        for (ReadVariantProfile& profile : profiles) {
            read_gap_cache_value(in, profile.read_id);
            read_gap_cache_value(in, profile.start_var_idx);
            read_gap_cache_value(in, profile.end_var_idx);
            read_gap_cache_vector(in, profile.alleles);
            read_gap_cache_vector(in, profile.alt_qi);
            read_gap_cache_vector(in, profile.graph_alleles);
            read_gap_cache_vector(in, profile.bam_alleles);
            read_gap_cache_vector(in, profile.bam_qi);
        }
        chunk.candidates = std::move(candidates);
        chunk.read_var_profile = std::move(profiles);
    }
    if (!in) throw std::runtime_error("truncated gap evidence cache: " + path);
}

static void populate_gap_msa_cache(std::vector<PhasingChunk>& chunks,
                                   const std::vector<PhaseGap>& gaps,
                                   const Options& opts) {
    std::atomic<size_t> next_chunk{0};
    std::exception_ptr first_error;
    std::mutex error_mutex;
    const size_t worker_count = std::min<size_t>(
        static_cast<size_t>(opts.threads), chunks.size());
    std::vector<std::thread> workers;
    workers.reserve(worker_count);
    for (size_t worker_i = 0; worker_i < worker_count; ++worker_i) {
        workers.emplace_back([&]() {
            try {
                while (true) {
                    const size_t chunk_i = next_chunk.fetch_add(1);
                    if (chunk_i >= chunks.size()) break;
                    PhasingChunk& chunk = chunks[chunk_i];
                    std::vector<Interval> intervals;
                    for (const PhaseGap& gap : gaps) {
                        if (gap.tid != chunk.region.tid ||
                            gap.right_beg - gap.left_end > kGapRecoveryMaxMsaSpan)
                            continue;
                        const hts_pos_t beg = std::max(
                            chunk.region.beg,
                            gap.left_end - kGapRecoveryMsaFlank);
                        const hts_pos_t end = std::min(
                            chunk.region.end,
                            gap.right_beg + kGapRecoveryMsaFlank);
                        if (beg <= end) intervals.push_back({beg, end, 0});
                    }
                    if (intervals.empty()) continue;
                    std::sort(intervals.begin(), intervals.end(),
                              [](const Interval& first, const Interval& second) {
                                  return first.beg < second.beg;
                              });
                    std::vector<Interval> merged;
                    for (const Interval& interval : intervals) {
                        if (merged.empty() || interval.beg > merged.back().end + 1)
                            merged.push_back(interval);
                        else
                            merged.back().end =
                                std::max(merged.back().end, interval.end);
                    }

                    const std::vector<int> saved_haps = chunk.haps;
                    const std::vector<hts_pos_t> saved_phase_sets = chunk.phase_sets;
                    struct ReadPhaseState {
                        int clean_agree;
                        int clean_conflict;
                        int bridge_agree;
                        int bridge_conflict;
                        int margin;
                        int scored;
                    };
                    std::vector<ReadPhaseState> read_states;
                    read_states.reserve(chunk.reads.size());
                    for (const ReadRecord& read : chunk.reads) {
                        read_states.push_back({
                            read.n_clean_agree_snps, read.n_clean_conflict_snps,
                            read.n_bridge_agree_snps, read.n_bridge_conflict_snps,
                            read.hap_score_margin, read.n_vars_scored});
                    }
                    struct CandidatePhaseState {
                        VariantKey key;
                        hts_pos_t phase_set;
                        int hap_alt;
                        int hap_ref;
                        std::array<std::vector<int>, 3> profiles;
                        std::array<int, 3> consensus;
                    };
                    std::vector<CandidatePhaseState> candidate_states;
                    candidate_states.reserve(chunk.candidates.size());
                    for (const CandidateVariant& candidate : chunk.candidates) {
                        candidate_states.push_back({
                            candidate.key, candidate.phase_set, candidate.hap_alt,
                            candidate.hap_ref, candidate.hap_to_alle_profile,
                            candidate.hap_to_cons_alle});
                    }
                    const std::vector<Interval> saved_noisy_regions =
                        chunk.noisy_regions;
                    for (const Interval& interval : merged) {
                        prepare_gap_msa_regions(chunk, interval.beg, interval.end);
                        run_gap_msa_tier(
                            chunk, opts, interval.beg, interval.end, true);
                        run_gap_msa_tier(
                            chunk, opts, interval.beg, interval.end, false);
                    }
                    chunk.haps = saved_haps;
                    chunk.phase_sets = saved_phase_sets;
                    for (size_t read_i = 0; read_i < chunk.reads.size(); ++read_i) {
                        ReadRecord& read = chunk.reads[read_i];
                        const ReadPhaseState& state = read_states[read_i];
                        read.n_clean_agree_snps = state.clean_agree;
                        read.n_clean_conflict_snps = state.clean_conflict;
                        read.n_bridge_agree_snps = state.bridge_agree;
                        read.n_bridge_conflict_snps = state.bridge_conflict;
                        read.hap_score_margin = state.margin;
                        read.n_vars_scored = state.scored;
                    }
                    for (const CandidatePhaseState& state : candidate_states) {
                        const auto candidate = std::lower_bound(
                            chunk.candidates.begin(), chunk.candidates.end(), state.key,
                            [](const CandidateVariant& value, const VariantKey& key) {
                                return exact_comp_var_site(&value.key, &key) < 0;
                            });
                        if (candidate == chunk.candidates.end() ||
                            exact_comp_var_site(&candidate->key, &state.key) != 0)
                            continue;
                        candidate->phase_set = state.phase_set;
                        candidate->hap_alt = state.hap_alt;
                        candidate->hap_ref = state.hap_ref;
                        candidate->hap_to_alle_profile = state.profiles;
                        candidate->hap_to_cons_alle = state.consensus;
                    }
                    VariantKeySet original_keys;
                    for (const auto& state : candidate_states) original_keys.insert(state.key);
                    for (auto& candidate : chunk.candidates) {
                        if (original_keys.count(candidate.key)) continue;
                        candidate.phase_set = -1;
                        candidate.hap_alt = candidate.hap_ref = 0;
                        candidate.hap_to_cons_alle = {-1, -1, -1};
                        for (auto& profile : candidate.hap_to_alle_profile) profile.clear();
                    }
                    chunk.noisy_regions = saved_noisy_regions;
                }
            } catch (...) {
                std::lock_guard<std::mutex> lock(error_mutex);
                if (!first_error) first_error = std::current_exception();
            }
        });
    }
    for (std::thread& worker : workers) worker.join();
    if (first_error) std::rethrow_exception(first_error);
}

static ReadRecord clone_cached_read(const ReadRecord& source) {
    ReadRecord read;
    read.tid = source.tid;
    read.input_index = source.input_index;
    read.beg = source.beg;
    read.end = source.end;
    read.reverse = source.reverse;
    read.nm = source.nm;
    read.mapq = source.mapq;
    read.qname = source.qname;
    if (source.alignment) read.alignment.reset(bam_dup1(source.alignment.get()));
    read.qual = source.qual;
    read.digars = source.digars;
    read.noisy_regions = source.noisy_regions;
    read.is_skipped = source.is_skipped;
    read.is_ont_palindrome = source.is_ont_palindrome;
    read.total_cand_events = source.total_cand_events;
    return read;
}

static void merge_cached_allele(int& dest, int source) {
    constexpr int kConflictingAllele = -3;
    if (source == -1 || dest == kConflictingAllele) return;
    if (source == kConflictingAllele || (dest >= 0 && source >= 0 && dest != source)) {
        dest = kConflictingAllele;
    } else if (dest < 0 || source >= 0) dest = source;
}

static int remap_cached_allele(int allele, const CandidateVariant& source,
                              const CandidateVariant& dest) {
    if (allele <= 0) return allele;
    std::string sequence;
    if (source.msa_insertion_alts.empty()) {
        if (allele != 1) return -3;
        sequence = source.key.alt;
    } else {
        if (static_cast<size_t>(allele) > source.msa_insertion_alts.size()) return -3;
        sequence = source.msa_insertion_alts[allele - 1];
    }
    if (dest.msa_insertion_alts.empty()) return sequence == dest.key.alt ? 1 : -3;
    const auto found = std::find(dest.msa_insertion_alts.begin(), dest.msa_insertion_alts.end(), sequence);
    return found == dest.msa_insertion_alts.end() ? -3 : static_cast<int>(found - dest.msa_insertion_alts.begin()) + 1;
}

static PhasingChunk build_cached_gap_proposal(
        const std::vector<PhasingChunk>& chunks, const RegionChunk& window) {
    PhasingChunk proposal;
    proposal.region = window;
    proposal.ref_beg = window.beg;
    proposal.ref_end = window.end;
    proposal.ref_seq.assign(
        static_cast<size_t>(window.end - window.beg + 1), 'N');

    std::vector<CandidateVariant> candidates;
    for (const PhasingChunk& source : chunks) {
        if (source.region.tid != window.tid || source.region.end < window.beg ||
            source.region.beg > window.end) {
            continue;
        }
        const hts_pos_t copy_beg = std::max(window.beg, source.ref_beg);
        const hts_pos_t copy_end = std::min(window.end, source.ref_end);
        if (copy_beg <= copy_end && !source.ref_seq.empty()) {
            proposal.ref_seq.replace(
                static_cast<size_t>(copy_beg - window.beg),
                static_cast<size_t>(copy_end - copy_beg + 1), source.ref_seq,
                static_cast<size_t>(copy_beg - source.ref_beg),
                static_cast<size_t>(copy_end - copy_beg + 1));
        }
        for (const Interval& region : source.low_complexity_regions)
            if (region.end >= window.beg && region.beg <= window.end)
                proposal.low_complexity_regions.push_back(region);
        for (const Interval& region : source.noisy_regions)
            if (region.end >= window.beg && region.beg <= window.end)
                proposal.noisy_regions.push_back(region);
        for (const CandidateVariant& candidate : source.candidates)
            if (candidate.key.sort_pos() >= window.beg &&
                candidate.key.sort_pos() <= window.end)
                candidates.push_back(candidate);
    }
    std::sort(candidates.begin(), candidates.end(),
              [](const CandidateVariant& first, const CandidateVariant& second) {
                  return exact_comp_var_site(&first.key, &second.key) < 0;
              });
    for (const CandidateVariant& candidate : candidates) {
        if (!proposal.candidates.empty() &&
            exact_comp_var_site(&proposal.candidates.back().key, &candidate.key) == 0) {
            const bool graph_site = proposal.candidates.back().graph_site || candidate.graph_site;
            if (!proposal.candidates.back().lcd_make_variants_region_pass &&
                candidate.lcd_make_variants_region_pass)
                proposal.candidates.back() = candidate;
            proposal.candidates.back().graph_site = graph_site;
        } else {
            proposal.candidates.push_back(candidate);
        }
    }

    using ReadKey = std::pair<int, std::string>;
    std::map<ReadKey, size_t> read_indices;
    for (const PhasingChunk& source : chunks) {
        if (source.region.tid != window.tid || source.region.end < window.beg ||
            source.region.beg > window.end) {
            continue;
        }
        for (size_t source_i = 0; source_i < source.reads.size(); ++source_i) {
            const ReadRecord& read = source.reads[source_i];
            if (read.tid != window.tid || read.end < window.beg ||
                read.beg > window.end) {
                continue;
            }
            const ReadKey key{read.input_index, read.qname};
            auto inserted = read_indices.emplace(key, proposal.reads.size());
            if (inserted.second) {
                proposal.reads.push_back(clone_cached_read(read));
                ReadVariantProfile profile;
                profile.read_id = static_cast<int>(proposal.reads.size() - 1);
                profile.start_var_idx = proposal.candidates.empty() ? -1 : 0;
                profile.end_var_idx = static_cast<int>(proposal.candidates.size()) - 1;
                profile.alleles.assign(proposal.candidates.size(), -1);
                profile.alt_qi.assign(proposal.candidates.size(), -1);
                profile.graph_alleles.assign(proposal.candidates.size(), -1);
                profile.bam_alleles.assign(proposal.candidates.size(), -1);
                profile.bam_qi.assign(proposal.candidates.size(), -1);
                proposal.read_var_profile.push_back(std::move(profile));
            }
            if (!inserted.second && proposal.reads[inserted.first->second].is_skipped && !read.is_skipped)
                proposal.reads[inserted.first->second] = clone_cached_read(read);
            if (source_i >= source.read_var_profile.size()) continue;
            const ReadVariantProfile& source_profile =
                source.read_var_profile[source_i];
            ReadVariantProfile& dest =
                proposal.read_var_profile[inserted.first->second];
            for (int source_vi = source_profile.start_var_idx;
                 source_vi <= source_profile.end_var_idx; ++source_vi) {
                if (source_vi < 0 ||
                    static_cast<size_t>(source_vi) >= source.candidates.size())
                    continue;
                const CandidateVariant& source_candidate =
                    source.candidates[static_cast<size_t>(source_vi)];
                const auto found = std::lower_bound(
                    proposal.candidates.begin(), proposal.candidates.end(),
                    source_candidate.key,
                    [](const CandidateVariant& candidate, const VariantKey& key) {
                        return exact_comp_var_site(&candidate.key, &key) < 0;
                    });
                if (found == proposal.candidates.end() ||
                    exact_comp_var_site(&found->key, &source_candidate.key) != 0)
                    continue;
                const size_t dest_vi = static_cast<size_t>(
                    found - proposal.candidates.begin());
                const size_t allele_i = static_cast<size_t>(
                    source_vi - source_profile.start_var_idx);
                if (allele_i < source_profile.bam_alleles.size()) {
                    merge_cached_allele(dest.bam_alleles[dest_vi], source_profile.bam_alleles[allele_i]);
                    if (allele_i < source_profile.bam_qi.size())
                        dest.bam_qi[dest_vi] = source_profile.bam_qi[allele_i];
                }
                if (allele_i < source_profile.graph_alleles.size() &&
                    source_profile.graph_alleles[allele_i] >= 0)
                    merge_cached_allele(dest.graph_alleles[dest_vi], source_profile.graph_alleles[allele_i]);
                if (allele_i >= source_profile.alleles.size() ||
                    source_profile.alleles[allele_i] == -1)
                    continue;
                merge_cached_allele(dest.alleles[dest_vi],
                    remap_cached_allele(source_profile.alleles[allele_i], source_candidate, *found));
                if (allele_i < source_profile.alt_qi.size())
                    dest.alt_qi[dest_vi] = source_profile.alt_qi[allele_i];
            }
        }
    }
    proposal.haps.assign(proposal.reads.size(), 0);
    proposal.phase_sets.assign(proposal.reads.size(), -1);
    proposal.ordered_read_ids.resize(proposal.reads.size());
    std::iota(proposal.ordered_read_ids.begin(), proposal.ordered_read_ids.end(), 0);
    std::sort(proposal.ordered_read_ids.begin(), proposal.ordered_read_ids.end(),
              [&](int first, int second) {
                  return proposal.reads[static_cast<size_t>(first)].beg <
                         proposal.reads[static_cast<size_t>(second)].beg;
              });
    cgranges_t* cr = cr_init();
    for (size_t read_i = 0; read_i < proposal.read_var_profile.size(); ++read_i) {
        const ReadVariantProfile& profile = proposal.read_var_profile[read_i];
        if (profile.start_var_idx < 0) continue;
        cr_add(cr, "cr", profile.start_var_idx, profile.end_var_idx + 1,
               static_cast<int32_t>(read_i));
    }
    cr_index(cr);
    proposal.read_var_cr.reset(cr);
    return proposal;
}

static bool gap_recovery_jobs_conflict(const PhaseGap& first,
                                       const PhaseGap& second) {
    if (first.tid != second.tid) return false;
    const hts_pos_t first_beg =
        std::max(first.region_beg, first.left_end - kGapRecoveryFlank);
    const hts_pos_t first_end =
        std::min(first.region_end, first.right_beg + kGapRecoveryFlank);
    const hts_pos_t second_beg =
        std::max(second.region_beg, second.left_end - kGapRecoveryFlank);
    const hts_pos_t second_end =
        std::min(second.region_end, second.right_beg + kGapRecoveryFlank);
    return first_beg <= second_end && second_beg <= first_end;
}

struct GapRecoveryJobResult {
    bool joined = false;
    bool has_edge = false;
    GapPhaseEdge edge{};
    std::string report_rows;
    // Final proposal state, kept only when !joined, for independent-block
    // consideration strictly after EVERY mechanism that can still claim or
    // relabel a block has finished -- not just this gap's own stitch decision
    // and apply_gap_phase_edges, but also the batch-level stitch_chunk_haps
    // call that runs again after recover_hybrid_gaps returns (it re-examines
    // block boundaries across the whole batch and does not know about an
    // independent block's provenance, so emitting before it ran let it
    // re-merge/relabel a handful of reads based on its own, separate voting
    // logic). Applying independent blocks only after that second stitch call
    // leaves nothing downstream able to touch them again.
    bool has_proposal = false;
    PhasingChunk proposal;
};

/// Gaps whose own bridge decision never reached `joined`, together with the
/// evidence needed to consider emitting each as an independent block later,
/// once every other mechanism that could still relabel a block has run.
struct PendingIndependentGaps {
    std::vector<std::pair<PhaseGap, PhasingChunk>> proposals;
    std::unique_ptr<GapReadIndex> read_index;
    int joined_this_round = 0;
};

static void write_gap_audit_input(std::ostream& out, const PhasingChunk& proposal,
                                  const PhaseGap& gap, const GapReadIndex& index) {
    out << "SCHEMA\t1\nGAP\t" << gap.tid << '\t' << gap.left_ps << '\t'
        << gap.right_ps << '\t' << gap.left_end << '\t' << gap.right_beg
        << '\t' << proposal.ref_beg << '\t' << proposal.ref_end << '\n';
    out << "REFERENCE\t" << proposal.ref_seq << '\n';
    for (size_t vi = 0; vi < proposal.candidates.size(); ++vi) {
        const auto& v = proposal.candidates[vi];
        out << "SITE\t" << vi << '\t' << v.key.pos << '\t'
            << static_cast<int>(v.key.type) << '\t' << v.key.ref_len << '\t'
            << v.key.alt << '\t' << static_cast<int>(v.ref_base) << '\t'
            << v.lcd_var_i_to_cate << '\t' << v.msa_verified << '\t'
            << v.is_homopolymer_indel << '\t' << v.phase_set << '\t'
            << v.hap_to_cons_alle[1] << '\t' << v.hap_to_cons_alle[2];
        for (const auto& allele : v.msa_insertion_alts) out << '\t' << allele;
        out << '\n';
    }
    for (size_t ri = 0; ri < proposal.reads.size(); ++ri) {
        const auto& r = proposal.reads[ri];
        const auto found = index.assignments.find({r.input_index, r.qname});
        const auto original = found == index.assignments.end()
            ? std::make_pair(0, static_cast<hts_pos_t>(-1)) : found->second;
        out << "READ\t" << ri << '\t' << r.input_index << '\t' << r.qname
            << '\t' << r.beg << '\t' << r.end << '\t' << r.mapq << '\t'
            << r.is_skipped << '\t' << original.first << '\t' << original.second
            << '\t' << (r.alignment ? r.alignment->core.flag : -1) << '\n';
        const auto& p = proposal.read_var_profile[ri];
        for (size_t pi = 0; pi < p.alleles.size(); ++pi) {
            const int vi = p.start_var_idx + static_cast<int>(pi);
            if (vi < 0 || static_cast<size_t>(vi) >= proposal.candidates.size()) continue;
            const int graph = pi < p.graph_alleles.size() ? p.graph_alleles[pi] : -1;
            if (p.alleles[pi] < 0 && graph < 0) continue;
            const int qi = pi < p.alt_qi.size() ? p.alt_qi[pi] : -1;
            const int quality = r.alignment && qi >= 0 && qi < r.alignment->core.l_qseq
                ? bam_get_qual(r.alignment.get())[qi] : -1;
            out << "OBS\t" << ri << '\t' << vi << '\t' << p.alleles[pi]
                << '\t' << graph << '\t' << qi << '\t' << quality << '\n';
        }
    }
}

static GapRecoveryJobResult recover_one_hybrid_gap(
        std::vector<PhasingChunk>& chunks, const Options& opts,
        const std::string& contig,
        const std::unordered_map<std::string, std::string>& chrom_remap,
        const BamAuthorityIntervals* bam_authority,
        const GapReadIndex& read_index, const PhaseGap& initial_gap,
        size_t gap_index, std::mutex& chunks_mutex, const GapEvidence& frozen, bool audit = false) {
    GapRecoveryJobResult job_result;
    PhaseGap gap = initial_gap;
    RegionChunk window;
    window.tid = gap.tid;
    window.beg = std::max(gap.region_beg, gap.left_end - kGapRecoveryFlank);
    window.end = std::min(gap.region_end, gap.right_beg + kGapRecoveryFlank);
    window.chunk_id = static_cast<int>(gap_index);
    Options local_opts = opts;
    local_opts.gap_recovery_beg = gap.left_end;
    local_opts.gap_recovery_end = gap.right_beg;
    // The recovery proposal is additive and confined to this gap window.
    // Enable its graph solver for every tier, including the clean-only first
    // pass; otherwise the clean block-bridge evidence is calculated only while
    // populating the cache and is never used to decide the recovered edge.
    local_opts.private_msa_admit_all_in_region = true;
    if (!opts.phase_matrix_dump_prefix.empty()) {
        local_opts.phase_matrix_dump_prefix =
            opts.phase_matrix_dump_prefix + ".tid" + std::to_string(gap.tid) +
            ".gap" + std::to_string(gap_index);
    }
    (void)chrom_remap;
    (void)bam_authority;
    std::vector<PhasingChunk> local;
    local.push_back(frozen.project(local_opts));
    auto& proposal = local.front();
    std::ofstream audit_input, audit_votes, audit_reads;
    if (audit) {
        const auto prefix = (std::filesystem::path(opts.gap_decision_audit) /
            ("tid" + std::to_string(gap.tid) + "." + std::to_string(gap.left_end) +
             "." + std::to_string(gap.right_beg))).string();
        frozen.write_audit(prefix);
        for (const auto* suffix : {".evidence.tsv", ".votes.tsv", ".reads.tsv"})
            if (std::filesystem::exists(prefix + suffix))
                throw std::runtime_error("gap audit files exist; use a new directory: " + prefix);
        audit_input.open(prefix + ".evidence.tsv");
        audit_votes.open(prefix + ".votes.tsv");
        audit_reads.open(prefix + ".reads.tsv");
        if (!audit_input || !audit_votes || !audit_reads)
            throw std::runtime_error("cannot write gap decision audit: " + prefix);
        write_gap_audit_input(audit_input, proposal, gap, read_index);
        audit_votes << "VIEW\tTIER\tPS\tL11\tL12\tL21\tL22\tR11\tR12\tR21\tR22\tJOINED\tFLIP\n";
        audit_reads << "VIEW\tTIER\tREAD\tHP\tPS\tSKIPPED\n";
    }
    int graph_observations = 0;
    int graph_conflicts = 0;
    for (const ReadVariantProfile& profile : proposal.read_var_profile) {
        const size_t count = std::min(profile.alleles.size(),
                                      profile.graph_alleles.size());
        for (size_t i = 0; i < count; ++i) {
            if (profile.graph_alleles[i] < 0) continue;
            ++graph_observations;
            if (profile.alleles[i] >= 0 &&
                profile.alleles[i] != profile.graph_alleles[i])
                ++graph_conflicts;
        }
    }
    std::ostringstream report_rows;
    const int recovery_passes = opts.graph_gap_bam ? 2 : 1;
    // Tried preserving pass 0's fully-solved proposal here (pass 1 reprojects
    // and overwrites `proposal` in place, and its graph-observation filter in
    // select_graph_gap_bam_reads excludes any read lacking a graph-channel
    // call, even one that formed a coherent local split during pass 0), so
    // the independent-block fallback below could pick whichever pass phased
    // more reads. Measured on chr20: a modest, unclear net change in
    // independent-block yield, but it also caused 112-116 reads that were
    // fine in a fresh flag-off baseline to become entirely unphased (not
    // corrupted -- no concordant read flipped to DISCORDANT -- but lost from
    // evaluation), through an interaction with the gap-link-by-alleles
    // one-sided partial-link attachment inside this same gap's own pass-0,
    // non-orientation-only stitch_gap_proposal tier attempts (not fully
    // root-caused before running out of session budget). Reverted: keep
    // whichever proposal the pass loop naturally ends on, unconditionally, as
    // before -- see CHECKPOINT.md, 2026-09-15, for the measurements.
    for (int recovery_pass = 0; recovery_pass < recovery_passes; ++recovery_pass) {
        int selected_graph_reads = 0;
        if (recovery_pass == 1) {
            // Reproject into a scratch view first. When no read qualifies, this
            // pass has nothing to solve and pass 0's solved split is the only
            // proposal the independent-block fallback below can still use --
            // overwriting it in place discarded that split for an unphased
            // reprojection.
            PhasingChunk bam_view = frozen.project(local_opts, true);
            selected_graph_reads = select_graph_gap_bam_reads(bam_view, gap, local_opts);
            if (selected_graph_reads == 0) break;
            proposal = std::move(bam_view);
        }
        assign_hap_based_on_germline_het_vars_kmeans(
            proposal, local_opts, kCandGermlineClean);
        std::vector<uint32_t> original_flags;
        original_flags.reserve(proposal.candidates.size());
        for (const CandidateVariant& candidate : proposal.candidates)
            original_flags.push_back(candidate.lcd_var_i_to_cate);
        constexpr int kGapHomopolymerTier = 4;
        for (int tier = 1; tier <= kGapHomopolymerTier; ++tier) {
            const size_t previous_sites = proposal.candidates.size();
            if (tier == kGapHomopolymerTier) {
                // Skipping this tier in the reprojected pass hid the only join some
                // gaps have. On chr20:61,738,233-61,757,551 the audit export, which
                // bypassed the skip, showed pass 1 reaching JOINED for proposal set
                // 61725696 with left support 5 and right support 4 -- clear of
                // min_block_link_reads and of the per-haplotype orientation
                // agreement -- while pass 0's differently anchored proposal linked
                // one side only. The tier's own guards are unchanged: it still needs
                // a homopolymer candidate in the gap, --link-by-alleles, that
                // orientation agreement, and the BAM validation below.
                // Last resort means the clean evidence plus homopolymer sites, not
                // everything at once. Restoring every flag also readmitted the
                // non-homopolymer verified indels tier 3 had just failed with, and
                // those are frequently the worst sites in the window: on
                // chr20:36,247,421-36,268,291 they segregate at 0.509 and 0.644
                // against read truth while the homopolymer deletion segregates at
                // 0.923, and carrying them into this tier holds the read partition
                // at 0.846 accuracy where clean plus the homopolymer site alone
                // reaches 1.000. Escalating to a weaker class should not re-admit a
                // class that already had its turn.
                for (size_t vi = 0; vi < proposal.candidates.size(); ++vi) {
                    CandidateVariant& candidate = proposal.candidates[vi];
                    candidate.lcd_var_i_to_cate = original_flags[vi];
                    if (candidate.lcd_var_i_to_cate != kCandNoisyCandHet) continue;
                    const bool allowed = candidate.msa_verified &&
                        (candidate.key.type == VariantType::Snp ||
                         candidate.is_homopolymer_indel);
                    if (!allowed)
                        candidate.lcd_var_i_to_cate &= ~kCandGermlineVarCate;
                }
                const bool has_hp = std::any_of(
                    proposal.candidates.begin(), proposal.candidates.end(),
                    [&](const CandidateVariant& v) {
                        return v.is_homopolymer_indel &&
                               v.lcd_var_i_to_cate == kCandNoisyCandHet &&
                               v.key.sort_pos() >= gap.left_end &&
                               v.key.sort_pos() <= gap.right_beg;
                    });
                if (!has_hp || !opts.link_by_alleles) break;
                local_opts.gap_hp_link_beg = gap.left_end;
                local_opts.gap_hp_link_end = gap.right_beg;
                assign_hap_based_on_germline_het_vars_kmeans(
                    proposal, local_opts, kCandGermlineVarCate);
            } else if (tier > 1) {
                // The chromosome pass has already discovered and MSA-verified
                // these sites. Select the requested evidence tier from the cache
                // instead of rerunning consensus alignment for every gap.
                for (size_t vi = 0; vi < proposal.candidates.size(); ++vi) {
                    CandidateVariant& candidate = proposal.candidates[vi];
                    candidate.lcd_var_i_to_cate = original_flags[vi];
                    if (candidate.lcd_var_i_to_cate != kCandNoisyCandHet) continue;
                    const bool allowed = candidate.msa_verified &&
                        !candidate.is_homopolymer_indel &&
                        (tier == 3 || candidate.key.type == VariantType::Snp);
                    if (!allowed)
                        candidate.lcd_var_i_to_cate &= ~kCandGermlineVarCate;
                }
                assign_hap_based_on_germline_het_vars_kmeans(
                    proposal, local_opts, kCandGermlineVarCate);
            }
            filter_hybrid_reads_by_margin(
                local, opts.min_read_hap_margin, tier > 1);
            filter_hybrid_small_phase_sets(local, opts.min_phase_set_reads);

            GapStitchResult result;
            // Set when the homopolymer tier's own stitch joined and only the
            // BAM-only confirmation revoked it.
            bool hp_vetoed = false;
            std::vector<GapLinkEvidence> evidence;
            {
                std::unique_lock<std::mutex> lock(chunks_mutex, std::defer_lock);
                if (!audit) lock.lock();
                result = stitch_gap_proposal(
                    chunks, proposal, gap, local_opts, &read_index, true,
                    audit || recovery_pass == 1 || tier == kGapHomopolymerTier,
                    audit ? &evidence : nullptr);
            }
            // A repeat-driven graph proposal must agree with a BAM-only solve
            // before its orientation can propagate into either trusted flank.
            // Both checks are read-only; apply accepted parity edges once below.
            if (audit) {
                for (const auto& e : evidence) {
                    audit_votes << recovery_pass << '\t' << tier << '\t' << e.proposal_ps;
                    for (const auto& side : e.votes)
                        for (const int count : side) audit_votes << '\t' << count;
                    audit_votes << '\t' << result.joined << '\t' << result.right_flip << '\n';
                }
                for (size_t ri = 0; ri < proposal.reads.size(); ++ri)
                    audit_reads << recovery_pass << '\t' << tier << '\t' << ri << '\t'
                        << proposal.haps[ri] << '\t' << proposal.phase_sets[ri] << '\t'
                        << proposal.reads[ri].is_skipped << '\n';
            }
            if (!audit && tier == kGapHomopolymerTier && result.joined) {
                std::vector<PhasingChunk> validation;
                validation.push_back(frozen.project(local_opts, true));
                auto& bam_proposal = validation.front();
                const int bam_reads = select_graph_gap_bam_reads(bam_proposal, gap, local_opts);
                if (bam_reads > 0) {
                    assign_hap_based_on_germline_het_vars_kmeans(bam_proposal, local_opts, kCandGermlineClean);
                    assign_hap_based_on_germline_het_vars_kmeans(bam_proposal, local_opts, kCandGermlineVarCate);
                    // No output filters on the validation view. The check asks
                    // whether an independent BAM-only solve reaches the same
                    // orientation, and stitch_gap_proposal counts a vote only
                    // from a proposal read that still holds a hap and phase set
                    // -- so filtering this view for reporting silenced the
                    // validator itself and vetoed joins that were right.
                    // Measured on chr20:36,247,421-36,268,291: the vetoed
                    // proposal scores 1.0000 against read truth over 95 reads
                    // and matches the frozen haplotype of all 26 left-flank and
                    // 69 right-flank reads it holds.
                }
                std::unique_lock<std::mutex> lock(chunks_mutex, std::defer_lock);
                if (!audit) lock.lock();
                const auto bam_result = stitch_gap_proposal(
                    chunks, bam_proposal, gap, local_opts, &read_index, true, true);
                if (bam_reads == 0 || !bam_result.joined || bam_result.right_flip != result.right_flip) {
                    result.joined = false;
                    hp_vetoed = true;
                }
            }
            int msa_snps = 0;
            int msa_indels = 0;
            for (const auto& v : proposal.candidates) {
                if (v.lcd_var_i_to_cate != kCandNoisyCandHet) continue;
                if (v.key.type == VariantType::Snp) ++msa_snps;
                else ++msa_indels;
            }
            report_rows << contig << '\t' << initial_gap.left_end << '\t'
                        << initial_gap.right_beg << '\t' << tier << '\t'
                        << window.beg << '\t' << window.end << '\t'
                        << proposal.candidates.size() - previous_sites << '\t'
                        << msa_snps << '\t' << msa_indels << '\t'
                        << result.left_linked << '\t' << result.right_linked << '\t'
                        << result.left_link_ps << '\t' << result.right_link_ps << '\t'
                        << result.reads_added << '\t' << result.right_flip << '\t'
                        << graph_observations << '\t' << graph_conflicts << '\t'
                        << recovery_pass << '\t' << selected_graph_reads << '\t'
                        << (hp_vetoed ? "vetoed"
                                : tier == kGapHomopolymerTier && !result.joined
                                      ? "rejected"
                                : result.joined
                                      ? "joined"
                                : result.left_linked && result.right_linked
                                      ? "split"
                                : result.left_linked || result.right_linked
                                      ? "partial"
                                      : "open");
            report_rows << '\n';
            if (result.joined && !audit) {
                job_result.joined = true;
                job_result.has_edge = true;
                job_result.edge = {gap.left_ps, gap.right_ps, result.right_flip};
                break;
            }
        }
        if (job_result.joined) break;
        local_opts.gap_hp_link_beg = -1;
        local_opts.gap_hp_link_end = -1;
    }
    if (audit) {
        audit_input.flush();
        audit_votes.flush();
        audit_reads.flush();
        if (!audit_input || !audit_votes || !audit_reads)
            throw std::runtime_error("failed writing gap decision audit: " + opts.gap_decision_audit);
    }
    // The gap could not be confidently bridged to either flank (stitch never
    // reached `joined`, so nothing above was applied to `chunks` for it).
    // Its own reads may still have converged to a coherent local split among
    // themselves; keep this proposal so a later, strictly sequential pass
    // (after every gap's own stitch decision is fully settled) can emit it as
    // a new, independent block instead of discarding it -- see the race-
    // avoidance note on GapRecoveryJobResult::has_proposal. Skipped during
    // audit, which must stay strictly read-only against the frozen
    // pre-recovery state.
    if (!job_result.joined && !audit && opts.gap_independent_min_reads > 0) {
        job_result.has_proposal = true;
        job_result.proposal = std::move(proposal);
    }
    job_result.report_rows = report_rows.str();
    return job_result;
}

static void recover_hybrid_gaps(std::vector<PhasingChunk>& chunks, const Options& opts,
                                const std::string& contig,
                                const std::unordered_map<std::string, std::string>& chrom_remap,
                                const BamAuthorityIntervals* bam_authority,
                                std::ostream* report,
                                PendingIndependentGaps* pending_out = nullptr,
                                const std::vector<PreFilterLabels>* prefilter = nullptr) {
    const auto initial_gaps = find_phase_gaps(chunks);
    if (initial_gaps.empty()) return;
    const auto cache_begin = std::chrono::steady_clock::now();
    const uint64_t cache_signature =
        gap_evidence_input_signature(chunks, initial_gaps, opts);
    bool loaded_cache = false;
    if (!opts.gap_evidence_cache.empty() &&
        std::filesystem::exists(opts.gap_evidence_cache)) {
        read_gap_evidence_cache(
            opts.gap_evidence_cache, chunks, cache_signature);
        loaded_cache = true;
    } else {
        populate_gap_msa_cache(chunks, initial_gaps, opts);
        if (!opts.gap_evidence_cache.empty())
            write_gap_evidence_cache(
                opts.gap_evidence_cache, chunks, cache_signature);
    }
    const auto cache_end = std::chrono::steady_clock::now();
    GapReadIndex read_index(chunks, prefilter);
    std::vector<size_t> pending(initial_gaps.size());
    std::iota(pending.begin(), pending.end(), 0);
    std::vector<GapRecoveryJobResult> results(initial_gaps.size());
    std::mutex chunks_mutex;
    // Freeze every gap before any proposal can extend a block. These records
    // retain source observations; later tiers only construct mutable views.
    std::vector<std::unique_ptr<GapEvidence>> evidence(initial_gaps.size());
    std::atomic<size_t> next_snapshot{0};
    std::exception_ptr snapshot_error;
    std::mutex snapshot_error_mutex;
    std::vector<std::thread> snapshot_workers;
    for (size_t wi = 0; wi < std::min<size_t>(opts.threads, initial_gaps.size()); ++wi) {
        snapshot_workers.emplace_back([&]() {
            try {
                while (true) {
                    const size_t gi = next_snapshot.fetch_add(1);
                    if (gi >= initial_gaps.size()) break;
                    const auto& gap = initial_gaps[gi];
                    RegionChunk window;
                    window.tid = gap.tid;
                    window.beg = std::max(gap.region_beg, gap.left_end - kGapRecoveryFlank);
                    window.end = std::min(gap.region_end, gap.right_beg + kGapRecoveryFlank);
                    window.chunk_id = static_cast<int>(gi);
                    evidence[gi] = std::make_unique<GapEvidence>(build_cached_gap_proposal(chunks, window), gap);
                }
            } catch (...) {
                std::lock_guard<std::mutex> lock(snapshot_error_mutex);
                if (!snapshot_error) snapshot_error = std::current_exception();
            }
        });
    }
    for (auto& worker : snapshot_workers) worker.join();
    if (snapshot_error) std::rethrow_exception(snapshot_error);
    // Audit every proposal against the same pre-recovery chunks. No production
    // recovery is allowed to mutate them until all audit workers have finished.
    if (!opts.gap_decision_audit.empty()) {
        std::filesystem::create_directories(opts.gap_decision_audit);
        const auto manifest_path = std::filesystem::path(opts.gap_decision_audit) /
            ("tid" + std::to_string(initial_gaps.front().tid) + ".manifest.tsv");
        if (std::filesystem::exists(manifest_path))
            throw std::runtime_error("gap audit snapshot exists; use a new directory: " + manifest_path.string());
        std::ofstream manifest(manifest_path);
        manifest << "SCHEMA\t1\nCONTIG\t" << contig << "\nCACHE_SIGNATURE\t" << cache_signature
                 << "\nGAPS\t" << initial_gaps.size() << "\nREFERENCE\t" << opts.ref_fasta
                 << "\nBAM\t" << opts.primary_bam_file() << "\nGRAPH_SITES\t" << opts.graph_sites_vcf
                 << "\nMIN_BLOCK_LINK_READS\t" << opts.min_block_link_reads
                 << "\nSTITCH_RULE\t" << opts.stitch_rule << "\nSTITCH_MARGIN\t" << opts.stitch_min_margin
                 << "\nMIN_READ_MARGIN\t" << opts.min_read_hap_margin
                 << "\nGRAPH_BAM\t" << opts.graph_gap_bam << '\n';
        manifest.close();
        if (!manifest) throw std::runtime_error("failed writing gap audit manifest: " + manifest_path.string());
        const auto block_path = std::filesystem::path(opts.gap_decision_audit) /
            ("tid" + std::to_string(initial_gaps.front().tid) + ".blocks.tsv");
        const auto member_path = std::filesystem::path(opts.gap_decision_audit) /
            ("tid" + std::to_string(initial_gaps.front().tid) + ".members.tsv");
        if (std::filesystem::exists(block_path) || std::filesystem::exists(member_path))
            throw std::runtime_error("gap block snapshot exists; use a new directory: " + opts.gap_decision_audit);
        std::ofstream blocks_out(block_path), members_out(member_path);
        members_out << "INPUT\tREAD\tHP\tPS\n";
        std::vector<GapReadIndex::Key> members;
        for (const auto& entry : read_index.assignments) members.push_back(entry.first);
        std::sort(members.begin(), members.end());
        std::map<hts_pos_t, std::array<int, 2>> block_counts;
        for (const auto& key : members) {
            const auto assignment = read_index.assignments.at(key);
            members_out << key.first << '\t' << key.second << '\t'
                        << assignment.first << '\t' << assignment.second << '\n';
            ++block_counts[assignment.second][assignment.first - 1];
        }
        std::map<hts_pos_t, std::pair<hts_pos_t, hts_pos_t>> block_bounds;
        for (const auto& chunk : chunks)
            for (const auto& v : chunk.candidates) {
                if (!block_counts.count(v.phase_set) || v.hap_to_cons_alle[1] < 0 ||
                    v.hap_to_cons_alle[2] < 0 || v.hap_to_cons_alle[1] == v.hap_to_cons_alle[2]) continue;
                const auto pos = v.key.sort_pos();
                auto inserted = block_bounds.emplace(v.phase_set, std::make_pair(pos, pos));
                inserted.first->second.first = std::min(inserted.first->second.first, pos);
                inserted.first->second.second = std::max(inserted.first->second.second, pos);
            }
        blocks_out << "PS\tBEGIN\tEND\tHP1_READS\tHP2_READS\n";
        for (const auto& [ps, counts] : block_counts) {
            const auto found = block_bounds.find(ps);
            const auto bounds = found == block_bounds.end()
                ? std::make_pair(static_cast<hts_pos_t>(-1), static_cast<hts_pos_t>(-1)) : found->second;
            blocks_out << ps << '\t' << bounds.first << '\t' << bounds.second
                       << '\t' << counts[0] << '\t' << counts[1] << '\n';
        }
        blocks_out.close();
        members_out.close();
        if (!blocks_out || !members_out)
            throw std::runtime_error("failed writing gap block snapshot: " + opts.gap_decision_audit);
        std::atomic<size_t> next_audit{0};
        std::exception_ptr audit_error;
        std::mutex audit_error_mutex;
        std::vector<std::thread> audit_workers;
        const size_t count = std::min<size_t>(opts.threads, initial_gaps.size());
        for (size_t wi = 0; wi < count; ++wi) {
            audit_workers.emplace_back([&]() {
                try {
                    while (true) {
                        const size_t gi = next_audit.fetch_add(1);
                        if (gi >= initial_gaps.size()) break;
                        recover_one_hybrid_gap(chunks, opts, contig, chrom_remap,
                            bam_authority, read_index, initial_gaps[gi], gi, chunks_mutex, *evidence[gi], true);
                    }
                } catch (...) {
                    std::lock_guard<std::mutex> lock(audit_error_mutex);
                    if (!audit_error) audit_error = std::current_exception();
                }
            });
        }
        for (auto& worker : audit_workers) worker.join();
        if (audit_error) std::rethrow_exception(audit_error);
    }
    size_t wave_count = 0;
    while (!pending.empty()) {
        std::vector<size_t> wave;
        std::vector<size_t> deferred;
        for (const size_t gap_index : pending) {
            const bool conflicts = std::any_of(
                wave.begin(), wave.end(), [&](const size_t selected_index) {
                    return gap_recovery_jobs_conflict(
                        initial_gaps[gap_index], initial_gaps[selected_index]);
                });
            if (conflicts) deferred.push_back(gap_index);
            else wave.push_back(gap_index);
        }
        pending = std::move(deferred);
        ++wave_count;

        const size_t worker_count = std::min<size_t>(
            static_cast<size_t>(opts.threads), wave.size());
        std::atomic<size_t> next_job{0};
        std::exception_ptr first_error;
        std::mutex error_mutex;
        std::vector<std::thread> workers;
        workers.reserve(worker_count);
        for (size_t worker_i = 0; worker_i < worker_count; ++worker_i) {
            workers.emplace_back([&]() {
                try {
                    while (true) {
                        const size_t job = next_job.fetch_add(1);
                        if (job >= wave.size()) break;
                        const size_t gap_index = wave[job];
                        results[gap_index] = recover_one_hybrid_gap(
                            chunks, opts, contig, chrom_remap, bam_authority,
                            read_index, initial_gaps[gap_index], gap_index,
                            chunks_mutex, *evidence[gap_index]);
                    }
                } catch (...) {
                    std::lock_guard<std::mutex> lock(error_mutex);
                    if (!first_error) first_error = std::current_exception();
                }
            });
        }
        for (std::thread& worker : workers) worker.join();
        if (first_error) std::rethrow_exception(first_error);
    }
    const auto recovery_end = std::chrono::steady_clock::now();
    int joined = 0;
    std::vector<GapPhaseEdge> edges;
    for (const auto& result : results) {
        joined += result.joined;
        if (result.has_edge) edges.push_back(result.edge);
        if (report) *report << result.report_rows;
    }
    const int edge_conflicts = apply_gap_phase_edges(chunks, edges);
    const int anchored_sites = anchor_orphan_msa_sites(chunks, opts, &initial_gaps);
    std::cerr << "Gap recovery: " << initial_gaps.size() << " initial gaps, "
              << joined << " joined in " << wave_count << " wave(s) on "
              << contig << ", " << edge_conflicts
              << " conflicting edge(s) rejected\n";
    std::cerr << "Gap recovery: attached " << anchored_sites
              << " orphan MSA sites to read-supported blocks\n";
    if (pending_out != nullptr) {
        // The caller applies these strictly after the batch's second
        // stitch_chunk_haps call; read_index must outlive that call, and it
        // borrows read names from `chunks`, which the caller keeps alive
        // (batch.chunks) across both calls.
        pending_out->joined_this_round = joined;
        pending_out->read_index = std::make_unique<GapReadIndex>(std::move(read_index));
        for (size_t gi = 0; gi < results.size(); ++gi) {
            if (!results[gi].has_proposal) continue;
            pending_out->proposals.emplace_back(initial_gaps[gi], std::move(results[gi].proposal));
        }
    }
    std::cerr << "Gap recovery timing: evidence cache "
              << (loaded_cache ? "load " : "build ")
              << std::chrono::duration<double>(cache_end - cache_begin).count()
              << " s, cached gap solves "
              << std::chrono::duration<double>(recovery_end - cache_end).count()
              << " s\n";
}

void run_collect_hybrid_variation(const Options& opts) {
    if (opts.graph_sites_vcf.empty())
        throw std::runtime_error("--graph-sites required for hybrid mode");
    if (opts.gaf_file.empty())
        throw std::runtime_error("--gaf required for hybrid mode");

    std::unique_ptr<PgbamSidecarData> pgbam_sidecar;
    if (!opts.pgbam_file.empty())
        pgbam_sidecar = std::make_unique<PgbamSidecarData>(
            load_pgbam_sidecar(opts.pgbam_file));

    std::unique_ptr<faidx_t, FaiDeleter> fai(load_reference_index(opts.ref_fasta));

    const std::vector<RegionChunk> chunks = load_region_chunks(opts);
    SamFile bam(opts.primary_bam_file(), 1, opts.ref_fasta);
    std::unique_ptr<bam_hdr_t, HeaderDeleter> header(sam_hdr_read(bam.get()));
    if (!header) throw std::runtime_error("failed to read BAM header");
    ReferenceCache ref(fai.get());

    VariantKeySet private_keys;
    const VariantKeySet* private_keys_ptr = nullptr;
    if (!opts.private_sites_vcf.empty()) {
        private_keys = load_private_variant_keys(opts.private_sites_vcf, header.get());
        private_keys_ptr = &private_keys;
        std::cerr << "Loaded " << private_keys.size()
                  << " private BAM candidate key(s) from "
                  << opts.private_sites_vcf << "\n";
    }

    BamAuthorityIntervals bam_authority;
    const BamAuthorityIntervals* bam_authority_ptr = nullptr;
    if (!opts.bam_authoritative_bed.empty()) {
        bam_authority = load_bam_authority_intervals(
            opts.bam_authoritative_bed, header.get());
        bam_authority_ptr = &bam_authority;
        size_t interval_count = 0;
        for (const auto& [tid, rows] : bam_authority) {
            (void)tid;
            interval_count += rows.size();
        }
        std::cerr << "Loaded " << interval_count
                  << " merged BAM-authoritative interval(s) from "
                  << opts.bam_authoritative_bed << "\n";
    }

    // Contig name resolution between BAM (e.g. "chr20") and pangenome paths
    // (e.g. "GRCh38#0#chr20") is handled by suffix matching in the tabix
    // query layers (append_graph_sites_tabix_filtered for VCF,
    // tbx_seq_tid_with_pangenome_fallback for GAF).  The VCF parser also
    // extracts ref_contig as the suffix after '#', so site coordinates use
    // BAM-compatible contig names.  No explicit remap table is needed.
    std::unordered_map<std::string, std::string> chrom_remap;

    // Open output files (same as BAM pipeline).
    std::ofstream variant_out(opts.output_tsv);
    if (!variant_out)
        throw std::runtime_error("failed to open output: " + opts.output_tsv);
    write_variants_tsv_header(variant_out);

    std::ofstream vcf_out;
    if (!opts.output_vcf.empty()) {
        vcf_out.open(opts.output_vcf);
        if (!vcf_out)
            throw std::runtime_error("failed to open VCF output: " + opts.output_vcf);
        write_variants_vcf_header(vcf_out, opts, header.get());
    }
    std::ofstream phased_vcf_out;
    if (!opts.output_phased_vcf.empty()) {
        phased_vcf_out.open(opts.output_phased_vcf);
        if (!phased_vcf_out)
            throw std::runtime_error("failed to open phased VCF: " + opts.output_phased_vcf);
        write_phased_variants_vcf_header(phased_vcf_out, opts, header.get());
    }

    std::unique_ptr<PhasedAlignmentWriter> phased_aln_writer;
    if (!opts.output_aln.empty())
        phased_aln_writer = std::make_unique<PhasedAlignmentWriter>(opts, header.get());

    std::ofstream recovery_report;
    if (!opts.gap_recovery_report.empty()) {
        recovery_report.open(opts.gap_recovery_report);
        if (!recovery_report) throw std::runtime_error("failed to open recovery report: " + opts.gap_recovery_report);
        recovery_report << "CHROM\tGAP_LEFT\tGAP_RIGHT\tTIER\tWINDOW_BEG\tWINDOW_END\tNEW_SITES\tMSA_HET_SNPS\tMSA_HET_INDELS\tLEFT_LINK\tRIGHT_LINK\tLEFT_LINK_PS\tRIGHT_LINK_PS\tREADS_ADDED\tRECOVERY_FLIP\tGRAPH_OBSERVATIONS\tGRAPH_CONFLICTS\tGRAPH_BAM_PASS\tSELECTED_GRAPH_READS\tSTATUS\n";
    }

    size_t n_variants = 0;
    size_t n_out_aln_reads = 0;
    size_t batch_begin = 0;
    while (batch_begin < chunks.size()) {
        size_t batch_end = batch_begin + 1;
        while (batch_end < chunks.size() &&
               chunks[batch_end].reg_chunk_i == chunks[batch_begin].reg_chunk_i) {
            ++batch_end;
        }

        // Determine the graph query contig for this batch.
        const std::string batch_contig =
            header->target_name[chunks[batch_begin].tid];
        // Suffix matching in the tabix query layers resolves BAM contig
        // names (e.g. "chr20") against pangenome paths in the VCF/GAF index.
        const std::string& graph_query_contig = batch_contig;

        ChunkBatchResult batch = collect_hybrid_chunk_batch_parallel(
            opts, chunks, batch_begin, batch_end,
            graph_query_contig, chrom_remap, private_keys_ptr,
            bam_authority_ptr);
        stitch_chunk_haps(batch.chunks, &opts, pgbam_sidecar.get());
        // Snapshot the labels the solve produced, before the output filters
        // erase the ones gap recovery's link votes are counted from. The gap
        // inventory itself is still derived from the filtered chunks, so this
        // changes what the votes can see and nothing about what is emitted.
        std::vector<PreFilterLabels> prefilter_labels;
        if (opts.recover_gaps) {
            prefilter_labels.reserve(batch.chunks.size());
            for (const PhasingChunk& chunk : batch.chunks)
                prefilter_labels.push_back({chunk.haps, chunk.phase_sets});
        }
        filter_hybrid_reads_by_margin(batch.chunks, opts.min_read_hap_margin,
                                      opts.private_msa_admit_all_in_region);
        filter_hybrid_small_phase_sets(batch.chunks, opts.min_phase_set_reads);
        if (opts.recover_gaps) {
            // Emitting an independent block, or bridging a gap, can both
            // create new, smaller gaps that did not exist in the original
            // inventory -- e.g. one newly-phased block now sits between an
            // existing flank and where a bridge previously had nothing to
            // reach. find_phase_gaps (inside recover_hybrid_gaps) discovers
            // the current gap inventory fresh each call, so repeating the
            // whole sequence lets the SAME, already-validated bridge and
            // independent-block logic reach those new opportunities too,
            // rather than requiring a separate general stitch solver. Bounded
            // and terminated on the first round with no progress at all.
            for (int round = 0; round < opts.gap_recovery_max_rounds; ++round) {
                // A cache is keyed by a signature over the current gap
                // inventory; every round after the first changes that
                // inventory (a newly-independent block splits an old gap
                // into two smaller ones), which read_gap_evidence_cache
                // rejects outright as an input mismatch. Only round 0 can
                // legitimately reuse a cache built for the original gaps.
                Options round_opts = opts;
                if (round > 0) round_opts.gap_evidence_cache.clear();
                PendingIndependentGaps pending;
                recover_hybrid_gaps(batch.chunks, round_opts, graph_query_contig, chrom_remap,
                                    bam_authority_ptr, recovery_report.is_open() ? &recovery_report : nullptr,
                                    &pending, &prefilter_labels);
                stitch_chunk_haps(batch.chunks, &opts, pgbam_sidecar.get());
                // Only after this second stitch call -- which re-examines
                // block boundaries across the whole batch and would
                // otherwise be free to re-merge or relabel a block it does
                // not know is independent -- can a gap-only local split
                // safely be applied as new, untouchable output. See
                // PendingIndependentGaps.
                int independent_gaps = 0, independent_reads = 0;
                // Index (not pointer) into pending.proposals, so the matching
                // proposal can be looked up cleanly for the bridge attempt
                // below without pending.proposals being touched in between.
                std::vector<std::pair<size_t, hts_pos_t>> emitted;
                // Adjacent gaps' windows overlap by construction, so two
                // gaps in this same sequential loop can independently
                // reconstruct the same leftmost-het phase-set id from
                // overlapping gap-only read pools. `pending.read_index` is a
                // whole-batch snapshot frozen before this loop starts, so it
                // cannot see an id a prior iteration of *this* loop already
                // claimed -- track that separately and thread it through.
                std::set<hts_pos_t> emitted_this_round;
                for (size_t pi = 0; pi < pending.proposals.size(); ++pi) {
                    auto& [gap, proposal] = pending.proposals[pi];
                    hts_pos_t chosen_ps = -1;
                    const int gained = emit_independent_gap_block(
                        batch.chunks, proposal, gap, opts, *pending.read_index,
                        opts.gap_independent_min_reads, &chosen_ps, &emitted_this_round);
                    if (gained > 0) {
                        ++independent_gaps;
                        independent_reads += gained;
                        emitted.emplace_back(pi, chosen_ps);
                    }
                }
                if (independent_gaps > 0 || opts.gap_independent_min_reads > 0)
                    std::cerr << "Gap recovery: " << independent_gaps
                              << " gap(s) emitted as new independent blocks, "
                              << independent_reads << " previously-unphased read(s) phased\n";
                // Cheap final stitch: try bridging each new block to its own
                // two flanks as two small sub-gaps, reusing the SAME proposal
                // (the gap's own local evidence already spans both flanks --
                // no re-extraction needed) against a read index rebuilt from
                // the just-updated chunks, so it correctly sees the new
                // block's reads as assigned to its own phase set. Bounded,
                // cheap alternative to a full extra recovery round: it can
                // only ever connect a block just emitted this round, so one
                // pass suffices and it never risks the unbounded, uncached
                // re-extraction cost of round > 0.
                int stitched_independent = 0;
                if (opts.gap_bridge_independent_blocks && !emitted.empty()) {
                    const GapReadIndex fresh_index(batch.chunks);
                    std::vector<GapPhaseEdge> independent_edges;
                    for (const auto& [pi, chosen_ps] : emitted) {
                        const PhaseGap& gap = pending.proposals[pi].first;
                        const PhasingChunk& proposal = pending.proposals[pi].second;
                        for (const bool right_side : {false, true}) {
                            PhaseGap sub = gap;
                            if (right_side) sub.left_ps = chosen_ps;
                            else sub.right_ps = chosen_ps;
                            // Deferred: a chain of gap.left_ps <-> chosen_ps
                            // <-> gap.right_ps must be resolved together by
                            // one union-find pass below, not by two
                            // independent immediate renames -- applying the
                            // left edge's rename first would leave the right
                            // edge relabelling gap.right_ps's reads to a
                            // chosen_ps that no longer names anything,
                            // splitting what should be one merged block into
                            // two surviving phase-set ids instead of one.
                            const auto result = stitch_gap_proposal(
                                batch.chunks, proposal, sub, opts, &fresh_index,
                                /*defer_phase_set_merge=*/true);
                            if (result.joined) {
                                independent_edges.push_back(
                                    {sub.left_ps, sub.right_ps, result.right_flip});
                            }
                        }
                    }
                    stitched_independent = static_cast<int>(independent_edges.size());
                    apply_gap_phase_edges(batch.chunks, independent_edges);
                }
                if (stitched_independent > 0)
                    std::cerr << "Gap recovery: " << stitched_independent
                              << " newly-independent block/flank edge(s) bridged\n";
                if (pending.joined_this_round == 0 && independent_gaps == 0 &&
                    stitched_independent == 0) break;
            }
            for (auto& chunk : batch.chunks) prune_not_candidate_variants(chunk);
        }
        CandidateTable variants = merge_chunk_candidates(batch.chunks);
        n_variants += variants.size();
        write_variants_tsv_records(variant_out, header.get(), ref, variants);
        if (!opts.output_vcf.empty())
            write_variants_vcf_records(vcf_out, opts, header.get(), ref, variants);
        if (!opts.output_phased_vcf.empty())
            write_phased_variants_vcf_records(
                phased_vcf_out, opts, header.get(), ref, variants);
        if (phased_aln_writer)
            n_out_aln_reads += static_cast<size_t>(
                phased_aln_writer->write_chunks(batch.chunks));

        batch_begin = batch_end;
    }

    std::cerr << "Hybrid: processed " << chunks.size() << " region chunks with "
              << opts.threads << " worker thread(s)\n";
    std::cerr << "Collected " << n_variants << " candidate variant sites into "
              << opts.output_tsv << "\n";
    if (!opts.output_vcf.empty())
        std::cerr << "Wrote candidate VCF to " << opts.output_vcf << "\n";
    if (!opts.output_phased_vcf.empty())
        std::cerr << "Wrote phased candidate VCF to " << opts.output_phased_vcf << "\n";
    if (!opts.output_aln.empty()) {
        phased_aln_writer.reset();
        std::cerr << "Output " << n_out_aln_reads << " reads to phased alignment\n";
        if (opts.refine_aln)
            coordinate_sort_refined_alignment_file_or_throw(opts);
        if (opts.output_aln_format == OutputAlignmentFormat::Bam ||
            opts.output_aln_format == OutputAlignmentFormat::Cram) {
            if (sam_index_build(opts.output_aln.c_str(), 0) != 0)
                throw std::runtime_error(
                    "failed to index output alignment: " + opts.output_aln);
        }
    }
}

} // namespace pgphase_collect

// ════════════════════════════════════════════════════════════════════════════
// CLI entry point
// ════════════════════════════════════════════════════════════════════════════

namespace pgphase_collect {

/**
 * @brief `getopt_long` option codes for flags without short aliases (values ≥ 1000).
 */
enum LongOption {
    kMinAltDepthOption = 1000,
    kMinAfOption,
    kMaxAfOption,
    kNoisyRegMergeDisOption,
    kMinSvLenOption,
    kChunkSizeOption,
    kHifiOption,
    kOntOption,
    kShortReadsOption,
    kStrandBiasPvalOption,
    kNoisyMaxXgapsOption,
    kMaxNoisyFracOption,
    kNoisySlideWinOption,
    kDebugSiteOption,
    kInputIsListOption,
    kRefineAlnOption,
    kPhasedVcfOutputOption,
    kPgbamFileOption,
    kPgbamPrimaryMarginOption,
    kPgbamPrimaryMinWinningOption,
    kNoPgbamCleanupPassOption,
    kPgbamCleanupMarginOption,
    kPgbamCleanupMinWinningOption,
    kNoPgbamRelaxedCleanupPassOption,
    kPgbamRelaxedCleanupMarginOption,
    kPgbamRelaxedCleanupMinWinningOption,
    kAmbBaseOption,
    kRefOption,
    kBamOption,
    kPhaseMatrixDumpOption,

};

/**
 * @brief Reads a newline-separated list of BAM/CRAM paths.
 *
 * Strips surrounding whitespace, skips blank lines and `#` comments. Requires at least one path.
 *
 * @param path Text file passed with `--input-is-list` / `-L`.
 * @return Non-empty list of alignment file paths.
 */
static std::vector<std::string> load_bam_list(const std::string& path) {
    std::ifstream in(path);
    if (!in) throw std::runtime_error("failed to open BAM/CRAM list: " + path);
    std::vector<std::string> files;
    std::string line;
    while (std::getline(in, line)) {
        const size_t first = line.find_first_not_of(" \t\r\n");
        if (first == std::string::npos) continue;
        if (line[first] == '#') continue;
        const size_t last = line.find_last_not_of(" \t\r\n");
        files.push_back(line.substr(first, last - first + 1));
    }
    if (files.empty()) throw std::runtime_error("BAM/CRAM list is empty: " + path);
    return files;
}

/**
 * @brief Prints usage and option summary for `collect-bam-variation` to stdout.
 */
static void print_collect_help() {
    std::cout
        << "Usage: pgphase collect-bam-variation [options]\n"
        << "\n"
        << "Required:\n"
        << "      --ref FILE                Reference FASTA (indexed)\n"
        << "      --bam FILE                Input BAM/CRAM file (or list with -L)\n"
        << "\n"
        << "Options:\n"
        << "  -L, --input-is-list          Treat --bam path as a list of BAM/CRAM files\n"
        << "  -X, --extra-bam FILE         Extra input BAM/CRAM file; may be repeated\n"
        << "  -t, --threads INT             Region worker threads [1]\n"
        << "  -q, --min-mapq INT            Minimum read mapping quality [30]\n"
        << "  -B, --min-bq INT              Minimum base quality for candidate sites [10]\n"
        << "  -D, --min-depth INT           Minimum total depth for clean candidates [5]\n"
        << "      --min-alt-depth INT       Minimum alternate depth for clean candidates [2]\n"
        << "      --min-af FLOAT            Minimum allele fraction for clean het candidates [0.20]\n"
        << "      --max-af FLOAT            Maximum allele fraction for clean het candidates [0.80]\n"
        << "  -r, --region STR              Optional region; may be repeated\n"
        << "      --region-file FILE        BED file of regions to process\n"
        << "      --autosome                Process chr1-22 / 1-22 only\n"
        << "  -j, --max-var-ratio FLOAT     Skip reads above this variant/ref-span ratio [0.05]\n"
        << "      --max-noisy-frac FLOAT    Skip reads with > this fraction in noisy regions [0.5]\n"
        << "      --include-filtered        Include QC-fail and duplicate reads\n"
        << "      --amb-base                Emit VCF rows with ambiguous (non-ACGT) REF/ALT bases\n"
        << "  -o, --output FILE             Output TSV file [output.tsv]\n"
        << "  -v, --vcf-output FILE         Optional VCF output for collected candidates\n"
        << "      --phased-vcf-out FILE      Optional phased VCF (GT:DP:AD:VAF:GQ:PS)\n"
        << "  -S/b/C --out-sam/bam/cram FILE\n"
        << "                                output phased SAM/BAM/CRAM file []\n"
        << "                                note: multiple input BAM/CRAM files will be merged in SAM/BAM/CRAM output\n"
        << "      --refine-aln              refine alignment in SAM/BAM/CRAM output;\n"
        << "                                coordinate-sorts with samtools sort before BAM/CRAM indexing (samtools on PATH)\n"

        << "      --pgbam-file FILE         Optional .pgbam sidecar for fallback chunk stitching when common-read signal is absent\n"
        << "      --pgbam-primary-margin INT         Thread polarity margin for primary .pgbam stitching [2]\n"
        << "      --pgbam-primary-min-winning INT    Winning shared polarized threads for primary .pgbam stitching [2]\n"
        << "      --no-pgbam-cleanup-pass            Disable final .pgbam cleanup pass with margin/min [2/1]\n"
        << "      --pgbam-cleanup-margin INT         Thread polarity margin for final cleanup pass [2]\n"
        << "      --pgbam-cleanup-min-winning INT    Winning shared polarized threads for final cleanup pass [1]\n"
        << "      --no-pgbam-relaxed-cleanup-pass    Disable relaxed .pgbam cleanup pass with margin/min [1/1]\n"
        << "      --pgbam-relaxed-cleanup-margin INT Thread polarity margin for relaxed cleanup pass [1]\n"
        << "      --pgbam-relaxed-cleanup-min-winning INT\n"
        << "                                      Winning shared polarized threads for relaxed cleanup pass [1]\n"

        << "      --chunk-size INT          Region chunk size in bp [500000]\n"
        << "      --noisy-merge-dis INT     Max distance (bp) to merge noisy/SV windows [500]\n"
        << "      --min-sv-len INT          min_sv_len for noisy-region cgranges merge [30]\n"
        << "      --noisy-slide-win INT     Slide window (bp) for per-read noisy regions [HiFi 100 / ONT/short reads 25]\n"
        << "      --hifi                    HiFi mode: 100 bp noisy window [default; no ONT Fisher strand test]\n"
        << "      --ont                     ONT mode: 25 bp window + Fisher exact test for alt strand bias\n"
        << "      --short-reads             Short-read mode: 25 bp noisy window (no ONT Fisher strand test)\n"
        << "      --strand-bias-pval FLOAT  max p-value for ONT strand filter [0.01]\n"
        << "      --noisy-max-xgaps INT     max indel len (bp) for STR/homopolymer flags [5]\n"
        << "  -V, --verbose INT            Verbosity level; 2 prints noisy-region diagnostic logs [0]\n"
        << "\n"
        << "Examples:\n"
        << "  pgphase collect-bam-variation \\\n"
        << "      --ref ref.fa \\\n"
        << "      --bam hifi.bam \\\n"
        << "      -o candidates.tsv \\\n"
        << "      --phased-vcf-out phased.vcf \\\n"
        << "      -t 8 \\\n"
        << "      -r chr11:1000-2000 \\\n"
        << "      -r chr12:1-500\n"
        << "\n"
        << "  pgphase collect-bam-variation \\\n"
        << "      --ref ref.fa \\\n"
        << "      --bam ont_reads.bam \\\n"
        << "      --ont \\\n"
        << "      --phased-vcf-out phased.vcf \\\n"
        << "      -t 16\n";
}

} // namespace pgphase_collect

/**
 * @brief CLI entry for the `collect-bam-variation` subcommand.
 *
 * Parses GNU long options into `Options`, validates numeric thresholds, resolves the BAM list,
 * then calls `run_collect_bam_variation`. Expects `argv` with the `collect-bam-variation` token
 * already removed by the caller.
 *
 * @param argc Argument count.
 * @param argv Argument vector (reference FASTA, input BAM or list, optional region strings).
 * @return 0 on success, 1 on usage or validation error, 1 if `run_collect_bam_variation` throws.
 */
int collect_bam_variation(int argc, char* argv[]) {
    using namespace pgphase_collect;
    Options opts;
    {
        std::ostringstream cmd;
        cmd << "pgphase collect-bam-variation";
        for (int i = 1; i < argc; ++i) {
            cmd << ' ' << argv[i];
        }
        opts.command_line = cmd.str();
    }
    std::vector<std::string> extra_bam_files;
    optind = 1;

    const struct option long_options[] = {
        {"threads",                   required_argument, nullptr, 't'},
        {"min-mapq",                  required_argument, nullptr, 'q'},
        {"min-bq",                    required_argument, nullptr, 'B'},
        {"min-depth",                 required_argument, nullptr, 'D'},
        {"min-alt-depth",             required_argument, nullptr, kMinAltDepthOption},
        {"min-af",                    required_argument, nullptr, kMinAfOption},
        {"max-af",                    required_argument, nullptr, kMaxAfOption},
        {"region",                    required_argument, nullptr, 'r'},
        {"region-file",               required_argument, nullptr, 'R'},
        {"autosome",                  no_argument,       nullptr, 'a'},
        {"max-var-ratio",             required_argument, nullptr, 'j'},
        {"max-noisy-frac",            required_argument, nullptr, kMaxNoisyFracOption},
        {"include-filtered",          no_argument,       nullptr, 'f'},
        {"amb-base",                  no_argument,       nullptr, kAmbBaseOption},
        {"output",                    required_argument, nullptr, 'o'},
        {"vcf-output",                required_argument, nullptr, 'v'},
        {"phased-vcf-out",            required_argument, nullptr, kPhasedVcfOutputOption},
        {"out-sam",                   required_argument, nullptr, 'S'},
        {"out-bam",                   required_argument, nullptr, 'b'},
        {"out-cram",                  required_argument, nullptr, 'C'},
        {"refine-aln",                no_argument,       nullptr, kRefineAlnOption},

        {"pgbam-file",                required_argument, nullptr, kPgbamFileOption},
        {"pgbam-primary-margin",      required_argument, nullptr, kPgbamPrimaryMarginOption},
        {"pgbam-primary-min-winning", required_argument, nullptr, kPgbamPrimaryMinWinningOption},
        {"no-pgbam-cleanup-pass",     no_argument,       nullptr, kNoPgbamCleanupPassOption},
        {"pgbam-cleanup-margin",      required_argument, nullptr, kPgbamCleanupMarginOption},
        {"pgbam-cleanup-min-winning", required_argument, nullptr, kPgbamCleanupMinWinningOption},
        {"no-pgbam-relaxed-cleanup-pass", no_argument,   nullptr, kNoPgbamRelaxedCleanupPassOption},
        {"pgbam-relaxed-cleanup-margin", required_argument, nullptr, kPgbamRelaxedCleanupMarginOption},
        {"pgbam-relaxed-cleanup-min-winning", required_argument, nullptr, kPgbamRelaxedCleanupMinWinningOption},

        {"chunk-size",                required_argument, nullptr, kChunkSizeOption},
        {"noisy-merge-dis",           required_argument, nullptr, kNoisyRegMergeDisOption},
        {"min-sv-len",                required_argument, nullptr, kMinSvLenOption},
        {"noisy-slide-win",           required_argument, nullptr, kNoisySlideWinOption},
        {"debug-site",                required_argument, nullptr, kDebugSiteOption},
        {"dump-phase-matrix",         required_argument, nullptr, kPhaseMatrixDumpOption},
        {"extra-bam",                 required_argument, nullptr, 'X'},
        {"input-is-list",             no_argument,       nullptr, 'L'},
        {"hifi",                      no_argument,       nullptr, kHifiOption},
        {"ont",                       no_argument,       nullptr, kOntOption},
        {"short-reads",               no_argument,       nullptr, kShortReadsOption},
        {"strand-bias-pval",          required_argument, nullptr, kStrandBiasPvalOption},
        {"noisy-max-xgaps",           required_argument, nullptr, kNoisyMaxXgapsOption},
        {"ref",                       required_argument, nullptr, kRefOption},
        {"bam",                       required_argument, nullptr, kBamOption},
        {"verbose",                  required_argument, nullptr, 'V'},
        {"help",                      no_argument,       nullptr, 'h'},
        {nullptr, 0, nullptr, 0}
    };

    int opt = 0;
    int long_index = 0;
    bool read_technology_was_set = false;
    bool read_technology_conflict = false;
    /** Records technology mode; flags conflict if more than one of --hifi/--ont/--short-reads is set. */
    const auto set_read_technology = [&](ReadTechnology tech) {
        if (read_technology_was_set && opts.read_technology != tech) {
            read_technology_conflict = true;
            return;
        }
        opts.read_technology = tech;
        read_technology_was_set = true;
    };

    while ((opt = getopt_long(argc, argv, "t:q:B:D:r:R:aj:o:v:S:b:C:hX:LV:", long_options, &long_index)) != -1) {
        switch (opt) {
            case 't': opts.threads = parse_int_arg(optarg, "--threads"); break;
            case 'q': opts.min_mapq = parse_int_arg(optarg, "--min-mapq"); break;
            case 'B': opts.min_bq = parse_int_arg(optarg, "--min-bq"); break;
            case 'D': opts.min_depth = parse_int_arg(optarg, "--min-depth"); break;
            case kMinAltDepthOption:    opts.min_alt_depth = parse_int_arg(optarg, "--min-alt-depth"); break;
            case kMinAfOption:          opts.min_af = parse_double_arg(optarg, "--min-af"); break;
            case kMaxAfOption:          opts.max_af = parse_double_arg(optarg, "--max-af"); break;
            case 'r': opts.regions.push_back(optarg); break;
            case 'R': opts.region_file = optarg; break;
            case 'a': opts.autosome = true; break;
            case 'j': opts.max_var_ratio_per_read = parse_double_arg(optarg, "--max-var-ratio"); break;
            case kMaxNoisyFracOption:   opts.max_noisy_frac_per_read = parse_double_arg(optarg, "--max-noisy-frac"); break;
            case 'f': opts.include_filtered = true; break;
            case kAmbBaseOption: opts.output_ambiguous_bases = true; break;
            case 'o': opts.output_tsv = optarg; break;
            case 'v': opts.output_vcf = optarg; break;
            case 'S':
                opts.output_aln = optarg;
                opts.output_aln_format = OutputAlignmentFormat::Sam;
                break;
            case 'b':
                opts.output_aln = optarg;
                opts.output_aln_format = OutputAlignmentFormat::Bam;
                break;
            case 'C':
                opts.output_aln = optarg;
                opts.output_aln_format = OutputAlignmentFormat::Cram;
                break;
            case kPhasedVcfOutputOption: opts.output_phased_vcf = optarg; break;
            case kRefineAlnOption:      opts.refine_aln = true; break;

            case kPgbamFileOption:      opts.pgbam_file = optarg; break;
            case kPgbamPrimaryMarginOption: opts.pgbam_primary_polarity_margin = parse_int_arg(optarg, "--pgbam-primary-margin"); break;
            case kPgbamPrimaryMinWinningOption: opts.pgbam_primary_min_winning_threads = parse_int_arg(optarg, "--pgbam-primary-min-winning"); break;
            case kNoPgbamCleanupPassOption: opts.pgbam_cleanup_pass = false; break;
            case kPgbamCleanupMarginOption: opts.pgbam_cleanup_polarity_margin = parse_int_arg(optarg, "--pgbam-cleanup-margin"); break;
            case kPgbamCleanupMinWinningOption: opts.pgbam_cleanup_min_winning_threads = parse_int_arg(optarg, "--pgbam-cleanup-min-winning"); break;
            case kNoPgbamRelaxedCleanupPassOption: opts.pgbam_relaxed_cleanup_pass = false; break;
            case kPgbamRelaxedCleanupMarginOption: opts.pgbam_relaxed_cleanup_polarity_margin = parse_int_arg(optarg, "--pgbam-relaxed-cleanup-margin"); break;
            case kPgbamRelaxedCleanupMinWinningOption: opts.pgbam_relaxed_cleanup_min_winning_threads = parse_int_arg(optarg, "--pgbam-relaxed-cleanup-min-winning"); break;
            case kChunkSizeOption:      opts.chunk_size = parse_ll_arg(optarg, "--chunk-size"); break;
            case kNoisyRegMergeDisOption: opts.noisy_reg_merge_dis = parse_int_arg(optarg, "--noisy-merge-dis"); break;
            case kMinSvLenOption:       opts.min_sv_len = parse_int_arg(optarg, "--min-sv-len"); break;
            case kNoisySlideWinOption:  opts.noisy_reg_slide_win = parse_int_arg(optarg, "--noisy-slide-win"); break;
            case kDebugSiteOption:      opts.debug_site = optarg; break;
            case kPhaseMatrixDumpOption: opts.phase_matrix_dump_prefix = optarg; break;
            case 'X': extra_bam_files.push_back(optarg); break;
            case 'L': opts.input_is_list = true; break;
            case kHifiOption:           set_read_technology(ReadTechnology::Hifi); break;
            case kOntOption:            set_read_technology(ReadTechnology::Ont); break;
            case kShortReadsOption:     set_read_technology(ReadTechnology::ShortReads); break;
            case kStrandBiasPvalOption: opts.strand_bias_pval = parse_double_arg(optarg, "--strand-bias-pval"); break;
            case kNoisyMaxXgapsOption:  opts.noisy_reg_max_xgaps = parse_int_arg(optarg, "--noisy-max-xgaps"); break;
            case kRefOption:            opts.ref_fasta = optarg; break;
            case kBamOption:            opts.bam_files.push_back(optarg); break;
            case 'V': opts.verbose = parse_int_arg(optarg, "--verbose"); break;
            case 'h': print_collect_help(); return 0;
            default:  print_collect_help(); return 1;
        }
    }

    if (opts.threads < 1 || opts.min_mapq < 0 || opts.min_bq < 0 || opts.chunk_size < 1 ||
        opts.min_depth < 0 || opts.min_alt_depth < 0 || opts.noisy_reg_merge_dis < 0 ||
        opts.min_sv_len < 0 || opts.min_af < 0.0 || opts.max_af < opts.min_af ||
        opts.strand_bias_pval < 0.0 || opts.strand_bias_pval > 1.0 ||
        opts.max_var_ratio_per_read < 0.0 || opts.max_noisy_frac_per_read < 0.0 ||
        opts.noisy_reg_max_xgaps < 0 ||
        opts.noisy_reg_slide_win < -1 || opts.verbose < 0 ||
        opts.pgbam_primary_polarity_margin < 1 || opts.pgbam_primary_min_winning_threads < 1 ||
        opts.pgbam_cleanup_polarity_margin < 1 || opts.pgbam_cleanup_min_winning_threads < 1 ||
        opts.pgbam_relaxed_cleanup_polarity_margin < 1 || opts.pgbam_relaxed_cleanup_min_winning_threads < 1) {
        std::cerr << "Error: numeric thresholds are invalid\n";
        return 1;
    }
    if (read_technology_conflict) {
        std::cerr << "Error: choose only one of --hifi, --ont, or --short-reads\n";
        return 1;
    }
    if (opts.input_is_list && !opts.bam_files.empty()) {
        const std::string list_path = opts.bam_files.front();
        opts.bam_files = load_bam_list(list_path);
    }
    opts.bam_files.insert(opts.bam_files.end(), extra_bam_files.begin(), extra_bam_files.end());

    if (optind < argc) {
        std::cerr << "Error: unexpected positional argument: " << argv[optind]
                  << "\n       Use --ref and --bam instead of positional arguments.\n";
        return 1;
    }
    if (opts.ref_fasta.empty()) {
        std::cerr << "Error: --ref is required\n";
        print_collect_help();
        return 1;
    }
    if (opts.bam_files.empty()) {
        std::cerr << "Error: --bam is required\n";
        print_collect_help();
        return 1;
    }
    opts.bam_file = opts.bam_files.front();

    try {
        run_collect_bam_variation(opts);
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << "\n";
        return 1;
    }

    return 0;
}
