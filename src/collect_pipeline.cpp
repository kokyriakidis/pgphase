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

    // Backfill allele counts on graph-only candidates from BAM profiles.
    // The BAM profile builder records alleles but doesn't update candidate
    // counts; this pass accumulates the missing ref/alt/total coverage.
    if (graph_authoritative)
        clear_bam_evidence_at_graph_candidates(chunk, graph_owned_cands);
    else
        backfill_graph_candidate_counts(chunk, graph_only_cands);

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
            if (credit_bridge_snps && read.n_bridge_agree_snps > 0) {
                const int bridge_margin =
                    (read.n_clean_agree_snps + read.n_bridge_agree_snps) -
                    (read.n_clean_conflict_snps + read.n_bridge_conflict_snps);
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
static constexpr uint32_t kGapEvidenceCacheVersion = 2;
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
            if (!proposal.candidates.back().lcd_make_variants_region_pass &&
                candidate.lcd_make_variants_region_pass)
                proposal.candidates.back() = candidate;
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
                proposal.read_var_profile.push_back(std::move(profile));
            }
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
                if (allele_i >= source_profile.alleles.size() ||
                    source_profile.alleles[allele_i] == -1)
                    continue;
                dest.alleles[dest_vi] = source_profile.alleles[allele_i];
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
};

static GapRecoveryJobResult recover_one_hybrid_gap(
        std::vector<PhasingChunk>& chunks, const Options& opts,
        const std::string& contig,
        const std::unordered_map<std::string, std::string>& chrom_remap,
        const BamAuthorityIntervals* bam_authority,
        const GapReadIndex& read_index, const PhaseGap& initial_gap,
        size_t gap_index, std::mutex& chunks_mutex) {
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
    if (!opts.phase_matrix_dump_prefix.empty()) {
        local_opts.phase_matrix_dump_prefix =
            opts.phase_matrix_dump_prefix + ".tid" + std::to_string(gap.tid) +
            ".gap" + std::to_string(gap_index);
    }
    (void)chrom_remap;
    (void)bam_authority;
    std::vector<PhasingChunk> local;
    {
        std::lock_guard<std::mutex> lock(chunks_mutex);
        local.push_back(build_cached_gap_proposal(chunks, window));
    }
    auto& proposal = local.front();
    assign_hap_based_on_germline_het_vars_kmeans(
        proposal, local_opts, kCandGermlineClean);
    std::vector<uint32_t> original_flags;
    original_flags.reserve(proposal.candidates.size());
    for (const CandidateVariant& candidate : proposal.candidates)
        original_flags.push_back(candidate.lcd_var_i_to_cate);
    constexpr int kGapHomopolymerTier = 4;
    std::ostringstream report_rows;
    for (int tier = 1; tier <= kGapHomopolymerTier; ++tier) {
        const size_t previous_sites = proposal.candidates.size();
        if (tier == kGapHomopolymerTier) {
            for (size_t vi = 0; vi < proposal.candidates.size(); ++vi)
                proposal.candidates[vi].lcd_var_i_to_cate = original_flags[vi];
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
            local_opts.private_msa_admit_all_in_region = true;
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
        {
            std::lock_guard<std::mutex> lock(chunks_mutex);
            result = stitch_gap_proposal(
                chunks, proposal, gap, local_opts, &read_index, true);
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
                    << result.reads_added << '\t'
                    << (tier == kGapHomopolymerTier && !result.joined
                            ? "rejected"
                            : result.joined
                                  ? "joined"
                                  : result.left_linked || result.right_linked
                                        ? "partial"
                                        : "open")
                    << '\n';
        if (result.joined) {
            job_result.joined = true;
            job_result.has_edge = true;
            job_result.edge = {gap.left_ps, gap.right_ps, result.right_flip};
            break;
        }
    }
    job_result.report_rows = report_rows.str();
    return job_result;
}

static void recover_hybrid_gaps(std::vector<PhasingChunk>& chunks, const Options& opts,
                                const std::string& contig,
                                const std::unordered_map<std::string, std::string>& chrom_remap,
                                const BamAuthorityIntervals* bam_authority,
                                std::ostream* report) {
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
    const GapReadIndex read_index(chunks);
    std::vector<size_t> pending(initial_gaps.size());
    std::iota(pending.begin(), pending.end(), 0);
    std::vector<GapRecoveryJobResult> results(initial_gaps.size());
    std::mutex chunks_mutex;
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
                            chunks_mutex);
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
    std::cerr << "Gap recovery: " << initial_gaps.size() << " initial gaps, "
              << joined << " joined in " << wave_count << " wave(s) on "
              << contig << ", " << edge_conflicts
              << " conflicting edge(s) rejected\n";
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
        recovery_report << "CHROM\tGAP_LEFT\tGAP_RIGHT\tTIER\tWINDOW_BEG\tWINDOW_END\tNEW_SITES\tMSA_HET_SNPS\tMSA_HET_INDELS\tLEFT_LINK\tRIGHT_LINK\tREADS_ADDED\tSTATUS\n";
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
        filter_hybrid_reads_by_margin(batch.chunks, opts.min_read_hap_margin,
                                      opts.private_msa_admit_all_in_region);
        filter_hybrid_small_phase_sets(batch.chunks, opts.min_phase_set_reads);
        if (opts.recover_gaps) {
            recover_hybrid_gaps(batch.chunks, opts, graph_query_contig, chrom_remap,
                                bam_authority_ptr, recovery_report.is_open() ? &recovery_report : nullptr);
            stitch_chunk_haps(batch.chunks, &opts, pgbam_sidecar.get());
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
