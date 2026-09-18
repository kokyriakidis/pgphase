/**
 * @file collect_pipeline.cpp
 * @brief Region chunking, parallel candidate collection, streaming writers, and CLI for collect-bam-variation.
 *
 * @details Coordinates are 1-based inclusive (BED is converted in `load_bed_regions`). Chunks are batched by
 * `reg_chunk_i` in `run_collect_bam_variation` so TSV/VCF can be streamed without holding
 * all candidates in memory.
 */

#include "collect_pipeline.hpp"
#include <array>

#include "arg_parse.hpp"
#include "bam_digar.hpp"
#include "collect_bam_output.hpp"
#include "collect_output.hpp"
#include "collect_phase.hpp"
#include "collect_phase_pgbam.hpp"
#include "collect_var.hpp"

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
///
/// Step 1 of --retry-unphased-with-bam, and the only step that decides WHERE the
/// retry applies. Its output becomes `retry_opts.retry_windows`, which
/// allele_depths_call_het then uses to confine the widened het admission; a
/// window this function does not report is a window the retry cannot touch.
///
/// Parameters, both exposed and both acting as floors on what counts as a
/// failure worth re-solving:
///   min_reads -- unphased reads that must pile up in the interval
///                (--retry-min-unphased-reads, default 5). Below this the
///                interval is a few stray reads, not a failed window.
///   min_bp    -- how long the interval must be
///                (--retry-min-window-bp, default 10000). Below this the solve
///                did not fail over a span worth re-solving.
/// Positions are accumulated in 1 kb bins, so both floors are applied to
/// contiguous runs of bins rather than to individual reads.
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

// Recover a window the first solve could not phase by running the ordinary
// alignment pipeline over that region alone, then stitching its answer in.
//
// The gap this exists for is a mapping-quality hole: at
// chr20:26,029,591-26,088,679 every read carries MAPQ 3, the default floor of
// 30 keeps them out of the chunk, and the pipeline discovers ZERO candidates
// across 50 kb while a competitor phases 120 heterozygotes there. Lowering the
// floor for the whole run closes the gap but also re-solves the flanks, where
// the pipeline was already right: the left flank falls from 100% to 93.3% read
// concordance. Waking those reads inside the parent chunk instead discovers the
// sites without hurting the flanks, but the interior fragments into 1-, 2- and
// 6-site blocks, because it is being solved inside a chunk whose read set is
// mostly asleep.
//
// So the region gets its own chunk, at its own floor, solved by the same code
// as any other chunk. Measured standalone over this window: one block of 393
// sites spanning 118.7 kb at 99.54% read concordance, both flanks consistent,
// sharing 95 reads with the parent's left block and 121 with its right.
//
// The stitch is the pipeline's own rule. Reads tagged in both solves vote an
// n11/n12/n21/n22 table per parent block and select_stitch_orientation decides
// -- the same function, and the same default net-margin standard, that joins
// adjacent chunks. Its refusal carries over too: a parent block sharing no read
// with the targeted solve cannot be merged, which is the guard a bespoke gap
// link does not have.
// One read length: beyond it no read reaches the window, so a wider region
// adds sites the stitch cannot use. Matches the seam span the block-link
// evidence saturates at.
constexpr hts_pos_t kTargetedSolveFlank = 30000;

/// Intervals between consecutive phase blocks -- the other kind of failure.
///
/// collect_unphased_windows finds where the solve left READS unphased. That is
/// one failure mode, and it is the one the alignment-driven pipeline produces:
/// a stretch with no usable site leaves its reads untagged. But a solve can
/// also place every read and still not join, leaving two blocks whose relative
/// phase is unknown. The reads there are phased, so no unphased-read window is
/// reported, and the targeted solve never sees the seam.
///
/// That distinction decides whether a graph-first hybrid can work. Measured on
/// chr20:5,309,406 with graph sites only outside the gap: three blocks with
/// 34.6 kb and 12.0 kb seams between them, the retry reporting 3 unphased-read
/// windows, none at a seam, and no targeted solve running. The blocks stayed
/// apart because nothing asked about the interval between them.
///
/// Seams are taken from the candidates rather than the reads because a block's
/// extent is first to last phased SITE; read starts understate it by up to a
/// read length at each end, a mistake made earlier in this work.
static std::vector<std::pair<hts_pos_t, hts_pos_t>> collect_block_seams(
        const PhasingChunk& chunk) {
    std::map<hts_pos_t, std::pair<hts_pos_t, hts_pos_t>> extent;
    for (const CandidateVariant& cand : chunk.candidates) {
        if (cand.phase_set == 0) continue;
        if (cand.hap_alt == 0 && cand.hap_ref == 0) continue;
        const hts_pos_t pos = cand.key.sort_pos();
        auto it = extent.find(cand.phase_set);
        if (it == extent.end()) extent.emplace(cand.phase_set, std::make_pair(pos, pos));
        else {
            it->second.first = std::min(it->second.first, pos);
            it->second.second = std::max(it->second.second, pos);
        }
    }
    std::vector<std::pair<hts_pos_t, hts_pos_t>> spans;
    spans.reserve(extent.size());
    for (const auto& kv : extent) spans.push_back(kv.second);
    std::sort(spans.begin(), spans.end());

    std::vector<std::pair<hts_pos_t, hts_pos_t>> seams;
    for (size_t i = 1; i < spans.size(); ++i) {
        const hts_pos_t beg = spans[i - 1].second;
        const hts_pos_t end = spans[i].first;
        if (end <= beg) continue;   // overlapping or touching blocks
        seams.emplace_back(beg, end);
    }
    return seams;
}

/// The region handed to process_chunk uses the chunk's own contig id, which is
/// the BAM's because the hybrid is the only caller. A caller whose header is not
/// the BAM's -- the graph pipeline built a synthetic one in reference-index
/// order -- must translate by contig NAME rather than pass its tid through;
/// passing it through targeted a different contig outright.
///
/// Every site the sub-solve phases is imported. A caller whose output table is
/// index-parallel to per-site metadata cannot accept that: appending candidates,
/// and the reorder that follows, decouples the arrays and mispairs records.
/// Sub-solve options: the recovery MAPQ floor, the alignment pipeline's own
/// noisy-k-means default, no recursion, and one thread because the caller
/// parallelises across regions instead.
static Options targeted_solve_options(const Options& opts) {
    Options sub = opts;
    sub.min_mapq = std::min(opts.min_mapq, opts.recovery_min_mapq);
    sub.skip_noisy_kmeans = false;
    sub.retry_unphased_with_bam = false;
    sub.threads = 1;
    sub.verbose = 0;
    return sub;
}

struct TargetedWindowGroup {
    RegionChunk region;
    std::vector<std::pair<hts_pos_t, hts_pos_t>> members;
};

/// Windows padded by kTargetedSolveFlank and merged where the padded regions
/// touch, so one region is one piece of work. Shared by the per-chunk path and
/// the batch prewarm so both derive identical region keys.
static std::vector<TargetedWindowGroup> build_targeted_groups(
        const std::vector<std::pair<hts_pos_t, hts_pos_t>>& windows, int solve_tid) {
    std::vector<std::pair<hts_pos_t, hts_pos_t>> window_list(windows.begin(), windows.end());
    std::sort(window_list.begin(), window_list.end());
    std::vector<TargetedWindowGroup> groups;
    groups.reserve(window_list.size());
    for (const auto& window : window_list) {
        const hts_pos_t beg = std::max<hts_pos_t>(1, window.first - kTargetedSolveFlank);
        const hts_pos_t end = window.second + kTargetedSolveFlank;
        if (!groups.empty() && beg <= groups.back().region.end) {
            groups.back().region.end = std::max(groups.back().region.end, end);
            groups.back().members.push_back(window);
            continue;
        }
        TargetedWindowGroup group;
        group.region.tid = solve_tid;
        group.region.beg = beg;
        group.region.end = end;
        group.region.chunk_id = -1;
        group.members.push_back(window);
        groups.push_back(std::move(group));
    }
    return groups;
}

/// Keep only what the apply phase reads, so a chromosome's worth of solved
/// regions can be held at once.
static TargetedSolveResult slim_targeted_result(PhasingChunk&& solved) {
    TargetedSolveResult out;
    if (solved.haps.size() != solved.reads.size()) return out;
    out.qnames.reserve(solved.reads.size());
    for (const auto& read : solved.reads) out.qnames.push_back(read.qname);
    out.haps = std::move(solved.haps);
    out.phase_sets = std::move(solved.phase_sets);
    out.candidates = std::move(solved.candidates);
    return out;
}

/// Solve the given regions in parallel, largest first.
///
/// Largest-first matters because the regions are very uneven -- a merged group
/// can hold twenty windows -- and taking them in coordinate order routinely
/// leaves the biggest one starting last with every other thread idle. Each
/// worker builds its own WorkerContext: htslib file handles are not shareable,
/// and this is what collect_chunk_batch_parallel does for the main chunk loop.
static std::vector<TargetedSolveResult> solve_targeted_regions(
        const std::vector<RegionChunk>& regions, const Options& sub, int threads) {
    std::vector<TargetedSolveResult> out(regions.size());
    if (regions.empty()) return out;

    std::vector<size_t> order(regions.size());
    for (size_t i = 0; i < order.size(); ++i) order[i] = i;
    std::sort(order.begin(), order.end(), [&](size_t a, size_t b) {
        return (regions[a].end - regions[a].beg) > (regions[b].end - regions[b].beg);
    });

    const size_t worker_count = std::min<size_t>(
        std::max<size_t>(1, static_cast<size_t>(threads)), regions.size());
    std::atomic<size_t> next_region{0};
    std::exception_ptr first_error;
    std::mutex error_mutex;
    std::vector<std::thread> workers;
    workers.reserve(worker_count);
    for (size_t worker_i = 0; worker_i < worker_count; ++worker_i) {
        workers.emplace_back([&]() {
            try {
                WorkerContext local_context(sub);
                while (true) {
                    const size_t k = next_region.fetch_add(1);
                    if (k >= order.size()) break;
                    const size_t ri = order[k];
                    out[ri] = slim_targeted_result(
                        process_chunk(regions[ri], sub, local_context));
                }
            } catch (...) {
                std::lock_guard<std::mutex> lock(error_mutex);
                if (!first_error) first_error = std::current_exception();
            }
        });
    }
    for (std::thread& worker : workers) worker.join();
    if (first_error) std::rethrow_exception(first_error);
    return out;
}

static size_t recover_windows_with_targeted_solve(
        PhasingChunk& chunk,
        const std::vector<std::pair<hts_pos_t, hts_pos_t>>& windows,
        const Options& opts,
        WorkerContext& context,
        int solve_tid,
        bool allow_import,
        TargetedSolveCache* cache) {
    if (windows.empty()) return 0;

    // Parent block extents, first to last phased site, so the import interval
    // can be the space between the blocks actually being joined.
    std::map<hts_pos_t, std::pair<hts_pos_t, hts_pos_t>> parent_extent;
    for (const CandidateVariant& cand : chunk.candidates) {
        if (cand.phase_set == 0) continue;
        if (cand.hap_alt == 0 && cand.hap_ref == 0) continue;
        const hts_pos_t pos = cand.key.sort_pos();
        auto it = parent_extent.find(cand.phase_set);
        if (it == parent_extent.end())
            parent_extent.emplace(cand.phase_set, std::make_pair(pos, pos));
        else {
            it->second.first = std::min(it->second.first, pos);
            it->second.second = std::max(it->second.second, pos);
        }
    }

    // Index the parent's tagged reads once.
    std::unordered_map<std::string, size_t> parent_by_name;
    parent_by_name.reserve(chunk.reads.size() * 2u);
    for (size_t i = 0; i < chunk.reads.size(); ++i)
        parent_by_name.emplace(chunk.reads[i].qname, i);

    size_t merged_total = 0;
    size_t imported_total = 0;
    size_t adopted_total = 0;
    // The sub-solve options do not depend on the window.
    const Options sub = targeted_solve_options(opts);
    std::vector<TargetedWindowGroup> groups = build_targeted_groups(windows, solve_tid);
    std::vector<RegionChunk> regions;
    regions.reserve(groups.size());
    for (const TargetedWindowGroup& group : groups) regions.push_back(group.region);
    if (opts.verbose > 0)
        fprintf(stderr, "[targeted] %zu window(s) -> %zu merged region(s)\n",
                windows.size(), groups.size());

    // Regions already solved by a batch prewarm are taken from the cache; the
    // rest are solved here, still in parallel. Either way the apply loop below
    // runs in region order, so the outcome does not depend on where a solve came
    // from or on which one finished first.
    std::vector<TargetedSolveResult> solved_results(regions.size());
    std::vector<RegionChunk> to_solve;
    std::vector<size_t> to_solve_idx;
    for (size_t ri = 0; ri < regions.size(); ++ri) {
        if (cache != nullptr) {
            auto it = cache->find(std::make_tuple(regions[ri].tid, regions[ri].beg,
                                                  regions[ri].end));
            if (it != cache->end()) {
                solved_results[ri] = std::move(it->second);
                continue;
            }
        }
        to_solve.push_back(regions[ri]);
        to_solve_idx.push_back(ri);
    }
    if (!to_solve.empty()) {
        std::vector<TargetedSolveResult> fresh =
            solve_targeted_regions(to_solve, sub, opts.threads);
        for (size_t k = 0; k < fresh.size(); ++k)
            solved_results[to_solve_idx[k]] = std::move(fresh[k]);
    }

    for (size_t window_i = 0; window_i < regions.size(); ++window_i) {
        const RegionChunk& region = regions[window_i];
        const std::vector<std::pair<hts_pos_t, hts_pos_t>>& group_windows =
            groups[window_i].members;
        const TargetedSolveResult& solved = solved_results[window_i];
        if (solved.haps.size() != solved.qnames.size()) continue;

        // votes[parent_ps][sub_ps] = n11, n12, n21, n22 over reads tagged in both
        std::map<std::pair<hts_pos_t, hts_pos_t>, std::array<int, 4>> votes;
        for (size_t j = 0; j < solved.qnames.size(); ++j) {
            const int sub_hap = solved.haps[j];
            if (sub_hap != 1 && sub_hap != 2) continue;
            auto it = parent_by_name.find(solved.qnames[j]);
            if (it == parent_by_name.end()) continue;
            const size_t i = it->second;
            if (i >= chunk.haps.size()) continue;
            const int par_hap = chunk.haps[i];
            if (par_hap != 1 && par_hap != 2) continue;
            const auto key = std::make_pair(chunk.phase_sets[i], solved.phase_sets[j]);
            votes[key][(par_hap == 1 ? 0 : 2) + (sub_hap == 1 ? 0 : 1)] += 1;
        }

        // A parent block joins the targeted solve when the standard says so.
        // Record, per targeted block, which parent blocks it carries and in
        // which orientation; two or more is a bridge.
        std::map<hts_pos_t, std::vector<std::pair<hts_pos_t, bool>>> bridged;
        for (const auto& kv : votes) {
            bool do_flip = false;
            if (!select_stitch_orientation(kv.second, &opts, do_flip)) continue;
            bridged[kv.first.second].emplace_back(kv.first.first, do_flip);
        }

        if (opts.verbose > 0) {
            size_t linked = 0;
            for (const auto& kv : bridged) linked += kv.second.size();
            fprintf(stderr,
                    "[targeted] %lld-%lld: %zu vote pair(s), %zu parent block(s) linked,"
                    " %zu targeted block(s) carrying links\n",
                    (long long)region.beg, (long long)region.end,
                    votes.size(), linked, bridged.size());
        }
        for (const auto& kv : bridged) {
            if (kv.second.size() < 2) continue;  // nothing bridged
            const hts_pos_t keep_ps = kv.second.front().first;
            const bool keep_flip = kv.second.front().second;
            for (size_t k = 1; k < kv.second.size(); ++k) {
                const hts_pos_t drop_ps = kv.second[k].first;
                // Relative orientation of the two parent blocks, via the
                // targeted solve they both agree with.
                const bool flip = (kv.second[k].second != keep_flip);
                for (size_t i = 0; i < chunk.phase_sets.size(); ++i) {
                    if (chunk.phase_sets[i] != drop_ps) continue;
                    chunk.phase_sets[i] = keep_ps;
                    if (flip && (chunk.haps[i] == 1 || chunk.haps[i] == 2))
                        chunk.haps[i] = 3 - chunk.haps[i];
                }
                for (CandidateVariant& cand : chunk.candidates) {
                    if (cand.phase_set != drop_ps) continue;
                    cand.phase_set = keep_ps;
                    if (flip) std::swap(cand.hap_alt, cand.hap_ref);
                }
                ++merged_total;
            }

            // Carry across the sites in the space being joined. The parent never
            // discovered them -- that is what a mapping-quality hole means -- so
            // without this the block spans an interval it reports nothing in and
            // the heterozygotes the competitor calls there stay invisible.
            //
            // The bound is the space between the blocks the vote just joined,
            // not the detected window. The detected window is where the solve
            // left READS unphased, or a seam, and it is routinely narrower than
            // the unreported span: on chr20:26,029,591 the narrow bound imported
            // 38 sites where the same solve had 282 to give, which cost 241
            // in-gap records. Outside the joined blocks the parent has its own,
            // better-supported calls, so the bound stays closed there.
            std::vector<std::pair<hts_pos_t, hts_pos_t>> import_spans;
            {
                std::vector<std::pair<hts_pos_t, hts_pos_t>> joined;
                for (const auto& entry : kv.second) {
                    auto pe = parent_extent.find(entry.first);
                    if (pe != parent_extent.end()) joined.push_back(pe->second);
                }
                std::sort(joined.begin(), joined.end());
                for (size_t j = 1; j < joined.size(); ++j)
                    if (joined[j].first > joined[j - 1].second)
                        import_spans.emplace_back(joined[j - 1].second, joined[j].first);
                if (import_spans.empty())
                    for (const auto& member : group_windows)
                        import_spans.emplace_back(member.first, member.second);
            }
            for (const CandidateVariant& src : solved.candidates) {
                if (src.phase_set != kv.first) continue;
                const hts_pos_t pos = src.key.sort_pos();
                bool inside = false;
                for (const auto& span : import_spans)
                    if (pos > span.first && pos < span.second) { inside = true; break; }
                if (!inside) continue;
                if (src.hap_alt == 0 && src.hap_ref == 0) continue;

                // A site the parent already holds but could not phase is the
                // common case here, not the exception: the graph catalog is
                // injected over the whole region, so inside a mapping-quality
                // hole the parent carries the sites with almost no read support
                // and leaves them unphased. Measured on chr20:26,029,591 -- of
                // 405 candidates the targeted solve phased in the joined
                // interval, 244 were already present. Skipping them as
                // duplicates left them unphased and unemitted, which is why the
                // interval reported 35 records instead of 276.
                //
                // So adopt rather than skip: take the targeted solve's phasing
                // and its counts, because the parent's counts are the starved
                // ones the hole produced and the solve's come from reads that
                // can actually see the site.
                CandidateVariant* existing = nullptr;
                for (CandidateVariant& have : chunk.candidates)
                    if (have.key.pos == src.key.pos && have.key.type == src.key.type &&
                        have.key.ref_len == src.key.ref_len && have.key.alt == src.key.alt) {
                        existing = &have; break;
                    }
                if (existing != nullptr) {
                    if (existing->phase_set != 0) continue;  // the parent's own call stands
                    // Take the targeted solve's COUNTS only where the parent's
                    // are unusable. That overwrite exists for the
                    // mapping-quality hole, where the parent holds injected
                    // catalog sites with almost no read support and its depths
                    // are starved -- there the solve's counts are the real ones.
                    //
                    // Applied unconditionally it destroys good calls. Measured
                    // when the graph pipeline drives the first pass, whose
                    // counts come from GAF evidence and are not starved:
                    // chr20:5,393,876 went from CLEAN_HET_SNP at DP 67, 27/40 to
                    // LOW_AF at DP 45, 1/44, and nine such candidates were then
                    // pruned outright -- 33 candidates down to 24, 23 phased
                    // heterozygotes down to 14. The phasing is always adopted;
                    // the evidence is not.
                    const bool parent_starved =
                        existing->counts.category == VariantCategory::LowCoverage ||
                        existing->counts.total_cov < opts.min_depth;
                    if (parent_starved) {
                        existing->counts = src.counts;
                        existing->lcd_var_i_to_cate = src.lcd_var_i_to_cate;
                    }
                    existing->hap_alt = keep_flip ? src.hap_ref : src.hap_alt;
                    existing->hap_ref = keep_flip ? src.hap_alt : src.hap_ref;
                    // The emitter derives the genotype from the CONSENSUS
                    // alleles, not from hap_alt/hap_ref: is_alt_genotype calls
                    // derive_hap_alt_ref_from_consensus, which reads
                    // hap_to_cons_alle[1] and [2]. Setting hap_alt/hap_ref alone
                    // left 244 adopted sites carrying a phase set, a clean
                    // category and DP near 80 that the VCF still dropped.
                    existing->hap_to_cons_alle = src.hap_to_cons_alle;
                    if (keep_flip)
                        std::swap(existing->hap_to_cons_alle[1],
                                  existing->hap_to_cons_alle[2]);
                    existing->phase_set = keep_ps;
                    existing->lcd_make_variants_region_pass = true;
                    ++adopted_total;
                    continue;
                }
                if (!allow_import) continue;
                CandidateVariant imported = src;
                imported.phase_set = keep_ps;
                if (keep_flip) std::swap(imported.hap_alt, imported.hap_ref);
                imported.lcd_make_variants_region_pass = true;
                // Relabel to the chunk's contig id. The solve ran in the BAM's
                // header, which is not the caller's when the graph pipeline is
                // driving -- its synthetic header is in reference-index order.
                imported.key.tid = chunk.region.tid;
                chunk.candidates.push_back(std::move(imported));
                ++imported_total;
            }
            if (opts.verbose > 0)
                fprintf(stderr,
                        "[targeted] %s:%lld-%lld bridged %zu block(s) into PS %lld\n",
                        context.primary_header() != nullptr &&
                                chunk.region.tid < context.primary_header()->n_targets
                            ? context.primary_header()->target_name[chunk.region.tid] : ".",
                        (long long)region.beg, (long long)region.end,
                        kv.second.size() - 1, (long long)keep_ps);
        }
    }
    if (adopted_total > 0 && opts.verbose > 0)
        fprintf(stderr, "[targeted] adopted phasing for %zu site(s) the parent held unphased\n",
                adopted_total);
    if (imported_total > 0) {
        // No reorder. merge_chunk_candidates sorts by variant key before the
        // records are written, so sorting here bought nothing, and keeping the
        // appended candidates at the tail is what lets a caller whose output is
        // index-parallel to per-site metadata extend those arrays.
        if (opts.verbose > 0)
            fprintf(stderr, "[targeted] imported %zu site(s) the parent never discovered\n",
                    imported_total);
    }
    return merged_total;
}

/// Recover what a solve could not phase, from alignment evidence.
///
/// Two window sources, because there are two failure modes and they do not
/// overlap: collect_unphased_windows reports where reads were left unphased --
/// what an alignment-driven solve produces -- and collect_block_seams reports
/// intervals between consecutive blocks, which is what a catalog-driven solve
/// produces, every read placed but the blocks unjoined.
///
/// `contig_name` names the chunk's contig so the BAM's own tid can be resolved.
/// Pass nullptr when the caller's header IS the BAM's; the graph pipeline builds
/// a synthetic header in reference-index order, so its tid for a contig is not
/// the BAM's and passing it through targeted a different contig outright.
///
/// `allow_import` adds sites the caller's pass never discovered. A caller whose
/// output table is index-parallel to per-site metadata must either extend those
/// arrays for the appended tail or pass false: appending candidates it has no
/// metadata for makes graph_chunks_to_candidate_table skip them.
size_t recover_unphased_windows_from_bam(PhasingChunk& chunk, const Options& opts,
                                        WorkerContext& context,
                                        const char* contig_name,
                                        bool allow_import,
                                        TargetedSolveCache* cache) {
    const int solve_tid =
        contig_name != nullptr
            ? sam_hdr_name2tid(context.primary_header(), contig_name)
            : chunk.region.tid;
    if (solve_tid < 0) return 0;
    std::vector<std::pair<hts_pos_t, hts_pos_t>> windows = collect_unphased_windows(
        chunk, opts.retry_min_unphased_reads, opts.retry_min_window_bp);
    for (const auto& seam : collect_block_seams(chunk))
        windows.push_back(seam);
    if (windows.empty()) return 0;
    std::sort(windows.begin(), windows.end());
    return recover_windows_with_targeted_solve(chunk, windows, opts, context, solve_tid,
                                              allow_import, cache);
}

void prewarm_targeted_solves(
        const std::vector<std::pair<PhasingChunk*, const char*>>& chunks,
        const Options& opts, WorkerContext& context, TargetedSolveCache& cache) {
    std::vector<RegionChunk> regions;
    for (const auto& entry : chunks) {
        PhasingChunk* chunk = entry.first;
        if (chunk == nullptr) continue;
        const int solve_tid = entry.second != nullptr
                                  ? sam_hdr_name2tid(context.primary_header(), entry.second)
                                  : chunk->region.tid;
        if (solve_tid < 0) continue;
        std::vector<std::pair<hts_pos_t, hts_pos_t>> windows = collect_unphased_windows(
            *chunk, opts.retry_min_unphased_reads, opts.retry_min_window_bp);
        for (const auto& seam : collect_block_seams(*chunk)) windows.push_back(seam);
        if (windows.empty()) continue;
        for (const TargetedWindowGroup& group : build_targeted_groups(windows, solve_tid))
            regions.push_back(group.region);
    }
    if (regions.empty()) return;

    // Regions from different chunks can still be duplicates of each other at a
    // chunk boundary, where both sides see the same seam. Solve each key once.
    std::sort(regions.begin(), regions.end(), [](const RegionChunk& a, const RegionChunk& b) {
        return std::make_tuple(a.tid, a.beg, a.end) < std::make_tuple(b.tid, b.beg, b.end);
    });
    regions.erase(std::unique(regions.begin(), regions.end(),
                              [](const RegionChunk& a, const RegionChunk& b) {
                                  return a.tid == b.tid && a.beg == b.beg && a.end == b.end;
                              }),
                  regions.end());

    if (opts.verbose > 0)
        fprintf(stderr, "[targeted] prewarming %zu region(s) across %d thread(s)\n",
                regions.size(), opts.threads);

    const Options sub = targeted_solve_options(opts);
    std::vector<TargetedSolveResult> solved = solve_targeted_regions(regions, sub, opts.threads);
    for (size_t ri = 0; ri < regions.size(); ++ri)
        cache.emplace(std::make_tuple(regions[ri].tid, regions[ri].beg, regions[ri].end),
                      std::move(solved[ri]));
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
    // Graph-first: the catalog's sites drive the first pass, so the alignment
    // channel's own discoveries do not enter it. Clearing here rather than
    // after injection is deliberate -- and necessary. Gating on ownership after
    // the claim pass does not work: measured over
    // chr20:25,979,591-26,138,679, injection claims 3,741 of 3,743 candidates,
    // whether tested by the graph_site flag or by its all_graph_cands index
    // set, so an ownership filter withholds 2 sites and the mode is inert. The
    // catalog's sites are added by inject_graph_sites immediately below, and
    // recover_windows_with_targeted_solve puts alignment sites back inside
    // every window the first pass could not phase.
    // The catalog's sites are the phasing anchors: the alignment channel's own
    // candidates do not enter the first solve. They are not wasted -- the
    // discovery that produced them also built chunk.noisy_regions, which scopes
    // the noisy-region MSA that supplies most of this chunk's phased sites.
    // Measured by skipping discovery outright on
    // chr20:25,979,591-26,138,679: 31% faster (4.05 s to 2.78 s) for a collapse
    // from 447 phased heterozygotes to 137. So discovery runs, its candidates
    // are withheld here, and the alignment returns in the gaps through
    // recover_windows_with_targeted_solve.
    {
        const size_t withheld = chunk.candidates.size();
        chunk.candidates.clear();
        if (opts.verbose > 0)
            fprintf(stderr,
                    "[hybrid] withheld %zu alignment-discovered candidate(s) from the"
                    " first solve; the catalog's sites phase it\n", withheld);
    }

    // Load graph sites and GAF reads for this region.
    // What the catalog load kept and what it discarded. A loader that drops a
    // record silently looks exactly like one that loaded the file correctly, so
    // the counts are reported rather than inferred from the downstream effect.
    GraphSiteCatalog chunk_catalog = load_sites_for_region(
        sites_handle, graph_query_contig, region.beg, region.end);
    if (opts.verbose >= 1) {
        std::cerr << "[graph-sites] " << graph_query_contig << ':' << region.beg << '-'
                  << region.end << ": " << chunk_catalog.stats.summary() << '\n';
    }
    for (GraphSite& s : chunk_catalog.sites) {
        auto it = chrom_remap.find(s.chrom);
        if (it != chrom_remap.end()) s.chrom = it->second;
        if (!s.ref_contig.empty()) {
            auto it2 = chrom_remap.find(s.ref_contig);
            if (it2 != chrom_remap.end()) s.ref_contig = it2->second;
        }
    }
    // Every catalog site in the chunk is injected. There used to be a BED that
    // excluded sites inside it so the alignment channel would "own" those
    // intervals; with the alignment's own candidates withheld from the first
    // solve, excluding the catalog's sites there left the interval with nothing
    // at all, so the option only deleted evidence.
    GraphSiteCatalogView chunk_view;
    chunk_view.source = &chunk_catalog.sites;
    chunk_view.indices.reserve(chunk_catalog.sites.size());
    for (size_t site_i = 0; site_i < chunk_catalog.sites.size(); ++site_i)
        chunk_view.indices.push_back(site_i);

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

    // The graph owns the candidates it contributed and the alignment channel
    // does not have. A catalog-MATCHED candidate is not owned: the alignment
    // channel measured it from reads, and a claim only promotes its category.
    //
    // An authoritative mode used to sit here, making every matched candidate
    // graph-owned and discarding the alignment channel's read evidence at each
    // one. Removed after measurement on chr20:5,309,406-5,345,085: it reduced
    // the hybrid to the graph channel's own phasing power -- the same 23 phased
    // heterozygotes collect-graph-variation reaches alone, against 36 by
    // default -- for 175 fewer reads tagged, 22 discordant against 7, and one
    // block fragmented into three. It discarded evidence and gained nothing the
    // graph could phase by itself. A whitelist (--private-sites) used to enable
    // it implicitly, which is what made that route look like a site-selection
    // change when it was an evidence change.
    const std::unordered_set<int>& graph_owned_cands = graph_only_cands;

    // Step 3.1: build BAM read profiles against augmented candidate table.
    collect_var_build_profiles(chunk, opts);

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

    // Steps 3.2-4: k-means plus the noisy-region MSA.
    collect_var_run_phasing(chunk, opts);

    // --retry-unphased-with-bam: the whole mechanism, in order.
    //
    // WHAT THE FIRST SOLVE LEFT BEHIND. It ran with the hybrid's own override
    // `skip_noisy_kmeans = true` (hybrid_collect.cpp), so only the CLEAN class
    // entered the k-means: clean het SNPs, clean het indels, clean hom. The
    // noisy-region MSA still ran and still built its candidates -- so the
    // window is NOT missing candidates. Measured inside
    // chr20:48,176,831-48,229,447 on stock defaults: 80 candidates, of which 6
    // are NOISY_CAND_HET and ZERO of those 6 carry a phase set. That is the
    // whole failure. The class that holds the interior evidence is present and
    // unoriented.
    //
    // And it is the only interior evidence there is. Across the six panel
    // windows the alignment channel places 266 CLEAN_HOM, 49 NOISY_CAND_HOM,
    // 24 NOISY_CAND_HET, zero clean het SNPs and exactly one clean het indel
    // strictly inside a gap; 17 of the 18 heterozygotes a competitor phases
    // inside these gaps are in our noisy class.
    //
    // WHAT THIS RETRY CHANGES, in order:
    //
    //   1. collect_unphased_windows -- find the intervals the solve failed on,
    //      measured from the result (no eligible het carries a phase set, or it
    //      carries one that does not continue the preceding block) rather than
    //      predicted from site spacing. These become `retry_windows`.
    //   2. collect_var_run_phasing  -- re-solve on a COPY of the options with
    //      three fields changed (see below). The chunk is mutated in place, so
    //      the second solve replaces the first one's labels.
    //   3. Inside that re-solve, collect_noisy_vars_step4 runs its k-means over
    //      the WIDER kCandGermlineVarCate mask, because skip_noisy_kmeans is now
    //      clear -- this is what orients the candidates the first solve left at
    //      phase_set 0. Same window after the retry: 82 candidates, 8
    //      NOISY_CAND_HET, and all 8 carry a phase set.
    //   4. force_noisy_msa additionally enables split_nested_msa_deletions in
    //      make_vars_from_msa_cons_aln, which is where those 2 extra candidates
    //      come from (6 -> 8).
    //
    // Measured on the six-window panel: the default spans 0 of 6 at 99.68% read
    // concordance; this arm spans 4 of 6 at 99.49% with ZERO reads moved from
    // concordant to discordant. Chromosome-wide it costs accuracy against the
    // default's 0.559% read Hamming, which is why it is a flag and not a
    // default -- it helps in the deficit windows and harms elsewhere.
    if (opts.retry_unphased_with_bam) {
        const auto windows = collect_unphased_windows(
            chunk, opts.retry_min_unphased_reads, opts.retry_min_window_bp);
        // Re-solve whenever a window failed. This used to be gated on a count of
        // categories restored from a recovery pass, which made the flag a no-op
        // on its own: with nothing zeroed there was nothing to restore, so the
        // re-solve below -- the part that actually admits the noisy class --
        // never ran.
        if (!windows.empty()) {
            // A COPY, deliberately: the three changes below must apply to this
            // re-solve and to nothing else in the chunk's remaining work.
            Options retry_opts = opts;
            // Asking for the MSA step BY NAME is deliberate. An earlier version
            // reached it by clearing a broader recovery flag, which also turned
            // off every other guard gated on that flag -- it silently disabled
            // split_nested_msa_deletions, which then emitted both nested forms
            // of one tandem-repeat deletion as independent heterozygotes.
            // What each field does, and why it is this field and not another:
            //
            //   force_noisy_msa    -- asks make_vars_from_msa_cons_aln for the
            //     noisy-region MSA by name, and also enables
            //     split_nested_msa_deletions inside it.
            //   skip_noisy_kmeans  -- clears the hybrid's own override, so
            //     collect_noisy_vars_step4 runs the kCandGermlineVarCate
            //     k-means (clean | noisy het | noisy hom) after recalling
            //     candidates. Without this the recalled sites exist but are
            //     never oriented, and nothing gains a phase set.
            //   retry_windows      -- the ONLY scoping. It confines
            //     allele_depths_call_het to these intervals, so the widened
            //     het admission applies where the solve failed and nowhere
            //     else. Admitting that class chunk-wide instead is a measured
            //     bad trade: chromosome-wide it roughly doubled the read
            //     Hamming error (0.878% -> 1.837%, 1,831 -> 3,415 discordant
            //     reads) while spanning 14 of 196 gaps.
            // Wake the reads the ordinary floor excluded, but only inside the
            // windows nothing could phase. A read below min_mapq was parsed and
            // marked skipped at load; here it becomes visible to discovery,
            // allele counting and the k-means for this re-solve. Reads skipped
            // for their variant load stay skipped -- skipped_for_mapq is what
            // distinguishes them.
            size_t woken = 0;
            for (ReadRecord& read : chunk.reads) {
                if (!read.skipped_for_mapq || !read.is_skipped) continue;
                if (read.mapq < opts.recovery_min_mapq) continue;
                bool in_window = false;
                for (const auto& w : windows)
                    if (read.end >= w.first && read.beg <= w.second) { in_window = true; break; }
                if (!in_window) continue;
                read.is_skipped = false;
                ++woken;
            }
            if (opts.verbose > 0 && woken > 0)
                fprintf(stderr, "[retry] woke %zu read(s) below the mapq floor in %zu window(s)\n",
                        woken, windows.size());

            retry_opts.force_noisy_msa = true;
            retry_opts.skip_noisy_kmeans = false;
            retry_opts.retry_windows = windows;
            collect_var_run_phasing(chunk, retry_opts);

            // Whatever the in-chunk re-solve still could not phase gets its own
            // chunk. Re-derive the windows first: the re-solve closes some, and
            // running a targeted solve over a window that is now phased would
            // re-litigate a settled answer.
            auto residual = collect_unphased_windows(
                chunk, opts.retry_min_unphased_reads, opts.retry_min_window_bp);
            // Both kinds of failure, not just one. A seam between two blocks is
            // a window the solve could not phase ACROSS even though it phased
            // the reads on either side, and no unphased-read window is ever
            // reported there.
            for (const auto& seam : collect_block_seams(chunk)) {
                bool covered = false;
                for (const auto& w : residual)
                    if (seam.first >= w.first && seam.second <= w.second) { covered = true; break; }
                if (!covered) residual.push_back(seam);
            }
            if (!residual.empty())
                recover_windows_with_targeted_solve(chunk, residual, opts, context,
                                                   chunk.region.tid, /*allow_import=*/true,
                                                   /*cache=*/nullptr);
        }
    }

    prune_not_candidate_variants(chunk);
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
        const std::unordered_map<std::string, std::string>& chrom_remap) {
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
                        graph_query_contig, chrom_remap);
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
            graph_query_contig, chrom_remap);
        stitch_chunk_haps(batch.chunks, &opts, pgbam_sidecar.get());
        // Snapshot the labels the solve produced, before the output filters
        // erase the ones gap recovery's link votes are counted from. The gap
        // inventory itself is still derived from the filtered chunks, so this
        // changes what the votes can see and nothing about what is emitted.
        // The second argument was the removed private-whitelist mode's
        // admit-all flag, default false.
        filter_hybrid_reads_by_margin(batch.chunks, opts.min_read_hap_margin,
                                      /*admit_all_in_region=*/false);
        filter_hybrid_small_phase_sets(batch.chunks, opts.min_phase_set_reads);
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
