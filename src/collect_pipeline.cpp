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
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <fstream>
#include <filesystem>
#include <getopt.h>
#include <iostream>
#include <iterator>
#include <limits>
#include <map>
#include <memory>
#include <mutex>
#include <numeric>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <thread>
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
    chunk.phase_sets.assign(chunk.reads.size(), kUnphasedReadPhaseSet);

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

// Return one PhasingChunk per input offset so stitching sees genomic order.
// Run process_chunk on chunks[batch_begin..batch_end) using a thread pool.
// Each worker opens its own BAM/FAI handles.  First exception is rethrown
// after all workers join.
static std::vector<PhasingChunk> collect_chunk_batch_parallel(const Options& opts,
                                                                const std::vector<RegionChunk>& chunks,
                                                                size_t batch_begin,
                                                                size_t batch_end) {
    if (batch_begin > batch_end || batch_end > chunks.size()) {
        throw std::runtime_error("invalid chunk batch range");
    }
    const size_t batch_size = batch_end - batch_begin;
    std::vector<PhasingChunk> result(batch_size);
    if (batch_size == 0) return result;

    const size_t worker_count = std::min<size_t>(static_cast<size_t>(opts.threads), batch_size);
    std::atomic<size_t> next_offset{0};
    std::exception_ptr first_error;
    std::mutex error_mutex;
    std::vector<std::thread> workers;
    workers.reserve(worker_count);

    for (size_t worker_i = 0; worker_i < worker_count; ++worker_i) {
        workers.emplace_back([&]() {
            try {
                WorkerContext context(opts);
                while (true) {
                    const size_t offset = next_offset.fetch_add(1);
                    if (offset >= batch_size) break;
                    result[offset] =
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

        // Stitch overlapping chunks before merging candidate rows and writing.
        std::vector<PhasingChunk> batch = collect_chunk_batch_parallel(
            opts, chunks, batch_begin, batch_end);
        stitch_chunk_haps(batch, &opts, pgbam_sidecar.get());
        CandidateTable variants = merge_chunk_candidates(batch);
        // Two records at one locus cannot assign different alleles to the same
        // haplotype. Make their genotypes complementary, then remove conflicts.
        if (opts.collapse_colocated_alleles) {
            make_colocated_alleles_complementary(variants, opts.min_alt_depth);
            drop_conflicting_haplotype_alleles(variants);
        }
        n_variants += variants.size();
        write_variants_tsv_records(variant_out, header.get(), ref, variants);
        if (!opts.output_vcf.empty()) {
            write_variants_vcf_records(vcf_out, opts, header.get(), ref, variants);
        }
        if (!opts.output_phased_vcf.empty()) {
            write_phased_variants_vcf_records(phased_vcf_out, opts, header.get(), ref, variants);
        }
        if (phased_aln_writer) {
            n_out_aln_reads += static_cast<size_t>(phased_aln_writer->write_chunks(batch));
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

#include "graph_sites.hpp"
#include "graph_query.hpp"

namespace pgphase_collect {

// Recover a seam by running the ordinary alignment pipeline over it, then
// merging the discovered sites and observations into the live graph chunk.
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

/// How far a targeted region reaches past its window, measured in PARENT PHASED
/// SITES rather than base pairs.
///
/// The extension exists so the sub-solve's phase set shares tagged reads with
/// the parent blocks on either side: `select_stitch_orientation` votes on reads
/// carrying BOTH a parent haplotype and a sub-solve haplotype, and the reads
/// inside a window are precisely the ones the parent left unphased. The main
/// chunk loop needs no such padding because a read straddling a chunk boundary
/// is already in both chunks (`initialize_chunk_overlap_state`); a window has no
/// equivalent, so it has to reach out to where the parent did phase.
///
/// A fixed 30 kb was wrong in both directions. Measured over the first 10 Mb of
/// chr20 (9,813 parent phased sites, 41 merged regions), reaching three parent
/// sites on both sides needs a median of 5.6 kb and a p90 of 16.9 kb -- so 30 kb
/// was 5x more than needed in the median case -- while ONE region needed 47.3 kb
/// and therefore had no parent site to vote against at all.
constexpr size_t kTargetedSolveFlankSites = 3;
constexpr hts_pos_t kTargetedSolveFlankMin = 2000;
constexpr hts_pos_t kTargetedSolveFlankMax = 60000;

/// A recovery boundary must be an oriented heterozygote assigned to a real
/// phase set. Both BAM and graph candidates use 0 until phasing assigns an
/// anchor; reads separately use -1 while unphased. Testing `> 0` encodes the
/// candidate contract directly, and unequal haplotype alleles exclude
/// homozygous rows.
static bool is_phase_set_anchor(const CandidateVariant& cand) {
    return cand.phase_set > 0 && cand.hap_alt != cand.hap_ref;
}

/// Return the parent's phased candidate positions in coordinate order.
///
/// Candidate tables are position-sorted by the graph builder and preserve that
/// invariant through recovery merges. Filtering the table preserves order, so
/// this stage needs neither a sort nor an associative container. Co-located
/// rows remain separate because the established recovery behavior measures
/// context in phased candidates rather than distinct coordinates.
static std::vector<hts_pos_t> parent_phased_positions(const PhasingChunk& chunk) {
    std::vector<hts_pos_t> positions;
    positions.reserve(chunk.candidates.size());
    for (const CandidateVariant& cand : chunk.candidates)
        if (is_phase_set_anchor(cand)) positions.push_back(cand.key.sort_pos());
    return positions;
}

/// Return positive-width gaps between neighboring phase-set extents.
///
/// Only oriented heterozygotes with `phase_set > 0` contribute coverage.
/// Unphased candidates use 0, and homozygous rows can carry allele values, so
/// testing the phase-set label or alleles alone would create false boundaries.
///
/// A phase-set label is the coordinate of the variant that began the block.
/// Each accepted candidate therefore covers [phase_set, candidate_position].
/// Candidates are position-sorted, so those interval ends arrive in order even
/// when a later candidate reconnects to an older phase set and moves the start
/// backward.
///
/// `covered` is a flat merge stack of disjoint coordinate extents. A new
/// interval absorbs every component it reaches, then appends once. Each
/// component is pushed and popped at most once, giving amortized O(C) time for C
/// candidates and O(K) contiguous storage for K covered components. A hash
/// table is unnecessary, and gaps cannot be emitted during the first pass
/// because a later interval may bridge a provisional gap.
///
/// Terminal and wholly unanchored regions are intentionally absent because
/// outside-in recovery needs established phase boundaries on both sides.
static std::vector<std::pair<hts_pos_t, hts_pos_t>> collect_phase_set_seams(
        const PhasingChunk& chunk) {
    struct CoveredExtent {
        hts_pos_t beg;
        hts_pos_t end;
    };

    std::vector<CoveredExtent> covered;
    for (const CandidateVariant& cand : chunk.candidates) {
        if (!is_phase_set_anchor(cand)) continue;

        CoveredExtent merged{cand.phase_set, cand.key.sort_pos()};
        while (!covered.empty() && merged.beg <= covered.back().end) {
            merged.beg = std::min(merged.beg, covered.back().beg);
            merged.end = std::max(merged.end, covered.back().end);
            covered.pop_back();
        }
        covered.push_back(merged);
    }

    std::vector<std::pair<hts_pos_t, hts_pos_t>> seams;
    if (covered.size() < 2) return seams;
    seams.reserve(covered.size() - 1);
    for (size_t i = 1; i < covered.size(); ++i)
        seams.emplace_back(covered[i - 1].end, covered[i].beg);
    return seams;
}

// Both the BAM command and exact-row recovery use these longcallD settings.
static void use_longcalld_bam_options(Options& opts) {
    opts.anchored_stage2 = false;
    opts.merge_colocated_msa_alleles = false;
    opts.refresh_msa_observations = false;
    opts.add_unplaced_msa_observations = false;
    opts.upstream_msa_insertion_hp = true;
    opts.phase_set_scoped_clean_rounds = false;
    opts.msa_sites_vote_without_gap_link = true;
    opts.infer_complement_at_multiallelic = true;
    opts.upstream_read_scoring = true;
    opts.upstream_assign_hap = true;
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
/// noisy-k-means default, no recursion, and one thread because the outer graph
/// worker pool already parallelises independent chunks.
static Options targeted_solve_options(const Options& opts) {
    Options sub = opts;
    sub.min_mapq = std::min(opts.min_mapq, opts.recovery_min_mapq);
    sub.skip_noisy_kmeans = false;
    sub.threads = 1;
    sub.verbose = 0;
    return sub;
}

struct TargetedWindowGroup {
    RegionChunk region;
    size_t first_window = 0;
    size_t past_last_window = 0;
};

static hts_pos_t clamp_targeted_flank(hts_pos_t distance) {
    return std::clamp(distance, kTargetedSolveFlankMin, kTargetedSolveFlankMax);
}

/// Add phased context to sorted seams and merge touching solve regions.
///
/// Both inputs are coordinate ordered. A monotone cursor visits each parent
/// candidate at most once, making target construction O(A + W) for A anchors
/// and W seams. Group membership is stored as a half-open index range into
/// `windows`; this avoids a vector allocation
/// and coordinate copies for every merged group.
///
/// The left flank reaches the configured number of phased candidates strictly
/// before the left boundary; the right flank reaches that number at or after
/// the right boundary. This preserves the validated recovery context exactly.
/// If a side has too few candidates, the maximum flank is used. Regions are
/// clamped before grouping so a sub-solve cannot pull reads from the neighboring
/// chunk.
static std::vector<TargetedWindowGroup> build_targeted_groups(
        const std::vector<std::pair<hts_pos_t, hts_pos_t>>& windows,
        int solve_tid,
        const std::vector<hts_pos_t>& parent_sites,
        hts_pos_t chunk_beg,
        hts_pos_t chunk_end) {
    std::vector<TargetedWindowGroup> groups;
    groups.reserve(windows.size());

    const size_t required_sites = kTargetedSolveFlankSites;
    size_t site_i = 0;
    for (size_t window_i = 0; window_i < windows.size(); ++window_i) {
        const auto& window = windows[window_i];

        while (site_i < parent_sites.size() && parent_sites[site_i] < window.first)
            ++site_i;
        const size_t left_end = site_i;

        hts_pos_t left_flank = kTargetedSolveFlankMax;
        if (left_end >= required_sites) {
            left_flank = clamp_targeted_flank(
                window.first - parent_sites[left_end - required_sites]);
        }

        while (site_i < parent_sites.size() && parent_sites[site_i] < window.second)
            ++site_i;
        hts_pos_t right_flank = kTargetedSolveFlankMax;
        if (parent_sites.size() - site_i >= required_sites) {
            right_flank = clamp_targeted_flank(
                parent_sites[site_i + required_sites - 1] - window.second);
        }

        const hts_pos_t beg = std::max(chunk_beg, window.first - left_flank);
        const hts_pos_t end = std::min(chunk_end, window.second + right_flank);
        if (beg >= end) continue;

        if (!groups.empty() && beg <= groups.back().region.end) {
            groups.back().region.end = std::max(groups.back().region.end, end);
            groups.back().past_last_window = window_i + 1;
            continue;
        }

        TargetedWindowGroup group;
        group.region.tid = solve_tid;
        group.region.beg = beg;
        group.region.end = end;
        group.region.chunk_id = -1;
        group.first_window = window_i;
        group.past_last_window = window_i + 1;
        groups.push_back(group);
    }
    return groups;
}

/// Return the seam containing pos, using the group's sorted half-open range.
///
/// Recovery intentionally uses strict seam bounds: the two boundary anchors
/// belong to the established blocks, while only sites inside the gap are new.
static const std::pair<hts_pos_t, hts_pos_t>* find_containing_window(
        const std::vector<std::pair<hts_pos_t, hts_pos_t>>& windows,
        const TargetedWindowGroup& group,
        hts_pos_t pos) {
    const auto first = windows.begin() + static_cast<std::ptrdiff_t>(group.first_window);
    const auto last = windows.begin() + static_cast<std::ptrdiff_t>(group.past_last_window);
    const auto after = std::upper_bound(
        first, last, pos,
        [](hts_pos_t value, const auto& window) { return value < window.first; });
    if (after == first) return nullptr;

    const auto& window = *std::prev(after);
    return pos > window.first && pos < window.second ? &window : nullptr;
}

namespace {

/// Identity of a candidate, for matching the alignment's calls against the
/// catalog's. Same fields exact_comp_var_site compares.
struct CandKey {
    hts_pos_t pos;
    int type;
    int ref_len;
    std::string alt;
    bool operator<(const CandKey& other) const {
        return std::tie(pos, type, ref_len, alt) <
               std::tie(other.pos, other.type, other.ref_len, other.alt);
    }
    bool operator==(const CandKey& other) const {
        return pos == other.pos && type == other.type && ref_len == other.ref_len &&
               alt == other.alt;
    }
};

CandKey cand_key_of(const CandidateVariant& cand) {
    return CandKey{cand.key.sort_pos(), static_cast<int>(cand.key.type), cand.key.ref_len,
                   cand.key.alt};
}

/// One read's allele at each merged candidate, keyed by read name.
using AlleleByCand = std::map<CandKey, std::pair<int, int>>;  // -> (allele, alt_qi)

using CandidateIndex = std::map<CandKey, size_t>;

struct ParentCandidateMatch {
    size_t index;
    bool is_raw;
};

static ParentCandidateMatch find_parent_candidate(
        const CandidateIndex& raw_index,
        const CandidateIndex& sequence_index,
        const CandKey& key,
        size_t missing_index) {
    const auto raw = raw_index.find(key);
    if (raw != raw_index.end()) return ParentCandidateMatch{raw->second, true};

    const auto sequence = sequence_index.find(key);
    if (sequence != sequence_index.end())
        return ParentCandidateMatch{sequence->second, false};
    return ParentCandidateMatch{missing_index, false};
}

struct TransferredCandidate {
    CandidateVariant candidate;
    bool orientation_decided;
    bool flip;
};

}  // namespace

void write_recovery_audit(const std::string& path,
                          const std::vector<RecoveredCandidate>& rows) {
    if (path.empty() || rows.empty()) return;
    static std::mutex audit_mu;
    std::lock_guard<std::mutex> lock(audit_mu);
    const bool fresh = !std::ifstream(path).good();
    std::ofstream out(path, std::ios::app);
    if (!out) return;
    if (fresh)
        out << "POS\tTYPE\tREF_LEN\tALT\tCATEGORY\tWIN_BEG\tWIN_END\tKNOWN_RAW\t"
               "KNOWN_TRANSLATED\tINSIDE_WINDOW\tCATEGORY_OK\tAPPENDED\tMETA_BUILT\t"
               "META_REF\tMETA_ALTS\tALN_VERIFIED\n";
    for (const RecoveredCandidate& r : rows)
        out << r.pos << '\t' << r.type << '\t' << r.ref_len << '\t'
            << (r.alt.empty() ? "." : r.alt) << '\t' << r.category << '\t'
            << r.win_beg << '\t' << r.win_end << '\t' << (r.known_raw ? 1 : 0) << '\t'
            << (r.known_translated ? 1 : 0) << '\t' << (r.inside_window ? 1 : 0) << '\t'
            << (r.category_admitted ? 1 : 0) << '\t' << (r.appended ? 1 : 0) << '\t'
            << (r.meta_built ? 1 : 0) << '\t'
            << (r.meta_ref.empty() ? "." : r.meta_ref) << '\t' << r.meta_alts << '\t'
            << (r.alignment_verified ? 1 : 0) << '\n';
}

// A BAM-verified repeat may replace a catalog context veto only when the BAM
// observations themselves link it to a heterozygote on each side. Use the same
// read count and purity gates as ordinary graph links; the k-means round still
// decides its phase after admission.
struct BamPairLinkCounts {
    int same = 0;
    int cross = 0;
    bool seen_ref = false;
    bool seen_alt = false;
};

static BamPairLinkCounts bam_pair_link_counts(const PhasingChunk& src,
                                               size_t a, size_t b,
                                               size_t tested) {
    BamPairLinkCounts counts;
    for (size_t ri = 0; ri < src.read_var_profile.size(); ++ri) {
        if (src.reads[ri].is_skipped) continue;
        const ReadVariantProfile& prof = src.read_var_profile[ri];
        if (prof.start_var_idx < 0 ||
            a < static_cast<size_t>(prof.start_var_idx) ||
            b < static_cast<size_t>(prof.start_var_idx) ||
            a > static_cast<size_t>(prof.end_var_idx) ||
            b > static_cast<size_t>(prof.end_var_idx)) continue;
        const int aa = prof.alleles[a - static_cast<size_t>(prof.start_var_idx)];
        const int bb = prof.alleles[b - static_cast<size_t>(prof.start_var_idx)];
        if (aa < 0 || aa > 1 || bb < 0 || bb > 1) continue;
        const int tested_allele = tested == a ? aa : bb;
        counts.seen_ref |= tested_allele == 0;
        counts.seen_alt |= tested_allele == 1;
        if (aa == bb) ++counts.same;
        else ++counts.cross;
    }
    return counts;
}

static bool bam_pair_has_pure_link(const PhasingChunk& src, size_t a, size_t b,
                                   size_t tested, const Options& opts) {
    const BamPairLinkCounts counts = bam_pair_link_counts(src, a, b, tested);
    const int total = counts.same + counts.cross;
    return counts.seen_ref && counts.seen_alt &&
           total >= opts.min_block_link_reads &&
           static_cast<double>(std::max(counts.same, counts.cross)) / total >=
               opts.link_earned_min_purity;
}

constexpr int kRecoveryFrontierMinMargin = 10;
constexpr hts_pos_t kRecoveryFrontierMaxStepBp = 10000;

struct RecoveryFrontierScore {
    int margin = 0;
    int total = 0;
    bool decisive = false;
};

static RecoveryFrontierScore recovery_frontier_score(const PhasingChunk& src,
                                                      size_t a, size_t b,
                                                      size_t tested,
                                                      const Options& opts) {
    const hts_pos_t distance =
        std::abs(src.candidates[a].key.sort_pos() - src.candidates[b].key.sort_pos());
    if (distance > kRecoveryFrontierMaxStepBp) return {};
    const BamPairLinkCounts counts = bam_pair_link_counts(src, a, b, tested);
    const int min_margin = std::max(opts.min_block_link_reads,
                                    kRecoveryFrontierMinMargin);
    const int margin = std::abs(counts.same - counts.cross);
    const int total = counts.same + counts.cross;
    return RecoveryFrontierScore{
        margin,
        total,
        counts.seen_ref && counts.seen_alt && margin >= min_margin};
}

static bool bam_site_has_pure_flank_links(const PhasingChunk& src, size_t ci,
                                          const Options& opts) {
    auto is_het = [](const CandidateVariant& cand) {
        const VariantCategory cat = cand.counts.category;
        return cat == VariantCategory::CleanHetSnp ||
               cat == VariantCategory::CleanHetIndel ||
               cat == VariantCategory::NoisyCandHet;
    };
    bool left = false, right = false;
    for (size_t j = ci; j > 0; ) {
        --j;
        if (!is_het(src.candidates[j])) continue;
        left = bam_pair_has_pure_link(src, j, ci, ci, opts);
        break;
    }
    for (size_t j = ci + 1; j < src.candidates.size(); ++j) {
        if (!is_het(src.candidates[j])) continue;
        right = bam_pair_has_pure_link(src, ci, j, ci, opts);
        break;
    }
    return left && right;
}

// Expand only from established phase-block boundaries. One call admits the
// best-supported frontier locus for each disconnected phase-set pair, then the
// caller re-runs phasing before another layer is considered.
size_t expand_recovery_frontiers_once(PhasingChunk& chunk, const Options& opts) {
    const auto is_anchor = [](const CandidateVariant& cand) {
        if (cand.phase_set <= 0 || cand.hap_alt == cand.hap_ref) return false;
        if (cand.bam_injected && !cand.alignment_verified) return false;
        if (cand.is_homopolymer_indel && !cand.gap_link_supported) return false;
        if (!cand.msa_insertion_alts.empty() && !cand.gap_link_supported) return false;
        return true;
    };
    const auto is_frontier_candidate = [](const CandidateVariant& cand) {
        if (!cand.alignment_verified || cand.gap_link_supported) return false;
        if (!cand.is_homopolymer_indel && cand.msa_insertion_alts.empty()) return false;
        return (cand.lcd_var_i_to_cate & kCandGermlineVarCate) != 0;
    };

    // One representative per phased locus is enough to define the two block
    // boundaries. Separate BAM rows remain separate candidates and are tested
    // individually when their locus reaches a frontier.
    std::vector<size_t> anchors;
    for (size_t ci = 0; ci < chunk.candidates.size(); ++ci) {
        if (!is_anchor(chunk.candidates[ci])) continue;
        if (!anchors.empty() &&
            chunk.candidates[anchors.back()].key.sort_pos() ==
                chunk.candidates[ci].key.sort_pos()) continue;
        anchors.push_back(ci);
    }

    struct RankedLocus {
        size_t first = 0;
        size_t last = 0;
        size_t boundary = 0;
        hts_pos_t pos = 0;
        hts_pos_t distance = 0;
        int margin = 0;
        int total = 0;
        bool found = false;
    };

    const auto choose_locus = [&](size_t first, size_t last, size_t boundary,
                                  RankedLocus& best) {
        size_t locus_first = first;
        while (locus_first < last) {
            const hts_pos_t pos = chunk.candidates[locus_first].key.sort_pos();
            size_t locus_last = locus_first + 1;
            while (locus_last < last &&
                   chunk.candidates[locus_last].key.sort_pos() == pos)
                ++locus_last;

            RecoveryFrontierScore locus_score;
            for (size_t ci = locus_first; ci < locus_last; ++ci) {
                if (!is_frontier_candidate(chunk.candidates[ci])) continue;
                // Evidence back to the phase set the row already belongs to
                // confirms that assignment; it does not extend a boundary.
                if (chunk.candidates[ci].phase_set > 0 &&
                    chunk.candidates[ci].phase_set ==
                        chunk.candidates[boundary].phase_set) continue;
                const size_t a = std::min(ci, boundary);
                const size_t b = std::max(ci, boundary);
                const RecoveryFrontierScore score =
                    recovery_frontier_score(chunk, a, b, ci, opts);
                if (!score.decisive) continue;
                if (!locus_score.decisive || score.margin > locus_score.margin ||
                    (score.margin == locus_score.margin && score.total > locus_score.total))
                    locus_score = score;
            }

            const hts_pos_t distance = std::abs(
                pos - chunk.candidates[boundary].key.sort_pos());
            // Rank a locus by its strongest separate BAM row. Summing rows at
            // one coordinate would count the same reads more than once.
            if (locus_score.decisive &&
                (!best.found || locus_score.margin > best.margin ||
                 (locus_score.margin == best.margin && locus_score.total > best.total) ||
                 (locus_score.margin == best.margin && locus_score.total == best.total &&
                  distance < best.distance) ||
                 (locus_score.margin == best.margin && locus_score.total == best.total &&
                  distance == best.distance && pos < best.pos))) {
                best = RankedLocus{locus_first, locus_last, boundary, pos, distance,
                                   locus_score.margin, locus_score.total, true};
            }
            locus_first = locus_last;
        }
    };

    std::vector<RankedLocus> winners;
    for (size_t ai = 1; ai < anchors.size(); ++ai) {
        const size_t left = anchors[ai - 1];
        const size_t right = anchors[ai];
        if (chunk.candidates[left].phase_set == chunk.candidates[right].phase_set)
            continue;
        RankedLocus best;

        // Only the next recovery locus on each side is a frontier. A weak
        // intermediate locus cannot be skipped to reach stronger evidence
        // deeper in the gap.
        size_t left_first = left + 1;
        while (left_first < right &&
               !is_frontier_candidate(chunk.candidates[left_first]))
            ++left_first;
        if (left_first < right) {
            size_t left_last = left_first + 1;
            const hts_pos_t pos = chunk.candidates[left_first].key.sort_pos();
            while (left_last < right &&
                   chunk.candidates[left_last].key.sort_pos() == pos)
                ++left_last;
            choose_locus(left_first, left_last, left, best);
        }

        size_t right_last = right;
        while (right_last > left + 1 &&
               !is_frontier_candidate(chunk.candidates[right_last - 1]))
            --right_last;
        if (right_last > left + 1) {
            size_t right_first = right_last - 1;
            const hts_pos_t pos = chunk.candidates[right_first].key.sort_pos();
            while (right_first > left + 1 &&
                   chunk.candidates[right_first - 1].key.sort_pos() == pos)
                --right_first;
            choose_locus(right_first, right_last, right, best);
        }
        if (best.found) winners.push_back(best);
    }

    if (winners.empty()) return 0;

    // Keep every BAM row at each selected locus separate. A row enters only on
    // its own decisive evidence; a winning row selects the coordinate but does
    // not lend support to its co-located neighbors.
    size_t admitted = 0;
    for (const RankedLocus& best : winners) {
        for (size_t ci = best.first; ci < best.last; ++ci) {
            if (!is_frontier_candidate(chunk.candidates[ci])) continue;
            if (chunk.candidates[ci].phase_set > 0 &&
                chunk.candidates[ci].phase_set ==
                    chunk.candidates[best.boundary].phase_set) continue;
            const size_t a = std::min(ci, best.boundary);
            const size_t b = std::max(ci, best.boundary);
            if (!recovery_frontier_score(chunk, a, b, ci, opts).decisive) continue;
            chunk.candidates[ci].gap_link_supported = true;
            ++admitted;
        }
    }
    return admitted;
}

size_t recover_phase_set_seams_in_place(GraphChunkBuildResult& graph_chunk,
                                         const Options& opts,
                                         WorkerContext& context,
                                         const char* contig_name) {
    PhasingChunk& chunk = graph_chunk.chunk;
    if (chunk.candidates.empty() || chunk.reads.empty()) return 0;

    const int solve_tid = contig_name != nullptr
                              ? sam_hdr_name2tid(context.primary_header(), contig_name)
                              : chunk.region.tid;
    if (solve_tid < 0) return 0;

    std::vector<std::pair<hts_pos_t, hts_pos_t>> windows =
        collect_phase_set_seams(chunk);
    if (windows.empty()) return 0;

    const std::vector<hts_pos_t> parent_sites = parent_phased_positions(chunk);
    const std::vector<TargetedWindowGroup> groups =
        build_targeted_groups(windows, solve_tid, parent_sites,
                              chunk.ref_beg, chunk.ref_end);
    if (groups.empty()) return 0;

    Options sub = targeted_solve_options(opts);
    // Scope the depth-based het escape to the windows actually being recovered.
    // allele_depths_call_het is the one reader of retry_windows, and nothing
    // has filled that list since the post-hoc retry path was removed -- so the
    // escape it guards (admitting a depth-clear homopolymer indel to the link
    // list, and keeping a genuine het from collapsing to 1|1) has been dead in
    // every run. Enabling it chromosome-wide was measured at 1.161% -> 4.502%
    // read hamming with 333 -> 526 blocks, so it is not a latent win: it was
    // tuned for exactly this scoped use, inside a window the first pass could
    // not phase.
    sub.retry_windows.clear();
    sub.retry_windows.reserve(windows.size());
    for (const TargetedWindowGroup& group : groups)
        for (size_t wi = group.first_window; wi < group.past_last_window; ++wi)
            sub.retry_windows.push_back(windows[wi]);
    // The parent re-solves over the merged sites and collapses them unless it
    // applies the same depth-based het test, which is scoped by this list.
    graph_chunk.recovery_windows = sub.retry_windows;
    std::vector<PhasingChunk> discovered;
    discovered.reserve(groups.size());
    for (const TargetedWindowGroup& group : groups)
        discovered.push_back(process_chunk(group.region, sub, context));

    // What the alignment found INSIDE the windows, plus every read's allele at
    // those sites and at the catalog sites the alignment also called. The second
    // part is what lets a read the catalog never saw link across the gap: it
    // carries alleles at sites on both sides, so one solve places it.
    // A catalog candidate is identified by its ALLELE WALK (key.alt is
    // ">114849551>114849554"), while the alignment identifies the same variant by
    // SEQUENCE ("T"). Matching the raw keys therefore never succeeds -- measured
    // on chr20:1-3,000,000, 0 of 198 sub-solve candidates matched a parent key,
    // including 0 of 152 SNPs -- which left the two channels' sites disjoint, so
    // consecutive sites had no read in common and the VCF came out in 4,680
    // blocks instead of 47. The sequence-level identity lives in the parallel
    // site_meta array, and vcf_to_variant_key applies the same anchor trimming
    // inject_graph_sites uses in the other direction.
    CandidateIndex parent_seq_index;
    const bool have_meta = graph_chunk.site_meta.size() == chunk.candidates.size();
    CandidateIndex parent_cand_index;
    for (size_t ci = 0; ci < chunk.candidates.size(); ++ci) {
        const CandKey key = cand_key_of(chunk.candidates[ci]);
        parent_cand_index.emplace(key, ci);
        if (!have_meta) continue;
        const GraphSiteMeta& meta = graph_chunk.site_meta[ci];
        if (meta.ref.empty()) continue;
        for (const std::string& alt : meta.alts) {
            if (alt.empty()) continue;
            const VariantKey translated =
                vcf_to_variant_key(solve_tid, meta.pos, meta.ref, alt);
            const CandKey key{translated.sort_pos(), static_cast<int>(translated.type),
                              translated.ref_len, translated.alt};
            parent_seq_index.emplace(key, ci);
        }
    }

    std::vector<std::set<hts_pos_t>> split_positions(discovered.size());
    for (size_t gi = 0; gi < discovered.size(); ++gi)
        for (const CandidateVariant& cand : discovered[gi].candidates) {
            const CandKey key = cand_key_of(cand);
            const ParentCandidateMatch parent =
                find_parent_candidate(parent_cand_index, parent_seq_index, key,
                                      chunk.candidates.size());
            if (cand.msa_insertion_alts.size() < 2 ||
                parent.index < chunk.candidates.size()) {
                continue;
            }
            if (find_containing_window(windows, groups[gi], key.pos) != nullptr)
                split_positions[gi].insert(key.pos);
        }

    Options exact_sub = sub;
    use_longcalld_bam_options(exact_sub);
    std::vector<PhasingChunk> exact_discovered(groups.size());
    for (size_t gi = 0; gi < groups.size(); ++gi)
        if (!split_positions[gi].empty())
            exact_discovered[gi] =
                process_chunk(groups[gi].region, exact_sub, context);

    size_t adopted = 0;
    // Every sub-solve candidate and what the merge decided about it, so a
    // site that is found and then silently dropped is visible in output
    // rather than only under a probe.
    std::vector<RecoveredCandidate> audit;
    std::map<CandKey, size_t> audit_of;
    std::map<CandKey, TransferredCandidate> new_cands;
    std::map<std::string, AlleleByCand> observed;
    std::map<std::string, int> observed_mapq;
    // Parent haplotype per read, for the orientation vote below. The parent
    // chunk was solved before recovery ran, so these labels are the gauge every
    // imported block has to be expressed in.
    std::map<std::string, int> parent_hap;
    if (opts.stitch_recovered && chunk.haps.size() == chunk.reads.size())
        for (size_t ri = 0; ri < chunk.reads.size(); ++ri)
            if (chunk.haps[ri] != 0) parent_hap.emplace(chunk.reads[ri].qname, chunk.haps[ri]);

    for (size_t gi = 0; gi < discovered.size(); ++gi) {
        const PhasingChunk& src = discovered[gi];
        if (src.read_var_profile.size() != src.reads.size()) continue;
        // Orient this sub-solve against the parent on the reads they share --
        // the same vote select_stitch_orientation runs between chunks. Without
        // it a carried consensus is an arbitrary orientation asserted as fact,
        // which is why pinning alone made things worse.
        bool flip_group = false;
        int vote_same = 0, vote_cross = 0;
        if (opts.stitch_recovered && src.haps.size() == src.reads.size()) {
            for (size_t ri = 0; ri < src.reads.size(); ++ri) {
                if (src.haps[ri] == 0) continue;
                auto ph = parent_hap.find(src.reads[ri].qname);
                if (ph == parent_hap.end()) continue;
                if (src.haps[ri] == ph->second) ++vote_same;
                else ++vote_cross;
            }
            flip_group = vote_cross > vote_same;
            if (getenv("PGPHASE_VOTE") != nullptr)
                fprintf(stderr, "VOTE group=%zu same=%d cross=%d flip=%d shared=%d\n",
                        gi, vote_same, vote_cross, (int)flip_group, vote_same + vote_cross);
        }
        const bool orient_ok = opts.stitch_recovered &&
                               (vote_same + vote_cross) >= 2 &&
                               vote_same != vote_cross;
        const auto containing_window = [&](hts_pos_t pos) {
            return find_containing_window(windows, groups[gi], pos);
        };
        for (size_t ci = 0; ci < src.candidates.size(); ++ci) {
            const CandidateVariant& cand = src.candidates[ci];
            const CandKey key = cand_key_of(cand);
            const ParentCandidateMatch parent =
                find_parent_candidate(parent_cand_index, parent_seq_index, key,
                                      chunk.candidates.size());

            const auto* member = containing_window(key.pos);
            if (!opts.recovery_audit_out.empty()) {
                RecoveredCandidate rec;
                rec.pos = key.pos;
                rec.type = key.type;
                rec.ref_len = key.ref_len;
                rec.alt = key.alt;
                rec.category = static_cast<int>(cand.counts.category);
                rec.known_raw =
                    parent.index < chunk.candidates.size() && parent.is_raw;
                rec.known_translated =
                    parent.index < chunk.candidates.size() && !parent.is_raw;
                rec.inside_window = member != nullptr;
                rec.category_admitted =
                    (cand.lcd_var_i_to_cate & kCandGermlineVarCate) != 0;
                if (member != nullptr) {
                    rec.win_beg = member->first;
                    rec.win_end = member->second;
                }
                const bool inserted = audit_of.emplace(key, audit.size()).second;
                if (inserted) audit.push_back(std::move(rec));
            }

            if (parent.index < chunk.candidates.size()) {
                // A sequence-identical graph row is still alignment recovered:
                // the merge adds BAM observations to the existing row instead
                // of appending a duplicate representation.
                const bool recovered_usable =
                    cand.counts.category == VariantCategory::NoisyCandHet ||
                    cand.counts.category == VariantCategory::CleanHetIndel ||
                    cand.counts.category == VariantCategory::CleanHetSnp;
                if (member != nullptr && recovered_usable) {
                    chunk.candidates[parent.index].alignment_verified =
                        cand.counts.category != VariantCategory::NoisyCandHet ||
                        cand.msa_verified;
                }

                // A reference-context demotion can be superseded only when the
                // BAM reconstruction links cleanly on both sides.
                if (member != nullptr &&
                    bam_site_has_pure_flank_links(src, ci, opts)) {
                    CandidateVariant& parent_cand = chunk.candidates[parent.index];
                    const bool parent_demoted =
                        parent_cand.counts.category == VariantCategory::RepeatHetIndel;
                    const bool sub_usable =
                        cand.counts.category == VariantCategory::NoisyCandHet ||
                        cand.counts.category == VariantCategory::CleanHetIndel ||
                        cand.counts.category == VariantCategory::CleanHetSnp;
                    if (parent_demoted && sub_usable &&
                        (cand.lcd_var_i_to_cate & kCandGermlineVarCate) != 0) {
                        // Keep the sub-solve's own category. Forcing clean
                        // indel measured worse on whole chr20 (2.369% versus
                        // 1.124%) without reducing fragmentation.
                        parent_cand.counts.category = cand.counts.category;
                        parent_cand.counts.candvarcate_initial =
                            cand.counts.category;
                        parent_cand.lcd_var_i_to_cate = cand.lcd_var_i_to_cate;
                        parent_cand.msa_verified = cand.msa_verified;
                        parent_cand.alignment_verified =
                            cand.counts.category != VariantCategory::NoisyCandHet ||
                            cand.msa_verified;
                        if (!cand.msa_insertion_alts.empty())
                            parent_cand.msa_insertion_alts =
                                cand.msa_insertion_alts;
                        ++adopted;
                    }
                }
                continue;
            }
            if (member == nullptr || split_positions[gi].count(key.pos) != 0)
                continue;
            // Admit what the solve itself would admit; the emitter's own
            // category gate runs later and independently.
            if ((cand.lcd_var_i_to_cate & kCandGermlineVarCate) == 0) continue;
            // ... but not a category the writer will never publish. A LOW_COV or
            // LOW_AF merged site still takes part in the solve, receives a phase
            // set and tags reads, while the emitter drops it -- which leaves
            // reads carrying a phase set no record describes. Measured on chr20:
            // all 17 candidates behind the five such phase sets were merged
            // sites, 16 of them LOW_COV or LOW_AF.
            // Deliberately NOT filtered here on what the writer will publish. A
            // site the writer calls LOW_COV or LOW_AF still carries evidence the
            // solve uses, and withholding those costs more than the tidiness it
            // buys: measured on chr20, excluding them removes 3 of the 4 phase
            // sets whose reads no record describes, and costs 11 blocks of
            // contiguity (324 -> 335) and 14 more misplaced reads (2,543 ->
            // 2,557, hamming 1.161% -> 1.168%). Contiguity and accuracy are the
            // deliverables; a read tagged with a phase set the VCF does not
            // describe is a cosmetic inconsistency.
            new_cands.emplace(
                key, TransferredCandidate{cand, orient_ok, flip_group});
            { auto ai = audit_of.find(key);
              if (ai != audit_of.end()) audit[ai->second].appended = true; }
        }
        for (size_t ri = 0; ri < src.reads.size(); ++ri) {
            const ReadVariantProfile& prof = src.read_var_profile[ri];
            if (prof.start_var_idx < 0) continue;
            AlleleByCand& per_read = observed[src.reads[ri].qname];
            observed_mapq[src.reads[ri].qname] = src.reads[ri].mapq;
            for (size_t k = 0; k < prof.alleles.size(); ++k) {
                const size_t ci = static_cast<size_t>(prof.start_var_idx) + k;
                if (ci >= src.candidates.size()) break;
                if (prof.alleles[k] < 0) continue;
                const CandKey observed_key = cand_key_of(src.candidates[ci]);
                const ParentCandidateMatch parent =
                    find_parent_candidate(parent_cand_index, parent_seq_index,
                                          observed_key, chunk.candidates.size());
                if (split_positions[gi].count(observed_key.pos) != 0 &&
                    parent.index == chunk.candidates.size()) {
                    continue;
                }
                per_read.emplace(observed_key,
                                 std::make_pair(prof.alleles[k],
                                                k < prof.alt_qi.size() ? prof.alt_qi[k] : 0));
            }
        }
    }
    // A graph-default MSA may combine co-located BAM alleles. Replace that
    // candidate and its observations with the separate rows produced by the
    // longcallD BAM settings; do not merge the rows again.
    for (size_t gi = 0; gi < exact_discovered.size(); ++gi) {
        const PhasingChunk& src = exact_discovered[gi];
        if (src.read_var_profile.size() != src.reads.size()) continue;
        for (const CandidateVariant& cand : src.candidates) {
            const CandKey key = cand_key_of(cand);
            const ParentCandidateMatch parent =
                find_parent_candidate(parent_cand_index, parent_seq_index, key,
                                      chunk.candidates.size());
            if (split_positions[gi].count(key.pos) == 0 ||
                parent.index < chunk.candidates.size() ||
                (cand.lcd_var_i_to_cate & kCandGermlineVarCate) == 0) {
                continue;
            }
            const auto* member =
                find_containing_window(windows, groups[gi], key.pos);
            if (member == nullptr) continue;
            new_cands.emplace(
                key, TransferredCandidate{cand, false, false});
            if (!opts.recovery_audit_out.empty()) {
                auto ai = audit_of.find(key);
                if (ai == audit_of.end()) {
                    RecoveredCandidate rec;
                    rec.pos = key.pos; rec.type = key.type;
                    rec.ref_len = key.ref_len; rec.alt = key.alt;
                    ai = audit_of.emplace(key, audit.size()).first;
                    audit.push_back(std::move(rec));
                }
                RecoveredCandidate& rec = audit[ai->second];
                rec.category = static_cast<int>(cand.counts.category);
                rec.category_admitted = true;
                rec.inside_window = true;
                rec.appended = true;
                rec.win_beg = member->first;
                rec.win_end = member->second;
            }
        }
        for (size_t ri = 0; ri < src.reads.size(); ++ri) {
            const ReadVariantProfile& prof = src.read_var_profile[ri];
            if (prof.start_var_idx < 0) continue;
            AlleleByCand& per_read = observed[src.reads[ri].qname];
            observed_mapq[src.reads[ri].qname] = src.reads[ri].mapq;
            for (size_t k = 0; k < prof.alleles.size(); ++k) {
                const size_t ci = static_cast<size_t>(prof.start_var_idx) + k;
                if (ci >= src.candidates.size()) break;
                const CandKey key = cand_key_of(src.candidates[ci]);
                if (prof.alleles[k] < 0 || split_positions[gi].count(key.pos) == 0 ||
                    new_cands.find(key) == new_cands.end()) continue;
                per_read.emplace(key, std::make_pair(prof.alleles[k],
                    k < prof.alt_qi.size() ? prof.alt_qi[k] : 0));
            }
        }
    }
    // A seam whose two sides are already called has nothing NEW between them --
    // and bailing here threw away the alignment's read evidence for the sites
    // that are already there, which is the evidence the linker was short of.
    // Measured at the 505 bp break chr20:26,624,496-26,625,001: the chunk holds
    // 13 observations at the left site and NONE at the right, while the BAM has
    // 114 and 100 reads across them. The transfer below already merges an
    // alignment allele onto a candidate the parent owns and adds reads the
    // parent never had; only this early return kept it from running.
    size_t refreshed = 0;
    for (const auto& per_read : observed)
        for (const auto& obs : per_read.second) {
            const ParentCandidateMatch parent =
                find_parent_candidate(parent_cand_index, parent_seq_index,
                                      obs.first, chunk.candidates.size());
            if (parent.index < chunk.candidates.size()) ++refreshed;
        }
    if (new_cands.empty() && refreshed == 0 && adopted == 0) return 0;

    // Re-index. The candidates must stay position-sorted: the solve's outward
    // sweep walks them in INDEX order, so appending in-gap sites at the tail
    // would make them adjacent to the chunk's last site instead of to their
    // positional neighbours. So the merged table is rebuilt in key order and
    // every index-parallel array is rebuilt with it -- the read profiles and the
    // graph arm's site_ids / site_meta / site_allele_orig_idx.
    struct Slot {
        CandKey key;
        long old_index;  // -1 for an alignment candidate
    };
    std::vector<Slot> slots;
    slots.reserve(chunk.candidates.size() + new_cands.size());
    for (size_t ci = 0; ci < chunk.candidates.size(); ++ci)
        slots.push_back(Slot{cand_key_of(chunk.candidates[ci]), static_cast<long>(ci)});
    for (const auto& entry : new_cands) slots.push_back(Slot{entry.first, -1});
    std::stable_sort(slots.begin(), slots.end(),
                     [](const Slot& a, const Slot& b) { return a.key < b.key; });

    CandidateTable merged_cands;
    std::vector<std::string> merged_ids;
    std::vector<GraphSiteMeta> merged_meta;
    std::vector<std::vector<int>> merged_orig;
    merged_cands.reserve(slots.size());
    merged_ids.reserve(slots.size());
    merged_meta.reserve(slots.size());
    merged_orig.reserve(slots.size());
    std::vector<long> old_to_new(chunk.candidates.size(), -1);
    std::map<CandKey, size_t> index_of;
    const bool had_meta = graph_chunk.site_meta.size() == chunk.candidates.size();
    std::map<size_t, CandKey> seq_key_of_parent;
    for (const auto& entry : parent_seq_index) seq_key_of_parent.emplace(entry.second, entry.first);
    for (const Slot& slot : slots) {
        index_of.emplace(slot.key, merged_cands.size());
        if (slot.old_index >= 0) {
            auto seq = seq_key_of_parent.find(static_cast<size_t>(slot.old_index));
            if (seq != seq_key_of_parent.end()) index_of.emplace(seq->second, merged_cands.size());
        }
        if (slot.old_index >= 0) {
            const size_t ci = static_cast<size_t>(slot.old_index);
            old_to_new[ci] = static_cast<long>(merged_cands.size());
            merged_cands.push_back(chunk.candidates[ci]);
            merged_ids.push_back(ci < graph_chunk.site_ids.size() ? graph_chunk.site_ids[ci]
                                                                 : std::string());
            merged_meta.push_back(had_meta ? graph_chunk.site_meta[ci] : GraphSiteMeta{});
            merged_orig.push_back(ci < graph_chunk.site_allele_orig_idx.size()
                                      ? graph_chunk.site_allele_orig_idx[ci]
                                      : std::vector<int>{});
            continue;
        }
        const TransferredCandidate& transferred = new_cands.at(slot.key);
        const CandidateVariant& cand = transferred.candidate;
        // The graph writer emits one record per entry of counts.alle_covs
        // (graph_chunks_to_candidate_table loops new_a = 1 .. alle_covs.size()),
        // but the alignment path reports depth as ref_cov/alt_cov and leaves
        // alle_covs empty on most candidates -- measured, 952 of 1,377 merged
        // sites over chr20:1-5,000,000. Those sites were phased and tagged reads
        // while emitting no record at all. A merged site is biallelic here (one
        // synthesized alt, orig index {0, 1}), so the two depths ARE the allele
        // depth vector.
        CandidateVariant merged_cand = cand;
        // Injected as discovered: the counts and the consensus below are the
        // sub-solve's own, and the flag keeps the parent from re-deriving them.
        merged_cand.bam_injected = true;
        merged_cand.alignment_verified =
            cand.counts.category == VariantCategory::CleanHetSnp ||
            cand.counts.category == VariantCategory::CleanHetIndel ||
            (cand.counts.category == VariantCategory::NoisyCandHet && cand.msa_verified);
        // Carry the sub-solve's consensus THROUGH the orientation, or drop it
        // and let the parent derive one. A consensus without a vote behind it
        // is worse than none.
        const bool o_ok = transferred.orientation_decided;
        const bool o_flip = transferred.flip;
        // Measured, and it did not work: carrying the sub-solve's per-site
        // consensus hurts even when the orientation vote is decisive and says
        // no flip is needed (same=157, cross=14 over the shared reads of
        // chr20:55,290,000-55,380,000), with emitted records falling 54 -> 46.
        // The consensus is defined against the sub-solve's own read set, so as
        // a fixed constraint in the parent it contradicts reads the sub-solve
        // never saw. Only the orientation itself is applied.
        if (o_ok && o_flip) {
            std::swap(merged_cand.hap_to_cons_alle[1], merged_cand.hap_to_cons_alle[2]);
            std::swap(merged_cand.hap_alt, merged_cand.hap_ref);
        }
        if (merged_cand.counts.alle_covs.size() < 2)
            merged_cand.counts.alle_covs = {merged_cand.counts.ref_cov,
                                            merged_cand.counts.alt_cov};
        merged_cands.push_back(std::move(merged_cand));
        merged_ids.push_back(std::string());
        // Synthesized metadata, so the emitter's index-parallel lookup finds an
        // entry for a site the catalog never held. Without it the merged site is
        // skipped at output (graph_collect.cpp:71) even though it phased.
        // VariantKey stores ref_len, not the reference bases, so REF comes from
        // the chunk's own reference slice.
        // The metadata must be in VCF form, because that is what the emitter
        // writes and what vcf_to_variant_key reads back. A VariantKey is the
        // TRIMMED form: an insertion carries ref_len = 0 with the anchor base at
        // key.pos - 1 and alt = the inserted bases only, and a pure deletion
        // carries alt = "". Pairing a reference slice at key.pos with key.alt
        // therefore emitted an insertion as a substitution -- measured over
        // chr20:1-5,000,000, 17 of 260 merged records had REF and ALT sharing no
        // anchor base (REF=G ALT=CTC for an insertion of CTC after G) -- and gave
        // a pure deletion an empty ALT.
        GraphSiteMeta meta;
        meta.chrom = contig_name != nullptr ? contig_name : std::string();
        // Reference bases come from the DISCOVERED alignment chunks, not from the
        // graph chunk: build_graph_chunk never fills chunk.ref_seq (it passes a
        // slice to the noise filter and keeps nothing), while process_chunk does
        // (bam_digar.cpp:1336). Reading the graph chunk's empty slice is why every
        // merged site ended up with an empty REF, which makes the writer skip it
        // (graph_collect.cpp:99) -- so until now the in-gap sites phased reads in
        // the BAM and never appeared in the VCF at all. That is the mechanism
        // behind the 54 phase sets carrying 3,259 tagged reads and no records.
        const auto ref_at = [&](hts_pos_t pos, int len) -> std::string {
            if (len <= 0) return std::string();
            for (const PhasingChunk& src : discovered) {
                if (src.ref_seq.empty()) continue;
                const hts_pos_t off = pos - src.ref_beg;
                if (off < 0) continue;
                if (static_cast<size_t>(off) + static_cast<size_t>(len) > src.ref_seq.size())
                    continue;
                return src.ref_seq.substr(static_cast<size_t>(off), static_cast<size_t>(len));
            }
            return std::string();
        };
        std::string vcf_ref, vcf_alt;
        hts_pos_t vcf_pos = cand.key.pos;
        // Anchoring depends only on the key, so REF and POS are shared by every
        // allele of the site; only the ALT string varies. Building it per allele
        // is what lets a multiallelic candidate survive the merge.
        const auto anchored = [&](const std::string& allele) -> std::string {
            if (cand.key.type == VariantType::Snp) return allele;
            if (cand.key.type == VariantType::Insertion && cand.key.ref_len == 0)
                return ref_at(cand.key.pos - 1, 1) + allele;
            if (cand.key.type == VariantType::Insertion) return allele;
            return ref_at(cand.key.pos - 1, 1) + allele;  // deletion: left anchor
        };
        if (cand.key.type == VariantType::Snp) {
            vcf_ref = ref_at(cand.key.pos, std::max(1, cand.key.ref_len));
        } else if (cand.key.type == VariantType::Insertion && cand.key.ref_len == 0) {
            vcf_pos = cand.key.pos - 1;
            vcf_ref = ref_at(vcf_pos, 1);
        } else if (cand.key.type == VariantType::Insertion) {
            vcf_ref = ref_at(cand.key.pos, cand.key.ref_len);
        } else {  // Deletion: anchor one base to the left so ALT is never empty.
            vcf_pos = cand.key.pos - 1;
            vcf_ref = ref_at(vcf_pos, 1) + ref_at(cand.key.pos, cand.key.ref_len);
        }
        vcf_alt = anchored(cand.key.alt);
        // Emit only what round-trips: if the VCF form does not convert back to
        // the key it came from, the record would describe a different variant
        // than the one that was phased, so the site is dropped instead.
        //
        // A site whose VCF form cannot be built, or that does not convert back to
        // the key it came from, gets an EMPTY meta rather than being skipped:
        // site_meta is addressed BY CANDIDATE INDEX
        // (graph_chunks_to_candidate_table, graph_collect.cpp:68-72), so dropping
        // an entry shifts every later site's metadata onto the wrong candidate.
        // An empty ref makes the writer skip that one record by itself
        // (graph_collect.cpp:99) while the arrays stay parallel.
        bool usable = !vcf_ref.empty() && !vcf_alt.empty();
        if (usable) {
            const VariantKey round_trip =
                vcf_to_variant_key(cand.key.tid, vcf_pos, vcf_ref, vcf_alt);
            usable = round_trip.type == cand.key.type && round_trip.pos == cand.key.pos &&
                     round_trip.ref_len == cand.key.ref_len && round_trip.alt == cand.key.alt;
        }
        // A candidate heterozygous between two ALTERNATE alleles -- hap_to_cons_alle
        // (1,2), no reference reads -- carries both allele strings in
        // msa_insertion_alts, the same ordered ALT list the alignment writer emits
        // (collect_output.cpp:117-123), indexed so entry i is consensus allele
        // i + 1 and matches counts.alle_covs[i + 1]. Keeping only key.alt flattened
        // those sites to one ALT, and the writer then found consensus index 2
        // pointing past meta.alts and dropped the record
        // (graph_collect.cpp:119). Measured in chr20:4,766,928-4,792,960: the two
        // het sites the gap needs -- 4,785,719 ('ATTTT' at 22 reads against a pure
        // 25 bp deletion at 43) and 4,791,668 (16 T at 34 against 17 T at 24) --
        // were both lost this way, leaving the gap with no phased record at all.
        std::vector<std::string> alts;
        if (cand.msa_insertion_alts.size() >= 2) {
            for (const std::string& allele : cand.msa_insertion_alts) {
                const std::string a = anchored(allele);
                if (a.empty()) { alts.clear(); break; }
                alts.push_back(a);
            }
        }
        if (alts.empty() && usable) alts.push_back(vcf_alt);

        if (usable && !alts.empty()) {
            meta.pos = vcf_pos;
            meta.ref = vcf_ref;
            meta.alts = alts;
        }
        merged_meta.push_back(std::move(meta));
        std::vector<int> orig;
        orig.reserve(alts.size() + 1);
        for (int i = 0; i <= static_cast<int>(alts.size()); ++i) orig.push_back(i);
        if (orig.size() < 2) orig = std::vector<int>{0, 1};
        merged_orig.push_back(std::move(orig));
    }

    // Reads the catalog never had. They are the ones that carry the gap: the
    // alignment sees them and the GAF does not, so without them the merged
    // in-gap sites have almost no observing read and each one starts its own
    // phase set. Measured on the first 10 Mb of chr20 with them omitted: 2,844
    // sites merged but the VCF came out in 4,680 blocks instead of 47.
    // Both inputs are qname-ordered. A monotone merge scan finds alignment-only
    // reads without building a second tree containing every existing qname.
    const size_t existing_read_count = chunk.reads.size();
    size_t existing_read_i = 0;
    for (const auto& entry : observed) {
        while (existing_read_i < existing_read_count &&
               chunk.reads[existing_read_i].qname < entry.first) {
            ++existing_read_i;
        }
        if (existing_read_i < existing_read_count &&
            chunk.reads[existing_read_i].qname == entry.first) {
            continue;
        }
        if (entry.second.empty()) continue;
        hts_pos_t beg = std::numeric_limits<hts_pos_t>::max();
        hts_pos_t end = 0;
        for (const auto& obs : entry.second) {
            beg = std::min(beg, obs.first.pos);
            end = std::max(end, obs.first.pos);
        }
        ReadRecord read;
        read.tid = chunk.region.tid;
        read.input_index = 0;
        read.qname = entry.first;
        auto mq = observed_mapq.find(entry.first);
        read.mapq = mq != observed_mapq.end() ? mq->second : opts.min_mapq;
        read.is_skipped = false;
        read.beg = beg;
        read.end = std::max(beg, end);
        chunk.reads.push_back(std::move(read));
        chunk.read_var_profile.push_back(ReadVariantProfile{});
    }

    // Rebuild each read's profile over the new index space, filling in the
    // alignment's alleles where it observed the read.
    std::vector<ReadVariantProfile> merged_profiles;
    merged_profiles.reserve(chunk.reads.size());
    for (size_t ri = 0; ri < chunk.reads.size(); ++ri) {
        const ReadVariantProfile& old_prof = chunk.read_var_profile[ri];
        std::map<size_t, std::pair<int, int>> alleles;
        if (old_prof.start_var_idx >= 0) {
            for (size_t k = 0; k < old_prof.alleles.size(); ++k) {
                const size_t ci = static_cast<size_t>(old_prof.start_var_idx) + k;
                if (ci >= old_to_new.size() || old_to_new[ci] < 0) continue;
                if (old_prof.alleles[k] < 0) continue;
                alleles.emplace(static_cast<size_t>(old_to_new[ci]),
                                std::make_pair(old_prof.alleles[k],
                                               k < old_prof.alt_qi.size() ? old_prof.alt_qi[k] : 0));
            }
        }
        auto it = observed.find(chunk.reads[ri].qname);
        if (it != observed.end())
            for (const auto& entry : it->second) {
                auto idx = index_of.find(entry.first);
                if (idx != index_of.end()) alleles.insert_or_assign(idx->second, entry.second);
            }
        ReadVariantProfile prof;
        prof.read_id = static_cast<int>(ri);
        if (alleles.empty()) {
            merged_profiles.push_back(std::move(prof));
            continue;
        }
        prof.start_var_idx = static_cast<int>(alleles.begin()->first);
        prof.end_var_idx = static_cast<int>(alleles.rbegin()->first);
        const size_t span = static_cast<size_t>(prof.end_var_idx - prof.start_var_idx + 1);
        prof.alleles.assign(span, -1);
        prof.alt_qi.assign(span, 0);
        for (const auto& entry : alleles) {
            const size_t off = entry.first - static_cast<size_t>(prof.start_var_idx);
            prof.alleles[off] = entry.second.first;
            prof.alt_qi[off] = entry.second.second;
        }
        merged_profiles.push_back(std::move(prof));
    }

    const size_t added = new_cands.size();
    chunk.candidates = std::move(merged_cands);
    chunk.read_var_profile = std::move(merged_profiles);
    graph_chunk.site_ids = std::move(merged_ids);
    graph_chunk.site_meta = std::move(merged_meta);
    if (!opts.recovery_audit_out.empty()) {
        for (size_t ci = 0; ci < chunk.candidates.size(); ++ci) {
            auto ai = audit_of.find(cand_key_of(chunk.candidates[ci]));
            if (ai == audit_of.end()) continue;
            RecoveredCandidate& r = audit[ai->second];
            if (ci < graph_chunk.site_meta.size()) {
                r.meta_built = !graph_chunk.site_meta[ci].ref.empty() &&
                               !graph_chunk.site_meta[ci].alts.empty();
                r.meta_ref = graph_chunk.site_meta[ci].ref;
                r.meta_alts = graph_chunk.site_meta[ci].alts.size();
                const CandidateVariant& mc = chunk.candidates[ci];
                r.alignment_verified =
                    mc.alignment_verified ||
                    mc.counts.category == VariantCategory::CleanHetSnp ||
                    mc.counts.category == VariantCategory::CleanHetIndel ||
                    (mc.counts.category == VariantCategory::NoisyCandHet && mc.msa_verified);
            }
        }
        write_recovery_audit(opts.recovery_audit_out, audit);
    }
    graph_chunk.site_allele_orig_idx = std::move(merged_orig);

    // Every read-indexed vector grows with the reads. Appending without this
    // leaves haps and phase_sets short -- caught by verify_chunk_invariants on
    // its first live run, "haps has 1950 entries for 1967 reads" -- and the
    // re-sort below then SKIPS reordering them, because it only moves vectors
    // whose length matches, so a read's label would no longer belong to that
    // read. The resetting solve happens to overwrite both, which is why this
    // stayed invisible; the anchored solve does not.
    if (chunk.haps.size() != chunk.reads.size()) chunk.haps.resize(chunk.reads.size(), 0);
    if (chunk.phase_sets.size() != chunk.reads.size())
        chunk.phase_sets.resize(chunk.reads.size(), kUnphasedReadPhaseSet);

    // Restore the qname ordering of chunk.reads. The cross-chunk stitch pairs
    // reads with a MERGE-JOIN over the two chunks' read vectors
    // (populate_graph_chunk_pair_overlap_impl, graph_bam_adapter.cpp:1112-1120),
    // which walks both with single advancing indices and therefore requires both
    // to be sorted by qname. Appending the alignment-only reads leaves exactly
    // one inversion at the junction -- measured, one out-of-order pair per chunk
    // -- and that is enough for the join to skip every read past it, costing the
    // stitch the overlap evidence it votes on. Everything indexed by read id
    // moves with the reads.
    {
        std::vector<size_t> order(chunk.reads.size());
        for (size_t i = 0; i < order.size(); ++i) order[i] = i;
        std::stable_sort(order.begin(), order.end(), [&](size_t a, size_t b) {
            return chunk.reads[a].qname < chunk.reads[b].qname;
        });
        bool already_sorted = true;
        for (size_t i = 0; i < order.size(); ++i)
            if (order[i] != i) { already_sorted = false; break; }
        if (!already_sorted) {
            std::vector<ReadRecord> reads_sorted;
            std::vector<ReadVariantProfile> profiles_sorted;
            reads_sorted.reserve(order.size());
            profiles_sorted.reserve(order.size());
            const bool have_haps = chunk.haps.size() == chunk.reads.size();
            const bool have_ps = chunk.phase_sets.size() == chunk.reads.size();
            std::vector<int> haps_sorted;
            std::vector<hts_pos_t> ps_sorted;
            if (have_haps) haps_sorted.reserve(order.size());
            if (have_ps) ps_sorted.reserve(order.size());
            for (size_t new_i = 0; new_i < order.size(); ++new_i) {
                const size_t old_i = order[new_i];
                reads_sorted.push_back(std::move(chunk.reads[old_i]));
                ReadVariantProfile prof = std::move(chunk.read_var_profile[old_i]);
                prof.read_id = static_cast<int>(new_i);
                profiles_sorted.push_back(std::move(prof));
                if (have_haps) haps_sorted.push_back(chunk.haps[old_i]);
                if (have_ps) ps_sorted.push_back(chunk.phase_sets[old_i]);
            }
            chunk.reads = std::move(reads_sorted);
            chunk.read_var_profile = std::move(profiles_sorted);
            if (have_haps) chunk.haps = std::move(haps_sorted);
            if (have_ps) chunk.phase_sets = std::move(ps_sorted);
        }
    }

    // The read<->variant interval tree is keyed by CANDIDATE INDEX, and the
    // solve looks every site's reads up through it (collect_phase.cpp:589, 983).
    // Re-indexing the candidates invalidates it, so it has to be rebuilt exactly
    // as build_graph_chunk and the injection path do. Without this the sweep
    // reads the wrong reads for every site and the chunk comes apart: 4,635
    // phase-set blocks over the first 10 Mb of chr20 against a baseline of 47.
    rebuild_read_var_cr(chunk);

    // Fail loudly here rather than as a block count three runs later.
    hts_pos_t region_lo = 0, region_hi = 0;
    for (const TargetedWindowGroup& group : groups) {
        region_lo = region_lo == 0
                        ? group.region.beg
                        : std::min(region_lo, group.region.beg);
        region_hi = std::max(region_hi, group.region.end);
    }
    verify_chunk_invariants(chunk, graph_chunk.site_ids.size(),
                            graph_chunk.site_meta.size(),
                            graph_chunk.site_allele_orig_idx.size(),
                            region_lo, region_hi);

    if (opts.verbose > 0)
        fprintf(stderr, "[in-pass] %zu window(s) -> %zu region(s), merged %zu alignment site(s)\n",
                windows.size(), groups.size(), added);
    // Refreshed evidence on a site the parent already owned is as much a reason
    // for the caller to re-solve as a new site: it is what lets the linker see
    // the alignment's reads at a seam whose two sides were already called.
    return added + (refreshed > 0 ? 1 : 0);
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

    // Use longcallD's BAM phasing and MSA behavior. These assignments override
    // shared defaults used by the graph path; see Options for each gate.
    use_longcalld_bam_options(opts);

    try {
        run_collect_bam_variation(opts);
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << "\n";
        return 1;
    }

    return 0;
}
