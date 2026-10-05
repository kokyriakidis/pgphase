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
#include "collect_phase_noisy.hpp"
#include "collect_phase_pgbam.hpp"
#include "collect_var.hpp"
#include "fisher_exact.hpp"

#include <algorithm>
#include <array>
#include <atomic>
#include <chrono>
#include <cmath>
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
#include <string_view>
#include <thread>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include <htslib/sam.h>

namespace pgphase_collect {

// A clean biallelic site creates a four-point haplotype separation. Require
// six for an unvalidated block so a single clean-site vote cannot emit a tag.
static constexpr int kIndependentBamReadMinHapScoreMargin = 6;

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

// Recover each bounded graph seam with a small standalone BAM solve. The solve
// owns its local HP/PS gauge; the merge imports that state and the graph pass
// later stitches it to the established blocks on read-backed allele evidence.
// Each side gets at least 50 kb and enough span to include three parent
// anchors, capped at 60 kb and clamped to the owning chunk.
constexpr size_t kTargetedSolveFlankSites = 3;
constexpr hts_pos_t kTargetedSolveFlankMin = 50000;
constexpr hts_pos_t kTargetedSolveFlankMax = 60000;

/// Return a candidate's canonical reference position.
///
/// Graph candidates retain the catalog site's internal key, whose position can
/// differ from the minimal VCF/BAM representation. Recovery windows and BAM
/// candidates must use one coordinate system, so catalog rows are translated
/// through their selected allele metadata. Injected BAM rows already use the
/// canonical key and take the fallback.
static hts_pos_t recovery_position(const GraphChunkBuildResult& graph_chunk,
                                   size_t candidate_index) {
    if (candidate_index < graph_chunk.site_meta.size() &&
        candidate_index < graph_chunk.site_allele_orig_idx.size()) {
        const GraphSiteMeta& meta = graph_chunk.site_meta[candidate_index];
        const std::vector<int>& orig =
            graph_chunk.site_allele_orig_idx[candidate_index];
        if (!meta.ref.empty() && orig.size() > 1) {
            const int alt_index = orig[1] - 1;
            if (alt_index >= 0 &&
                alt_index < static_cast<int>(meta.alts.size()) &&
                !meta.alts[static_cast<size_t>(alt_index)].empty()) {
                const std::string& alt =
                    meta.alts[static_cast<size_t>(alt_index)];
                if (meta.ref.size() == alt.size()) {
                    size_t first_change = 0;
                    while (first_change < meta.ref.size() &&
                           std::toupper(static_cast<unsigned char>(
                               meta.ref[first_change])) ==
                           std::toupper(static_cast<unsigned char>(
                               alt[first_change])))
                        ++first_change;
                    if (first_change < meta.ref.size())
                        return meta.pos +
                            static_cast<hts_pos_t>(first_change);
                }
                return vcf_to_variant_key(
                    graph_chunk.chunk.region.tid, meta.pos, meta.ref,
                    alt).sort_pos();
            }
        }
    }
    return graph_chunk.chunk.candidates[candidate_index].key.sort_pos();
}

/// Return the parent's phased candidate positions in reference-coordinate order.
///
/// Catalog rows normally arrive in reference order. Minimal allele
/// normalization can move a row within a repeat, so the uncommon out-of-order
/// case is sorted explicitly before the monotone target-building scan.
static std::vector<hts_pos_t> parent_phased_positions(
        const GraphChunkBuildResult& graph_chunk) {
    const PhasingChunk& chunk = graph_chunk.chunk;
    std::vector<hts_pos_t> positions;
    positions.reserve(chunk.candidates.size());
    for (size_t ci = 0; ci < chunk.candidates.size(); ++ci)
        if (is_phase_set_anchor(chunk.candidates[ci]))
            positions.push_back(recovery_position(graph_chunk, ci));
    if (!std::is_sorted(positions.begin(), positions.end()))
        std::sort(positions.begin(), positions.end());
    return positions;
}

/// Return positive-width gaps between neighboring phased anchor loci.
///
/// Only oriented heterozygotes with `phase_set > 0` contribute boundaries.
/// Unphased candidates use 0, and homozygous rows can carry allele values, so
/// testing the phase-set label or alleles alone would create false boundaries.
///
/// Recovery needs the actual interval between the last heterozygous anchor of
/// one block and the first heterozygous anchor of the next. A phase-set label is
/// only the block's nominal start and can lie inside that interval. Using it as
/// a boundary trimmed useful BAM sites from the solve. The graph candidate key
/// can likewise retain a padded catalog coordinate, so positions are translated
/// to the same minimal reference representation used by BAM candidates.
///
/// Co-located graph rows are treated as one locus. No seam is emitted when two
/// neighboring loci share any phase set, because that common block already
/// connects them. The usual path is a linear scan; normalization displacement
/// triggers one sort to preserve coordinate order.
///
/// Terminal and wholly unanchored regions are intentionally absent because
/// outside-in recovery needs established phase boundaries on both sides.
std::vector<RecoverySeam> collect_phase_set_seams(
        const GraphChunkBuildResult& graph_chunk) {
    struct Anchor {
        hts_pos_t pos;
        hts_pos_t phase_set;
    };

    const PhasingChunk& chunk = graph_chunk.chunk;
    std::vector<Anchor> anchors;
    anchors.reserve(chunk.candidates.size());
    for (size_t ci = 0; ci < chunk.candidates.size(); ++ci) {
        const CandidateVariant& cand = chunk.candidates[ci];
        if (is_phase_set_anchor(cand))
            anchors.push_back(Anchor{recovery_position(graph_chunk, ci),
                                     cand.phase_set});
    }
    if (!std::is_sorted(anchors.begin(), anchors.end(),
                        [](const Anchor& a, const Anchor& b) {
                            return a.pos < b.pos;
                        })) {
        std::stable_sort(anchors.begin(), anchors.end(),
                         [](const Anchor& a, const Anchor& b) {
                             return a.pos < b.pos;
                         });
    }

    std::vector<RecoverySeam> seams;
    std::vector<hts_pos_t> previous_phase_sets;
    hts_pos_t previous_pos = -1;

    size_t ai = 0;
    while (ai < anchors.size()) {
        const hts_pos_t pos = anchors[ai].pos;
        std::vector<hts_pos_t> phase_sets;
        do {
            const hts_pos_t phase_set = anchors[ai].phase_set;
            if (std::find(phase_sets.begin(), phase_sets.end(), phase_set) ==
                phase_sets.end()) {
                phase_sets.push_back(phase_set);
            }
            ++ai;
        } while (ai < anchors.size() && anchors[ai].pos == pos);

        if (!previous_phase_sets.empty()) {
            bool connected = false;
            for (const hts_pos_t phase_set : phase_sets) {
                if (std::find(previous_phase_sets.begin(), previous_phase_sets.end(),
                              phase_set) != previous_phase_sets.end()) {
                    connected = true;
                    break;
                }
            }
            if (!connected && previous_pos < pos)
                seams.push_back(RecoverySeam{previous_pos, pos,
                                              previous_phase_sets.front(),
                                              phase_sets.front()});
        }
        previous_pos = pos;
        previous_phase_sets = std::move(phase_sets);
    }
    return seams;
}

// Both the BAM command and targeted graph recovery use these longcallD settings.
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
/// Sub-solve options: the longcallD candidate and observation behavior, the
/// recovery MAPQ floor, noisy k-means enabled, no recursion, and one thread
/// because the outer graph worker pool already parallelises independent chunks.
static Options targeted_solve_options(const Options& opts) {
    Options sub = opts;
    use_longcalld_bam_options(sub);
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
    bool focused_retry = false;
    // A newly admitted conflicting-majority MSA retry can resolve a compound
    // flank allele. Transfer its complete certified source path, including
    // private rows outside the seam that define that allele's orientation.
    bool preserve_source_flanks = false;
    // A validated focused solve owns this one seam; the original matrix serves
    // the other seams without copying its reads or running a second broad MSA.
    std::optional<size_t> omitted_window;
};

static bool group_owns_window(const TargetedWindowGroup& group, size_t wi) {
    return group.first_window <= wi && wi < group.past_last_window &&
        group.omitted_window != wi;
}

static bool group_omits_position(const std::vector<RecoverySeam>& windows,
                                 const TargetedWindowGroup& group, hts_pos_t pos) {
    return group.omitted_window &&
        windows[*group.omitted_window].beg <= pos &&
        pos <= windows[*group.omitted_window].end;
}

static bool source_seam_has_msa_conflict(const PhasingChunk& source,
                                        const RecoverySeam& seam,
                                        const Options& opts) {
    std::vector<size_t> sites;
    for (size_t ci = 0; ci < source.candidates.size(); ++ci) {
        const auto& candidate = source.candidates[ci];
        const hts_pos_t pos = candidate.key.sort_pos();
        if (pos < seam.beg || pos > seam.end + 1 ||
            !is_phase_set_anchor(candidate)) continue;
        sites.push_back(ci);
    }
    // Insertion anchors precede their key position. Sort indices rather than
    // rows, so the profiles retain their original candidate offsets.
    std::stable_sort(sites.begin(), sites.end(), [&source](size_t a, size_t b) {
        return source.candidates[a].key.sort_pos() < source.candidates[b].key.sort_pos();
    });
    for (size_t si = 1; si < sites.size(); ++si)
        if (msa_source_conflict_is_supported(source, {sites[si - 1], sites[si]}, opts))
            return true;
    return false;
}

// A second MSA solve is admitted for a weak internal edge or two-block boundary with
// missing, dropped-out, or conflicting calls on crossing molecules. Before
// CIGAR backfill, the first pass measures original MSA dropout; after it,
// the ordinary pass measures the remaining paired-call defects.
// Sparse conflicting pairs and deep indel dropout name one seam for a focused
// solve; other signals retain the grouped solve and its source-row guard.
static bool source_seam_needs_unplaced_msa(
        const PhasingChunk& source, const TargetedWindowGroup& group,
        const std::vector<RecoverySeam>& windows, const Options& opts,
        bool& preserve_source_rows, bool& ordinary_retry,
        std::optional<size_t>& sparse_window,
        bool only_original_dropout = false,
        std::optional<size_t> only_window = std::nullopt,
        bool allow_complementary = false,
        bool* conflicting_majority = nullptr,
        bool allow_nonunique_focused_retry = true,
        bool* internal_conflict_retry = nullptr) {
    if (conflicting_majority != nullptr) *conflicting_majority = false;
    if (internal_conflict_retry != nullptr) *internal_conflict_retry = false;
    preserve_source_rows = false;
    ordinary_retry = false;
    sparse_window.reset();
    if (source.read_var_profile.size() != source.reads.size()) return false;
    constexpr int kAdmissionMinMapq = 30;
    constexpr int kUnknownMapq = 255;
    constexpr int kAdmissionMinCrossingReads = 20;
    constexpr double kMissingPairMaxP = 0.01;
    constexpr int kDropoutMinPairedCalls = 6;
    constexpr int kDropoutMinOtherAlleleCalls = 2;
    constexpr double kDropoutMaxP = 0.05;
    constexpr int kSparseMinPairedCalls = 2;
    constexpr double kDecisivePairMaxP = 0.05;
    constexpr int kMixedParityMinPairedCalls = 6;
    constexpr int kMixedParityMinOpposingCalls = 2;
    constexpr double kMixedParityErrorRate = 0.05;
    constexpr double kMixedParityMaxP = 0.01;
    struct SourceRun {
        hts_pos_t phase_set;
        size_t first;
        size_t last;
    };
    const size_t first = only_window.value_or(group.first_window);
    const size_t past_last = only_window ? *only_window + 1 : group.past_last_window;
    for (size_t wi = first; wi < past_last; ++wi) {
        if (!group_owns_window(group, wi)) continue;
        const RecoverySeam& seam = windows[wi];
        if (source_seam_has_msa_conflict(source, seam, opts)) {
            preserve_source_rows = true;
            sparse_window = wi;
            if (internal_conflict_retry != nullptr) *internal_conflict_retry = true;
            return true;
        }
        std::vector<SourceRun> runs;
        for (size_t ci = 0; ci < source.candidates.size(); ++ci) {
            const CandidateVariant& candidate = source.candidates[ci];
            const hts_pos_t pos = candidate.key.sort_pos();
            // Include both graph boundaries. A BAM insertion key can lie
            // one base beyond the right boundary's VCF anchor.
            if (pos < seam.beg || pos > seam.end + 1 ||
                !is_phase_set_anchor(candidate))
                continue;
            if (runs.empty() || runs.back().phase_set != candidate.phase_set)
                runs.push_back({candidate.phase_set, ci, ci});
            else
                runs.back().last = ci;
            if (runs.size() > 2) break;
        }
        if (runs.size() != 2) continue;
        const size_t left_i = runs.front().last;
        const size_t right_i = runs.back().first;
        const hts_pos_t left_pos = source.candidates[left_i].key.sort_pos();
        const hts_pos_t right_pos = source.candidates[right_i].key.sort_pos();
        if (left_pos >= right_pos) continue;
        // Two co-located, complementary phased rows can request a focused
        // retry without merging their allele representations. Grouped retries
        // still require unique boundaries, and the focused retry must retain
        // every original phased row before it replaces the source solve.
        const std::array<size_t, 2> boundaries{left_i, right_i};
        const bool mixed_msa_pair = only_original_dropout &&
            msa_boundary_dropout_is_supported(source, boundaries, opts,
                                               allow_complementary);
        const auto unique_boundary = [&source, mixed_msa_pair](hts_pos_t pos) {
            return std::count_if(source.candidates.begin(),
                                 source.candidates.end(),
                                 [pos, mixed_msa_pair](const CandidateVariant& candidate) {
                                     return (!mixed_msa_pair || is_phase_set_anchor(candidate)) &&
                                         candidate.key.sort_pos() == pos;
                                 }) == 1;
        };
        const bool left_unique = unique_boundary(left_pos);
        const bool right_unique = unique_boundary(right_pos);
        const bool unique_pair = left_unique && right_unique;
        // A single seam needs the same full adjacent-block context as an
        // outer seam of a group. Its complementary rows must not prevent the
        // focused retry merely because no other padded seam touches it.
        const bool can_focus = group.focused_retry || mixed_msa_pair ||
            wi == group.first_window || wi + 1 == group.past_last_window;
        const auto complementary_boundary = [&](size_t selected) {
            const CandidateVariant& first = source.candidates[selected];
            const CandidateVariant* other = nullptr;
            for (size_t ci = 0; ci < source.candidates.size(); ++ci) {
                if (ci == selected ||
                    (mixed_msa_pair && !is_phase_set_anchor(source.candidates[ci])) ||
                    source.candidates[ci].key.sort_pos() !=
                        first.key.sort_pos())
                    continue;
                if (other != nullptr) return false;
                other = &source.candidates[ci];
            }
            return other != nullptr &&
                first.phase_set == other->phase_set &&
                (first.key.type == other->key.type ||
                 (mixed_msa_pair && first.msa_verified && other->msa_verified &&
                  first.key.type != VariantType::Snp &&
                  other->key.type != VariantType::Snp)) &&
                first.hap_to_cons_alle[1] >= 0 &&
                first.hap_to_cons_alle[1] <= 1 &&
                first.hap_to_cons_alle[2] >= 0 &&
                first.hap_to_cons_alle[2] <= 1 &&
                first.hap_to_cons_alle[1] != first.hap_to_cons_alle[2] &&
                first.hap_to_cons_alle[1] == other->hap_to_cons_alle[2] &&
                first.hap_to_cons_alle[2] == other->hap_to_cons_alle[1];
        };
        if (!unique_pair &&
            (!allow_nonunique_focused_retry || !can_focus ||
             (!left_unique && !complementary_boundary(left_i)) ||
             (!right_unique && !complementary_boundary(right_i))))
            continue;
        const bool left_msa_snp = source.candidates[left_i].msa_verified &&
            source.candidates[left_i].key.type == VariantType::Snp;
        const bool right_msa_snp = source.candidates[right_i].msa_verified &&
            source.candidates[right_i].key.type == VariantType::Snp;
        if (only_original_dropout && !source.candidates[left_i].msa_verified &&
            !source.candidates[right_i].msa_verified) continue;
        int snp_calls[2] = {0, 0};
        int spanning = 0;
        int callable = 0;
        int left_alleles[2] = {0, 0};
        int right_alleles[2] = {0, 0};
        int same_parity = 0;
        int cross_parity = 0;
        for (size_t ri = 0; ri < source.reads.size(); ++ri) {
            const ReadRecord& read = source.reads[ri];
            if (read.is_skipped || read.mapq < kAdmissionMinMapq ||
                (only_original_dropout && read.mapq == kUnknownMapq) ||
                read.beg > left_pos || read.end < right_pos)
                continue;
            ++spanning;
            const ReadVariantProfile& profile = source.read_var_profile[ri];
            if (only_original_dropout && profile.start_var_idx >= 0) {
                for (size_t side = 0; side < boundaries.size(); ++side) {
                    if (boundaries[side] <
                        static_cast<size_t>(profile.start_var_idx))
                        continue;
                    const size_t offset = boundaries[side] -
                        static_cast<size_t>(profile.start_var_idx);
                    // Any known allele is a call, including a non-reference
                    // allele at a multiallelic SNP. Only absence is dropout.
                    if (offset < profile.alleles.size() &&
                        profile.alleles[offset] >= 0)
                        ++snp_calls[side];
                }
            }
            if (profile.start_var_idx < 0 ||
                left_i < static_cast<size_t>(profile.start_var_idx) ||
                right_i > static_cast<size_t>(profile.end_var_idx))
                continue;
            const size_t left_offset =
                left_i - static_cast<size_t>(profile.start_var_idx);
            const size_t right_offset =
                right_i - static_cast<size_t>(profile.start_var_idx);
            if (right_offset < profile.alleles.size() &&
                profile.alleles[left_offset] >= 0 &&
                profile.alleles[right_offset] >= 0) {
                ++callable;
                const int left_allele = profile.alleles[left_offset];
                const int right_allele = profile.alleles[right_offset];
                if (left_allele <= 1 && right_allele <= 1) {
                    ++left_alleles[left_allele];
                    ++right_alleles[right_allele];
                    if (left_allele == right_allele)
                        ++same_parity;
                    else
                        ++cross_parity;
                }
            }
        }
        if (only_original_dropout) {
            if (mixed_msa_pair && can_focus) {
                preserve_source_rows = true;
                sparse_window = wi;
                return true;
            }
            // Backfill can supply these calls only after the source HP/PS solve.
            // Test the original SNP matrix so a repaired call cannot conceal
            // an under-supported consensus. Correct for the two endpoints
            // when both boundaries are MSA SNPs.
            const int missing = std::max(left_msa_snp ? spanning - snp_calls[0] : 0,
                                         right_msa_snp ? spanning - snp_calls[1] : 0);
            if (unique_pair &&
                spanning >= std::max(opts.min_depth, kAdmissionMinCrossingReads) &&
                missing > spanning - missing &&
                (static_cast<int>(left_msa_snp) + static_cast<int>(right_msa_snp)) *
                    binomial_upper_tail(spanning, missing, 0.5) <= kMissingPairMaxP) {
                preserve_source_rows = true;
                ordinary_retry = true;
                sparse_window = wi;
                return true;
            }
            continue;
        }
        // A monomorphic indel call on molecules carrying both alleles of
        // the other boundary indicates allele dropout in the BAM profiles.
        // Use the exact two-sided tail for a balanced heterozygote to avoid
        // retrying on a small random sample of one haplotype.
        const bool left_indel =
            source.candidates[left_i].key.type != VariantType::Snp;
        const bool right_indel =
            source.candidates[right_i].key.type != VariantType::Snp;
        const int paired_biallelic = left_alleles[0] + left_alleles[1];
        const auto dropout = [&](const int affected[2], const int other[2]) {
            return std::min(affected[0], affected[1]) == 0 &&
                std::min(other[0], other[1]) >= kDropoutMinOtherAlleleCalls;
        };
        if (unique_pair &&
            paired_biallelic >= std::max(opts.min_depth, kDropoutMinPairedCalls) &&
            std::ldexp(2.0, -paired_biallelic) <= kDropoutMaxP &&
            ((left_indel && dropout(left_alleles, right_alleles)) ||
             (right_indel && dropout(right_alleles, left_alleles)))) {
            preserve_source_rows = true;
            ordinary_retry = true;
        }
        // These source blocks are still split. A decisive majority cannot
        // conceal conflicting indel calls that exceed the sequencing-error
        // model: MSA may repair the allele representation before stitching.
        // Admission supplies no phase parity; transfer validates the new path.
        const int opposing = std::min(same_parity, cross_parity);
        // Test the observed majority, not the probability of unanimous calls:
        // 6:1 votes have a 0.0625 tail under 50:50, while 7:0 has 0.0078125.
        // Retry just this seam if a grouped MSA would disturb other sites;
        // unanimously oriented sparse calls do not need a retry.
        if ((left_indel || right_indel) && same_parity > 0 &&
            cross_parity > 0 &&
            paired_biallelic >= kSparseMinPairedCalls &&
            binomial_upper_tail(paired_biallelic,
                                std::max(same_parity, cross_parity), 0.5) >
                kDecisivePairMaxP) {
            preserve_source_rows = true;
            if (!sparse_window) sparse_window = wi;
        }
        // A strong majority is an existing link signal. Revisit it only for
        // a compound flank representation: a verified phased MSA deletion
        // covers the graph anchor, but the BAM has no SNP row for that base.
        // Unplaced consensus reads can recover that missing allele context;
        // ordinary noisy indel/SNP majorities need no new MSA solve.
        const auto compound_flank = [&](hts_pos_t pos) {
            if (std::any_of(source.candidates.begin(), source.candidates.end(),
                    [pos](const CandidateVariant& candidate) {
                        return candidate.key.type == VariantType::Snp && candidate.key.pos == pos;
                    })) return false;
            return std::any_of(source.candidates.begin(), source.candidates.end(),
                [pos](const CandidateVariant& candidate) {
                    return candidate.msa_verified && is_phase_set_anchor(candidate) &&
                        candidate.key.type == VariantType::Deletion &&
                        candidate.key.pos <= pos &&
                        pos < candidate.key.pos + candidate.key.ref_len;
                });
        };
        const bool decisive_majority = binomial_upper_tail(paired_biallelic,
            std::max(same_parity, cross_parity), 0.5) <= kDecisivePairMaxP;
        if (unique_pair && (left_indel || right_indel) &&
            (!decisive_majority || compound_flank(seam.beg) || compound_flank(seam.end)) &&
            paired_biallelic >=
                std::max(opts.min_depth, kMixedParityMinPairedCalls) &&
            opposing >= kMixedParityMinOpposingCalls &&
            binomial_upper_tail(paired_biallelic, opposing,
                                kMixedParityErrorRate) <= kMixedParityMaxP) {
            preserve_source_rows = true;
            ordinary_retry = true;
            if (conflicting_majority != nullptr &&
                decisive_majority)
                *conflicting_majority = true;
        }
        if (!unique_pair) continue;
        if (spanning < std::max(opts.min_depth, kAdmissionMinCrossingReads) ||
            spanning - callable <= callable)
            continue;
        // Exact upper binomial tail for missing calls under a 50:50 null.
        // This is an admission trigger, not evidence for joining the blocks.
        const int missing = spanning - callable;
        const double tail = binomial_upper_tail(spanning, missing, 0.5);
        if (tail <= kMissingPairMaxP) {
            if (!opts.phase_matrix_dump_prefix.empty())
                std::fprintf(stderr,
                             "[recovery-msa-admit] %" PRId64 "-%" PRId64
                             " pair=%" PRId64 "-%" PRId64
                             " spanning=%d callable=%d p=%.4g\n",
                             static_cast<int64_t>(seam.beg),
                             static_cast<int64_t>(seam.end),
                             static_cast<int64_t>(left_pos),
                             static_cast<int64_t>(right_pos),
                             spanning, callable, tail);
            preserve_source_rows = left_pos <= seam.beg ||
                                   right_pos >= seam.end;
            ordinary_retry = true;
            return true;
        }
    }
    return ordinary_retry || sparse_window.has_value();
}

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
/// the right boundary. Each side is then clamped to 50--60 kb. If a side has
/// too few candidates, the maximum flank is used. Regions are
/// clamped before grouping so a sub-solve cannot pull reads from the neighboring
/// chunk.
static std::vector<TargetedWindowGroup> build_targeted_groups(
        const std::vector<RecoverySeam>& windows,
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

        while (site_i < parent_sites.size() && parent_sites[site_i] < window.beg)
            ++site_i;
        const size_t left_end = site_i;

        hts_pos_t left_flank = kTargetedSolveFlankMax;
        if (left_end >= required_sites) {
            left_flank = clamp_targeted_flank(
                window.beg - parent_sites[left_end - required_sites]);
        }

        while (site_i < parent_sites.size() && parent_sites[site_i] < window.end)
            ++site_i;
        hts_pos_t right_flank = kTargetedSolveFlankMax;
        if (parent_sites.size() - site_i >= required_sites) {
            right_flank = clamp_targeted_flank(
                parent_sites[site_i + required_sites - 1] - window.end);
        }

        const hts_pos_t beg = std::max(chunk_beg, window.beg - left_flank);
        const hts_pos_t end = std::min(chunk_end, window.end + right_flank);
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
static const RecoverySeam* find_containing_window(
        const std::vector<RecoverySeam>& windows,
        const TargetedWindowGroup& group,
        hts_pos_t pos) {
    const auto first = windows.begin() + static_cast<std::ptrdiff_t>(group.first_window);
    const auto last = windows.begin() + static_cast<std::ptrdiff_t>(group.past_last_window);
    const auto after = std::upper_bound(
        first, last, pos,
        [](hts_pos_t value, const RecoverySeam& window) { return value < window.beg; });
    if (after == first) return nullptr;

    const auto selected = std::prev(after);
    if (!group_owns_window(group, static_cast<size_t>(selected - windows.begin())))
        return nullptr;
    const auto& window = *selected;
    return pos > window.beg && pos < window.end ? &window : nullptr;
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

// A targeted BAM homozygote can veto a graph-only seam anchor only when
// molecules assigned to both BAM haplotypes independently carry that allele.
static bool supported_bam_homozygous_snp(
        const PhasingChunk& source, size_t candidate_i) {
    const CandidateVariant& cand = source.candidates[candidate_i];
    const bool verified_homozygote =
        cand.counts.category == VariantCategory::CleanHom ||
        (cand.counts.category == VariantCategory::NoisyCandHom &&
         cand.msa_verified);
    if (cand.key.type != VariantType::Snp || !verified_homozygote ||
        cand.hap_to_cons_alle[1] < 0 ||
        cand.hap_to_cons_alle[1] != cand.hap_to_cons_alle[2])
        return false;
    const int allele = cand.hap_to_cons_alle[1];
    constexpr int kHomozygousAnchorMinMapq = 30;
    std::array<int, 2> support{};
    for (size_t ri = 0; ri < source.reads.size(); ++ri) {
        if (ri >= source.haps.size() ||
            ri >= source.read_var_profile.size() ||
            source.reads[ri].mapq < kHomozygousAnchorMinMapq ||
            source.haps[ri] < 1 || source.haps[ri] > 2)
            continue;
        const ReadVariantProfile& profile = source.read_var_profile[ri];
        if (profile.start_var_idx < 0 ||
            static_cast<int>(candidate_i) < profile.start_var_idx ||
            static_cast<int>(candidate_i) > profile.end_var_idx)
            continue;
        const size_t offset = candidate_i -
            static_cast<size_t>(profile.start_var_idx);
        if (offset < profile.alleles.size() &&
            profile.alleles[offset] == allele)
            ++support[static_cast<size_t>(source.haps[ri] - 1)];
    }
    // Under a true diploid heterozygote, the chance of observing only one
    // allele this many times is two-sided Binomial(n, 0.5).
    constexpr double kHomozygousAnchorPValue = 0.01;
    return support[0] > 0 && support[1] > 0 &&
        std::ldexp(2.0, -(support[0] + support[1])) <=
            kHomozygousAnchorPValue;
}

struct TransferredCandidate {
    CandidateVariant candidate;
    size_t source_id = 0;
    hts_pos_t source_phase_set = 0;
};

struct TransferredReadPhase {
    int hap = 0;
    hts_pos_t phase_set = 0;
    size_t source_id = 0;
};

} // namespace

namespace {

static int physical_snp_call(const ReadRecord& read, hts_pos_t pos,
                             char ref_base, char alt_base,
                             int* base_quality = nullptr) {
    return read.is_skipped ? -1 : pgphase_collect::physical_snp_call(
        read.alignment.get(), pos, ref_base, alt_base, base_quality);
}

// A sparse graph seam can use a physical SNP molecule only when the same
// molecule also confirms the phase of a neighboring SNP in either block.
// Otherwise one incorrectly oriented endpoint can flip an entire block.
static std::optional<bool> physical_graph_snp_bridge(
        const PhasingChunk& source, size_t left_i, size_t right_i,
        int left_graph_hap1, int right_graph_hap1,
        const std::vector<std::pair<size_t, int>>& left_sites,
        const std::vector<std::pair<size_t, int>>& right_sites,
        double max_wrong_parity_probability) {
    constexpr int kMinBridgeMapq = 30;
    constexpr int kUnknownMapq = 255;
    constexpr int kMinBridgeBaseq = 30;
    constexpr int kMissingBaseQuality = 255;
    constexpr hts_pos_t kMinCorroboratingSnpSpacing = 100;
    const CandidateVariant& left = source.candidates[left_i];
    const CandidateVariant& right = source.candidates[right_i];
    if (left.key.type != VariantType::Snp ||
        right.key.type != VariantType::Snp ||
        left.counts.category != VariantCategory::CleanHetSnp ||
        right.counts.category != VariantCategory::CleanHetSnp ||
        left.key.alt.size() != 1 || right.key.alt.size() != 1 ||
        left.key.pos >= right.key.pos ||
        left.key.pos < source.ref_beg || right.key.pos < source.ref_beg ||
        left_graph_hap1 < 0 || left_graph_hap1 > 1 ||
        right_graph_hap1 < 0 || right_graph_hap1 > 1)
        return std::nullopt;
    const size_t left_offset = static_cast<size_t>(
        left.key.pos - source.ref_beg);
    const size_t right_offset = static_cast<size_t>(
        right.key.pos - source.ref_beg);
    if (left_offset >= source.ref_seq.size() ||
        right_offset >= source.ref_seq.size())
        return std::nullopt;
    const char left_ref = source.ref_seq[left_offset];
    const char right_ref = source.ref_seq[right_offset];
    std::optional<bool> flip;
    std::optional<bool> observed_relation;
    std::unordered_set<std::string> counted_reads;
    double wrong_parity_bound = 1.0;
    for (const ReadRecord& read : source.reads) {
        if (read.is_skipped || read.mapq < kMinBridgeMapq ||
            read.mapq == kUnknownMapq ||
            read.beg > left.key.pos || read.end < right.key.pos)
            continue;
        int left_quality = kMissingBaseQuality;
        int right_quality = kMissingBaseQuality;
        const int left_call = physical_snp_call(
            read, left.key.pos, left_ref, left.key.alt[0], &left_quality);
        const int right_call = physical_snp_call(
            read, right.key.pos, right_ref, right.key.alt[0], &right_quality);
        if ((left_call != 0 && left_call != 2) ||
            (right_call != 0 && right_call != 2) ||
            left_quality < kMinBridgeBaseq ||
            right_quality < kMinBridgeBaseq ||
            left_quality == kMissingBaseQuality ||
            right_quality == kMissingBaseQuality)
            continue;
        const bool left_hap1 =
            (left_call == 2) == (left_graph_hap1 == 1);
        const bool right_hap1 =
            (right_call == 2) == (right_graph_hap1 == 1);
        const bool relation = left_hap1 != right_hap1;
        // A callable contradiction vetoes the edge even when that read lacks
        // enough neighboring SNPs to supply positive block corroboration.
        if (observed_relation && *observed_relation != relation)
            return std::nullopt;
        observed_relation = relation;
        const auto corroborates = [&](const std::vector<std::pair<size_t, int>>& sites,
                                      size_t selected_i,
                                      bool expected_hap1) -> std::optional<bool> {
            bool extra = false;
            const hts_pos_t selected_pos =
                source.candidates[selected_i].key.pos;
            for (const auto& [candidate_i, graph_hap1] : sites) {
                if (candidate_i == selected_i) continue;
                const CandidateVariant& site = source.candidates[candidate_i];
                if (site.key.alt.size() != 1 ||
                    site.key.pos < source.ref_beg ||
                    read.beg > site.key.pos || read.end < site.key.pos)
                    continue;
                const size_t offset = static_cast<size_t>(
                    site.key.pos - source.ref_beg);
                if (offset >= source.ref_seq.size()) continue;
                int quality = kMissingBaseQuality;
                const int call = physical_snp_call(
                    read, site.key.pos, source.ref_seq[offset],
                    site.key.alt[0], &quality);
                if ((call != 0 && call != 2) ||
                    quality < kMinBridgeBaseq ||
                    quality == kMissingBaseQuality)
                    continue;
                if (((call == 2) == (graph_hap1 == 1)) != expected_hap1)
                    return std::nullopt;
                if (std::llabs(site.key.pos - selected_pos) >=
                    kMinCorroboratingSnpSpacing)
                    extra = true;
            }
            return extra;
        };
        const std::optional<bool> left_extra =
            corroborates(left_sites, left_i, left_hap1);
        const std::optional<bool> right_extra =
            corroborates(right_sites, right_i, right_hap1);
        if (!left_extra || !right_extra) return std::nullopt;
        if (!*left_extra && !*right_extra) continue;
        if (!counted_reads.insert(read.qname).second) continue;
        flip = relation;
        const auto error_probability = [](int quality) {
            return std::pow(10.0, -static_cast<double>(quality) / 10.0);
        };
        wrong_parity_bound *=
            error_probability(left_quality) +
            error_probability(right_quality) +
            error_probability(read.mapq);
    }
    return flip && wrong_parity_bound <= max_wrong_parity_probability
        ? flip : std::nullopt;
}

struct SourcePathEvidence {
    size_t site_count = 0;
    std::vector<hts_pos_t> weak_cuts;
    std::vector<hts_pos_t> quality_supported_cuts;
};

// Two independently sequenced molecules can validate a clean-SNP cut even
// when both carry the same haplotype. Bound the chance that both observed
// parities are wrong by the product of their base and mapping error bounds.
static bool high_quality_snp_cut_support(const PhasingChunk& source,
                                         size_t left_i, size_t right_i,
                                         hts_pos_t source_ps) {
    constexpr int kMinBridgeMapq = 30;
    constexpr int kUnknownMapq = 255;
    constexpr int kMinBridgeBaseq = 20;
    constexpr int kMissingBaseQuality = 255;
    constexpr int kMinIndependentMolecules = 2;
    constexpr double kMaxWrongParityProbability = 0.01;
    const CandidateVariant& left = source.candidates[left_i];
    const CandidateVariant& right = source.candidates[right_i];
    if (left.key.type != VariantType::Snp ||
        right.key.type != VariantType::Snp ||
        left.counts.category != VariantCategory::CleanHetSnp ||
        right.counts.category != VariantCategory::CleanHetSnp ||
        left.key.alt.size() != 1 || right.key.alt.size() != 1 ||
        left.key.pos < source.ref_beg || right.key.pos < source.ref_beg)
        return false;
    const size_t left_offset = static_cast<size_t>(
        left.key.pos - source.ref_beg);
    const size_t right_offset = static_cast<size_t>(
        right.key.pos - source.ref_beg);
    if (left_offset >= source.ref_seq.size() ||
        right_offset >= source.ref_seq.size())
        return false;
    const char left_ref = source.ref_seq[left_offset];
    const char right_ref = source.ref_seq[right_offset];
    std::unordered_set<std::string> counted_reads;
    int consistent = 0;
    double wrong_parity_bound = 1.0;
    for (size_t ri = 0; ri < source.reads.size(); ++ri) {
        const ReadRecord& read = source.reads[ri];
        if (read.is_skipped || read.mapq < kMinBridgeMapq ||
            read.mapq == kUnknownMapq ||
            read.beg > left.key.pos || read.end < right.key.pos)
            continue;
        int left_quality = kMissingBaseQuality;
        int right_quality = kMissingBaseQuality;
        const int left_call = physical_snp_call(
            read, left.key.pos, left_ref, left.key.alt[0], &left_quality);
        const int right_call = physical_snp_call(
            read, right.key.pos, right_ref, right.key.alt[0], &right_quality);
        if (left_call != 0 && left_call != 2) continue;
        if (right_call != 0 && right_call != 2) continue;
        if (left_quality < kMinBridgeBaseq ||
            right_quality < kMinBridgeBaseq ||
            left_quality == kMissingBaseQuality ||
            right_quality == kMissingBaseQuality)
            continue;
        const int left_allele = left_call == 2 ? 1 : 0;
        const int right_allele = right_call == 2 ? 1 : 0;
        const int left_hap = left_allele == left.hap_to_cons_alle[1] ? 1 : 2;
        const int right_hap = right_allele == right.hap_to_cons_alle[1] ? 1 : 2;
        // A physical contradiction vetoes the bridge even if that molecule
        // was unassigned by the BAM sub-solve. Only source-assigned reads may
        // supply independent positive support.
        if (left_hap != right_hap) return false;
        if (ri >= source.phase_sets.size() ||
            source.phase_sets[ri] != source_ps ||
            ri >= source.haps.size() ||
            (source.haps[ri] != 1 && source.haps[ri] != 2))
            continue;
        if (source.haps[ri] != left_hap) return false;
        if (!counted_reads.insert(read.qname).second) continue;
        ++consistent;
        const auto error_probability = [](int quality) {
            return std::pow(10.0, -static_cast<double>(quality) / 10.0);
        };
        wrong_parity_bound *=
            error_probability(left_quality) +
            error_probability(right_quality) +
            error_probability(read.mapq);
    }
    return consistent >= kMinIndependentMolecules &&
           wrong_parity_bound <= kMaxWrongParityProbability;
}

// A repeat deletion may be verified by MSA but absent from the sparse
// profiles of reads carrying the reference haplotype. Validate that missing
// SNP-to-deletion path directly against the original alignment. Only exact
// reference bases qualify: a shifted repeat deletion cannot impersonate this
// candidate's ALT, and a physical contradiction vetoes the whole cut.
static bool high_quality_snp_repeat_ref_cut_support(
        const PhasingChunk& source, size_t left_i, size_t right_i) {
    constexpr int kMinBridgeMapq = 30;
    constexpr int kMinBridgeBaseq = 10;
    constexpr int kHighQualityBaseq = 30;
    constexpr int kMinHighQualityMolecules = 2;
    constexpr int kMissingBaseQuality = 255;
    constexpr int kUnknownMapq = 255;
    constexpr int kMinIndependentMolecules = 5;
    constexpr double kMaxRandomParityP = 0.05;
    const CandidateVariant& left = source.candidates[left_i];
    const CandidateVariant& right = source.candidates[right_i];
    if (left.key.type != VariantType::Snp ||
        left.counts.category != VariantCategory::CleanHetSnp ||
        left.key.alt.size() != 1 ||
        right.key.type != VariantType::Deletion ||
        !right.msa_verified || !right.is_homopolymer_indel ||
        right.counts.category != VariantCategory::NoisyCandHet ||
        left.key.pos < source.ref_beg ||
        left.key.pos - source.ref_beg >=
            static_cast<hts_pos_t>(source.ref_seq.size()) ||
        right.key.pos < source.ref_beg ||
        right.key.pos + right.key.ref_len - source.ref_beg >
            static_cast<hts_pos_t>(source.ref_seq.size()))
        return false;
    const char left_ref = source.ref_seq[
        static_cast<size_t>(left.key.pos - source.ref_beg)];
    std::unordered_set<std::string> counted_reads;
    std::unordered_set<std::string> high_quality_reads;
    for (const ReadRecord& read : source.reads) {
        if (read.is_skipped || !read.alignment ||
            read.mapq < kMinBridgeMapq || read.mapq == kUnknownMapq ||
            read.beg > left.key.pos || read.end < right.key.pos)
            continue;
        int left_quality = kMissingBaseQuality;
        const int left_call = physical_snp_call(
            read, left.key.pos, left_ref, left.key.alt[0], &left_quality);
        if ((left_call != 0 && left_call != 2) ||
            left_quality < kMinBridgeBaseq ||
            left_quality == kMissingBaseQuality)
            continue;
        int indel_qi = -1;
        const int right_call = bam_exact_indel_allele(
            read.alignment.get(), right, kMinBridgeBaseq, &indel_qi);
        if (right_call != 0) continue;
        bool exact_reference = true;
        for (hts_pos_t pos = right.key.pos;
             pos < right.key.pos + right.key.ref_len; ++pos) {
            exact_reference &= physical_snp_call(
                read, pos, source.ref_seq[
                    static_cast<size_t>(pos - source.ref_beg)], 'N') == 0;
        }
        if (!exact_reference) continue;
        const int left_allele = left_call == 2 ? 1 : 0;
        const int left_hap =
            left_allele == left.hap_to_cons_alle[1] ? 1 : 2;
        const int right_hap =
            right.hap_to_cons_alle[1] == 0 ? 1 : 2;
        if (left_hap != right_hap) return false;
        counted_reads.insert(read.qname);
        if (left_quality >= kHighQualityBaseq &&
            bam_get_qual(read.alignment.get())[indel_qi] >=
                kHighQualityBaseq)
            high_quality_reads.insert(read.qname);
    }
    return counted_reads.size() >= kMinIndependentMolecules &&
        high_quality_reads.size() >= kMinHighQualityMolecules &&
        std::ldexp(1.0, -static_cast<int>(counted_reads.size())) <=
            kMaxRandomParityP;
}

static SourcePathEvidence source_phase_set_path_evidence(
        const PhasingChunk& source, hts_pos_t source_ps) {
    std::vector<size_t> sites;
    for (size_t ci = 0; ci < source.candidates.size(); ++ci) {
        const CandidateVariant& cand = source.candidates[ci];
        if (cand.phase_set == source_ps && is_phase_set_anchor(cand))
            sites.push_back(ci);
    }
    SourcePathEvidence evidence;
    evidence.site_count = sites.size();
    if (sites.size() < 2) return evidence;
    std::sort(sites.begin(), sites.end(),
              [&source](size_t a, size_t b) {
                  return source.candidates[a].key.sort_pos() <
                         source.candidates[b].key.sort_pos();
              });
    std::vector<int> local_index(source.candidates.size(), -1);
    for (size_t i = 0; i < sites.size(); ++i)
        local_index[sites[i]] = static_cast<int>(i);
    std::vector<std::array<int, 3>> delta(sites.size() + 1);
    std::vector<std::array<int, 3>> other_phase_delta(sites.size() + 1);
    constexpr int kMinPathMapq = 30;
    constexpr int kUnknownMapq = 255;
    for (size_t ri = 0; ri < source.reads.size(); ++ri) {
        if (ri >= source.phase_sets.size() ||
            ri >= source.read_var_profile.size())
            continue;
        const bool assigned_to_source = source.phase_sets[ri] == source_ps;
        const ReadRecord& read = source.reads[ri];
        if (!assigned_to_source &&
            (source.phase_sets[ri] <= 0 || read.is_skipped ||
             read.mapq < kMinPathMapq || read.mapq == kUnknownMapq))
            continue;
        auto& path_delta = assigned_to_source ? delta : other_phase_delta;
        const ReadVariantProfile& profile = source.read_var_profile[ri];
        if (profile.start_var_idx < 0) continue;
        size_t previous = sites.size();
        int previous_hap = 0;
        for (size_t offset = 0; offset < profile.alleles.size(); ++offset) {
            const size_t ci = static_cast<size_t>(profile.start_var_idx) + offset;
            if (ci >= local_index.size()) break;
            const int local = local_index[ci];
            const int allele = profile.alleles[offset];
            if (local < 0 || allele < 0) continue;
            const CandidateVariant& cand = source.candidates[ci];
            const int hap = allele == cand.hap_to_cons_alle[1] ? 1
                          : allele == cand.hap_to_cons_alle[2] ? 2 : 0;
            const size_t current = static_cast<size_t>(local);
            if (previous != sites.size() && previous < current) {
                const size_t vote = hap != 0 && hap == previous_hap
                                        ? static_cast<size_t>(hap - 1) : 2;
                ++path_delta[previous][vote];
                --path_delta[current][vote];
            }
            previous = current;
            previous_hap = hap;
        }
    }
    constexpr double kSourceOneHapBridgeMaxP = 0.05;
    std::array<int, 3> crossing{};
    std::array<int, 3> other_phase_crossing{};
    for (size_t i = 0; i + 1 < sites.size(); ++i) {
        for (size_t vote = 0; vote < crossing.size(); ++vote) {
            crossing[vote] += delta[i][vote];
            other_phase_crossing[vote] += other_phase_delta[i][vote];
        }
        const hts_pos_t left = source.candidates[sites[i]].key.sort_pos();
        const hts_pos_t right = source.candidates[sites[i + 1]].key.sort_pos();
        // A conflict-free one-haplotype bridge is informative when its exact
        // one-sided binomial probability under random polarity is small. The
        // second haplotype need not physically span every source-site pair.
        const int consistent = crossing[0] + crossing[1];
        const bool both_haps = crossing[0] > 0 && crossing[1] > 0;
        const bool single_hap_confident = !both_haps && crossing[2] == 0 &&
            consistent > 0 &&
            std::ldexp(1.0, -consistent) <= kSourceOneHapBridgeMaxP;
        bool supported = consistent > crossing[2] &&
            (both_haps || single_hap_confident);
        if (!supported && crossing[2] == 0)
            supported = high_quality_snp_repeat_ref_cut_support(
                source, sites[i], sites[i + 1]);
        // Complementary MSA indels can start a new candidate PS while reads
        // calling both indel rows retain the preceding read PS. Those reads
        // repair the resulting false path cut only when high-MAPQ observations
        // support both haplotypes without a contradictory pair. Clean-SNP
        // source cuts keep the original read-label rule.
        const int extra_hap1 = crossing[0] + other_phase_crossing[0];
        const int extra_hap2 = crossing[1] + other_phase_crossing[1];
        const int extra_conflicts = crossing[2] + other_phase_crossing[2];
        constexpr int kMinExtraPerHaplotype = 2;
        constexpr double kExtraPathMaxP = 1e-6;
        const CandidateVariant& cut_site = source.candidates[sites[i]];
        const auto complementary_indel = [&](size_t other_i) {
            const CandidateVariant& other = source.candidates[sites[other_i]];
            return cut_site.key.type != VariantType::Snp &&
                   other.key.type != VariantType::Snp &&
                   cut_site.msa_verified && other.msa_verified &&
                   cut_site.key.sort_pos() == other.key.sort_pos() &&
                   cut_site.hap_to_cons_alle[1] != other.hap_to_cons_alle[1] &&
                   cut_site.hap_to_cons_alle[2] != other.hap_to_cons_alle[2];
        };
        const bool paired_indel_locus =
            (i > 0 && complementary_indel(i - 1)) ||
            (i + 1 < sites.size() && complementary_indel(i + 1));
        const bool extra_supported =
            paired_indel_locus &&
            extra_hap1 >= kMinExtraPerHaplotype &&
            extra_hap2 >= kMinExtraPerHaplotype &&
            extra_conflicts == 0 &&
            std::ldexp(1.0, -(extra_hap1 + extra_hap2)) <= kExtraPathMaxP;
        if (left < right && !supported && !extra_supported) {
            evidence.weak_cuts.push_back(left);
            if (high_quality_snp_cut_support(source, sites[i], sites[i + 1],
                                             source_ps))
                evidence.quality_supported_cuts.push_back(left);
        }
    }
    return evidence;
}

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

// The BAM solver can split two complete source paths at an MSA insertion even
// when the original alignments give a decisive SNP-to-insertion relation.
// Require a consistent exact-SNP gauge for each graph phase set before
// offering that physical edge to the normal stitcher.
static std::optional<bool> physical_graph_snp_insertion_bridge(
        const GraphChunkBuildResult& graph_chunk, const RecoverySeam& seam,
        const PhasingChunk& source) {
    constexpr int kMinBridgeMapq = 30;
    constexpr int kMinBridgeBaseq = 30;
    constexpr int kUnknownMapq = 255;
    constexpr int kMissingBaseq = 255;
    constexpr int kMinSupportPerHap = 2;
    constexpr double kMaxRandomParityP = 0.01;
    struct Anchor {
        size_t source_index = 0;
        hts_pos_t pos = 0;
        int graph_hap1 = -1;
    };
    struct Flank {
        std::vector<Anchor> anchors;
    };
    std::map<CandKey, size_t> source_by_key;
    for (size_t ci = 0; ci < source.candidates.size(); ++ci)
        source_by_key.emplace(cand_key_of(source.candidates[ci]), ci);
    std::array<Flank, 2> flanks;
    for (size_t gi = 0; gi < graph_chunk.chunk.candidates.size(); ++gi) {
        const CandidateVariant& graph_site = graph_chunk.chunk.candidates[gi];
        const size_t side = graph_site.phase_set == seam.left_phase_set ? 0 :
            graph_site.phase_set == seam.right_phase_set ? 1 : 2;
        if (side == 2 || graph_site.key.type != VariantType::Snp ||
            graph_site.counts.category != VariantCategory::CleanHetSnp ||
            graph_site.hap_to_cons_alle[1] < 0 ||
            graph_site.hap_to_cons_alle[1] > 1 ||
            graph_site.hap_to_cons_alle[2] !=
                1 - graph_site.hap_to_cons_alle[1] ||
            gi >= graph_chunk.site_meta.size())
            continue;
        const std::string* alt = selected_graph_candidate_alt(graph_chunk, gi);
        const GraphSiteMeta& meta = graph_chunk.site_meta[gi];
        if (alt == nullptr) continue;
        const VariantKey key = vcf_to_variant_key(
            source.region.tid, meta.pos, meta.ref, *alt);
        if (key.type != VariantType::Snp || key.ref_len != 1 ||
            key.alt.size() != 1)
            continue;
        const hts_pos_t pos = key.sort_pos();
        Flank& flank = flanks[side];
        const auto hit = source_by_key.find(CandKey{
            pos, static_cast<int>(key.type), key.ref_len, key.alt});
        if (hit == source_by_key.end()) continue;
        const CandidateVariant& bam_site = source.candidates[hit->second];
        if (bam_site.phase_set <= 0 ||
            (side == 1 &&
             (bam_site.hap_to_cons_alle[1] < 0 ||
              bam_site.hap_to_cons_alle[1] > 1 ||
              bam_site.hap_to_cons_alle[2] !=
                  1 - bam_site.hap_to_cons_alle[1])))
            continue;
        flank.anchors.push_back(Anchor{
            hit->second, pos, graph_site.hap_to_cons_alle[1]});
    }
    for (Flank& flank : flanks) {
        if (flank.anchors.size() < 2) return std::nullopt;
        std::sort(flank.anchors.begin(), flank.anchors.end(),
                  [](const Anchor& a, const Anchor& b) {
                      return a.pos < b.pos;
                  });
    }
    const Anchor& left = flanks[0].anchors.back();
    if (left.pos > seam.beg || left.pos >= seam.end)
        return std::nullopt;
    const CandidateVariant& left_bam = source.candidates[left.source_index];
    const hts_pos_t left_ps = left_bam.phase_set;
    const SourcePathEvidence left_path =
        source_phase_set_path_evidence(source, left_ps);
    if (left_path.site_count < 2 || !left_path.weak_cuts.empty())
        return std::nullopt;

    // The first verified insertion after the preselected left SNP is the
    // candidate bridge. Its exact CIGAR calls below supply physical validation;
    // do not search later alleles for a favorable vote.
    std::optional<size_t> insertion_index;
    for (size_t ci = 0; ci < source.candidates.size(); ++ci) {
        const CandidateVariant& site = source.candidates[ci];
        if (site.key.sort_pos() <= left.pos ||
            site.key.sort_pos() >= seam.end ||
            site.key.type != VariantType::Insertion ||
            !site.msa_verified ||
            site.phase_set <= 0 || site.phase_set == left_ps ||
            site.hap_to_cons_alle[1] < 0 ||
            site.hap_to_cons_alle[1] > 1 ||
            site.hap_to_cons_alle[2] != 1 - site.hap_to_cons_alle[1])
            continue;
        insertion_index = ci;
        break;
    }
    if (!insertion_index) return std::nullopt;
    const CandidateVariant& insertion = source.candidates[*insertion_index];
    const hts_pos_t right_ps = insertion.phase_set;
    const SourcePathEvidence right_path =
        source_phase_set_path_evidence(source, right_ps);
    if (right_path.site_count < 2 || !right_path.weak_cuts.empty())
        return std::nullopt;

    std::array<size_t, 2> matched{};
    std::array<int, 2> parity{-1, -1};
    for (size_t side = 0; side < 2; ++side) {
        const hts_pos_t required_ps = side == 0 ? left_ps : right_ps;
        for (const Anchor& anchor : flanks[side].anchors) {
            const CandidateVariant& bam_site =
                source.candidates[anchor.source_index];
            if (bam_site.phase_set != required_ps ||
                bam_site.counts.category != VariantCategory::CleanHetSnp)
                continue;
            const int this_parity =
                bam_site.hap_to_cons_alle[1] != anchor.graph_hap1;
            if (parity[side] >= 0 && parity[side] != this_parity)
                return std::nullopt;
            parity[side] = this_parity;
            ++matched[side];
        }
        if (matched[side] < 2)
            return std::nullopt;
    }
    // The BAM caller may phase a graph-biallelic SNP as alleles 1/2 after
    // local MSA. Its exact physical REF/ALT calls still orient the graph row;
    // the other shared clean SNPs certify the left block's source gauge.
    if (left_bam.hap_to_cons_alle[1] >= 0 &&
        left_bam.hap_to_cons_alle[1] <= 1 &&
        left_bam.hap_to_cons_alle[2] == 1 - left_bam.hap_to_cons_alle[1] &&
        left_bam.hap_to_cons_alle[1] !=
            (parity[0] == 0 ? left.graph_hap1 : 1 - left.graph_hap1))
        return std::nullopt;
    const int insertion_graph_hap1 = parity[1] == 0
        ? insertion.hap_to_cons_alle[1]
        : insertion.hap_to_cons_alle[2];
    if (left.pos < source.ref_beg ||
        left.pos - source.ref_beg >=
            static_cast<hts_pos_t>(source.ref_seq.size()) ||
        insertion.key.pos <= source.ref_beg ||
        insertion.key.pos - source.ref_beg >
            static_cast<hts_pos_t>(source.ref_seq.size()) ||
        left_bam.key.alt.size() != 1)
        return std::nullopt;
    const char left_ref = source.ref_seq[
        static_cast<size_t>(left.pos - source.ref_beg)];
    const char insertion_anchor_ref = source.ref_seq[
        static_cast<size_t>(insertion.key.pos - 1 - source.ref_beg)];
    std::array<int, 2> votes{};
    std::array<int, 2> support_by_left_hap{};
    std::unordered_set<std::string> counted;
    for (const ReadRecord& read : source.reads) {
        if (read.is_skipped || !read.alignment ||
            read.mapq < kMinBridgeMapq || read.mapq == kUnknownMapq ||
            read.beg > left.pos || read.end < insertion.key.pos ||
            counted.count(read.qname) != 0)
            continue;
        int left_quality = kMissingBaseq;
        int anchor_quality = kMissingBaseq;
        const int left_call = physical_snp_call(
            read, left.pos, left_ref, left_bam.key.alt[0], &left_quality);
        const int anchor_call = physical_snp_call(
            read, insertion.key.pos - 1, insertion_anchor_ref, 'N',
            &anchor_quality);
        if ((left_call != 0 && left_call != 2) || anchor_call != 0 ||
            left_quality < kMinBridgeBaseq ||
            anchor_quality < kMinBridgeBaseq ||
            left_quality == kMissingBaseq ||
            anchor_quality == kMissingBaseq)
            continue;
        int query_index = -1;
        const int insertion_call = bam_exact_indel_allele(
            read.alignment.get(), insertion, kMinBridgeBaseq,
            &query_index);
        if (insertion_call != 0 && insertion_call != 1) continue;
        counted.insert(read.qname);
        const bool left_hap1 =
            (left_call == 2 ? 1 : 0) == left.graph_hap1;
        const bool right_hap1 = insertion_call == insertion_graph_hap1;
        ++votes[left_hap1 != right_hap1 ? 1 : 0];
        ++support_by_left_hap[left_hap1 ? 0 : 1];
    }
    const int total = votes[0] + votes[1];
    const int winner = std::max(votes[0], votes[1]);
    if (votes[0] == votes[1] ||
        support_by_left_hap[0] < kMinSupportPerHap ||
        support_by_left_hap[1] < kMinSupportPerHap ||
        binomial_upper_tail(total, winner, 0.5) > kMaxRandomParityP)
        return std::nullopt;
    return votes[1] > votes[0];
}

// Targeted MSA can verify a SNP hidden from the ordinary clean category.
// Its physical bridge allele is still checked against the original BAM base
// and base quality before it can orient two blocks.
static bool verified_bam_bridge_snp(const CandidateVariant& site) {
    return is_phase_set_anchor(site) && site.key.type == VariantType::Snp &&
        (site.counts.category == VariantCategory::CleanHetSnp ||
         (site.counts.category == VariantCategory::NoisyCandHet &&
          site.msa_verified));
}

// A downstream MSA deletion may precede the first shared graph SNP. Its
// existing source allele and exact CIGAR must agree later; deletion REF also
// needs matching reference bases. Keep the established upstream anchor rule.
static bool verified_bam_bridge_site(const CandidateVariant& site,
                                     bool allow_msa_deletion = false) {
    return verified_bam_bridge_snp(site) ||
        (is_phase_set_anchor(site) &&
         site.key.type != VariantType::Snp &&
         (site.counts.category == VariantCategory::CleanHetIndel ||
          (allow_msa_deletion && site.key.type == VariantType::Deletion &&
           site.counts.category == VariantCategory::NoisyCandHet &&
           site.msa_verified)));
}

// A CIGAR REF call for a deletion only proves that its reference span is
// aligned. Check the bases too, so a different nearby allele cannot serve as
// reference evidence for the graph/BAM deletion.
static bool physical_deletion_reference_matches(
        const ReadRecord& read, const CandidateVariant& deletion,
        const PhasingChunk& source) {
    for (hts_pos_t base = deletion.key.pos;
         base < deletion.key.pos + deletion.key.ref_len; ++base) {
        if (base < source.ref_beg ||
            base - source.ref_beg >=
                static_cast<hts_pos_t>(source.ref_seq.size()) ||
            physical_snp_call(read, base, source.ref_seq[
                static_cast<size_t>(base - source.ref_beg)], 'N') != 0)
            return false;
    }
    return true;
}

constexpr hts_pos_t kSingletonValidationFlank = 5000;

struct ValidatedPhysicalBridge {
    bool flip = false;
    hts_pos_t pre_attach_source_ps = 0;
    bool multiallelic_cohort_corroborated = false;
};

// A deletion-to-SNP link needs both haplotypes and a significant allele
// association. A second deletion length at the same repeat locus can then
// check the first row's polarity without collapsing their representations.
static std::optional<bool> significant_deletion_cohort_flip(
        const PhasingChunk& source, size_t deletion_i, size_t right_i,
        int min_mapq) {
    const CandidateVariant& deletion = source.candidates[deletion_i];
    const CandidateVariant& right = source.candidates[right_i];
    std::array<int, 2> same_by_hap{};
    std::array<int, 2> cross_by_hap{};
    for (size_t ri = 0; ri < source.reads.size(); ++ri) {
        const ReadRecord& read = source.reads[ri];
        if (read.is_skipped || read.mapq < min_mapq || read.mapq == 255 ||
            ri >= source.read_var_profile.size())
            continue;
        const ReadVariantProfile& profile = source.read_var_profile[ri];
        if (profile.start_var_idx < 0 ||
            static_cast<int>(deletion_i) < profile.start_var_idx ||
            static_cast<int>(right_i) > profile.end_var_idx)
            continue;
        const size_t deletion_offset =
            deletion_i - static_cast<size_t>(profile.start_var_idx);
        const size_t right_offset =
            right_i - static_cast<size_t>(profile.start_var_idx);
        if (deletion_offset >= profile.alleles.size() ||
            right_offset >= profile.alleles.size())
            continue;
        const int deletion_allele = profile.alleles[deletion_offset];
        const int right_allele = profile.alleles[right_offset];
        if ((deletion_allele != 0 && deletion_allele != 1) ||
            (right_allele != 0 && right_allele != 1))
            continue;
        const int deletion_hap =
            deletion_allele == deletion.hap_to_cons_alle[1] ? 1 : 2;
        const int right_hap =
            right_allele == right.hap_to_cons_alle[1] ? 1 : 2;
        ++(deletion_hap == right_hap
            ? same_by_hap[static_cast<size_t>(deletion_hap - 1)]
            : cross_by_hap[static_cast<size_t>(deletion_hap - 1)]);
    }
    const int same = same_by_hap[0] + same_by_hap[1];
    const int cross = cross_by_hap[0] + cross_by_hap[1];
    const bool flip = cross > same;
    const auto& supporting = flip ? cross_by_hap : same_by_hap;
    constexpr double kMaxIndelCohortP = 0.01;
    if (same == cross || supporting[0] == 0 || supporting[1] == 0 ||
        binomial_upper_tail(same + cross, std::max(same, cross), 0.5) >
            kMaxIndelCohortP)
        return std::nullopt;
    return flip;
}

// A physical bridge can orient two blocks only after the full BAM solve confirms
// that each graph block keeps one allele gauge away from the seam. At 62.7 Mb
// a Q40 boundary read is locally correct, but its right graph block switches
// internally; a local-only check would join hundreds of misplaced reads.
static std::optional<ValidatedPhysicalBridge> validated_physical_bridge(
        const GraphChunkBuildResult& graph_chunk,
        const RecoverySeam& seam,
        const PhasingChunk& source,
        const Options& opts, bool allow_msa_deletion) {
    const PhasingChunk& graph = graph_chunk.chunk;
    std::map<CandKey, size_t> source_by_key;
    for (size_t ci = 0; ci < source.candidates.size(); ++ci)
        source_by_key.emplace(cand_key_of(source.candidates[ci]), ci);

    struct BlockEvidence {
        size_t graph_sites = 0;
        size_t matched_sites = 0;
        hts_pos_t first_graph_site = std::numeric_limits<hts_pos_t>::max();
        hts_pos_t last_graph_site = 0;
        hts_pos_t first_matched_site = std::numeric_limits<hts_pos_t>::max();
        hts_pos_t last_matched_site = 0;
        hts_pos_t source_ps = 0;
        int parity = -1;
        size_t boundary_graph = 0;
        size_t boundary_source = 0;
        hts_pos_t boundary_pos = 0;
        bool has_boundary = false;
        bool inconsistent = false;
    };
    std::array<BlockEvidence, 2> blocks;
    const std::array<hts_pos_t, 2> graph_ps{
        seam.left_phase_set, seam.right_phase_set};
    for (size_t gi = 0; gi < graph.candidates.size(); ++gi) {
        const CandidateVariant& candidate = graph.candidates[gi];
        const size_t side = candidate.phase_set == graph_ps[0] ? 0 :
                            candidate.phase_set == graph_ps[1] ? 1 : 2;
        if (side == 2 || !is_phase_set_anchor(candidate)) continue;
        const std::string* alt = selected_graph_candidate_alt(graph_chunk, gi);
        if (alt == nullptr) continue;
        const GraphSiteMeta& meta = graph_chunk.site_meta[gi];
        const VariantKey key = vcf_to_variant_key(
            source.region.tid, meta.pos, meta.ref, *alt);
        const CandKey normalized{key.sort_pos(), static_cast<int>(key.type),
                                 key.ref_len, key.alt};
        const auto matched = source_by_key.find(normalized);
        const hts_pos_t pos = key.sort_pos();
        BlockEvidence& block = blocks[side];
        const bool clean_snp = candidate.key.type == VariantType::Snp &&
            candidate.counts.category == VariantCategory::CleanHetSnp;
        const bool clean_indel = candidate.key.type != VariantType::Snp &&
            candidate.counts.category == VariantCategory::CleanHetIndel;
        if (clean_snp || clean_indel) {
            ++block.graph_sites;
            block.first_graph_site = std::min(block.first_graph_site, pos);
            block.last_graph_site = std::max(block.last_graph_site, pos);
        }
        if (matched == source_by_key.end()) continue;
        const CandidateVariant& bam_site = source.candidates[matched->second];
        if (!is_phase_set_anchor(bam_site) ||
            (clean_snp && (bam_site.key.type != VariantType::Snp ||
                           bam_site.counts.category != VariantCategory::CleanHetSnp)) ||
            (clean_indel && (bam_site.key.type != candidate.key.type ||
                             bam_site.counts.category != VariantCategory::CleanHetIndel)) ||
            (!clean_snp && !clean_indel))
            continue;
        const int parity = candidate.hap_to_cons_alle[1] ==
                           bam_site.hap_to_cons_alle[1] ? 0 : 1;
        // A graph block can cross several independently numbered BAM blocks.
        // Their numeric PS labels do not imply a graph polarity conflict;
        // compare the allele gauges, and use the nearest shared site to pick
        // the BAM source that actually reaches this seam.
        if (block.parity >= 0 && block.parity != parity)
            block.inconsistent = true;
        const bool closer = !block.has_boundary ||
            (side == 0 ? pos > block.boundary_pos : pos < block.boundary_pos);
        if (closer) {
            block.boundary_graph = gi;
            block.boundary_source = matched->second;
            block.boundary_pos = pos;
            block.source_ps = bam_site.phase_set;
            block.parity = parity;
            block.has_boundary = true;
        }
        ++block.matched_sites;
        block.first_matched_site = std::min(block.first_matched_site, pos);
        block.last_matched_site = std::max(block.last_matched_site, pos);
    }
    constexpr double kMinMatchedGraphSiteFraction = 0.5;
    for (const BlockEvidence& block : blocks) {
        // Exact source sites at both ends can validate a block even when its
        // repeat-rich middle has graph sites that the BAM caller demoted.
        // The full-block parity check above still vetoes conflicting matches.
        const bool spans_graph_block = block.matched_sites >= 2 &&
            block.first_matched_site < block.last_matched_site &&
            block.first_matched_site <=
                block.first_graph_site + kSingletonValidationFlank &&
            block.last_matched_site + kSingletonValidationFlank >=
                block.last_graph_site;
        if (block.inconsistent || !block.has_boundary ||
            block.source_ps <= 0 || block.matched_sites == 0 ||
            (static_cast<double>(block.matched_sites) <
                 kMinMatchedGraphSiteFraction * block.graph_sites &&
             !spans_graph_block))
            return std::nullopt;
    }
    // The first BAM site of the right source block can lie inside the graph
    // gap. A molecule need not reach the first catalog site of that block if
    // the complete BAM path already connects those two source sites.
    for (size_t side = 0; side < 2; ++side) {
        BlockEvidence& block = blocks[side];
        const size_t matched_source = block.boundary_source;
        const hts_pos_t matched_pos = block.boundary_pos;
        bool found = false;
        for (size_t ci = 0; ci < source.candidates.size(); ++ci) {
            const CandidateVariant& site = source.candidates[ci];
            if (site.phase_set != block.source_ps ||
                !verified_bam_bridge_site(site, allow_msa_deletion && side == 1))
                continue;
            const hts_pos_t pos = site.key.sort_pos();
            if ((side == 0 && pos >= seam.end) ||
                (side == 1 && pos <= seam.beg))
                continue;
            if (!found || (side == 0 ? pos > block.boundary_pos :
                                      pos < block.boundary_pos)) {
                block.boundary_source = ci;
                block.boundary_pos = pos;
                found = true;
            }
        }
        if (!found) return std::nullopt;
        const SourcePathEvidence path = source_phase_set_path_evidence(
            source, block.source_ps);
        const hts_pos_t graph_pos = recovery_position(
            graph_chunk, block.boundary_graph);
        const hts_pos_t lo = std::min(block.boundary_pos, graph_pos);
        const hts_pos_t hi = std::max(block.boundary_pos, graph_pos);
        if (std::any_of(path.weak_cuts.begin(), path.weak_cuts.end(),
                        [lo, hi](hts_pos_t cut) {
                            return lo <= cut && cut < hi;
                        })) {
            // A nearer private BAM site must not inherit the graph anchor's
            // gauge across an unsupported source cut. The exact shared site
            // remains a valid boundary and may have its own direct link.
            block.boundary_source = matched_source;
            block.boundary_pos = matched_pos;
        }
    }
    if (blocks[0].boundary_pos >= blocks[1].boundary_pos)
        return std::nullopt;
    // A single BAM phase set can still contain a weak internal cut. The
    // selected bridge must not silently turn that cut into a strong link.
    if (blocks[0].source_ps == blocks[1].source_ps) {
        const SourcePathEvidence path = source_phase_set_path_evidence(
            source, blocks[0].source_ps);
        if (std::any_of(path.weak_cuts.begin(), path.weak_cuts.end(),
                        [&](hts_pos_t cut) {
                            return blocks[0].boundary_pos <= cut &&
                                   cut < blocks[1].boundary_pos;
                        }))
            return std::nullopt;
    }

    constexpr int kMinBridgeMapq = 30;
    int paired_reads = 0;
    std::array<int, 2> paired_haplotypes{};
    std::optional<bool> flip;
    for (size_t ri = 0; ri < source.reads.size(); ++ri) {
        const ReadRecord& read = source.reads[ri];
        if (read.is_skipped || !read.alignment ||
            read.mapq < kMinBridgeMapq || read.mapq == 255 ||
            read.beg > blocks[0].boundary_pos ||
            read.end < blocks[1].boundary_pos ||
            ri >= source.read_var_profile.size())
            continue;
        const ReadVariantProfile& profile = source.read_var_profile[ri];
        std::array<int, 2> hap{};
        for (size_t side = 0; side < 2; ++side) {
            const BlockEvidence& block = blocks[side];
            const int ci = static_cast<int>(block.boundary_source);
            const int offset = ci - profile.start_var_idx;
            const CandidateVariant& site = source.candidates[block.boundary_source];
            int allele = profile.start_var_idx >= 0 &&
                         ci <= profile.end_var_idx && offset >= 0 &&
                         static_cast<size_t>(offset) < profile.alleles.size()
                ? profile.alleles[static_cast<size_t>(offset)] : -1;
            if (site.key.type != VariantType::Snp) {
                int query_index = -1;
                const int physical = bam_exact_indel_allele(
                    read.alignment.get(), site, opts.min_bq, &query_index);
                // An indel cannot gain a bridge allele from a missing sparse
                // profile call: require independent agreement with the BAM
                // solve's existing observation.
                if (physical < 0 || allele != physical ||
                    (physical == 0 && site.key.type == VariantType::Deletion &&
                     !physical_deletion_reference_matches(read, site, source)))
                    continue;
            } else {
                if (site.key.ref_len != 1 || site.key.alt.size() != 1 ||
                    site.key.pos < source.ref_beg ||
                    site.key.pos - source.ref_beg >=
                        static_cast<hts_pos_t>(source.ref_seq.size()))
                    continue;
                const char ref_base = source.ref_seq[
                    static_cast<size_t>(site.key.pos - source.ref_beg)];
                // The BAM profile can mask a clean SNP inside a noisy read
                // interval. For this one physical bridge, use its original
                // aligned base when the profile omitted the allele.
                if (allele < 0 && site.key.alt.size() == 1 &&
                    site.key.pos >= source.ref_beg &&
                    site.key.pos - source.ref_beg <
                        static_cast<hts_pos_t>(source.ref_seq.size())) {
                    const char ref = source.ref_seq[
                        static_cast<size_t>(site.key.pos - source.ref_beg)];
                    const int physical = physical_snp_call(
                        read, site.key.pos, ref, site.key.alt[0]);
                    allele = physical == 0 ? 0 : physical == 2 ? 1 : -1;
                }
                const uint8_t quality = bam_snp_observation_quality(
                    read.alignment.get(), site.key.pos, ref_base, site.key.alt[0], allele);
                if (quality == 0 || quality < opts.min_bq) continue;
            }
            if (allele != 0 && allele != 1) continue;
            const int source_hap = allele == site.hap_to_cons_alle[1] ? 1 :
                                   allele == site.hap_to_cons_alle[2] ? 2 : 0;
            if (source_hap == 0) continue;
            hap[side] = block.parity == 0 ? source_hap : 3 - source_hap;
        }
        if (hap[0] == 0 || hap[1] == 0) continue;
        ++paired_reads;
        ++paired_haplotypes[static_cast<size_t>(hap[0] - 1)];
        const bool read_flip = hap[0] != hap[1];
        if (flip && *flip != read_flip) return std::nullopt;
        flip = read_flip;
    }
    // Ordinary multi-read SNP links keep their statistical stitch path. A
    // private downstream MSA deletion can instead carry the boundary when
    // its exact physical calls cover both haplotypes and all pairs agree.
    // Missing or contradictory pairs never acquire an allele from read HP.
    const CandidateVariant& right_boundary =
        source.candidates[blocks[1].boundary_source];
    const bool msa_deletion_boundary =
        right_boundary.key.type == VariantType::Deletion &&
        right_boundary.msa_verified &&
        right_boundary.counts.category == VariantCategory::NoisyCandHet;
    if ((!allow_msa_deletion && paired_reads == 1) ||
        (allow_msa_deletion && paired_reads > 1 && msa_deletion_boundary &&
         paired_haplotypes[0] > 0 && paired_haplotypes[1] > 0))
        return ValidatedPhysicalBridge{*flip, 0};
    if (paired_reads != 0) return std::nullopt;

    constexpr int kMinIndelBridgeBaseq = 10;
    std::optional<bool> indel_flip;
    bool multiallelic_cohort_corroborated = false;
    for (size_t di = 0; di < source.candidates.size(); ++di) {
        const CandidateVariant& deletion = source.candidates[di];
        const hts_pos_t pos = deletion.key.sort_pos();
        if (deletion.key.type != VariantType::Deletion ||
            deletion.counts.category != VariantCategory::NoisyCandHet ||
            !deletion.msa_verified || !is_phase_set_anchor(deletion) ||
            deletion.phase_set != blocks[0].source_ps ||
            pos <= blocks[0].boundary_pos || pos >= blocks[1].boundary_pos)
            continue;
        const std::optional<bool> cohort_flip =
            significant_deletion_cohort_flip(
                source, di, blocks[1].boundary_source, kMinBridgeMapq);
        if (!cohort_flip) continue;
        bool complementary_cohort = false;
        for (size_t other_i = 0; other_i < source.candidates.size(); ++other_i) {
            const CandidateVariant& other = source.candidates[other_i];
            if (other_i == di || other.key.type != VariantType::Deletion ||
                other.key.pos != deletion.key.pos ||
                other.key.ref_len == deletion.key.ref_len ||
                other.counts.category != VariantCategory::NoisyCandHet ||
                !other.msa_verified || !is_phase_set_anchor(other) ||
                other.phase_set != deletion.phase_set ||
                other.hap_to_cons_alle[1] == deletion.hap_to_cons_alle[1])
                continue;
            const std::optional<bool> other_flip =
                significant_deletion_cohort_flip(
                    source, other_i, blocks[1].boundary_source,
                    kMinBridgeMapq);
            if (other_flip && *other_flip == *cohort_flip) {
                complementary_cohort = true;
                break;
            }
        }

        int physical_pairs = 0;
        std::optional<bool> physical_flip;
        const CandidateVariant& left =
            source.candidates[blocks[0].boundary_source];
        for (size_t ri = 0; ri < source.reads.size(); ++ri) {
            const ReadRecord& read = source.reads[ri];
            if (read.is_skipped || !read.alignment ||
                read.mapq < kMinBridgeMapq || read.mapq == 255 ||
                read.beg > blocks[0].boundary_pos || read.end < pos ||
                ri >= source.read_var_profile.size())
                continue;
            const ReadVariantProfile& profile = source.read_var_profile[ri];
            const int offset = static_cast<int>(blocks[0].boundary_source) -
                               profile.start_var_idx;
            if (profile.start_var_idx < 0 || offset < 0 ||
                static_cast<size_t>(offset) >= profile.alleles.size() ||
                left.key.ref_len != 1 || left.key.alt.size() != 1 ||
                left.key.pos < source.ref_beg ||
                left.key.pos - source.ref_beg >=
                    static_cast<hts_pos_t>(source.ref_seq.size()))
                continue;
            const int left_allele = profile.alleles[static_cast<size_t>(offset)];
            if (left_allele != 0 && left_allele != 1) continue;
            const uint8_t quality = bam_snp_observation_quality(
                read.alignment.get(), left.key.pos,
                source.ref_seq[static_cast<size_t>(left.key.pos - source.ref_beg)],
                left.key.alt[0], left_allele);
            if (quality == 0 || quality < opts.min_bq) continue;
            int query_index = -1;
            const int deletion_allele = bam_exact_indel_allele(
                read.alignment.get(), deletion,
                std::max(opts.min_bq, kMinIndelBridgeBaseq), &query_index);
            if (deletion_allele < 0) continue;
            if (deletion_allele == 0 &&
                !physical_deletion_reference_matches(read, deletion, source))
                continue;
            const int profile_offset = static_cast<int>(di) -
                                       profile.start_var_idx;
            if (profile_offset >= 0 &&
                static_cast<size_t>(profile_offset) < profile.alleles.size()) {
                const int called = profile.alleles[
                    static_cast<size_t>(profile_offset)];
                if (called >= 0 && called != deletion_allele)
                    return std::nullopt;
            }
            if (ri >= source.phase_sets.size() ||
                source.phase_sets[ri] != blocks[0].source_ps ||
                ri >= source.haps.size() ||
                (source.haps[ri] != 1 && source.haps[ri] != 2))
                continue;
            const int left_source_hap =
                left_allele == left.hap_to_cons_alle[1] ? 1 : 2;
            const int deletion_source_hap =
                deletion_allele == deletion.hap_to_cons_alle[1] ? 1 : 2;
            if (source.haps[ri] != left_source_hap) continue;
            const int left_graph_hap = blocks[0].parity == 0
                ? left_source_hap : 3 - left_source_hap;
            const int right_source_hap = *cohort_flip
                ? 3 - deletion_source_hap : deletion_source_hap;
            const int right_graph_hap = blocks[1].parity == 0
                ? right_source_hap : 3 - right_source_hap;
            const bool bridge_flip = left_graph_hap != right_graph_hap;
            if (physical_flip && *physical_flip != bridge_flip)
                return std::nullopt;
            physical_flip = bridge_flip;
            ++physical_pairs;
        }
        if (physical_pairs != 1) continue;
        if (indel_flip && *indel_flip != *physical_flip)
            return std::nullopt;
        indel_flip = physical_flip;
        multiallelic_cohort_corroborated |= complementary_cohort;
    }
    return indel_flip
        ? std::optional<ValidatedPhysicalBridge>(ValidatedPhysicalBridge{
              *indel_flip, blocks[0].source_ps,
              multiallelic_cohort_corroborated})
        : std::nullopt;
}


static void apply_selected_msa_observations(
        PhasingChunk& source, const GraphChunkBuildResult& graph_chunk,
        const Options& opts, std::map<hts_pos_t, SourcePathEvidence>& original_paths) {
    if (source.pending_msa_observations.empty()) return;
    if (!opts.phase_matrix_dump_prefix.empty())
        dump_recovery_phase_state(source, opts, "msa-transfer-pending");
    CandidateIndex source_index;
    for (size_t ci = 0; ci < source.candidates.size(); ++ci)
        source_index.emplace(cand_key_of(source.candidates[ci]), ci);
    using Gauge = std::pair<hts_pos_t, bool>;
    std::map<hts_pos_t, std::optional<Gauge>> gauges;
    std::set<hts_pos_t> graph_owned_source_phase_sets;
    for (size_t gi = 0; gi < graph_chunk.chunk.candidates.size(); ++gi) {
        const auto& graph = graph_chunk.chunk.candidates[gi];
        if (!graph.graph_site || graph.bam_injected) continue;
        const std::string* alt = selected_graph_candidate_alt(graph_chunk, gi);
        if (alt == nullptr) continue;
        const auto& meta = graph_chunk.site_meta[gi];
        const VariantKey key = vcf_to_variant_key(source.region.tid, meta.pos, meta.ref, *alt);
        const auto found = source_index.find(CandKey{
            key.sort_pos(), static_cast<int>(key.type), key.ref_len, key.alt});
        if (found == source_index.end()) continue;
        const auto& bam = source.candidates[found->second];
        if (!is_phase_set_anchor(bam)) continue;
        graph_owned_source_phase_sets.insert(bam.phase_set);
        if (!is_phase_set_anchor(graph) ||
            graph.counts.category != VariantCategory::CleanHetSnp ||
            bam.counts.category != VariantCategory::CleanHetSnp) continue;
        const Gauge gauge{graph.phase_set,
            (graph.hap_to_cons_alle[1] == 1) != (bam.hap_to_cons_alle[1] == 1)};
        const auto [it, inserted] = gauges.try_emplace(bam.phase_set, gauge);
        if (!inserted && it->second != gauge) it->second.reset();
    }
    // A local diploid recall may extend one established graph gauge. A BAM
    // block touching different graph phase sets needs independent stitching
    // evidence; its inherited PS cannot certify their relative orientation.
    auto& pending = source.pending_msa_observations;
    std::set<hts_pos_t> recalled_phase_sets;
    for (const auto& observation : pending) {
        if (!observation.update_counts) continue;
        const auto found = source_index.find(CandKey{
            observation.key.sort_pos(), static_cast<int>(observation.key.type),
            observation.key.ref_len, observation.key.alt});
        if (found != source_index.end())
            recalled_phase_sets.insert(source.candidates[found->second].phase_set);
    }
    std::ofstream admission;
    if (!opts.phase_matrix_dump_prefix.empty()) {
        const std::string path = opts.phase_matrix_dump_prefix + ".msa-admission.tsv";
        admission.open(path);
        if (!admission) throw std::runtime_error("cannot write MSA admission trace: " + path);
        admission << "qname\tpos\ttype\tref_len\talt\tallele\tupdate_counts\tstatus\n";
    }
    pending.erase(std::remove_if(pending.begin(), pending.end(),
        [&](const DeferredMsaObservation& observation) {
            const auto found = source_index.find(CandKey{
                observation.key.sort_pos(), static_cast<int>(observation.key.type),
                observation.key.ref_len, observation.key.alt});
            // Physical corrections support the newly recalled local diploid
            // contrast. Without that fixed-consensus context, keep the source's
            // existing insertion projection and read-rescue behavior.
            std::string_view status = "eligible";
            if (found == source_index.end()) status = "source_site_missing";
            else {
                const hts_pos_t ps = source.candidates[found->second].phase_set;
                const auto gauge = gauges.find(ps);
                if (observation.update_counts) {
                    // Keep the verified allele independently of graph ownership.
                    // Extra calls must not erase a weak cut in the original
                    // connection certificate of an unresolved graph/BAM gauge.
                    const bool unresolved_graph_gauge = gauge == gauges.end()
                        ? graph_owned_source_phase_sets.count(ps) != 0 : !gauge->second;
                    if (unresolved_graph_gauge && original_paths.count(ps) == 0)
                        original_paths.emplace(ps, source_phase_set_path_evidence(source, ps));
                } else if (recalled_phase_sets.count(ps) == 0)
                    status = "no_fixed_consensus_context";
                else if (gauge == gauges.end()) {
                    status = "no_shared_clean_snp";
                } else if (!gauge->second) status = "inconsistent_graph_gauge";
            }
            if (admission.is_open())
                admission << source.reads[observation.read_id].qname << '\t'
                          << observation.key.pos << '\t' << static_cast<int>(observation.key.type)
                          << '\t' << observation.key.ref_len << '\t' << observation.key.alt
                          << '\t' << observation.allele << '\t' << observation.update_counts
                          << '\t' << status << '\n';
            return status != "eligible";
        }), pending.end());
    apply_pending_msa_observations(source);
    if (!opts.phase_matrix_dump_prefix.empty())
        dump_recovery_phase_state(source, opts, "msa-transfer-selected");
}

bool recover_phase_set_seams_in_place(GraphChunkBuildResult& graph_chunk,
                                      const Options& opts,
                                      WorkerContext& context,
                                      const char* contig_name,
                                      const std::vector<RecoverySeam>* completed_seams) {
    PhasingChunk& chunk = graph_chunk.chunk;
    if (chunk.candidates.empty() || chunk.reads.empty()) return false;

    const int solve_tid = contig_name != nullptr
                              ? sam_hdr_name2tid(context.primary_header(), contig_name)
                              : chunk.region.tid;
    if (solve_tid < 0) return false;

    std::vector<RecoverySeam> windows = collect_phase_set_seams(graph_chunk);
    std::vector<RecoveryPhysicalSnpBridge> retry_bridges;
    if (completed_seams != nullptr) {
        // A BAM block imported by the first pass may split an old graph seam.
        // Retry its new phase-set pair only when physical SNP calls already
        // support one orientation. Overlap alone does not identify an old pair.
        windows.erase(std::remove_if(windows.begin(), windows.end(),
            [&](const RecoverySeam& seam) {
                const bool already_solved = std::any_of(
                    completed_seams->begin(), completed_seams->end(),
                    [&seam](const RecoverySeam& completed) {
                        return seam.left_phase_set == completed.left_phase_set &&
                               seam.right_phase_set == completed.right_phase_set;
                    });
                if (already_solved) return true;
                const bool overlaps = std::any_of(
                    completed_seams->begin(), completed_seams->end(),
                    [&seam](const RecoverySeam& completed) {
                        return seam.beg < completed.end &&
                               completed.beg < seam.end;
                    });
                if (!overlaps) return false;
                // Adjacent anchors leave no unphased reference base for the
                // second BAM solve. Re-solving them can retag reads without
                // closing an actual gap.
                if (seam.end <= seam.beg + 1) return true;
                std::optional<bool> graph_parity;
                if (!has_direct_snp_parity_for_retry(
                        graph_chunk, seam, context, solve_tid, &graph_parity))
                    return true;
                // Admission proved a graph SNP orientation, not just coverage.
                // Carry that relation into the stitch so an indel edge cannot
                // silently reverse the blocks that admitted this solve.
                if (graph_parity)
                    retry_bridges.push_back(RecoveryPhysicalSnpBridge{
                        seam.left_phase_set, seam.right_phase_set, *graph_parity});
                return false;
            }), windows.end());
    }
    if (windows.empty()) return false;
    if (!opts.phase_matrix_dump_prefix.empty()) {
        for (const RecoverySeam& seam : windows) {
            std::fprintf(stderr,
                         "[recovery-seam] %" PRId64 "-%" PRId64
                         " left_ps=%" PRId64 " right_ps=%" PRId64 "\n",
                         static_cast<int64_t>(seam.beg),
                         static_cast<int64_t>(seam.end),
                         static_cast<int64_t>(seam.left_phase_set),
                         static_cast<int64_t>(seam.right_phase_set));
        }
    }

    const std::vector<hts_pos_t> parent_sites = parent_phased_positions(graph_chunk);
    // The ordinary BAM solve groups touching seams. A rejected MSA retry can
    // instead solve one seam over the complete adjacent graph phase sets.
    // These extents are from the original graph solve, before any BAM transfer.
    std::map<hts_pos_t, std::pair<hts_pos_t, hts_pos_t>> phase_set_extents;
    for (size_t ci = 0; ci < chunk.candidates.size(); ++ci) {
        const CandidateVariant& candidate = chunk.candidates[ci];
        if (!is_phase_set_anchor(candidate)) continue;
        const hts_pos_t pos = recovery_position(graph_chunk, ci);
        const auto [it, inserted] = phase_set_extents.try_emplace(
            candidate.phase_set, pos, pos);
        if (!inserted) {
            it->second.first = std::min(it->second.first, pos);
            it->second.second = std::max(it->second.second, pos);
        }
    }
    const std::vector<TargetedWindowGroup> initial_groups =
        build_targeted_groups(windows, solve_tid, parent_sites,
                              chunk.ref_beg, chunk.ref_end);
    if (initial_groups.empty()) return false;

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
    for (const TargetedWindowGroup& group : initial_groups)
        for (size_t wi = group.first_window; wi < group.past_last_window; ++wi)
            sub.retry_windows.emplace_back(windows[wi].beg, windows[wi].end);
    // Keep depth-based heterozygote repair local to the BAM sub-solve.
    // The same windows become the final left-to-right stitch targets.
    graph_chunk.recovery_windows = windows;
    std::vector<TargetedWindowGroup> groups;
    std::vector<PhasingChunk> discovered;
    std::vector<std::map<hts_pos_t, SourcePathEvidence>> source_paths;
    groups.reserve(initial_groups.size());
    discovered.reserve(initial_groups.size());
    const auto backfill = [&](PhasingChunk& target,
                              const TargetedWindowGroup& target_group,
                              const Options& target_opts) {
        for (size_t wi = target_group.first_window;
             wi < target_group.past_last_window; ++wi) {
            if (!group_owns_window(target_group, wi)) continue;
            backfill_msa_observations(target, target_opts,
                                      windows[wi].beg, windows[wi].end);
        }
    };
    const auto preserves_hets = [](const PhasingChunk& before,
                                   const PhasingChunk& after,
                                   bool preserve_clean_gauge = false) {
        std::map<CandKey, const CandidateVariant*> after_hets;
        for (const CandidateVariant& candidate : after.candidates)
            if (is_phase_set_anchor(candidate))
                after_hets.emplace(cand_key_of(candidate), &candidate);
        std::map<hts_pos_t, std::pair<hts_pos_t, bool>> clean_gauges;
        return std::all_of(
            before.candidates.begin(), before.candidates.end(),
            [&](const CandidateVariant& candidate) {
                if (!is_phase_set_anchor(candidate)) return true;
                const auto found = after_hets.find(cand_key_of(candidate));
                if (found == after_hets.end()) return false;
                if (!preserve_clean_gauge ||
                    candidate.counts.category != VariantCategory::CleanHetSnp)
                    return true;
                const CandidateVariant& next = *found->second;
                const bool flip = candidate.hap_to_cons_alle[1] == next.hap_to_cons_alle[2];
                if (!flip && candidate.hap_to_cons_alle[1] != next.hap_to_cons_alle[1])
                    return false;
                const auto [gauge, inserted] = clean_gauges.emplace(candidate.phase_set,
                    std::make_pair(next.phase_set, flip));
                return inserted || gauge->second == std::make_pair(next.phase_set, flip);
            });
    };
    const auto improves_source_paths = [](const PhasingChunk& before,
                                          const PhasingChunk& after,
                                          const RecoverySeam& seam) {
        std::map<CandKey, size_t> after_index;
        std::map<hts_pos_t, std::vector<size_t>> before_sites;
        for (size_t ci = 0; ci < after.candidates.size(); ++ci)
            if (is_phase_set_anchor(after.candidates[ci]))
                after_index.emplace(cand_key_of(after.candidates[ci]), ci);
        for (size_t ci = 0; ci < before.candidates.size(); ++ci)
            if (is_phase_set_anchor(before.candidates[ci]))
                before_sites[before.candidates[ci].phase_set].push_back(ci);
        std::map<hts_pos_t, SourcePathEvidence> after_paths;
        std::map<hts_pos_t, hts_pos_t> after_owners;
        const auto has_cut = [](const SourcePathEvidence& path,
                                hts_pos_t left, hts_pos_t right) {
            const auto cut = std::lower_bound(path.weak_cuts.begin(), path.weak_cuts.end(), left);
            return cut != path.weak_cuts.end() && *cut < right;
        };
        bool improved = false;
        for (auto& [ps, sites] : before_sites) {
            // Partial improvements cannot merge independent source gauges.
            // Cross-block connections still belong to the validated stitch.
            for (const size_t ci : sites) {
                const auto found = after_index.find(cand_key_of(before.candidates[ci]));
                if (found == after_index.end()) return false;
                const hts_pos_t next_ps = after.candidates[found->second].phase_set;
                const auto [owner, inserted] = after_owners.emplace(next_ps, ps);
                if (!inserted && owner->second != ps) return false;
            }
            const SourcePathEvidence old_path = source_phase_set_path_evidence(before, ps);
            std::sort(sites.begin(), sites.end(), [&before](size_t a, size_t b) {
                return before.candidates[a].key.sort_pos() < before.candidates[b].key.sort_pos();
            });
            for (size_t si = 1; si < sites.size(); ++si) {
                const auto& left = before.candidates[sites[si - 1]];
                const auto& right = before.candidates[sites[si]];
                const hts_pos_t beg = left.key.sort_pos(), end = right.key.sort_pos();
                if (beg >= end) continue;
                const auto ai = after_index.find(cand_key_of(left));
                const auto bi = after_index.find(cand_key_of(right));
                if (ai == after_index.end() || bi == after_index.end()) return false;
                const auto& next_left = after.candidates[ai->second];
                const auto& next_right = after.candidates[bi->second];
                bool connected = false;
                if (next_left.phase_set == next_right.phase_set) {
                    auto path = after_paths.find(next_left.phase_set);
                    if (path == after_paths.end())
                        path = after_paths.emplace(next_left.phase_set,
                            source_phase_set_path_evidence(after, next_left.phase_set)).first;
                    connected = !has_cut(path->second, beg, end);
                }
                const bool was_connected = !has_cut(old_path, beg, end);
                if (was_connected && !connected) return false;
                if (!was_connected && connected && beg >= seam.beg && end <= seam.end + 1)
                    improved = true;
            }
        }
        return improved;
    };
    for (const TargetedWindowGroup& group : initial_groups) {
        // All BAM sub-solves use chunk_id=-1. Scope their diagnostic output to
        // this graph chunk and window so a later seam cannot overwrite it.
        const std::string dump_prefix = opts.phase_matrix_dump_prefix.empty() ?
            std::string() :
            opts.phase_matrix_dump_prefix + ".recovery.chunk" +
                std::to_string(chunk.region.chunk_id) + ".window" +
                std::to_string(group.first_window);
        Options source_opts = sub;
        source_opts.recall_unplaced_msa_insertions = true;
        if (!dump_prefix.empty())
            source_opts.phase_matrix_dump_prefix = dump_prefix + ".initial";
        PhasingChunk source = process_chunk(group.region, source_opts, context);
        // Measure SNP dropout in the matrix used by the source solve before
        // exact-CIGAR backfill hides it. A newly admitted retry retains every
        // source heterozygote; filling calls after solving does not rephase them.
        bool raw_preserve = false;
        bool raw_ordinary = false;
        bool raw_internal_retry = false;
        std::optional<size_t> raw_sparse;
        const bool raw_needs_retry = source_seam_needs_unplaced_msa(
            source, group, windows, opts, raw_preserve, raw_ordinary,
            raw_sparse, true, std::nullopt, false, nullptr, true,
            &raw_internal_retry);
        // Scan every requested seam before backfill hides its original dropout.
        // A rejected earlier request must not suppress a later independent one.
        std::vector<size_t> focused_windows;
        std::vector<size_t> focused_repair_windows;
        for (size_t wi = group.first_window; wi < group.past_last_window; ++wi) {
            bool preserve = false, ordinary = false;
            std::optional<size_t> requested;
            if (source_seam_needs_unplaced_msa(source, group, windows, opts,
                    preserve, ordinary, requested, true, wi, true) && !ordinary) {
                focused_repair_windows.push_back(wi);
            }
        }
        backfill(source, group, source_opts);
        // Physical homopolymer calls can expose a source conflict, but are
        // insufficient alone to tag reads. Use them only to select a validated
        // MSA retry, then restore the source projection before any transfer.
        auto source_profiles = source.read_var_profile;
        auto source_index = std::move(source.read_var_cr);
        // A singleton graph flank cannot supply independent orientation context.
        // Diagnose its isolated MSA deletion dropout before requesting the
        // padded retry; ordinary block boundaries retain the paired-row rule.
        const auto has_singleton_flank = [&](size_t wi) {
            const RecoverySeam& seam = windows[wi];
            const auto left = phase_set_extents.find(seam.left_phase_set);
            const auto right = phase_set_extents.find(seam.right_phase_set);
            return (left != phase_set_extents.end() &&
                    left->second.first == left->second.second) ||
                   (right != phase_set_extents.end() &&
                    right->second.first == right->second.second);
        };
        for (size_t wi = group.first_window; wi < group.past_last_window; ++wi) {
            if (group_owns_window(group, wi))
                backfill_msa_retry_deletions(
                    source, source_opts, windows[wi].beg, windows[wi].end, false);
        }
        // Backfill can reveal an internal conflict in a nominally connected
        // source PS. Keep every focused request, including later group seams.
        for (size_t wi = group.first_window; wi < group.past_last_window; ++wi) {
            if (!group_owns_window(group, wi) ||
                !source_seam_has_msa_conflict(source, windows[wi], opts)) continue;
            const auto found = std::lower_bound(focused_repair_windows.begin(),
                                                focused_repair_windows.end(), wi);
            if (found == focused_repair_windows.end() || *found != wi)
                focused_repair_windows.insert(found, wi);
        }
        TargetedWindowGroup remainder = group;
        bool focused_accepted = false;
        bool preserve_source_rows = false;
        bool ordinary_retry = false;
        std::optional<size_t> sparse_window;
        bool conflicting_majority = false;
        bool filled_internal_retry = false;
        bool filled_needs_retry = source_seam_needs_unplaced_msa(
            source, group, windows, opts, preserve_source_rows,
            ordinary_retry, sparse_window, false, std::nullopt, false,
            &conflicting_majority, true, &filled_internal_retry);
        // Isolated homopolymer calls may corroborate an internal conflict,
        // but must not change the ordinary retry/dropout decision above.
        // Existing pair conflicts retain their established selection/context.
        std::vector<size_t> isolated_conflict_windows;
        for (size_t wi = group.first_window; wi < group.past_last_window; ++wi) {
            if (!group_owns_window(group, wi) || !has_singleton_flank(wi)) continue;
            if (source_seam_has_msa_conflict(source, windows[wi], opts)) continue;
            if (backfill_msa_retry_deletions(
                    source, source_opts, windows[wi].beg, windows[wi].end, true) == 0 ||
                !source_seam_has_msa_conflict(source, windows[wi], opts)) continue;
            isolated_conflict_windows.push_back(wi);
            const auto found = std::lower_bound(focused_repair_windows.begin(),
                                                focused_repair_windows.end(), wi);
            if (found == focused_repair_windows.end() || *found != wi)
                focused_repair_windows.insert(found, wi);
        }
        if (!isolated_conflict_windows.empty()) {
            sparse_window = isolated_conflict_windows.front();
            filled_needs_retry = filled_internal_retry = true;
            ordinary_retry = false;
        }
        source.read_var_profile = std::move(source_profiles);
        source.read_var_cr = std::move(source_index);
        if (!filled_needs_retry && raw_needs_retry) {
            preserve_source_rows = true;
            ordinary_retry = raw_ordinary;
            sparse_window = raw_sparse;
        }
        const bool needs_retry = filled_needs_retry || raw_needs_retry;
        const bool internal_retry = filled_needs_retry
            ? filled_internal_retry : raw_internal_retry;
        // Try the established selection first. If it fails certification,
        // later dropout and internal-conflict requests get independent attempts.
        // Internal conflicts authorize only a focused trial, never the broad
        // fallback controlled by needs_retry above.
        // One accepted source owns the seam; importing overlapping focused
        // context from multiple solves would duplicate source evidence.
        if (sparse_window && !ordinary_retry) focused_windows.push_back(*sparse_window);
        for (const size_t wi : focused_repair_windows)
            if (!sparse_window || wi != *sparse_window || ordinary_retry)
                focused_windows.push_back(wi);
        for (const size_t wi : focused_windows) {
            const bool fallback = !sparse_window || ordinary_retry || wi != *sparse_window;
            // Newly exposed internal conflicts need the same solve padding as
            // raw MSA dropout. A singleton flank otherwise ends at the seam,
            // making the strict containment check silently skip its retry.
            const bool focused_repair = fallback
                ? std::binary_search(focused_repair_windows.begin(),
                                     focused_repair_windows.end(), wi)
                : (raw_needs_retry && !raw_ordinary) ||
                    std::binary_search(isolated_conflict_windows.begin(),
                                       isolated_conflict_windows.end(), wi);
            const std::string focused_prefix = dump_prefix.empty() ? std::string() :
                dump_prefix + ".focused.seam" + std::to_string(wi);
            const RecoverySeam& seam = windows[wi];
            const auto left = phase_set_extents.find(seam.left_phase_set);
            const auto right = phase_set_extents.find(seam.right_phase_set);
            if (left != phase_set_extents.end() &&
                right != phase_set_extents.end()) {
                TargetedWindowGroup isolated = group;
                isolated.first_window = wi;
                isolated.past_last_window = wi + 1;
                isolated.focused_retry = true;
                isolated.region.beg = std::max(chunk.ref_beg, left->second.first);
                isolated.region.end = std::min(chunk.ref_end, right->second.second);
                // Singleton graph blocks need BAM solve context as well as
                // their one anchor. Use the existing flank, bounded by the chunk.
                if (focused_repair && isolated.region.beg == seam.beg)
                    isolated.region.beg = std::max(chunk.ref_beg,
                        seam.beg - kTargetedSolveFlankMin);
                if (focused_repair && isolated.region.end == seam.end)
                    isolated.region.end = std::min(chunk.ref_end,
                        seam.end + kTargetedSolveFlankMin);
                if (isolated.region.beg < seam.beg &&
                    isolated.region.end > seam.end) {
                    Options isolated_opts = sub;
                    isolated_opts.retry_windows = {{seam.beg, seam.end}};
                    if (!dump_prefix.empty())
                        isolated_opts.phase_matrix_dump_prefix =
                            focused_prefix;
                    PhasingChunk local = process_chunk(
                        isolated.region, isolated_opts, context);
                    bool local_raw_preserve = false, local_raw_ordinary = false;
                    std::optional<size_t> local_raw_sparse;
                    const bool local_raw_retry = source_seam_needs_unplaced_msa(
                        local, isolated, windows, opts, local_raw_preserve,
                        local_raw_ordinary, local_raw_sparse, true, std::nullopt,
                        fallback);
                    backfill(local, isolated, isolated_opts);
                    auto local_profiles = local.read_var_profile;
                    auto local_index = std::move(local.read_var_cr);
                    backfill_msa_retry_deletions(
                        local, isolated_opts, seam.beg, seam.end, false);
                    bool local_preserve = false;
                    bool local_ordinary = false;
                    std::optional<size_t> local_sparse;
                    const bool local_filled_retry = source_seam_needs_unplaced_msa(
                        local, isolated, windows, opts, local_preserve,
                        local_ordinary, local_sparse);
                    if (std::binary_search(isolated_conflict_windows.begin(),
                                           isolated_conflict_windows.end(), wi))
                        backfill_msa_retry_deletions(
                            local, isolated_opts, seam.beg, seam.end, true);
                    const bool internal_conflict =
                        source_seam_has_msa_conflict(local, seam, opts);
                    local.read_var_profile = std::move(local_profiles);
                    local.read_var_cr = std::move(local_index);
                    if (local_filled_retry ||
                        (local_raw_retry && !local_raw_ordinary) || internal_conflict) {
                        isolated_opts.add_unplaced_msa_observations = true;
                        if ((local_raw_retry && !local_raw_ordinary) ||
                            internal_conflict) {
                            isolated_opts.joint_het_orientation = true;
                            isolated_opts.link_by_alleles = true;
                        }
                        if (!dump_prefix.empty())
                            isolated_opts.phase_matrix_dump_prefix =
                                focused_prefix + ".msa";
                        PhasingChunk local_retry = process_chunk(
                            isolated.region, isolated_opts, context);
                        if (isolated_opts.joint_het_orientation) {
                            for (CandidateVariant& candidate : local_retry.candidates) {
                                if (candidate.msa_verified &&
                                    candidate.key.type != VariantType::Snp &&
                                    candidate.lcd_var_i_to_cate == kCandNoisyCandHet)
                                    candidate.read_rescue_requires_validation = true;
                            }
                        }
                        // Certify the matrix that transfer will actually use.
                        // Post-solve CIGAR calls can introduce a conflicting
                        // cut even when the raw MSA source path was complete.
                        backfill(local_retry, isolated, isolated_opts);
                        if (!dump_prefix.empty())
                            dump_recovery_phase_state(local_retry, isolated_opts, "trial-source");
                        // A focused replacement must actually carry a
                        // read-supported BAM path across both graph flanks.
                        // Retaining row keys alone does not prevent an MSA
                        // retry from shifting another established block.
                        std::map<hts_pos_t, std::pair<hts_pos_t, hts_pos_t>>
                            retry_extents;
                        for (const CandidateVariant& candidate :
                             local_retry.candidates) {
                            if (!is_phase_set_anchor(candidate)) continue;
                            const hts_pos_t pos = candidate.key.sort_pos();
                            const auto [it, inserted] = retry_extents.try_emplace(
                                candidate.phase_set, pos, pos);
                            if (!inserted) {
                                it->second.first =
                                    std::min(it->second.first, pos);
                                it->second.second =
                                    std::max(it->second.second, pos);
                            }
                        }
                        bool complete_source_span = false;
                        for (const auto& [source_ps, extent] : retry_extents) {
                            if (extent.first > seam.beg + 1 ||
                                extent.second < seam.end)
                                continue;
                            const SourcePathEvidence path =
                                source_phase_set_path_evidence(
                                    local_retry, source_ps);
                            if (path.site_count >= 2 &&
                                path.weak_cuts.empty()) {
                                if (fallback) {
                                    // A newly selected source can retain old
                                    // keys yet orient its added MSA rows poorly.
                                    // Independently certify the flank relation
                                    // with the existing physical SNP cut check.
                                    std::optional<size_t> previous_snp;
                                    bool reaches_right = false;
                                    for (size_t ci = 0; ci < local_retry.candidates.size(); ++ci) {
                                        const CandidateVariant& site = local_retry.candidates[ci];
                                        if (site.phase_set != source_ps ||
                                            !is_phase_set_anchor(site) ||
                                            site.counts.category != VariantCategory::CleanHetSnp)
                                            continue;
                                        const hts_pos_t pos = site.key.sort_pos();
                                        if (pos <= seam.beg) {
                                            previous_snp = ci;
                                            continue;
                                        }
                                        if (!previous_snp ||
                                            !high_quality_snp_cut_support(local_retry,
                                                *previous_snp, ci, source_ps))
                                            break;
                                        previous_snp = ci;
                                        if (pos >= seam.end) {
                                            reaches_right = true;
                                            break;
                                        }
                                    }
                                    if (!reaches_right) continue;
                                }
                                complete_source_span = true;
                                break;
                            }
                        }
                        // Independent BAM blocks may repair part of a graph
                        // seam. Retain every old supported edge and clean SNP
                        // gauge; at least one weak internal edge must improve.
                        const bool partial_source_improvement =
                            !complete_source_span &&
                            internal_conflict &&
                            preserves_hets(local, local_retry, true) &&
                            improves_source_paths(local, local_retry, seam);
                        // A newly exposed conflict beside a singleton cannot
                        // authorize rearranging an established clean-SNP gauge
                        // elsewhere in the solve. Key retention alone misses
                        // an internal reversal or a split of that old block.
                        const bool preserve_clean_gauge =
                            std::binary_search(isolated_conflict_windows.begin(),
                                               isolated_conflict_windows.end(), wi);
                        if ((complete_source_span || partial_source_improvement) &&
                            preserves_hets(local, local_retry, preserve_clean_gauge)) {
                            local = std::move(local_retry);
                            // The focused source owns this seam's observations
                            // and gauge. The original solve keeps its other
                            // seams, without repeating the broad MSA.
                            groups.push_back(isolated);
                            source_paths.emplace_back();
                            apply_selected_msa_observations(local, graph_chunk,
                                isolated_opts, source_paths.back());
                            discovered.push_back(std::move(local));
                            if (wi == group.first_window)
                                remainder.first_window = wi + 1;
                            else if (wi + 1 == group.past_last_window)
                                remainder.past_last_window = wi;
                            else
                                remainder.omitted_window = wi;
                            focused_accepted = true;
                            if (!dump_prefix.empty()) {
                                Options dump_opts = sub;
                                dump_opts.phase_matrix_dump_prefix = focused_prefix;
                                dump_recovery_phase_state(discovered.back(), dump_opts,
                                                          "recovery-source");
                            }
                            break;
                        }
                    }
                }
            }
        }
        // Nonunique boundaries were admitted for a focused, certified solve.
        // If it fails, they cannot also authorize the broad fallback, which
        // lacks the focused path certificate. Reapply grouped admission with
        // its original unique-boundary requirement before replacing source.
        if (!focused_accepted && !internal_retry &&
            (ordinary_retry ||
             (needs_retry && group.past_last_window - group.first_window == 1 &&
              source_seam_needs_unplaced_msa(source, group, windows, opts,
                  preserve_source_rows, ordinary_retry, sparse_window,
                  false, std::nullopt, false, &conflicting_majority, false)))) {
            source_opts.add_unplaced_msa_observations = true;
            if (!dump_prefix.empty())
                source_opts.phase_matrix_dump_prefix = dump_prefix + ".msa";
            PhasingChunk retried = process_chunk(group.region, source_opts, context);
            if (conflicting_majority) {
                for (CandidateVariant& candidate : retried.candidates)
                    if (candidate.msa_verified && is_phase_set_anchor(candidate))
                        candidate.read_rescue_requires_validation = true;
            }
            if (!preserve_source_rows || preserves_hets(source, retried, conflicting_majority)) {
                source = std::move(retried);
                remainder.preserve_source_flanks = conflicting_majority;
                backfill(source, group, source_opts);
            }
        }
        if (remainder.first_window == remainder.past_last_window) continue;
        // Keep new fixed-consensus evidence out of retry admission. A repaired
        // earlier seam must not displace the established solve for a later one.
        // The selected source's original missing calls were queued before any
        // CIGAR backfill, so their verified MSA alleles now take precedence.
        source_paths.emplace_back();
        apply_selected_msa_observations(source, graph_chunk, source_opts, source_paths.back());
        groups.push_back(std::move(remainder));
        discovered.push_back(std::move(source));
        if (!dump_prefix.empty()) {
            Options dump_opts = sub;
            dump_opts.phase_matrix_dump_prefix = dump_prefix;
            dump_recovery_phase_state(discovered.back(), dump_opts,
                                      "recovery-source");
        }
    }

    // Match graph alleles by their selected reference sequence, not their
    // graph walk. This index is also used below for candidate transfer.
    CandidateIndex parent_seq_index;
    CandidateIndex parent_cand_index;
    std::vector<uint8_t> suffix_padded_parents(chunk.candidates.size(), 0);
    for (size_t ci = 0; ci < chunk.candidates.size(); ++ci) {
        parent_cand_index.emplace(cand_key_of(chunk.candidates[ci]), ci);
        const std::string* alt = selected_graph_candidate_alt(graph_chunk, ci);
        if (alt == nullptr) continue;
        const GraphSiteMeta& meta = graph_chunk.site_meta[ci];
        const VariantKey translated =
            vcf_to_variant_key(solve_tid, meta.pos, meta.ref, *alt);
        const CandKey key{translated.sort_pos(), static_cast<int>(translated.type),
                          translated.ref_len, translated.alt};
        parent_seq_index.emplace(key, ci);
        if (meta.ref.size() != alt->size()) {
            size_t prefix = 0;
            while (prefix < std::min(meta.ref.size(), alt->size()) &&
                   meta.ref[prefix] == (*alt)[prefix])
                ++prefix;
            const size_t remaining = meta.ref.size() + alt->size() - 2 * prefix;
            if (static_cast<size_t>(translated.ref_len) + translated.alt.size() < remaining)
                suffix_padded_parents[ci] = 1;
        }
    }

    // A BAM homozygote observed on both source haplotypes is not a usable
    // heterozygous graph seam anchor. Keep its genotype row, but exclude it
    // from the phase-set boundaries and recompute those boundaries only when
    // the existing targeted solve still covers the expanded seam.
    std::vector<size_t> disputed_anchors;
    for (const PhasingChunk& source : discovered) {
        for (size_t ci = 0; ci < source.candidates.size(); ++ci) {
            if (!supported_bam_homozygous_snp(source, ci)) continue;
            const ParentCandidateMatch match = find_parent_candidate(
                parent_cand_index, parent_seq_index,
                cand_key_of(source.candidates[ci]), chunk.candidates.size());
            if (match.index >= chunk.candidates.size()) continue;
            const CandidateVariant& graph_site = chunk.candidates[match.index];
            if (graph_site.counts.category != VariantCategory::CleanHetSnp ||
                !is_phase_set_anchor(graph_site))
                continue;
            const hts_pos_t pos = recovery_position(graph_chunk, match.index);
            const bool is_right_anchor = std::any_of(
                windows.begin(), windows.end(), [pos](const RecoverySeam& seam) {
                    return seam.end == pos;
                });
            if (is_right_anchor &&
                std::find(disputed_anchors.begin(), disputed_anchors.end(),
                          match.index) == disputed_anchors.end())
                disputed_anchors.push_back(match.index);
        }
    }
    if (!disputed_anchors.empty()) {
        std::vector<hts_pos_t> old_phase_sets;
        old_phase_sets.reserve(disputed_anchors.size());
        for (const size_t ci : disputed_anchors) {
            old_phase_sets.push_back(chunk.candidates[ci].phase_set);
            chunk.candidates[ci].phase_set = kUnsetCandidatePhaseSet;
        }
        std::vector<RecoverySeam> expanded = collect_phase_set_seams(graph_chunk);
        bool covered = expanded.size() == windows.size();
        for (size_t wi = 0; covered && wi < expanded.size(); ++wi) {
            const RecoverySeam& before = windows[wi];
            const RecoverySeam& after = expanded[wi];
            const auto group = std::find_if(
                groups.begin(), groups.end(), [wi](const TargetedWindowGroup& value) {
                    return group_owns_window(value, wi);
                });
            covered = group != groups.end() &&
                before.left_phase_set == after.left_phase_set &&
                before.right_phase_set == after.right_phase_set &&
                group->region.beg <= after.beg &&
                after.end <= group->region.end;
        }
        if (covered) {
            windows = std::move(expanded);
            graph_chunk.recovery_windows = windows;
        } else {
            for (size_t i = 0; i < disputed_anchors.size(); ++i)
                chunk.candidates[disputed_anchors[i]].phase_set = old_phase_sets[i];
        }
    }

    // Each targeted solve has its own numeric HP gauge. BAM chunks retain
    // coordinate order, while graph reads use qname order, so a two-pointer
    // merge silently loses matches. Index the stable parent qname strings once
    // and reuse the lookup for every targeted solve in this graph chunk.
    std::unordered_map<std::string_view, size_t> parent_read_by_qname;
    parent_read_by_qname.reserve(chunk.reads.size());
    for (size_t read_i = 0; read_i < chunk.reads.size(); ++read_i)
        parent_read_by_qname.try_emplace(chunk.reads[read_i].qname, read_i);

    using SourceBlockGauge = std::pair<hts_pos_t, hts_pos_t>;
    std::vector<std::map<hts_pos_t, std::pair<int, int>>> gauge_votes(
        discovered.size());
    using DiploidGaugeCounts = std::array<std::array<int, 2>, 2>;
    std::vector<std::map<SourceBlockGauge, DiploidGaugeCounts>>
        block_gauge_votes(discovered.size());
    struct SharedCandidateVotes {
        int same = 0;
        int cross = 0;
        hts_pos_t first_pos = 0;
        bool has_first = false;
        bool has_distinct_loci = false;
    };
    std::vector<std::map<SourceBlockGauge, SharedCandidateVotes>>
        shared_candidate_votes(discovered.size());
    for (size_t gi = 0; gi < discovered.size(); ++gi) {
        const PhasingChunk& src = discovered[gi];
        for (size_t source_i = 0; source_i < src.reads.size(); ++source_i) {
            const auto parent =
                parent_read_by_qname.find(src.reads[source_i].qname);
            if (parent == parent_read_by_qname.end()) continue;
            const size_t parent_i = parent->second;
            if (source_i >= src.haps.size() ||
                source_i >= src.phase_sets.size() ||
                (src.haps[source_i] != 1 && src.haps[source_i] != 2) ||
                src.phase_sets[source_i] <= 0) {
                continue;
            }

            // Preserve the ordinary whole-read graph/BAM vote when both
            // solvers emitted one.
            hts_pos_t tagged_graph_phase_set = kUnphasedReadPhaseSet;
            if (parent_i < chunk.haps.size() &&
                parent_i < chunk.phase_sets.size() &&
                chunk.phase_sets[parent_i] > 0 &&
                (chunk.haps[parent_i] == 1 ||
                 chunk.haps[parent_i] == 2)) {
                tagged_graph_phase_set = chunk.phase_sets[parent_i];
                auto& vote = gauge_votes[gi][tagged_graph_phase_set];
                if (chunk.haps[parent_i] == src.haps[source_i])
                    ++vote.first;
                else
                    ++vote.second;
                auto& block_vote =
                    block_gauge_votes[gi][SourceBlockGauge{
                        tagged_graph_phase_set,
                        src.phase_sets[source_i]}];
                ++block_vote[
                    static_cast<size_t>(chunk.haps[parent_i] - 1)]
                    [static_cast<size_t>(src.haps[source_i] - 1)];
            }

        }
    }

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

    // Every sub-solve candidate and what the merge decided about it, so a
    // site that is found and then silently dropped is visible in output
    // rather than only under a probe.
    std::vector<RecoveredCandidate> audit;
    std::map<CandKey, size_t> audit_of;
    std::map<CandKey, TransferredCandidate> new_cands;
    std::map<size_t, TransferredCandidate> shared_msa_genotypes;
    std::set<size_t> ambiguous_shared_msa_genotypes;
    std::map<std::string, AlleleByCand> observed;
    std::map<std::string, std::set<CandKey>> owned_calls;
    std::map<std::string, std::map<CandKey, uint8_t>> observed_snp_qualities;
    std::map<std::string, int> observed_mapq;
    std::map<std::string, TransferredReadPhase> observed_phase;
    std::map<std::string, TransferredReadPhase> overlay_phase;
    std::set<std::string> ambiguous_overlay_reads;
    struct SourceSite {
        CandKey key;
        size_t solve_id;
        size_t candidate_index;
        hts_pos_t phase_set;
        int hap1_allele;
        int hap2_allele;
        hts_pos_t graph_phase_set;
        int graph_hap1_allele;
        bool graph_clean_snp;
        bool graph_clean_indel;
        bool can_adopt;
    };
    std::vector<SourceSite> source_sites;

    // Select complete BAM phase sets that contribute a phased site to a seam.
    std::vector<std::set<hts_pos_t>> selected_source_phase_sets(discovered.size());
    for (size_t gi = 0; gi < discovered.size(); ++gi) {
        for (const CandidateVariant& cand : discovered[gi].candidates) {
            if (!is_phase_set_anchor(cand))
                continue;
            if (find_containing_window(
                    windows, groups[gi], cand.key.sort_pos()) != nullptr)
                selected_source_phase_sets[gi].insert(cand.phase_set);
        }
    }

    // A validated focused replacement or conflicting-majority retry must
    // carry its whole selected source path through transfer. Dropping private flank rows truncates evidence
    // and can remove sites that a later recovery pass would otherwise find.
    // Cache each path once for admission here and stitch metadata below.
    for (size_t gi = 0; gi < discovered.size(); ++gi) {
        for (const hts_pos_t ps : selected_source_phase_sets[gi])
            if (source_paths[gi].count(ps) == 0)
                source_paths[gi].emplace(ps, source_phase_set_path_evidence(discovered[gi], ps));
    }

    for (size_t gi = 0; gi < discovered.size(); ++gi) {
        PhasingChunk& src = discovered[gi];
        if (src.read_var_profile.size() != src.reads.size()) continue;
        const size_t source_id = gi;
        const auto containing_window = [&](hts_pos_t pos) {
            return find_containing_window(windows, groups[gi], pos);
        };
        for (size_t ci = 0; ci < src.candidates.size(); ++ci) {
            const CandidateVariant& cand = src.candidates[ci];
            const CandKey key = cand_key_of(cand);
            if (group_omits_position(windows, groups[gi], key.pos)) continue;
            const ParentCandidateMatch parent =
                find_parent_candidate(parent_cand_index, parent_seq_index, key,
                                      chunk.candidates.size());

            const auto* member = containing_window(key.pos);
            if (is_phase_set_anchor(cand)) {
                const CandidateVariant* matched =
                    parent.index < chunk.candidates.size()
                        ? &chunk.candidates[parent.index] : nullptr;
                const bool comparable = matched != nullptr &&
                    matched->phase_set > 0 &&
                    matched->counts.n_uniq_alles == 2 &&
                    cand.counts.n_uniq_alles == 2 &&
                    matched->hap_to_cons_alle[1] >= 0 &&
                    matched->hap_to_cons_alle[2] >= 0 &&
                    matched->hap_to_cons_alle[1] != matched->hap_to_cons_alle[2];
                const bool can_adopt = cand.counts.n_uniq_alles == 2 &&
                    cand.hap_to_cons_alle[1] < 2 &&
                    cand.hap_to_cons_alle[2] < 2 &&
                    (matched == nullptr || matched->counts.n_uniq_alles == 2);
                const bool graph_clean_snp = comparable &&
                    cand.key.type == VariantType::Snp &&
                    cand.counts.category == VariantCategory::CleanHetSnp &&
                    matched->key.type == VariantType::Snp &&
                    matched->counts.category == VariantCategory::CleanHetSnp;
                const bool graph_clean_indel = comparable &&
                    cand.key.type != VariantType::Snp &&
                    cand.counts.category == VariantCategory::CleanHetIndel &&
                    matched->key.type == cand.key.type &&
                    matched->counts.category == VariantCategory::CleanHetIndel;
                if (selected_source_phase_sets[gi].count(cand.phase_set) != 0 ||
                    graph_clean_snp || graph_clean_indel)
                    source_sites.push_back(SourceSite{key, gi, ci, cand.phase_set,
                        cand.hap_to_cons_alle[1], cand.hap_to_cons_alle[2],
                        comparable ? matched->phase_set : 0,
                        comparable ? matched->hap_to_cons_alle[1] : -1,
                        graph_clean_snp, graph_clean_indel, can_adopt});
            }
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
                    rec.win_beg = member->beg;
                    rec.win_end = member->end;
                }
                const bool inserted = audit_of.emplace(key, audit.size()).second;
                if (inserted) audit.push_back(std::move(rec));
            }

            if (parent.index < chunk.candidates.size()) {
                // The same clean heterozygote in both representations is a
                // direct gauge anchor: allele 0 and allele 1 have identical
                // sequence meaning after find_parent_candidate translated the
                // graph walk. Preserve that relationship before the merge
                // keeps the graph row and discards the BAM row's consensus.
                //
                // This is a consensus anchor, not another independent read, so
                // keep it separate from the molecule counts below.
                const CandidateVariant& graph_cand =
                    chunk.candidates[parent.index];
                const auto is_clean_het = [](const CandidateVariant& value) {
                    return (value.lcd_var_i_to_cate == kCandCleanHetSnp ||
                            value.lcd_var_i_to_cate == kCandCleanHetIndel) &&
                           value.phase_set > 0 &&
                           value.hap_to_cons_alle[1] >= 0 &&
                           value.hap_to_cons_alle[2] >= 0 &&
                           value.hap_to_cons_alle[1] !=
                               value.hap_to_cons_alle[2];
                };
                if (is_clean_het(graph_cand) && is_clean_het(cand)) {
                    auto& vote = shared_candidate_votes[gi][SourceBlockGauge{
                        graph_cand.phase_set, cand.phase_set}];
                    int* matched = nullptr;
                    if (graph_cand.hap_to_cons_alle[1] ==
                            cand.hap_to_cons_alle[1] &&
                        graph_cand.hap_to_cons_alle[2] ==
                            cand.hap_to_cons_alle[2]) {
                        matched = &vote.same;
                    } else if (graph_cand.hap_to_cons_alle[1] ==
                                   cand.hap_to_cons_alle[2] &&
                               graph_cand.hap_to_cons_alle[2] ==
                                   cand.hap_to_cons_alle[1]) {
                        matched = &vote.cross;
                    }
                    if (matched != nullptr) {
                        ++*matched;
                        if (vote.has_first)
                            vote.has_distinct_loci |= vote.first_pos != key.pos;
                        else {
                            vote.first_pos = key.pos;
                            vote.has_first = true;
                        }
                    }
                }
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
                const auto path = source_paths[gi].find(cand.phase_set);
                const bool complete_source_anchor =
                    (groups[gi].focused_retry || groups[gi].preserve_source_flanks) &&
                    is_phase_set_anchor(cand) && path != source_paths[gi].end() &&
                    path->second.site_count >= 2 && path->second.weak_cuts.empty();
                const bool isolated_source_site = member != nullptr &&
                    path != source_paths[gi].end() &&
                    bam_site_has_only_weak_links(src.candidates, ci, path->second.weak_cuts);
                if ((member != nullptr || complete_source_anchor) &&
                    (suffix_padded_parents[parent.index] != 0 || isolated_source_site) &&
                    graph_cand.graph_site && !graph_cand.bam_injected &&
                    graph_cand.phase_set <= 0 &&
                    graph_cand.counts.category == VariantCategory::RepeatHetIndel &&
                    graph_cand.counts.n_uniq_alles == 2 &&
                    cand.counts.category == VariantCategory::NoisyCandHet &&
                    cand.msa_verified && is_phase_set_anchor(cand) &&
                    cand.counts.n_uniq_alles == 2 &&
                    selected_source_phase_sets[gi].count(cand.phase_set) != 0) {
                    const auto [it, inserted] = shared_msa_genotypes.emplace(
                        parent.index, TransferredCandidate{cand, gi, cand.phase_set});
                    if (!inserted &&
                        (it->second.source_id != gi ||
                         it->second.source_phase_set != cand.phase_set ||
                         it->second.candidate.hap_to_cons_alle != cand.hap_to_cons_alle))
                        ambiguous_shared_msa_genotypes.insert(parent.index);
                }

                continue;
            }
            const auto path = source_paths[gi].find(cand.phase_set);
            const bool complete_source_anchor =
                (groups[gi].focused_retry || groups[gi].preserve_source_flanks) &&
                is_phase_set_anchor(cand) && path != source_paths[gi].end() &&
                path->second.site_count >= 2 && path->second.weak_cuts.empty();
            if (member == nullptr && !complete_source_anchor) continue;
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
                key, TransferredCandidate{cand, source_id, cand.phase_set});
            { auto ai = audit_of.find(key);
              if (ai != audit_of.end()) audit[ai->second].appended = true; }
        }
        for (size_t ri = 0; ri < src.reads.size(); ++ri) {
            ReadVariantProfile& prof = src.read_var_profile[ri];
            if (prof.start_var_idx < 0) continue;
            prof.bam_base_qualities.assign(prof.alleles.size(), 0);
            AlleleByCand& per_read = observed[src.reads[ri].qname];
            observed_mapq[src.reads[ri].qname] = src.reads[ri].mapq;
            if (ri < src.haps.size() && ri < src.phase_sets.size() &&
                src.haps[ri] > 0 && src.phase_sets[ri] > 0) {
                observed_phase.emplace(
                    src.reads[ri].qname,
                    TransferredReadPhase{src.haps[ri], src.phase_sets[ri], source_id});
                if (selected_source_phase_sets[gi].count(src.phase_sets[ri]) != 0) {
                    const auto [it, inserted] = overlay_phase.emplace(
                        src.reads[ri].qname,
                        TransferredReadPhase{src.haps[ri], src.phase_sets[ri], source_id});
                    if (!inserted &&
                        (it->second.hap != src.haps[ri] ||
                         it->second.phase_set != src.phase_sets[ri] ||
                         it->second.source_id != source_id))
                        ambiguous_overlay_reads.insert(src.reads[ri].qname);
                }
            }
            for (size_t k = 0; k < prof.alleles.size(); ++k) {
                const size_t ci = static_cast<size_t>(prof.start_var_idx) + k;
                if (ci >= src.candidates.size()) break;
                if (prof.alleles[k] < 0) continue;
                const VariantKey& key = src.candidates[ci].key;
                uint8_t quality = 0;
                if (key.type == VariantType::Snp && key.ref_len == 1 &&
                    key.alt.size() == 1 && key.pos >= src.ref_beg &&
                    key.pos - src.ref_beg < static_cast<hts_pos_t>(src.ref_seq.size())) {
                    quality = bam_snp_observation_quality(src.reads[ri].alignment.get(),
                        key.pos, src.ref_seq[static_cast<size_t>(key.pos - src.ref_beg)],
                        key.alt[0], prof.alleles[k]);
                    prof.bam_base_qualities[k] = quality;
                }
                const CandKey observed_key = cand_key_of(src.candidates[ci]);
                if (group_omits_position(windows, groups[gi], observed_key.pos)) continue;
                const auto owner = new_cands.find(observed_key);
                const bool owns_call = owner != new_cands.end() &&
                    owner->second.source_id == gi;
                const auto read_owned_calls = owned_calls.find(src.reads[ri].qname);
                const bool already_owned = read_owned_calls != owned_calls.end() &&
                    read_owned_calls->second.count(observed_key) != 0;
                const auto call = std::make_pair(prof.alleles[k],
                    k < prof.alt_qi.size() ? prof.alt_qi[k] : 0);
                auto retained_observation = per_read.emplace(observed_key, call);
                if (owns_call && !already_owned) {
                    // The selected source owns this row's genotype and gauge.
                    // Transfer its call as part of that complete BAM block.
                    retained_observation.first->second = call;
                    owned_calls[src.reads[ri].qname].insert(observed_key);
                    observed_snp_qualities[src.reads[ri].qname].erase(observed_key);
                } else if ((!already_owned || owns_call) && !retained_observation.second &&
                    retained_observation.first->second.first != prof.alleles[k]) {
                    // Both source matrices retain their own calls. The merged
                    // slot abstains permanently instead of choosing by order.
                    retained_observation.first->second = {kConflictingBamAllele, -1};
                    observed_snp_qualities[src.reads[ri].qname].erase(observed_key);
                }
                if (quality > 0 && retained_observation.first->second.first >= 0 &&
                    (!already_owned || owns_call))
                    observed_snp_qualities[src.reads[ri].qname].emplace(observed_key, quality);
            }
        }
        for (const auto& rejection : src.rejected_msa_observations) {
            const CandKey key{rejection.key.sort_pos(), static_cast<int>(rejection.key.type),
                              rejection.key.ref_len, rejection.key.alt};
            const std::string& qname = src.reads[rejection.read_id].qname;
            const auto owner = new_cands.find(key);
            if (owner != new_cands.end() && owner->second.source_id != gi &&
                owned_calls[qname].count(key) != 0) continue;
            if (!group_omits_position(windows, groups[gi], key.pos))
                observed[qname].insert_or_assign(
                    key, std::make_pair(kConflictingBamAllele, -1));
            observed_snp_qualities[qname].erase(key);
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
    if (new_cands.empty() && refreshed == 0) return false;

    // A sub-solve's PS coordinate is meaningful only inside that solve. It can
    // numerically collide with a graph PS while using the opposite HP gauge.
    // Give each imported local block an unused anchor based on its first
    // transferred candidate; candidates and reads share this source-scoped map.
    using SourcePhaseSet = std::pair<size_t, hts_pos_t>;
    std::set<hts_pos_t> used_phase_sets;
    for (const CandidateVariant& candidate : chunk.candidates)
        if (candidate.phase_set > 0) used_phase_sets.insert(candidate.phase_set);
    for (const hts_pos_t phase_set : chunk.phase_sets)
        if (phase_set > 0) used_phase_sets.insert(phase_set);
    std::map<SourcePhaseSet, hts_pos_t> phase_set_remap;
    std::set<hts_pos_t> imported_phase_sets;
    for (size_t gi = 0; gi < selected_source_phase_sets.size(); ++gi) {
        for (const hts_pos_t source_ps : selected_source_phase_sets[gi]) {
            hts_pos_t first_pos = std::numeric_limits<hts_pos_t>::max();
            for (const SourceSite& site : source_sites)
                if (site.solve_id == gi && site.phase_set == source_ps)
                    first_pos = std::min(first_pos, site.key.pos);
            if (first_pos == std::numeric_limits<hts_pos_t>::max()) continue;
            hts_pos_t mapped = first_pos;
            while (mapped <= 0 || used_phase_sets.count(mapped) != 0)
                ++mapped;
            used_phase_sets.insert(mapped);
            phase_set_remap.emplace(SourcePhaseSet{gi, source_ps}, mapped);
        }
    }
    graph_chunk.recovery_source_path_supported.clear();
    graph_chunk.recovery_source_weak_cuts.clear();
    graph_chunk.recovery_source_quality_cuts.clear();
    for (const auto& [source, mapped] : phase_set_remap) {
        if (selected_source_phase_sets[source.first].count(source.second) == 0)
            continue;
        SourcePathEvidence path = std::move(source_paths[source.first].at(source.second));
        graph_chunk.recovery_source_path_supported.emplace(
            mapped, path.site_count >= 2 && path.weak_cuts.empty());
        graph_chunk.recovery_source_weak_cuts.emplace(
            mapped, std::move(path.weak_cuts));
        graph_chunk.recovery_source_quality_cuts.emplace(
            mapped, std::move(path.quality_supported_cuts));
    }
    for (auto& [key, transferred] : new_cands) {
        (void)key;
        if (transferred.source_phase_set <= 0) {
            transferred.candidate.phase_set = kUnsetCandidatePhaseSet;
            continue;
        }
        const SourcePhaseSet source{transferred.source_id,
                                    transferred.source_phase_set};
        auto mapped = phase_set_remap.find(source);
        if (mapped == phase_set_remap.end()) {
            hts_pos_t phase_set = transferred.candidate.key.sort_pos();
            while (phase_set <= 0 || used_phase_sets.count(phase_set) != 0)
                ++phase_set;
            used_phase_sets.insert(phase_set);
            mapped = phase_set_remap.emplace(source, phase_set).first;
            imported_phase_sets.insert(phase_set);
        }
        transferred.candidate.phase_set = mapped->second;
        imported_phase_sets.insert(mapped->second);
    }

    // Sequence matching must preserve the independently phased BAM genotype,
    // not merely add observations to an unphased catalog repeat. Adoption is
    // source-scoped and deferred until the unused PS labels exist. Multiple
    // independent claims stay unphased rather than picking a gauge by order.
    std::set<size_t> adopted_msa_indices;
    std::set<hts_pos_t> independent_msa_phase_sets;
    for (const auto& entry : shared_msa_genotypes) {
        if (ambiguous_shared_msa_genotypes.count(entry.first) != 0) continue;
        const TransferredCandidate& transferred = entry.second;
        const auto mapped = phase_set_remap.find(
            SourcePhaseSet{transferred.source_id, transferred.source_phase_set});
        if (mapped == phase_set_remap.end()) continue;
        hts_pos_t target_phase_set = mapped->second;
        if (suffix_padded_parents[entry.first] == 0) {
            // A verified allele is not proof of its source block's orientation.
            // Newly retained shared rows start independently; paired reads must
            // establish their connection before any source gauge is inherited.
            target_phase_set = transferred.candidate.key.sort_pos();
            while (target_phase_set <= 0 || used_phase_sets.count(target_phase_set) != 0)
                ++target_phase_set;
            used_phase_sets.insert(target_phase_set);
        }
        if (adopt_unphased_graph_allele_from_bam(
                chunk.candidates[entry.first], transferred.candidate, target_phase_set)) {
            if (suffix_padded_parents[entry.first] == 0) {
                chunk.candidates[entry.first].read_rescue_requires_validation = true;
                chunk.candidates[entry.first].bam_independent_genotype = true;
                independent_msa_phase_sets.insert(target_phase_set);
            }
            imported_phase_sets.insert(target_phase_set);
            adopted_msa_indices.insert(entry.first);
        }
    }

    graph_chunk.recovery_phase_gauges.clear();
    graph_chunk.recovery_phase_gauges.reserve(groups.size());
    std::ofstream source_evidence;
    if (!opts.phase_matrix_dump_prefix.empty()) {
        const std::string path = opts.phase_matrix_dump_prefix + ".chunk" +
            std::to_string(chunk.region.chunk_id) +
            (completed_seams == nullptr ? ".bam-source-evidence.tsv" : ".bam-source-evidence-retry.tsv");
        source_evidence.open(path);
        if (!source_evidence) throw std::runtime_error("cannot write BAM source evidence: " + path);
        source_evidence << "solve\tpos\ttype\tref_len\talt\tqname\tallele\tquality\tkind\n";
    }
    for (size_t gi = 0; gi < groups.size(); ++gi) {
        RecoveryPhaseGauge gauge;
        gauge.focused_retry = groups[gi].focused_retry;
        gauge.beg = groups[gi].region.beg;
        gauge.end = groups[gi].region.end;
        for (const auto& [source, mapped_phase_set] : phase_set_remap)
            if (source.first == gi)
                gauge.imported_phase_sets.push_back(mapped_phase_set);
        std::vector<std::pair<hts_pos_t, hts_pos_t>> source_labels;
        for (const auto& [source, mapped_phase_set] : phase_set_remap)
            if (source.first == gi)
                source_labels.emplace_back(source.second, mapped_phase_set);
        // The strict seam selection controls ownership transfer only. Keep
        // flank-only source blocks too: a block starting at the right boundary
        // can carry spanning reads even though none of its rows is injected.
        retain_recovery_bam_evidence(discovered[gi], source_labels, gauge);
        if (source_evidence.is_open()) {
            for (const RecoveryBamRead& read : gauge.bam_reads)
                for (size_t oi = 0; oi < read.observations.size(); ++oi) {
                    const auto& observation = read.observations[oi];
                    const VariantKey& key = gauge.bam_sites[observation.first].key;
                    source_evidence << gi << '\t' << key.sort_pos() << '\t'
                        << static_cast<int>(key.type) << '\t' << key.ref_len << '\t'
                        << key.alt << '\t' << read.qname << '\t' << observation.second
                        << '\t' << static_cast<int>(read.base_qualities[oi]) << "\tcall\n";
                }
            for (const RecoveryBamRecall& alternative : gauge.conflicting_recalls)
                source_evidence << gi << '\t' << alternative.key.sort_pos() << '\t'
                    << static_cast<int>(alternative.key.type) << '\t' << alternative.key.ref_len
                    << '\t' << alternative.key.alt << '\t' << alternative.qname << '\t'
                    << alternative.allele << "\t0\t"
                    << (alternative.fixed_consensus ? "alternative_msa" : "alternative_physical") << '\n';
        }
        for (const auto& [phase_set, vote] : gauge_votes[gi])
            gauge.graph_votes.push_back(
                PhaseSetGaugeVote{phase_set, vote.first, vote.second});
        for (const auto& [source, vote] : block_gauge_votes[gi]) {
            const auto mapped = phase_set_remap.find(
                SourcePhaseSet{gi, source.second});
            if (mapped == phase_set_remap.end()) continue;
            gauge.block_votes.push_back(RecoveryBlockGaugeVote{
                source.first, mapped->second, vote});
        }
        for (const auto& [source, vote] : shared_candidate_votes[gi]) {
            const auto mapped = phase_set_remap.find(
                SourcePhaseSet{gi, source.second});
            if (mapped == phase_set_remap.end()) continue;
            auto existing = std::find_if(
                gauge.block_votes.begin(), gauge.block_votes.end(),
                [&](const RecoveryBlockGaugeVote& candidate) {
                    return candidate.graph_phase_set == source.first &&
                           candidate.bam_phase_set == mapped->second;
                });
            if (existing == gauge.block_votes.end()) {
                RecoveryBlockGaugeVote candidate_vote;
                candidate_vote.graph_phase_set = source.first;
                candidate_vote.bam_phase_set = mapped->second;
                candidate_vote.shared_candidate_same = vote.same;
                candidate_vote.shared_candidate_cross = vote.cross;
                candidate_vote.has_distinct_shared_loci = vote.has_distinct_loci;
                gauge.block_votes.push_back(std::move(candidate_vote));
            } else {
                existing->shared_candidate_same += vote.same;
                existing->shared_candidate_cross += vote.cross;
                existing->has_distinct_shared_loci |= vote.has_distinct_loci;
            }
        }
        constexpr hts_pos_t kPhysicalBridgeFlankContext = 2000;
        constexpr size_t kMaxPhysicalBridgeSitesPerFlank = 3;
        constexpr double kMaxPhysicalBridgeError = 0.01;
        for (size_t wi = groups[gi].first_window;
             wi < groups[gi].past_last_window; ++wi) {
            if (!group_owns_window(groups[gi], wi)) continue;
            const RecoverySeam& seam = windows[wi];
            std::vector<const SourceSite*> left_candidates;
            std::vector<const SourceSite*> right_candidates;
            std::vector<std::pair<size_t, int>> left_sites;
            std::vector<std::pair<size_t, int>> right_sites;
            for (const SourceSite& site : source_sites) {
                if (site.solve_id != gi || !site.graph_clean_snp)
                    continue;
                const std::pair<size_t, int> anchor{
                    site.candidate_index, site.graph_hap1_allele};
                if (site.graph_phase_set == seam.left_phase_set) {
                    left_sites.push_back(anchor);
                    if (site.key.pos >= seam.beg - kPhysicalBridgeFlankContext &&
                        site.key.pos <= seam.end)
                        left_candidates.push_back(&site);
                } else if (site.graph_phase_set == seam.right_phase_set) {
                    right_sites.push_back(anchor);
                    if (site.key.pos >= seam.beg &&
                        site.key.pos <= seam.end + kPhysicalBridgeFlankContext)
                        right_candidates.push_back(&site);
                }
            }
            std::sort(left_candidates.begin(), left_candidates.end(),
                      [](const SourceSite* a, const SourceSite* b) {
                          return a->key.pos > b->key.pos;
                      });
            std::sort(right_candidates.begin(), right_candidates.end(),
                      [](const SourceSite* a, const SourceSite* b) {
                          return a->key.pos < b->key.pos;
                      });
            if (left_candidates.size() > kMaxPhysicalBridgeSitesPerFlank)
                left_candidates.resize(kMaxPhysicalBridgeSitesPerFlank);
            if (right_candidates.size() > kMaxPhysicalBridgeSitesPerFlank)
                right_candidates.resize(kMaxPhysicalBridgeSitesPerFlank);
            if (left_candidates.empty() || right_candidates.empty()) {
                const std::optional<bool> insertion_flip =
                    physical_graph_snp_insertion_bridge(
                        graph_chunk, seam, discovered[gi]);
                if (insertion_flip)
                    gauge.physical_snp_bridges.push_back(
                        RecoveryPhysicalSnpBridge{
                            seam.left_phase_set, seam.right_phase_set,
                            *insertion_flip});
                continue;
            }
            const double max_error = kMaxPhysicalBridgeError /
                static_cast<double>(left_candidates.size() *
                                    right_candidates.size());
            std::optional<bool> supported_flip;
            bool conflict = false;
            for (const SourceSite* left : left_candidates) {
                for (const SourceSite* right : right_candidates) {
                    if (left->key.pos >= right->key.pos ||
                        left->phase_set == right->phase_set)
                        continue;
                    const std::optional<bool> flip =
                        physical_graph_snp_bridge(
                            discovered[gi], left->candidate_index,
                            right->candidate_index,
                            left->graph_hap1_allele,
                            right->graph_hap1_allele,
                            left_sites, right_sites, max_error);
                    if (!flip) continue;
                    if (supported_flip && *supported_flip != *flip) {
                        conflict = true;
                        break;
                    }
                    supported_flip = flip;
                }
                if (conflict) break;
            }
            if (!conflict && supported_flip)
                gauge.physical_snp_bridges.push_back(
                    RecoveryPhysicalSnpBridge{
                        seam.left_phase_set, seam.right_phase_set,
                        *supported_flip});
        }
        for (size_t wi = groups[gi].first_window;
             wi < groups[gi].past_last_window; ++wi) {
            if (!group_owns_window(groups[gi], wi)) continue;
            const RecoverySeam& seam = windows[wi];
            for (const RecoveryPhysicalSnpBridge& bridge : retry_bridges)
                if (bridge.left_phase_set == seam.left_phase_set &&
                    bridge.right_phase_set == seam.right_phase_set)
                    gauge.physical_snp_bridges.push_back(bridge);
        }
        graph_chunk.recovery_phase_gauges.push_back(std::move(gauge));
    }

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
        // Keep the sub-solve consensus in its independent local gauge. The
        // left-to-right stitch decides parity only after the merge has rebuilt
        // all candidate and read observations.
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
        std::map<size_t, int> graph_observations;
        std::map<size_t, std::pair<int, int>> bam_observations;
        std::map<size_t, uint8_t> bam_snp_qualities;
        if (old_prof.start_var_idx >= 0) {
            for (size_t k = 0; k < old_prof.alleles.size(); ++k) {
                const size_t ci = static_cast<size_t>(old_prof.start_var_idx) + k;
                if (ci >= old_to_new.size() || old_to_new[ci] < 0) continue;
                // BAM calls establish the adopted genotype. Newly retained
                // rows also keep independent graph calls for read rescue;
                // those calls must never become source MSA observations.
                if (adopted_msa_indices.count(ci) != 0 &&
                    suffix_padded_parents[ci] != 0) continue;
                const size_t final_i = static_cast<size_t>(old_to_new[ci]);
                const bool bam_owned = adopted_msa_indices.count(ci) != 0;
                if (!bam_owned && old_prof.alleles[k] >= 0)
                    alleles.emplace(final_i,
                                    std::make_pair(old_prof.alleles[k],
                                        k < old_prof.alt_qi.size() ? old_prof.alt_qi[k] : 0));
                const int graph_allele = old_prof.graph_alleles.empty()
                    ? old_prof.alleles[k]
                    : (k < old_prof.graph_alleles.size()
                        ? old_prof.graph_alleles[k] : -1);
                if (graph_allele >= 0) {
                    graph_observations.emplace(final_i, graph_allele);
                    // Retain the independent graph channel without declaring
                    // it a BAM/MSA call. Include its unknown primary slot in
                    // the profile extent, even on graph-only reads.
                    if (bam_owned)
                        alleles.emplace(final_i, std::make_pair(-1, 0));
                }
                if (bam_owned) continue;
                if (k < old_prof.bam_alleles.size() &&
                    (old_prof.bam_alleles[k] >= 0 ||
                     old_prof.bam_alleles[k] == kConflictingBamAllele))
                    bam_observations.emplace(final_i,
                        std::make_pair(old_prof.bam_alleles[k],
                            k < old_prof.bam_qi.size() ? old_prof.bam_qi[k] : 0));
                if (k < old_prof.bam_base_qualities.size() &&
                    old_prof.bam_base_qualities[k] > 0)
                    bam_snp_qualities.emplace(final_i,
                        old_prof.bam_base_qualities[k]);
            }
        }
        auto it = observed.find(chunk.reads[ri].qname);
        const auto read_qualities =
            observed_snp_qualities.find(chunk.reads[ri].qname);
        if (it != observed.end())
            for (const auto& entry : it->second) {
                auto idx = index_of.find(entry.first);
                if (idx == index_of.end()) continue;
                const auto prior_bam = bam_observations.find(idx->second);
                // A retry is another independent solve, not a resolution of
                // the earlier contradictory calls. Preserve that abstention.
                if (prior_bam != bam_observations.end() &&
                    prior_bam->second.first == kConflictingBamAllele)
                    continue;
                const auto graph_call = graph_observations.find(idx->second);
                if (entry.second.first == kConflictingBamAllele && graph_call != graph_observations.end())
                    // A BAM disagreement does not erase an independent graph
                    // call. Keep it in the working profile and expose the BAM
                    // ambiguity separately, without borrowing its quality.
                    alleles.insert_or_assign(idx->second, std::make_pair(graph_call->second, 0));
                else
                    alleles.insert_or_assign(idx->second, entry.second);
                bam_observations.insert_or_assign(idx->second, entry.second);
                // This replay replaces the BAM call, so its certificate must
                // replace the previous one too, including an absent quality.
                bam_snp_qualities.erase(idx->second);
                if (entry.second.first >= 0 &&
                    read_qualities != observed_snp_qualities.end()) {
                    const auto quality = read_qualities->second.find(entry.first);
                    if (quality != read_qualities->second.end())
                        bam_snp_qualities.insert_or_assign(idx->second, quality->second);
                }
            }
        // Channel calls can outlive a missing primary allele. Their complete
        // union defines the common extent; unknown primary slots stay unknown.
        for (const auto& entry : graph_observations)
            alleles.try_emplace(entry.first, std::make_pair(-1, 0));
        for (const auto& entry : bam_observations)
            alleles.try_emplace(entry.first, std::make_pair(-1, 0));
        for (const auto& entry : bam_snp_qualities)
            alleles.try_emplace(entry.first, std::make_pair(-1, 0));
        ReadVariantProfile prof;
        prof.read_id = static_cast<int>(ri);
        const auto source_mapq = observed_mapq.find(chunk.reads[ri].qname);
        prof.bam_mapq = source_mapq != observed_mapq.end()
            ? source_mapq->second : old_prof.bam_mapq;
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
        // Preserve the independent graph and BAM calls. A disagreeing BAM
        // allele may drive the primary recovery profile, but the graph call
        // remains available as conflict evidence for site reconciliation.
        if (!graph_observations.empty()) {
            prof.graph_alleles.assign(span, -1);
            for (const auto& [candidate_i, allele] : graph_observations)
                prof.graph_alleles[candidate_i -
                    static_cast<size_t>(prof.start_var_idx)] = allele;
        }
        if (!bam_observations.empty()) {
            prof.bam_alleles.assign(span, -1);
            prof.bam_qi.assign(span, 0);
            if (!bam_snp_qualities.empty()) {
                prof.bam_base_qualities.assign(span, 0);
                for (const auto& [candidate_i, quality] : bam_snp_qualities)
                    prof.bam_base_qualities[candidate_i -
                        static_cast<size_t>(prof.start_var_idx)] = quality;
            }
            for (const auto& [candidate_i, observation] : bam_observations) {
                const size_t off = candidate_i -
                    static_cast<size_t>(prof.start_var_idx);
                prof.bam_alleles[off] = observation.first;
                prof.bam_qi[off] = observation.second;
            }
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
    graph_chunk.recovery_source_sites.clear();
    for (const SourceSite& site : source_sites) {
        const auto mapped = phase_set_remap.find(
            SourcePhaseSet{site.solve_id, site.phase_set});
        const auto final_index = index_of.find(site.key);
        if (mapped == phase_set_remap.end() || final_index == index_of.end())
            continue;
        graph_chunk.recovery_source_sites.push_back(RecoverySourceSite{
            final_index->second, mapped->second,
            site.hap1_allele, site.hap2_allele,
            site.graph_phase_set, site.graph_hap1_allele,
            site.graph_clean_snp, site.can_adopt &&
                independent_msa_phase_sets.count(
                    chunk.candidates[final_index->second].phase_set) == 0});
    }

    // Padded BAM solves have independent HP gauges. A clean SNP represented
    // identically in both solves and the graph can anchor each source after
    // their separate allele/read votes. Other duplicated rows have no safe
    // shared gauge and cannot be adopted by iteration order.
    std::unordered_map<size_t, size_t> source_claims;
    for (const RecoverySourceSite& site : graph_chunk.recovery_source_sites)
        ++source_claims[site.candidate_index];
    for (RecoverySourceSite& site : graph_chunk.recovery_source_sites)
        if (source_claims[site.candidate_index] != 1 &&
            !site.clean_shared_snp)
            site.can_adopt = false;

    // Every read-indexed vector grows with the reads. Appending without this
    // leaves haps and phase_sets short; the qname re-sort below would then be
    // unable to move those labels with their reads.
    if (chunk.haps.size() != chunk.reads.size()) chunk.haps.resize(chunk.reads.size(), 0);
    if (chunk.phase_sets.size() != chunk.reads.size())
        chunk.phase_sets.resize(chunk.reads.size(), kUnphasedReadPhaseSet);

    // Keep established graph assignments. New and previously unphased reads
    // inherit the independent gauge produced by the BAM sub-solve only when a
    // transferred candidate represents that local phase set.
    for (size_t ri = 0; ri < chunk.reads.size(); ++ri) {
        if (chunk.haps[ri] != 0 || chunk.phase_sets[ri] > 0) continue;
        const auto recovered = observed_phase.find(chunk.reads[ri].qname);
        if (recovered == observed_phase.end()) continue;
        const SourcePhaseSet source{recovered->second.source_id,
                                    recovered->second.phase_set};
        const auto mapped = phase_set_remap.find(source);
        if (mapped == phase_set_remap.end() ||
            imported_phase_sets.count(mapped->second) == 0)
            continue;
        chunk.haps[ri] = recovered->second.hap;
        chunk.phase_sets[ri] = mapped->second;
    }

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

    graph_chunk.recovery_source_reads.clear();
    std::unordered_map<std::string_view, size_t> final_read_by_qname;
    final_read_by_qname.reserve(chunk.reads.size());
    for (size_t ri = 0; ri < chunk.reads.size(); ++ri)
        final_read_by_qname.try_emplace(chunk.reads[ri].qname, ri);
    if (!opts.phase_matrix_dump_prefix.empty()) {
        // Exact source keys resolve through index_of, including graph-walk
        // aliases. Record selected-source calls separately from phase-label
        // adoption so lost observations cannot hide behind a successful join.
        const std::string path = opts.phase_matrix_dump_prefix + ".chunk" +
            std::to_string(chunk.region.chunk_id) +
            (completed_seams == nullptr ? ".transfer.tsv" : ".transfer-retry.tsv");
        std::ofstream trace(path);
        if (!trace) throw std::runtime_error("cannot write recovery transfer trace: " + path);
        trace << "solve\tpos\ttype\tref_len\talt\tflags\tmsa_verified\tqname"
              << "\tsource_allele\tdestination_index\tbam_allele\tprimary_allele\tstatus"
              << "\tsource_phase_set\tphase_anchor\n";
        for (size_t gi = 0; gi < discovered.size(); ++gi) {
            const PhasingChunk& src = discovered[gi];
            for (size_t ri = 0; ri < src.reads.size(); ++ri) {
                const ReadVariantProfile& profile = src.read_var_profile[ri];
                if (profile.start_var_idx < 0) continue;
                const auto read = final_read_by_qname.find(src.reads[ri].qname);
                for (size_t oi = 0; oi < profile.alleles.size(); ++oi) {
                    const size_t ci = static_cast<size_t>(profile.start_var_idx) + oi;
                    if (ci >= src.candidates.size()) break;
                    const CandidateVariant& site = src.candidates[ci];
                    if (profile.alleles[oi] < 0 ||
                        (!site.msa_verified && site.counts.category != VariantCategory::CleanHetSnp &&
                         site.counts.category != VariantCategory::CleanHetIndel)) continue;
                    const CandKey key = cand_key_of(site);
                    const auto destination = index_of.find(key);
                    const bool omitted = group_omits_position(windows, groups[gi], key.pos);
                    int bam_allele = -1, primary_allele = -1;
                    if (destination != index_of.end() && read != final_read_by_qname.end()) {
                        const auto& target = chunk.read_var_profile[read->second];
                        const int offset = static_cast<int>(destination->second) - target.start_var_idx;
                        if (target.start_var_idx >= 0 && offset >= 0) {
                            const size_t k = static_cast<size_t>(offset);
                            if (k < target.bam_alleles.size()) bam_allele = target.bam_alleles[k];
                            if (k < target.alleles.size()) primary_allele = target.alleles[k];
                        }
                    }
                    trace << gi << '\t' << key.pos << '\t' << key.type << '\t' << key.ref_len
                          << '\t' << key.alt << '\t' << site.lcd_var_i_to_cate << '\t'
                          << site.msa_verified << '\t' << src.reads[ri].qname << '\t'
                          << profile.alleles[oi] << '\t'
                          << (destination == index_of.end() ? -1 : static_cast<long>(destination->second))
                          << '\t' << bam_allele << '\t' << primary_allele << '\t'
                          << (omitted ? "other_seam_owner" : destination == index_of.end()
                              ? "candidate_not_transferred" : read == final_read_by_qname.end()
                              ? "read_not_transferred" : "mapped") << '\t'
                          << site.phase_set << '\t' << is_phase_set_anchor(site) << '\n';
                }
            }
        }
    }
    // A child snarl can be callable only on one parent path. In that case its
    // graph REF and ALT calls may both come from deletion-bearing reads, while
    // reads carrying the physical reference base bypass the child entirely.
    // Recover the allele gauge only when a phased BAM deletion spans the SNP,
    // the exact CIGAR observations establish both sides, and the graph ALT is
    // enriched on the deletion side. Only a terminal graph SNP may move to the
    // overlapping BAM block; the rest of an established graph block stays put.
    std::map<hts_pos_t, std::vector<size_t>> terminal_graph_snps;
    std::map<size_t, char> snp_alt_base;
    std::map<hts_pos_t, size_t> graph_ps_site_count;
    std::map<hts_pos_t, hts_pos_t> graph_ps_last_pos;
    for (size_t ci = 0; ci < chunk.candidates.size(); ++ci) {
        const CandidateVariant& candidate = chunk.candidates[ci];
        if (graph_chunk.site_ids[ci].empty() || candidate.phase_set <= 0)
            continue;
        ++graph_ps_site_count[candidate.phase_set];
        auto& last_pos = graph_ps_last_pos[candidate.phase_set];
        last_pos = std::max(last_pos, candidate.key.sort_pos());
    }
    for (size_t ci = 0; ci < chunk.candidates.size(); ++ci) {
        const CandidateVariant& candidate = chunk.candidates[ci];
        if (graph_chunk.site_ids[ci].empty() ||
            candidate.counts.category != VariantCategory::CleanHetSnp ||
            candidate.phase_set <= 0 ||
            candidate.key.sort_pos() !=
                graph_ps_last_pos[candidate.phase_set])
            continue;
        const std::string* alt = selected_graph_candidate_alt(graph_chunk, ci);
        if (alt == nullptr) continue;
        const GraphSiteMeta& meta = graph_chunk.site_meta[ci];
        const VariantKey key =
            vcf_to_variant_key(solve_tid, meta.pos, meta.ref, *alt);
        if (key.type == VariantType::Snp && key.alt.size() == 1) {
            terminal_graph_snps[key.pos].push_back(ci);
            snp_alt_base.emplace(ci, key.alt[0]);
        }
    }
    constexpr double kOverlapAlleleMaxP = 0.01;
    constexpr double kOverlapHapMaxP = 0.05;
    std::set<size_t> projected_snps;
    for (size_t gi = 0; gi < discovered.size(); ++gi) {
        const PhasingChunk& source = discovered[gi];
        std::map<hts_pos_t, std::vector<hts_pos_t>> source_positions;
        for (const CandidateVariant& candidate : source.candidates)
            if (candidate.phase_set > 0 &&
                candidate.hap_to_cons_alle[1] >= 0 &&
                candidate.hap_to_cons_alle[2] >= 0 &&
                candidate.hap_to_cons_alle[1] !=
                    candidate.hap_to_cons_alle[2])
                source_positions[candidate.phase_set].push_back(
                    candidate.key.sort_pos());
        for (auto& [phase_set, positions] : source_positions) {
            (void)phase_set;
            std::sort(positions.begin(), positions.end());
            positions.erase(std::unique(positions.begin(), positions.end()),
                            positions.end());
        }
        for (size_t source_ci = 0; source_ci < source.candidates.size(); ++source_ci) {
            const CandidateVariant& deletion = source.candidates[source_ci];
            if (deletion.key.type != VariantType::Deletion ||
                deletion.key.ref_len <= 0 || deletion.phase_set <= 0 ||
                deletion.hap_to_cons_alle[1] < 0 ||
                deletion.hap_to_cons_alle[2] < 0 ||
                deletion.hap_to_cons_alle[1] == deletion.hap_to_cons_alle[2])
                continue;
            const auto mapped = phase_set_remap.find(
                SourcePhaseSet{gi, deletion.phase_set});
            if (mapped == phase_set_remap.end()) continue;
            const hts_pos_t del_beg = deletion.key.pos;
            const hts_pos_t del_end = del_beg + deletion.key.ref_len - 1;
            const auto& positions = source_positions.at(deletion.phase_set);
            const auto next = std::upper_bound(positions.begin(),
                                               positions.end(), del_beg);
            const hts_pos_t next_source_site =
                next == positions.end() ? 0 : *next;
            const auto cuts = graph_chunk.recovery_source_weak_cuts.find(
                mapped->second);
            if (next_source_site == 0 ||
                (cuts != graph_chunk.recovery_source_weak_cuts.end() &&
                 std::any_of(cuts->second.begin(), cuts->second.end(),
                             [del_beg, next_source_site](hts_pos_t cut) {
                                 return del_beg <= cut &&
                                        cut < next_source_site;
                             })))
                continue;
            for (auto snps = terminal_graph_snps.lower_bound(del_beg);
                 snps != terminal_graph_snps.end() && snps->first <= del_end;
                 ++snps) {
                if (snps->second.size() != 1 ||
                    next_source_site <= snps->first)
                    continue;
                const size_t ci = snps->second.front();
                if (projected_snps.count(ci) != 0) continue;
                const hts_pos_t pos = snps->first;
                const int old_graph_hap = chunk.candidates[ci].hap_to_cons_alle[1];
                const hts_pos_t old_graph_ps = chunk.candidates[ci].phase_set;
                int ref_by_hap[2] = {0, 0};
                int del_by_hap[2] = {0, 0};
                int graph_del_alt = 0, graph_del_ref = 0;
                int graph_ref_alt = 0, graph_ref_ref = 0;
                int physical_alt = 0;
                std::array<std::array<int, 2>, 2> source_event_calls{};
                int source_ref = 0, source_del = 0;
                // -1 = uncallable, 0 = exact reference, 1 = deletion,
                // 2 = physical SNP ALT. Cache these calls for read transfer.
                std::vector<int8_t> physical_calls(source.reads.size(), -1);
                for (size_t ri = 0; ri < source.reads.size(); ++ri) {
                    const ReadRecord& read = source.reads[ri];
                    if (read.is_skipped) continue;
                    const GraphSiteMeta& meta = graph_chunk.site_meta[ci];
                    const size_t ref_offset =
                        static_cast<size_t>(pos - meta.pos);
                    if (ref_offset >= meta.ref.size()) continue;
                    const int physical = physical_snp_call(
                        read, pos, meta.ref[ref_offset], snp_alt_base.at(ci));
                    if (physical < 0) continue;
                    physical_calls[ri] = static_cast<int8_t>(physical);
                    if (physical == 2) { ++physical_alt; continue; }
                    if (ri < source.phase_sets.size() &&
                        ri < source.haps.size() &&
                        source.phase_sets[ri] == deletion.phase_set &&
                        (source.haps[ri] == 1 || source.haps[ri] == 2)) {
                        int& count = physical == 0
                            ? ref_by_hap[source.haps[ri] - 1]
                            : del_by_hap[source.haps[ri] - 1];
                        ++count;
                    }
                    if (ri < source.read_var_profile.size()) {
                        const ReadVariantProfile& profile =
                            source.read_var_profile[ri];
                        const int offset = static_cast<int>(source_ci) -
                                           profile.start_var_idx;
                        if (profile.start_var_idx >= 0 && offset >= 0 &&
                            static_cast<size_t>(offset) < profile.alleles.size()) {
                            const int allele = profile.alleles[offset];
                            if (allele == 0 || allele == 1) {
                                ++source_event_calls[physical][allele];
                                if (physical == 0) ++source_ref;
                                else ++source_del;
                            }
                        }
                    }
                    const auto parent = final_read_by_qname.find(read.qname);
                    if (parent == final_read_by_qname.end()) continue;
                    const ReadVariantProfile& profile =
                        chunk.read_var_profile[parent->second];
                    const int offset = static_cast<int>(ci) -
                                       profile.start_var_idx;
                    if (profile.start_var_idx < 0 || offset < 0 ||
                        static_cast<size_t>(offset) >= profile.graph_alleles.size()) {
                        if (physical == 0) ++graph_ref_ref;
                        continue;
                    }
                    const int graph_allele = profile.graph_alleles[offset];
                    if (physical == 1) {
                        if (graph_allele == 1) ++graph_del_alt;
                        else if (graph_allele == 0) ++graph_del_ref;
                    } else {
                        if (graph_allele == 1) ++graph_ref_alt;
                        else ++graph_ref_ref;
                    }
                }
                const int deletion_hap =
                    deletion.hap_to_cons_alle[1] == 1 ? 0 : 1;
                const int reference_hap = 1 - deletion_hap;
                // Independent physical events must agree with the MSA alleles
                // on both sides of the diploid contrast. Test their association
                // with the same Fisher threshold as the graph contrast; one
                // sequencing error must not veto a well-supported projection.
                const bool source_events_agree =
                    source_event_calls[0][0] > source_event_calls[0][1] &&
                    source_event_calls[1][1] > source_event_calls[1][0] &&
                    fisher_exact_two_tail(source_event_calls[0][0], source_event_calls[0][1],
                        source_event_calls[1][0], source_event_calls[1][1]) <= kOverlapAlleleMaxP;
                if (physical_alt != 0 || !source_events_agree ||
                    source_ref == 0 || source_del == 0 ||
                    graph_del_alt == 0 || graph_ref_ref == 0 ||
                    graph_ref_alt != 0 ||
                    fisher_exact_two_tail(graph_del_alt, graph_del_ref,
                                          graph_ref_alt, graph_ref_ref) >
                        kOverlapAlleleMaxP ||
                    ref_by_hap[reference_hap] <= ref_by_hap[deletion_hap] ||
                    del_by_hap[deletion_hap] <= del_by_hap[reference_hap] ||
                    fisher_exact_two_tail(
                        ref_by_hap[reference_hap],
                        ref_by_hap[deletion_hap],
                        del_by_hap[reference_hap],
                        del_by_hap[deletion_hap]) > kOverlapHapMaxP)
                    continue;
                CandidateVariant& snp = chunk.candidates[ci];
                snp.phase_set = mapped->second;
                snp.hap_to_cons_alle[1] = deletion.hap_to_cons_alle[1];
                snp.hap_to_cons_alle[2] = deletion.hap_to_cons_alle[2];
                if (old_graph_hap != snp.hap_to_cons_alle[1]) {
                    std::swap(snp.hap_to_alle_profile[1],
                              snp.hap_to_alle_profile[2]);
                    std::swap(snp.hap_alt, snp.hap_ref);
                }
                // A one-site graph block has no other gauge to preserve.
                // Its HP labels were based solely on the child walk, so exact
                // BAM events can replace them. In a larger graph block, keep
                // its read labels anchored by the other graph sites.
                if (graph_ps_site_count[old_graph_ps] == 1) {
                    for (size_t ri = 0; ri < source.reads.size(); ++ri) {
                        const int allele = physical_calls[ri];
                        if (allele != 0 && allele != 1) continue;
                        const auto parent = final_read_by_qname.find(
                            source.reads[ri].qname);
                        if (parent == final_read_by_qname.end() ||
                            chunk.phase_sets[parent->second] != old_graph_ps)
                            continue;
                        chunk.haps[parent->second] =
                            deletion.hap_to_cons_alle[1] == allele ? 1 : 2;
                        chunk.phase_sets[parent->second] = mapped->second;
                    }
                }
                projected_snps.insert(ci);
            }
        }
    }
    for (const auto& [qname, source_read] : overlay_phase) {
        if (ambiguous_overlay_reads.count(qname) != 0) continue;
        const auto mapped = phase_set_remap.find(
            SourcePhaseSet{source_read.source_id, source_read.phase_set});
        const auto read = final_read_by_qname.find(qname);
        if (mapped == phase_set_remap.end() || read == final_read_by_qname.end())
            continue;
        graph_chunk.recovery_source_reads.push_back(RecoverySourceRead{
            read->second, mapped->second, source_read.hap});
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

    // Screen physical coverage before the expensive whole-flank solve. The
    // singleton route retains its geometry; a downstream MSA deletion can
    // also validate multiple molecules through the checks above.
    for (size_t wi = 0; wi < windows.size(); ++wi) {
        const RecoverySeam& seam = windows[wi];
        const auto left_extent = phase_set_extents.find(seam.left_phase_set);
        const auto right_extent = phase_set_extents.find(seam.right_phase_set);
        if (left_extent == phase_set_extents.end() ||
            right_extent == phase_set_extents.end())
            continue;
        constexpr hts_pos_t kMinSingletonGap = 10000;
        if (seam.end - seam.beg < kMinSingletonGap) continue;
        const auto group = std::find_if(
            groups.begin(), groups.end(), [wi](const TargetedWindowGroup& g) {
                return group_owns_window(g, wi);
            });
        if (group == groups.end()) continue;
        const size_t source_id = static_cast<size_t>(group - groups.begin());
        const PhasingChunk& local = discovered[source_id];
        std::array<std::map<hts_pos_t, int>, 2> source_votes;
        for (const SourceSite& site : source_sites) {
            if (site.solve_id != source_id ||
                (!site.graph_clean_snp && !site.graph_clean_indel))
                continue;
            if (site.graph_phase_set == seam.left_phase_set)
                ++source_votes[0][site.phase_set];
            else if (site.graph_phase_set == seam.right_phase_set)
                ++source_votes[1][site.phase_set];
        }
        if (source_votes[0].empty() || source_votes[1].empty()) continue;
        std::array<hts_pos_t, 2> local_source_ps{};
        for (size_t side = 0; side < 2; ++side)
            local_source_ps[side] = std::max_element(
                source_votes[side].begin(), source_votes[side].end(),
                [](const auto& a, const auto& b) {
                    return a.second < b.second;
                })->first;
        const auto has_geometry = [&](bool allow_msa_deletion) {
            std::array<std::optional<hts_pos_t>, 2> local_boundary;
            for (const CandidateVariant& site : local.candidates) {
                if (!verified_bam_bridge_site(site) &&
                    !(allow_msa_deletion && site.phase_set == local_source_ps[1] &&
                      verified_bam_bridge_site(site, true)))
                    continue;
                const hts_pos_t pos = site.key.sort_pos();
                if (site.phase_set == local_source_ps[0] && pos < seam.end &&
                    (!local_boundary[0] || pos > *local_boundary[0]))
                    local_boundary[0] = pos;
                if (site.phase_set == local_source_ps[1] && pos > seam.beg &&
                    (!local_boundary[1] || pos < *local_boundary[1]))
                    local_boundary[1] = pos;
            }
            if (!local_boundary[0] || !local_boundary[1] ||
                *local_boundary[0] >= *local_boundary[1])
                return false;
            constexpr int kMinScreenMapq = 30;
            int spanning = 0;
            for (const ReadRecord& read : local.reads) {
                if (!read.is_skipped && read.mapq >= kMinScreenMapq &&
                    read.mapq != 255 && read.beg <= *local_boundary[0] &&
                    read.end >= *local_boundary[1])
                    ++spanning;
                if (spanning > 1) break;
            }
            if (spanning == 0) {
                // A read may end at an MSA deletion before the next graph SNP.
                // Screen that geometry before the full solve tests its cohort link.
                for (const CandidateVariant& site : local.candidates) {
                    const hts_pos_t pos = site.key.sort_pos();
                    if (site.key.type != VariantType::Deletion ||
                        site.counts.category != VariantCategory::NoisyCandHet ||
                        !site.msa_verified || !is_phase_set_anchor(site) ||
                        site.phase_set != local_source_ps[0] ||
                        pos <= *local_boundary[0] || pos >= *local_boundary[1])
                        continue;
                    int crossing = 0;
                    for (const ReadRecord& read : local.reads) {
                        if (!read.is_skipped && read.mapq >= kMinScreenMapq &&
                            read.mapq != 255 &&
                            read.beg <= *local_boundary[0] && read.end >= pos)
                            ++crossing;
                        if (crossing > 1) break;
                    }
                    if (crossing == 1) {
                        spanning = 1;
                        break;
                    }
                }
            }
            if (!allow_msa_deletion) return spanning == 1;
            const bool msa_deletion_boundary = std::any_of(
                local.candidates.begin(), local.candidates.end(),
                [&](const CandidateVariant& site) {
                    const hts_pos_t pos = site.key.sort_pos();
                    return site.phase_set == local_source_ps[1] &&
                           pos == *local_boundary[1] &&
                           site.key.type == VariantType::Deletion && site.msa_verified &&
                           site.counts.category == VariantCategory::NoisyCandHet;
                });
            return spanning > 1 && msa_deletion_boundary;
        };
        const bool original_geometry = has_geometry(false);
        const bool deletion_geometry = !original_geometry && has_geometry(true);
        if (!original_geometry && !deletion_geometry) continue;
        RegionChunk validation = group->region;
        validation.beg = std::max(chunk.ref_beg,
            std::min(left_extent->second.first, right_extent->second.first) -
                kSingletonValidationFlank);
        validation.end = std::min(chunk.ref_end,
            std::max(left_extent->second.second, right_extent->second.second) +
                kSingletonValidationFlank);
        Options validation_opts = targeted_solve_options(opts);
        validation_opts.retry_windows = {{seam.beg, seam.end}};
        if (!opts.phase_matrix_dump_prefix.empty())
            validation_opts.phase_matrix_dump_prefix =
                opts.phase_matrix_dump_prefix + ".validation.chunk" +
                std::to_string(chunk.region.chunk_id) + ".window" +
                std::to_string(wi);
        PhasingChunk full_source = process_chunk(validation, validation_opts, context);
        // Preserve the established certificate before extending its boundary.
        // A nearer recalled deletion must not displace a valid clean-site join.
        std::optional<ValidatedPhysicalBridge> bridge = original_geometry
            ? validated_physical_bridge(graph_chunk, seam, full_source, opts, false)
            : std::nullopt;
        if (!bridge &&
            (deletion_geometry || (original_geometry && has_geometry(true)))) {
            // Partial MSA reads can omit a boundary allele in this larger solve.
            // Restore targeted recovery's linkage calls, preserving discovery
            // counts, genotypes and source orientation, then validate physically.
            backfill_msa_observations(full_source, validation_opts, seam.beg, seam.end);
            bridge = validated_physical_bridge(
                graph_chunk, seam, full_source, opts, true);
        }
        // Keep the original certificate when it exists. A larger validation
        // solve can omit the verified insertion calls that connected the
        // targeted source to its first graph SNP; retry that observation mode
        // independently, then apply the same whole-block and physical checks.
        bool insertion_validation_bridge = false;
        if (!bridge && original_geometry) {
            const bool insertion_boundary = std::any_of(
                local.candidates.begin(), local.candidates.end(),
                [&](const CandidateVariant& site) {
                    return site.phase_set == local_source_ps[1] &&
                           site.key.type == VariantType::Insertion &&
                           site.counts.category == VariantCategory::NoisyCandHet &&
                           site.msa_verified && is_phase_set_anchor(site) &&
                           site.key.sort_pos() > seam.beg &&
                           site.key.sort_pos() < seam.end;
                });
            if (insertion_boundary) {
                Options insertion_opts = validation_opts;
                insertion_opts.add_unplaced_msa_observations = true;
                insertion_opts.recall_unplaced_msa_insertions = true;
                if (!insertion_opts.phase_matrix_dump_prefix.empty())
                    insertion_opts.phase_matrix_dump_prefix += ".insertion";
                PhasingChunk insertion_source = process_chunk(
                    validation, insertion_opts, context);
                bridge = validated_physical_bridge(
                    graph_chunk, seam, insertion_source, opts, false);
                insertion_validation_bridge = bridge.has_value();
            }
        }
        if (!bridge) continue;
        auto gauge = std::find_if(
            graph_chunk.recovery_phase_gauges.begin(),
            graph_chunk.recovery_phase_gauges.end(),
            [&](const RecoveryPhaseGauge& g) {
                return g.beg <= seam.beg && g.end >= seam.end;
            });
        if (gauge == graph_chunk.recovery_phase_gauges.end() ||
            std::any_of(gauge->physical_snp_bridges.begin(),
                        gauge->physical_snp_bridges.end(),
                        [&](const RecoveryPhysicalSnpBridge& existing) {
                            return existing.left_phase_set == seam.left_phase_set &&
                                   existing.right_phase_set == seam.right_phase_set;
                        }))
            continue;
        hts_pos_t pre_attach_source_ps = 0;
        if (bridge->pre_attach_source_ps > 0) {
            const auto mapped = phase_set_remap.find(SourcePhaseSet{
                source_id, bridge->pre_attach_source_ps});
            // The full-block validation solve can give this path a different
            // numeric PS from the targeted injection solve. Its physical
            // graph-flank bridge is still valid; only pre-attachment requires
            // an imported source block with the same PS in this solve.
            if (mapped != phase_set_remap.end()) {
                pre_attach_source_ps = mapped->second;
            } else if (!bridge->multiallelic_cohort_corroborated) {
                // One physical homopolymer allele is not enough to transfer a
                // whole graph-block orientation across distinct BAM solves.
                continue;
            }
        }
        if (insertion_validation_bridge) {
            // Use the normal stitcher's outer-allele and gauge conflict vetoes
            // on a disposable state; only its certified core relation is kept.
            PhasingChunk proof;
            proof.region = chunk.region;
            proof.ref_beg = chunk.ref_beg;
            proof.ref_end = chunk.ref_end;
            proof.ref_seq = chunk.ref_seq;
            proof.candidates = chunk.candidates;
            proof.reads.resize(chunk.reads.size());
            for (size_t ri = 0; ri < chunk.reads.size(); ++ri) {
                proof.reads[ri].beg = chunk.reads[ri].beg;
                proof.reads[ri].end = chunk.reads[ri].end;
                proof.reads[ri].mapq = chunk.reads[ri].mapq;
                proof.reads[ri].qname = chunk.reads[ri].qname;
                proof.reads[ri].is_skipped = chunk.reads[ri].is_skipped;
            }
            proof.read_var_profile = chunk.read_var_profile;
            proof.haps = chunk.haps;
            proof.phase_sets = chunk.phase_sets;
            RecoveryPhaseGauge proof_gauge = *gauge;
            proof_gauge.physical_snp_bridges.push_back(RecoveryPhysicalSnpBridge{
                seam.left_phase_set, seam.right_phase_set, bridge->flip, 0});
            Options proof_opts = opts;
            constexpr int kRecoveryLinkWindow = 128;
            proof_opts.block_link_window = kRecoveryLinkWindow;
            proof_opts.min_block_link_reads = 1;
            proof_opts.link_by_alleles = true;
            if (stitch_recovery_phase_sets_left_to_right(
                    proof, {seam}, {proof_gauge}, proof_opts,
                    &graph_chunk.recovery_source_path_supported,
                    &graph_chunk.recovery_source_weak_cuts,
                    &graph_chunk.recovery_source_quality_cuts) == 0)
                continue;
            std::optional<hts_pos_t> proof_ps;
            std::array<std::optional<bool>, 2> proof_flips;
            bool consistent = true;
            for (size_t ci = 0; ci < chunk.candidates.size(); ++ci) {
                const CandidateVariant& site = chunk.candidates[ci];
                if (!is_phase_set_anchor(site)) continue;
                const size_t side = site.phase_set == seam.left_phase_set ? 0 :
                    site.phase_set == seam.right_phase_set ? 1 : 2;
                if (side == 2) continue;
                const CandidateVariant& joined = proof.candidates[ci];
                const bool flip = site.hap_to_cons_alle[1] !=
                                  joined.hap_to_cons_alle[1];
                if ((proof_ps && *proof_ps != joined.phase_set) ||
                    (proof_flips[side] && *proof_flips[side] != flip))
                    consistent = false;
                proof_ps = joined.phase_set;
                proof_flips[side] = flip;
            }
            if (!consistent || !proof_flips[0] || !proof_flips[1] ||
                (*proof_flips[0] != *proof_flips[1]) != bridge->flip)
                continue;
            DeferredPhysicalBridge deferred;
            deferred.flip = bridge->flip;
            for (const CandidateVariant& site : chunk.candidates) {
                if (!is_phase_set_anchor(site)) continue;
                if (site.phase_set == seam.left_phase_set)
                    deferred.left_anchors.emplace_back(
                        site.key, site.hap_to_cons_alle[1]);
                if (site.phase_set == seam.right_phase_set)
                    deferred.right_anchors.emplace_back(
                        site.key, site.hap_to_cons_alle[1]);
            }
            graph_chunk.deferred_physical_bridges.push_back(std::move(deferred));
            continue;
        }
        gauge->physical_snp_bridges.push_back(RecoveryPhysicalSnpBridge{
            seam.left_phase_set, seam.right_phase_set, bridge->flip,
            pre_attach_source_ps});
    }

    if (opts.verbose > 0)
        fprintf(stderr, "[in-pass] %zu window(s) -> %zu region(s), merged %zu alignment site(s)\n",
                windows.size(), groups.size(), added);
    // Refreshed evidence alone is enough to run the final stitch: it can
    // supply the missing allele edge between two existing graph blocks.
    return true;
}

static void add_bam_observation_to_graph_profile(
        ReadVariantProfile& profile, size_t candidate_i, int allele, int alt_qi) {
    if (allele != 0 && allele != 1) return;
    if (profile.start_var_idx >= 0 &&
        candidate_i >= static_cast<size_t>(profile.start_var_idx)) {
        const size_t offset = candidate_i - static_cast<size_t>(profile.start_var_idx);
        // Targeted recovery already chose this BAM observation. A separate
        // whole-chunk solve may use different MSA clusters; its fallback call
        // must not replace the selected allele or borrow its base quality.
        if (offset < profile.bam_alleles.size() &&
            (profile.bam_alleles[offset] >= 0 ||
             profile.bam_alleles[offset] == kConflictingBamAllele))
            return;
    }

    const int old_start = profile.start_var_idx;
    const int old_end = profile.end_var_idx;
    const int new_start =
        old_start < 0 ? static_cast<int>(candidate_i)
                      : std::min(old_start, static_cast<int>(candidate_i));
    const int new_end =
        old_end < 0 ? static_cast<int>(candidate_i)
                    : std::max(old_end, static_cast<int>(candidate_i));
    const size_t new_size =
        static_cast<size_t>(new_end - new_start + 1);
    if (old_start < 0 || new_start != old_start || new_end != old_end) {
        const size_t shift =
            old_start < 0 ? 0 : static_cast<size_t>(old_start - new_start);
        const auto expand = [&](std::vector<int>& values, int fill) {
            std::vector<int> expanded(new_size, fill);
            if (old_start >= 0) {
                std::copy(values.begin(), values.end(),
                          expanded.begin() +
                              static_cast<std::ptrdiff_t>(shift));
            }
            values = std::move(expanded);
        };

        expand(profile.alleles, -1);
        expand(profile.alt_qi, 0);
        if (!profile.graph_alleles.empty()) expand(profile.graph_alleles, -1);
        if (!profile.bam_alleles.empty()) expand(profile.bam_alleles, -1);
        if (!profile.bam_qi.empty()) expand(profile.bam_qi, 0);
        if (!profile.bam_base_qualities.empty()) {
            std::vector<uint8_t> expanded(new_size, 0);
            if (old_start >= 0)
                std::copy(profile.bam_base_qualities.begin(),
                          profile.bam_base_qualities.end(),
                          expanded.begin() + static_cast<std::ptrdiff_t>(shift));
            profile.bam_base_qualities = std::move(expanded);
        }
        profile.start_var_idx = new_start;
        profile.end_var_idx = new_end;
    }
    const size_t offset =
        candidate_i - static_cast<size_t>(new_start);

    // The graph observation remains authoritative when both channels called
    // different alleles. BAM fills only an absent graph observation; retaining
    // the conflict in bam_alleles keeps it available for diagnostics without
    // turning disagreement into a phasing vote.
    if (profile.bam_alleles.empty())
        profile.bam_alleles.assign(new_size, -1);
    if (profile.bam_qi.empty())
        profile.bam_qi.assign(new_size, 0);
    profile.bam_alleles[offset] = allele;
    profile.bam_qi[offset] = alt_qi;
}

/// Add exact sequence-matched BAM alleles to existing graph candidates.
///
/// The whole-chunk BAM solve used for output fallback sees reads and private
/// indels that the GAF profile can omit. Candidate state and graph phasing stay
/// unchanged; this function only fills missing per-read observations so the
/// post-stitch statistical rescue can evaluate them in the final graph gauge.
static void attach_bam_observations_to_graph_profiles(
        GraphChunkBuildResult& graph_chunk,
        const PhasingChunk& bam,
        int solve_tid,
        const std::unordered_map<std::string_view, size_t>& graph_read_by_qname) {
    PhasingChunk& graph = graph_chunk.chunk;
    CandidateIndex raw_index;
    CandidateIndex sequence_index;
    for (size_t ci = 0; ci < graph.candidates.size(); ++ci) {
        const CandidateVariant& candidate = graph.candidates[ci];
        if (candidate.counts.n_uniq_alles != 2) continue;
        raw_index.emplace(cand_key_of(candidate), ci);
        const std::string* alt = selected_graph_candidate_alt(graph_chunk, ci);
        if (alt == nullptr) continue;
        const GraphSiteMeta& meta = graph_chunk.site_meta[ci];
        const VariantKey translated =
            vcf_to_variant_key(solve_tid, meta.pos, meta.ref, *alt);
        sequence_index.emplace(
            CandKey{translated.sort_pos(), static_cast<int>(translated.type),
                    translated.ref_len, translated.alt},
            ci);
    }

    const size_t missing_candidate = graph.candidates.size();
    std::vector<size_t> graph_candidate_by_bam_candidate(
        bam.candidates.size(), missing_candidate);
    for (size_t bam_candidate_i = 0;
         bam_candidate_i < bam.candidates.size(); ++bam_candidate_i) {
        const CandidateVariant& candidate = bam.candidates[bam_candidate_i];
        if (candidate.counts.n_uniq_alles > 2) continue;
        const ParentCandidateMatch parent = find_parent_candidate(
            raw_index, sequence_index, cand_key_of(candidate),
            missing_candidate);
        graph_candidate_by_bam_candidate[bam_candidate_i] = parent.index;
    }

    for (size_t bam_read_i = 0;
         bam_read_i < bam.reads.size() &&
         bam_read_i < bam.read_var_profile.size(); ++bam_read_i) {
        const auto graph_read =
            graph_read_by_qname.find(bam.reads[bam_read_i].qname);
        if (graph_read == graph_read_by_qname.end()) continue;
        const size_t graph_read_i = graph_read->second;
        if (graph_read_i >= graph.read_var_profile.size()) continue;

        const ReadVariantProfile& source = bam.read_var_profile[bam_read_i];
        if (source.start_var_idx < 0) continue;
        for (size_t offset = 0; offset < source.alleles.size(); ++offset) {
            const size_t bam_candidate_i =
                static_cast<size_t>(source.start_var_idx) + offset;
            if (bam_candidate_i >= bam.candidates.size()) break;
            const int allele = source.alleles[offset];
            if (allele != 0 && allele != 1) continue;
            const size_t graph_candidate_i =
                graph_candidate_by_bam_candidate[bam_candidate_i];
            if (graph_candidate_i == missing_candidate) continue;

            add_bam_observation_to_graph_profile(
                graph.read_var_profile[graph_read_i],
                graph_candidate_i, allele,
                offset < source.alt_qi.size() ? source.alt_qi[offset] : 0);
        }
    }
}

size_t recover_independent_bam_read_blocks_in_place(
        GraphChunkBuildResult& graph_chunk,
        const Options& opts,
        WorkerContext& context,
        const char* contig_name) {
    PhasingChunk& graph = graph_chunk.chunk;
    if (graph.reads.empty()) return 0;

    const int solve_tid = contig_name != nullptr
                              ? sam_hdr_name2tid(context.primary_header(), contig_name)
                              : graph.region.tid;
    if (solve_tid < 0) return 0;

    RegionChunk region = graph.region;
    region.tid = solve_tid;
    region.beg = graph.ref_beg;
    region.end = graph.ref_end;

    Options sub = targeted_solve_options(opts);
    sub.retry_windows.clear();
    if (!opts.phase_matrix_dump_prefix.empty())
        sub.phase_matrix_dump_prefix =
            opts.phase_matrix_dump_prefix + ".whole-bam";
    PhasingChunk bam = process_chunk(region, sub, context);

    std::unordered_map<std::string_view, size_t> graph_read_by_qname;
    graph_read_by_qname.reserve(graph.reads.size());
    for (size_t read_i = 0; read_i < graph.reads.size(); ++read_i)
        graph_read_by_qname.try_emplace(graph.reads[read_i].qname, read_i);

    // Preserve the whole-chunk BAM solve's exact allele observations for
    // sequence-identical graph candidates. They are consumed only by the
    // post-stitch read rescue and cannot alter graph candidate phase state.
    if (!opts.phase_matrix_dump_prefix.empty())
        dump_recovery_phase_state(graph, opts, "bam-overlay-input");
    attach_bam_observations_to_graph_profiles(
        graph_chunk, bam, solve_tid, graph_read_by_qname);
    if (!opts.phase_matrix_dump_prefix.empty())
        dump_recovery_phase_state(graph, opts, "bam-overlay-output");

    struct BamBlockEvidence {
        std::unordered_map<hts_pos_t, size_t> link_by_graph_phase_set;
        std::vector<IndependentBamBlockLink> links;
    };
    std::unordered_map<hts_pos_t, BamBlockEvidence> evidence_by_bam_phase_set;

    // Validate each BAM block against graph-assigned reads. Each graph phase
    // set has its own arbitrary HP orientation, so its 2x2 table remains
    // separate until the validator chooses that link's better orientation.
    for (size_t bam_read_i = 0; bam_read_i < bam.reads.size(); ++bam_read_i) {
        if (bam_read_i >= bam.haps.size() ||
            bam_read_i >= bam.phase_sets.size() ||
            (bam.haps[bam_read_i] != 1 && bam.haps[bam_read_i] != 2) ||
            bam.phase_sets[bam_read_i] <= 0) {
            continue;
        }
        const auto found = graph_read_by_qname.find(bam.reads[bam_read_i].qname);
        if (found == graph_read_by_qname.end()) continue;
        const size_t graph_read_i = found->second;
        if (graph_read_i >= graph.haps.size() ||
            graph_read_i >= graph.phase_sets.size() ||
            (graph.haps[graph_read_i] != 1 && graph.haps[graph_read_i] != 2) ||
            graph.phase_sets[graph_read_i] <= 0) {
            continue;
        }

        BamBlockEvidence& block =
            evidence_by_bam_phase_set[bam.phase_sets[bam_read_i]];
        const auto [link_it, inserted] = block.link_by_graph_phase_set.try_emplace(
            graph.phase_sets[graph_read_i], block.links.size());
        if (inserted) block.links.emplace_back();
        ++block.links[link_it->second]
              .counts[static_cast<size_t>(bam.haps[bam_read_i] - 1)]
                     [static_cast<size_t>(graph.haps[graph_read_i] - 1)];
    }

    std::unordered_set<hts_pos_t> supported_bam_phase_sets;
    supported_bam_phase_sets.reserve(evidence_by_bam_phase_set.size());
    for (const auto& entry : evidence_by_bam_phase_set) {
        if (independent_bam_block_is_supported(entry.second.links))
            supported_bam_phase_sets.insert(entry.first);
    }

    graph.haps.resize(graph.reads.size(), 0);
    graph.phase_sets.resize(graph.reads.size(), kUnphasedReadPhaseSet);
    graph.bam_fallback_haps.assign(graph.reads.size(), 0);
    graph.bam_fallback_phase_sets.assign(
        graph.reads.size(), kUnphasedReadPhaseSet);
    graph.bam_output_fallback_reads.clear();

    size_t recovered = 0;
    for (size_t bam_read_i = 0; bam_read_i < bam.reads.size(); ++bam_read_i) {
        if (bam_read_i >= bam.haps.size() ||
            bam_read_i >= bam.phase_sets.size() ||
            (bam.haps[bam_read_i] != 1 && bam.haps[bam_read_i] != 2) ||
            bam.phase_sets[bam_read_i] <= 0 ||
            (supported_bam_phase_sets.count(bam.phase_sets[bam_read_i]) == 0 &&
             bam.reads[bam_read_i].hap_score_margin <
                 kIndependentBamReadMinHapScoreMargin) ||
            bam.phase_sets[bam_read_i] >
                std::numeric_limits<int32_t>::max() - kBamFallbackPsOffset) {
            continue;
        }
        const hts_pos_t fallback_phase_set =
            bam.phase_sets[bam_read_i] + kBamFallbackPsOffset;
        const auto found = graph_read_by_qname.find(bam.reads[bam_read_i].qname);
        if (found == graph_read_by_qname.end()) {
            // A read without a catalog-site observation has no graph profile.
            // Preserve its independent BAM assignment for phased-BAM output;
            // it must not become graph evidence.
            graph.bam_output_fallback_reads.push_back(
                {bam.reads[bam_read_i].qname, bam.haps[bam_read_i],
                 fallback_phase_set});
            ++recovered;
            continue;
        }
        const size_t graph_read_i = found->second;
        const bool graph_primary =
            graph_read_i < graph.haps.size() &&
            graph_read_i < graph.phase_sets.size() &&
            (graph.haps[graph_read_i] == 1 || graph.haps[graph_read_i] == 2) &&
            graph.phase_sets[graph_read_i] > 0;
        if (graph_primary ||
            graph.bam_fallback_haps[graph_read_i] == 1 ||
            graph.bam_fallback_haps[graph_read_i] == 2) {
            continue;
        }

        graph.bam_fallback_haps[graph_read_i] = bam.haps[bam_read_i];
        graph.bam_fallback_phase_sets[graph_read_i] = fallback_phase_set;
        ++recovered;
    }
    return recovered;
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
