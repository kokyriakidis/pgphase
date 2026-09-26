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

/// A recovery boundary must be an oriented heterozygote assigned to a real
/// phase set. Both BAM and graph candidates use 0 until phasing assigns an
/// anchor; reads separately use -1 while unphased. Testing `> 0` encodes the
/// candidate contract directly, and unequal haplotype alleles exclude
/// homozygous rows.
static bool is_phase_set_anchor(const CandidateVariant& cand) {
    if (cand.phase_set <= 0 || cand.hap_to_cons_alle.size() <= 2) return false;
    const int hap1 = cand.hap_to_cons_alle[1];
    const int hap2 = cand.hap_to_cons_alle[2];
    // Graph candidates can retain an original multiallelic index here. The
    // graph writer projects the selected biallelic row by testing allele 1, so
    // recovery must use that same emitted genotype when deciding whether a row
    // is an anchor. BAM candidates are already 0/1 and follow the same rule.
    return (hap1 == 1) != (hap2 == 1);
}

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
                return vcf_to_variant_key(
                    graph_chunk.chunk.region.tid, meta.pos, meta.ref,
                    meta.alts[static_cast<size_t>(alt_index)]).sort_pos();
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
static std::vector<RecoverySeam> collect_phase_set_seams(
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
};

// Exact binomial upper tail for local MSA retry admission.
static double binomial_upper_tail(int n, int first, double p) {
    const double log_term =
        std::lgamma(static_cast<double>(n + 1)) -
        std::lgamma(static_cast<double>(first + 1)) -
        std::lgamma(static_cast<double>(n - first + 1)) +
        static_cast<double>(first) * std::log(p) +
        static_cast<double>(n - first) * std::log1p(-p);
    double term = std::exp(log_term);
    double tail = term;
    for (int k = first; k < n; ++k) {
        term *= static_cast<double>(n - k) /
                static_cast<double>(k + 1) * p / (1.0 - p);
        tail += term;
    }
    return std::min(1.0, tail);
}

// A second MSA solve is admitted only for a two-block BAM boundary with
// missing, dropped-out, or conflicting indel calls on crossing molecules.
// Sparse conflicting pairs name one seam for a focused solve; other admission
// signals retain the existing grouped solve and its source-row guard.
static bool source_seam_needs_unplaced_msa(
        const PhasingChunk& source, const TargetedWindowGroup& group,
        const std::vector<RecoverySeam>& windows, const Options& opts,
        bool& preserve_source_rows, bool& ordinary_retry,
        std::optional<size_t>& sparse_window) {
    preserve_source_rows = false;
    ordinary_retry = false;
    sparse_window.reset();
    if (source.read_var_profile.size() != source.reads.size()) return false;
    constexpr int kAdmissionMinMapq = 30;
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
    constexpr double kDecisiveParityMaxP = 0.05;
    struct SourceRun {
        hts_pos_t phase_set;
        size_t first;
        size_t last;
    };
    for (size_t wi = group.first_window; wi < group.past_last_window; ++wi) {
        const RecoverySeam& seam = windows[wi];
        std::vector<SourceRun> runs;
        for (size_t ci = 0; ci < source.candidates.size(); ++ci) {
            const CandidateVariant& candidate = source.candidates[ci];
            const hts_pos_t pos = candidate.key.sort_pos();
            // Include both graph boundaries. A BAM insertion key can lie
            // one base beyond the right boundary's VCF anchor.
            if (pos < seam.beg || pos > seam.end + 1 ||
                candidate.phase_set <= 0 ||
                candidate.hap_to_cons_alle[1] < 0 ||
                candidate.hap_to_cons_alle[2] < 0 ||
                candidate.hap_to_cons_alle[1] ==
                    candidate.hap_to_cons_alle[2])
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
        const auto unique_boundary = [&source](hts_pos_t pos) {
            return std::count_if(source.candidates.begin(),
                                 source.candidates.end(),
                                 [pos](const CandidateVariant& candidate) {
                                     return candidate.key.sort_pos() == pos;
                                 }) == 1;
        };
        const bool left_unique = unique_boundary(left_pos);
        const bool right_unique = unique_boundary(right_pos);
        const bool unique_pair = left_unique && right_unique;
        const bool can_focus = group.focused_retry ||
            (group.past_last_window - group.first_window > 1 &&
             (wi == group.first_window ||
              wi + 1 == group.past_last_window));
        const auto complementary_boundary = [&](size_t selected) {
            const CandidateVariant& first = source.candidates[selected];
            const CandidateVariant* other = nullptr;
            for (size_t ci = 0; ci < source.candidates.size(); ++ci) {
                if (ci == selected ||
                    source.candidates[ci].key.sort_pos() !=
                        first.key.sort_pos())
                    continue;
                if (other != nullptr) return false;
                other = &source.candidates[ci];
            }
            return other != nullptr &&
                first.phase_set == other->phase_set &&
                first.key.type == other->key.type &&
                first.hap_to_cons_alle[1] >= 0 &&
                first.hap_to_cons_alle[1] <= 1 &&
                first.hap_to_cons_alle[2] >= 0 &&
                first.hap_to_cons_alle[2] <= 1 &&
                first.hap_to_cons_alle[1] != first.hap_to_cons_alle[2] &&
                first.hap_to_cons_alle[1] == other->hap_to_cons_alle[2] &&
                first.hap_to_cons_alle[2] == other->hap_to_cons_alle[1];
        };
        if (!unique_pair &&
            (!can_focus ||
             (!left_unique && !complementary_boundary(left_i)) ||
             (!right_unique && !complementary_boundary(right_i))))
            continue;
        int spanning = 0;
        int callable = 0;
        int left_alleles[2] = {0, 0};
        int right_alleles[2] = {0, 0};
        int same_parity = 0;
        int cross_parity = 0;
        for (size_t ri = 0; ri < source.reads.size(); ++ri) {
            const ReadRecord& read = source.reads[ri];
            if (read.is_skipped || read.mapq < kAdmissionMinMapq ||
                read.beg > left_pos || read.end < right_pos)
                continue;
            ++spanning;
            const ReadVariantProfile& profile = source.read_var_profile[ri];
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
        // A true diploid link has one allele parity. Retry MSA only when
        // opposite calls exceed the error model AND neither parity has a
        // decisive majority; an already supported link needs no new solve.
        const int opposing = std::min(same_parity, cross_parity);
        // Conflicting sparse indel calls have no statistically decisive
        // parity. Retry just this seam if a grouped MSA would disturb sites
        // elsewhere; unanimously oriented sparse calls do not need a retry.
        if ((left_indel || right_indel) && same_parity > 0 &&
            cross_parity > 0 &&
            paired_biallelic >= kSparseMinPairedCalls &&
            std::ldexp(1.0, -paired_biallelic) > kDecisivePairMaxP) {
            preserve_source_rows = true;
            if (!sparse_window) sparse_window = wi;
        }
        if (unique_pair && (left_indel || right_indel) &&
            paired_biallelic >=
                std::max(opts.min_depth, kMixedParityMinPairedCalls) &&
            opposing >= kMixedParityMinOpposingCalls &&
            binomial_upper_tail(paired_biallelic, opposing,
                                kMixedParityErrorRate) <= kMixedParityMaxP &&
            binomial_upper_tail(paired_biallelic,
                                std::max(same_parity, cross_parity), 0.5) >
                kDecisiveParityMaxP) {
            preserve_source_rows = true;
            ordinary_retry = true;
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

    const auto& window = *std::prev(after);
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

// Check the BAM solver's own phase path before translating its sites. A numeric
// phase-set label is insufficient: a single-haplotype bridge can hide a switch
// between two locally pure parts of that phase set.
// Read one physical base from the retained original alignment. The BAM solve
// has already freed its bulky parsed Digar vector, but keeps the bam1_t record.
// Returns 0 for reference, 1 for a deletion, 2 for SNP ALT, -1 otherwise.
static int physical_snp_call(const ReadRecord& read, hts_pos_t pos,
                             char ref_base, char alt_base,
                             int* base_quality = nullptr) {
    const bam1_t* aln = read.alignment.get();
    if (aln == nullptr || read.is_skipped) return -1;
    hts_pos_t ref_pos = aln->core.pos + 1;
    int query_pos = 0;
    const uint32_t* cigar = bam_get_cigar(aln);
    for (uint32_t i = 0; i < aln->core.n_cigar; ++i) {
        const int op = bam_cigar_op(cigar[i]);
        const int len = bam_cigar_oplen(cigar[i]);
        const int consumed = bam_cigar_type(op);
        if ((consumed & 2) != 0 && ref_pos <= pos && pos < ref_pos + len) {
            if (op == BAM_CDEL) return 1;
            if (op != BAM_CMATCH && op != BAM_CEQUAL && op != BAM_CDIFF)
                return -1;
            const int qi = query_pos + static_cast<int>(pos - ref_pos);
            if (qi < 0 || qi >= aln->core.l_qseq) return -1;
            const char base = seq_nt16_str[bam_seqi(bam_get_seq(aln), qi)];
            if (base_quality != nullptr)
                *base_quality = bam_get_qual(aln)[qi];
            if (base == ref_base) return 0;
            if (base == alt_base) return 2;
            return -1;
        }
        if ((consumed & 2) != 0) ref_pos += len;
        if ((consumed & 1) != 0) query_pos += len;
        if (ref_pos > pos) break;
    }
    return -1;
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
    constexpr int kMinBridgeBaseq = 30;
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

static SourcePathEvidence source_phase_set_path_evidence(
        const PhasingChunk& source, hts_pos_t source_ps) {
    std::vector<size_t> sites;
    for (size_t ci = 0; ci < source.candidates.size(); ++ci) {
        const CandidateVariant& cand = source.candidates[ci];
        if (cand.phase_set == source_ps &&
            cand.hap_to_cons_alle[1] >= 0 &&
            cand.hap_to_cons_alle[2] >= 0 &&
            cand.hap_to_cons_alle[1] != cand.hap_to_cons_alle[2])
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
    for (size_t ri = 0; ri < source.reads.size(); ++ri) {
        if (ri >= source.phase_sets.size() ||
            source.phase_sets[ri] != source_ps ||
            ri >= source.read_var_profile.size())
            continue;
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
                ++delta[previous][vote];
                --delta[current][vote];
            }
            previous = current;
            previous_hap = hap;
        }
    }
    constexpr double kSourceOneHapBridgeMaxP = 0.05;
    std::array<int, 3> crossing{};
    for (size_t i = 0; i + 1 < sites.size(); ++i) {
        for (size_t vote = 0; vote < crossing.size(); ++vote)
            crossing[vote] += delta[i][vote];
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
        const bool supported = consistent > crossing[2] &&
            (both_haps || single_hap_confident);
        if (left < right && !supported) {
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

// A REF SNP has no alt_qi in the longcallD profile. Recover its actual BAM
// base quality from the alignment before the targeted source is released.
static uint8_t bam_snp_base_quality(const bam1_t* bam, hts_pos_t pos) {
    if (bam == nullptr) return 0;
    hts_pos_t ref_pos = bam->core.pos + 1;
    int query_pos = 0;
    const uint32_t* cigar = bam_get_cigar(bam);
    for (uint32_t i = 0; i < bam->core.n_cigar; ++i) {
        const int length = bam_cigar_oplen(cigar[i]);
        const int consumption = bam_cigar_type(bam_cigar_op(cigar[i]));
        if ((consumption & 3) == 3 && pos >= ref_pos &&
            pos < ref_pos + length) {
            const uint8_t quality = bam_get_qual(bam)[
                query_pos + static_cast<int>(pos - ref_pos)];
            return quality == 255 ? 0 : quality;
        }
        if ((consumption & 1) != 0) query_pos += length;
        if ((consumption & 2) != 0) ref_pos += length;
    }
    return 0;
}

bool recover_phase_set_seams_in_place(GraphChunkBuildResult& graph_chunk,
                                      const Options& opts,
                                      WorkerContext& context,
                                      const char* contig_name) {
    PhasingChunk& chunk = graph_chunk.chunk;
    if (chunk.candidates.empty() || chunk.reads.empty()) return false;

    const int solve_tid = contig_name != nullptr
                              ? sam_hdr_name2tid(context.primary_header(), contig_name)
                              : chunk.region.tid;
    if (solve_tid < 0) return false;

    std::vector<RecoverySeam> windows = collect_phase_set_seams(graph_chunk);
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
    groups.reserve(initial_groups.size());
    discovered.reserve(initial_groups.size());
    const auto backfill = [&](PhasingChunk& target,
                              const TargetedWindowGroup& target_group,
                              const Options& target_opts) {
        for (size_t wi = target_group.first_window;
             wi < target_group.past_last_window; ++wi)
            backfill_msa_observations(target, target_opts,
                                      windows[wi].beg, windows[wi].end);
    };
    const auto preserves_hets = [](const PhasingChunk& before,
                                   const PhasingChunk& after) {
        std::set<CandKey> after_hets;
        for (const CandidateVariant& candidate : after.candidates)
            if (is_phase_set_anchor(candidate))
                after_hets.insert(cand_key_of(candidate));
        return std::all_of(
            before.candidates.begin(), before.candidates.end(),
            [&](const CandidateVariant& candidate) {
                return !is_phase_set_anchor(candidate) ||
                    after_hets.count(cand_key_of(candidate)) != 0;
            });
    };
    for (const TargetedWindowGroup& group : initial_groups) {
        Options source_opts = sub;
        PhasingChunk source = process_chunk(group.region, source_opts, context);
        // Exact-CIGAR backfill is part of the ordinary source evidence.
        backfill(source, group, source_opts);
        bool preserve_source_rows = false;
        bool ordinary_retry = false;
        std::optional<size_t> sparse_window;
        const bool needs_retry = source_seam_needs_unplaced_msa(
            source, group, windows, opts, preserve_source_rows,
            ordinary_retry, sparse_window);
        const bool isolated_edge = sparse_window && !ordinary_retry &&
            group.past_last_window - group.first_window > 1 &&
            (*sparse_window == group.first_window ||
             *sparse_window + 1 == group.past_last_window);
        if (isolated_edge) {
            const size_t wi = *sparse_window;
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
                if (isolated.region.beg < seam.beg &&
                    isolated.region.end > seam.end) {
                    Options isolated_opts = sub;
                    isolated_opts.retry_windows = {{seam.beg, seam.end}};
                    PhasingChunk local = process_chunk(
                        isolated.region, isolated_opts, context);
                    backfill(local, isolated, isolated_opts);
                    bool local_preserve = false;
                    bool local_ordinary = false;
                    std::optional<size_t> local_sparse;
                    if (source_seam_needs_unplaced_msa(
                            local, isolated, windows, opts, local_preserve,
                            local_ordinary, local_sparse)) {
                        isolated_opts.add_unplaced_msa_observations = true;
                        PhasingChunk local_retry = process_chunk(
                            isolated.region, isolated_opts, context);
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
                                complete_source_span = true;
                                break;
                            }
                        }
                        if (complete_source_span &&
                            preserves_hets(local, local_retry)) {
                            local = std::move(local_retry);
                            backfill(local, isolated, isolated_opts);
                            // The focused source owns this seam's observations
                            // and gauge. The original solve keeps its other
                            // seams, without repeating the broad MSA.
                            groups.push_back(isolated);
                            discovered.push_back(std::move(local));
                            TargetedWindowGroup remainder = group;
                            if (wi == group.first_window)
                                remainder.first_window = wi + 1;
                            else
                                remainder.past_last_window = wi;
                            groups.push_back(remainder);
                            discovered.push_back(std::move(source));
                            continue;
                        }
                    }
                }
            }
        }
        if (ordinary_retry ||
            (needs_retry && group.past_last_window - group.first_window == 1)) {
            source_opts.add_unplaced_msa_observations = true;
            PhasingChunk retried = process_chunk(group.region, source_opts, context);
            if (!preserve_source_rows || preserves_hets(source, retried)) {
                source = std::move(retried);
                backfill(source, group, source_opts);
            }
        }
        groups.push_back(group);
        discovered.push_back(std::move(source));
        dump_recovery_phase_state(discovered.back(), source_opts,
                                  "recovery-source");
    }

    // Match graph alleles by their selected reference sequence, not their
    // graph walk. This index is also used below for candidate transfer.
    CandidateIndex parent_seq_index;
    CandidateIndex parent_cand_index;
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
                    return value.first_window <= wi && wi < value.past_last_window;
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
    std::vector<std::map<SourceBlockGauge, std::pair<int, int>>>
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
    std::map<std::string, AlleleByCand> observed;
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
        bool can_adopt;
    };
    std::vector<SourceSite> source_sites;

    // Select complete BAM phase sets that contribute a phased site to a seam.
    std::vector<std::set<hts_pos_t>> selected_source_phase_sets(discovered.size());
    for (size_t gi = 0; gi < discovered.size(); ++gi) {
        for (const CandidateVariant& cand : discovered[gi].candidates) {
            if (cand.phase_set <= 0 ||
                cand.hap_to_cons_alle[1] < 0 ||
                cand.hap_to_cons_alle[2] < 0 ||
                cand.hap_to_cons_alle[1] == cand.hap_to_cons_alle[2])
                continue;
            if (find_containing_window(
                    windows, groups[gi], cand.key.sort_pos()) != nullptr)
                selected_source_phase_sets[gi].insert(cand.phase_set);
        }
    }

    for (size_t gi = 0; gi < discovered.size(); ++gi) {
        const PhasingChunk& src = discovered[gi];
        if (src.read_var_profile.size() != src.reads.size()) continue;
        const size_t source_id = gi;
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
            if (cand.phase_set > 0 &&
                cand.hap_to_cons_alle[1] >= 0 &&
                cand.hap_to_cons_alle[2] >= 0 &&
                cand.hap_to_cons_alle[1] != cand.hap_to_cons_alle[2]) {
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
                if (selected_source_phase_sets[gi].count(cand.phase_set) != 0 ||
                    graph_clean_snp)
                    source_sites.push_back(SourceSite{key, gi, ci, cand.phase_set,
                        cand.hap_to_cons_alle[1], cand.hap_to_cons_alle[2],
                        comparable ? matched->phase_set : 0,
                        comparable ? matched->hap_to_cons_alle[1] : -1,
                        graph_clean_snp, can_adopt});
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
                    if (graph_cand.hap_to_cons_alle[1] ==
                            cand.hap_to_cons_alle[1] &&
                        graph_cand.hap_to_cons_alle[2] ==
                            cand.hap_to_cons_alle[2]) {
                        ++vote.first;
                    } else if (graph_cand.hap_to_cons_alle[1] ==
                                   cand.hap_to_cons_alle[2] &&
                               graph_cand.hap_to_cons_alle[2] ==
                                   cand.hap_to_cons_alle[1]) {
                        ++vote.second;
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

                continue;
            }
            if (member == nullptr) continue;
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
            const ReadVariantProfile& prof = src.read_var_profile[ri];
            if (prof.start_var_idx < 0) continue;
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
                const CandKey observed_key = cand_key_of(src.candidates[ci]);
                per_read.emplace(observed_key,
                                 std::make_pair(prof.alleles[k],
                                                k < prof.alt_qi.size() ? prof.alt_qi[k] : 0));
                if (src.candidates[ci].key.type == VariantType::Snp) {
                    const uint8_t quality = bam_snp_base_quality(
                        src.reads[ri].alignment.get(), src.candidates[ci].key.sort_pos());
                    if (quality > 0)
                        observed_snp_qualities[src.reads[ri].qname].emplace(
                            observed_key, quality);
                }
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
        SourcePathEvidence path = source_phase_set_path_evidence(
            discovered[source.first], source.second);
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

    graph_chunk.recovery_phase_gauges.clear();
    graph_chunk.recovery_phase_gauges.reserve(groups.size());
    for (size_t gi = 0; gi < groups.size(); ++gi) {
        RecoveryPhaseGauge gauge;
        gauge.focused_retry = groups[gi].focused_retry;
        gauge.beg = groups[gi].region.beg;
        gauge.end = groups[gi].region.end;
        for (const auto& [source, mapped_phase_set] : phase_set_remap)
            if (source.first == gi)
                gauge.imported_phase_sets.push_back(mapped_phase_set);
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
                candidate_vote.shared_candidate_same = vote.first;
                candidate_vote.shared_candidate_cross = vote.second;
                gauge.block_votes.push_back(std::move(candidate_vote));
            } else {
                existing->shared_candidate_same += vote.first;
                existing->shared_candidate_cross += vote.second;
            }
        }
        constexpr hts_pos_t kPhysicalBridgeFlankContext = 2000;
        constexpr size_t kMaxPhysicalBridgeSitesPerFlank = 3;
        constexpr double kMaxPhysicalBridgeError = 0.01;
        for (size_t wi = groups[gi].first_window;
             wi < groups[gi].past_last_window; ++wi) {
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
            if (left_candidates.empty() || right_candidates.empty())
                continue;
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
                const size_t final_i = static_cast<size_t>(old_to_new[ci]);
                if (old_prof.alleles[k] >= 0)
                    alleles.emplace(final_i,
                                    std::make_pair(old_prof.alleles[k],
                                        k < old_prof.alt_qi.size() ? old_prof.alt_qi[k] : 0));
                const int graph_allele = old_prof.graph_alleles.empty()
                    ? old_prof.alleles[k]
                    : (k < old_prof.graph_alleles.size()
                        ? old_prof.graph_alleles[k] : -1);
                if (graph_allele >= 0)
                    graph_observations.emplace(final_i, graph_allele);
                if (k < old_prof.bam_alleles.size() &&
                    old_prof.bam_alleles[k] >= 0)
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
                alleles.insert_or_assign(idx->second, entry.second);
                bam_observations.insert_or_assign(idx->second, entry.second);
                if (read_qualities != observed_snp_qualities.end()) {
                    const auto quality = read_qualities->second.find(entry.first);
                    if (quality != read_qualities->second.end())
                        bam_snp_qualities.insert_or_assign(idx->second, quality->second);
                }
            }
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
            site.graph_phase_set, site.graph_hap1_allele, site.can_adopt});
    }

    // A shared graph row can occur in two padded BAM solves. Their HP gauges
    // are independent, so neither source may claim the row by iteration order.
    std::unordered_map<size_t, size_t> source_claims;
    for (const RecoverySourceSite& site : graph_chunk.recovery_source_sites)
        ++source_claims[site.candidate_index];
    for (RecoverySourceSite& site : graph_chunk.recovery_source_sites)
        if (source_claims[site.candidate_index] != 1)
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
                int physical_alt = 0, source_conflicts = 0;
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
                                if (allele != physical) ++source_conflicts;
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
                if (physical_alt != 0 || source_conflicts != 0 ||
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
    PhasingChunk bam = process_chunk(region, sub, context);

    std::unordered_map<std::string_view, size_t> graph_read_by_qname;
    graph_read_by_qname.reserve(graph.reads.size());
    for (size_t read_i = 0; read_i < graph.reads.size(); ++read_i)
        graph_read_by_qname.try_emplace(graph.reads[read_i].qname, read_i);

    // Preserve the whole-chunk BAM solve's exact allele observations for
    // sequence-identical graph candidates. They are consumed only by the
    // post-stitch read rescue and cannot alter graph candidate phase state.
    attach_bam_observations_to_graph_profiles(
        graph_chunk, bam, solve_tid, graph_read_by_qname);

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
