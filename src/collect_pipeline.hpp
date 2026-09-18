#ifndef PGPHASE_COLLECT_PIPELINE_HPP
#define PGPHASE_COLLECT_PIPELINE_HPP

/**
 * @file collect_pipeline.hpp
 * @brief Public declarations for region chunking and collect-bam-variation orchestration.
 */

#include "collect_types.hpp"

#include <map>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

namespace pgphase_collect {

/**
 * @brief Public API for region chunking and parallel collect-bam-variation orchestration.
 *
 * Region partitioning into RegionChunk units, worker dispatch,
 * and merged candidate tables.
 */

/**
 * @brief Parses `chr`, `chr:pos`, or `chr:start-end` (commas allowed) into `RegionFilter`.
 * @param region Region literal from CLI.
 * @return Filter; empty string yields `enabled == false`.
 */
RegionFilter parse_region(const std::string& region);

/**
 * @brief Loads 3-column BED regions as inclusion filters (0-based BED → 1-based inclusive).
 * @param path BED file path.
 * @return List of filters; throws on I/O or parse errors.
 */
std::vector<RegionFilter> load_bed_regions(const std::string& path);

/**
 * @brief Tiles the genome (or user filters) into chunks of `opts.chunk_size` and annotates neighbors.
 * @param opts Chunk size and region inputs.
 * @param header BAM header.
 * @param fai Reference index (contig presence checked).
 * @return Annotated chunk list.
 */
std::vector<RegionChunk> build_region_chunks(const Options& opts,
                                             const bam_hdr_t* header,
                                             const faidx_t* fai);

// Overload that uses pre-built filters instead of re-parsing opts.regions.
// Used by the graph path to pass contig-resolved filters.
std::vector<RegionChunk> build_region_chunks(const Options& opts,
                                             const bam_hdr_t* header,
                                             const faidx_t* fai,
                                             const std::vector<RegionFilter>& filters);

/**
 * @brief Opens primary BAM + reference index and returns chunks from `build_region_chunks`.
 * @param opts Must set `primary_bam_file()`, `ref_fasta`, and region fields.
 * @return Chunk list; throws if BAM/index/FAI cannot be opened.
 */
std::vector<RegionChunk> load_region_chunks(const Options& opts);

/**
 * @brief Streaming driver: batch by `reg_chunk_i`, write TSV/VCF incrementally.
 * @param opts Output paths, reference, and BAM list.
 */
void run_collect_bam_variation(const Options& opts);

/**
 * @brief Hybrid BAM+graph driver: BAM variant calling augmented with graph
 *        snarl site observations before k-means phasing.
 *
 * Requires opts.graph_sites_vcf and opts.gaf_file (or gbz_db+gaf_db) in
 * addition to the standard BAM pipeline inputs.
 */
void run_collect_hybrid_variation(const Options& opts);


/// Recover what a solve could not phase from alignment evidence: both the
/// windows where reads were left unphased and the seams between blocks are
/// re-solved as their own chunks and stitched in on shared reads. See the
/// definition for `contig_name` and `allow_import`.
/// What the recovery's apply phase needs from a targeted sub-solve: read names
/// with their haplotype and phase set, and the phased candidates. Sequences and
/// alignments are deliberately dropped -- caching whole chunks for a
/// chromosome's worth of regions costs gigabytes and nothing downstream reads
/// them.
struct TargetedSolveResult {
    std::vector<std::string> qnames;
    std::vector<int> haps;
    std::vector<hts_pos_t> phase_sets;
    std::vector<CandidateVariant> candidates;
};

/// Region key (tid, beg, end) -> solved result.
using TargetedSolveCache =
    std::map<std::tuple<int, hts_pos_t, hts_pos_t>, TargetedSolveResult>;

/// Solve the recovery regions of EVERY listed chunk in one parallel batch.
///
/// Without this, recovery runs once per chunk and each call parallelises only
/// its own windows: measured on the graph arm over the first 10 Mb of chr20, 18
/// invocations with 1-4 merged regions each, a mean parallel width of 2.3 of 16
/// threads, run one after another. The regions are disjoint after merging, so
/// there is no ordering constraint between them. Results land in `cache`, and
/// the per-chunk entry point below finds its regions already solved; the regions
/// and the order they are applied in are unchanged, so this only widens the
/// parallelism.
void prewarm_targeted_solves(
        const std::vector<std::pair<PhasingChunk*, const char*>>& chunks,
        const Options& opts, WorkerContext& context, TargetedSolveCache& cache);

size_t recover_unphased_windows_from_bam(PhasingChunk& chunk, const Options& opts,
                                        WorkerContext& context,
                                        const char* contig_name = nullptr,
                                        bool allow_import = true,
                                        TargetedSolveCache* cache = nullptr);

} // namespace pgphase_collect

/**
 * @brief CLI entry for `collect-bam-variation` (argv without subcommand name).
 * @param argc Argument count.
 * @param argv Arguments after subcommand removal.
 * @return Exit code (0 success, 1 error).
 */
int collect_bam_variation(int argc, char* argv[]);

#endif
