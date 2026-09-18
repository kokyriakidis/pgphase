#ifndef PGPHASE_COLLECT_PIPELINE_HPP
#define PGPHASE_COLLECT_PIPELINE_HPP

/**
 * @file collect_pipeline.hpp
 * @brief Public declarations for region chunking and collect-bam-variation orchestration.
 */

#include "collect_types.hpp"

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


/// Recover what a solve could not phase, from the alignment channel: both the
/// windows where reads were left unphased and the seams between blocks are
/// re-solved as their own chunks at the recovery mapq floor and stitched in on
/// shared reads. Returns the number of parent blocks bridged. Requires a
/// WorkerContext, so the caller must have supplied a BAM.
/// `contig_name` names the chunk's contig so the BAM's own tid can be resolved;
/// pass nullptr when the caller's header IS the BAM's.
/// `allow_import` adds sites the caller's pass never discovered. It must be
/// FALSE for the graph pipeline: that path's output table is index-parallel to
/// per-site metadata (GraphChunkBuildResult::site_meta, site_ids,
/// site_allele_orig_idx), so appending candidates -- and the reorder that
/// follows -- decouples the arrays and silently mispairs or drops records.
/// Measured with it on: 33 candidates down to 24 and 23 phased heterozygotes
/// down to 14 on chr20:5,309,406. Adoption alone is index-safe: it mutates
/// candidates in place and adds none.
size_t recover_unphased_windows_from_bam(PhasingChunk& chunk, const Options& opts,
                                        WorkerContext& context,
                                        const char* contig_name = nullptr,
                                        bool allow_import = true);

} // namespace pgphase_collect

/**
 * @brief CLI entry for `collect-bam-variation` (argv without subcommand name).
 * @param argc Argument count.
 * @param argv Arguments after subcommand removal.
 * @return Exit code (0 success, 1 error).
 */
int collect_bam_variation(int argc, char* argv[]);

#endif
