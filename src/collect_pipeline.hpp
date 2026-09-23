#ifndef PGPHASE_COLLECT_PIPELINE_HPP
#define PGPHASE_COLLECT_PIPELINE_HPP

/**
 * @file collect_pipeline.hpp
 * @brief Public declarations for region chunking and collect-bam-variation orchestration.
 */

#include "collect_types.hpp"
#include "graph_bam_adapter.hpp"

#include <string>
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

/// Recover bounded seams between neighboring phase sets inside a graph chunk.
///
/// The alignment caller supplies candidate rows, local candidate/read HP/PS,
/// and per-read observations. The merge keeps the local gauge under a
/// collision-free PS label and leaves established graph assignments unchanged.
///
/// Returns true when new candidates or refreshed graph-site evidence was merged.
bool recover_phase_set_seams_in_place(GraphChunkBuildResult& graph_chunk,
                                       const Options& opts,
                                       WorkerContext& context,
                                       const char* contig_name);

/// Phase graph-unassigned reads from an independent whole-chunk BAM solve.
///
/// Only fallback HP/PS vectors are updated. Graph candidates, observations,
/// primary assignments, and cross-chunk stitching inputs remain unchanged.
size_t recover_independent_bam_read_blocks_in_place(
    GraphChunkBuildResult& graph_chunk,
    const Options& opts,
    WorkerContext& context,
    const char* contig_name);

/// One candidate the recovery sub-solve found inside a recovery window, with
/// every decision the merge made about it.
///
/// The merge translates keys, matches them against the parent, appends the new
/// ones and synthesises per-site metadata -- and the writer skips any candidate
/// whose metadata is missing, silently. Measured on chr20:55,336,460, where the
/// sub-solve finds the two candidates carrying the window's only informative
/// signal, both are merged into the chunk as NoisyCandHet, and both are dropped
/// at emission because site_meta does not extend to their indices. This record
/// makes each step observable, so completeness is checked rather than probed.
struct RecoveredCandidate {
    hts_pos_t pos = 0;
    int type = 0;
    int ref_len = 0;
    std::string alt;
    int category = 0;
    hts_pos_t win_beg = 0;
    hts_pos_t win_end = 0;
    bool known_raw = false;          ///< parent holds this exact key
    bool known_translated = false;   ///< parent holds it via its VCF form
    bool inside_window = false;
    bool category_admitted = false;
    bool appended = false;           ///< reached the merged candidate table
    bool meta_built = false;
    /// The alignment path vouched for this candidate: clean, or noisy with the
    /// MSA's verification. Recorded so a demoted locus that came back usable is
    /// visible in the audit rather than only inferable from its new category.
    bool alignment_verified = false;         ///< per-site metadata synthesised for it
    size_t meta_alts = 0;
    std::string meta_ref;
};

/// Append one TSV row per recovered candidate. Thread-safe.
void write_recovery_audit(const std::string& path,
                          const std::vector<RecoveredCandidate>& rows);

} // namespace pgphase_collect

/**
 * @brief CLI entry for `collect-bam-variation` (argv without subcommand name).
 * @param argc Argument count.
 * @param argv Arguments after subcommand removal.
 * @return Exit code (0 success, 1 error).
 */
int collect_bam_variation(int argc, char* argv[]);

#endif
