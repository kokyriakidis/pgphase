#ifndef PGPHASE_GRAPH_BAM_ADAPTER_HPP
#define PGPHASE_GRAPH_BAM_ADAPTER_HPP

// Bridge between the graph pipeline (snarl sites + allele walk matching) and
// the shared phasing pipeline (k-means clustering, chunk stitching, VCF
// output).  build_graph_chunk converts per-read allele observations from
// graph_query into a PhasingChunk with CandidateVariant + ReadVariantProfile,
// applying parent-snarl gating, multi-allelic→biallelic decomposition, and
// three-phase depth/AF filtering along the way.

#include "collect_types.hpp"
#include "graph_query.hpp"
#include "graph_sites.hpp"

#include <htslib/sam.h>

#include <array>
#include <iosfwd>
#include <optional>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

namespace pgphase_collect {

// Per-read phasing result accumulated across chunks during streaming output.
struct PhaseReadOutputRow {
    std::string read_name;
    int chunk_id = -1;
    int hap = 0;                // 1 or 2 when phased, 0 when unphased
    hts_pos_t phase_set = kUnphasedReadPhaseSet;  // positive when phased
    bool has_phased_assignment = false;
    bool has_primary_assignment = false;
    int copies = 0;             // number of chunks that observed this read
    std::unordered_map<std::string, int> allele_by_site;  // site_id → allele
};

// Record for a site that was dropped during depth/AF filtering (Phase 1-2).
struct FilteredGraphSite {
    std::string site_id;
    std::string chrom;          // output contig, for the --filtered-sites-out dump
    hts_pos_t pos = 0;          // VCF POS, so drops can be stratified by region
    int ref_cov = 0;
    int alt_cov = 0;
    int total_cov = 0;
    double allele_fraction = 0.0;
    std::string filter_reason;  // e.g. "ref_only", "low_depth", "high_af"
};

// VCF-level metadata for a snarl site, stored per candidate so the output
// converter doesn't need a global site lookup table.
struct GraphSiteMeta {
    std::string chrom;       // output contig (ref_contig or chrom)
    hts_pos_t pos = 0;
    std::string ref;
    std::vector<std::string> alts;
};

// Output of build_graph_chunk: a PhasingChunk ready for k-means phasing,
// plus graph-specific bookkeeping for VCF output and diagnostics.
struct RecoverySourceSite {
    size_t candidate_index = 0;
    hts_pos_t phase_set = 0;
    int hap1_allele = -1;
    int hap2_allele = -1;
    hts_pos_t graph_phase_set = 0;
    int graph_hap1_allele = -1;
    // Exact clean SNP in both the graph and this independent BAM solve.
    bool clean_shared_snp = false;
    bool can_adopt = false;
};

struct RecoverySourceRead {
    size_t read_index = 0;
    hts_pos_t phase_set = 0;
    int hap = 0;
};

struct GraphChunkBuildResult {
    /// Complete selected BAM phase-set membership, including shared graph rows.
    std::vector<RecoverySourceSite> recovery_source_sites;
    std::vector<RecoverySourceRead> recovery_source_reads;
    std::unordered_map<hts_pos_t, bool> recovery_source_path_supported;
    /// Coordinate cuts where the original BAM solve lacks a two-haplotype path.
    std::unordered_map<hts_pos_t, std::vector<hts_pos_t>>
        recovery_source_weak_cuts;
    /// Weak cuts with independent high-quality physical SNP calls on both sides.
    /// Only the seam containing such a cut may use it to validate a BAM bridge.
    std::unordered_map<hts_pos_t, std::vector<hts_pos_t>>
        recovery_source_quality_cuts;
    /// Bounded graph seams solved by BAM recovery. Each entry retains the
    /// canonical boundaries and exact adjacent graph PS identities for the
    /// final left-to-right stitch.
    std::vector<RecoverySeam> recovery_windows;
    /// Read-assignment gauge for each targeted BAM solve. The final stitch uses
    /// it when separate allele rows do not share a callable observation.
    std::vector<RecoveryPhaseGauge> recovery_phase_gauges;
    /// Core insertion joins certified on the final recovery matrix. Candidate
    /// indices stay stable through chunk stitching; apply after read rescue so
    /// output-only groups retain their independently assigned HP/PS gauges.
    std::vector<std::pair<size_t, size_t>> equivalent_insertion_joins;
    PhasingChunk chunk;
    // Snarl site ID per candidate (parallel to chunk.candidates).
    std::vector<std::string> site_ids;
    // VCF-level metadata per candidate (parallel to chunk.candidates).
    std::vector<GraphSiteMeta> site_meta;
    // Maps final biallelic allele indices back to original allele-walk indices
    // in GraphSite (0 = ref walk, 1 = first alt walk, etc.).
    std::vector<std::vector<int>> site_allele_orig_idx;
    // Sites that survived Phase 1 but were dropped in Phase 2.
    std::vector<FilteredGraphSite> filtered_sites;
};

/// Retain a verified binary BAM genotype on an unphased padded catalog repeat.
/// The caller must establish exact allele identity, suffix padding and an unused
/// source phase-set label; parallel graph metadata continues to describe it.
bool adopt_unphased_graph_allele_from_bam(CandidateVariant& graph,
                                         const CandidateVariant& source,
                                         hts_pos_t phase_set);

// Convert graph-space allele observations into a PhasingChunk for phasing.
// Applies parent-snarl gating, multi-allelic→biallelic decomposition,
// and three-phase depth/AF filtering.
/// Rebuild the read<->variant interval tree from chunk.read_var_profile. The
/// tree is keyed by CANDIDATE INDEX, so anything that inserts, removes or
/// reorders candidates must call this before the chunk is solved again.
void rebuild_read_var_cr(PhasingChunk& chunk);

/// Return the VCF ALT represented by a biallelic graph candidate.
/// Split candidates retain the full source ALT list; allele index 1 selects
/// exactly one original ALT. Whole multiallelic candidates have no binary match.
const std::string* selected_graph_candidate_alt(
    const GraphChunkBuildResult& graph_chunk, size_t candidate_index);

/// Fill missing phased SNP observations from unique, identical catalog branches.
/// Keeps the original candidate counts, genotypes, PS labels and read gauges;
/// conflicting source observations, repeated branches and gated children abstain.
/// Supplemental depth must preserve heterozygosity under the original AF limits.
/// Returns the number of observations added and rebuilds the profile index once.
size_t supplement_phased_snp_branches(const GraphSiteCatalogView& catalog,
    const std::vector<GraphReadAllele>& rows, GraphChunkBuildResult& graph_chunk,
    const Options& opts);

/// Read one physical BAM base at a 1-based SNP coordinate. Returns 0 for REF,
/// 2 for ALT, 1 for a deletion, and -1 when the base is not callable.
int physical_snp_call(const bam1_t* alignment, hts_pos_t pos,
                      char ref_base, char alt_base,
                      int* base_quality = nullptr);

/// Quality of an existing SNP observation only when the original aligned BAM
/// base calls that same REF/ALT allele. Missing quality, a third base, an indel
/// or a contradictory MSA-derived allele supplies no physical certificate.
uint8_t bam_snp_observation_quality(const bam1_t* alignment, hts_pos_t pos,
                                  char ref_base, char alt_base, int allele);

/// Admit an overlapping second recovery solve only with callable nearest SNPs.
/// Selected graph alleles are normalized to reference bases. A projected graph
/// anchor needs independent paired molecules supporting both haplotypes and
/// rejecting random parity at p <= 0.001, with consistent quality odds. BAM-only
/// pairs retain the established source admission rule. If requested, graph_parity
/// retains the admitted graph relation for stitching and is otherwise cleared.
/// This does not stitch blocks.
bool has_direct_snp_parity_for_retry(
    const GraphChunkBuildResult& graph_chunk, const RecoverySeam& seam,
    WorkerContext& context, int tid,
    std::optional<bool>* graph_parity = nullptr);

/// Certify a candidate pair's parity with a one-sided binomial test and
/// agreeing read halves. Pure BAM pairs use BAM observations and BAM MAPQ;
/// other pairs retain the working matrix and graph read eligibility.
/// Returns whether the target's current haplotype orientation must flip.
std::optional<bool> local_run_boundary_flip(
    const PhasingChunk& chunk, size_t source_index, size_t target_index,
    int min_mapq = 0, double max_p = 0.05);

/// Find the BAM prefix of a complementary deletion pair after its rightward
/// bridge is independently certified. A nearest clean SNP must have unanimous
/// exclusive-ALT pairs on both allele classes. Every moved anchor must share
/// that source gauge and lie in the same audited weak-cut-free run. Catalog
/// anchors up to graph_end stay independent; returns the first BAM coordinate.
std::optional<hts_pos_t> bam_prefix_before_deletion_pair(
    const GraphChunkBuildResult& graph_chunk, size_t first_index,
    size_t second_index, hts_pos_t graph_end);

/// Certify a detached BAM run's source gauge and weak-cut-free extent.
/// Oriented heterozygotes must belong to one source with a consistent flip;
/// homozygous rows contribute neither a gauge nor a run boundary. When an
/// independently called bridge admits shared graph rows, those rows need the
/// same exact source provenance. Only cuts inside the current extent matter.
bool bam_source_run_supported(const GraphChunkBuildResult& graph_chunk,
                              hts_pos_t phase_set,
                              bool include_shared_graph = false);

/// Detach private BAM islands attached across an audited weak source cut.
/// Eligible owners in the preceding source component are checked independently.
/// Catalog anchors and observed read bridges preserve additional owners' joins.
/// Detached sites and exclusive source reads recover
/// their original BAM gauge and an independent phase-set label.
void detach_bam_sites_across_weak_cuts(GraphChunkBuildResult& graph_chunk);

/// Certify an imported site's original BAM path and its current block gauge.
/// A relabeled site cannot borrow the destination block's source certificate;
/// another oriented site from its own cut-free source must confirm the gauge.
bool bam_source_site_path_supported(const GraphChunkBuildResult& graph_chunk,
                                    size_t candidate_index);

GraphChunkBuildResult build_graph_chunk(const GraphSiteCatalogView& catalog,
                                               const std::vector<GraphReadAllele>& rows,
                                               const std::string& contig,
                                               hts_pos_t beg,
                                               hts_pos_t end,
                                               int chunk_id,
                                               const Options& opts);

// Reclassify indel candidates that sit in homopolymer, tandem-repeat, or
// low-complexity reference contexts as RepeatHetIndel.  Mirrors the BAM
// pipeline's classify_variant_initial noise gate.
//
// @param result  Build result whose chunk.candidates are updated in place.
// @param ref_seq Reference sequence slice covering the chunk region.
// @param ref_beg 1-based start of ref_seq on the chromosome.
// @param ref_end 1-based end of ref_seq (inclusive).
// @param max_xgaps Maximum indel span to check (opts.noisy_reg_max_xgaps).
/// Re-admit repeat-context het indels that agree with a trusted neighbour.
///
/// `apply_graph_noise_filter` demotes every het indel in a homopolymer/STR
/// context, judging the locus by its reference sequence and never by whether
/// the reads separate there. This pass gives such a site one way back: take the
/// reads that observe both it and a nearby clean het SNP, and admit it when
/// they agree. Returns the number of sites promoted.
size_t promote_link_supported_repeat_indels(GraphChunkBuildResult& result,
                                            const Options& opts);

void apply_graph_noise_filter(GraphChunkBuildResult& result,
                              const std::string& ref_seq,
                              hts_pos_t ref_beg,
                              hts_pos_t ref_end,
                              int max_xgaps);

// Find reads shared between adjacent chunks (by name merge-intersect)
// and record them in PhasingChunk overlap bookkeeping for stitching.
void populate_graph_chunk_overlaps(std::vector<GraphChunkBuildResult>& graph_chunks);

// K-means hap assignment for every chunk, then overlap detection and stitching.
void phase_graph_chunks(std::vector<GraphChunkBuildResult>& graph_chunks,
                            const Options& opts);

/// Haplotag reads that only observe sites excluded from the clean graph solve.
///
/// Already phased reads orient an excluded biallelic site within one existing
/// phase set. A directly phased locus, two independently inferred loci, or one
/// inferred locus with a statistically bounded primary-read error rate may
/// assign an unphased read to a separate read-only block. This pass never
/// changes candidate phasing or joins phase sets.
size_t rescue_unphased_graph_reads(
    PhasingChunk& chunk, const std::vector<RecoverySeam>& recovery_windows = {});

/// Refresh inherited read HP using agreeing Q30, MAPQ30 physical SNP calls in
/// its final PS. Clean graph SNPs need two distinct loci at least 100 bases apart;
/// an imported MSA SNP may certify alone after physical links to clean SNPs
/// validate its gauge. Contradictions abstain; candidate gauges and PS stay fixed.
size_t refresh_recovered_read_haps_from_bam_snps(PhasingChunk& chunk);

/// Read links between one independently solved BAM block and one graph block.
/// Rows are BAM haplotypes and columns are graph haplotypes.
struct IndependentBamBlockLink {
    std::array<std::array<int, 2>, 2> counts{};
};

/// Return true when every informative graph link supports a clean BAM block.
///
/// Each link must contain both haplotypes and reject random association. The
/// combined discordance rate must also have a 95% Wilson upper bound at or
/// below 10%, so a high-depth but weak association cannot pass merely because
/// its p-value is small.
bool independent_bam_block_is_supported(
    const std::vector<IndependentBamBlockLink>& links);

/// Apply statistically validated BAM assignments only to reads that remain
/// unassigned after graph stitching and excluded-site rescue.
size_t apply_independent_bam_read_blocks(PhasingChunk& chunk);

// Fold one stitched chunk's per-read hap/PS assignments into a running map.
void merge_graph_chunk_into_read_rows(
    std::unordered_map<std::string, PhaseReadOutputRow>& rows_by_read,
    const GraphChunkBuildResult& gc,
    int min_read_hap_margin = 0);

// Write phased-BAM records for reads whose chunks are fully stitched.
// Reads still needed by the next chunk (next_chunk_qnames) are held back;
// pass nullptr to flush all remaining reads (end of contig).
void flush_graph_phase_bam_after_merge(
    samFile* phase_bam_out,
    sam_hdr_t* phase_bam_hdr,
    std::unordered_map<std::string, PhaseReadOutputRow>& rows_by_read,
    const std::unordered_set<std::string>* next_chunk_qnames,
    std::unordered_set<std::string>& emitted_read_names);

// Phase-sites TSV output (used by tests).
void write_graph_phase_sites_tsv_header(std::ostream& out);
void write_graph_phase_sites_tsv_rows(std::ostream& out,
                                          const GraphChunkBuildResult& gc);
void write_graph_phase_sites_tsv(std::ostream& out,
                                     const std::vector<GraphChunkBuildResult>& graph_chunks);

} // namespace pgphase_collect

#endif
