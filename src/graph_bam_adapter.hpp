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
    bool bam_alt_deletion_no_ref = false;
    bool bam_low_fraction_snp = false;
    std::vector<std::array<int, 2>> physical_allele_strands = {};
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
    // Preserve the independent consensus edit when a catalog row keeps its
    // graph representation and therefore cannot retain the source MSA flag.
    std::optional<VariantKey> msa_key = std::nullopt;
};

struct RecoverySourceRead {
    size_t read_index = 0;
    hts_pos_t phase_set = 0;
    int hap = 0;
};

/// A partial rescue transfer must not split a cohort with greater local read
/// coverage than its destination core. Coordinates are half-open intervals.
bool core_dominates_rescue_coverage(const PhasingChunk& chunk, hts_pos_t phase_set);

/// Return the deletion-ALT haplotype when both observation channels call ALT,
/// opposite-haplotype MSA SNPs are masked by the deletion, and no phased call
/// contradicts that haplotype. Physical CIGAR and source support are separate.
int masked_snp_deletion_haplotype(const PhasingChunk& chunk,
                                const ReadVariantProfile& profile, size_t deletion_i);

struct DeferredPhysicalBridge {
    std::vector<std::pair<VariantKey, int>> left_anchors;
    std::vector<std::pair<VariantKey, int>> right_anchors;
    bool flip = false;
    // A calibrated shared deletion also certifies its retained rescue marker.
    std::optional<VariantKey> shared_deletion;
    std::optional<VariantKey> calibrated_insertion;
    bool calibrated_insertion_repeat = false;
    bool calibrated_insertion_source_prefix = false;
    // Complete uncut sources with disjoint physical deletion calibration.
    std::optional<VariantKey> calibrated_source_deletion;
    // Boundary certificates retain the downstream owning chunk until rescue ends.
    std::optional<size_t> right_chunk_index;
    // Independently called molecules, in the certificate's left-anchor gauge.
    std::unordered_map<std::string, int> calibrated_read_haps;
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
    /// Independently certified graph edges after physical switch repair.
    /// The boolean records the relative allele gauge, invariant under block flips.
    std::map<std::pair<std::string, std::string>, bool> recovery_physical_graph_edges;
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
    /// Whole-block physical certificates applied after output-only read rescue.
    std::vector<DeferredPhysicalBridge> deferred_physical_bridges;
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

/// True only when every neighboring source anchor is separated by a weak cut.
/// Co-located rows retain their source contrast; a source singleton has no
/// internal edge to classify. `weak_cuts` must be coordinate ordered.
bool bam_site_has_only_weak_links(const std::vector<CandidateVariant>& candidates,
                                  size_t site_index,
                                  const std::vector<hts_pos_t>& weak_cuts);

/// Retain a verified binary BAM genotype on an unphased catalog repeat.
/// The caller establishes exact allele identity and a safe independent PS label;
/// parallel graph metadata continues to describe the selected catalog allele.
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

/// Certify a shared recalled SNP with Q30 physical calls and two independent
/// clean source loci per molecule, on both haplotypes of a complete source.
bool retained_source_snp_path_anchor_supported(
    const GraphChunkBuildResult& graph_chunk, size_t candidate_index);

/// Recover the physical catalog deletion behind a shared, alignment-verified
/// MSA source row. A graph representation may retain repeat classification.
/// Read evidence must separately certify its current haplotype gauge.
std::optional<VariantKey> retained_shared_deletion_key(
    const GraphChunkBuildResult& graph_chunk, size_t candidate_index);

/// Check a staged BAM read against an independently calibrated clean SNP.
/// Physically deleted REF calls may abstain; every other phased call must agree.
int masked_bam_snp_haplotype(const PhasingChunk& chunk,
                            const ReadVariantProfile& profile, size_t witness_index,
                            const std::vector<size_t>& physically_masked);

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

/// Certify a source run from a BAM marker to its first exact shared graph SNP.
/// A complete graph path must separately certify the rest of the live block.
bool bam_source_prefix_to_graph_supported(const GraphChunkBuildResult& graph_chunk,
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

/// Require a diploid indel gauge and unanimous physical pairs. Bridge
/// error stays <=5%; calibration Wilson error, <=1% physical gauge error and
/// bridge error together stay <=20%.
std::optional<bool> calibrated_indel_bridge_flip(
    const IndependentBamBlockLink& gauge, const std::array<int, 2>& parity,
    double log_odds);

/// Require unanimous, diploid physical calibration and unanimous bridge calls.
/// The union bound across the actual calibration and bridge bases must be <=20%.
std::optional<bool> calibrated_source_deletion_bridge_flip(
    const IndependentBamBlockLink& gauge, const std::array<int, 2>& parity,
    double quality_error_bound);

/// Separate non-reference repeat lengths; zero length and equidistant calls abstain.
int complementary_insertion_length_class(int observed, int first, int second);

/// Require independently calibrated length classes and unanimous physical
/// bridges from both classes. Bound a coherent wrong union using both disjoint
/// class calibrations and their measured molecule errors.
std::optional<bool> calibrated_repeat_insertion_bridge_flip(
    const IndependentBamBlockLink& gauge,
    const std::array<std::array<int, 2>, 2>& parity,
    const std::array<double, 2>& call_errors);

/// Orient a physical repeat bridge after independent read gauges calibrate
/// both flanks. The joint calibration and molecule error must stay below 20%.
std::optional<bool> calibrated_repeat_snp_bridge_flip(
    const std::array<IndependentBamBlockLink, 2>& gauges,
    const std::array<int, 2>& parity, double wrong_parity_bound);

/// Join a physically validated repeat chain. The established left gauge must
/// pass its existing calibration; both allele classes must cross the interior
/// edge, and independent right SNP molecules must agree at least 80%.
std::optional<bool> calibrated_repeat_chain_flip(
    const IndependentBamBlockLink& gauge,
    const std::array<std::array<int, 2>, 2>& bridge,
    const std::array<int, 2>& right_parity, double right_log_odds);

/// Join two repeat contrasts calibrated on disjoint clean-SNP molecules.
std::optional<bool> calibrated_deletion_chain_flip(
        const std::array<IndependentBamBlockLink, 2>& gauges,
        const std::array<int, 2>& parity);

/// Orient a terminal insertion from two disjoint primary-molecule cohorts.
/// Both cohorts must independently support the same diploid SNP gauge.
/// Discordance plus call error is capped at 20% per cohort and 10% combined.
/// Confirm an already verified source insertion in independent diploid cohorts.
/// Both cohorts reject random association and agree on orientation; the pooled
/// Wilson discordance bound plus 1% physical call error cannot exceed 20%.
std::optional<int> calibrated_verified_insertion_hap1(
    const std::array<IndependentBamBlockLink, 2>& cohorts);

std::optional<int> calibrated_terminal_insertion_hap1(
    const std::array<IndependentBamBlockLink, 2>& cohorts);

/// Retire validated HOM ALT and REF-absent substitution/deletion anchors after the
/// initial solve, retaining the
/// gauge and read assignments of phase sets with surviving heterozygotes.
void reclassify_physically_validated_graph_snps(
    const GraphSiteCatalogView& catalog, GraphChunkBuildResult& graph_chunk);

/// Reject a graph heterozygote with significant absence of physical REF.
/// Homozygous ALT needs no deletion; mixed ALT/deletion retains the substantial
/// deletion requirement. Any REF or third base vetoes exclusion.
bool graph_snp_ref_absence_supported(int ref_count, int alt_count,
    int deletion_count, int other_count, size_t tested_sites);

/// A padded SNP may retain a significant substitution/deletion contrast with
/// two deletions. Its phase gauge must be calibrated separately before use.
bool graph_snp_padded_deletion_supported(int ref_count, int alt_count,
    int deletion_count, int other_count, size_t tested_sites);

/// Reject a repeat substitution whose physical ALT fraction fails the caller
/// floor and a familywise diploid balance test.
bool graph_snp_low_alt_fraction_supported(int ref_count, int alt_count,
    int deletion_count, int other_count, double min_af, size_t tested_sites);

/// Orient a deletion from unanimous independent, quality-weighted witnesses.
std::optional<int> physical_deletion_gauge_haplotype(
    const std::array<int, 2>& hap_counts, double wrong_gauge_bound);

struct GraphSnpReferenceEvidence {
    int reference_class_reads = 0;
    int alternate_class_reference_reads = 0;
    int alternate_class_deletions = 0;
    double wrong_alternate_bound = 1.0;
};

/// Exclude two SNPs observed on the same graph molecules from anchoring when
/// one is physically deleted and the other is reference-only. Missing evidence
/// (including any credible physical ALT) and weak contradiction abstain.
bool graph_snp_cohort_is_physically_contradicted(
    const std::optional<GraphSnpReferenceEvidence>& terminal,
    const std::optional<GraphSnpReferenceEvidence>& partner);

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
