#ifndef PGPHASE_GRAPH_BAM_ADAPTER_HPP
#define PGPHASE_GRAPH_BAM_ADAPTER_HPP

// Bridge between the graph pipeline (snarl sites + allele walk matching) and
// the shared phasing pipeline (k-means clustering, chunk stitching, VCF
// output).  build_graph_chunk converts per-read allele observations from
// graph_query into a PhasingChunk with CandidateVariant + ReadVariantProfile,
// applying parent-snarl gating, multi-allelic→biallelic decomposition, and
// three-phase depth/AF filtering along the way.

#include "collect_types.hpp"
#include "allele_identity.hpp"
#include "graph_query.hpp"
#include "graph_sites.hpp"

#include <htslib/sam.h>

#include <array>
#include <iosfwd>
#include <optional>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <tuple>
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
    // Decomposed graph zero means "not selected ALT", not literal REF.
    bool non_selected_alt_class = false;
};

struct CandidateAlleleContrast {
    AlleleContrastKey key;
    // Canonical sequence order mapped to the original profile's allele indices.
    std::array<int, 2> local_alleles;
};

struct GraphChunkBuildResult {
    /// Reference-equivalent descriptions; original rows retain recovery/output
    /// contracts. Rebuild this index after changing candidate indices.
    std::vector<std::pair<CandidateIdentityKey, std::vector<size_t>>> joint_candidate_loci;
    std::vector<std::optional<CandidateAlleleContrast>> candidate_allele_contrasts;
    std::vector<std::pair<AlleleContrastKey, std::vector<size_t>>> joint_allele_contrasts;
    /// Union gap phasing: last-resort read labels applied after chunk stitching,
    /// as (read, anchor read of the same phase set, same haplotype as the anchor).
    std::vector<std::tuple<size_t, size_t, bool>> deferred_read_labels;
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

/// Group complete biallelic REF/ALT descriptions without combining their contracts.
void rebuild_joint_candidate_loci(GraphChunkBuildResult& graph_chunk,
    const std::function<char(hts_pos_t)>& reference_base);

/// Project compatible aliases to one molecule vote in the initial clean solve,
/// then restore original profiles/metadata and share the resulting locus phase.
void phase_joint_graph_candidates(GraphChunkBuildResult& graph_chunk, const Options& opts,
    const std::function<char(hts_pos_t)>& reference_base);

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

/// Haplotag reads that only observe sites excluded from the clean graph solve.
///
/// Already phased reads orient an excluded biallelic site within one existing
/// phase set. A directly phased locus, two independently inferred loci, or one
/// inferred locus with a statistically bounded primary-read error rate may
/// assign an unphased read to a separate read-only block. This pass never
/// changes candidate phasing or joins phase sets.
size_t rescue_unphased_graph_reads(PhasingChunk& chunk);

/// Find a shared pure insertion motif up to eight bases; the first allele may
/// be reference. Return zero for a compound or incompatible sequence pair.
size_t tandem_insertion_motif_length(const std::string& first, const std::string& second);

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

} // namespace pgphase_collect

#endif
