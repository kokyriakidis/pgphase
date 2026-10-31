#ifndef PGPHASE_GRAPH_BAM_ADAPTER_HPP
#define PGPHASE_GRAPH_BAM_ADAPTER_HPP

// Bridge between the graph pipeline (snarl sites + allele walk matching) and
// the shared phasing pipeline (k-means clustering, chunk stitching, VCF
// output).  build_graph_chunk converts per-read allele observations from
// graph_query into a PhasingChunk with CandidateVariant + ReadVariantProfile,
// applying parent-snarl gating, multi-allelic→biallelic decomposition, and
// three-phase depth/AF filtering along the way.

#include "refined_alignment.hpp"
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
    std::string trace;          // PGPHASE_TRACE only: each chunk's label and how it came about (reads.tsv)
    hts_pos_t beg = 0, end = 0;  // PGPHASE_TRACE only: the span of the sites the read was called at
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
    // Written no record; its alleles still serve the reads' calls (unlike an empty REF).
    bool no_record = false;
    // Written 1/1: both consensus sequences of the MSA region it lies in carry its allele (union_select.cpp).
    bool msa_hom = false;
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
    // --refine-aln: the reads' alignments the alignment solve's MSA rewrote (aligned output).
    std::vector<RefinedAlignment> refined_alignments;
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

/// A read's haplotype label in its chunk: the solve's (haps, phase_sets) where it phased the read, else a fill
/// label's (gap_haps, gap_phase_sets: rescue markers, the EM's island blocks) in its block's gauge, without
/// kGapFillPsOffset. hap 0 where the read has neither.
struct ChunkReadLabel {
    int hap = 0;
    hts_pos_t phase_set = kUnphasedReadPhaseSet;
    bool primary = false;  // the solve's own label
};
ChunkReadLabel read_label(const PhasingChunk& chunk, size_t ri);

/// A phased het alone in its phase set phases nothing on its own. It is
/// written in the phase set holding most of its labelled reads (at least
/// 3) when they split there by allele on both haplotypes (80% each),
/// with those reads' orientation; the phaser's own labels are untouched.
/// Such a site typically sits between two blocks that were not joined,
/// its reads labelled in either: its phase against one of them comes
/// from that block's reads alone. Over chr12/18/20 +2,599 phased truth
/// hets (+366 when 80% of the reads had to sit in one phase set),
/// switches, flips and reads unchanged.
/// Two het SNPs at one position with different ALTs, both left unphased (a 1|2 site), are placed the same way from
/// the labels of the reads carrying each ALT, on opposite haplotypes of one phase set. Neither changes a read label or
/// a join; the VCF writer applies them (graph_collect.cpp).
struct UnlinkedHetPlacement {
    hts_pos_t phase_set = 0;
    int hap1_allele = -1, hap2_allele = -1;
};
std::unordered_map<size_t, UnlinkedHetPlacement> place_unlinked_hets(
    const PhasingChunk& chunk, std::unordered_map<size_t, std::string>* trace_votes = nullptr);

/// The sequence of a row's allele `k` (its profile's allele index) over its site: the site's REF for the reference
/// walk, else that walk's ALT (site_allele_orig_idx). Catalog alleles carry the site's full flanks, so each replaces
/// all of REF. Null where the row has no allele k or no meta (a row the writer skips has an empty REF).
const std::string* row_allele_seq(const GraphChunkBuildResult& graph_chunk, size_t ci, size_t k);

/// A row's change as its site holds it: the site's REF to the row's allele 1. A row's allele 0 is its site's REF
/// (Phase 2 splits every site into REF/ALT rows; an injected row's alleles are 0, 1[, 2]).
struct RowAllele {
    hts_pos_t pos = 0;  // 1-based, of ref[0]
    const std::string* ref = nullptr;  // the site's whole REF: a 13-allele site's is 753 bp (chr20:23,847,472)
    const std::string* alt = nullptr;
};
std::optional<RowAllele> row_allele(const GraphChunkBuildResult& graph_chunk, size_t ci);

/// The row's change in minimal VCF form (upper case, trim_to_minimal_vcf): the reference the row's allele itself
/// replaces, [pos - 1, pos - 1 + ref.size()) 0-based. Of that 753 bp REF one row's allele replaces 702 bp, another's
/// one base.
struct VcfAllele {
    hts_pos_t pos = 0;  // 1-based
    std::string ref, alt;
};
std::optional<VcfAllele> row_vcf_allele(const GraphChunkBuildResult& graph_chunk, size_t ci);

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
    const GraphChunkBuildResult& gc);

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
