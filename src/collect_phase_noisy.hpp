#ifndef PGPHASE_COLLECT_PHASE_NOISY_HPP
#define PGPHASE_COLLECT_PHASE_NOISY_HPP

/**
 * @file collect_phase_noisy.hpp
 * @brief Step 4: iterative noisy-region MSA variant calling and re-phasing.
 *
 * Step 4: noisy-region MSA variant recall and re-phasing.
 *   sort_noisy_regs → collect_noisy_vars1 loop → assign_hap k-means (kCandGermlineVarCate).
 *
 * Alignment-heavy functions (WFA2-lib / abPOA) live in align.hpp / align.cpp.
 * This file covers the outer loop, MSA variant extraction, and merge into the chunk.
 */

#include "align.hpp"
#include "collect_var.hpp"

#include <array>
#include <cstdint>
#include <string_view>
#include <utility>
#include <vector>

namespace pgphase_collect {

/// Inclusive insertion positions producing the same reference edit, with the
/// inserted motif rotated at each step. Reference boundaries and ambiguous
/// bases stop the interval; no alignment or distance cutoff is involved.
std::pair<hts_pos_t, hts_pos_t> insertion_equivalent_positions(
    hts_pos_t pos, const std::string& alt, const PhasingChunk& chunk);
std::pair<hts_pos_t, hts_pos_t> insertion_equivalent_positions(
    hts_pos_t pos, const std::string& alt, ReferenceCache& reference,
    int tid, const bam_hdr_t* header);

/// Compare two pure insertion edits on the same supplied reference background.
/// Positions are 1-based insertion coordinates; reference_between covers
/// [left_pos, right_pos). The caller may supply verified common SNP alleles.
/// Ambiguous bases, reversed positions and incomplete background abstain.
bool insertion_edits_are_equivalent(
    hts_pos_t left_pos, std::string_view left_alt,
    hts_pos_t right_pos, std::string_view right_alt,
    std::string_view reference_between);

/// Query index of a shifted multi-base insertion ALT whose reference edit and
/// full crossed path match the candidate with known base qualities. Return -1
/// for exact placement, missing, different or compound events. This only calls
/// an ALT; it never changes source observations, assignments or candidate rows.
int bam_shifted_repeat_insertion_query_index(
    const bam1_t* read, const CandidateVariant& insertion,
    const PhasingChunk& chunk, int min_baseq);
int bam_shifted_repeat_insertion_query_index(
    const bam1_t* read, const CandidateVariant& insertion,
    ReferenceCache& reference, int tid, const bam_hdr_t* header, int min_baseq);

/// Call a site only when both consensus alignment paths agree with exact local flanks.
int call_msa_site_allele(const std::array<AlnStr, 2>& alignments,
                         const VariantKey& key, hts_pos_t ref_beg,
                         const std::array<AlnStr, 2>* consensuses = nullptr);


/// Extend MSA het profiles only if the expanded observations pass the existing AF gate.
void add_msa_site_observations(const Options& opts,
                                const std::vector<UnassignedMsaRead>& reads,
                                hts_pos_t ref_beg,
                                std::vector<CandidateVariant>& vars,
                                std::vector<ReadVariantProfile>& profiles,
                                const std::array<AlnStr, 2>* consensuses = nullptr);

/// Remove consensus-derived calls outside each physical BAM alignment and
/// recount depths in O(observations + sites). Bounds are one-based, inclusive;
/// indels require both surviving flanks. Profiles must align with reads.
/// Used by unplaced-read recovery; ordinary upstream BAM MSA is unchanged.
void restrict_msa_observations_to_read_coverage(
    const std::vector<ReadRecord>& reads,
    std::vector<CandidateVariant>& vars,
    std::vector<ReadVariantProfile>& profiles);

/// Fill the strand tallies of candidates whose counts the MSA path built,
/// derived in one sweep from the read profiles. update_variant_depth_fields
/// derives ref_cov/alt_cov from alle_covs and leaves the strand fields at zero,
/// which exempts those candidates from the strand-bias screen. Only a record
/// whose derivation reproduces its own ref_cov/alt_cov is filled.
void derive_msa_candidate_strand_counts(PhasingChunk& chunk);

/// nt4 code (A=0, C=1, G=2, T/U=3) for an ASCII base, case-insensitive; 4 for
/// anything else, including 'N'. Declared here so the predicate tests can
/// exercise it directly.
uint8_t base_to_nt4(char base);

/// Is this indel in a homopolymer context? The BAM MSA path can use
/// longcallD's raw-reference versus nt4 insertion comparison; other callers
/// compare bases in nt4. `alt` is ignored for deletions.
bool var_is_homopolymer_indel(const PhasingChunk& chunk,
                              hts_pos_t ref_pos,
                              VariantType type,
                              int ref_len,
                              const std::string& alt,
                              bool upstream_reference_bytes = false);

/// Call an exact BAM indel CIGAR allele with quality-checked flanks.
/// Returns 0 for REF, 1 for ALT, and -1 for an ambiguous alignment.
int bam_exact_indel_allele(const bam1_t* bam, const CandidateVariant& var,
                           int min_bq, int* alt_qi);

/// Return 1 for an exact or sequence-equivalent deletion, 0 for verified REF,
/// and -1 for missing, low-quality, compound, or different-allele observations.
/// Candidate coordinates and allele rows are preserved. Reference belongs to
/// the calling worker; this check performs no alignment.
int bam_equivalent_deletion_allele(const bam1_t* read,
                                    const CandidateVariant& deletion,
                                    ReferenceCache& reference, int tid,
                                    const bam_hdr_t* header, int min_baseq);

/// Certify only an exact deletion ALT sequence between surviving query anchors.
/// Shifted CIGAR deletions with compensating mismatches may match; false means
/// unverified, never REF. No read labels, candidate rows or alignments change.
bool bam_matches_deletion_sequence(const bam1_t* read,
                                    const CandidateVariant& deletion,
                                    ReferenceCache& reference, int tid,
                                    const bam_hdr_t* header, int min_baseq);

/// Recover missing ALT calls for shifted, single-base MSA insertions in targeted
/// recovery windows, inclusive of their VCF anchors. Require MAPQ30 and Q30
/// sequence-equivalent CIGAR evidence;
/// Preserve existing calls and complex alleles. Reindex profiles before phasing.
/// Returns the number of added calls; ordinary BAM runs with no windows are inert.
int backfill_shifted_msa_insertions(PhasingChunk& chunk, const Options& opts);

/// Update the primary allele, growing the sparse range and retaining the site
/// offsets of populated provenance channels. New sites have no BAM base quality.
void update_read_var_profile_with_allele(int var_idx, int allele, int alt_qi,
                                         ReadVariantProfile& profile);

/// Temporarily fill missing MSA deletion calls for source retry admission.
/// `include_isolated` is reserved for singleton graph flanks; the resulting
/// retry must preserve clean-SNP gauges. Homopolymer calls require MAPQ30/Q30.
/// Complementary co-located rows require a jointly missing pair and one verified
/// ALT; the other row records ALT absence. Third alleles, compound events,
/// ambiguous loci and existing observations are not projected.
/// Callers restore profiles and the read index before validation or transfer:
/// physical calls diagnose a missing/reversed source edge, never tag reads.
int backfill_msa_retry_deletions(PhasingChunk& chunk, const Options& opts,
                                 hts_pos_t beg, hts_pos_t end, bool include_isolated);

/// Fill missing observations at admitted MSA sites from overlapping BAM reads.
/// `beg` and `end` are inclusive VCF anchors; CIGAR calls retain internal keys.
/// Missing ALT at a simple phased deletion can use Q30/MAPQ30 edit equivalence
/// with a clean SNP confirming the source gauge, including when exact-position
/// REF misses an equivalent shifted ALT. Missing SNP evidence retains the
/// source ALT-absence contrast; contradictory SNP evidence leaves it unknown.
/// Separate allele rows and existing MSA calls are preserved. For a missing insertion
/// call, a Q30/MAPQ30 equivalent shifted edit takes precedence over an exact-
/// coordinate REF call only with the same independent SNP gauge check; no
/// candidate key, genotype or phase label is changed. For an unpaired binary
/// insertion, a verified shifted ALT contradicted by an independent clean SNP
/// stays unknown, never literal REF. Missing SNP evidence and complementary
/// MSA contrasts retain their existing source projection. Site depth fields
/// remain the discovery genotype census, not the enlarged recovery matrix.
int backfill_msa_observations(PhasingChunk& chunk, const Options& opts,
                              hts_pos_t beg, hts_pos_t end);

/// Commit allele recalls after choosing the recovery source. Only sites
/// missing in its original MSA projection are queued; original calls and source
/// genotypes/phase labels remain authoritative. Supplementary CIGAR calls at
/// those queued sites may be corrected or explicitly rejected by local evidence.
bool apply_pending_msa_observations(PhasingChunk& chunk);

// ════════════════════════════════════════════════════════════════════════════
// Step 4 top-level entry
// ════════════════════════════════════════════════════════════════════════════

/**
 * @brief Attribute each co-located deletion read to one candidate.
 *
 * A longer deletion may also satisfy a shorter candidate's window. Keep the
 * longest supported allele and change the other observations to reference.
 */
void make_colocated_deletions_exclusive(std::vector<CandidateVariant>& vars,
                                        std::vector<ReadVariantProfile>& profiles);

/// Process smaller noisy regions first and re-phase when MSA admits candidates.
void collect_noisy_vars_step4(PhasingChunk& chunk, const Options& opts,
                              const VariantKeySet* site_whitelist = nullptr);

/**
 * @brief Align reads in one noisy region and merge MSA candidates into the chunk.
 *
 * Returns the number of admitted candidates, or -1 if MSA could not resolve
 * the region. A nonnegative result marks the region done in the outer loop.
 */
int collect_noisy_vars1(PhasingChunk& chunk, const Options& opts, int noisy_reg_i,
                        const VariantKeySet* site_whitelist = nullptr);

// ════════════════════════════════════════════════════════════════════════════
// Sorting
// ════════════════════════════════════════════════════════════════════════════

/**
 * @brief Return sorted noisy-region indices: label asc, then length asc.
 *
 * Sort noisy regions by variant count (label) ascending, then by length.
 * `label` = variant count recorded during noisy-window detection.
 */
std::vector<int> sort_noisy_regs(const PhasingChunk& chunk);

// ════════════════════════════════════════════════════════════════════════════
// Read and reference collection
// ════════════════════════════════════════════════════════════════════════════

/**
 * @brief Collect indices of non-skipped reads overlapping `[beg, end]`.
 *
 * Collect reads overlapping a noisy region.  Iterates in
 * `chunk.ordered_read_ids` order.
 * Overlap test: `r.beg > end || r.end <= beg` → skip (beg/end from
 * `chunk->digars[read_i].beg/end`; pgPhase equivalent: `ReadRecord::beg/end`).
 */
std::vector<int> collect_noisy_reg_reads(const PhasingChunk& chunk,
                                         hts_pos_t beg, hts_pos_t end);

/**
 * @brief Extract the reference slice `[beg, end]` as 2-bit encoded bytes.
 *
 * Extract reference subsequence for a noisy region.  Encoding follows
 * nt4 encoding: A=0, C=1, G=2, T/U=3, other=4.  `beg` and `end`
 * are clipped to `[chunk.ref_beg, chunk.ref_end]` in-place.
 *
 * @param beg in/out — clipped to chunk.ref_beg if smaller.
 * @param end in/out — clipped to chunk.ref_end if larger.
 * @return Vector of 2-bit bases; empty when the clipped interval is degenerate
 *         or `chunk.ref_seq` is empty.
 */
std::vector<uint8_t> collect_reg_ref_bseq(const PhasingChunk& chunk,
                                           hts_pos_t& beg, hts_pos_t& end);

// ════════════════════════════════════════════════════════════════════════════
// Variant extraction from MSA
// ════════════════════════════════════════════════════════════════════════════

/**
 * @brief Extract `NOISY_CAND_HET` / `NOISY_CAND_HOM` candidates from MSA strings.
 *
 * Extract variant candidates from MSA consensus alignment strings.  Walks
 * consensus-vs-reference and read-vs-consensus alignment strings to locate
 * positions where the consensus differs from the reference, producing new
 * `CandidateVariant` objects and per-read allele profiles for the noisy region.
 *
 * @param noisy_vars     [out] Newly discovered variants.
 * @param noisy_var_cate [out] Category for each new variant.
 * @param noisy_rvp      [out] Per-read allele profiles covering the new variants.
 * @return Number of new variants produced.
 */
int make_vars_from_msa_cons_aln(
    const Options& opts, PhasingChunk& chunk,
    int n_noisy_reads, const std::vector<int>& read_ids,
    hts_pos_t noisy_reg_beg,
    int n_cons,
    const std::array<int, 2>& clu_n_seqs,
    const std::array<std::vector<int>, 2>& clu_read_ids,
    const std::array<std::vector<AlnStr>, 2>& aln_strs,
    std::vector<CandidateVariant>& noisy_vars,
    std::vector<VariantCategory>& noisy_var_cate,
    std::vector<ReadVariantProfile>& noisy_rvp);

// ════════════════════════════════════════════════════════════════════════════
// Merging into chunk
// ════════════════════════════════════════════════════════════════════════════

/**
 * @brief Merge newly discovered noisy variants into `chunk.candidates` and
 *        `chunk.read_var_profile`, then rebuild `chunk.read_var_cr`.
 *
 * Merge MSA-recalled variants into the chunk's candidate table.  Inserts new
 * `CandidateVariant` objects into the sorted candidate table, re-indexes all
 * existing `ReadVariantProfile` entries to account for shifted variant indices,
 * merges `noisy_rvp` allele observations into the affected reads' profiles, and
 * rebuilds the `chunk.read_var_cr` interval tree so that the subsequent
 * k-means call (`kCandGermlineVarCate`) sees the complete, merged profile.
 * When `site_whitelist` is set, unlisted calls are discarded and an exact
 * whitelisted collision may replace an existing repeat-indel row.
 *
 * `replace_sites` explicitly permits replacing untrusted existing rows and
 * their observations when adopting a recovery proposal.
 *
 * @return Number of MSA candidates admitted into the merged table.
 */
int merge_var_profile(PhasingChunk& chunk,
                      const std::vector<CandidateVariant>& noisy_vars,
                      const std::vector<VariantCategory>& noisy_var_cate,
                      const std::vector<ReadVariantProfile>& noisy_rvp,
                      const VariantKeySet* site_whitelist = nullptr,
                      bool admit_all_in_region = false,
                      const VariantKeySet* replace_sites = nullptr);

} // namespace pgphase_collect

#endif
