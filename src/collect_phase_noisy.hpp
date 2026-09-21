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
#include <vector>

namespace pgphase_collect {

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

/// Fill missing observations at admitted MSA sites from every overlapping BAM read.
int backfill_msa_observations(PhasingChunk& chunk, const Options& opts,
                              hts_pos_t beg, hts_pos_t end);

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
