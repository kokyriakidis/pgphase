#ifndef PGPHASE_COLLECT_PHASE_HPP
#define PGPHASE_COLLECT_PHASE_HPP

// k-means read-haplotype clustering: iterative read→haplotype assignment
// driven by per-variant consensus allele profiles (hap_to_cons_alle),
// up to 10 refinement rounds, plus cross-chunk stitching.

#include "collect_types.hpp"

#include <cstdint>
#include <optional>
#include <set>
#include <unordered_map>
#include <utility>
#include <vector>

namespace pgphase_collect {

/// An oriented heterozygote with a positive PS can anchor recovery. Match the
/// selected graph row's emitted biallelic genotype; unknown and homozygous
/// candidates cannot define a boundary or an observation-path cut.
bool is_phase_set_anchor(const CandidateVariant& candidate);


// ════════════════════════════════════════════════════════════════════════════
// Candidate-category bitmask flags
// ════════════════════════════════════════════════════════════════════════════

// Candidate-category bitmask flags for lcd_var_i_to_cate and k-means target selection.
// Individual category bits:
constexpr uint32_t kLongcalldLowCovVar      = 0x001u;  // low coverage / low alt depth
constexpr uint32_t kLongcalldStrandBiasVar   = 0x002u;  // strand bias
constexpr uint32_t kCandCleanHetSnp          = 0x004u;  // clean het SNP
constexpr uint32_t kCandCleanHetIndel        = 0x008u;  // clean het indel
constexpr uint32_t kLongcalldRepHetVar       = 0x010u;  // repeat/homopolymer het indel
constexpr uint32_t kLongcalldCandSomaticVar  = 0x040u;  // candidate somatic variant
constexpr uint32_t kCandCleanHom             = 0x080u;  // clean homozygous
constexpr uint32_t kCandNoisyCandHet         = 0x100u;  // noisy-region MSA het
constexpr uint32_t kCandNoisyCandHom         = 0x200u;  // noisy-region MSA hom
constexpr uint32_t kLongcalldLowAfVar        = 0x400u;  // low allele fraction
constexpr uint32_t kLongcalldNonVar          = 0x800u;  // non-variant placeholder
// A het that is real enough to call and phase, but whose allele fraction is far
// enough from 0.5 that letting it anchor k-means risks separating paralogs
// instead of haplotypes (see --anchor-af-margin).  Kept out of the anchor mask,
// kept in the germline mask so it is still emitted and phased.
constexpr uint32_t kCandNonAnchorHet         = 0x1000u;

// Composite masks:
// Categories excluded from noisy-region containment checks.
constexpr uint32_t kLongcalldNotCandVarCate =
    kLongcalldNonVar | kLongcalldLowCovVar | kLongcalldStrandBiasVar;
// All het candidate categories.
constexpr uint32_t kCandHetVarCate =
    kCandCleanHetSnp | kCandCleanHetIndel | kCandNoisyCandHet;
// Categories allowed to anchor k-means read assignment.  Narrower than the
// germline mask on purpose: anchoring drives haplotype assignment, so a site
// that is merely callable must not automatically get a vote.
constexpr uint32_t kCandAnchorClean = kCandCleanHetSnp | kCandCleanHetIndel | kCandCleanHom;
// Clean germline categories (used for VCF INFO CLEAN flag and output gating).
constexpr uint32_t kCandGermlineClean =
    kCandCleanHetSnp | kCandCleanHetIndel | kCandCleanHom | kCandNonAnchorHet;
// All germline categories (clean + noisy-recalled).
constexpr uint32_t kCandGermlineVarCate =
    kCandGermlineClean | kCandNoisyCandHet | kCandNoisyCandHom;

// ════════════════════════════════════════════════════════════════════════════
// Public interface
// ════════════════════════════════════════════════════════════════════════════

// Maps a VariantCategory enum to its bitmask flag (0 for unmapped categories).
uint32_t category_to_flag(VariantCategory c);

void stitch_chunk_haps(std::vector<PhasingChunk>& chunks,
                       const Options* opts = nullptr,
                       const PgbamSidecarData* pgbam_sidecar = nullptr);

/// Apply the normal overlap-read stitching rule to 11/12/21/22 haplotype votes.
bool select_stitch_orientation(const std::array<int, 4>& votes,
                               const Options* opts, bool& do_flip);

/// Update block links and orient candidate alleles for one k-means iteration.
int iter_update_var_hap_cons_phase_set(PhasingChunk& chunk,
                                      const std::vector<int>& valid_var_idx,
                                      const Options& opts);

// Assign haplotypes and phase sets to reads via iterative k-means clustering.
//
// Phase 1: initial sweep from the highest-confidence pivot variant outward,
//   assigning each read to hap 1 or 2 based on per-variant consensus allele
//   profiles (hap_to_cons_alle), updating profiles after each assignment.
//
// Phase 2: up to 10 k-means iterations — detect phase-set breaks, flip
//   consensus alleles when conflict reads outnumber agreement, re-assign all
//   reads and rebuild profiles.  Stops early on convergence.
//
// Phase 3: assign chunk.phase_sets[read_i] from the first phased het variant
//   each read overlaps.
//
// Phase 4: fill hap_alt/hap_ref on every CandidateVariant from the finalized
//   consensus alleles (0=ref, 1=alt, -1=unresolved).
//
// `flags` is a bitmask of kCand* constants selecting which candidate categories
// participate.  Masking uses lcd_var_i_to_cate, not counts.category.
// read_var_profile and read_var_cr must be populated before calling.
/// @param anchored Keep the incoming read labels and the consensus alleles a
///                 previous round decided, pinning them through the iteration,
///                 so this round refines that solution instead of resetting and
///                 re-solving. Off by default: the shipped second round discards
///                 the first one.
void assign_hap_based_on_germline_het_vars_kmeans(PhasingChunk& chunk, const Options& opts,
                                                  uint32_t flags, bool anchored = false);

/// True when a read is mapped confidently enough to carry a haplotype call
/// into the output. Reads that fail this still supply allele evidence, so
/// lowering `min_mapq` below `min_assign_mapq` admits their sites without
/// letting an ambiguously placed read own an HP/PS tag.
bool read_carries_phase_tags(int mapq, const Options& opts);


// ── Ported from longcallD, name-for-name ────────────────────────────────────
// These are also used by the direct upstream C differential test. For the BAM
// path, pass preserve_decided=false, msa_sites_vote_without_gap_link=true,
// infer_complement_at_multiallelic=true, and upstream_read_scoring=true.
void var_init_hap_profile_cons_allele(bool is_ont, CandidateTable& variants,
                                     const std::vector<int>& valid_var_idx,
                                     bool preserve_decided = false);
void update_var_hap_to_cons_alle(bool is_ont, CandidateVariant& var, int hap);
int read_to_cons_allele_score(CandidateVariant& var, int hap, int allele_i,
                              bool msa_sites_vote_without_gap_link,
                              bool infer_complement_at_multiallelic,
                              bool upstream_read_scoring);

// These functions carry their upstream names so parity can be checked against
// the original C. Each comment gives the upstream definition.

/// longcallD assign_hap.c:307 -- do two variants agree for one read, under the
/// haplotype it is being tested against? Our `check_agree_alleles` is a
/// SEPARATE helper that consults chunk.haps instead; it is not this function.
int check_agree_haps(const PhasingChunk& chunk, int read_i, int hap, int var1, int var2);

/// longcallD assign_hap.c:151 -- pick the haplotype whose consensus a read's
/// alleles match; -1 when the read carries no usable allele.
int init_assign_read_hap_based_on_cons_alle(PhasingChunk& chunk, int read_i, uint32_t flags,
                                           std::optional<hts_pos_t> phase_set = std::nullopt,
                                           bool msa_sites_vote = false,
                                           bool infer_complement = false,
                                           bool upstream_read_scoring = false);

/// longcallD assign_hap.c:292 -- add this read's alleles to the per-haplotype
/// profile of EVERY variant it covers, under one hap for the whole read.
void update_var_hap_profile_based_on_read_hap(PhasingChunk& chunk, int read_i, int hap,
                                             uint32_t flags,
                                             std::optional<hts_pos_t> phase_set = std::nullopt);

/// longcallD assign_hap.c:270 -- same, against each variant's consensus allele.
void update_var_hap_profile_cons_alle_based_on_read_hap(PhasingChunk& chunk, bool is_ont,
                                                        int read_i, int hap, uint32_t flags);

} // namespace pgphase_collect

#endif
