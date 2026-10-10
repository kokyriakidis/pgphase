#ifndef PGPHASE_PHASING_TYPES_HPP
#define PGPHASE_PHASING_TYPES_HPP

// Pipeline-neutral types shared by both the BAM and graph pipelines:
// constants, enums, Options, region/chunk structs, candidate/variant types,
// read profiles, and PhasingChunk.

#include <array>
#include <cstdint>
#include <memory>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

#include <htslib/sam.h>

extern "C" {
#include "cgranges.h"
}

namespace pgphase_collect {


// ════════════════════════════════════════════════════════════════════════════
// Constants
// ════════════════════════════════════════════════════════════════════════════

// Keep the two longcallD phase-set sentinels distinct. Candidate rows start at
// 0 until phasing assigns a genomic anchor; read rows use -1 while unphased.
// A real phase-set anchor is therefore always tested with `phase_set > 0`.
constexpr hts_pos_t kUnsetCandidatePhaseSet = 0;
constexpr hts_pos_t kUnphasedReadPhaseSet = -1;
// Independent BAM solves disagree on this working slot. Unlike an ordinary
// missing call (-1), a later fallback must not silently choose an allele.
constexpr int kConflictingBamAllele = -2;

// Default thresholds and window sizes for variant calling and phasing.
constexpr int kDefaultMinMapq = 30;
constexpr int kDefaultMinBaseq = 10;
constexpr int kMinSvLen = 30;
constexpr hts_pos_t kDefaultChunkSize = 500000;
constexpr double kDefaultMaxVarRatioPerRead = 0.05;
constexpr double kDefaultMaxNoisyFracPerRead = 0.5;
constexpr int kDefaultMinDepth = 5;
constexpr int kDefaultMinAltDepth = 2;
constexpr double kDefaultMinAf = 0.20;
constexpr double kDefaultMaxAf = 0.80;
constexpr int kDefaultNoisyRegSlideWinHifi = 100;
constexpr int kDefaultNoisyRegSlideWinOnt = 25;
constexpr int kDefaultNoisyRegSlideWinShortReads = 25;
constexpr int kLongClipLength = 30;
constexpr int kClipFlank = 100;
constexpr hts_pos_t kReferenceFlank = 50000;
constexpr int kSdustThreshold = 5;
constexpr int kSdustWindow = 20;
constexpr int kNoisyRegMergeDis = 500;
constexpr int kNoisyRegFlankLen = 10;
constexpr int kLongcalldMinSvLen = 30;
constexpr double kDefaultStrandBiasPvalOnt = 0.01;
constexpr int kDefaultNoisyRegMaxXgaps = 5;
// Chunk-stitch abstain margin: adjacent chunks merge only when the absolute
// flip-vote score strictly exceeds this margin.  0 reproduces the original
// behavior (merge on any non-zero vote).  Higher values abstain on weakly
// supported boundaries to avoid over-merging discordant blocks.
constexpr int kDefaultStitchMinMargin = 0;

// Chunk-stitch decision rule (experimental).  Selects how flip_chunk_hap
// decides whether to merge two adjacent chunks:
//   0 = net-margin (default): merge when |flip_votes - noflip_votes| > margin.
//   1 = both-strands-bridged (Reading A): merge when the winning orientation
//       has >=1 read confirming EACH of its two haplotype links (e.g. for a
//       no-flip merge, >=1 read on pre1->cur1 AND >=1 on pre2->cur2).  Guards
//       against merging on evidence from a single haplotype.
//   2 = literal (Reading B): merge when >=2 reads with >=1 supporting flip and
//       >=1 supporting no-flip (stitches on contested seams; diagnostic only).
//   3 = both-strands + net margin (Reading C): Reading A AND the net vote still
//       wins by more than the configured margin.
constexpr int kStitchRuleNetMargin = 0;
constexpr int kStitchRuleBothStrands = 1;
constexpr int kStitchRuleLiteral = 2;
constexpr int kStitchRuleBothStrandsMargin = 3;
// BAM/graph pipelines default to the original net-margin rule.
constexpr int kDefaultStitchRule = kStitchRuleNetMargin;

// scratch-buffer noisy k-means, so their phase sets occupy a disjoint namespace
// well above any genomic-coordinate PS id and never collide with core blocks.
constexpr hts_pos_t kGapFillPsOffset = 1000000000;

// Graph-only het-indel anchor gating (hybrid pipeline only).  Graph het indels
// added to k-means as CleanHetIndel can mis-orient reads when the genotype is
// unreliable.  Keep an indel anchor only when its allele fraction sits within
// kDefaultGraphIndelAfMargin of 0.5 and it has at least
// kDefaultGraphIndelMinAlt supporting alt reads.  Defaults reproduce prior
// behavior (AF window 0.2-0.8, no extra alt floor beyond min_alt_depth).
constexpr double kDefaultGraphIndelAfMargin = 0.11;
constexpr int kDefaultGraphIndelMinAlt = 0;

// ════════════════════════════════════════════════════════════════════════════
// Enumerations
// ════════════════════════════════════════════════════════════════════════════

/**
 * @brief Variant class aligned with BAM CIGAR op types used in digars.
 */
enum class VariantType : uint8_t {
    Snp = 8,       // BAM_CDIFF
    Insertion = 1, // BAM_CINS
    Deletion = 2   // BAM_CDEL
};

enum class VariantCategory : uint8_t {
    LowCoverage,        // below min_depth or min_alt_depth
    LowAlleleFraction,  // below min_af; folded to LowCoverage in final category
    StrandBias,         // ONT strand-bias filter (Fisher exact test)
    CleanHetSnp,        // het SNP passing all filters
    CleanHetIndel,      // het indel passing all filters
    CleanHom,           // homozygous (AF > max_af)
    NoisyCandHet,       // het recalled by noisy-region MSA (one haplotype consensus)
    NoisyCandHom,       // hom recalled by noisy-region MSA (both haplotype consensuses)
    NoisyResolved,      // noisy-region variant resolved after MSA refinement
    RepeatHetIndel,     // het indel in homopolymer / STR context
    NonVariant          // non-variant site (used as placeholder)
};

/**
 * @brief Per-read alignment event type in the digar stream (match, SNP, indel, clip, skip).
 */
enum class DigarType : uint8_t {
    Equal,
    Snp,
    Insertion,
    Deletion,
    SoftClip,
    HardClip,
    RefSkip
};

/**
 * @brief Sequencing preset: affects noisy window size and ONT strand-bias testing.
 */
enum class ReadTechnology : uint8_t {
    Hifi,
    Ont,
    ShortReads
};

enum class OutputAlignmentFormat : uint8_t {
    None,
    Sam,
    Bam,
    Cram
};

// ════════════════════════════════════════════════════════════════════════════
// Options & region input
// ════════════════════════════════════════════════════════════════════════════

/**
 * @brief Parsed CLI state for collect-bam-variation (threads, thresholds, paths, regions).
 */
struct Options {
    int threads = 1;
    int min_mapq = kDefaultMinMapq;
    // Reads below min_mapq never enter a phasing chunk, so that floor governs
    // candidate discovery and block linking as well as haplotype assignment.
    // This floor governs only the assignment: an admitted read below it still
    // contributes its alleles, but is left without HP/PS in the output because
    // its placement is too ambiguous to own a haplotype call. Equal values
    // reproduce single-floor behavior exactly, which is the default.
    int min_assign_mapq = kDefaultMinMapq;
    /// Floor for reads the RECOVERY may use. Reads at or above this but below
    /// min_mapq are excluded from the graph solve but admitted by a targeted
    /// BAM sub-solve over a seam the graph left unphased.
    ///
    /// The reason this floor exists separately: at chr20:26,029,591-26,088,679
    /// the reads carry MAPQ 3 and the default floor of 30 excludes them, so the
    /// pipeline discovers ZERO candidates across 50 kb while a competitor
    /// phases 120 heterozygotes there, 39 of 40 sampled segregating cleanly
    /// against read truth. Admitting them everywhere closes the gap but drops
    /// the left flank from 100% to 93.3% read concordance. Admitting them only
    /// where nothing could be phased is the point of the separate floor.
    int recovery_min_mapq = 1;
    int min_bq = kDefaultMinBaseq;
    int min_depth = kDefaultMinDepth;
    int min_alt_depth = kDefaultMinAltDepth;
    hts_pos_t chunk_size = kDefaultChunkSize;
    double min_af = kDefaultMinAf;
    double max_af = kDefaultMaxAf;
    double max_var_ratio_per_read = kDefaultMaxVarRatioPerRead;
    double max_noisy_frac_per_read = kDefaultMaxNoisyFracPerRead;
    int noisy_reg_slide_win = -1; // if <0, pick default based on read_technology
    // Max gap between adjacent noisy intervals before they are merged.
    int noisy_reg_merge_dis = kNoisyRegMergeDis;
    int min_sv_len = kLongcalldMinSvLen;
    ReadTechnology read_technology = ReadTechnology::Hifi;
    double strand_bias_pval = kDefaultStrandBiasPvalOnt;
    int noisy_reg_max_xgaps = kDefaultNoisyRegMaxXgaps;
    /// Merge co-located MSA alleles into one multiallelic record.
    ///
    /// OFF by default, because longcallD does not do it: measured on its own
    /// chr20 output, upstream emits 0 records with a comma in ALT and 2,401
    /// positions carrying two biallelic rows, and the alignment arm is a port.
    /// With the merge on, our alignment arm emits 1,664 comma-ALT records and
    /// 1,657 GT 1|2, which accounts for 1,532 positions where upstream writes
    /// two rows and we write one.
    ///
    /// The merge exists because splitting measures each allele against a
    /// reference no read carries, so allele fraction runs to 1 and both halves
    /// can classify homozygous. Keeping it available, off by default, makes
    /// that a measurable choice rather than a silent divergence.
    /// ON by default, and the reason is measured rather than preferred.
    ///
    /// longcallD does not merge: its own chr20 output has 0 records with a
    /// comma in ALT and 2,401 positions carrying two biallelic rows. But its
    /// split is COHERENT -- counting positions where one haplotype claims two
    /// different ALT alleles, upstream scores 0. Ours does not: with the merge
    /// off our alignment arm scores 412 such positions. So parity with upstream
    /// is not "stop merging"; it is "make the split complementary the way
    /// upstream's is", and until that is done the merged form is the correct
    /// one to ship.
    bool merge_colocated_msa_alleles = true;
    /// Collapse two co-located haplotype-specific alleles into one record.
    ///
    /// Where both haplotypes carry an allele at one position, this project
    /// reduces the pair to a single row -- as a homozygous call when the two
    /// alleles are equal, or by keeping one when they differ. longcallD emits
    /// one biallelic row per haplotype instead: at chr20:3,863,176 it writes
    /// `C>CAAAAAAAAA 1|0` beside `C>CAAAAAAAAA 0|1` where we write `1|1`, and at
    /// 7,829,788 `C>CA 1|0` beside `C>CAAA 0|1` where we keep only `C>CA`.
    /// That accounts for 31 of the 38 positions where the two disagree on
    /// alleles. Off for the ported path.
    bool collapse_colocated_alleles = true;
    /// Re-derive every two-cluster MSA candidate's counts and per-read profile
    /// from `call_local_msa_allele` after `update_cand_var_profile_from_cons_
    /// aln_str2` has already built them.
    ///
    /// This step has NO counterpart in longcallD: `update_cand_var_profile_
    /// from_cons_aln_str21` (collect_var.c:2178) is upstream's only producer of
    /// these counts, and upstream never revisits them. It is this project's own
    /// addition, and it is a second, disagreeing classifier layered over a
    /// faithful port.
    ///
    /// What it overwrites is load-bearing. In the port, a read belonging to the
    /// OTHER cluster than the candidate is recorded as REFERENCE outright
    /// (`allele_i = 0`, collect_var.c:2205) -- that forced complementarity is
    /// exactly how upstream's two rows at one locus come out `1|0` and `0|1`
    /// "by construction". `call_local_msa_allele` re-measures those reads and
    /// calls many of them ALT at a nested shorter allele too, collapsing the
    /// complementary pair into one wrongly-homozygous row. Measured on chr20:
    /// of the 751 upstream-only records at two-row loci, 664 come back when
    /// this is off; the median NOISY_CAND_HET depth rises from 52 to its true
    /// 61; record identity goes 99.13% -> 99.60%.
    bool refresh_msa_observations = true;
    /// Add observations, from reads the MSA could not place, to candidates that
    /// have none for them. No counterpart upstream: longcallD's depth at a noisy
    /// candidate is exactly the reads its two cluster alignments cover.
    bool add_unplaced_msa_observations = true;
    /// Match longcallD's raw-reference versus nt4 comparison for MSA insertions.
    bool upstream_msa_insertion_hp = false;
    /// In the clean rounds, re-score a read once per phase set it spans and
    /// update each block's allele profile with that block's own verdict.
    ///
    /// No counterpart upstream: `iter_update_var_hap_to_cons_alle`
    /// (assign_hap.c:437) computes ONE hap per read and applies it to every
    /// variant the read covers, with no phase-set scoping anywhere.
    bool phase_set_scoped_clean_rounds = true;
    /// Let a two-cluster MSA candidate vote on which haplotype a read belongs
    /// to.
    ///
    /// `read_to_cons_allele_score` returns 0 -- no vote -- for any candidate
    /// with `msa_insertion_alts`, which is the whole NOISY_CAND_HET class. longcallD's counterpart (`assign_hap.c:127-147`)
    /// has no such condition: every candidate in the mask votes. The effect is
    /// not a small one, because it is self-reinforcing -- in a region whose only
    /// candidates come from the noisy MSA, no read scores, `n_vars_used` stays
    /// 0, `init_assign_read_hap_based_on_cons_alle` returns -1, and the
    /// unlabelled reads are then counted toward BOTH haplotype profiles, so the
    /// argmax calls reference twice and the site is dropped at output for
    /// carrying no ALT.
    ///
    /// Measured at chr20:3,870,827 (DP 72, 43/29, phased): our profiles are
    /// [25,17] and [31,27], consensus [0,0], nothing emitted, while upstream on
    /// the identical depth emits 0|1. Over 20 sampled sites of that class we
    /// labelled 20% of spanning reads against upstream's 45%.
    bool msa_sites_vote_without_gap_link = false;
    /// Infer the unknown haplotype's consensus allele as the complement at any
    /// site, not only a biallelic one.
    ///
    /// longcallD does this unconditionally (`assign_hap.c:141-142`): if one
    /// haplotype's consensus is -1 and the other is known, the missing one
    /// becomes `1 - other`, in place, so the site scores for every later read.
    /// We restricted it to sites with exactly two allele slots on the grounds
    /// that a multiallelic site has no unique complement -- true in itself, but
    /// it means a site with three slots and one unknown haplotype returns 0 for
    /// every read instead of voting.
    ///
    /// That is what withholds the labels behind the 254 records upstream emits
    /// and we do not: in chr20:3.86-3.88 Mb the only het candidates are noisy
    /// MSA sites, and reads there span up to 15 usable candidates yet end with
    /// n_vars_used == 0.
    bool infer_complement_at_multiallelic = false;
    /// Score reads exactly as longcallD does, dropping four restrictions this
    /// project added to `init_assign_read_hap_based_on_cons_alle`,
    /// `read_to_cons_allele_score` and `update_read_phase_set`:
    ///
    ///   - a homopolymer indel is skipped unconditionally (`assign_hap.c:166`);
    ///   - the clean agree/conflict tallies count any clean SNP, heterozygous or
    ///     homozygous (`assign_hap.c:174`), where we counted only a clean het SNP
    ///     with at most two alleles;
    ///   - a clean het SNP or indel weighs 2 regardless of allele count
    ///     (`assign_hap.c:130-131`), where we required at most two alleles;
    ///   - a read takes the phase set of the first heterozygous site it covers
    ///     (`assign_hap.c:328-336`), where we additionally skipped homopolymer,
    ///     noisy-hom and ungap-linked sites and required the read's own allele to
    ///     match one of the two consensus alleles.
    bool upstream_read_scoring = false;
    bool link_earned_repeat_indels = false;
    int link_earned_min_reads = 15;
    double link_earned_min_purity = 0.90;
    // Chunk-stitch abstain margin (see kDefaultStitchMinMargin).  Adjacent
    // chunks merge only when |flip_hap_score| > stitch_min_margin.
    int stitch_min_margin = kDefaultStitchMinMargin;
    // Chunk-stitch decision rule (see kStitchRule* constants above).
    int stitch_rule = kDefaultStitchRule;
    // Graph-only het-indel anchor gates (hybrid pipeline).  See constants above.
    double graph_indel_af_margin = kDefaultGraphIndelAfMargin;
    /** Max |AF-0.5| for ANY graph site to vote in k-means.  Default 0.5 is
        off: min_af/max_af already bound AF to [0.20, 0.80], and comparing
        against 0.30 would reject AF exactly 0.20 on floating-point rounding. */
    double anchor_af_margin = 0.5;
    /** Compute a graph site's allele fraction against total site depth rather
        than against the ref+alt pair only.  Identical for biallelic sites;
        differs only where a site has several observed alleles. */
    bool af_vs_site_depth = false;
    int graph_indel_min_alt = kDefaultGraphIndelMinAlt;
    // When true, keep step-4 noisy-region MSA variant recall (so noisy variants
    // still appear in the output VCF) but skip the kCandGermlineVarCate k-means
    // re-run that re-orients reads using those noisy candidates. On the hybrid
    // pipeline this re-orientation phased ~8k extra reads at ~65% error and
    // poisoned the BAM-shared core (Hamming 0.71% -> 3.11%); skipping it brings
    // hybrid accuracy to graph-pipeline level (99.18%) at higher contiguity.
    // Defaults true for hybrid (set in collect_hybrid_variation), false for the
    // BAM pipeline. See CHECKPOINT.md "Hybrid step-4 re-orientation".
    bool skip_noisy_kmeans = false;

    /// Anchored stage 2: the second k-means KEEPS stage 1's read labels and the
    /// consensus alleles stage 1 decided, and refines within that gauge instead
    /// of resetting and re-solving over the wider site set.
    ///
    /// As shipped, stage 2 discards stage 1 entirely -- collect_phase.cpp clears
    /// every read's haplotype and phase set and re-initialises every
    /// participating site's consensus, then sweeps outward from a pivot chosen
    /// over the NEW site set. A noisy site can therefore overturn the parity the
    /// clean sites established. longcallD resets the same way (assign_hap.c:491)
    /// while its own comment (collect_var.c:2939) says a round should use the
    /// previously obtained phasing as initialization, so anchoring is a knowing
    /// divergence from upstream in the direction upstream documented.
    bool anchored_stage2 = true;
    // Minimum WFA score gap between the two haplotype consensuses before an
    // excluded read is admitted as bridge evidence.  1 is "strictly better
    // wins", which measured +139 discordant reads on chr20; 24 (4x the
    // mismatch cost) was the only setting that beat the baseline on both
    // accuracy and contiguity.  See evaluations/.../chr20_15019294_15130077.md.
    //
    // Renamed when the private-whitelist mode was removed: the margin was never
    // specific to it, and governs the MSA consensus rescoring in align.cpp for
    // every noisy region.
    int msa_ambiguity_margin = 24;

    // Step 4: noisy-region MSA options.
    int max_noisy_reg_len = 50000; // skip regions longer than this bp
    int max_noisy_reg_cov = 1000;  // skip regions with more overlapping reads
    int min_hap_full_reads = 1;    // min full-cover reads per hap for hap-aware MSA
    int min_hap_reads = 2;         // min reads per hap total
    int min_noisy_reg_size_to_sample_reads = 10000; // subsample reads in regions >= this bp
    int noisy_reg_flank_len = kNoisyRegFlankLen;     // flank for scoring boundary deletions
    // WFA2 alignment scoring for noisy-region MSA.
    int match      = 2;
    int mismatch   = 6;
    int gap_open1  = 6;
    int gap_ext1   = 2;
    int gap_open2  = 24;
    int gap_ext2   = 1;
    int gap_aln    = 1;   // gap alignment direction: 1=left-most
    double partial_aln_ratio = 1.1; // max longer/shorter ratio for partial alignment
    int verbose = 0;
    // If false (default), skip VCF rows with non-ACGT REF/ALT bases.
    bool output_ambiguous_bases = false;
    bool include_filtered = false;
    bool autosome = false;
    bool input_is_list = false;
    std::string region_file;
    std::vector<std::string> regions;
    std::string ref_fasta;
    std::string bam_file;
    std::vector<std::string> bam_files;
    std::string output_tsv = "output.tsv";
    std::string output_vcf;
    std::string output_phased_vcf;
    std::string output_phased_bam;
    /** If non-empty, write a TSV of sites dropped during graph candidate collection. */
    std::string output_filtered_sites;
    /** If non-empty, write retained graph site candidates with source SITE_ID. */
    std::string output_phase_sites;
    /** If non-empty, write a per-read TSV of phasing evidence (observations, agree/conflict). */
    std::string output_phase_reads;
    /** Minimum (clean-SNP agree - conflict) margin to commit a read to a haplotype. 0 disables. */
    int min_read_hap_margin = 0;
    /** Spanning reads required to carry a phase block across two adjacent het
        variants.  Below this the block is cut.  Lowering it trades a higher
        chance of a wrong link for longer blocks; 2 is the historical value. */
    int min_block_link_reads = 2;
    // How many preceding het variants to consider when looking for a spanning-read
    // link.  1 reproduces the original adjacent-pair chain; higher values let a
    // block survive a single weakly covered variant by linking across it.
    int block_link_window = 1;
    /// Follow assign_hap.c's phasing rules on the BAM arm.
    /// The graph arm retains its wider link and final read-label handling.
    bool upstream_assign_hap = false;
    // Link variants by the allele pattern a read carries rather than by agreement
    // with the read's assigned haplotype, so untagged reads still contribute
    // connectivity.
    bool link_by_alleles = false;
    // Emit and phase hets that fail the anchor AF margin instead of discarding
    // them.  They still never anchor k-means; this only restores their output.
    bool emit_nonanchor_hets = false;
    // Widen only the GAF read query around each chunk.  Sites are loaded on the
    // unpadded window, but a snarl near a chunk edge is traversed by reads whose
    // alignment starts outside it; without padding those sites see zero evidence
    // and are dropped as if uncovered.
    int gaf_pad = 0;
    // Treat a multi-allelic snarl as one locus: for each alt, contrast reads
    // carrying it against reads carrying any other allele, instead of against
    // the graph reference alone.
    bool snarl_allele_phasing = false;
    // With snarl_allele_phasing, keep multi-allelic snarls as single n-allelic
    // anchors instead of scoring each alt separately.
    bool snarl_keep_whole = false;
    // Minimum share of a snarl's reads that must fall on its two best-supported
    // alleles for it to be used as a phasing anchor.
    double snarl_top2_frac = 0.9;
    // If non-empty, output phased SAM/BAM/CRAM with HP/PS tags.
    std::string output_aln;
    // Output alignment format selected by -S/-b/-C.
    OutputAlignmentFormat output_aln_format = OutputAlignmentFormat::None;
    /** If true, refine phased read alignments from per-read digars before writing output alignment. */
    bool refine_aln = false;

    /** Optional pgbam sidecar file used for fallback chunk stitching when overlap reads have no signal. */
    std::string pgbam_file;
    int pgbam_primary_polarity_margin = 2;
    int pgbam_primary_min_winning_threads = 2;
    bool pgbam_cleanup_pass = true;
    int pgbam_cleanup_polarity_margin = 2;
    int pgbam_cleanup_min_winning_threads = 1;
    bool pgbam_relaxed_cleanup_pass = true;
    int pgbam_relaxed_cleanup_polarity_margin = 1;
    int pgbam_relaxed_cleanup_min_winning_threads = 1;
    /** Reference sample name for GBZ interval queries (e.g. "CHM13"; auto-derived from FASTA if empty). */
    std::string graph_sample;
    /** Optional GBZ-base database for future graph-native snarl/read queries. */
    std::string gbz_db;
    /** Optional raw GAF alignments for future graph-native read traversal support. */
    std::string gaf_file;
    /** Optional GAF-base database/cache for graph-native ReadSet queries. */
    std::string gaf_db;
    /** Optional precomputed vg deconstruct VCF used as a development/debug graph-site catalog. */
    std::string graph_sites_vcf;
    /** Optional whitelist of BAM-derived sites added to graph+private joint phasing. */
    /** Keep GAF as the sole evidence source at every graph-represented site. */
    /** BED intervals where clean BAM candidates replace graph/GAF evidence. */
    /** Suppress output assignments from phase sets with fewer phased reads. */
    std::string debug_site; // CHR:POS, emits per-read digar hits to stderr
    // If non-empty, dump the per-read x per-variant allele matrix consumed by
    // k-means to "{prefix}.chunk{id}.flags{flags}.tsv" for offline optimizer
    // experiments (MEC/EM vs the production greedy k-means). Debug-only.
    std::string phase_matrix_dump_prefix;
    /** CLI command string used for PG:CL header field in phased SAM/BAM/CRAM output. */
    std::string command_line;

    /**
     * @brief First BAM path: `bam_files.front()` when multi-input, else legacy `bam_file`.
     */
    const std::string& primary_bam_file() const {
        return bam_files.empty() ? bam_file : bam_files.front();
    }

    /** @brief True when `--ont` mode (Fisher strand bias, shorter noisy window default). */
    bool is_ont() const {
        return read_technology == ReadTechnology::Ont;
    }

    /** @brief True when `--short-reads` mode. */
    bool is_short_reads() const {
        return read_technology == ReadTechnology::ShortReads;
    }
};

// ════════════════════════════════════════════════════════════════════════════
// Region & read data
// ════════════════════════════════════════════════════════════════════════════

/**
 * @brief User-specified inclusion interval on one contig (1-based inclusive coordinates).
 */
struct RegionFilter {
    bool enabled = false;
    std::string chrom;
    hts_pos_t beg = 1;  // 1-based inclusive
    hts_pos_t end = -1; // 1-based inclusive, -1 means contig end
};

/**
 * @brief One tiled window on a contig plus neighbor metadata for chunk-boundary overlap logic.
 */
struct RegionChunk {
    int chunk_id = -1;
    int tid = -1;
    hts_pos_t beg = 1; // 1-based inclusive
    hts_pos_t end = 0; // 1-based inclusive
    int reg_chunk_i = -1;
    int reg_i = -1;
    int prev_chunk_id = -1;
    int prev_tid = -1;
    hts_pos_t prev_beg = 0;
    hts_pos_t prev_end = 0;
    int next_chunk_id = -1;
    int next_tid = -1;
    hts_pos_t next_beg = 0;
    hts_pos_t next_end = 0;

    /** @brief True if `prev_chunk_id` points to the preceding chunk on the same contig. */
    bool has_prev_region() const { return prev_chunk_id >= 0; }
    /** @brief True if `next_chunk_id` points to the following chunk on the same contig. */
    bool has_next_region() const { return next_chunk_id >= 0; }
};

/**
 * @brief Candidate identity: contig, 1-based position, type, ref length, alternate sequence.
 */
struct VariantKey {
    int tid = -1;
    hts_pos_t pos = 0; // 1-based. Insertions are between pos-1 and pos.
    VariantType type = VariantType::Snp;
    int ref_len = 0;
    std::string alt;

    // Position key for sorting; indels use pos-1 so they sort before the anchor base.
    hts_pos_t sort_pos() const {
        return type == VariantType::Snp ? pos : pos - 1;
    }
};

/**
 * @brief Closed interval on the reference with optional integer label (noisy size, merge metadata).
 */
struct Interval {
    hts_pos_t beg = 0; // 1-based inclusive
    hts_pos_t end = 0; // 1-based inclusive
    int label = 0;
};

/**
 * @brief One CIGAR/MD/CS-derived event on a read (position, type, length, bases, quality flag).
 */
struct DigarOp {
    hts_pos_t pos = 0; // 1-based reference coordinate
    DigarType type = DigarType::Equal;
    int len = 0;
    int qi = 0;
    bool low_quality = false;
    std::string alt;
};


/** @brief `std::unique_ptr` deleter calling `bam_destroy1`. */
struct AlignmentDeleter {
    void operator()(bam1_t* p) const { bam_destroy1(p); }
};

/** @brief `std::unique_ptr` deleter calling `cr_destroy`. */
struct CgrangesDeleter {
    void operator()(cgranges_t* p) const { if (p) cr_destroy(p); }
};



/** @brief An output-only HP/PS assignment for a read absent from graph profiles. */

/**
 * @brief One input read after parsing: coordinates, digars, qualities, noisy subregions.
 */
struct ReadRecord {
    int tid = -1;
    int input_index = 0;
    hts_pos_t beg = 0;
    hts_pos_t end = 0;
    bool reverse = false;
    int nm = 0;
    int mapq = 0;
    std::string qname;
    // Kept for future phasing/BAM-output consumers; candidate collection itself works from digars.
    std::unique_ptr<bam1_t, AlignmentDeleter> alignment;
    // Per-base qualities are copied out because allele quality/depth accounting needs them after parsing.
    std::vector<uint8_t> qual;
    // Ordered per-read alignment events: =/X/I/D/clip/refskip operations.
    std::vector<DigarOp> digars;
    // Dense-error intervals discovered while building digars; merged into chunk-level noisy regions.
    std::vector<Interval> noisy_regions;
    bool is_skipped = false;
    /// Skipped solely because mapq < min_mapq, so the recovery may un-skip it.
    /// A read skipped for its variant load is NOT eligible.
    bool skipped_for_mapq = false;
    bool is_ont_palindrome = false;
    int n_clean_agree_snps = 0;    // populated during phasing (Step 2)
    int n_clean_conflict_snps = 0; // populated during phasing (Step 2)
    // Same agree/conflict counting, but restricted to admitted MSA SNP calls
    // (NOISY_CAND_HET, biallelic, key.type==Snp) rather than pre-existing
    // clean candidates. A read whose only informative sites in a stretch are
    // MSA-recovered bridge SNPs has zero clean-SNP margin and would otherwise
    // be silently stripped of its HP tag at output even when its haplotype
    // assignment (which already uses this evidence via hap_scores) is correct.
    // Kept separate from the clean counters so the output-margin gate can
    // require this evidence to come specifically from admitted bridge SNPs,
    // not from ordinary noisy-region recall run without --private-msa.
    int n_bridge_agree_snps = 0;
    int n_bridge_conflict_snps = 0;
    // |hap_scores[1] - hap_scores[2]| from the last init_assign_read_hap_based_on_cons_alle call,
    // and the number of informative variants behind the winning haplotype.
    // The clean-SNP agree/conflict counts above see only germline clean SNPs;
    // these see every category that votes, so they separate a confidently
    // assigned read from a marginal one more finely.
    int hap_score_margin = 0;
    int n_vars_scored = 0;
    int total_cand_events = 0; // total candidate variant events (includes long-clip noisy windows)
};

/**
 * @brief Counts of reads skipped at chunk boundaries when loading from a BAM slice.
 */
struct OverlapSkipCounts {
    int upstream = 0;
    int downstream = 0;
};

/**
 * @brief Parsed pgbam sidecar map: set_id -> graph thread IDs.
 */
struct PgbamSidecarData {
    std::unordered_map<uint32_t, std::vector<uint64_t>> set_to_threads;
};

// ════════════════════════════════════════════════════════════════════════════
// Variant data
// ════════════════════════════════════════════════════════════════════════════

/**
 * @brief Per-candidate depth, strand tallies, categories, and allele fraction.
 */
struct VariantCounts {
    int total_cov = 0;
    int ref_cov = 0;   // reference allele depth (alle_covs[0])
    int alt_cov = 0;   // alternate allele depth (alle_covs[1])
    int low_qual_cov = 0;
    int forward_ref = 0;  // strand_to_alle_covs[0][0]
    int reverse_ref = 0;  // strand_to_alle_covs[1][0]
    int forward_alt = 0;  // strand_to_alle_covs[0][1]
    int reverse_alt = 0;  // strand_to_alle_covs[1][1]
    // Number of distinct alleles (width of alle_covs / hap_to_alle_profile).
    int n_uniq_alles = 2;
    // Per-allele coverage counts.  When empty, only ref_cov/alt_cov are used.
    // When non-empty, k-means init picks the highest-coverage allele (lower index wins ties).
    std::vector<int> alle_covs;
    // Final category after all classification passes (noisy overlap, AF→LowCov rewrite, etc.).
    VariantCategory category = VariantCategory::LowCoverage;
    // Category from the first classification pass only (before noisy-region adjustments).
    VariantCategory candvarcate_initial = VariantCategory::LowCoverage;
    double allele_fraction = 0.0;
};

/**
 * @brief Key plus counts and reference-encoded bases for classification.
 */
struct CandidateVariant {
    VariantKey key;
    VariantCounts counts;
    // Alternate insertion sequences for one MSA site; index 0 remains genomic reference.
    std::vector<std::string> msa_insertion_alts = {};
    uint8_t ref_base = 4;     // 0-3=ACGT, 4=unknown; SNPs only
    uint8_t alt_ref_base = 4; // INS/DEL anchor base: 0-3=ACGT consensus, 4=use ref, >3=gap (skip VCF)
    hts_pos_t phase_set = kUnsetCandidatePhaseSet;
    int hap_alt = 0;
    int hap_ref = 0;
    // True for indels in homopolymer context (set by MSA gap analysis, not by classification).
    bool is_homopolymer_indel = false;
    // True when the site was independently recovered from an MSA consensus.
    bool msa_verified = false;
    /// Discovered by the recovery sub-solve from the alignment and injected
    /// into this chunk, carrying the allele depths and per-read alleles the
    /// solve that saw the reads measured.
    ///
    /// Its haplotype consensus is deliberately NOT carried: hap_to_cons_alle is
    /// expressed in the sub-solve's own haplotype labels, and the parent resets
    /// reads to its own gauge, so pinning it asserts an arbitrary orientation
    /// as fact. Measured on chr20:55,290,000-55,380,000: pinning dropped the
    /// emitted records from 56 to 46 and collapsed the window to one site.
    bool bam_injected = false;
    /// The alignment path vouched for this site: it came back from a solve over
    /// the reads as a clean heterozygote, or as a noisy candidate the MSA
    /// reconstructed and verified.
    ///
    /// This is the distinction the graph arm was missing. In the alignment
    /// pipeline RepeatHetIndel is a POINTER to a noisy region, not a verdict --
    /// classify_cand_vars_pgphase adds the locus to noisy_var_cr, the MSA
    /// rebuilds it, and it re-enters as NoisyCandHet. In the graph arm the same
    /// label is terminal: over chr20:22,930,000-23,060,000 the alignment arm
    /// has 17 NoisyCandHet and 0 RepeatHetIndel where the graph arm has 0 and
    /// 17, the same loci. A site carrying this tag has been through the
    /// alignment's own verification and is not screened out again on reference
    /// context alone.
    bool alignment_verified = false;
    // Site identity from the graph catalog, independent of observed coverage.
    bool graph_site = false;
    // True when the variant's VCF POS falls inside the chunk's active region.
    // Used during tiling-overlap dedup: prefer the copy that passes this gate.
    bool lcd_make_variants_region_pass = true;
    // Category bitmask for VCF output filtering.  Composed from kCand* flag
    // constants; VCF INFO CLEAN is set when this matches kCandGermlineClean.
    uint32_t lcd_var_i_to_cate = 0;

    // Per-haplotype allele count profiles for k-means phasing.
    // Indexed [hap 0–2][allele i]; hap 0 unused, haps 1–2 are diploid.
    // Length matches counts.n_uniq_alles.  Reset each phasing call.
    std::array<std::vector<int>, 3> hap_to_alle_profile{};
    // Consensus allele per haplotype after k-means.
    // -1 = unknown; 0 = ref; 1 = alt.
    std::array<int, 3> hap_to_cons_alle{-1, -1, -1};
};

/** @brief Ordered list of candidates for one chunk or merged batch. */
using CandidateTable = std::vector<CandidateVariant>;

struct ReadVariantProfile {
    int read_id = -1;
    int start_var_idx = -1;
    int end_var_idx = -1;
    std::vector<int> alleles;
    std::vector<int> alt_qi;
    // Parallel GAF observation channel. It preserves graph/BAM conflicts for
    // recovery without changing the BAM allele used by the normal phaser.
    std::vector<int> graph_alleles;
    // Original BAM observations, retained before GAF injection/MSA replacement.
    std::vector<int> bam_alleles;
    std::vector<int> bam_qi;
    // Physical SNP quality certifies agreement with the original CIGAR base.
    // Zero means absent quality or no matching original BAM REF/ALT call.
    std::vector<uint8_t> bam_base_qualities;
    int bam_mapq = -1;  // Mapping quality of this channel's BAM alignment.
};

// ════════════════════════════════════════════════════════════════════════════
// Chunk data
// ════════════════════════════════════════════════════════════════════════════

/**
 * @brief Working state for one region chunk: reads, reference slice, noisy intervals, candidates.
 */
struct PhasingChunk {
    RegionChunk region;
    hts_pos_t ref_beg = 0;
    hts_pos_t ref_end = 0;
    std::string ref_seq;
    std::vector<Interval> low_complexity_regions;
    std::vector<int> ordered_read_ids;
    std::vector<Interval> noisy_regions;
    std::vector<ReadRecord> reads;
    std::vector<std::vector<int>> up_ovlp_read_i;
    std::vector<std::vector<int>> down_ovlp_read_i;
    std::vector<int> n_up_ovlp_skip_reads;
    std::vector<int> n_down_ovlp_skip_reads;
    std::vector<int> haps;
    std::vector<hts_pos_t> phase_sets;
    // only for reads the clean core left unphased (haps[i]==0) and recovered by
    // the scratch-buffer noisy k-means; empty otherwise. Kept separate from
    // `haps` so they never influence cross-chunk stitching (which inspects
    // `haps` and skips hap==0 reads). Applied only at BAM-write time. PS values
    // already carry kGapFillPsOffset.
    std::vector<int> gap_haps;
    std::vector<hts_pos_t> gap_phase_sets;
    // True only for assignments created in the second rescue pass from an
    // exact BAM observation missing in the graph profile. These are fill-only
    // when overlapping chunks already supplied any phased assignment.
    std::vector<bool> gap_from_bam_observation;
    std::vector<ReadVariantProfile> read_var_profile;
    /** Interval tree [start_var_idx, end_var_idx+1) → read_i, built by `collect_read_var_profile`. */
    std::unique_ptr<cgranges_t, CgrangesDeleter> read_var_cr;
    CandidateTable candidates;
    // Interval-tree cache for noisy-read ratio computation, built lazily before classification.
    std::unique_ptr<cgranges_t, CgrangesDeleter> var_noisy_read_cov_cr; // read coverage spans
    std::unique_ptr<cgranges_t, CgrangesDeleter> var_noisy_read_err_cr; // merged XID intervals per read
    std::vector<int> var_noisy_read_marks; // per-read dedup mark (size = reads.size())
    int var_noisy_read_mark_id = 0;
    int chunk_min_qual = 0;
    int chunk_first_quar_qual = 0;
    int chunk_median_qual = 0;
    int chunk_third_quar_qual = 0;
    int chunk_max_qual = 0;
    /**
     * When `erase_non_variant_candidates_fully_in_noisy_spans` removes a contained `NON_VAR` row
     * that had a clean first-pass category, we stash `candvarcate_initial` here so step-4 MSA
     * merges can restore the initial category on the replacement row.
     */
};

} // namespace pgphase_collect

#endif
