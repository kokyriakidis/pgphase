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
constexpr int kStrTractMaxPeriod = 6;       // tandem-repeat tracts a noisy region extends over
constexpr int kStrTractMinLen = 10;
constexpr int kStrTractMaxLen = 500;
constexpr int kNoisySplitFlank = 30;           // a long noisy region's cut point: this far from any candidate
constexpr double kMsaMinConfidentShare = 0.5;  // below this share of confident reads, an MSA takes those alone
constexpr int kMsaMinConfidentReads = 6;         // ... and needs this many of them
constexpr int kLongcalldMinSvLen = 30;
constexpr double kDefaultStrandBiasPvalOnt = 0.01;
constexpr int kDefaultNoisyRegMaxXgaps = 5;
// scratch-buffer noisy k-means, so their phase sets occupy a disjoint namespace
// well above any genomic-coordinate PS id and never collide with core blocks.
constexpr hts_pos_t kGapFillPsOffset = 1000000000;

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
    /// Noisy-region MSA: "plegma", "abpoa" (longcallD's), or "" (plegma when built with it, else abPOA).
    std::string msa;
    /// plegma's homopolymer model: "" the built-in one for read_technology, "none", or a shasta2 csv.
    std::string msa_model;
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

    /// A noisy region longer than max_noisy_reg_len is never resolved by MSA;
    /// when set, it does not suppress the candidates inside it either (see
    /// apply_noisy_containment_filter). Union gap phasing's alignment solve.
    bool keep_calls_in_unresolved_noisy = false;

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
    // accuracy and contiguity.
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
    /// Follow assign_hap.c's phasing rules on the BAM arm.
    /// The graph arm retains its wider link and final read-label handling.
    bool upstream_assign_hap = false;
    // If non-empty, output phased SAM/BAM/CRAM with HP/PS tags.
    std::string output_aln;
    // Output alignment format selected by -S/-b/-C.
    OutputAlignmentFormat output_aln_format = OutputAlignmentFormat::None;
    /** If true, refine phased read alignments from per-read digars before writing output alignment. */
    bool refine_aln = false;
    /** The noisy-region MSA rewrites each clustered read's digars from its alignment to the region's consensus, and
     *  phasing keeps the digars for the output: refine_aln with an aligned output (output_aln, or the union's). */
    bool refine_msa_digars = false;
    /** collect-graph-variation's aligned output (--out-bam/--out-cram/--out-sam): the --bam records of the regions
     *  processed with their reads' union labels, and with refine_aln their refined alignments. Not output_aln, which
     *  also changes how stitching labels reads (stitch_chunk_haps). */
    std::string output_union_aln;
    OutputAlignmentFormat output_union_aln_format = OutputAlignmentFormat::None;
    /** In that output, every record's insertions and deletions left-aligned (--left-align). */
    bool left_align_aln = false;

    /** Coordinate-indexed GAF of the reads' graph alignments (collect-graph-variation). */
    std::string gaf_file;
    /** Site catalog (.pgsites folder or file, or sites VCF) from build-snarl-catalog (collect-graph-variation). */
    std::string graph_sites;
    std::string debug_site; // CHR:POS, emits per-read digar hits to stderr
    // If non-empty, dump the per-read x per-variant allele matrix consumed by
    // k-means to "{prefix}.chunk{id}.flags{flags}.tsv" for offline optimizer
    // experiments (MEC/EM vs the production greedy k-means). Debug-only.
    std::string phase_matrix_dump_prefix;
    /// Diagnostic, off by default: when set (PGPHASE_DIAG, union alignment
    /// solve), collect_var_classify writes every candidate before pruning and
    /// the noisy regions to this directory (aln.TID.BEG.tsv, noisy.TID.BEG.tsv).
    std::string diag_dir;
    /** CLI command string used for PG:CL header field in phased SAM/BAM/CRAM output. */
    std::string command_line;

    /**
     * @brief The first BAM path; empty when there is none.
     */
    const std::string& primary_bam_file() const {
        static const std::string none;
        return bam_files.empty() ? none : bam_files.front();
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
    bool piece = false;  // a noisy region: a piece of one longer than the MSA takes (split_long_noisy_regions)
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
    // The consensus letter before each of msa_insertion_alts (alt_ref_base's encoding), each allele's own: a merged
    // allele can follow a SNP the other does not have (msa_allele_anchor).
    std::vector<uint8_t> msa_insertion_anchors = {};
    uint8_t ref_base = 4;     // 0-3=ACGT, 4=unknown; SNPs only
    uint8_t alt_ref_base = 4; // INS/DEL anchor base: 0-3=ACGT consensus, 4=use ref, >3=gap (skip VCF)
    hts_pos_t phase_set = kUnsetCandidatePhaseSet;
    int hap_alt = 0;
    int hap_ref = 0;
    // True for indels in homopolymer context (set by MSA gap analysis, not by classification).
    bool is_homopolymer_indel = false;
    // True when the site was independently recovered from an MSA consensus.
    bool msa_verified = false;
    // An MSA-verified site whose het genotype the MSA establishes alone (plegma with a homopolymer model, which
    // decides run lengths from both haplotypes' reads): written unphased even where its reads are too noisy to
    // phase it (union EM). abPOA's hets need the reads: written so, they added 420 false hets for 150 true
    // (HG002 HiFi chr12/18/20).
    bool msa_genotyped = false;
    // An MSA call: the reads its noisy region saw (spanning it or not), and those its consensus clusters took (the
    // reads it was genotyped from); 0 for any other candidate.
    int msa_region_reads = 0;
    int msa_cluster_reads = 0;
    // An MSA call's allele on each of the two consensus sequences (1: among that consensus's own variants), -1 for a
    // merged two-allele record or a region with one consensus.
    std::array<int8_t, 2> msa_cons_allele = {-1, -1};
    // An MSA-verified het the union injection keeps out of the solve (too few
    // confidently mapped reads, or the depth gate): it enters only as an
    // island site, never as an EM site (no germline bits, no VCF record).
    bool island_only = false;
    // Phased for the reads but written unphased in the VCF (union EM: an
    // unverified call, or a learned error above kEmVcfMaxSnpError /
    // kEmVcfMaxIndelError).
    bool vcf_unphased = false;
    // A catalog repeat indel (RepeatHetIndel, never an EM site) its allele window phased (union EM pass 2): written
    // with the window's genotype and block.
    bool window_phased = false;
    // One ALT of a multi-allelic site written as its own record (graph output rebuild): its depth counts the site's
    // other alleles, which its REF count does not.
    bool split_allele = false;
    // An injected MSA het (union): the index of the het MSA region it lies in (PhasingChunk::msa_regions of its
    // alignment solve) and the allele that solve put on its haplotype 1, -1 when none. Two calls of one region relate
    // through these as the MSA's consensus sequences do: cis when the same side carries both ALTs.
    int msa_region = -1;
    int msa_hap1_allele = -1;
    // An injected call the alignment solve phased: its phase set there and the allele on its haplotype 1 (-1 none).
    hts_pos_t aln_phase_set = 0;
    int aln_hap1_allele = -1;
    // A clean het the union injection admits although its alleles lack
    // confidently mapped reads (the confident-read exemption): phased like an
    // unverified call (unresolved_noisy) -- it joins no blocks and labels
    // reads only as their last resort.
    bool low_mapq_only = false;
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
    /// A pileup call kept inside a noisy region too long for MSA (see
    /// apply_noisy_containment_filter): never verified, so union phasing
    /// treats it as weak -- it phases and labels, but cannot carry a join.
    bool unresolved_noisy = false;
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

/// The anchor letter of a multi-allelic MSA site's allele i (msa_insertion_alts[i]): its consensus's, else the
/// site's alt_ref_base. 0-3 = ACGT, 4 = the reference's.
inline uint8_t msa_allele_anchor(const CandidateVariant& c, size_t i) {
    return c.msa_insertion_anchors.size() == c.msa_insertion_alts.size() && i < c.msa_insertion_anchors.size()
               ? c.msa_insertion_anchors[i]
               : c.alt_ref_base;
}

/// Whether an insertion or deletion is written as a replacement with no padding base: REF the bases it consumes and
/// ALT its own sequence (key.alt), both from key.pos. Every other indel leads with a padding base at key.pos - 1.
inline bool written_unpadded(const CandidateVariant& c) {
    return c.key.type != VariantType::Snp && !c.key.alt.empty() && c.alt_ref_base == 4 &&
           c.msa_insertion_alts.empty() && (c.key.type == VariantType::Deletion || c.key.ref_len > 0);
}

/// POS of the candidate's VCF record (collect_output.cpp build_vcf_record_core). Not key.sort_pos(): a replacement
/// sorts a base before the first base it writes.
inline hts_pos_t vcf_record_pos(const CandidateVariant& c) {
    if (c.key.type == VariantType::Snp || written_unpadded(c)) return c.key.pos;
    return c.key.pos > 1 ? c.key.pos - 1 : 1;
}

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
    // Island blocks from the union EM: reads it leaves unlabelled, phased by
    // their own call at one site that cannot vote, in a block of that site
    // alone that is never joined. Recorded by the EM, applied after the rescue
    // to reads still unphased (apply_island_labels), through the gap_* channel.
    struct IslandSite { size_t candidate; hts_pos_t phase_set; int hap1; int hap2; };
    std::vector<int> island_haps;
    std::vector<hts_pos_t> island_phase_sets;
    std::vector<IslandSite> island_sites;
    // Trace only (PGPHASE_TRACE; empty otherwise): per read of this chunk, how its calls and its label came about,
    // written to DIR/reads.tsv with the read's phased BAM row (pgphase_trace::read_trace).
    struct ReadTrace {
        bool alignment_only = false;  // added from the alignment solve (add_alignment_only_reads)
        int realigned = 0;            // calls realignment set at injected indels (realign_indel_observations)
        int graded = 0;               // MSA calls realignment graded there
        int filled = 0;               // calls given to an alignment-only read (fill_missing_observations)
        int snp_replaced = 0;         // graph SNP calls the alignment's base replaced (confirm_catalog_snp_calls)
        std::string label;            // the EM's label or why none, then deferred / island labels (union_em.cpp)
        std::string fill;             // a fill label applied after stitching (rescue, deferred or island)
        int aln_hap = -1;             // the alignment solve's k-means label (-1: not in that solve)
        hts_pos_t aln_ps = 0;         // and its phase set there
    };
    std::vector<ReadTrace> read_trace;
    std::vector<ReadVariantProfile> read_var_profile;
    /// A noisy region the MSA resolved (plegma), and where its consensus sequences differ from the reference
    /// (1-based, inclusive; an insertion as the two letters either side of it). `spelled` holds each consensus (two
    /// for a het region, one for a hom one) along the reference: at [p - beg] the letter it aligns to p (none for a
    /// deletion), then any it inserts before p + 1.
    struct MsaRegion {
        hts_pos_t beg;
        hts_pos_t end;
        std::vector<Interval> differs;
        std::vector<std::vector<std::string>> spelled = {};
        // The reads of each consensus sequence's cluster (chunk read indices).
        std::array<std::vector<int>, 2> cluster_reads = {};

        // Two consensus sequences, not the same: the region's haplotypes differ.
        bool het() const { return spelled.size() == 2 && spelled[0] != spelled[1]; }
        // Consensus sequences carrying `alt` (one letter) at p, a base strictly inside the region.
        int carriers_of_snp(hts_pos_t p, char alt) const {
            if (p <= beg || p >= end) return -1;
            int n = 0;
            for (const auto& cons : spelled) {
                const std::string& at = cons[static_cast<size_t>(p - beg)];
                n += !at.empty() && at[0] == alt;
            }
            return n;
        }
    };
    std::vector<MsaRegion> msa_regions;
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
