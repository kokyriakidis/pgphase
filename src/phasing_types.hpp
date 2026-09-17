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

// Hybrid pipeline default stitch margin.  Graph-augmented chunks carry extra
// overlap reads, so abstaining on weak seams reduces over-merging without
// losing contiguity.  Matched-eval (HG002 chr20) shows margin 10 lowers
// Hamming and raises genome covered / perfect phase sets versus margin 0, at
// no switch/flip cost.  The BAM pipeline keeps kDefaultStitchMinMargin (0).
constexpr int kHybridDefaultStitchMinMargin = 10;

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
// Hybrid pipeline defaults to both-strands-bridged: it improves contiguity at
// no accuracy cost on HG002 chr20 (auN 11.91M -> 12.00M, +163 reads phased,
// Hamming flat).  --stitch-rule overrides it.
constexpr int kHybridDefaultStitchRule = kStitchRuleBothStrands;

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
    // Experimental graph+private mode: run noisy-region MSA only around
    // whitelisted private sites, admit exact-key MSA calls, and re-phase with
    // those calls plus the clean graph core. Off by default pending a stronger
    // bridge-support gate; the first chr20 test preserved N50 but added errors.
    bool private_msa = false;
    // Minimum WFA score gap between the two haplotype consensuses before an
    // excluded read is admitted as bridge evidence.  1 is "strictly better
    // wins", which measured +139 discordant reads on chr20; 24 (4x the
    // mismatch cost) was the only setting that beat the baseline on both
    // accuracy and contiguity.  See evaluations/.../chr20_15019294_15130077.md.
    int private_msa_margin = 24;
    // Trust the whole noisy region once it overlaps a whitelisted window,
    // instead of requiring the MSA call to land on an exact whitelist key.
    // Use this when the whitelist marks WHERE to look (e.g. a phase-block
    // junction) rather than WHAT to expect -- the point of running MSA there
    // is to discover a position no one already knows.
    bool private_msa_admit_all_in_region = false;
    // Try phasing a junction with MSA SNP calls alone before admitting MSA
    // indel calls in the same region. Only meaningful together with
    // private_msa_admit_all_in_region.
    bool private_msa_snp_first = false;
    // Internal last-resort trial: verified homopolymer links only inside this
    // unresolved gap. Negative bounds disable the trial in ordinary rounds.
    // Internal bounds for non-repeat MSA sites admitted for one gap retry.
    // Minimum gap-only reads (no original flank/block assignment at all)
    // sharing one locally-derived phase set before emitting it as a new,
    // independent block. Purely additive: never overwrites an existing
    // assignment. 0 disables independent-block emission entirely.
    // Bound on how many times the bridge/independent-block recovery sequence
    // repeats per batch. Each round can create new, smaller gaps (a fresh
    // independent block now sits between an existing flank and a bridge that
    // previously had nothing to reach); repeating lets the same validated
    // logic reach those, terminating itself on the first round with no
    // progress. 1 reproduces the original single-pass behavior.
    //
    // is keyed to the gap inventory, which every later round changes), so
    // round >=1 rebuilds raw evidence uncached. On chr20 this was observed to
    // grow resident memory past 47 GB within 10 minutes with no sign of
    // bounding -- not yet safe to enable by default. Raise this only with
    // memory monitored, on a pipeline that has evidence-building costs under
    // control for repeated, partial re-invocation.
    // Attempt to bridge each newly-independent block to its own two flanks,
    // reusing its proposal (no re-extraction) against a read index rebuilt
    // post-emission. Off by default: measured on chr20 it added only 6 reads
    // chromosome-wide with 50% accuracy on that specific population -- these
    // gaps already failed the identical vote test once (that is why they
    // could not join originally), and routing through the new block does not
    // change the flank-side vote count, so any success here rides on
    // marginal evidence. Never corrupts an existing read (same guards as
    // emit_independent_gap_block), but not worth enabling until the bridge
    // criterion itself is strengthened for this specific, already-once-
    // rejected evidence.
    // Additively attach reads with no committed pre-recovery haplotype/phase-
    // set at all (GapReadIndex::assignments misses them -- e.g. they sat
    // right at a block boundary and the main pass left them unassigned) to a
    // flank they individually agree with by allele
    // (CandidateVariant::hap_to_cons_alle), same principle as
    // --link-by-alleles / check_agree_alleles for the main phasing pass.
    //
    // Two earlier versions fed this allele agreement directly into
    // stitch_gap_proposal's `votes` (i.e. let it help decide which side a
    // gap joins to, and whether it flips) and both measurably corrupted
    // existing phase sets on chr20 -- the allele-derived signal is self-
    // consistent per read (confirmed by requiring agreement across every one
    // of a read's overlapping-tile-chunk locations, which changed nothing)
    // but simply not reliable enough, concentrated in already-hard regions,
    // to trust for the join/orientation decision itself. See CHECKPOINT.md,
    // 2026-09-15, for both experiments' numbers.
    //
    // The current implementation instead computes `votes`/orientation
    // selection (which side wins, whether it flips) from committed
    // (first_assignment) reads only -- byte-identical to this flag being
    // off -- and consults allele agreement only *after* a join is already
    // accepted on committed evidence alone, purely to attach more reads to
    // it (the same additive contract emit_independent_gap_block uses).
    // Verified on chr20, whole chromosome, against
    // ../pgphase-eval-data/truth/chr20/diplinator_merged.bam: gaps bridged
    // (112) and their orientation are unchanged, byte-for-byte, from this
    // flag being off. Zero of the 184,974 previously-evaluated reads
    // regressed (no concordant -> DISCORDANT transition anywhere); 22,676
    // additional reads got evaluated, 96.4% of them concordant (in line with
    // emit_independent_gap_block's 96.2% on its own, much smaller,
    // population); contiguity improved slightly (176 phase sets vs 179,
    // fewer/larger blocks) because many reads that previously only qualified
    // for a same-PS independent block now attach directly to the correct
    // flank instead.
    /// catalog is too sparse to phase, nothing is assigned and the region becomes
    /// a gap for a later pass to recover. Admitting the BAM sites there and
    /// solving again keeps the work in the normal pass, and confines the cost:
    /// admitting them chromosome-wide doubled the read Hamming error on chr20,
    /// 0.878% -> 1.837%, while spanning only 14 of the 196 gaps.
    /// Give a biallelic candidate the joint two-haplotype orientation whenever
    /// its own allele depths call it heterozygous, not only inside a retry
    /// window. iter_update_var_hap_to_cons_alle otherwise recomputes each
    /// haplotype's consensus independently by majority, so both haplotypes can
    /// select the deeper allele and the site is emitted homozygous -- the case
    /// the multi-allele branch beside it already guards, with the comment
    /// "Independent haplotype majorities can select the same allele twice."
    /// Measured on chr20:48,204,383 (DEL, DP 71, 30 ref / 41 alt, AF 0.577,
    /// NoisyCandHet): emitted 1|1, where hiphase calls 0|1 and uses it to cross
    /// a 41.8 kb interval that otherwise carries no phased heterozygote, so our
    /// chain links across it with no spanning read and picks an arbitrary
    /// orientation.
    bool joint_het_orientation = false;
    /// ON BY DEFAULT. After the first solve, find windows it left unphased and
    /// re-solve each with the noisy class admitted (force_noisy_msa,
    /// skip_noisy_kmeans = false, scoped by retry_windows). Measured on the six
    /// panel windows: 4 of 6 gaps span against 0 without it, 23 heterozygotes
    /// phased strictly inside gaps against 0, at 99.49% read concordance
    /// against 99.67%. Disable with --no-retry-unphased-with-bam.
    bool retry_unphased_with_bam = true;
    /// split_nested_msa_deletions, which then emitted both nested forms of one
    /// tandem-repeat deletion as independent hets -- the exact false bridge that
    /// function exists to prevent.
    bool force_noisy_msa = false;
    /// Windows the retry is re-solving, in reference coordinates. Carried so the
    /// het-seeding repair in iter_update_var_hap_cons_phase_set can be confined
    /// to them: outside a failed window the provisional-label collapse is not
    /// this pass's business to repair.
    std::vector<std::pair<hts_pos_t, hts_pos_t>> retry_windows;
    /// Reads overlapping a window with no phase set before it counts as failed.
    int retry_min_unphased_reads = 5;
    /// Width a failed window must reach before the retry admits sites in it.
    hts_pos_t retry_min_window_bp = 10000;
    // When true (gap recovery only), an MSA-verified private het SNP inside a
    // gap may act as a block-bridge anchor, on the same terms as a recovered
    // MSA indel: it must pass the anchor allele-segregation test that sets
    // gap_link_supported, and a read's observation of it must be confirmed
    // directly in the BAM. The bridge anchor set otherwise admits only exactly
    // clean het SNPs and non-homopolymer MSA indels, which is the evidence a
    // graph gap interior is least likely to contain.
    // When true (hybrid + skip_noisy_kmeans only), recover the reads that
    // skip_noisy_kmeans leaves unphased: re-run the kCandGermlineVarCate k-means
    // into a scratch buffer and adopt its haplotype for reads the clean core
    // left unphased (hap==0), shifting their phase-set ids into a disjoint
    // namespace (PS += kGapFillPsOffset) so no core phase set is renumbered,
    // merged, or re-oriented. Additive: the clean core is never modified. This
    // is the in-binary form of scripts/gapfill.py (1-source / BAM-pipeline gap).

    // Trim graph-only catalog alleles to minimal VCF form before the hybrid
    // indel noise filter, matching apply_graph_noise_filter. Graph catalog
    // alleles are non-minimal (full repeat run on both flanks); trimming shrinks
    // the derived indel length so max_xgaps-bounded repeat detection matches the
    // standalone graph pipeline. Recovers ~342 over-demoted het indels as k-means
    // anchors (hamming 0.819->0.768%, switch 319->291). Defaults true for hybrid
    // (set in collect_hybrid_variation); --no-hybrid-trim disables it.
    bool exp_hybrid_trim = false;
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
    std::string private_sites_vcf;
    /** Keep GAF as the sole evidence source at every graph-represented site. */
    bool graph_authoritative = false;
    /** BED intervals where clean BAM candidates replace graph/GAF evidence. */
    std::string bam_authoritative_bed;
    /** Suppress output assignments from phase sets with fewer phased reads. */
    int min_phase_set_reads = 0;
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
    // Observations of a homopolymer indel the homopolymer tier admitted as a
    // gap's last-resort evidence (CandidateVariant::hp_gap_scorable). Kept
    // separate from the bridge-SNP counters so the margin filter can credit
    // them only where that tier is what phased the read.
    int n_hp_gap_agree = 0;
    int n_hp_gap_conflict = 0;
    // |hap_scores[1] - hap_scores[2]| from the last init_assign_read_hap call,
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
    hts_pos_t phase_set = 0;
    int hap_alt = 0;
    int hap_ref = 0;
    // True for indels in homopolymer context (set by MSA gap analysis, not by classification).
    bool is_homopolymer_indel = false;
    // Recomputed from clean-site evidence for MSA indel bridges in recovery gaps.
    bool gap_link_supported = false;
    // Set only by select_gap_link_sites, only for an MSA-verified homopolymer
    // indel inside the homopolymer tier's gap window. Lets that one site score
    // reads in the gap it was admitted for without relaxing the global rule
    // that homopolymer indels do not contribute read scores.
    bool hp_gap_scorable = false;
    // True when the site was independently recovered from an MSA consensus.
    bool msa_verified = false;
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

// Per-read allele observations across a contiguous range of candidates.
// alleles[i] corresponds to candidate start_var_idx + i.
// Values: 0=reference, 1=alternate, -1=non-informative, -2=low-quality alt.
// alt_qi normally stores a nonnegative BAM query index.  This sentinel records
// an allele independently confirmed by the same read's graph walk.
constexpr int kGraphConfirmedAltQi = -2;
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
    std::vector<std::pair<VariantKey, VariantCategory>> erased_clean_signal_initial;
};

} // namespace pgphase_collect

#endif
