# Gap read attrition and recovery parity audit

Date: 2026-09-22

## Question

Do pgphase read filters remove the molecules that HiPhase uses to cross the 34
noncentromeric chr20 gaps where the frozen HiPhase output is locally accurate?
If not, where does the evidence disappear?

## Reproduction

`analyze.py` reads the 34-row target table, the parental read truth map, the raw
alignment BAM, and the phased BAM outputs from graph recovery, standalone BAM,
and HiPhase. It counts a qname as bridging only when accepted records reach both
gap flanks. pgphase eligibility uses the recovery floor MAPQ 1 and primary
alignments. HiPhase eligibility uses its default MAPQ 5 and permits
supplementary alignments, combining records by qname.

The committed output is `per_gap.tsv`. The target input is
`evaluations/2026-09-21-short-gap-recovery/remaining_hiphase_correct_gaps.tsv`
and the truth map is `test_data/derived/chr20_truth_hap.tsv`.

## Read attrition result

Across the 34 gaps:

| stage | bridging read names |
|---|---:|
| primary alignments physically spanning both flanks | 860 |
| eligible for pgphase recovery | 859 |
| eligible under HiPhase defaults | 788 |
| gained by HiPhase supplementary handling | 0 |
| tagged by graph pgphase | 772 |
| tagged by standalone BAM pgphase | 649 |
| tagged by both pgphase solves | 593 |
| tagged by HiPhase | 780 |

Nineteen of 34 gaps have fewer than eight bridging reads tagged by both pgphase
solves. BAM recovery admission is not the limiting filter: it admits more
bridging reads than HiPhase, and supplementary records add no unique bridge in
this panel. The graph solve was different. Its former MAPQ 30 GAF floor removed
reads before graph-site clustering and phase-block construction, so the later
MAPQ 1 BAM sub-solve could not restore their effect on the adjacent graph blocks.

HiPhase starts from supplied heterozygous VCF sites, globally optimizes their
diploid orientations with its A-star solver, defines blocks from variants
connected by reads, and assigns one-allele reads after solving. Its advantage in
these gaps comes from a dense set of usable supplied sites and a joint solve,
not from a more permissive read filter.

## Correctness bug

Recovery matched parent graph reads to targeted BAM reads with a monotone
two-pointer scan. Graph reads are qname ordered, while the BAM sub-solve retains
coordinate order. At 4.78 Mb the 493 BAM reads contain 242 qname-order
inversions. The scan reduced complete graph/BAM block votes to 5--0; a hash
lookup over all parent qnames measures 180--1 on the left and 197--0 on the
right. The implementation now builds one `string_view` qname index per graph
chunk and looks up every source read.

Using every recovered vote without chain validation is unsafe. It confidently
joined internally inconsistent BAM blocks and produced wrong graph joins. The
rejected full-chr20 trials were:

| policy | tagged reads | discordant | discordance |
|---|---:|---:|---:|
| previous retained output | 219,579 | 4,759 | 2.17% |
| corrected qname votes, permissive internal replay | 219,396 | 5,206 | 2.37% |
| corrected qname votes, independent failed blocks | 219,396 | 5,098 | 2.33% |
| direct outer relation required everywhere | 219,396 | 4,641 | 2.12% |
| final statistically anchored fallback | 219,396 | 4,641 | 2.12% |

Read-count significance alone did not separate correct and wrong BAM block
joins. Wrong edges had 84--100% apparent read support, while the correct 55.381
Mb internal edge had 86.5%. These observations are correlated through the same
local solve, so a binomial read test overstates independent evidence.

## Retained design

Graph/BAM attachment requires all of the following:

1. a source-specific 2x2 read HP vote with one-sided exact binomial `p <= 0.01`;
2. a sequence-identical clean heterozygote with the same local orientation;
3. an atomic graph-to-graph chain whose outer parity is independently known.

Direct outer-candidate evidence is preferred. A multi-block chain may use the
whole-window gauge only when that gauge passes `p <= 0.01` and each outer
boundary has a consistent shared-candidate vote at `p <= 0.05`. Candidate and
gauge outer tests must agree when both exist. Any failed edge restores the
whole seam snapshot; internal BAM joins are not replayed. Imported BAM phase
blocks and their read assignments remain independent.

The final chr20 output has 59,891 phased heterozygotes, 671 VCF phase blocks,
N50 412,113 bp, 219,396 tagged reads, and 4,641 discordant of 219,233 evaluated
reads: 97.88% parental truth accuracy. The tracked panel spans 12/48 coordinate
cases, including the duplicate coordinate descriptions of the 4.78 Mb gap. No
spanned tracked gap is classified as a haplotype switch by the local truth
scorer. The prior output spanned 14/48 but had 116 more discordant reads and its
two extra physical joins cannot be justified independently from the correlated
BAM block evidence.


## Graph MAPQ experiment and retained policy

Lowering only `collect-graph-variation --min-mapq` from 30 to 5 increased the
full-chr20 candidate count from 77,549 to 80,160 and phased heterozygotes from
59,891 to 61,630. Tagged reads increased from 219,396 to 225,005. Tracked gap
spans increased from 12/48 to 17/48: six current misses closed and one historical
span was lost, for a net gain of five. The local flank truth scorer reports no
`SWITCH` among the 17 spans; the additional cases are `UNCERTAIN` where one
flank lacks enough separately placed truth reads.

Truth accuracy changed from 4,641 discordant among 219,233 evaluated reads
(97.88%) to 6,112 among 224,943 (97.28%). This accepted tradeoff makes MAPQ 5
the graph command default. The BAM/longcallD command keeps MAPQ 30, and the
command-line `--min-mapq` override remains available.

A narrower experiment kept the graph solve at MAPQ 30 and let otherwise
unassigned recovery BAM reads vote against already oriented graph sites. It
left the panel at 12/48 and changed full-chr20 truth discordance by one read.
This confirms that the lower-MAPQ reads must be present while the neighboring
graph phase sets are constructed; adding orientation votes after those blocks
are fixed is too late. That special path was removed.
