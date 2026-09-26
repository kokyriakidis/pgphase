# Graph evidence between recovered BAM phase blocks

Date: 2026-09-23

## Question and method

After targeted BAM recovery, can retained graph sites between two still-separate
BAM-derived phase blocks bridge them correctly?

The current `pgphase` binary was run on four noncentromeric HG002 chr20 regions
using the CHM13 graph-site catalog, indexed GAF, annotated BAM, and reference.
Each run used `collect-graph-variation --bam` with `--phase-matrix-dump`,
`--phase-sites-out`, and `--recovery-audit-out`; all other phasing settings were
defaults. Regions (1-based inclusive) were 15,006,025-15,121,132,
21,773,066-21,894,359, 33,742,638-33,847,964, and
38,233,000-38,343,000. The regional runs each
completed successfully in approximately 1-2 seconds.

For every pair of consecutive, distinct BAM-injected phase sets, the input
matrix was scanned from the left block's last phased candidate to the right
block's first. `bam_injected=0` marks an original graph candidate; it does **not**
necessarily mean BAM never called the same allele. Matrix observations were
compared with the source BAM block's read HP labels. The BAM sub-solve's numeric
PS may differ from the imported PS because recovery allocates a collision-free
label; the comparison maps these labels via the same phased candidate.
Graph-site allele separation was checked against
`test_data/derived/chr20_truth_hap.tsv`, maximizing polarity because 0/1 has
no inherent maternal/paternal meaning. Truth was used only for this audit.

Across the four runs, all **12** inter-BAM-block intervals remained open.
They contain **33** original graph candidates: **27** repeat-het indels and
**six** clean SNP rows (five heterozygous, one homozygous). This sample is
small and selected for multiple recovered BAM blocks; it is not a
chromosome-wide rate.

## Representative gaps

| Region between BAM blocks | Graph evidence | Result |
|---|---|---|
| 15,023,124-15,039,543 | Three unphased graph repeat indels. The best-covered candidate at 15,038,246 separates truth reads only 37/54 (68.5%). It has 46 source-BAM reads in the right block but none in the left block. | No safe bridge. |
| 21,823,067-21,844,360 | The graph deletion at 21,823,068 separates truth reads 65/68 (95.6%) and agrees with the left BAM block on 39 versus 7 of 46 shared reads (one-sided binomial p < 0.000001). It has no right-block shared read. The next two graph indels have near-tied allele links (19:15 and 33:32). | Good one-sided evidence, no supported chain to the right block. It remains unphased as a repeat-het candidate. |
| 33,792,638-33,819,343 | Six graph rows. The graph SNP at 33,797,964 has 55/55 truth-consistent allele observations and its PS joins the left BAM block correctly, but it shares no read with the right BAM block. The next graph SNP at 33,812,283 separates truth reads only 26/50 (52%). An indel at 33,819,337 separates 50/55 (90.9%) but has no left-block shared read. Allele links across the intervening graph rows are weak or tied. | The clean SNP helps attach the left side; no reliable connection reaches the right BAM block. |
| 38,253,203-38,283,562 | Two clean graph SNPs at 38,259,278 and 38,279,170 share graph PS 38,259,278 between two BAM blocks. Their truth allele separation is 68/71 (95.8%) and 55/56 (98.2%). Source-specific graph/BAM read gauges favor the same HP orientation 42:2 on the left and 55:0 on the right; their one-sided binomial p-values are 5.6e-11 and 2.8e-17. All three blocks have the same parental orientation. | **A plausible missed graph-mediated bridge.** The final matrix retains three separate PS labels. The left graph/BAM edge lacks a sequence-identical shared candidate, which the direct stitch gate requires; the fallback also did not close it. A diagnostic same-gauge union of 209 unique truth-scored reads has 202 correct and 7 discordant, with two conflicting HP assignments among shared reads. This is not a production output score. |
| 33,762,265-33,771,787 | One original graph repeat insertion at 33,771,786 shows striking allele/HP association to both BAM blocks (29:1 left; 53:0 right), yet its graph allele separates parental truth only 43/83 (51.8%). A nearby BAM insertion uses a different internal key; biological equivalence was not established. | A strong BAM-label association here is not independent evidence that the graph allele distinguishes haplotypes. Do not join on this candidate alone. |

The two clean graph SNPs in the 33.8 Mb region were already phased before BAM
transfer. Recovery can use established graph phase sets for attachment, and
that is visible for the 33,797,964 SNP. The 17 repeat indels in the first three
regions remained unoriented in the final matrix; a graph site with strong
evidence on only one side cannot establish parity between BAM blocks.

The 38.26-38.28 Mb gap is different. The graph already phased two clean SNPs
into one block between the BAM blocks. The left edge has 44 source-specific
shared-read votes but no sequence-identical graph/BAM candidate, so the direct
`block_gauge_flip` gate refuses it. The right edge has 55 same-gauge votes and
one shared clean candidate. The current atomic seam and fallback logic leaves
the blocks separate. Truth supports a same-gauge join here, but a production
rule must be evaluated on more windows before changing the safety gate.

## Conclusion

Graph sites **can** supply a strong bridge between BAM-derived phase blocks.
The 38.26-38.28 Mb example is a concrete missed join: two already phased
graph SNPs link the blocks in the same parental orientation, while the current
shared-candidate requirement leaves them separate. Joining them would extend
phase-set continuity; the graph SNP block already phases the middle reads, so
this measurement does not establish a gain in tagged-read count. In other
sampled gaps,
repeat indels are noisy or informative only on one side. A rule that uses
independent clean graph anchors to validate the source-specific read gauge
could recover the positive case; treating all interblock graph observations
as bridge evidence would admit the noisy examples too.

## Follow-up experiment (2026-09-24)

The 38.2 Mb graph bridge now joins in the guarded graph/BAM fallback. The
first graph/BAM join had previously marked its BAM block as used, preventing a
second independently supported join. Lifting that restriction unconditionally
joined a 59 kb panel window at 26.03 Mb with no read crossing its intervening
BAM block. The retained version permits reuse only when reads observe both end
sites of the reused block in the expected phase on both haplotypes. It keeps
the 38.2 Mb join (4:1 concordant/conflicting end-to-end observations across
the left BAM block; 2:0 across the middle graph SNP block, one on each haplotype)
and leaves the 26.03 Mb panel window open. A tied full-block MEC result may
use a clean shared candidate plus decisive source-specific read gauge only on
the second attachment of an already validated, end-to-end-supported graph
block; broader use failed existing synthetic atomicity tests.

The graph SNP reported at 38,259,278 has a 19 bp snarl representation. Its
sequence-normalized VCF position is 38,259,286, matching a BAM SNP there.
The existing translated-key lookup already finds that BAM row, but the BAM
subsolve classifies it as a repeat heterozygote, so the clean-candidate gauge
gate excludes it. The two distinct centered SNPs directly pass the full and
both disjoint read-fold parity tests and validate that left edge. The right
edge has a sequence-identical clean candidate and a decisive source-specific
read gauge even though its full-block MEC parity ties.

On the exact 38,233,000-38,343,000 run, a primary block now contains 209
truth-scored reads, 200 concordant and 9 discordant. The same reads and scores
were present in separate blocks before; 124 PS labels change and no HP label
changes. The three other audited regions are output-identical. The chr20 gap
panel passes all 232 prior assertions without changing expectations, and a
focused test checks the new VCF phase-set bridge.

A matched full-chr20 A/B with only the stitcher toggled confirms the local
result at chromosome scale. Truth-scored phased reads move from 236,378 to
236,379; correct from 227,799 to 227,800; discordant remain 8,579. Read phase
sets fall from 871 to 861. No shared read changes truth correctness, and the
single newly phased read is correct. The final VCF adds two variant keys and
reverses the numeric phase of 277 existing heterozygous GTs; it changes no
existing site's heterozygosity. This comparison uses the same current site
representation and VCF writer in both arms.
