# HiPhase connectivity versus graph/BAM recovery (2026-09-24)

## Reproducible case

The chr20 graph+BAM production run uses the reference, graph catalog, GAF,
and annotated BAM in `test_data/` with `--link-by-alleles`,
`--block-link-window 8`, `--chunk-size 1000000`, eight threads, and the
CHM13 chr20 path. `--dump-phase-matrix PREFIX` writes the pre-stitch hybrid
`recovery-input` and the BAM subsolve `recovery-source` matrices. The
parental-read map is `test_data/derived/chr20_truth_hap.tsv`.

The retained baseline after the next-seam guard has 236,520 phased reads,
228,184 truth-correct and 8,336 discordant (96.476% within independently
oriented read phase sets); its VCF has 457 blocks and 456,233 bp N50. At the
51 Mb seam it leaves the clean SNP at 51,262,081 in PS 51,262,081, the SNPs
at 51,270,774 and 51,286,463 in PS 51,270,774, and the SNP at 51,287,372
in PS 51,287,372. It correctly separates 51,235,063 from 51,262,081 because
no source read bridges that 27-kb interval.

Running HiPhase on **pgphase's own VCF** puts 51,262,081, 51,270,774,
51,286,463 and 51,287,372 in one PS 51,262,081, while splitting
51,235,063. The frozen DeepVariant/HiPhase output makes the same split and
join. This establishes that a different variant catalog is not required to
connect this particular chain.

## How HiPhase makes the connection

Inspected HiPhase source revision `ec2ccdfe21f76b67d39f56b119217d3d56a94d2e`:

- `src/block_gen.rs::get_longest_multispan` and the block generation loop
  admit consecutive heterozygous sites when a read mapping spans them
  (or a qualifying supplementary mapping connects them). The default
  `--min-connecting-reads` is one. With no mapping over the 51,235,063 to
  51,262,081 pair, a new phase problem starts.
- `src/astar_phaser.rs` solves the alleles of each phase problem jointly.
  `src/phaser.rs::get_phase_block_ids` then traverses a graph of *valid
  read allele observations* over the solved heterozygous sites. It assigns
  one PS to a connected component. A read need not span the two outermost
  sites if intermediate site links connect them.

This explanation is specific to the inspected implementation; it does not
assume that HiPhase's A* search always finds a global optimum.

## Signal actually available to pgphase

In the graph+BAM `recovery-input` matrix, all four sites and the following
valid (0/1) paired observations are present:

| Adjacent SNP pair | Informative observations by allele pair | Total |
|---|---|---:|
| 51,262,081 -> 51,270,774 | 0/1: 17; 1/0: 13 | 30 |
| 51,270,774 -> 51,286,463 | 0/0: 2; 1/1: 2 | 4 |
| 51,286,463 -> 51,287,372 | 0/1: 23; 1/0: 26 | 49 |

The BAM `recovery-source` matrix has the same pairs. Its subsolve assigns all
four sites to one PS 51,262,081. The hybrid retains shared graph rows at
51,262,081 and 51,287,372 as graph candidates, while the imported BAM-only
rows retain PS 51,270,774. The imported subset itself has two disconnected
observation components because intervening shared sites were left with the
graph blocks. Treating that subset as one indivisible block loses the
transitive site path from the BAM solve.

The outer stitch transaction for seam 51,262,081-51,287,372 initially joins
both local graph/BAM edges. It then rolls them back because no **single** read
observes the two outermost graph anchors. This rejects the valid 30 -> 4 ->
49 read chain. Its full-block fallback exceeds the 20-variable exact-search
bound; the scoped fallback requires a candidate-count anchor test, which the
single left graph anchor cannot pass even though the direct allele pair has
30 consistent observations. This is a block transfer/stitch limitation, not
missing BAM information.

## Rejected stitch experiments

A local-window end-support relaxation joined the 51,262,081 graph block to
imported BAM sites in a 1-Mb replay. That replay gained 378 phased reads
(360 truth-correct, 18 discordant) and N50 rose from 631,930 to 735,806 bp.
On matched full chr20, however, it yielded 236,521 phased, 227,906 correct,
8,615 discordant and N50 460,310 bp. It wrongly absorbed a pure 279-read
PS 19,414,719 into the large PS 18,999,993: the extra 279 discordant reads
are a whole-block polarity error. This relaxation was reverted.

A narrower experiment allowed the existing two-half direct-SNP bridge when
the whole-block exact problem exceeded its variable bound. On full chr20 it
changed only 28 VCF rows around the target seam (five toward left PS, 23
toward right PS), retained 236,520 phased / 228,184 correct / 8,336
discordant reads, and left N50 at 456,233 bp. It did not solve the complete
chain, so it was reverted too. No experimental stitch change is retained.

## Implication for the next implementation

The replacement must preserve a site-level path across graph and BAM rows,
including shared rows, and cut it at unsupported read-observation gaps. A
whole-block flip of a separately solved graph PS is unsafe: direct boundary
reads only validate its local end, as the 19 Mb trial and the earlier
51,235,063/51,262,081 counterexample show. A focused regression must assert
both the supported 51,262,081 -> 51,287,372 connection and the unsupported
51,235,063 -> 51,262,081 split, plus a full-chr20 truth check to catch
distant switches.
