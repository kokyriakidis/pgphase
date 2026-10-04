# Retain fixed-consensus tandem insertion evidence after source selection

## Target and defect

The target is chr20:41,879,449–41,880,908, outside the centromeric exclusion.
The current chromosome audit nominates 36 distinct competitor-spanned gaps;
these nominations require matching parental orientation on disjoint flanks
and at least 95% local and flank purity. They are not an exhaustive count of
all missing connections. The saved same-BAM HiPhase runs are from October 1;
no new competitor timing or updated competitor input catalog is claimed.

The BAM source keeps two complementary insertion rows at 41,880,908:
`T→TCACACACACACACA` and `T→TCACACACACACACACACA`. Original MSA discovery
records only 11 observations per row despite many more physically covering
reads. The local fixed-consensus recall gate admitted homopolymers but rejected
this CA repeat. Supplementary CIGAR projection then recorded ALT absence for
many missing calls at the MSA anchors, weakening the path from the upstream
`A→AGA` insertion at 41,879,449.

## Retained design

Permit two different non-reference insertion lengths sharing a sequence
prefix at the same reference anchor. Keep the existing disjoint one-base
error neighborhoods, surviving flanks, physical coverage and agreement of
both composed fixed-consensus paths. The two candidate rows remain separate;
a third sequence called REF to both rows supplies no diploid observation.

Homopolymer recall keeps its existing timing. Queue newly eligible calls only
at original missing observations, after ordinary source discovery and phasing.
Commit them after recovery chooses its source solve, so the added calls cannot
displace another seam's validated retry. A queued consensus call may replace
supplementary CIGAR backfill but never an original MSA call. Coalesce duplicates,
abstain on conflicting recalls, preserve an assigned read's source membership
and allele gauge, and consume the queue once. Update the MSA observation census
and interval index without changing source genotype, key or phase label.

Require the source block to have one coherent graph gauge: every exact shared
clean SNP must map to one graph PS with one constant allele relation, and at
least one shared clean SNP must exist. Multiple owners, inconsistent allele
orientations and missing clean anchors retain the original source projection.
A local repeat recall cannot certify a connection between independent graph
blocks. Existing transfer and stitch path checks still decide the connection.
Production uses no truth, competitor calls or fixture coordinates. It extends
the existing fixed-consensus observation pass without adding an aligner or
merging variant representations.

The first guarded closure lost eight owning-chunk read assignments, failing the
existing 4,154-read floor. Keep that floor. Exact shifted-insertion ALT corrections
are staged with the same coherent graph-owner guard, only in source blocks
with newly recalled complementary fixed-consensus insertion observations,
when the per-read SNP gauge is absent; contradictory per-read SNPs still abstain. These physical calls never
change discovery counts. Third repeat lengths are explicitly rejected on both
complementary rows; transfer clears older primary and BAM calls rather than
treating the rejection as an ordinary missing slot. Original MSA calls are
untouched. Finally, equally supported independent read-rescue blocks may prefer a directly phased
SNP whose orientation passes the existing primary-read singleton confidence
check. Inferred excluded-site associations cannot break the tie. Preserve equal
confidence ties and never merge blocks.
All 4,154 original owning-chunk assignments now remain phased: 17 old errors
become correct, with no old correct read becoming incorrect.

## Target measurements

The owning 41–42 Mb replay retains all 965 VCF keys and all old phased SNPs.
The repaired 341,860 bp block joins the old 340,401 bp graph block with its
private complementary insertion rows. The downstream 41.900 Mb connection
remains intact, and the outer graph blocks remain separate. VCF blocks fall
from six to five with constant original SNP allele gauges.

| Owning-chunk read metric | Before | After |
|---|---:|---:|
| Truth-scored phased reads | 4,154 | 4,154 |
| Correct | 4,109 | 4,126 |
| Discordant | 45 | 28 |
| Accuracy | 98.9167% | 99.3260% |

Of the 81 parental-truth reads physically overlapping the gap, pgphase changes
from 81 phased / 62 correct / 19 discordant to 81 / 79 / 2. Saved HiPhase with
pgphase calls has 80 / 80 / 0. The VCF connection is correctly oriented; local
read purity still trails HiPhase. Do not claim this single closure
matches all of its read assignments. Independent pgphase read phase sets still
partition these local overlaps.

## Rejected trials

- Internal dropout detection exposes more MSA calls but does not certify a
  complete connection. Restore the original retry admission and scheduling.
- Accepting a shifted binary insertion ALT before source selection without a
  clean SNP source gauge changes the retry matrix. Keep these physical
  corrections pending and require a coherent graph-block gauge before admission.
- Committing the tandem calls before retry selection closes this gap but loses
  20 owning-chunk keys and reopens the downstream connection. Stage the calls.
- Applying physical corrections outside newly recalled blocks drops an old
  5.31 Mb read assignment. Unconditional SNP tie preference adds 20 owning
  34 Mb assignments, ten discordant, from inferred excluded sites; requiring
  a primary-read certificate alone is insufficient. Require a directly phased
  SNP as well. The existing 34 and 35 Mb correctness floors then pass.
  `rejected-broad-physical-parity.json` retains the rejected chromosome trial.
- Staging the calls without the coherent graph-owner guard preserves keys but
  wrongly joins the 19.374 Mb blocks. Full chr20 converts 368 previously correct
  assignments to incorrect and raises discordance from 7,120 to 7,446.
  Reject; `rejected-multigauge-parity.json` and its transfer audit retain the
  counterexample. Keep all existing split expectations.

## Regression

Add the new coordinate case to the committed panel with exact `spans=1`, an
in-gap heterozygote floor and measured accuracy/separation floors. Keep the
upstream insertion in the required-site manifest. Expand the existing owning
41 Mb test to require both separate repeat rows, their opposite ALT orientation,
the upstream insertion's supported genotype relation, the downstream connection
and whole-chunk correct/error limits. The starting executable fails six checks.
The retained executable passes it and the 19.374 Mb conflict guard.

Standalone units and all 1,465 predicate assertions in 47 cases pass. HiFi/ONT
TSV and VCF goldens pass; HiFi one/four-thread outputs are identical. The abPOA
unused SIMD helper warning is pre-existing; no new warnings are introduced.

## Final whole chr20

| Metric | Before | After |
|---|---:|---:|
| Truth-scored phased reads | 237,225 | 237,308 |
| Truth-correct assignments | 230,105 | 230,211 |
| Discordant assignments | 7,120 | 7,097 |
| Conditional read accuracy | 96.998630% | 97.009372% |
| Truth-scored read phase sets | 654 | 654 |
| VCF variant keys | 64,188 | 64,188 |
| VCF phase blocks | 333 | 332 |
| Span N50 | 806,449 bp | 806,449 bp |
| Span NG50 | 684,798 bp | 684,798 bp |
| Coordinate panel spans | 89/103 | 90/103 |

The expanded denominator includes the new coordinate case, open in the starting
output. No existing connection reopens. The sole new VCF union is the target
block, with no lost SNP, mixed old SNP gauge, lost key or new key. All 230,105 old correct
assignments remain correct. Twenty-one previously incorrect reads become
correct and two become unphased. All 85 newly scored reads are correct.
Net +106 correct / −23 discordant / +83 truth-scored. Extra correct recovery outside the target adds
read assignments without merging other VCF blocks. The unchanged N50/NG50
reflect that this 1.459 kb extension does not cross either cumulative threshold.

The full run takes 369.89 s at eight threads while native regression requests
run concurrently; this is validation timing, not a controlled runtime comparison.
Starting outputs: `test_data/tmp_gap_fix45/full-scoped/`.
Final outputs: `test_data/tmp_gap_fix46/full-provenance/`.
Starting binary SHA256:
`772e57fe587c75f7acdfb1b312cf791573194e31ab1a1732bc667a0ea10655a0`.
Final binary SHA256:
`9035673c60ec6798a356e2651b44880ba3b87f2967fb104ee4deed4c3d4b9477`.

`chromosome-parity.json` records exact read truth transitions and variant parity.
`block-transfer.json` records the sole new block connection and all panel spans.
`span-metrics.json`, `owning41-parity.json` and `local-reads.json` retain the
whole-chromosome, owning-chunk and competitor-local measurements respectively.

## Final regression validation

The final executable passes all 7,236 window assertions in 73 Catch2 cases,
covering the 103-coordinate panel and owning-chunk regressions. All old split,
coverage, correct-read and discordant-read requirements remain unchanged.
The new owning41 checks strengthen the measured floors to 4,154 scored /
4,126 correct / at most 28 discordant. Standalone units and all 1,465
predicate assertions in 47 cases pass, including weak and inferred SNP ties
and the co-located certificate counterexample.

Generate all 106 exact native argument vectors with the final executable
(488.02 s with four concurrent requests), then rescore those outputs using
`replay_cached_panel.py`. `native-cache.json` pins the final executable SHA256,
all input sizes/timestamps, every normalized CLI argument and output hashes;
an unknown request runs the native executable. This validates native outputs,
not baseline replay data. The cached scoring pass is separate from the native
pipeline timing. Evaluation cache paths require the retained local fixtures.

Final suite log: `test_data/tmp_gap_fix46/window-provenance-verified.log`.
