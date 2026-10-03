# Repeat representation and discovery-census guards

## Outcome

No new gap closure is retained in this round. Restore every production option
and algorithm to the accepted binary, SHA256
`b36dd6d8ea32b1f1707b64efc51878b2d8641344ef1b0f009f08eee94b1bac0b`.
The rebuilt final binary has exactly that hash. Retain comments explaining the
backfill contract and two synthetic regression sections that protect it.
No expectation is lowered or refreshed, and no candidate row is merged.

The accepted chr20 state remains 237,199 truth-scored phased reads, 230,072
correct and 7,127 discordant (96.995350%), 661 read phase sets, 63,630 VCF
keys, 333 VCF blocks and 774,189 bp span N50. The coordinate panel remains
86/100 spanned: ten historical competitor targets and four controls are open.
The separate 7.28/24.12 Mb owning boundary cases remain split too. This is a
tracked-panel inventory, not a claim that every chromosome gap is enumerated.
No new accepted chromosome run or runtime comparison is claimed: the final
executable is identical to the previously validated accepted executable.

## Corrected diagnosis at 19.377 Mb

The existing source MSA emits complementary insertion rows with 10-base
`ATATATAGAG` and 14-base `ATATATATAGAGAG` ALTs. Recovery deliberately uses
the port's cluster-derived source observations without enabling the project's
additional local refresh. This preserves the longcallD source behavior.

Instrument a diagnostic copy of the existing refresh, leaving its results
uncommitted. Join the first region's read indices to the initial native source
matrix, whose READ rows are emitted in read-index order. Both consensuses and
all 58 assigned reads have full alignment bounds; every read also physically
covers the variant. The independent local caller gives a valid contrast for
42 reads and unknown for 16. Across the two rows, all 84 callable tuples agree
with their existing source calls: 42 REF/REF and 42 ALT/ALT. The remaining 32
tuples are unknown locally but retain source-cluster alleles. There is no
callable 0-to-1 or 1-to-0 disagreement in this comparison.

This corrects the earlier interpretation that the disabled local refresh
was flipping exact locally recognizable alleles. Literal BAM CIGAR edits and
composed MSA site slices can describe the repeat differently; a CIGAR length
by itself does not prove the equivalent complete MSA allele. The port's
contrast separates parental reads poorly, but enabling local refresh does
not establish the missing phase relation. The auxiliary composition audit
also finds no mismatch between the reference/consensus query sequence and
the consensus/read target sequence anywhere in the owning replay.

DeepVariant labels the 14-base insertion at 19,377,345 homozygous, whereas
the two source MSA consensuses create a heterozygous 10/14 contrast. This is
an unresolved genotype/context difference, not evidence sufficient to demote
one row in production. Local saved DV-HiPhase scoring from the preceding
audit is 111 correct / three discordant overlaps in one block; accepted
pgphase has 94 correct / 31 discordant overlaps in three blocks. The competitor
numbers are saved evaluation results, not a fresh HiPhase benchmark this round.

`msa-read-audit.txt` and `local-versus-source-calls.json` retain the diagnostic
slices, coverage and source-call comparison. The diagnostic source/binary and
native outputs are under `test_data/tmp_gap_next35/`; none is production input.
Truth and competitor calls are used only in evaluation.

## Rejected experiments

| Trial | Useful local effect | Reason rejected |
|---|---|---|
| Enable existing MSA refresh in all targeted solves | Owning 19 Mb errors 41 -> 33 | Correct reads 4,031 -> 4,027; other chunks lose keys, tags and continuity |
| Refresh only exact complementary insertion contrasts | Owning 19 Mb correct reads +1 | Owning 2 Mb errors 62 -> 111; 58 correct reads become wrong |
| Recount supplementary backfill depths immediately | Corrected matrix/AD/AF/strand bookkeeping | Full chromosome produces 13 correct-to-wrong tags and no new tracked span |
| Defer recount until source selection finishes | Four checked owning replays preserve tags | Same full-chromosome regressions as immediate recount; the timing hypothesis is insufficient |
| Enable anchored noisy source solving | Retains some initial SNP orientations | Reopens 1.18 Mb, loses sites elsewhere; owning 19 Mb exchanges 26 correct/wrong pairs without fixing the gap |

Both recount implementations give the identical chromosome totals below.
Each receives a fresh native chromosome run and 105 fresh native panel
requests. They are rejected outputs, not results of the restored binary.

| Metric | Accepted baseline | Rejected recount |
|---|---:|---:|
| Truth-scored phased reads | 237,199 | 237,244 |
| Correct reads | 230,072 | 230,088 |
| Discordant reads | 7,127 | 7,156 |
| Conditional read accuracy | 96.995350% | 96.983696% |
| VCF keys | 63,630 | 63,631 |
| Tracked coordinate spans | 86 | 86 |

The recount adds 32 correct and 18 discordant formerly unphased reads;
13 formerly correct tags become wrong, three correct tags abstain, and two
former errors abstain. One established SNP at 41 Mb changes its relative
orientation. A 33 Mb block connection also appears, but the full trial fails
the accuracy requirement and is not retained as a new closure. No regression
expectation is changed to accept it.

Immediate/deferred native chromosome times are 292.33/265.37 seconds at eight
threads. Their 105-request panel times are 390.22/352.01 seconds with four
workers. These concurrent validation timings are not controlled benchmarks.
Raw native inputs, manifests, outputs and source snapshots are in tmp_gap_next35;
the accompanying JSON files preserve counts and individual truth transitions.

## Retained contract and verification

Post-solve CIGAR backfill extends linkage observations while the depth and AF
fields retain the discovery genotype census. Replacing that census with the
expanded projection changes subsequent retry/transfer/stitch eligibility.
It cannot be treated as an output-only bookkeeping correction. A reliable
future genotype repair must preserve complete allele identities and validate
both local assignments and the established block connection.

Add synthetic insertion and SNP sections checking this contract: calls reach
the matrix, discovery depth/AD/AF/strands and genotype/PS/HP remain fixed, and
repeated backfill is inert. These guards use neither truth nor coordinates
from chr20. The rejected recount implementation fails 15 assertions in the
updated test; the restored implementation passes 1,020 assertions in 42 cases.

The build and all standalone units pass. HiFi/ONT TSV and VCF goldens pass,
including HiFi one/four-thread determinism. The full window suite passes
7,060 assertions in 69 cases. It rescores the prior native cache after checking
its exact accepted-binary hash, all input stats, CLI arguments and output
hashes. Its unmatched owning request executes natively. This verified reuse
is explicitly not called a fresh 105-request accepted run. Production remains
truth agnostic, and all existing span, parental-orientation, read-accuracy,
coverage and key-preservation expectations are unchanged.
