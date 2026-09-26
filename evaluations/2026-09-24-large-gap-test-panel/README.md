# Larger chr20 gaps added to the window regression panel

The current full-chromosome pgphase graph+BAM VCF was compared with the saved
HiPhase 1.6.0 VCF and BAM under `/tmp/hiphase-on-pgphase-final-chr20/`.
HiPhase used the same original annotated BAM and pgphase-derived variant calls.
The saved HiPhase input predates the final 528,827 representation correction,
so every selected window was also checked against that input VCF: the complete
set of `(POS, REF, ALT)` records from 50 kb left of the gap through 50 kb
right of it is identical to the current pgphase VCF. The two otherwise
qualifying 23.4-Mb windows with a local call-set difference were excluded.

Selection used consecutive nonoverlapping pgphase phase-set extents with a gap
larger than 10 kb. The exact left and right boundary records had to share a
HiPhase phase set. Gaps intersecting chr20:26.0–29.5 Mb were excluded. That
screen found 94 larger HiPhase-joined gaps with at least 20 truth-labeled
HiPhase reads; 40 had at least 98% local truth purity, at least 95% purity on
each 10-kb flank, and the same majority parental orientation on both sides.

At panel creation, sixteen high-confidence gaps across the chromosome were
replayed with the committed window harness. Pgphase already closed 8,638,940–8,662,670 in that
local context, so the other **15 still-open gaps** were added to
`evaluations/2026-09-16-test-panel/panel.tsv` and their measured floors to
`src/test_gap_windows_expect.tsv`. They range from **18,671 to 26,615 bp**.
Each has at least one primary BAM read physically crossing both boundaries.
The exact boundary alleles, HiPhase phase set, local and flank truth evidence,
read coverage, and both tools' separated-read scores are in `selection.tsv`.

At panel creation, for these 15 windows HiPhase correctly placed **73.17–90.12%** of all
truth-scorable reads overlapping the gap into one phase set. Pgphase's local
single-block baseline floors are **31–50%**, and none of the 15 isolated
replays spanned both boundaries. The test preserves those open states as exact
expectations until a proposed join is checked for molecule support and parental
orientation. It also checks phased sites inside each gap, read concordance,
separated-read coverage, absence of confident flank switches, and absence of
contradictory haplotype alleles at one position.

The panel's pre-existing rows and accuracy floors were retained. Only the 15
new rows and the TOTAL in-gap-heterozygote floor changed. The `competitor`
label `hiphase_same_sites` distinguishes these scores from older panel rows
that came from a different HiPhase call set. No production phasing code changed.

## 2026-09-25 recovery follow-up

The 47,636,774–47,659,910 window now spans in the isolated graph+BAM replay.
Two distinct MAPQ-60 molecules make Q40 physical calls at the clean SNPs
flanking the BAM source's weak cut (47,636,774 to 47,659,910). Both support
the same haplotype connection. Recovery retains that cut as weak for the
source phase set as a whole and uses the two-read evidence only for the graph
seam containing it. Both graph/BAM flank gauges must also have significant
read votes and unanimous exact shared-site orientation. The separated-read
score in the replay rose from 62/147 (42.18%) to 120/147 (81.63%), matching
the saved HiPhase score, with 2 discordant of 422 tagged truth reads and no
confident parental flank switch.

The first trial promoted every quality-validated cut into a complete source
path. A matched whole-chr20 run lost 909 truth-correct reads because a
distant 47-Mb graph seam consumed the same source path. Scoping the cut to
its own seam still lost 356: the next seam reused graph gauge votes saved
before its left block had flipped. Translating that saved gauge into current
candidate orientation fixed the erroneous join. The final matched run has:

| Run | Phased truth reads | Correct | Discordant | Purity | VCF blocks | N50 |
|---|---:|---:|---:|---:|---:|---:|
| Previous | 236,834 | 228,914 | 7,920 | 96.6559% | 392 | 491,812 bp |
| Seam-scoped cut plus current-gauge stitch | 236,834 | 228,927 | 7,907 | 96.6614% | 389 | 506,103 bp |

The exact 47–48 Mb chunk regression checks that the newly joined block
has at least 2,000 truth-labeled reads and at least 99% single-block parental
purity. The known 19-Mb wrong-join control remains split. Of the 15 larger
tracked gaps, one now joins in the whole-chromosome VCF. The other 14 remain
open; most have at most one exact callable boundary-spanning read or have
ambiguous indel calls, so this clean-SNP exception does not apply to them.

### 56.323-Mb residual seam audit

The isolated 56,323,427–56,343,002 SNP gap has three physical read pairs
with the same allele relation. Its recovered BAM source phase set, however,
has weak cuts at 56,313,636 and 56,323,427; only the latter passes the
clean-SNP quality check. The 56,323,427–56,350,725 graph seam needs that
latter cut, but the current whole-source block also contains the earlier,
unsupported cut. Joining the entire source phase set would propagate an
unproved orientation outside this seam. The safe seam-scoped rule leaves it
split until source runs can be represented as separate phase blocks at their
weak cuts. The diagnostic replay used `--phase-matrix-dump` on
`CHM13#0#chr20:56273427-56393002`; temporary instrumentation was removed.

## 2026-09-25 sparse SNP bridge follow-up

A source block with an unsupported cut outside the 56.323-Mb seam had a
quality-supported cut inside it. The seam stitch connected the graph flanks,
but a later one-sided source-run pass moved the left boundary out of that
phase set. Keeping such a locally bridged source block out of both later
source-run passes preserves the correct 56,323,427–56,343,002 join. Its
regional truth score remains 493/495; the separated-read score becomes
157/183, matching HiPhase.

A first direct-SNP trial let one Q40/MAPQ-60 molecule join whole graph blocks.
It closed the 24.357- and 48.022-Mb gaps locally, but full chr20 discordance
rose from 7,907 to 8,615. One false 62.623–62.645-Mb join caused 727
misplaced reads: its sole spanning molecule called only one SNP in each
block. The retained rule requires that molecule to confirm its block
orientation at another clean SNP at least 100 bp away, rejects every
contradictory high-quality pair, and divides its error bound by the number of
nearby pairs tested. It also considers up to three SNPs on each side within
2 kb of the seam, so the clean right-flank SNP beyond 20.875 Mb can testify
when the nominal right endpoint has low base quality.

The locally verified new joins are 20,875,468–20,895,876 (513/514 correct
reads, 120/144 separated), 24,357,615–24,379,633 (433/433 correct,
109/131 separated), and 48,022,132–48,043,584 (486/488 correct,
111/139 separated). The 19-Mb and 62-Mb false-join controls remain split.
The full-chromosome score and final panel result are recorded in CHECKPOINT.md.

The original 21.435-Mb 50-kb-padded panel replay omitted both boundary
variants. The test now uses its full 21–22-Mb owning chunk; its concordance
floor reflects that corrected replay, while its gap remains split.

## 2026-09-25 indel MSA retry follow-up

The 24.581-Mb indel pair is joined by admitting unplaced reads to the existing
BAM MSA sub-solve when paired boundary calls show statistically significant
allele dropout. Its right BAM insertion is keyed one base after the graph VCF
anchor, so the detector must include that boundary. The 63.945-Mb pair is
joined when paired indel calls show statistically significant mixed parity and
neither parity dominates. Both joins pass local parental truth (582/583 and
469/471) and the full chr20 window regression panel (810 assertions).

The accepted full-chr20 output has 7/15 selected larger pairs in one phase set,
236,807 truth-evaluable phased reads, 228,967 correct, 7,840 discordant,
382 VCF blocks, and 517,052-bp N50. Eight selected pairs remain split. The
unprotected MSA trial closes three more locally but misphases 3.573 Mb badly;
those candidate-replacing joins have not been promoted.


## 2026-09-25 full-chunk 56.064-Mb follow-up

The local 120-kb replay was misleading because the production BAM solve merged
this seam with later seams. A grouped MSA retry changed phased rows in those
later seams and was rejected. The focused retry now solves only the complete
adjacent graph phase-set extents, while full-chunk graph read assignments still
orient both flanks. The source block's supported site-to-site path replaces an
impossible direct first-to-last read requirement when only one private BAM row
was injected. This alternative applies only to a focused retry; both
attachments also require direct SNP bridges and matching MEC/read-gauge
parity. The unscoped trial changed 650 truth assignments near the centromeric
26-Mb shoulder (529 improved, 121 worsened) and was rejected. The scoped
output has no truth-assignment changes there.

In the full chr20 output, boundary positions 56,064,697 and 56,083,708 and
the next right-flank SNP at 56,107,742 share PS 56,040,613. Tracked larger-gap
joins rise from 7/15 to 8/15. Truth-evaluable reads: 236,824 phased,
228,992 correct, 7,832 discordant (96.6929% purity), versus 236,807,
228,967, and 7,840 (96.6893%) in the preceding accepted run. N50 remains
517,052 bp and VCF blocks remain 382. Eighteen previously gap-fill-tagged
56-Mb reads become unphased when the graph gap closes; the chromosome still
has a net 17 more phased reads. The full owning-chunk 56-Mb regression and
the 63.945-Mb mixed-parity control pass. The full panel passes 823
assertions across 22 test cases.


## 2026-09-25 complete-source and complementary-row follow-up

The full 48–49-Mb chunk closes 48,929,511–48,950,388 by reusing a focused BAM
phase path with 1,122 sites and zero weak cuts. Both graph/BAM flank gauges
agree, and the boundary-scoped exact MEC solve selects their parity. Its
owning-chunk truth count stays 3,844/3,901. The narrow 120-kb replay remains
split because its retry loses an established flanking source row; a dedicated
full-chunk regression gates the join.

The 41,900,800–41,919,471 seam closes when two complementary BAM deletion
rows are admitted to a focused MSA retry without merging or dropping either.
The retried source has 107 sites and no weak cut. The preceding 41.866-Mb seam
remains split. The local graph panel now expects this span, and an owning-chunk
truth regression checks it. A similar retry at 11.6 Mb had one weak source cut
and harmed read truth; the new complete-path admission guard rejects it, with a
dedicated negative-control regression.

The guarded full chr20 output closes 10/15 selected larger gaps versus 8/15
previously. It has 236,821 truth-evaluable phased reads, 228,988 correct,
7,833 discordant (96.6924% purity), 381 VCF blocks, and 517,052-bp N50.
The previous output had 236,824 / 228,992 / 7,832 (96.6929%), 382 blocks,
and the same N50. The five remaining selected gaps each have one primary BAM
read spanning both exact boundaries; none has independent clean SNP calls on
both flanks from that read. The full panel passes 864 assertions in 25 cases;
unit tests and the main build pass.
