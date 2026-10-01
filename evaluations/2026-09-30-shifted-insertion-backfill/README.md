# Shifted insertion backfill at 19.373 Mb: confirmed observation bug, rejected broad fix

All trials used the same chr20 HiFi BAM, graph catalog, GAF, and parental truth
only for evaluation. No production phasing rule from this experiment was kept.

## The missing allele information

The 19–20 Mb owning-chunk BAM source has two independent MSA-verified insertion
rows at internal position 19,373,923: four-base `TTCC` and eight-base
`TTCCTTCC`. They are mutually exclusive in the initial MSA matrix. In the
post-backfill source matrix, **47 reads call REF on both rows**. Inspection of
the original BAM shows that 42 of these 47 carry a nearby insertion: 18
four-base, 23 eight-base, and one distinct five-base insertion. The 4- and
8-base events appear at CIGAR positions 19,373,927 or 19,373,933 rather
than the candidate's 19,373,923. Comparing edited reference strings assigns
each four-base event uniquely to `TTCC` and each eight-base event uniquely to
`TTCCTTCC`; a shifted event must not be called REF at both rows. Parental
truth independently shows 20 four-base reads maternal and 26 eight-base
reads paternal in this local BAM screen, but did not enter allele calling.

The existing `backfill_msa_observations` uses `bam_exact_indel_allele`, which
recognizes only the candidate's exact CIGAR position. It calls a shifted
insertion REF if the candidate position itself aligns. That explains the
false REF observations. The initial MSA matrix has seven mutually exclusive
calls at this locus; the later broad backfill floods the source with
position-dependent REF calls.

## Reversible trials

An ALT-only indel-backfill trial removed many valid REF observations across the
19–20 Mb chunk. Truth-scored read accuracy fell from **4,031/4,071** to
**3,743/4,072**, and it did not close the 19.373 Mb gap. It was reverted.

A sequence-equivalent insertion caller validated the original CIGAR event,
inserted bases, reference-edit equality, and base quality before backfill.
The two shifted-read regression observations became their correct ALT rows.
Local truth-scored reads improved to **4,035/4,074** (39 discordant versus
40 in the accepted local run), but the 19.373 Mb gap remained split: only two
one-sided read pairs reach the clean right SNP beyond the final intermediate
insertion. HiPhase leaves that intermediate insertion unphased.

On complete chr20, this broad equivalence trial kept the same **62,361**
variant keys but changed 194 shared-key phased GT fields. Truth-scored tagged reads
changed **236,880 -> 236,919**, correct reads **229,143 -> 229,275**, and
discordant reads **7,737 -> 7,644**. However, it newly joined the protected
36.332 Mb gap at only **140/179 = 78.21%** correctly separated reads, below
its existing 80% floor, and split the established 41.900 Mb complementary
boundary join. A narrower trial that merely suppressed shifted-insertion REF
backfill still caused both regressions, so the problem is not the new ALT calls
alone. Neither trial closed one of the high-confidence tracked open gaps.

Both production trials and the temporary regression were reverted. The
accepted binary again passes the 36.332 Mb and 41.900 Mb focused tests.
Artifacts: `/tmp/pgphase-19373-alt-only-trial/`,
`/tmp/pgphase-19373-equiv-insertion-trial/`,
`/tmp/pgphase-equiv-insertion-full/`, and
`/tmp/pgphase-19373-guard-ref-trial/`.

A durable fix needs the source MSA and backfill to share one
representation-aware allele model, then must revalidate phase-set continuity
using the corrected observations. Correcting backfill after the BAM sub-solve
leaves the source HP/PS labels unchanged, but changes the later graph transfer
and stitch evidence. That mismatch can break previously supported joins.

## Transfer-only counterexample at 4.766 Mb

I next left the BAM source observations and HP/PS labels intact and corrected
only the alleles copied into the graph recovery matrix. A broad guard against
nearby CIGAR insertions improved the 19–20 Mb local read score from
4,031/4,071 to 4,033/4,072 but created an unsupported 36.332 Mb join
(138/179 correctly separated). Requiring a true reference-edit-equivalent
insertion avoided that 36.332 Mb regression and preserved the 41.900 Mb join,
but the full window panel found another regression at 4,766,928–4,792,960:
the accepted `spans=1` became `spans=0`, local separation fell from the
committed >=61% floor to 73/175 (41.7%), and concordance fell below 98%.
Calling matching shifted reads ALT and the other row UNKNOWN, instead of just
suppressing false REF, still split this gap.

The 4–5 Mb source has one BAM PS `4713175` on the left and another `4785719`
on the right. In the accepted graph `recovery-input`, the 4,792,960 graph SNP
is independent. After stitching, both BAM PS and that SNP join under the left
PS. With sequence-equivalent REF suppression, the right BAM PS still attaches
to the graph SNP, but its join to the left BAM PS fails. The source gauge votes
and the 94 shared candidate votes on the right graph flank are unchanged.
The lost edge is the aggregate read-allele connection between the two BAM
source blocks, across a repeat boundary with complementary deletion rows at
4,785,720. This is a real dependency on the insertion REF observations;
removing them without providing a replacement independent bridge is unsafe.

No trial was kept. Artifacts are
`/tmp/pgphase-19373-transfer-guard/`,
`/tmp/pgphase-19373-equiv-transfer/`,
`/tmp/pgphase-4766-equiv-transfer/`, and
`/tmp/pgphase-4766-baseline-audit/`.
