# Calibrated deletion chain at the 56 Mb boundary

The next qualifying uncovered interval belongs to HiPhase's rank-43 block,
chr20:55,815,775–56,293,581 (477,807 bp). All larger uncovered intervals in
`screened-seams.json` have fewer than 80% correct HiPhase reads, counting
unphased original primary reads. This evaluation closes its largest seam,
55,999,194–56,040,612 (41,418 bp). The entire HiPhase block remains split.

| Same original overlapping primary reads | Correct | Discordant | Unphased | Connected core correct |
|---|---:|---:|---:|---:|
| Before | 131 | 0 | 115 | 68 |
| pgphase | 214 | 0 | 32 | 214 |
| HiPhase | 213 | 1 | 32 | 213 |

The denominator is 246: pgphase correctly phases 87.0%, exceeds HiPhase by
one read both in total and connected core, and retains every earlier correct
assignment. The resulting pgphase block is 55,882,616–56,156,527 (273,912 bp).

## Cause and fix

There is no original read spanning both endpoint SNPs. HiPhase uses the
interior repeat variants. Pgphase has the same original BAM and discovers
complementary deletions of four/six bases at 56,007,501 and one/two bases at
56,027,379. Its short source projections miss many six-base ALT observations
when CIGAR places the same deletion at 56,007,532, 31 bases away. Existing
physical bridges cannot compose both complementary deletion contrasts across
this SNP desert, and strict unanimous bridge admission rejects the one
repeat-length error among five original crossing molecules.

The repair reconstructs both alleles through each complete bounded repeat
with external flanks. Separate SNP-bearing molecules calibrate both allele
classes on either side; crossing molecules are disjoint from calibration.
The accepted chain has calibration counts `12,3,0,8` and `9,1,2,9` after
normalizing their SNP gauges, and bridge counts `1,4`. Smoothed empirical
repeat-length errors supply the bridge likelihood, rather than treating base
qualities as repeat-genotype certainty. Calibration and crossing agreement
must each reach 80%, diploid calibration must reject random association, and
bridge likelihood odds must reach 4:1. A poorly calibrated nearby deletion
pair at 56,026,465 is rejected.

New read labels need their own physical allele. Conflicting two-marker calls
abstain. Exact complete-repeat calls can retain the shorter class through one
base of length slippage, including a zero-base CIGAR deletion. Lower-quality
external anchors can restore reads only after the connection is certified;
they never calibrate or bridge blocks. Independently called BAM molecules
without GAF catalog profiles are retained in the connected core. Chunk-boundary
replay carries its certified calls along with the checked endpoint gauges.
Variant descriptions, counts and categories stay intact.

A broad two-consensus retry was tested and rejected: it lost 949 previously
correct owning-region assignments. The production repair leaves the discovery
and source-solve rules unchanged. Parental truth is used only for evaluation.

## Preservation and remaining work

Full chromosome: 250 blocks, N50 **955,496 bp**, largest block **3,012,193 bp**.
Correct reads increase 230,811→230,894; discordant stay 6,644; unphased decrease
19,155→19,074. Two previously absent assignment records now carry correct
labels, explaining why the increase in correct records is two larger than
the decrease in explicitly unphased records. All earlier correct/phased
assignments, block extents and 82,287 variant keys/evidence are preserved.
No accuracy changes occur outside original reads overlapping this seam's
neighborhood. One- and four-thread owners have identical candidate TSVs,
VCF records and BAM assignments.

The full 477,807 bp target is **not closed**. Its overlapping prefix blocks
remain separate, along with seams 56,156,527–56,175,630 and
56,175,630–56,192,736. HiPhase reaches only 128/170 and 104/160 correct reads
at those latter seams. Across the whole target, pgphase has 1,966 correct,
57 discordant and 197 unphased reads, with 1,160 correct in its largest core;
HiPhase has 2,002 correct, 13 discordant and 205 unphased, with 2,002 core
correct. The remaining whole-target deficits are 36 total and 842 core reads.
Those joins are not certified by this change.

## Verification and reproduction

The committed panel, certification manifest, measured expectations, endpoint
requirements and owning-region override include this seam. The native fixture
checks 214/246 total/core correct, separate parental flanks, restored low-quality
and GAF-absent molecules, and abstention on the contradictory crossing read.
The task-start binary reproduces 11 failures; the corrected native fixture
passes 666 assertions. In-memory cases cover both parities, sparse or missing
bridges, poor calibration, single-class calibration and sub-80% agreement.

Run from the repository root:

```bash
LD_LIBRARY_PATH=/usr/lib/x86_64-linux-gnu make gap-dev-check
LD_LIBRARY_PATH=/usr/lib/x86_64-linux-gnu make gap-owner-check GAP='56 Mb boundary'
LD_LIBRARY_PATH=/usr/lib/x86_64-linux-gnu make unit-tests check window-tests
bash evaluations/2026-10-06-next-split-repeat-block/replay.sh
```

Use the pysam Python environment to run `audit_block.py`, `audit_owner.py`,
`audit_preservation.py`, `audit_panel.py` and `audit_physical.py` in this
directory. `results.json`, `owner-results.json` and `full-preservation.json`
contain the full measurements. Ranking and screening use the task-start
full-chromosome output, keeping the target-selection evidence reproducible.
The HiPhase geometry audit uses its indexed BAM over the target and reuses
previously validated chromosome-wide tag orientations.

Final build, all units, 1,551 phasing-predicate assertions and HiFi/ONT goldens
pass. The complete gap suite passes 17,602 assertions across four cases. All
241 task-start native output labels / 127 independent replay requests are
unchanged against the final binary. The new cached regression takes 0.48 s.
The four replay-cache helper tests pass. Full validation identities and counts
are in `validation.json`; native preservation is in `panel-audit.json`.
