# Retained source reads in the overlapping 36.259 Mb block

HiPhase's rank-8 block is chr20:36,137,653–37,393,302 (1,255,650 bp).
The baseline is `dc64635`, full output `test_data/tmp_gap_fix86/final/0`.
Ranking only gaps between phase-set extents missed the separate BAM orphan
36,259,922–36,261,301 inside an already spanning graph phase set.
`screen.py` also examines adjacent phased-site PS transitions. Among the
remaining extent-uncovered HiPhase blocks, this is the largest with a
qualifying missed read connection; larger remaining seams have HiPhase
correctness below 80% of all original truth-scorable primary overlaps.

| Interval | Original reads | Before correct / core | After correct / core | HiPhase correct / core |
|---|---:|---:|---:|---:|
| 36,252,907–36,259,922 | 90 | 74 / 65 | 81 / 80 | 80 / 80 |
| 36,261,301–36,268,291 | 100 | 93 / 84 | 97 / 96 | 96 / 96 |

Correct fractions include unphased original reads: 90% and 97% after repair.
The dominant connected core now matches HiPhase at both intervals, and total
correct reads exceed HiPhase by one. The comparison uses the same original
alignment geometry, sequence, base qualities, mapping qualities and flags;
`audit_block.py` verifies these against HiPhase's BAM.

## Defect and repair

A later focused recovery replaced the earlier BAM source metadata while
retaining its variants and read observations. The retained D1 marker at
36,252,908 (VCF 36,252,907 CA>C) consequently lost its original shared-SNP
provenance. Its source cut leaves the marker: it does not invalidate the
already anchored suffix from the previous shared SNP at 36,247,421.

Retries now retain a complete old source cohort only when every exact member
is unique, still present and unclaimed by the new source. Source IDs and
candidate/read indices are remapped without changing the source allele gauge;
weak and quality cuts remain intact. Missing, refreshed or ambiguous members
reject retention of the whole original cohort.

A bounded read assignment uses the retained D1 evidence only before a separate
BAM orphan inside an already spanning core. It requires two earlier shared
clean SNPs with a consistent gauge, an uncut suffix from the nearest previous
shared SNP, the original exit cut, and a following same-core graph SNP within
20 kb. A verified MSA anchor within 32 bp vetoes the binary length call.
Original primary MAPQ30 reads need complete bounded repeat sequence and
length support; contrary clean phased SNP observations veto assignment.
This assigns reads to an existing core and does not certify a new variant join.

HiPhase does not create the same I5/I6 and I1 heterozygous repeat orphan: its
repeat calls are homozygous. It can use the D1 and flanking SNP evidence
without that artificial contrast. This repair retains pgphase's independently
verified source gauge and uses each original read's physical D1 contrast.
The orphan variant calls themselves remain separate. The target's leading
36,137,653–36,140,083 interval also remains open: HiPhase has only 40/82
correct original reads there. The entire HiPhase block is not declared closed.

## Preservation and validation

The native 36–37 Mb owner gains seven correct reads (five previously unphased,
two previously discordant), loses no previously correct or phased labels,
and preserves all variant keys, evidence and counts. Fifteen HP/PS labels
change, including eight previously correct rescue reads attached to the core.
One- and four-thread candidate TSVs, VCF records and BAM assignments agree.
A broader trial lost 11 correct labels elsewhere and was rejected; the final
orphan and neighboring-anchor gates avoid those unrelated changes.

The old binary fails 17 assertions in the new native regression; the repair
passes all 699. Both intervals are added to the committed panel, certified
manifest, HiPhase reference and expectations. The owning replay retains the
original phase-set context. Tests verify restored molecules, both parental
flanks, all-original denominators and connected core counts.

Full chr20 has 250 blocks, N50 955,496 bp and largest block 3,012,193 bp,
unchanged. Correct reads increase 230,894→230,901, discordant 6,644→6,642,
and unphased 19,074→19,069. No previously correct or phased read is lost,
all 82,287 variant keys/evidence and previous block extents survive, and no
read-status change occurs outside the target neighborhood. Thirty-two HP/PS
labels change across the chromosome: the 15 owner labels and 17 already
correct rescue labels at 56.671 Mb attached to their existing core after
source provenance is retained. All VCF records, including GT/PS, are unchanged.
Whole-target reads improve 4,841→4,848 correct and 4,744→4,759 correct in the
main core; HiPhase has 4,906 for both, leaving deficits of 58/147.

The full gap suite passes 17,715 assertions across four cases, and all four
cache-helper tests pass. The panel audit checks 244 native output labels
(including new baseline controls) against 128 current independent requests.
Thirteen labels change only their BAM read tags: no old correct/phased label
or VCF record is lost or changed. Build, units, 1,551 phasing-predicate
assertions and HiFi/ONT golden/determinism gates pass. The warm targeted
regression takes 0.66 seconds. Detailed checks are in `validation.json`.

## Reproduce

Run from the repository root, with pysam available for the audit scripts:

```bash
bash evaluations/2026-10-06-overlapping-source-block/replay.sh
LD_LIBRARY_PATH=/usr/lib/x86_64-linux-gnu make gap-owner-check GAP=36.259
LD_LIBRARY_PATH=/usr/lib/x86_64-linux-gnu make window-tests
python evaluations/2026-10-06-overlapping-source-block/audit_block.py
python evaluations/2026-10-06-overlapping-source-block/audit_preservation.py
python evaluations/2026-10-06-overlapping-source-block/audit_owner.py
python evaluations/2026-10-06-overlapping-source-block/audit_panel.py
```

Original inputs and the cached competitor are unchanged. Full target results,
N50 and whole-block deficits are in `results.json`; preservation, native-panel
comparison and required checks are recorded alongside it.
