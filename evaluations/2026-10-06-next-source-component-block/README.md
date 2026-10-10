# Retained source components at the 46.403 Mb overlapping orphan

The refreshed scan of the 50 largest HiPhase blocks finds the next qualifying
miss in rank 30: chr20:45,866,904–46,636,707, 769,804 bp. The starting executable
is `098e6b4895707ecc6a7d2af558b1c4aa9688741f9f742b807e339851596d117e`, with
full output `test_data/tmp_gap_fix88/final/0`. Pgphase already covers the target's
coordinate span; the missing product is correctly connected read coverage.
All larger remaining seams either meet HiPhase's total/core correct counts or
have HiPhase below 80% correctness on original primary overlaps. Discovery
includes overlapping PS transitions (`ranked-blocks.json`, `screened-seams.json`).

| Connection | Original reads | Before correct / core | After correct / core | HiPhase correct / core |
|---|---:|---:|---:|---:|
| 46,389,864–46,402,979 | 118 | 104 / 88 | 107 / 107 | 106 / 106 |
| 46,402,979–46,405,048 | 75 | 70 / 62 | 72 / 72 | 72 / 72 |

Correct fractions include all unphased original primary overlaps: 90.68% and
96.00%. The first retains two discordant and nine unphased reads; the second
retains one discordant and two unphased reads. Both total and dominant connected
correctness meet or exceed HiPhase. `audit_block.py` verifies identical read
names, alignment geometry, flags, mapping quality, sequences and base qualities
against HiPhase's BAM.

## Defect and repair

HiPhase phases a 14 bp deletion at 46,388,261 and nearby one-base repeat edits.
Pgphase retains the same deletion as an alignment-verified graph catalog row,
including the original MSA source edit and shared-SNP gauge. The single-base
source deletion at 46,402,980 remains an overlapping orphan inside the main
core; the source cut is at 46,389,864. The two source components on its flanks
already have multiple consistently oriented clean SNPs in that core.

The old read promotion requires a new shared-deletion bridge certificate and
an entirely uncut source, and its coverage guard vetoes the whole core whenever
rescue coverage exceeds core coverage anywhere in the chunk. Those requirements
apply to a block transfer but exclude independent read evidence inside an
already spanning block. Original complete-sequence controls find 30 maternal
REF / 32 paternal ALT calls at the 14 bp deletion, with independent upstream
SNP calls 17 maternal / 22 paternal. The one-base repeat calls are noisier;
using a single repeat's allele to reorient every read is rejected.

The repair verifies the marker's original uncut component on at least two
shared clean graph SNP loci. Every retained source anchor in that component
must keep a single consistent current gauge. Weak and quality cuts delimit
its evidence. A nearby BAM-only island must already have clean SNPs in the
same core on both sides. This certifies local reads without a new block join;
other ordinary core assignments remain protected.

A shared pure deletion supplies each read's physical REF/ALT call and can
restore abstentions or correct a wrong rescue label. Insertions promote only
existing agreeing rescues, using the complete repeat, 16 bp external flanks,
bounded length slippage and agreeing unique sequence distance. An insertion
cannot overwrite an established core label, preserving stronger deletion
assignments when a noisy repeat disagrees. Calls require primary nonduplicate
alignments, MAPQ at least 20, known base qualities and <=5% combined mapping
and physical error. Contrary clean phased SNPs veto assignment. Production
uses no parental truth or HiPhase calls.

The native 46–47 Mb owner improves from 3,905 correct / 17 discordant / 120
unphased to 3,908 / 16 / 118. Twenty HP/PS labels change: two abstentions and one
wrong label become correct, 16 already correct rescues enter the core, and one
old discordant rescue keeps its haplotype while entering the core. All old
correct/phased reads survive. Every candidate TSV byte and complete VCF record
is identical, and one-/four-thread candidates, VCF records and BAM assignments
match (`owner-results.json`).

## Regression and reproduction

Both connections enter the committed panel, certified manifest, identical-input
HiPhase reference, replay map, required sites and measured expectations. Their
native owner regression checks original denominators, total/core floors,
parental flanks, restored molecules and unchanged whole-owner accuracy.
The starting executable fails 15 of 720 assertions; the fix passes all 724.
The warm focused run takes 0.76 seconds. Sixteen new in-memory checks exercise
component boundaries, quality cuts, source provenance, distinct SNP loci,
contrary gauges and coherent reversals.

The first broad preservation audit rejected extending the local exception to
unshared short source deletions: a 53 Mb replay moved 13 rescue labels and left
six labels with the opposite parental majority, causing two unchanged reads
to lose their measured orientation. Unshared deletions now retain their old
whole-cohort guards; the local deletion exception is restricted to the shared
catalog deletion. The extra unit checks explicitly require this distinction.
The original acceptance thresholds and all existing expectations are retained.

The final executable is
`12beb76e831986814ebb50a6d55eacbfd024f517224f7ca11039c80d06272058`.
Build, all units, 1,551 phasing-predicate assertions in 47 cases, HiFi/ONT
goldens and HiFi thread determinism pass. The complete gap suite passes
17,954 assertions in five cases; four cache-helper tests pass. The native-panel
audit covers 248 output labels, 128 current independent requests and 129
before/after comparison pairs. Seven output labels change read tags, with no
old correct/phased assignment lost and every complete VCF record preserved.
The 53 Mb control retains every HP/PS label and VCF record exactly.

The final full chr20 replay has 230,925 correct / 6,620 discordant / 19,067
unphased primary reads, versus 230,922 / 6,621 / 19,069 before. All 20 changed
tags overlap the target neighborhood. Every prior correct/phased read, old
block extent, all 64,483 complete VCF records and every byte of the 82,287-site
candidate TSV survive. N50 remains **955,496 bp** and the largest block remains
3,012,193 bp because this repair connects reads within an existing block.

These two connections meet the acceptance gates, but the **whole 769,804 bp
HiPhase target still has a read deficit**: pgphase has 3,038 total / 3,020
connected-core correct reads, versus HiPhase's 3,042 / 3,042 on 3,129 original
primary overlaps. The remaining deficits are four overall and 22 in the core;
whole-block parity is not claimed. `after-screened-seams.json` records the
refreshed transition scan. `validation.json`, `results.json`,
`full-preservation.json` and `panel-audit.json` identify the final results.

```bash
LD_LIBRARY_PATH=/usr/lib/x86_64-linux-gnu make gap-dev-check
LD_LIBRARY_PATH=/usr/lib/x86_64-linux-gnu make gap-owner-check GAP=46.403
LD_LIBRARY_PATH=/usr/lib/x86_64-linux-gnu make window-tests
bash evaluations/2026-10-06-next-source-component-block/replay.sh
python evaluations/2026-10-06-next-source-component-block/audit_block.py
python evaluations/2026-10-06-next-source-component-block/audit_preservation.py
python evaluations/2026-10-06-next-source-component-block/audit_owner.py
python evaluations/2026-10-06-next-source-component-block/audit_panel.py
```

Run from the repository root with pysam available for Python audits. Full
chromosome, whole-target and panel preservation results are in the JSON reports.
