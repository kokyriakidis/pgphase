# Calibrated complementary deletion reads at the 4.637 Mb orphan

The refreshed audit of the 50 largest HiPhase blocks identifies rank 34,
chr20:4,570,122–5,196,043 (625,922 bp), as the next qualifying missed connection.
The starting binary is
`12beb76e831986814ebb50a6d55eacbfd024f517224f7ca11039c80d06272058`, with full
output `test_data/tmp_gap_fix89/final/0`. Pgphase already covers the coordinates;
its missing product is correctly connected original-read coverage. Discovery
includes overlapping phase-set transitions (`ranked-blocks.json`,
`screened-seams.json`). Larger remaining transitions either meet HiPhase's
counts or have HiPhase below 80% correctness on original primary overlaps.

| Connection | Original primary reads | Before correct / core | After correct / core | HiPhase correct / core |
|---|---:|---:|---:|---:|
| 4,625,182–4,637,298 | 127 | 95 / 95 | 113 / 113 | 113 / 113 |
| 4,637,298–4,642,430 | 87 | 69 / 69 | 76 / 76 | 76 / 76 |

The final correct fractions are 88.98% and 87.36%, including 14 and 11
unphased original reads respectively. Neither repaired connection has a
pgphase discordant read. HiPhase has one discordant read at the first.
The comparison uses identical primary alignments, read geometry, flags,
mapping quality, sequences and base qualities (`audit_block.py`).

## Why HiPhase phases these reads

HiPhase uses a heterozygous 4 bp / 8 bp AGGA-repeat deletion contrast at
4,625,192. Pgphase retains equivalent normalized source deletions at
4,625,183, both MSA- and alignment-verified, with complementary gauges already
in the connected core. Its source-read promotion treats each deletion as a
separate binary REF/ALT marker and rejects a nearby compound contrast.
Existing complete-deletion repeat calling also only recognizes one- or
two-base motifs. Eighteen reads inside the SNP desert therefore stay unphased
although their complete original sequence distinguishes the retained alleles.

Comparing the whole AGGA repeat and 16 bp external flanks yields 10 paternal
short-deletion and eight maternal long-deletion calls among those abstentions.
All 18 agree with parental truth. Independent original Q30 SNP calls at
4,618,685 calibrate the retained gauge in deterministic molecule halves:
`[[0,8],[11,0]]` and `[[1,7],[8,0]]`, with rows the deletion class and columns
SNP haplotype. Both halves observe both classes and support the same gauge.
Production never uses parental truth or HiPhase's variants/tags.

The nearby 2 bp / 3 bp homopolymer deletion at 4,637,299 is much less reliable:
complete sequence calls mix both parents in each length class. That orphan
has no clean SNP in its own phase set, so it cannot certify new reads. Using
its noisier length to overwrite existing read labels is unnecessary.

## Repair and regressions

Complete-repeat deletion calling now uses the existing bounded fundamental
motif detector (up to eight bases), with both net deletion length and unique
whole-sequence distance required. New read filling requires complementary
pure source deletions, verified in the same existing core, and independent
diploid calibration against the nearest clean graph SNP within 20 kb. Both
molecule halves must agree with the retained gauge, pass association and the
20% gauge-error bound. Calibration calls have <=1% combined physical/mapping/
base error. Assignment requires primary nonduplicate known MAPQ30 alignments,
known external base qualities and <=5% physical plus mapping error. Contrary
clean phased SNP observations veto assignment. Only reads without any core
or rescue haplotype are eligible.

The 4–5 Mb native owner changes from 3,844 correct / 26 discordant / 256
unphased to 3,862 / 26 / 238. Exactly 18 abstentions become correct core reads.
Every previous phased HP/PS label is identical, all old correct reads remain
correct, every candidate TSV byte and complete VCF record is unchanged.
One- and four-thread owner outputs are audited in `owner-results.json`.

Both connections enter the panel, certified manifest, identical-input HiPhase
reference, owner replay map, required sites and measured expectations. The
native regression checks exact original denominators, HiPhase total/core
floors, four restored molecules, independent parental SNP flanks and native
accuracy. The starting binary fails 15 of 734 assertions; the repaired fixture
passes all 738 and takes 0.88 seconds warm. Four additional in-memory checks
cover the four-base motif, measured diploid cohorts, contradictory gauges and
noisy homopolymer abstention. Existing thresholds and floors are retained.

The final binary is
`f836c348e7cb5e7f7970e8cf08d506f8f1e3a064eda15ec171e502988a4aba9a`.
Build, all units, 1,551 predicates in 47 cases, HiFi/ONT golden output and
HiFi determinism gates pass. The complete gap suite passes 18,067 assertions
in five cases; all four cache-helper tests pass.
The native-panel audit covers 250 output labels, 128 current independent
requests and 129 before/after comparison pairs. Thirty-seven labels improve
read tags, with every old correct/phased read and complete VCF record retained
(`panel-audit.json`, `validation.json`).

The final full chr20 run has 230,965 correct / 6,620 discordant / 19,027
unphased primary reads, versus 230,925 / 6,620 / 19,067 before. All 40 changed
labels were unphased and become correct core labels. Eighteen overlap the
target neighborhood; the other 22 lie in the 3, 11, 19, 41 and 48 Mb chunks.
Every old phased HP/PS label and correct assignment is preserved, as are all
old block extents, all 64,483 complete VCF records and every byte of the
82,287-site candidate TSV. `tag-changes.json` locates all 40 assignments in
their original alignments. N50 remains **955,496 bp** and the largest block
remains 3,012,193 bp; the repair assigns reads within existing blocks.

Whole-block read parity remains open. On 2,622 original primary overlaps,
pgphase has 2,511 total / 2,426 dominant-core correct reads versus HiPhase's
2,515 / 2,515. The remaining deficits are four overall and 89 in the core.
The separate 4,874,129–4,884,130 transition retains 54 total / 31 core correct
reads versus HiPhase's 53 / 53 on 64 original overlaps; this rescue-to-core
connection is not changed or certified by the present repair.
`after-screened-seams.json` records the remaining transitions.

```bash
LD_LIBRARY_PATH=/usr/lib/x86_64-linux-gnu make gap-dev-check
LD_LIBRARY_PATH=/usr/lib/x86_64-linux-gnu make gap-owner-check GAP=4.637
LD_LIBRARY_PATH=/usr/lib/x86_64-linux-gnu make window-tests
bash evaluations/2026-10-06-complementary-deletion-reads/replay.sh
python evaluations/2026-10-06-complementary-deletion-reads/audit_block.py
python evaluations/2026-10-06-complementary-deletion-reads/audit_owner.py
python evaluations/2026-10-06-complementary-deletion-reads/audit_preservation.py
python evaluations/2026-10-06-complementary-deletion-reads/audit_panel.py
```

Run from the repository root with pysam available. Whole-target and final
preservation measurements belong in the final JSON reports; the two accepted
connections do not by themselves establish whole-block read parity.
