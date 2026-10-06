# Shared graph/BAM source-run SNP bridge

The genuine open gap at chr20:49,720,311–49,742,973 has a callable physical
SNP pair, but neither flank qualified under the existing graph-only path
rules. The left block contains shared graph rows and the right block is a
local run before an original source cut. The existing
`bam_source_run_supported(..., true)` validator certifies both: every phased
heterozygous anchor has exactly one adoptable claim in a single independent
BAM source, its current gauge is consistent, and no weak or quality cut
crosses the local anchor interval. Homozygotes do not supply orientation.

The original physical evidence is one primary MAPQ-60 molecule,
`m84031_231217_062403_s3/200282771/ccs`, calling ALT at both SNPs at BQ40.
Its two-base/two-mapping error bound is 0.000202, and the same-gauge log odds
are −8.50704, exceeding the existing 0.001 wrong-parity threshold. Both
endpoints have independently called biallelic coverage (37/37 and 30/34).
No confidence threshold or candidate classification changes.
See `physical-pairs.json`. Native HiPhase gives both exact SNP descriptions
GT 0|1 / PS 49716722 (`hiphase-endpoints.json`); pgphase retains its existing
GT 1|0 gauge and joins into PS 49700407.

The new route accepts two complete source runs, including shared graph rows,
when the ordinary path check fails. It records every anchor in the existing
deferred-bridge representation. Application after rescue checks the current
block gauges again and merges only core variants/reads. It does not admit
previously unphased reads or use parental truth in production.

## Matched HiPhase comparison

All 162 primary input overlaps have parental truth. The comparator is native
HiPhase/DV at `test_data/tmp_gap_fix48/competitor/hiphase_dv/phased.bam`.
Its read names, coordinates, CIGAR and sequences match all 162 original BAM
alignments. Unphased reads remain in the accuracy denominator; rescue PS
labels are excluded from the dominant core count.

| Output | Correct | Wrong | Unphased | Correct in one core |
|---|---:|---:|---:|---:|
| Before | 137 | 0 | 25 | 74 |
| After | 137 | 0 | 25 | 137 |
| HiPhase | 137 | 0 | 25 | 137 |

Thus the closure matches HiPhase at 137/162 (84.57%) total and core accuracy.
The owning 49–50 Mb audit preserves all 3,489 scored read outcomes, all
variant keys, genotype alleles and nonphase fields, all old block gauges,
and every output-only rescue tag. Disjoint flanks retain the same parental
orientation (93/0 and 72/7 votes). See `owner-audit.json` and
`owner-hiphase-comparison.json`.

An initial eager implementation merged the core blocks before rescue. It
combined independently inferred rescue cohorts, reclassifying 21 formerly
correct and 12 formerly discordant owning-region reads. That implementation
was rejected. Deferring the union preserves every rescue cohort and all
previous read outcomes.

## Reproduction

Run from the repository root with the bench-phasers Python environment and
`LD_LIBRARY_PATH=/usr/lib/x86_64-linux-gnu`. Inputs and frozen binaries are
local evaluation artifacts under `test_data/tmp_gap_fix63/`.

```bash
python evaluations/2026-10-05-shared-source-snp-bridge/evaluate_gap.py \
  --before test_data/tmp_gap_fix62/full-certified \
  --after test_data/tmp_gap_fix63/full-deferred-source-runs \
  --output evaluations/2026-10-05-shared-source-snp-bridge/full-audit.json
python evaluations/2026-10-05-shared-source-snp-bridge/compare_hiphase.py
```

The new owning-region native case fails on the frozen preceding binary:
span=false, separated=74/162, and unequal endpoint phase sets. It passes
with the new binary (698 assertions). The committed panel adds this gap
with `spans=1`, a 0.84 separation floor and both endpoint witnesses; owning
context is fixed at 49,000,001–50,000,000. Existing expectations and required
sites are retained.

The full chromosome audit closes exactly this gap (98/111 tracked spans),
with no read outcome change: 230,566 correct and 6,768 discordant remain
unchanged across 237,334 scored reads. All rescue tags are byte-for-byte
unchanged. Only 133 core read phase labels and six variant phase labels
change; all 64,188 variant keys, genotype alleles and nonphase fields are
preserved. VCF blocks fall 325→324 and read phase sets 645→644; N50 remains
856,770. See `full-audit.json`, `full-parity.json` and
`hiphase-comparison.json`. Frozen binary SHA256:
`6824c6287134b7178f12e6690bcfca026b0a32850a861b36c041842e0144e29e`.

Final required checks pass: build without warnings, all unit tests, 47 phase
predicate cases / 1,536 assertions, deterministic HiFi/ONT validation gates,
and the complete native suite of 85 cases / 11,497 assertions over all 111
panel windows. The four native shards use the frozen final production and
test binaries; case coverage is complete with no duplicates. See
`validation.txt`, `native-validation.json` and `manifest.json`.
