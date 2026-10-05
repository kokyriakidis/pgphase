# Physical graph SNP switch repair at 62 Mb

The graph arm now closes chr20 **62,623,253–62,642,316** in the owning
62,000,001–63,000,000 chunk and the full chromosome. The accepted pre-change
baseline is `test_data/tmp_gap_fix55/full-final`, binary SHA256
`bdce74cda1d0ec407fe7a0cfaf9e30d2c77a60b605d6f193ee92e2cb16503933`.

The right graph flank carries an internal switch between C>T SNPs at
62,718,395 and 62,722,021. All 38 independent primary BAM pairs at MAPQ/BQ
>=30 contradict its old relative orientation: 22 and 16 represent the two
haplotypes, none agree. Reversal log odds are 318.814 and the two-sided
unanimous-count probability is 7.28e-12. Neighboring SNP edges retain their
orientation with 49 and 56 unanimous pairs. `inspect_switch.py` and
`physical-evidence.json` reproduce these measurements from the baseline VCF
and original BAM; parental truth is not used to nominate or apply the repair.

The implementation repairs a single failed internal edge on a recovery seam's
right flank. It requires unique normalized SNP rows, primary Q30/MQ30 physical
calls, unanimous reversal, at least two molecules on each haplotype, count
p<=0.01 and quality-weighted wrong-parity probability <=0.001. It flips the
candidate suffix and records an exact graph-site edge certificate with its
relative REF/ALT gauge, comparing selected ALT presence rather than raw
allele IDs. The complete remaining graph SNP path must pass or the
repair is rolled back. Read HP gauges predominantly supported by phased heterozygous anchors in
that suffix flip with it; homozygous observations do not vote; prefix reads retain their assignments. The existing physical
seam bridge then validates and joins the two flanks. No chromosome coordinates
or parental labels occur in the production rule.

## Measurements

| Scope | Before | After |
|---|---:|---:|
| Owning-chunk truth-scored reads | 3,750 | 3,750 |
| Owning-chunk correct / discordant | 3,394 / 356 | 3,722 / 28 |
| Full-chromosome correct / discordant | 230,235 / 7,095 | 230,563 / 6,767 |
| Full-chromosome VCF blocks | 328 | 327 |
| Existing tracked gaps connected | 94/107 | 94/107 |

The full run has exactly one newly spanning old gap and no reopened tracked
gap. All 64,188 variant keys, every old phased key, genotype alleles, counts,
categories and other nonphase fields are preserved. The old right block
PS=62,637,077 deliberately changes its internal gauge to repair the switch;
all other old blocks retain uniform gauges. Span N50 remains 856,770 bp.
Of the 119 truth-scorable gap overlaps, 116 remain phased and 113 correct.
The correctly placed largest phase group increases from 53 to 95 reads
(separated fraction 0.445378→0.798319) because the core flanks join. Disjoint 10 kb flanks choose the same parent, 31/31 left and 36/36 right.

This is a correction of an existing internal switch, so it does change read
assignments: **329 previously discordant reads become correct; one previously
correct read becomes discordant**, a net reduction of 328 errors. The remaining
read is `m84031_231217_062403_s3/248908890/ccs`: its only high-quality phased
physical SNP is REF A at 62,645,168 (BQ40), matching its retained HP against the
repaired genotype. Its label is not special-cased using parental truth.
Twenty output-only rescued reads retain their exact HP and group membership;
their independent PS label changes from 1,062,637,077 to 1,062,432,427 as the
associated core label changes. There is no union with another rescue group.
`evaluate.py` checks all these properties explicitly, including the rescue
membership/HP bijection; `full-audit.json` and `full-parity.json` retain the
complete comparisons.

The panel adds the measured `spans=1`, 0.99 owning-read concordance and 0.79
separation floors; no existing window floor is changed. The graph TOTAL span
floor rises from 92 to 93 by adding this measured window. The dedicated owning
regression pins the corrected SNP relationship, parental orientation, at least
3,722 correct of 3,750 scored reads and at most 28 discordant. It passes 1,240
assertions; the old baseline fails four connection/accuracy/switch assertions.
The new panel window and required boundary/switch alleles also pass. The old
lone-boundary-SNP regression explicitly required these blocks to stay split.
It now permits a join only when the internal SNP relationship is repaired and
the whole owning chunk passes the same parental/error bounds; its 22 assertions
pass. It continues to reject an uncertified one-molecule join.

Build, unit tests, 1,536 predicate assertions, HiFi/ONT golden outputs and HiFi
thread determinism pass. The complete final native suite passes **9,757 assertions in 82 cases**,
including all 108 panel windows, in four disjoint batches against the frozen
`test_data/tmp_gap_fix58/pgphase-projected-final`. `native-summary.json`,
`native-cases.json` and the four `native-shard*.log` files record full coverage.
The final projected-gauge binary has exactly unchanged read tags and VCF rows
relative to the audited candidate (`projected-gauge-parity.json`).

## Rejected experiments

An isolated false terminal SNP at 57,481,589 has 63 high-quality REF A and no
ALT G observations, including paired support on both preceding haplotypes.
Clearing that anchor and admitting mixed complementary insertion/deletion
boundaries does not close any gap in a full-chromosome trial. Those changes
are reverted. A diagnostic 62 Mb candidate-only suffix flip closes the gap
but worsens owning read concordance to 71.12%; correcting the corresponding
read gauges is necessary. No diagnostic coordinate-specific rule is retained.
Trial snapshots and intermediate matrices remain under `test_data/tmp_gap_fix58/`.

Run the full comparison from the repository root with the reference, annotated
BAM, striped sites VCF and coordinate GAF using `collect-graph-variation -t 8`,
no region restriction, and phased VCF/BAM outputs. The frozen candidate,
full-run outputs and complete test workdirs remain under `test_data/tmp_gap_fix58/`.
Use the bench-phasers Python environment with pysam and
`LD_LIBRARY_PATH=/usr/lib/x86_64-linux-gnu` for the independent audit scripts.
