# Close the 8.978–9.014 Mb gap across owning chunks

Connect chr20:8,977,829–9,014,032 across the 9 Mb chunk boundary. The accepted
baseline is `test_data/tmp_gap_fix54/full-gated`, binary SHA256
`3a227a265b28ca21b15e4e59e272f0dc2cc06e7feb2046d50108365e481152c5`.
The retained binary is `test_data/tmp_gap_fix55/pgphase-final`, SHA256
`bdce74cda1d0ec407fe7a0cfaf9e30d2c77a60b605d6f193ee92e2cb16503933`.

## Certificate and implementation

The bounded 8,927,829–9,064,032 replay contains two complementary BAM insertions
at VCF position 8,998,972: `T>TAC` and `T>TACAC`. Both are MSA- and
alignment-verified, opposite phased heterozygotes in one independently supported
BAM source run. Each distinct ALT has one exact physical molecule spanning the
left `T>A` graph SNP. Both reads have MAPQ 60 and SNP Q40; the insertion and
aligned-flank minimum qualities are Q35 for AC and Q40 for ACAC. The two
molecules unanimously orient the relation, with log odds 16.286107 and
wrong-parity probability below 0.001. Literal REF of either insertion never
substitutes for the other insertion ALT.

Both current graph SNP paths and the original BAM source path must pass their
existing checks. Detection is confined to an already targeted recovery window.
The certificate retains every original phased anchor and is applied after
read rescue. Raw-key lookup must ignore unphased multiallelic shadow rows: at
9,048,832 only the third of five rows sharing one graph key is phased. The
uniform live phase-set and allele-gauge checks remain mandatory.

For the owning-chunk union, require both exact boundary SNPs in the same replay
phase set and at least two exact shared clean SNPs per original owning block,
with uniform replay orientation. Save both blocks' original anchors and their
owning chunks, then apply the union after independent read rescue finishes.
Private replay genotypes, depths, observations and candidate rows are not
transferred. Independently rescued HP/PS labels remain unchanged.

## Exact audits

| Scope | Correct | Discordant | VCF blocks before → after | Read PS before → after | Span N50 before → after |
|---|---:|---:|---:|---:|---:|
| Bounded seam | 530 | 7 | 2 → 1 | 3 → 2 | 60,405 → 112,849 |
| Paired owning 8–10 Mb | 8,517 | 8 | 6 → 5 | 10 → 9 | 985,788 → 1,053,649 |
| Full chr20 | 230,235 | 7,095 | 329 → 328 | 651 → 650 | 856,770 → 856,770 |

Every individual parental truth status is unchanged in all three scopes. The
full chromosome retains all 256,610 output reads, 237,330 scored reads and
64,188 variant keys. Exactly this one old gap closes; no tracked gap reopens.
94/107 committed panel gaps now connect. All old phased sites, genotype
alleles and non-phase fields remain. Every absorbed old block moves uniformly.
The full output has 5,634 core tag relabels and 1,272 VCF phase changes, with
zero output-only HP/PS changes. Paired output has 4,229 core tag relabels and
988 VCF phase changes, also with zero output-only changes.

The gap has 223 truth-scorable input reads. 180 are phased and 177 correct
before and after; three phase groups become two, with the largest correctly
placed group increasing from 64 to 126 reads (126/223 = 0.565022). The other
independent rescue group remains intact. Disjoint 10 kb flanks select the same
parent, 35/35 left and 45/45 right. The new panel row uses measured `spans=1`
and a 0.56 separated-fraction floor; no existing expectation is relaxed.

The panel comparator is native HiPhase on the same BAM with contig names
normalized. `competitor.json` checks identical sequences, qualities, CIGARs,
coordinates and MAPQ for all 223 overlapping truth-scorable reads: 180 phased,
179 correct in one group (0.802691 separated).

## Verification and reproduction

Build, standalone unit tests, all 1,536 predicate assertions, BAM HiFi/ONT
goldens and thread determinism pass. The focused owning-chunk regression passes
24 assertions; the accepted baseline fails its connection/orientation checks.
The complete final native suite passes 8,234 assertions in all 81 cases,
using 112 fresh pipeline replays. Results are recorded in `native-tests.log`.

`{bounded,paired,full}-{parity,audit}.json` preserve the exact output and
parental comparisons. `physical-evidence.json` audits both ALT molecules and
the multiallelic shadow directly from the input BAM and final bounded matrix.
`manifest.json` records binary, source, output and input provenance.

Run `reproduce.sh BEFORE_BINARY AFTER_BINARY OUTPUT_DIR` from the repository
with `PYTHON` pointing to a pysam-enabled Python. It performs fresh bounded,
paired and full replays and all audits. The final native suite uses the frozen
binary through `PGPHASE_BIN`; its working directory is
`test_data/tmp_gap_fix55/native-counts-final`.
