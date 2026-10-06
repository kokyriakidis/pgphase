# Recover the seventh HiPhase block’s terminal deletion

The selected boundary is **40,633,644–40,636,353 (2,709 bp)** at the end of
HiPhase’s seventh-largest chr20 block, **39,273,205–40,636,353 (1,363,149 bp)**.
The fresh ranking uses the preceding cut-free insertion-prefix output. All
larger uncovered boundaries fail HiPhase’s 80% original-overlap rule, including
terminal differences; see `selection.json` and `screened-seams.json`.

## Cause and repair

The graph catalog describes `TG>TT,GT` at 40,636,353, selecting the padded
`G>T` SNP. Its graph classes contain 9 REF and 28 ALT observations, but original
primary alignments contain **no G**: they substitute T or delete G. Validation
previously admitted only literal one-base, single-ALT catalog SNPs, so this
padded selection with an unused ALT escaped physical validation. The false
singleton overrode the upstream haplotypes and misassigned 13 gap reads.

BAM MSA recalled the substitution as homozygous ALT and did not retain the
terminal deletion. Simply dropping the graph SNP leaves that deletion
unsupported and lets a later source block capture some tail reads. That trial
failed both the connected-core floor and preservation checks and was rejected.

Validation now inspects the minimal selected SNP when only one candidate
survives from its catalog site. A padded substitution/deletion contrast needs
zero physical REF/other, corrected 1% significance, at least two deletions and
20% deletion support. It is retired **after** the initial solve so surviving
anchors keep their gauges. Existing unpadded exclusion gates remain unchanged.

After recovery and rescue, Q30 primary deletion/SNP pairs calibrate the terminal
deletion against the nearest upstream MSA- and alignment-verified SNP. At least
two pairs must agree, with joint wrong-gauge bound <=0.001; contrary pairs veto
reconstruction. The compound genotype retains substitution and deletion as
separate nonreference alternatives. Its two phased VCF rows are complementary,
with depths and strand counts derived from exact original primary calls.
Original graph/BAM observations remain available; primary calls at the invalid
SNP are masked rather than reused with the new allele IDs.

Assignments use the calibrated deletion or the upstream physical SNP with
MAPQ/BQ >=20 and summed error <=0.01. All local MSA observations must agree.
An uncertain or unobserved terminal base does not cancel a valid upstream call;
reads without either informative witness keep their existing assignments.
Parental truth is used only for evaluation.

## Native results

All denominators include unphased original primary truth-scorable overlaps.
Rescue phase sets do not count toward connected core coverage.

| Tool | Correct | Discordant | Unphased | Scorable | Dominant correct core |
|---|---:|---:|---:|---:|---:|
| Before | 32 | 13 | 0 | 45 | 23 |
| Fixed | **45** | **0** | **0** | **45** | **41** |
| HiPhase | 40 | 0 | 5 | 45 | 40 |

The owner has 4,061 primary reads: **3,848→3,863 correct**, **86→73
incorrect**, and **127→125 unphased**. Every previously correct and tagged
read is preserved. The changed evidence is the corrected terminal `G>T` and
new `TG>T` deletion; unrelated variant alleles/counts/filters are unchanged.
One- and four-thread candidates, VCF records and BAM assignments are identical
(`owner-results.json`). Disjoint upstream and terminal-only parental cohorts
verify orientation independently of the gap overlap counts.

The new panel row, certified manifest, exact `spans=1` expectation and native
40–41 Mb owning regression enforce total correctness and core parity. The
saved old binary fails 11 assertions in that regression. The optimized build
has zero warnings; units, 1,551 predicate assertions and HiFi/ONT goldens and
determinism pass. HiPhase’s new row is independently regenerated on identical
alignments by `make gap-benchmark`.

## Full chromosome and preservation

The main block now spans **39,273,205–40,636,354 (1,363,150 bp)**, covering
HiPhase’s entire seventh span and extending one base beyond its endpoint.
No other HiPhase block becomes newly covered. **N50 stays 944,186 bp**, block
count stays 258, and the largest block stays 3,012,193 bp.

Across all 256,610 primary reads, correct **230,651→230,666**, discordant
**6,705→6,692**, unphased **19,254→19,252**. Every previous correct/tagged
assignment and block extent is preserved (`full-preservation.json`).

Full-block read parity is separate from terminal-gap acceptance: over all
5,763 reads overlapping this 1.36 Mb block, pgphase has **5,667 correct/core
5,606**, versus HiPhase **5,671 correct/core 5,671**. The remaining deficit
is four total correct and 65 connected correct reads. Covering the complete
span does not claim those separate read deficits are fixed.

One other phased catalog description, **60,171,648 T>C**, is withdrawn by
the same physical validation: original MAPQ30/BQ20 alignments have **37 C,
34 deletion, zero T**. Its deletion spans neighboring bases, so exact
one-base reconstruction abstains. Its read assignments and block extents are
preserved. Apart from the corrected terminal G>T, all remaining old variant
alleles/counts/filters are unchanged. Adding TG>T and withdrawing this row
leave the full VCF at 64,483 records.

The full suite passes **16,185 assertions in four cases**, plus four cache
helper tests. `panel-audit.json` covers **224 previous labels/129 current
cache requests** across system/local htslib runtimes. No correct read is lost;
only two 60 Mb labels lose the invalid T>C description, with identical read
tags/counts. The saved new owner state reruns 583 assertions in **0.47 s**.

Final optimized binary SHA256:
`bf4ab0d3c70cb98483199663c2308bedbf35223d43fa3cfb3adf4ba4286a4a35`.

## Reproduction

```sh
make -j8
make gap-dev-check
make gap-owner-check GAP=terminal-deletion
make unit-tests predicate-tests check
LD_LIBRARY_PATH=/usr/lib/x86_64-linux-gnu make window-tests
bash evaluations/2026-10-06-terminal-source-extension/replay.sh
/home/kokyriakidis/micromamba/envs/bench-phasers/bin/python evaluations/2026-10-06-terminal-source-extension/audit_block.py
/home/kokyriakidis/micromamba/envs/bench-phasers/bin/python evaluations/2026-10-06-terminal-source-extension/audit_owner.py
/home/kokyriakidis/micromamba/envs/bench-phasers/bin/python evaluations/2026-10-06-terminal-source-extension/audit_preservation.py
```
