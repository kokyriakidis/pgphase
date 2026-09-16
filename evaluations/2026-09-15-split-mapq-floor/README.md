# Split MAPQ floor (`--min-assign-mapq`): implemented, validated on 1 of 3 gaps

`min_mapq` gated both jobs at once: reads below it never enter a phasing chunk,
so they were lost to candidate discovery and block linking as well as to
haplotype assignment. `--min-assign-mapq` separates the second job. A read
admitted above `-q` but below `--min-assign-mapq` now supplies allele evidence
and leaves the output without HP/PS
(`read_carries_phase_tags`, `src/collect_phase.cpp`, consulted in
`PhasedAlignmentWriter::write_chunks`). Equal floors reproduce the previous
behavior, which is the default: the default-arm `phased.bam` is byte-identical
to the pre-change build.

## The single gap the idea came from

`chr20:25,834,662-25,883,079`, 148 kb window. The in-pipeline result reproduced
the post-hoc estimate exactly on all three arms — reads / blocks / discordant /
N50:

| arm | in-pipeline | predicted |
|---|---|---|
| `-q 30` (default) | 395 / 2 / 0 / 166,590 | 395 / 2 / 0 / 166,590 |
| `-q 1 --min-assign-mapq 30` | 399 / 1 / 0 / 270,669 | 399 / 1 / 0 / 270,669 |
| `-q 1 --min-assign-mapq 5` | 505 / 1 / 0 / 275,646 | 505 / 1 / 0 / 275,646 |

## It does not generalize to the other two trustworthy gaps

Same arms on the other two gaps whose sub-floor reads were shown to carry
correct phase (`results.tsv`, `gate.tsv`):

| gap | arm | phased | blocks | discordant | concordant -> DISCORDANT | N50 |
|---|---|---:|---:|---:|---:|---:|
| 25,834,662-25,883,079 | default | 395 | 2 | 0 | — | 166,590 |
| | `-q 1` assign 30 | 399 | **1** | 0 | **0** | 270,669 |
| | `-q 1` assign 5 | 505 | **1** | 0 | **0** | 275,646 |
| 25,944,471-25,986,123 | default | 311 | 2 | 0 | — | 166,590 |
| | `-q 1` assign 30 | 311 | 2 | 6 | **6** | 166,590 |
| | `-q 1` assign 5 | 389 | 2 | 19 | **6** | 168,727 |
| 26,029,591-26,088,679 | default | 292 | 2 | 1 | — | 124,172 |
| | `-q 1` assign 30 | 266 | **1** | 6 | **6** | 300,159 |
| | `-q 1` assign 5 | 376 | **1** | 20 | **6** | 343,437 |

The regression gate fails on two of the three: 12 previously concordant reads
become discordant, and **the count is identical at assignment floor 30 and 5**.
That is the informative part — the corrupted reads are not the low-MAPQ reads
being tagged, since at assign 30 no read below MAPQ 30 carries a tag at all.
Opening *discovery and linking* to MAPQ 1 is what moves them: the sub-floor
reads' alleles and link votes change the solution for confidently mapped reads.
Coverage is not monotonic either — gap 26.03 Mb joins its flanks at assign 30
while losing 26 phased reads.

So competitor evidence that a gap's low-MAPQ reads are phaseable (96-100% in all
three of these gaps) does not imply our mechanism phases them correctly, and the
one gap the design was derived from was not a representative sample.

## Status

The option ships opt-in with defaults byte-identical, useful as a diagnostic
lever, and is **not** a default candidate. The next step is selective admission
rather than wholesale: gate the sub-floor evidence per site (the balance signal
from `../2026-09-15-mapq-starved-gaps`, median |VAF-0.5| 0.052 across the
trustworthy gaps against 0.115 across the noise-verdict ones) or down-weight sub-floor
link votes, then re-measure the gate on all three gaps before any whole-chr20 run.

## Files

- `run.sh` — three gaps x three arms, with per-gap truth subsetting by read name.
- `results.tsv` — coverage, blocks, accuracy and N50 per gap and arm.
- `gate.tsv` — read-level transitions against each gap's default arm.
