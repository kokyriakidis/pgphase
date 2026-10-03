# Complementary insertion observations in gap recovery

## Reproduced defect and repair

Two separate MSA rows describe the four-T and eight-T haplotypes. An exact
read matching the eight-T consensus previously became unknown on the four-T
row. That row uses zero for absence of its selected ALT, so the correct
projection is four-T `{1,0}` and eight-T `{0,1}` across the two rows.

The caller now accepts an exact other-consensus match only when both fixed
consensuses cover the same reference footprint, their query flanks agree, and
exactly one carries the selected ALT. Both composed read alignments must call
the same allele with supported flanks. Incomplete or conflicting alignments,
ambiguous contexts and unrelated third insertions abstain. Candidate rows
remain separate; neither a new alignment nor an allele merge is introduced.
Production code receives no truth, competitor calls or fixture coordinates.

Restoring exact recognition also lets the existing one-error caller use two
valid, distinct consensus alleles. That caller is unchanged: common flanks,
a strictly closer consensus within one edit, and agreement between both
independently composed read paths remain required. A seven-T read can therefore
support the eight-T haplotype while a six-T read remains ambiguous. This is
existing sequence comparison, not a new alignment or guessed repeat allele.

The committed synthetic regression covers exact complementary calls,
unavailable or incomplete consensuses, unequal flanks, third alleles,
conflicting read paths, separate row identities and rebuilt observation counts.
Seven-T reads use the existing one-error fallback; six-T reads stay unknown.
The pre-fix linked test failed five checks (`representation-before.log`).

## Rejected recovery changes

Missing callable pairs can hide an unsupported edge inside one nominal BAM
phase set: the existing internal-conflict detector cannot count absent pairs.
An experimental trigger used the existing MSA dropout statistic to request a
focused fresh-observation retry. Simply adding this retry produced a wrong
24.12 Mb flank relation. Adding a joint exact MEC solve against complete clean
SNP flanks rejected that local trial, but the broader chromosome experiment
still regressed existing assignments.

The first broad joint guard changed 57 previously correct reads to incorrect
and lost 301 variant keys. A narrower dropout-only version preserved all old
keys but changed 979 previously correct reads to incorrect, largely through
an incorrect block connection around 17.8 Mb. Neither trigger nor new joint
wrapper remains in production. Established retry admission and stitch rules
are retained.

An exact-only restriction on the complementary insertion caller phased 21
more reads but lost 50 truth-correct assignments and increased discordance by
71, with two fewer VCF blocks. It discarded usable error-tolerant evidence and
was removed. Its results are in `rejected-exact-only-parity.json`. The broad
and dropout-only joint trials are in the other two `rejected-*-parity.json`
files. None is an accepted benchmark result.

## Targeted validation

The owning 24–25 Mb replay retains all 1,234 variant keys and every read tag:
3,895 phased/scored, 3,893 correct, two discordant. In-gap truth-scored reads
remain 60/60 correct. The 24.121713–24.131707 Mb gap remains open. Exact
insertion projection repairs one demonstrable information loss; it does not
make the current boundary matrix a certified bridge. The corresponding
7.264321–7.280346 Mb gap also remains open.

A permanent owning-chunk regression retains both insertion descriptions,
checks opposite insertion orientations when co-phased, rejects a joined SNP
flank with the wrong parental relation, and protects owning/local read scores.
The two owning-chunk tests protect representation and orientation independently of read-placement accuracy.


## Read labels after the final stitch

A second reproduced problem is an inherited read HP that contradicts independent
physical SNPs in its final PS. The new graph-recovery post-stitch pass requires
known MAPQ30, Q30 single-base clean heterozygotes, agreement between physical BAM
and primary observations, and at least 100 bases between witnesses. Any eligible
SNP disagreement vetoes correction; co-located rows count once. It changes only
an already assigned read's HP and scoring counters, never candidate genotypes,
PS labels, joins or unassigned reads. Both graph worker paths call it.

Synthetic adapter tests cover the repair, fixed candidate gauge, low/unknown
MAPQ, low base quality, physical disagreement, channel disagreement, duplicate
coordinates, adjacent witnesses, noisy SNPs, another PS and unassigned reads.
Compared with insertion recall alone, 36 previously incorrect reads become
correct and one previously correct read becomes incorrect. The net gain is 35
correct reads, with exactly identical VCF rows and phased-read counts. A
primary-graph-only version lacks the physical-quality certificate and undoes two
correct assignments in the 3 Mb replay; it was rejected. The old file names
`rejected-snp-correction-4.json` and `rejected-primary-markers-3.json` describe
intermediate probes, not accepted chromosome results.

## Accepted chromosome comparison

All runs use the same chr20 reference, graph catalog, GAF, BAM and parental truth
map. Truth is loaded by evaluation only. The starting binary is
`cb3b4ad62fb671abb1f375e02a48a44b1eac2553158b7ab39bc6c9755d94fc06`;
the final binary is
`2f0513f66ddcbc2a06f0a1b37d3674b9acb4f6df6c735540c7c8df8af0a3a03e`.
Full native outputs are `test_data/tmp_gap_next27/accepted/` and
`test_data/tmp_gap_next27/insertion-bam-refresh/`.

| Metric | Starting version | Both fixes |
|---|---:|---:|
| Phased, truth-scored reads | 237,134 | 237,172 |
| Truth-correct reads | 229,981 | 230,022 |
| Discordant reads | 7,153 | 7,150 |
| Read concordance | 96.983562% | 96.985310% |
| Read phase sets | 665 | 661 |
| VCF variant keys | 63,492 | 63,630 |
| VCF phase blocks | 333 | 333 |
| Variant block span N50 | 756,878 bp | 774,189 bp |

This is a net improvement, not unchanged accuracy for every read: 54 previously
correct reads become incorrect, 67 incorrect become correct, 45 unphased become
correct, ten unphased become incorrect and 17 correct become unphased. All old
variant keys remain. All 57,769 previously phased SNP keys retain a constant
gauge within each old SNP block, with no old block splitting or internal
reorientation. Three pairs of old SNP blocks merge; this does not imply that
every intervening indel boundary joins. The SNP refresh changes no VCF rows, so
`old-snp-block-gauges.json` applies unchanged to both fixes.

The full final run took 251.48 seconds at eight threads. This is recorded wall
time, not a controlled runtime comparison with competitors.

## The 17.50 Mb connection and local costs

The exact pair 17,502,614 `CT>C` and 17,521,341 `A>G` now shares PS 16,715,242
and the same genotype orientation in both full chr20 and the owning 16–18 Mb
replay. The left block starts in the preceding chunk. A permanent two-chunk
regression asserts the exact connection, retained representations, local read
scores and agreement of independent parental majorities on the two flanks.
No in-gap heterozygote is required when reads span the gap directly.

Local reads do worsen: 136 phased overlaps previously had 131 correct and five
incorrect assignments in separate blocks; the joined block has 109 correct and
27 incorrect (80.147%). Saved same-BAM HiPhase DV output has 144 phased overlaps,
138 correct and six incorrect (95.833%). This gap is connected in the correct
parental gauge but does **not** yet match HiPhase's local read accuracy.

The 4.866–4.874 Mb exact boundary pair remains split despite a joined core SNP
path. In its owning replay coverage increases 3,844→3,870 and correct reads
3,828→3,844, with errors 16→26; local overlaps increase 40/40 correct to 51/52.
Its read-score regression records that explicit cost and extra coverage, without
changing a `spans=0` expectation. Two existing 17 Mb owning-score floors are
adjusted from 99% to the measured 98.65%/98.71% floors caused by the earlier
17.50 Mb connection; their exact suffix connections, split boundaries and
no-switch assertions are unchanged. The 3 Mb 98% floor is retained after the
physical-SNP correction restores 3,771/3,847 correct.

The first insertion-only panel run failed six quality assertions in four cases;
`initial-window-tests.log` preserves that evidence. Only the above measured
tradeoffs are accepted. Expectations are not wholesale refreshed. The panel now
contains 100 coordinates with 86 required connections; the new 17.50 Mb case has
`spans=1` and a parental-orientation check. The 24 Mb anti-inversion regression
retains its full chunk and separate four-T/eight-T allele descriptions.

## Verification

Unit tests and 906 predicate assertions across 40 cases pass. The HiFi and ONT
BAM TSV/VCF golden gates pass, including HiFi one/four-thread determinism.
The expanded gap regression suite uses fresh native outputs from the final
binary. `replay_cached_panel.py` checks binary SHA, exact argument vectors, input
stats and hashes of every output before reuse, and executes missing requests
natively. It never substitutes starting-version output. The final suite passes 7,031 assertions in 68 test cases, covering 100
coordinates and 86 required connections (`final-window-tests.log`). The 105
fresh native replay requests completed in 344.63 seconds at four concurrent
workers; the scoring pass reuses only verified final-binary outputs.
