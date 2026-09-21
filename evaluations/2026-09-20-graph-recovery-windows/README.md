# Noncentromeric graph recovery windows (2026-09-20)

The five intervals in `panel.tsv` are outside the low-MAPQ chr20:27.5–29 Mb
centromeric region. HiPhase spans each at 99.38–100% parental-truth read
concordance (`evaluations/2026-09-16-stock-deficit/deficit_scored.tsv` for
5.31 Mb; `evaluations/2026-09-16-current-deficit/deficit_scored.tsv` for the
others). Before the exact-row transfer, graph+BAM recovery split all five. `current.tsv`
records that baseline; `accepted.tsv` records the transferred result, which
spans 5.31 Mb and remains split at the other four. The separated fraction is
the share of scorable window reads one graph block places correctly.

| Gap (Mb) | HiPhase read concordance | Graph spans | Graph in-gap phased hets | Graph separated |
|---|---:|---:|---:|---:|
| 5.309–5.345 | 99.55% | no | 3 | 0.53 |
| 11.235–11.262 | 100% | no | 2 | 0.34 |
| 24.104–24.142 | 100% | no | 3 | 0.32 |
| 48.177–48.229 | 100% | no | 5 | 0.35 |
| 60.033–60.058 | 99.38% | no | 1 | 0.40 |

To rerun the panel after building `pgphase` and `test_gap_windows` and generating
`test_data/derived/chr20_truth_hap.tsv` with `scripts/make_truth_hap_map.sh`:

```bash
PGPHASE_PANEL=evaluations/2026-09-20-graph-recovery-windows/panel.tsv \
PGPHASE_EXPECT=evaluations/2026-09-20-graph-recovery-windows/current.tsv \
./test_gap_windows '[gap][windows]'
```

## Exact BAM rows at 11.23 Mb

Standalone BAM emits two complementary deletion records at 11,255,369. In its
second phasing matrix, 47 reads observe each row and the clean right SNP at
11,262,361. The first row's alleles agree with that SNP on 44/47 reads; the
second row's complementary orientation agrees on 45/47. Graph recovery instead
runs the BAM `process_chunk` with graph defaults, including
`merge_colocated_msa_alleles=true`. It produces one three-allele candidate
before transfer; allele 2 appears with right-SNP allele 0 on 21 reads and with
allele 1 on 9. The graph merge preserves those source observations, so changing
the candidate-table copy alone cannot restore the separate BAM rows.

Disabling only that MSA merge lost in-gap phased heterozygotes in committed
regressions (17 to 4 at 3.85 Mb, 6 to 1 at 4.77 Mb). Applying all BAM port
options produced separate rows but lowered graph read concordance to 0.917 at
3.85 Mb and 0.912 at 60.03 Mb, and lowered separated reads to 0.323 at 5.31 Mb.
All trials were reverted. A graph-style admission gate closed 5.31 Mb at 0.947
separated reads and 99.6% concordance, but its source still merged BAM rows, so
it was also reverted after the no-merge requirement was clarified.

## Raw catalog equivalence

Minimal `CandKey` inequality does not prove a BAM indel is private. In the
3.85 Mb recovery audit, 32 of 39 alignment-verified heterozygotes called
unmatched describe the same edited local reference as a raw catalog ALT at a
different repeat anchor. `classify_private.py` reproduces the 32/39 count from
`--recovery-audit-out`, the indexed catalog VCF, and CHM13 FASTA. The example
at 3,855,872 is already a multiallelic `CA`-anchored catalog site but appears
as `C`-anchored BAM insertion rows.

A trial excluding raw-catalog-equivalent rows removed most apparent private
sites in these windows and failed the existing phased-site floors. This shows
that sites absent from the active graph candidate table and sites absent from
the raw graph catalog are different sets. Candidate admission and exact BAM
site transfer need separate handling; no experimental behavior from this audit
was retained.


## Additional transfer experiments

The raw catalog VCF already contains both 11,255,369 deletion alleles:
`GAAA>GAA` and `GAAA>GA`. They are absent from the active graph candidate
table, not from the raw catalog. A private-site check therefore needs both
identities; a key mismatch alone would reintroduce a second description of a
catalog locus. With the deletion rows separate, the graph's allele linker
observed 42 agreeing and 5 conflicting reads between the two rows, then 45
agreeing and 2 conflicting reads to the right SNP. It observed no allele link
from the left boundary to the first deletion. The output consequently remained
split between the left boundary and the deletion. Preserving the rows is
necessary but cannot by itself create the missing left link.

A recovery solve with MSA merging, refresh, and unplaced observations all off
kept separate rows and retained the 24.1 and 48.2 Mb phased-site counts, but
the 5.31 Mb separated fraction fell from 0.53 to 0.32. Keeping unplaced
observations and adding read-supported admission closed 5.31 Mb at 0.95
separated, but lost two phased sites at 24.1 Mb. Neither passed the committed
window suite. An exact standalone BAM solve used only for new rows, alongside
graph-style observations at catalog sites, still reduced 60.03 Mb read
concordance from about 0.99 to 0.91. Even filtering those new rows to
alignment-verified heterozygotes did not restore it; the changed candidate
evidence, including records not emitted to VCF, affects the graph read solve.

The retained implementation replaces newly combined MSA candidates with
the standalone BAM rows and admits read-supported catalog sites. It closes
5.31 Mb at 0.947 separated and 0.996 read concordance while preserving the
other four focus-window floors. The 3.85 Mb window falls from about 0.99 to
0.909 read concordance (45 discordant reads). The user accepted that local
change and the full chr20 tradeoff below.

| Full chr20 | Evaluated reads | Discordant | Error rate |
|---|---:|---:|---:|
| Before | 219,058 | 2,527 | 1.154% |
| Separate BAM MSA rows | 220,640 | 3,302 | 1.497% |

Both rows use the same diplinator truth BAM and
`scripts/evaluate_phase_accuracy.py` settings. Among 219,011 reads scored in
both, 542 become discordant and 326 become concordant; the new run also tags
additional reads. The largest new-error cluster is chr20:19.3–19.4 Mb (350
newly discordant shared reads). This is a measured accuracy tradeoff, not a
parity claim. The output at 11,255,369 contains separate `GA>G` and
`GAA>G` records, as standalone BAM does.

The graph phaser stores `bam_injected` on imported candidates, but
`allele_depths_call_het` checks for a retry window before that flag, while
the graph re-solve has no retry windows. Allowing injected sites through that
guard changed the committed 22.98 Mb window's separated fraction from at
least 0.36 to 0.308, without closing the focus gaps. Restricting it to
`alignment_verified` candidates still failed. This guard requires
read-link validation before its scope can safely be widened. That guard
experiment was reverted.

## Direct whole-chunk injection trial

A second screen of the noncentromeric competitor deficits found four windows
where HiPhase is at least 95.5% concordant and the standalone BAM pipeline also
has one phase set across the endpoints: 4.85, 47.74, 60.95, and 64.91 Mb.
`bam_recoverable.tsv` records them.

Using the standalone BAM solve as the only recovery source did inject its rows
verbatim and then rerun both normal graph chunk rounds. It closed none of the
four. It also reopened 5.31 Mb (separated fraction 0.947 to 0.323), reduced the
3.85 Mb in-gap site count 18 to 17, and lowered 60.03 Mb read concordance to
0.912. Applying the longcallD phasing options to the whole graph chunk still
closed none and lowered 3.85 Mb concordance to 0.887. Those trials were
reverted.

The reason is observable in the emitted rows. At 4.85 Mb the exact BAM
`ATTT>A` row joins the left graph block, but the sequence-equivalent catalog
site at 4,878,943 remains a separate block. At 47.74 Mb the complementary BAM
deletions join the left block, while the right insertion starts a new block.
Promoting every sequence-equivalent graph site from the BAM verdict closes
64.91 Mb, but also makes an unsupported switch across 60.03 Mb: concordance
falls to 59.2% (144 discordant reads). Requiring pure BAM links on both flanks
prevents that false join but closes none of the four.

The retained union therefore keeps graph observations, replaces only a
co-located MSA candidate that graph compressed with the standalone BAM rows,
and re-runs the ordinary graph chunk rounds. Shared repeat sites still require
two pure BAM flank links. This preserves the exact 11.23 Mb deletion
representation and the accepted 5.31 Mb closure without importing BAM phase
labels or making unsupported joins.

## Graph-only gaps shorter than 10 kb

The recovery-disabled graph VCF has 74 internal phase-span gaps below 10 kb,
totalling 281,849 bp. The first BAM comparison used a noon chr20 output that
predated the final longcallD parity fixes. It reported 13 spans and must not be
used. `compare_short_gaps_bam.py` reruns the measurement from the current
parity-fixed output: standalone BAM spans 36 gaps (119,456 bp). All 36 have
truth-scorable HP-tagged reads physically crossing both graph flanks; 1,135 of
1,235 read-gap observations agree with the phase set's whole-block truth
orientation (91.90%), and 19 gaps are at least 98% concordant.

`short_gap_hiphase_comparison.tsv` checks the same 74 flank pairs against the
frozen HiPhase 1.6 output made from the shared DeepVariant calls. HiPhase places
50 gaps (195,219 bp) inside one VCF phase set. Forty-nine have truth-scorable
crossing reads. Their whole-block-oriented concordance is 1,631/1,867 (87.36%):
21 gaps reach 98%. Allowing each gap to flip locally gives 1,727/1,867 (92.50%)
and 24 gaps at 98%, but the extra three lie in a 3.02 Mb phase set whose local
orientation has switched relative to its global truth assignment.

### Why BAM misses correct HiPhase spans

HiPhase has seven at-least-98%-concordant spans that current BAM splits. One,
26,917,467–26,925,879, is on the centromeric low-MAPQ shoulder and is excluded
from recovery work: only 3 of 58 alignments overlapping its 18 HiPhase SNPs meet
BAM's default MAPQ 30 filter.

The six noncentromeric misses are linkage failures, not missing BAM evidence.
BAM emits the relevant variants in every interval. Four breaks arise because
the exact longcallD linker excludes homopolymer indels from its heterozygous
link list. The apparent gap is therefore measured from a much older eligible
anchor: 21.0–45.5 kb rather than 3.7–7.5 kb. The right candidate then has 0 or 1
supporting links, below longcallD's two-read minimum. This explains 6.58,
11.26, 12.27, and 48.23 Mb. At 48.23 Mb there is an additional genotype error:
BAM projects the HiPhase-heterozygous 48,204,383 deletion as homozygous.

The 4.78 and 34.10 Mb breaks are representation/order failures. HiPhase phases
one multiallelic record at each locus. Exact BAM output preserves two
complementary rows, and upstream-compatible phase linking checks only the
immediately preceding heterozygote. At 4.78 Mb the first row has 0/0 votes back
to the old block and starts a new phase set; the second has 44/0 votes only to
the first row. At 34.10 Mb the corresponding votes are 1/0 and 63/14. Strong
support on the second row cannot repair the break introduced by the first.

`hiphase_correct_bam_misses.tsv` records the seven diagnoses and their link
votes. `compare_short_gaps_bam.py` and `compare_short_gaps_hiphase.py` reproduce
the two span tables from a phased VCF/BAM and the existing per-read truth
results.


## Outside-in recovery frontier

Recovery now injects the BAM rows, runs the normal graph solve, and expands the
neighboring phase blocks one locus at a time. Only the closest
alignment-verified recovery locus on each side is exposed. The exposed loci
for each disconnected phase-set pair are ranked by net same-versus-cross
support, total paired reads, distance, and coordinate. One locus per pair
enters before the chunk is solved again; co-located BAM rows remain separate
and each row must qualify independently. Keeping the ranking per pair prevents
stronger unrelated gaps in the same production chunk from consuming both
recovery waves.

A row already carrying the tested boundary's phase set is not an expansion.
This check removes a false chr20:15.10 Mb join where a 29/6 vote merely
confirmed the deletion's existing block, then made that deletion eligible to
join the next block. With two waves, the minimum needed for 48.23 Mb, full
chr20 has 3,595 discordant of 220,610 truth-scored reads (1.629%). Allowing the
same-block evidence produced 4,691/220,610 (2.126%); unrestricted convergence
produced 4,682/220,610 (2.122%). The retained exact-row result before frontier
growth was 3,302/220,640 (1.497%). The 15,095,642–15,101,261 interval is in
the committed panel as a no-span control; the retained path leaves it split at
99.2% read concordance. The full chr20 VCF produced with normal chunking spans
all four target gaps; `full_chr20_target_spans.tsv` records the covering phase
set and extent for each.

`frontier_short_gaps.tsv` records the six noncentromeric HiPhase-correct BAM
misses. The frontier closes 6.58, 11.26, 12.27, and 48.23 Mb. Their window read
concordance is 99–100%. The 4.78 and 34.10 Mb representation cases stay open
because their first external edge is below the two-read threshold. The original
six-window regression panel remains green, including the deliberately open
3.85 Mb interval where HiPhase's join is only 54.3% correct.

## Seam-only targeting

Recovery targeting was reduced to bounded gaps between neighboring phase sets.
The unphased-read bin scan and graph-seeded noisy regions were removed: the seam
already contains the unphased reads between its two anchors, and a recovery
site without decisive boundary support remains excluded by the frontier gate.
Terminal and wholly unanchored regions are intentionally outside this extension
model because they do not have phase boundaries on both sides.

The rebuilt seam-only binary keeps every intended internal closure. The 15.10
Mb no-span control remains split and measures 98.71% local concordance, so its
concordance floor moved from 0.99 to 0.98 while the exact no-span assertion stays
unchanged. The resulting suite passes all 229 assertions.

Whole chr20 phases 217,543 of 245,053 emitted reads in 320 phase sets. Of 217,529
truth-evaluated reads, 3,174 are discordant (1.459%; 98.54% accuracy). Compared
with the preceding frontier build, seam-only targeting removes 3,103 phased
reads and 35 phase sets while reducing discordant reads by 421 and improving
accuracy from 98.37%. The lost coverage came from recovery regions without two
neighboring phase-set anchors, which are intentionally outside the extension
model. `--graph-noisy-msa` and its seeding helper were removed with the only
consumer of those seeded regions.

The seam detector was then audited independently. The graph adapter had
initialized unphased candidates to `-1`, while longcallD initializes candidates
to `0` and reserves `-1` for unphased reads. Graph candidate and read sentinels
now match that contract. One shared anchor predicate requires a positive
phase-set label and different haplotype alleles, so neither an unassigned
candidate nor a homozygous row can define a seam. Extents are accumulated in
candidate coordinate order through an unordered label-to-index map, then
overlapping or touching extents are coalesced before gaps are emitted. This
changes expected detector cost from O(C log B + B log B) to O(C+B). A full
chr20 rerun retains 77,484 candidate rows and produces byte-identical phased
VCF and BAM output. Exactly 17,276 TSV rows differ, only in `PHASE_SET`, and
each difference is the intended `-1` to `0` candidate-sentinel normalization.

The detector and target construction were subsequently converted to flat,
ordered scans. Each candidate contributes `[phase_set, position]` to a vector
merge stack, so seam detection is amortized O(C) without a hash table. Target
construction scans the ordered anchors and seams once, stores group membership
as index ranges, and uses binary search for later candidate membership.

A distinct-locus, symmetric-flank trial was rejected: it produced 3,201
discordant reads and evaluated 54 fewer reads on full chr20. Retaining the
validated flank semantics yields 217,543 phased reads in 320 phase sets and
3,170/217,529 discordant truth-scored reads (1.457%; 98.54% accuracy). All 229
window assertions pass. The output contains 77,448 candidate rows.

The recovery transfer itself was then audited without changing its decisions.
Parent membership now comes from the raw and translated indexes directly,
candidate orientation is stored with the transferred candidate, disabled audit
output performs no audit-row string copies, and missing-read detection is a
linear qname merge scan. The full chr20 TSV, phased VCF and phased BAM are
byte-identical to the preceding flat-detector output.
