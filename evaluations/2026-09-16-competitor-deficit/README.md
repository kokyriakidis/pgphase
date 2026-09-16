# Gaps the competitors close and we do not

Our unresolved gaps are only a competitive deficit where a competitor actually
phases across them, so `find_deficit_gaps.py` tests every gap in the pass-1
inventory against each competitor's phase sets (same CHM13 coordinates, so the
bounds are used directly; a spanning set must start at or before the gap and end
at or after it).

## The deficit set

| class | gaps | span |
|---|---:|---:|
| spanned by all 3 tools | 114 | 2.28 Mb |
| spanned by 2 of 3 | 54 | 1.72 Mb |
| spanned by 1 of 3 | 19 | 0.74 Mb |
| **deficit total (>= 1 tool)** | **187** | **4.74 Mb** |
| spanned by none -- nobody phases these | 124 | 13.07 Mb |

So of 311 gaps and 17.81 Mb, the genuine deficit is **187 gaps and 4.74 Mb**;
the remaining 13.07 Mb is not a deficit at all -- no tool in the panel phases it.
That reframes the earlier drilldowns: `35,919,404` and `13,429,829`, the two
largest gaps in the eight-window panel, are in the nobody-spans class (checked
directly: HiPhase, LongPhase and WhatsHap all break at the same coordinates we
do), which is consistent with their interiors having steps no read crosses.

Deficit gaps are small. Mean size is 25 kb against 105 kb for the
nobody-spans class, so this is a many-small-gaps problem.

One qualification on the nobody-spans class: "no tool spans the whole gap" does
not mean no tool phases *inside* it. In `35,919,404` HiPhase forms a 131.8 kb
block (`35,981,805-36,113,633`) covering the gap's middle, which we split into
pieces -- it bridges the second of our three zero-spanning-read deserts using
interior sites. Of its 23 phased sites in that block we call 13 at the same
position and most of the rest are indels with a call within a base or two
(`35,994,697`/`35,994,698`, `36,040,406`/`36,040,407`, `36,077,134`/`36,077,135`,
`36,107,334`/`36,107,325`), with 4 genuinely absent. So part of the 13.07 Mb is
recoverable as sub-blocks even though the full gap is not closable by anyone;
the full-span test used here measures joins, not interior phasing.

## Gap 1: `chr20:61,732,321-61,810,469` (78.1 kb, all 3 tools span it)

The largest deficit gap. Hybrid with recovery already covers most of it -- blocks
`61,690,751-61,738,233` and `61,757,551-61,858,547` -- so the real unphased break
is **19.3 kb, `61,738,233-61,757,551`**. The audit and normal runs produce
identical final blocks -- but they are not equivalent in what they evaluate, and
that turns out to matter (see the tier-4 section below): `--gap-decision-audit`
makes the homopolymer tier run in recovery pass 1, which a normal run skips.

Not a discovery problem. Of HiPhase's 17 phased sites inside the gap we hold 8 at
the same position and 6 more within 15 bp (indel anchor placement differs), and
only 3 are genuinely absent -- all multi-allelic repeat expansions
(`TTTTTTT > T,TTTTTTTTTTTTTTTTT,...`, `ACACACA > A,ACA`, `C > CA`).

Inside the break the evidence is good and the chain is complete:

| pos | type | msa_verified | homopolymer | our category | reads | segregation vs truth | admitted by |
|---|---|---:|---:|---|---:|---:|---|
| 61,747,508 | DEL | 1 | **1** | NoisyMsaHet | 56 | **0.929** | tier 4 only |
| 61,747,510 | DEL | 0 | 0 | **RepHetIndel** | 45 | **0.956** | never |
| 61,755,064 | INS | 1 | 0 | NoisyMsaHet | 51 | **1.000** | tiers 3, 4 |
| 61,757,551 | SNP | — | 0 | CleanHetSnp | 49 | 1.000 | all (endpoint) |

Reads spanning each consecutive step: 60, 114, 42, 80. Nothing is coverage-limited
here. And HiPhase is right: on the 138 reads overlapping our break it scores
**132 concordant, 1 discordant (99.25%)**.

So the blocker is the indel-class exclusions, and here they are excluding correct
signal. Tier 1 admits only clean germline, tier 2 only MSA-verified SNPs (both
MSA sites here are indels), tier 3 admits `61,755,064` but excludes `61,747,508`
on `is_homopolymer_indel`, and `61,747,510` -- the 0.956-segregating repeat-demoted
indel -- is admitted by no tier at all. Tier 4 exists for exactly this case,
admits the homopolymer indel, and still returns `rejected`.

Two distinctions worth keeping: this is a different mechanism from the two gaps
already diagnosed (`36,332,599` was the private-SNP flag gate, `35,919,404` is
read coverage), and it does not contradict the earlier class-level finding that
repeat-demoted indels are only 23% informative -- these two specific sites score
0.929 and 0.956 and happen to be load-bearing, which is an argument for
per-site evidence rather than for admitting the class.

Next question for this gap: **why does tier 4 return `rejected`** when it admits
`61,747,508` and the reads chain. Candidates are the bridge-vote requirement that
link support be established against a clean anchor (`kCandGermlineClean`, and the
only clean site in the break is the right endpoint 9.9 kb away) and the
`min_block_link_reads` margin.

## Gap 2: `chr20:36,217,274-36,268,291` (51.0 kb, all 3 tools span it)

Same protocol (`diagnose_gap.py`). Our blocks are `36,168,059-36,247,421` and
`36,268,291-36,268,558`, so the residual break is **20.9 kb,
`36,247,421-36,268,291`**. HiPhase holds 8 sites in the gap; we hold 5 at the
same position and 3 are absent, all multi-allelic repeat expansions
(`C > CT,CTT`, `T > TTT,TTTT`, `T > TT`).

| pos | type | msa | homopolymer | our category | reads | segregation | admitted by |
|---|---|---:|---:|---|---:|---:|---|
| 36,247,421 | SNP | — | 0 | CleanHetSnp | 64 | 0.906 | all (endpoint) |
| 36,252,908 | DEL | 1 | **1** | NoisyMsaHet | 65 | **0.923** | tier 4 only |
| 36,259,923 | DEL | 1 | 1 | NoisyMsaHet | 61 | 0.607 | tier 4 only |
| 36,261,164 | INS | 1 | 0 | NoisyMsaHet | 57 | **0.509** | tiers 3-4 |
| 36,261,302 | INS | 1 | 0 | NoisyMsaHet | 59 | 0.644 | tiers 3-4 |
| 36,261,311 | INS | 0 | 0 | RepHetIndel | 59 | 0.644 | never |
| 36,268,291 | SNP | — | 0 | CleanHetSnp | 71 | 0.986 | all (endpoint) |

Reads span every step (80, 78, 114, 116, 116, 58), so again not coverage-limited.
HiPhase joins it but less cleanly than gap 1: **147 concordant, 7 discordant
(95.45%)** on the 157 reads overlapping our break.

The contrast with gap 1 matters. Here the sites tiers 3-4 *do* admit are the
uninformative ones -- 0.509 is chance, 0.607 and 0.644 little better -- and the
one good interior site (`36,252,908`, 0.923) is homopolymer-flagged, so again
only tier 4 can use it. **MSA verification is not a proxy for informativeness**:
all four `NoisyMsaHet` sites here carry `msa_verified = 1` and they range from
0.509 to 0.923.

## Why tier 4 does not deliver: two different reasons

Reading the tier-4 vote matrices against `min_block_link_reads = 2`:

- **Gap 1 has an available join that is never evaluated.** In the audit export,
  tier 4 in **recovery pass 1** reports `JOINED = 1` for proposal set `61725696`
  with left support 5 and right support 4 -- both clear of the threshold, and
  clear of tier 4's extra requirement that each original haplotype independently
  favour the same orientation (`gap_recovery.cpp:349-355`). A normal run never
  sees it: `collect_pipeline.cpp:1944` skips the homopolymer tier whenever
  `recovery_pass == 1 && !audit`. In pass 0, where tier 4 does run, the proposal
  is anchored differently (left support 53, right support 0) and cannot join.
- **Gap 2 has no two-sided link in any tier or pass.** The left block carries
  left 13 / right 0, the right block left 0 / right 18; no proposal set ever
  holds support on both flanks, so `joined` cannot become true however admission
  is widened.

A third constraint sits behind both: `stitch_gap_proposal` is called
orientation-only whenever `tier == kGapHomopolymerTier` or `recovery_pass == 1`
(`collect_pipeline.cpp:1991`), so a tier-4 link is never applied in place and can
only act through the validation and independent-block path that follows.

Next step for gap 1, as a measurement rather than a fix: let the homopolymer tier
run in pass 1 and see whether the join materialises and what it costs on the
read-level gate (concordant -> discordant, plus tags lost). Gap 2 needs something
else -- its problem is the proposal's block structure, not admission.

## Fix tested: let the homopolymer tier run in the reprojected pass

One-line change in `collect_pipeline.cpp`: drop `if (recovery_pass == 1 &&
!audit) break;` from the `kGapHomopolymerTier` branch. The tier's own guards are
untouched -- it still requires a homopolymer candidate in the gap,
`--link-by-alleles`, the per-haplotype orientation agreement
(`gap_recovery.cpp:349-355`), and the BAM-only validation that can still veto it.
The join propagates as a parity edge (`job_result.edge`), so the orientation-only
stitch call is not an obstacle. Clean build, no new warnings, all five unit-test
binaries pass.

Verified on gap 1 (`chr20:61,732,321-61,810,469`), flag off, normal run:

| | before | after |
|---|---|---|
| tier 4, pass 1 | not evaluated | `joined`, `LEFT_LINK_PS = RIGHT_LINK_PS = 61725696` |
| blocks over the gap | `61,690,751-61,738,233` (8 sites) + `61,757,551-61,858,547` (48) | **`61,690,751-61,858,547` (56 sites)** |
| phase sets in window | 2 | **1** |
| reads tagged | 575 | 575 |
| concordant / scored | 567 / 575 (98.61%) | 567 / 575 (98.61%) |
| gate concordant -> discordant | — | **0** |
| concordant tags lost | — | **0** |

So the 78.1 kb gap closes and **not one read changes label** -- every read keeps
its haplotype and its concordance, and only the phase-set identity merges. That
is a cleaner result than the private-SNP gate fix, which closed its gap at a cost
of 23 correct read tags.

Gap 2 is the control and behaves as predicted: tier 4 now runs in pass 1, reports
`rejected` because no proposal block holds support on both flanks, and the output
is byte-identical (2 phase sets, gate 0, 358 reads tagged before and after).

Pre-existing observation, unchanged by this patch and not chased here: gap 2's
window scores only 77.65% per-read against truth in both arms, meaning one of its
two blocks is internally inconsistent. That is a separate defect from the gap not
closing.

Not promoted to default on this evidence: one gap closed, one control unchanged.
The chromosome-wide gate run is what decides it, since the change affects every
gap with a homopolymer candidate, and this tier is the designed last resort.

### Panel regression: inert everywhere else

Re-ran the eight-window panel with the patched build against the matched
pre-change runs (`compare_panel.py`, `hp_pass1_panel.tsv`). Nothing moves:

- newly joined: **0**; newly spanning: **0**
- gate concordant -> discordant: **0**; concordant tags lost: **0**; newly tagged: **0**
- tagged reads 3807 -> 3807; concordant 3785 -> 3785

`36,332,599` still joins (the private-SNP fix), and the rest still do not close.
So this change is not a broad win -- on the panel it changes nothing at all --
but it costs nothing and it closes gap 1, whose 19.3 kb break had complete,
informative, read-backed evidence and a competitor joining it at 99.25%.

One reporting nuance it introduces: five panel gaps whose status read `split` or
`partial` now read `rejected`, because tier 4 is evaluated in pass 1 and its
label (`tier == kGapHomopolymerTier && !joined`) is the last report row. The
underlying state is unchanged -- same blocks, same reads, same tags -- but any
consumer reading the final tier row loses the `split` versus `partial`
distinction. Worth fixing in the report rather than in the tier.

## Deficit gap 2: diagnosed to the line, not yet closed

`chr20:36,217,274-36,268,291`. The evidence is sufficient and the correct join
exists; what blocks it is an ordering flaw plus a veto, both upstream of the
recovery tiers.

### The evidence is sufficient

An offline solver over the audit's per-read observations (orient each site by
pairwise agreement, assign reads by majority vote at margin 2) shows what the
window supports:

| site set | gap 1 | gap 2 |
|---|---|---|
| clean only | no partition | no partition |
| + verified indels, non-repeat (tier 3 today) | 38 reads, 1.000 | 26 reads, **0.538** |
| + homopolymer, cumulative (tier 4 today) | 73 reads, 1.000 | 65 reads, 0.846 |
| clean + homopolymer only | 11 reads, 1.000 | **52 reads, 1.000** |

Gap 2's poison is the *non-repeat* verified indels (truth segregation 0.509 and
0.644); its one good site is the homopolymer deletion at 0.923. Because the
tiers were cumulative, the junk admitted at tier 3 was never removed when tier 4
added the good site. Class hierarchy does not predict per-site quality here --
though it does on average, which is why the ladder is still clean SNPs first:
measured earlier on this chromosome, verified SNPs segregate at a median 0.988
against 0.895 for verified indels.

### Two truth-free screens fail, one works

Do not reach for these again as they stand:

- **Pairwise agreement with a clean anchor**: the 0.923 site scored 0.744 while a
  0.644 site scored 0.818. Agreement between two noisy sites compounds their
  errors, so a good site against a 0.906 anchor expects only ~0.84.
- **Leave-one-out consistency against the consensus partition**: the junk sites
  scored *higher* (0.820, 0.857) than the good ones (0.754, 0.737).

The reason is positional. `36,261,164`, `36,261,302` and `36,261,311` are one
repeat event called three times inside 150 bp; they agree with each other and
form a self-consistent clique that outvotes the true signal -- the same failure
shape as the earlier column-discovery result, where self-consistency read 1.000
at chance-level truth concordance. **Collapsing candidates within 300 bp to the
best-covered one** lifts gap 2 from 0.846 to 0.981 and leaves gap 1 at 1.000.
That is the screen worth implementing.

### What actually blocks the join

`stitch_gap_proposal` counts a link vote only from a read that still holds a hap
and one of the flanks' phase sets. Two read-tagging filters run *before*
`recover_hybrid_gaps` in the batch loop, so the labels the stitch needs are
already erased:

| reads overlapping | n | frozen phase set |
|---|---:|---|
| left endpoint `36,247,421` | 64 | **62 unphased**, 2 in the flank block |
| right endpoint `36,268,291` | 71 | 68 in the flank block |

Those 62 reads observe exactly one clean het SNP, and the margin counts only
`kCandCleanHetSnp` (its rescue credits bridge SNPs, never indels), so 1 < 2
strips them. No tier configuration can reach them, and no evidence-crediting fix
helps because the break contains zero MSA-verified SNPs.

**The correct join does exist.** Running recovery on unfiltered labels, both
flanks link to the same proposal phase set (`36247421` on both sides), the
tier-4 proposal scores **1.0000** against read truth over 95 reads, and it
agrees with the frozen haplotype of **all 26** left-flank and **all 69**
right-flank reads it holds.

### Lowering the margin is not the fix

Ten-window panel, current build, read-level gate matched by read name:

| | margin 2 | margin 1 | margin 0 (8 windows) |
|---|---:|---:|---:|
| reads tagged | 4,687 | 5,740 | 3,699 -> 5,034 |
| accuracy | 97.63% | **95.44%** | 99.38% -> **96.48%** |
| concordant -> discordant | -- | **99** | **30** |

Gap 2 does span at margin 1 -- by flipping 70 previously concordant reads, so
that join is wrong. `--min-read-margin` is also not a BAM-side remnant: it is
the global `min_read_hap_margin`, parsed by both subcommands and consumed by the
chromosome pass, gap recovery's proposal and validation, and the graph path's own
read gate. Note the built-in default is 0 while this project's canonical scripts
pass 2.

### Reordering the filters is not the implementation either

Moving both filters after recovery does fix the starvation -- all 64 left-boundary
reads keep a frozen phase set -- but it then makes the output filter judge reads
on counters recomputed over a narrow gap window, and tagging collapses: 358 -> 177
reads in this gap and 575 -> 398 in gap 1, losing 105 and 169 concordant tags with
0 concordant -> discordant. Crediting bridge and last-resort evidence at that call
does not restore them. Reverted.

The implementation that follows from the diagnosis is a **snapshot**: build the
stitch's `read_index` from the haps and phase sets as they stand before the
output filters, and leave output filtering exactly where it is. That preserves
current tagging by construction while giving the votes an unstarved view.

### And the veto needs its own look

With the starvation removed, tier 4 reports `vetoed`, not `joined`. The
homopolymer tier requires an independent BAM-only solve to reach the same
orientation, which is a guard worth keeping -- margin 1 proved a wrong join is
readily available in this window. But it fires here against a proposal measured
1.0000 accurate, with `SELECTED_GRAPH_READS=157` and `GRAPH_BAM_PASS=1` in pass 1,
so either the validation solve does not join or its orientation differs. That is
the next thing to instrument.

### Kept from this session

Three changes, each measured, all opt-in-by-default-behaviour:

1. Tier 4 is now a real last resort -- clean plus verified SNPs plus verified
   homopolymer indels -- instead of restoring every flag and readmitting the
   non-repeat indels tier 3 had just failed with.
2. `CandidateVariant::hp_gap_scorable` lets the one site the homopolymer tier
   admitted contribute read scores. `init_assign_read_hap` had skipped *every*
   homopolymer indel, so the tier could earn link support and still not phase the
   reads it was reached for. Measured effect on gap 2: the in-gap block grew
   69 -> 95 reads and extended ~17 kb leftward, reads reaching the left endpoint
   1 -> 26. `gap_link_supported` could not be reused for this: line 1188 grants it
   from `msa_insertion_alts` alone, also outside a homopolymer gap.
3. `ReadRecord::n_hp_gap_agree/conflict`, credited in the recovery margin
   filter's rescue for the same reason bridge SNPs are.

Gap 1 still joins with all three applied; gap 2 does not yet.

## The starvation bug is fixed; gap 2's last blocker is the validator

The diagnosis above named the fix as a snapshot, and that is what was built --
but it took three attempts to scope, because the filtered view turned out to be
load-bearing in three places, not one.

`GapReadIndex` gains `vote_assignments`: the same read-to-(hap, phase set) map,
falling back to the label a read held **before** the output filters ran wherever
its current one is empty. `PreFilterLabels` is captured in the batch loop just
before `filter_hybrid_reads_by_margin` and threaded into `recover_hybrid_gaps`.
The gap inventory is still derived from the filtered chunks, so nothing about
what is emitted changes -- only what the link votes can see.

**Scoping was the hard part.** Each widening of the filtered view cost read
coverage, because three consumers treat "no committed assignment" as permission
to act:

| widened | effect |
|---|---|
| `assignments` (whole index) | `emit_independent_gap_block` stops treating filter-zeroed reads as gap-only and re-phasing them -- 105 and 169 concordant tags lost in the two gaps |
| `first_assignment` | also builds `supported_phase_sets`, which gates emission -- same loss |
| `original` in `stitch_gap_proposal` | also means "already in a block, do not attach" -- same loss |

So the final form keeps `assignments`, `first_assignment` and `original`
filter-accurate, and adds one wide map, `vote_original`, used by the vote loop
alone. The validation view's own output filters were also removed, for the same
reason they were wrong on the flanks: `stitch_gap_proposal` counts a vote only
from a proposal read still holding a hap and phase set, so filtering the
validator for reporting silenced the validator.

### Measured

Left-boundary voters at `chr20:36,247,421` went from **2 of 64** reads to
**64 of 64**, and both flanks now link to the same proposal phase set
(`36247421` on both sides) -- the join condition. Ten-window panel against the
matched pre-change runs (`votesnapshot_panel.tsv`):

- gate concordant -> discordant: **0**
- tagged reads 4,687 -> 4,714; concordant 4,576 -> **4,606** (net +30)
- nine of ten windows byte-identical; all movement is in `61,732,321`, which
  gains 78 concordant and 2 discordant tags while 48 concordant tags move into
  the merged block, and its window accuracy rises 98.61% -> 99.17%

### Gap 2 is still not closed, and the remaining blocker is now isolated

Tier 4 reports `vetoed`, and the veto's own counters say why:

```
gap=36247421-36268291  bam_reads=157  bam_joined=0  bam_flip=0  graph_flip=0
```

Not an orientation disagreement -- the orientations agree. The homopolymer
tier's validator has 157 reads and its independent BAM-only solve simply does
not link the flanks, so it cannot confirm a join that scores 1.0000 against read
truth over 95 reads. The guard is left in place: margin 1 proved a wrong join is
readily available in this window, so a veto that fires too often is the right
failure direction. What needs work is the confirmation mechanism --
`select_graph_gap_bam_reads` restricts the validator to reads carrying a
graph-channel call, so it is asked to reproduce a join from a smaller read set
than the proposal it judges. That is the next measurement, and it is a design
change to the validator rather than another scoping fix.

## Why the validator fails on a proposal that is 1.0000 accurate

A temporary probe in the veto branch (removed again) dumped the validator's own
proposal. It is not disagreeing about orientation -- it never forms a block that
reaches both flanks:

```
gap=36247421-36268291  left_ps=36209945  right_ps=36268291
selected=157  total=510  skipped=353  phased=151
  ps=36247421   reads=63   left-flank reads=63   right-flank reads=0
  ps=36259922   reads=41   left-flank reads=0    right-flank reads=41
  ps=36261163   reads=47   left-flank reads=0    right-flank reads=29
```

Three blocks, each touching exactly one flank, against the graph proposal's
single 95-read block holding 26 left-flank and 70 right-flank reads. Two things
cause it:

1. **The validator solves on a subset.** `select_graph_gap_bam_reads` admits
   only reads carrying a graph-channel call -- 157 of the 510 reads in the
   window, with 353 skipped -- so it is asked to reproduce a join from
   materially thinner coverage than the proposal it judges.
2. **It fragments exactly at the repeat cluster.** Its block boundaries are
   `36,259,922` and `36,261,163`, which are the junk sites `36,259,923` and
   `36,261,164` -- the same one-repeat-event-called-three-times that drags the
   full solve from 1.000 down to 0.846. With 157 reads there is not enough
   overlap to chain across them, so the partition breaks there.

Both point the same way, and neither is fixed by weakening the veto. The
root-cause fix is the positional screen the site analysis already identified:
collapse candidates within 300 bp to the best-covered one, which lifts the full
solve 0.846 -> 0.981 and would remove the boundaries the validator breaks at.
Widening the validator's read set to match the proposal it judges is the second
candidate, and the weaker one -- it makes the check easier to pass without making
the underlying evidence any cleaner.

## Applying the default pipeline's own merge rule to recovery: tested, wrong

The default pipeline merges blocks in three steps, and only two mechanisms:

| mechanism | scope | evidence |
|---|---|---|
| `flip_chunk_hap` -> `select_stitch_orientation` -> `apply_chunk_flip_and_merge` | adjacent **chunks** | overlap reads present in both chunks, counted as an n11/n12/n21/n22 table |
| `stitch_phase_blocks_with_pgbam` / `stitch_adjacent_chunks_with_pgbam` | blocks **within** a chunk, and chunk seams where the read vote fails | GBWT haplotype threads, requires `--pgbam-file` |

`select_stitch_orientation` is the whole standard: under the default
`kStitchRuleNetMargin` it merges iff `|(n12+n21) - (n11+n22)| > stitch_min_margin`
and takes the orientation from the sign. After a merge,
`apply_chunk_flip_and_merge` rewrites the downstream phase set to the upstream id
and flips hap labels if needed, and
`propagate_overlap_read_phase_to_output_owner` copies HP/PS onto overlap reads
that were unphased -- but only for pairs that merged, because otherwise the
relative phase is unknown.

Recovery's `stitch_gap_proposal` already builds the same table, per flank, and
`GapStitchResult::right_flip` is the same orientation bit. So the parity question
is narrow: recovery adds the homopolymer tier's independent BAM-only re-solve on
top, and nothing in the default pipeline does that.

**Tested: skip the re-solve when the table is unambiguous.** Gap 2's tier-4 vote
table is 26-0 on the left and 69-0 on the right -- under the default rule, scores
of -26 and -69, both unambiguous merges, three times the support that joins
deficit gap 1 (19 and 10). Making confirmation conditional on that produced the
join:

| | before | after |
|---|---:|---:|
| spans the gap | no | **yes** (`36,168,059-36,268,558`) |
| reads tagged | 358 | 358 |
| window accuracy | 77.65% | **58.38%** |
| concordant -> discordant | -- | **70** |

**The join is wrong, and the veto was right.** All 70 broken reads were in the
right flank's phase set `36268291`, all land in the merged block, and all sit
inside the break: the merge inverted the right flank relative to the left. Change
reverted.

### Why the same rule gives a wrong answer here

The rule is not the problem; its inputs are. In the default pipeline the voting
reads are phased **independently, twice** -- once by each chunk's own
full-coverage solve -- and a chunk seam is an artifact of chunking, so the two
solves really are independent estimates of the same reads' haplotypes. In
recovery the voters at the left flank's boundary are, measured, **62 of 64 reads
that the output filter had zeroed**, each phased from a single clean het SNP.
Admitting them makes the table unanimous, but they share one point of failure, so
26 votes carry about one vote's worth of independent information. Unanimity among
correlated weak voters is not confirmation, which is exactly the hole the
BAM-only re-solve was covering.

That reframes what "parity with the default pipeline" should mean for recovery.
Not a looser confirmation, but voters that resemble the default pipeline's: reads
whose own phase is independently supported. Two candidates, in order:

1. **Collapse the repeat cluster** (`36,259,923`/`36,261,164`/`36,261,311` are one
   event called three times inside 150 bp). It is the root cause of both the
   solve's 0.846 ceiling and the validator's fragmentation, and it is what makes
   the boundary reads' single-site phase unreliable in the first place.
2. **Weight votes by supporting evidence** -- admit a filter-zeroed read as a
   voter only when more than one site supports its hap. On this gap that would
   return the left flank to 2 voters and correctly refuse the join.

Note the interaction with the committed vote-snapshot change: with the veto in
place nothing regresses (gap 2 stays `vetoed`, panel gate clean), but that change
is what makes this wrong join *reachable* if the confirmation is ever loosened.
The two must not be relaxed together.

## Removing the repeats: measured, and it fixes the orientation

The wrong join above means the merge inverted one flank, i.e. the proposal's own
polarity differs between its left-touching and right-touching reads. That is
testable from the audit observations with truth alone, no flank labels needed:
build the partition, pick the single global orientation truth prefers, then score
the reads that reach each endpoint separately. A proposal that is internally
consistent scores the same both ways; one carrying a switch does not.

| site set | sites | left reads | left acc | right reads | right acc |
|---|---:|---:|---:|---:|---:|
| all in-break sites | 7 | 30 | 93.3% | 22 | **81.8%** |
| drop the `REP` class | 6 | 28 | 100.0% | 20 | **85.0%** |
| collapse to one site per 300 bp | 5 | 28 | **100.0%** | 17 | **100.0%** |
| clean + homopolymer only | 4 | 27 | **100.0%** | 18 | **100.0%** |
| clean only | 2 | -- | no partition | -- | -- |
| oracle, informative sites only | 3 | 29 | 100.0% | 4 | 100.0% |

So yes -- removing the repeats fixes it. With all seven sites the two halves
disagree under one orientation (93.3% against 81.8%), which is exactly the switch
that inverts the right flank when the blocks merge. Collapse the cluster and both
halves read 100.0% while still holding 28 left-reaching and 17 right-reaching
reads, so the merge orientation would be correct and the 70-read flip would not
occur.

**The removal has to be positional, not by class.** Dropping the `REP` category
leaves the right half at 85.0%, because only one of the three correlated calls is
`REP` -- `36,261,164` and `36,261,302` are `NoisyMsaHet` with `msa_verified = 1`.
Class hierarchy cannot see that three calls inside 150 bp are one event; only
position can. This is the same conclusion the site-quality work reached from the
other direction, where leave-one-out consistency ranked that clique *above* the
genuinely informative sites.

Also note `clean only` yields no partition at all: two sites 20.9 kb apart cannot
chain. So the answer is not "trust only clean SNPs" -- it is clean SNPs plus one
representative of each repeat event, which is what the 300 bp collapse produces
and what reaches 100% on both halves.

### Caveat: the offline result did not transfer to the pipeline

Implemented the collapse where it belongs -- a positional demotion applied to
`original_flags` right after they are captured, so every tier inherits it, with
redundant cluster members set to `kLongcalldRepHetVar` (which
`kCandGermlineVarCate` excludes, so they drop out of the solve as well as out of
link eligibility). Paired with the default-rule confirmation, gap 2 still joins
**and still flips the same 70 reads**, byte-identical to the run without the
collapse: 358 tags, 58.38%.

So one of two things is true, and the next diagnostic separates them: either the
demotion is not reaching the candidates the solve actually uses (the `in_gap`
selection or the observation counts are wrong), or the switch inside the
pipeline's k-means has a cause other than the three correlated calls, and the
offline agreement at 100%/100% was a property of my solver rather than of the
evidence. The measurement to run is the pipeline analogue of the halves test:
export the tier-4 proposal and score its left-reaching and right-reaching reads
separately under one orientation, with and without the collapse.

Both experimental changes are reverted; `src/` is back to the pushed commit,
gap 2 `vetoed`, nothing regressed.

## The phasing inside both gaps, scored on its own (no join)

Fresh audits from the committed build. For each recovery pass and tier, the
largest proposal block is scored against read truth under a single global
orientation, and the reads reaching each endpoint are scored separately -- so an
internally inconsistent partition (a switch) shows up as one half disagreeing
with the other.

### Gap 1, `chr20:61,732,321-61,810,469`, break `61,738,233-61,757,551`

Interior sites: `61,747,508` MSA DEL homopolymer seg 0.929; `61,747,510` REP DEL
seg 0.956; `61,755,064` MSA INS seg 1.000; `61,757,551` clean SNP seg 1.000.

| pass/tier | blocks | reads | accuracy | left half | right half |
|---|---:|---:|---:|---|---|
| pass 0, tiers 1-2 | 1 | 115 | 100.0% | 11/11 | -- |
| pass 0, tier 3 | 1 | 109 | 100.0% | 9/9 | -- |
| pass 0, tier 4 | 2 | 109 | 100.0% | 9/9 | -- |
| pass 1, tiers 1-3 | 2 | 20 | 100.0% | 20/20 | -- |
| **pass 1, tier 4** | 1 | **45** | **100.0%** | **27/27** | **20/20** |

### Gap 2, `chr20:36,217,274-36,268,291`, break `36,247,421-36,268,291`

Interior sites: clean SNPs at the edges (`36,247,421` seg 0.906, `36,268,291` seg
0.986), `36,252,908` MSA DEL homopolymer seg 0.923, then the correlated cluster
`36,259,923` seg 0.607, `36,261,164` seg 0.509, `36,261,302` seg 0.644,
`36,261,311` seg 0.644.

| pass/tier | blocks | reads | accuracy | left half | right half |
|---|---:|---:|---:|---|---|
| pass 0/1, tiers 1-2 | 2 / 1 | 69 | 100.0% | 1/1 | 69/69 |
| pass 0, tier 3 | 2 | 41 | **97.6%** | -- | -- |
| pass 1, tier 3 | 1 | 41 | 100.0% | -- | 41/41 |
| **pass 0/1, tier 4** | 2 / 1 | **95** | **100.0%** | **26/26** | **70/70** |

### Two things this settles

**The phasing is not the problem in either gap.** At tier 4 both gaps produce a
single block that reaches both endpoints and is perfectly concordant with truth
-- gap 1 at 45 reads (27 left, 20 right), gap 2 at 95 reads (26 left, 70 right).
Every read scored, in both halves, in both gaps. Tier 3 is the only tier that
degrades (gap 2, pass 0: 41 reads at 97.6%), and it is the tier that admits the
non-repeat verified indels, consistent with the escalation-ladder result.

**Correction: gap 2's proposal carries no internal switch.** An earlier offline
solve of the same observations put its left-reaching reads at 93.3% and its
right-reaching reads at 81.8% and attributed the wrong merge to that
disagreement. The pipeline's own tier-4 proposal reads 26/26 and 70/70 -- no
switch at all. That 93.3/81.8 split was a property of the offline solver, which
voted over all seven interior sites with a plain majority instead of the tier's
filtered set, and it should not have been carried into an explanation of the
pipeline's behaviour.

So the 70-read flip seen when the join is forced does **not** originate in the
proposal. It has to come from the flank side of the comparison: the left flank
contributes only 26 voters, all of them reads the output filter had zeroed and
each phased from a single clean het SNP -- and the nearest such SNP,
`36,247,421`, itself segregates at only 0.906. The next diagnostic is therefore
whether those 26 boundary reads are correctly phased *relative to their own
block*, not whether the proposal is self-consistent.

## Why the join cannot be correct: an earlier join in the same window is wrong

Gap 2's phasing is perfect and its votes are sound, so the failure had to be on
the other side of the comparison. It is.

### The voters are fine

Probed every read that casts a link vote for proposal `36247421` (temporary
probe, removed). Left flank `36209945`: 64 voters, proposal labels 90.6% correct
against truth, flank labels 90.6% correct, both `hap1=PAT`. Right flank
`36268291`: 69 voters, 100.0% both, also `hap1=PAT`. Each flank's voters agree
with the proposal *and* with truth, so the local comparison is sound on both
sides.

### The block they point at carries a switch

The emitted left block `36168059` spans `36,144,761-36,267,248` with 287 reads at
only **72.5%** truth consistency, and the opposite-polarity reads are not
scattered -- they appear abruptly:

| 10 kb bin | reads | opposite |
|---|---:|---:|
| 36,140,000-36,190,000 | 125 | **0%** |
| 36,200,000 | 41 | 2% |
| 36,210,000 | 39 | 31% |
| 36,220,000 | 23 | 57% |
| 36,230,000 | 36 | **92%** |
| 36,240,000 | 23 | 87% |

Everything before ~36,210,000 is perfectly consistent; everything after is
inverted. That is a switch error, and `36,209,945` -- the phase set recovery
links gap 2's left side to -- is exactly where it starts.

### The switch is a gap-recovery join (culprit corrected below)

The same run's recovery report holds an upstream gap:

```
gap 36,172,778-36,209,945  tier 1  L=1 R=1  leftPS=36168059  rightPS=36168059
                                   reads_added=81  STATUS=joined
```

Re-running the identical window **without** `--recover-gaps` settles it:

| | block `36168059` | window total |
|---|---|---|
| recovery off | 169 reads, span `36,148,502-36,233,980`, **99.4%** | 241/242 = **99.59%** |
| recovery on | 287 reads, span `36,144,761-36,247,400`, **72.5%** | 278/358 = **77.65%** |

Recovery is what breaks this window. The tier-1 join across
`36,172,778-36,209,945` attaches 81 reads and inverts everything downstream of
the seam, taking the block from 99.4% to 72.5% and the window from 99.59% to
77.65%. It buys 116 more scored reads and 37 more concordant ones while creating
79 discordant ones.

### What this means

Gap 2 is not a gap-recovery *admission* problem and never was. Its own solve is
perfect, its votes are correct, and it correctly aligns with the segment next to
it -- but that segment is already inverted relative to the rest of its own block,
so any join propagates our correct local orientation into a block whose global
polarity disagrees with it, and the 70 reads we attach score discordant. **The
veto has been preventing us from compounding an error we had already made.**

It also explains three earlier results that looked unrelated: the window's 77.65%
baseline (noted at the time as a pre-existing inconsistency, now root-caused),
why `--min-read-margin 1` and the unanimity skip both flipped exactly the same 70
reads, and why removing the repeat cluster changed nothing -- the defect was never
inside gap 2.

The priority therefore inverts. The bug worth fixing is the **tier-1** join at
`36,172,778-36,209,945` -- tier 1 being the clean-SNP-only tier, the one we trust
most -- not gap 2's confirmation rule. And the measurement that sizes it is the
read-level gate for `--recover-gaps` on versus off across the whole chromosome,
since a single wrong join costs more concordant reads than several correct ones
gain.

## Root cause: allele attachment on a gap that never joined

The section above pinned the corruption on recovery and on the tier-1 join in the
same window. Recovery is right; the tier-1 join is not. Bisecting the pipeline
settles it -- read labels dumped and scored at each stage:

| stage | window | key block |
|---|---|---|
| pre-recovery | 157/161 (**97.52%**) | `36209945`: 49 reads, 91.8% |
| after the gap threads | 278/360 (**77.22%**) | `36209945`: **164 reads, 50.6%** |
| after the parity edges | 278/360 (77.22%) | merged, 72.0% |

So the damage happens inside the gap threads, before any edge is applied, and the
join is not doing it: the joined gap's target block goes 44 -> 125 reads and stays
at **100.0%**. `apply_gap_phase_edges` is a parity union-find that relabels
uniformly per phase set and already guards self-edges, so it cannot corrupt half a
block either.

Attributing every write single-threaded (parallel threads interleave stderr and
lose lines -- an earlier count of 44 was an artifact of that) gives 199 new labels,
of which **155 come from `implied_assignment`**, the allele-agreement path behind
`--link-by-alleles`:

| gap | verdict | target | reads attached | accuracy vs truth |
|---|---|---|---:|---:|
| 36,172,778-36,209,945 | **joined** | 36168059 | 38 | **100.0%** |
| 36,247,421-36,268,291 | **never joins** | 36209945 | **115** | **67.0%** |
| 36,247,421-36,268,291 | never joins | 36268291 | 2 | 50.0% |

`implied_assignment` derives a read's haplotype from its own allele agreement with
the flank's live candidates and commits it with no flip, and it is a pure
unanimity test with **no minimum count**: a read that observes a single site of
that phase set and matches it gets a haplotype. On a joined gap that is sound,
because the flank's orientation has been confirmed from both sides. On gap 2,
which never joins, nothing confirms it -- and the 115 reads it wrote took flank
36209945 from 49 reads at 91.8% to 164 at 50.6%, chance.

**That is why the gap looked unjoinable.** Its own phasing is perfect (95 reads,
100.0%, both halves), its votes are sound, and every attempt to force the join --
`--min-read-margin 1`, skipping the confirmation, collapsing the repeat cluster --
flipped the same ~70 reads, because the block it was joining to had already been
filled with coin-flip haplotypes by the same gap's own one-sided attachment.

### The fix, and why it is opt-in

`--gap-allele-attach-join-only` requires a gap to have joined before allele
agreement may attach reads to its flanks. On this window: **77.65% -> 98.51%**
(265 concordant of 269 tagged, against 241/242 with recovery off, so recovery now
adds 24 correct reads instead of 79 wrong ones). Default behaviour is unchanged
byte-for-byte.

It is not a default because across the ten-window panel the restriction is a bad
trade overall (`implied_attach_panel.tsv`):

| | before | after |
|---|---:|---:|
| tagged | 4,687 | 3,576 |
| concordant | 4,576 | 3,555 |
| discordant | 111 | **21** |
| accuracy | 97.63% | **99.41%** |
| gate concordant -> discordant | -- | 2 |

It removes 90 discordant reads and 1,021 concordant ones. The allele path is a
net gain in most windows and harmful only where a one-sided link is unvalidated.

Read agreement count does **not** separate the two cases, which rules out the
obvious threshold: single-site attachments are 98.1% correct at chr20:7,073,919
and 74.8% at 36,217,274. The discriminator has to be a property of the flank being
attached to, not of the read. That is the next measurement: score attachments
against the polarity of their target block across the panel, split by whether the
flank carries committed read support of its own.

## Deficit gap 3: cannot be joined, but its interior was being thrown away

`chr20:35,919,404-36,156,319`, 236.9 kb. A different failure from gap 2: every
tier links **both** flanks (`L=1 R=1`) but to different proposal phase sets
(`35873360` and `36156319`), and joining requires one proposal block to hold both.

### Why no single block spans it

The proposal fragments into four to six blocks, all internally accurate. At tier 1
pass 0:

| proposal ps | reads | span | accuracy | reaches left / right edge |
|---|---:|---|---:|---|
| 35873360 | 299 | 35,849,159-35,961,779 | 100.0% | 63 / 0 |
| 36011708 | 58 | 35,998,813-36,041,729 | 100.0% | 0 / 0 |
| 36083127 | 30 | 36,076,496-36,123,378 | 96.7% | 0 / 0 |
| 36156319 | 114 | 36,134,775-36,194,144 | 100.0% | 0 / 71 |

The interior sites cover the gap edge to edge (58 sites with >= 10 observations:
11 clean SNPs, 2 clean indels, 36 MSA, 9 repeat), so this is not a site-discovery
failure. One spacing breaks the chain: **35,959,001 -> 35,981,762, 22.8 kb, with
zero reads covering both sites.** Every other large spacing carries 8-20 spanning
reads. Coverage is not the problem either -- 65-71x throughout, 168 reads inside
the interval -- the longest read in the region is 29.1 kb and none of them spans
from one flanking het site to the other.

So the gap genuinely cannot be closed by read linkage, and abstaining is correct.
No read-based phaser can cross that point.

### What was being lost

The output kept only the two flanks and discarded every interior block:

```
ps=35873360  379 reads  35,849,159-35,958,948   97.9%
ps=36156319  144 reads  36,131,191-36,172,639  100.0%
```

`emit_independent_gap_block` exists for exactly this case and was not firing. Its
own diagnostic (`--verbose 2`) shows why:

```
GapIndependentBlock  locally_phased=522  has_original=134
                     gap_only_ps_groups=5  chosen_ps=35902410  chosen_size=111
GapIndependentBlockApplied  chosen_size=111  applied=5
```

Five gap-only groups qualify and the function takes **only the largest**. That
choice is also the wrong one: the largest gap-only group is the flank-adjacent
component, whose reads this same recovery round has already attached, so it
applies 5 reads while four real interior blocks are discarded. Raising
`--gap-independent-min-reads` does not help (default 3; at 20 the window moves by
4 reads) because the limit was never the threshold.

### Fix: emit every group that clears min_reads

`min_reads` already screens the noise fragments the single-group rule was
guarding against. On this window:

| | before | after |
|---|---:|---:|
| tagged | 527 | 801 |
| concordant | 518 | **788** |
| accuracy | 98.29% | 98.38% |
| gate concordant -> discordant | -- | **0** |

Four interior blocks now emit: `35981805` 41 reads at 100.0%, `36011708` 91 at
98.9%, `36056601` 50 at 100.0%, `36077134` 92 at 96.7%.

Across the ten-window panel (`independent_blocks_panel.tsv`) it is a net gain
rather than a trade, which is why this one is a default and the allele-attachment
restriction is not:

| | before | after |
|---|---:|---:|
| blocks | 21 | 26 |
| tagged | 4,687 | 5,037 |
| concordant | 4,576 | **4,922** |
| discordant | 111 | 115 |
| accuracy | 97.63% | 97.72% |
| gate concordant -> discordant | -- | **0** |

394 reads newly tagged concordant against 9 discordant. Two windows carry the
change: this one and `13,429,829` (+49 tags, +46 concordant, one extra block).
`61,732,321` churns -- 48 concordant tags lost against a larger gain, net +30
concordant -- because emitting more interior blocks changes which reads its
homopolymer-tier join claims first; no read there flips concordant to discordant.

`emit_independent_gap_block` had no unit test at all; it now has one pinning both
halves of the contract (every qualifying group emitted, groups below `min_reads`
still screened).
