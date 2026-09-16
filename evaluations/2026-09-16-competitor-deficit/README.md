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
