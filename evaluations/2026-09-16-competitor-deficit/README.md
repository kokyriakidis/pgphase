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
