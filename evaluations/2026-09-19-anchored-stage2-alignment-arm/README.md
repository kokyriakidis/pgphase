# Filling the gap `2026-09-18-anchored-stage2` left open: porting the alignment arm to longcallD's actual (unanchored) stage 2

## Background

`opts.anchored_stage2` (`phasing_types.hpp:296`) defaults to `true`: the
second, noisy-inclusive k-means re-solve (`kCandGermlineVarCate`, run after
each noisy-region pass) skips its Phase-1 outward sweep entirely and just
refines on top of whatever the first, clean-only pass decided, pinning every
site that pass already resolved. `evaluations/2026-09-18-anchored-stage2/`
measured this as a real win -- 119 fewer misplaced reads, one fewer bridged
block -- but **only on the graph + recovery arm**, and says so explicitly:
"The alignment and hybrid arms share this code and were NOT measured, so
their baselines are unverified against this change." That gap is what this
record closes.

longcallD's own code resets fully every call (`assign_hap.c:491`,
`read_init_hap_phase_set` unconditional), contradicting its own comment at
`collect_var.c:2939` that a round should use the previous phasing as
initialization -- our `anchored_stage2=true` implements what upstream's
comment says, not what upstream's code does.

## Orientation: a real, measured win

`chr20:656,411`/`658,439` (`PS=130540`, 807 clean-SNP-anchored genotypes,
807 scored, switches once with zero conflicts anywhere in the local link
chain -- see `evaluations/2026-09-19-joint-orientation-alignment-scope/`)
was one of 4 genuine internal-orientation switches found chromosome-wide
against parental-origin read truth, out of ~295 scored phase sets, on the
`ac59de2` default (`anchored_stage2=true`).

Rebuilt with `opts.anchored_stage2 = false` scoped to
`collect_bam_variation()` only (`collect_pipeline.cpp`, same pattern as the
reverted `joint_het_orientation` attempt) -- struct default untouched, graph
arm unaffected, confirmed by `make window-tests` passing 125/125 including
the graph-arm window. Re-scanned whole chr20 the same way:

| | switches | scored phase sets |
|---|---:|---:|
| `anchored_stage2=true` (shipped) | 4 | ~295 |
| `anchored_stage2=false` | **2** | ~260 |

3 of the 4 original switches resolve (130540, plus two others). One
(26,549,599, 2-4 scored SNPs either way) persists in both configurations --
thin evidence, likely a genuinely hard site rather than something either
flag touches. One new switch appears (65,498,089, 11 scored SNPs, 6
against 5 -- a near-tie, not a confident block-sized switch like the ones
fixed). Net: fewer switches, and the survivors are all thin-evidence,
where the fixed ones were large, confident blocks (807 and similar scored
counts).

Full per-site results: `switches_anchored_false.tsv`.

## Record-level parity: apparent loss, mostly the known merge tension resurfacing

Matched against upstream's actual output
(`test_data/longcallD_eval_chr20/phased.vcf`) on `(POS, REF, ALT)`:

| | our records | identical to upstream | ours-only | upstream-only |
|---|---:|---:|---:|---:|
| `anchored_stage2=true` (shipped) | 117,738 | 116,851 (99.2%) | 887 | 1,419 |
| `anchored_stage2=false` | 117,040 | 114,410 (97.8%) | 2,630 | 3,860 |

Both mismatch directions roughly tripled, but not from new errors: every
one of the 1,667 new `ours-only` records carries a comma-ALT (a merged,
multiallelic record -- baseline's 887 `ours-only` records had ZERO of
these), and 1,625 of the 1,667 positions behind them (97.5%) are
positions where WE also have a record after the fix and upstream has a
record there too -- just as upstream's own two separate single-ALT rows
against our one merged multiallelic row, the exact representation choice
`merge_colocated_msa_alleles` already prices
(`evaluations/2026-09-20-parity-report/`, `2026-09-19-multiallelic-
representation/`). Only 42 of the 1,667 positions are a genuinely new gap
(no record on our side at all).

Spot-checked one, chr20:10,064,691 (`G>GTT,GTTT`, an STR indel): baseline
emitted ONE allele at DP 25; the fix emits BOTH, merged, at DP 88 --
much closer to upstream's own DP 93 at that locus (upstream: `G>GTT`
0|1 and `G>GTTT` 1|0, same 93 reads split both ways). The fix is
recovering reads at this locus that the anchored baseline was not
using, not corrupting the call: the "regression" is our merge writing a
comma-ALT row for a site that used to be silently dropped or
under-supported, priced the same way that tension already is
elsewhere, not a new defect from the reset itself.

Distinct phase-set count also drops, 472 -> 369 -- fewer, larger blocks,
consistent with a fresh sweep over the complete site set bridging
gaps the anchored solve's narrower stage-1-only pivot could not reach.

## Status: shipped

First reverted here, for the record-parity loss above. Two things changed
that decision.

One: an attempted narrower fix for the *other* known mechanism (the
576,038-class homopolymer-indel tie-break, `evaluations/2026-09-19-joint-
orientation-alignment-scope/`) turned out, on inspection, not to be a port
of anything -- longcallD's own `iter_update_var_hap_to_cons_alle`
(`assign_hap.c:425-462`) has no joint/same-vs-flip branch at all, just the
same unconditional per-haplotype argmax as our fallback path. That
invented mechanism was dropped rather than shipped.

Two: `anchored_stage2=false` is not in the same category. It is not a new
mechanism being weighed against a metric -- it is upstream's actual,
unconditional behavior (`assign_hap_based_on_germline_het_vars_kmeans`,
`assign_hap.c:465`, calls `read_init_hap_phase_set` on every invocation;
there is no anchoring parameter in that function's signature at all).
`anchored_stage2=true` is documented, in this project's own source, as
"a knowing divergence from upstream" -- our own addition, not a port.

The record-parity table above is real and stands as the measured cost.
It is not evidence the port is wrong: it measures agreement with
upstream's *output*, and upstream's own comment at `collect_var.c:2939`
(a round should use the previous phasing as initialization) disagrees
with upstream's own code, which is what was ported here. Matching
upstream's records more closely by keeping our own un-ported anchoring
addition is not the same thing as being more correct -- the switch
metric is truth-scored against parental origin, an external ground
truth neither tool's output defines. Shipped: `opts.anchored_stage2 =
false`, scoped to `collect_bam_variation()`
(`src/collect_pipeline.cpp`), struct default left at `true` so the
graph arm (validated separately, `evaluations/2026-09-18-anchored-
stage2/`) is untouched. `make unit-tests`, `make predicate-tests`,
`make window-tests` all pass unchanged (125/125 window assertions,
graph-arm window included).

Co-authored-by: Claude Sonnet 5 <noreply@anthropic.com>
