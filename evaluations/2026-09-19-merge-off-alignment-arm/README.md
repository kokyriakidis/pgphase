# merge_colocated_msa_alleles off, on top of the anchored_stage2 port: the old accuracy cost was mostly the same unported divergence

## Background

`merge_colocated_msa_alleles` (`phasing_types.hpp:243`, default `true`) is
a real, acknowledged non-port: longcallD never merges a locus where both
haplotypes carry a different ALT into one multiallelic record
(`collect_var.c:1329-1336` keeps the older candidate and frees the MSA's;
there is no merge step on the upstream side). `evaluations/2026-09-20-
parity-report/`'s own closing caveat: "it exists so the comparison...can
be run, not as a configuration knob, and by the standing instruction it
should be deleted rather than kept once the representation question is
settled." That report measured turning it off as roughly doubling
misplaced reads (0.699% -> 1.385%) and kept it on for that reason.

That measurement predates this session's `anchored_stage2` port
(`evaluations/2026-09-19-anchored-stage2-alignment-arm/`). Re-measured
with the port already in place, the picture changes.

## Record-level parity, with the port fix already shipped

Matched on `(POS, REF, ALT)` against longcallD's actual output
(`test_data/longcallD_eval_chr20/phased.vcf`):

| | our records | identical | ours-only | upstream-only |
|---|---:|---:|---:|---:|
| anchored_stage2=false, merge on (prior state) | 117,040 | 114,410 (97.8%) | 2,630 | 3,860 |
| **anchored_stage2=false, merge off** | 117,896 | **116,868 (99.1%)** | 1,028 | 1,402 |
| original baseline (merge on, anchored on) | 117,738 | 116,851 (99.2%) | 887 | 1,419 |

Turning merge off recovers almost all the record-identity the anchored
fix alone had cost, and lands slightly ahead of even the pre-session
baseline.

## Read accuracy, re-measured against the diplinator truth BAM

`scripts/evaluate_phase_accuracy.py`, whole chr20, both builds with
`anchored_stage2=false`:

| | accuracy | discordance | phase sets |
|---|---:|---:|---:|
| merge on | 99.34% | 0.66% (1,425/216,942) | 323 |
| merge off | 99.28% | 0.72% (1,570/217,638) | 361 |

Not a doubling. A 0.06 percentage point move, 145 more discordant reads
for 1,458 more records matching upstream exactly. The original 0.699% ->
1.385% doubling was measured against the `anchored_stage2=true` baseline;
most of that cost was the same read-partition instability the anchoring
port already fixes elsewhere in this session's work, not the merge
representation choice itself.

## Status: shipped

`opts.merge_colocated_msa_alleles = false`, scoped to
`collect_bam_variation()` (`src/collect_pipeline.cpp`), alongside
`anchored_stage2 = false`. Struct default left at `true` -- the graph arm
never sets this field itself and was previously untouched; confirmed by
`make window-tests` (125/125, including the graph-arm windows) after
fixing a testing-only mistake (the struct default was left flipped from
an earlier manual experiment and briefly leaked into a build before this
fix's own scoped-only version was rebuilt and re-verified).

Unit, predicate and window suites all pass unchanged.

Co-authored-by: Claude Sonnet 5 <noreply@anthropic.com>
