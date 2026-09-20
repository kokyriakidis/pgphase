# The noisy-MSA branch gate: a real divergence, measured, and deliberately not shipped

## The bug, verified on both sides

`collect_noisy_reg_aln_strs` picks one of two paths per noisy region
(`align.cpp:1942-1952`):

```
if (ps > 0)                       -> phase-set-guided consensus
else if (n_full_reads >= min_depth) -> unguided MSA
                                   -> otherwise the region is abandoned
```

Upstream's unguided path bails only when there are NO fully-covering reads:
`wfa_collect_noisy_aln_str_no_ps_hap` returns early on `n_full_reads == 0`
(`align.c:1175`) and then on the longest read reaching `max_noisy_reg_len`
(`align.c:1177`). The second guard is already ours at `align.cpp:1626`; the
first is not -- we demand `min_depth` (5) fully-covering reads instead of 1.

That scales the requirement with region LENGTH, because a longer region is
spanned end to end by fewer reads. Measured on the 200 kb span around the
cluster:

| region | length | noisy reads | fully covering | `ps` | branch |
|---|---:|---:|---:|---|---|
| 30,802,901-30,804,117 | 1,217 bp | 16 | **3** | -1 | **abandoned** |
| 30,804,329-30,804,872 | 544 bp | 16 | 14 | 30,782,343 | guided, resolves |

Branch census over that span: 42 guided, 1 unguided, **20 abandoned**.

Everything around the gate is a faithful port and was checked line for line:
`collect_phase_set_with_both_haps` against `align.c:1230-1272` (same
full/partial accounting, same minor-haplotype checks, same thresholds 1 and 2),
the abPOA setup (`wb -1`, `inc_path_score 1`, `cons_algrm ABPOA_MF`,
`min_freq = min_af`), the scoring defaults (match 2, mismatch 6, gaps 6/2/24/1),
and the chunk size (500 kb). So `ps = -1` here on both sides, and upstream
therefore takes the unguided path with 3 reads where we abandon the region.

## What fixing it does

Changing the gate to `n_full_reads > 0`:

| | records | identical to upstream | upstream-only | ours-only |
|---|---:|---:|---:|---:|
| before | 116,896 | 114,407 | 3,863 | 2,489 |
| after | 121,800 | 114,433 | 3,837 | **7,367** |

In the target cluster it works as intended -- 17 records become 59, of which 42
are upstream's (25 recovered), leaving 26 of upstream's 68 still missing. But
chromosome-wide it buys 26 upstream-matching records for 4,878 that upstream
does not have, and read placement gets worse on both arms:

| arm | before | after |
|---|---|---|
| alignment | 216,962 tagged, 349 blocks, 1,516 misplaced, 0.699% | 216,967 tagged, 341 blocks, 1,524 misplaced, 0.702% |
| graph | 219,055 tagged, 326 blocks, 2,526 misplaced, 1.153% | 219,069 tagged, 327 blocks, **2,696** misplaced, **1.231%** |

## What decided it: the GIAB benchmark

Scored against `HG002_CHM13v2.0_v5.0q_smvar.vcf.gz` inside its own confident
BED (648 chr20 intervals, 58.9 Mb, 174,492 truth variants), exact
`(POS, REF, ALT)` matching:

| arm | TP | FP | FN | precision | recall | F1 |
|---|---:|---:|---:|---:|---:|---:|
| ours, before the gate fix | 91,894 | 2,964 | 3,166 | 0.9688 | 0.9667 | **0.9677** |
| ours, after the gate fix | 91,894 | 2,964 | 3,166 | 0.9688 | 0.9667 | **0.9677** |
| upstream longcallD | 91,478 | 3,458 | 3,582 | 0.9636 | 0.9623 | 0.9629 |

The fix is **exactly neutral** on the benchmark: identical TP, FP and FN, which
means all 4,904 extra records fall OUTSIDE the confident regions -- they are in
the low-mapping-quality and segmental-duplication territory the benchmark
excludes. So the change costs read-level accuracy and buys nothing measurable.

**It stays reverted**, recorded here as a known, deliberate divergence from
upstream with its cost priced.

Caveat: exact-match scoring without allele normalisation undercounts
representation-equivalent calls, so the absolute numbers are conservative. All
three arms were scored identically, so the comparison between them holds.

## The finding that matters more

On this benchmark we are already **better than upstream** -- F1 0.9677 against
0.9629, with 416 more true positives, 494 fewer false positives and 416 fewer
false negatives. Record-level parity with upstream is therefore the wrong
objective wherever it conflicts with accuracy, and this gate is precisely such a
place: matching upstream's rule here would add thousands of calls outside the
confident regions and misplace 170 more reads.
