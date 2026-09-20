# The MSA observation refresh was overwriting the one mechanism that makes upstream's co-located rows complementary

## The class this started from

Of the ~2,290 record-level mismatches against longcallD on chr20, the single
largest systematic class was two-row loci: 751 upstream-only records sat at
positions where longcallD emits **two** rows (both haplotypes differ from
reference, complementary `1|0` and `0|1`) and we emit one, or none.

chr20:155,439, upstream:

```
155439  CT   C   1|0:61:37,24:0.393:60:130540
155439  CTT  C   0|1:61:38,23:0.377:60:130540
```

ours, before this change:

```
155439  CTT  C   1|1:34:3,31:0.912:12
```

One row instead of two, wrongly homozygous, at DP 34 against upstream's 61.

## Ruling out the port

Every function on the path was diffed against upstream and is faithful:

- `is_match_aln_str_del` — identical walk, identical `started_check_del` /
  `n_non_del` logic, identical returns. Upstream's extra `float sim_thres`
  parameter is not referenced anywhere in its own body.
- `update_cand_var_profile_from_cons_aln_str21` (collect_var.c:2178) — line
  for line, including the `delta_ref_alt` running adjustment and the
  `var_beg_in_ref_str - delta_ref_alt` argument.
- `update_cand_var_profile_from_cons_aln_str2` (collect_var.c:2225) — the
  hap1/hap2 merge, `var_from_cons_idx` 1/2/3, and the
  `NOISY_CAND_HET` / `NOISY_CAND_HOM` split.
- `exact_comp_var_site` — ours calls it directly where upstream calls
  `exact_comp_cand_var`, which is a two-line wrapper around the same function.
  The `VariantType` enum is defined as the BAM opcodes (`Snp = 8`,
  `Insertion = 1`, `Deletion = 2`), so the type ordering matches too.

## Where upstream's complementarity actually comes from

`update_cand_var_profile_from_cons_aln_str21` scores a read against a candidate
in one of two ways, chosen by whether the candidate came from that read's own
cluster:

```c
if (var_from_cons_idx[i] & clu_idx) {          // candidate is from THIS cluster
    allele_i = get_var_allele_i_from_cons_aln_str(...);   // measure it
} else {                                       // candidate is the other haplotype's
    ... full_cover = get_full_cover_...(...);  // only check coverage
    allele_i = 0;                              // and record REFERENCE outright
}
```

That `allele_i = 0` is the whole mechanism. A read in cluster 2 is *declared*
reference at every candidate discovered from cluster 1, and vice versa, so the
two records at one locus come out anti-correlated by construction — which is
why longcallD scores **zero** loci on chr20 where one haplotype carries two
different alleles, and why its AD fields at 155,439 read 37,24 and 38,23
rather than both counting all 47 deletion-bearing reads.

## The bug

`refresh_assigned_msa_observations` (`collect_phase_noisy.cpp`, called from
`make_vars_from_msa_cons_aln`'s two-cluster branch) re-derives every
`NoisyCandHet` candidate's per-read profile and aggregate counts from a
second, independent classifier, `call_local_msa_allele`.

Neither function exists upstream. Grepping longcallD's whole source for
`refresh`, `local_msa_allele`, `reclassif` returns nothing:
`update_cand_var_profile_from_cons_aln_str21` is upstream's only producer of
these observations and upstream never revisits them.

`call_local_msa_allele` does not know about cluster membership. It re-measures
the other cluster's reads against the candidate and, at a nested deletion,
calls many of them ALT at the shorter allele too — because the read's deletion
genuinely does span that reference position. So it overwrites the deliberate
`allele_i = 0` with a measurement, the two records stop being complementary,
and the pair collapses into one wrongly-homozygous row.

## Which half does the damage

The refresh has two separable halves: a per-read profile rewrite and an
aggregate DP/AD rebuild. They were measured apart, whole chr20:

| | record identity | two-row class | 155,439 |
|---|---:|---:|---|
| shipped baseline (refresh on) | 99.13% | 751 | one row, `1\|1`, DP 34 |
| aggregate rebuild off only | 99.11% | 705 | one row, `1\|1`, DP 69 |
| **whole refresh off** | **99.60%** | **87** | **both rows, `1\|0` + `0\|1`, PS 130540** |

Dropping only the count rebuild fixes the depth (34 → 69, the real coverage)
and changes essentially nothing else — 705 of 751 still missing. The parity
damage is entirely in the per-read rewrite, which is the half that overwrites
`allele_i = 0`. Both halves therefore come out together.

With it off, 155,439 matches upstream on POS, REF, ALT, GT **and phase set**
(130540) for both rows.

## Whole chr20

Record-level parity, matched on (POS, REF, ALT) against
`test_data/longcallD_eval_chr20/phased.vcf`:

| | our records | identical | ours-only | upstream-only |
|---|---:|---:|---:|---:|
| shipped baseline | 118,038 | 117,009 (99.13%) | 1,029 | 1,261 |
| **refresh off** | 118,152 | **117,677 (99.60%)** | **475** | **593** |

Both mismatch buckets more than halve; total mismatches 2,290 → 1,068.

Two secondary signals confirm the mechanism rather than just the outcome:

- Median `NOISY_CAND_HET` depth over all 31,853 such candidates rises from
  **52 to 61**. The refresh was suppressing depth on every two-cluster MSA
  candidate chromosome-wide, not only at nested loci, because the rebuild
  counts only reads `call_local_msa_allele` can attribute.
- The other dominant ours-only class — SNPs upstream never calls — falls from
  845 to 398. Those spurious records had median DP **11** against 52 for
  `NOISY_CAND_HET` as a whole: they were precisely the badly-undercounted ones.

## The read-accuracy move is convergence, not regression

Against the diplinator truth BAM, whole chr20:

| | accuracy | discordant | phased reads |
|---|---:|---:|---:|
| shipped baseline | 99.27% | 1,595 | 217,741 |
| **refresh off** | **98.34%** | 3,658 | 219,851 |
| **longcallD itself** | **97.03%** | 6,506 | 219,090 |

Read in isolation this looks like the same trap that reverted
`evaluations/2026-09-20-refresh-observations-bug/`'s fix 2 and the nested
deletion split. It is not. longcallD's own score on the same truth BAM is
97.03% — the refresh was making this arm *better than the tool it ports*,
through a mechanism that tool does not have, while breaking the record parity
the arm exists to hold. Turning it off moves accuracy toward upstream and
still leaves us 1.3 points ahead of it, and moves phased-read count toward
upstream's too (217,741 → 219,851 against upstream's 219,090).

The loss is also concentrated, not diffuse: four phase sets carry 1,708 of the
3,658 discordant reads, against a worst-case of 185 in the baseline. It is not
over-merging — `PS=30782288` has the *same* 3,545 reads over the same span in
both runs, and the VCF block boundaries are unchanged and still match
upstream's own break at 31,736,302. It is an internal orientation flip inside
an identical block, which is a separate problem from record identity and is
the natural next thing to trace.

## Status: shipped

`opts.refresh_msa_observations = false`, scoped to `collect_bam_variation()`
(`src/collect_pipeline.cpp`), alongside `anchored_stage2 = false` and
`merge_colocated_msa_alleles = false`. Struct default left at `true`; the graph
arm never sets the field and is unaffected, confirmed by `make window-tests`
(125/125, graph-arm windows included). Unit (ALL PASS) and predicate (151/151)
suites pass unchanged. The final build was re-run on whole chr20 and produces a
byte-identical VCF to the measurement above.

## Noted while auditing, not acted on

Diffing every function in `collect_phase_noisy.cpp` against longcallD's source
turned up the other live mechanisms with no upstream counterpart:
`add_msa_site_observations` (adds observations from reads the MSA could not
place — the likely reason our DP at 155,439 is 80/69 where upstream's is 61/61,
since upstream counts only the two clusters' reads) and
`make_colocated_deletions_exclusive` (forces a read that is ALT at two nested
deletion records to be reference at the shorter one). Both are unconditional in
the alignment arm. `backfill_msa_observations` is declared and defined but
never called from anywhere — dead code, the same pattern as
`split_nested_msa_deletions` after `--recover-gaps` was removed.

Co-authored-by: Claude Opus 5 (1M context) <noreply@anthropic.com>
