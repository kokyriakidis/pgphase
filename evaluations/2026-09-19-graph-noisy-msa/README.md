# The graph arm never ran stage 2, and why wiring it up is not enough

## The question

Whether the repeat-context het indels the graph arm discards are rescued by the
noisy MSA in the default BAM pipeline -- i.e. whether `msa_verified` reaches
them after all.

## Answer: yes in the alignment arm, never in the graph arm

`classify_cand_vars_pgphase` does not treat `RepeatHetIndel` as a verdict. It
treats the locus as a POINTER TO A NOISY REGION (`collect_var.cpp:1690-1694`):

```
if (c == VariantCategory::RepeatHetIndel) {
    if (pos in region) cr_add_var_to_noisy_cr(noisy_var_cr, low_comp, key, chunk, false, opts);
    continue;
}
```

The MSA then reconstructs the locus and the rebuilt site enters round 2 as
`NoisyCandHet`. Upstream longcallD does the same (`collect_var.c:2917-2972`:
`pre_process_noisy_regs` -> k-means on clean -> `collect_noisy_vars1` per noisy
region -> k-means over `LONGCALLD_CAND_GERMLINE_VAR_CATE`). Fetched
`origin/main` at 491f055 -- 5 commits ahead of the local checkout, all README /
makefile / usage; the pipeline is unchanged.

Measured on `chr20:22,930,000-23,060,000`, the same 130 kb through both arms:

| category | alignment arm | graph arm |
|---|---:|---:|
| CLEAN_HET_SNP | 53 | 55 |
| CLEAN_HET_INDEL | 4 | 8 |
| **NOISY_CAND_HET** | **17** | **0** |
| **REP_HET_INDEL** | **0** | **17** |

Exactly the same 17 loci. In the alignment arm `REP_HET_INDEL` is a transient
label that routes a site INTO the MSA; in the graph arm it is terminal. The
site blocking `chr20:22,980,600-23,008,891` is one of them: the alignment arm
emits `23,007,537 TA>T` as `NOISY_CAND_HET`, phased at `PS=22,930,405`, while
the graph arm leaves it `REP_HET_INDEL` at `PS=-1`.

This also corrects a claim made earlier today. "`msa_verified` cannot reach
these sites" was measured on graph-arm candidates only (all 808 carry 0) and
generalised wrongly: in the alignment arm the very same loci are rebuilt BY the
MSA and do carry it.

## Wiring stage 2 into the graph arm: three pieces, and a wall

Behind `--graph-noisy-msa` (default off):

1. **Reference slice.** `build_graph_chunk` never fills `chunk.ref_seq`, which
   the MSA reads. `seed_graph_noisy_regions` fills it from the slice the caller
   already fetched for the noise filter.
2. **Noisy regions.** Seeded from the demoted candidates with the same widening
   `cr_add_var_to_noisy_cr` applies (a locus inside a low-complexity tract takes
   the whole tract). 17 candidates -> 15 merged regions on the test window.
3. **Round 2 and metadata.** `collect_noisy_vars_step4`, then
   `assign_hap_based_on_germline_het_vars_kmeans(kCandGermlineVarCate)`, then
   `synthesize_meta_for_appended_candidates` -- the writer looks metadata up BY
   CANDIDATE INDEX and silently skips anything past the end of `site_meta`, so
   a reconstructed site would phase and never be emitted.

**It reconstructs nothing.** With 15 regions seeded and 73 candidates, the main
solve returns `noisy_het = 0`. The reason is one line in
`collect_noisy_reg_reads` (`collect_phase_noisy.cpp:1042`):

```
if (r.digars.empty()) return;  // graph-only reads have no sequences
```

A graph chunk's reads carry no per-read alignment stream, so every read is
skipped, the region has zero reads, and the MSA has nothing to align. The
in-chunk recovery's sub-solve DOES produce 4 reconstructed sites in the same
chunk -- because it re-reads the BAM through `process_chunk`, which builds
digars.

## What that implies for the design

Stage 2 cannot run on a graph chunk as it stands; it needs reads with
alignment detail. Two ways:

- build digars for graph reads at chunk build time -- chunk-wide cost for a
  pass that only touches repeat loci;
- build them only for reads overlapping the seeded noisy regions, which is what
  the in-chunk recovery already does for unphased windows, and it already
  yields `NoisyCandHet` sites there.

The second is the same shape as the machinery that exists, so the next step is
to feed the seeded noisy regions into the recovery's window list rather than to
duplicate digar construction.

Default path byte-identical with the flag off; unit 3/3, predicate 151/151,
window 125/125.

## A build lesson from this change

`0ac38e2` was committed with `make unit-tests` broken and the failure misread
as passing: the check counts `ALL PASS` lines, and a suite that fails to LINK
prints none, so "ALL PASS=0/3" looks like the same shape as a pass count. The
two new functions lived in `graph_bam_adapter.cpp` and call
`vcf_to_variant_key`, `populate_low_complexity_intervals` and
`variant_genomic_span`, which live in `collect_var.o` -- an object
`test_graph_bam_adapter` does not link. Moved both definitions to
`graph_collect.cpp`, the only caller and a translation unit that links the full
set. A zero from a counting check has to be read as a failure, not a number.
