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

## Routing the seeded regions through the recovery: it works, and it regresses

The MSA cannot run on the graph chunk, but the recovery's sub-solve re-reads the
BAM through `process_chunk`, which builds digars -- and `process_chunk` consults
NO catalog (no `load_sites_for_region`, no `inject_graph_sites`), so everything
it finds in a gap window is BAM-discovered. Routing the seeded noisy regions
into `retry_unphased_windows_in_place`'s window list therefore gives exactly the
alignment pipeline's stage 2, on BAM evidence only.

Two fixes were needed to make it land:

1. **The sub-solve's verdict was discarded for loci the parent already held.**
   The merge skips any key in `parent_keys`, transferring only read
   observations, so a parent `REP_HET_INDEL` kept a label derived from the
   reference context while the sub-solve had just answered the same question
   with reads. It now adopts the sub-solve's category for demoted loci.
2. **The adoption lookup silently found nothing.** A graph-derived candidate
   stores the catalog's ALLELE WALK in `key.alt`
   (`">115859261>115859263"`), not a sequence, so the alignment's key never
   matches the raw parent key -- only the translated VCF form in `site_meta`,
   which `parent_seq_index` records. Looking only at raw keys matched zero
   times.

On the motivating window this is exactly the intended effect: `23,007,537 TA>T`
is admitted, phased, and pulls `23,008,891` and `23,009,196` into a block with
it; repeat-class candidates drop 17 -> 12 and records rise 63 -> 70.

Whole chr20, against the shipped arm at 1.161% / 333 VCF blocks:

| arm | tagged | misplaced | hamming | VCF blocks |
|---|---:|---:|---:|---:|
| off (shipped) | 219,061 | 2,543 | **1.161%** | 333 |
| stage 2, sub-solve category adopted | 225,650 | 9,471 | 4.197% | 343 |
| stage 2, forced CLEAN adoption | 221,702 | 5,253 | 2.369% | 1,689 |
| stage 2, **no second round anywhere** | 219,533 | 2,468 | **1.124%** | **1,655** |

The third row was initially mislabelled as "without the added round 2". It is
not: the ablation guard matched SIX k-means call sites, including the two the
pipeline already shipped, so that arm disables the second round entirely. It is
the only arm that improves accuracy (75 fewer misplaced reads, +472 tagged) and
it fragments the VCF fivefold while the READ blocks barely move (323 -> 371) --
the admitted sites keep the sub-solve's phase-set ids and are never reconciled
with the parent's blocks.

So every configuration that lets these loci into a chunk-wide solve costs 2-4x
the read accuracy, exactly as the historical `--retry-unphased-with-bam`
measurement did (0.559% -> 2.723% for +6,669 reads; here +6,589 reads for
+6,928 misplaced). The recovery mechanism is now correct and the sites are
reachable; what is missing is that admitting them re-solves the chunk.

Off by default, flag-off byte-identical over whole chr20.
