# The BAM pipeline as a faithful port: five divergences removed, and what each cost

The alignment path is a port of longcallD, so every place it behaves
differently is a defect unless it is a knowing, measured choice. Five such
places existed. All five are now off for `collect_bam_variation`, each behind a
named option so the graph arm -- which is ours, not a port -- keeps its own
behaviour.

| # | option | what we did | what longcallD does |
|---|---|---|---|
| 1 | `merge_colocated_msa_alleles` | merge a locus whose two haplotypes carry different ALTs into one multiallelic record | never merges; emits two complementary biallelic rows (`collect_var.c:1523-1559`, `vcf_utils.c:156-159`) |
| 2 | `refresh_msa_observations` | re-derive a two-cluster MSA candidate's counts and per-read profile with `call_local_msa_allele` after the MSA already built them | never revisits them; a read from the other cluster is recorded as reference outright (`collect_var.c:2205`), which is what makes its two rows complementary by construction |
| 3 | `add_unplaced_msa_observations` | add observations from reads the MSA could not place to candidates that have none | depth at a noisy candidate is exactly the reads its two cluster alignments cover |
| 4 | `phase_set_scoped_clean_rounds` | in the clean rounds, re-score a read once per phase set it spans and update each block with that block's own verdict, taking `chunk.haps` from the first | one hap per read, applied to every variant it covers (`assign_hap.c:437-447`); no phase-set scoping anywhere |
| 5 | `anchored_stage2` | stage 2 pins stage 1's read labels and per-site consensus | `assign_hap_based_on_germline_het_vars_kmeans` has no anchoring parameter; stage 2 resets and re-sweeps from a fresh pivot |

## What each is worth

Whole chr20, alignment arm, starting from all five faithful and restoring one
divergence at a time. Identity is against upstream's own `native.vcf`
(`##source=longcallD version=0.0.11-23e369d`) on `(POS, REF, ALT)` over phased
records; read metrics are truth-scored per phase set; F1 is GIAB
`HG002_CHM13v2.0_v5.0q_smvar` inside its confident BED (95,060 chr20 truth
variants, exact-match).

| configuration | records | identity | misplaced | hamming | blocks | GIAB F1 |
|---|---:|---:|---:|---:|---:|---:|
| **all faithful (now shipped)** | 118,273 | **99.69%** | 3,370 | 1.559% | 381 | 0.9632 |
| + `refresh_msa_observations` | 118,002 | 99.11% | 1,513 | 0.705% | 405 | 0.9628 |
| + `add_unplaced_msa_observations` | 118,114 | 99.55% | 3,685 | 1.676% | 325 | 0.9637 |
| + `phase_set_scoped_clean_rounds` | 118,311 | 99.64% | 3,343 | 1.546% | 380 | 0.9632 |
| + `anchored_stage2` | 117,209 | 98.80% | 2,092 | 0.968% | 387 | 0.9632 |
| all ours (previous default) | 116,896 | 96.73% | 1,516 | 0.699% | 349 | 0.9677 |
| upstream longcallD | 118,270 | -- | -- | -- | -- | 0.9629 |

Upstream-only records fall from 3,863 to 370 and ours-only from 2,489 to 373.

## The cost, stated plainly

Parity is not free on this data. The faithful configuration gives up
**0.0045 GIAB F1** (0.9677 -> 0.9632, 415 true positives and 451 false
positives) and **1,854 read placements** (1,516 -> 3,370 misplaced, 0.699% ->
1.559%), and it lands within 0.0003 F1 of upstream itself. In other words the
five divergences were the entire margin by which this pipeline beat the tool it
ports, and removing them removes that margin along with the divergence.

Two further readings from the table, both useful later:
- `refresh_msa_observations` is the load-bearing one for read placement: alone
  it takes misplaced reads 3,370 -> 1,513 while still holding 99.11% identity.
- `phase_set_scoped_clean_rounds` is nearly inert either way (1.546% against
  1.559%, identity 99.64% against 99.69%), so removing it costs nothing and
  simplifies the clean rounds.

## One fix that is not a divergence removal

`refresh_msa_observations` had a real defect of its own, fixed here rather than
deleted with the flag, because the graph arm still uses that path:
`call_local_msa_allele` attributes a read only when BOTH haplotype consensuses
yield a consistent allele, which fails when a homopolymer is unstable enough
that the two consensuses disagree on event TYPE (one insertion, one deletion at
the same locus). Every read then failed classification for that variant and its
counts were wiped. A pre-pass now records, per candidate, whether any read
produced a valid call, and the refresh skips candidates with none. Graph arm,
whole chr20: 62,458 -> 62,488 records, hamming 1.153% -> 1.154%.

## Scope

The five overrides sit in `collect_bam_variation` (`collect_pipeline.cpp`), so
`collect-graph-variation` is unaffected except through the pre-pass above. The
option defaults in `phasing_types.hpp` remain our own behaviour, which is what
the graph arm reads.
