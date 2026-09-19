# A predicate that was dead everywhere, and what it costs to wake

## The question

Whether BAM discovery supplies the sites needed to phase like the competitor.
On `chr20:22,980,600-23,008,891` it does -- and they are still not used.

## The site is present, observed, categorised and phased

The alignment arm holds `23,007,537 TA>T` as `NOISY_CAND_HET`, DP 67, 27 ref /
40 alt, AF 0.597, with a phase set. The pipeline's own phase-matrix dump shows
it co-observed with the next clean SNP `23,008,891` on **59 reads: 7 agree, 52
conflict** -- almost exactly the read-truth table (8 / 49). The evidence for the
join is in the pipeline, not merely in the BAM.

The linker never uses it. With `--verbose 2` the walk reports `23008891 0 0`:
no link found at all, while its neighbour 305 bp further on reports `57 0`.

## Why: the site is excluded from the LINK LIST

`iter_update_var_hap_cons_phase_set` (`collect_phase.cpp:561-563`) keeps a
homopolymer indel out of the link list unless its own allele depths call it a
clear heterozygote:

```
const bool hp_indel_blocks_link =
    var.is_homopolymer_indel && !allele_depths_call_het(var, opts);
```

Probing every clause of the admission gate for this site:

```
HETGATE pos=23007537 cons=[0,1] hp=1 depths_het=0 msa_alts=0 gap_link=0 two_allele=0
```

`cons=[0,1]` -- both haplotype consensus alleles set and different. The only
failing clause is `depths_het`, on a site at 27/40 and AF 0.597.

## Why that: the predicate is dead in every run

`allele_depths_call_het` opens with

```
if (opts.retry_windows.empty() && !opts.joint_het_orientation) return false;
```

and **nothing in the tree sets either**. `joint_het_orientation` is declared and
never assigned; `retry_windows` was filled by the post-hoc retry path, removed
in the one-arm purge. So the predicate returns false on every call, which
silently disabled both of its consumers:

- the link-list admission above, so EVERY homopolymer indel is excluded from
  linking regardless of its depths;
- the orientation branch at `collect_phase.cpp:833`, which keeps a genuine het
  from collapsing to `1|1` -- the failure its comment records at
  chr20:48,204,383 (AT>A, 30 ref / 41 alt, the only het in a 41.8 kb stretch
  and the site hiphase bridges with).

Both were written with measurements attached. Both were inert.

## Waking it chromosome-wide is wrong; scoping it is right

| arm | tagged | misplaced | hamming | VCF blocks | N50 |
|---|---:|---:|---:|---:|---:|
| HEAD default | 219,061 | 2,543 | 1.161% | 333 | 486 kb |
| escape revived chromosome-wide | 219,154 | 9,866 | **4.502%** | 526 | 466 kb |
| **escape scoped to recovery windows** | 219,055 | **2,526** | **1.153%** | 339 | 486 kb |
| scoped + `--graph-noisy-msa` | 225,645 | 9,451 | 4.188% | 350 | 503 kb |

Removing the gate outright costs 3.3 pp of read accuracy and 193 blocks, and
the window suite catches it. The predicate was tuned for the scoped use its
name suggests -- inside a window the first pass could not phase -- so the
recovery now fills `sub.retry_windows` with the windows it is solving. That is
a small net win on the default path: 17 fewer misplaced reads for 6 more
blocks, N50 unchanged.

The lesson for the removal that caused it: deleting the option that FILLS a
list silently changes the behaviour of every predicate that READS it. A grep
for the remaining readers would have shown two live call sites going quiet.
