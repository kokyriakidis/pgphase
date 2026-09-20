# Traced a targeted mismatch to a genuine mechanism; the fix is correct locally and a net loss chromosome-wide

## The trace

Targeted one of the ~1,000 remaining record mismatches: chr20:13,920,218.
We call `CTT>C` (a 2bp deletion) `1|1` (DP 54: 24 ref / 30 alt). longcallD
calls `CTTT>C` (a 3bp deletion at the same anchor) `1|0` (DP 66: 33/33).
Two red flags together: a "homozygous" call with 44% of reads showing
REF, and a different deletion length than upstream at the same anchor.

Instrumented `hap1_vars`/`hap2_vars` directly (the two haplotype
consensuses' own independently-discovered candidates, before they get
merged): hap1's own consensus shows a 2bp deletion at this position,
hap2's shows a 3bp deletion at the *same* position. `exact_comp_var_site`
compares `(pos, type, ref_len)`; a 2bp and 3bp deletion at the same pos
don't match, so they should become two independent het records -- but the
actual output shows one, homozygous. The mechanism: this project already
has code for exactly this case, `split_nested_msa_deletions` (decomposes
two same-anchor, different-length deletions into a common prefix shared
by both haplotypes plus a residual private to the longer one), gated on
`opts.force_noisy_msa`. Grepped: `force_noisy_msa` is never assigned
`true` anywhere in the current source -- it was wired only to
`--recover-gaps`'s retry path, which this codebase no longer has. The
split has been completely unreachable, in every configuration, since
that removal.

## The fix, verified locally

Two changes, mirroring the pattern of the two previously-shipped
CLI-scoped fixes: `opts.force_noisy_msa = true` in
`collect_bam_variation()` (alignment path only), keeping the existing
`if (opts.force_noisy_msa)` gate in `collect_phase_noisy.cpp` rather than
removing it -- removing it broke a graph-arm window test that expects the
*undecomposed* nested forms (`GTTAT>G` / `GTTATTTAT>G`) to exist as
candidates, because graph-site injection matches specific allele forms
against the graph's own snarl catalog. Scoping to the alignment path
fixed that regression; `make window-tests` passes 125/125 with the fix
scoped this way.

At chr20:13,920,219 the fix works exactly as designed: the locus splits
into a common 2bp deletion (DP 66, 3 ref / 63 alt -- correctly
homozygous) plus a private 1bp residual (DP 66, 33 ref / 33 alt --
correctly the heterozygous signal, matching longcallD's AD field exactly:
33,33). The spurious "1|1 with 44% ref" call is gone.

## Whole chr20: a net loss on both metrics that matter

| | record identity | read accuracy |
|---|---:|---:|
| before (3 shipped fixes) | 99.13% (117,009/118,038) | 99.27% (0.73% disc.) |
| **+ this fix** | **98.80%** (116,582/118,400) | **99.21%** (0.79% disc.) |

Both got worse, not just one. Two separate reasons:

1. **Representation, not just genotype.** Even where the split fixes the
   genotype correctly (as at 13,920,219), it produces `(pos, ref, alt)`
   text that does not match longcallD's own record -- longcallD emits the
   whole 3bp deletion as one record on one haplotype; we now emit a
   decomposed common+residual pair. Neither of the two decomposed records
   is byte-identical to upstream's single one. Fixing the genotype and
   matching the record text turned out to be different problems here.
2. **The split fires far more often than the one locus checked.** Both
   `ours-only` and `upstream-only` mismatch counts rose (789 and 427
   more respectively) -- consistent with many other nested-deletion loci
   where the decomposition is not as cleanly beneficial as the one
   example traced by hand, and some of it evidently introduces new
   disagreement rather than resolving it.

## Status: reverted

`git diff` shows only the three previously-shipped fixes
(`anchored_stage2`, `merge_colocated_msa_alleles`, the per-candidate
`refresh_assigned_msa_observations` fix). Same pattern as
`2026-09-20-refresh-observations-bug/`'s first attempt: a real,
root-caused mechanism, individually verified correct at the locus that
motivated it, net-negative when measured chromosome-wide. The residual
half of the split (the private 1bp deletion at 13,920,221) also never
resolves a genotype on its own -- `is_homopolymer_indel` excludes it from
scoring the same way chr20:411,654 was excluded before that fix, and its
`HAP_ALT`/`HAP_REF` stayed `0/0` even with the split enabled -- so even
the one locus this was measured against is only partially fixed by this
change alone.

Co-authored-by: Claude Sonnet 5 <noreply@anthropic.com>
