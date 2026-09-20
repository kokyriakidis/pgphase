# Function-by-function audit against longcallD: what was removed, and what cannot be

Every function that carries upstream's name was extracted from both sources with
brace-balanced bounds and compared body-for-body -- 26 pairs with an extractable
definition on both sides. The audit reported 47 differences across 18 of them.
Each was then verified by hand against the upstream line before anything moved,
because this kind of sweep produces false positives: it flagged the
"both consensus alleles unset -> return 0" early return in
`read_to_cons_allele_score` as an addition when upstream has it at
`assign_hap.c:139`.

## Removed: four behavioural divergences in read scoring

All four are restrictions this project added, now off for the BAM path behind
`Options::upstream_read_scoring`:

| divergence | upstream |
|---|---|
| a homopolymer indel was spared the skip when `hp_gap_scorable` was set | skipped unconditionally, `assign_hap.c:166` |
| the clean agree/conflict tallies counted only a clean het SNP with at most two alleles | any clean SNP, heterozygous or homozygous, `assign_hap.c:174` |
| a clean het SNP or indel weighed 2 only with at most two alleles | weighs 2 regardless, `assign_hap.c:130-131` |
| a read took a phase set only from a site that was not homopolymer, not noisy-hom, not ungap-linked, and whose own allele matched a consensus | the first heterozygous site it covers, `assign_hap.c:328-336` |

Whole chr20, alignment arm: records 118,273 and identity against upstream's own
output unchanged at 117,900 of 118,270 (99.69%), upstream-only 370, ours-only
373, GIAB F1 unchanged at 0.9632. Read placement improves slightly -- misplaced
3,370 -> 3,226, hamming 1.559% -> 1.523% -- on 4,367 fewer tagged reads
(216,208 -> 211,841) across 389 blocks rather than 381.

## Restored: two sanity checks upstream carries

`wfa_trim_aln_str` has `assert(query_end <= target_end)` and
`assert(query_start >= target_start)`. Both are back, and neither fires on whole
chr20 in either arm.

## Attempted and reverted: the guards are not inventions

The audit's largest group was defensive guards absent upstream, which looked
like obvious deletions. Two were tried and both broke:

**`get_var_init_max_cov_allele`'s `!alle_covs.empty()` branch.** Upstream loops
`n_uniq_alles` over `alle_covs` with no guard because it always fills that
array. We do not: over `chr20:5-6 Mb`, **662 candidates carry
`n_uniq_alles = 2` with `alle_covs` EMPTY**, every one of them `CLEAN_HET_SNP`,
because a BAM-discovered candidate carries its depths in `ref_cov`/`alt_cov`
instead. Removing the branch indexes `alle_covs[0]` on an empty vector and the
run ends in `Segmentation fault (core dumped)`. The second branch -- deriving
the allele from `ref_cov`/`alt_cov` -- is not an invention layered on the port;
it is how those candidates are represented here at all.

**The `start_var_idx < 0` guards** in `init_assign_read_hap_based_on_cons_alle`
and `update_var_hap_profile_based_on_read_hap`. Removing them failed four graph
windows (`test_gap_windows.cpp:689`, windows 3,852,321 / 4,766,928 / 5,309,406 /
22,980,600). Graph-only reads legitimately carry no variant profile, and these
functions are shared with that arm, so the guard is what keeps a shared code
path correct for a caller upstream does not have. Upstream's own
`update_read_phase_set` checks the same thing at `assign_hap.c:326`.

So the guards stay, and the reason is on the record: they encode a difference in
the data model, not a difference in judgement. Exact textual parity there would
mean crashing where upstream reads adjacent memory.

One consequence for the test suite: the parity case asserting that a read with
no profile returns -1 pins OUR guard rather than upstream behaviour. It is kept,
because the guard is kept, but it should be read as pinning a deliberate
divergence.

## Not pursued, with reasons

- Dropped `LONGCALLD_VERBOSE` debug blocks in `wfa_end2end_aln`,
  `wfa_collect_aln_str` and `edlibAlignmentToXGAPS`: cannot affect output.
- `wfa_end2end_aln`'s CIGAR collection: we never consume a CIGAR from it, and
  its return value carries the alignment length instead of a constant 0.
- `main`'s dispatch: upstream dispatches `call`, we dispatch three `collect`
  subcommands. Not a port question.
- `get_var_site_start` starting its binary search at 0 rather than at the
  caller's `cur_site_i`: same result, more work.
- `update_read_var_profile_with_allele`'s backward extension and
  `graph_alleles` maintenance, and `push_digar_alt_seq` copying `alt` where
  upstream nulls it: both have upstream counterparts (`bam_utils.c:251`,
  `align.c`) and both are still open. Not attempted here.
