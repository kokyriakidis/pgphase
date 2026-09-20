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

## The two open items, closed

### `push_digar_alt_seq` now nulls `alt` where upstream does

`bam_utils.c:573-577` sets `alt_seq = NULL` on the new-element path and copies
it only for `BAM_CDIFF` and `BAM_CINS`. We pushed the whole op, so a deletion,
match or clip kept whatever alt string it arrived with. Now cleared for every
type except `Snp` and `Insertion`.

The neighbouring `push_digar0` was checked at the same time and needs no change:
upstream assigns `alt_seq = d.alt_seq` there unconditionally, with no type test
(`bam_utils.c:600`), which is what we already do. The nulling is specific to the
`_alt_seq` variant.

Whole chr20: both arms **byte-identical**. So the field was never read for those
types on this data -- the divergence was latent, and is now closed rather than
left to surface on data that does read it.

### `update_read_var_profile_with_allele`: one branch was dead, the other is caller-driven

Upstream writes into a PRE-SIZED profile at `var_i - start_var_idx` and sets
`end_var_idx = var_i` unconditionally (`bam_utils.c:250-254`). Ours grows the
vectors, extends `end_var_idx` only upward, and had a branch that grew the
profile leftward when a call arrived below `start_var_idx`. Instrumenting every
call over `chr20:5-9 Mb`:

| arm | calls | below `end_var_idx` | below `start_var_idx` |
|---|---:|---:|---:|
| alignment | 64,000,000 | **0** | **0** |
| graph | 2,000,000 | 1,490 | **0** |

Two conclusions, each acted on differently.

The leftward branch was **dead** -- no call on either arm arrives below
`start_var_idx` -- so it is removed, and whole chr20 is byte-identical in both
arms.

The `end_var_idx` difference is **not observable on the ported path**: zero of
64 million alignment-arm calls arrive below the current end, so upstream's
unconditional assignment and our upward-only extension cannot be told apart
there. They differ only for the graph arm's merge paths, which renumber
candidates and do call below the end 1,490 times -- a caller upstream does not
have, and one for which upstream's version would shrink the span and truncate
the profile. So the conditional stays, now with the number that says why.

## Is the residual worth closing? Scored against GIAB truth

Every record in both residual sets was classified against
`HG002_CHM13v2.0_v5.0q_smvar` inside the benchmark BED:

| | upstream emits, we do not (370) | we emit, upstream does not (373) |
|---|---:|---:|
| outside the benchmark region | 299 (81%) | 344 (92%) |
| **TRUE variant, exact allele** | **12** | **13** |
| **FALSE positive, no truth variant** | **25** | **9** |
| truth has the position, other allele | 34 | 7 |

Matching upstream exactly on these records would **gain 12 true variants and
lose 13, while gaining 25 false positives and shedding 9** -- net one fewer true
variant and sixteen more false ones.

For the 254 held-but-unemitted class specifically: 25 false positives against 12
true variants inside the benchmark. Our abstention there is right about twice as
often as it is wrong, which is the same thing the labelling comparison said from
the other direction -- upstream labels 45% of spanning reads at 60.2%
truth-consistency where we label 20% at 65.6%. It emits more because it commits
more, and most of what it commits at these loci is not in the truth set.

So the residual is not an accuracy target. It is only worth closing if the goal
is bit-parity with the reference implementation, which is a different goal and
should be chosen deliberately.

Caveat on the size of the signal: 81-92% of both sets fall outside the benchmark
region, so the verdict rests on 37 and 22 assessable records respectively. The
direction is consistent with the independent labelling measurement, but these are
small numbers and should not be quoted as precise rates.

## The real cost is already known, and it is reversible

The nine divergences removed today were suppressing false positives. Their
removal took GIAB F1 from 0.9677 to 0.9632 and misplaced reads from 1,516 to
3,226. Every one is behind a named option, defaulted to this project's behaviour
and set faithful only in `collect_bam_variation`, so the choice between "matches
upstream" and "scores better" is now one flag per divergence rather than a fork.
