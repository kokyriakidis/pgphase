# A catalog claim survives the repeat screen

## Why claim-driven promotions never reached a gap

Over the six panel windows the hybrid holds 37 sites that are `CLEAN_HET` only
because the graph claimed them -- and every one is in a **flank**, none inside a
gap. The reason is a demotion, traced at `chr20:5,315,591` with an env-gated
probe on the candidate's category at three points in one run:

| point | category | counts |
|---|---|---|
| after `classify_graph_only_candidates` | `0x008` `CLEAN_HET_INDEL`, `graph_site = 1` | DP 68, 37/31 |
| before `collect_noisy_vars_step4` | `0x010` `REP_HET_INDEL` | DP 68, 37/31 |
| after it | `0x100` `NOISY_CAND_HET` | DP 68, 37/31 |

The counts never move; only the class does. `apply_hybrid_noise_filter`
(`hybrid_inject.cpp:887`) asks `is_noisy_site` -- does this locus sit in a repeat
tract -- and on a yes overwrites the category **and** `candvarcate_initial`. That
question is true of most real indels in these windows and says nothing about
whether this one segregates: read truth puts `5,315,591` at **1.000 over 69
reads**, and the competitor bridges this gap with it.

Once the site is `REP_HET_INDEL` it is screened out of k-means by construction,
and `merge_var_profile`'s `replace_repeat` rule then lets the MSA's own
`NOISY_CAND_HET` call displace it -- a class the hybrid excludes from its solve.
So a gap, being repeat-dense, undoes exactly the promotions it needs.

## The change

The demotion is skipped when the candidate is `graph_site` **and** already
classified `CLEAN_HET_SNP` or `CLEAN_HET_INDEL` -- that is, when the catalog
claims the locus and the alignment channel's own depths call it a clean het. It
is not a blanket admission of the repeat class: measured earlier against read
truth, `REP_HET_INDEL` as a class is 23% informative and 62% phantom, and
admitting it wholesale injects about three phantoms per real site.

## Panel, stock defaults

| | before | after |
|---|---:|---:|
| in-gap phased hets | 12 | **14** |
| reads tagged | 2,794 | **2,831** |
| read concordance | 99.68% | **99.72%** |
| concordant -> discordant | | **0** |

39 newly concordant reads, 0 newly discordant, and one concordant tag lost: a
9.7 kb read at `48,169,957-48,179,694` (truth PATERNAL) in window 48,183,976,
which was concordant in block `48,162,480` and is now untagged. That window
gains nothing, so its `min_tagged` floor drops 393 -> 392; the retry floors for
window 5,309,406 and the retry panel total drop by the same one read.

## The target gap: what closes it, and what does not

`chr20:5,309,406-5,345,085`. Two sites inside it are informative against truth:
`5,315,591` (`C>CT`, 1.000) which the catalog claims, and `5,331,265`
(`TAGAC>T`, 1.000) which it does **not** -- so the claim lever reaches the first
and cannot reach the second.

On stock defaults the fix extends the left block ~6 kb into the gap on a
42-to-0 link, and the gap stays open:

| site | class | phase set | link |
|---|---|---|---|
| 5,309,406 `C>T` | `CLEAN_HET_SNP` | 5,272,413 | 11 / 0 |
| **5,315,591 `C>CT`** | **`CLEAN_HET_INDEL`** | **5,272,413** | **42 / 0** |
| 5,331,266 `TAGAC>T` | `NOISY_CAND_HET` | 0 | -- |
| 5,345,085 `G>A` | `CLEAN_HET_SNP` | 5,345,085 | 0 / 0 |

With `--retry-unphased-with-bam` the window-scoped admission reaches
`5,331,266` as well and the gap **closes**: one block `5,263,741-5,394,961`
over 37 sites, chained `5,315,592` (42/0) -> `5,331,266` (6/0) ->
`5,339,369` (31/0).

The orientation is correct, which is the part a read-level accuracy number
cannot show on its own: scoring the merged block's reads separately at each end
puts the same parent on haplotype 1 on both sides (left 0.968, right 0.994), so
there is no switch across the gap.

| arm | blocks | tagged | concordance | across the gap |
|---|---:|---:|---:|---|
| claim fix only | 3 | 466 | 99.79% | one side each |
| claim fix + retry | **1** | **544** | 98.71% | **consistent** |

The cost is visible and worth stating: 78 more reads carry a tag and the
discordant count goes from about 1 to about 7. None of those are reads that were
correct and became wrong -- the panel gate reports 0 concordant-to-discordant,
348 newly concordant and 7 newly discordant -- so the trade is newly tagged
reads, a minority of which are wrong, against coverage.

Not measured: chromosome-wide effect. Per instruction, verification for this
change is the window and the panel only.

## The window tests no longer assert read counts

A `min_tagged` floor fails on changes that are improvements. This fix loses one
tagged read in window 48,183,976 -- a 9.7 kb read at 48,169,957-48,179,694,
truth PATERNAL, concordant in block 48,162,480, now untagged -- while adding 39
concordant reads elsewhere, and that tripped three floors by exactly one read.
Lowering them would have laundered the loss; keeping them would have blocked the
gain. Neither is the property worth asserting.

The expectations file is now two columns, `spans` and `min_in_gap_hets`, and the
assertions that matter do not live in a file at all and cannot drift:

| assertion | why |
|---|---|
| `spans` equality per window | a span where none is expected is a join across an interval no read crosses |
| in-gap phased hets, floor | a count of SITES phased inside the gap |
| **no clean het inside a gap without a phase set** | a site we hold and are allowed to use, left unused |
| **no block switches across a gap** | the hazard read-level accuracy cannot see: with no read crossing, an inverted join still scores 100% and only the two ends disagree |
| concordance >= 0.95 | a loose floor that fires on a collapse, not on a decision that moves a handful of reads |

The switch check places each end of a block independently -- at least five
scored reads and at least 90% agreement -- and fails when two confidently placed
ends disagree.

## The closure is now guarded by the sites it rests on

`spans` alone can pass by finding some other way across. So every heterozygote a
closing arm phases **strictly inside** the gap is recorded in
`src/test_gap_windows_required.tsv` -- 18 rows over the four windows that close
-- and each is asserted twice in that arm:

| assertion | what its failure means |
|---|---|
| RETRIEVED: a candidate exists at the position (+/-2 bp) | a discovery or injection regression -- we stopped calling the site |
| USED: that candidate carries a non-zero `PHASE_SET` | an admission regression -- we call it and exclude it from the solve |

For the target window the two rows are `5,315,591 C>CT` (`CLEAN_HET_INDEL`, the
site the claim fix recovers) and `5,331,265 TAGAC>T` (`NOISY_CAND_HET`, the one
only the window-scoped retry reaches). Either going missing or going unused
fails the suite even if the gap still closes by some other route.

The file is written by the test binary during a refresh, from the same run that
writes the expectations, so it cannot drift from a second implementation of
"which sites close this gap". The +/-2 bp match is needed because an insertion
anchors its candidate one base past the position the VCF reports.

Validated by injecting both faults: an invented position reports
`NOT RETRIEVED`, and requiring `5,331,265` in the **default** arm -- where it is
retrieved but unphased -- reports `RETRIEVED BUT UNUSED`. Both fail the suite;
the real file passes at 178 assertions.
