# Graph-first, as a mode rather than an invocation

Until now "graph-first" meant a hand-made pair of files: a header-only
`--private-sites` VCF to drop the alignment channel's candidates and a per-gap
`--bam-authoritative-bed` to let them back in. That is not a mode, and it was
not testable as an arm.

`--graph-first` now does it directly: the alignment channel's own candidates are
withheld from the first pass, the catalog's sites are injected as usual, and
`recover_windows_with_targeted_solve` puts alignment evidence back inside every
window the first pass could not phase -- a full alignment solve over that
region at the recovery mapq floor, stitched by shared reads.

## Two false starts, both instructive

**Gating on ownership after the claim pass does not work.** The obvious place to
withhold is after `inject_graph_sites`, keeping what the catalog owns. Measured
over `chr20:25,979,591-26,138,679`, injection claims **3,741 of 3,743**
candidates -- the same count whether tested by the `graph_site` flag or by the
`all_graph_cands` index set it returns. An ownership filter there withholds two
sites and the mode is inert.

**And erasing entries there would have been a bug.** `inject_graph_sites`
returns index sets into `chunk.candidates`, and every later stage -- including
`backfill_graph_candidate_counts` -- addresses candidates by those indices.
Erasing renumbers them silently. The first implementation did exactly that; the
second zeroed categories instead, which is index-safe, before measurement showed
the whole approach was inert anyway.

So the withholding happens **before** injection, which is also where the
`--private-sites` route acts.

## What each flag actually does, separated

| arm | blocks | in-gap | reads tagged | read concordance | discordant |
|---|---:|---:|---:|---:|---:|
| default (union) | 2 | 278 | 293 | 99.66% | 1 |
| `--graph-first` alone | 2 | 278 | 293 | 99.66% | 1 |
| `--graph-authoritative` alone | 2 | 278 | **254** | 99.21% | 2 |
| both | **1** | 278 | 254 | 99.21% | 2 |
| the old BED route | 1 | 278 | 254 | 99.21% | 2 |

**This corrects the attribution recorded earlier.** The BED arm's behaviour was
credited to phasing on graph sites only. It is not: restricting the site set is
**inert on its own**. What `--graph-first` withholds is small -- measured,
`[graph-first] withheld 102 alignment-discovered candidate(s) from the first
pass` over this region, against a final candidate table of 451 rows -- and the
catalog's own sites carry the first pass either way. The read-tagging difference
comes from `--graph-authoritative`, which replaces the alignment channel's read
evidence with the graph's at every claimed site. The block merge needs both.

(The 3,741-of-3,743 figure above is a different quantity at a different stage:
candidates claimed by injection at the point the ownership filter was tried,
before pruning. It is not the catalog's share of the emitted table, and an
earlier version of this document mixed the two.)

Which means the hybrid's first pass is already graph-driven in the site sense.
What "graph-first" adds beyond that is evidence ownership.

## The arm is now under test, and it fails honestly

`graphfirst` = `--graph-first --graph-authoritative` is a window-test arm, so
regressions are caught while it is iterated on. Its first run failed the
concordance gate:

```
graphfirst / window 5309406
  CHECK( got.concordance() >= 0.95 )  ->  0.9403794038 >= 0.95
  spans=yes in_gap_hets=3 blocks=3 tagged=369 discordant=22
```

That is the read-placement cost, measured: a third of the clean het anchors
outside a gap are alignment-only, including every clean het indel in the regions
examined, so fewer reads find a site to sit on.

The floor is now per arm, emitted by the test binary and rounded **down** --
`default` 0.98, `noretry` 0.99, `graphfirst` 0.94. Recording graph-first's
measured floor states the deficit in the baseline, where a loose global bound
would hide it and a fixed 0.95 would block the suite. Any further drop fails.

It remains a rate, not a read count: how many reads get tagged still has no
expectation, for the reason recorded earlier.

## Next

Read placement is the open problem, and it is now the only one separating
graph-first from the default: both report 278 in-gap records at the same site
accuracy. Recovering placement outside gaps without re-admitting the full
alignment site set is the work.

Not measured: chromosome-wide. Those runs are on hold by instruction.
