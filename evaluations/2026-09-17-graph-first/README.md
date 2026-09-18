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

## Graph-first is now the default

`graph_first = true`, with `--no-graph-first` to get the union first pass back.
The architecture is single: the catalog drives the first pass, and the alignment
channel is a recovery mechanism that fixes what the first pass could not phase,
per window, at the recovery mapq floor.

The committed baseline, `src/test_gap_windows_expect.tsv`. Its columns are
`spans` (a flag, not a block count), `min_in_gap_hets` and `min_concordance`:

| arm | 5,309,406 | 26,029,591 |
|---|---|---|
| **default (graph-first)** | spans, 2 in-gap, floor 0.98 | spans, 278 in-gap, floor 0.99 |
| `nographfirst` (union) | spans, 2 in-gap, floor 0.98 | spans, 278 in-gap, floor 0.99 |
| `noretry` | no span, 1 in-gap, floor 0.99 | no span, 0 in-gap, floor 0.99 |
| `graphauth` | spans, 3 in-gap, floor **0.94** | spans, 278 in-gap, floor 0.99 |

An earlier version of this table printed the `spans` flag as a block count, so
every row read "1 blk". Block counts are not in the expectations file; where one
is quoted below it is measured directly.

**The flip costs nothing on the panel and is not a behaviour change there** --
the default and the union arm agree on every recorded quantity, because the
catalog's sites carry the first pass either way and the targeted solve supplies
the rest. What changed is which channel is primary, and therefore what happens
in a region the catalog covers but the alignment channel would also have called:
the catalog's call stands.

**Evidence ownership stays opt-in and separate.** `--graph-authoritative`
replaces the alignment channel's read evidence with the graph's at every claimed
site. It is the axis that costs accuracy. Measured directly on
chr20:5,309,406, default against `--graph-authoritative`:

| arm | blocks | spans | in-gap hets | reads tagged | read concordance | discordant |
|---|---:|---|---:|---:|---:|---:|
| default | **1** | yes | 2 | 544 | 98.71% | 7 |
| `graphauth` | **3** | yes | 3 | 369 | **94.04%** | 22 |

It spans in both arms, so `spans` alone does not separate them -- the block
count, the reads tagged and the concordance do. Hence an arm with its own
recorded floor rather than part of the default.

Arms are now `default` (graph-first), `noretry`, `nographfirst` and `graphauth`;
132 assertions.

Not measured: chromosome-wide. Those runs are on hold by instruction, and this
changes the default path, so it is the number to take before calling the flip
settled.

## Does the graph channel phase what graphauth hands it?

`--graph-authoritative` discards the alignment channel's read evidence at every
claimed site and keeps the graph's. Whether that is a trade or a loss depends on
whether the graph's own evidence phases those sites. Measured on
chr20:5,309,406-5,345,085, all three over the same region and the same reads
(639 GAF records, 639 BAM reads):

| arm | phased hets | blocks | spans | in-gap hets |
|---|---:|---:|---|---:|
| graph channel alone (`collect-graph-variation`) | **23** | 3 | **no** | **0** |
| hybrid + `--graph-authoritative` | **23** | 3 | yes | 3 |
| hybrid default | **36** | 1 | yes | 2 |

The graph channel alone does **not** phase this gap. Its largest block is
5,272,413-5,309,406, stopping exactly on the gap's left bound, and it phases
nothing inside.

And the hybrid under `--graph-authoritative` lands on the same 23 phased hets as
the graph channel alone. That is the flag's effect stated as a number: at claimed
sites it reduces the hybrid to the graph channel's own phasing power. Its span
comes entirely from the targeted alignment recovery, not from the graph
evidence it substituted in.

So the evidence swap is a loss here, not a trade -- it discards 13 phased
heterozygotes' worth of alignment evidence and gains nothing the graph could
phase on its own. That is the measurement behind keeping it out of the default.

## --graph-authoritative is removed

The mechanism is gone, not merely defaulted off: the `Options` field, the CLI
flag, `clear_bam_evidence_at_graph_candidates` and its declaration, the
`graph_owned_cands` switch, and the unit-test assertions that covered the
clearing. `graph_owned_cands` is now always `graph_only_cands` -- the candidates
the graph contributed and the alignment channel does not have -- and
`backfill_graph_candidate_counts` is no longer conditional.

A whitelist (`--private-sites`) used to enable the mode implicitly unless
`--private-msa-admit-all-in-region` was also given. That is what made the
whitelist route look like a site-selection change when it was an evidence
change, and it is the mis-attribution corrected earlier in this document. A
whitelist now scopes retention only.

Default unchanged: `spans`, 2 in-gap hets, floor 0.98 on chr20:5,309,406 and
`spans`, 278 in-gap, floor 0.99 on chr20:26,029,591. Arms are `default`,
`noretry`, `nographfirst`; 99 assertions, unit 4/4.

## Is the hybrid's first pass already the graph channel's result?

That was the goal stated for the architecture, so it is worth measuring rather
than assuming. On chr20:5,309,406-5,345,085:

| arm | phased hets | blocks | spans |
|---|---:|---:|---|
| graph channel alone (`collect-graph-variation`) | 23 | 3 | no |
| hybrid, first pass only (`--no-retry-unphased-with-bam`) | **26** | 3 | no |
| hybrid, default (first pass + recovery) | **36** | **1** | **yes** |

As position sets, the first pass against the graph channel: **22 shared, 4 the
hybrid phases and the graph does not, 1 the graph phases and the hybrid does
not.** All 4 of the extra sites are in the graph catalog -- so the first pass is
working the graph's site set, and simply phases more of it, because the
alignment reads carry the evidence at those sites rather than the GAF records.

That is the shape asked for, with one qualification worth stating: the first
pass is not *identical* to the graph channel and should not be. Forcing identity
means substituting the graph's read evidence at every claimed site, which is the
mode just removed -- it costs 13 phased heterozygotes and 175 reads here.
"Same sites, better evidence" is the achievable version.

**The one divergence is a disagreement, not a miss.** At 5,293,282 the graph
channel phases `C>A`; the hybrid holds a different event at that locus, an
insertion `C>ACACAA` at DP 80 and 40/40, classed noisy and so unphased in the
first pass. Same position, two representations.

**Then recovery does exactly what it should:** 26 phased hets to 36, three
blocks to one, and the gap spans. The added sites are alignment-derived and
arrive only inside the windows the first pass could not phase.
