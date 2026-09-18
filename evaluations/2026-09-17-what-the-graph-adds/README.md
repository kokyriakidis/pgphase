# What do the graph sites add?

Asked directly: why do we need them at all? Measured by running the hybrid with
an **empty catalog** -- a header-only sites VCF, bgzipped and indexed, so the
only thing removed is the graph's sites while every other code path, including
the targeted recovery, is untouched. The loader confirms it:
`[graph-sites] ...: 0 data lines, 0 sites parsed, 0 eligible`.

## The architecture, first, because the question rests on it

The hybrid is the **alignment pipeline with the graph injected into it**, not the
graph with the BAM injected. `process_chunk_hybrid` runs
`collect_var_classify` (alignment discovery and classification) first, then
Phase A `inject_graph_sites`, then profiles, then Phase B
`inject_graph_reads`. Every injection-side function is named `graph`; there is
no BAM-injection path, because the alignment candidates are the table the graph
is merged into.

So "phase with the graph first, then add the BAM" does not describe the code.

## chr20:26,029,591-26,088,679 -- the catalog is actively limiting

| arm | hets | blocks | largest | spans | in-gap records | reads tagged | read concordance | discordant |
|---|---:|---:|---:|---|---:|---:|---:|---:|
| alignment channel alone | 151 | 3 | 43.5 kb | no | 0 | 293 | 99.66% | 1 |
| hybrid, real catalog | 198 | 2 | 151.1 kb | YES | 35 | 293 | 99.66% | 1 |
| **hybrid, empty catalog** | **430** | 2 | 150.6 kb | YES | **276** | 293 | 99.66% | 1 |

Site-level quality of the in-gap calls, scored against read truth:

| arm | SNPs scored | clean | wrong | error |
|---|---:|---:|---:|---:|
| real catalog | 17 | 15 | 2 | 12% |
| **empty catalog** | **246** | **243** | 3 | **1%** |
| longphase | 120 | 119 | 1 | 1% |

With no catalog the pipeline reports **twice the competitor's scorable SNP count
at the same 1% error**, spans the gap, and tags reads exactly as accurately.

The cause is not that catalog sites are bad calls. It is that the catalog's own
phased sites near the gap **narrow the unphased window** `collect_unphased_windows`
reports, and the targeted solve's import is confined to that window: 38 sites
imported with the catalog, 282 without. The 37-against-120 shortfall recorded
earlier is therefore a scoping artifact of the import, not a limit of the
recovery.

## chr20:5,309,406-5,345,085 -- equivalent, and a correction

| arm | blocks | spans | in-gap | tagged | concordance | discordant |
|---|---:|---|---:|---:|---:|---:|
| real catalog | 1 | YES | 2 | 544 | 98.71% | 7 |
| empty catalog | 1 | YES | 2 | 545 | 98.53% | 8 |

**This corrects an earlier attribution.** That window's closure was credited to
the catalog claim at 5,315,591 surviving the repeat screen. That was true when
it was measured, before the targeted solve existed. With the targeted solve in
place the window closes with no catalog at all, one extra tagged read and one
extra discordant one.

## What this does and does not establish

It does not establish that the graph sites are unnecessary. Both windows are
**gaps**, which is to say places the graph-only pass already failed -- the
selection is biased to exactly where the graph is weakest, and the gap set
itself was derived from a graph-only baseline.

What it establishes is narrower and still worth having:

- on these two windows the catalog adds nothing to phasing, and on the larger
  one it costs 241 in-gap records;
- the import scoping, not the recovery, is what limits in-gap recovery;
- the graph channel alone phases 103 hets here in two blocks that stop exactly
  on the gap bounds, so it does not reach into the hole either.

The test that would answer the general question is a chromosome-wide
empty-catalog arm against the standing 212,320 tagged / 452 blocks / 0.559%
baseline. It has not been run -- chromosome-wide runs are on hold by
instruction.

## Can we use the graph sites and inject BAM sites only in the gap?

It needs no new architecture -- the mode already exists. `--private-sites FILE`
restricts the table to graph sites plus a whitelist, and
`--bam-authoritative-bed FILE` lets clean BAM sites back in inside given
intervals. Passing a header-only whitelist and a BED holding just the gap is
exactly "graph sites everywhere, BAM sites in the gap".

Measured on both windows, with the targeted recovery active in every arm.

### chr20:26,029,591-26,088,679 -- graph-first wins on structure and sites

| arm | blocks | spans | hets | in-gap | in-gap SNPs | wrong | tagged | concordance | discordant |
|---|---:|---|---:|---:|---:|---:|---:|---:|---:|
| **graph + BAM in gap only** | **1** | YES | 338 | **276** | 246 | 3 | 254 | 99.21% | 2 |
| current default (union) | 2 | YES | 198 | 35 | 17 | 2 | 293 | 99.66% | 1 |
| empty catalog (BAM only) | 2 | YES | 430 | 276 | 246 | 3 | 293 | 99.66% | 1 |

One block instead of two, and the full in-gap set at the competitor's 1% error
rate -- but 39 fewer reads tagged and one more discordant.

### chr20:5,309,406-5,345,085 -- graph-first loses clearly

| arm | blocks | spans | hets | in-gap | tagged | concordance | discordant |
|---|---:|---|---:|---:|---:|---:|---:|
| graph + BAM in gap only | 3 | YES | 21 | 2 | 437 | **96.80%** | **14** |
| current default (union) | **1** | YES | 37 | 2 | **544** | **98.71%** | 7 |

Three blocks instead of one, 107 fewer reads tagged, and **twice the discordant
reads**. The cause is the provenance measurement above: outside the gap a third
of the clean het anchors are alignment-only, including every clean het indel in
these regions, so excluding them costs read placement -- which is what read
concordance measures.

### Conclusion

The architecture is implementable today and is not a clear win: it buys block
structure and in-gap site recovery on one window and costs read-level accuracy
on both, badly on one.

The better route is to take what graph-first bought without its cost. The
in-gap recovery it achieved -- 276 records, 246 scorable SNPs at 1% error -- is
the same figure the empty-catalog arm reached **while keeping the union first
pass and its 99.66% / 1 discordant**. So the gain is not coming from removing
BAM sites in the flanks; it is coming from the recovery importing the whole gap
rather than the narrow detected window. Widening that import scope gets the
benefit with the flanks untouched, and it is a change to one predicate rather
than to the pipeline's shape.

## Why graph-first fragmented, and a correction to the verdict above

The verdict "graph-first loses clearly" on chr20:5,309,406 needs qualifying.

Graph-first produced three blocks with 34.6 kb and 12.0 kb seams. The recovery
never touched them, and the reason is a gap in the trigger, not in the solve:
`collect_unphased_windows` reports where the solve left READS unphased. At a
seam the reads are phased -- into different blocks -- so no window is reported
and the targeted solve never sees it. The alignment-driven pipeline mostly
produces the first failure mode; a graph-first pipeline produces the second.

`collect_block_seams` now adds every interval between consecutive blocks as a
recovery window, taken from candidate positions rather than read starts because
a block's extent is first to last phased site. With it the fragmented window
runs seven targeted solves instead of two, and the new per-window report shows
what each one found.

**And the seam turns out to be unbridgeable.** For the 34.6 kb seam:

```
[targeted] 5315085-5421703: 2 vote pair(s), 2 parent block(s) linked,
                            2 targeted block(s) carrying links
```

Each parent block linked to a *different* targeted block -- the targeted solve
breaks at the same place. Directly measured: **0 reads span 5,345,085-5,379,675**.
No read-based method can join it, so the refusal is correct.

The union arm *does* join across it, in one block of 37 sites, and by truth the
join is right: left PAT-on-1 at 0.981 over 323 reads, right PAT-on-1 at 1.000
over 54. But with no read crossing the interval that orientation is a coin flip
that landed correctly -- the same pattern recorded for chr20:48,176,830 earlier
in this work, where an identical-looking join was wrong and no read-level gate
could see it.

So part of graph-first's higher block count is **honesty, not weakness**: it
refuses a join the union arm makes on no evidence. What remains a genuine cost
is read placement -- 437 reads tagged against 544, and 14 discordant against 7 --
which follows from a third of the clean het anchors outside the gap being
alignment-only. That is a separate axis from the seam question and is where the
work would go next for a graph-first hybrid.

Seam recovery is inert on the current default on both windows -- identical
blocks, sites, tags and concordance -- because the union pipeline's seams are
either already joined or, as here, unbridgeable.
