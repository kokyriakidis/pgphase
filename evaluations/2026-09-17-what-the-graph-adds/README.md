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
