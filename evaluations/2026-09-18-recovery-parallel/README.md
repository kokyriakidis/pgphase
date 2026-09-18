# The recovery step ran on one core

`recover_windows_with_targeted_solve` looped over its windows serially, and each
sub-solve set `sub.threads = 1`. So the whole step used a single core while the
parent pipeline had `opts.threads` workers -- the parent's own chunk loop opens a
`WorkerContext` per worker (`collect_pipeline.cpp:540`), and the recovery did
not.

Measured on `CHM13#0#chr20:1-10000000` at `-t 16`:

| arm | wall |
|---|---:|
| hybrid, `--no-retry-unphased-with-bam` | 31.9 s |
| hybrid with recovery, serial loop | 246.7 s |
| hybrid with recovery, parallel loop | **88.7 s** |

So the recovery portion went 214.8 s -> 56.8 s, **3.8x**, and the phased VCF is
**byte-identical** (19,618 records both ways); the `[targeted]` log lines are
identical as a set, 206 windows either way.

## What changed

The loop is now two phases. Phase 1 solves every window in parallel, each
sub-solve still single-threaded, each worker constructing its own
`WorkerContext` because htslib file handles are not shareable. Phase 2 applies
the results in the original window order, so the parent mutations happen in the
same sequence as before and the outcome does not depend on completion order.

## What is still wasteful, and needs a decision rather than a patch

Over the 196 of 206 logged windows whose region parses cleanly:

| | |
|---|---:|
| total bp re-solved | 15.06 Mb (1.5x the 10 Mb region itself) |
| union of those regions | 6.83 Mb in 38 disjoint groups |
| overlap redundancy | **2.20x** |
| actual gap/seam span | 3.30 Mb |
| flank overhead (`kTargetedSolveFlank` = 30 kb each side) | **4.6x** |

Merging the 196 regions into their 38 disjoint groups would cut the work by
2.2x, and shrinking the flank would cut it further -- but neither is
behaviour-preserving. A merged region is solved as one chunk, so its k-means
sees a different site set and can reach a different orientation; the flank is
what gives the stitch enough overlapping reads to vote. Both therefore need
measuring against the arm baseline, unlike this change.
