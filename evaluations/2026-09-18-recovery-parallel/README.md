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

## Round two: one solve per merged region

Each window is padded by `kTargetedSolveFlank` = 30 kb per side, and that
padding is what makes neighbouring windows overlap. Measured on the graph +
recovery arm over `CHM13#0#chr20:1-10000000`:

| | |
|---|---:|
| windows | 188 |
| raw windows (padding removed), median span | 18.0 kb |
| raw windows -> disjoint groups | 188 -> 134, so **54 overlap before any padding** |
| exact duplicate raw windows | 0 |
| disjoint raw pairs closer than 60 kb | 98 (overlap created purely by the flank) |
| bp of raw windows alone | 4.25 Mb |
| bp actually solved, padded | 15.53 Mb |
| after merging | 6.24 Mb, 41 regions |

The 54 pre-padding overlaps come from the window list being the concatenation of
`collect_unphased_windows` and `collect_block_seams`, which describe some of the
same neighbourhoods. Windows whose padded regions touch are one piece of work,
so they are now merged and solved once, with the import-span fallback covering
every member window of the group.

| 10 Mb, graph + recovery | wall | bridged | tagged | blocks | discordant | read hamming |
|---|---:|---:|---:|---:|---:|---:|
| per-window, serial recovery | 155.4 s | 38 | 36,645 | 43 | 97 | 0.265% |
| grouped + parallel | **36.8 s** | 36 | 36,645 | 44 | 97 | 0.265% |

**4.2x end to end with read accuracy unchanged** -- same tagged count, same 97
discordant reads, same 0.265%. It is not free: two bridges are lost and the
block count goes 43 -> 44, so the cost is contiguity rather than misplacement.
A merged region is solved as one chunk, so its k-means sees every site in the
group instead of one window's worth.

## Round three: one batch for the whole chromosome

Recovery was called once per chunk, from a serial loop over the chunks
(`graph_collect.cpp:865`), and each call parallelised only its own windows. Over
the first 10 Mb of chr20 that was **18 sequential invocations** with 1-4 merged
regions each -- a mean parallel width of **2.3 of 16 threads**.

The merged regions are disjoint, so nothing orders them. `prewarm_targeted_solves`
now collects every chunk's windows, merges them with the same
`build_targeted_groups` the per-chunk path uses, deduplicates region keys (two
chunks can see the same seam at their shared boundary), and solves them all in
one batch across the full thread pool, largest region first so the biggest does
not start last. The per-chunk entry point then finds its regions already solved
and applies them in region order exactly as before.

Two things keep this a pure speedup rather than a behaviour change: the regions
and their application order are unchanged, and the cache holds a
`TargetedSolveResult` -- read names, haplotypes, phase sets, candidates -- rather
than whole `PhasingChunk`s, so a chromosome's worth of solves fits in memory
without keeping sequences and alignments nothing downstream reads.

| graph + recovery | 10 Mb | whole chr20 |
|---|---:|---:|
| per-window, serial recovery | 155.4 s | -- |
| per-window, parallel | -- | 443 s |
| grouped, parallel per chunk | 36.8 s | -- |
| grouped, one global batch | **16.8 s** | **137 s** |

The 10 Mb output is **byte-identical** to the grouped per-chunk run, and whole
chr20 emits the same 56,032 records. Bridged blocks chromosome-wide go 132 ->
126, which is the merged-region cost recorded above, not the batching.

Suites: unit 4/4, window 66/66, predicate 130/130.

## Round four: the flank measured in sites, not base pairs

Why a flank at all, given the main pipeline stitches chunks without one: a chunk
boundary is straddled by reads, so `initialize_chunk_overlap_state` already has
the same read in both chunks and the stitch can vote on it. A window has no
equivalent -- the reads inside it are exactly the ones the parent left unphased,
so they carry no parent haplotype. `select_stitch_orientation` needs reads
carrying BOTH a parent haplotype and a sub-solve haplotype, so the region has to
reach out to where the parent did phase.

How far, measured over the first 10 Mb of chr20 (9,813 parent phased sites, 41
merged regions) -- the extension needed to reach parent sites on **both** sides:

| sites reached | median | p90 | max |
|---|---:|---:|---:|
| 1 | 0.8 kb | 7.3 kb | 46.2 kb |
| 2 | 3.4 kb | 16.7 kb | 46.8 kb |
| 3 | 5.6 kb | 16.9 kb | 47.3 kb |

A fixed 30 kb is therefore wrong in both directions: 5x more than needed in the
median case, and not enough for one of the 41 regions, which had no parent site
to vote against at all.

The flank now reaches `kTargetedSolveFlankSites` = 3 parent phased sites per
side, clamped to [2 kb, 60 kb], with the clamp maximum used when there are not
three sites that way -- nothing to anchor against nearby, so reach as far as
allowed rather than solve a region that cannot vote.

| chr20 | wall | regions | bp solved | bridged | tagged | blocks | discordant | hamming |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| fixed 30 kb | 137 s | 236 | 39.55 Mb | 126 | 203,751 | 283 | 2,851 | 1.399% |
| site-based | **120 s** | 356 | **29.78 Mb** | **127** | 203,751 | 281 | 2,851 | 1.399% |

25% less sequence solved, one more block bridged -- the region that could not
reach a parent site before -- two fewer blocks, and read accuracy identical to
the read: same 203,751 tagged, same 2,851 discordant, same 1.399%. The same 10 Mb
comparison gives 6.38 -> 3.98 Mb solved at 36 bridged and 0.265% both ways.

Region count rises (236 -> 356) because smaller extensions overlap less, so
fewer windows merge -- more solves, each much smaller.

## A bug the invariant check found: a merge flip never reached the VCF

Asked whether the recovery leaves already-computed phase sets alone, the answer
from the code is yes -- adoption is guarded on `existing->phase_set != 0` ("the
parent's own call stands"), import only adds, and a block merge rewrites a whole
block uniformly. Checking it against output found something else.

Comparing the arm against itself with the recovery loop skipped (an env-gated
probe, same binary and inputs, so only the recovery differs) over the first
10 Mb of chr20: all **9,801** parent-phased records came out with an **identical**
genotype, 0 flipped. But a probe on the merge itself reports **6 of 18 merges
carry `flip = 1`**, covering 1,987 candidates.

The cause: the merge path swapped `hap_alt`/`hap_ref`, and the emitter does not
read that pair. `is_alt_genotype` calls `derive_hap_alt_ref_from_consensus`
(`collect_output.cpp:383`, `:574`), which reads `hap_to_cons_alle[1]` and `[2]`.
So a flipped merge relabelled the absorbed block into the kept phase set while
leaving its emitted genotypes in the original orientation -- the two halves came
out anti-phased -- and the phased BAM disagreed with the VCF, because
`chunk.haps` above does flip. The adoption path 90 lines below already handled
this, with a comment recording that setting the pair alone left 244 adopted sites
dropped by the VCF.

Parity of each `flip = 1` merge against read truth, hap1's parent per half:

| merge | before | after |
|---|---|---|
| 130540 <- 195305 | ANTI-PHASED | CONSISTENT |
| 195305 <- 545002 | ANTI-PHASED | (separate blocks, see below) |
| 1690618 <- 1907095 | ANTI-PHASED | CONSISTENT |
| 2342132 <- 2550439 | ANTI-PHASED | CONSISTENT |
| 5239260 <- 5272413 | ANTI-PHASED | CONSISTENT |
| 8106314 <- 8176552 | ANTI-PHASED | CONSISTENT |

All six were anti-phased before; five are consistent after. Whole chr20: 1,808
of 55,907 phased records re-oriented (3.2%), with bridged blocks, record count,
block count and read accuracy all unchanged -- 127 bridged, 56,032 records, 281
blocks, 203,751 tagged, 2,851 discordant, 1.399% -- because the reads were always
flipped correctly and only the VCF was wrong.

## Still open: a merge can chain onto a stale phase-set label

The sixth row above is a different defect. Recovery runs per chunk, so a merge
relabels only the candidates in the chunk being processed, and a block spanning a
chunk boundary is relabelled in one chunk and not the other:

```
parent PS=130540  -> final {130540: 6}
parent PS=195305  -> final {130540: 558, 195305: 11}    <- split
parent PS=545002  -> final {195305: 379}
```

`195305` lost 558 sites to `130540` and kept 11, and a later merge then joined
`545002` onto that stale `195305` label instead of onto `130540`. The cost is
contiguity -- a block reported as two, and a join that lands on the wrong label --
not orientation, since each final label is internally consistent. The fix is
alias resolution: follow `keep_ps` to its current label before merging, and
relabel across chunks rather than within one. Not attempted here.
