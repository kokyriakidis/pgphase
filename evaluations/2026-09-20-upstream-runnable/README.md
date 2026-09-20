# longcallD now builds and runs here, and the first experiment reframes the parity gap

## The capability

`~/Downloads/longcallD` builds clean with `make -j20` (submodules already
checked out: WFA2-lib v2.3.5-9, abPOA v1.5.5-3, htslib 1.21-14). Binary at
`bin/longcallD`. Parity questions are now direct experiments instead of source
readings.

Its CLI differs from ours in three ways that each cost a wasted run:

| intent | upstream |
|---|---|
| region | **positional argument**, after ref and BAM: `longcallD call ref.fa in.bam "CHM13#0#chr20:30795000-30815000"` |
| `-r` | `--ref-idx` (the `.fai`), **not** region |
| `-v` | prints the version and exits; verbosity is `-V INT` |
| output | `-o FILE`, default stdout |

## Experiment 1: on identical region inputs, the two tools agree exactly

Same reference, same BAM, same region `CHM13#0#chr20:30,795,000-30,815,000`,
both tools run directly:

| | phased records | positions |
|---|---:|---:|
| upstream longcallD | 426 | 424 |
| pgphase alignment arm | 426 | 424 |
| identical `(POS, REF, ALT)` | **426** | -- |
| upstream-only | **0** | |
| ours-only | **0** | |

Perfect record-level parity on this window. That is a much stronger statement
than the chromosome-wide 98.8%, and it localizes the whole remaining gap to
something other than the per-region calling logic.

## Experiment 2: upstream's own output moves with the run span

Counting records inside the SAME fixed inner window (30,795,000-30,815,000)
while widening the span the tool is asked to process:

| span processed | upstream | ours |
|---|---:|---:|
| 20 kb (the window itself) | 426 | 426 |
| 200 kb | **477** | 426 |
| 1 Mb | 426 | 426 |
| whole chr20 | **477** | 426 |

Not thread-dependent: `-t 1` and `-t 8` give the same counts, and the 426 set is
always a strict subset of the 477 set. Upstream's extra 51 records appear only
when the locus is interior to one of its 500 kb chunks with full read context,
and disappear when a chunk boundary lands near it -- in the 1 Mb run its chunks
start at 30,300,000, putting a boundary at 30,800,000, 2.9 kb before the
cluster.

## What those 51 records are

All in one tight cluster, `30,802,900-30,804,730` -- 1.8 kb carrying 46 SNPs,
2 deletions and 3 insertions, one variant every 39 bp, which upstream calls
`0|1` 17 times, `1|0` 18 times and `1|1` 16 times. A hyper-divergent locus.

We form and keep noisy regions there -- `30,802,910-30,804,107` (1,198 bp) and
`30,804,338-30,804,843` (506 bp), both at noisy/total = 16/16, ratio 1.000 --
and reconstruct **17** candidates across `30,802,900-30,804,730`, all 17 of
which we emit. So the region is found, kept and processed; the yield differs.

Upstream's verbose output lists several smaller regions in the same span (a
1,232 bp one ending at 30,802,459, then 417, 321 and 297 bp entries). Which
pipeline stage each of those lines belongs to was NOT established, so the
granularity comparison is a lead, not a result.

## Why this changes the target

31% of the 162 never-discovered positions (51 of them) are records upstream
itself produces in some chunkings and not others, at identical settings. They
are not a stable parity target: matching them would mean matching a
chunk-boundary artifact.

The remaining work on this bucket should therefore start by re-deriving it from
region-scoped runs of both tools -- where, in this window, the two tools agree
exactly -- rather than from the whole-chromosome outputs that carry each tool's
chunking with them.
