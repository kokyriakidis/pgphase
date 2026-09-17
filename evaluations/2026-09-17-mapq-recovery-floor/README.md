# A gap we discover nothing in, and why

## The window

`chr20:26,029,591-26,088,679`, 59.1 kb, the largest gap in the deficit set.
LongPhase spans it at **99.68% over 313 reads** and phases **120 heterozygotes
inside it**; 39 of the first 40 scored segregate cleanly against read truth, so
the evidence is real.

We phase none of them. Candidate density across the run region, per 10 kb:

```
26,000,000  28   ############################
26,010,000  30   ##############################
26,020,000  33   #################################
26,030,000   0   <-- gap
26,040,000   0   <-- gap
26,050,000   0   <-- gap
26,060,000   0   <-- gap
26,070,000   0   <-- gap
26,080,000   2   ## <-- gap
26,090,000  14   ##############
```

Zero candidates across 50 kb with a sharp boundary. This is not the usual
admission failure -- there is nothing to admit.

## What it is not

| hypothesis | test | result |
|---|---|---|
| assembly gap or masked sequence | reference composition | 0% N, GC 45.6%, ordinary sequence |
| broken alignments | CIGAR inspection | 124 reads, **0 soft-clipped**, 0 with a >= 1 kb indel, median length 15.4 kb |
| reads dropped as too noisy | set `max_var_ratio_per_read` and `max_noisy_frac_per_read` to 0 and rerun | **no change** -- still 2 candidates in the gap |
| the noisy-region machinery swallowing it | list every noisy region in the run | none overlaps the gap |
| variants found but filtered | probe counting digars per 10 kb | no digars there to filter |

## What it is

A probe on `collect_candidate_sites_from_records` counting read START positions
per 10 kb shows the answer: **zero reads start between 26,010,000 and
26,080,000** in the pipeline's view, while samtools reports 124 reads
overlapping the middle of that interval. The reads exist; the pipeline never
loads them.

Their mapping qualities are the reason:

```
reads starting inside 26,010,000-26,080,000: all primary, no secondary,
supplementary, duplicate or QC-fail
MAPQ 3  x236     MAPQ 12 x11     MAPQ 4  x10     MAPQ 9  x4     MAPQ 21 x4
```

`kDefaultMinMapq = 30`, applied in `read_passes_filters`, so none of them enters
a phasing chunk. The gap is a **mapping-quality hole**, and the competitor uses
those reads.

## What admitting them everywhere costs

| arm | in-gap hets | phased in gap | spans | read concordance | across the gap |
|---|---:|---:|---|---:|---|
| default (`-q 30`) | 0 | 0 | no | 99.66% | -- |
| `-q 1` | 293 | 294 | **YES** | 97.61% | consistent (L 0.933 n=90, R 0.995 n=193) |

So `-q 1` closes the gap and the orientation is **correct** -- but the left
flank falls from 100% to 93.3%, because the low-MAPQ reads also re-solve the
flanks, where the pipeline was already right. Coverage bought at the cost of
places that did not need help.

## The recovery floor

`--recovery-min-mapq` (default 1). Reads at or above it but below `--min-mapq`
are parsed into the chunk and immediately marked skipped, so every stage behaves
as if they were never loaded. The re-solve then un-skips those overlapping a
window the first solve left unphased -- `skipped_for_mapq` distinguishes them
from reads skipped for their variant load, which stay skipped.

**One latent bug had to be fixed to make loading them safe at all.**
`collect_bam_output` re-iterates the BAM and walks the chunk's read vector **in
parallel by index**, applying its own `min_mapq` test. With the load floor lower
than that test, the two lists desynchronise and every read's tags are written
onto a different read: measured read concordance fell to **52%**, chance. The
output test now uses the same floor as the load. Today `min_mapq` is always the
load floor, so the bug was unreachable; it is an invariant made explicit.

## Status: the floor is in place and does not yet close the gap

On the closed window `5,309,406` it is effectively inert -- 544 tagged, 1 block,
98.71% concordance, 7 discordant, identical with the floor at 1 or 30, differing
in 2 candidate rows.

On the target window the re-solve reports `woke 305 read(s) below the mapq floor
in 2 window(s)` and **nothing changes**: 0 in-gap hets, no span. The reason is
the same ordering constraint that governs injection and the MSA -- the retry
re-enters at `collect_var_run_phasing`, which is *after*
`collect_var_classify`. Discovery has already run, so the woken reads' variants
are never turned into candidates, and the noisy-region MSA cannot help because
no noisy region covers the gap.

Closing this gap therefore needs discovery to run again over the woken reads,
scoped to the window. That is the next step, and it is not a one-line change:
`collect_candidate_sites_from_records` clears the whole table, so a scoped
append is required rather than a re-run.

## The fix: run the pipeline on the gap as its own chunk

Waking the recovery reads inside the parent chunk discovered the sites without
hurting the flanks -- 0 to 292 in-gap candidates at unchanged 99.66% read
concordance -- but the interior fragmented into 1-, 2- and 6-site blocks,
because it was being solved inside a chunk whose read set is mostly asleep.
That approach was reverted.

Instead the window gets **its own chunk**, at its own floor, solved by
`process_chunk` -- the same function that solves every other chunk -- with
`min_mapq = recovery_min_mapq`. Then the answer is stitched in using the
pipeline's own rule: reads tagged in both solves vote an n11/n12/n21/n22 table
per parent block and `select_stitch_orientation` decides, the same function and
the same default net-margin standard that joins adjacent chunks. Its refusal
carries over: a parent block sharing no read with the targeted solve cannot be
merged, which is the guard a bespoke gap link lacks.

Sites strictly inside the window are imported with the merged phase set, flipped
if the vote says so. Without that the block spans an interval it reports nothing
in.

| arm | spans | in-gap hets | tagged | read concordance | left flank | right flank |
|---|---|---:|---:|---:|---:|---:|
| default before | no | 0 | 293 | 99.66% | -- | -- |
| `-q 1` everywhere | YES | 293 | 293 | 97.61% | 0.933 | 0.995 |
| **targeted solve** | **YES** | **37** | 293 | **99.66%** | **1.000** | 0.995 |

The gap closes with **no read-level cost**: 1 discordant read, the same as
before, against 7 for the global floor. One block now spans 151.1 kb, correctly
oriented. Runtime for the window is 3.6 s.

## What is not yet right

**37 in-gap heterozygotes, against the competitor's 120.** The import is
confined to the detected unphased window and to the one targeted block the vote
bridged; the gap is wider than that window. The interval actually reported runs
26,029,671-26,086,807.

**2 of the 18 scorable in-gap SNP calls do not segregate against read truth**
(16 do; 19 more are indels or too shallow to score). The competitor's calls in
the same interval were 39 of 40 clean, so the interior calls are noticeably
worse than its, even though the block they sit in is correctly oriented and the
reads are tagged as accurately as before.

Neither is a reason to hold the change -- the alternative is reporting nothing
across 59 kb -- but both are the next work, and the in-gap call quality matters
more than the count.
