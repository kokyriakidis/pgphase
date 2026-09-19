# Cross-chunk seams have no headroom; within-chunk breaks have all of it

Deleting the post-hoc pass left cross-chunk seams with no recovery, and that was
recorded as headroom. It is not. Measured on the single arm's own chr20 output
(commit 2c8e99b, 333 phase blocks, 330 breaks between consecutive block
extents), classifying a break as cross-chunk when a 500 kb chunk boundary falls
strictly inside it:

| break class | count | bp | bridgeable (>= 3 reads spanning the break) |
|---|---:|---:|---:|
| cross-chunk | 34 | 5.02 Mb | **0** |
| within-chunk | 296 | 7.19 Mb | **106** |

Not one cross-chunk seam has three reads crossing it. A seam-only recovery pass
would have nothing to work with, so it is not worth building: the chunking is
not costing contiguity, because chunk boundaries fall where reads do not reach
anyway. (An earlier count of 43 "seam" breaks came from testing the MIDPOINT of
a break for proximity to a boundary, which labels a 682 bp break 4.7 kb from a
boundary as cross-chunk when both its sides sit in the same chunk. The strict
test is whether a boundary lies between the two blocks.)

## Where the fragmentation actually is

106 within-chunk breaks have reads crossing them, many trivially:

| break | width | spanning reads |
|---|---:|---:|
| 26,624,496 - 26,625,001 | 505 bp | 100 |
| 11,157,014 - 11,157,086 | 72 bp | 90 |
| 45,329,938 - 45,330,583 | 645 bp | 80 |
| 22,565,006 - 22,565,036 | **30 bp** | 71 |
| 8,946,171 - 8,946,528 | 357 bp | 69 |

A 30 bp break with 71 reads across it is a defect, not a limitation. The
candidate table says what is happening: the site after each break opens a new
phase set at its own position while the site before it belongs to a long
running one, and the graph channel's depth at these sites is far below the read
depth -- DP 13 and 10 across a break 100 BAM reads span, DP 7/9/23 across one
that 71 span. The link between two sites is decided on reads carrying an allele
at both, and in the graph channel there are too few.

These breaks are NOT a missing-site problem, which is what recovery fixes: the
flanking sites exist and are phased, and a 30 bp gap has nothing to merge into
it. What is missing is a link, and the evidence for it is in the alignment.

So the next piece is not a seam pass over chunk boundaries but a bridge over
within-chunk block seams that scores the two flanking sites against the reads
that carry both -- the same decision the deleted post-hoc pass made with a vote,
applied where reads actually exist. Upper bound if every bridgeable break
joined: 333 blocks -> 227.

## The bug in the merge, and why it was not the cause

`retry_unphased_windows_in_place` collects the sub-solve's read observations for
EVERY site it sees, including sites the parent chunk already owns, and the
transfer below it already merges an alignment allele onto a parent candidate and
adds reads the parent never had. But the function returned early --
`if (new_cands.empty()) return 0;` -- so a seam whose two sides are already
called, with nothing new between them, threw that evidence away. Fixed: the
early return now also requires that no parent site gained an observation.

Measured on whole chr20: **VCF byte-identical**, 219,061 tagged against 219,059,
same 323 read blocks, same 2,543 misplaced, 1.161%, 120 s. One additional chunk
runs the merge (125 against 124). So the bug is real and the fix is right, but
it is not what holds the 106 breaks open.

## What the 106 breaks are NOT

Two hypotheses tested and rejected:

- *The blocks start at unsupported catalog sites.* No: median DP at the 333
  block-start sites is 62, against 63 for all 76,621 candidates, and 62.5 at the
  starts following a bridgeable break. Low depth is mildly enriched at block
  starts (9% below DP 10 against 2% overall) but describes a minority.
- *The alignment channel would stitch them.* No: running the alignment arm over
  four of these breaks reproduces the break at the same position in two cases
  (its own previous block ends 29 kb earlier at 11,157,086) and does not call
  the left site at all in the other two.

What the link probe does show is that the linker is working from far fewer reads
than overlap the locus: every link around chr20:26,624,4xx is decided on
`overlap_reads=13`, and no link into 26,625,001 is attempted at all, while the
BAM has 114 and 100 reads across those positions. DP is not that number -- a read
counts for a link only if it carries a CALLED ALLELE at both sites, and in the
graph channel most reads covering a site do not.
