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
