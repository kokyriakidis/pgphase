# Guarded local exact-MEC retry

> Superseded control flow: the retained result is unchanged, but boundary scope
> is now attempted only when the full problem exceeds the exact variable bound.
> A tied or contradictory full solve no longer triggers the narrower solve. See
> `evaluations/2026-09-23-single-flow-recovery/`.

Date: 2026-09-23

## Reproduced bug

The full 500 kb chunk around `chr20:14,264,549-14,272,741` carries decisive,
consistent orientation evidence:

- the two independent deletion rows relate to the right SNP at 35/35 reads;
- the graph/BAM block gauge is 275/275 for the cross orientation;
- 59 sequence-identical candidates independently select that orientation.

The trusted MEC solver nevertheless abstained. It expanded the exact problem
from the 8 kb seam over the complete neighboring phase block. Unrelated
unphased sites then raised the problem above the 20-variable exact-search bound.
The gauge translation was correct.

## Retained design

The ordinary exact solve still covers the read-connected extent of both atomic
blocks. If it cannot solve that scope, a second attempt restricts the problem to
the selected boundary-anchor interval. This local retry is available only when
sequence-identical candidate votes select one orientation at one-sided exact
binomial `p <= 0.05`. Full reads, both deterministic halves, the source-specific
read gauge, and the candidate anchors must all agree. BAM/BAM joins and second
attachments of an imported block remain prohibited.

The guard matters at `chr20:51.27 Mb`. Its local edge and 56 shared-read votes
look decisive, but it has only one shared candidate and joins graph blocks with
opposite parental orientations. An unguarded local solve made 54 correct reads
discordant in that block. The statistical candidate guard leaves it separate.

## Full chr20 result

The comparison used 500 kb chunks, 16 workers, graph MAPQ 5, the same chr20 BAM,
GAF, graph catalog, and diplinator truth BAM as the retained baseline.

The guarded retry closes the 14.264 Mb target, loses no tracked closure, and
reduces VCF phase blocks from 452 to 419. VCF N50 remains 482,085 bp and all
225,005 reads remain tagged. Direct comparison against parental truth finds no
read changing from correct to incorrect or from incorrect to correct. The formal
evaluator scores 6,258 discordant among 224,985 evaluated reads (97.2185%);
nine additional reads become eligible through the merged phase sets.

The unguarded arm is retained only as a rejection result. It closes the same
tracked target but adds 49 net discordant reads, concentrated at 51.27 Mb.
