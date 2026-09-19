# Letting a repeat-context indel earn admission, and why it still regresses

## The gap that motivated it

`chr20:22,980,600-23,008,891`. Hiphase spans it with ONE block, 169 scored
reads, **99.4% correct** -- a join truth supports, not a competitor artefact.
We split it in two.

Hiphase uses 5 het sites inside the window. We phase 3 and reject 2, and both
rejections are `REP_HET_INDEL`:

| site | hiphase | us |
|---|---|---|
| 22,980,600 `G>T` | 1\|0 | phased |
| 22,985,065 `T>TT` | 1\|0 | **rejected REP_HET_INDEL** |
| 22,989,976 `CAAAAAAAA>C,CAAAAA` | 0\|1 | phased |
| **23,007,537 `TA>T`** | 0\|1 | **rejected REP_HET_INDEL** |
| 23,008,891 `T>G` | 1\|0 | phased, in the NEXT block |

23,007,537 is discovered from the alignment at DP 67, 40 alt against 27 ref,
AF 0.60, and demoted purely because it sits in an A-run. By read truth it
segregates 0.875 -- its reference allele is 25/25 paternal -- and it links into
the next block on 49 reads (28 maternal / 21 paternal, 8 conflicting). The site
we DO keep two steps earlier, 23,001,780, segregates 0.660 and does not
separate the parents at all. The rule that rejected it reads the reference
context and never the reads.

## The criterion is sound at site level

Scoring 218 `REP_HET_INDEL` candidates at DP >= 20 on chr20:20-25 Mb against
read truth: **103 are informative** (segregation >= 0.90). Nearly half the class
is signal.

| gate | admitted | informative | precision | recall |
|---|---:|---:|---:|---:|
| reject all (shipped) | 0 | 0 | -- | 0% |
| read-concentration (top-2 net-length classes >= 90%) | 81 | 53 | 65% | -- |
| **link purity >= 0.90 over >= 15 co-observing reads** | 100 | 95 | **95%** | **92%** |

The control behaves: `CLEAN_HET_INDEL` sites at the same threshold are 98%
informative, so the statistic is not an artefact.

## Delivering it as a category promotion does NOT work

`promote_link_supported_repeat_indels` (`graph_bam_adapter.cpp`, behind
`--link-earned-repeat-indels`) re-labels an earning site `CLEAN_HET_INDEL`
before the first k-means. On chr20:20-25 Mb that looks like a win -- read
hamming 0.570% -> 0.533% at purity 0.90 and 0.497% at 0.85, with more reads
tagged. Chromosome-wide it is a regression at every threshold:

| arm | tagged | read blocks | misplaced | hamming | VCF blocks | N50 |
|---|---:|---:|---:|---:|---:|---:|
| off (shipped) | 219,061 | 323 | 2,543 | **1.161%** | 333 | 486 kb |
| purity 0.85 | 219,752 | 314 | 4,145 | 1.886% | 322 | 506 kb |
| purity 0.90 | 219,695 | 321 | 3,152 | 1.435% | 330 | 506 kb |
| purity 0.95 | 219,446 | 322 | 4,473 | 2.038% | 332 | 479 kb |

**The response is not monotone in the threshold** -- 0.95 admits the fewest
sites and is the worst arm. That rules out "we admitted some phantoms" as the
explanation. A promoted site becomes a full participant in the chunk's k-means,
so it does not merely contribute its own link: it changes the clustering
trajectory for every read in the chunk, and the sign of that effect is
essentially arbitrary. The 5 Mb region was not representative, which is the
trap this project has fallen into before.

Admitting the site also did not close the motivating gap. At purity 0.85 it is
emitted and joins the RIGHT block (the boundary moves 1.4 kb left), but the
27 kb step from 22,980,600 remains: no read covers both ends of it, so the join
needs the intermediate sites, and 22,985,065 (DP 76, 43/33) is still rejected
while 22,989,976 and 23,001,780 are not candidates at all under this chunking.

## What to do instead

The criterion earns a LINK, so it should buy a link and nothing else: use an
earning site to join or extend blocks across a break, without letting it into
the clustering that decides every read's haplotype. That is a change to the
block-linking path rather than to the category, and it can be gated on the same
evidence this pass already computes.

Shipped OFF by default; flag-off output is byte-identical to the previous
commit (0 body differences over whole chr20).
