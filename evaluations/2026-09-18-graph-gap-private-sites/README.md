# Can a targeted alignment run recover a graph-only gap? One gap, measured

The graph pass on whole chr20 leaves **373 seams between consecutive blocks,
19.71 Mb in total**. Test gap: `chr20:32,024,200-32,206,530` -- 182.3 kb, 297
reads, a 486-site block on its left, and **283 heterozygotes phased inside it by
hiphase**, so the interval is genuinely heterozygous.

## What the alignment finds there

`collect-bam-variation` over the seam plus 30 kb flanks, classified against the
catalog by `(POS, REF, ALT)` and by position:

| category inside the seam | private | in catalog |
|---|---:|---:|
| `CLEAN_HET_SNP` | **10** | 10 |
| `NOISY_CAND_HET` | **32** | 44 |
| `CLEAN_HOM` | 4 | 7 |
| `NOISY_CAND_HOM` | 11 | 20 |

**42 private heterozygotes, of which 27 segregate at >= 0.90 against read
truth** and 15 do not. So yes: the alignment finds phaseable sites the catalog
does not have, and the answer to the question as asked is yes.

## But private sites are not the reason the gap exists

The seam contains **5,512 catalog records**, every one carrying an `AT` walk
annotation, and **334 GAF records** overlap it -- more than the 179 over the
28 kb flank that phases fine. Yet the graph pass carries **2 candidates** there,
both `REP_HET_INDEL`. Reproduced by running the graph pass on the interval
alone: 2 candidates from 5,512 records, against 119 candidates (104
`CLEAN_HET_SNP`, all phased) in the adjacent flank.

The sites are **absent from the candidate table**, not demoted within it, so
they are lost inside `build_graph_chunk` before classification -- its per-alt
`min_alt_depth` cut, the biallelic decomposition's depth/AF filters, or the
parent-gated observation counting. The interval is **100% soft-masked** (182,331
bp, 0% N, GC 40.4%), which is consistent with reads' walks through repeat snarls
not matching any catalog allele walk exactly, leaving the sites with no
observations at all. That last step is narrowed, not proven: it needs a count of
the `GraphReadAllele` rows the query returns for the seam against the flank.

So the recovery this gap needs is not primarily private-site discovery. It is
re-measuring the catalog's own sites from the reads, which the alignment channel
does successfully in the same interval -- it phases **97 heterozygotes inside
the seam** where the graph pass phases 2.

## And it would not bridge this gap

The alignment-only solve over the seam does not span it either: three blocks of
28-30 kb (`31,994,598-32,022,665`, `32,174,519-32,204,115`,
`32,204,477-32,234,748`). A recovery here would add interior blocks and interior
sites, not a bridge. Whether those interior blocks are worth having depends on
the phasing metric: they add phased sites and tagged reads without joining the
flanks.

## What this says about a recovery mode

Worth building, with the scope set by the measurement rather than by the name:
the targeted run should admit **both** the catalog's sites re-measured from
reads and the alignment's private sites, because in a repeat-masked interval the
catalog's sites are the larger loss (5,512 unusable against 42 private found).
The blocker recorded for the graph path still applies -- its output table is
index-parallel to per-site metadata, so appended candidates need a synthesized
`GraphSiteMeta` each and no reordering.

Not measured: whether the same pattern holds across the other 372 seams, and
whether the 19.71 Mb of seam is mostly soft-masked.

## A target with the competitor's span verified correct

The first gap was picked on competitor *presence*. This one is picked on
competitor *correctness*: of the 373 seams, 84 have at least 20 phased sites on
both flanks and span 20-400 kb, and **4 of those have a competitor block that
spans the seam with both flanks scored independently at >= 0.95 against read
truth and agreeing with each other**.

| seam | kb | tool | accuracy | scored reads | left flank | right flank |
|---|---:|---|---:|---:|---|---|
| 10,706,319-10,760,419 | 54.1 | hiphase | 98.26% | 574 | MAT 1.000 | MAT 1.000 |
| **12,411,075-12,458,901** | **47.8** | **hiphase** | **100.00%** | **625** | **PAT 1.000** | **PAT 1.000** |
| 10,296,487-10,343,595 | 47.1 | hiphase | 98.17% | 602 | MAT 1.000 | MAT 1.000 |
| 58,703,912-58,747,200 | 43.3 | hiphase | 99.61% | 509 | MAT 0.993 | MAT 1.000 |

## What recovery gets on that target

`chr20:12,411,075-12,458,901`, 47,827 bp, 50.9% soft-masked, 286 reads.

| what is there | count |
|---|---:|
| catalog records inside the seam | 686 |
| graph-only candidates inside the seam | **2** |
| alignment candidates inside the seam | 13 |

Same collapse as the first gap -- 2 of 686 catalog records reach the graph pass's
table -- and this seam is only half soft-masked, so the pattern is not confined
to fully masked intervals.

Of the 13 the alignment finds: 4 `CLEAN_HET_SNP` (**2 private**, 2 catalog), 2
`NOISY_CAND_HET` (1 private, 1 catalog), 7 `CLEAN_HOM`. Scored against read
truth, **2 of 3 private hets and 2 of 3 catalog hets are informative**.

**The alignment solve spans the seam, correctly:** one block
`12,382,782-12,488,586`, 105.8 kb, 143 sites, and per-block truth scoring gives
**100.00% over 508 tagged reads, 0 discordant**, with the left flank MAT 1.000
(n=116) and the right flank MAT 1.000 (n=117) -- **consistent across the seam**,
so the join is earned rather than a coin flip. Hiphase reaches the same verdict
independently on the same interval.

**And it phases inside:** 6 heterozygotes within the seam, which cuts the
largest unphased hole from **47.8 kb to 17.0 kb**.

So on this target recovery delivers all three things: the seam is bridged with a
verified-correct orientation, the unphased span shrinks by 64%, and six sites
that were unphased become phased -- four of them informative by truth, two of
them private to the alignment.

Not measured: whether the other three verified seams behave the same way, and
whether the bridge survives being stitched into the graph pass's blocks rather
than standing alone.

## Why a seam with every site present still inverts

`chr20:55,843,827-55,889,113`, which graph+recovery merges at **76.06%** while
hiphase spans it at **99.83%** (both flanks consistent, 0.99/1.00).

**Site supply is not the problem.** Hiphase phases 5 heterozygotes inside the
seam; we hold a candidate for **5 of 5** and our sub-solve phases **5 of 5** --
3 from the catalog, 2 private to the alignment. Only 2 of the 5 are informative
by read truth (the two endpoints, both at 1.000).

**The failure is one link.** Scoring each site's parity against truth, our
output agrees with hiphase at the left endpoint (hap1 = PAT) and at all three
interior sites, then disagrees at the right endpoint `55,889,113` -- hiphase
`1|0` (hap1 = PAT), ours `0|1` (hap1 = MAT). That site is the right flanking
block's terminal site, so the whole right parent block goes with it, which is
the ~24% of misplaced reads. Read hamming counts reads, so one bad edge at a
block boundary costs a block's worth of reads while the site count looks
perfect.

**And the link's read evidence is unanimous the other way.** 37 reads observe
both `55,883,019` and `55,889,113`:

| 55,883,019 | 55,889,113 | reads | parent |
|---|---|---:|---|
| REF | REF | 15 | all MATERNAL |
| REF | ALT | 22 | all PATERNAL |

Zero contradicting reads. So the orientation *rule* is not at fault: the
decision was not made on this evidence at all.

**The cause is the representation of `55,883,019`.** The record is
`AATA>A,AA` -- ALT1 is net -3, ALT2 is net -2. The reads carry neither:

| net length | reads | MAT | PAT |
|---:|---:|---:|---:|
| **-6** | 31 (46.3%) | 26 | 5 |
| **-4** | 28 (41.8%) | 1 | 27 |
| -8 | 4 | 4 | 0 |
| -2 | 1 | 0 | 1 |

**No read matches either declared allele**, so every read scores REF there and
the site carries no phase information -- which is why both ALTs give an
identical link table above. The locus is a tandem-repeat deletion whose two
haplotypes are -6 (maternal, 26/31) and -4 (paternal, 27/28): informative if
represented correctly, and both non-reference, which is the multiallelic class
whose representation was the subject of the earlier session work. Represented as
-3/-2 against REF it is inert, the chain has nothing at that step, and the
parity across the last 6.1 kb is then decided by something weaker.

So this seam is a **representation** defect, not a k-means or stitching defect.
The two measurements that would close it: emit the locus with its real alleles
(-6 and -4) and re-check that the link resolves to the flip truth requires, and
confirm that the sub-solve's parity at `55,889,113` follows.
