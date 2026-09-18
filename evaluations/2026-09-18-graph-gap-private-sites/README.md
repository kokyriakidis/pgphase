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

**RETRACTED: the cause is not the representation of `55,883,019`.** The
paragraph that stood here said the record declares net -3 and -2 while the reads
carry -6 and -4, so no read matches either declared allele and the site is
inert. That was wrong, and wrong because of my own display: the strings were
printed through `ref[:4]+'>'+alt[:4]`, which turned `AATATAT>AAT,A` into
`AATA>A,AA`. The record's alleles are **net -4 (`AAT`) and net -6 (`A`)** -- 
exactly the two events the reads carry -- and it is emitted `2|1` with
`AD=0,29,33`. The claim was also false on its own terms even for the misread
alleles: one read does carry net -2, so "zero reads match either declared
allele" was never true, and the table as published dropped the -7 and -5 rows.

The full distribution at the locus, reads covering the window:

| net length | reads | MAT | PAT |
|---:|---:|---:|---:|
| **-6** | 31 (46.3%) | 26 | 5 |
| **-4** | 28 (41.8%) | 1 | 27 |
| -8 | 4 (6.0%) | 4 | 0 |
| -7 | 2 (3.0%) | 2 | 0 |
| -2 | 1 (1.5%) | 0 | 1 |
| -5 | 1 (1.5%) | 0 | 1 |

So the site is well represented and informative: -6 is maternal (26/31) and -4
is paternal (27/28).

**The real cause is an unlinked step.** Re-walking the parity with full allele
strings -- resolving each tool's GT to the actual allele on each haplotype
rather than to "the first ALT", which is a different allele in the two tools --
puts the divergence at `55,883,019` and it persists to the endpoint:

| site | hiphase hap1 | parent | ours hap1 | parent | |
|---|---|---|---|---|---|
| 55,843,827 | G (ref) | PAT 1.000 | G | PAT 1.000 | agree |
| 55,862,239 | T | PAT 1.000 | T | PAT 1.000 | agree |
| 55,862,269 | TCAGTAAAT | PAT 1.000 | TCAGTAAAT | PAT 1.000 | agree |
| 55,883,019 | AAT (-4) | PAT 0.964 | A (-6) | MAT 0.839 | diverge |
| 55,889,113 | T | PAT 1.000 | A (ref) | MAT 1.000 | diverge |

And the chain's consecutive links:

| step | kb | reads observing both | verdict |
|---|---:|---:|---|
| 55,862,269 -> 55,883,019 | 20.8 | **2** | no link |
| 55,883,019 -> 55,889,113 | 6.1 | 37 | no flip (30 same / 7 cross) |

**Both tools agree on the 6.1 kb link.** The divergence is entirely across the
20.75 kb step that only 2 reads observe -- median read length there is 17.7 kb,
so almost nothing spans it. Those 2 reads say no flip, which is hiphase's
answer; we chose flip. (The earlier "FLIP, 15/22, zero contradictions" table for
that 6.1 kb link was computed against the misread allele; with the record's real
ALT1 it is 30 same / 7 cross, no flip.)

So the defect is not the orientation rule -- it had 2 reads -- and not the
representation. It is that one phase set is emitted across a step with no read
linkage, where an unsupported parity costs a whole block's reads while splitting
would cost only contiguity. That is the rule this session established for
zero-spanning-read intervals, applied one level down: to consecutive sites
within a solve rather than to a seam between blocks.

**A minimum-link split does not separate the cases, though.** Minimum
consecutive-link support along each sub-solve's chain across the seam:

| seam | verdict | sites | min link | weakest step |
|---|---|---:|---:|---|
| 45,347,748 | BAD 52% | 8 | 9 | 13.2 kb |
| 55,843,827 | BAD 76% | 11 | 5 | 17.4 kb |
| 12,717,796 | ok 100% | 16 | **3** | 18.9 kb |
| 41,866,917 | ok 100% | 13 | 11 | 13.6 kb |
| 12,411,075 | ok 100% | 11 | 5 | 17.0 kb |
| 22,625,517 | ok 96% | 8 | 27 | 9.7 kb |

A threshold that rejects the 76% seam at 5 also rejects two correct joins at 3
and 5, and the 52% seam passes at 9. So minimum link support is not the gate
either, and the 52% seam must fail for a different reason than the 76% one --
its weakest link has 9 reads, so its sub-solve being at chance throughout
(0.54/0.68) is not explained by linkage at all.

