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
