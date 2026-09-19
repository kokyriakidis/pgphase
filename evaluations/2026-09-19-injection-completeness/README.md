# Does every alignment-path site in the gap reach the output, and do we hold every competitor site?

Two windows, each run twice: the alignment arm (`collect-bam-variation`, BAM
only) and the graph arm with recovery routing (`--graph-noisy-msa
--stitch-recovered --recovery-audit-out`), plus hiphase over the same span.
Only heterozygous phased records strictly inside the gap are compared.

## 1. Injection is complete except for one record

Matching on the full record `(POS, REF, ALT)` -- not position alone, because a
multiallelic locus emits two rows at one position:

| window | alignment het records | reaching the graph arm | missing |
|---|---:|---:|---:|
| `chr20:22,980,600-23,008,891` | 3 | 3 | **0** |
| `chr20:55,313,902-55,358,363` | 6 | 5 | **1** |

The single loss is `55,336,460 CGTGT>C 1|0` (AD 3,26) -- the deletion half of
the complementary pair whose insertion half `C>CGT 0|1` IS emitted. That is the
open defect already traced: nothing removes it, and it is lost inside the
writer's per-allele loop.

Everything else arrives with the same alleles AND the same haplotype
assignment. Position-only matching reported "0 missing" for both windows and
was wrong -- it matched the deletion against the insertion at the same
position.

## 2. We hold every site the competitor uses

Three hiphase heterozygotes looked absent on a position match and all three are
present, re-represented:

| hiphase record | ours |
|---|---|
| `22,985,065 T>TT` | `22,985,062 G>GT` -- same 1 bp insertion, anchored 3 bp earlier |
| `55,331,021 T>TATATATATATATATA 1\|2` | `55,331,014 G>GTATA...(21)` + `G>GTATA...(27)`, the same two alleles as two rows |
| `55,336,493 GTGTG>G 1\|0` | the `55,336,460` pair -- insertion emitted, deletion the loss above |

So competitor site supply is not a deficit in either window. Apparent absences
were anchor offsets of a few bp and multiallelic-versus-split representation.

## 3. The graph arm never emits a multiallelic record

Systematic, and quantified over whole chr20:

| arm | phased records | records with a comma in ALT | positions carrying 2+ het rows |
|---|---:|---:|---:|
| graph, default | 62,352 | **0** | 400 (800 records) |
| graph, routed | 69,906 | **0** | 653 (1,306 records) |
| alignment, w55 window | 277 | 8 | 7 |
| alignment, w22 window | 96 | 2 | 0 |

Where the alignment arm and hiphase emit one record with two ALTs and `GT=1|2`,
the graph arm emits two biallelic rows with `1|0` and `0|1`. Verified
content-equivalent on full strings:

```
alignment  55331014  G > GTATATATATATATATATATA,GTATATATATATATATATATATATATA  1|2  AD 0,15,15
graph      55331014  G > GTATATATATATATATATATA                              1|0  AD 0,16
graph      55331014  G > GTATATATATATATATATATATATATA                        0|1  AD 0,15
```

Same alleles, same haplotype for each, consistent depths. No information is
lost -- but a record-level comparison against the alignment arm or a competitor
shows a mismatch at every one of these loci, and that decomposition is exactly
what made two comparisons in this session miscount before full strings were
read.
