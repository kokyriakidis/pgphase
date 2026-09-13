# Generated Phasing Benchmark Report

Baseline frozen: `2026-09-13T06:30:55.909110+00:00`  
Panel: `HG002 CHM13 chr12/chr18/chr20`  
Competitor lock: `96647bf88fba101df40b9daf7d6be5a3ca8458e86bc714a1185664eef5fc426c`

## Pooled Results

| method | assessed pairs | variant Hamming | phased reads | read Hamming | median chr NGC50 (kb) |
|---|---|---|---|---|---|
| graph | 251641 | 0.034% | 811382 | 0.110% | 251 |
| bam | 253738 | 0.515% | 998736 | 1.890% | 278 |
| hybrid | 251679 | 0.290% | 814172 | 0.461% | 295 |
| graph_lock | 251886 | 0.035% | 822021 | 0.123% | 252 |
| whatshap | 261464 | 4.230% | 1028246 | 5.412% | 359 |
| whatshap_opt | 259356 | 2.945% | 1026344 | 3.694% | 502 |
| hiphase | 260676 | 3.514% | 1064369 | 4.377% | 644 |
| longphase | 225783 | 0.805% | 1027359 | 2.339% | 454 |

The frozen LongPhase baseline is the measured `--pb` SNP mode; it does
not include LongPhase's optional `--indels` mode.

## Correct Competitor Bridges

Blocks require >=99% read-truth accuracy, >=50 evaluated reads, and no
variant switch-error interval across the pgphase gap.

| tool | pgphase break reason | bridges | gap bp | competitor sites | graph sites | private sites |
|---|---|---|---|---|---|---|
| hiphase | block_stitching | 86 | 4116394 | 424 | 125 | 22 |
| hiphase | catalog_site_not_candidate | 64 | 873582 | 37 | 0 | 3 |
| hiphase | clean_candidate_unphased | 34 | 829988 | 177 | 0 | 6 |
| hiphase | repeat_indels_excluded | 250 | 5692132 | 405 | 0 | 15 |
| longphase | block_stitching | 44 | 1650334 | 99 | 49 | 18 |
| longphase | catalog_site_not_candidate | 72 | 994019 | 25 | 0 | 3 |
| longphase | clean_candidate_unphased | 30 | 654985 | 101 | 0 | 6 |
| longphase | repeat_indels_excluded | 164 | 3049996 | 85 | 0 | 20 |
| whatshap | block_stitching | 58 | 2738765 | 289 | 84 | 17 |
| whatshap | catalog_site_not_candidate | 52 | 691350 | 8 | 0 | 1 |
| whatshap | clean_candidate_unphased | 23 | 551461 | 54 | 0 | 5 |
| whatshap | repeat_indels_excluded | 167 | 3452221 | 157 | 0 | 12 |
| whatshap_opt | block_stitching | 67 | 3207822 | 317 | 102 | 19 |
| whatshap_opt | catalog_site_not_candidate | 61 | 855457 | 6 | 0 | 2 |
| whatshap_opt | clean_candidate_unphased | 29 | 676174 | 69 | 0 | 5 |
| whatshap_opt | repeat_indels_excluded | 195 | 4021848 | 161 | 0 | 15 |

## Why HiPhase Is More Contiguous

HiPhase correctly bridges 434 graph breaks under the strict accuracy rule. Their median width is 22,575 bp. At these breaks HiPhase phases 1,043 linking sites while pgphase phases 125.

The dominant mechanism is repeat-indel exclusion: 250 bridges (57.6%) contain repeat heterozygous indel candidates but no graph-phased anchor. Another 86 (19.8%) already contain graph-phased sites and fail at block stitching. The remaining breaks are clean candidates left unphased or catalog sites that never become candidates.

### Largest Correct HiPhase Bridges

| region | gap bp | read accuracy | reason | shared hets | HiPhase sites | graph sites | repeat indels |
|---|---|---|---|---|---|---|---|
| chr18:46185393-46302851 | 117458 | 99.83% | block_stitching | 16 | 16 | 2 | 24 |
| chr12:46725702-46838918 | 113216 | 99.21% | block_stitching | 17 | 17 | 2 | 7 |
| chr20:15019294-15130077 | 110783 | 99.71% | block_stitching | 17 | 15 | 3 | 18 |
| chr18:33319321-33429051 | 109730 | 99.32% | block_stitching | 17 | 17 | 4 | 13 |
| chr12:63758735-63866834 | 108099 | 99.85% | block_stitching | 27 | 25 | 3 | 28 |
| chr18:71439208-71538123 | 98915 | 99.88% | block_stitching | 7 | 7 | 6 | 12 |
| chr18:52053305-52147183 | 93878 | 99.97% | block_stitching | 6 | 6 | 2 | 2 |
| chr20:47671540-47762233 | 90693 | 99.80% | block_stitching | 11 | 10 | 3 | 19 |
| chr12:68757844-68841697 | 83853 | 99.82% | block_stitching | 6 | 6 | 4 | 16 |
| chr12:78062641-78138755 | 76114 | 99.59% | block_stitching | 9 | 6 | 2 | 5 |
| chr12:15723767-15798039 | 74272 | 99.26% | block_stitching | 4 | 4 | 1 | 11 |
| chr12:87661162-87735057 | 73895 | 99.37% | block_stitching | 11 | 11 | 4 | 16 |
| chr12:95750003-95822767 | 72764 | 99.90% | block_stitching | 12 | 11 | 2 | 20 |
| chr18:6325477-6398154 | 72677 | 99.82% | block_stitching | 7 | 7 | 2 | 7 |
| chr20:61737962-61810469 | 72507 | 99.15% | block_stitching | 18 | 14 | 2 | 18 |

## Reproduction

```bash
python3 scripts/benchmark_panel.py verify-competitors
python3 scripts/benchmark_panel.py run
python3 scripts/benchmark_panel.py report
```

Normal `run` mode verifies and reuses frozen competitors. It cannot
invoke competitor phasers unless the lower-level driver is explicitly
called with `RUN_COMPETITORS=1`.
