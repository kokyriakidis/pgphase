# One arm, iterated

`arm.sh` is the only configuration. It is the alignment machinery run over the
union of the alignment's own candidates and the graph's catalog sites, one solve
per chunk, with nothing added on top: no `--recover-gaps` (the second pass zeroes
every non-graph candidate's category before phasing and skips the noisy-region
MSA), no retry, no margin filter (`min_read_hap_margin` defaults to 0 and the
hybrid subcommand does not override it), `--keep-noisy-kmeans` because
`skip_noisy_kmeans = true` is a hybrid-specific override the alignment pipeline
does not apply, and `-q 1`.

```sh
OUT=/tmp/arm REGION='CHM13#0#chr20:48126830-48279446' \
  ./evaluations/2026-09-16-single-arm/arm.sh
```

## Baseline on chr20:48,176,830-48,229,446

2 het blocks, no span, 6 in-gap phased hets, 393 reads tagged, **100.00%**
concordant against read truth. Blocks `48,149,548-48,229,226` (8 sites) and
`48,229,446-48,279,445` (49 sites).

Both site defects diagnosed under the recovery arm are **absent** here:
`48,177,780` is emitted `0|1` (the recovery arm emitted `NOISY_CAND_HOM 1|1`),
and the two sites carrying no haplotype information are not emitted at all.

## Iteration 1: the crossing site was genotyped homozygous

The block boundary sits where the chain has nothing to link across. Inside the
41.8 kb interval `48,183,976-48,225,787` the arm holds:

| source | in the interval |
|---|---|
| graph catalog | **500 sites** |
| arm candidates | 59 `CLEAN_HOM`, 1 `CLEAN_HET_INDEL`, 2 `NOISY_CAND_HET`, 6 `NOISY_CAND_HOM` |
| arm candidate at `48,204,384` | DEL, DP 71, **30 ref / 41 alt, AF 0.5774**, `NOISY_CAND_HET` |
| arm emitted record at `48,204,383` | `AT>A` **`GT=1|1`** |
| hiphase record at `48,204,383` | `AT>A` **`GT=0|1`** |

So the signal is not missing from discovery. The site hiphase crosses on is
found, classified a noisy het, and measured at a textbook heterozygous depth --
and then emitted homozygous, which is why the chain has no phased heterozygote
inside the interval and links across it with no spanning read.

`iter_update_var_hap_to_cons_alle` recomputes each haplotype's consensus
**independently by majority**, so both haplotypes select the deeper allele. The
multi-allele branch beside it already guards this case, with the comment
*"Independent haplotype majorities can select the same allele twice."* The
predicate that applies the joint orientation existed but was scoped to retry
windows; `--joint-het-orientation` applies it wherever a biallelic candidate's
own allele depths call it heterozygous, under the same conditions the predicate
already tested (category het, both haplotype profiles present, both depths above
`min_alt_depth`, allele fraction within `[min_af, max_af]`).

| | baseline | `--joint-het-orientation` |
|---|---|---|
| `48,204,383` | `1\|1` | **`0\|1`** (hiphase agrees) |
| `48,202,056` | absent | `0\|1` |
| in-gap phased hets | 6 | **9** |
| blocks / spans | 2 / no | 2 / no |
| reads tagged | 393 | 393 |
| concordant | **100.00%** | **100.00%** |

Per-site orientation shows the interval is now crossed **correctly**: the chain
runs `48,183,976` (1.000) to `48,204,383` (0.958), both PATERNAL-on-hap1 like the
rest of the body, over links of 18.1 kb / 7 reads, 2.3 kb / 58 reads and 21.4 kb
/ 7 reads. Off by default; measured on this window only.

## Iteration 2, the next target

One switch remains, and it is now isolated: `48,204,383` is PATERNAL at 0.958
and `48,225,786` (`CAAAA>C`, consistency **1.000**) is MATERNAL, so the
orientation flips across the 21.4 kb link carried by 7 reads. Everything from
`48,225,786` onward, including the right block, is MATERNAL and mutually
consistent.

Separately, the 220 bp boundary `48,229,226 -> 48,229,446` reports chain link
evidence `agree=0 conflict=0` despite **60** reads spanning it, because both of
the left block's terminal sites are homopolymer indels and those are excluded
from read scoring. That exclusion must not be lifted before the switch above is
fixed: doing so earlier produced a spanning block at 59.29% accuracy.
