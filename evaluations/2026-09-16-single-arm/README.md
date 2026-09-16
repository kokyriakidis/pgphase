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

## Recheck: do we hold all of hiphase's signal, and encode it usably?

Every phased heterozygote hiphase emits in `48,145,000-48,240,000`, matched to
the arm's candidate and emitted record at the same locus (indel anchors differ by
a base or two, so matching is by locus, not by exact position), with each site
genotyped from the alignment and scored against read truth.

**Hiphase emits 19, all in a single phase set `48149548` spanning
`48,149,548-48,235,309`. The arm holds 15 of them and splits the window into
two blocks.**

| verdict | n |
|---|---:|
| used, alleles identical | 12 |
| used, re-represented | 3 |
| candidate only, screened out | 2 |
| absent from our candidates | 2 |

The four we do not use all carry real signal:

| site | hiphase allele | our status | n | segregation |
|---|---|---|---:|---:|
| 48,149,567 | `A>AGAG` | **absent** | 64 | **0.984** |
| 48,173,317 | `TGGG>T` | **absent** | 66 | **1.000** |
| 48,177,725 | `T>TA` | candidate, `REP_HET_INDEL` | 74 | 0.865 |
| 48,234,100 | `C>CCT` | candidate, `REP_HET_INDEL` | 54 | **0.981** |

So two informative sites are never discovered, and two more are discovered and
then screened out by repeat demotion while hiphase phases both.

The three re-represented sites are encoded differently but remain usable:
`48,177,780` is multi-allelic for hiphase (`GAGA>G,GA`, `2|1`) and biallelic for
us (`0|1`, segregation 0.946), so one allele is dropped; `48,225,788`
(`AAAA>A`, `1|0`) is our nested-deletion split into a common `CA>C` emitted
`1|1`, which carries nothing, plus a residual `CAAAA>C` emitted `0|1` at
segregation **1.000**; `48,149,548` differs in anchor and allele only.

### This does not explain the switch

**None of the four unused sites lies inside the 21.4 kb link where our
orientation flips, and hiphase has no heterozygote there either.** It crosses
`48,204,383 -> 48,225,786` on the same 7 spanning reads we have, between two
sites we also have and which are both informative in our own data (0.958 and
**1.000**), and it gets the orientation right where we do not.

So the two findings are separate, and in this order:

1. **The switch is a decision, not a data deficit** -- same sites, same reads,
   different outcome. That is iteration 2.
2. **We do discard usable signal**, but recovering it lengthens and strengthens
   the chain rather than fixing the switch: two sites to discover
   (`48,149,567`, `48,173,317`) and two the repeat screen removes at 0.865 and
   0.981 (`48,177,725`, `48,234,100`).

## Iteration 2: a site excluded from the link list still inherited a phase set

The switch was not a decision made on thin evidence -- it was no decision at
all. `48,225,786` (`CAAAA>C`, in an A run) is `is_homopolymer_indel`, and the
het **link** list excludes such a site unless a recovery homopolymer window is
active:

```cpp
if (... && (var.msa_insertion_alts.empty() || var.gap_link_supported) &&
    (!var.is_homopolymer_indel || gap_hp_link) && !unsupported_gap_indel) {
    is_het[_vi] = true;
    het_var_idx.push_back(_vi);
}
```

The arm has no such window, so the site never entered the list and
`het_rank[_vi]` stayed `-1`. But the emit loop still hands it the running phase
set:

```cpp
const int hi = het_rank[_vi];
if (hi >= 0) {
    phase_set = het_ps[hi];
    if (parity[hi] == 1) std::swap(var.hap_to_cons_alle[1], var.hap_to_cons_alle[2]);
}
var.phase_set = phase_set;      // reached with hi < 0 too
```

So it joined the block carrying whatever orientation its own consensus produced,
with `parity` never applied -- never reconciled against the block it was labelled
into. That is why the switch was invisible in read space: the reads covering the
site belong to the next block, so per-block read concordance stayed at 100.00%
while the site itself sat on the opposite haplotype.

The fix admits a site to the link list when its own allele depths call it a clear
heterozygote (`allele_depths_call_het`, the predicate iteration 1 generalised),
so its orientation is decided by spanning reads like any other het instead of
inherited.

| | before | after |
|---|---|---|
| switches within a block | **1** | **0** |
| blocks over the window | 2 | 3 |
| reads tagged | 393 | 393 |
| concordant | **100.00%** | **100.00%** |
| `48,229,446` onward | own block `PS=48229446` | **merged into `PS=48225786`** |

Blocks after the fix, each internally consistent against read truth:
`48,147,225-48,202,056` (8 sites, all paternal-on-hap1 at 0.946-1.000),
`48,204,383` alone, and `48,225,786-48,279,445` (52 sites, all maternal-on-hap1
at 0.933-1.000). The 220 bp boundary that no merge could close is now joined by
the ordinary chain, because the site on its left finally has a link.

## Iteration 3, the next target

The window now fragments where it used to switch. `48,204,383` is a singleton:
its links to `48,202,056` (2.3 kb, 58 spanning reads) and to `48,225,786`
(21.4 kb, 7 spanning reads) are both refused. `48,202,056` -- the site carrying
no haplotype information, segregation 0.507 -- sits between it and the left
block, so the chain's only route from the body to `48,204,383` runs through a
site whose alleles are noise. Removing that site, or refusing to link through it,
is the next step; hiphase holds no het there at all.

Measured on this window only; `--joint-het-orientation` remains off by default.
