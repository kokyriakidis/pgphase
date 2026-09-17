# A fixed panel of windows we leave unphased and a competitor phases correctly

Every iteration before this re-derived its own target, so results were not
comparable across attempts. This is the standing panel: six windows drawn from
the stock-default deficit set where a competitor spans the gap at **100.00%**
read-level accuracy over at least 150 scored reads, spread across chr20 and
ranging 30-45 kb.

```sh
OUT=/tmp/panel-stock ./evaluations/2026-09-16-test-panel/run_panel.sh
OUT=/tmp/panel-x FLAGS="--some-flag" ./evaluations/2026-09-16-test-panel/run_panel.sh
python3 evaluations/2026-09-16-test-panel/score_panel.py \
  --panel evaluations/2026-09-16-test-panel/panel.tsv \
  --arm /tmp/panel-x --baseline /tmp/panel-stock \
  --truth-map /tmp/truth_hap.tsv --output /tmp/arm_x.tsv
```

`FLAGS` empty is the stock pipeline. The scorer reports, per window, whether one
block spans the gap, how many heterozygous sites are phased inside it, reads
tagged, read concordance against the diplinator truth, and -- against a baseline
arm -- the gate that decides admissibility: a read that was concordant and is now
discordant is a regression no coverage gain offsets. Concordance is computed per
phase set with each block oriented by its own majority, so a switch inside a
block lowers the score instead of hiding.

## Baseline: stock defaults

| gap_left | gap_bp | blocks | spans | in-gap hets | tagged | concordance |
|---|---:|---:|---|---:|---:|---|
| 48,183,976 | 45,470 | 2 | no | **2** | 393 | 100.00% |
| 55,843,827 | 45,286 | 2 | no | **2** | 413 | 99.76% |
| 24,105,188 | 37,099 | 3 | no | **2** | 495 | 100.00% |
| 5,309,406 | 35,679 | 3 | no | **2** | 442 | 99.55% |
| 12,717,796 | 34,495 | 2 | no | **2** | 529 | 98.87% |
| 39,838,293 | 30,518 | 2 | no | **2** | 522 | 100.00% |

**Every window has exactly two in-gap phased heterozygotes**, and on the two
windows examined in detail those two are the gap's own boundary sites. So on
stock defaults the pipeline phases **no interior site at all** in any of the six,
while the competitor crosses each one. That uniformity is the panel's value: it
is one signature, not six unrelated failures.

## Iteration 1: `--retry-unphased-with-bam`

| gap_left | in-gap hets | tagged | concordance | c->d | new conc / disc |
|---|---:|---:|---|---:|---:|
| 48,183,976 | 2 -> **8** | 393 | 100.00% | 0 | 0 / 0 |
| 55,843,827 | 2 -> **10** | 413 -> **501** | 99.76% -> 97.80% | **1** | 79 / 9 |
| 24,105,188 | 2 -> 2 | 495 | 100.00% | 0 | 0 / 0 |
| 5,309,406 | 2 -> **4** | 442 | 99.55% | 0 | 0 / 0 |
| 12,717,796 | 2 -> **3** | 529 | 98.87% | 0 | 0 / 0 |
| 39,838,293 | 2 -> **9** | 522 | 100.00% | 0 | 0 / 0 |

Panel: in-gap hets **12 -> 36**, tagged 2,794 -> 2,882, concordance 99.68% ->
99.34%, gate 1 concordant->discordant and 79 new concordant against 9 new
discordant. **Still zero windows spanned.**

Three things this says that the single-window work could not:

1. **Admitting the sites works and generalizes** -- five of six windows gain
   interior heterozygotes, 12 -> 36 across the panel.
2. **It is nowhere near sufficient.** No window spans, and in four of the six the
   read tags do not move at all: sites become phased, no read gains a tag. So the
   missing step is linkage, not discovery, in most of the panel.
3. **The panel's own cost is negligible** -- nine new discordant reads and one
   flip -- while chromosome-wide the same flag takes read hamming from 0.559% to
   2.723%. That is 0.559% of 212,320 against 2.723% of 218,989, i.e. about
   1,190 discordant reads against about 5,960: roughly **4,800 extra**
   discordant reads, not the 5,300 an earlier version of this line stated.
   Six windows cannot produce that,
   so the chromosome-wide harm comes from windows the detector opens *outside* the
   deficit set. That relocates the next fix from the admission gate to the
   detector's scope, which the single window could not have shown.

Any default change is measured here first and against the 0.559% chromosome-wide
baseline second; the panel is cheap enough to run every iteration, the whole-chr20
arm is not.
