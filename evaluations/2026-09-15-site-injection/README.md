# Injecting BAM-derived sites into the graph solve

The competitors find sites in the BAM and phase those, so the natural move is to
add BAM-derived sites to the graph's site set and let one solve span the gap.
That mechanism already exists -- it is what `collect-hybrid-variation` does via
`hybrid_inject.cpp` -- so this is a measurement, not a design.

Two arms per window, on the same eight gap windows used for column discovery and
phase transfer, at `-q 1` with a 50 kb flank: the site union alone, and the union
plus the gap-recovery tiers.

## The sites really are there

Counting candidates inside the *true* residual break intervals (first-to-last
phased-variant spans, not read starts):

| source | clean het SNPs | noisy candidate hets | repeat het indels |
|---|---:|---:|---:|
| BAM channel | 23 | 79 | 0 |
| graph channel | 9 | — | 74 |

The BAM retrieves 2.5x the clean het SNPs the graph has where phasing breaks.
Unevenly, though: four of the six mid-size gaps get only 0-2 extra clean het
SNPs in their break, while the two 200 kb+ gaps get 6 and 12.

## But no single phase set spans any gap

**0 of 8 gaps spanned, in both arms.** Every window still ends with 2-7 phase
sets and none covering the gap interval end to end. The blocks that do form are
not wrong -- the target block scores 97.9-100.0% against truth across the sixteen
window-arms (the union arm is 100.0% in all eight; the floor is 97.9% at gap
35,919,404 with recovery) -- so this is a linkage failure, not a corruption.

What injection does do is narrow the break. Largest uncovered stretch inside each
gap, by approach:

| gap | gap bp | graph-only | phase transfer | injection | injection + recovery |
|---|---:|---:|---:|---:|---:|
| 7,073,919 | 44,697 | 44,697 | 44,698 | 44,697 | **21,262** |
| 7,163,303 | 43,347 | 43,347 | **20,038** | 43,347 | **20,038** |
| 10,296,487 | 47,108 | 47,108 | **12,622** | 25,260 | 19,723 |
| 13,429,829 | 201,996 | 201,996 | 189,935 | **57,744** | 137,367 |
| 16,673,759 | 41,482 | 41,482 | **22,624** | 41,482 | **22,624** |
| 35,919,404 | 236,915 | 236,915 | 227,126 | **43,482** | 56,364 |
| 36,332,599 | 48,420 | 48,420 | **8,041** | 48,420 | 10,898 |
| 60,084,884 | 47,398 | 47,398 | **16,416** | 24,472 | **16,416** |
| **total** | **711,362** | **711,362** | 541,481 | 328,882 | **304,692** |

Injection transforms the two large gaps -- 202 kb down to 58 kb and 237 kb down
to 43 kb, where transfer barely touched them -- which tracks their 6 and 12 extra
clean het SNPs. On the mid-size gaps injection alone changes nothing in four of
six, because their breaks contain 0-2 usable sites: the union has nothing there
to link with. Union plus recovery is the best single configuration overall
(711 -> 305 kb, a 57% reduction) but still closes no gap.

Anomaly worth chasing: `--recover-gaps` makes the two large gaps *worse* than the
plain union (57,744 -> 137,367 bp and 43,482 -> 56,364 bp) while improving the
mid-size ones. Recovery appears to fragment coverage it should be extending.

## Where this leaves the approach

Site retrieval from the BAM is necessary and already in place; it is not
sufficient. The whole-chromosome evidence says the same thing -- the hybrid path
is this union, and it still leaves 164 unresolved gaps spanning 15.2 Mb.

The mid-size breaks are genuine site deserts, so the remaining options are the
two that do not need a called site: the GAF `cs:Z:` private-variant channel
(present on every record, keyed on node and offset), or the pgbam haplotype-thread
stitching the CLI already exposes -- `--pgbam-file` with the
`--pgbam-*-min-winning` thresholds, documented as the fallback "when common-read
signal is absent", using the `hs`/`hb`/`he` GBWT thread tags. The latter is the
pipeline's own designed answer for exactly this case and has not been tested here.
