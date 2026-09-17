# Gather the sites the reads in the gap participate in

The proposal: instead of taking every site in a gap corridor, take the reads
that cross the gap and keep the alignment-derived sites those reads participate
in. Measured on the six panel windows.

## No read spans any gap

| gap | gap | longest read | reads overlapping |
|---|---:|---:|---:|
| 48,183,976 | 45.5 kb | 29.9 kb | 269 |
| 55,843,827 | 45.3 kb | 31.1 kb | 252 |
| 24,105,188 | 37.1 kb | 30.2 kb | 229 |
| 5,309,406 | 35.7 kb | 30.4 kb | 226 |
| 12,717,796 | 34.5 kb | 29.3 kb | 219 |
| 39,838,293 | 30.5 kb | 30.0 kb | 194 |

So "reads that span the gap" is an empty set in every window, and the gap can
only be crossed by a chain of overlapping reads. The usable question is the
sites the *overlapping* reads participate in.

## What that selection yields

| gap | sites in gap | give a read partition | informative (seg >= 0.90) |
|---|---:|---:|---:|
| 48,183,976 | 69 | 13 | 3 |
| 55,843,827 | 48 | 13 | 3 |
| 24,105,188 | 19 | 4 | 2 |
| 5,309,406 | 150 | 11 | 3 |
| 12,717,796 | 56 | 9 | 2 |
| 39,838,293 | 9 | 8 | 7 |

## Those sites alone do not chain -- but the missing ones are not missing

Chaining flank to flank through the informative sites only, the largest step in
five of six windows is crossed by **zero** reads:

| gap | largest step | reads crossing it |
|---|---:|---:|
| 48,183,976 | 41.8 kb | 0 |
| 55,843,827 | 26.8 kb | 0 |
| 24,105,188 | 37.1 kb | 0 |
| 5,309,406 | 21.9 kb | 0 |
| 12,717,796 | 34.5 kb | 0 |
| 39,838,293 | 8.5 kb | 27 |

Because "informative" is scored against read truth, that is an oracle bound: no
site-selection rule over this set can cross those five steps.

The competitor does cross them, so the set is incomplete -- and every site it
uses there is already in our candidate table. Of its 14 phased heterozygotes
inside the five steps, **0 are absent from our candidates**, 10 are informative
and 4 are phantoms:

| gap | competitor site | segregation | our category |
|---|---|---:|---|
| 48,183,976 | 48,183,976 `CCT>C` | 1.000 | `CLEAN_HET_INDEL` |
| 48,183,976 | 48,204,383 `AT>A` | 0.606 | `NOISY_CAND_HET` (phantom) |
| 55,843,827 | 55,883,019 `AATAT>A` | 0.597 | `NOISY_CAND_HET` (phantom) |
| 55,843,827 | 55,889,113 `A>T` | 1.000 | `CLEAN_HET_SNP` |
| 24,105,188 | 24,121,713 `C>CTTTT` | 0.757 | `NOISY_CAND_HET` (phantom) |
| 5,309,406 | **5,315,591 `C>CT`** | **1.000** | **`NOISY_CAND_HET`** |
| 5,309,406 | **5,331,265 `TAGAC>T`** | **1.000** | **`NOISY_CAND_HET`** |
| 12,717,796 | **12,721,112 `A>AT`** | **0.945** | **`NOISY_CAND_HET`** |
| 12,717,796 | 12,735,894 `TA>T` | 0.520 | `NOISY_CAND_HOM` (phantom) |

(The four flank-anchor rows at 24,105,188 / 24,142,287 / 5,309,406 / 12,717,796
/ 12,752,291 are clean SNPs we already use.)

So the bridging sites are **present, measured, and excluded by class**: the
hybrid sets `skip_noisy_kmeans`, so `NOISY_CAND_HET` never enters its solve. The
constraint is admission, not discovery and not read length.

One nuance worth keeping: these were partitioned using the COMPETITOR's REF/ALT.
At 5,315,591 our own record's alleles are not identical, and a site can
partition cleanly under one allele representation and not under another, so
"present in our table" is not the same as "usable exactly as we hold it".

## What this does and does not justify

It does not justify admitting the noisy class chunk-wide: that is the measured
bad trade, 0.559% to 2.723% read hamming. It does say that the window-scoped
admission already shipped -- `--retry-unphased-with-bam` -- is working on the
right material, and that its 4-of-6 panel spans at 0 concordant-to-discordant
reads are limited by which sites it admits inside the window rather than by any
shortage of evidence.
