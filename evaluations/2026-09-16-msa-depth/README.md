# Losing sites to a classification made on under-counted depth

Worked on panel window `chr20:39,838,293-39,868,811` (30.5 kb; hiphase spans it
at 100.00% over 192 reads; our reads there are already 100% correct, so nothing
is corrupted and the question is purely why we do not link).

## We hold every site the competitor crosses on

| position | hiphase | us |
|---|---|---|
| 39,838,293 `A>G` | `0\|1` | `0\|1` |
| 39,846,791 `T>TATATATATATATATA` | `1\|2` (multiallelic) | **two records**, `1\|0` and `0\|1` |
| 39,848,373 `ACATTATATATAAT>A` | `1\|0` | `1\|0` |
| 39,848,887 `TTATATAATACACA>T` | `0\|1` | `0\|1` |
| 39,856,144 `C>T` | `0\|1` | `1\|0` |
| 39,860,491 `C>CACAC` | `0\|1` | `39,860,487 T>TACAC` (same event) |
| 39,868,811 `T>C` | `1\|0` | `0\|1` |

So discovery is not the problem. Our blocks are `39,788,782-39,848,887` (82
sites) and `39,856,144-39,913,010` (122 sites): the chain breaks between them
over **7.3 kb**, and the vote there is `agree=0 conflict=0` although **36 reads
span** `39,848,887 -> 39,856,144`.

## The bug: MSA sites were counted over the clustered reads only

The candidate table gives the break's two sites depths that cannot be right:

| position | type | reported DP | real coverage | ratio | initial class |
|---|---|---:|---:|---:|---|
| 39,845,582 | DEL | 34 | 56 | 0.61 | LOW_COV |
| 39,846,792 | INS | 18 | 68 | 0.26 | LOW_COV |
| 39,848,374 | DEL | 22 | 37 | 0.59 | LOW_COV |
| 39,848,888 | DEL | 19 | 32 | 0.59 | LOW_COV |
| **39,849,434** | INS | **61** | 63 | **0.97** | **CLEAN_HOM** |
| **39,856,144** | **SNP** | **16** | **65** | **0.25** | **LOW_COV** |
| 39,860,488 | INS | 21 | 59 | 0.36 | LOW_COV |

A clean heterozygous SNP at 65x is recorded at DP 16 with 5 ref / 11 alt, while
the one site in the neighbourhood that is not MSA-derived reads 0.97 of its true
coverage. Base quality is not the cause: 53 of the 70 reads have BQ >= 30 at that
base, `LOW_QUAL_COUNT` is 0, and the pileup is a clean ~50/50 split.

The cause is the rescue channel measured on the other panel window. A read the
MSA clustering does not place gets no observation at any site in the region, and
a site's counts are accumulated from the same set, so both the depth and the
per-read observations miss every unplaced read. The composition path that would
recover them was gated behind `--private-msa`.

**Fixed** by composing every unplaced read whenever the caller asks for them, and
by having the caller always ask (it previously asked only under `recover_gaps`).
The consumer, `add_msa_site_observations`, only adds an observation where there is
none and accepts an allele only where two independently composed paths agree, so
it can raise a site's depth toward its true coverage but cannot overturn an
existing observation.

| position | before | after | real |
|---|---|---|---:|
| 39,856,144 `C>T` | DP 16, 5/11, AF 0.688 | **DP 70, 34/36, AF 0.514** | 65 |
| 39,848,888 | DP 19, 12/7, AF 0.368 | **DP 63, 29/34, AF 0.540** | 32 |
| 39,848,374 | DP 22, ratio 0.59 | **DP 65, 37/28, AF 0.431** | 37 |
| 39,845,582 | DP 34, ratio 0.61 | **DP 66, 25/41, AF 0.621** | 56 |

The `before` column above reports **DP / covering reads**, the ratio the
discovery pass printed (34/56 and 22/37), not an allele fraction: no ref/alt
split was computed for these two records before the fix, so they have no
measured `AF` to compare against. An earlier version of this table labelled both
columns `AF`, which put two different quantities under one header. The `after`
column's values are allele fractions.

Panel result: **byte-identical to stock in every window** -- 0
concordant->discordant, 0 tags lost, 0 new tags, same blocks, same in-gap hets.
So this is a data-correctness fix with no behavioural change, and it is safe as a
default.

## Why it does not yet recover the site: the verdict is stale

`INIT_CAT` stays `LOW_COV` and the category stays `NOISY_CAND_HET` even at
DP 70 with 34/36. Classification runs in `collect_var_classify` ->
`classify_chunk_candidates` (`collect_var.cpp:2010`) **before**
`collect_noisy_vars_step4` (`collect_var.cpp:2113`) creates and refreshes the
MSA sites, and nothing reclassifies afterwards. So the site is judged at DP 16,
filed `LOW_COV`, and keeps that verdict once its counts are corrected -- which is
precisely how the site is lost: as a noisy candidate it cannot be a clean anchor,
and on stock defaults the noisy class never enters phasing at all
(`hybrid_collect.cpp:128` sets `skip_noisy_kmeans = true`).

The next change is to re-run classification on the corrected counts. It is not a
one-liner: `classify_variant_initial` (`collect_var.cpp:1462`) needs only chunk
data and recomputes the allele fraction itself, but it sets
`counts.candvarcate_initial`, whereas the mask phasing reads
(`lcd_var_i_to_cate`) is assigned by the later passes of
`classify_cand_vars_pgphase` (`collect_var.cpp:1518`). Re-running that whole
function after step 4 would also re-screen the MSA sites, which is the behaviour
we want for the repeat indels and the risk to measure: a 34/36 SNP at DP 70
should become `CLEAN_HET_SNP`, but several neighbouring indels sit at AF ~0.5 with
corrected depth in the 60s while segregating at chance, and promoting those would
be a regression. The panel plus the 0.559% chromosome-wide baseline are the gates.

Also visible here and not yet addressed: `39,846,791` is emitted as **two
records with identical alleles and opposite genotypes**, each at AF 1.000 with no
reference reads, where hiphase emits one multiallelic `1|2`. That is the same
representation defect found at `48,225,787`, in a second window.

## Reclassifying on the corrected counts: tried, destructive, reverted

The obvious follow-up is to re-run classification once step 4 has corrected the
counts. Two attempts, both instructive.

**Attempt 1 guarded on the candidate count and was inert.** `collect_noisy_vars_step4`
*inserts* its own candidates, so the vector grows and a `size()` comparison
skipped the whole block. Nothing changed, which is why `INIT_CAT` still read
`LOW_COV` at DP 70.

**Attempt 2 compared by variant key instead, so it ran -- and deleted the sites.**
After `classify_cand_vars_pgphase` was re-run, the entire neighbourhood
`39,845,000-39,862,000` was reduced to the single non-MSA site
(`39,849,434 CLEAN_HOM`): every MSA-derived candidate, including the
34/36-at-DP-70 SNP, disappeared from the output.

The reason is in that function's own comment: it is the initial-discovery pass,
and it assumes `cand.alt_ref_base` holds a reference base, while noisy candidates
built by `make_cand_vars_from_baln0` carry consensus column bytes there. Re-run
over MSA candidates it therefore misreads them, assigns a category that
`prune_not_candidate_variants` deletes, and the site is lost outright -- a worse
outcome than the stale verdict it was meant to fix.

Reverted; the panel is byte-identical to the baseline again and the suite passes.

So reclassification cannot be done by re-running that function. The narrower
route is to call `classify_variant_initial` alone on the corrected counts -- it
takes only the key, the counts and the reference slice, and recomputes the allele
fraction itself -- and then set the phasing mask directly for the clean verdicts,
rather than re-deriving every candidate's mask. That needs the
category-to-`lcd_var_i_to_cate` mapping the later passes apply, which is the part
to read next.

Still open in this window, and unchanged by any of the above:
`39,846,791` is emitted as two records with identical-length but different
insertion sequences and opposite genotypes, each at AF 1.000 with **zero
reference reads** (DP 27 and 23 after the depth fix, 0 ref in both). Each record
is internally consistent over its own read subset; together they are one
multiallelic locus that hiphase emits as a single `1|2`. Emitting it as two
biallelic hets is what lets both claim AF 1.000.
