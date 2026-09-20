# The co-located allele pairs: a dead rule, and what fixing it bought

Of the 38 positions where our alleles differ from upstream's, 31 had upstream
emitting two rows where we emitted one. At `chr20:3,863,176` upstream writes
`C>CAAAAAAAAA 1|0` beside `C>CAAAAAAAAAAA 0|1`; we wrote `C>CAAAAAAAAA 1|1` and
nothing for the longer allele.

## We held both alleles all along

From the candidate table at that locus:

| allele | DP | ref/alt | HAP_ALT | outcome |
|---|---:|---:|---:|---|
| `C>AAAAAAAAA` (9 A) | 75 | 29/46 | 3 (both haplotypes) | emitted `1\|1` |
| `C>AAAAAAAAAAA` (11 A) | 75 | 51/24 | 0 (neither) | dropped, carries no ALT |

So this was never a representation or discovery problem. The consensus gave one
allele both haplotypes and the other none.

## The rule meant to fix this was dead on the alignment arm

`make_colocated_alleles_complementary` exists for exactly this, and two things
stopped it:

1. Its depth test read `c.counts.alle_covs[1]`, and a BAM-discovered candidate
   leaves `alle_covs` EMPTY, carrying its depths in `ref_cov`/`alt_cov` -- 662
   such candidates per megabase. So `alt_depth` returned 0 and the
   `min_alt_support` test rejected every pair. Now falls back to `alt_cov`.
2. It only handled the case where both candidates claim allele 1 on the SAME
   haplotype. The case here -- one claims BOTH, the other NEITHER -- fell
   through the `continue` labelled "already complementary", which it was not.

The new branch gives each allele one haplotype, choosing the orientation from
the dropped allele's own per-haplotype profile and firing only when that profile
is present and asymmetric, so the assignment comes from read evidence rather
than from candidate order.

## What it bought, and what it cost

Whole chr20, alignment arm:

| | before | after |
|---|---:|---:|
| records | 118,273 | 118,290 |
| identical to upstream | 117,900 | **117,911** |
| upstream-only | 370 | **359** |
| ours-only | 373 | 379 |
| allele-difference positions | 38 | **33** |
| GIAB TP / FP | 91,479 / 3,415 | 91,479 / **3,431** |
| GIAB F1 | 0.9632 | 0.9631 |
| misplaced reads | 3,226 | 3,226 |

Of the 17 newly emitted records, 15 fall at positions where the truth set has a
variant with a different allele representation, 1 is a false positive and 1 is
outside the benchmark. Under exact-allele matching those 15 score as false
positives, which is where the +16 FP comes from; a representation-normalising
comparison would likely count most of them as matches. No true positive was
lost and read placement is unchanged.

## The 33 that remain

| why | count |
|---|---:|
| both/none pair the new branch rejected (profile absent or symmetric) | 16 |
| two candidates, neither claiming both nor neither -- a different shape | 15 |
| we hold only one candidate where upstream has two alleles (DP 6 SNPs) | 2 |

The 16 are the tractable next step and need one probe: why the per-haplotype
profile is empty or symmetric at those loci when the rule runs on the final
table. The 2 are a genuine discovery difference at depth 6, not a
representation one.

## A correction made along the way

I first set `collapse_colocated_alleles = false` on the ported path, reasoning
that collapsing was our invention and upstream splits. That was backwards:
`make_colocated_alleles_complementary` is what PRODUCES upstream's
one-row-per-haplotype shape, so disabling it moves away from parity. Measured
inert either way in that configuration, and left enabled. The option and its
call-site gates remain, documented.
