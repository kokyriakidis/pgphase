# chr20 55.9 Mb: correct the bridge rather than reject its read noise

The shared SNPs at 55999561 and 56040612 are both 0|1 in truth. HiPhase
phases both as 0|1 in PS 55815775. WhatsHap, optimized WhatsHap, and LongPhase
have separate PS values at these endpoints; only HiPhase bridges this case.

The original gap-specific margin trial joined two blocks in the wrong
orientation. All 39 original right-block reads became discordant and all 269
left-block reads stayed concordant; no new reads were tagged. The shared SNPs
became 1|0 and 0|1 in the same PS, independently confirming the wrong stitch.

## Root causes

MSA found a four-base and a six-base deletion at 56007501, and one-base and
two-base deletions at 56027379. The latter is a genuine multiallelic genotype:
truth and DeepVariant both contain 56027378 TAA→T,TA with GT 1/2 (orientation
aside). The former is also 4/6 in truth, while DeepVariant calls only the
four-base event. The old MSA representation made each deletion an independent
het against reference. The site-level flank rule allowed a longer deletion
to count as ALT for the shorter event as well. All five observed spanning
reads had ALT=1 at both four-base and six-base rows. Four BAM alignments have
a six-base deletion shifted right within the tandem repeat; the fifth is noisy.

Splitting the shared deletion from the length difference makes the common
four-base event homozygous and places a heterozygous two-base difference at
56007505. The right site similarly becomes a shared deletion plus a one-base
difference at 56027380. The two biological repeat lengths now yield opposite
allele observations, covered by a failing-before/passing-after regression test.

A second problem remained: provisional HP assignments can give both haplotype
profiles the same majority allele at a balanced, MSA-verified repeat. The core
then excludes the site from its connectivity graph before its edges can resolve
the arbitrary orientation of the neighboring blocks. During the gap-local HP
fallback, verified sites passing existing depth/AF checks now retain a
heterozygous genotype for that solve. Ordinary clean phasing is unchanged by
this genotype preservation. Repeat edges use the existing net-margin threshold;
read-level disagreement is allowed.

## Measured result

`chr20_nested_seed` and `chr20_observed_owner`: two blocks→one, 308 assessed
reads, zero discordance and zero switch/flip errors. Shared SNPs are both 1|0
in PS 55949998, equivalent to truth after one global HP-label flip. The 54.5 Mb,
57.1 Mb, 62.4 Mb and 46.4 Mb positive controls remain joined correctly locally.

The full 114-case genotype-preservation trial also introduced bad joins in
other regions; it is not a chromosome-wide safety result. The subsequent
missing-observation ownership fix is evaluated separately. Frozen binaries,
commands and matrices are under /tmp/pgphase-chr20_nested_seed/ and
/tmp/pgphase-5599-root-cause/, with the representation-only diagnostic matrix
under /tmp/pgphase-5599-root-fixed/.

## Direct candidate-site validation (current local result)

The final direct-site check distinguishes the true residual deletion from the
competing noisy site using the same clean SNP at 56040612:

| MSA site (one-based event position) | Clean allele 0: ref / alt | Clean allele 1: ref / alt | Eligible |
| --- | --- | --- | --- |
| 56027380 | 9 / 2 | 1 / 9 | yes (both margins >= 2) |
| 56029103 | 10 / 1 | 5 / 6 | no (second margin only 1) |

The true site's leftward bridge has 3 supporting and 1 conflicting observations;
its rightward bridge has 18 supporting and 3 conflicting observations. Individual
read conflicts are retained, while the competing poorly separating site is not
admitted as a repeat link. The resulting chr20_direct_site_anchor local run joins
the target with all 308 truth-assessed reads concordant and no reversed original
block. This requires no competitor or truth labels in the algorithm. Broad panel
validation is still pending.
