# Blocks phased at chance, and 27% of the reported error that is not error

Started from "what made the scoped escape improve things", which turned out to
be the wrong question: it did not improve anything.

## 1. The scoped escape is churn, not a gain

Whole chr20, HEAD default vs the scoped escape (752f06a): 510 genotype changes,
786 phase-set changes, 388 records lost and 367 gained. Per read, against the
parental truth map:

| | reads |
|---|---:|
| wrong -> right | 87 |
| right -> WRONG | 72 |
| **net** | **+15 of 219,055 (0.007%)** |

All of it inside 7 blocks in one region, and the same blocks appear on both
sides. The commit's "1.161% -> 1.153%" is arithmetically true and within churn.
It should not be read as an improvement the change earned.

## 2. 15% of the reported error is unscorable, not wrong

18 blocks holding 949 reads are **more than 90% one parent**. With one parent's
reads present, any haplotype split scores about 50% by construction, so these
blocks cannot be scored by read hamming at all. They carry **371 of 2,526
misplaced reads**.

chr20 read hamming, scoped-escape arm:

| measurement | value |
|---|---:|
| as reported all session | 1.153% |
| excluding single-parent blocks | **0.988%** |

`PS=26896221` was the worst-looking block on the chromosome at 48.6% minority.
Its reads are 68 maternal against 4 paternal. It is not a defect.

Arm-to-arm comparisons in earlier records stand, because the artifact is common
to every arm, but the absolute figures are inflated.

## 3. A real defect: blocks emitted with no separating evidence

Of 301 blocks with at least 20 scored reads, **13 have no site that separates
the parents** (sampled up to 12 records each, purity >= 0.90 against read
truth). They are not single-parent blocks -- both parents are present:

| PS | reads | records | informative | minority |
|---|---:|---:|---:|---:|
| 45,391,776 | 132 | 12 | 0 | 47.0% |
| 44,300,910 | 74 | 6 | 0 | 47.3% |
| 45,859,664 | 74 | 5 | 0 | 45.9% |
| 44,779,026 | 70 | 3 | 0 | 48.6% |
| 45,330,583 | 70 | 10 | 0 | 45.7% |

Mean minority in this class is **34.2%**, against **4.77%** where at least one
informative site exists. They tag 854 reads and carry **315 of 2,526 misplaced
reads (12%)**. `45,391,776` is about 50/50 correct in every 20 kb positional
bucket, so it is not a switch -- it is a coin flip over 132 reads.

Together with section 2, roughly **27% of the reported read error is metric
artifact or coin-flip tagging rather than mis-phasing**.

## 4. What those regions are, ruled in and out

For `chr20:45,391,776-45,431,550`:

- **Not a competitor advantage.** Hiphase spans it with one block and scores
  **49.0% minority over 100 reads** -- it fails identically. Its 3 phased hets
  against our 12: its only scorable site separates at 0.552, the same event and
  the same numbers as our 45,392,720.
- **Not a discovery gap.** The one site hiphase has that we lack,
  `45,410,343 TCCTGTTGTGGACC>T`, IS discovered by the alignment arm
  (`45,410,345 CTGTTGTGGACCAG>.`) and called homozygous.
- **Not a MAPQ discard.** All 53 reads overlapping that locus are MAPQ 60 and
  `samtools depth` is identical at every floor from 0 to 60. The depth of 5 is
  real: **53 of 53 reads carry a >= 10 bp deletion spanning the locus**, so
  almost no read has an aligned base at the anchor. Net length over 51 fully
  covering reads is **-19 on 41 of them, 21 maternal and 20 paternal** -- the
  deletion is on BOTH haplotypes. Our 1|1 is right and hiphase's 1|0 is wrong.
- **Not a truth-map failure.** Truth covers 100% of reads there, and mean
  parental purity of balanced clean het SNPs is 1.000 in essentially every 2 Mb
  bin of the chromosome.
- **Not a collapsed duplication.** Median record depth in these blocks is 61
  against 64 in well-phased blocks.

What is left is the simple explanation: these are **runs of homozygosity**. In
`45,405,000-45,415,000` the alignment arm emits about 36 records of which only
4 are heterozygous; the rest are 1|1. There is nothing to phase, and the few
"hets" present do not segregate.

## The actionable bug

The pipeline still opens a phase set there and tags every read in it. A block
whose sites carry no read-separating evidence should not tag reads: the split
is arbitrary, it is indistinguishable from phasing in the output, and it costs
315 misplaced reads chromosome-wide. The truth-free form of the test is whether
any site in the block partitions the assigned reads by allele -- which is
computable from the phased BAM alone.

Not yet implemented; no code change in this record.
