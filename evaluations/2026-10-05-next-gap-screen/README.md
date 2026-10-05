# Unresolved gap screen after the 41.881 Mb closure

Baseline production binary SHA256:
`69d998cab0b2319f6875e3514a7be4eaa1dc9b514419831614aeb012a13c4386`.
This investigation closes no new gap. All experimental production edits were
reverted; panel expectations and production behavior are unchanged.

## Physical evidence

At 11,235,277–11,255,369, three primary MAPQ-60 reads span both boundaries.
All have a four-base deletion near the left boundary. Two have a two-base
deletion near the right boundary; the third has a one-base deletion. Therefore
the original alignments do not give a unanimous diploid boundary relation.
`inspect_deletion_spanners.py` reproduces these events without using read HP.

At 7,264,321–7,280,346, the blocking graph SNP pair is 7,280,346 C>T and
7,280,356 C>T. Its graph path votes are 3 same / 3 cross, all same votes on
one haplotype. Original MAPQ-30-or-higher primary alignments give C/C=37,
C/T=25 and T/C=3; requiring Q30 at both bases gives 23, 20 and 2 respectively.
These physical calls do not independently certify the graph block's existing
opposite SNP orientations. This is an observation/representation discrepancy,
not evidence for dropping the statistical reversal veto.
`inspect_snp_pair.py` reproduces the counts.

## Rejected trials

1. Report statistically weak GAF reversals as candidate cuts for the existing
   independent physical SNP check. A second trial applies that check to the
   right graph path, requiring both haplotypes and the existing 0.001 count
   threshold. Both leave the 7 Mb output unchanged (3,926 correct / 165
   discordant reads, eight VCF blocks); neither closes the nomination.
2. Permit a complete BAM source to anchor its first missing graph edge at a
   clean boundary, retaining the exact-site, orientation and cut checks. The
   actual blocking edge is not missing: it has reversal observations. This
   does not address the diagnosed evidence defect and was reverted.
3. Let the full-block physical validator inspect ordinary multi-read geometry
   as well as singletons. Accept multi-read pairs only when all physical
   relations agree, both haplotypes occur, and the two-sided random-parity
   probability is at most 0.001. All existing full-block gauge and source-path
   checks remain in place. Six fresh owning-chunk comparisons are unchanged:

| Owning Mb | Correct reads | Discordant reads | Changed tags / VCF rows | New nominated closures |
| --- | ---: | ---: | ---: | ---: |
| 7 | 3,926 | 165 | 0 / 0 | 0 |
| 15 | 4,272 | 16 | 0 / 0 | 0 |
| 21 | 3,111 | 210 | 0 / 0 | 0 |
| 34 | 3,539 | 28 | 0 / 0 | 0 |
| 48 | 3,844 | 57 | 0 / 0 | 0 |
| 56 | 3,914 | 74 | 0 / 0 | 0 |

The committed parity JSON files contain exact read-tag, VCF-row, key and
parental-truth comparisons. Fresh replay inputs use the standard chr20 BAM,
catalog, coordinate GAF and reference, four threads, and owning regions
`CHM13#0#chr20:{Mb*1000000+1}-{(Mb+1)*1000000}`. Native dumps and build logs
are under `test_data/tmp_gap_fix53/`. No full-chromosome improvement or new
regression expectation is claimed.
