# A physically certified graph path closes the 47 Mb gap

Target: chr20:47,003,897–47,713,869 in the 23-target HiPhase-correct panel.
This is a 710 kb gap in phase-set labels, not a 710 kb stretch lacking sites:
the large graph/BAM block between its endpoints is already phased. The
six-site MSA-verified BAM island at its left boundary retained a separate
phase-set label but no read tags. The right SNP block began at 47,713,869.

## Diagnosis and retained checks

* The source island has a complete BAM path with no nonempty weak or quality
  cut. A clean graph SNP at 47,002,481 and the island insertion have direct
  physical same-haplotype evidence. The graph SNP's internal allele ID is a
  graph walk, so the physical caller now uses its normalized VCF REF/ALT base.
* A weak graph SNP edge at 47,671,540–47,689,418 has two independent MAPQ-60
  reads calling both SNPs at Q35/40. Their relation agrees with the existing
  block. This independently certifies the edge before the source island or
  right block can inherit its orientation.
* The MSA insertion at VCF 47,694,119 has equivalent CIGAR placements nearby.
  One bridge read has an inserted base at Q10 but clean aligned flanks and
  right SNP. The final bridge allows that base and uses its actual Q10 error
  probability in the parity likelihood. It observes both insertion alleles on
  two reads each, with decisive orientation at a wrong-parity bound of 0.001.
  A direct clean SNP-to-insertion check certifies its link to the graph suffix.

The earlier readless block is attached as candidates only; no read tags are
moved. The graph block and right SNP then join as a whole after both internal
links and the outer insertion-to-SNP bridge pass. The ordinary high-quality
insertion caller remains unchanged for other call sites.

## Measurements

| Metric | Before | After |
|---|---:|---:|
| Tracked gaps closed | 19/23 | 20/23 |
| VCF `(CHROM,POS,REF,ALT)` keys | 62,154 | 62,154 |
| Truth-scored phased reads | 236,854 | 236,854 |
| Correct reads | 229,098 | 229,097 |
| Discordant reads | 7,756 | 7,757 |
| Read phase sets | 702 | 700 |

The full runs are `/tmp/pgphase-gap59825-full/` and
`/tmp/pgphase-gap47003-final-full/`. Accuracy changes 96.72541% ->
96.72499%. Exactly 11 genotype phase strings change, only in the merged 47 Mb
blocks; allele dosages stay the same. The local owning-chunk replay has 3,993
truth-scored phased reads and no parental switch; one read changes from
correct to discordant. Its exact island, suffix insertion, and right SNP rows
share a phase set, with opposite phase genotypes on the boundary alleles.

The three tracked targets still open are 21,377,985–21,395,286;
23,421,003–23,445,252; and 32,234,664–32,246,127. Their available
boundary evidence is weaker than this bridge; they were left split.

The panel now replays the owning 47–48 Mb chunk for this target and asserts
its exact boundary join and local truth orientation. The graph total span
expectation rises 48 -> 49.

Validation: `make -j4`, `make unit-tests`, `make check`, and
`make window-tests` pass. The full panel has 2,870 assertions in 40 cases.
The prior 47–48 Mb graph-gauge check now observes 2,967/2,998 reads on the
majority parental orientation (98.966%); its floor is 98.9%, with exact
site-gauge and opposite-genotype checks retained.
