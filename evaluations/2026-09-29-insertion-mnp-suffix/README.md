# BAM insertion to graph MNP gap, 2026-09-29

The tracked chr20:17,865,146–17,883,198 gap has an MSA-verified BAM insertion
(`T→TTGTG`) on the left and a graph MNP (`AC→GT`) on the right. Five reference
bases separate the insertion's VCF position from the equivalent CIGAR placement
used by some spanning reads. The old left-to-right stitch demanded an insertion
at one exact CIGAR coordinate, and its physical boundary picker accepted only
single-base SNPs. The seam detector also placed an equal-length MNP one base
before its first changed base, so the targeted interval ended before the right
boundary.

The repaired bridge verifies both changed MNP bases and admits a shifted
insertion only when the two edits produce the same local sequence with clean
flanks. At MAPQ/base quality at least 30, the outer insertion-to-MNP edge has
three REF and one ALT insertion calls, with quality-weighted log odds -26.8793
for the inherited phase. The left BAM source has a weak cut immediately before
the insertion. Moving that insertion alone closes the VCF gap, but mixes
oppositely oriented rescued reads in one fallback phase set: its owning-chunk
truth falls from 4,049/4,082 to 4,058/4,106 and 20 previously correct reads
become discordant.

The retained transfer uses the nearest clean BAM SNP at 17,852,024 before the
last weak cut. Direct high-quality molecules call that SNP and the insertion
on both alleles (five reference-insertion and two alternate-insertion calls)
and certify their inherited same-haplotype relation. The earlier source prefix
remains independent; the SNP, both complementary deletion rows at 17,856,666,
and the insertion transfer together to the right MNP's phase set. Crossing
read tags are cleared. In the owning 17–18 Mb replay the gap closes and truth
improves from 4,049/4,082 to 4,052/4,082 correct reads, with no mixed fallback
phase set. The dedicated regression requires the exact boundary rows to share
one PS, keeps the upstream BAM SNP separate, and checks local read truth.

An exploratory change that included every graph MNP in the clean-SNP path
validator reopened two previously closed tracked gaps at 1.508 and 3.529 Mb.
The retained code leaves that path rule unchanged while allowing an MNP as a
physical boundary for the targeted allele bridge. A broader lookup for
BAM-only SNPs in all physical seams was also rejected after it showed the same
regression in a full-chromosome trial.

| Full chr20 | Before | After |
|---|---:|---:|
| Tracked HiPhase-correct gaps closed | 17/23 | 18/23 |
| VCF keys | 62,154 | 62,154 |
| Truth-scored phased reads | 236,847 | 236,854 |
| Correct / discordant | 229,088 / 7,759 | 229,098 / 7,756 |
| Truth-scored read phase sets | 708 | 704 |
| Truth accuracy | 96.7240% | 96.7254% |

All VCF genotype strings and dosages remain unchanged. Four common reads become
truth-correct, one becomes discordant, seven previously unphased reads gain
truth-scored tags, and none lose tags. The full comparison uses
`/tmp/pgphase-gap1383-full/` and `/tmp/pgphase-gap17865-final-full/`.

The graph panel's 17.865 Mb span changes from open to closed, raising its
expected total from 46 to 47. The owning-chunk test also checks the separate
upstream PS and local truth. Final validation: `make -j4`, `make unit-tests`,
`make check`, and `make window-tests` pass; the window suite reports 1,848
assertions in 38 test cases. `git diff --check` is clean.
