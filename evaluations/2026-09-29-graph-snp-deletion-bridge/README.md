# Graph SNP to complementary BAM deletion bridge, 2026-09-29

The tracked chr20:59,825,454–59,842,960 gap separated a graph SNP block from
an MSA-verified BAM block with two overlapping, opposite-haplotype deletions:
`GTGTT→G` and `GTGTTTGTT→G`. The three primary BAM reads spanning the two
boundaries have left SNP base qualities 10, 22, and 40. A physical SNP-pair
bridge therefore had no qualifying paired calls, although all three reads
already had a graph haplotype and an observed graph allele at the left SNP.

Two spanning maternal reads support the short deletion. Their CIGARs place its
four-base deletion at 59,842,961 and 59,842,969. The two edits yield the same
reference sequence in the local repeat, with clean aligned bases through the
verified span. The older deletion caller rejected the exact-placement read
because a separate one-base insertion at 59,842,939 fell within its broad
search window, although it was outside the allele's verified span. A third
spanning read has a different three-base deletion and abstains.

The repaired caller checks whether other indels intersect the bases needed to
verify a sequence-equivalent allele. The new stitch uses each read's established
graph HP and graph SNP observation on the left, and a clean, physically
verified ALT call for exactly one of the two deletion rows on the right. It
requires two independent MAPQ-30 reads with one consistent orientation, a
wrong-parity bound of 0.001, and supported graph SNP paths on both blocks.
Overlapping deletion REF calls remain ambiguous and do not vote. The BAM
candidate rows retain their original representation and coordinates.

The 59–60 Mb owning-chunk replay moves the left SNP and both right deletion
rows into PS 59,778,440. Its truth-scored read result is unchanged at
3,654/3,668 correct (14 discordant); the separate before/after phase-set
orientation votes agree on the same parent. The short 117.5 kb panel replay
also closes. The exact-row owning-chunk regression checks their phase and
opposite deletion haplotypes and rejects a parental switch.

| Full chr20 | Before | After |
|---|---:|---:|
| Tracked HiPhase-correct gaps closed | 18/23 | 19/23 |
| VCF keys | 62,154 | 62,154 |
| Truth-scored phased reads | 236,854 | 236,854 |
| Correct / discordant | 229,098 / 7,756 | 229,098 / 7,756 |
| Truth-scored read phase sets | 704 | 702 |
| Truth accuracy | 96.7254% | 96.7254% |

All VCF genotype strings and dosages are unchanged. No truth-scored read gains
or loses a tag or changes correctness. Full comparison:
`/tmp/pgphase-gap17865-final-full/` versus `/tmp/pgphase-gap59825-full/`.
The graph panel's 59.825 Mb span and total expected count move from 0 to 1
and from 47 to 48, respectively.

Final validation: `make -j4`, `make unit-tests`, `make check`, and
`make window-tests` pass. The window suite reports 2,853 assertions in 39
test cases; `git diff --check` is clean.
