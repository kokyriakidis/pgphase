# Complementary deletion suffix stitch (2026-09-30)

At chr20:64,138,752–64,140,314, the original seam list contains the
64,140,314–64,144,256 seam. Moving the 64,140,314 BAM deletion into the right
block exposes the earlier seam, but the list was a snapshot and the physical
stitch did not revisit it. The left boundary deletion alone has 11 matching
and one conflicting BAM matrix allele pair to the right deletion. An earlier
complementary pair at 64,134,226 supplies 11 unambiguous paired reads with
both left ALT classes and no conflicts (five maternal, six paternal by the
truth-only audit). Its original BAM source has no weak or quality cut between
the pair and the boundary. The graph SNP at 64,118,182 has only one paired
observation to the BAM run, so the graph prefix must remain separate.

The stitch now revisits a seam exposed by a boundary move and uses only
MSA-verified BAM rows and observed BAM-channel alleles. Its quality guard reads
the BAM channel's MAPQ, not the GAF read's MAPQ; otherwise valid BAM bridges
are discarded when the graph alignment has lower confidence. It requires at least
two distinct reads per left allele class, unanimous phase relation, an exact
one-sided binomial tail at most 0.001, exactly two BAM rows at the left locus,
and no source cut inside the transferred run. It transfers the suffix beginning at
the complementary pair while retaining the earlier rows. One trial passed
the graph candidate's padded key position as a physical split coordinate;
that incorrectly absorbed the graph prefix and gave 86.79% local read
concordance. Using the selected graph allele's physical position fixed the
coordinate error. The 64 Mb owning-chunk regression then closed the tracked
gap and passed 37 assertions, including parental orientation, the retained
earlier BAM SNP at 64,128,828, and 99% read concordance for both adjacent
windows.

The gap was already in the committed panel. Its graph-arm `spans` expectation
is updated from zero to one, and the dedicated owning-chunk test now checks
the new deletion connection and unchanged graph-prefix split. The test replay
cache introduced in the same work session shares identical owning-chunk runs
within one Catch2 process; it does not change pipeline output.

## Full chr20 comparison

The final pair-only implementation was run over the entire chr20 input in
`/tmp/pgphase-comp-left-bam-mapq-final/` and compared with the immediately preceding
`/tmp/pgphase-gap61757-ordered-full/` baseline. Both emit the same 62,361 VCF keys and the same genotypes. Only three VCF
sample fields change, all PS labels on the two complementary deletion rows
and the following boundary deletion. Both runs tag the same 236,880
parental-truth reads. Both score 229,144 correct and
7,736 discordant reads. Twenty reads change PS label and three change HP
label, but no read changes parental correctness. Read phase-set count rises
692 to 693 because the unsupported prefix remains a separate block. The
64,138,752 and 64,140,314 deletion records now share PS 64,905,482 in the
full VCF; the earlier graph SNP and BAM SNP remain in PS 63,896,021.

The one-sided significance and BAM-MAPQ guards were added after the initial
full run. The final VCF and BAM are byte-identical to that earlier trial;
the guards do not change this chr20 fixture's emitted calls or read labels.

## Certified read-label transfer

The generic suffix splitter clears reads with calls on both sides of the cut.
Six of those reads also call the complementary left deletion and the right
deletion, so their right-block haplotype is already established by the same
11-read stitch certificate. Retaining only those six reads in the right PS
raises the target window's correctly separated reads from 32/53 (60.38%) to
38/53 (71.70%); 45/53 truth-scorable reads remain tagged. The full chr20
trial at `/tmp/pgphase-comp-left-read-trial/` has a VCF identical to the
pair-only baseline and changes six read PS labels, no HP labels. All 236,880
truth-scored tagged reads retain their correctness status: 229,144 correct,
7,736 discordant, and 693 read phase sets. The 0.71 separated floor is now
asserted in both the window expectation and the owning-chunk regression.

Validation: `make -j20`, `make unit-tests`, and `make check` passed. The full
gap suite on the final runtime passed 3,731 assertions across 45 cases. The
strengthened owning-chunk case then passed 38 assertions, and the isolated
one-window panel with the 0.71 floor and required-site row passed 20.
