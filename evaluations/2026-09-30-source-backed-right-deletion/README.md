# Source-backed right-deletion seam (2026-09-30)

The graph/BAM run left chr20:8,166,027–8,172,072 split after recovery. Its
right boundary is an MSA-verified BAM deletion. The final physical stitch
retried newly exposed MNP and long-insertion boundaries, but skipped this
right deletion. The nearest callable clean SNPs lie outside the narrow seam:
8,165,210 on the left and 8,176,552 on the right. Five distinct MAPQ/base-Q
30 BAM reads call both SNPs; four carry left ALT/right REF and one carries
left REF/right ALT. All five imply the same phase orientation, with a
quality-weighted log odds of 41.08.

The left block mixes catalog and BAM sites. The BAM source has a complete
path without weak or quality cuts, but its 8,165,210 SNP was unphased in the
original graph and adopted from the source after recovery. Requiring that SNP
to have had an original graph PS rejected the exact source allele. A direct
read edge to the preceding clean graph SNP confirms the adopted orientation.
The older graph edge at 8,132,545–8,151,005 has no direct read pair in the
recovery matrix, but both endpoints are exact shared SNPs in the same complete
BAM source. The stitch now accepts that source certificate for an absent or sparse
agreeing edge; any directly reversed pair vetoes the whole-block join. The right graph
path remains an independent requirement.

The owning 8–9 Mb run now puts both complementary insertion alleles at
8,166,027 and the deletion at 8,172,072 in PS 8,106,315. The two flanking
clean SNP ALT alleles at 8,165,210 and 8,176,552 occupy opposite haplotypes.
The new panel case spans, has two phased heterozygotes strictly inside its
bounds, at least 99% local read concordance and 81% correctly separated reads.
The dedicated owning-chunk test checks all five exact VCF alleles and parental
read orientation; it passed 29 assertions. Same-callset
HiPhase leaves both left insertion alleles unphased, so this case is labeled
`none_spanning` in the panel.

Full chr20 at `/tmp/pgphase-8166-final-full/` retains all 62,361 variant keys and
unordered genotypes and all 236,880 truth-scored tagged reads. The previous
run at `/tmp/pgphase-comp-left-read-trial/` scored 229,144 correct and 7,736
discordant reads in 693 read phase sets; the new run scores 229,143 correct and
7,737 discordant in 691 read phase sets. The two principal read blocks have 299/299 and 2,513/2,513 parentally
correct reads before the join; the merged block has 2,812/2,812 after it.
Exactly one truth-scored maternal read in an auxiliary read PS changes
correctness when its right-block PS is absorbed; no other read changes
correctness. The 650 changed VCF sample fields are the expected
right-block PS rename and haplotype orientation flip over 8,172,072–8,748,147;
no variant genotype changes.

The final direct-reversal veto was rerun over chr20; its VCF and phased BAM are
byte-identical to the preceding trial, and the owning-chunk regression still
passes.

Validation on the final build: `make -j8 pgphase`, `make unit-tests`,
`make check`, and `make window-tests` passed. The complete gap suite passed
3,788 assertions across 46 test cases; `git diff --check` passed.
