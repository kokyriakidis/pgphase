# Defer a clean-SNP stitch past a recovered insertion

The existing HiPhase same-callset run phases the clean chr20 SNPs at
54,483,506 T>G and 54,490,229 G>C in one phase set and places all 93
truth-scored reads overlapping their 6,723-bp gap correctly. The accepted
pgphase graph/BAM output split the SNPs. Thirty-two independent MAPQ-30
primary reads call both SNPs at base quality 30 or higher. Thirteen carry
left ALT/right REF and 19 carry left REF/right ALT; none contradict this
relative phase. Both allele classes are represented.

The prior exact-boundary audit missed this SNP pair because the right
pgphase block begins at an insertion, not at its first clean SNP.

The original graph recovery seam ends at a recovered insertion at
54,489,927. The physical stitch tries that insertion, then limits its
clean-SNP search to the old seam endpoint. It never considers the clean
right SNP 302 bp later. The fix lets the existing physical-SNP likelihood
and block-path checks use the first clean SNP in the same right phase set
when no other phased block intervenes. It does not reinterpret the repeat
allele at the boundary.

Trying this extra pair immediately joined the intended SNPs, but changed
the right block before the already supported 54,547,514–54,569,072 seam
was processed. The resulting full chr20 VCF lost five rows and the BAM lost
one truth-correct tagged read. That ordering was rejected. Deferring the
extra SNP pair until ordinary seams finish preserves the downstream join
and its path certification. In the owning 54–55 Mb replay, the new and
existing joins share one phase set with the truth-derived ALT relationships:
left versus right opposite; right versus 54,547,514 same; 54,547,514 versus
54,569,072 opposite. The 93 overlapping reads are all in one truth-correct
block instead of 67 in the dominant block. Owning-chunk read accuracy stays
4,306/4,315, while read phase sets fall seven to six.

| Full chr20 measure | Before | Deferred SNP stitch |
|---|---:|---:|
| VCF variant keys | 62,361 | 62,361 |
| Truth-scored tagged reads | 236,866 | 236,866 |
| Truth-correct / discordant reads | 229,131 / 7,735 | 229,131 / 7,735 |
| Read phase sets | 695 | 694 |
| VCF phase sets | 345 | 344 |

A final guard also stops at an intervening phased indel, using the graph
allele's minimal reference coordinate. Its full chr20 VCF, candidate TSV,
and tagged BAM are byte-identical to the deferred trial.

Exactly 12 VCF sample fields change, all in the 54 Mb chunk. The previously
separate 160/160-correct and 875/880-correct read blocks combine into one
1,035/1,040-correct block, exactly the sum of their correct and discordant
counts. The panel has a new exact-SNP span and parental-orientation check;
the existing downstream gap remains a span. Runtime uses no parental truth.

Reproduction inputs are the shared chr20 fixture. Before output:
`/tmp/pgphase-orphan-full.vcf` and `.bam`. Final output:
`/tmp/pgphase-next-snp-final-full/phased.vcf` and `.bam`. The rejected
immediate trial is `/tmp/pgphase-next-snp-full/`.

Validation on the final binary: `make -j8 pgphase`, `make unit-tests`,
`make check`, and `make window-tests` pass. The clean-build window suite
reports 3,632 assertions in 45 test cases; the focused new window passes
33 assertions. `git diff --check` is clean.
