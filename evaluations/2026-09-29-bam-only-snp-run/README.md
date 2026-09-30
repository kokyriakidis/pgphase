# Certify BAM-only source runs and prioritize clean SNPs

Two adjacent noncentromeric chr20 gaps, 19,395,544–19,403,172 and
19,403,172–19,414,720, are one correctly oriented HiPhase block. HiPhase's
dominant read block places 91/95 and 110/113 truth-scorable local reads.
The previous pgphase run placed 59/95 and 67/113, respectively.

At the right gap, 19 independent MAPQ/base-quality >=30 primary reads
call both boundary SNPs with decisive phase parity. Pgphase rejected their
edge because the left phase set held no graph SNPs: it consists of two BAM
insertion rows at 19,397,607 and a BAM SNP at 19,403,172. All three belong
to BAM source phase set 19,272,142. That source has weak cuts at 19,365,840,
19,373,922, and 19,377,345, before the local three-row run. The new path
check certifies only that run, requiring one source, consistent allele
orientation, and no internal weak cut. It joins the right gap without
absorbing the earlier source prefix.

The join exposes the upstream clean SNP at 19,395,544. Twenty-six Q30
primary reads call both that SNP and the recovered 19,403,172 SNP; all
26 support the same VCF ALT relation. An existing SNP-to-repeat-insertion
stitch instead gave a conflicting length vote (9 REF-like, 14 ALT-like;
log odds -17.72) and attached the block with the wrong haplotype connection.
The clean SNP pair passes the existing split-half profile test at p<=0.001.
It now takes priority only when its recovered SNP shares the insertion's
BAM source phase set and no source weak cut intervenes. All three emitted
SNP ALTs then lie on the same haplotype, matching HiPhase.

The owning 19–20 Mb replay improves dominant-block read separation to
79/95 at the upstream gap and 98/113 at the downstream gap. Two exact-site
window regressions pin the shared PS and haplotype orientation.

A broad retry initially changed chr20:36,620,864 G>A from singleton PS
36,620,864 into PS 36,549,255. Its local truth labels were mixed (31/54
correct in its original block), HiPhase leaves that genotype unphased, and
the trial lost seven truth-scored phased reads. A retry now requires a
callable clean SNP at the newly adjacent left boundary; the 36.62 Mb
indel-to-singleton seam stays split. An owning 36–37 Mb control asserts
that nonspan.

| Full chr20 | Accepted baseline | Guarded fix |
|---|---:|---:|
| Variant keys | 62,154 | 62,154 |
| Truth-scored phased reads | 236,859 | 236,859 |
| Truth-correct / discordant | 229,108 / 7,751 | 229,108 / 7,751 |
| Read phase sets | 698 | 697 |

Exactly 32 VCF sample fields change, all between 19,397,607 and
19,464,219. Candidate keys and genotypes are unchanged. The two new gap
cases and the anti-join control are in the window panel; the earlier gap
expectations retain their previous floors.

Validation: `make -j4`, `make check`, `make unit-tests`, and the complete
`make window-tests` suite passed (3,396 assertions in 42 test cases).
