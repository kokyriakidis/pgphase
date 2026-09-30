# Chr20 SNP-to-deletion gap stitch

At chr20:23,421,003–23,445,252, recovery preserves a clean left SNP and an
MSA-verified right 9 bp BAM deletion in separate phase sets. The physical
stitch previously handled SNP-to-insertion boundaries but had no
SNP-to-deletion case. The recovery window ends at 23,445,261, nine bases after
the right deletion's VCF anchor, so the deletion is selected by membership in
the window rather than equality with its end.

Only one primary BAM read spans both exact boundary alleles. It has MAPQ 60,
left ALT at Q17, and the complete right deletion REF at Q40. Conservative Q30
for the deletion gives parity log odds 3.84422, or about 2.1% probability of
the opposite relation under the read-error model. The right BAM source has a
complete path with no weak or quality cuts, and both graph SNP paths are
supported. The bridge rejects opposite-haplotype overlapping deletion rows.
The owning 23–24 Mb chunk joins the exact VCF alleles. Its 116 locally tagged
reads retain 114 correct parental assignments in one phase set, compared with
the same 114 correct reads split across two phase sets before. Four additional
chunk reads receive truth-correct tags.

The first chromosome trial allowed a complete BAM source without validating
its right graph path. It joined a different seam at 1.152–1.197 Mb, where two
reads had log odds -10.65 but the right graph SNP path was unsupported. That
early join caused the previously closed 1.508 Mb target to reopen; the broad
trial closed only 21/23 tracked gaps and changed 720 VCF sample fields. The
right graph-path check rejects this counterexample and retains the target
23.421 Mb join.

Full chr20, before → guarded result:

| Measure | Before | Guarded result |
|---|---:|---:|
| VCF variant keys | 62,154 | 62,154 |
| Truth-scored phased reads | 236,855 | 236,859 |
| Truth-correct reads | 229,097 | 229,101 |
| Discordant reads | 7,758 | 7,758 |
| Read phase sets | 699 | 699 |
| Tracked HiPhase-correct gaps closed | 21/23 | 22/23 |

Only three VCF sample fields change, all in the right 23 Mb block: its
phased genotypes reverse to match the left block, while allele dosages stay
unchanged. The only tracked open gap is now
chr20:21,377,985–21,395,286. Its two physical crossing reads each have an
uncertain allele at a different boundary, so the same stitch rule does not
apply.

The graph panel's 23.421 Mb window replays the whole 23–24 Mb owning chunk.
It scores 3,495/3,601 phased reads correctly (97.06%) across that chunk and
places 114/147 local truth-scorable reads correctly in one phase set (77.55%).
The previous 99% concordance floor came from a shorter padded replay; the
owning-chunk floor is 97% and the separated-read floor rises from 38% to 77%.

A standalone `collect-bam-variation` run over chr20:21–22 Mb also keeps the
remaining 21.378 Mb insertion pair in separate phase sets (21,334,973 and
21,395,286). The graph recovery transfer preserves those same rows and labels,
so the remaining split is already present in the BAM solve. The two spanning
reads have complementary quality limitations: one carries a distinct nearby
2 bp insertion and has Q3 at the right REF flank; the other has Q10/Q22 at
the left REF flank and Q40 at the right. A stronger connection needs a
representation-aware likelihood or an independent graph allele observation,
not a phase-set transfer change.
