# Complementary BAM deletion rows bridge a chr20 graph seam

The tracked chr20:60,453,499–60,467,115 gap has a clean right graph SNP and
two separate MSA-verified BAM deletion rows attached to the left graph block:
8 bp at 60,453,494 and 2 bp at 60,453,500. Their reference spans overlap,
and their ALT alleles occupy opposite phased haplotypes. A REF call at one row
cannot represent the other haplotype, because the other deletion removes part
of that row's reference span. The previous SNP-only physical bridge had no
clean left SNP in this seam and therefore could not use either deletion row.

Ten MAPQ >= 30 primary molecules cross the deletion/SNP boundary. Six make an
exact, quality-checked 8 bp CIGAR deletion call and right SNP ALT call; two
make an exact 2 bp deletion call and right SNP REF call. The other two carry
different nearby deletions and do not vote. The eight accepted votes agree on
one phase relation. Both original graph blocks pass the existing continuous
clean-SNP path check, and the BAM deletion rows were already attached to the
left block by the source-path/gauge checks. The stitch uses both rows as they
are; it does not merge or recode them.

A second bug appeared while auditing the physical calls: `physical_snp_call`
compared uppercase BAM bases to raw reference characters, so soft-masked
lowercase reference bases could never produce a REF call. The base comparison
now normalizes case. The new deletion-flank check uses the existing worker
reference cache, which returns uppercase bases. Exact CIGAR calls and matching
reference flanks are both required before a read votes.

The 60–61 Mb owning chunk changes its read phase sets from 18 to 17, with
3,290 correct and 36 discordant of 3,326 truth-scored reads in both runs.
The tracked gap now has a single phase-set ID and opposite boundary ALT
haplotypes. Full chr20 tracked HiPhase-correct gap closures rise 11/23 to
12/23; read phase sets fall 716 to 715. All 62,154 VCF keys and genotypes are
identical to the preceding build. The same 236,847 truth-scored reads retain
their individual correctness outcomes: 229,090 correct and 7,757 discordant
(96.7249%). The measured graph-window span count rises 39 to 40, and this
window's separated fraction rises from the 0.54 floor to 0.89.

Reproduce with `collect-graph-variation` over
`CHM13#0#chr20:60000001-61000000` and `CHM13#0#chr20:1-66210255`, using
the chr20 reference, surjected BAM, striped site VCF and coordinate GAF in
`test_data/`. Score with `/tmp/score_pgphase_weakcut.py` or run the owning
chunk and full window regressions.

Validation: `make -j8`, `make unit-tests`, `make check` and `make window-tests`
pass. The full window suite reports 1,754 assertions across 34 cases, including
the owning 60–61 Mb chunk and the unsafe 62.6 Mb control.
