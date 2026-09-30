# Exact BAM insertion closes the chr20 19 Mb chunk-edge seam

The chr20:18,983,414–18,999,993 gap has a left MSA-verified BAM insertion
(`G>GTGTGTA`) attached to the earlier graph block and a clean right graph SNP
(`A>G`). The physical seam stitch previously required two clean SNPs or two
complementary deletion rows. Recovery retained the insertion but could not
orient it to the right SNP. In the owning 18–19 Mb graph chunk, the right phase
set contains exactly one graph candidate, seven bases before the chunk edge.
The two-site graph-path check had no internal edge to evaluate.

Nine MAPQ-60 primary molecules cross both sites. Three have an exact six-base
CIGAR insertion and right SNP ALT; four have insertion REF and right SNP REF.
The eighth callable pair is insertion REF with right SNP ALT and carries a
different nearby four-base insertion. One other read has a different insertion
at the target and abstains. All accepted bases have quality at least 30.
Parental truth places the exact insertion and right SNP ALT on the paternal
haplotype. The new bridge preselects the nearest eligible MSA-verified BAM
insertion, checks its exact CIGAR allele and both reference flanks, requires
both insertion alleles and the existing 0.001 wrong-parity likelihood bound,
and certifies the complete left graph SNP path. It accepts a singleton right
graph block only when it has exactly one candidate in the chunk, so no
internal right path can be reversed. Existing multi-site right blocks still
require the full graph path check.

The exact 18–19 Mb owning chunk joins the two boundary rows with ALT on the
same haplotype. It retains 4,254 correct and 41 discordant among 4,295
truth-scored reads. The short panel replay joins too, with 81/151 local
truth-scorable reads correctly separated in one block (0.536); its required
normalized insertion site is 18,983,415. Full chr20 tracked closures rise
12/23 to 13/23 and read phase sets fall 715 to 713. All 62,154 variant keys
and genotypes remain the same. The same 236,847 truth-scored read names keep
their individual correctness outcomes: 229,090 correct and 7,757 discordant
(96.7249%). The graph panel span total rises 40 to 41.

Reproduce using `collect-graph-variation` on
`CHM13#0#chr20:18000001-19000000` and `CHM13#0#chr20:1-66210255`, with the
chr20 reference, surjected BAM, striped site VCF and coordinate GAF in
`test_data/`. The owning-chunk regression is `physical bridges survive their
owning graph chunks` in `src/test_gap_windows.cpp`.
