# chr20 35.9 Mb: verified indel lost between recovery tiers

The native gap is 35919404–35939110. HiPhase, WhatsHap, WhatsHap-opt, and LongPhase all phase both endpoint SNPs
with the same relative orientation in one block. Its only phased heterozygous site strictly inside that
gap is an A insertion at VCF position 35920927. The graph catalog contains the
corresponding multiallelic repeat site at 35920926; truth left-normalizes the
insertion to 35920913 (G→GA). MSA independently verifies it at internal event
position 35920914.

## Root cause and correction

The SNP-only tier computes this indel but does not admit it, as intended.
It then extends the left flank to a recovered SNP near 35927397. The indel
tier incorrectly narrows its MSA search using that new boundary and omits the
previously examined insertion window (35920904–35920944). The fallback can then
join through a misleading T deletion near 35929194 instead, reversing all 48
truth-assessed reads in the original right block.

Recovery now keeps the original gap's MSA window across tiers. Gap boundaries
may still move for stitching, but tier escalation does not discard previously
examined evidence. The insertion is admitted in tier 3 with 37 reference and
22 alternate observations. The gap joins before the repeat fallback, with both
flanking SNPs 1|0 in PS 35873360, matching truth up to a block-wide label swap.

A separate allele-representation fix restores the non-deleted allele of the
nearby residual TTCC deletion at 35927394 when its fixed consensus contains
an overlapping C→T substitution. Literal-reference matching previously lost
these observations (1 reference / 6 alternate); verified consensus matching
recovers 19 reference / 7 alternate. Unsupported sequences and conflicting
alignment paths remain unknown. This fix alone did not correct the join.

## Validation

`chr20_fixed_tier_window`: one block, all 293 truth-assessed reads retained,
48→0 discordant reads, 3→0 switch/flip errors. The 18-case regression panel
adds no discordance or switch/flip errors relative to `chr20_best_site_anchor`.
The expanded 114-case panel also introduces no discordance relative to the
prior full `chr20_direct_site_anchor` trial. Known wrong joins at 46.7 and
60.1 Mb remain. Build and all unit tests pass; the overlapping-SNP allele
regression fails before the correction and passes afterward.
