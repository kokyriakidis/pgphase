# chr20 60.033 Mb: shifted insertion bridge

The tracked `60,033,052 C>G` to `60,048,237 A>AA` gap remained split in the
full chr20 and owning 60–61 Mb chunk. The right insertion is MSA-verified in
the recovered BAM source. Several primary reads put a one-base insertion at
60,048,231 or 60,048,234 in the same A run. Inserting A there produces the
same local sequence as the candidate insertion at 60,048,238. The old final
stitch required a clean SNP on the right and never tested this edge.

With MAPQ and base quality at least 30, seven read pairs call insertion REF
and two call sequence-equivalent insertion ALT. Their quality-weighted
relation favors a flip of the right block. The left phase set contains one
original graph SNP at 60,005,003 and attached BAM sites through 60,033,052.
The old path check rejected it because it required two original graph SNPs,
although there is no internal graph SNP edge to validate. The right graph
path passes.

The new bridge uses edited-reference equivalence within 16 bases, checks all
inserted and flanking bases, rejects other nearby indels, and requires at
least two independent reads for each insertion allele. Its left singleton
certificate also requires exactly one original graph candidate in that phase
set; the right path uses the ordinary graph SNP check. SNP pairs keep priority
and no candidate site row is rewritten.

The short regression window changes from split, 0.46 truth-separated reads
to joined, 0.71. The full owning chunk keeps all 3,326 truth-scored reads with
3,290 correct and 36 discordant before and after. Its read phase sets fall
17 to 15. The initial full chr20 trial raises tracked closures 14/23 to
15/23 and reduces read phase sets 711 to 709, with all 236,847 individual
truth outcomes unchanged (229,090 correct, 7,757 discordant). All 62,154 VCF
keys are unchanged. Eight phase genotypes invert together within the joined
block, with no dosage changes. The final singleton guard narrows the scope of
that trial; its validation is recorded below.

The older 60,033,052–60,058,235 panel window also closes. It has 0.75
truth-separated reads versus 0.40 before, and its interior insertion at
60,048,238 is recorded as a required site. Both 60.033 Mb windows raise the
graph panel's span count from 42 to 44.

Final validation after the narrowed singleton guard reproduces the initial
full-chr20 result: 15/23 tracked closures, 709 read phase sets, 62,154
identical VCF keys, eight phase-only genotype inversions, and no individual
read-truth changes among 236,847 scored reads (229,090 correct; 7,757
discordant). `make unit-tests` and `make check` pass. The full
`make window-tests` suite passes 1,793 assertions across 34 cases.
