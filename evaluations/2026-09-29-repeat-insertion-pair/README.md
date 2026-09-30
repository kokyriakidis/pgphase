# Repeat insertion pair at chr20:21.378–21.395 Mb

The last open member of the 23 HiPhase-correct tracked gaps has two
MSA-verified BAM insertion rows at its boundaries: 21,377,985 C>CTGTGTG
and 21,395,286 A>AT. The ordinary physical stitch only visited graph SNP
boundaries, so it did not evaluate this BAM-only pair. Two distinct MAPQ-60
primary reads span both insertions. Their local CIGAR length classes agree:
one calls the left +2 bp slippage and right zero length, and one calls zero
length at both. The left source has an earlier weak edge at SNP 21,343,641;
two independent MAPQ-60 reads call that SNP and the next insertion REF,
certifying the inherited left orientation. The right source has its weak
edge at the singleton insertion.

The new seam pass uses those primary BAM length observations to orient the
two insertion blocks. It requires two agreeing molecules, zero opposing
votes, an exact zero-length call at each boundary, and a repeated 1–2 bp
motif in the reference. Multiple nearby indel events abstain. Only right
source reads with a callable primary BAM boundary allele enter the joined
phase set; seven other reads retain their HP in a separate read-only phase
set. Eleven reads that appear to start after the boundary in graph
projection actually cross it in the primary BAM and correctly stay joined.

An unguarded trial also joined a different 23.792–23.806 Mb seam. Its right
boundary has both an insertion and a complementary deletion at the same
locus. An insertion-only length vote cannot represent that diploid allele
choice. The final stitch vetoes a boundary with another phased allele at
the locus or an overlapping deletion. The guarded full run preserves the
23.8 Mb rows and their original phase sets.

| Full chr20 measure | Previous accepted run | Guarded run |
|---|---:|---:|
| VCF variant keys | 62,154 | 62,154 |
| Truth-scored phased reads | 236,859 | 236,859 |
| Truth-correct phased reads | 229,101 | 229,108 |
| Discordant phased reads | 7,758 | 7,751 |
| Read phase sets | 699 | 699 |
| Tracked HiPhase-correct gaps closed | 22/23 | 23/23 |

Only the 21,395,286 A>AT VCF sample field changes: 0|1 with
PS=21,395,286 becomes 1|0 with PS=21,334,973, matching the left
insertion's ALT haplotype. The owning-chunk replay has 122/144
truth-scorable reads in its dominant parental orientation (84.72%);
HiPhase has 131/144 (90.97%) there. The join closes the tracked gap
without a local parental switch, although the local read purity remains
below HiPhase.

The panel now includes the exact 21.378 Mb owning-chunk bridge and a
23.8 Mb multiallelic counterexample. The former asserts shared phase
set and parental orientation; the latter asserts that the insertion and
complementary deletion remain paired while the earlier block stays split.
Validation: make -j4, make unit-tests, make check, targeted owning-chunk
regressions, and the full make window-tests panel (2,939 assertions in
41 test cases) all pass.
