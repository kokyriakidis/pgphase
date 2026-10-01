# Why the 50.548 Mb deletion seam still abstains

The accepted full chr20 graph/BAM output leaves `50548245:CA>C` and
`50548245:CAA>C` in PS 50,491,033 and `50562066:A>G` in PS 50,562,066.
HiPhase joins the boundary rows on the same callset; the committed gap panel
records its 101/119 local truth-correct separation. This audit used the
50–51 Mb owning chunk and the same input BAM. Parental truth was used only
to score output, never to select a stitch.

Twelve MAPQ-at-least-30 primary reads span the deletion and right SNP. Their
BAM source and graph read labels agree 6/6 on both haplotypes, but this is
not a direct allele bridge: all twelve have allele `-1` for **both** deletion
rows in the BAM subsolve and transferred recovery matrices. Most begin after
the last informative left site, so their inherited BAM source HP cannot
certify the deletion-to-SNP relationship. The source matrix and owning-chunk
outputs are under `/tmp/pgphase-50548-current/` while available.

CIGAR inspection found one- and two-base deletions shifted two bases into
the A repeat in several spanning reads. The complementary-deletion stitch
currently asks for exact CIGAR placement. A trial reused the existing
reference-equivalent deletion caller and removed a redundant exact-coordinate
flank check. Under its Q30 reference-base checks, only four independent
paired calls survived: one for one deletion allele, three for the other.
Their relative phase votes were 1 versus 3 (quality-weighted log odds 13.62),
so the current rule's unanimity check correctly abstained. The right graph
SNP path also failed both its ordinary certificate and its existing complete
BAM source fallback. Temporarily bypassing that path check still did not
join because the four allele votes conflict. The trial was discarded; the
source and binary were restored to the previous validated implementation.

A scan of current open panel gaps found no decisive Q30 clean-SNP pair
within 15 kb of both boundaries. At 21.159–21.172 Mb, the strongest such
pair calls the left SNP in both allele classes but calls the right SNP REF
on all 16 reads (9 and 7), so it cannot orient the right block. This is a
prioritization aid, not a catalog of every indel or longer-range chain.

A future solution needs local diploid allele evidence across the repeat and
an independently certified path through the right block. Whole-read source
labels or one-sided CIGAR matches alone are insufficient to justify merging
the two established phase sets.
