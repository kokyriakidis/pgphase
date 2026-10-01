# Graph indel as an independent two-hop bridge (2026-09-30)

The tracked HiPhase-correct chr20:35,328,965–35,347,817 gap remained open in
the 35–36 Mb owning chunk. Both boundaries are imported BAM phase blocks inside
the larger graph seam 35,079,659–35,490,918. The existing outer seam
transaction rolled back, and its fallback intentionally left adjacent imported
BAM blocks independent. The exact BAM-to-BAM allele pair was inconclusive.

At the intermediate graph catalog site 35,342,612, the selected padded REF and
ALT strings give the same edited reference sequence as the BAM `G>GAGATAGAT`
insertion at 35,342,608. These remain two distinct candidate rows. The graph
site's read alleles form a consistent chain: 15 left-boundary pairs are 9
ALT/ALT and 6 REF/REF; 11 pairs to the right clean SNP at 35,357,391 are
9 ALT/ALT and 2 REF/REF. All these reads are MAPQ 60. The BAM MSA insertion
row calls 28 of its graph-ALT reads REF, explaining why an exact BAM-row
stitch missed this path.

The new fallback uses an unphased graph indel only as an allele bridge between
two adjacent imported BAM blocks. Each hop requires a one-sided exact binomial parity p-value at most 0.001 on
MAPQ-30 calls, both allele classes represented at least twice, no opposing
pairs, and
at most 25 kb to the graph site. The left BAM anchor must be clean or
MSA-verified; the right must be a clean SNP. Multiple qualifying paths must
agree. A significant direct BAM-block vote of opposite parity vetoes the join.
The graph candidate stays unphased and the BAM sites stay independent.

| Measure | Previous full chr20 | Two-hop trial |
| --- | ---: | ---: |
| VCF variant keys | 62,361 | 62,361 |
| VCF PS fields changed | — | 7, all 35.342–35.370 Mb |
| Truth-scored tagged reads | 236,880 | 236,878 |
| Correct | 229,590 | 229,592 |
| Discordant | 7,290 | 7,286 |
| Read phase sets | 691 | 689 |

In the owning 35–36 Mb chunk, correct reads changed 3,189 to 3,191 and
discordant reads 174 to 170. The target's VCF boundary rows now share PS
35,305,645. The owning-chunk gap test checks the span, parental orientation,
read concordance, and the required interior BAM insertion. Its previous 100%
concordance expectation came from a different 100 kb replay; the old owning
chunk scored 94.826% and the joined one 94.942%, so the new 94.9% floor is
stricter than the comparable baseline.

The first normalization trials changed a distant 23 Mb block and closed no
tracked gap; they were reverted. The new rule uses graph observations as a
bridge without changing site keys or candidate ownership.
