# Open chr20 gap bridge screen (2026-09-30)

The graph plus BAM recovery pipeline was replayed on owning 1 Mb chunks with
`--phase-matrix-dump`. The same aligned BAM and graph catalog used by the
committed gap panel were used; no parental truth entered the phasing run.
Numbers below count paired, nonnegative matrix allele observations unless
specified otherwise. These are diagnostics, not sufficient evidence to join
whole phase sets.

| Boundary positions (chr20) | Paired evidence | Finding |
| --- | --- | --- |
| 19,373,922–19,395,544 | 1 pair per complementary left insertion row to the SNP | The shared BAM source PS is not an independent bridge; its direct boundary link is weak. |
| 64,138,752–64,140,314 | 12 paired deletion calls: 11 support one relation, 1 opposes | The BAM source label has the opposite orientation across the cut. The existing regression requires the split because the downstream block has a separate weak edge. |
| 35,498,368–35,516,845 | 0 matrix pairs; 3 MAPQ 60 reads physically cross | One read has an insertion at the deletion locus. The other two are deletion REF/SNP ALT, but fail Q30 flank checks at Q10 and Q22. Lowering the flank gate and using their measured quality did not close the gap with existing path checks. |
| 15,351,845–15,367,755 | 6 pairs, all left REF | No support for the left ALT haplotype. |
| 21,130,205–21,148,338 | 3 pairs, all left REF | Two right REF, one right ALT. |
| 21,594,343–21,612,458 | 1 pair | Insufficient direct bridge. |
| 21,823,066–21,844,359 | 0 pairs | No callable boundary bridge. |
| 23,792,418–23,806,565 | 4 pairs, all left REF | The right insertion has both calls; the left ALT is unobserved. |
| 49,720,311–49,742,973 | 1 clean SNP pair | Insufficient direct bridge. |
| 58,366,458–58,385,702 | 0 pairs | The original BAM source PS spans the interval, but recovery has no callable boundary edge. |
| 3,963,987–3,981,464 | 0 BAM-deletion/SNP matrix pairs; 7 physical spanning reads | Four reads have a nearby one-base deletion, three do not; SNP bases in both groups are mixed. The graph's six-base deletion row is a different allele and must not substitute for the BAM one-base row. |
| 52,414,285–52,422,209 | 30 insertion-row/SNP pairs | Complementary insertion rows have mixed SNP votes; right SNP ALT depth is 14 of 63. |

Two reversible trials changed `stitch_msa_indel_to_snp`: admitting one allele
class with a decisive likelihood, and then also lowering deletion-flank base
quality to Q10 while weighting by observed quality. Both left the 35.498 Mb
boundary split. The source was restored and rebuilt after each trial. A
statistical graph-deletion vote change was considered but not applied because
the examined 60.706 Mb boundaries are BAM-derived and have zero callable pairs;
it would not address the tracked case.

This screen does not claim the remaining gaps are impossible. It rules out
simple whole-block joins from these named boundary pairs. Further work needs
a validated intermediate allele chain or a read-based gauge that certifies
both sides of a weak internal edge.

## Exhaustive short-boundary screen

I screened the 40 remaining gap coordinates in
`/tmp/pgphase-current-open-gaps.tsv` against primary BAM alignments with
MAPQ >=30. The diagnostic caller used Q20 SNP bases and exact one-event
CIGAR matches within 3 bp for short indels. It deliberately did not attempt
repeat normalization, so these counts are a prioritization screen, not stitch
certificates. The most significant raw paired counts were 17 cross versus 3
same at 17,839,398–17,839,399 and 26 same versus 9 cross at
60,706,792–60,711,400. Both gaps have opposite parental flank orientations
in the truth-only audit; accepting a raw vote majority would risk the wrong
whole-block connection. The high-purity 10.727 Mb gap has only two matrix
pairs to its right insertion, both REF on that side. The 57.854 Mb gap has
no callable matrix pair from its BAM deletion rows to the right SNP.

Read-tag comparison found 23 truth-labeled reads in the 57.854 Mb window
that HiPhase tags and pgphase leaves untagged; ten of those HiPhase tags
disagree with the majority parental orientation of its local phase set. At
13.593 Mb, twelve HiPhase-only truth reads all call REF at both competing
right-insertion rows in the pgphase matrix. At 23.792 Mb, nine HiPhase-only
truth reads likewise call REF at both right-insertion rows, and the interval
has no heterozygous interior SNP. Assigning these reads from either insertion
row would use an allele they do not carry.

The gap test helper now shares identical owning-chunk replays within a test
process. Two neighboring 64 Mb windows passed 40 assertions in 17.08 s,
versus 16.80 s for one window alone; their VCFs were hard links to the same
replay artifact. This changes test execution only.
