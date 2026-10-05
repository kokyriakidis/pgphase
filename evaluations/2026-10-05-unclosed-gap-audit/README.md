# New-gap search: no accepted production change

The starting and restored binaries have SHA256
`bdce74cda1d0ec407fe7a0cfaf9e30d2c77a60b605d6f193ee92e2cb16503933`.
The full-chromosome baseline is `test_data/tmp_gap_fix55/full-final`.
All trial changes were removed; this investigation does not close a new gap.

The corrected search uses the maximum covered endpoint of all preceding phase
blocks. Adjacent blocks in sorted-start order are insufficient: an earlier,
longer block can already cover their apparent gap. The screens cover 55 simple
SNP pairs, 139 simple SNP/indel or indel pairs, and 12 SNP/complementary-insertion
pairs at actual uncovered boundaries no longer than 50 kb. None passes its
unanimous, two-haplotype physical-call nomination criterion. The complementary
screen accepts sequence-equivalent placement within 16 bp. All screens exclude
secondary, supplementary, duplicate and QC-failed reads and require MAPQ and
called base qualities of at least 30. The exact-indel screen is intentionally
conservative and is not an exhaustive sequence-equivalence proof.

A discarded trial connected the BAM-only insertion component at 2,311,483 to
its neighboring graph block at 2,298,732 using six consistent physical calls
and the existing supported bypass of one weak graph SNP. This was not a new
gap: the original 1,684,336 phase block already spans through 2,316,265 in the
full chromosome. Its owning replay likewise already spans the proposed gap.
The independent parental audit also rejects the insertion-side flank (12
matching and four discordant votes). Exact variant keys, genotype alleles,
non-phase fields, all old block gauges, previously correct read assignments
and output-only rescue tags were unchanged in that trial. It was still rejected.

Broadening the independent insertion-observation retry to either flank also
produced identical VCFs in the owning 12, 19 and 52 Mb replays. Direct checks
of the 15.023--15.040 Mb deletion/SNP boundary and the 17.839--17.852 Mb SNP
boundary found conflicting physical allele combinations. Repeat-insertion
probes in the 22, 33 and 44 Mb chunks likewise found both parities. The 34,
38 and 56 Mb owning-boundary replays did not supply an outer-SNP physical link.

The three Python scripts reproduce the final nomination screens using pysam;
run them from the repository root with the derived test data available. They
write their JSON into `test_data/tmp_gap_fix56/`. Saved JSON and the rejected
owning audit are adjacent to this document. No expected span or parental floor
was weakened and no window was added for the rejected candidate.
