# Verified deletion bridge regression

The recovery block-bridge validator previously supported MSA insertions only:
every deletion observation returned false. At chr20:53,926,537–53,962,801,
the missing signal is the MSA-verified deletion at native event position
53,943,646 (VCF anchor 53,943,645, deleted sequence GT).
Read `m84031_231217_034919_s2/94048220/ccs` observes the reference allele there
and reaches the right-block SNP with Q40 bases. Both original blocks are
perfect against parental read truth and require opposite HP orientations.

Recovery now validates the exact BAM deletion or complete aligned reference
allele, requiring Q30 flanking bases and, for REF, Q30 bases throughout the
event. Reference skips, shifted/overlapping unmatched events, and homopolymer
indels are not accepted by this path. A validated deletion can anchor a
singleton bridge to a strong clean-SNP component; contradictory component
votes still veto it. The existing phasing and deferred stitching apply the
result. No competitor calls or parental labels enter phasing.

## Results

`before_after.tsv` compares the previous `graph_bam_v2` run with
`deletion_bridge_v1`, using corrected endpoint coordinates in both runs.
The two arms within `deletion_bridge_v1` share this core fix; their difference
is only whether graph-selected BAM fallback is enabled.

- Local 53.9 Mb: split → joined, 486 evaluated reads, zero discordant reads
  and zero switch/flip errors before and after.
- All 11 targets: 1 joined, 9 split, 1 unphased endpoint. No read-error increases.
- Full chr20: 122 → 123 of 280 initial gaps joined; no previous join lost.
  Read discordance remains 2,731/189,001; switches 353 and flips 413 are unchanged.
  Cached solve time: 197.722 s, compared with 195.14 s previously.
- Whole-chromosome output: `/tmp/pgphase-deletion-bridge`; regional runs:
  `/tmp/pgphase-gap-trials/runs/deletion_bridge_v1`.

The trial endpoint checker had compared VCF anchors with native insertion
positions. Native indel positions require subtraction of one. Thus the prior
three unresolved endpoints were two phased insertions that were missed by the
checker (12,721,112 and 37,984,529), plus one genuinely unphased insertion
(23,480,815). The corrected before-state is 10 splits and 1 unphased endpoint,
not 8 splits and 3 absent candidates. No nearby-position approximation is used.

## Remaining evidence and rejected trials

- 1.08 Mb has a single primary spanning read with MAPQ60 but Q10 at the only
  observed right-block heterozygous SNP (1,110,921).
- At 37.98 Mb, read `m84031_231217_034919_s2/158272567/ccs` has confident
  observations at the verified insertion (37,984,529 anchor), deletion
  (38,005,399 anchor), and SNPs 38,005,401 and 38,008,026. Its right-block
  observations imply opposing haplotypes. The component unanimity test drops
  that anchor; simply allowing biallelic insertion singleton bridges does not
  resolve it. This conflict remains unresolved.
- Increasing the link window from 8 to 64 did not close additional targets.
- Allowing validated biallelic insertion singleton anchors, alone or with
  reseeding supported biallelic indels, did not close additional targets.
  Both experimental changes were removed.

`--verbose 2` now reports `GapBlockVotes` (block anchors, two orientation
counts, strong counts, clean counts); `--verbose 3` additionally reports
`GapAnchorObservation` (read, anchor, observed allele, HP1 allele, HP2 allele,
BAM-confidence flag). Matrix dumps include comment metadata with MSA,
homopolymer and link flags; the dump precedes link selection.

Build and unit tests pass. Added regression cases cover REF/ALT orientation,
low quality, reference skips, contradictory cached alleles, opposing bridge
reads, homopolymer exclusion, and indel endpoint coordinate matching.
