# chr20 19 Mb: a missing allele is not a phase-block anchor

`update_read_phase_set` previously selected the first phased heterozygous
variant in a read profile's index interval without inspecting that read's
allele at the variant. Profiles are sparse: -1 means no observation, not a
reference call. Reads inferred from right-block variants could therefore
inherit a left-block PS, supplying false overlap votes to gap stitching.

A four-read/two-disconnected-site unit test reproduces this bug. Each read
observes only one site, but its profile spans both indexes. Before the fix,
right-site reads inherit the left site's PS. The corrected owner requires an
observed allele matching one of the site's two phased alleles. The test fails
before and passes after the change. This fixes the common core, including
initial phasing, rather than adding a recovery-only confidence threshold.

In this region, nine primary MAPQ60 reads span MSA insertions 18983382/18983415
and clean SNP 18999993. Eight have -1 at both MSA sites in the old matrix and
an observed right SNP; only one has observations at both sides. The false
join's 208 discordant reads disappear with the ownership fix, while all 489
truth-assessed reads are retained. The region correctly remains split pending
recovery of actual bridge observations; it is not yet phased to competitor
contiguity. Whole-MSA windows cover the sites (18983305–18983451), so the missing
observations require alignment/representation diagnosis, not a claim that no
spanning reads exist.

## Subsequent alignment fix and local validation

Composing read-to-consensus and consensus-to-reference alignments left
sequence-equivalent repeat insertions at different coordinates. Left-normalizing
the composed alignment recovers the missing direct allele observations. The
18983382 insertion now links to SNP 18999993 with 3 reads supporting the correct
orientation and 1 conflicting read, instead of a single usable observation.
With observed-allele phase-set ownership and this normalization, the
`chr20_local_alleles` and `chr20_unassigned_local` local runs join the target:
489 truth-assessed reads, zero discordant reads, and no reversed original block.
This is a local result; chromosome-wide corrected recovery remains unvalidated.
