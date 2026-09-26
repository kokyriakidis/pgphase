# Five open chr20 gaps: HiPhase signal audit

The comparison uses the saved HiPhase 1.6.0 run in
`/tmp/hiphase-on-pgphase-final-chr20/`, which phased pgphase-derived variant
calls against the original annotated BAM. The current pgphase result is
`/tmp/pgphase-complete-focus-chr20/`. The original alignment is
`test_data/HG002_chr20_hifi_mapped_to_CHM13_chr20_annotated.bam`.

For every gap below, the HiPhase input VCF contains exactly the two boundary
heterozygotes and **no internal heterozygote**. HiPhase gives the two boundaries
one phase-set ID, and its output BAM tags the sole primary read spanning both.
Its run log records `Minimum spanning reads: 1`, `Minimum connecting reads: 1`,
`Minimum call quality: 0`, and `Minimum mapping quality: 5`. Thus the same-sites
HiPhase join is explainable by one molecule; there is no hidden chain of
internal HiPhase variant calls. Fetching all primary, secondary, and
supplementary alignments at both endpoints finds the same one shared read
name in each gap, so split mappings add no other endpoint bridge. The primary
BAM read counts count reads intersecting each interval, including reads that
reach only one boundary.

| Boundary positions | Primary BAM reads in interval | Reads spanning both | Spanning-read boundary BQ | BAM source observations | HiPhase spanning-read HP |
|---|---:|---:|---|---|---:|
| 1,086,625–1,110,921 | 170 | 1, MAPQ 60 | 40 / 10 | ALT / REF | 2 |
| 3,573,979–3,596,663 | 161 | 1, MAPQ 60 | 40 / 10 at deleted-base REF | ALT / missing | 1 |
| 21,435,750–21,456,826 | 157 | 1, MAPQ 60 | 40 / 40 | REF / REF | 1 |
| 21,736,539–21,757,517 | 143 | 1, MAPQ 60 | 40 / 40 | REF / REF | 2 |
| 61,664,136–61,690,751 | 205 | 1, MAPQ 60 | 40 / 17 | REF / ALT | 2 |

The 3.573-Mb source matrix lacks the right deletion observation on the bridge
read. Its alignment carries the complete reference allele, but the deleted
base has BQ10. The ordinary reference-call path in
`collect_read_var_profile` does not test reference base quality. The missing
observation could instead come from an overlapping read event or a read-level
noisy interval; that branch has not yet been traced. HiPhase's saved outputs
do not expose its per-read realigned allele call.
The other four bridge molecules retain both boundary observations in pgphase's
saved BAM source or local phasing matrix. Pgphase's graph seam still abstains:
its direct physical-SNP path requires Q30 calls and an independent clean SNP
at least 100 bp from the selected boundary on one flank. The 21.736-Mb left
boundary is an indel and cannot enter that SNP-only path. A previous trial
that accepted one Q40 molecule without corroboration made a wrong whole-block
join near 62 Mb and misplaced 727 reads; see the larger-gap panel README.

## Graph-only deletion at 21,742,440

The full 21-Mb graph recovery matrix is
`/tmp/pgphase-complementary-chunk21/matrix.chunk0.recovery-final.tsv`.
Its graph-only deletion is a repeat-indel candidate (category mask 16),
with 43 REF and 14 ALT observations; it is absent from the HiPhase input VCF.
The matrix has 34 callable left-boundary/deletion read pairs and 6
callable deletion/right-boundary pairs. Those are close to the 38 and 8
primary BAM reads physically crossing the corresponding coordinate pairs.
Four and two physical crossers, respectively, lack a call on that **deletion
row**. All six have an indexed GAF path through the distinct `CT→CTT` allele
(`>115804235>115804237>115804238>115804239`), whereas the deletion row
represents `CT→C` (`>115804235>115804239`). Marking those reads unknown on
this biallelic row is correct; counting them as deletion REF would be wrong.
The earlier claim of zero paired observations was incorrect.

The deletion does not currently give a clean diploid bridge. Among the 34
left/deletion pairs, graph allele counts are (left REF, deletion REF)=11,
(REF, ALT)=7, (left ALT, deletion REF)=14, and (ALT, ALT)=2. Thus both left
haplotypes predominantly call the deletion REF. The only full bridge read has graph observations REF/ALT/REF. Its BAM alignment
has a one-base T deletion at 21,742,441, matching the graph's `CT→C` allele:
the `CT` REF includes the anchoring C and does **not** describe a two-base
deletion. Its BAM source has only the boundary REF/REF calls. The earlier
claim that the bridge's CIGAR and graph deletion differed in length was wrong.

The reference has a 14-T run at 21,742,441–21,742,454. Among 64 primary BAM
reads covering both flanks of the run, direct read-sequence lengths compared
with graph calls on the deletion row are:

| Read T-run length | Graph deletion REF | Graph deletion ALT | Absent from row |
|---:|---:|---:|---:|
| 14 (reference) | 43 | 0 | 0 |
| 13 (one-T deletion) | 0 | 13 | 0 |
| 15 (one-T insertion) | 0 | 0 | 7 |
| 12 (two-T deletion) | 0 | 1 | 0 |

Thus the graph's allele representation agrees with the read sequence for
the 13 one-T deletion reads and all 43 reference reads. One 12-T read is
coarsely assigned to the one-T deletion walk. The 15-T reads belong to the
separate insertion allele and must stay unknown on the deletion row. These
are graph-path observations, **not MSA-verified BAM calls**: candidate 361 has
`msa_verified=0` in the graph matrix, and the saved BAM recovery-source
matrix has no candidate at this run. `vcf_to_variant_key` would normalize
`CT→C` to a one-base deletion at 21,742,441 if a matching BAM candidate
were available to transfer. The graph site's weak segregation with the left
boundary, not a two-base representation mismatch, remains the measured
limit on its use as a bridge.

The immediate deficit is the **stitch decision**, with one missing BAM
allele at 3.573 Mb. That site is MSA-verified but flagged as a homopolymer
indel, which the targeted exact-CIGAR backfill currently excludes. The
original BAM does not contain more physical boundary-spanning reads than
pgphase sees. A future singleton join rule needs an allele-likelihood check
and the 62-Mb wrong-join control before it can be used on whole phase blocks.
The graph deletion's present calls do not segregate well enough to justify
using it as a bridge.

## Follow-up on the suspected dropped calls

The six GAF paths above all traverse the catalog's third allele, `CT→CTT`.
The default graph projection keeps that allele unknown on the separate
`CT→C` row. A unit regression now checks this distinction. These six rows
were not lost coverage and must not be filled as REF.

The 3.573-Mb right boundary is an MSA-verified deletion in a 19-base A run.
The targeted CIGAR backfill deliberately excludes homopolymer indels. A trial
required one read to match the entire reference run and both distinct flanks
before adding REF. It restored the bridge read's missing observation, but did
not close the gap. The full gap panel then showed new concordance regressions
at 6.578, 12.269, 3.573, and 21.736 Mb, so the trial was removed. The 3.573-Mb
focused window returned to its accepted 20-assertion baseline. Broadening
homopolymer REF backfill is not a safe fix for these gaps.

A full 21–22-Mb owning-chunk replay with `--snarl-allele-phasing` lets the
other ALT provide contrast to the deletion row, but the 21,736,539 and
21,757,517 boundaries still have different phase sets. The saved trial is
`/tmp/pgphase-multiallelic-21736/`; the default rule was retained.

## Follow-up: 21.435 Mb joined

The 21,435,750–21,456,826 pair now joins in its 21–22-Mb owning chunk.
Recovery retains both boundary Q40 BAM SNP calls and the bridge read's
MSA-verified deletion at 21,433,911. The two adjacent BAM source paths have
no weak cut, and the third site corroborates the bridge read's left-block
orientation. The transfer now carries SNP base quality; prior transfer kept
only BAM alleles and MAPQ. The 62.623-Mb single-SNP false join remains split.
Full chr20 changes only three PS labels, with 228,988/236,821 truth-correct
reads and 7,833 discordant reads in both baseline and updated runs. VCF block
count improves from 381 to 380; N50 remains 517,052 bp.
