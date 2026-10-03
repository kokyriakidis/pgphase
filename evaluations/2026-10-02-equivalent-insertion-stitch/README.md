# Equivalent insertion stitch on a shared SNP background

## Result

Close the noncentromeric chr20 VCF seam **1,196,894–1,196,967** and join its
core read blocks. Both pgphase and native HiPhase place **78/78** truth-scorable
reads overlapping the seam correctly in one block. This change improves block
continuity; it phases no additional reads.

| Full chr20 metric | Before | After |
|---|---:|---:|
| Output read tags | 256,601 | 256,601 |
| Phased / truth-scorable reads | 237,127 | 237,127 |
| Truth-correct reads | 229,974 | 229,974 |
| Discordant reads | 7,153 | 7,153 |
| Read concordance | 96.983473% | 96.983473% |
| Read phase sets | 674 | 673 |
| VCF keys | 63,492 | 63,492 |
| VCF blocks | 338 | 337 |
| VCF span N50 | 756,878 bp | 756,878 bp |

Every read HP is unchanged. Exactly 2,030 core read PS labels change from
1450106 to 1142089; output-only rescue and BAM fallback assignments are
unchanged. Exactly 712 VCF rows change PS, with no changes to keys, GT, depth,
or any other sample field. The prior long-insertion run remains connected.
The full run used eight threads and took 277.71 seconds; the prior run took
281.99 seconds. These runs were not isolated timing benchmarks.

The existing Oct 1 native HiPhase/DV output uses the same primary BAM and
reference. Its local 78/78 result is the comparison here; HiPhase was not
rerun for this change. The saved pgphase-callset HiPhase arm in the wider
audit uses the **Oct 1 callset**, not literally the current VCF.

## Root cause and proof

The BAM insertion and graph insertion both add 155 bases, at different repeat
placements. Comparing their complete edited intervals on genomic FASTA gives
two mismatches, at edited-string offsets 55 and 210. BAM MSA also calls a
common homozygous **1196950 C→T** SNP, with 2 REF and 76 ALT observations.
Applying that common background makes the two insertion edits exactly equal.
The rows stay separate; this is allele-identity evidence, not an allele merge
or another alignment.

The final recovery matrix contains 76 distinct callable, high-MAPQ pairs:
**36 REF/REF, 40 ALT/ALT, zero conflicts**. Both BAM and graph MAPQ must be
known and at least 30. Both allele classes need independent molecules, and
association must pass the existing p<=0.001 test. An ALT/ALT graph genotype
cannot use genomic REF as its binary baseline. Both complete graph SNP paths
are audited, and the graph insertion must independently agree with its next
clean physical SNP. The existing right-path audit can bypass one weak site
only through an independently significant two-haplotype edge; here the
1349327→1351500 bypass has 19/31 haplotype supports and no opposing call.
The exact graph SNP 1197105 T→C retains the same ALT haplotype as both
insertion descriptions.

## Why the merge is deferred

Core stitching does not certify output-only rescued read groups. Merging
before excluded-site rescue rebuilt their groups under a shared parent PS
and pooled eight reads from a different read-only gauge with the downstream
cohort. It also admitted one additional read. The owning chunk changed from
4044/4034/10 phased/correct/discordant to 4045/4027/18.

Record the core join after both recovery solves finalize their candidate
indices, then apply it after cross-chunk stitching, read rescue and independent
BAM assignment. Candidate allele orientation is read again at application time,
so intervening flips are handled. Update the entire downstream core PS across
all batch chunks, including its tail. Keep independently assigned output-only
HP/PS untouched. The accepted owning replay remains **4044/4034/10**.

The trace also found single-locus calls in complex MSA contexts. Blanket
removal of CIGAR backfill lost older supported connections. Site-wide singleton
restrictions lost correct reads and still did not solve all the uncertainty;
one broad trial increased 58 Mb errors. Those trials are rejected. This change
does **not** claim to repair those individual complex-context allele calls.

## Regressions and remaining work

- Add the 73 bp seam to the committed panel: **95 coordinate cases, 81 spans**.
- Pin both insertion descriptions, the adjacent clean SNP, their allele gauge,
  and the parental orientation of all 78 local reads.
- Retain the owning-chunk ceiling of 10 discordant reads and every prior floor.
- Add pure edit-identity tests for common background, rotations, lowercase,
  missing sequence, ambiguous bases and invalid positions; check both flipped
  and unflipped core joins preserve output-only gauges.
- Accepted baseline binary fails the new window's span assertion (12/13 other
  assertions pass); final binary passes it.
- 105 native command replays are freshly generated from the final binary;
  their manifest checks its SHA256 before reusing outputs for Catch2. Separate
  fresh 1 Mb and 57 Mb matrix replays supply diagnostic evidence.
- Window tests: **61 cases / 4,730 assertions**, passing. Predicate tests:
  **37 cases / 822 assertions**, passing. Standalone units pass.

A fresh full-chromosome extent audit finds **316 gaps** and **42 distinct
competitor-supported coordinate targets** (53 competitor records), down from
317/43/54. Exclude 25–30 Mb. Require at least 95% local purity on ten reads,
and at least 95% purity on five reads in each disjoint 10 kb flank with agreeing
parental orientation. These are nominations, not proof of exact allele joins.
The historical panel's ten open targets and four controls are not an exhaustive
list of the current nominations.

Truth, fixture coordinates and competitor output are used only in evaluation
and regression tests. Production decisions use allele strings, candidate
classifications, observed reads and existing confidence/path checks.

## Reproduce

Build with `make -j8 pgphase`, then run:

```bash
./pgphase collect-graph-variation \
  --ref test_data/chm13v2.0.chr20.renamed.fa \
  --bam test_data/HG002_chr20_hifi_mapped_to_CHM13_chr20_annotated.bam \
  --sites test_data/chr20.sites.striped.vcf.gz \
  --gaf test_data/HG002.chr20.annotated.coord.gaf.gz \
  -r CHM13#0#chr20 -t 8 \
  -o /tmp/pgphase-gap-next21/deferred-full/candidates.tsv \
  --phased-vcf-out /tmp/pgphase-gap-next21/deferred-full/phased.vcf \
  --phased-bam-out /tmp/pgphase-gap-next21/deferred-full/phased.bam
```

For the evidence matrix, use the same inputs with `-r
CHM13#0#chr20:1000001-2000000 -t 1 --phase-matrix-dump /tmp/matrix`.
`measure_insertion_stitch.py --help` gives the exact edit/output comparison
arguments; `audit_competitor_gaps.py --help` gives the competitor audit inputs.

Before binary SHA256:
`9e50548ea3690f94d4922571c3bd2958f287cb7f5118bcc2dbc36c9360933d35`.
Final binary SHA256:
`860544d71c534e13ca37e0b23d806b0d62ef58c4a8c6d6bd2995f04cb76381e6`.
