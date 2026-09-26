# Atomic graph attachment and local BAM source runs (2026-09-24)

## Why the short gaps remained

The 35 tracked gaps are the noncentromeric chr20 gaps under 10 kb that a full
HiPhase replay joined at at least 98% local truth purity with at least 20
truth-scored reads, while the previous pgphase full-chr20 output left them open.
The fixed input list is
`evaluations/2026-09-24-hiphase-on-pgphase-sites/short_gap_replay.tsv`.

Three transfer problems were verified:

- At 11.357-11.360 Mb, a BAM source block moved only the boundary row of an
  existing graph phase set, splitting a graph connection that was present
  before recovery.
- At 14.577-14.584 Mb, the BAM block had 145/145 conflict-free haplotype
  votes and 42 agreeing shared-site votes against the left graph block, but
  source transfer required the *same* block to pass the right graph flank.
  HiPhase joins the left SNP to the BAM SNP and leaves the next right flank
  independent. An earlier graph stitch had also changed the left block's
  phase-set ID, so testing only its original ID missed the attachment.
- At 8.632-8.639 Mb, the BAM source path was sound across the gap, but its
  global path check failed at a later 8.638940-8.662670 Mb cut (2 supporting
  reads from one haplotype, none from the other, one conflict). Rejecting the
  entire source block lost the valid prefix.

The accepted transfer moves an approved graph block atomically into the BAM
source gauge, allows one supported flank, and uses the graph block's current
phase-set label and orientation. When the complete BAM source path fails, it
uses the original BAM solve's weak-cut positions to attach a supported run
to exactly one approved graph block. Sites and reads beyond the cut keep their
source phase set. A source read observing sites from another run does not
inherit the attached graph block's label.

A competing experiment split every BAM source phase set at weak cuts before
transfer. It closed 8.63 Mb but incorrectly joined the truth-opposite 19 Mb
control and broke the valid 55 Mb bridge, so it was removed.

## Matched full chr20 run

Command (substitute the desired output directory):

```bash
./pgphase collect-graph-variation \
  --ref test_data/chm13v2.0.chr20.renamed.fa \
  --bam test_data/HG002_chr20_hifi_mapped_to_CHM13_chr20_annotated.bam \
  --sites test_data/chr20.sites.striped.vcf.gz \
  --gaf test_data/HG002.chr20.annotated.coord.gaf.gz \
  -r 'CHM13#0#chr20' -t 8 \
  -o OUTPUT/candidates.tsv \
  --phased-vcf-out OUTPUT/native.vcf \
  --phased-bam-out OUTPUT/phased.bam
```

The previous output is `/tmp/pgphase-try5538-final`; the new output is
`/tmp/pgphase-hi-gap-source-cuts`. The final two-site run under
`/tmp/pgphase-hi-gap-two-site` has byte-identical VCF and BAM output.
Read truth was scored per output phase set
against `test_data/derived/chr20_truth_hap.tsv`. N50 and block counts come
from `scripts/phase_block_stats.py` on the phased VCF.

| Measure | Previous | New |
|---|---:|---:|
| Truth-evaluable phased reads | 236,527 | 236,675 |
| Truth-correct reads | 228,190 | 228,700 |
| Discordant reads | 8,337 | 7,975 |
| Read truth purity | 96.48% | 96.63% |
| Phased heterozygotes | 61,646 | 61,747 |
| VCF phase blocks | 463 | 411 |
| VCF block N50 | 456,233 bp | 486,139 bp |
| Tracked high-purity HiPhase-only short gaps closed | 0/35 | 27/35 |

The new local-run pass adds 8.632-8.639, 47.118-47.122 and
56.348-56.351 Mb beyond atomic, one-sided whole-block transfer. Their
boundary truth votes are respectively 90/90 and 79/79, 73/73 and 90/90,
and 71/71 and 63/63 in a consistent parental orientation. Across all 27
closed tracked gaps, no flank pair has a confident parental switch. The
19 Mb wrong-join control remains split; the 51 Mb and 55 Mb positive
controls remain joined as intended.

The test panel also gained a 26.030-26.089 Mb span. Its thinnest interval
has 40 physical crossing reads, the two boundary SNPs have 81/81 and 71/71
consistent parental assignments, and its window truth score is 312/313.
That specific expected span and separated-read floor were raised.

Eight tracked gaps remain open: 0.529, 0.542-0.545, 14.236-14.240,
17.617-17.626, 32.490, 34.844, 36.016-36.019, and
37.462-37.467 Mb. Several boundaries are repeat indels or an exceptionally
long structural variant; the transfer fix does not establish a safe allele
path through every such representation. HiPhase still has higher block N50
on the updated identical call set, as measured below.

The exact VCF boundary types of the eight still-open tracked gaps are:

| Gap (chr20, 1-based) | Left boundary | Right boundary |
|---|---|---|
| 528,827-528,828 | SNP | 40 bp insertion |
| 542,052-545,002 | 1 bp deletion | SNP |
| 14,235,594-14,239,930 | 14 bp deletion | SNP |
| 17,616,778-17,625,527 | 1 bp deletion | 1 bp insertion |
| 32,490,058-32,490,150 | 5,555 bp deletion | two-base substitution |
| 34,844,194-34,844,579 | 1 bp insertion | SNP |
| 36,016,138-36,018,966 | 2 bp deletion | SNP |
| 37,461,999-37,466,820 | two overlapping deletion rows | SNP |

Every remaining boundary pair includes an indel or structural variant. At
34.844 Mb, 60 reads in the recovery matrix observe both the injected insertion
and right SNP, but their direct alleles split 35 same to 25 alternate. HiPhase
joins those exact VCF rows with 117/117 locally truth-consistent reads. This
suggests a broader, representation-aware chain is needed; forcing that one
weak pair is not justified. Relaxing the source path rule globally already
failed the 19 Mb safety control.

## HiPhase replay on the updated pgphase VCF

HiPhase `1.6.0-ac3f399` was run on the **new** unmodified pgphase VCF
(`/tmp/pgphase-hi-gap-two-site/native.vcf`), compressed and tabix-indexed
without changing records. It used the same annotated BAM, CHM13 reference,
`--threads 8 --ignore-read-groups`, and default phasing thresholds. The
command was:

```bash
hiphase --bam test_data/HG002_chr20_hifi_mapped_to_CHM13_chr20_annotated.bam \
  --vcf /tmp/hiphase-on-pgphase-source-cuts/input.vcf.gz \
  --output-vcf /tmp/hiphase-on-pgphase-source-cuts/phased.vcf.gz \
  --output-bam /tmp/hiphase-on-pgphase-source-cuts/phased.bam \
  --reference test_data/chm13v2.0.chr20.renamed.fa \
  --threads 8 --ignore-read-groups
```

| Measure | pgphase graph + BAM | HiPhase on pgphase VCF |
|---|---:|---:|
| Truth-evaluable phased reads | 236,675 | 228,362 |
| Truth-correct reads | 228,700 | 216,279 |
| Discordant reads | 7,975 | 12,083 |
| Read truth purity | 96.63% | 94.71% |
| Phased heterozygotes | 61,747 | 59,456 |
| VCF blocks | 411 | 211 |
| VCF block N50 | 486,139 bp | 919,153 bp |

HiPhase still puts the exact boundary variants of all eight open tracked gaps
in one PS on this updated VCF. Their local truth purities are, in table order
above: 98.63%, 98.80%, 98.59%, 99.40%, 100%, 100%, 98.13%, and 99.30%.
Thus the remaining continuity deficit is real; the data here do not justify
joining these indel/SV boundaries under the current pgphase evidence model.

## One-haplotype source path at 34.844 Mb

A later targeted audit found that source PS 34,811,748 had exactly one weak
cut at source SNP 34,818,318. Five assigned BAM molecules consistently
observed that site and the next heterozygote at 34,835,167 on one haplotype;
none conflicted. The original path predicate rejected the cut solely because
the other haplotype had no spanning molecule. The exact one-sided binomial
probability of five consistent votes under random polarity is 1/32 (0.03125).

The source path now admits a conflict-free one-haplotype cut at `p <= 0.05`.
The existing graph/BAM transfer still requires exact shared sites and
conflict-free, two-haplotype read votes at both graph flanks. The 34,844,194
insertion and 34,844,579 SNP now share a phase set. They are both `0|1` in
pgphase and both `1|0` in the HiPhase replay: the relative parental connection
is identical after the arbitrary block-label flip. The dedicated regional
regression reports no parental switch and at least 98% truth concordance. The
19 Mb wrong-join control remains split, and its existing regression passes.

On matched full chr20, this closes one more tracked high-purity HiPhase gap:
28/35 rather than 27/35. VCF blocks decrease from 411 to 410. Phased hets
stay at 61,747. Read truth remains exactly 236,675 phased, 228,700 correct,
and 7,975 discordant (96.63%). Block N50 remains 486,139 bp. The run is
`/tmp/pgphase-onehap-chr20`; its comparator is
`/tmp/pgphase-hi-gap-two-site`.
