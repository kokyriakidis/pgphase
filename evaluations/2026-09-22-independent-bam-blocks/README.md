# Independent BAM recovery blocks

Date: 2026-09-22

## Question

Does preserving every targeted BAM sub-solve phase set as an independent block
reduce wrong graph-recovery joins, and can a final direct-evidence stitch recover
safe continuity?

## Compared implementations

1. **Prior recovery:** generic strongest-site, aggregate, DP, and MEC paths may
   absorb imported blocks into one chain.
2. **Rejected aggregate:** imported blocks begin independently, but graph/BAM
   and BAM/BAM pairs can stitch from aggregate allele evidence.
3. **Retained:** graph/BAM attachment requires a source-specific 2x2
   read-haplotype vote. Adjacent BAM/BAM blocks may stitch from their exact
   candidate memberships. Both tests use a one-sided exact binomial parity test
   at `p <= 0.01`.

## Full chr20 command

```bash
./pgphase collect-graph-variation \
  --ref test_data/chm13v2.0.chr20.renamed.fa \
  --bam test_data/HG002_chr20_hifi_mapped_to_CHM13_chr20_annotated.bam \
  --sites test_data/chr20.sites.striped.vcf.gz \
  --gaf test_data/HG002.chr20.annotated.coord.gaf.gz \
  -r 'CHM13#0#chr20' \
  --chunk-size 500000 \
  --threads 16 \
  -o candidates.tsv \
  --phased-vcf-out native.vcf \
  --phased-bam-out phased.bam
```

The fast integration panel was emitted without changing committed expectations:

```bash
PGPHASE_EMIT_EXPECTATIONS=/tmp/gauge-only-expect.tsv \
PGPHASE_REQUIRED=/tmp/gauge-only-required.tsv \
PGPHASE_TEST_WORKDIR=/tmp/pgphase-gauge-only \
./test_gap_windows '[gap][windows]'
```

## Results

| metric | retained | rejected aggregate | prior recovery |
|---|---:|---:|---:|
| phased heterozygotes | 59,892 | 59,891 | 59,906 |
| phase sets | 657 | 338 | 187 |
| block span | 54.32 Mb | 54.67 Mb | 56.71 Mb |
| N50 | 412.1 kb | 490.0 kb | 1,136.4 kb |
| largest block | 1,415.3 kb | 2,041.6 kb | 3,004.3 kb |
| tagged reads | 219,579 | 219,621 | 221,836 |
| truth-discordant reads | 4,762 | 9,140 | 36,220 |
| truth discordance | 2.169% | 4.162% | 16.327% |
| target gaps spanned | 14/48 | 28/48 | 31/48 |

The aggregate graph/BAM fallback gains only 42 tagged reads and 14 target
spans. It also creates truth-switched joins at 22.98 Mb (51.1% concordance) and
48.23 Mb (68.3%). The source-specific rule removes both. The retained
48-target result has no `SWITCH`.

## Decision

Retain independent imported blocks, source-specific graph/BAM attachment, and
pair-specific BAM/BAM stitching. This phases 98.98% as many reads as the prior
recovery while reducing truth-discordant reads by 86.9%. Unsupported boundaries
remain separate phase sets.
