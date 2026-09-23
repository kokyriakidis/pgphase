# Statistically validated independent BAM read fallback

Date: 2026-09-23

## Question

Can the graph pipeline use additional read assignments from the BAM signal
without importing candidates, changing graph phase sets, or attaching a BAM
block to the wrong graph haplotype?

## Design

With `--bam`, each graph chunk receives one ordinary BAM `process_chunk` solve.
BAM blocks remain independent. Shared reads already assigned by the graph form
one 2x2 BAM-HP by graph-HP table per graph phase set. A table contributes only
when both haplotypes occur on both sides. Every contributing table must reject
a random 50:50 association at exact two-sided `p <= 0.01`. After independently
choosing each graph block's arbitrary HP orientation, the combined disagreement
rate must have a one-sided 95% Wilson upper bound at or below 10%.

Assignments from a passing BAM block are staged separately and applied only to
reads still unphased after cross-chunk stitching and excluded-site rescue. They
use `PS + 1,500,000,000`, separate from graph PS and the excluded-site rescue
namespace. No BAM candidate, allele observation, consensus, or stitching vote
is transferred by this pass.

A broader assignment-only experiment without block validation added 4,452
reads but those additions were only 69.18% truth-correct. A rule accepting a
BAM block when any single graph overlap supported it also admitted blocks with
contradictory overlaps elsewhere. Both were rejected.

## Reproduction

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

The baseline is commit `62f25bc` with its excluded-site rescue enabled.
Read truth is `test_data/derived/chr20_truth_hap.tsv`; each PS is oriented by
its majority parental label before reads are scored.

## Full chr20 result

| metric | baseline | validated BAM fallback | change |
|---|---:|---:|---:|
| tagged shared reads | 230,072 | 230,468 | +396 |
| truth-correct reads | 222,801 | 223,190 | +389 |
| truth-discordant reads | 7,271 | 7,278 | +7 |
| truth accuracy | 96.839685% | 96.842078% | +0.002393 pp |
| lead over HiPhase on shared reads | 275 | 671 | +396 |
| phased heterozygotes | 61,644 | 61,644 | 0 |
| VCF phase sets | 419 | 419 | 0 |

All 230,072 prior HP/PS assignments are byte-for-byte unchanged, no read is
lost, and the phased VCF is byte-identical. Of the 396 added reads, 389 (98.23%)
match parental truth. The pass recovers 393 of the former 5,523 HiPhase-only
reads; 387/393 (98.47%) match truth. The remaining shared HiPhase-only set is
5,130 reads. pgphase now has 5,801 shared reads HiPhase leaves unphased, versus
5,130 in the other direction.

The full run completed in 81.24 seconds at 1,268% CPU and 26.4 GB maximum RSS on
the evaluation host. `make unit-tests` passes, and `make window-tests` passes
all 232 assertions without changing expectations.

## Decision

Retain the statistically validated assignment-only fallback by default whenever
`collect-graph-variation` receives `--bam`. It increases coverage and accuracy
while leaving the graph variant result and every existing read assignment
unchanged.
