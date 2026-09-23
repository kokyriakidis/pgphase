# Independent BAM read fallback

Date: 2026-09-23

## Question

Can the graph pipeline use read assignments supported by sample-private BAM
variation, including reads with no graph-catalog observation, while leaving the
graph phase result unchanged and retaining higher truth accuracy than HiPhase?

## Design

With `--bam`, each graph chunk receives one ordinary BAM `process_chunk` solve.
The imported BAM phase sets remain independent. No BAM candidate, observation,
consensus, or stitch vote is transferred by this pass.

Shared reads already assigned by the graph form one 2x2 BAM-HP by graph-HP
table per graph phase set. A table contributes only when both haplotypes occur
on both sides. Every contributing table must reject random 50:50 association at
exact two-sided `p <= 0.01`. After independently choosing each graph block
orientation, the combined disagreement rate must have a one-sided 95% Wilson
upper bound at or below 10%. Every assignment from a passing BAM block may be
staged.

A BAM block without sufficient or consistent graph overlap remains an
independent output phase set. An individual read from such a block is staged
only when its BAM haplotype score margin is at least 4. A clean biallelic SNP or
indel contributes +2 to one haplotype and -2 to the other, so 4 is the first
complete clean-site separation. Reads with graph profiles fill only rows still
unphased after graph stitching and excluded-site rescue. Reads absent from the
graph profiles are carried in a separate output-only list. A primary graph
assignment from any overlapping chunk always wins.

Both fallback paths use `PS + 1,500,000,000`. They cannot alter graph candidate
GT/PS, join a graph block, or affect cross-chunk stitching.

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

Read truth is `test_data/derived/chr20_truth_hap.tsv`. Each output phase set is
oriented independently by its majority parental label before scoring.

## Threshold experiment

The conservative validated-block implementation at commit `7b98f30` phases
230,468 reads. An unvalidated all-block arm phases 236,863, but narrows the
accuracy lead over HiPhase to 0.101 percentage points. The per-read margin
sweep gives:

| minimum margin | phased reads | coverage | correct | discordant | accuracy |
|---:|---:|---:|---:|---:|---:|
| **4** | **234,787** | **86.3137%** | **226,223** | **8,564** | **96.3524%** |
| 5 | 232,749 | 85.5645% | 224,914 | 7,835 | 96.6337% |
| 6 | 232,749 | 85.5645% | 224,914 | 7,835 | 96.6337% |

Margins 5 and 6 phase 604 fewer reads than HiPhase. Margin 4 is the only tested
score boundary that exceeds HiPhase coverage while retaining higher accuracy.
Compared with the validated-block baseline, it adds 4,319 phased reads. The
selected output phases 2,133 BAM reads absent from the graph/GAF population and
2,186 additional graph-observed reads.

All pre-existing phased HP/PS assignments are unchanged. The selected phased
VCF is byte-identical to the validated-block VCF: 61,644 phased heterozygotes,
419 phase sets, and a 482,085 bp span N50.

## Three-tool whole-chr20 comparison

All three tools consume the same underlying 272,016 annotated `vg giraffe` BAM
records. pgphase and longcalld use the annotated BAM directly. HiPhase uses a
lossless contig-header reheader of that BAM, and its DeepVariant VCF was called
from the same reheadered records and reference. pgphase additionally receives
the graph catalog and GAF; longcalld discovers variants internally.

| tool | phased reads | coverage | truth correct | truth discordant | truth accuracy | correct yield |
|---|---:|---:|---:|---:|---:|---:|
| **pgphase graph+recovery** | **234,787** | **86.3137%** | **226,223** | **8,564** | **96.3524%** | **83.1653%** |
| HiPhase | 233,353 | 85.7865% | 223,800 | 9,553 | 95.9062% | 82.2746% |
| longcalld | 219,090 | 80.5431% | 212,584 | **6,506** | **97.0304%** | 78.1513% |

pgphase phases 1,434 more reads than HiPhase, produces 2,423 more correct
assignments, produces 989 fewer discordant assignments, and has 0.4462
percentage points higher conditional accuracy. It therefore dominates HiPhase
on all four read-level metrics in this fixture.

Longcalld remains the higher-accuracy, lower-coverage operating point. Its
conditional accuracy is 0.6780 percentage points above pgphase, while pgphase
phases 15,697 more reads and produces 13,639 more correct assignments.

The 252,292-qname graph-observed subset remains useful for diagnosing the graph
phaser, but it is not the primary end-to-end comparison:

| tool | phased graph-observed reads | coverage | truth accuracy |
|---|---:|---:|---:|
| pgphase graph+recovery | 232,654 | 92.2162% | 96.6293% |
| HiPhase | 229,797 | 91.0837% | 96.4321% |
| longcalld | 217,619 | 86.2568% | 97.2581% |

On the 212,405 reads phased by all three, independently orienting phase sets on
that common subset gives 98.2755% accuracy for pgphase, 97.1470% for HiPhase,
and 98.2538% for longcalld.

The longcalld totals include 315 reads carrying HP with `PS=0`, matching its
existing evaluator convention and published 219,090-read total. Variant counts
are reported separately because the tools start from different callsets:

| tool | phased heterozygotes | VCF phase sets | span N50 |
|---|---:|---:|---:|
| pgphase graph+recovery | 61,644 | 419 | 482,085 bp |
| HiPhase | 77,123 | 196 | 1,005,183 bp |
| longcalld | 83,013 | 431 | 416,121 bp |

These VCF counts measure output quantity and continuity, not variant accuracy.
The machine-readable primary summary is in `three_tool_chr20.tsv`.

## Decision

Retain the margin-4 independent assignment fallback by default whenever
`collect-graph-variation` receives `--bam`. It is the most accurate tested
configuration that also exceeds HiPhase read coverage, and it does so without
changing the graph VCF or any prior graph read assignment.
