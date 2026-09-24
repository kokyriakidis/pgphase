# Independent BAM read fallback and graph context

Date: 2026-09-23

## Question

Can the graph pipeline use read assignments supported by sample-private BAM
variation, including reads with no graph-catalog observation, while exceeding
HiPhase read coverage and retaining higher truth accuracy?

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
only when its BAM haplotype score margin is at least 6. One clean biallelic SNP
or indel creates a four-point separation, so the retained threshold requires
additional supporting evidence. Reads with graph profiles fill only rows still
unphased after graph stitching and excluded-site rescue. Reads absent from the
graph profiles are carried in a separate output-only list. A primary graph
assignment from any overlapping chunk always wins.

Both fallback paths use `PS + 1,500,000,000`. They cannot alter graph candidate
GT/PS, join a graph block, or affect cross-chunk stitching.

With `--bam`, the graph command now defaults to 1 Mb chunks. Graph-only and
standalone BAM runs keep their 500 kb default, and an explicit `--chunk-size`
always wins. Graph and independent BAM solves receive the same chunk context.

## Reproduction

```bash
./pgphase collect-graph-variation \
  --ref test_data/chm13v2.0.chr20.renamed.fa \
  --bam test_data/HG002_chr20_hifi_mapped_to_CHM13_chr20_annotated.bam \
  --sites test_data/chr20.sites.striped.vcf.gz \
  --gaf test_data/HG002.chr20.annotated.coord.gaf.gz \
  -r 'CHM13#0#chr20' \
  --threads 16 \
  -o candidates.tsv \
  --phased-vcf-out native.vcf \
  --phased-bam-out phased.bam
```

Read truth is `test_data/derived/chr20_truth_hap.tsv`. Each output phase set is
oriented independently by its majority parental label before scoring.

## Experiments

Accepting every assignment from an unvalidated BAM block is not selective
enough. The 2,076 reads added beyond the margin-4 result were only 57.47%
truth-correct. Most had margin 2; MAPQ, number of observations, and rejected
graph-link summaries did not isolate a large high-accuracy subset.

Of the 6,301 reads HiPhase phased but the broad BAM arm did not, none received
BAM HP/PS in any profiled pgphase chunk. Of these, 2,735 were in the excluded
centromeric interval. The standalone BAM writer MAPQ gate was not responsible:
the recovery pass consumes the in-memory phasing state before that output gate.

Extending only the independent BAM solve by 250 kb on each side was rejected.
It produced 234,708 phased reads, 226,055 correct and 8,653 discordant. The graph
and BAM solves need the same context.

The joint context and confidence sweep was:

| graph chunk | BAM margin | phased | coverage | correct | discordant | accuracy |
|---:|---:|---:|---:|---:|---:|---:|
| 500 kb | 4 | 234,787 | 86.3137% | 226,223 | 8,564 | 96.3524% |
| 500 kb | 6 | 232,749 | 85.5645% | 224,914 | 7,835 | 96.6337% |
| 1 Mb | 4 | 236,918 | 87.0971% | 228,010 | 8,908 | 96.2400% |
| **1 Mb** | **6** | **235,835** | **86.6989%** | **227,293** | **8,542** | **96.3780%** |
| 2 Mb | 4 | 239,129 | 87.9099% | 228,903 | 10,226 | 95.7236% |

The selected 1 Mb, margin-6 arm strictly improves the previous 500 kb,
margin-4 result: 1,048 more phased reads, 1,070 more correct assignments, 22
fewer discordant assignments, and 0.0255 percentage points higher accuracy.
The 2 Mb arm phases more reads but falls below HiPhase accuracy and is rejected.

Changing the graph chunk context changes graph block boundaries. The selected
VCF contains 61,726 phased heterozygotes, 452 phase sets and a 460,310 bp span
N50, compared with 61,644, 419 and 482,085 bp at 500 kb. These structural counts
are reported separately from read accuracy.

## Three-tool whole-chr20 comparison

All three tools consume the same underlying 272,016 annotated `vg giraffe` BAM
records. pgphase and longcalld use the annotated BAM directly. HiPhase uses a
lossless contig-header reheader of that BAM, and its DeepVariant VCF was called
from the same reheadered records and reference. pgphase additionally receives
the graph catalog and GAF; longcalld discovers variants internally.

| tool | phased reads | coverage | truth correct | truth discordant | truth accuracy | correct yield |
|---|---:|---:|---:|---:|---:|---:|
| **pgphase graph+recovery** | **235,835** | **86.6989%** | **227,293** | **8,542** | **96.3780%** | **83.5587%** |
| HiPhase | 233,353 | 85.7865% | 223,800 | 9,553 | 95.9062% | 82.2746% |
| longcalld | 219,090 | 80.5431% | 212,584 | **6,506** | **97.0304%** | 78.1513% |

pgphase phases 2,482 more reads than HiPhase, produces 3,493 more correct
assignments, produces 1,011 fewer discordant assignments, and has 0.4718
percentage points higher conditional accuracy. It dominates HiPhase on all four
read-level metrics in this fixture.

Longcalld remains the higher-accuracy, lower-coverage operating point. Its
conditional accuracy is 0.6525 percentage points above pgphase, while pgphase
phases 16,745 more reads and produces 14,709 more correct assignments.

The 252,292-qname graph-observed subset is a phaser diagnostic rather than the
primary end-to-end comparison:

| tool | phased graph-observed reads | coverage | truth accuracy |
|---|---:|---:|---:|
| pgphase graph+recovery | 232,964 | 92.3396% | 96.7802% |
| HiPhase | 229,797 | 91.0837% | 96.4321% |
| longcalld | 217,619 | 86.2568% | 97.2581% |

The output-only BAM channel phases 2,871 reads absent from graph/GAF profiles.
On the 212,076 reads phased by all three, independently orienting phase sets on
that common subset gives 98.5760% accuracy for pgphase, 97.1779% for HiPhase,
and 98.2709% for longcalld.

The longcalld totals include 315 reads carrying HP with `PS=0`, matching its
existing evaluator convention and published 219,090-read total. Variant counts
are reported separately because the tools start from different callsets:

| tool | phased heterozygotes | VCF phase sets | span N50 |
|---|---:|---:|---:|
| pgphase graph+recovery | 61,726 | 452 | 460,310 bp |
| HiPhase | 77,123 | 196 | 1,005,183 bp |
| longcalld | 83,013 | 431 | 416,121 bp |

These VCF counts measure output quantity and continuity, not variant accuracy.
The machine-readable primary summary is in `three_tool_chr20.tsv`.

## Decision

Use 1 Mb chunks by default for `collect-graph-variation --bam` and require
margin 6 for individual reads from unvalidated BAM blocks. This is the tested setting
that improves both read coverage and read accuracy over the previous pgphase
configuration while retaining the accuracy lead over HiPhase.
