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

## Three-tool whole-chr20 comparison

The final pgphase result was compared with the frozen upstream longcalld and
HiPhase results. All three consume the same underlying 272,016 BAM records.
pgphase and longcalld use
`test_data/HG002_chr20_hifi_mapped_to_CHM13_chr20_annotated.bam`; HiPhase uses
`shared_calls/chr20/HG002.chr20.normalized.bam`, which `samtools reheader`
created from that exact BAM without changing records, coordinates, CIGARs, or
qnames. HiPhase's DeepVariant VCF was called from that same normalized BAM and
reference. pgphase additionally uses the graph catalog and GAF; longcalld
discovers its own variants from the shared BAM.

The primary comparison counts every BAM qname. A read absent from a tool's
output is unphased. Each emitted phase set is independently oriented to maximize
agreement with `chr20_truth_hap.tsv`, then every tagged read is scored.

| tool | phased reads | coverage | truth correct | truth discordant | truth accuracy | correct yield |
|---|---:|---:|---:|---:|---:|---:|
| pgphase graph+recovery | 230,468 | 84.7259% | 223,190 | 7,278 | 96.8421% | 82.0503% |
| **HiPhase** | **233,353** | **85.7865%** | **223,800** | 9,553 | 95.9062% | **82.2746%** |
| longcalld | 219,090 | 80.5431% | 212,584 | **6,506** | **97.0304%** | 78.1513% |

No tool dominates all read metrics. HiPhase phases 2,885 more reads than
pgphase and produces 610 more correct assignments, but also 2,275 more wrong
assignments and has 0.936-point lower conditional accuracy. Longcalld has the
highest conditional accuracy, 0.188 points above pgphase, while pgphase phases
11,378 more reads and produces 10,606 more correct assignments. pgphase is the
middle operating point: substantially more coverage than longcalld and
substantially fewer errors than HiPhase.

### Improvement targets from this baseline

| objective | delta from current pgphase |
|---|---:|
| match HiPhase coverage | +2,885 phased reads |
| match HiPhase correct yield | +610 truth-correct reads |
| match longcalld accuracy at the current 230,468 phased reads | at most 6,843 discordant reads, 435 fewer |

Future changes should report all three quantities. Coverage-only gains can move
toward HiPhase while adding too many errors, and accuracy-only filtering can
move toward longcalld by abstaining. A strict improvement increases correct
yield without losing existing correct assignments; the longer-term Pareto goal
is to exceed HiPhase's coverage while retaining at least longcalld's conditional
accuracy.

The graph/GAF path contains 252,292 of the 272,016 BAM qnames. Restricting all
three tools to those graph-observed reads is useful for diagnosing the phaser,
but is not the primary end-to-end coverage comparison:

| tool | phased graph-observed reads | coverage | truth accuracy |
|---|---:|---:|---:|
| pgphase graph+recovery | 230,468 | 91.3497% | 96.8421% |
| HiPhase | 229,797 | 91.0837% | 96.4321% |
| longcalld | 217,619 | 86.2568% | 97.2581% |

On the 210,832 reads phased by all three, accuracy is 98.3276% for pgphase,
97.2495% for HiPhase, and 98.3655% for longcalld. This separates assignment
quality from each tool's decision to abstain and leaves pgphase and longcalld
0.038 percentage points apart on the identical phased subset.

The longcalld totals include 315 reads carrying HP with `PS=0`, matching its
existing evaluator's convention and published 219,090-read total. Excluding
those rows gives 218,775 phased reads on the full BAM; the comparison retains
them because longcalld emitted an HP assignment.

VCF structure is reported separately because the tools do not start from the
same variant callset: pgphase uses the graph catalog, HiPhase uses DeepVariant
calls generated from the shared BAM, and longcalld discovers alignment variants
itself.

| tool | phased heterozygotes | VCF phase sets | span N50 |
|---|---:|---:|---:|
| pgphase graph+recovery | 61,644 | 419 | 482,085 bp |
| HiPhase | 77,123 | 196 | 1,005,183 bp |
| longcalld | 83,013 | 431 | 416,121 bp |

These VCF counts measure output quantity and continuity, not variant accuracy;
a variant-truth comparison would be required to rank the three callsets. The
machine-readable primary summary is in `three_tool_chr20.tsv`.
