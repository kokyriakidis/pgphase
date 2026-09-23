# HiPhase-style trusted-site recovery experiment

Date: 2026-09-22

## Question

Can graph recovery use the same joint binary MEC model as HiPhase without
realigning reads, using pgphase's existing graph and BAM allele observations?

## HiPhase behavior inspected

HiPhase builds read-by-variant allele costs, jointly chooses diploid site
orientations with bounded A-star search, forms phase blocks from read/variant
connectivity, and haplotags reads against the solved haplotypes. Its supplied
DeepVariant sites also exclude apparent heterozygotes that the graph catalog can
retain at strongly off-center allele fractions.

The implementation here is independent; no HiPhase source was copied. It uses
only allele observations pgphase already holds and performs no realignment.

## Retained model

The production fallback runs only after the existing atomic graph-to-graph
recovery transaction abstains:

1. Use biallelic sites with both alleles observed and `|AF - 0.5| <= 0.12`.
2. Build and solve the SNP-only read-connected component first. Admit centered
   indels only if SNPs cannot connect the two boundary blocks.
3. Minimize the exact diploid MEC objective. SNP cost is lexicographically
   prior to indel cost. Problems above 20 internal variables abstain.
4. Require the full read set and two deterministic FNV-split read halves to
   choose the same unique parity. A tie or split disagreement abstains.
5. For graph/BAM attachment, also require the existing source-specific 2x2
   block gauge and sequence-identical candidate anchor to exist and agree.
6. Keep imported BAM blocks independent. A block can participate in at most one
   trusted fallback attachment per chunk, including overlapping recovery
   windows, which prevents several locally stable edges from hiding an older
   internal polarity change.

The diagnostic phase-matrix dump now records ref depth, alt depth, and allele
fraction so the trusted-site decision can be replayed exactly.

## Experiments rejected while developing the gate

- Uniform-cost all-site MEC chose the wrong parity at 33.79 Mb because noisy
  indels and false heterozygous SNPs outweighed a 30-to-2 direct SNP link.
- Edge-local graph/BAM attachment without whole-block safeguards closed 28/48
  targets but reduced chr20 truth accuracy to 95.31%.
- Allowing a recovered block to attach through several overlapping windows
  produced near-50% merged blocks despite each local edge being stable.
- Allowing trusted BAM/BAM fallback joins changed reads in regions where no
  graph seam closed. Those blocks remain independent in the retained design.

## Final chr20 result

Same 500 kb chunks, 16 workers, graph catalog, GAF, recovery BAM, and parental
truth map in both arms:

| metric | MAPQ 5 baseline | trusted MEC | change |
|---|---:|---:|---:|
| tracked spans | 17/48 | 26/48 | +9 |
| current HiPhase-correct spans | 9/34 | 18/34 | +9 |
| historical spans | 8/14 | 8/14 | 0 |
| phase blocks | 659 | 467 | -192 |
| N50 | 412,113 bp | 482,085 bp | +69,972 bp |
| tagged reads | 225,005 | 225,005 | 0 |
| truth discordant | 6,116 | 6,259 | +143 |
| truth accuracy | 97.282% | 97.218% | -0.064 points |

New tracked closures, with none lost:

- 0,882,277-0,890,261
- 11,495,992-11,497,008
- 14,719,378-14,729,749
- 30,813,841-30,814,005
- 33,792,638-33,797,964
- 46,895,038-46,896,191
- 58,834,248-58,836,037
- 60,971,452-60,980,031
- 61,373,859-61,375,344

The per-block truth audit found no new large polarity reversal. The 3.85 Mb
integration control remains split and moves from 91% to 89.3% concordance; its
expectation floor is 89%. HiPhase is only 54.3% in that window.
