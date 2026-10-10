# Union gap phasing on chr20 (HG002 HiFi, CHM13)

`--union-gap-phasing` as implemented in `src/union_phase.cpp` (behaviour:
docs/IMPLEMENTATION.md, "Union gap phasing"). It is compared with the default
arm (graph + seam recovery), graph-only and the pinned competitors. Parental
truth is used only to score.

## Setup

- **Inputs.** The `test_data/` chr20 BAM, GAF, site catalog and reference.
  The command is the one in `iterate.sh`, with `--union-gap-phasing` and
  nothing else (`--min-mapq` defaults to 1 in this mode). 10 threads. The
  binary hash is in `binary.sha256`.
- **Read scores** (`audit_hybrid.py`, `region_reads.py`).
  - Denominator: the 272,016 truth-labelled primary reads. Unphased reads
    count.
  - Each phase set takes the majority parental orientation of its own reads.
  - Arms: alignment midpoint outside [26, 32) Mb, 246,372 reads. Cen: inside
    it, 25,644 reads.
- **Variant switches** (`switch_points.py`): adjacent phased hets against the
  GIAB HG002 CHM13v2.0 v5.0q small-variant truth.
- **Graph links** (`graph_contig.py`, `link_truth.py`): consecutive phased hets
  in one graph-only phase set. "Broken" means the run has no block spanning
  both sites.
- **Competitors.** Pinned scores in `pgphase-eval-data/results/pinned_chr20_scores/`.

## Results

Whole chr20:

| tool | correct | discordant | unphased | switches (arms / cen) |
|---|---:|---:|---:|---:|
| **union** | **232,591** | **6,231** | **33,194** | **38** (7 / 31) |
| legacy hybrid (default `--bam`) | 230,965 | 6,620 | 34,431 | 39 (11 / 28) |
| HiPhase | 223,800 | 9,553 | 38,663 | 106 (16 / 90) |
| longphase | 221,375 | 6,758 | 43,883 | 98 (8 / 90) |
| WhatsHap (optimised) | 215,911 | 8,059 | 48,046 | |
| longcallD | 212,407 | 6,368 | 53,241 | |
| WhatsHap | 211,420 | 12,864 | 47,732 | 143 (33 / 110) |
| graph only | 210,326 | 2,754 | 58,936 | |

Arms only (246,372 reads), which is the current focus:

| tool | correct | discordant | unphased | discordant / phased |
|---|---:|---:|---:|---:|
| **union** | **218,973** | 4,703 | **22,696** | 2.10% |
| legacy hybrid | 217,496 | 5,176 | 23,700 | 2.32% |
| HiPhase | 210,696 | 7,178 | 28,498 | 3.29% |
| WhatsHap (optimised) | 205,228 | 7,564 | 33,580 | 3.55% |
| longphase | 204,666 | **1,527** | 40,179 | **0.74%** |
| longcallD | 203,501 | 5,224 | 37,647 | 2.50% |
| WhatsHap | 200,750 | 12,355 | 33,267 | 5.80% |
| graph only | 198,417 | 2,122 | 45,833 | 1.06% |

**Graph contiguity on the arms.** 1 of the graph-only links is broken, and it
is a link the truth cannot judge. There are 12 relative-phase flips: 11 cannot
be judged and 1 is a graph-only error that union corrects. Legacy breaks 0 and
flips 21 on the arms.

**Gap panel (tier 1, 358 graph gaps).**

- Targets: 16 joined correctly, 2 switched, 62 open.
- hybrid_closed controls: 150 correct, 0 switched, 3 undetermined, 31 open.

**Runtime and determinism.** 90 s wall for full chr20 on 10 threads. Two runs
of the same tree give byte-identical VCFs.

## Where the remaining arm errors and unphased reads are

`unphased_reads.py` and `nohet_reads.py`, with the CHM13 v2.0 SD track
(`pgphase-eval-data/annotations/chm13v2.0_SD.bed`, sha256 alongside).

- **Reads that cannot be phased.** 25,491 arm reads (10.3%) span no truth
  heterozygote. The two haplotypes are identical over the whole read. These
  reads lie in 1,107 stretches covering 16.3 Mb, including a run of
  homozygosity from 43.7 to 45.8 Mb.
- **Union labels some of them anyway.** It labels 6,163 of these reads at 55%
  accuracy: 3,403 correct and 2,760 discordant. That is 59% of union's arm
  discordant reads. Other tools on the same reads:

  | tool | labelled | accuracy |
  |---|---:|---:|
  | legacy | 5,910 | 55% |
  | HiPhase | 1,500 | 54% |
  | longphase | 241 | 90% |

- **The sites behind those labels.** Phased records that sit inside these
  reads are mostly sites the truth does not hold: 691 for union, 109 for
  graph-only. The extra ones are injected clean SNPs (339 against 102) and
  MSA noisy indels (280). A false het in a homozygous stretch splits reads at
  random, and one isolated site cannot reveal that to the EM.
- **Unphased arm reads by cause.** Of 22,696:
  - 19,489 span no truth het;
  - 306 are at least 50% in segmental duplications;
  - 297 have MAPQ below 20 (outside SDs);
  - 2,604 are phaseable but unphased. For 1,400 of these the run has phased
    hets in the span but the read's posterior stayed below 0.9. For 1,204 the
    run emits no het in the span. HiPhase phases 1,727 of the 2,604.
- **Phaseable reads.** Of the 220,881 arm reads spanning at least one truth
  het, union labels 98.5%.

## Configuration choices measured during cleanup

| variant | correct | discordant | unphased | arm links broken |
|---|---:|---:|---:|---:|
| dev best (env toggles) | 232,652 | 6,210 | 33,154 | 29 |
| clean, no window MAPQ filter | 232,662 | 6,184 | 33,170 | |
| clean, window MAPQ ≥ 5 | 232,668 | 6,178 | 33,170 | 29 |
| + catalog blocks locked in the EM | 227,178 | 11,668 | 33,170 | 7 |
| **+ block cuts judged over every read (kept)** | **232,591** | **6,231** | **33,194** | **1** |

- **Two dev pieces left out with no loss.** Demotion of non-voting graph
  repeat rows describing an injected allele, and the union-only linker changes
  in `collect_phase.cpp`.
- **Why links broke in the dev best.** Reads below MAPQ 20 stopped teaching
  the model, so they also stopped supporting block continuity. Locking the
  graph's blocks outright made the EM's repairs into switch errors, so it was
  rejected. Letting every read decide cuts restores arm contiguity, at a cost
  of 77 correct reads chromosome-wide.

## Files

- `iterate.sh`: the one-command benchmark (tier-1 panel, full chr20,
  switches, arms/cen split, pinned rows).
- `region_reads.py`, `graph_contig.py`, `link_truth.py`, `unphased_reads.py`,
  `nohet_reads.py`, `false_het_sites.py`: scorers (evaluation only).
- `regions_union.tsv`, `regions_pinned.tsv`, `clean5.log`: outputs for the final binary
  and the pinned rows.
- `binary.sha256`: hash of the scored binary.
