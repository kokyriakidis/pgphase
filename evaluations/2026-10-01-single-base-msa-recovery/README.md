# Shifted single-base MSA insertion recovery (2026-10-01)

## Bug and retained change

A missing MSA call can become REF during recovery's post-solve exact-CIGAR
backfill when the same insertion is aligned at a nearby position. At internal
chr20:23,792,419, the BAM candidate inserts `T`; shifted original CIGAR events
represent that same allele, but the source solve initially sees only three REF
and five ALT observations. Fixing the transferred calls alone leaves BAM HP/PS
labels based on the incomplete source matrix.

The retained correction adds missing, sequence-equivalent **single-base ALT**
calls before the noisy-candidate BAM k-means pass, inside targeted recovery
windows only. It requires MAPQ >=30, known Q30 inserted and flanking bases,
one insertion within 32 bp, no interfering deletion, and exact reference-edit
equivalence. It preserves existing MSA calls, coordinates, independent rows,
and ordinary BAM behavior. Newly extended profiles are indexed before phasing.
Truth is used only below for evaluation; no coordinate, parental label, or
competitor result enters the repair.

Broader trials were rejected. Rephasing after all SNP/indel backfill created
large unrelated switches in the 3–4 Mb chunk. Restoring shifted multi-base
insertion calls also joined the protected 41,866,917–41,898,323 boundary.
Restricting repair to high-quality single-base events preserves that split and
the established complementary deletion connection at 41.900 Mb.

## Matched full chr20

Both runs use the annotated HG002 chr20 BAM, stripped graph catalog, coordinate
GAF, and CHM13 reference in `test_data/`, with default recovery and eight worker
threads. Baseline is commit `74bc786`; output directories are
`/tmp/pgphase-gap-next/baseline-full/` and `fixed-full/`. Each read phase set is
oriented independently against the existing parental read truth map.

| Measure | Before | After |
| --- | ---: | ---: |
| Output reads | 256,570 | 256,570 |
| Truth-scored phased reads | 236,878 | 237,070 |
| Phased-read fraction of output reads | 92.3249% | 92.3997% |
| Correct assignments | 229,592 | 229,890 |
| Discordant assignments | 7,286 | 7,180 |
| Read concordance | 96.9242% | 96.9714% |
| Read phase sets | 689 | 692 |
| VCF variant keys | 62,361 | 62,368 |
| VCF phase blocks | 341 | 341 |
| VCF span N50 | 643,699 bp | 672,998 bp |
| Wall time | 198.61 s | 193.41 s |

All shared variant keys retain their unordered allele genotypes; seven additional
single-base insertion rows become emit-eligible. There are 1,642 shared-key
GT/PS changes, including block gauge flips. Thus this is a source phasing change,
not only an output PS relabeling.

Four previously separated exact VCF boundaries now share a phase set:

| Left | Right | Primary BAM spanners | Local concordance before → after | Correctly separated fraction before → after |
| ---: | ---: | ---: | ---: | ---: |
| 14,446,295 | 14,458,292 | 12 | 96.62% → 98.25% | 56.76% → 76.58% |
| 23,792,418 | 23,806,565 | 4 | 83.05% → 92.90% | 22.94% → 53.21% |
| 53,945,086 | 53,947,946 | 43 | 97.37% → 98.18% | 37.50% → 26.56% |
| 53,947,946 | 53,962,801 | 9 | 97.61% → 98.36% | 59.66% → 60.50% |

Local scoring uses input-BAM primary read spans from 50 kb before the left
boundary through 50 kb after the right boundary, and the full-run read tags.
The 23.792 Mb joined block's two flank orientation votes are 111:13 and 185:0:
the majorities agree, while the left flank still contains read errors. The other
three connections have at least 98% concordant flanks in the joined block.
The first 53.945 Mb gap loses some correctly tagged interior reads despite
improved continuity and local concordance; the coverage cost is recorded here.

The original short 23.792 Mb replay now spans at 382/457 = 83.59% read
concordance, versus its committed 82% floor, and correctly separates 58/109
reads in one block. Its expectation is raised to `spans=1`, concordance >=83%,
and separation >=53%. The owning 23–24 Mb regression also requires both right
complementary rows, their relative allele orientation, and parental checks.

## Reproduction

`commands.sh` runs both supplied binaries over the same full chromosome.
The scoring scripts and detailed temporary matrices are under
`/tmp/pgphase-gap-next/`; small measured summaries are committed beside this
report. Regenerate the parental truth map with
`scripts/make_truth_hap_map.sh` if needed.

## Regression coverage and context

The permanent panel gains the 14.446 Mb case and both 53.945 Mb cases;
23.792 Mb changes its exact span expectation from zero to one. All existing
floors outside that intended closure remain unchanged. New full-chunk panel
floors are >=99% concordance, with measured separation floors of 76%, 26%,
and 60%. Complementary allele rows and recovered single-base insertions remain
required where the chain uses them.

The 15.351 Mb short replay regressed to 96.42% concordance when deprived of
its full source phase-set context. The previous and repaired owning 15–16 Mb
chunks instead score 4,257/4,279 and 4,259/4,280 (99.4859% and 99.5093%),
with identical shared GT/PS fields. That test now uses its owning chunk,
keeping its original 99% concordance, 52% separation, and one-site floors.
The new 14 Mb and 53 Mb controls also use their owning chunks so the regression
checks reproduce production evidence rather than truncated source gauges.

Thirteen previously tracked HiPhase-positive gaps remain open after this change.
All physical primary spanners in that tracked panel already have MAPQ >=30;
lowering the recovery read threshold does not supply a new molecule bridge.
The original graph observations were screened separately from merged BAM calls:
none of the remaining gaps has a conflict-free, two-haplotype graph two-hop
path under the existing proof rule. Complex shifted indels remain a distinct
source-consistency problem, as the rejected broad trials demonstrate.

## Validation

- `make -j8 pgphase`: passes with no new warnings.
- `make unit-tests predicate-tests parity-tests upstream-parity-tests`: passes;
  60 new shifted-insertion assertions and 168,696 upstream assertions pass.
- `make window-tests`: passes, 3,868 assertions in 48 Catch2 cases over the
  expanded 89-window panel, including the protected 41 Mb split.
- `git diff --check` and reproduction-script syntax: pass.
