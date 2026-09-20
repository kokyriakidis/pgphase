# Parity report: the alignment path against longcallD

State at `e855265`, measured on whole chr20 (HG002 HiFi, CHM13 chr20) against
upstream's own output in the eval data, `##source=longcallD version=0.0.11-23e369d`.

## Method

Two independent comparisons, because either alone misleads.

1. **Record level.** Match on the full record `(POS, REF, ALT)`, phased records
   only. Matching on position alone silently pairs an insertion with a deletion
   at the same anchor and reports false agreement -- it produced two wrong
   counts earlier in this work.
2. **Rule level.** Read the upstream source for each rule and compare it with
   ours, rather than inferring the rule from output. Genotype comparison is
   deliberately NOT used as a parity metric: a block's orientation is a free
   gauge, so 19,669 of the shared records differ in literal GT with no defect
   behind it.

## Where the alignment path stands

| arm | phased records | identical to upstream | upstream-only | ours-only | comma-ALT | two-row positions | two ALTs on one hap |
|---|---:|---:|---:|---:|---:|---:|---:|
| default, merge on (ships) | 116,896 | 114,407 | 3,863 | 2,489 | 1,667 | 828 | **0** |
| **parity config, merge off** | **117,738** | **116,851** | **1,419** | **887** | **0** | **1,671** | **0** |
| upstream longcallD | 118,270 | -- | -- | -- | 0 | 2,401 | 0 |

In the parity configuration we reproduce **116,851 of upstream's 118,270
phased records exactly -- 98.8%** -- with its representation (no multiallelic
records, split rows at a locus) and its coherence (no haplotype ever carrying
two alleles).

## The cost, which is the real finding

| arm | reads tagged | read blocks | misplaced | read hamming |
|---|---:|---:|---:|---:|
| default, merge on | 216,962 | 349 | 1,516 | **0.699%** |
| parity config, merge off | 217,641 | 390 | 3,014 | **1.385%** |

Representational parity is not free and not cosmetic. The co-located merge runs
BEFORE the solve, so a merged locus enters k-means as one site carrying two
alleles instead of two competing single-allele sites. Turning it off doubles
misplaced reads on this arm and fragments blocks 349 -> 390.

So the two goals are in tension on measured evidence: upstream's
representation, or our read-level accuracy. That is the decision to make, and
it is now priced.

## Fixed this round

| bug | measured effect |
|---|---|
| The VCF writer re-applied the DETECTION thresholds at emission (`passes_vcf_depth_gates`), against counts the MSA and merge had rewritten. Upstream's `make_variants` filters on category and active region only. | +726 records, identical-to-upstream 114,111 -> 114,407, read placement unchanged |
| `drop_conflicting_haplotype_alleles` was wired into the graph writer only, so the alignment arm emitted positions with two ALTs on one haplotype -- 8 with the merge on, 412 with it off. | 412 -> 0, read placement unchanged |
| Dropping was the wrong resolution where both alleles have support: `make_colocated_alleles_complementary` moves the weaker allele to the other haplotype, as upstream's structure does by construction. | recovers records that dropping cost; 882,277 now emits `A>ATC 1\|0` beside `A>ATCTC 0\|1` |

## Ruled out as divergences, with the evidence

- **Every threshold matches.** All seventeen: MAPQ 30, depth 5, alt depth 2, AF
  0.20/0.80, noisy region length 50,000 and coverage 1,000, X/gaps 5, slide
  window 100/25, noisy fraction per read 0.5, clip 30 with flank 100, merge
  distance 500, flank 10, sampling size 10,000, and the minimum reads
  supporting a noisy region -- which upstream leaves commented out and derives
  from `min_alt_dp`, exactly as we do.
- **The allele-fraction rule never touches an MSA candidate.** Classification
  runs before the noisy pass (`collect_var.cpp:2171`); MSA candidates take the
  two-haplotype verdict (`collect_phase_noisy.cpp:535-563`), as upstream does.
- **The merge collision rule is a divergence but inert here.** Upstream always
  keeps the old candidate and frees the MSA's (`collect_var.c:1329-1336`); our
  four `replace_*` paths are additions. Disabling all four leaves this arm at
  the same record count. They matter to the graph arm only.
- **The consensus rules are faithful ports.** The argmax update
  (`collect_phase.cpp:223-236` against `assign_hap.c:244-268`) and the
  complement inference (`:263-268` against `:139-142`) match line for line,
  including the prefer-reference tie-break and the ONT guard.

## Not a bug: the `is_alt_genotype` class

1,149 phased candidates are suppressed because their haplotype consensus is
`[0,0]`. Chased to the read level at `26,591,362`: both tools sample the SAME 8
reads and 7 of 8 labels agree. The eighth decides the record -- a paternal read
we place on HP2 and upstream places on HP1. Our labelling is 8/8 consistent
with parental truth, upstream's is 7/8, and its `0|1` rests on the misplaced
read. Of 13 evaluable positions in the class, 8 are real heterozygotes we lose
and 5 are mixed-parent artifacts upstream should not call.

`is_alt_genotype` itself is correct and stays: upstream emits no `0|0` records
either (chr20 is exactly 41,705 `1|0`, 41,308 `0|1`, 35,257 `1|1`).

## What remains, sized

Of the 1,419 upstream records we lack in the parity configuration:

| bucket | count | what is known |
|---|---:|---|
| positions we never emit, but hold a phased candidate for | 461 | mostly the `[0,0]` consensus class above |
| positions with no candidate at all | **162** | clustered in 12 windows of 10 kb, 146 SNPs; largest 50 at 30.80 Mb and 89 across 26.68-26.71 Mb, all low-mapping-quality (21 of 28 reads at MAPQ 1-19) |
| allele differs at a position we hold | ~796 | 148 indels and 24 SNPs characterised so far |

The 162 are the open candidate defect: not a threshold, not emission, so
digar-level read filtering or how a noisy region is cut in low-MAPQ regions.

## One caveat on the option

`merge_colocated_msa_alleles` has no command-line surface -- it is a source
default only. It exists so the comparison above can be run, not as a
configuration knob, and by the standing instruction it should be deleted rather
than kept once the representation question is settled.
