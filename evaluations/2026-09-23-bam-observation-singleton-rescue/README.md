# Exact BAM-observation singleton rescue

Date: 2026-09-23

## Question

Why does HiPhase correctly tag reads that the final graph+BAM pgphase output
leaves unphased, and can pgphase recover a truth-supported subset without
changing existing HP assignments or variant phase sets?

## Audit

Read names were compared exactly on the complete 272,016-read chr20 fixture.
Every phase set was independently oriented against
`test_data/derived/chr20_truth_hap.tsv`.

HiPhase tags 6,617 reads absent from the prior pgphase output. Of these, 5,009
agree with parental truth; 3,770 are outside chr20:26-30 Mb. Every one of the
3,770 overlaps a phased HiPhase heterozygote, and 2,871 overlap exactly one.
The evidence is predominantly indel based: 3,204 overlap only HiPhase indels,
320 only SNPs, and 246 both.

A one-base representation-aware coordinate audit classified the 3,770 reads:

| cause | reads |
|---|---:|
| corresponding pgphase `REP_HET_INDEL` excluded from the clean solve | 1,807 |
| no pgphase candidate within 1 bp | 1,145 |
| corresponding pgphase candidate phased, but read has no eligible graph observation | 810 |
| clean SNP present but tied/unassigned | 8 |

The graph read diagnostic independently found 3,364 reads with observations but
no eligible vote, 223 without a graph profile, 147 without a called graph
allele, 33 tied, and 3 with a positive margin but no final assignment.

## Retained design

The existing whole-chunk BAM sub-solve now retains allele observations at
sequence-identical biallelic graph candidates. The selected graph ALT is
normalized with `vcf_to_variant_key`, so anchored VCF indels match the internal
BAM representation. Matching is precomputed once per BAM candidate; read loops
use an integer candidate map.

Read rescue remains post-stitch and read-only:

1. Complete the original graph-observation fixed point.
2. Run a second fixed point that may fill a missing graph allele from the exact
   BAM observation.
3. Orient an excluded site only when both haplotypes and alleles occur in
   exactly one graph phase set and the exact one-sided binomial test is
   `p <= 0.01`.
4. Permit one inferred site to tag a read only when primary graph assignments
   alone also have a one-sided 95% Wilson discordance upper bound at or below
   15%. Read-only rescued assignments cannot validate this condition.
5. Preserve called graph alleles on BAM disagreement. The BAM channel fills
   missing alleles only.
6. Preserve primary graph tags, earlier graph-only rescues, and staged
   margin-6 BAM assignments. A BAM-observation rescue fills only an empty
   output assignment.

Candidates, genotypes, phase sets, and stitching votes are unchanged.

## Experiments

| arm | phased | correct | discordant | accuracy |
|---|---:|---:|---:|---:|
| prior 1 Mb / margin-6 baseline | 235,835 | 227,293 | 8,542 | 96.377976% |
| observation attachment only | 235,828 | 227,285 | 8,543 | 96.377445% |
| singleton Wilson upper <=10% | 236,171 | 227,603 | 8,568 | 96.372120% |
| **singleton Wilson upper <=15%** | **236,338** | **227,754** | **8,584** | **96.367914%** |
| HiPhase | 233,353 | 223,800 | 9,553 | 95.906202% |

The retained 15% arm adds 503 reads over baseline: 465 correct and 38
discordant (92.45% truth accuracy). It loses no read, changes no existing HP,
and changes 22 existing PS labels from an independent BAM namespace to a
graph-oriented rescue namespace while preserving HP. The phased VCF is byte
identical.

It recovers 464 of the 5,009 truth-correct HiPhase-only reads, including 462
outside chr20:26-30 Mb. The noncentromeric recovered set contains 403 reads at
excluded repeat indels, 31 at represented phased candidates, 26 without a
candidate within the audit's one-base proxy, and 2 clean-SNP cases.

Against HiPhase, the retained pgphase output phases 2,985 more reads, produces
3,954 more correct assignments, produces 969 fewer discordant assignments, and
remains 0.4617 percentage points more accurate. Coverage is 86.8839% and
correct yield is 83.7282%.

## Verification

- `make -j16 pgphase`
- `make unit-tests`
- `make window-tests`: 232 assertions passed
- optimized full chr20 run: 111.04 seconds, 29,859,688 KiB peak RSS
- optimized candidate-index implementation produced byte-identical phased VCF
  and read tags to the selected 15% experiment
