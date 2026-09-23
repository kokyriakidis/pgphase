# Chr20 unphased-read audit

Date: 2026-09-23

## Inputs

The current graph+BAM recovery output was compared by qname with the saved
HiPhase 1.6.0 output. Both use the HG002 chr20 read population. HiPhase used the
DeepVariant VCF recorded in its BAM header; pgphase used the graph catalog plus
targeted BAM recovery. Parental truth is
`test_data/derived/chr20_truth_hap.tsv`. Phase-set orientation was selected by
the majority truth relation within each block before scoring reads.

The pgphase phased BAM contains 252,292 graph-associated reads. The original
linear BAM contains 272,016 reads, so the comparison of pgphase's unphased reads
uses the 252,292 shared qnames. HiPhase tags 3,556 additional reads among the
19,724 that are absent from the graph output.

## Why reads remain unphased

Before the output fix, pgphase left 27,287 shared reads unphased. HiPhase tagged
9,183 of those, while both tools left 18,104 unphased. The graph
`--phase-reads-out` diagnostic classified pgphase's misses as:

| Cause | All pgphase-unphased | HiPhase tagged |
|---|---:|---:|
| Graph alleles observed, but no eligible site produced a haplotype score | 24,795 | 8,127 |
| No graph allele observation | 1,807 | 396 |
| Eligible evidence tied | 174 | 156 |
| Positive score margin but solver left the read unassigned | 14 | 9 |
| Internally assigned, then erased by an unphased downstream chunk | 497 | 495 |

A coordinate-level catalog audit found 19,544/27,287 reads overlap no phased
pgphase heterozygote. HiPhase phases 4,479 of those using its larger variant
set. HiPhase emits 77,123 phased heterozygotes versus pgphase's 61,644; 5,592 of
its 9,183 extra assignments rely on exactly one HiPhase site. The extra subset
is difficult: 8,071/9,183 (87.89%) HiPhase assignments agree with parental
truth, compared with 95.91% across all HiPhase-tagged reads. Only 505/9,183 have
BAM MAPQ below 30, and none is below HiPhase's MAPQ 5 floor, so MAPQ is not the
main deficit.

## Fixed output loss

A read seen in adjacent chunks used the last chunk's assignment unconditionally.
An unphased downstream visit therefore erased a valid upstream HP/PS tag. The
merge now lets a later phased assignment replace an earlier one but preserves a
valid assignment across a later unphased visit.

This adds 512 reads with no lost or changed existing assignments. Parental truth
supports 500/512 (97.66%). The phased VCF is byte-identical, so the fix changes
read output only. The current output has 225,517 phased and 26,775 unphased
shared reads. HiPhase phases 8,682 of the remainder; 7,612/8,682 (87.68%) agree
with parental truth.


## Post-solve excluded-site rescue

After final chunk stitching, phased reads now orient biallelic candidates that
were retained in the chunk but excluded from the clean graph solve. An exact
one-sided binomial test at p<=0.01 is required, evidence from different phase
sets is never pooled, SNP votes take precedence over indel votes, and
co-located conflicting rows abstain. One directly phased site may tag a read;
an indirectly oriented excluded site requires confirmation from a second
independent locus. Accepted assignments live in the existing read-only fallback
vectors under `PS + kGapFillPsOffset`; they do not change candidate GT/PS, join
phase sets, or override a primary read assignment. Fixed-point layers let a
supported site extend the independent read block toward the middle of an
observation gap.

On the shared chr20 population this adds 4,555 tags, loses no tags, changes no
existing HP/PS assignment, and leaves the phased VCF byte-identical. The output
now phases 230,072/252,292 reads, 275 more than HiPhase's 229,797 shared-read
tags. It recovers 3,159 of the particular 8,682 HiPhase-only baseline reads;
2,703/3,159 (85.57%) agree with parental truth, close to HiPhase's 87.68% on the
full target set. Across all added reads, 3,555/4,555 (78.05%) agree with truth;
whole-output accuracy is 222,801/230,072 (96.84%), compared with HiPhase's
223,800/233,353 (95.91%) on its full BAM population.

The graph window regression moves five read-concordance floors to their measured
values (0.88, 0.98, 0.99, 0.85 and 0.94). Candidate phasing and the phased VCF
are unchanged; these rows record the intentional extra read assignments. All
232 window assertions pass.

A controlled HiPhase 1.7.0 run on pgphase's own emitted graph VCF tagged only
225,978 shared reads and recovered 3,444/8,682 of the target reads. This shows that swapping in its A* solver alone cannot close the remaining
coverage gap. HiPhase's original
DeepVariant result contains 25,465 phased alleles absent from pgphase's phased
VCF, including 11,312 SNPs. Only 6,344 of those alleles, including 1,818 SNPs,
exist exactly in the full graph catalog. The remaining exact-read difference
therefore requires sample-private site discovery or representation recovery;
it cannot be recovered by changing the graph phaser alone.
