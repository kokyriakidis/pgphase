# Direct clean-SNP support for omitted graph indels (2026-09-24)

## Question and inputs

Could unphased graph indels inside bounded BAM-recovery gaps tag additional
reads without joining phase sets? All runs used the same annotated HG002 chr20
BAM, CHM13 reference, graph catalog and GAF as the existing chr20 evaluation.
Read truth came from `test_data/derived/chr20_truth_hap.tsv` and was used only
after phasing. The A/B baseline is the full-chr20 guarded graph/BAM stitch run
in `/tmp/pgphase-graph-bridge-chr20/` (236,379 phased truth-matched reads).

The initial rescue already infers an excluded site's allele-to-HP relation from
reads assigned to one phase set. A singleton normally needs a one-sided 95%
Wilson upper discordance bound of 15%. That leaves some graph indels unused
even when their alleles directly track a nearby clean SNP. On regional
matrices, the graph deletion at 21,823,068 has 44 reads co-observing the clean
SNP at 21,817,199, with 38 consistent and 6 conflicting allele pairs. Its 22
still-unphased observing reads are all correctly classified by that allele.
The graph deletion at 23,007,537 has 59/59 consistent direct observations and
eight still-unphased observing reads, all truth-correct in the regional replay.
The tempting insertion at 38,294,482 has 22/27 consistent direct observations,
but one independent read half is only 5/8 consistent and is not admitted.

## Trials

A broad trial accepted an excluded graph indel when its nearest phased clean
SNP agreed with the proposed phase-set orientation in the full read set and
both deterministic read halves at one-sided binomial p<=0.05. It added 207
truth-scored reads, 187 correct and 20 discordant; one prior tag was lost and
two existing reads became incorrect. Several false sites had strong allele
links only on a small subset of their observations. At 48,147,226, a graph
deletion linked perfectly on 15 reads but had 61 observations and separated
truth only 34/61. At 7,918,863 the direct link was 33:6, yet only 8/12
remaining unphased reads were classified correctly. Thus the broad rule was
rejected.

The retained rule applies only to one graph repeat indel inside a bounded
recovery seam, with no competing unphased graph alternative at that coordinate.
The closest phased clean SNP is chosen before testing the pair. In addition to
full and half-fold p<=0.05, both alleles and both SNP haplotypes must be
observed, the one-sided 95% Wilson upper bound on pair discordance must be
<=25%, and the upper bound on site observations missing that SNP must be <=50%.
These bounds require evidence that the link is reliable and represents the
site's reads, not just a small selected subset. The indel remains unphased as
a variant; eligible reads receive a read-only assignment in the existing
phase-set gauge. No candidate consensus or phase set is changed.

## Full chr20 result

| Arm | Truth-scored phased reads | Correct | Discordant | Read phase sets |
|---|---:|---:|---:|---:|
| Guarded stitch baseline | 236,379 | 227,800 | 8,579 | 861 |
| Direct-SNP rescue | 236,409 | 227,830 | 8,579 | 862 |

The retained rule adds **30 reads, all truth-correct**, loses no tag, and
changes no shared read's truth correctness. One shared numeric HP and one PS
label change, both truth-neutral. The phased VCF is byte-identical. The run
took 127.12 s wall time with maximum RSS 21,175,884 kB on eight threads.
The 21.823 Mb site accounts for 22 of the 30 new reads in the full run. The
7.919 Mb regional control adds none under the tightened rule.

`make unit-tests`, `make window-tests`, and `make check` gate this result.
