# Normalize graph SNPs and calibrate second recovery admission

## Retained correction

The guard for a newly exposed seam inside a completed recovery window checked
`candidate.key.alt.size() == 1`. Catalog SNP candidates store a graph walk in
that field, so valid graph SNPs never reached the physical BAM evidence query.

Project the selected original graph REF/ALT with the existing
`selected_graph_candidate_alt` and `vcf_to_variant_key` functions before finding
the nearest flanking SNPs. Use the normalized position and ALT for the BAM
query. Decomposed rows select their actual original ALT. A projected MNP or
indel, whole multiallelic row, or absent mapping is not a physical SNP anchor.
Injected BAM SNP keys stay unchanged.

Normalization alone exposed a second defect: high base-quality odds allowed
small or one-haplotype cohorts to admit another solve. Keep the original MAPQ30,
Q30, four-molecule, 75% majority, and signed quality-odds requirements, and also
require the exact one-sided binomial upper tail against random parity to be
at most 0.001 when either anchor is graph-derived. Both haplotypes must supply
at least two winning-parity calls, and that parity must dominate within each
haplotype. A large cohort cannot hide a reversal on the other haplotype.

BAM-only anchor pairs retain the existing admission rule. Their independent
BAM source solve carries additional phase evidence. Applying the new graph
count screen to those pairs regresses existing protected closures; that trial
is rejected rather than changing the panel's expectations.

The original guard also discarded the orientation that admitted the retry.
Retain that graph SNP relation as an existing `RecoveryPhysicalSnpBridge` in
its group's recovery gauge. The stitch compares this relation with aggregate
and source-gauge evidence; competing validated physical relations must agree.
An indel edge cannot silently reverse the SNP relation that admitted the solve.
Existing source-path and orientation checks still decide whether blocks join.
Runtime decisions use reference alleles and observed reads only. Truth and
competitor outputs do not enter the predicate. No candidate merge, new
realignment, fixture coordinate rule, or modified expectation floor is added.

## Reproduced defects

Move the production retry predicate and its existing physical SNP caller to
`graph_bam_adapter.cpp`, where the adapter tests can exercise an indexed
synthetic BAM and the actual query/normalization code. The existing exact
binomial routine moves unchanged into the shared phasing core. The adapter
and noise test targets now link the existing BAM normalization dependencies;
there is no new source object or alternate implementation.

There are **47 new direct checks**, including fixture creation/indexing. They
cover padded walk SNPs, original ALT2 selection, MNP/indel/multiallelic
exclusion, noisy SNP verification, primary/duplicate/read-name filtering,
missing/low mapping and base qualities, same/swapped phase relations,
nonunanimous decisive support, and weak or one-haplotype/reversing cohorts.
The actual preceding predicate fails **7** checks for graph exclusion and
missing parity output. A normalization-only implementation fails **10** count/
allele/haplotype checks. The original stitch's first-row choice fails **2**
checks for conflicting physical relations, one per row order. Exact
reproduction logs accompany this report.
BAM-only admission and clearing a prior parity output are explicitly tested.

At 61.747 Mb, ten Q30 physical SNP pairs (six REF/REF and four ALT/ALT)
call 61,757,551 and 61,773,799 and agree with the stored allele observations.
Admission succeeds, but the unconstrained primary stitch reversed the right
block. Read concordance in the existing short replay fell to 79.95%. Retaining
the admitted SNP parity joins without reversing the right alleles and restores
**368/379 truth-correct reads (97.0976%)**, passing the unchanged 97% floor.
The existing scored window now also requires both exact SNPs to share their PS
and retain the supported allele relation. The neighboring 61.757-Mb catalog
gap and owning-chunk test remain in the panel.

A new owning-chunk case reuses the existing 15-Mb window replay. It protects
15,056,025–15,071,132 and the adjacent deletion run: new graph admission checks
must preserve the BAM-only connection, both upstream SNPs' gauge, and the
downstream deletion connection. This supplements the existing protected
closure and parental-orientation tests. No panel floor is weakened.

## Rejected experiments

`retry_calls.tsv` records the 14 accepted nominations in the normalization-only
full-chromosome trial, including existing BAM-only nominations. That trial
changes 340 VCF blocks to 330, but phased reads decrease from 237,101 to 237,099
and discordant reads increase from 7,181 to 7,272. It closes none of the 14
remaining tracked targets. It is rejected.

Applying the new count screen to all anchor pairs reduces full-chromosome
errors by two but opens the protected 15.056-Mb gap. In the owning context,
15.056-Mb parental separation falls below its floor, and the 61.747-Mb replay
falls from its 97% concordance floor to 79.95%. Five panel assertions fail.
Keep those gates unchanged. Scope the added count screen to graph-derived
anchors and preserve the established BAM-only path. This restores 15.056 Mb;
the additional retained SNP parity fixes the remaining 61.747-Mb reversal.

## Final production validation

Final default chr20, eight threads, 194.61 seconds:

| Measure | Before | Retained |
|---|---:|---:|
| Phased/scored reads | 237,101 | 237,101 |
| Truth-correct reads | 229,920 | 229,921 |
| Discordant reads | 7,181 | 7,180 |
| Read concordance | 96.971333% | 96.971755% |
| Read phase sets | 686 | 686 |
| VCF blocks | 340 | 340 |
| VCF records | 62,798 | 62,798 |
| Span N50 | 739,888 bp | 739,888 bp |
| Panel spans in chromosome output | 74/92 | 74/92 |

All primary read names and every nonheader VCF row are retained. Exactly one
read changes HP/PS and becomes truth-correct. There is no additional default
chr20 panel closure: 14 tracked competitor targets and four controls remain
open. The short 61.747-Mb replay now retains its supported graph/BAM SNP join
without a wrong haplotype reversal.

All **4,041 assertions / 53 cases** pass across the unchanged 92-window panel
and its owning contexts. Every one of the 105 distinct native input commands
runs on the final production binary (four workers, 249.41 seconds); the stock
Catch runner scores those fresh outputs with exact-command and binary-hash
checks. Shared units, 398 predicate assertions, 27 port-parity assertions, and
168,696 original upstream assertions pass. The build introduces no warnings.
No expectation, required site, accuracy floor or existing span is weakened.

## Reproduction

```sh
make -j8 pgphase
make unit-tests predicate-tests parity-tests upstream-parity-tests
make window-tests

./pgphase collect-graph-variation \
  --ref test_data/chm13v2.0.chr20.renamed.fa \
  --bam test_data/HG002_chr20_hifi_mapped_to_CHM13_chr20_annotated.bam \
  --sites test_data/chr20.sites.striped.vcf.gz \
  --gaf test_data/HG002.chr20.annotated.coord.gaf.gz \
  -t 8 -r 'CHM13#0#chr20' -o /tmp/chr20-candidates.tsv \
  --phased-vcf-out /tmp/chr20-phased.vcf \
  --phased-bam-out /tmp/chr20-phased.bam
```

Final measurements and validation logs accompany this report. Competitors are
not rerun for this predicate correction.
