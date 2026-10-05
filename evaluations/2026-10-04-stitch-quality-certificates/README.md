# Keep physical stitch certificates attached to their exact allele

## Retained implementation fixes

1. `complete_recovery_block_flip` used the working graph allele but read its
   quality from the BAM channel without checking that the BAM allele was
   present and equal. A Q40 belonging to the opposite BAM allele could therefore
   certify the graph observation. Require equality of the two called alleles
   before admitting that physical certificate; ordinary graph evidence remains
   available without the certificate.
2. `corroborated_bam_block_flip` accepted BAM MAPQ 255 as high confidence.
   That value means unknown mapping quality and now abstains.
3. The same predicate accepted base quality 255 as Q255. Missing base quality
   now abstains, consistent with the other physical stitch predicates.

The two graph-channel regressions fail on the accepted pre-change predicate
(`certificate-before.log`). The missing-base-quality and unknown-MAPQ cases
fail after the first fix but before the sentinel fixes
(`missing-quality-before.log`). Positive stitch fixtures now carry actual BAM
alleles alongside their physical qualities, including when candidate indices
shift. Four added negative assertions guard the bugs. No stitch threshold,
source genotype, allele representation, or original BAM behavior changes.

## Gap experiments

Work directory: `test_data/tmp_gap_fix50`. Fresh owning-Mb probes cover 3, 21,
35, 57 and 58 Mb. The two 21 Mb source audits find a boundary SNP with no
physical ALT at 21,172,501 and one paired MSA-indel molecule at
21,594,343–21,612,458. They do not establish additional missing transferred
verified calls. Protected wrong-join controls remain split.

At 35,498,368–35,516,845, test the existing exact MEC kernel with both full
flank memberships fixed and every usable binary graph/BAM gap row as a free
variable. With five additional rows, uniform scores are 156 for the wrong
parental connection and 159 for the correct connection. SNP-priority scores
are 56,845 versus 56,852; one deterministic read half is tied. With centered
extra sites, the complete matrix still prefers the wrong connection. Truth is
used only after the solve to score its choice. `joint-matrix-35.json` retains
all candidates and scores. For the 57,854,341–57,866,713 case, centered extras
are disconnected; all extras tie under uniform weights and narrowly prefer
wrong polarity with SNP priority (`joint-matrix-57.json`). No forced MEC join
is retained.

Fresh HiPhase 1.6.0 default runs use the same owning-probe pgphase VCF sites
and complete original BAM alignments overlapping 35.445–35.560 Mb and
57.810–57.900 Mb. Genotypes are unphased before input. The named read segments
in serial trace output preserve each allele vector's explicit read identity
and site interval; do not zip anonymous global vectors with BAM fetch order.
Global alignment can skip records, so that proposed mapping fails the count
check. `audit_hiphase.py` checks the VCF/block variant count, vector extent and
molecule presence before recording boundary calls.

HiPhase has three MAPQ60 molecules calling left deletion REF and right SNP ALT
at 35.5 Mb. All three left calls are unknown in pgphase's initial BAM source,
before transfer. A trial extends deferred fixed-consensus recall to solitary
MSA deletions inside recovery seams and admits the necessary unplaced reads
at consensus collection. It keeps genotypes and existing HP/PS fixed during
observation admission. The target deletion gains 31 calls (9/5 REF/ALT becomes
32/13), restoring all three missing bridge calls. The final gap remains split.
The owning-Mb result adds 124 scored reads but increases discordance from
170 to 219, losing six previously correct assignments (five become unphased).
This broad recall trial is rejected and removed. The partial trial changing
site eligibility alone had no effect because consensus collection still
admitted only insertion contrasts. `deletion-trial-calls.json` and
`rejected-deletion-recall.json` preserve the results.

These targeted HiPhase solves connect both intervals, but do not meet the
accuracy of pgphase's existing independent read groups on these inputs:

| Gap | pgphase scored / correct | HiPhase scored / correct |
| --- | ---: | ---: |
| 35,498,368–35,516,845 | 83 / 80 (96.39%) | 100 / 93 (93.00%) |
| 57,854,341–57,866,713 | 81 / 81 (100%) | 104 / 94 (90.38%) |

Score unique primary molecules physically overlapping the gap, using the
original BAM spans and parental truth; each emitted PS has its arbitrary
orientation chosen independently. The diagnostic uses current pgphase sites,
not the previous DeepVariant competitor input or a full-chromosome HiPhase
run. It does not revise that benchmark. Raw boundary calls and local scoring
are in `hiphase-boundary-calls.json` and `local-read-truth.json`.

## Default chr20 validation

Accepted before: `test_data/tmp_gap_fix49/full-final` (SHA
`54e67e0911beb6402e2558b44ac5dbf31d04da8da8f26b93a66e9acd8a91348f`).
Retained after: `test_data/tmp_gap_fix50/full-final` (SHA
`9bda79a94b69241f8b9083ee902dceea6d0c4560f583ff09cfce2543afcb45df`).
Same reference, annotated BAM, catalog, GAF, defaults and eight threads.
Every unique primary HP/PS pair and every nonheader VCF row is identical.

| Metric | Before and after |
| --- | ---: |
| Unique primary output reads | 256,610 |
| Truth-scored phased reads | 237,328 |
| Correct / discordant | 230,232 / 7,096 |
| Read concordance | 97.010045% |
| Scored read phase sets | 654 |
| VCF keys / blocks | 64,188 / 331 |
| VCF span N50 | 806,449 bp |
| Tracked coordinate spans | 91 / 104 |

Ten competitor nominations and three split controls remain open. No new gap
is claimed or expectation changed. The eight-thread run takes 323.69 s while
other evaluation/build jobs run; this is not an isolated runtime benchmark.

Build passes with only the existing vendored abPOA unused SIMD helper warning
on rebuilding `align.cpp`. Standalone unit tests pass, predicates pass 1,536
assertions in 47 cases, and targeted native window tests pass 164 assertions
in seven cases, including source admission, private-block transfer, late
fallback preservation, quality provenance, full-block orientation, complex
left-flank closure and rejected retry scheduling. HiFi/ONT TSV/VCF goldens
and HiFi one-/four-thread determinism pass. The full native panel suite is not
rerun this turn; its spans are independently checked on the new full chr20
output, and its preceding 8,060-assertion run is recorded in the prior audit.
