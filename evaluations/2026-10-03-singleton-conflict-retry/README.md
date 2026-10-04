# Isolated MSA deletion conflict beside a singleton graph block

Target: chr20:61,738,239–61,747,506 (9,267 bp), outside the centromere.
The surrounding graph seam is 61,732,321–61,757,551; the right graph phase set
initially contains one site.

## Two defects

The source BAM SNP and homopolymer deletion shared a nominal phase set with
opposite ALT orientations. Only two paired calls remained in the source
matrix, both supporting ALT/ALT. All 30 high-MAPQ physical spanners had been
retained, but ordinary backfill excluded the homopolymer deletion. Missing
calls hid the contradictory source edge from the statistical retry detector.

The existing focused retry also omitted solve padding when its primary
request came from a newly exposed internal conflict. At the singleton graph
flank, the solve ended exactly at the seam. Its strict containment check then
silently skipped the retry.

## Fix

Only seams with a singleton graph flank additionally diagnose an isolated,
phased biallelic MSA homopolymer deletion. Known-MAPQ30/Q30 original CIGAR
calls use the existing exact/edit-equivalent deletion caller. Existing calls,
third deletion lengths, replacement alleles, compound events and co-located
MSA indel alternatives abstain. Ordinary graph block boundaries retain the
complementary-row diagnostic rule.

The temporary calls expose seven ALT/ALT and one REF/REF pair, all opposing
the frozen source relation. Six belong to the source-assigned cohort, enough
for the existing conflict statistic. No threshold or source-read membership
rule is weakened. Profiles and their read index are restored before retry
validation or transfer; physical calls alone never become rescue evidence.

Isolated calls participate only in the internal-conflict statistic, after
the ordinary retry/dropout decisions have been measured. They cannot request
a broad retry. An already supported complementary-row conflict retains its
original selection and context. Only a newly diagnosed isolated conflict adds
chunk-bounded singleton flank context to the existing focused BAM/MSA retry.
Its replacement must preserve old clean-SNP memberships and allele relations
as well as the source keys. Only its accepted, validated matrix is imported;
source path and atomic stitch checks remain in force. Production uses no truth, competitor
output or coordinates and adds no new alignment stage or merged row.

## Rejected broad diagnostic

Admitting isolated deletion diagnostics at every graph boundary also closed
the target, but reopened the protected 17.503 and 22.981 Mb gaps. Full chr20
accuracy fell from 96.995055% to 96.610648%, with 960 previously correct reads
becoming discordant. This broader trial is rejected; its parity and block
audits are retained as counterevidence.

A singleton-only trial still incorrectly let isolated physical calls influence
the generic dropout decision at 17.503 Mb, and padded an existing paired-row
conflict at 19.374 Mb. The latter retained clean-SNP source gauges yet raised
owning-chunk discordance by 341 reads. Those are separate defects: temporary
isolated calls must be confined to the internal-conflict statistic, and newly
admitted requests must not alter established retry selections. The final
routing fixes both, instead of weakening a score threshold. The singleton-only
full trial's 96.852973% accuracy and audits are retained as rejected evidence.

## Permanent regression

The panel grows from 101 to 102 intervals, with an explicit `spans=1`
requirement for this gap and a complete 61–62 Mb owning-chunk replay. The
existing dedicated owning test now requires the closure, equal clean-SNP
allele relations, equal MSA SNP/deletion relations, at least 3,468 scored and
3,447 correct reads, and at most 21 discordant reads. Its previous-binary
control fails six of 34 assertions; the accepted candidate passes all 34.
Existing panel floors, ceilings, required rows and span expectations are
preserved.

New unit tests cover exact/shifted ALT, verified literal REF, third alleles,
compound edits, replacement keys, co-located MSA alternatives, existing
observations, missing qualities, low/unknown MAPQ, unphased rows, seam bounds,
ordinary-flank abstention, idempotence and source genotype/count preservation.
They reproduce the hidden conflict without modifying source read assignments.
The previous paired-only helper fails eight of the original 75 new checks.
Protected owning-chunk comparisons at 17 and 19 Mb now have exactly identical
variant rows and read tags (`protected-owner17-parity.json` and
`protected-owner19-parity.json`). All standalone units and the final 1,404 predicate assertions/46 cases pass.

## Owning-chunk and local read evidence

The owning chunk improves from 3,452 scored / 3,423 correct / 29 discordant
reads to 3,468 / 3,447 / 21. All 16 new assignments are correct, eight previous
errors become correct, and no previously correct read worsens. The same 496
variant keys remain. VCF blocks fall 8→7 and span N50 rises 250,203→510,398 bp.

For the 96 parental-truth reads physically overlapping the gap, pgphase
improves from 83 phased / 79 correct / four discordant to 94 / 94 / zero.
Saved same-BAM HiPhase outputs, with both DeepVariant and pgphase calls, phase
92 / 91 / one. These use whole output-block parental gauges. pgphase still
has two read phase sets: its dominant block contains 85 correct overlapping
reads, versus HiPhase's 91 in one block. The VCF gap is closed; this result
does not claim all local reads are consolidated into one phase set.

## Reproduction

```bash
make -j4
make unit-tests predicate-tests
make window-tests
```

The native full chr20 run uses the default graph recovery, the same annotated
BAM/GAF, reference and graph catalog as the accepted baseline, and eight
threads. Cached regression replays verify the final executable hash, every
normalized argument, input metadata and output hashes; missing requests run
natively. Runtime with concurrent panel validation is not a competitor timing
comparison.

## Full chr20 result

| Metric | Before | Accepted fix |
|---|---:|---:|
| Truth-scored phased reads | 237,209 | 237,225 |
| Truth-correct reads | 230,081 | 230,105 |
| Discordant reads | 7,128 | 7,120 |
| Read accuracy | 96.995055% | 96.998630% |
| Truth-scored read phase sets | 656 | 654 |
| VCF phase blocks | 331 | 330 |
| Block span N50 | 790,093 bp | 806,449 bp |
| NG50 | 643,699 bp | 684,798 bp |
| Tracked connected intervals | 88/102 | 89/102 |

The target is the only newly connected tracked interval and the only merge
of old VCF blocks. Every old phased SNP retains its original allele gauge;
all 64,185 variant keys remain, with no additions or losses. All 16 new
truth-scored assignments are correct, eight previous errors become correct,
and no previously correct read worsens. This is the complete chromosome
output, beyond the shorter owning-chunk regression. Evidence is in
`full-parity.json`, `full-block-audit.json` and `ng50.json`.

The native eight-thread run took 300.26 seconds with panel validation running
concurrently. Final executable SHA256:
`336996e6338b3abc4a55b662231e71e86d378f4222615c4e39aa6cc94a645ce4`.

Final validation passes **7,176 window assertions in 72 cases**, including all
102 panel intervals, the owning context and parental-orientation checks.
All standalone units and 1,404 predicate assertions/46 cases pass. The final
cache contains 106 fresh native normalized requests; uncached requests execute
natively. Main build and `git diff --check` pass, with no new warnings. Details
and negative controls are recorded in `validation.json` and the saved logs.
