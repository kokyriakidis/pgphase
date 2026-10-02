# Exact MSA insertion contrast at a deletion row

## Reproduced defect

A recovered deletion row discarded a read matching the other verified MSA
haplotype when that haplotype inserted bases inside the deletion footprint.
The classifier accepted literal reference, the target deletion, and a verified
same-length substitution, but not the longer insertion sequence. For this
separate deletion row the insertion is absence of the selected deletion ALT.
It must not be inferred from cluster membership or relabeled as literal REF.

The synthetic reproducer uses `ACGTGCA` with one complete consensus deleting
`T` and the other inserting `T` at that footprint. The read's footprint is
`TT`. The accepted baseline returns unknown and omits the read from the
row's counts. The initial reproducer fails three assertions before correction.

## Retained behavior

Accept the exact insertion contrast only when both complete consensuses cover
the same reference footprint, exactly one deletes it, and the other has a
longer query sequence. Both independently composed read alignments must agree,
including supported flanks. Keep the candidate rows separate and update their
own observations and counts. A different deletion length stays under its
existing handling. A mixed insertion/deletion context cannot enable the
one-error local fallback: exact supported sequence is required.

The implementation has no truth, competitor, contig, coordinate, or read-name
condition and adds no realignment stage. The ordinary BAM port's cluster-based
observations remain unchanged. Existing source transfer, gauge, and stitch
certificates remain in force.

The new regression checks exact opposite-allele calls, missing or contradictory
consensuses, incomplete and ambiguous flanks, different deletion lengths,
separate insertion/deletion rows and counts, and rejection of fuzzy mixed edits.
It contains 37 assertions. The fuzzy-call case fails two assertions before its
guard; the different-deletion-length scope check fails one before narrowing.

## Targeted 15.351 Mb investigation

The original left noisy region has 65 physically covering reads but only eight
in its MSA clusters. The eight HiPhase boundary-paired reads have no callable
left/right pair in the original source matrix. A temporary focused MSA retry
admits its 57 unassigned reads and recovers deletion calls on all 65 reads.
Its eight paired calls match the diagnostic HiPhase vector: six cross versus
two same. They are not unanimous. The source solve then collapses the deletion
to homozygous and does not provide the complete certified path across the graph
seam. Observation recovery alone does not establish a whole-block join.

Temporary mixed-event retry admission, interior-seam replacement and singleton
context changes are removed. They did not close this target. The original
source ownership, detector, flank context and stitching implementation are kept.

## Rejected broader classifier

Allowing any different-length nonempty opposite consensus also changed partial
versus complete deletion rows. Across chr20 it added 91 phased reads, with
76 additional truth-correct reads and 15 additional discordant reads; accuracy
fell from 96.971755% to 96.966592%. It newly spanned the 15.095642-15.101261 Mb
control in its owning chunk. Its 15-Mb replay added 38 correct reads without
new errors and passed the disjoint parental-flank check at 99.513664% overall
concordance. However, the 34-Mb replay added 43 scored reads and raised errors
from 34 to 48. Removing the fuzzy fallback alone did not remove that cohort.
The broader classifier is not retained. Its apparent new closure is not added
to the accepted panel or claimed as progress.

The final insertion-only scope restores the accepted 15- and 34-Mb replay
outputs exactly, including read HP/PS and VCF keys. No panel expectations or
required sites are refreshed. Final chromosome measurements and validation
follow below.

## Final matched chr20

Binary SHA256:
`d1e05dfe2252cf89c422dfd87db06f240eda6d7f859bc8a300a840008cfa1dfb`.
The eight-thread default run takes 224.46 s concurrent with the native panel;
this is a functional comparison, not an isolated timing benchmark.

| Measurement | Accepted baseline | Retained correction |
|---|---:|---:|
| Primary reads | 256,586 | 256,586 |
| Phased / truth-scored | 237,101 | 237,101 |
| Correct / discordant | 229,921 / 7,180 | 229,921 / 7,180 |
| Read concordance | 96.971755% | 96.971755% |
| Read phase sets | 684 | 684 |
| VCF rows / blocks | 62,798 / 339 | 62,798 / 339 |
| VCF span N50 | 739,888 bp | 739,888 bp |
| Spanned panel windows | 75 / 93 | 75 / 93 |

Read names, HP/PS assignments, VCF keys, genotypes and phase-set labels all match
exactly. **No new tracked gap closes. Fourteen HiPhase target gaps remain open.**
The representation bug is reproduced and corrected, but its affected reads are
not admitted by the current default retry path in this fixture. The temporary
broader retry still did not yield a certified source join. It is not retained.

Standalone units pass; predicate tests pass 636 assertions / 29 cases;
port parity passes 27 / 9; original upstream phasing parity passes 168,696 / 7.
All 105 distinct final native commands run fresh with four workers in 276.94 s.
Exact-command and binary-SHA verified scoring passes 4,082 assertions / 53 cases
on the unchanged 93-window panel, including owning-chunk and parental-orientation
controls. No compiler warnings or whitespace errors are added.
