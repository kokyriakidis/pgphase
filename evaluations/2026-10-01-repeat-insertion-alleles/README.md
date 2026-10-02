# Repeat-shifted insertion allele audit and correction

## Finding

At the remaining 15,351,845–15,367,755 gap, the right BAM candidate is
`A→AATCT`, internal insertion position 15,367,756. Five of the eight paired
HiPhase reads carry the insertion in their original CIGAR at 15,367,800:
44 bases to the right, across eleven copies of `ATCT`. Four inserted sequences
have Q40 bases; one has Q22. The two placements produce the identical
reference edit. The physical graph stitch caller searches only 16 bases for
insertions longer than two bases, so it calls the high-quality shifted event
REF. The pre-solve source repair only handles single-base insertions within
32 bases; it does not repair this multi-base site.

The retained correction is an ALT-only physical-call fallback. Keep the existing
16-base / 64-base nearby-indel guards and callable ALT decisions. When a
multi-base insertion would otherwise be REF with no nearby indel, search its
complete reference-equivalent placement interval. Require one matching CIGAR
insertion, identical edited reference, and known qualities at the aligned-flank
threshold across every inserted and crossed reference base. Ambiguous reference,
different alleles and compound events cannot supply a new ALT. Use the actual
minimum inserted quality in the stitch likelihood. Keep the original BAM
profiles, HP/PS gauge, weak cuts, candidate keys and stitch gates.

This uses the existing CIGAR and reference, performs no realignment, and has no
truth, competitor, chromosome, read-name or fixture-coordinate condition.

## Rejected experiments

1. Inferring a missing homozygous BAM SNP row from balanced source-HP reads
   demoted the false graph SNP at 15,367,800 but did not close the target.
   A new 52.4 Mb join had worse local read accuracy than the archived HiPhase
   result. The experiment is removed.
2. Broadening physical interference windows over the entire repeat and adding
   missing multi-base ALT calls before the BAM source solve restored 24 right
   `ATCT` observations. Whole-chromosome discordance fell from 7,180 to 7,149,
   but the protected 47 Mb join split and the 17.865 Mb suffix joined the wrong
   source prefix. The 53-case suite failed five assertions. The experiment is
   removed; its chromosome-wide gain does not certify its joins.
3. Adding the multi-base calls in the ordinary postsolve CIGAR backfill still
   changed subsequent retry and source-path decisions. The same protected joins
   failed. This experiment is also removed. Applying a change after one solve
   does not guarantee it preserves the later source certificate.

## HiPhase boundary evidence

A diagnostic build of HiPhase **v1.6.0**, upstream commit
`74dbc7eaa46d2f0678cfafe07a3a22b9a4a8c121`, logged its per-read allele vectors
without changing its algorithm. Inputs are the accepted pgphase VCF rows,
unphased and renamed to the BAM sample, the same annotated primary BAM and
reference. Compact BAMs keep whole alignments overlapping the target context.
This is a targeted allele audit, not a new competitor chromosome benchmark.

| Boundary rows | REF/REF | REF/ALT | ALT/REF | ALT/ALT |
|---|---:|---:|---:|---:|
| 15,351,845 CT→C / 15,367,755 A→AATCT | 1 | 4 | 2 | 1 |
| 21,159,070 C→T / 21,172,487 CA→C | 10 | 1 | 7 | 0 |
| 57,854,341 CAA→C / 57,866,713 C→A | 0 | 1 | 1 | 4 |

Full row identities and read names are in `hiphase-boundary-calls.tsv`. The
57,854,341 locus also contains `CA→C`; the audit identifies `CAA→C` explicitly
and does not collapse the two deletion representations by coordinate.

At 15.351 Mb the direct HiPhase boundary votes are six cross versus two same,
not unanimous. pgphase's left deletion has no callable paired source observations
on these spanning reads; the missing calls at that boundary remain the blocker.
The graph SNP at 15,367,800 G→A is also contradicted by the primary BAM: all
82 aligned bases are G, with 80 at Q30 or above. HiPhase leaves this row
unphased with the same input callset. Correcting only the right insertion does
not restore the left boundary or certify a whole-block orientation.

## Final matched chr20 run

The final binary SHA256 is
`7eb3a9b57570972d2a64332f6790892a33ddfc7139dfc8e06756f84a26db2091`.
The eight-thread run takes 215.24 s while the native panel also runs with four
workers. This is a functional comparison, not an isolated runtime benchmark.

| Measurement | Accepted baseline | Corrected run |
|---|---:|---:|
| Primary reads | 256,586 | 256,586 |
| Phased / truth-scored | 237,101 | 237,101 |
| Correct / discordant | 229,921 / 7,180 | 229,921 / 7,180 |
| Read concordance | 96.971755% | 96.971755% |
| Read phase sets | 684 | 684 |
| VCF records / blocks | 62,798 / 339 | 62,798 / 339 |
| VCF span N50 | 739,888 bp | 739,888 bp |
| Spanned panel windows | 75 / 93 | 75 / 93 |

Every read name and HP/PS assignment, and every VCF key, genotype and phase-set
label matches the accepted baseline. **No new tracked gap closes; fourteen
competitor targets remain open.** No panel, required-site entry or expectation
is refreshed. The physical allele representation bug is fixed, but it is not
sufficient to close the audited gap.

## Regression verification

The new predicate case reproduces the 44-base shifted insertion and checks
motif rotation, shifts in both directions, original source labels and rows,
known qualities, reference mismatches and ambiguity, different allele lengths,
compound insertion/deletion paths, reference boundaries and the worker's
reference-cache channel. It passes 108 assertions. The two previously regressed
owning controls pass 39 assertions before the full panel.

Standalone units pass. Predicate tests pass 599 assertions / 28 cases;
port parity passes 27 / 9; upstream parity passes 168,696 / 7.
All 105 distinct native commands run fresh on the final binary in 267.45 s
with four workers. Binary-SHA and exact-command verified scoring passes
4,082 assertions / 53 cases over the unchanged 93-window panel and owning
orientation controls. No expectations are weakened.
No new compiler warnings or whitespace errors.
