# Close the 23.461–23.481 Mb gap without changing read rescue

The retained change connects chr20:23,460,963–23,480,815 and preserves the
preceding 23,421,003–23,445,252 connection. A fresh owning 23 Mb replay has
125/125 correctly oriented gap reads before and after; their core phase sets
become one. Disjoint 10 kb flanks agree with the same parental orientation:
29/29 left and 31/32 right. Both insertion descriptions at 23,480,815 retain
complementary genotypes and unchanged allele depths.

## Mechanism

The ordinary whole-flank BAM validation solve omits the targeted source's
verified insertion observations. When the original physical geometry and a
verified noisy right-source insertion exist, retry validation independently
with unplaced MSA observations and insertion recall enabled. The ordinary
whole-block shared-site gauges, source-cut checks and exact physical
certificate remain required. No discovery counts, genotype rows or read calls
from this supplemental solve are transferred.

Evaluate the proposed relation on a disposable stitch state using the existing
outer-allele and gauge conflict vetoes. Retain complete original anchor keys and
orientations on both sides. Apply the union after read rescue and equivalent
insertion joins, requiring every original anchor on each side to retain one
uniform live phase set and gauge. Recompute parity in those live gauges and
update core candidates and primary reads in all chunks carrying the absorbed
phase set. Output-only rescued read labels retain their independently assigned
HP/PS, as they already do for equivalent insertion unions.

## Owning replays

| Replay | Scored before → after | Correct before → after | Discordant before → after | VCF blocks before → after | Span N50 before → after |
|---|---:|---:|---:|---:|---:|
| 23–24 Mb | 3,602 → 3,602 | 3,542 → 3,542 | 60 → 60 | 7 → 6 | 208,433 → 680,357 |
| 5–6 Mb control | 4,024 → 4,024 | 4,001 → 4,001 | 23 → 23 | 3 → 3 | 605,720 → 605,720 |

The 5 Mb output has zero changed tags or VCF rows. All owning 23 Mb truth
statuses are invariant; 857 tag changes are uniform core phase-set/gauge
relabels. All 933 variant keys, genotype alleles and non-phase fields remain.
The source audit checks 23,912 independent calls and 40,669 overlay calls,
with zero source disagreements.

## Full chromosome and paired context

The final chromosome run preserves all 256,610 output reads, 237,330 scored
reads, 230,235 correct assignments and 7,095 discordant assignments. Every
individual truth status is unchanged. All 64,188 variant keys and non-phase
fields remain; every old phase block moves uniformly. Exactly one formerly
open gap closes, no tracked gap reopens, and 93/106 tracked gaps are connected.
VCF blocks fall from 330 to 329, read phase sets from 652 to 651, and span N50
rises from 806,449 to 856,770 bp. There are 857 core tag relabels and 256 VCF
phase changes, with no output-only read gauge changes.

The paired 22–24 Mb replay preserves all 7,428 scored reads: 7,307 correct and
121 discordant, with every individual truth status unchanged. Its VCF blocks
fall from 13 to 11 and N50 rises from 461,736 to 1,322,284 bp. The extra short
replay boundaries joined there are already connected in the accepted full
chromosome; they are not additional whole-chromosome gap closures.

## Rejected alternatives

`rejected-alias-parity.json` records the discarded source-alias shortcut: it
lost four correct reads and reopened the preceding connection. The first
successful core closure applied its supplemental bridge before read rescue.
Its full output (`rejected-read-cohort-{parity,audit}.json`) lost five correct
5 Mb tags and changed a 23 Mb read-only marker cohort, making three previously
correct reads discordant and two previously discordant reads correct. Its
variants had exactly the intended closure, but its read behavior is rejected.

Separating independent BAM validator cohorts and refreshing inherited BAM SNP
read labels earlier did not repair these changes. The raw whole-BAM solve
matrices were byte-identical. The actual cause is excluded-site marker rescue
consuming altered core phase-set cohorts; deferring the union preserves those
inputs. No original graph path or physical quality threshold is relaxed.

## Validation

Build, all standalone units and 1,536 predicate assertions in 47 cases pass.
BAM HiFi and ONT TSV/VCF goldens and thread determinism pass. The new owning
and paired regression cases pass 72 assertions. The original binary fails the
five connection checks; the discarded early-join binary fails all three exact
rescued-PS retention checks. No existing expectation is relaxed. The broad native run passes 8,168
assertions across 79 cases using 110 fresh final-binary pipeline requests.
The additional paired case and strengthened rescue-label checks are verified
in the final focused run (72 assertions across two cases), covering all 80
distinct cases in the final test source.

## Reproduction and provenance

Run `reproduce.sh BEFORE_BINARY AFTER_BINARY OUTPUT_DIR` from the repository,
with a Python interpreter providing pysam selected through `PYTHON`. It replays
the owning 5 Mb and 23 Mb chunks, paired 22–24 Mb, and full chr20, then measures
parity and runs the strict new-closure audit. `evaluate.py` asserts no correct
read becomes discordant or loses its tag, no old block splits or mixes gauges,
no allele/count changes, no key loss, no reopened tracked gap, and the new
closure with consistent parental orientation on disjoint flanks.

Baseline binary SHA256:
`69d998cab0b2319f6875e3514a7be4eaa1dc9b514419831614aeb012a13c4386`.
Final binary SHA256:
`3a227a265b28ca21b15e4e59e272f0dc2cc06e7feb2046d50108365e481152c5`.
Fresh output roots are `test_data/tmp_gap_fix54/deferred-gated/{5,23}`,
`test_data/tmp_gap_fix54/focused-gated/`, and
`test_data/tmp_gap_fix54/full-gated/`; the accepted full baseline is
`test_data/tmp_gap_fix52/full-current/`.
