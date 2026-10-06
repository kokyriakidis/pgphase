# Deletion-masked SNP read gap at chr20:1,194,233–1,196,894

The interval already spanned in the variant output, but its connected read core
was below HiPhase. On the same 91 truth-scorable original primary alignments:

| output | correct / all overlaps | connected core correct | discordant | unphased |
|---|---:|---:|---:|---:|
| pgphase before | 90/91 | 90 | 0 | 1 |
| pgphase after | 91/91 | 91 | 0 | 0 |
| HiPhase | 91/91 | 91 | 0 | 0 |

The previously unassigned maternal read
`m84031_231217_062403_s3/204411306/ccs` has a primary MAPQ60 alignment and an
exact Q40 six-base deletion of CACCAC at 1,194,228. It deletes the verified
opposite-haplotype SNP positions 1,194,230 and 1,194,233. Its transferred
primary and BAM profiles retain deletion ALT but cannot call those SNPs.
Ordinary rescue counts six primary graph donors, three per haplotype: even
6/6 agreement has a one-sided binomial tail of 0.015625, above its 1% gate.

The new fill uses the physical alleles of already phased reads in this existing
core, including staged BAM-only output rows. There are 72 callable Q30-or-better
votes: 37 deletion ALT / haplotype 1 and 35 REF / haplotype 2, no disagreements.
It requires both allele/haplotype classes, the existing two-sided association
p <= 0.01 and one-sided 95% Wilson discordance <= 10% gates, and the marker's
current orientation. It fills only entirely unassigned reads with exact CIGAR
ALT, matching surviving sequence, both profile channels calling ALT, masked
verified opposite-haplotype SNPs, and no contrary phased call in either channel.
MAPQ30/Q30 and pure MSA-verified deletions <=32 bases are mandatory. Its reference
and BAM handles belong to the worker; BAM contig IDs come from site metadata.

This local assignment does not require a complete independent BAM source path.
That condition protects whole-block joins; here the physical allele-to-HP gauge
is measured directly against existing core reads. Candidate genotypes, phase-set
connections, prior core assignments and prior rescues remain unchanged.

The certified owning replay is 1,000,001–2,000,000. Its audit compares 4,197
primary output reads and 1,528 variant records with the task-start binary:
exactly one new correct assignment, all 4,044 prior assigned tags unchanged,
10 discordant reads unchanged. The target is added to `src/test_gap_certified.tsv`
with its existing `spans=1` expectation and parental orientation assertions.

Final binary SHA256:
`6a8057fc5d3875b7fb3a608a579efe47a2f4d1ffce99aac84bdd0dcec6321f04`.
Task-start binary:
`f9efac33e864064418fe09e70c4502b4db43dae48f6825bc4e9556d699e6ed1b`.

Fast fixtures run in 0.03 seconds; a cached real owner check takes 0.34 seconds.
The changed production binary requires fresh replay state for final regressions.

`compare_hiphase.py` verifies all 91 competitor alignments against original
start/end, CIGAR and sequence, including abstentions in the denominator.
`audit_reads.py` verifies the owner; `audit_panel.py` compares matching saved
native outputs across the complete suite. Measurements live in their JSON files.

Final validation passes all 84 registered gap checks: 83 mechanisms and all
112 panel windows, exercised in 195 selected invocations with 82,416 assertions.
The three focused gap unit cases, 47 predicate cases (1,551 assertions), unit
binaries, four cache-helper tests and standard validation gates also pass.
The panel retains all 99 spans and improves full-contract passes from 33 to 34.
The saved-output audit compares 329 matching native replays. All variant records
and existing read tags are preserved in every replay; nine overlapping outputs
contain the same single new correct read assignment. No new discordant read
appears in the audited replays. This validation uses panel and owning replays;
a separate full-chromosome replay was not repeated for this patch.
