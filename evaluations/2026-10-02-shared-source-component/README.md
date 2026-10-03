# Preserve original BAM component identity across graph attachment

## Bug and correction

At **chr20:34,046,350–34,055,920**, the MSA-verified 32-base insertion has
20 independent usable physical BAM links to the right catalog SNP: eight REF
and twelve ALT molecules. The existing likelihood test accepts the connection.
Its left catalog/BAM component, however, has been relabeled from original
BAM PS 34,018,169 to catalog PS 34,018,168. The complete original BAM source
has a later weak cut at 34,058,504; that cut is outside the retained prefix.
The old path checks neither certify the graph-owned shared rows from their
exact BAM provenance nor distinguish that original identity from the attached
label.

A second transfer defect strands the private deletion at internal coordinate
34,094,605 under the earlier catalog label, across the source cut. The nearest
row before the cut belongs to another owner. Inspecting only that nearest
owner leaves the earlier owner's orphan island attached and stretches its
apparent block extent across the cut. This overlap also hides the real break
from an extent-only gap audit.

The correction checks eligible additional owners with an exact shared clean
SNP in the source component immediately before the cut. An observed read
bridge or matching far-side catalog anchor preserves their existing joins.
Only unlinked private islands are detached; their allele gauge and exclusive
source reads are restored together. The subsequent existing component-aware
attachment handles the island's actual neighbor. The original BAM label
participates in the closest-row lookup even though it needs no restoration.

After detachment, an independently certified physical insertion-to-SNP bridge
can use a cut-free attached source component as its left path certificate.
Every phased anchor, including catalog-owned rows, must have exact adoptable
provenance in one original BAM source and the same source/current allele gauge.
Weak and quality cuts inside the retained heterozygous extent veto the path.
Both allele classes and the existing physical likelihood gate remain required.
There is no truth input, competitor input, coordinate exception, allele merge,
new realignment, or lowered evidence threshold in production.

The detachment implementation is moved into the adapter so synthetic tests
exercise the actual transfer code and read-gauge restoration. Exploratory
instrumentation in `source-trace.log` is absent from production.

## Scope experiments

A trial checking all historical owners improved aggregate read accuracy but
split protected blocks: VCF blocks 336→340 and span N50 756,878→739,888 bp.
Restricting owners merely to the preceding component was insufficient. It
also wrongly omitted an original-source row from the closest-row search.
Those trials are rejected. The final eligibility and physical-link checks
preserve the accepted 53–54 Mb VCF and every read HP/PS tuple exactly
(`protected-53-comparison.json`).

## Local evidence

All tools are scored on the same 123 gap-overlapping input read names. Local
accuracy is the majority parental orientation per output block; the new window
regression separately checks both flank orientations and exact allele parity.

| Tool | Phased / correct / discordant | Local blocks | Correct reads in dominant block |
|---|---:|---:|---:|
| Accepted pgphase before | 123 / 118 / 5 | 2 | 68 |
| pgphase corrected | 123 / 118 / 5 | 1 | 118 |
| HiPhase, native DV calls | 123 / 118 / 5 | 1 | 118 |
| HiPhase, saved pgphase calls | 123 / 118 / 5 | 1 | 118 |

The saved pgphase-callset HiPhase run uses Oct 1 calls, not this new binary's
calls. Local concordance is **95.93496%** for all three outputs. The owning
34–35 Mb pgphase solve improves from 3,560 phased / 3,532 correct / 28 discordant
to **3,567 / 3,539 / 28**: seven additional correct reads, no additional errors.
The three right SNPs move to the left PS with reversed GT, which is the
supported parental connection. No VCF key is lost or added.

## Regression coverage

The committed panel grows to **97 coordinate cases / 83 required spans**.
The new case uses owning-chunk context and requires the unchanged insertion,
left SNP/insertion agreement, three right SNPs in the same block with opposite
ALT orientation, parental flank agreement, read separation at least 118/123,
and the owning ceiling of 28 errors. Existing read floors are preserved.
The accepted pre-fix binary fails eight of the 29 new assertions. Synthetic
units cover source identity, shared-row provenance, mixed gauges, cut extent,
restored private-row/read orientation, preserved physical and catalog links,
component scope, and the original-source closest-row lookup.

## Reproduction

Build with `make -j8`, run `make unit-tests`, `./test_phase_predicates` and
`make window-tests` with the chr20 fixture and derived parental map present.
The comparison scripts are evaluation-only; their `--help` describes inputs.
Fresh native requests are checked against a frozen binary SHA256 before and
after execution. Reusing their outputs in the Catch runner avoids rerunning
identical commands; matrix regressions run natively with the same binary.

## Whole chromosome

The final frozen build has SHA256
`06b32560fe18b80332994314fa9bbcaa3962bbbc445bb4108f05cf27e7f3602c`.
A fresh eight-thread chr20 run takes 295.76 s while native window jobs are also
running; this is a validation runtime, not a controlled competitor benchmark.

| Measurement | Accepted before | Corrected |
|---|---:|---:|
| Phased/scored reads | 237,127 | 237,134 |
| Correct / discordant reads | 229,974 / 7,153 | 229,981 / 7,153 |
| Read concordance | 96.983473% | 96.983562% |
| Read phase sets | 670 | 669 |
| VCF blocks | 336 | 335 |
| VCF keys | 63,492 | 63,492 |
| VCF span N50 | 756,878 bp | 756,878 bp |

All 229,974 previously correct reads remain correct, and all 7,153 discordant
reads retain their individual classification. Exactly seven formerly unphased
reads become correct. There are 179 changed HP/PS tuples and three changed
VCF rows; zero keys are lost or gained. `full-read-comparison.json` records
aggregate counts and individual-state transitions. The change has no effect
outside the owning 34 Mb chunk in this chromosome run.

A fresh noncentromeric extent audit still nominates 41 distinct coordinate
pairs (51 competitor records) across 315 VCF extent gaps. These are candidates,
not 41 proven allele bridges. This newly fixed break was masked by the orphan
row's overlapping pre-fix extent, so the nomination count does not fall.
`remaining-gaps.json` retains the evidence for subsequent investigations.

Final validation: the build introduces no warnings; all standalone unit tests
pass; phase predicates pass 822 assertions / 37 cases; window tests pass
4,806 assertions / 62 cases. All 105 fresh native replay requests complete in
393.85 s against the frozen SHA256; the additional matrix regression runs
natively. `native-validation.json`, `native-panel.log`, and the gate logs record
these checks. The exact new span is required, rather than silently accepting
any increase in span count.
