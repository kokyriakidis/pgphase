# Ignore homozygotes in a detached BAM run's phase certificate

## Retained correction

`bam_source_run_supported` used every candidate carrying the run's PS as a
phase anchor. BAM homozygotes can retain that PS with identical haplotype
alleles. Such a row cannot orient the run, but previously rejected the whole
certificate. Including its position would also extend the weak-cut interval
beyond the actual heterozygous anchors.

Ignore known homozygotes both when indexing source provenance and when
checking the run's gauge and extent. Keep the existing requirements for every
heterozygote: an injected BAM row, unique adoptable provenance, one original
source, a consistent same/swapped allele gauge, at least two distinct anchor
positions, and no original source weak cut inside their half-open extent.
Unknown allele rows still veto the certificate. Homozygotes alone cannot
certify a path. No read or candidate genotype is changed by this predicate.

Move the existing certificate from `graph_collect.cpp` to the already-linked
`graph_bam_adapter.cpp` so the actual production function is directly tested.
Its caller and all stitch likelihood/count rules remain unchanged. No new
object dependency or alignment step is added. Runtime decisions contain no
truth, comparator outputs, or fixture coordinates.

## Regression and production evidence

The adapter test adds 15 checks covering valid and globally flipped runs,
REF/ALT homozygotes, homozygous provenance duplicates, extents, malformed
rows, conflicting heterozygous gauges, interior cuts, absent metadata,
independent sources, duplicate heterozygous provenance, graph-owned anchors,
unadoptable gauges, singleton runs, and homozygotes alone. Before the fix,
the new source-labelled-homozygote and extent checks fail. All current checks
pass. `source-run-before.log` records those two reproduced failures.

A diagnostic whole-chr20 run reaches source-labelled homozygotes in BAM runs
PS 32,168,699 and 49,700,407. Removing that false veto does not establish
an additional accepted stitch on this fixture: other certificate/stitch
requirements still apply. The trial's chromosome metrics are unchanged.
The final production build is independently rerun and checked below.

## Additional experiments, not retained

- At 57,841,772–57,866,713, the initial source has complementary deletion
  rows at canonical 57,854,341 and 11 MAPQ-30 spanning reads, with zero
  callable pairs to the right SNP. The existing retry requires an outer
  seam of a grouped solve. A trial allows focused complementary-row dropout
  admission for an interior seam, retaining ownership of the other seams.
  It restores observations but fails both the complete-source-path and
  heterozygote-preservation certificates. It closes no tracked target and
  is not retained. This does not show that all recovered observations are
  erroneous; it shows that this replacement cannot safely satisfy the
  existing source contract.
- Checking singleton graph/BAM gauge consistency within each independently
  numbered BAM PS makes no change in the replayed remaining-target chunks.
  No relaxation of the full-block validator is retained.
- Allowing the existing BAM-only run certificate in complementary-deletion
  and single-indel stitch helpers makes no chromosome improvement. The
  production indel routing stays unchanged.

The accepted baseline is `/tmp/pgphase-gap-next6/verified-full`; current work
and native replays are under `/tmp/pgphase-gap-next7`. Native HiPhase and
Longcalld comparisons are not rerun. The existing 92-window panel and all its
floors are retained; no newly closed gap needs a new panel row this turn.

## Final result and validation

Fresh default graph+recovery chr20 run, same committed BAM, reference, catalog
and GAF, eight threads, 185.20 seconds:

| Metric | Before | After |
| --- | ---: | ---: |
| Primary read names | 256,586 | 256,586 |
| Phased/truth-scored reads | 237,101 | 237,101 |
| Correct / discordant reads | 229,920 / 7,181 | 229,920 / 7,181 |
| Read concordance | 96.971333% | 96.971333% |
| Read phase sets | 686 | 686 |
| VCF records / blocks | 62,798 / 340 | 62,798 / 340 |
| VCF span N50 | 739,888 bp | 739,888 bp |
| Previously established panel spans | 74/92 | 74/92 |

Every primary qname and HP/PS pair and every nonheader VCF row is identical
to the accepted baseline (`identity.json`). No additional tracked gap closes;
fourteen comparator targets and four intentional controls remain open.

Warning-free build; shared unit tests, 398 predicate assertions, 27 port-parity
assertions, and 168,696 original-upstream assertions pass (`gates.log`).
All 105 distinct window pipeline commands execute on the final production
binary using four concurrent workers in 238.21 seconds (`native-replays.log`).
The unchanged Catch runner then scores those fresh outputs using a wrapper
that verifies the production SHA256 and exact canonical command. No prior
native output is substituted. The input command inventory came from the
previous unchanged panel; that inventory contains inputs, not test results.
Temporary reproduction tools are `parallel-panel.py` and
`use-native-outputs.py` under `/tmp/pgphase-gap-next7`; ordinary
`make window-tests` retains its existing sequential behavior.

All **4,023 window assertions in 52 Catch cases** pass across the unchanged
92-window panel (`panel.log`). `git diff --check` passes. No files are staged,
committed, or pushed this turn.
