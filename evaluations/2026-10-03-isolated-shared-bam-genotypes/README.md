# Retain isolated verified BAM genotypes on shared graph rows

## Defect and repair

An exact sequence match could suppress a verified BAM heterozygote when its
catalog row was an unphased repeat. The existing ownership rule retained a
BAM genotype only for suffix-padded matches. At 21.831 Mb this kept the BAM
observations but discarded the `CT→C` genotype and its MSA proof from the
live phased candidate.

Preserve the existing suffix-padded ownership rule. For other exact binary
matches, retain the BAM genotype when every neighboring heterozygote in its
source phase set is separated by a weak source cut. Co-located source rows,
source singletons and sites with a supported neighboring edge retain the
existing rule. The predicate handles candidate order and indel anchor
coordinates explicitly and uses the already computed source path evidence.

The adopted row keeps its physical BAM key, genotype, depths and MSA proof,
with an independent unused PS. Its raw source PS cannot attach it to a graph
block. Source adoption is disabled for this independent row; ordinary paired
read evidence can establish a later connection. The original selected
catalog metadata still describes the same sequence allele.

Keep graph calls in their independent observation channel. The primary and
BAM channels use source calls, with unknown primary slots left unknown.
Graph rescue for these rows infers orientation from primary read associations
to one established block, rather than borrowing the independent BAM gauge.
The normal singleton confidence rule applies; two supported block gauges
cannot be pooled. Scope this behavior with `bam_independent_genotype` so
existing MSA rescue markers keep their established behavior.

Channel arrays must span the union of primary, graph, BAM and quality indices.
Sizing them from primary calls alone can write a retained channel call outside
the allocated array. The new channel-retention experiment exposed a heap
allocation failure; the union extent removes that invalid indexing.

Production uses source candidates, allele observations and weak-cut evidence.
Truth and competitor outputs are evaluation inputs only. No source-coordinate
admission, new alignment stage, genotype threshold or row merge is introduced.

## Full chr20

| Metric | Accepted baseline | Fixed |
|---|---:|---:|
| Truth-scored phased reads | 237,225 | 237,225 |
| Correct assignments | 230,105 | 230,105 |
| Discordant assignments | 7,120 | 7,120 |
| Conditional accuracy | 96.998630% | 96.998630% |
| Read phase sets with truth-scored reads | 654 | 654 |
| VCF keys | 64,185 | 64,188 |
| VCF phase blocks | 330 | 333 |
| Span N50, bp | 806,449 | 806,449 |
| Span NG50, bp | 684,798 | 684,798 |
| Coordinate panel connections | 89/102 | 89/102 |

Every read's HP/PS tuple and truth-correctness state is unchanged. Every old
VCF row survives unchanged, no phased SNP disappears, no old SNP block changes
gauge and no old block merges. The three additional records are independent
single-site blocks:

| VCF position | Alleles | Source AD | Final PS |
|---|---|---|---|
| 21,294,113 | `TA→T` | 3,2 | 21,294,113 |
| 21,831,480 | `CT→C` | 11,4 | 21,831,480 |
| 44,083,059 | `C→CT` | 16,54 | 44,083,060 |

The 21,823,066–21,844,359 gap remains open. Its missing BAM site now survives
transfer; this result does not claim a new supported flank connection.

## Regression and validation

Add an owning 21–22 Mb regression requiring the exact recovered deletion,
source depths, positive independent PS, separate flanking blocks and the
existing read coverage/accuracy floors. It fails on the starting executable
because the middle row is missing (two rows rather than three). The final
executable passes all 26 assertions. Synthetic checks cover either supported
neighbor, reversed source order, co-located contrasts, distinct source PS,
separate graph/BAM evidence, ambiguous graph gauges and genotype preservation.

- Main build, `make unit-tests`, and `make predicate-tests` pass without new
  warnings; predicates have **1,404 assertions / 46 cases**.
- Complete window suite: **7,198 assertions / 73 cases**, all pass. No
  expectation floor or required connection is weakened.
- HiFi/ONT TSV and VCF goldens pass; HiFi one/four-thread output is identical.
- Run **106 fresh native panel requests**, then score outputs with executable,
  exact normalized argument, input-stat and output-hash checks. The new owning
  test reuses the identical owning-chunk request; uncached requests run natively.
- Full chr20: **280.60 seconds**, eight threads. Concurrent native panel:
  **395.14 seconds**. These are validation timings, not a competitor runtime
  comparison.

`chromosome-parity.json`, `block-transfer.json`, `owning21-parity.json` and
`span-metrics.json` hold the measurements. The owning chunk retains exactly
3,321 scored / 3,110 correct / 211 discordant reads and adds two records.

Starting executable SHA256:
`336996e6338b3abc4a55b662231e71e86d378f4222615c4e39aa6cc94a645ce4`.
Final executable SHA256:
`772e57fe587c75f7acdfb1b312cf791573194e31ab1a1732bc667a0ea10655a0`.
Starting chromosome: `test_data/tmp_gap_fix43/full-accepted/`.
Final chromosome: `test_data/tmp_gap_fix45/full-scoped/`.

## Rejected broader trials

Allowing all newly matched verified genotypes to start independent blocks
changes existing rescue evidence, drops correct read coverage and leaves
unsupported source components exposed. Restrict adoption to isolated source
sites and preserve their graph channel. An unscoped rescue-marker fallback
also removes 23 read assignments outside the owning chunk; the explicit
independent-genotype state restores all existing assignments. The earlier
whole-source ownership trial is documented separately in
`evaluations/2026-10-03-shared-bam-genotype-trust/`.
