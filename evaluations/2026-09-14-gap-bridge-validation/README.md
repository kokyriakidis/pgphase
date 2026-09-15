# Validate the allele evidence crossing recovered gaps

The original recovery result had 2,427 discordant reads and ten incorrect
joins among original blocks with clear parental majority labels. The fresh
recovery-disabled run had 586 discordant reads. These are read-origin metrics;
variant switch/Hamming error was not measured here. The input is graph/GAF plus
surjected BAM/MSA evidence. Parental truth is evaluation-only; no DeepVariant
calls enter phasing or acceptance.

## Fixes

1. Honor `gap_link_supported` for ordinary biallelic MSA indels inside the
   recovery gap. The validation was computed but graph-node admission checked
   it only for multiallelic insertions and homopolymer links. A rejected ordinary
   indel could still connect flanks using abundant reference observations.
2. Evaluate an MSA allele edge as a 2x2 table. Both alleles at both endpoints must independently
   favor the same relative orientation (including after transposition), and the combined net margin must meet
   the existing `min_block_link_reads`. Previously, pooled winning counts could
   accept an edge even when one haplotype favored the opposite orientation.
   This uses the observed alleles across the edge, not proposal read HP labels
   or the number of reads subsequently attached to the flanks.

These changes apply to recovery's allele graph. They introduce no new
threshold or blanket ban on homopolymer context. Existing single-molecule
bridges retain their separate BAM/graph anchor-validation path.

## Experiments

All full runs reuse `/tmp/chr20-evidence-v5.gapev`; each output directory
contains its exact `command.json`, log and parental read evaluation.

| Arm | Output directory | Discordant reads | Audited wrong joins |
| --- | --- | ---: | ---: |
| Before | `/tmp/pgphase-confirm-original` | 2427 | 10 |
| Indel validation only | `/tmp/pgphase-validated-gap-indels` | 1662 | 5 |
| Rejected HOM-node exclusion | `/tmp/pgphase-validated-gap-genotypes` | 1934 | 6 |
| Row-only allele check | `/tmp/pgphase-validated-gap-alleles` | 807 | 0 |
| Final symmetric allele check | `/tmp/pgphase-validated-gap-symmetric` | 760 | 0 |

The HOM-node exclusion trial was rejected: although it enforced a plausible
category invariant, changing the available sites exposed an incorrect
alternative orientation at 56.15 Mb, adding 272 discordant reads relative to
the indel-only arm. That exclusion is not in the final code. Merely removing
sites cannot certify the remaining bridging evidence.

## Validation

`make -j8 unit-tests benchmark-tests` passes. Two new C++ regressions cover:

- An unsupported ordinary MSA indel must not join separate flank read groups;
  a validated indel still can.
- A pooled 9:3 edge must be rejected when its second haplotype opposes the
  orientation 3:1. A consistent control joins, and swapping one site's HP gauge
  does not alter acceptance.

The first regression failed before its fix. A separate mutant object with only
the new per-haplotype edge gate removed fails the second regression in both HP
gauges. Existing singleton insertion/deletion and independent-SNP bridge tests
continue to pass. Build logs are in `/tmp/gap-balanced-build.log` and
`/tmp/gap-balanced-all-tests.log`.

`compare.py` checks original blocks transform uniformly, measures common-read
status changes, and scores original-gap endpoints against parental majority
orientation. Blocks with >10% original discordance are marked uncertain rather
than used as clear edge labels. The artifacts cover all 280 original gaps, not
just successful competitor cases. This is chr20 regression evidence, not an
untouched validation set or a calibrated false-join probability.

## Final result and limits

The final code has 86 joins: 84 have correct relative orientation under the
original-block truth audit, two have uncertain original-block truth. All ten
previously identified incorrect joins are now split. Relative to the prior 125
joins, 27 previously correct joins and two uncertain joins are also withheld.
The 57.84 Mb control remains joined; 62.41 Mb is withheld. This is a correctness
safeguard with a contiguity cost, not recovery of all competitor-closed gaps.

Read discordance is 760/188714, versus 2427/189042 before. Corrected read-order
switch/flip diagnostics are 138/208, versus 205/267 before (the historical
351/414 used mixed parental coordinates). Relative to the no-recovery cohort,
40 originally concordant reads become discordant, versus 1452 before. Relative
to the previous recovery output, 28 common reads become discordant and 1561
become concordant; 331 formerly evaluated reads are no longer evaluated and
three are newly evaluated. These counts do not support a claim of zero read
regressions or zero errors. All 278 compared original blocks retain a uniform
HP/PS transformation.

The final run took 290.6 seconds including read-truth evaluation; its log reports
206.3 seconds for cached gap solves. Final full output and exact invocation are
under `/tmp/pgphase-validated-gap-symmetric`. Regression tests also cover
transposing endpoint observations. Source hashes and all ablation summaries
are retained beside this document.
