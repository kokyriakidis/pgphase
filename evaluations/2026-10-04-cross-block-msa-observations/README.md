# Preserve verified MSA observations across independent source blocks

## Bug and correction

`apply_deferred_msa_observations` discarded a fixed-consensus allele call
when its read already belonged to another source phase set. The operation
adds an observation; it does not move the read's HP or PS. Numerical HP
values from independent blocks cannot decide whether that call conflicts.
Rejecting it removed the evidence needed to compare those blocks later.

Retain the verified call while preserving both block gauges. The existing
same-phase-set allele contradiction veto remains. Discovery, retry selection,
candidate representation, genotype and stitch criteria remain unchanged.
Calls are still coalesced by exact key/read ID, conflicting recalls abstain,
and the selected source's graph-anchor checks still apply. No new CIGAR
recall or alignment is introduced. Production does not read parental truth
or competitor output.

## Measurements

Baseline is commit `0f44897874558f0cdab74738f3a4190bbe7bd447`, full native
chr20 output `test_data/tmp_gap_fix47/full-alt-only`. The trial output is
`test_data/tmp_gap_fix48/full-cross-block`, using the same BAM, graph catalog,
GAF and reference, with eight threads.

Twenty-five additional complementary MSA call pairs survive at VCF anchors
11,586,531 (8 reads), 13,619,764 (1), 15,101,262 (10), and 41,880,908 (6).
The eight separate insertion rows gain depth; all variant keys, genotypes
and phase labels are unchanged. Six owning-chunk probes (3, 10, 15, 21, 34
and 35 Mb) confirm the 15 Mb depth gain without a new connection.

Full chromosome read tags are identical: 237,308 truth-scored phased reads,
230,211 correct and 7,097 discordant, 97.009372% concordance. There are still
331 VCF blocks and 653 scored read phase sets. All 91 connected panel
coordinates remain connected; the 13 open panel coordinates remain open.
This corrects observation loss, not a newly closed gap.

See `read-parity.json` and `block-audit.json` for exact output comparison.
Parental truth is used only for these evaluations and regression assertions.

## Verification

- `make -j8`, no new compiler warnings.
- `make unit-tests`: passes.
- `./test_phase_predicates`: 1,517 assertions in 47 cases pass.
- `./test_gap_windows 'complementary insertion*'`: 1,310 assertions in
  two owning-context orientation regressions pass; expectations unchanged.
- `make check`: HiFi/ONT goldens and HiFi thread determinism pass.
- Full native chr20 read-tag parity and all committed panel connection checks.

The hermetic MSA test exercises both numerical HP values in another phase
set and verifies that observations survive while all read and candidate
HP/PS labels remain fixed. A separate section retains the contradictory
same-phase-set rejection.
