# Representation repair: step 1

## Scope

Remove order-dependent candidate routing before introducing reference-normalized
aliases or merging diploid loci. The starting executable is
`c6bcf86a43b5301cd90960eeb8b525d14944494dfb8381f9b676c5641952e674`.
Its source and binary snapshot remain in `test_data/tmp_representation_step1/`.
The evidence motivating this change is the preceding
`2026-10-06-joint-evidence-investigation/` representation audit.

`CandidateIdentityIndex` in `src/allele_identity.*` retains a unique candidate
index or a permanent ambiguous entry for each exact identity. Adding the same
candidate again is idempotent. Adding another candidate does not select either
one, including when the input order is reversed. Raw and selected-ALT sequence
indexes must agree if both contain the queried identity. Distinct alternatives
sharing a topology key can still match independently through sequence metadata.

The seam-recovery parent indexes, rebuilt transfer index and whole-chunk BAM
overlay now use that rule. A BAM candidate that matches an ambiguous existing
identity cannot be appended as a private new locus. The rebuilt transfer index
includes every parent's selected sequence alias, including duplicate aliases,
so reindexing cannot accidentally make an ambiguous parent unique. Observations,
base-quality certificates and retained source gauges are transferred only to
unique destinations. Existing graph rows and their graph observations remain
intact; this step does not deduplicate candidate rows or combine read counts.

The index compares exact position/type/reference-length/ALT keys within one
contig. It does not yet left-align shifted repeats or express complete ALT1/ALT2
genotypes. Candidate identity is separate from MSA provenance and source PS.

## Verification

The standalone `test_allele_identity` checks duplicate-input idempotence,
permanent ambiguity in both orders, distinct insertion and replacement alleles,
missing sites, selected alternatives sharing a graph key, duplicate sequence
identities, and agreement/conflict between the two matching indexes. A temporary
mutant restoring first-entry-wins fails five checks; the fixed implementation
passes. One hundred warm executions average about **1.3 ms**, including process
startup. `identity-checks.json` records the measured result and failing checks.
The target is included in `make unit-tests` and `make gap-dev-check`, and needs
no BAM, graph, truth map or pipeline replay.

Build and unit/predicate checks pass without new warnings. The 47 existing
predicate cases pass 1,551 assertions. HiFi TSV/VCF golden outputs and t1/t4
determinism pass; ONT TSV/VCF golden outputs pass. All four replay-cache tests
pass. Remeasurement on the original input and frozen competitor BAM validates
identical primary alignments for all 141 benchmark windows; the resulting
`hiphase-verified.tsv` is byte-identical to `src/test_gap_hiphase.tsv`.

The full owning/window verification runs every registered `all gaps` section
with its existing assertions and expectations. `check_windows.py` selects exact
sections and runs three independent processes, sharing the replay cache and
using separate output directories. This retains parental orientation checks,
the >=80% all-primary denominator, and both total/core HiPhase parity checks
for certified closures. Historical open gaps keep their measured floors.
All **107 registered checks pass**. Cold parallel verification takes 666.4 s;
the subsequent standard `make window-tests` passes **18,138 assertions in five
test cases** using completed replay states. All 141 native panel measurements
are unchanged from the preceding audit, including the historical deficits.
`window-checks.json`, `panel-after.tsv` and `panel-comparison.json` preserve
the results. The verified executable is
`088b5c56e17d04828af0ccbfede26b9ca3ff65476a62f796066af7a6ac63eb9d`.

The before comparison uses the frozen pre-change full output and preceding
native contract, not the initial suite launched during preparation: its later
replays overlapped the executable rebuild. The independently selected 107
checks and subsequent standard suite use the unchanged verified executable.
The full-chromosome replay also passes the preservation audit:

- 82,287 candidate rows are byte-identical; all 64,483 VCF data records match.
- All 256,612 primary read HP/PS pairs and parental classifications match.
- Parental totals remain 230,965 correct, 6,620 discordant and 19,027 unphased.
- The 250 variant blocks retain N50 955,496 bp and largest block 3,012,193 bp.

`audit_full.py` and `full-comparison.json` preserve the comparison. Thus the
new collision rule has no observed biological-output effect on this chromosome;
the first-entry-wins defects are exercised by the focused synthetic states.
This is a verified foundation for normalized aliases, not evidence that the
ten repeat-shifted duplicate pairs or remaining phasing deficits are repaired.

## Reproduction

From the repository root:

```bash
make -j$(nproc)
make test_allele_identity
./test_allele_identity
make unit-tests gap-dev-check
make check
python3 scripts/test_cache_gap_replay.py
python3 evaluations/2026-10-06-representation-step1/check_windows.py \
  --binary ./pgphase --out test_data/tmp_representation_step1/windows-step1
bash evaluations/2026-10-06-shared-insertion-source/replay.sh \
  test_data/tmp_representation_step1/full-current
/home/kokyriakidis/micromamba/envs/bench-phasers/bin/python \
  evaluations/2026-10-06-representation-step1/audit_full.py
```

## Following steps

1. Add a shared reference-anchored allele identity, tested against a normalization
   oracle and complete local haplotypes. Cover repeat shifts, common context,
   compound replacements and slice boundaries before transferring observations
   across normalized aliases.
2. Merge equivalent descriptions while retaining provenance and one observation
   per original molecule/locus. Duplicate descriptions must leave evidence
   weight unchanged; disagreements must remain visible and uncertain.
3. Keep distinct sample alternatives in a complete diploid locus instead of
   independent REF/ALT questions. Preserve the MSA consensus window and read
   confidence, and derive genotype/orientation jointly rather than importing
   the source HP as independent truth.

Each step first exercises immutable in-memory states, then its affected owning
contexts, then the unchanged full regression/parity gates. This step claims no
new gap closure or N50 improvement.
