# Representation step 13: compose neighbors through mandatory matched subpaths

Add a mapped fallback for compound parent/neighbor overlaps, using the verified
all-optimal REF-to-parent-ALT matches from step 12. Preserve parent changes,
source/ALT provenance and existing anchored behavior. Export a separate diagnostic
catalog and keep physical fitting and production phasing unchanged. No new gap
closure is claimed.

## Composition contract

`AlleleReferenceMap` now retains its validated uppercase REF and ALT strings,
binding cached matches to their exact inputs. `compose_allele_on_subpaths` takes
a returned map, the complete context's one-based start and explicit neighbor
edits. Callers prepare each parent once and reuse its map for all neighbor ALTs.
No source HP/PS, truth, category, channel or MSA flag selects an alignment or ALT.

The function first attempts the existing reference-validated anchored composer.
Valid and unsupported results retain their status and sequence. Only an overlap
can use the mapped fallback, which requires a complete map and independently
valid neighbor composition on the original reference. No-op entries impose no
genotype constraint, including where the parent differs from REF. Invalid DNA,
REF, bounds and length remain unsupported.

Every changed neighbor's entire raw REF must fit inside one mandatory matched
subpath. Equal-length replacements can touch its edges. Length-changing edits
conservatively require a matched base on both sides of the raw REF. A matched
base alone does not determine insertion order beside a parent insertion: for
example REF `AC`, parent `TAC` and neighbor `AC -> GAC` must not silently choose
`TGAC` over `GTAC`. The focused regression remains unresolved.

Translate accepted raw edits using certified ALT offsets, then run the existing
composer on the complete parent sequence. This revalidates their REF and keeps
its anchored duplicate, ambiguous-padding, shifted-alias, conflicting insertion
and overlap checks. Multi-edit fixtures check permutation, exact/padded duplicate
aliases, mixed SNP/indel composition and contradictory edits. The primitive
supports explicit neighbor sets; this diagnostic extension enumerates one
neighbor ALT per attempt. The existing multi-neighbor physical-site choices
and atomic source-alias coupling are untouched.

The existing 1,048,576-cell mapping and 4,096-base DNA bounds remain. A limited
map may retain an already supported anchored result but cannot authorize the
mapped fallback. No partial map is treated as a placement certificate. Cache
state is local to each context; no new shared mutable reference state is added.
These mathematical sequence hypotheses are not calibrated genotype or phase
constraints.

| Table | Meaning |
|---|---|
| `.nested-compositions.tsv` | All original parent/neighbor/source-ALT attempts with mapped-fallback status and complete sequence |
| `.nested-alleles.tsv` | Sorted unique union of the entire existing path catalog and valid nested sequences |

Old matrix, composition, path, mapping, cohort, cost, fit and source tables retain
their previous behavior and bytes. Original fitted allele indices remain separate
from catalog indices. Normal execution does not invoke this diagnostic extension.

## Independent and fast verification

All **14,713** context/composition/path/cohort checks pass. New fixtures cover
compound islands, substitutions at either matched edge, guarded I/D, mixed
SNP/indel permutations, duplicate and padded aliases, no-op non-constraints,
contradictory/masked edits, invalid raw inputs, parent ALT offset shifts, limited
mapping, same-boundary insertions and ambiguous repeat SNPs. Five hundred random
compound-parent substitutions are compared with independent exhaustive
optimal-traceback intersections. All **12,069** genotype checks, existing unit
tests and **47 predicate cases / 1,551 assertions** pass without new warnings.
HiFi/ONT golden outputs and thread-determinism gates pass.

`audit_nested.py` independently recounts all optimal alignments using arbitrary
precision, validates every mandatory span, and uses an independent anchored-edit
oracle for every legacy and nested attempt. It checks every previous owning
matrix against frozen step-12 hashes and bytes. The audit takes **20.0 seconds**.
All **16,176 attempts**, **734 original parents** and **70 physical contexts** agree.
All previously valid compositions retain their exact output.

| Owner | Attempts | Newly valid nested compositions | Added unique sequences | Final hypotheses |
|---|---:|---:|---:|---:|
| Default 4–5 Mb | 1,980 | 7 | 0 | 516 |
| Default 5–6 Mb | 68 | 0 | 0 | 91 |
| Default 65–66 Mb | 13,951 | 0 | 0 | 12,509 |
| Whole-snarl 4–5 Mb | 177 | 0 | 0 | 234 |

Seven of step 12's **30** whole-REF-contained overlap candidates are resolved;
**23** remain unresolved under the conservative length-change guard. The seven
are physical context 18's parents 0/1/2/4/5/6/7 with candidate 1094 at
**4,760,662, A -> C**, an MSA-verified injected BAM SNP. Source labels annotate
these attempts without influencing admission. The verified substitution touches
the beginning of a mandatory matched subpath; no extra indel boundary rule is
needed for its fixed-length replacement. All parent changes are retained.

The resulting sequences already exist through different parent/neighbor paths
in step 11's catalog. For example, the sequence for parent 0 plus candidate 1094
also occurs through reference parent 3 plus candidates 1091/1094 or 1091/1093.
New matched-parent provenance does not create another molecule, independent vote
or unique sequence. This step resolves a representation limitation, not a new
sequence-catalog deficit or production gap.

`audit_catalogs.py` proves each nested catalog equals its existing path catalog
byte-for-byte, including indices, and verifies all original cohorts/costs/fits
against prior manifests. Deterministic edit costs on identical sequences and
queries cannot change. Previously verified per-read path residuals remain
**712 / 38 / 6,455 / 323** over **1,796 / 413 / 425 / 1,337** read/context rows.
No redundant catalog scoring is required. The 65 Mb context 3's length deficit
from step 12 remains: 4,327 of 4,383 residual edits are unavoidable with its
55-base maximum hypotheses. Missing read-supported long sequences still need
representation and calibrated joint inference.

`test_allele_context --nested INPUT_FOLDER OUTPUT_PREFIX` rebuilds all ten
composition/path/map/nested tables from five raw input tables, without BAM or
precomputed outputs. Every output matches pipeline bytes. Missing raw site state
is rejected. Over 100 runs under acceptance load, reconstruction/export averages
**39.6 / 8.5 / 376.2 / 16.2 ms** across the four owners; focused fixtures average
**140.3 ms**. Parent alignment is computed once per map rather than per attempt.

Eight actual production mutations are caught: losing parent ALT offsets, ignoring
either boundary guard, ignoring map completeness, losing parent changes, treating
no-ops as REF constraints, selecting arbitrary repeat matches and discarding
supported anchored results. See `nested-checks.json`, `boundary-checks.json`,
`catalog-checks.json`, `mutation-checks.json`, `raw-state-checks.json` and
`fast-checks.json`.

## Reproduction

```bash
make -j10 pgphase unit-tests gap-dev-check
make check
python3 evaluations/2026-10-07-representation-step13/audit_nested.py
python3 evaluations/2026-10-07-representation-step13/audit_boundaries.py
python3 evaluations/2026-10-07-representation-step13/audit_catalogs.py
python3 evaluations/2026-10-07-representation-step13/verify_raw_state.py
python3 evaluations/2026-10-07-representation-step13/verify_fast.py
python3 evaluations/2026-10-06-representation-step1/check_windows.py \
  --binary test_data/tmp_representation_step13/candidate-pgphase \
  --out test_data/tmp_representation_step13/accepted-final/windows
bash evaluations/2026-10-07-representation-step13/replay_full.sh \
  test_data/tmp_representation_step13/accepted-final/full-frozen \
  test_data/tmp_representation_step13/candidate-pgphase
PGPHASE_BIN=test_data/tmp_representation_step13/candidate-pgphase make window-tests
```

Owning scripts remain step-5 `replay_matrix.sh`/`replay_whole_snarl.sh` with the
explicit frozen binary, original regions and fresh outputs. Full/owner/panel
parental audits use bench-phasers Python with pysam. The HiPhase helper rechecks
input/index/truth/panel/helper identities before reusing previous measurements.

## Production and HiPhase acceptance

Final frozen binary SHA-256:
`5459c5d8b7a755403919f9361f4e8b119dbe9c162cf1d11c83db61d74d767484`.
All **107/107** registered checks pass in **717.8 seconds** of cold
acceptance. The official suite using the same frozen executable and replay cache
passes **18,138 assertions in five cases**, plus all **four** cache tests. Source
and binary were not changed during acceptance. Cold pipeline runs belong to final
acceptance; development uses the measured raw-state replay timings.

All four owning candidate tables, VCF rows, source/read-channel state, primary
HP/PS and parental classifications match step 12. The explicit-binary chromosome
run preserves **82,287** candidates, **64,483** VCF rows and all **256,612** primary
read assignments and parental statuses, including connected core and rescues.
Correct/discordant/unphased remain **230,965/6,620/19,027**. Its **250** blocks retain
N50 **955,496 bp** and largest **3,012,193 bp**.

All **141** native/full panel metrics and HiPhase classifications match step 12.
Of **118** HiPhase >=80%-correct windows, **59** still pass and **59** still fail
the full closure contract. HiPhase input/index/truth/panel/helper fingerprints
are rechecked before prior measurements are reused. No floor, expectation,
panel or closure certificate is refreshed; no new production closure is claimed.
New closures still require >=80% correct original primary truth-scorable reads,
with abstentions in the denominator, and both total and connected-core correct
counts at least HiPhase; rescue phase sets do not count as connected core.

`manifest.json` freezes source, inputs and accepted evidence hashes. All original
matrix tables are compared against hash-protected step-12 evidence; catalog/
cohort/cost invariance also checks the earlier residual verification hashes.
