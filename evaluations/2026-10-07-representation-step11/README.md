# Representation step 11: physical multi-neighbor sequence paths

Add a bounded multi-neighbor path builder over the reference-validated composer,
with one ALT choice per identical physical neighbor table. Retain all source
aliases and original parent indices as provenance. The pipeline exports the
expanded catalog in matrix diagnostics; it does not fit it or change read phase.
No new gap closure is claimed.

## Path and physical-site contract

Each original parent sequence is combined with omitted-site or selected-ALT
choices at the scoped neighbor sites already saved in step 10. Omission means
no edit, not a literal REF genotype. Raw neighbor POS/REF and the sorted unique
complete ALT table define physical identity. Source channel, candidate order,
ALT order and redundant ALT entries do not create a second physical choice.
Different bounds, REF strings or full alternative tables remain separate.
Duplicate candidate IDs are unsupported.

A selected physical ALT is applied once and maps back to the first matching
original ALT index of every source alias. All raw entries remain in the original
site/provenance tables. Paths preserve parent index and sorted source candidate/
ALT selections. Source labels are not extra independent observations or votes.
Final sequence hypotheses are sorted and deduplicated separately from paths.

The depth-first search prunes an unresolved-overlap prefix because adding edits
cannot remove its existing conflict under the conservative anchored composer.
Unsupported prefixes continue: an insertion exceeding the sequence length bound
can become supported after a later deletion. Every terminal unsupported path is
counted. The search visits at most **65,536 prefixes** and recurses over at most
**64 physical neighbor sites** per context. Hitting a bound returns `limited`
and clears every partial path. The existing step-10 catalog is retained;
truncation cannot become a complete genotype hypothesis table.

`complete` means the conservative edit search finished. It does not certify
biological haplotypes, graph-walk compatibility, read confidence or resolved
internal alignment of compound parents. Ambiguous repeat placement and edits
inside an undecomposed parent replacement can still remain unresolved.

| New table | Meaning |
|---|---|
| `.path-status.tsv` | Complete/limited/unsupported, visited and overlap-pruned prefixes, unsupported terminal paths and emitted path count |
| `.sequence-paths.tsv` | Original parent index, source candidate/ALT selections and full sequence |
| `.path-alleles.tsv` | Entire step-10 catalog plus unique valid paths from complete searches |

Original physical allele indices, costs, full/held-out fits and all previous
matrix tables stay unchanged. Path-catalog indices are separate from fitted
original allele indices. No parental truth, source HP/PS, genotype or MSA flag
admits a neighbor or chooses its ALT.

## Duplicate-neighbor regression found during this step

The first builder enforced one ALT per candidate description. The 65 Mb context
4 had **7,286,424 raw neighbor combinations per parent**, reached its work limit
and produced no complete path catalog. Its five raw sites included two exact
physical alias pairs:

* Candidates 1752/1753: POS 65,769,539, REF length 5, 28 ALTs each.
* Candidates 1757/1758: POS 65,769,573, REF length 1, 37 ALTs each.

Both pairs have identical complete raw tables. Coupling them reduces five source
choices to three physical choices and **6,612 combinations per parent**. The
corrected search completes in 20,930 visited prefixes and emits 19,837 paths.

This also prevents a correctness defect, not just redundant work. Two aliases
of `AA -> {CA, AC}` could independently select different ALTs and invent `CC`.
The regression reproduced three failing assertions before the fix. The corrected
catalog has three choices (omitted, CA, AC) and maps each selected sequence back
to both original ALT-table orders. Different full tables at the same raw bounds
are not silently merged. See `alias-preflight-checks.log`, the focused fixtures
and the physical-identity mutation checks.

Pre-dedup owning, full and registered runs are preliminary evidence only. The
accepted final outputs live in `test_data/tmp_representation_step11/accepted-final/`.
An inherited full replay script also ignored an extra binary argument; that
incorrect-binary run is excluded. `replay_full.sh` requires an explicit binary,
and final acceptance uses the frozen `candidate-pgphase` executable. Original
worker outputs were allowed to finish and their directories were not reused.

## Independent complete-path verification

`audit_paths.py` checks all prior matrix tables against immutable step-10 hashes
and exact bytes, then independently groups the raw complete site tables. Its
Python oracle enumerates maximal shared-padding trims, overlap corridors,
duplicate anchored edits and isolated repeat aliases before rendering in original
reference coordinates. It independently reproduces the search counters and every
path's ordered provenance and sequence.

For all completed searches it additionally enumerates **every full assignment
without pruning**, proving no supported path disappeared through an intermediate
conflict or length failure. It checks both extreme maximal padding placements
of every selected edit on every valid path. All **62,577 full assignments** and
**24,260 emitted paths** agree. The independent audit takes about **4.2 seconds**.
Every raw site, previous hypothesis, physical cohort, cost, fit and source state
is unchanged from hash-protected step 10; original input metadata are verified.

| Owner | Complete contexts | Emitted paths | Paths selecting >=2 physical neighbors | Added unique sequences | Final hypotheses |
|---|---:|---:|---:|---:|---:|
| Default 4–5 Mb | 31/31 | 711 | 182 | 89 | 516 |
| Default 5–6 Mb | 6/6 | 101 | 0 | 0 | 91 |
| Default 65–66 Mb | 10/10 | 23,203 | 20,960 | 10,929 | 12,509 |
| Whole-snarl 4–5 Mb | 23/23 | 245 | 0 | 0 | 234 |

All **70/70** searches complete. The catalog adds **11,018 unique sequences** to
step 10's 2,332 hypotheses, yielding **13,350**. These counts span overlapping
comparison arms and contexts; they are not independent genome-wide loci, reads
or accepted phase constraints. Multiple source selections may describe one
physical edit, so the table counts physical neighbor groups rather than alias
labels.

`test_allele_context --paths INPUT_FOLDER OUTPUT_PREFIX` rebuilds both the legacy
one-neighbor tables and all path tables from the five raw input tables: complete
reference contexts, raw global sites, parent sequences, physical membership and
parent candidate metadata. All six replay outputs equal production bytes.
Raw-only replay has no precomputed path outputs and rejects a missing raw site
catalog. No BAM or pipeline run is needed. Over 100 runs during acceptance load,
reconstruction and export average **39.3/14.5/388.4/15.4 ms** for the four owners;
all focused fixtures average **230.1 ms**. See `raw-state-checks.json` and
`fast-checks.json`.

All **11,613 context/composition/path/cohort checks** pass, including 100 randomized
four-site cases, exhaustive raw edit application and permutation comparisons,
atomic alias choices, original ALT index mapping, different full tables,
parent provenance, conflict pruning, exact/insufficient/zero budgets, depth
limits and a compensating deletion after an overlong prefix. All **12,069 genotype
checks**, existing unit tests and 47 predicate cases/1,551 assertions pass.
Twelve mutations of the actual production builder are detected, including
independent alias choices, ignoring the full ALT table, unsupported-prefix
pruning, publishing truncated paths, lost backtracking/provenance and unstable
source ordering. See `mutation-checks.json`.

## Same-cohort residual measurements and limits

The evaluation-only scorer uses original physical query slices. Every original
cost vector equals the independently verified baseline. The path catalog must
retain the one-neighbor catalog, so its per-read minimum cannot increase.

| Owner | Physical read/context records | Step-10 best-sequence residual | Path residual | Improved records |
|---|---:|---:|---:|---:|
| Default 4–5 Mb | 1,796 | 803 | 712 | 66 |
| Default 5–6 Mb | 413 | 38 | 38 | 0 |
| Default 65–66 Mb | 425 | 6,670 | 6,455 | 39 |
| Whole-snarl 4–5 Mb | 1,337 | 323 | 323 | 0 |

This saves **306 edit mismatches across 105 records** beyond one-neighbor
composition. Per-read best-sequence costs are unconstrained lower bounds: reads
may choose different haplotypes, and a larger catalog can fit sequencing errors.
These are not diploid fits, calibrated confidence or phase assignments.
The difficult 65 Mb context 3 still retains its high residual and original
4,404 diploid cost. Its single neighbor and unresolved internal parent/repeat
placement are not resolved by enumerating more site combinations.

Next verify internal REF-to-ALT matched subpaths and repeat placement, then use
physical reads and graph paths to constrain joint inference and calibrate errors.
A large raw catalog must not be passed blindly to a quadratic diploid fitter.
Source aliases, expanded hypotheses and overlapping contexts cannot count as
additional independent evidence. No truth-derived threshold is installed.

## Production and HiPhase acceptance

Final binary:
`b94bc35b4fba83829f7c4f5f14512347b4fbc0fd61d18c328e985065c9b93527`.
Final build and fast suites pass without new warnings. HiFi/ONT golden and
thread-determinism checks pass. Four final owning outputs preserve candidates,
VCFs, source matrices, every primary HP/PS and parental status. Identical original
physical cohorts retain the same 1,683 distinct molecules whose input alignments
were verified against HiPhase in earlier steps. HiPhase panel fingerprints and
committed measurements remain unchanged.

The final explicit-binary chromosome run preserves 82,287 candidates, 64,483
VCF rows and all **256,612** primary HP/PS assignments and parental statuses,
including connected-core and rescue tags. Correct/discordant/unphased remain
**230,965/6,620/19,027**. Its 250 blocks retain N50 **955,496 bp** and largest
**3,012,193 bp**. See `full-checks.json`.
All **107/107** registered checks pass in **788.3 seconds** in final cold
acceptance. All **141** native/full panel metrics and HiPhase classifications
match step 10. Of 118 windows with HiPhase >=80% correct original scorable reads,
59 still pass and 59 still fail the full closure contract. The official cached
suite with `PGPHASE_BIN` set to the same frozen executable passes **18,138
assertions in five cases**, plus all four replay-cache tests. See
`panel-comparison.json` and `standard-window-tests-final.log`.
Replay cache signatures include the literal executable token, so an identical
binary at `./pgphase` misses the frozen-path cache. An extra default-path cold
suite was stopped after the frozen-path suite passed; its incomplete outputs
are excluded from acceptance.
No floor, expectation, panel or closure certificate is refreshed.
New production closures still require >=80% correct original primary scorable
reads with abstentions in the denominator, and both total and connected-core
correct counts at least HiPhase.

## Reproduction

```bash
make -j10 pgphase unit-tests gap-dev-check
make check
./test_allele_context --paths \
  test_data/tmp_representation_step11/accepted-final/matrix4-final \
  test_data/tmp_representation_step11/cached/matrix4-final
# Repeat for matrix5-final, matrix65-final and whole-final.
python3 evaluations/2026-10-07-representation-step11/audit_paths.py
python3 evaluations/2026-10-07-representation-step11/verify_raw_state.py
python3 evaluations/2026-10-07-representation-step11/verify_fast.py
python3 evaluations/2026-10-06-representation-step1/check_windows.py \
  --binary test_data/tmp_representation_step11/candidate-pgphase \
  --out test_data/tmp_representation_step11/accepted-final/windows
bash evaluations/2026-10-07-representation-step11/replay_full.sh \
  test_data/tmp_representation_step11/accepted-final/full-frozen \
  test_data/tmp_representation_step11/candidate-pgphase
PGPHASE_BIN=test_data/tmp_representation_step11/candidate-pgphase make window-tests
```

Owning scripts are the step-5 `replay_matrix.sh`/`replay_whole_snarl.sh`; pass
an explicit binary and the same original 4–5, 5–6 and 65–66 Mb regions. Use fresh
output directories for pipeline reruns. Compile `score_catalog.cpp` against
`src/allele_context.o`, `src/allele_identity.o`, `src/edlib.o` and `-lhts`, score
each final owner to `cached/OWNER.residuals.tsv`, and run `audit_residuals.py`.
Parental/full/panel audits need the bench-phasers Python with pysam installed;
run `audit_panel.py` and `compare_panel.py` after all registered checks finish.
`manifest.json` freezes final source, input identities and accepted evidence.
