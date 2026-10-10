# Representation step 9: one model per identical physical context

Descriptions with identical reference bounds and complete allele sequences now
share one original-BAM cohort, one full allele-cost table and one diploid/held-out
fit. This removes the duplicated physical evidence found in step 8 while
preserving every original contrast, source observation and diagnostic projection.
The stage remains diagnostic; no candidate genotype, read phase or gap closure
is changed or claimed.

## Physical identity and provenance

`PhysicalAlleleContextKey` contains one-based inclusive reference bounds and the
sorted unique full sequence table. Keys are scoped to the owning contig/chunk.
Selected REF/ALT versus ALT/ALT pair, source channel, description order and
redundant allele-table entries do not define another physical observation.
Coordinates and complete alternative tables do: overlapping contexts or contexts
with different alternatives are not merged. Unsupported contexts have no
physical identity. This is exact deduplication, not neighborhood composition.

The pipeline builds all complete contexts, groups them in deterministic physical
key order, collects original BAM slices once per group, scores each molecule
against the complete table once, and fits one full and held-out model. The
existing eligible-primary/unique-name/coverage/quality-preservation rules are
unchanged. No source allele, source phase gauge or parental truth admits reads
or selects a fit.

New tables:

| Table | Meaning |
|---|---|
| `.physical-contexts.tsv` | Physical identity and complete reference bounds |
| `.physical-members.tsv` | Original contrast and selected full-table allele indices |
| `.physical-alleles.tsv` | Canonical complete sequence table |
| `.physical-cohort.tsv` | One raw sequence, quality, query-bound and MAPQ row per physical context/molecule |
| `.physical-costs.tsv` | One full distance vector per physical context/molecule |
| `.physical-genotypes.tsv` | One full fit, including ties and runner-up cost |
| `.physical-heldout.tsv` | One leave-one-molecule-out fit and allele decision per physical context/molecule |

The `.joint-*` tables are compatibility/provenance projections of shared state,
not additional independent evidence units. They remain **byte-identical** to
step 8. Each selected contrast retains its own full-table mapping and original
source annotations, so sharing a fit cannot import the wrong allele gauge.

## Independent verification

| Owner | Descriptions | Physical contexts | Previous records | Physical records | Removed repeats |
|---|---:|---:|---:|---:|---:|
| Default 4–5 Mb | 31 | 31 | 1,796 | 1,796 | 0 |
| Default 5–6 Mb | 7 | 6 | 484 | 413 | 71 |
| Default 65–66 Mb | 11 | 10 | 483 | 425 | 58 |
| Whole-snarl 4–5 Mb | 23 | 23 | 1,337 | 1,337 | 0 |

There are **70 physical contexts**, **3,971 physical molecule/context records**,
and **129 removed repeats**. Default 5 Mb descriptions 5/6 share all 71 molecules;
default 65 Mb descriptions 7/8 share all 58. Physical scoring saves **9,459**
repeated individual allele distances (40,568 retained rather than 50,027), as
well as the repeated full and held-out solves. These are counts of actual
physical evidence units; they are not additional reads or closures.

`audit_physical.py` independently constructs identities and member mappings in
Python, verifies that every unique cohort and cost vector equals the original
BAM state independently audited in step 8, and recomputes all **70 full** and
**3,971 held-out** unordered-pair fits with fresh NumPy reductions. Complete
genotypes, ties, runner-up costs, full-table nearest-allele membership and held-out
stability all match. It checks both physical and compatibility saved-state
replay tables byte-for-byte. The audit takes **2.1 seconds**. Immutable baseline
file hashes and original input metadata are verified before reuse; old scalar
DP and BAM extraction are not needlessly repeated.

**1,170 context/group/cohort checks** and **12,069 genotype checks** pass, plus
all unit tests, 47 predicate cases/1,551 assertions and HiFi/ONT golden and
thread-determinism gates. Grouping fixtures include 100 description/allele-order
permutations, REF/ALT versus ALT/ALT views, redundant hypotheses, unsupported
contexts, coordinate changes and different full alternative sets.

Six mutations to the actual production identity code are detected: ignore left
bound (one failed check), ignore right bound (two), ignore full alternatives
(two), use selected pair as identity (one then missing-key failure), lose alias
membership (77), and make every description independent (743). All seven existing
production fitter mutants are still detected. See `group-mutation-checks.json`
and `mutation-checks.json`.

`test_allele_genotype --state INPUT_FOLDER OUTPUT_PREFIX` consumes saved raw
physical cohorts, reconstructs grouping from complete contexts, fits each group
once and exports both physical and compatibility state. It requires the physical
cohort and explicitly rejects a legacy-only input. Over 100 runs during broader
verification, full scoring, full/held-out fitting and both table exports average
**37.1/13.5/52.5/19.8 ms** for the four owners; focused genotype fixtures average
**8.6 ms**. No pipeline or BAM read is required. See `fast-checks.json`.

## HiPhase and parental evaluation

No new molecules, alignments, physical hypotheses or genotype decisions enter
this stage. The physical union still contains the same **1,683 distinct
molecules** whose original CIGAR, SEQ, qualities, MAPQ, flags and coordinates
were verified against HiPhase in step 8. Original inputs and that verified
state remain unchanged. HiPhase panel measurements also pass their complete
input/index/truth/panel/helper fingerprint and equal the committed table.

`audit_parental.py` projects the cached independent step-8 held-out parental
separability onto unique physical units and checks that aliases agree. Exact
raw cohorts and full/held-out fits are verified first by `audit_physical.py`.
This is evaluation-only parental orientation using other molecules; these
predictions are not production HP/PS, calibrated confidences or closure claims.
Physical contexts and the whole-snarl comparison arm are not independent folds.

| Owner | Correct | Discordant | Unassigned | HiPhase correct on identical physical records |
|---|---:|---:|---:|---:|
| Default 4–5 Mb | 1,134 | 105 | 557 | 1,748 |
| Default 5–6 Mb | 82 | 5 | 326 | 334 |
| Default 65–66 Mb | 269 | 16 | 140 | 26 |
| Whole-snarl 4–5 Mb | 1,023 | 85 | 229 | 1,299 |

The numerical differences from step 8 remove repeated measurements; they do not
indicate lost phasing. The identical 65 Mb context retains its two-way fit tie;
selected source pairs cannot force another physical genotype.

The next step is explicit representation of compatible neighboring edits in
complete haplotype sequences. Exact physical deduplication does not yet explain
neighboring variation or establish a phase constraint. For example, the 65 Mb
context with residual cost 4,404 over 39 molecules remains unchanged. Preserve
incompatible alternatives and source provenance, then validate the composed
sequences before consuming links or calibrating errors. No truth-derived
threshold is installed.

## Production regression verification

The frozen final binary is
`a7dd51f5d76b93c4189e356d9dad89d312ebd53c8f2d00a766f6f80e1f453c51`;
the starting step-8 binary is
`d69176d69537b5f6030710ef0db97705c4b3f2bc783d9e1e8e626861abfb87d8`.

All four owning outputs preserve candidate rows, VCF calls, all original
`.joint-*`/source tables, every primary HP/PS and parental status. The full chromosome run also preserves 82,287 candidates, 64,483 VCF rows,
and all **256,612** primary HP/PS assignments. Correct/discordant/unphased remain
**230,965/6,620/19,027**; connected-core and rescue tags are unchanged. The
**250 blocks** retain N50 **955,496 bp**, largest **3,012,193 bp**. See
`full-checks.json`. All **107/107** registered checks pass in **716.6 seconds** for the one cold
acceptance run. The standard cached suite passes **18,138 assertions in five
cases** plus all four replay-cache tests. All **141** native/full panel metrics
and HiPhase classifications are identical to step 8. Of 118 windows where
HiPhase correctly phases at least 80% of the original primary scorable reads,
59 still pass and 59 still fail the full closure contract. See
`window-checks.json`, `standard-window-tests.log` and `panel-comparison.json`.
No expectation, floor, panel or closure certificate is refreshed.

## Reproduction

```bash
make -j10 pgphase unit-tests gap-dev-check
make check
```

Replay the original reference/BAM/sites/GAF with
`evaluations/2026-10-07-representation-step5/replay_matrix.sh`, using
`CHM13#0#chr20:4000001-5000000`, `5000001-6000000`, `65000001-66000000`, and
`replay_whole_snarl.sh` for the whole-snarl 4–5 Mb comparison. Accepted outputs
are in `test_data/tmp_representation_step9/accepted/` under `matrix4-final`,
`matrix5-final`, `matrix65-final` and `whole-final`. Use fresh output directories
for new replays.

```bash
./test_allele_genotype --state \
  test_data/tmp_representation_step9/accepted/matrix4-final \
  test_data/tmp_representation_step9/cached/matrix4-final
# Repeat cached replay for the other three owners.
/home/kokyriakidis/micromamba/envs/bench-phasers/bin/python \
  evaluations/2026-10-07-representation-step9/audit_physical.py
/home/kokyriakidis/micromamba/envs/bench-phasers/bin/python \
  evaluations/2026-10-07-representation-step9/audit_parental.py
python3 evaluations/2026-10-07-representation-step9/verify_group_mutants.py
python3 evaluations/2026-10-07-representation-step9/verify_fast.py
python3 evaluations/2026-10-06-representation-step1/check_windows.py \
  --binary ./pgphase --out test_data/tmp_representation_step9/accepted/windows
bash evaluations/2026-10-06-shared-insertion-source/replay.sh \
  test_data/tmp_representation_step9/accepted/full-current
make window-tests
```

`audit_owners.py`, `audit_full.py`, `audit_panel.py` and `compare_panel.py` verify
unchanged owning/full phases and the complete panel/HiPhase contract against
step 8. Run the panel audit after the registered suite finishes its contract.
`manifest.json` records the final source, input identities and accepted evidence.
