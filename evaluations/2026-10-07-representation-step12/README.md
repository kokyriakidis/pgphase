# Representation step 12: all-optimal parent/reference matched subpaths

Add a bounded REF-to-parent-ALT mapper before allowing nested edits inside
compound parent replacements. Export its evidence in matrix diagnostics and
rebuild it from saved raw state. The anchored composer, hypothesis catalogs,
diploid fits and production phasing retain their previous behavior. No new gap
closure is claimed.

## Contract and implementation

`map_allele_reference` computes global unit edit distance over uppercase A/C/G/T
strings of at most 4,096 bases each. Empty strings are valid. Its backward cost
matrix and rolling forward row identify every optimal edge consuming each
reference base. Each complete alignment consumes that base once; only a unique
optimal consuming edge whose bases match certifies its ALT offset. Consecutive
certified edges in both strings are coalesced into exact matched subpaths.

A fixed subpath survives ambiguity elsewhere. Substitutions, deletions and
matches with alternative repeat placements return no map for those positions.
The result does not select an arbitrary optimal traceback or assume that shared
prefix/suffix padding fixes repeat placement. Mathematical alignment invariance
is not biological compatibility, a genotype call or calibrated confidence.

The default bound is **1,048,576 cells**, including boundary row and column.
Insufficient budget returns `limited`; invalid DNA or excessive sequence length
returns `unsupported`. Both retain distance -1 and no subpaths. `complete` means
all-optimal analysis finished, even when no position has a certified match. The
length cap bounds arithmetic; the cell cap bounds quadratic work and memory.
All state is local to the call; no shared mutable reference handle is introduced.

| Table | Meaning |
|---|---|
| `.parent-maps.tsv` | Original physical context/parent indices, map status, distance, matched bases and subpaths |
| `.matched-subpaths.tsv` | Zero-based REF/ALT offsets and subpath length, in complete context strings |

No parental truth, source HP/PS, channel, candidate category or MSA flag chooses
an alignment. Original allele indices are retained. Normal execution does not
compute these maps; they are emitted only with matrix diagnostics.

## Independent verification and fast saved-state checks

All **14,184** context/composition/path/cohort checks pass, including 500 random
short sequence pairs checked against exhaustive optimal-traceback intersections.
Targeted fixtures cover internal islands in compound parents, both repeat
insertion/deletion directions, equal-cost mismatch/gap alternatives, empty
strings, case normalization, invalid DNA and exact/insufficient/zero work budgets.
All existing unit tests, **12,069** genotype checks and **47** predicate cases/
**1,551** assertions pass. Builds introduce no warnings. HiFi/ONT golden outputs
and thread-determinism gates pass.

`audit_maps.py` uses an independent full scalar DP and arbitrary-precision counts
of all optimal alignments. A matching edge is mandatory only if the count of
complete alignments through it equals the total. Every emitted subpath and cost
agrees for all original parents in four owning replays. Earlier matrices,
catalogs, cohorts, costs, full/held-out fits and source annotations remain byte
identical to hash-verified step 11. The independent audit takes **18.9 seconds**.

| Owner | Contexts | Parents, all complete | Parents with multiple optimal alignments | Certified matched bases | Subpaths |
|---|---:|---:|---:|---:|---:|
| Default 4–5 Mb | 31 | 269 | 232 | 7,189 | 400 |
| Default 5–6 Mb | 6 | 91 | 79 | 1,809 | 156 |
| Default 65–66 Mb | 10 | 243 | 220 | 13,422 | 326 |
| Whole-snarl 4–5 Mb | 23 | 131 | 98 | 3,597 | 180 |

All **734** maps complete over **70** contexts. Of these parents, **629** have
multiple optimal alignments. The mapper certifies **26,017** of 53,679 reference
bases in **1,062** subpaths; the remainder includes altered/deleted bases and
ambiguous matches. Counts include overlapping comparison arms, not independent
loci or evidence units.

`test_allele_context --maps INPUT_FOLDER OUTPUT_PREFIX` reconstructs eight output
tables from the same five raw inputs used in step 11: reference contexts, raw
sites, parent sequences, parent membership and parent metadata. Replay does not
read precomputed compositions, paths or maps. All eight outputs equal pipeline
bytes, and a missing raw site catalog is rejected. Over 100 executions, complete
reconstruction/export averages **34.9/9.1/311.3/17.0 ms** for the four owners.
Focused fixtures average **137.7 ms** under acceptance load.

Eight mutations of actual production mapping code are detected: admitting
nonmandatory matches, ignoring deletion alternatives, mapping mismatches,
shifting ALT offsets, losing subpath coalescing, ignoring the budget, losing case
normalization and corrupting insertion/deletion boundary costs. See
`mutation-checks.json`, `fast-checks.json` and `raw-state-checks.json`.

## What this establishes for the next composition step

Among previously unresolved parent/neighbor attempts, **30** in default 4–5 Mb
have the neighbor's whole raw REF inside one mandatory matched subpath. These
are concrete candidates for certified nested composition; they are not accepted
sequence hypotheses yet. All old composition/path tables deliberately remain
unchanged. Whole-REF containment is conservative and does not certify a new ALT,
its boundary placement, graph-walk compatibility or diploid support.

The default 65 Mb context 3 still has no formerly overlapping neighbor wholly
inside a certified subpath. Its residual requires checking hypothesis coverage
as well as alignment ambiguity. Every parent and path hypothesis is at most
**55 bases**, while **21 of 39** original query slices are longer, up to **407
bases**. Independent aligned-pair extraction reproduces every query boundary
and sequence; internal CIGAR insertions/deletions exactly explain each length.
All 39 primary input alignments also match HiPhase. No query slicing artifact
is needed to explain the observed length imbalance.

The unavoidable length-only lower bound is **4,327** of **4,383** best-path residual
edits (**98.7%**), independently bounded per read and compared with hash-verified
step-11 residuals. The original diploid cost is still **4,404**. All 11 scoped
neighbor ALTs are shorter than their 23-base REF, so nested composition among
these existing tables cannot supply the missing repeat expansion. The observed
long sequences require read-supported hypothesis construction with graph and
MSA evidence, followed by diploid/error calibration; catalog growth alone must
not become a phase constraint. This measurement establishes missing sequence
length in the catalog, not a biological expansion genotype. See
`length-deficit-checks.json`.

## Reproduction

```bash
make -j10 pgphase unit-tests gap-dev-check
make check
python3 evaluations/2026-10-07-representation-step12/audit_maps.py
python3 evaluations/2026-10-07-representation-step12/verify_raw_state.py
python3 evaluations/2026-10-07-representation-step12/verify_fast.py
/home/kokyriakidis/micromamba/envs/bench-phasers/bin/python \
  evaluations/2026-10-07-representation-step12/audit_length_deficit.py
python3 evaluations/2026-10-06-representation-step1/check_windows.py \
  --binary test_data/tmp_representation_step12/candidate-pgphase \
  --out test_data/tmp_representation_step12/accepted-final/windows
bash evaluations/2026-10-07-representation-step12/replay_full.sh \
  test_data/tmp_representation_step12/accepted-final/full-frozen \
  test_data/tmp_representation_step12/candidate-pgphase
PGPHASE_BIN=test_data/tmp_representation_step12/candidate-pgphase make window-tests
```

Owning scripts remain step-5 `replay_matrix.sh`/`replay_whole_snarl.sh`; use fresh
outputs, the explicit frozen executable and original 4–5, 5–6 and 65–66 Mb regions.
Parental/full/panel audits use the bench-phasers Python with pysam installed.
The HiPhase measurement helper verifies the saved signature against current
inputs, indices, truth, panel and helper identities before reusing measurements.

## Production and HiPhase acceptance

Final frozen binary SHA-256:
`b93796a43194c7c9d37c1ce39295bd89eca95b7e0c72a68b14341b1bdb6d526e`.
All **107/107** registered checks pass in **722.0 seconds** of cold
acceptance. The official suite using the same frozen `PGPHASE_BIN` and shared
replay cache passes **18,138 assertions in five cases**, plus all four cache tests.
Each saved-state assertion still runs. Cold pipeline replay belongs to acceptance,
while the mapper/composer development loop uses the measured raw-state timings.

All four owning runs preserve candidates, VCF rows, source/read-channel states,
primary HP/PS and parental classifications. The explicit-binary full chr20 run
preserves **82,287** candidates, **64,483** VCF rows and all **256,612** primary
HP/PS assignments and parental statuses, including connected core and rescues.
Correct/discordant/unphased remain **230,965/6,620/19,027**. Its **250** blocks retain
N50 **955,496 bp** and largest **3,012,193 bp**.

All **141** native/full panel metrics and HiPhase classifications match step 11.
Of **118** HiPhase >=80%-correct windows, **59** still meet the full closure
contract and **59** still fail it. HiPhase input/index/truth/panel/helper identities
are reverified before cached measurements are reused; the 65 Mb diagnosis also
checks all 39 relevant primary alignments directly against HiPhase. No floor,
panel, expectation or closure certificate is refreshed. No new closure is claimed.
Production closure still requires >=80% correct original primary truth-scorable
reads, with abstentions in the denominator, plus total and connected-core correct
counts at least HiPhase; rescue phase sets are excluded from connected core.

`manifest.json` freezes final source, inputs and accepted output hashes;
`baseline-evidence.json` records the step-11 evidence reverified before comparison.
No worker output or binary was replaced during acceptance.
