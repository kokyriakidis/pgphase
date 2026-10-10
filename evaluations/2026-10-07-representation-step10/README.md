# Representation step 10: anchored parent/neighbor sequence composition

Add a reference-validated edit composer and a diagnostic catalog of complete
parent-plus-one-neighbor sequences. The catalog uses the same original BAM
cohorts as step 9. It retains all original hypotheses, all scoped raw neighbor
ALTs and each attempted composition's provenance and status. The added catalog
is not yet used by genotype fits or production phasing. No gap is claimed closed.

## Composition contract and repeat pitfalls

An input edit is raw one-based POS/REF/ALT within a complete original-reference
window. Validate every REF, including no-op inputs, and every DNA base before
composition. Bound windows and rendered sequences to 4,096 bases. A REF/no-op
input imposes no edit, rather than a literal reference-genotype constraint.

Trim common padding and apply edits in original reference coordinates. Exact
repeated anchored inputs and unambiguous padded duplicates count once.
Different insertions at one boundary, reference consumed by another edit,
shifted descriptions with the same isolated effect, and edits touching an
ambiguous padding corridor stay unresolved. Differently padded ambiguous edits
cannot be collapsed solely because one trimming placement yields equal keys.
No sequence is emitted for an unresolved overlap. This conservative result is
not proof that no compatible haplotype exists.

Two preflight failures shaped the contract:

* Independently normalizing/left-aligning neighbors can shift an insertion or
  deletion through a repeat and make a raw conflicting combination look valid.
  The composer now retains anchored coordinates; normalization of isolated
  variants is insufficient to establish their joint path.
* The independent owning audit found a parent and nested `C(T)n` insertion whose
  suffix-first trimmed keys matched. Prefix-first trimming could place the two
  edits separately and produce a longer sequence. Treating them as a duplicate
  would assert an unsupported path. Exact raw duplicates are distinguished from
  differently padded ambiguous descriptions, which remain unresolved.

For each physical parent context, enumerate every ALT of every other candidate
whose complete raw REF is strictly inside the aligned outer anchors. Exclude
candidate IDs belonging to the parent descriptions. Include all channels,
categories and ALT indices without source genotype, phase, parental truth or
MSA-verification admission. Render every parent hypothesis with each neighbor
ALT separately. This is a **one-neighbor catalog**, not an exhaustive combination
of all neighbor genotypes or an overlap-aware local assembly. Compound parents
remain one trimmed replacement; internal matched subpaths are not inferred.
Some compatible nested neighbors therefore remain unresolved until a verified
REF-to-ALT path can place them.

| Table | Meaning |
|---|---|
| `.site-catalog.tsv` | Complete raw candidate ALT metadata |
| `.composition-contexts.tsv` | Physical parent bounds and original reference |
| `.composition-sites.tsv` | Scoped neighbor ALTs and source annotations |
| `.compositions.tsv` | Every parent/neighbor attempt, selected ALT, status and valid sequence |
| `.composed-alleles.tsv` | Sorted unique union of old hypotheses and valid compositions |

Physical IDs and original allele indices are unchanged. The expanded sorted
catalog has its own indices; they must not be interpreted as old fitted allele
indices. No expanded catalog is scored or fitted in the pipeline in this step.

## Independent sequence and saved-state verification

`audit_compositions.py` verifies raw reference against FASTA, complete raw source
catalogs against original candidate/parent tables, scoped neighbors against an
independently constructed census, and the complete Cartesian attempt list.
Every valid sequence is independently rendered in Python using both extreme
maximal padding trims of each input. It checks that these placements agree and
that every old hypothesis is retained. Unresolved overlaps may abstain
conservatively. All original `.joint-*` and `.physical-*` files remain
byte-identical to immutable, hash-verified step-9 evidence; original input
metadata also match. The audit completes in about 0.44 seconds.

| Owner | Contexts | Neighbor ALTs | Attempts | Valid | Unresolved | Added unique sequences |
|---|---:|---:|---:|---:|---:|---:|
| Default 4–5 Mb | 31 | 106 | 1,980 | 260 | 1,720 | 158 |
| Default 5–6 Mb | 6 | 10 | 68 | 10 | 58 | 0 |
| Default 65–66 Mb | 10 | 349 | 13,951 | 2,506 | 11,445 | 1,337 |
| Whole-snarl 4–5 Mb | 23 | 15 | 177 | 114 | 63 | 103 |

Across these comparison arms, 16,176 attempts yield 2,890 valid compositions,
13,286 unresolved overlaps and **1,598 added unique sequences**. Original
734 hypotheses become 2,332 hypotheses. These are physical-context catalog
counts, not independent genome-wide loci or additional phase evidence.

`test_allele_context --compose INPUT_FOLDER OUTPUT_PREFIX` rebuilds neighbors,
compositions and final hypotheses from the raw saved site catalog, parent
membership and reference. All three replay outputs equal production bytes.
It does not read BAM or rerun the pipeline. Over 100 runs during acceptance,
owner replays average **15.7/7.9/60.9/12.6 ms**; all focused fixtures average
**15.2 ms**. A separate raw-only replay with just the five input tables reproduces
all output bytes and rejects a missing raw site catalog. See `fast-checks.json`
and `raw-state-checks.json`.

All **3,191 context/composition/cohort checks** pass, including 1,000 independent
raw mixed-edit applications and their permutation/duplicate comparisons.
Fixtures cover reference mismatches, unknown DNA, bounds, simultaneous SNPs,
insertions and deletions, insertion boundaries, masked SNPs, repeat ambiguity
and both preflight failures. All **12,069 genotype checks**, existing unit tests
and 47 predicate cases/1,551 assertions pass. Eight mutations of the actual
production composer are caught: ignore REF validation, ignore padding ambiguity,
collapse ambiguous nested duplicates, count exact raw duplicates, count padded
duplicates, count shifted aliases, permit consumed-reference overlaps and fail
to consume deleted reference. See `mutation-checks.json`.

## Representation benefit and remaining inference work

The evaluation-only `score_catalog.cpp` uses the shared sequence-distance helper
on the same immutable physical cohorts. `audit_residuals.py` verifies every
original distance vector against the independently audited step-9 state and
checks the expanded minimum cannot exceed the old minimum.

| Owner | Molecule/context records | Old best-sequence residual | Expanded residual | Improved records |
|---|---:|---:|---:|---:|
| Default 4–5 Mb | 1,796 | 1,260 | 803 | 183 |
| Default 5–6 Mb | 413 | 38 | 38 | 0 |
| Default 65–66 Mb | 425 | 7,331 | 6,670 | 104 |
| Whole-snarl 4–5 Mb | 1,337 | 397 | 323 | 15 |

Best-sequence residual decreases by **1,192 edits over 302 records**. This is an
unconstrained per-read lower bound: each read may select a different hypothesis.
A larger catalog can also fit sequencing errors. It does not certify a diploid
model, confidence, genotype or parental phase. The difficult 65 Mb physical
context 3 still has 4,383 residual edits over 39 molecules (formerly 4,384 in
this unconstrained metric); its original diploid fit remains 4,404. One-neighbor
composition alone does not solve that representation/inference deficit.

Next construct compatible paths involving multiple neighbors while preserving
ambiguous overlaps and source provenance, then independently verify those paths
before fitting or consuming phase links. Resolve repeat placement with physical
sequence/path evidence rather than arbitrary isolated normalization. No
truth-derived threshold is installed.

## Production and HiPhase regression verification

Final production binary SHA-256:
`bfe9a26e258cf295fd5e6af9a8c5771f1401d88150ed9ec3811733f2d378903e`.
The signedness warning cleanup was checked to produce an exactly byte-identical
executable before finishing acceptance; no running comparison changed binaries.
Final build/fast tests introduce no warnings. HiFi/ONT golden and thread gates
pass. All four owning outputs preserve candidates, VCFs, source matrices, primary
HP/PS and parental status. The same 1,683 distinct molecules and unchanged raw
physical cohorts retain their previously verified identical input alignments
with HiPhase. HiPhase panel input/index/truth/panel/helper fingerprints pass and
the measurements equal the committed competitor table.

The full chromosome preserves 82,287 candidates, 64,483 VCF rows and all
256,612 primary HP/PS assignments, including connected-core and rescue tags.
Correct/discordant/unphased remain **230,965/6,620/19,027**. Its 250 blocks retain
N50 **955,496 bp** and largest **3,012,193 bp**. See `full-checks.json`.
All **107/107** registered checks pass in **718.5 seconds** for the one cold
acceptance run. All **141** native/full panel metrics and HiPhase classifications
match step 9. Of 118 windows where HiPhase correctly phases at least 80% of the
original primary scorable reads, 59 still pass and 59 still fail the full closure
contract. See `window-checks.json` and `panel-comparison.json`.
The standard cached suite passes **18,138 assertions in five cases** plus all
four replay-cache tests. See `standard-window-tests.log`.
No floor, panel, expectation or closure certificate is refreshed; future
production closures still require >=80% correct original scorable reads including
abstentions and both total and connected-core HiPhase parity.

## Reproduction

```bash
make -j10 pgphase unit-tests gap-dev-check
make check
```

Replay the four owners with the step-5 `replay_matrix.sh`/`replay_whole_snarl.sh`
scripts using 4–5, 5–6 and 65–66 Mb as documented in step 9. Accepted output is
`test_data/tmp_representation_step10/accepted/{matrix4-final,matrix5-final,matrix65-final,whole-final}`.
Use new output directories for new pipeline replays.

```bash
./test_allele_context --compose \
  test_data/tmp_representation_step10/accepted/matrix4-final \
  test_data/tmp_representation_step10/cached/matrix4-final
# Repeat saved-state composition for the other three owners.
/home/kokyriakidis/micromamba/envs/bench-phasers/bin/python \
  evaluations/2026-10-07-representation-step10/audit_compositions.py
python3 evaluations/2026-10-07-representation-step10/verify_fast.py
python3 evaluations/2026-10-07-representation-step10/verify_raw_state.py
python3 evaluations/2026-10-06-representation-step1/check_windows.py \
  --binary ./pgphase --out test_data/tmp_representation_step10/accepted/windows
bash evaluations/2026-10-06-shared-insertion-source/replay.sh \
  test_data/tmp_representation_step10/accepted/full-current
make window-tests
```

Build `score_catalog.cpp` against `src/allele_context.o`, `src/allele_identity.o`,
`src/edlib.o` and `-lhts`, then run it on each owner, saving output under
`cached/OWNER.residuals.tsv`; `audit_residuals.py` verifies those tables.
`audit_owners.py`, `audit_full.py`, `audit_panel.py` and `compare_panel.py` verify
production and full closure contracts against step 9. Run panel audits only
after the registered suite finishes. `manifest.json` freezes accepted source,
input identities and evidence.
