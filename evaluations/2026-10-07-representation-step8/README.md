# Representation step 8: collect the complete physical genotype cohort

The diagnostic genotype fit now receives every eligible, uniquely identified
original BAM molecule with complete parent-context coverage, independently of
existing graph/BAM allele calls. This fixes the censored input identified in
step 7. The change remains a representation/input stage: candidate genotypes,
phase joins and HP/PS are unchanged, and no new gap closure is claimed.

## Admission and saved state

`collect_allele_read_slices` consumes the original alignments already retained by
the owning whole-BAM solve. Unmapped, secondary and supplementary records are
excluded. The recovery MAPQ floor and configured QC/duplicate policy match BAM
loading. Multiple eligible primary records with the same molecule name are
ambiguous before physical coverage is checked; an uncovered second primary
cannot allow its covering counterpart to masquerade as a unique alignment.

Complete aligned outside anchors are required. Internal insertions/deletions,
original reference-oriented sequence, raw base qualities (including 255), MAPQ
and query bounds survive unchanged. Reference skips, incomplete anchors, unknown
bases and query slices exceeding 4096 bases abstain. Solver `is_skipped` state,
graph membership, source allele labels and parental truth do not admit molecules.

The new `.joint-cohort.tsv` uses the sequence-state schema, with a source allele
as annotation or `.` when there is no retained source observation. The old
`.joint-sequences.tsv` preserves contrast-observed provenance byte-for-byte.
Full allele costs, diploid fits and exact leave-one-molecule-out fits consume the
new cohort. The genotype objective is unchanged. `test_allele_genotype --state`
requires the new table; a legacy censored state fails explicitly rather than
silently becoming its input.

## Independent census and genotype checks

| Owning run | Previous records | Complete records | Added | Affected contexts | Changed fitted pairs |
|---|---:|---:|---:|---:|---:|
| Default 4–5 Mb | 1,714 | 1,796 | 82 | 4 | 2 |
| Default 5–6 Mb | 412 | 484 | 72 | 3 | 1 |
| Default 65–66 Mb | 386 | 483 | 97 | 8 | 5 |
| Whole-snarl 4–5 Mb | 1,337 | 1,337 | 0 | 0 | 0 |

These are molecule/locus records, not distinct molecules or independent folds.
`audit_cohort.py` independently fetches original primary BAM records and walks
aligned query/reference pairs to verify complete cohort membership and every
saved sequence, quality byte, query bound and MAPQ. All **4,100** records match;
there are no missing or ineligible extras. The **251** additions agree with the
pre-change census. Original scored rows and their source annotations are equal.

The same physical context and full hypothesis table must yield identical costs
and genotypes even when its original selected contrast differs. This now holds
for default 5 Mb loci 5/6 (both 71 molecules) and default 65 Mb loci 7/8 (both 58).
The latter previously used 18 and 39 molecules and selected different unique
pairs; both now correctly retain a two-way tie. Default 4 Mb locus 1 expands from
12 to 49 molecules and changes its unique pair. Three previously unique fits
become ambiguous. See `changed-cohort-fits.json` and `cohort-checks.json`.

`audit_genotypes.py` independently solves every full and held-out unordered pair
with fresh NumPy reductions, including homozygous pairs and ties. It verifies
all **11,510** added individual allele distances with independent scalar dynamic
programming. The **38,517** retained distances must equal the independently
verified step-7 vectors; they are not expensively recomputed. All four saved
allele/cost/genotype/held-out tables are byte-identical to production output.
The independent cached audit takes 15.3 seconds for this acceptance run.

**667** context/cohort checks and **12,069** genotype checks pass. Six mutations
to production cohort code fail: secondary admission, ignored MAPQ, ignored
filtered policy, duplicate revival, coverage before uniqueness and reference
skip admission (102/3/1/53/3/3 failing checks). All seven existing production
fitter mutants also fail. These mutate the actual implementation, not test-only
copies of its logic. See `cohort-mutation-checks.json` and `mutation-checks.json`.

For 100 executions under concurrent broad verification, saved full scoring,
fitting, held-out fitting and writing average **41.5/20.4/65.4/20.6 ms** for the
four owners. Genotype fixtures average **15.7 ms**. These loops require no BAM
read, pipeline rerun or parental truth. See `fast-checks.json`.

## HiPhase and remaining modeling work

`prepare_hiphase.py` freshly verifies identical original primary alignments for
all **1,683 distinct cohort molecules** in the frozen HiPhase output, including
CIGAR, SEQ, quality, MAPQ, flags, coordinates and reference lengths. Only the
existing `CHM13#0#chr20` / `chr20` alias differs. A single global parental
orientation is used; the limited older panel tag cache is not used for this
cohort. See `hiphase-molecule-checks.json`.

The following is **evaluation-only parental separability**, using other reads
for a held-out genotype's parental orientation. It is not production phasing,
a closure certificate, or a calibrated confidence estimate. Repeated
molecule/locus records and the whole-snarl arm are not independent observations.

| Owning run | Correct | Discordant | Unassigned | HiPhase correct on identical records |
|---|---:|---:|---:|---:|
| Default 4–5 Mb | 1,134 | 105 | 557 | 1,748 |
| Default 5–6 Mb | 123 | 5 | 356 | 405 |
| Default 65–66 Mb | 269 | 16 | 198 | 26 |
| Whole-snarl 4–5 Mb | 1,023 | 85 | 229 | 1,299 |

Uncensoring removes incompatible source-dependent fitted pairs, but cannot make
single-parent edit distances a general phasing model. For example, default
65 Mb locus 3 still has cost 4,404 across 39 molecules and a runner-up margin of
one. Neighboring variants and alternative local haplotypes must be represented
jointly; identical parent contexts must not become independent phase votes.
The next implementation stage should deduplicate these physical contexts and
construct explicit compatible neighboring edit sequences before adding phase
constraints or calibrating errors. No truth-derived threshold is installed.

## Production regression verification

The final production SHA-256 is
`d69176d69537b5f6030710ef0db97705c4b3f2bc783d9e1e8e626861abfb87d8`;
the starting step-7 binary is
`5746188fb9e6e737be9017bf39cefb11fd4ea1f88c1a72b6325490b720fd7bbe`.
All four final owning runs retain identical candidate rows, VCF calls, original
sequence/provenance tables and every primary HP/PS and parental status.

The full chr20 run also retains identical 82,287 candidates, 64,483 VCF rows,
and all **256,612** primary HP/PS assignments. Correct/discordant/unphased remain
**230,965/6,620/19,027**. All existing connected-core and rescue tags survive.
There are **250 blocks**, N50 **955,496 bp**, largest **3,012,193 bp**.
See `owner-checks.json` and `full-checks.json`.

All **107/107** registered window/mechanism checks pass in **722.0 seconds**
for the one cold acceptance run. The standard cached `make window-tests` passes
**18,138 assertions in five cases**, plus all four replay-cache tests. All 141
native/full panel contracts and HiPhase classifications equal step 7. Of the
118 HiPhase-at-least-80%-correct windows, 59 pass and 59 still fail the full
closure contract; this representation step does not claim a new closure.
HiPhase panel measurements reuse a signature-verified unchanged original/BAM,
competitor/index, truth, panel and helper state and equal the committed table.
See `window-checks.json`, `panel-comparison.json`, `hiphase-validation.log` and
`standard-window-tests.log`.
Floor, panel and original input identities are unchanged; no expectation is
refreshed. `baseline-checks.json` verifies input metadata and the immutable,
independently scored baseline used to avoid repeating scalar DP on old records.
`manifest.json` freezes final source, evidence and accepted states.

## Reproduction

Build and run focused checks with:

```bash
make -j10 pgphase unit-tests gap-dev-check
make check
python3 evaluations/2026-10-07-representation-step8/verify_cohort_mutants.py
python3 evaluations/2026-10-07-representation-step8/verify_fast.py
```

Owning inputs are the original chr20 BAM/reference/sites/GAF. Replays use
`evaluations/2026-10-07-representation-step5/replay_matrix.sh` at
`CHM13#0#chr20:4000001-5000000`, `5000001-6000000` and `65000001-66000000`,
and `replay_whole_snarl.sh` for 4–5 Mb, under
`test_data/tmp_representation_step8/accepted/{matrix4-final,matrix5-final,matrix65-final,whole-final}`.
Use fresh directories for a rerun.

```bash
./test_allele_genotype --state \
  test_data/tmp_representation_step8/accepted/matrix4-final \
  test_data/tmp_representation_step8/cached/matrix4-final
# Repeat saved-state replay for the other three owners.
/home/kokyriakidis/micromamba/envs/bench-phasers/bin/python \
  evaluations/2026-10-07-representation-step8/audit_cohort.py
/home/kokyriakidis/micromamba/envs/bench-phasers/bin/python \
  evaluations/2026-10-07-representation-step8/prepare_hiphase.py
/home/kokyriakidis/micromamba/envs/bench-phasers/bin/python \
  evaluations/2026-10-07-representation-step8/audit_genotypes.py --scalar-costs
/home/kokyriakidis/micromamba/envs/bench-phasers/bin/python \
  evaluations/2026-10-07-representation-step8/audit_genotypes.py \
  --production --work test_data/tmp_representation_step8/accepted
python3 evaluations/2026-10-06-representation-step1/check_windows.py \
  --binary ./pgphase --out test_data/tmp_representation_step8/accepted/windows
bash evaluations/2026-10-06-shared-insertion-source/replay.sh \
  test_data/tmp_representation_step8/accepted/full-current
make window-tests
```

`audit_owners.py`, `audit_full.py`, `audit_panel.py` and `compare_panel.py` compare
owning/full outputs and the complete HiPhase contract with step 7. Run the
panel audit only after the registered runner has finished its complete contract.
