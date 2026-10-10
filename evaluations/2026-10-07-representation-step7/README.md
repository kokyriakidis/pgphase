# Representation repair: step 7, full-allele cohort fitting and held-out evidence

The diagnostic sequence model now retains **every allele's cost for every
molecule in the saved contrast-observed cohort**, rather than collapsing unselected parent alleles to their
minimum distance. Full sequences have a sorted, unique table independent of
source ALT order and the originally selected contrast.

`fit_diploid_alleles` evaluates all unordered pairs, including homozygous pairs,
minimizing the sum of each molecule's edit distance to the nearer allele.
Repeated identical molecule vectors count once; conflicting vectors for the
same molecule remain unknown. Cohort sums are 64-bit. The model records minimum
and runner-up costs and the number of tied pairs. Ties retain no selected
pair; a missing second allele is not guessed from the first ALT in a list.
This objective is **not a likelihood model or calibrated genotype confidence**.

Every original molecule also has a leave-one-molecule-out fit. Pair totals are
computed once; subtracting that molecule's contribution gives exactly the fit
on all other molecules without realignment or repeated training solves.
An empty training cohort abstains. Held-out membership requires a unique nearest
sequence across **all** full hypotheses and membership in the independently
fitted pair. The record separately states whether that pair matches the full
cohort's pair. A lone supporting molecule cannot certify its own second allele.
Unselected competitors, sequence ties and unstable pairs retain uncertainty.

`--phase-matrix-dump` adds `.chunkN.joint-alleles.tsv`, `.joint-costs.tsv`,
`.joint-genotypes.tsv` and `.joint-heldout.tsv`. Original complete contexts,
query sequences, qualities, physical bounds, source reductions and channel
provenance remain byte-identical to step 6. Production fitting sees no parental
truth or source HP/PS labels. With diagnostics disabled, the fitter does no work.
The candidate's consensus genotype, admission, original graph/BAM profiles,
component joins, read tags and output remain unchanged.

## What fitting the original molecules shows

| Replay | Contexts | Different unique pairs | Ambiguous full fits | Stable for every held-out molecule | Inferred ALT/ALT fits |
|---|---:|---:|---:|---:|---:|
| Default 4–5 Mb | 31 | 10 | 0 | 27 | 8 |
| Default 5–6 Mb | 7 | 2 | 1 | 6 | 0 |
| Default 65–66 Mb | 11 | 6 | 1 | 9 | 2 |
| Whole-snarl 4–5 Mb | 23 | 2 | 0 | 23 | 2 |

These are diagnostic fits, **not newly accepted genotypes or gap closures**.
Across 72 contexts there are 70 unique full fits; 20 differ from the original
contrast, and 65 are unchanged by holding out every molecule individually.
The fitted pair can include complete parent alleles absent from the original
REF/ALT contrast, including ALT/ALT. For example, default 4–5 Mb locus 4 reduces
cohort cost from 53 to 15 with a different pair, stable over all 68 held-out
molecules. Default 5–6 Mb loci 5 and 6 select the same full pair, despite their
different original selected ALTs. Another 5–6 Mb context has two equally good
pairs and remains ambiguous.

Among 111 source conflicts in the two default repeat-heavy owners, 67 are
uniquely closest members of a stable held-out pair; 44 remain outside that pair
or sequence-tied. Whole-snarl mode's 51 conflicts split 33/18. None becomes a
phasing vote at this stage.

## Why best fit and stability are insufficient confidence

Parental truth is used only in the evaluation, to measure whether a fitted split
separates parental molecules. Each evaluated read is held out of both genotype
fitting and the evaluation's allele-orientation vote. A parental orientation is
available only when the other molecules give distinct, untied parental-majority
labels for both alleles. This is an **evaluation-only separability experiment**,
not the accuracy of a production phaser or an externally known allele genotype.

The four owners have experimental correct/discordant/unassigned counts of
1,094/105/515, 130/25/257, 257/6/123 and 1,023/85/229 over their scored
molecule/locus records. Reads repeated across loci or modes are repeated records;
these totals are not counts of distinct genome-wide reads. Corresponding HiPhase
correct counts on those same records are 1,666, 333, 26 and 1,299. The 65–66 Mb
HiPhase count is measured on the actual cached locus molecules, not extrapolated
from its other blocks or from the gap panel.

Even exact full-sequence matches with minimum BQ >=40 and MAPQ >=40 have 16
experimental discordances among 439 stable-member records (413 correct, ten
unassigned). Larger nearest-allele margins correlate with better separability:
margin >=2 has 1,265 correct, one discordant and 17 unassigned records. These are
observational bins, not a chosen production threshold. Default and whole-snarl
modes share reads, so they cannot serve as independent calibration partitions.

Some loci have a stable best pair but a very large absolute residual, indicating
unmodeled sequence differences in the window. For example, the 65–66 Mb locus 3
has cost 4,331 over 25 molecules. Other sites fit exactly but still have mixed
parental groups. Source preference, MSA provenance, a minimum whole-slice BQ,
or held-out cost stability alone cannot establish phase confidence. The next
consumer first needs an uncensored molecule census, then joint neighboring-site/
graph constraints and a validated error model rather than converting these fits directly into read HP or block joins.
No truth-derived threshold has been installed in production.

An independent original-BAM coverage audit also finds a representation input
limitation. Current cached states admit molecules through the original contrast's
retained source observations. In the three default owners, **15 contexts omit
251 eligible, fully covered molecule/locus records** (82/72/97). Whole-snarl mode
omits none in this audited state. The audit applies the whole-BAM recovery MAPQ
floor, primary/flag filters, aligned outer anchors, no reference skip, known DNA
and the existing 4,096-base slice bound; it is not just a bounding-box count.
For example, default 4–5 Mb locus 19 retains 26 of 69 eligible molecules, and
5–6 Mb loci 5/6 retain 49/38 of the same 71 eligible molecules. Thus even the
same complete sequence universe can be fitted on different source-conditioned
cohorts. These fits cannot certify complete sample genotypes. The next input
step must obtain costs for every eligible fully covered molecule independently
of original allele-call admission, including originally unknown calls, before
using the model for phase inference. See `cohort-coverage-checks.json` and
`audit_cohort_coverage.py`.

See `heldout-parental-audit.tsv`, `calibration-checks.json`,
`cached-genotype-checks.json` and `production-genotype-checks.json`.

## Verification and fast development

The focused binary passes **12,069 checks**, including 1,000 independent exhaustive
matrix oracles, fresh held-out training fits, allele-index/read-order permutation,
complete ALT/ALT recovery, homozygous hypotheses, ties, empty training cohorts,
unselected competitors, duplicate/conflicting molecules and 64-bit sums.
The existing 545 physical-context checks and 47 phasing predicate cases
(1,551 assertions) still pass.

An independent scalar dynamic-programming oracle verifies all **38,517 individual
allele distances** across **3,849 original molecule/locus slices**. The one-time
audit takes 117.1 seconds. Independent NumPy exhaustive fits verify every full
and held-out model. The production dump matches cached scorer output byte for
byte. A replay serializer ordering loci as text was caught by this comparison;
it now uses the same numeric locus order as production.

Seven actual production mutations fail: retaining the first tied pair, training
on the held-out read, counting duplicate molecules, overwriting a conflicting
molecule, forcing heterozygous pairs, overflowing cohort sums and ignoring
other read hypotheses. See `mutation-checks.json`.

Over 100 executions including startup and concurrent broad-test load, focused
fixtures average **8.28 ms**. Complete cached scoring, fitting, held-out validation
and writing the four state tables average **28.87/13.38/49.89/18.28 ms** for
4–5/5–6/65–66 Mb and whole-snarl 4–5 Mb. No BAM access or pipeline replay occurs.
The final timing state and binary identity are in `fast-checks.json`.

Fresh HiPhase auditing verifies **1,568 distinct original primary molecules**
against its BAM, including alignment bounds, CIGAR, sequence, quality bytes,
MAPQ and flags. Its legacy header shortens `CHM13#0#chr20` to `chr20`; all
corresponding contig lengths and header order are equal. Its allele gauges are
oriented globally on the complete competitor output, not separately by locus.
The committed 141-window competitor measurements also pass input/index/truth/
panel/competitor/helper identity validation. See `hiphase-molecule-checks.json`
and `hiphase-verified.tsv(.json)`.

All four owning candidate/read-channel/quality matrices, original sequence
contexts/provenance tables, candidate/VCF rows, primary tags and parental states
remain identical to step 6. Default correct/discordant/unphased counts remain
3,862/26/238, 4,001/23/232, 2,962/6/3; whole-snarl 3,866/27/235 under its own
unchanged flags. See `owner-checks.json` and `contrast-checks.json`.

Starting production SHA256:
`e6a9aae89bdbe011be7d7759f67efcfe1f29b9c7c77d2c12298a867f37fdde68`.
Final production SHA256:
`5746188fb9e6e737be9017bf39cefb11fd4ea1f88c1a72b6325490b720fd7bbe`.

## Reproduction

Run from the repository root; the benchmark Python contains pysam and NumPy.

```sh
make -j10
make unit-tests gap-dev-check
./test_allele_genotype --state \
  test_data/tmp_representation_step6/accepted/matrix4-final \
  test_data/tmp_representation_step7/cached/matrix4-final
python3 evaluations/2026-10-07-representation-step7/verify_fast.py
/home/kokyriakidis/micromamba/envs/bench-phasers/bin/python \
  evaluations/2026-10-07-representation-step7/prepare_hiphase.py
/home/kokyriakidis/micromamba/envs/bench-phasers/bin/python \
  evaluations/2026-10-07-representation-step7/audit_genotypes.py --scalar-costs
bash evaluations/2026-10-07-representation-step5/replay_matrix.sh \
  test_data/tmp_representation_step7/accepted/matrix4-final
bash evaluations/2026-10-07-representation-step5/replay_matrix.sh \
  test_data/tmp_representation_step7/accepted/matrix5-final ./pgphase CHM13#0#chr20:5000001-6000000
bash evaluations/2026-10-07-representation-step5/replay_matrix.sh \
  test_data/tmp_representation_step7/accepted/matrix65-final ./pgphase CHM13#0#chr20:65000001-66000000
bash evaluations/2026-10-07-representation-step5/replay_whole_snarl.sh \
  test_data/tmp_representation_step7/accepted/whole-final
/home/kokyriakidis/micromamba/envs/bench-phasers/bin/python \
  evaluations/2026-10-07-representation-step7/audit_genotypes.py \
  --production --work test_data/tmp_representation_step7/accepted
make check
python3 evaluations/2026-10-06-representation-step1/check_windows.py \
  --binary ./pgphase --out test_data/tmp_representation_step7/accepted/windows
make window-tests
bash evaluations/2026-10-06-shared-insertion-source/replay.sh \
  test_data/tmp_representation_step7/accepted/full-current
```

The four cached owner states must exist before auditing production; generate each
with `test_allele_genotype --state` as above. After completed replays, run the local
`audit_owners.py`, `audit_contrasts.py`, `audit_full.py`, `audit_panel.py` and
`compare_panel.py`. `manifest.json` records original input metadata, extraction/
fitting sources, frozen evidence states and unchanged floors and panel hashes.

## Final regression and chromosome checks

Build, units, phasing predicates, HiFi/ONT TSV/VCF goldens and HiFi t1/t4
determinism pass without new compiler warnings. All **107/107** registered
window/mechanism checks pass in **716.4 seconds** cold. The standard cached
`make window-tests` passes **18,138 assertions in five cases** and all four
replay-cache tests.

The final full chr20 replay preserves all 82,287 candidate rows, 64,483 VCF
records and 256,612 primary HP/PS and parental states: 230,965 correct, 6,620
discordant, 19,027 unphased. N50 remains **955,496 bp**, largest block
**3,012,193 bp**, across 250 variant blocks. All **141 native/full panel
measurements and HiPhase classifications** match step 6. Among 118 HiPhase
>=80% windows, 59 still pass and 59 still fail the full closure contract.
The original primary denominator includes abstentions; total correct and
dominant connected-core correct must meet HiPhase, with rescue PS excluded
from core. No expectation, read floor, closure certificate or panel changed.
No new gap closure or biological accuracy gain is claimed for this diagnostic
fit stage. See `full-checks.json`, `panel-comparison.json` and `manifest.json`.
