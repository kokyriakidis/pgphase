# Representation step 14: molecule-supported complete-context hypotheses

The missing long sequence in the 65 Mb physical context is now represented as
a diagnostic hypothesis. Six distinct original primary molecules share its
exact 312-base sequence. Whole-cohort abPOA independently returns that sequence
and a 51-base short sequence, including after reversing input order and after
removing each of the six exact long supporters. This step adds hypotheses;
it does not yet change genotype fitting, read confidence or phase joins.

## Implementation

`build_read_allele_catalog` groups case-normalized complete-context read slices
by physical molecule name. At least two distinct names must share the exact
sequence. Duplicate source descriptions count once. A conflicting or invalid
name abstains entirely. Only nonempty A/C/G/T sequences up to 4,096 bases are
supported. Excluding a name removes all its records before validation and
discovery. The 64-hypothesis budget returns `limited` with no partial hypotheses
when exceeded. Exact support admits a hypothesis; it is not biological confidence.

Matrix diagnostics now preserve five additional tables: `read-catalog-status`,
`read-hypotheses`, `read-support`, `read-exclusions`, and `read-alleles`. The last
retains the entire previous nested catalog. No original candidate, source
observation, raw read slice, quality, physical context, original cost or fit is
replaced. Production discovery uses no MSA result, source preference, category,
genotype, parental truth, HP or PS. The optional MSA verifier is evaluation-only.

`test_allele_context --read-catalogs INPUT_FOLDER OUTPUT_PREFIX` reconstructs all
15 diagnostic tables from six raw tables, without BAM access or precomputed
outputs. The four owners have 70 contexts and 3,971 physical read/context rows
representing 1,683 distinct read names. They supply 255 exact-supported
hypotheses, including 63 sequences absent from the old catalog. All catalogs
complete; every old diagnostic table remains byte-identical.

## Residual measurement

The old path minimum is reused only after verifying its frozen hash and exact
catalog/query identity. Step 13 added no sequence to that catalog. Each scored
molecule is excluded from **new hypothesis discovery**, while the old catalog
remains fixed. This is not cross-validation of upstream discovery or genotype
fitting. Each novel sequence cost is independently checked by scalar dynamic
programming (2,917 checks). The reported minimum ranges over all hypotheses;
it is a representation lower bound, not a diploid fit or read-correctness rate.

| Owner | New unique sequences | Rows | Old residual | Full-cohort residual | Discovery-held-out residual | Improved held-out rows |
|---|---:|---:|---:|---:|---:|---:|
| 4–5 Mb | 25 | 1,796 | 712 | 398 | 436 | 86 |
| 5–6 Mb | 3 | 413 | 38 | 20 | 22 | 16 |
| 65–66 Mb | 21 | 425 | 6,455 | 1,258 | 1,360 | 176 |
| Whole-snarl 4–5 Mb | 14 | 1,337 | 323 | 203 | 213 | 64 |

Physical context 3 (`65,769,387–65,769,441`) contains 39 original molecules. Its
previous hypotheses were at most 55 bases long, while 21 query slices exceeded
55 bases. The four repeated exact sequences have lengths/support 51/12, 53/2,
55/2 and 312/6. Three are novel. Its edit residual falls from 4,383 to 711 when
all reads participate in discovery, and to 715 when each scored name is excluded
before discovery. The latter improves 28 of 39 rows. The six-read long group
still has five independent exact supporters in every excluded-name replay.

The independent whole-cohort MSA partitions the 39 molecules into a 16-read
long cluster and 23-read short cluster, with consensuses exactly matching the
312- and 51-base hypotheses. Both memberships and consensuses are identical
under reversed input order. Each long-supporter exclusion yields the same two
sequences and removes only that molecule from its cluster. The verifier also
checks that every input MSA row ungaps to the original query slice. The old
catalog lacks both MSA consensuses. This establishes a missing-hypothesis defect;
it does not establish that the GBZ has no walk carrying the sequence.

Post-hoc parental truth checks classify all six exact long supporters as maternal
and all twelve exact short supporters as paternal. The broader MSA clusters,
however, contain 11 maternal/5 paternal long reads and 4 maternal/19 paternal
short reads. Giving each cluster its majority parental label would yield only
30/39 correct (76.9%) in this context, below the requested 80% threshold.
Consensus-sequence agreement therefore cannot certify its cluster memberships
as phase calls. Truth is used only for this audit, never for discovery.

## Verification

Focused fixtures: **14,933 checks, zero failures**, including aliases, case,
invalid/conflicting names, singleton rejection, exact long sequences, the
length/work bounds, exclusion before discovery and 200 input permutations.
All eight deliberate production-rule mutations are rejected. Fixtures average
166 ms over 100 runs. Raw diagnostic replay averages 47 / 10 / 428 / 19 ms for
the four owners over 100 runs; the largest owner still enumerates the previous
12,509-path catalog. These measurements use the same unchanged test binary.

The four owning pipeline replays preserve candidate rows, VCF genotypes/PS,
primary HP/PS, parental status and original channel evidence. All unit tests,
the build, and HiFi/ONT goldens pass; HiFi t1/t4 output remains deterministic.
The full chromosome replay retains 82,287 candidate rows, 64,483 VCF rows and all
256,612 primary HP/PS assignments and parental statuses, including connected
core and previously phased reads. Correct/discordant/unphased remain
230,965/6,620/19,027. Its 250 blocks retain N50 955,496 bp and largest block
3,012,193 bp. Regression floors and the committed panel are unchanged.

All **107/107** registered checks pass in 720.6 seconds of final cold acceptance.
The official suite passes **18,138 assertions in five cases**, plus all four
cache tests. All **141** native/full panel metrics and HiPhase classifications
match step 13. Of 118 HiPhase >=80%-correct windows, 59 still pass and 59 still
fail the complete closure contract. Prior competitor measurements are reused
only after checking input/index/truth/panel/helper fingerprints. No floor,
expectation or closure certificate is refreshed. Source and binary remain
unchanged throughout final acceptance.

## Next step

Fit diploid models using the expanded hypotheses and rebuild discovery without
each scored molecule. Test residuals, pair stability and competing hypotheses
before calibrating confidence or allowing phase joins. Exact two-read support
can recur through systematic repeat/alignment errors and does not replace that
validation. Any eventual gap closure still needs at least 80% correct original
truth-scorable overlapping reads, with abstentions in the denominator, and
correct-read/connected-core performance on par with or better than HiPhase.

## Evidence and reproducibility

Fast development checks use saved state:

```bash
make gap-dev-check
./test_allele_context --read-catalogs \
  test_data/tmp_representation_step14/accepted-final/matrix65-final \
  /tmp/read-catalog
python3 evaluations/2026-10-07-representation-step14/audit_read_catalogs.py
python3 evaluations/2026-10-07-representation-step14/audit_residuals.py
python3 evaluations/2026-10-07-representation-step14/audit_msa.py
```

Owning pipeline replay uses step 5's `replay_matrix.sh` and
`replay_whole_snarl.sh` with the explicit frozen binary. Final acceptance uses
step 1's `check_windows.py` with `--binary` and a fresh `--out`, this directory's
`replay_full.sh OUTPUT_DIRECTORY BINARY`, and `make window-tests` with the same
`PGPHASE_BIN` and replay cache. Bench-phasers Python runs `audit_owners.py`,
`audit_full.py` and `audit_panel.py`; `compare_panel.py` checks every panel metric
against step 13. `make_manifest.py` freezes acceptance and
`verify_manifest.py` independently hashes it.

- [Independent catalog audit](read-catalog-checks.json)
- [Residual checks](residual-checks.json) and [scorer](score_read_catalog.cpp)
- [Independent MSA checks](msa-checks.json) and [verifier](msa_verify.cpp)
- [Raw-state replay](raw-state-checks.json)
- [Mutation checks](mutation-checks.json) and [timings](fast-checks.json)
- [Owning output comparison](owner-checks.json)
- [Full chromosome comparison](full-checks.json)
- [Complete panel comparison](panel-comparison.json)

The candidate binary SHA256 is
`5a17f417187a0a29a893781a9dc40b1c4072470c02d314a7a07fbc4ed57827bd`.
The unchanged prior accepted binary is
`5459c5d8b7a755403919f9361f4e8b119dbe9c162cf1d11c83db61d74d767484`.
`manifest.json` freezes the final code, inputs, baseline and all accepted evidence.
