# One recovery path

The hybrid carried three ways to recover a window the graph's sites could not
phase. Two are removed; the restructured one is kept.

## What is kept

`--retry-unphased-with-bam`. `collect_unphased_windows` finds a window where the
solve left reads unphased, and the window is re-solved by asking for the
noisy-region MSA **by name** (`force_noisy_msa`) with the noisy k-means enabled
inside `retry_windows`. Nothing else in the pipeline is disturbed: asking by name
rather than clearing `recover_gaps` is what keeps the other guards in place.

## What was removed, and why each removal is safe

Every gate removed here took the branch the **default already took**, because
each was off by default.

| removed | lines | note |
|---|---:|---|
| `gap_recovery.cpp` / `.hpp`, `gap_evidence.cpp` / `.hpp` | 1,507 | the `--recover-gaps` pass |
| `test_phase_block_stitch.cpp` | 2,372 | tested only that pass |
| pipeline driver: `recover_one_hybrid_gap`, `recover_hybrid_gaps`, `write_gap_audit_input`, `populate_gap_msa_cache`, `gap_recovery_jobs_conflict`, the whole `gap_cache_*` family, the recovery rounds block | 1,316 | |
| the gap-link machinery in `collect_phase.cpp` | 408 | see below |
| CLI surface and `Options` fields for 13 flags, plus three validation blocks whose only job was rejecting combinations of them | 117 | |

**The gap-link machinery could not survive the flag.** `recovery_graph` was
`recover_gaps && link_by_alleles && private_msa_admit_all_in_region`, so with
`recover_gaps` gone every branch it guarded is unreachable: the multiallelic
seeding branch, `gap_hp_link`, `unsupported_gap_indel`, the allele-vote
accumulation, the MSA bridge margin check, the edge list with its union-find
pass, and `select_gap_link_sites`.

**The readmission loop inside the retry went too.** It restored categories that
recovery had zeroed, and -- as the comment beside it recorded -- could only ever
fire when `recover_gaps` was also set. The re-solve below it is the retry's
actual mechanism and is untouched.

**Two gates become unconditional**, because their false branch was the default:
`collect_noisy_vars_step4` in `collect_var_run_phasing`, which ran whenever
`recover_gaps` was off, and `prune_not_candidate_variants` at the end of
`process_chunk_hybrid`, which was skipped only because recovery consumed the
candidate-indexed profiles. `collect_phase_noisy`'s gate narrows from
`recover_gaps || force_noisy_msa` to `force_noisy_msa`.

## Verification

| check | result |
|---|---|
| **chromosome-wide, stock defaults** | **212,320 tagged / 452 blocks / 0.559% -- exactly the standing baseline** |
| panel, stock defaults | records identical to the pre-removal build: 0 lost, 0 gained, 0 concordant-to-discordant |
| panel, `--retry-unphased-with-bam` | identical: 0 lost, 0 gained; 4 of 6 spanned at 99.49% with 0 concordant-to-discordant |
| injection verifier | 1,575 candidates: dropped 0, duplicated 0, missing 0, attributes 0, same single verdict item |
| `-Wall -Wextra` on every touched unit | no unused function or variable |
| unit suite | 4/4 (the fifth binary was the recovery-stitch test) |

Net: **5,787 lines removed, 35 added.**

## What this leaves to iterate on

Two arms, and nothing else that recovers a window:

- **default** -- 212,320 tagged / 452 blocks / 0.559%, panel 0 of 6 spanned at 99.68%.
- **`--retry-unphased-with-bam`** -- panel 4 of 6 spanned at 99.49%, 0 concordant-to-discordant.

The historical evaluation directories still name the removed flags. They are
records of what was run at the time and are deliberately left as they are; the
measurements that justified each removal live beside them, in
`evaluations/2026-09-17-gap-bam-only/` and `evaluations/2026-09-17-two-stage-refine/`.
