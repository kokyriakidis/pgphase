# Representation repair: step 3, source molecule reduction

The whole-chunk BAM overlay now reduces equivalent source descriptions once per
read and unique graph destination. Previously, normalized source collisions
were discarded and exact source matches were transferred one at a time, so the
first eligible call could win. The new `MoleculeAlleleEvidence` reducer accepts
agreeing source calls as one observation and records opposing calls as
`kConflictingBamAllele`, which supplies no BAM vote. Repeating a call is
idempotent. Original source indices, alleles and query coordinates remain in the
reducer during transfer; the source chunk is not rewritten. Differing query
coordinates become zero (unknown), rather than being added or treated as
confidence scores.

The destination still must satisfy the step-1/2 identity rules. Multiple graph
rows claiming the same identity remain ambiguous; source aggregation cannot
select one of them. Every newly normalized source call still requires coverage
of both physical edits and their outside flanks, plus an independent graph
observation. Existing targeted-recovery BAM calls, conflicts and measured SNP
qualities retain precedence. No source HP/PS is imported into a graph gauge.

This is the source-evidence part of molecule merging. Candidate rows retain
their original graph/MSA representation, genotype/classification and recovery
path state. In particular, marking a graph-key row `bam_injected` would make its
topology identifier look like physical ALT to the writer; discarding an injected
row also discards recovery predicates that require that provenance. A later
joint-locus representation must preserve those contracts before compacting rows.
The existing 19 normalized duplicate output pairs and split sample ALT contrasts
are not claimed as repaired here. No new gap closure is claimed.

## Verification

The frozen step-2 executable is
`5c47ca57f979344cce86d167b18225b20d9ffc7accc47a99dd602c7236694ee7`.
Its source snapshot and executable are retained in
`test_data/tmp_representation_step3/baseline_src/` and `pgphase.before`.
The step-3 production executable is
`fe0ae0862c1f3596921218fee6b3f4b5f6cdcef9d6ee674df1cd687cc3d69b50`.

Fast fixtures cover 42 order permutations, repeated descriptions, permanent
conflicts (including opposing calls from the same source), unknown calls,
query-coordinate disagreement, distinct equal-position ALT sequences, and
different reads. A normalized repeat-shifted cohort exercises matching and
reduction together and verifies that destination ambiguity still abstains.
The 3,328 independent complete-haplotype normalization fixtures remain enabled.
One hundred runs average 3.7 ms including process startup. Mutants that ignore
conflicts, borrow query coordinates, or retain repeated provenance fail 25, 21,
and 48 checks respectively (`verify_fast.py`, `mutation-checks.json`).

The 5–6 Mb owning matrix is byte-identical to step 2: 1,825 candidate states,
4,256 read states, independent graph calls and measured SNP qualities are
unchanged. This chunk has no duplicate source cohort to reduce. Its emitted
candidate/VCF rows and all primary HP/PS and parental statuses are unchanged:
4,001 correct, 23 discordant and 232 unphased. See `audit_matrix.py`,
`audit_owner.py`, and their JSON results. The two owning regressions that rejected
broader normalized admission in step 2 pass their unchanged checks: the 33 Mb
insertion source (1,007 assertions) and 47 Mb independent SNP pairs (333).

Build, unit tests, the 47 phasing predicate cases (1,551 assertions), HiFi/ONT
TSV/VCF goldens, HiFi t1/t4 determinism, and all four replay-cache tests pass.
`bcftools norm` agrees with the shared identity helper for all 64,483 frozen
chromosome records, including 504 realigned records. The validated HiPhase panel
state reuses the unchanged original alignment, truth, panel and competitor
identities; its measurements are unchanged. No expectation or floor was relaxed.

The full chromosome replay has the same 82,287 candidate rows and 64,483 VCF
records. Candidate TSV and phased BAM are byte-identical; parsed VCF records
are identical. All 256,612 primary HP/PS assignments and parental statuses are
unchanged: 230,965 correct, 6,620 discordant and 19,027 unphased. There are still
250 variant blocks, N50 955,496 bp, and largest block 3,012,193 bp. See
`full-comparison.json`. These results verify preservation; they do not establish
an improvement in biological phasing from this source-reduction stage.

All 107 named window/mechanism checks pass (707.3 seconds cold). All 141 native
and full-chromosome panel measurements are identical to step 2, including the
80% denominator, total-correct and dominant-core HiPhase comparisons. Among the
118 windows where HiPhase reaches 80%, 59 still pass the full closure contract
and 59 still fail it. These historical deficits remain visible rather than
being reclassified as closures. See `window-checks.json`, `panel-audit.json` and
`panel-comparison.json`.

The standard cached `make window-tests` also passes all 18,138 assertions in
five test cases. The final build retains the evaluated executable checksum.

## Reproduction

Run from the repository root. Use
`/home/kokyriakidis/micromamba/envs/bench-phasers/bin/python` for audit scripts
that import `pysam`.

```sh
make -j10
make unit-tests gap-dev-check
python3 evaluations/2026-10-07-representation-step3/verify_fast.py
make gap-owner-check GAP="a recovered short insertion retains"
make gap-owner-check GAP="independent BAM SNP pairs certify"
make check
python3 evaluations/2026-10-06-representation-step1/check_windows.py \
  --binary ./pgphase --out test_data/tmp_representation_step3/windows
make window-tests
bash evaluations/2026-10-06-shared-insertion-source/replay.sh \
  test_data/tmp_representation_step3/full-current
```

The matrix replay uses the same graph arguments as `replay.sh`, region
`CHM13#0#chr20:5000001-6000000`, one thread, output directory
`test_data/tmp_representation_step3/matrix5`, and
`--phase-matrix-dump test_data/tmp_representation_step3/matrix5/matrix`.
Run `audit_matrix.py`, `audit_owner.py`, `audit_full.py` and `audit_panel.py`
after their corresponding replays. The full panel audit scores the original
overlapping primary alignments, including abstentions, against parental truth
and HiPhase, with rescue phase sets excluded from the connected-core count.
