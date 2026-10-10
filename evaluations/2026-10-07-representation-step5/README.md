# Representation repair: step 5, explicit allele pairs and molecule evidence

The joint view now retains **both selected sequence alleles**, including ALT/ALT,
with a canonical-to-original allele-index mapping. Each full allele is validated
against reference and normalized independently. Literal REF has an explicit
identity; an unknown call has no identity. Reversing the ALT list or source
haplotype order changes the mapping, not the unordered pair. Two REF/ALT rows
sharing one ALT cannot masquerade as an ALT/ALT genotype.

`candidate_allele_contrasts` preserves each candidate's pair and original indices;
`joint_allele_contrasts` groups equal pairs. Binary candidates retain their existing
pair. Multiallelic candidates require an already selected distinct consensus pair;
this view does not infer a genotype from depth. Metadata, original classes, counts,
source phase gauges, recovery indices, physical contexts and quality certificates
remain on their original rows. The old REF/ALT locus view keeps its original 0/1
contracts.

Graph decomposition now records `GraphSiteMeta::non_selected_alt_class` when zero
means another non-selected ALT. Such a collapsed class cannot certify a literal
sequence contrast. `joint_contrast_molecule_evidence` maps original graph/BAM
observations to a complete pair's canonical gauge. Dependent descriptions provide
one call. Opposing calls and pre-existing BAM conflicts remain conflicts; calls
outside the selected pair abstain, including genomic REF at an ALT/ALT site.
Original graph/BAM channels are inspected together so primary-channel preference
cannot erase a disagreement. Query indices survive only when agreeing and are
never added as quality. Source observations remain traceable to original candidates.

With `--phase-matrix-dump`, `.chunkN.allele-contrasts.tsv` records pair identities,
original indices, classification, source PS, MSA provenance and collapsed-class
flags. `.chunkN.joint-molecules.tsv` records one canonical reduction per
molecule/duplicate contrast, including all retained observations. Raw sequence
contexts, original channels and physical qualities remain in the full matrix and
candidate metadata. Dumps are diagnostic evidence, not additional phasing votes.

This stage supplies the explicit-pair evidence adapter. It does not widen clean
admission, certify a primitive REF call as the full REF of a parent snarl, select
a winning channel by quality, or replace the recovery/component solver. Those
consumers still need calibrated full-context observation support. The step-4
coverage and independent-parent guards remain intact. Parental truth is used
only in evaluations.

## Verification

Starting executable SHA256:
`7bccceb86f200cb11ae7526d215d368c9ddc47255945747101fad008a024758f`.
Production SHA256:
`cf849f5fa6fac264f2b5d7b01f6c768590da9cf8ac578e478264014ea53b1750`.
The baseline binary, source snapshots and initial file hashes are preserved under
`test_data/tmp_representation_step5/`.

The identity tests verify **5,664 independent complete haplotype-pair fixtures**,
including reference, ALT/ALT, motif shifts, padded edits and reversed allele order.
The production adapter checks remapped graph walk indices, reordered MSA ALT
lists, matching binary/multiallelic profile gauges, alias-only sparse coverage,
conflict retention, unknown/outside-pair abstention, invalid consensus indices,
ambiguous unselected genotypes and unchanged source metadata/quality channels.
The real graph builder's collapsed-class flag is checked as well.

Eight mutations of the actual production implementation fail: dropping the first
allele, ignoring order, inventing REF membership, preferring the primary channel,
accepting collapsed classes, erasing an existing conflict, omitting the graph
class flag and ignoring original walk indices. Failure counts are 9, 7, 1, 3, 1,
1, 1 and 2. Over 100 executions, the focused adapter averages **2.04 ms** and the
complete-pair oracle fixtures **8.07 ms**, including startup. See `verify_fast.py`,
`fast-checks.json` and `mutation-checks.json`.

Build, unit tests, all 47 phasing predicate cases (1,551 assertions), HiFi/ONT
TSV/VCF goldens and HiFi t1/t4 determinism pass. The edit normalizer still agrees
with bcftools 1.24 on all **64,483** frozen records, including 504 realignments.
HiPhase was measured afresh on **141 windows with identical primary alignments**;
the measurements match the existing committed competitor state.

All **107/107 registered window/mechanism checks** pass in 714.8 seconds cold.
The cached standard `make window-tests` passes **18,138 assertions in five cases**,
including all four replay-cache checks.
All **141 native and full panel measurements** and their HiPhase classifications
match step 4. Among 118 windows where HiPhase correctly phases at least 80% of
original primary truth-scorable reads, the full closure contract still has 59
passing and 59 failing windows. No expectation or threshold was changed.
See `window-checks.json`, `panel-audit.json` and `panel-comparison.json`.


## Real observation audit

| Owning region | Complete duplicate pairs | Molecule/locus records | Conflicts | Distinct conflicting reads |
| --- | ---: | ---: | ---: | ---: |
| 4–5 Mb | 31 | 1,749 | 104 | 101 |
| 5–6 Mb | 7 | 413 | 7 | 7 |
| 65–66 Mb | 11 | 391 | 0 | 0 |

All original saved matrices, candidate/VCF rows, primary HP/PS and parental
statuses match step 4 in all three owners. The audit independently verifies
canonical call membership, coverage of all retained graph/BAM calls, reduction,
query-index agreement and traceability to original channels.

In 4–5 Mb, **94** conflicting molecule/loci have graph repeat ALT versus BAM
noisy/MSA REF, **nine** the reverse, and **one** opposing graph descriptions.
At 5–6 Mb the seven conflicts split six/one in those two graph/BAM directions.
None arise from opposing retained graph/BAM calls at the same candidate; they
cross descriptions of the same normalized contrast. Identity therefore does
not resolve source confidence. These observations remain explicit conflicts,
with original physical contexts and classifications available for calibration.
See `audit_contrasts.py` and `contrast-checks.json`.

These default owning replays contain only REF/ALT candidates; they do not
establish recovered ALT/ALT genotypes. ALT/ALT mapping is checked by the
production-adapter fixtures and the separate existing whole-snarl mode replay.

The 4–5 Mb whole-snarl replay (`--snarl-allele-phasing --snarl-keep-whole`) has
1,614 candidates: **1,431 explicit pairs**, including **three real ALT/ALT pairs**.
It retains **99 collapsed other-ALT class flags** and 84 further rows without an
explicit selected/valid pair. There are 23 duplicate complete pairs and 1,370
molecule/locus records, with 51 conflicts across 50 reads. All retained calls
pass the same independent original-channel audit. Its full original matrix,
emitted rows and primary HP/PS/parental classifications match the starting binary
under those same options: 3,866 correct/27 discordant/235 unphased. This is a
representation check of the existing mode, not a new closure acceptance.


The accepted full chr20 replay preserves all **82,287 candidate rows**, **64,483
VCF records** and **256,612 primary HP/PS assignments and parental states**:
230,965 correct, 6,620 discordant, 19,027 unphased. N50 remains **955,496 bp**,
largest block **3,012,193 bp**, across 250 variant blocks. No new closure or
biological accuracy gain is claimed for this adapter stage. See `full-checks.json`.

## Reproduction

Run from the repository root. Output audits use the benchmark Python containing
pysam: `/home/kokyriakidis/micromamba/envs/bench-phasers/bin/python`.

```sh
make -j10
make unit-tests gap-dev-check
./test_allele_identity --contrasts
./test_graph_bam_adapter --joint-contrasts
python3 evaluations/2026-10-07-representation-step5/verify_fast.py
bash evaluations/2026-10-07-representation-step5/replay_matrix.sh \
  test_data/tmp_representation_step5/matrix4-final
bash evaluations/2026-10-07-representation-step5/replay_matrix.sh \
  test_data/tmp_representation_step5/matrix5-final ./pgphase CHM13#0#chr20:5000001-6000000
bash evaluations/2026-10-07-representation-step5/replay_matrix.sh \
  test_data/tmp_representation_step5/matrix65-final ./pgphase CHM13#0#chr20:65000001-66000000
bash evaluations/2026-10-07-representation-step5/replay_whole_snarl.sh \
  test_data/tmp_representation_step5/whole-before test_data/tmp_representation_step5/pgphase.before
bash evaluations/2026-10-07-representation-step5/replay_whole_snarl.sh \
  test_data/tmp_representation_step5/whole-final
make check
python3 evaluations/2026-10-06-representation-step1/check_windows.py \
  --binary ./pgphase --out test_data/tmp_representation_step5/accepted/windows
make window-tests
bash evaluations/2026-10-06-shared-insertion-source/replay.sh \
  test_data/tmp_representation_step5/accepted/full-current
```

After the corresponding complete replays, run `audit_owners.py`,
`audit_contrasts.py`, `audit_full.py`, `audit_panel.py` and `compare_panel.py`.
The gap contract retains the original primary-read denominator including
abstentions, total correct >= HiPhase and dominant connected-core correct >=
HiPhase. No expectation, read floor or closure manifest is relaxed.
