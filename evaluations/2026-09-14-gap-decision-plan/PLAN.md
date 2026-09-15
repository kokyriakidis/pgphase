# Experiment plan: joint evidence for gap orientation

## Objective and limits

Increase correctly connected phase blocks while controlling whole-component
orientation errors. Robust handling includes an explicit unresolved outcome:
identical available reads do not guarantee that a trustworthy orientation can
be inferred by this implementation. No promise to close every gap.

Use only reference, graph/GAF and surjected BAM/MSA observations for inference.
Competitor calls select diagnostic cases; parental truth is evaluation-only.
DeepVariant calls must not enter discovery, weighting or stitching.

Freeze the accepted `/tmp/pgphase-hp-confirm-only` build/output as the baseline:
123/280 recovered gaps, 2,427 discordant / 189,042 evaluated reads,
351 switches and 414 flips. Record binary and source hashes: HEAD alone does
not identify this uncommitted implementation.

## What the current implementation loses

- Tier processing exits on the first join, so later evidence is unexamined.
- `GapStitchResult` retains booleans, not the two orientation scores or the
  difference between missing support and positive contradiction.
- Graph/BAM/MSA observations from one molecule are correlated. A second solve
  from these reads is a sensitivity check, not independent replication.
- Sparse opposite observations can veto a bridge without regard to reliability.
- Frozen original read assignments already exist, but proposal construction
  reads mutable chunks across recovery waves. A mutex prevents concurrent
  access; it does not prove order independence or eliminate evidence feedback.
- Small-region rephasing changes the anchor blocks. It cannot by itself validate
  a chromosome-wide connection, as the 23.48 Mb cascade demonstrated.

## Phase 0 — Reproducible, immutable evidence replay

First retain the current decision rule and export chromosome-context inputs
before any gap mutates reads or candidates. Keep the existing chromosome MSA
cache; add a compact per-gap view, rather than another copy of all BAMs.

Each export identifies input/cache schema hashes, options, original block IDs
and orientations, complete original block bounds, local sequence, candidate
allele sequences, original read memberships, and per-read graph/BAM/MSA
observations. Preserve missing and contradictory observations explicitly;
missing is never reference. Preserve the BAM query index and source provenance
separately from graph confirmation. Deduplicate by input + molecule identity
across chunks and alignments. Keep untagged molecules as evidence.

Export original 2x2 flank-vote matrices and candidate proposals for every tier
and evidence view; do not short-circuit when one joins. These are alternative
hypotheses, not independent votes. Initially generate proposals with the
existing k-means. Score without writing BAM or rerunning MSA/k-means. If all
proposals fail to express a connection supported by the observation graph,
label that a proposal-generation failure rather than weakening the scorer.

Verification: old decisions reproduced from frozen chromosome evidence;
reversing job order and changing worker count cannot change decisions;
scoring changes cannot change export hashes; truth never appears in exports.
The export must retain the 23.42→23.48 Mb relationship, not crop it away.

## Phase 1 — Comparable allele observations

Group overlapping SNP/indel descriptions into local sequence events. Confirm
sequence equivalence using the cached reference and observed read sequence;
coordinate proximity alone is not equivalence. Preserve distinct alternatives
including GT 1|2. An ambiguous repeat-length observation stays ambiguous.

For each molecule/event retain allele evidence, base and mapping quality,
sequence agreement, alignment/repeat context, source agreement/disagreement,
and usable flank anchors. Do not sum a SNP and overlapping deletion as two
independent observations, or treat MSA discovery and its contributing reads as
independent confirmations. Do not promote a HOM call from allele balance alone.

Verification fixtures: shifted-equivalent indels; distinct repeat alleles;
correlated boundary SNPs; unknown qualities; duplicate observations; missing
reference support; multi-allelic insertions/deletions; untagged bridging reads.
Every transformation must be auditable back to the original observations.

## Phase 2 — Compare SAME, FLIP and unresolved

Assess both relative orientations against exactly the same evidence. Existing
block labels define a fixed gauge, with read membership/allele uncertainty
retained; numeric HP 1/2 is not parental truth. Reads are latent diploid
assignments. A supported chain through gap sites can connect the blocks even
if no individual read spans both ends.

Start with the simplest molecule-level score that can express the known cases.
For each event, compare support for SAME and FLIP while allowing uncertain or
erroneous observations. Count a molecule once per evidence contribution and
account for correlated local events. Qualities inform reliability but must not
be treated as calibrated repeat-alignment error probabilities. An uncalibrated
score is not a posterior probability. Expose support, opposition, missingness,
and the decisive molecules/events rather than only a total score.

Keep outcomes SAME, FLIP, INSUFFICIENT, CONFLICTING and UNSUPPORTED_REPRESENTATION.
No evidence under one source cannot count as evidence for the opposite parity.
A clean site can outweigh noisy sites only through an explicit, validated
reliability model, not by a region-specific exception.

Evaluate influence by leaving out a molecule or correlated event group.
Distinguish losing all support from reversing orientation. A single-event
bridge is not automatically rejected, but must be evaluated/calibrated as a
separate sparse-evidence class. Grouped read holdouts test circularity where
coverage permits; full MSA leave-out rebuilding is a diagnostic for decisive
borderline cases, not part of every scoring replay.

## Phase 3 — Controlled ablations before threshold tuning

Use the same immutable input and apply these changes one at a time:

| Arm | Question |
| --- | --- |
| A: existing rule | Does frozen replay reproduce the reference? |
| B: all tier proposals, three-way evidence status | Does missing-versus-conflicting explain the withheld correct links? |
| C: molecule/event deduplication and representation handling | Were counts inflated or useful alleles mismatched? |
| D: reliability-weighted SAME/FLIP score | Can strong evidence survive weak disagreement without admitting wrong joins? |
| E: influence checks | Can we identify orientation decisions dominated by an unstable molecule/event? |
| F: final component consistency | Do accepted local relationships compose without new cascades? |

Reject added complexity that does not improve held-out correctness/coverage.
Choose acceptance thresholds on development data only. Plot accepted-edge
error versus accepted coverage; stratify sparse bridges, repeat contexts,
source disagreement and block sizes. Report sample sizes and uncertainty;
zero errors on a handful of sites is not a validated error probability.

## Panels and leakage controls

Use the already studied regions as development/regression cases:

- Correct rescues: 0.86, 19.37, 23.48, 37.98 and 53.9 Mb.
- Incorrect upstream connections: 23.42 and 37.64 Mb.
- Withheld previously correct connections: 57.84 and 62.41 Mb.
- Remaining splits: 1.08, 12.72, 17.61, 47.67 and 56.00 Mb.
- Explicit reporting case: insertion endpoint 23,480,815 (GT 1|2).

Add ordinary already-correct boundaries and ambiguous/no-bridge examples; a
panel containing only competitor successes cannot estimate false-join risk.
Artificial splits inside well-supported blocks are useful controls, but do not
substitute for real gaps. Label old blocks from independent parental evidence;
mark internally mixed or weakly labelled blocks as uncertain.

Group train/evaluation splits by shared molecules, connected blocks and
nearby windows, not by target row. All previously examined chr20/chr12/chr18
regions are development data. Reserve documented untouched blocks for locked
validation before scoring them; after inspecting an error there, that block
becomes development data and a new validation set is required. Without a
sufficient untouched set, report regression performance only.

## Phase 4 — Component-safe application

Treat original phase sets as nodes and accepted SAME/FLIP relationships as
edges. Check contradictory parity cycles using all available nonadjacent
links. Do not silently discard contradictory cycle evidence because a greedy
union already connected its endpoints. Cycle consistency cannot certify a
wrong bridge in a tree; local confidence and influence checks remain required.

Apply accepted relationships once to reads and candidates. Attach gap-only
sites after orientation resolution. Initially do not retag existing reads;
adding new tags is a separate measured change. No worker may alter evidence
used to decide another gap. Unit tests must verify global HP-label invariance,
edge ordering, thread-count invariance, duplicate inputs, inconsistent cycles,
and failed proposals leaving both flanks unchanged.

## Phase 5 — Promotion gates and efficient execution

Run tests in increasing cost order:
1. Synthetic representation/parity/invariance regressions.
2. In-memory replay of development panels, with no BAM output or realignment.
3. Replay all 280 chr20 gaps from the chromosome snapshot.
4. One full output/evaluation for the best unchanged configuration.
5. Locked validation on untouched blocks/chromosomes before default adoption.

Freeze original read cohorts and evaluate each new edge's parity plus final
component orientation. Measure previously concordant→discordant transitions,
reads newly tagged or lost, switches, and the number/span of reads or sites
affected by an erroneous edge. Per-block majority orientation can hide effects
of splits, so report split/merge changes and fixed-cohort transitions together.
The known 23.42 and 37.64 Mb erroneous edges must stay rejected; recovering
57.84/62.41 Mb is an explicit development objective, not a threshold override.

Retain the current regression gate of no new previously-concordant→discordant
reads on the fixed audited cohort. Report uncertain truth separately without
silently excluding it. On locked validation report false-join rate with
uncertainty and coverage, not an all-cases-correct claim. Measure variant-level
switch/Hamming and NGC50 on an independently validated comparable site set;
if unavailable, explicitly state they were not measured. A larger read-based
N50 alone does not pass. Competitor/DV variants, if used to define an evaluation
site set, remain outside phasing and calibration inputs.

Profile export time, replay time, peak memory and final-output time separately.
Aim for a panel scoring sweep in seconds; this is an engineering target, not
a measured result. Load evidence once per process, reuse the sparse observation
structure, run independent gaps with a bounded worker pool, and keep worker
results immutable. Never launch jobs×threads beyond the intended CPU budget.
Only rebuild MSA for a demonstrated missing observation/representation, once,
with explicit cache versioning. Threshold changes must not invalidate raw
observation caches. Preserve manifests, hashes, scores, reasons, and metrics
for every arm; do not repeatedly rerun chromosome BAM output for score tuning.

## Immediate implementation order

1. Frozen chromosome-context exporter and read-only replay reproducing A.
2. SAME/FLIP score diagnostics and explicit insufficient/conflicting states.
3. Representation/provenance fixes and the B–E ablation loop.
4. Component checks and one final application pass.
5. Full chr20 plus locked validation; only then replace default acceptance.

This plan does not enable a new production scorer. The initial experiment below
checks reproducibility of the existing regional runner; it is not evidence that
the proposed score or frozen chromosome exporter already exists or is correct.

## Initial experiment completed

Replayed all 11 target rows with the current binary and warm evidence caches:
`robust_plan_replay_v1` used 2 jobs × 2 threads; `robust_plan_serial_v1` used
1 job × 1 thread. Both ran the existing baseline and graph-BAM arms, giving
22 paired comparisons (44 regional executions total).

All 22 comparisons have identical read identities/alignment locations/flags,
HP/PS tags, native VCF data records and endpoint statuses. Each arm reproduces
six joined target rows and five splits. Median per-arm execution was 0.802 s
in the parallel run and 0.918 s in the serial run. These overlapping single
runs do not establish a speedup. The results show the existing regional loop
is usable for fast preliminary checks; they do not validate a new score,
prove chromosome-wide order independence, or add new truth measurements.

Artifacts: `replay_comparison.json`, both manifests/result tables,
`source_manifest.json`, and `compare_replays.py`. The latter compares semantic
read assignments and VCF data rather than compression/header byte differences.
Reproduce the comparison with:

```sh
python3 evaluations/2026-09-14-gap-decision-plan/compare_replays.py \
  --left /tmp/pgphase-gap-trials/runs/robust_plan_replay_v1 \
  --right /tmp/pgphase-gap-trials/runs/robust_plan_serial_v1 \
  --output /tmp/gap-replay-comparison.json
```
