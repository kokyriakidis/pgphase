# Evaluation Records

Each dated directory is a durable experiment record. Keep enough information
in the repository to reproduce a result without relying on session scratch
space.

Every evaluation directory should contain:

- `README.md`: hypothesis, input inventory, result table, interpretation, and
  follow-up decision;
- `commands.sh`: executable commands with every non-default phasing threshold;
- `results.tsv`: one stable, machine-readable row per evaluated configuration;
- per-configuration `summary.json` and `summary.txt` from
  `scripts/evaluate_phase_accuracy.py`;
- `worst_phase_sets.tsv` and `bad_phase_regions.bed` when read truth is used.

The maintained chr12/18/20 panel is controlled by
`scripts/benchmark_panel.py`. Competitor results are immutable artifacts in a
checked lockfile; normal development runs execute pgphase only. Use
`make benchmark-tests` for framework tests and `make benchmark-report` to
rebuild presentation tables from the frozen baseline.

Large BAM, GAF, VCF, matrix, and per-read files stay outside git. Their paths
and derivation commands must be recorded in the experiment README/script.

Current evaluations:

- [`2026-09-12-chr12-18-20-comparison`](2026-09-12-chr12-18-20-comparison/README.md):
  HiPhase-style shared-call comparison across chr12, chr18, and chr20, with
  variant truth, read truth, chromosome NGC50, and competitor-gap diagnosis.
- [`2026-09-12-native-bam-graph-lock`](2026-09-12-native-bam-graph-lock/README.md):
  chr12/chr18 native-BAM and chr20 DeepVariant-private validation of
  graph-locked private-gap recovery.
