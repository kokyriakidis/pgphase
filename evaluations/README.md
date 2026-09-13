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

Large BAM, GAF, VCF, matrix, and per-read files stay outside git. Their paths
and derivation commands must be recorded in the experiment README/script.

Current evaluations:

- [`2026-09-12-native-bam-graph-lock`](2026-09-12-native-bam-graph-lock/README.md):
  chr12/chr18 validation of graph-locked native-BAM private-gap recovery.
