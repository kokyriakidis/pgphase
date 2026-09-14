# Combined scoped-k-means and MSA anchor diagnostic

Not accepted. This combines the MSA experiment archived in remaining_msa_site_anchors_experiment.patch with per-phase-set iterative updates and output scoring. Full local panel: 87 joined, 25 split, two endpoint-unphased. It fixes both 1.9 Mb endpoint targets and adds three MSA joins, but loses four existing joins. No original-block majority reversals; the ambiguous 26.6 Mb case adds one read error and one switch/flip. See baseline_comparison.tsv and summary.json.

Subsequent experiments separate the k-means bug from MSA allele/eligibility changes. Competitor and read-truth data are evaluation inputs only.
