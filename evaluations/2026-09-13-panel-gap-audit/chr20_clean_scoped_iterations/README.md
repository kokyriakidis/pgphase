# Clean-round phase-set scoping: accepted local-panel fix

Completed all 114 local chr20 cases: **88 joined, 24 split, two endpoint-unphased**. All 86 confident_clean_bridge joins remain. The two repaired targets are overlapping windows at chr20 1.9 Mb: 1912570_1913054 and 1912570_1913912.

The SNP observations already separate the truth haplotypes cleanly. K-means pooled arbitrary HP labels across disconnected phase sets, contaminating read assignment and allele-profile updates. Clean-candidate iterations now score and update each read separately within each observed phase set. Output HP is recomputed within its reported PS using the final consensus. MSA-round iterative updates retain their existing path because the broader trials exposed unresolved repeat-link regressions.

In each repaired window, 354 reads are assessed: discordant reads decrease **30 to 3**, and read switch/flips **10 to 3**. WhatsHap 2.8 variant comparison assesses 203 instead of 171 variants in the first window (204 instead of 172 in the second), with blockwise hamming **5 to 1**. Across all 114 windows there are **no variant-hamming-count increases and no original-block majority reversals**.

The ambiguous chr20_26602087_26651086 endpoint remains unphased. Its discordant-read count increases 34 to 35 and read switch/flips 15 to 16, with unchanged variant-truth metrics. Existing recovery read errors at 7.2 and 61.7 Mb remain. This is not a claim of zero individual read errors or whole-chromosome NGC50/hamming validation.

Build and all unit tests pass. The new fixture checks internal allele-profile partitioning and output HP/PS consistency for reads crossing disconnected blocks; it fails against the original core. The frozen binary hash is in binary.sha256. The final source binary matches it.

Reproduce the panel and comparisons from the repository root:

```bash
python3 evaluations/2026-09-13-panel-gap-audit/run_link_panel.py --chromosomes chr20 --label chr20_clean_scoped_iterations --all-concordant --workers 8
python3 evaluations/2026-09-13-panel-gap-audit/summarize_link_panel.py --label chr20_clean_scoped_iterations --baseline-label chr20_confident_clean_bridge
python3 evaluations/2026-09-13-panel-gap-audit/compare_panel_variant_truth.py --label chr20_clean_scoped_iterations --baseline-label chr20_confident_clean_bridge --whatshap /home/kokyriakidis/micromamba/envs/bench-phasers/bin/whatshap --workers 8
```

A label freezes its binary on first use. Use a new label to evaluate changed code. The case manifest and command files record the local regions and inputs; truth and competitor records are used for evaluation only. The remaining inventory is ../chr20.remaining.clean_scoped_iterations.tsv. Rejected MSA trials and their source patches are documented in CHECKPOINT.md and neighboring experiment reports.
