# Gap-local homopolymer recovery experiment

The ordinary clean → clean + MSA SNP → clean + MSA SNP + MSA indel sequence
runs first. Only an unresolved gap can reach tier 4. This tier reuses the same
observations and k-means, enabling MSA-verified homopolymer link nodes only
between the current flanks. It does not collect another set of sites.

A repeat edge must have net allele support at least `min_block_link_reads`.
The proposal must connect both original flanks under the normal stitching
rule; on each flank, each original haplotype must independently favor the
chosen orientation by that same margin. Otherwise the trial is rejected
without changing the original blocks or adding reads. Earlier successful
tiers never reach this trial. These gates use input reads, not competitor
phase or truth labels.

Reproduce the local panel:

```sh
python3 evaluations/2026-09-13-panel-gap-audit/run_link_panel.py --chromosomes chr20 --label chr20_adaptive_hp --all-concordant --workers 4
python3 evaluations/2026-09-13-panel-gap-audit/summarize_link_panel.py --label chr20_adaptive_hp --baseline-label chr20_graph_support2
```

The label freezes its binary in `/tmp/pgphase-chr20_adaptive_hp/pgphase`;
use a fresh label after changing source. `binary.sha256` records provenance.
`manifest.json` lists the same 114 truth-concordant competitor endpoint pairs
used in previous local trials. Some windows overlap and represent the same
physical break; sums of their read counts are not chromosome-wide metrics.

Initial checks: chr20 62.4 Mb joins with 385 assessed reads and zero discordance.
At 24.1 Mb and 36.0 Mb ordinary indel recovery already succeeds and remains
unchanged, avoiding the 53 and 5 errors from unrestricted homopolymer admission.
The pre-existing incorrect ordinary-tier join at 19 Mb remains unresolved.

This is a local validation of evidence selection. The separate full-chromosome
`graph_support2` run excludes this fallback and reveals worse blockwise hamming
than nearest-link recovery despite local improvements. Its metrics are in
`../chr20.metrics.tsv`; no full-chromosome accuracy claim follows from this panel.


## Completed result: superseded acceptance rule

The first fallback trial retained all 85 previous target joins and added four.
However, chr20 55.9 Mb introduced 39 discordant reads (0→39 of 308), so its
net-support acceptance rule is insufficient and is superseded. Diagnostics
showed opposing votes on both homopolymer bridge edges (4:1 and 18:3), despite
perfectly consistent original-block anchor votes. Anchor agreement alone does
not validate the bridge orientation. The successor experiment is
`../chr20_adaptive_unanimous/`, which requires unambiguous repeat edges and
retains the existing minimum support. The original frozen results here are
preserved, including the regression.
