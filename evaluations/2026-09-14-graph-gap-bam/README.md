# Graph-selected BAM gap recovery trials

Follow-up: the endpoint-coordinate checker was corrected and a verified
deletion bridge was fixed. See [the follow-up evaluation](../2026-09-14-gap-deletion-bridge/README.md).
The original results below preserve the earlier run; its three unresolved
endpoints comprise two coordinate-matching errors and one unphased insertion.

The trial runner compares the existing recovery (`--no-graph-gap-bam`) with
the graph-selected BAM fallback. Both arms use the same binary, inputs and
per-region MSA evidence cache. Truth is passed only to the evaluator.

```bash
python3 scripts/trial_graph_gap_bam.py \
  --targets evaluations/2026-09-14-graph-gap-bam/targets.tsv \
  --ref test_data/chm13v2.0.chr20.renamed.fa \
  --bam test_data/HG002_chr20_hifi_mapped_to_CHM13_chr20_annotated.bam \
  --graph-sites test_data/chr20.sites.striped.vcf.gz \
  --gaf test_data/HG002.chr20.annotated.coord.gaf.gz \
  --truth-bam ../pgphase-eval-data/truth/chr20/diplinator_merged.bam \
  --output /tmp/pgphase-gap-trials --run-name graph_bam_v2 \
  --threads 2 --jobs 2
```

Use a new `--run-name` after each code change. The runner refuses to overwrite
an existing run. Evidence caches persist under `--output/cache`; mismatched
inputs are rejected by pgphase. The manifest records the exact target panel,
arguments and executable SHA256. Commands, logs and parental read evaluations
are retained under each region. A failed region produces a nonzero exit status
and is recorded in `failures.json` without cancelling other regions.

Audit selection requires `split_block`, a concordant competitor endpoint pair,
and competitor block accuracy >= 0.99. The current audit yields 11 unique gaps.
A custom fixed panel can instead supply `chromosome`, `left`, and `right` TSV
columns. Coordinates are one-based. Default input chromosome names are `chr20`;
`--contig-prefix` supplies the BAM/graph contig prefix.

Review `results.tsv` for remaining splits, unresolved endpoints, new joins,
read Hamming errors, switch/flip errors and runtime. An unresolved endpoint is
not a successful join. Cold-cache baseline timing includes MSA discovery;
compare warm runs when measuring performance. Local windows may phase differently
from whole chromosomes, so confirm accepted changes on the complete chr20 run.
Do not automatically accept a join merely because contiguity increases or
reject one solely for a few extra discordant reads: check the orientation of
the original blocks against parental truth.

The final chromosome trial joined 122 of 280 gaps (previously 121). The new
13,894,047–13,941,253 relationship was correctly unflipped against parental
read truth. Read discordance remained 2,731 / 189,001. Cached recovery took
195.14 seconds, versus approximately 124 seconds before the additional pass.
These are read-truth results, not newly measured shared-VCF NGC50 statistics.

The 11-region paired panel completed without failures or additional discordant
reads. Both arms have eight split targets and three unresolved endpoints. In
the latter cases, the exact audited endpoint is absent from the native candidate
table (12,721,112; 23,480,815; 37,984,529). These need candidate/representation
analysis before interpreting them as stitching failures. Warm fallback phasing
takes 0.78–1.25 seconds per region, excluding truth evaluation. The panel is a
baseline for subsequent fixes; it does not demonstrate recovery of these 11
remaining targets. Results are preserved in `results.tsv`.
