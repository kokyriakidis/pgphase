# The pipeline

One mode: **the catalog's sites phase, the alignment recovers the gaps.** There
is no configuration that selects a different architecture.

```
pgphase collect-hybrid-variation \
  --ref chm13.fa --bam reads.bam \
  --graph-sites sites.vcf.gz --gaf reads.coord.gaf.gz \
  -r 'CHM13#0#chr20' -o candidates.tsv --phased-vcf-out phased.vcf -b phased.bam
```

## Per chunk, in order

| step | call | what it contributes |
|---|---|---|
| 1 | `load_and_prepare_chunk` | reads and their digars. Loads down to `min(min_mapq, recovery_min_mapq)`; anything below `min_mapq` is parsed and immediately marked skipped, so every stage behaves as if it were absent |
| 2 | `collect_var_classify` | alignment discovery, allele counts, the **noisy-region model**, classification |
| 3 | *(withhold)* | `chunk.candidates.clear()` -- the alignment channel's own candidates do not enter the first solve |
| 4 | `load_sites_for_region` → `inject_graph_sites` | the catalog's sites become the candidate table. These are the phasing anchors |
| 5 | `collect_var_build_profiles` | each read's allele at each candidate it overlaps |
| 6 | `inject_graph_reads` | GAF-derived observations for reads the alignment does not carry |
| 7 | `backfill_graph_candidate_counts` | derives every count on a graph-only candidate from the final profile state. The single writer, so double counting is inexpressible |
| 8 | `classify_graph_only_candidates` → `apply_hybrid_noise_filter` | category per injected candidate, then the indel noise screen |
| 9 | `collect_var_run_phasing` | the clean k-means plus the noisy-region MSA. The hybrid's own `skip_noisy_kmeans` keeps the noisy class out of this solve |
| 10 | **recovery** | below |
| 11 | `prune_not_candidate_variants` | drop what no longer qualifies |

Then `stitch_chunk_haps` joins adjacent chunks on shared reads, the read filters
run, and the outputs are written.

## Why step 2 survives step 3

Step 3 discards the candidates step 2 produced, which makes the discovery look
like dead cost. It is not: `chunk.noisy_regions` is derived from those
candidates, that model scopes the noisy-region MSA, and the MSA supplies most of
a chunk's phased sites. Skipping discovery outright on
`chr20:25,979,591-26,138,679` is 31% faster (4.05 s to 2.78 s) and collapses
phased heterozygotes from **447 to 137**.

## Recovery (step 10), on by default

Two failure modes, because they do not overlap:

- `collect_unphased_windows` -- intervals where the solve left **reads**
  unphased, which is what an alignment-driven solve produces;
- `collect_block_seams` -- intervals between consecutive blocks, every read
  placed but the blocks unjoined, which is what a catalog-driven solve produces.

Each window is treated twice, in order:

1. **In place.** The chunk is re-solved with `force_noisy_msa` (asks for the
   noisy-region MSA by name, and enables `split_nested_msa_deletions` with it),
   `skip_noisy_kmeans = false` (so the recalled sites are actually oriented), and
   `retry_windows` as the **only** scoping -- it confines the widened het
   admission to the failed intervals. Admitting that class chunk-wide instead
   roughly doubled the chromosome-wide read Hamming error (0.878% → 1.837%).
2. **As its own chunk.** Whatever is still unphased, plus the seams, goes to
   `recover_windows_with_targeted_solve`: `process_chunk` over the window plus
   one read length at the recovery mapq floor, with the noisy class admitted and
   no further recursion. The result is stitched to each adjacent parent block by
   `select_stitch_orientation` over the reads tagged in both -- which refuses
   when no read is shared, so an interval no read crosses yields two honest
   blocks instead of a coin flip. Sites the parent holds unphased **adopt** the
   sub-solve's phasing; sites it never discovered are **imported**.

## The other two subcommands

`collect-bam-variation` (alignment only) and `collect-graph-variation` (catalog
sites from GAF evidence only) are separate tools, not modes of the hybrid. They
are what the hybrid is measured against; neither takes the other's input.

## What is not a mode

`--no-retry-unphased-with-bam` disables recovery and exists for regression
attribution, not as a supported configuration. Everything else on the hybrid is
a threshold or an output path. Removed as modes: `--recover-gaps`,
`--msa-verified-refine`, `--gap-bam-only`, `--graph-first`/`--no-graph-first`,
`--graph-authoritative`, `--private-sites`, `--bam-authoritative-bed`, and
`collect-graph-variation --bam`.
