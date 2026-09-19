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
`--graph-authoritative`, `--private-sites` and `--bam-authoritative-bed`.

## The graph arm with `--bam`

`collect-graph-variation --bam` is the arm under active work, and it is not the
hybrid: the catalog's sites are phased from GAF evidence, and the alignment is
consulted ONLY to recover what those sites could not join.

```
pgphase collect-graph-variation \
  --ref chm13.fa --sites sites.vcf.gz --gaf reads.coord.gaf.gz --bam reads.bam \
  -r 'CHM13#0#chr20' -o candidates.tsv \
  --phased-vcf-out phased.vcf --phased-bam-out phased.bam
```

Two recovery placements, and they are a real choice:

| | where | flag |
|---|---|---|
| post-hoc (default) | after the pass, as its own sub-solve grafted onto the parent blocks | -- |
| in-chunk | inside each chunk, before the stitch | `--in-chunk-recovery` |

In-chunk recovery merges the alignment's in-gap candidates into the LIVE chunk
and re-runs both solve rounds over the union, so the recovered sites are ordinary
members of the chunk's own solve rather than a graft whose internal parity
nothing checks. `--no-anchored-stage2` makes stage 2 reset and re-solve over the
wider site set instead of refining stage 1; on this path the resetting form
measures better, which is the reverse of the alignment arm and is not explained.

chr20 at `-t 16`, scored against read-level parental truth:

| | wall | tagged | read blocks | discordant | read hamming | VCF records | phased | hom |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| post-hoc | -- | 203,751 | 281 | 2,732 | 1.341% | 56,032 | 55,907 | 125 |
| in-chunk, `--no-anchored-stage2` | 114 s | 219,059 | 323 | **2,543** | **1.161%** | 61,789 | **61,650** | 139 |

### What the merge must preserve

Mutating a live chunk means restoring every property the rest of the pass
assumes, and six defects shipped because nothing checked them.
`verify_chunk_invariants` (`collect_phase.cpp`) runs at the end of the merge and
throws, naming the first violation:

- the three per-site arrays (`site_ids`, `site_meta`, `site_allele_orig_idx`) as
  long as `candidates` -- the writer addresses metadata BY CANDIDATE INDEX;
- `read_var_profile`, `haps` and `phase_sets` as long as `reads`;
- `candidates` position-sorted, and `read_var_profile[i].read_id == i`;
- `reads` qname-sorted -- the cross-chunk stitch pairs them with a merge-join,
  so one inversion makes it skip everything past that point;
- the re-solved region inside `[chunk.ref_beg, chunk.ref_end]` -- discovering
  outside pulls in the neighbour's reads, which then enter the stitch's vote.

### What a merged site contributes

A merged site takes part in the solve like any other, but it is WRITTEN only
when this writer's own classification (`graph_collect.cpp:196-223`, which
reclassifies from depth) calls it a het. Two reasons: the alignment's in-gap
discovery also calls homozygous variants, and a merged candidate often carries
`ref_cov = 0`, which that reclassification reads as homozygous. Emitting either
would add calls that appear only inside recovery windows -- a biased subset of
the genome -- to a VCF whose contract is the catalog's sites plus what recovery
phased.

Four phase sets chromosome-wide still tag reads without a record describing
them. Their sites are merged candidates the writer classifies LOW_COV or LOW_AF;
withholding those at admission removes three of the four and costs 11 blocks of
contiguity and 14 more misplaced reads, so they are left in deliberately.

## Test gates

`make unit-tests` (4 binaries), `make window-tests` (the committed chr20 gap
windows) and `make predicate-tests` (the phasing predicates and the chunk
invariants). All three must pass before a commit.

There is no injection suite: `src/test_bam_site_injection.cpp` was deleted in
c092785 when the tests were narrowed to the window under work. A compiled binary
outlived it in some working trees and kept printing pass/fail counts from
retired expectations; it is not a build target and its output means nothing.
