# The graph pipeline, with alignment recovery in the gaps

`collect-graph-variation --bam FILE` runs the graph pipeline exactly as before --
the catalog's sites phased from GAF evidence, no alignment pass -- and then uses
the alignment only for what it could not phase.

```
pgphase collect-graph-variation \
  --ref ... --sites chr20.sites.striped.vcf.gz --gaf ...coord.gaf.gz \
  --bam HG002...bam \
  -r 'CHM13#0#chr20:5259406-5395085' \
  -o candidates.tsv --phased-vcf-out native.vcf --phased-bam-out phased.bam
```

The recovery is the hybrid's own, reached through one new entry point,
`recover_unphased_windows_from_bam`. It takes both kinds of first-pass failure:
`collect_unphased_windows` for intervals where reads were left unphased, and
`collect_block_seams` for intervals between consecutive blocks -- the mode a
graph-driven pass produces, every read placed but the blocks unjoined. Each is
re-solved as its own chunk from the alignment at the recovery mapq floor and
stitched in on shared reads.

## Result

| window | arm | phased hets | blocks | spans | in-gap | reads tagged | concordance |
|---|---|---:|---:|---|---:|---:|---:|
| 5,309,406 | graph alone | 23 | 3 | no | 0 | 422 | 99.53% |
| | **+ recovery** | 23 | **2** | **YES** | 0 | 422 | 99.53% |
| 26,029,591 | graph alone | 103 | 2 | no | 0 | 293 | 99.66% |
| | **+ recovery** | 103 | **1** | **YES** | 0 | 293 | 99.66% |

Both gaps span where the graph pipeline alone left them open, no phased
heterozygote is lost, and read concordance is unchanged.

## Two integration bugs found, and one limitation left

**The contig id is not shared across the boundary.** The graph pipeline builds a
synthetic BAM header from the reference index, so its tid for a contig is not the
BAM's. Passing the chunk's tid straight into the targeted solve made it fetch the
wrong contig -- `failed to fetch reference contig: CHM13#0#chr10` where tid 0 was
chr20 in one header and chr10 in the other. The entry point now resolves the
contig by NAME in the BAM's header, and candidates coming back are relabelled to
the chunk's id.

**Adopting the solve's counts destroys good calls when the parent is not
starved.** The count overwrite exists for the mapping-quality hole, where the
parent holds injected catalog sites with almost no read support. Applied to a
GAF-derived call it is destructive: `chr20:5,393,876` went from `CLEAN_HET_SNP`
at DP 67, 27/40 to `LOW_AF` at DP 45, 1/44. It is now conditional on the parent's
counts actually being unusable -- `LowCoverage`, or total coverage below
`min_depth`. The phasing is always adopted; the evidence is not.

**Limitation: in-gap sites are not emitted yet, so `allow_import = false` here.**
The graph path's output table is index-parallel to per-site metadata --
`GraphChunkBuildResult::site_meta`, `site_ids`, `site_allele_orig_idx`, all
indexed by candidate position in `chunk.candidates` and skipped outright at
`graph_collect.cpp:71` when the index runs past `site_meta`. Recovery appends
candidates and then reorders the table, which decouples the arrays: measured with
import on, 33 candidates fell to 24 and 23 phased heterozygotes to 14 on
chr20:5,309,406, with records mispaired rather than merely dropped.

Adoption alone is index-safe -- it mutates candidates in place and adds none --
which is why the arm bridges blocks and loses nothing. Emitting the alignment's
in-gap sites needs the parallel arrays kept in step: a synthesized `GraphSiteMeta`
per appended candidate and no reordering. That is the next step for this arm.

Not measured: chromosome-wide.

## Against the hybrid

Same two windows, same reads, scored the same way. Note the graph path writes an
**unaligned** BAM, so every record reads as unmapped: filtering on that drops all
of them and reports 0 tagged. Reads carrying `HP` are counted regardless.

| window | arm | phased hets | blocks | spans | in-gap | reads tagged | concordance | discordant |
|---|---|---:|---:|---|---:|---:|---:|---:|
| 5,309,406 | graph alone | 23 | 3 | no | 0 | 422 | 99.53% | 2 |
| | graph + recovery | 23 | **2** | **YES** | 0 | 422 | 99.53% | 2 |
| | hybrid (default) | **36** | **1** | YES | **2** | **544** | 98.71% | 7 |
| 26,029,591 | graph alone | 103 | 2 | no | 0 | 293 | 99.66% | 1 |
| | graph + recovery | 103 | **1** | **YES** | 0 | 293 | 99.66% | 1 |
| | hybrid (default) | **436** | 2 | YES | **278** | 293 | 99.66% | 1 |

**The recovery mechanism performs well and is free.** Both gaps go from open to
spanned and every other column is unchanged from graph alone -- same phased
heterozygotes, same reads tagged, same concordance, same discordant count. That
is the mechanism doing exactly what it should.

**The graph-alone base is the weaker part.** It phases 23 against the hybrid's
36, and 103 against 436 -- four times fewer on the larger window. On
chr20:26,029,591 the hybrid reaches those 436 with *identical* read numbers (293
tagged, 99.66%, 1 discordant), so its extra 333 phased heterozygotes and 278
in-gap records cost nothing at the read level.

**One column favours the graph arm.** On chr20:5,309,406 it reads 99.53% with 2
discordant against the hybrid's 98.71% with 7, for 122 fewer reads tagged --
fewer reads placed, more purely. Whether that trade is worth taking is a choice,
not a defect.

**And the in-gap zeros are the arm's limitation, not the architecture's.** They
are the index-coupling blocker above, with `allow_import = false`. Until a
`GraphSiteMeta` is synthesized per appended candidate, this arm can bridge a gap
but cannot report a variant inside one -- which is most of what closing a gap is
for.
