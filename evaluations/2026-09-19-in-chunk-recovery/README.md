# Recovery inside the chunk, before the stitch

Chunks do disjoint work, so recovery belongs in the chunk worker rather than in a
pass bolted on after everything is finished. `--in-chunk-recovery` does that: each
worker, having solved its chunk, collects its own failed windows, discovers the
alignment's in-gap sites over them, merges those sites into the chunk, and
re-solves it -- all before `populate_graph_chunk_overlaps` and
`stitch_chunk_haps` run.

Nothing is grafted onto finished state, so the reconciliation layer disappears:
no block merge, no orientation vote, no phase-set relabel. Both of 2026-09-18's
defects lived in that layer.

## Result, graph + recovery arm, chr20 at -t 16

| | wall | tagged | read blocks | discordant | read hamming | VCF records | VCF blocks |
|---|---:|---:|---:|---:|---:|---:|---:|
| post-hoc (default) | 118 s | 203,751 | 281 | 2,732 | 1.341% | 55,907 | 285 |
| in-chunk | 138 s | **213,905** | 312 | **2,482** | **1.160%** | 55,863 | **267** |

10,154 more reads phased, 250 fewer misplaced, and a more contiguous VCF -- 267
blocks against 285 -- for 20 s. 124 of the chunks recovered something. The
default path is byte-identical with the flag off.

On `chr20:1-10,000,000`: 37,145 tagged / 39 read blocks / 101 discordant /
0.272% / 38 VCF blocks, against the post-hoc path's 36,645 / 44 / 97 / 0.265% /
47.

## Three bugs found getting here, in order

**1. The insertion point was not the problem.** Moving the merge from after the
stitch to inside the worker changed nothing: 4,624 blocks became 4,635. That
killed the explanation and forced the real one.

**2. The stale read<->variant index was.** `chunk.read_var_cr` is a cgranges tree
keyed by CANDIDATE INDEX, and the solve looks up each site's reads through it
(`collect_phase.cpp:589`, `:983`). Merging re-indexes every candidate, so the
tree was stale and the sweep read the wrong reads for every site. `build_graph_chunk`
and the injection path both rebuild it; the merge did not. One call to
`rebuild_read_var_cr` took the 10 Mb slice from **4,635 blocks to 84**.

**3. The merged sites were being ignored.** They are `NOISY_CAND_HET`, and the
graph arm's k-means runs over `kCandGermlineClean` only, so the sites were merged
and then excluded from the solve. Running the alignment pipeline's two rounds --
clean, then `kCandGermlineVarCate` -- admits them: 84 blocks -> 38, and tagged
reads 36,732 -> 37,145.

## The anchored second round must be OFF on this path

Measured on the 10 Mb slice, in-chunk, second round over `kCandGermlineVarCate`:

| second round | tagged | read blocks | discordant | read hamming |
|---|---:|---:|---:|---:|
| anchored | 37,125 | 51 | 1,412 | **3.803%** |
| reset | 37,145 | 39 | 101 | **0.272%** |

This is the opposite of the post-hoc measurement, where anchoring took chr20 from
1.399% to 1.341% and is now the default. The difference is what stage 1 knows: on
the post-hoc path it has already seen every site it will ever see, so pinning its
gauge protects a good answer. In-chunk, stage 1 runs BEFORE the in-gap sites and
the reads that carry them exist, so pinning freezes a gauge derived without them
and the newly added reads are labelled against it. Whatever the mechanism in
detail, the measurement is unambiguous and the arm runs with
`--no-anchored-stage2`.

Not yet explained: read blocks rise 281 -> 312 chromosome-wide while VCF blocks
fall 285 -> 267, and records fall by 24.
