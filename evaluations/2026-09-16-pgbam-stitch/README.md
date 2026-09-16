# Can the existing machinery make multiple blocks and stitch them?

Yes to the first half, and there is a dedicated mechanism for the second -- but on
this dataset it stitches at chance, and the read-based stitch it would complement
is not wired for the case the deficit gaps fall into.

## What exists

`stitch_chunk_haps` (`collect_phase.cpp:1658`) is a four-stage design:

1. `stitch_phase_blocks_with_pgbam(chunk, ...)` per chunk -- merges blocks
   **within** a chunk. Gated on `--pgbam-file`.
2. `flip_chunk_hap(prev, cur, opts)` for each adjacent chunk pair -- the read-vote
   stitch, via `select_stitch_orientation`.
3. `stitch_adjacent_chunks_with_pgbam(...)` -- rescues a chunk seam that the read
   vote failed to decide. Gated on `--pgbam-file`.
4. Two more `stitch_phase_blocks_with_pgbam` passes over all chunks, at cleanup
   and relaxed-cleanup thresholds. Both default **on**, gated on `--pgbam-file`.

So the pipeline already produces many blocks (238 VCF phase sets, 191 read blocks
on chr20) and already has a within-chunk block stitcher. Every run in this
investigation omitted `--pgbam-file`, so stages 1, 3 and 4 were inactive, while
the sidecar for this exact BAM has been in `test_data/` all along
(`HG002_chr20_hifi_mapped_to_CHM13_chr20_annotated.pgbam`, 126 MB).

## What the pgbam stitch links on

`decide_phase_block_concordance` intersects the sets of **GBWT haplotype threads**
polarized to haplotype 1 and 2 in each block (`s11`, `s12`, `s21`, `s22`) and
merges when one orientation wins and the winner clears `min_winning_threads`.
That is evidence of a different kind from reads: a graph thread continues through
a 27 kb read-linkage break that no read can cross. It is the right shape of
evidence for these gaps.

## Measured: it is not usable here

Whole chr20, otherwise identical settings (`-q 1`, recovery on), scored per block
against the diplinator read truth under each block's own best orientation.

| arm | VCF blocks | phased Mb | gaps | gap Mb | read blocks | read Hamming |
|---|---:|---:|---:|---:|---:|---:|
| no pgbam (committed default) | 238 | 56.02 | 196 | 11.88 | 191 | **0.878%** |
| `--pgbam-file`, defaults | 2 | 66.21 | 0 | 0.00 | 1 | **47.511%** |
| primary pass only, win/margin 2 | 5 | 63.34 | 2 | 2.86 | 5 | 46.853% |
| primary pass only, win/margin 5 | 19 | 62.75 | 11 | 3.45 | 11 | 41.087% |
| primary pass only, win/margin 20 | 133 | 59.31 | 89 | 7.27 | 92 | **24.050%** |

At its defaults it merges chr20 into a single read block at 47.5% error -- chance.
Tightening the thresholds trades merging for accuracy monotonically and never
approaches the baseline while still merging: even at 20 winning threads with a
margin of 20 the error is 24.05% against 0.878%. The thread sets in this sidecar
are not discriminative for this sample's haplotypes, so intersecting them yields
near-random polarity. (The sidecar parses cleanly -- correct magic and version --
but it was generated months before the current pipeline, so a stale or mismatched
sidecar is not excluded; the mechanism is not condemned by this measurement, only
its usability on this input.)

That also explains the runtime: the default pgbam arm finishes in 51 s against
10.5 min, because once everything is merged into one block there are no gaps left
and gap recovery does nothing (1,313 tier rows -> 1).

## The gap this leaves

`select_stitch_orientation` -- the project's read-vote standard -- is called in
exactly three places: the adjacent-chunk seam in `flip_chunk_hap`, and the two
gap-recovery link votes. **There is no read-vote stitch between two blocks inside
one chunk.** Chunks are 500 kb by default and the deficit gaps run 3-72 kb, so
they sit inside a chunk essentially always: the read-based stitch never applies
to them, and the only mechanism that would is the pgbam path measured above.

This is why gap recovery has its own flank-linking at all, and why it is
all-or-nothing: it asks whether one proposal block links both flanks, rather than
composing a chain of block-to-block links. The generic missing piece is a
within-chunk, read-vote block stitch at the same `select_stitch_orientation`
standard, of which recovery's flank link is a special case.

Whether that generic stitch would close the 31 recoverable deficit gaps is open
and measurable: `chain_proposal_blocks.py` prototypes exactly this vote on the
recovery proposals and found the seams thin -- adjacent pairs sharing 0-12 reads
that observe consensus sites in both, with one 12-read pair splitting 6/6. It
should be re-prototyped on the main solve's blocks, where the site consensus is
built from hundreds of reads rather than the proposal's dozens, before any C++.

## Change made

The hybrid subcommand accepted `--pgbam-file` but none of the eight threshold
options that `collect-bam-variation` exposes, so the path could only be run at
defaults that produce 47.5% error, with no way to tighten it. Those options are
now exposed in `collect-hybrid-variation` too (`--pgbam-primary-margin`,
`--pgbam-primary-min-winning`, `--no-pgbam-cleanup-pass`,
`--pgbam-cleanup-margin`, `--pgbam-cleanup-min-winning`,
`--no-pgbam-relaxed-cleanup-pass`, `--pgbam-relaxed-cleanup-margin`,
`--pgbam-relaxed-cleanup-min-winning`). Defaults unchanged; the sweep above is
what they make possible.
