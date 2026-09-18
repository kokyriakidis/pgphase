# Recovering sites where reads stay unphased: what the evidence constrains, and the design it points to

Everything below is measured on the one arm of record, `collect-graph-variation
--sites --gaf --bam` on chr20 at `-t 16`, against read-level truth.

## The baseline any change must not damage

| | whole chr20 | chr20:1-10,000,000 |
|---|---:|---:|
| wall | 118 s | 16.3 s |
| phased VCF records | 56,032 | 9,813 in 47 blocks |
| reads tagged | 203,751 | 36,645 |
| read blocks | 281 | 44 |
| discordant reads | 2,732 | 97 |
| read hamming | **1.341%** | **0.265%** |
| blocks bridged from the alignment | 126 | 36 |
| targeted regions / bp solved | 356 / 29.78 Mb | 72 / 3.98 Mb |

## Six things today's measurements settled

1. **A window's own reads cannot vote.** `select_stitch_orientation` needs reads
   carrying BOTH a parent haplotype and a sub-solve haplotype, and the reads
   inside a window are exactly the ones the parent left unphased. The main chunk
   loop needs no padding because a read straddling a chunk boundary is already in
   both chunks (`initialize_chunk_overlap_state`); a window has no equivalent.
   So a recovery region must reach out to where the parent did phase.

2. **The reach should be counted in sites, not base pairs.** Extension needed to
   reach parent phased sites on both sides, over 9,813 parent sites in 10 Mb:
   one site median 0.8 kb / p90 7.3 kb, three sites median 5.6 kb / p90 16.9 kb,
   max 47.3 kb. A fixed 30 kb was ~5x too large in the median case and too small
   for one of 41 regions, which therefore had no anchor at all.

3. **The recovery is not slow per bp; it solves too many bp.** 3.2 s/Mb for the
   first pass against 3.7 s/Mb for the recovery, but 4.25 Mb of actual gap and
   seam became 29.78 Mb solved. Cost is span, not speed.

4. **The reconciliation layer is where the defects are.** Grafting a detached
   sub-solve onto finalised state needs a block merge, an orientation vote and a
   phase-set relabel. Two defects came out of exactly that machinery today: a
   merge flip that never reached the emitted genotypes (6,168 chr20 records were
   anti-phased, fixed in `cc0dc05`), and a merge chaining onto a phase-set label
   another chunk had already retired (open -- 6 parent blocks split in 10 Mb).

5. **Merging the sites into the chunk is nearly sound; re-solving it naively is
   not.** With the alignment's in-gap sites merged into the graph chunk and NO
   re-solve, the output went 47 -> 59 blocks. With a bare
   `assign_hap_based_on_germline_het_vars_kmeans` after the merge it went to
   4,624 blocks. The fragmentation is in the re-entry, not the merge.

6. **The two channels name sites differently.** A catalog candidate's `key.alt`
   is an allele walk (`>114849551>114849554`); the alignment's is sequence
   (`T`). Raw key matching found 0 of 198 sub-solve candidates, including 0 of
   152 SNPs. The bridge is the parallel `site_meta` array plus
   `vcf_to_variant_key`, which applies the same anchor trimming
   `inject_graph_sites` uses in the other direction. With it, only ~64 of ~2,844
   merged sites were duplicates -- in-gap sites really are new, not renamed.

## The design these point to

Recover in the chunk, retry the chunk's own phasing, and keep the detached path
only for what genuinely cannot be done in-chunk.

### Phase 1 -- discovery stays where it is, adoption moves in-chunk

Keep using the existing targeted machinery to DISCOVER: `process_chunk` over the
site-anchored region gives in-gap candidates and, importantly, per-read alleles
at both the in-gap sites and the flanking catalog sites it re-derives. Discard
its phase sets. Nothing about its solve is used, so the sub-solve's own gauge --
the thing that had to be stitched -- never enters the pipeline.

### Phase 2 -- merge, in position order, with the key translation

Insert the discovered candidates into `chunk.candidates` in key order (NOT
appended: the solve's outward sweep walks candidates in INDEX order, so tail
appends make an in-gap site adjacent to the chunk's last site rather than to its
positional neighbours), rebuild `read_var_profile` over the new index space, and
extend `site_ids` / `site_meta` / `site_allele_orig_idx` in the same order --
synthesizing metadata per new site, without which the emitter's index-parallel
lookup skips it even though it phased. Match against the catalog through
`vcf_to_variant_key` on `site_meta`, so an alignment read's observation at a
catalog site lands on the catalog candidate. That shared observation is what
replaces the stitch vote: the read carries alleles on both sides of the gap, so
one solve places it.

Also add the reads the catalog never had, with their alleles at both the in-gap
and the matched catalog sites. They are the ones that cross the gap.

### Phase 3 -- retry the chunk's phasing the way the first pass ran it

This is the step the reverted attempt got wrong. The graph arm's blocks are not
produced by the k-means alone: `populate_graph_chunk_overlaps` and
`stitch_chunk_haps` run after it, and the phase-set numbering comes from that
sequence. Re-entry must replay the same sequence the first pass used for that
chunk, not call the k-means directly -- 4,624 blocks against 47 is what calling
it directly costs.

Concretely: re-run the chunk's phase step, then let the normal overlap/stitch
path run as it already does for every chunk. The retry then differs from the
first pass in exactly one respect -- more sites are admissible inside the
windows -- which is the whole intent.

### Phase 4 -- admission, escalated

Admit in two tiers, the boundary being repeat context rather than category:

- **Tier 1**: everything except repeat-context indels. The homopolymer detector
  is trustworthy again as of `510f865` (its insertion branch had never fired;
  fixing it took chr20 read hamming from 3.665% to 0.658% on the alignment arm),
  so `is_homopolymer_indel` is a usable tier boundary now in a way it was not
  this morning.
- **Tier 2**, only where the junction is still open after tier 1: the repeat
  indels too. The machinery for this already exists in `collect_noisy_vars_step4`
  -- a conservative pass, `noisy_region_still_broken` as the gate, and a second
  pass admitting the weaker source -- and is unreachable only because its
  predicate was derived from options removed with the whitelist mode and
  hardcodes `false`.

Reviving that predicate is a smaller change than it looks and is independent of
the rest of this plan, so it can be measured on its own.

### Phase 5 -- what stays on the detached path

Seams that cross a chunk boundary. `collect_block_seams` derives them from
solved candidates, and a cross-chunk seam does not exist until stitching has
happened, so no in-chunk retry can see it. That residual set keeps the current
sub-solve path, which means the open stale-label defect still needs its fix:
resolve `keep_ps` through to its current label before merging, and relabel across
chunks rather than within one. The prewarm log already reports how many regions
that residual set is, so its size is known before the work starts.

## Order of work, and the gate at each step

Each step is measured on `chr20:1-10,000,000` first (16.3 s, 36 bridged, 44
blocks, 97 discordant, 0.265%), then whole chr20 (118 s, 126 bridged, 281
blocks, 2,732 discordant, 1.341%). The criterion is the one this project has
used all along: read accuracy must not get worse, and yield must not fall.

1. Revive the escalation predicate (tier 1 / tier 2). Independent, smallest, and
   testable on the existing path.
2. Merge + replayed chunk phasing on ONE chunk, behind a flag, with the block
   count as the first thing to look at -- 47 blocks is the number that says the
   re-entry is right.
3. Extend to all chunks, measure, and only then consider making it default.
4. Fix the stale-label chaining on the residual cross-chunk path.

## What would falsify this

If the replayed chunk phasing still fragments, the assumption that phase-set
numbering is a property of the chunk's phase sequence is wrong, and the in-chunk
route is dead -- at which point the honest fallback is the current detached path
plus the stale-label fix, and the effort moves to admission instead.

If merging sites in-chunk improves nothing on accuracy, the value was never in
the reconciliation; it is in which sites are admitted, and Phase 4 is the whole
project.
