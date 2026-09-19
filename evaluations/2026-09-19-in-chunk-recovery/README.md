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

## Two follow-ups: a broken precondition, and the read-block anomaly explained

**The merge violated an invariant the cross-chunk stitch depends on.** The stitch
pairs reads between adjacent chunks with a MERGE-JOIN over their read vectors
(`populate_graph_chunk_pair_overlap_impl`, `graph_bam_adapter.cpp:1112-1120`):
two indices advancing on `qname.compare`, which is correct only if both vectors
are sorted by qname. `build_graph_chunk` produces them that way; appending the
alignment-only reads left exactly **one** inversion per chunk (measured), at the
junction between the original block and the appended one. A merge-join does not
fail on that -- it silently stops matching past the inversion, so those reads
never pair and the stitch loses the overlap evidence it votes on.

The merge now restores the qname ordering before rebuilding the index, moving
`read_var_profile` (and each `read_id`), `haps` and `phase_sets` with the reads.

chr20: discordant 2,482 -> 2,456, read hamming 1.160% -> 1.148%, VCF blocks
267 -> 265, tagged 213,905 -> 213,877, wall 138 s -> 126 s.

This is the second invariant of the same family, so they are worth listing
together. Merging into a LIVE chunk means restoring everything derived from it:

| invariant | depended on by | how it failed |
|---|---|---|
| candidates sorted by position | the solve's sweep walks by index | handled in the merge |
| `site_ids` / `site_meta` / `site_allele_orig_idx` parallel to candidates | the VCF emitter | handled in the merge |
| `read_var_cr` keyed by candidate index | the solve finds each site's reads through it | 4,635 blocks; 84 after the rebuild |
| `reads` sorted by qname | the cross-chunk stitch's merge-join | lost overlap votes, silently |

None are asserted or documented where a chunk would be mutated, so each was found
by its symptom.

**The read-block rise is not a defect.** Read blocks go 281 -> 310 while VCF
blocks fall 285 -> 265, and the ordering fix does not account for it (312 -> 310).
Counting phase sets that carry tagged reads against those carrying phased VCF
records:

| | read blocks | VCF blocks | blocks with tagged reads but NO VCF record |
|---|---:|---:|---:|
| post-hoc | 281 | 285 | 0 (0 reads) |
| in-chunk | 310 | 265 | **54 (3,259 reads)** |

The recovered in-gap sites tag reads in the phased BAM without emitting phased
records into the catalog-shaped VCF. So the extra read blocks are recovered
coverage the VCF does not describe, not fragmentation. Whether those sites should
also be emitted is a separate question: they are alignment-discovered and absent
from the catalog the graph arm's VCF is built from.

## Audit of the merge: three more bugs, one of them older than this work

Asked to check what the recovery fires on, what it changes, and whether the
merged sites are represented and used correctly, with the requirement that the
cross-chunk stitch stay invariant. Probing the merge before and after, per chunk,
over `chr20:1-5,000,000`:

| quantity | before merge | after (as shipped in bb34e9b) |
|---|---:|---:|
| candidates | 6,208 | 7,584 |
| positions carrying more than one candidate | 246 | 260 |
| reads whose span lies entirely outside the chunk | **0** | **693** |
| merged sites with usable metadata | -- | **0 of 1,377** |

**1. Recovery regions were not clamped to the chunk.** A group's extension can
reach past the chunk's own span, and discovering there pulls in reads belonging
to the neighbour -- 693 of them over 5 Mb, from zero. Those reads then enter the
cross-chunk stitch's merge-join as though they straddled the boundary, changing
the evidence it votes on. Regions are now clamped to `[chunk.ref_beg,
chunk.ref_end]`; the neighbour recovers its own territory. Reads outside: 0.

**2. The merged sites never reached the VCF, and that predates this work.**
`build_graph_chunk` never fills `chunk.ref_seq` -- it passes a reference slice to
the noise filter and keeps nothing -- while the alignment path does
(`bam_digar.cpp:1336`). The metadata synthesis read the graph chunk's empty
slice, so every merged site got an empty REF, and an empty REF makes the writer
skip that record (`graph_collect.cpp:99`). All 1,377 merged sites over 5 Mb were
affected. This is the mechanism behind the 54 phase sets carrying 3,259 tagged
reads and no records reported above. Reference bases now come from the discovered
alignment chunks, which carry them.

**3. The synthesized alleles were not in VCF form.** A `VariantKey` is the
TRIMMED form: an insertion has `ref_len = 0` with the anchor base at `key.pos - 1`
and `alt` = the inserted bases only; a pure deletion has `alt = ""`. Pairing a
reference slice at `key.pos` with `key.alt` describes an insertion as a
substitution and gives a deletion an empty ALT. The synthesis now builds the
anchored VCF form per type and keeps only what round-trips back through
`vcf_to_variant_key` to the key it came from. Merged sites with malformed indel
alleles: 0 of 1,377.

A fourth issue was introduced and caught during the audit: skipping a site that
failed the round-trip desynchronised `site_meta` from `candidates` (6,208 against
7,584), and since the writer addresses metadata BY CANDIDATE INDEX that shifts
every later site's metadata onto the wrong candidate. An unusable site now gets
an EMPTY metadata entry instead, which the writer skips on its own while the
arrays stay parallel.

## Result, chr20 at -t 16

| | tagged | read blocks | discordant | read hamming | VCF records (phased / blocks) | read-only blocks |
|---|---:|---:|---:|---:|---|---:|
| post-hoc | 203,751 | 281 | 2,732 | 1.341% | 56,032 (55,907 / 285) | 0 |
| in-chunk, pre-audit | 213,877 | 310 | 2,456 | 1.148% | 56,006 (55,867 / 265) | 54 |
| in-chunk, audited | 219,059 | 323 | 2,543 | **1.161%** | **63,897 (59,948 / 316)** | 12 |

Against the post-hoc default: 15,308 more reads phased, 189 fewer misplaced,
read hamming 1.341% -> 1.161%, and 4,041 more phased records -- the in-gap sites
appearing in the VCF for the first time. The VCF block count rises 285 -> 316
because those in-gap phase sets are now visible as records rather than existing
only as read tags; read-only blocks fall 54 -> 12.

Against the pre-audit in-chunk run the accuracy is slightly worse (1.148% ->
1.161%) and that is the clamp: confining discovery to the chunk removes evidence
the unclamped version was using. It is kept because the alternative is feeding
the cross-chunk stitch reads that never straddled the boundary, which is the
invariant the audit was asked to protect.

On `chr20:1-10,000,000`: 37,939 tagged, 40 read blocks, 102 discordant, 0.269%,
10,322 records (9,932 phased / 42 blocks), 0 read-only blocks. The default path
remains byte-identical with the flag off. Unit 4/4, window 66/66, predicate
130/130.

Still open: 12 read-only blocks chromosome-wide, and the 14 extra position
collisions the merge introduces (246 -> 260 over 5 Mb), which are positions where
a catalog site and a merged alignment site coexist.

## Second audit round: two more bugs, one non-bug

**4. Most merged sites emitted nothing, because the writer reads a depth vector
they did not have.** `graph_chunks_to_candidate_table` emits one record per entry
of `counts.alle_covs` (`for new_a = 1 .. alle_covs.size()`), which is the graph
path's per-allele depth vector. The alignment path reports depth as
`ref_cov`/`alt_cov` and leaves `alle_covs` empty on most candidates: **952 of
1,377** merged sites over `chr20:1-5,000,000`. Those sites were phased, tagged
reads, and produced no record. A merged site is biallelic here -- one synthesized
alt, orig index `{0, 1}` -- so the two depths are exactly that vector, and it is
now filled in. All 1,377 have it; chromosome-wide phased records rose 59,948 ->
61,659 and read-only blocks fell 12 -> 4.

**5. Emitting every merged site turned the arm into a variant caller.** With the
depths filled in, chr20 produced **74,630 records against 56,032** -- and the
extra 12,846 were `1|1` homozygous calls with no phase set (125 such records
before, 12,971 after). They come from the alignment's in-gap discovery, which
calls hom variants too. Writing them is worse than not: they would appear ONLY
inside recovery windows, a biased subset of the genome, in a VCF whose contract
is the catalog's sites plus what recovery phased. A merged site is now written
only when it carries phase (`hap_to_cons_alle[1]` and `[2]` both set and
different). Records 74,630 -> 62,197, hom 12,971 -> 547, phased records and read
accuracy unchanged.

**Not a bug: the duplicate candidates.** Identical-key duplicates -- same
position, type, ref_len and alt -- measured before and after the merge in each
chunk: **166 -> 166** over 5 Mb. The merge introduces none; they are a
pre-existing property of the catalog's own decomposition. The position-level
collisions counted in the first audit round (246 -> 260) are distinct variants
sharing a position, not duplicate representations.

## Final state, chr20 at -t 16

| | tagged | read blk | discordant | read hamming | VCF records | phased | blocks | hom | read-only blk |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| post-hoc | 203,751 | 281 | 2,732 | 1.341% | 56,032 | 55,907 | 285 | 125 | 0 |
| in-chunk | 219,059 | 323 | 2,543 | **1.161%** | 62,197 | **61,650** | 324 | 547 | 4 |

15,308 more reads phased, 189 fewer misplaced, 5,743 more phased records. The
default path stays byte-identical with the flag off. Unit 4/4, window 66/66,
predicate 130/130.

Still open: 4 read-only blocks, and 547 hom records against the default's 125 --
422 more than expected once merged hom sites are excluded, which is unexplained.

## Third round: the two items left open

**The 422 extra hom records were manufactured by the depth synthesis.** The
writer reclassifies every record from depth (`graph_collect.cpp:196`):
`is_hom_alt = (ref_cov == 0 && alt_cov >= min_alt_depth)` -> `CleanHom` -> `1|1`.
An alignment candidate merged from a gap frequently carries `ref_cov = 0`
(measured: `ref_cov 0, alt_cov 57`), so filling `alle_covs = {ref_cov, alt_cov}`
made the writer call it homozygous. One wrong theory was tested first and discarded: that the sites failed an
allele index check, since `hap_to_cons_alle` can hold an index of 2 on a site
collapsed to biallelic. Gating on index validity produced byte-identical
output, and the GT construction at `graph_collect.cpp:240-245` turns out to
handle an out-of-range consensus index correctly, so that was not the cause.
The attribution was then settled by an env-gated probe suppressing every
merged site at emission: hom fell 547 -> 139, against 125 for the default,
which locates the records in the merged sites and leaves no contribution from
catalog candidates at the same positions.

A merged site is now written only when the writer's own classification calls it
a het. chr20: records 62,197 -> 61,789, hom **547 -> 139**, phased records
(61,650), tagged reads, blocks and read accuracy all unchanged.

**The read-only phase sets are accepted, with the cost of removing them
measured.** All 17 candidates behind them are merged sites, 16 of them LOW_COV or
LOW_AF -- sites the writer never publishes, which still take a phase set and tag
reads. Excluding them at admission (predicting the writer's own depth rule, since
the category carried at merge time is not the one the writer uses) does remove
them, and costs more than it buys:

| | read-only blocks | VCF blocks | discordant | read hamming |
|---|---:|---:|---:|---:|
| merged sites admitted | 4 | 324 | 2,543 | 1.161% |
| unpublishable ones withheld | 1 | 335 | 2,557 | 1.168% |

Three fewer read-only phase sets for 11 blocks of contiguity and 14 more
misplaced reads. Contiguity and accuracy are the deliverables and a read tagged
with a phase set no record describes is a cosmetic inconsistency, so the sites
stay in and the four blocks stay. The filter is left in the source as a comment
recording the measurement.

## State at the end of the audit, chr20 at -t 16

| | tagged | read blk | discordant | read hamming | VCF records | phased | blocks | hom |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| post-hoc | 203,751 | 281 | 2,732 | 1.341% | 56,032 | 55,907 | 285 | 125 |
| in-chunk | 219,059 | 323 | **2,543** | **1.161%** | 61,789 | **61,650** | 324 | 139 |

Six bugs found and fixed across three rounds (unclamped regions, missing
reference source, wrong allele form, absent allele depths, unscoped emission,
manufactured hom calls), plus one introduced and caught (metadata desync). One
item accepted with its cost measured (4 read-only phase sets). Default path
byte-identical with the flag off; unit 4/4, window 66/66, predicate 130/130.
