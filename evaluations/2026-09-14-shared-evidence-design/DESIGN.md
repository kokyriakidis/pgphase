# Shared graph/BAM/MSA evidence boundary

## Recommendation

Introduce one immutable, source-preserving event-and-molecule representation
before constructing `PhasingChunk`. Freeze the original graph-derived phase
blocks and their gaps first. Only private surjected-BAM events inside those gaps
become new sites in the phasing observation graph. Each gap owns its evidence
view; existing graph sites in its adjacent blocks provide fixed anchors. Private events do not require
editing the GBZ. Keep the existing phaser behind a checked adapter initially;
changing representation and changing the optimizer together would obscure
whether a result improved because information was preserved or because the
algorithm changed.

Implementation status (2026-09-14): the gap-owned evidence snapshot and checked
adapter are implemented in `src/gap_evidence.hpp/cpp`. Separate cached BAM and
graph observations survive injection and MSA replacement; each original gap is
projected independently into the existing phaser. This implements the simple-event
boundary below, not the proposed general complex-event equivalence resolver,
and does not establish that all phasing errors are eliminated.

## Scope: one evidence set per graph phase gap

The user explicitly restricts new BAM candidates to gaps between the original
graph-derived phase blocks. A BAM-only phase gap elsewhere does not trigger
recovery. Gap definitions and anchor memberships are frozen before recovery;
newly recovered blocks do not recursively expand the BAM-discovery territory.

A `GapEvidence` view contains:

- A stable gap ID, reference interval, and the two original graph block IDs.
- Existing graph anchor events/read memberships from the adjacent blocks.
- Clean BAM SNPs/indels and MSA-verified SNPs/indels owned by this gap.
- Per-molecule observations and source provenance at those events, including
  untagged reads. Already tagged reads are not a prerequisite for evidence.
- Context references needed for alignment/MSA, separately from admitted events.

Read retrieval and MSA windows may extend into the flanks to align an event or
identify which graph block a molecule supports. Such context does not authorize
private BAM candidates outside the gap to enter its solve. Existing graph
anchors remain fixed while the gap's evidence is assessed.

Use canonical reference footprints for ownership, not whichever VCF anchor
position a source happened to choose. Simple events wholly inside the gap are
owned by it. Events spanning a gap boundary require explicit mapping against
the existing graph anchor event; initially retain them as boundary/unsupported
records instead of silently replacing a flank or cropping their sequence.
Shared alignment context may be referenced by several gap views, but each
molecule/event observation has one canonical identity and is not multiplied
when views overlap.

Reuse the BAM variation already collected: query/index its sites by these gap
intervals, retaining only the required per-gap candidate sets. Fill a genuinely
missing observation once and cache it. There is no reason to rerun whole-BAM
variant discovery for each gap or stitch trial.

The data flow is:

```text
original graph phase blocks -> frozen gap inventory
                                      |
existing BAM site/observation cache -> gap-owned evidence views
                                      |
                         checked legacy-phaser adapter
                                      |
                   local block-orientation proposal + evidence
                                      |
                    validate and apply accepted block parities
```

Keep discovery tier as provenance (clean / MSA SNP / MSA indel). All available
tiers can live in the same gap record; choosing a subset or weighting does not
require recollection, change ownership, or turn multiple views of a molecule
into independent supporting reads.

## Concrete problems in the current boundary

- `CandidateVariant` (`src/phasing_types.hpp`) combines allele identity,
  aggregate counts, discovery/classification flags, MSA status, current PS/HP,
  and mutable k-means consensus. The gap cache serializes these together
  (`write_gap_cache_candidate` in `src/collect_pipeline.cpp`). Reusing a site
  therefore also reuses derived state unless every caller explicitly resets it.
- `ReadVariantProfile::alt_qi` stores either a BAM query position or a graph
  confirmation sentinel. `inject_graph_reads` in `src/hybrid_inject.cpp` replaces
  the query coordinate even when an existing BAM allele agrees with GAF. Source
  agreement should not erase the coordinate used to retrieve base quality.
- `alleles` is a working BAM/graph/MSA channel; `graph_alleles` preserves GAF
  observations separately, but there is no equivalent independent immutable
  BAM and MSA history. `merge_var_profile` chooses old/new site records and
  remaps parallel observation vectors; it is not a source-preserving evidence
  merge. `select_graph_gap_bam_reads` consequently rescans BAM CIGARs to recover
  SNP observations for a BAM-only view.
- The allele vocabulary is inconsistent: comments describe 0/1 observations,
  general count vectors admit more alleles, and only MSA insertions have an
  explicit list of additional alternate sequences. An integer allele has no
  safe meaning without its owning event's sequence dictionary.
- `vcf_to_variant_key` classifies a length-increasing replacement as insertion:
  AC -> ATG becomes C -> TG with ref_len=1 and type=Insertion. In contrast,
  `slice_msa_site` treats insertion as consuming zero reference bases. A generic
  replacement event needs to retain both its reference span and alternate
  sequence; SNP/INS/DEL is a derived description, not its identity. This is a
  concrete incompatibility in the conventions, not a measured explanation for
  a particular chr20 error.
- Site validation and link eligibility can diverge. The current recovery fix
  closes one example: ordinary biallelic MSA indels previously bypassed a
  computed `gap_link_supported` result. A per-site `msa_verified` flag also does
  not mean every read observation at that site was individually validated.

`provenance_audit.json` measures this limitation in four frozen gap windows.
All matching graph/working observations there carry the query-position
sentinel and cannot export a query base quality directly. Some are graph-only
observations that never had a BAM query position; these counts must not be
interpreted as the number of BAM qualities actually overwritten.

## Minimal data contract

| Record | Required meaning |
| --- | --- |
| Event | Reference identity/version, contig, 0-based half-open replacement interval, explicit reference sequence and alternate sequences, stable event ID |
| Source allele mapping | Source site ID plus source allele ID -> canonical event allele ID; preserve unmappable and ambiguous mappings |
| Molecule | Stable input namespace plus read name; alignment records are children, not additional independent molecules |
| Observation | Molecule ID, event ID, allele ID or explicit unknown/ambiguous state, source, alignment/window provenance, query interval, available quality and mapping quality |
| Derived phasing view | Selected molecule/event observations, recomputed counts, genotype/eligibility result and reasons; separate from raw evidence |
| Phasing result | Event-to-HP orientation, molecule HP/PS and block parity edges; never written back into evidence |

Allele 0 is the canonical reference sequence, even when no read carries it.
A site with two alternate alleles remains genotype 1/2. Allele IDs are local to
an event. A missing observation is not reference; a read spanning a coordinate
is not sufficient evidence for an indel's reference allele. Low quality,
not covered, incompatible alignment and ambiguous repeat placement remain
explicit states with their reasons.

Store query coordinates with an explicit alignment strand/convention; graph and
BAM observations must map to the same original molecule orientation.

Use half-open coordinates inside the new representation and checked conversion
at legacy one-based boundaries. Insertions have an empty replaced interval;
replacements have a nonempty interval and alternate sequence, regardless of
whether the alternate is longer or shorter.

## Encoding and merging

1. Normalize source events against the same reference. Compare spelled local
   haplotype sequences to reconcile repeat-shifted or differently decomposed
   representations. Trimming/left alignment helps simple events but does not
   establish equivalence for all complex graph/BAM descriptions.
2. Preserve aliases to original graph sites and BAM/MSA candidates. Merge only
   proven equivalent allele events. Overlapping events remain related evidence;
   coordinate overlap alone does not justify merging them or counting them as
   independent signals. Retain unsupported complex mappings explicitly.
3. Record BAM, GAF and MSA observations separately for the same molecule/event.
   Their agreement is corroboration on one molecule, not three independent
   reads. Their disagreement must survive selection of a phasing view.
4. MSA discovery produces candidate allele sequences. Each contributing read
   needs an observation supported by its aligned sequence; membership of an
   MSA cluster is not itself an allele call. Retain the MSA window and alignment
   provenance to identify correlated artifacts.
5. Build a deterministic projection into the current candidate/profile format.
   Counts and eligibility derive from exactly the projected observations.
   Rejected or unrepresentable events stay in the evidence store with reasons;
   they do not acquire an approximate legacy encoding silently.

Do not attach a supposedly calibrated probability to graph/BAM/MSA agreement.
Preserve raw quality and context first; establish weights with held-out tests.
Separate genotype support, observation reliability, and relative-block linkage.
None of these three decisions can substitute for the others.

## Efficient layout and recovery

Store allele strings once per event. Keep observations in compact sparse
vectors, with region/event and molecule indexes; gap windows reference records
rather than copying and repeatedly reconstructing them. Retain the cached
indel/window sequence evidence needed for ambiguous events, not necessarily
whole read sequences in every record.

Cache discovery and read observations with schema, reference/input hashes and
normalization/extractor versions. Keep original block membership in a separate
snapshot and decisions in separate outputs. Changing a stitch threshold should
not invalidate raw evidence. Re-extract only when inputs or extraction rules
change; unsupported sequence comparisons may require new work once, not on
every trial.

This keeps graph traversal and BAM/MSA discovery out of ordinary recovery
scoring. A gap query produces an evidence view, the existing phaser produces a
proposal, an allele-based validator evaluates its bridge, and accepted parity
edges are applied to blocks separately.

## Implementation order and acceptance tests

1. Add the immutable records and lossless adapters alongside the current path.
   Start with SNPs and simple indels; represent complex/unmapped events without
   pretending the old phaser can consume them. Preserve exact provenance.
2. Add a round-trip audit: event/allele identity, molecule counts, source
   observations and qualities must survive merge, sort, cache and projection.
   Matched ordinary events must reproduce the current core's input exactly.
3. Move gap recovery to cached views of these records. Compare source-only and
   combined views without re-reading the BAM or silently replacing evidence.
4. Then evaluate alternate scoring/solver strategies on the same frozen input.
   Preserve the current accuracy fixes as regression controls during migration.

Required fixtures include:

- SNPs, pure insertion/deletion, C -> TG replacement and multi-base substitution;
- repeated-sequence shifts that are equivalent, and nearby events that are not;
- two alternate alleles with no reference support and permuted source allele IDs;
- reference support versus no coverage, deletion versus ref-skip, partial indels;
- GAF/BAM agreement preserving query position and quality, and explicit conflict;
- one molecule observed in two chunks or three evidence sources counting once;
- read/candidate/source ordering, strand, chunk partition and cache round-trip invariance;
- private BAM candidates remain gap-owned even when extraction windows overlap;
- context outside a gap does not become a newly admitted BAM candidate;
- immutable evidence hashes across rephasing and HP-label changes;
- counts matching the chosen observation view and eligibility reasons matching
  the actual graph-node admission decision.

Use all original chr20 gaps plus the known bad 46.84/52.39 Mb bridges and good
57.84/62.41 Mb controls. Evaluate lost correct joins as well as wrong joins.
Keep truth and competitor calls out of event normalization and phasing inputs.
Require a documented untouched evaluation set before describing the new model
as generally calibrated or robust across samples.
