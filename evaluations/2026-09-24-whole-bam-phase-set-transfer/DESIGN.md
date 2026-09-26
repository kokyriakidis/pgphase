# Whole BAM phase-set transfer into graph recovery

Status: design proposal, 2026-09-24. A conservative seam-local subset is now
implemented and measured in `RESULTS.md`. The broader site-level reconciliation
below remains a proposal; `docs/IMPLEMENTATION.md` describes current behavior.

## Goal and evidence

A targeted BAM solve already phases a connected set of sites. Recovery must
carry that *site-to-phase-set relationship* through shared graph/BAM sites,
not just transfer private candidates and per-read alleles. It must preserve
graph-only sites and reads, and it must be able to cut an existing graph phase
set when a local source block starts or the graph phase relation conflicts
with observed molecules.

The concrete chr20 case is:

```text
coordinate       51,235,063   51,262,081   51,270,774   51,286,463   51,287,372
BAM source PS     earlier PS   └─────────────── PS 51,262,081 ────────────────┘
read pairs                    0             30            4             49
current hybrid   graph PS ─────────┘  graph PS | imported PS ───────| graph PS
wanted output    earlier PS  |  └──────────── one connected PS ──────────────┘
```

The BAM source block has 155 oriented sites from 51.262 to 51.414 Mb. Its
read-observation graph is connected. The transferred hybrid contains the
observations at all four shown sites, but its current PS labels separate them.
HiPhase also connects those four sites when given pgphase's own VCF. See
`evaluations/2026-09-24-hiphase-connectivity-audit/README.md` for the matched
source/hybrid matrices and full-chr20 counterexamples.

## Why the present transfer loses the block

`recover_phase_set_seams_in_place` calls the BAM solver over padded seam groups.
It matches BAM candidates to graph candidates by exact sequence-normalized
identity. For a match, it records a shared-site orientation vote and refreshes
read observations, then keeps the graph row and its graph PS. It appends a BAM
candidate only when it is absent from the graph and lies strictly inside a
seam. The source PS is remapped only for those appended candidates and eligible
previously unphased reads. Thus a source PS can span several graph PSs while
only scattered private rows carry its membership.

The read-profile rebuild also copies only the primary allele and query-index
arrays. At a shared row, BAM's observation replaces the graph observation in
that array; both channels and their disagreement are not retained there. A
later whole-chunk BAM observation pass can populate `bam_alleles`, but it does
not reconstruct the missing source PS membership or the original graph
observation at this recovery boundary.

## Representation

Keep one output candidate for each *exactly equivalent* allele. Retain every
BAM-private candidate row and its allele representation as produced by the
BAM solve, including complementary indel rows. The BAM result also supplies a
recovery-local overlay, indexed by the final candidate table:

```cpp
struct SourceSiteMembership {
    size_t candidate_index;          // final graph/BAM candidate index
    size_t solve_id;                  // identifies the targeted BAM solve
    hts_pos_t source_phase_set;       // local to solve_id
    int source_hap1_allele;
    int source_hap2_allele;
    // Exact allele-index translation if this is a shared graph row.
};
```

The exact field layout can follow existing containers, but the information
must remain separate from `CandidateVariant::phase_set`: one shared site can
legitimately have a graph PS and a BAM source PS before reconciliation.
Source PS identity is `(solve_id, source_phase_set)`, never the numeric PS
alone. Keep the BAM read's source HP/PS in the same overlay. Build the overlay
for **every phased source site in a selected BAM PS within the active graph
chunk**, including the padded flanks. Select a BAM PS if any of its sites
intersects a recovery seam; do not truncate its membership to the strict gap.
When a source PS crosses a chunk boundary, represent its in-chunk part and
let the existing cross-chunk logic decide the later connection.

At a shared site, translate the graph allele walk through `GraphSiteMeta` to
reference allele sequence and match an exact normalized `(POS, REF, ALT)`
identity and allele-index mapping. One output row retains both provenances.
Do not infer equivalence for overlapping but different indel descriptions.
Keep those rows separate and mark their locus as a representation conflict so
they cannot become two independent votes for the same molecular event.
Multi-allelic graph alleles are mapped individually; an unmapped allele is
unknown rather than silently converted to a biallelic REF/ALT observation.

Preserve BAM and GAF observations in parallel arrays through the candidate
re-index. A qname is one molecule: agreement between its BAM and GAF calls at
one site supplies one observation; disagreement is recorded as a conflict and
cannot vote for a join. Keep source MAPQ, graph MAPQ, BAM base quality and site
category available for evidence selection. New BAM-only reads need their BAM
quality metadata, not just an inferred coordinate range.

## Reconciliation algorithm

1. **Build site connectivity.** Use phased heterozygous loci in coordinate
   order. For each unique molecule, connect successive loci at which it has
   valid, usable allele observations. This is a sparse incidence graph, not an
   all-pairs graph. A source PS is an orientation proposal, never a synthetic
   read edge. A zero-read cut, such as 51,235,063 to 51,262,081, separates
   components even if the old graph PS crosses it. Co-located alternative
   rows are one locus for cut placement.

2. **Use the complete BAM block as the local scaffold.** Within each
   read-connected source PS component, copy the source site's diploid allele
   orientation through the exact allele mapping. This preserves the BAM
   solver's already computed 155-site chain without a second whole-chunk
   k-means run or a 20-variable exact-search bound. Source sites keep their
   relative orientation; the common HP gauge can be flipped freely.

3. **Segment graph phase sets at site level.** Compare each overlapping graph
   site's orientation with the source orientation. Split an old graph PS at a
   zero-read cut, a source PS boundary, or a change in its relative parity to
   the source along independently observed shared sites. Validate each
   resulting segment with the actual molecules that cover its sites. A local
   boundary vote can orient that segment; it cannot flip an unrelated distant
   part of the original graph PS. Ambiguous or conflicting segments remain
   independent. This is the guard missing from the rejected 19 Mb whole-block
   join, which made 279 previously correct reads discordant.

4. **Attach graph-only information.** Preserve graph-only candidate rows and
   their GAF-only molecules. A graph-only segment joins a source component
   when its allele observations give a stable orientation against that
   component; otherwise it remains a separate PS. Unphased graph-only sites
   can use the existing local phaser with the established source sites held as
   anchors. The source BAM consensus itself is not rederived from graph read
   labels.

5. **Emit final labels and reads.** Assign one PS to each accepted, connected
   component. Orient a graph segment by its accepted local allele mapping,
   updating only candidates and reads belonging to that segment. Re-score
   affected reads once from their final sites, with each molecule and each
   locus counted once; leave unaffected graph assignments alone. BAM-only
   reads retain the source block's HP after its gauge is mapped. Output one
   VCF description per exact shared allele and preserve each distinct BAM
   indel row.

For a proposed graph-to-source attachment, compare the read likelihood under
both relative orientations using unique molecules and their allele quality.
Clean SNPs supply the strongest evidence; verified indels can contribute at a
calibrated lower reliability. Shared-site consensus establishes an allele
mapping but is not counted as an extra read. Require a stable winning
orientation in deterministic disjoint molecule subsets. If the data do not
separate the orientations, keep the blocks independent. This rule uses a
likelihood/error target rather than a fixed net-read count; its calibration
must be measured on held-out joins and checked against parental truth, with
no truth used by the runtime algorithm.

The source PS need not pass a second stitch test at every adjacent source
site: it is the BAM pipeline's existing phasing result. Its raw read links are
still checked for physical connectivity and conflicting representation. If
there is a true zero-read cut or a clear parity conflict, split the imported
source component. This permits the 4-read middle link while refusing the
0-read cut to the left.

## Invariants and cost

- Every source site in a selected source PS maps to exactly one final
  candidate or remains an explicit unmatched source-only row; no silent drop.
- Existing graph and BAM observations survive with provenance. A read name
  contributes at most one molecular vote per locus.
- A final PS has a path of real allele observations between its sites.
  Neither a pre-existing graph PS label nor a BAM PS number creates a link.
- An accepted graph attachment changes only its validated site segment.
  A split never flips or relabels the distant remainder of a graph block.
- Candidate, read-profile and graph metadata arrays remain index-parallel;
  output records and read PS labels describe the same final component.
- The normal path is one candidate sort plus linear scans over candidates
  and sparse read observations. No chromosome-wide rerun or combinatorial
  search bound controls whether an existing BAM source PS can transfer.

## Implementation and verification sequence

1. Add the source-site/read overlay to the transfer result and retain all
   selected source PS members. Unit-test exact shared matching, multi-allelic
   allele maps, private-site completeness, distinct complementary deletion
   rows, source-ID collisions and chunk clipping.
2. Preserve both observation channels through re-indexing and add an audit
   assertion that the source PS's site/read incidence survives transfer.
   Test one read present in both BAM and GAF, including disagreeing alleles.
3. Add read-connected site components and graph-PS segmentation. Make a
   matrix-only replay test, without BAM/GAF parsing, for the 0/30/4/49 chain.
   Assert 51,235,063 remains separate and 51,262,081 through 51,287,372
   share one PS with the source genotype orientation.
4. Add a counterexample where the graph PS switches internally and only its
   local end agrees with the source. Assert the distant 279-read segment is
   not absorbed. Add exact and alternative indel representation cases.
5. Run the existing chr20 window panel and matched full chr20, scoring both
   read truth and site switches per final PS. The retained baseline is
   236,520 phased / 228,184 correct / 8,336 discordant reads and 456,233-bp
   VCF N50. Require the target chain to close without reopening the 51.235-
   51.262 Mb false connection; inspect every newly joined block for internal
   polarity changes. Then test the same logic on the available chr12/chr18
   fixtures and profile runtime/memory before replacing the old stitch path.

`docs/IMPLEMENTATION.md` should be changed only when this behavior is actually
implemented and validated.
