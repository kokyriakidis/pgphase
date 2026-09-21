# Planning graph recovery with private BAM variants

## Current path

`collect-graph-variation` loads catalog sites and GAF allele observations into
`GraphChunkBuildResult` (`build_graph_chunk`), applies the graph noise filter,
and phases clean catalog candidates. With `--bam`, each worker then calls
`run_in_chunk_recovery` before cross-chunk stitching. That calls
`retry_unphased_windows_in_place`, which finds unphased-read intervals and
phase-block seams, pads them to nearby phased parent sites, and runs the same
BAM `process_chunk` used by `collect-bam-variation` over those windows at the
recovery MAPQ floor.

The merge matches catalog alleles to BAM sequence alleles through
`GraphSiteMeta` and `vcf_to_variant_key`, imports missing candidates, transfers
BAM observations by read name, builds VCF metadata for imported candidates,
then re-runs clean and noisy-inclusive k-means over the full parent chunk.
`--recovery-audit-out` reports whether each BAM candidate was known, inside a
window, appended, and given usable metadata. `--graph-noisy-msa` additionally
routes graph-demoted repeat loci through the BAM sub-solve; it is off by default.

## A bounded chr20 check

Input: `CHM13#0#chr20:5,259,406-5,395,085`, containing the committed gap
`5,309,406-5,345,085`. The panel records a HiPhase spanning block at 97.7%
read separation; the default graph arm does not span it.

In the default recovery audit, 148 BAM candidates are strictly inside the
recovery window. One matches a catalog allele after sequence translation;
147 are absent from the graph candidate set and all 147 are appended. The
unmatched candidates comprise 127 `CleanHom`, 18 `NoisyCandHom`, and two
`NoisyCandHet`; only two are marked alignment verified. The BAM-derived deletion
at 5,331,265 is imported and emitted, but the graph VCF still has distinct
phase sets at the left and right gap endpoints (5,272,413 and 5,331,265).
The standalone BAM solve also splits there. Site import alone cannot give this
gap a phase relationship that its reads and solver have not established.

On this same window, `--stitch-recovered` leaves the gap split.
`--graph-noisy-msa` produces one phase set spanning both endpoints. Against
`test_data/derived/chr20_truth_hap.tsv`, the default graph arm places 495 of
497 scored tagged reads concordantly; the flagged arm places 516 of 518.
Matched by read name and oriented per phase set, the flagged arm gains 21
concordant tags, loses none, and changes zero concordant reads to discordant.
This is a local result. The earlier whole-chr20 experiment in
`evaluations/2026-09-19-graph-noisy-msa/` reported read Hamming error rising
from 1.161% to 4.197% when repeat-locus recovery was broadly admitted.
The flag cannot simply become the default.

## A direct private-site proof case

Screening the 11 gaps in `evaluations/2026-09-14-graph-gap-bam/targets.tsv`
with the current parity-tested BAM solver found one where a single BAM phase
set spans the endpoints: `chr20:19,358,995-19,395,544`. The graph run over
`19,308,995-19,445,544` remains split. Its recovery audit shows eight
BAM candidates unmatched to active graph candidates strictly inside the gap,
including five MSA-verified
`NoisyCandHet` candidates. All eight are appended with usable metadata.
Four of the verified loci appear as phased heterozygous VCF positions inside
the gap (two have complementary allele rows), in the left graph phase set.
The right catalog SNP at 19,395,544 matches the BAM call after allele
translation, yet stays in a second graph phase set. The BAM VCF puts that
same SNP and the unmatched BAM sites in one set spanning 19,310,197-19,444,196.

These are not yet proven absent from the RAW catalog: the catalog VCF has
indels within 1-11 bp of several BAM calls in these repeat tracts, and
`vcf_to_variant_key` trims shared flanks without left-aligning equivalent
indels across repeats. Position or minimal-key inequality is therefore only
"unmatched to the active candidate table," not proof of a private allele.
The first implementation step must compare local haplotype sequence before
injecting a second description of such a locus.

This is the first implementation target for the requested design. Discovery
and import reached the unmatched sites; the unresolved question
is whether BAM observations at the right catalog SNP and the unmatched BAM sites
form the same read-allele chain in the graph chunk just before k-means. The
current candidate audit cannot answer that, because it does not export the
post-merge observation matrix or shared-read counts at the boundary. That
matrix should be checked before changing the solver.

The current `retry_unphased_windows_in_place` returns immediately when a graph
chunk has no candidates or no GAF reads. Such a chunk cannot currently trigger
BAM discovery through the usual window detector, even if the BAM covers it.
This is a separate targeting blind spot; it does not explain the measured 5.3 Mb
gap, which has graph candidates and reads.

## Intended design: inject, then use the ordinary phaser

The graph solve made before recovery is only a way to locate missing-evidence
windows. It must not become a second authority on the recovered phase. In each
selected window, use the parity-tested BAM `process_chunk` to discover candidates
and genotype reads. Match the BAM candidates to the catalog by normalized
sequence allele, add candidates absent from the catalog, and add their
per-read allele observations to the graph profiles. Then clear provisional
HP/PS and run the normal clean and noisy-inclusive phasing rounds on that union.
No private-site phase set or sub-solve consensus is pinned, and no special
block-stitch rule decides the result. The ordinary solver decides naturally
from the augmented site/read matrix, as in the BAM arm.

This is close to the code already present: `retry_unphased_windows_in_place`
imports BAM sites and observations, and `run_in_chunk_recovery` runs both shared
k-means rounds again. The remaining work is to make the selection and transfer
complete, explain each rejected private call, and test whether the union has a
real read-allele chain across a target gap. Broadly turning on
`--graph-noisy-msa` is a different policy: it also revisits demoted catalog
repeat loci, and its whole-chr20 result regressed badly.

The 5.3 Mb gap is a counterexample to claiming that private-site injection by
itself closes every competitor gap. The current code already injects 147
candidates unmatched to active graph candidates there and leaves two phase sets. `--graph-noisy-msa` closes
it because it additionally recovers a repeat-demoted catalog locus. If the
scope is strictly candidates absent from the catalog, that window is an audit
control, not a promised win. The five unspanned graph
windows in the committed six-window panel checked with the standalone BAM
solver also remain split; the 19.36 Mb proof case comes from the separate
11-gap evaluation panel.

## Implementation sequence

1. **Start with the 19.36 Mb proof case.** First classify the nearby
   catalog repeat alleles by local haplotype sequence, so only truly absent
   BAM alleles are called private. Export the exact post-merge read-by-site
   matrix for the five verified unmatched loci and the right
   catalog SNP, then compare its co-observations to the BAM sub-solve. If the
   right-SNP observations are absent, fix the sequence-key/read-name transfer;
   if they are present, trace where the ordinary k-means separates that SNP.
   For each other graph gap or seam,
   record the normalized BAM candidates absent from the catalog, their original
   BAM categories, and how many reads co-observe each private site and a phased
   graph flank. Separate gaps where a private candidate supplies a continuous
   read-allele chain from gaps where BAM itself splits or no read bridges a cut.
   The existing `--recovery-audit-out` provides most candidate decisions; add
   the missing chain and final-PS fields before changing admission.
2. **Complete the candidate/read union.** Reuse the BAM pipeline's candidate
   discovery and allele profiles without transcribing another caller or
   inventing genotypes. Match against catalog `GraphSiteMeta` after canonical
   VCF normalization and local repeat-equivalence checks, insert only candidates
   with no matching catalog allele,
   transfer BAM observations by read name, and preserve catalog/GAF evidence
   at shared sites. Keep the index-parallel metadata arrays synchronized and
   retain actual BAM alignment spans for reads added to the graph chunk.
   Cover the current early return when a graph chunk has no candidates or GAF
   reads but BAM reads exist, so private-only sequence is not invisible.
3. **Run the existing phasing stages once on the union.** Use the preliminary
   graph solve only to identify target windows. After injection, discard its
   provisional read labels and call the shared clean and noisy-inclusive
   k-means stages in their normal order. Let the normal category masks decide
   which imported BAM candidates participate. Do not copy the BAM sub-solve's
   phase-set IDs or hand-orient the private variants. Keep a gap open if the
   augmented read/site matrix lacks a supported chain.
4. **Verify on real gaps before expanding scope.** Add unit tests for exact
   sequence matching, no duplicate catalog alleles, read observations,
   metadata order, and catalog-empty chunks. Run the committed gap panel and
   the 11 unresolved targets from
   `evaluations/2026-09-14-graph-gap-bam/`, with parental-truth scoring of
   spans, in-gap phased hets, read concordance, and concordant-to-discordant
   transitions. Then run whole chr20, comparing read Hamming error, blocks,
   output records, and runtime to the shipped graph arm. The graph-only arm
   and parity-tested BAM arm must remain unchanged.

If private BAM calls are present but the ordinary union solve still splits,
report that as a linkage or solver problem. The 5.3 Mb control demonstrates why
copying the calls again or forcing their phase labels would not solve it.

## 2026-09-20 no-merge audit

The requested design requires the BAM channel to inject the caller's separate
site rows exactly. The current recovery sub-solve violates that requirement because it inherits the
graph default `merge_colocated_msa_alleles=true`. At 11.23 Mb, standalone BAM
emits two complementary deletion rows at 11,255,369; the recovery sub-solve
emits one three-allele candidate. Both descriptions reach the graph merge
unchanged, so the defect occurs in the sub-solve options, not in reindexing.

Setting the BAM port's exact options in recovery restores separate rows and
per-read MSA calls, but it regresses the committed graph windows. Disabling
only the merge loses phased heterozygotes because the graph-specific MSA refresh
then changes the split rows' genotypes. Neither trial was retained. This is
why the next implementation must separate exact BAM discovery/profile transfer
from graph catalog evidence without turning a failed local test into a shipped
regression. No experimental code from this audit is in the current behavior.

A second bug affects the meaning of “private.” In the 3.85 Mb window, 32/39
BAM heterozygotes audited as unmatched are local-haplotype equivalents of raw
catalog ALTs with different anchors. Checking only `CandKey` or VCF position
would inject duplicates. The same issue appears in the 11.23 and 24.10 Mb
windows. Matching against all raw catalog alleles removed duplicates in a
trial, but also removed most in-gap phasing evidence and failed the committed
regression suite. The record-level distinction between absent from raw catalog
and absent from retained phasing candidates needs to be resolved before the
private-site policy can be finalized.

The chromosome-wide empty-catalog fallback was tried and reverted. It added
2,271 calls across the low-MAPQ 27.5–29 Mb region and increased read discordance
from 1.154% to 1.633%. The noncentromeric five-window screen and matrix counts
are recorded in `evaluations/2026-09-20-graph-recovery-windows/`.
