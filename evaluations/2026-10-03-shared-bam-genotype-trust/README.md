# Preserve BAM alleles without assuming their whole-block orientation

## Question

Can targeted recovery simply trust verified BAM sites inside an unphased
graph gap? A gap can reflect missing sites, but can also reflect missing
allele observations or an unsupported phase connection. MSA verification
establishes an allele representation; it does not independently certify all
site orientations in a source phase set.

## Observed transfer loss

In the owning 21–22 Mb replay, the MSA-verified BAM deletion at key position
21,831,481 has genotype `0|1` in source PS 21,817,199 and 11 REF / four ALT
observations. Sequence projection matches an existing, unphased catalog repeat
row. Transfer marks it alignment verified and retains BAM-channel calls, but
the graph row keeps its repeat demotion and unset genotype: the current
genotype ownership rule applies to suffix-padded matches only.

This is evidence to improve shared-row transfer, not proof of a safe stitch.
The source deletion at 21,823,067 shares that PS with 21,831,481, yet no
retained source read has callable alleles at both rows. The latter deletion
and the right insertion at 21,844,360 have 15 paired observations, with raw
allele combinations `(0,0)=10`, `(1,0)=3`, `(0,1)=1`, `(1,1)=1`.
Graph observations at the shared middle row are also ambiguous. Merely
inheriting the source PS would conceal an unsupported edge.

These observations are from native matrix dumps under
`test_data/tmp_gap_fix44/probe/21/`. Matrix `VAR` column five is the phasing
**weight**, not the number of distinct alleles; the header omits several
actual columns. Allele-count eligibility must be checked against the code,
not inferred from this weight.

## Rejected chromosome experiment

Broaden shared-row genotype adoption beyond suffix-padded matches when the
complete source block has at least two phased heterozygotes and **no weak
read-supported path cut**. Retain the other binary, verified MSA,
unphased-repeat and independent-source conflict checks. Adoption copies the
BAM key, genotype, depths and remapped PS together and uses the source
observations rather than graph calls. This was a constrained trial, not a
test that trusts every BAM label unconditionally.

| Full chr20 metric | Accepted baseline | Rejected trial |
|---|---:|---:|
| Truth-scored phased reads | 237,225 | 237,178 |
| Correct assignments | 230,105 | 230,063 |
| Discordant assignments | 7,120 | 7,115 |
| VCF keys | 64,185 | 64,485 |
| VCF phase blocks | 330 | 332 |
| Span N50, bp | 806,449 | 790,093 |
| Required coordinate connections | 89/102 | 89/102 |

Aggregate accuracy improves slightly, but 26 previously correct assignments
become incorrect and 62 previously correct reads become unphased. Only 19
previously unphased reads become correct. No old key disappears, but two old
SNP blocks acquire split gauges. The target at 21.823 Mb remains open.
Therefore the trial is rejected; no production behavior from it is retained.

`chromosome-parity.json` records all read-state transitions;
`block-transfer.json` records old SNP gauges and every panel span.
Truth is evaluation only and is never read by the phasing implementation.

Starting executable SHA256:
`336996e6338b3abc4a55b662231e71e86d378f4222615c4e39aa6cc94a645ce4`.
Trial executable SHA256:
`32d529f27a5b3ba1f2b208dd3febe27a25de400bd62602f17e94debf7a897e3d`.
The trial used the same annotated BAM, reference, catalog and GAF as the
accepted baseline, with eight threads over the whole chromosome; elapsed
time was 242.28 seconds, not a controlled competitor runtime comparison.
Outputs: `test_data/tmp_gap_fix44/full-path-adoption/`.

Preserving individual verified BAM sites remains appropriate. A safe next
repair must preserve their independently supported component orientations
and the existing clean SNP gauges, rather than treating a source PS label
as a certificate for every internal edge or for a graph-block join.
