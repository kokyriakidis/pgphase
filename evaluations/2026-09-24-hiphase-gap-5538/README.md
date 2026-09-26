# Why HiPhase closes chr20:55.381 Mb and graph recovery does not

Date: 2026-09-24. This is a noncentromeric, current full-chr20 miss. The
comparison uses the same annotated HG002 BAM, CHM13 reference and pgphase
variant rows. Frozen HiPhase also spans the interval using its DeepVariant
calls. The fresh HiPhase replay below uses `--ignore-read-groups` because
pgphase's VCF sample is `SAMPLE` while the BAM read group names HG002.

## Current output and actual molecules

The current full-chr20 pgphase VCF has a left PS 55,360,776 at the repeat
sites 55,360,776/55,360,780 and a right PS 55,373,606 containing the clean
SNP at 55,382,729. Frozen HiPhase puts the left sites and 55,382,729 into
PS 55,360,776. Fresh HiPhase run on pgphase's *own local VCF* does the same;
it leaves pgphase's 55,373,606 G>T row heterozygous but unphased.

Two primary reads physically span 55,360,780 to 55,382,729; both are maternal
and HiPhase tags both in PS 55,360,776. One of them,
`m84031_231217_034919_s2/120324954/ccs`, has a valid C at the left SNP and
A at the right SNP, both BQ40 and MAPQ60. The other has no valid right-SNP
base because its alignment deletes that coordinate. HiPhase's block builder
accepts one spanning mapping by default, and `get_phase_block_ids` traverses
valid allele observations with one connecting read by default. The fresh
HiPhase replay demonstrates that this site path is enough for its joint solve.

The middle indel pair has much stronger physical coverage: 57 primary reads
span 55,381,221 to 55,381,472 in frozen HiPhase output, and all 57 carry
PS 55,360,776 with truth-consistent HP labels (32 maternal HP1 and 25
paternal HP2). In pgphase's targeted BAM source matrix, corresponding
55,381,222 deletion and 55,381,473 insertion observations have 44
same-gauge allele pairs versus six conflicting pairs, with one uncalled.
The source BAM solver nevertheless assigns the two sites to distinct local
phase sets, 55,360,776 and 55,381,472.

## First loss: a false graph heterozygote shortens the seam

Pgphase's graph row at 55,373,606 is called a clean G>T heterozygote from
13 G versus 50 T graph observations and becomes the first site of the right
PS. The matched primary BAM has 49/49 T bases at that coordinate, 48 at
MAPQ >=30; the targeted BAM caller labels it homozygous T (20 T, zero G in
its local candidate). The ten graph-reference observations still present in
the hybrid matrix belong to reads with no BAM base at that coordinate and
have no useful parental separation (seven maternal, three paternal).
HiPhase's DeepVariant VCF has no heterozygous call there; when supplied with
pgphase's own VCF, HiPhase leaves that row unphased.

Recovery detects a seam from graph PS 55,331,033 to 55,373,606. Its strict
candidate transfer ends at 55,373,606, so the informative BAM indels at
55,381,222 and 55,381,473 lie outside the seam and are not injected. This is
a concrete false-anchor problem, not a shortage of spanning molecules.

## Second loss: an outer transaction vetoes an independent inner block

A controlled local graph catalog containing the same rows except
55,373,606 moves the seam's right flank to 55,382,729. Pgphase then injects
the two BAM indels, but leaves them in PS 55,360,776 and PS 55,381,472:
the supported internal link is still open. The left *outer graph* PS
55,331,033 has a nearly tied whole-window BAM gauge (26 same, 21 cross).
The current recovery stitch requires a complete graph-to-graph transaction;
when it cannot prove that outer attachment, it rolls back internal BAM/BAM
joins too. The later fallback explicitly leaves two imported BAM blocks
separate. HiPhase starts a new read-connected block at 55,360,776 and does
not need to attach the unrelated 55,331,033 flank.

The local pgphase A/B is also informative: with the false row present it
phases 1,243 truth-evaluable reads, 1,196 correct and 47 discordant. Removing
only that row gives 1,222 phased, 1,195 correct and 27 discordant. Thus
masking improves local purity but loses 21 phased reads and still does not
close the target; it is not a complete production fix. Fresh HiPhase on
both local pgphase VCF arms connects 55,360,776 to 55,382,729 and leaves
the disputed G>T row unphased in the unmasked arm.

## Implication

A fix needs two operations with separate evidence: reconcile graph
heterozygotes with exact BAM homozygous calls before using them as seam
anchors, then build the maximally supported *independent* phase component
inside a seam from read-connected BAM and graph sites. Failure to attach an
outer graph flank must not erase a locally valid component. The prior
permissive internal-BAM replay raised full-chr20 discordance (recorded in
`CHECKPOINT.md`), so the 44:6 pair alone is not a safe general join rule.
The two long reads, local site path, representation agreement and parental
truth make this interval a useful positive control for a guarded redesign.

Local pgphase diagnostics and the temporary one-site catalog A/B are under
`/tmp/pgphase-miss-5538/` and `/tmp/pgphase-miss-5538-narrow/`. The matched
full-chr20 baseline is `/tmp/pgphase-wholeps-final2-chr20/`; frozen HiPhase
is under `~/Downloads/pgphase-eval-data/results/chr12-18-20-comparison/chr20/hiphase/`.

## Guarded implementation and whole-chromosome A/B

The implementation first compares a graph seam's right clean-SNP anchor with
an exact matching SNP in the targeted BAM solve. A clean or MSA-verified BAM
homozygote can remove the anchor's phase label only when MAPQ >= 30 molecules
from both BAM haplotypes support the same allele at two-sided binomial
`p <= 0.01`, and the expanded seam remains within the BAM region already
solved. The graph row retains its original heterozygous genotype unphased.
At 55,373,606 this admits the previously omitted 55,381,221 deletion and
55,381,472 insertion into recovery.

After the outer graph-to-graph stitch transaction fails, a narrow independent
inner-component path can join the last two imported BAM blocks to the right
graph block. It requires a significant allele vote between the two BAM blocks,
the second BAM block already attached to the right graph block, and at least
one conflict-free direct MAPQ >= 30 BAM molecule with clean or MSA-verified
heterozygous SNP calls on the first BAM block and the right graph block.
The direct check uses the BAM mapping quality stored with its allele channel,
not the graph alignment's mapping quality. The unrelated left graph block
stays separate.

A first version of the anchor check counted low-MAPQ BAM observations and
incorrectly unphased the true 30,794,399 A>C graph SNP. Its primary BAM bases
are 31 C and nine A, mostly below MAPQ 30; parental truth supports both
alleles. The MAPQ guard restores its original phased call. At 38,177,812,
by comparison, all 37 primary BAM bases are A (20 maternal, 17 paternal),
so leaving the graph heterozygote unphased is consistent with the alignment
and truth evidence.

Matched full chr20 with the same graph catalog, GAF, BAM and command options:

| Metric | Previous | Guarded implementation |
|---|---:|---:|
| Truth-evaluable phased reads | 236,526 | 236,527 |
| Truth-correct phased reads | 228,189 | 228,190 |
| Discordant reads | 8,337 | 8,337 |
| VCF phase blocks | 464 | 463 |
| VCF block N50 | 456,233 bp | 456,233 bp |

The 55,360,776 left BAM block, both middle indels and 55,382,729 right SNP
now share PS 55,360,776. The false graph SNP at 55,373,606 remains
heterozygous but unphased. The separate 55,331,033 graph PS stays separate;
the known 19 Mb wrong-join control also remains split. The integration
regression checks this site path, the truth-scored absence of a local switch,
and preservation of the true low-MAPQ graph SNP. The final run is under
`/tmp/pgphase-try5538-final/`, compared with
`/tmp/pgphase-wholeps-final2-chr20/`.
