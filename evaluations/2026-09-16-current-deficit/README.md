# The deficit, measured against the current pipeline

The earlier deficit set (187 of 311 gaps, 4.74 Mb) was derived from the
**graph-only pass-1** gap inventory, which does not describe what the pipeline
leaves open. `chr20:30,794,962-30,814,005` makes the point: all three competitors
close it, it was on the deficit list, and the hybrid path already covers it with a
single 76.1 kb phase set of 305 sites at 99.6% read accuracy -- gap recovery is
never even asked about it, and 116 clean het SNPs sit inside.

`run.sh` runs the current full pipeline on chr20 (graph sites + BAM channel,
`-q 1`, recovery on; 10.5 min on 16 threads), `find_current_deficit.py` derives
gaps from its phased VCF and intersects them with competitor spans, and
`score_deficit.py` scores each competitor span against the diplinator read truth.

## What is actually open

| | gaps | span |
|---|---:|---:|
| our phase blocks | 238 | 56.02 Mb phased |
| gaps between consecutive blocks | 196 | 11.88 Mb |
| spanned by at least one competitor | 47 | 1.20 Mb |
| spanned by nobody | 149 | 10.67 Mb |

So the deficit is 1.20 Mb, not 4.74 Mb, and nine tenths of the remaining gap span
is closed by no tool at all.

## A span is only an opportunity if it is right

| tool | spans (>=20 scored reads) | read concordance | >= 98% accurate | < 90% accurate |
|---|---:|---:|---|---|
| hiphase | 45 gaps | 91.11% | 31 gaps, 0.80 Mb | 8 gaps, 0.25 Mb |
| whatshap_opt | 27 gaps | 92.78% | 17 gaps, 0.43 Mb | 5 gaps, 0.11 Mb |
| longphase | 7 gaps | 97.26% | 5 gaps, 0.16 Mb | 1 gap |

Several spans are worse than chance -- hiphase reads 31.8% at
`36,059,715-36,099,955` and 35.9% at `58,994,554-59,046,421`, whatshap 43.9% at
`60,033,052-60,058,235` -- so in those gaps our abstention is better than their
join. **The genuinely recoverable deficit is 31 gaps, 0.81 Mb**, where some tool
spans at >= 98% over >= 50 scored reads.

## In those gaps, what do we hold at their sites?

`site_deficit.py` takes every phased het record hiphase places inside the 30 gaps
it spans at >= 98% (0.79 Mb) and asks what our own candidate table says at that
position. `genotype_sites.py` then genotypes each site from the alignment and
scores its allele partition against read truth -- substitutions from the base at
the position, indels from the net insert-minus-delete length over a window, with
the classes we already trust as the control.

| their site | our category | n | median segregation | informative | phantom |
|---|---|---:|---:|---:|---:|
| SNP | `CLEAN_HET_SNP` | 44 | 1.000 | 100% | 0% |
| DEL/INS | `CLEAN_HET_INDEL` | 6 | 1.000 | 100% | 0% |
| DEL/INS | `NOISY_CAND_HET` | 33 | 0.904 | 58% | 18% |
| DEL | `NOISY_CAND_HOM` | 7 | 0.939 | 57% | 29% |
| DEL/INS | absent from our table | 5 | 0.953 | 80% | 0% |
| DEL/INS | `REP_HET_INDEL` | 3 | 0.800 | 33% | 0% |

Three results follow.

**We already have every SNP they use** -- all 44, as clean het SNPs. SNP discovery
is not the deficit.

**The withheld class is majority-informative here.** `NOISY_CAND_HET` indels
inside these gaps are 58% informative against the 23% measured for
repeat-demoted indels chromosome-wide, so the global figure understates them
badly in exactly the windows that matter.

**Seven of their het sites we call homozygous.** `NOISY_CAND_HOM` at a site that
segregates at a median 0.939 is a genotyping error, not a screening decision: a
site called homozygous can never link.

The controls are what make this readable. Clean het SNPs and clean het indels both
come out at median 1.000 and 100% informative; an earlier run of the same script
scored the SNP control at 0.542 with 0% informative, which was a one-base
coordinate error in the genotyper (SAM POS and VCF POS are both 1-based and the
reference walk starts at SAM POS), not a property of the sites.

## MSA verification is not the blocker

On `chr20:48,176,830-48,229,446` (52.6 kb, hiphase spans at 100.0% over 252
reads), hiphase crosses on 8 het records, 7 of them indels in homopolymer or
tandem context (`AAC>A`, `T>TA`, `GAGAAAGA>G`, `AAAAAAA>A`). The audit shows our
pipeline **does** verify them: every `msa_verified = 1` site there segregates at
0.89-1.00 against truth. The evidence is present, correctly discovered and
correctly verified.

What blocks the join is structural. Every tier links both flanks but to different
proposal phase sets, so no single block holds both, and recovery requires one that
does. The proposal fragments at the holes between site clusters, and the holes are
not equivalent:

| hole | width | reads covering both flanking sites |
|---|---:|---:|
| 48,096,582 -> 48,123,657 | 27.1 kb | **0** |
| 48,123,657 -> 48,147,230 | 23.6 kb | **0** |
| 48,204,384 -> 48,225,788 | 21.4 kb | 7 |
| 48,183,977 -> 48,202,057 | 18.1 kb | 7 |
| 48,162,480 -> 48,176,831 | 14.4 kb | 15 |
| 48,149,549 -> 48,162,480 | 12.9 kb | 19 |

Coverage at the two hard breaks is 72.5x and 66.9x with longest reads of 29.8 and
29.5 kb, so they are linkage breaks, not coverage holes: no read covers both
flanking het sites, and no read-based phaser crosses them.

**The two hard breaks are outside the deficit interval.** Recovery's gap for this
region is 134.7 kb and swallows them, while the interval a competitor actually
spans is the right 52.6 kb, whose widest holes carry 7 spanning reads each.
Because recovery is all-or-nothing over a whole gap, the bridgeable part is
abandoned together with the unbridgeable part. That is the method limit these gaps
share, and it is separate from site admission.

## Gaps that are correct abstentions

Two gaps were diagnosed and should not be pursued.

`chr20:35,919,404-36,156,319` (236.9 kb): 58 interior sites cover it edge to edge,
but the spacing `35,959,001 -> 35,981,762` (22.8 kb) has zero reads covering both
sites where every other large spacing has 8-20, at 65-71x coverage. Its interior
blocks were being discarded, which is fixed (commit a302c67), but the gap itself
cannot be joined.

`chr20:13,429,829-13,631,825` (202.0 kb): no competitor spans it either --
hiphase places 11 phased hets inside in 4 phase sets, longphase 7, whatshap 10 --
and it holds 8 clean het SNPs in 202 kb at 75x depth across all MAPQ. Its
proposal blocks cannot be chained: `chain_proposal_blocks.py` finds adjacent pairs
sharing 0-12 reads that observe consensus sites in both, and the one pair with 12
votes splits exactly 6/6. A heterozygosity-poor region where abstaining is right
for everyone.
