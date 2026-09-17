# The alignment pipeline's two-stage solve, applied in the hybrid

The alignment pipeline solves in two rounds and the machinery is already there:

| stage | call | mask |
|---|---|---|
| 1 | `collect_var.cpp:2105` `assign_hap_based_on_germline_het_vars_kmeans(chunk, opts, kCandGermlineClean)` | clean SNPs, indels, clean hom |
| 2 | `collect_var.cpp:2065` same k-means with `kCandGermlineVarCate` | `kCandGermlineClean \| kCandNoisyCandHet \| kCandNoisyCandHom` |

But stage 2 as it stands is **read recovery, not refinement**.
`gap_fill_unphased_reads` runs it, harvests a haplotype only for reads the core
left at `hap == 0` into a disjoint phase-set namespace, and then restores every
candidate field -- `phase_set`, `hap_alt`, `hap_ref`, `hap_to_cons_alle` -- so no
noisy site ever gains a phase set and none can extend or link a block.

`--msa-verified-refine` runs the same second round and **adopts** it, for reads
and candidates, over the clean sites plus the **MSA-verified** noisy ones (an
unverified noisy candidate is hidden from this k-means alone and keeps its
category for every other consumer).

## It does what it should, per site

On `chr20:24,105,188-24,142,287` the hybrid now reaches the alignment channel's
own state, which it never did before:

| site | hybrid before | hybrid after | alignment channel |
|---|---|---|---|
| `24,121,714` (multiallelic, informative) | `PS = 0` | **`PS = 24,103,779`** | `PS = 24,103,779` |
| `24,131,708` (phantom) | `PS = 0` | `PS = 24,103,779`, no haplotype | same |

and the blocks grow with it: the left block goes from `24,103,779-24,105,188`
(2 sites) to `24,103,779-24,121,713` (3 sites), the right from 64 to 77 sites.

## And it costs more than it buys

### First measurement, before the bridge bug was found

| arm | spanned | in-gap hets | concordance | c -> d |
|---|---:|---:|---:|---:|
| stock defaults | 0 of 6 | 12 | 99.68% | 0 |
| `--retry-unphased-with-bam` (shipped) | 4 of 6 | 31 | **99.49%** | **0** |
| `--msa-verified-refine` | 4 of 6 | 25 | 92.19% | **152** |
| `--msa-verified-refine --retry-unphased-with-bam` | **5 of 6** | 30 | 97.46% | **64** |
| ... `--min-block-link-reads 3` | 3 of 6 | 30 | 97.18% | 65 |
| ... `--min-block-link-reads 5` | 2 of 6 | 30 | 99.21% | 1 |

The refine spans **more** gaps than the shipped path -- five against four -- so
the trade is spans against flips, not a strict loss. It is the flips that make it
unusable as written, and they turned out to be a bug rather than the price of the
extra span.

The damage is one window per arm and it is the same failure mode both times: the
window spans and is **switched**. Alone it is `55,843,827` at 53.89%; with the
retry it is `24,105,188` at 87.07%, whose join rests on a single usable read --
`24,121,714 -> 24,142,287` is 20.6 kb with two spanning reads, of which one gets
an allele at the site and the other is three edits from both.

Three gates were tried and none separates the good joins from that one:

- **Verification** does not: all eight sites the retry admits at
  `55,843,827-55,889,113` carry `msa_verified = 1`, yet one segregates at or
  above 0.90 against read truth and five sit below 0.70.
- **A homopolymer/depth screen** does not, and is exactly inert here -- identical
  numbers with and without it -- because the phantom `24,131,708` is a 1 bp
  deletion at 39/29, allele fraction 0.427: a textbook heterozygote by depth
  whose net lengths (`+0` at 17/15, `-1` at 14/8, `-2` at 6/2) separate nothing.
- **Minimum link support** trades spans for flips without recovering the win: at
  5 reads the panel is down to 2 spans for 1 flip, where the retry alone already
  gives 4 spans for 0.

## The flips were a bug: a phantom bridging the k-means

The second round's damage is not the price of its extra span. Tracing
`24,105,188` in the refine+retry arm, the chain evidence reads:

```
24,105,188   agree=53  conflict=0     the left boundary
24,121,714   agree=10  conflict=1     the informative multiallelic site
24,131,708   agree=9   conflict=5     the phantom -- and the bridge
24,142,287   agree=0   conflict=0     absorbed with NO evidence
```

and by truth the whole right flank -- thirteen sites each segregating at 1.000 --
was inverted relative to the left block. The linker itself is not at fault: at
`collect_phase.cpp:860` a het with no supported link correctly starts a **new**
phase set. The merge happens inside the k-means, which is a global clustering
with no pairwise evidence test, so any admitted site can act as a bridge. The
phantom `24,131,708` supplied the path on a 9-against-5 vote -- a net margin of 4
that the link loop's repeat rule would refuse.

The screen let it through because of its depth escape: the phantom is a 1 bp
deletion at 39/29, allele fraction 0.427, a textbook heterozygote by depth that
segregates at 0.636. Removing the escape -- a homopolymer indel is admitted only
when it carries both of the locus' alleles -- fixes it:

| arm | spanned | in-gap hets | concordance | c -> d |
|---|---:|---:|---:|---:|
| stock defaults | 0 of 6 | 12 | 99.68% | **0** (byte-identical) |
| `--retry-unphased-with-bam` | **4 of 6** | 31 | **99.49%** | **0** |
| `--msa-verified-refine --retry-unphased-with-bam` | **4 of 6** | 29 | **99.49%** | **0** |

`24,105,188` now declines to span rather than spanning inverted, which is the
honest outcome for a window whose only crossing evidence is one usable read.

## Still open: the refine-alone arm, and where the real gate belongs

`--msa-verified-refine` without the retry still loses `55,843,827` at 53.89%
(152 concordant-to-discordant), and the bridges there are **not** homopolymers,
so no class screen reaches them: `55,871,837` glues on 11 agree against 7, and
`55,883,020` on **9 against 8**. Both segregate weakly against truth (0.778 and
0.839) while the flank they invert segregates at 1.000.

That is the general statement: a bridge is identified by its **link evidence**
being near-even, not by its site class, and the check that encodes exactly that
standard already exists at `collect_phase.cpp:811` -- both-haplotype support plus
`|agree - conflict| >= min_block_link_reads` -- but it is gated on
`recovery_graph`, so it never runs in the normal pipeline.

Ungating it was tried and reverted: it is too blunt as written, costing the
shipped retry arm three of its four spans (4 -> 1) and introducing a flip where
there had been none. The right version has to apply the margin to the evidence
that glues the k-means' components without discarding the well-supported joins,
which is a change to where the check sits, not to the standard it applies.
