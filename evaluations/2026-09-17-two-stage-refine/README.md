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

| arm | spanned | in-gap hets | concordance | c -> d |
|---|---:|---:|---:|---:|
| stock defaults | 0 of 6 | 12 | 99.68% | 0 |
| `--retry-unphased-with-bam` (shipped) | **4 of 6** | 31 | **99.49%** | **0** |
| `--msa-verified-refine` | 4 of 6 | 25 | 92.19% | **152** |
| `--msa-verified-refine --retry-unphased-with-bam` | **5 of 6** | 30 | 97.46% | **64** |
| ... `--min-block-link-reads 3` | 3 of 6 | 30 | 97.18% | 65 |
| ... `--min-block-link-reads 5` | 2 of 6 | 30 | 99.21% | 1 |

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

So the second round ships **off by default**, with its cost measured. The
shipped path dominates it on every axis the panel measures, and what would make
the refine worth enabling is not a better site screen but a join that refuses to
merge across an interval carrying one usable read -- the same conclusion the
chain search reached from the other direction.
