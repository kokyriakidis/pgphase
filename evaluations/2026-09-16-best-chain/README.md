# Is there a chain through the union of BAM and graph sites?

`chain_search.py` takes every site available in a panel window -- the pipeline's
own candidates (alignment-derived plus injected catalog sites) **and** every graph
catalog site, most of which never become candidates -- genotypes each one directly
from the alignment, and searches for a path of sites linking the left flank to the
right flank. Sites are partitioned without truth (a substitution splits on the
base, an indel on net length across its tract) and the path is chosen without
truth; truth is applied only afterwards, so a chain reported here is one the
pipeline could in principle find.

```sh
python3 evaluations/2026-09-16-best-chain/chain_search.py \
  --panel evaluations/2026-09-16-test-panel/panel.tsv --arm /tmp/panel-noisy \
  --catalog test_data/chr20.sites.striped.vcf.gz --bam test_data/HG002...bam \
  --truth-map /tmp/truth_hap.tsv --min-shared 2 --corridor 5000 \
  --outdir evaluations/2026-09-16-best-chain
```

## The union is not short of sites

Inside each gap corridor the union yields **113-292** sites with a usable
two-group read partition, against the 2-3 the pipeline phases. The catalog alone
holds 575 sites inside `24,105,188-24,142,287`, where the pipeline carries 20
candidates.

## Three truth-free selection rules, all of which fail

| rule | result |
|---|---|
| max-bottleneck path, links >= 5 shared reads | chain in 2 of 6 windows, both correct -- and both are the windows that already span |
| max-bottleneck path, links >= 2 shared reads, free routing | chain in 6 of 6, but only **3 correct**; bottleneck does not separate them (correct 14/6/3 reads, switched 4/4/2) nor does link consistency (correct 1.000/0.977/0.900, switched 0.929/0.919/0.900) |
| same, but the path must cross the gap (5 kb corridor) | chain in 3 of 6; the one new chain runs on a 2-read link over 21.3 kb at consistency **1.000** and is **switched** -- two reads agreeing perfectly is not confirmation |
| keep only sites agreeing with a majority of their neighbours at >= 0.95 | selects the **phantoms**: for `5,309,406` it keeps 78 sites, every one of them chance-level. The phantoms are the majority and they agree with each other |
| aggregate: every site votes, reads split by the leading eigenvector of the read-read agreement matrix, no per-link gate | 53.67 / 57.75 / 68.02 / 54.02 / 63.64% -- chance -- and 88.07% only in the window that already spans |

## Why they fail: the per-read allele call, not the site set

Per-site segregation against read truth, inside the panel's gaps:

| partition made on | sites | >= 0.90 | < 0.70 | median |
|---|---:|---:|---:|---:|
| substitution (base) | 13 | **12 (92%)** | 1 | **1.000** |
| indel (net length) | 879 | 120 (14%) | 737 (84%) | **0.573** |

So the hundreds of extra sites the union contributes are noise *as we call them*,
and the reliable material is **1-3 substitutions per gap**. Those are spaced
**22-67 kb apart with zero reads covering two consecutive ones**, which is why
every chain rule above fails: gating links rejects everything, and aggregating
lets 84% noise swamp the three real sites.

Where each chain dies is the same story. The frontier site is itself a phantom in
every failing window -- `48,203,445 AC>A` at 0.535 over 71 reads, `24,125,663
GG>GC` at 0.588 over 68, `12,755,317 AT>A` at 0.562 over 48 -- and in
`12,717,796` two substitutions 0.6 kb ahead segregate **1.000** over 47 reads yet
link to that frontier at only 0.600. Plenty of reads (40-66 span the failing
jumps of 0.6-6.0 kb); no usable allele.

## What the competitor actually does

Site-level orientation across each gap, at the flanking substitutions truth
resolves -- the honest test, since read-level accuracy is structurally blind
across an interval no read spans:

| window | hiphase | longphase | whatshap |
|---|---|---|---|
| 55,843,827 | CORRECT | no join claimed | CORRECT |
| 24,105,188 | CORRECT | no join claimed | no join claimed |
| 5,309,406 | CORRECT | no join claimed | CORRECT |
| 12,717,796 | CORRECT | no join claimed | **SWITCHED** |
| 39,838,293 | CORRECT | CORRECT | CORRECT |

Hiphase is right **5 for 5** where no read spans the flanking pair, which is not a
coin flip. And its chain is small: **3-7 phased hets per gap, 1-3 of them
substitutions and the rest indels, with every consecutive step read-linked** --
widest step 8.5-21.4 kb, carried by 2-27 spanning reads.

So the chain that closes these gaps is about five sites long, it bridges the wide
substitution spacings **through the indels in between**, and we already hold every
one of those sites. The binding constraint is the **per-read allele call at those
indels**: exact-allele matching admits 13 of 74 reads at `24,121,713`, and the
net-length substitute tested here segregates at chance for 84% of indel sites.

Chain selection is therefore not the bottleneck and no amount of extra sites
helps. The work is upstream: call each read's allele at an indel by comparing it
against both candidate haplotype sequences rather than by exact match or by a
crude net-length split, and represent the locus multiallelically so both
haplotypes' alleles exist to be matched.
