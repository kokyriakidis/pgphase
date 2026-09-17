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

## Do we retrieve the gap sites correctly? Allele, depth and assignment audited

`audit_retrieval.py` asks the three questions separately for every candidate
inside the panel's gaps, measured against the alignment directly: is the emitted
allele a real mode of the read distribution and is every supported mode emitted;
does DP match the reads that actually cover the locus; and of the reads a record
does count, are they parentally pure -- the call being right, as distinct from
complete.

Conventions were calibrated against the emitted VCF at shared loci: an indel
record's POS is the VCF POS + 1 with the anchor base dropped, so an insertion's
delta is `+len(ALT)` and a deletion's is `-len(REF)`, while a **substitution's POS
is the VCF POS itself**. Reading substitutions one base to the left scored the
wrong column and returned chance-level purity, which is how that was caught.

| site class | n | DP / coverage | reads matching an emitted allele | loci with an unemitted mode (>= 5 reads) | assignment purity >= 0.95 |
|---|---:|---:|---:|---:|---:|
| **het substitutions** | 12 | **0.997** | **99.2%** | **0 of 12** | **12 of 12** |
| **het indels** | 24 | 0.884 | **77.3%** | **15 of 24** | **11 of 23** |
| clean hom substitutions | 272 | 0.992 | 99.8% | 0 | n/a (hom) |
| clean hom indels | 30 | 0.983 | 84.2% | 13 | n/a (hom) |
| noisy hom indels | 24 | **0.545** | 68.1% | 19 | n/a (hom) |

**Substitutions are retrieved correctly on all three counts.** Depth is 0.997 of
true coverage, essentially every covering read matches an emitted allele, no locus
hides a supported mode, and every one of the twelve het substitutions assigns its
reads to a single parent at >= 0.95.

**Het indels fail all three.** The worst in-gap loci:

| gap | position | covering reads | DP | reads matching an emitted allele | assignment purity | biggest unemitted mode |
|---|---|---:|---:|---:|---:|---|
| 5,309,406 | 5,339,369 | 65 | 36 | **0%** | unmeasurable | **+34 on 21 reads** |
| 55,843,827 | 55,862,240 | 61 | 60 | 67% | 1.000 | **+14 on 18 reads** |
| 48,183,976 | 48,204,384 | 71 | 71 | 79% | 0.931 | -2 on 12 reads |
| 24,105,188 | 24,121,714 | 74 | **8** | 39% | 1.000 | +3 on 11 reads |
| 39,838,293 | 39,845,582 | 71 | 66 | 73% | **0.600** | -2 on 11 reads |
| 55,843,827 | 55,882,617 | 68 | 66 | 75% | **0.522** | -2 on 10 reads |

Three distinct defects, and they are not the same defect:

1. **Incomplete allele sets.** 15 of 24 het indel loci have a read-length mode
   carried by at least five reads that no record represents -- 21 reads at
   `5,339,369`, 18 at `55,862,240`. At `5,339,369` **not one** of the 65 covering
   reads matches anything we emitted. Exactly one locus per window is emitted
   multiallelically, so the second haplotype's allele is usually simply absent.
2. **Short depth at individual loci**, even though the class mean looks close.
   `24,121,714` records DP 8 against 74 covering reads; `5,339,369` records 36 of
   65. Noisy hom indels as a class sit at 0.545 of coverage.
3. **Wrong assignments, not merely missing ones.** Purity is 0.600 at
   `39,845,582` and 0.522 at `55,882,617`: of the reads assigned to an allele
   there, the parental split is near even, so those calls are incorrect rather
   than absent. Elsewhere -- `55,862,240`, `24,121,714` -- purity is 1.000, so the
   reads that are counted are counted correctly and only coverage is missing.

So the answer is: **yes for substitutions, no for indels.** The allele identity we
discover is often right -- at `24,121,714` we do emit both `+4` and `+8`, the two
true modes -- but with DP 6 and 2 against 13 and 16 reads actually carrying them,
each as a separate biallelic record claiming AF 1.000 with no reference read.
Every gap's chain depends on indels for its widest steps, which is why this
matters more than the substitution result suggests.

## Does the injection work correctly, and inject in the right way?

Compared the hybrid arm against `collect-bam-variation` run standalone over the
same six windows, so the only differences are the catalog claim and what follows
from it. The hybrid arm is the `--keep-noisy-kmeans` one, which matches the
alignment channel's noisy handling; without that flag the comparison would also
carry the hybrid-only `skip_noisy_kmeans` override and would not isolate
injection.

**Faithful in what it adds.** Injection adds 6-16 candidates per window
(152->158, 245->252, 254->263, 302->311, 350->366, 242->253) and no locus the
alignment already had is dropped.

**Complete with respect to real heterozygotes.** Each gap holds 636-896 catalog
sites and only 2-130 are injected -- but of the catalog SNPs inside the gaps that
read truth confirms as heterozygous, **none is missing** in any of the six
windows. The catalog's remaining hundreds of in-gap sites are not heterozygous in
this sample, and declining them is correct rather than a loss.

**Its verdicts are right where the allele representation is right.** Injection
changes the record at 128 shared loci. **106 are category changes and 22 are not**
-- in those 22 the `CATEGORY` field is identical and only `INIT_CAT` moves, 17 of
them `LOW_COV -> REP_HET_INDEL` with the counts untouched and 5 `LOW_COV ->
LOW_COV` where only the counts move. Of the 106 category changes, 105 are
promotions and one is the wrong verdict dissected below:

| transition | loci |
|---|---:|
| `NOISY_CAND_HOM` -> `CLEAN_HOM` | 68 |
| `NOISY_CAND_HET` -> `CLEAN_HET_SNP` | 36 |
| `NOISY_CAND_HET` -> `NOISY_CAND_HET` (not a promotion: `INIT_CAT` `LOW_COV` -> `REP_HET_INDEL`) | 17 |
| `NOISY_CAND_HOM` -> `NOISY_CAND_HOM` (not a promotion: counts only) | 5 |
| `NOISY_CAND_HET` -> `CLEAN_HET_INDEL` | 1 |
| `NOISY_CAND_HET` -> `CLEAN_HOM` | **1** |

**All 36 het promotions are true heterozygotes** -- purity of the ALT-carrying
reads is 1.000 at 35 of them and 0.967 at `24,169,753`. And where injection
changes DP it usually corrects it: of 80 DP changes, **58 move closer to true
coverage** and several dramatically so (`5,325,566` 7 -> 69 against 71 covering
reads, `5,323,528` 12 -> 68 against 70, `5,319,294` 27 -> 70 against 70,
`5,335,995` 30 -> 74 against 74), 21 move farther by 1-2 reads, and one is
unchanged in distance.

**The single wrong verdict is the multiallelic defect again, with a worse
consequence.** At `55,896,396` (`TTTTTTTTTTTTTT>.`, a 14 bp deletion in a T
homopolymer) the alignment channel says `DP 30, 0/30, NOISY_CAND_HET` and the
hybrid says `DP 37, 2/35, AF 0.946, CLEAN_HOM`. By truth the locus is one of the
cleanest heterozygotes in that window:

| net length | MAT | PAT |
|---:|---:|---:|
| -16 / -17 | 0 | **32** |
| -3 / -4 | **32** | 0 |

Both haplotypes are non-reference, so a single `-14` record sees 35 alt against 2
reference reads and is promoted to homozygous -- and a homozygote is dropped from
phasing entirely. The promotion logic is behaving correctly on an allele set that
misrepresents the locus.

So: injection is faithful, complete with respect to real heterozygotes, and its
promotions and depth backfills are correct -- with the two qualifications above,
that 22 of its 128 record changes promote nothing (they relabel `INIT_CAT` or move
counts only) and that one of the 106 category changes is wrong. It is not the defect. But its verdicts
are only as good as the allele representation handed to it, and where both
haplotypes are non-reference that representation converts an informative het into
a homozygote.
