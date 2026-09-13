# How the phasing field evaluates, and how we should

Reference analysis of the HiPhase benchmark (the most recent and most directly
comparable), what it does and does not specify, and what we must add to be read
as comparable. Written after reproducing their metric on our data.

> **Correction (2026-09-12): shared-call-set pgphase rows below are retracted.**
> The DeepVariant input already contained 67,127 phased heterozygotes, and the
> old `phase_vcf_from_hp.py` preserved inherited phase at sites without pgphase
> support. The corrected script clears phase/PS from every heterozygote first.
> Correct chr20 graph-transfer results are 59,015 assessed, 276 blocks,
> N50/NGC50 215/204 kb, 15 switchflips, and Hamming 16. The best experimental
> all-linear-site plus conservative PS-merge result is 59,447 assessed, 247 blocks,
> N50/NGC50 401/374 kb, 39 switchflips, and Hamming 79. Competitor rows are
> unaffected. Historical experiment sections below remain diagnostic only; use
> `docs/HANDOFF.md` for the current comparison.

Sources: [paper](https://academic.oup.com/bioinformatics/article/40/2/btae042/7588891)
(Holt et al., *Bioinformatics* 2024, btae042) ·
[performance.md](https://github.com/PacificBiosciences/HiPhase/blob/main/docs/performance.md) ·
[repo](https://github.com/PacificBiosciences/HiPhase)

---

## 1. What HiPhase actually did

### Experimental design

| element | HiPhase benchmark |
|---|---|
| samples | HG001, HG002, HG005 (docs); 3 replicates of HG002 on Revio (paper) |
| coverage | ~30x |
| platform | Sequel II (docs), Revio (paper) |
| reference | GRCh38 |
| scope | **full autosomes** |
| alignment | pbmm2 |
| variant input | DeepVariant small variants; pbsv SVs; TRGT tandem repeats |
| truth | NIST/GIAB **v4.2.1 phased variant sets** |
| error assessment | **`whatshap compare`** |
| competitors | WhatsHap v1.4, twice: default `--indel`, and "optimised" `--indel --distrust-genotypes` |

Every tool receives **the same BAM and the same VCF**. That is the fairness
convention in this field and we must match it.

### Metrics they report

1. **NG50** — block length at which blocks >= that length cover 50% of the
   reference. Genome-denominated (3.1 Gb).
2. **NGC50** — "corrected NG50": blocks are *split at every switchflip* before
   recomputing. Their headline contiguity number, and the one that resists the
   "long but wrong" criticism.
3. **Switchflips** — switches plus flips, summed.
4. **Hamming distance** — summed, blockwise.
5. **Phased variants** — het variants emitted with `|`.
6. **Genes fully phased (%)** — genes whose het variants all share one phase set.
7. **Wall clock** — noting WhatsHap single-threaded vs HiPhase 16 threads.

### Their results (mean/sum over 3 samples)

| metric | WhatsHap | WhatsHap opt. | HiPhase | HiPhase+SV |
|---|---|---|---|---|
| NG50 | 248,730 | 238,626 | 310,547 | **312,620** |
| Switchflips | 4,196 | **1,856** | 2,888 | 2,768 |
| NGC50 | 215,392 | 225,435 | 285,874 | **287,129** |
| Hamming | 110,379 | **71,497** | 122,866 | 125,239 |
| Genes fully phased | 68.35% | 67.65% | 71.78% | **71.89%** |
| Wall clock (s) | 5,580 | 7,403 | **5,461** | 12,178 |

**Worth noting for our framing:** HiPhase does *not* sweep its own table.
WhatsHap-optimised beats it on switchflips (1,856 vs 2,768) and on Hamming
(71,497 vs 122,866). HiPhase wins contiguity and gene phasing. A method that
wins switches, switchflips *and* Hamming simultaneously would be a stronger
claim than theirs.

### What is not reproducible from their materials

The repo contains source and docs only — **no benchmarking scripts**. Absent:

- exact `whatshap compare` invocation and filters
- how phase blocks were extracted (presumably `whatshap stats --block-list`)
- the NGC50 implementation
- the gene annotation used for "genes fully phased"
- exact data accessions

So we cannot replicate their numbers; we can only replicate their *method*. Our
`scripts/compute_ngc50.py` already reimplements NGC50 from the description
(blocks -> BED, subtract switch-error BED, NG50 over 3.1 Gb).

---

## 2. Where our evaluation differs

| | HiPhase benchmark | ours (to date) |
|---|---|---|
| unit of truth | phased **variants** (GIAB VCF) | phased **reads** (diploid assembly) |
| truth source | NIST v4.2.1 | Q100/v5.0q assembly-derived |
| scope | full autosomes, 3 samples | chr20/18/12, HG002 only |
| contiguity | NG50 + NGC50 | N50, auN, median span |
| error metric | switchflips, Hamming (variant) | switch, flip, Hamming (read) |
| gene phasing | reported | script exists, not run |

Both units are legitimate and they answer different questions. Read-level
accuracy is what matters for haplotagging and for supplying HP tags to a caller;
variant-level is the field's convention. **We should report both**, and lead
with whichever the venue expects.

### Their exact methodology (from the supplementary, section 3.2-3.3)

The supplementary publishes the command templates the repo omits:

```
whatshap compare --tsv-pairwise {tsv} --switch-error-bed {error_bed} {truth} {vcf}
whatshap stats   --tsv {tsv} --block-list {blocks} {vcf}
bedtools subtract -a {phase_block_bed} -b {error_bed} > {corrected_bed}   # -> NGC50
bedtools intersect -a {refseq_gene_bed} -b {phase_block_bed} -f 1.0 -wa -A # -> genes
```

- **switchflips = switches + flips**, parsed from `all_switchflips` and summed
  over chromosomes. **Hamming** is `blockwise_hamming`, likewise summed.
- **NGC50**: phase blocks -> BED, `bedtools subtract` the switch-error BED to
  split blocks at every error, then recompute NG50 on the surviving sub-blocks.
- **Genes fully phased**: RefSeq GRCh38 GFF3 (NCBI Annotation Release 110),
  keeping `gene`/`pseudogene` from BestRefSeq/RefSeq/Gnomon/Curated Genomic on
  primary chromosomes. Phase blocks are first **extended** in both directions
  until the next het variant, to absorb homozygous stretches inside genes, then
  a gene counts as phased if one extended block covers it entirely (`-f 1.0`).
- They note `whatshap stats` under-reports HiPhase because it mishandles
  multi-allelic sites, so HiPhase's own `--blocks-file`/`--summary-file` were
  used for its block metrics. Any comparison must use one block source per tool
  consistently and say which.

### Reproducing their metric on our data

`whatshap compare` against the v5.0q CHM13 phased truth, HG002 chr20, all tools
given the same reads:

| method | assessed | switches | switchflips | sf rate | Hamming | ham rate |
|---|---|---|---|---|---|---|
| whatshap 2.8 | 62,509 | 475 | 143/166 | 0.004943 | 3,196 | 0.050975 |
| LongPhase 2.0.2 | 54,679 | 188 | 98/45 | 0.002615 | 605 | 0.011021 |
| HiPhase 1.6.0 | 62,264 | 238 | 106/66 | 0.002762 | 2,279 | 0.036509 |
| **graph (margin 2)** | 47,894 | **24** | **8/8** | **0.000334** | **19** | **0.000394** |

7.8x fewer switches, 7.8x lower switchflip rate, 32x lower Hamming than the best
competitor — on their metric, with their tool.

**But that table is not apples-to-apples**, and correcting it changes the
picture. Above, each tool phased a different variant set: competitors phased
DeepVariant's calls, ours phased its own. `scripts/phase_vcf_from_hp.py`
transfers our read HP tags onto the *same* DeepVariant call set (per-site
majority vote of tagged reads), so all four phase identical variants:

| method | phased | switchflips | Hamming | N50 kb | NGC50 kb |
|---|---|---|---|---|---|
| whatshap 2.8 | 76,567 | 309 | 3,196 | 673 | 338 |
| LongPhase 2.0.2 | 62,287 | **143** | 605 | 475 | 400 |
| HiPhase 1.6.0 | 75,060 | 172 | 2,279 | **856** | **644** |
| graph HP -> DV calls | 71,504 | 145 | **235** | 329 | 312 |

Honest reading of the corrected table:

- **Hamming: we win clearly** — 235 against LongPhase's 605 (2.6x) and
  HiPhase's 2,279 (10x). Hamming counts variants placed on the wrong haplotype,
  which is the quantity that matters for haplotagging and for supplying HP tags
  downstream.
- **Switchflips: we tie LongPhase** (145 vs 143), both well ahead of HiPhase
  (172) and whatshap (309).
- **Contiguity: we lose.** NGC50 312 kb against HiPhase's 644. Transferring onto
  a larger call set exposes it: 10,108 DV het sites had too little
  haplotype-tagged support and 2,592 were mixed, and each unphased site can
  break a block.

So the defensible claim is **accuracy, not contiguity** — and specifically
Hamming, where the margin is large and consistent. Claiming a contiguity win
would not survive this table.

(NGC50 here is chr20-denominated, 66 Mb. HiPhase's published ~310 kb is
genome-denominated over 3.1 Gb of autosomes. The two are not comparable; only
the within-table ordering is.)

---

## 3. Gaps to close before this is publishable as comparable

### Blocking

1. **Identical variant sets.** Competitors phased DeepVariant's calls; we phased
   our own, so we assessed 47,894 het pairs against their 55-62k. Rates
   normalise for this but a reviewer will still object. **Fix:** transfer our
   read HP tags onto the same DeepVariant VCF (majority vote of tagged reads per
   het site) so all four phase an identical call set. Small script, removes the
   objection entirely.
2. **NG50/NGC50 need whole-genome runs.** Their NG50 ~310 kb is denominated on
   3.1 Gb across full autosomes. Single-chromosome runs give 0 — which is
   exactly what our summary reports. Either run full autosomes or report
   per-chromosome N50 and state the difference explicitly. Do not put our
   chr20 N50 next to their genome-wide NG50.
3. **NGC50 must be reported alongside N50.** It is their headline and it is the
   metric that neutralises "long but wrong" blocks. Given our switch counts are
   ~10x lower, NGC50 should favour us far more than raw N50 does — this is
   likely our strongest contiguity argument, and we are currently not making it.

### Expected by reviewers

4. **Multiple samples.** They use three. We have one (HG002). HG003 or HG005.
5. **Genes fully phased.** `scripts/gene_phasing_completeness.py` exists; needs
   a gene BED (`scripts/prepare_refseq_gene_bed.sh` exists too).
6. **Wall clock and threads.** They report it; cheap for us to add, and our
   runtime (~36 s/chromosome at 20 threads) is likely favourable.
7. **WhatsHap in both modes.** Default *and* `--distrust-genotypes`. The latter
   is their strongest competitor on error metrics and omitting it looks like
   cherry-picking.

### Optional but differentiating

8. **SV and tandem-repeat phasing.** HiPhase's actual novelty. We do not do it
   and should say so plainly rather than let a reviewer notice.
9. **Read-level accuracy as our own contribution.** No competitor reports
   haplotagging accuracy against an assembly. It is both a fairer measure of
   what phasing is used for and a place where our margin is largest.

---

## 4. Recommended table for the paper

One table, all tools, same reads, same variant set, both units:

```
                     variant level (whatshap compare)      read level (assembly truth)    contiguity
tool        version  assessed  switches  switchflips  ham   reads  discordant  hamming   N50   NGC50
whatshap    2.8      ...
whatshap    2.8 -dg  ...
LongPhase   2.0.2    ...
HiPhase     1.6.0    ...
longcallD   ...
ours        ...
```

with median block span beside N50, since all methods produce ~870 kb typical
blocks and the N50 spread comes from a handful of long blocks.

---

## 5. Why we lose contiguity, and what is actually fixable

The corrected table in §2 leaves one clear deficit: NGC50 312 kb against
HiPhase's 644. This section takes that apart. Three candidate causes were
tested; two were eliminated and one was a bug in our own harness.

### 5.1 Chunking is not the cause (eliminated)

`collect-graph-variation` phases in chunks and stitches. If a block could not
outlive a chunk, contiguity would be capped by construction. Raising the chunk
size 20x does nothing:

| chunk size | blocks | N50 kb | max block kb | switchflips | Hamming |
|---|---|---|---|---|---|
| 500 kb (default) | 366 | 217 | 1,235 | — | — |
| 2 Mb | 369 | 212 | 1,199 | 15 | 17 |
| 10 Mb | 371 | 206 | 1,199 | 12 | 12 |

Blocks already span 2.5 default chunks, so stitching works. Blocks break for a
local reason, not at chunk seams. Chunk size stays at 500 kb.

### 5.2 The block-link threshold is a real lever with a bad exchange rate

The break rule in `iter_update_var_hap_cons_phase_set` starts a new phase set
when fewer than `min_block_link_reads` reads either agree or conflict across a
site. It was hardcoded at 2; it is now `--min-block-link-reads`.

| setting | blocks | N50 kb | NGC50 kb | switchflips | Hamming |
|---|---|---|---|---|---|
| 2 (default) | 366 | 217 | 198 | 16 | 19 |
| 1 | 313 | 302 | 280 | 18 | **247** |

Contiguity improves 41% for a 13x worse Hamming. Since Hamming is the metric we
actually win on, this is not a trade worth taking, and the default stays at 2.
The flag is exposed so the exchange rate can be reported rather than assumed.

### 5.3 Anatomy of our block breaks

Of our 365 native breaks, 231 fall strictly inside a HiPhase block. But the
median gap is **30 kb** — wider than a HiFi read, so HiPhase is not linking
those by read evidence either. Counting what lives in each gap:

| gap width | n | truth hets | our hets | HiPhase hets |
|---|---|---|---|---|
| 0–1 kb | 4 | 4 | 1 | 0 |
| 1–10 kb | 3 | 0 | 1 | 7 |
| 10–50 kb | 273 | 3 | 1 | 3 |
| >50 kb | 84 | 7 | 1 | 7 |

**This table is wrong and is kept only to mark the error.** Phase blocks nest
and overlap, so consecutive-block differences are not the uncovered regions; the
counts above are diluted by variants that other blocks do cover. §7 recomputes
the gaps as the complement of the *merged* blocks and reaches the opposite
conclusion: the unphased regions are het-**rich**, at 1.8x the chromosome
average. Any gap analysis must merge blocks first.

The variant-sensitivity reading survives the correction; the "het-poor"
explanation does not.

### 5.4 longcallD settles it: we lose *coverage*, not contiguity

Adding longcallD to the native-call-set table is what makes the cause
unambiguous. All tools, chr20 HG002, each on its own calls, against v5.0q:

| method | assessed | blocks | N50 kb | NGC50 kb | switches | switchflips | Hamming |
|---|---|---|---|---|---|---|---|
| whatshap 2.8 | 62,510 | 228 | 673 | 338 | 475 | 309 | 3,196 |
| HiPhase 1.6.0 | 62,265 | 194 | **856** | **644** | 238 | 172 | 2,279 |
| LongPhase 2.0.2 | 54,680 | 226 | 475 | 400 | 188 | 143 | 605 |
| longcallD (default) | 65,750 | 431 | 319 | 273 | 117 | 83 | 884 |
| **graph (margin 2)** | 47,895 | 366 | 217 | 198 | **24** | **16** | **19** |
| graph (link-reads 1) | 47,946 | 313 | 302 | 280 | 27 | 18 | 247 |

longcallD is the second-best tool on accuracy — ahead of HiPhase and whatshap on
every error metric — so it is the right yardstick. We beat it 5.2x on
switchflips and 46x on Hamming.

Now the block geometry, which is the decisive panel:

| method | blocks | total span Mb | % of chr20 | mean kb | median kb | longest kb |
|---|---|---|---|---|---|---|
| whatshap 2.8 | 228 | 56.6 | 85.4 | 248 | 57 | 3,012 |
| HiPhase 1.6.0 | 194 | 57.7 | 87.2 | 298 | 67 | 3,024 |
| LongPhase 2.0.2 | 226 | 54.1 | 81.7 | 239 | 91 | 1,721 |
| longcallD | 431 | 55.2 | **83.3** | 128 | 41 | 1,236 |
| **graph (margin 2)** | 366 | 46.5 | **70.2** | 127 | 39 | 1,235 |

Our block-size distribution is **indistinguishable from longcallD's** — mean 127
vs 128 kb, median 39 vs 41 kb, longest 1,235 vs 1,236 kb (the same locus, capped
by the same feature of the chromosome). We do not build shorter blocks, and
longcallD makes *more* of them than we do (431 vs 366) while still scoring a
higher NGC50.

The entire difference is that our blocks cover **70.2% of chr20 against 83.3%** —
about 8.7 Mb of chromosome that we never place in any block. NGC50 is
denominated over the whole chromosome, so uncovered span depresses it directly.

The cause is upstream of phasing: we phase **47,895** het sites where longcallD
phases 65,750 and HiPhase 62,265. `--min-read-margin` cannot help — it gates
which reads get an HP tag at BAM-write time and leaves the VCF untouched
(margins 0, 1 and 2 give byte-identical block structure: 366 blocks, 70.2%,
NGC50 198). This is **site sensitivity in the snarl catalog**, the same
limitation documented in `publication_plan.md` C1.5: het-poor catalog regions,
graph-blind satellite, and alt-vs-alt hets scored AF = 1.0 and dropped.

**Consequence for the paper.** The contiguity gap is not a block-linking
deficit and should not be chased with block-linking thresholds (§5.2 shows what
that costs). It is a *call-set* deficit. But it is worse than a contiguity
problem alone -- see §5.5, which shows the skipped sites also carry the errors.

### 5.5 The error advantage was site selection (claim retracted)

The tables above compare each tool on its own call set, so we phase 47,895 sites
where longcallD phases 65,750. That invites an obvious objection: are we winning
because we phase only the easy 73%? The answer is yes.

Restricting every tool to the 50,877 het sites we phase, then rescoring:

| method | assessed | switches | switchflips | Hamming |
|---|---|---|---|---|
| whatshap 2.8 | 47,998 | 37 | 28 | 1,981 |
| HiPhase 1.6.0 | 47,986 | 9 | 8 | 1,276 |
| LongPhase 2.0.2 | 46,446 | 10 | **7** | **15** |
| longcallD | 48,233 | **4** | **3** | 488 |
| **graph (margin 2)** | 47,895 | 24 | 16 | 19 |

Thinning alone does not explain this. A size-matched *random* subset of each
tool's own het sites is the control:

| method | full | random subset | our subset |
|---|---|---|---|
| whatshap 2.8 | 309 sf / 3,196 ham | 242 / 2,038 | **28 / 1,981** |
| HiPhase 1.6.0 | 172 / 2,279 | 156 / 1,429 | **8 / 1,276** |
| LongPhase 2.0.2 | 143 / 605 | 131 / 565 | **7 / 15** |
| longcallD | 83 / 884 | 66 / 791 | **3 / 488** |

Random removal of 23% of sites removes ~20% of errors, as expected. Removing
*our* skipped 23% removes **~90%** of them. The sites we do not phase are
precisely where every tool makes its mistakes.

Two conclusions, both of which change the paper:

1. **The switchflip claim is retracted.** On an equal site set we are *last*:
   16 switchflips against longcallD's 3, LongPhase's 7, HiPhase's 8. The
   published-looking 5.2x win over longcallD was an artifact of comparing
   different site sets.
2. **The Hamming lead is real but partly structural.** We beat longcallD 19 vs
   488 and HiPhase 19 vs 1,276 on identical sites, but lose to LongPhase
   (15). Blockwise Hamming is computed per block and penalises long blocks, and
   ours are the shortest of any tool (§5.4), so some of the margin is block
   geometry rather than accuracy. It should not be quoted without NGC50 beside
   it.

**The target changes accordingly.** Skipping hard sites is not a strategy -- it
inflates every error metric while depressing every contiguity metric at the same
time. The goal is to phase the *full* call set and still be more accurate on it,
using the pangenome alignment as the source of the extra evidence. That is a
harder claim, and the only one worth making.

**Scope of the retraction.** It applies to comparisons on our *native* call set,
which is too sparse (47,895 sites, 70.2% of chr20). It does **not** apply once
every tool is given the same dense call set -- see §5.6, where we assess 98% of
what HiPhase assesses and the accuracy win survives intact. The lesson is that
our own variant emission, not our phasing, was producing the unfair comparison.

### 5.6 On a shared call set the accuracy claim holds

HiPhase, LongPhase and whatshap are phasers, not callers: all three phase the
same DeepVariant VCF. Transferring our read haplotypes onto that VCF puts all
four on one call set. chr20 HG002 vs v5.0q:

| method | assessed | blocks | % chr20 | N50 kb | NGC50 kb | switches | switchflips | Hamming |
|---|---|---|---|---|---|---|---|---|
| whatshap 2.8 | 62,510 | 228 | 85.4 | 673 | 338 | 475 | 309 | 3,196 |
| HiPhase 1.6.0 | 62,265 | 194 | 87.2 | **856** | **644** | 238 | 172 | 2,279 |
| LongPhase 2.0.2 | 54,680 | 226 | 81.7 | 475 | 400 | 188 | 143 | 605 |
| **graph -> DV (m0)** | 60,616 | 549 | **92.1** | 338 | 319 | **186** | **138** | **224** |

And the §5.5 control, re-run on this set -- restricting every phaser to the sites
we phase -- is now almost a no-op, which is the point:

| method | switchflips (full) | switchflips (our sites) | Hamming (full) | Hamming (our sites) |
|---|---|---|---|---|
| HiPhase 1.6.0 | 172 | 172 | 2,279 | 2,270 |
| LongPhase 2.0.2 | 143 | 143 | 605 | 605 |
| whatshap 2.8 | 309 | 309 | 3,196 | 3,185 |

Our site set now covers essentially all of theirs, so there is no site-selection
advantage left to strip out. On equal footing we win switchflips (138 vs 143 /
172 / 309) and Hamming (224 vs 605 / 2,279 / 3,196, a 2.7x margin over the best
competitor), while covering more of the chromosome than any of them.

Loosening the transfer gates barely matters, which is worth knowing:

| config | phased | unphased (support / mixed) | NGC50 | switchflips | Hamming |
|---|---|---|---|---|---|
| margin 2, min-reads 2 | 83.1% | 10,290 / 2,472 | 319 | 142 | 233 |
| margin 0, min-reads 2 | 83.8% | 9,429 / 2,660 | 319 | 138 | 224 |
| margin 0, min-reads 1 | 84.7% | 8,276 / 2,893 | 319 | 141 | 228 |

Tagging 95% of reads instead of 82% recovers only 0.7 pp of sites, and dropping
the support floor to one read trades `support` failures for `mixed` failures
one-for-one. The ~15% of het sites we cannot phase are not being withheld by a
gate -- they have no usable haplotype-tagged evidence at all.

### 5.7 What is actually left: fragmentation, not coverage

The one metric where the gap survives a fair comparison is contiguity: **549
blocks against HiPhase's 194**, NGC50 319 vs 644, *despite* our blocks covering
92.1% of chr20 against their 87.2%.

This is the opposite of the §5.4 diagnosis, and both are correct in their own
regime. On our sparse native call set the deficit was coverage (70.2%). On a
dense call set coverage is no longer the problem and the deficit is purely that
we emit ~2.8x as many phase sets. That is a property of the block-break rule
(§5.2), not of the call set, and it is the remaining open problem.

### 5.8 A bug in our own comparison harness (fixed)

Transferring our HP tags onto the DeepVariant call set produced **548** blocks —
worse than our native 366 — and 16% of them *overlapped* their neighbour, which
no real phaser does (HiPhase and LongPhase: 0%). Overlapping blocks are not a
phasing outcome; they are an artifact.

Cause: `phase_vcf_from_hp.py` pooled covering reads by `HP` alone. **An HP tag
is only meaningful within a phase set** — HP=1 in one block and HP=1 in the next
are unrelated, because each block's haplotype labelling is arbitrary. Every site
near a block boundary is covered by reads from two phase sets, whose votes then
cancel or scramble. Of 2,318 adjacent phase-set changes in the output, **1,682
were returns to a phase set already left** — flapping, not real breaks. Only 636
were genuine.

Fix: group covering reads by `(PS, HP)` and decide the site within the phase set
contributing the most haplotype-informative reads. This is not a tuning choice —
pooling across phase sets is simply incorrect, and it inflated both our block
count and our "mixed" (unphaseable) site count.

Any benchmark that transfers read haplotypes onto a shared call set has this
hazard, so it is worth stating explicitly in the paper's methods.

---

## 6. How HiPhase achieves its contiguity

On the shared DeepVariant call set HiPhase emits 194 phase blocks to our 549
(§5.7). This section works out why, from measurement rather than from its paper.

### 6.1 Its extra bridges are correct, not confident

Of our 466 block junctions, a single HiPhase block spans **295**. If HiPhase
were simply bridging on weak evidence, those junctions would be where it makes
its mistakes. They are not:

| | switch-error rate |
|---|---|
| HiPhase, inside the junctions we break at | **5.4%** (16/295) |
| HiPhase, size-matched random windows inside its own blocks | 5.0% |

Identical to its own baseline. So this is not a risk-appetite difference that we
could match by loosening a threshold — HiPhase genuinely links across these
junctions. §5.2 already showed that simply lowering `--min-block-link-reads`
buys contiguity at 13x the Hamming cost, which is what "bridging on weak
evidence" actually looks like. HiPhase is not doing that.

### 6.2 Our break rule is a chain, not a graph

`iter_update_var_hap_cons_phase_set` (`src/collect_phase.cpp:457-470`) counts
linking reads for exactly one pair of variants: `het_var_idx[hi-1]` and
`het_var_idx[hi]`, the two *adjacent* het variants. If that single pairwise link
is thin, a new phase set starts:

```cpp
if (n_agree[_vi] < opts.min_block_link_reads &&
    n_conflict[_vi] < opts.min_block_link_reads) {
    phase_set = var.key.sort_pos();     // break
}
```

A variant is therefore connected only to its immediate predecessor. Reads
linking it to the variant two or three back contribute nothing, so one weakly
supported variant severs a block that is otherwise densely spanned. Phasers that
build a connectivity *graph* over variants — where any read covering two
variants is an edge, and blocks are connected components — do not have this
failure mode.

Further, `check_agree_haps` consults `chunk.haps[read_i]`: a read only links if
it has already been assigned a haplotype. Our pipeline is **read-first**
(k-means assigns reads to haplotypes, then derives variants), so an unassigned
read is invisible as linking evidence even when it spans both variants cleanly.

### 6.3 The junctions are where reads are untagged — self-reinforcingly

Measured at the 245 narrow (<50 kb) junctions HiPhase bridges:

| junction width | n | spanning reads (median) | of those, HP-tagged | % |
|---|---|---|---|---|
| 0–5 kb | 47 | 34 | **0** | 0.0 |
| 5–20 kb | 108 | 11 | 5 | 45.5 |
| 20–50 kb | 90 | 1 | 0 | 0.0 |
| **overall** | 245 | 3,817 total | 1,208 | **31.6** |

Chromosome-wide at margin 0 about 95% of reads carry HP. At these junctions it
is 31.6%, and in the narrowest bucket 34 reads span end to end and *none* is
tagged. The reads are present; they are simply not usable as evidence under a
read-first design.

Partly this is mapping quality — junction-spanning reads are MAPQ-depleted
relative to block interiors:

| | median MAPQ | MAPQ = 0 | MAPQ < 20 | tagged |
|---|---|---|---|---|
| spanning our junctions | 60 | 5.4% | **39.0%** | 31.6% |
| inside our blocks | 60 | 1.9% | 19.3% | 68.5% |

But MAPQ is not the whole story: 61% of junction reads have MAPQ >= 20 while only
32% are tagged. The rest cover no *clean candidate site*, so k-means never scores
them.

This closes a loop: no tagged reads across the junction -> no pairwise link ->
block breaks -> the reads there are never brought into a phase set. The
fragmentation and the untagged reads cause each other.

### 6.4 What HiPhase does differently, and the implication

HiPhase also phases indels as linking evidence: of the 13,328 het variants it
phased inside our junction regions, **2,869 (21.5%) are indels**.

The architectural contrast is the point. Ours is read-first: assign reads to
haplotypes, then variants follow, so a read must already be confidently placed
to contribute anything. HiPhase is variant-first: variants are nodes, reads are
edges, and allele assignment is solved over the resulting component — a read
with no prior haplotype still contributes connectivity.

**This is the right place to attack, and it is the one that does not trade
against accuracy.** Loosening the break threshold (§5.2) costs 13x Hamming
because it makes the same fragile pairwise test more permissive. Replacing the
adjacent-pair chain with transitive connectivity over a window of preceding het
variants adds evidence rather than lowering a bar, which is consistent with
HiPhase bridging these junctions at its ordinary error rate.

---

## 7. Attempting the fix: two mechanisms, both negative

§6 argued the break rule was the place to attack. Both implementations work and
neither helps, and the measurement that explains why also relocates the problem.

### 7.1 What was built

Two flags, both defaulting to the previous behaviour:

- **`--block-link-window INT`** (default 1). The break rule linked only the
  immediately preceding het variant, so one weakly covered variant severed a
  block reads otherwise spanned. The window searches the previous N hets and
  takes the nearest link that clears `--min-block-link-reads`, chaining flip
  parity through whichever variant it actually linked to. At window 1 it
  reproduces the old chain exactly (verified: 366 blocks, NGC50 198, 16
  switchflips, Hamming 19 — identical to baseline).
- **`--link-by-alleles`** (default off). `check_agree_haps` consults
  `chunk.haps[read_i]` and so ignores untagged reads — and §6.3 showed junction
  reads are mostly untagged. `check_agree_alleles` links two variants by the
  allele pattern a read carries, no haplotype required.

### 7.2 Results: neither works

chr20, native call set, margin 2 + af 0.12:

| config | assessed | blocks | % chr20 | N50 kb | NGC50 | switchflips | Hamming |
|---|---|---|---|---|---|---|---|
| window 1, hap-linked (baseline) | 47,895 | 366 | 70.2 | 217 | 198 | 16 | **19** |
| window 10, hap-linked | 47,900 | 360 | 70.6 | 241 | 206 | 16 | 23 |
| window 1, allele-linked | 47,898 | 363 | 70.2 | 217 | 203 | 17 | 246 |
| window 10, allele-linked | 47,901 | 359 | 70.6 | 241 | 206 | 17 | 250 |

Widening the window moves 366 blocks to 360. Allele-linking adds nothing and
costs 13x Hamming — the same exchange rate as loosening the threshold (§5.2),
which is the signature of admitting noise rather than evidence.

### 7.3 Why: 84% of breaks cannot be spanned by any read

| gap width | n | median reads spanning | breaks with 0 spanning |
|---|---|---|---|
| 0–5 kb | 6 | 36 | 0 (0%) |
| 5–15 kb | 4 | 18 | 1 (25%) |
| 15–25 kb | 82 | 1 | 36 (44%) |
| >25 kb | 272 | 0 | **270 (99%)** |

**307 of 364 breaks (84%) have no read spanning end to end**, median gap 32.7 kb
against a 15–20 kb HiFi read. No linking rule can repair a junction no read
crosses. The flags are kept (behaviour-preserving at their defaults, and the
window is a prerequisite for any future graph-based linking) but they are not
the fix, and `--link-by-alleles` should stay off.

### 7.4 The unphased regions are het-RICH, which relocates the problem

Recomputing the gaps properly — merge nested/overlapping blocks, take the
complement:

| | |
|---|---|
| merged blocks | 365, covering 46.5 Mb (70.2%) |
| genuine gaps | 364, totalling **19.6 Mb** (29.8%) |
| truth hets inside them | **30,035**, i.e. **1,530/Mb** |
| chr20 average het density | ~852/Mb |

The regions we fail to phase carry **1.8x the chromosome's average
heterozygosity**. This reverses §5.3 (whose gap arithmetic did not merge nested
blocks) and it is the single most important correction in this document: we are
not failing on barren stretches, we are failing on the most polymorphic ones —
precisely where a pangenome ought to be strongest, and precisely where 30,035
available het sites would subdivide unspannable 33 kb jumps into read-spannable
steps.

Where those 30,035 sites go:

| | count | share |
|---|---|---|
| never emitted as candidates | 19,833 | 66% |
| emitted but filtered | 10,202 | 34% |

and the filtered ones, by reason:

| filter reason | sites in gaps | true hets | precision |
|---|---|---|---|
| `ref_only` | 507,248 | 9,363 | 1.85% |
| `high_af` | 14,708 | 420 | 2.86% |
| `low_af` | 1,688 | 78 | 4.62% |
| `low_depth` | 1,527 | 341 | 22.33% |

**These per-reason counts are wrong and §8 replaces them.** They are inflated by
a reporting artifact: chunks overlap, so a site near a boundary is built twice
and sees reads in only one copy, and the empty copy was also labelled
`ref_only`. 860,440 of the 969,625 zero-coverage `ref_only` site_ids also appear
elsewhere *with* coverage. The table above therefore counts phantoms, not losses.
§8 fixes the label and redoes the attribution.

---

## 8. Where the missing sites actually go

§7.4 blamed the filters. Chasing a single site through the pipeline showed that
attribution was built on a reporting artifact, and the corrected answer points
somewhere else entirely.

### 8.1 The trace

`chr20:863406` (LV1, truth het C>A) sits in a gap and was dumped as `ref_only`
with `REF_COV = 0, TOTAL_COV = 0`. Every step was checked:

| suspected cause | measured | verdict |
|---|---|---|
| parent gating drops nested observations | 0 gated observations chromosome-wide — the catalog has **no `PA` field**, only `LV`/`PS`/`AT`, so `conditional_parent_alleles` is always empty | not it (dead code here) |
| GAF rows dropped (unknown site / allele out of range) | 57,810,015 rows, **100% kept** | not it |
| per-read conflict sentinel (tandem repeats re-entering a snarl) | 57,810,015 observations, **0 conflicted** | not it |
| low graph MAPQ | all 17 reads at the site are MAPQ 58–60, default `-q` is 30 | not it |
| walk matching fails | `PGPHASE_DEBUG_SITE=863406` shows all 17 reads matched, 11 ref + 6 alt, orientation handled | not it |

The matcher works. What the compaction trace then showed:

```
[compact 863406] min_alt_depth=2 counts before=[12,5,0] after=[12,5]   <- real evidence
[compact 863406] min_alt_depth=2 counts before=[0,0,0]  after=[0]      <- empty duplicate
```

**Chunks overlap, so the site is built in two chunks and sees reads in only
one.** The empty copy was labelled `ref_only` with zero coverage — and that
label, not any real filter, is what dominated the §7.4 table. Across chr20,
**860,440 of 969,625** zero-coverage `ref_only` site_ids also appear elsewhere
with coverage.

### 8.2 Fix: name the artifact

`src/graph_bam_adapter.cpp` now distinguishes the two cases when a site has no
surviving alt:

```cpp
const char* reason = (ac[0] == 0) ? "no_reads_in_chunk" : "ref_only";
```

`ref_only` now means what it says (ref reads present, no alt cleared
`min_alt_depth`); `no_reads_in_chunk` marks the chunking artifact so the dump
can be deduplicated. Reason counts before and after:

| reason | before | after |
|---|---|---|
| `ref_only` | 1,829,985 | 826,647 |
| `no_reads_in_chunk` | — | 1,003,338 |
| `high_af` | 42,755 | 42,755 |
| `low_af` | 7,046 | 7,046 |
| `low_depth` | 3,371 | 3,371 |

This is diagnostics only — no phasing behaviour changes — but it is the
difference between a dump that points at the filters and one that points at the
cause.

### 8.3 The corrected attribution

Deduplicating by site id and keeping the most informative copy, for the 30,035
truth hets inside genuine gaps:

| reason (deduplicated) | sites in gaps | true hets | precision |
|---|---|---|---|
| `ref_only` | 218,502 | 2,027 | 0.93% |
| `no_reads_in_chunk` | 29,705 | 2,735 | 9.21% |
| `high_af` | 14,703 | 420 | 2.86% |
| `low_af` | 1,688 | 78 | 4.62% |
| `low_depth` | 1,526 | 341 | 22.35% |

| outcome | hets | share |
|---|---|---|
| **never emitted as a candidate at all** | **24,434** | **81%** |
| attributable to a real filter | 2,866 | 10% |
| lost only to the chunking artifact | 2,735 | 9% |

And the sites that were never emitted are not a position-representation
mismatch. Matching truth hets against *any* catalog record with a tolerance
window:

| tolerance | truth hets in gaps near a catalog record |
|---|---|
| exact | 5,423 / 30,035 (18.1%) |
| ±5 bp | 9,405 (31.3%) |
| ±25 bp | 14,142 (47.1%) |

Even allowing 25 bp of slack, **53% of the het sites we fail to phase have no
corresponding record in the snarl catalog at all.** The pangenome, as
deconstructed, carries no bubble there for this sample's variation.

### 8.4 What this means

The chain of causation, corrected end to end:

1. Contiguity is limited by unspannable gaps (§7.3) — 84% of breaks have no read
   crossing them.
2. Those gaps are het-**rich**, 1.8x chromosome average (§7.4).
3. What would subdivide them is the 30,035 het sites inside them.
4. **81% of those were never candidates**, and over half are absent from the
   catalog entirely — not filtered, not mis-scored, not lost in phasing.

So the binding constraint on this pipeline is **snarl-catalog completeness**,
not the phaser, not the filters, and not the linking rule. Every phasing-side
lever tried here (§5.2, §7.1) moved contiguity by a few percent at best, which
is consistent: they were all operating on the 19% of the problem.

This also retrospectively justifies the pantree direction
(`docs/pantree_catalog_experiment.md`): a denser, more accurate catalog is the
lever with real headroom. The two recoverable engineering items are small by
comparison — the 9% lost to the chunking artifact, and `low_depth` at 22%
precision, the only filter whose losses are worth revisiting.

---

## 9. One fix landed: the anchor AF gate was also suppressing output

Tracing `chr20:863406` to the end found a discrepancy between what the code says
it does and what it does.

### 9.1 The bug

`--anchor-af-margin` exists to stop paralog-collapsed sites (AF near 0.25/0.75)
from anchoring k-means, and its comment in `src/graph_bam_adapter.cpp` states the
intent plainly:

> Only lcd_var_i_to_cate changes: the site is still emitted as a call, it just
> stops voting.

It was not emitted. Failing the margin set `lcd_var_i_to_cate =
kLongcalldLowAfVar`, which sits outside `kCandGermlineClean` and outside the
k-means target mask, so the site received no consensus alleles and no phase set
and was dropped from the VCF entirely. At margin 0.12 this silently removed
**4,922 het records** (50,950 -> 55,872 with the gate off) — a completeness cost
paid for an accuracy mechanism that only needed to change voting.

`chr20:863406` is a clean example: 17 reads, 12 ref / 5 alt, AF 0.294, MAPQ
58–60, walk matching correct. It was discarded solely for being off-centre.

### 9.2 The change

A new category expresses the documented intent:

```cpp
constexpr uint32_t kCandNonAnchorHet = 0x1000u;   // callable, phaseable, never anchors
constexpr uint32_t kCandAnchorClean  = kCandCleanHetSnp | kCandCleanHetIndel | kCandCleanHom;
constexpr uint32_t kCandGermlineClean = kCandAnchorClean | kCandNonAnchorHet;
```

Non-anchor hets take part in consensus and phase-set assignment but return
weight 0 from `read_to_cons_allele_score` and `phase_matrix_var_weight`, so they
can never influence which haplotype a read is assigned to. Gated behind
`--emit-nonanchor-hets` (default off). All six unit-test binaries pass and the
default path is unchanged.

### 9.3 Measured, chr20 native call set

| config | assessed | blocks | % chr20 | N50 kb | NGC50 | switchflips | Hamming |
|---|---|---|---|---|---|---|---|
| margin 0.12 (current default) | 47,895 | 366 | 70.2 | 217 | 198 | **16** | **19** |
| margin 0.50 (gate fully off) | 51,492 | 368 | 71.6 | **280** | **246** | 26 | 49 |
| margin 0.12 + `--emit-nonanchor-hets` | 51,407 | **352** | 71.0 | 260 | 214 | 27 | 168 |

Both routes recover the same ~3,500 sites (+7.3%). The honest reading is that
**the principled fix is not the better one**: letting the off-centre sites
anchor as well (margin 0.50) gives a better Hamming (49 vs 168) than admitting
them as non-voting passengers. Their consensus is established but never
cross-checked by read scoring, so their genotypes are less reliable.

So `--emit-nonanchor-hets` ships off, and the operating-point question is simply
where to set `--anchor-af-margin`. Against the alternatives that matter — every
competitor's Hamming on its own call set is 605 (LongPhase) to 3,196 (whatshap),
and longcallD's is 884 — a move from 19 to 49 buys +24% NGC50 while staying an
order of magnitude ahead. That is likely the right trade for a completeness
claim, but it should be confirmed on the shared call set (§5.6) before the
default changes.

### 9.4 What this does not fix

+7.3% of sites. The §8.3 breakdown still stands: 81% of the het sites inside our
gaps were never candidates, and over half have no catalog record within 25 bp.
This fix recovers part of the 10% "filtered" slice, not the 81%.

---

## 10. Exactly which het sites competitors use that we miss

Direct question, directly measured: take the truth het sites inside *our* phase-
block gaps, and ask which of them each phaser actually phased.

### 10.1 Everyone struggles here; they just struggle less

| phaser | truth hets it phased inside our gaps | share of the 30,035 |
|---|---|---|
| HiPhase | 2,833 | 9.4% |
| whatshap | 2,718 | 9.0% |
| longcallD | 2,504 | 8.3% |
| LongPhase | 2,216 | 7.4% |
| **ours** | **328** | **1.1%** |

No tool phases most of what is there — HiPhase manages under 10%. But it manages
8.6x more than we do, and that is what lets it bridge junctions we break at.

**4,061** of these sites are phased by at least one competitor; **3,736** of those
are phased by none of ours. That set is the target.

### 10.2 What those 3,736 sites are

| | |
|---|---|
| SNVs | 2,855 (76.4%) |
| deletions | 452 (12.1%) |
| insertions | 425 (11.4%) |
| homopolymer run < 4 bp | 3,389 (90.7%) |
| short tandem repeat | 68 (1.8%) |

This is the opposite of the expected profile. They are **not** indels in repeats —
they are ordinary SNVs in clean sequence. Whatever is losing them is not a
repeat-resolution problem, and `--gaf-pad`/read-level fixes will not touch them.
All 3,736 are absent from our VCF entirely: not mis-genotyped, not left
unphased — never emitted.

### 10.3 Why each one is lost

Matching each site to the catalog and to the filtered-site dump with a ±25 bp
window (exact-position matching is unsafe: dump POS is the *snarl* position):

| cause | n | share | median REF / ALT cov |
|---|---|---|---|
| `ref_only` — site covered, every read matched the REF walk | 1,216 | 32.5% | 21 / 0 |
| no catalog record within 25 bp | 1,079 | 28.9% | — |
| `no_reads_in_chunk` — chunking artifact | 722 | 19.3% | 0 / 0 |
| `high_af` — AF = 1.0, graph REF allele not carried by this sample | 428 | 11.5% | 0 / 22 |
| `low_depth` | 171 | 4.6% | 0 / 2 |
| `low_af` | 120 | 3.2% | 31 / 3 |

Nothing falls outside these buckets — there is no residual "passed all filters
but silently vanished" class once positions are matched correctly.

The `ref_only` row is the important one. **ALT_COV is 0 for all 1,216, and
REF_COV >= 10 for 63% of them.** These sites are well covered and *every* read
matched the reference walk. That cannot mean the sample is homozygous — truth
says het and four other tools phase it. It means **the alternate haplotype's walk
is not among the snarl's enumerated alleles**. The catalog has the site but not
the allele.

Together with the 1,079 sites that have no catalog record at all, **2,295 of
3,736 (61%) are catalog incompleteness** — 29% at site level, 32% at allele
level. Add `high_af` (428, where the graph's REF allele is one this sample does
not carry, the known alt-vs-alt representation case) and it is **73%**.

### 10.4 What was tried and what it recovers

| fix | recovers of the 3,736 |
|---|---|
| `--anchor-af-margin 0.50` | 312 (8.4%) |
| `--emit-nonanchor-hets` | 236 (6.3%) |
| `--gaf-pad 50000` | ~0 |

`--gaf-pad` widens the per-chunk GAF read query so a snarl near a chunk edge
still sees reads. It is a correct thing to have, but it is not the cause here:
`no_reads_in_chunk` fell only 1,003,338 -> 1,002,736 and the phasing metrics are
unchanged (NGC50 198, Hamming 19 -> 21). Those sites are long snarl records that
tabix returns into chunks far from where their nodes actually sit, not reads
missed at a boundary. Ships default 0.

### 10.5 The recovery plan, in priority order

| target | sites | share | what it needs |
|---|---|---|---|
| missing **alleles** in existing snarls | 1,216 | 32.5% | catalog rebuild, or genotype reads against a locally assembled alt |
| missing **sites** entirely | 1,079 | 28.9% | catalog rebuild (pantree) or local reassembly in gaps |
| alt-vs-alt (`high_af`) | 428 | 11.5% | representation fix — `--af-vs-site-depth` exists but degrades phasing; needs rework |
| chunk-placement artifact | 722 | 19.3% | place a site in the chunk containing its *nodes*, not its VCF POS |
| thin/under-observed | 291 | 7.8% | depth/AF thresholds — smallest and least reliable slice |

So "recover all the het sites we miss" is **61–73% a catalog problem**, and the
remaining ~27% splits into one genuine engineering bug (chunk placement, 19%)
and a thin tail. No phasing-side change reaches the majority of it, which is
consistent with every phasing lever tried in §5.2, §7.1 and §9 moving the
needle by single-digit percentages.

---

## 11. Are the missed sites in the graph at all? Mostly yes — we discard them

§10 concluded "catalog incompleteness". Testing that claim directly shows it was
half wrong, and the correction is more actionable.

### 11.1 The walk matcher is not the problem

Instrumenting every read/site traversal (`[walk-match]`, `--verbose 1`):

| | count | share of spanned |
|---|---|---|
| spanned a snarl and matched an allele | 57,810,015 | **100.0%** |
| spanned but matched no enumerated allele | 11,817 | 0.02% |
| touched a boundary without spanning | 509,540 | — |

So "the alternate haplotype's walk is absent from the catalog" is **not** a
general phenomenon: essentially every read that crosses a snarl matches one of
its enumerated alleles. §10.3's interpretation of the `ref_only` bucket was
wrong.

### 11.2 Where the other half of the reads actually go

`chr20:145067` (truth het G>C) is representative. The catalog carries it at
**145061** — a 6 bp offset inside a TG microsatellite — as an LV1 snarl
`>114853298>114853301` nested inside a parent repeat snarl at 145056 with **6
alt alleles**.

Tracing it: 51 reads span the child snarl and **all 51 match allele 0**; none
matches allele 1. BAM depth at the site is 91. The missing 40 reads never span
the child at all — they traverse a *different allele of the parent*, so the
child's boundary nodes are not on their path.

The site is genuinely hom-ref *conditional on that parent allele*. The
heterozygosity is real but it lives one level up, in **which parent allele each
read carries**. We record that allele identity — and then discard it.

### 11.3 Biallelic decomposition is the loss channel

`build_graph_chunk` Phase 2 decomposes every snarl into ref-vs-alt_i pairs
"mirroring the BAM path where every candidate is a single ref/alt pair", with:

> A read observing ref contributes allele 0 to every pair from that site; a read
> observing alt_i contributes allele 1 to only the pair for alt_i.

A read carrying alt_1 therefore contributes **nothing** to the alt_2 pair. When
the sample is alt_1/alt_2 het — carrying no reference allele at all — every pair
sees only its own alt and reads AF = 1.0, which is exactly the `high_af` bucket
(428 missed sites, median REF_COV 0 / ALT_COV 22). The locus is cleanly
heterozygous in the graph and unusable after decomposition.

The enrichment confirms the mechanism:

| | catalog | missed sites |
|---|---|---|
| multi-allelic snarls (>1 ALT) | 80,373 of 977,275 (**8.2%**) | 1,256 of 3,736 (**33.6%**) |

**Multi-allelic snarls are 8.2% of the catalog but account for 33.6% of what we
miss — a 4x enrichment.** Of those, 683 sit under snarls with 4+ alts, and 597
are nested (LV >= 1), matching the §11.2 pattern.

### 11.4 Corrected accounting, and what to build

Of the 3,736 sites competitors phase and we do not, by richest catalog record
within 25 bp:

| | sites | share | recoverable from data we already have? |
|---|---|---|---|
| multi-allelic snarl nearby | 1,256 | 33.6% | **yes** — allele identity is recorded, then decomposed away |
| only biallelic records nearby | 1,401 | 37.5% | partly — includes offset representation (§11.2) and the anchor gate (§9) |
| no catalog record | 1,079 | 28.9% | no — genuine catalog gap |

So roughly **70% of the sites we miss are represented in the graph**, and the
dominant reason we do not use them is that the pipeline flattens a multi-allelic,
hierarchically nested structure into independent biallelic ref/alt pairs to match
the BAM code path. The pangenome's richest signal — *which* allele of *which*
snarl a read carries — is computed correctly (§11.1) and then thrown away.

**The fix is to phase on snarl allele identity.** The matcher already returns an
exact allele index per read per snarl with a 100% match rate. k-means needs to
treat a snarl as a categorical multi-allelic locus (read carries allele k of n)
rather than as a set of independent binary ref/alt sites. Two reads agree if they
carry the same allele index, disagree otherwise — which is well defined whether
or not either allele is the reference, and needs no AF-vs-reference reasoning at
all. That single change addresses the `high_af` bucket, the alt-vs-alt case that
`--af-vs-site-depth` failed to fix, and the nested-divergence case in §11.2
together.

That is a real change to the phasing core and is not attempted here. This
section is the evidence that it is the right one.

---

## 12. Implementing the §11 fix: two attempts, one honest failure

§11 argued the fix was to phase on snarl allele identity rather than on
decomposed biallelic pairs. Both forms of that were implemented and measured.

### 12.1 `--snarl-allele-phasing` — alt vs *other*, not alt vs ref

For a multi-allelic snarl, pair_i now contrasts "carries alt_i" against "carries
any other allele", instead of against the graph reference alone. A read on alt_j
appears in pair_i as allele 0 rather than being omitted — it is evidence
*against* alt_i, not a missing observation. Phase 2 recomputes the pair's
reference count as `site_total - alt_c`; Phase 3 emits the allele-0 observations
that were previously dropped.

This is what `--af-vs-site-depth` should have been: that flag fixed only the AF
denominator and left each pair's read profile monomorphic, which is why it
corrected the allele fractions and still degraded phasing.

It works mechanically — `high_af` falls 42,806 -> 35,173, `low_af` rises
7,053 -> 11,818 — and changes almost nothing:

| config | assessed | blocks | % chr20 | NGC50 | switchflips | Hamming |
|---|---|---|---|---|---|---|
| baseline | 47,895 | 366 | 70.2 | 198 | 16 | 19 |
| `--snarl-allele-phasing` | 47,986 | 375 | 70.8 | 198 | 17 | 20 |

**9** of the 3,736 target sites recovered.

### 12.2 Why: the missed sites live in high-multiplicity snarls

The sites that moved out of `high_af` landed in `low_af` with median AF 0.103
(median ALT 6, REF 58). They are not balanced alt_1/alt_2 heterozygotes. Allele
multiplicity of the richest snarl covering each missed site:

| alts | missed sites | whole catalog |
|---|---|---|
| 1 (biallelic) | 52.7% | 91.8% |
| 2–3 | 21.6% | 6.2% |
| 4–7 | 9.7% | 1.2% |
| **8+** | **16.0%** | **0.8%** |

**Snarls with 8+ alts are 0.8% of the catalog and 16% of what we miss — a 20x
enrichment**; 4+ alts is 2.0% vs 26%, a 13x enrichment. In a repeat snarl with
eight alleles, no single alt reaches an informative allele fraction, so *no*
binary contrast recovers the site. The representation fix is correct and
insufficient.

### 12.3 `--snarl-keep-whole` — one n-allelic anchor per snarl

The k-means core turns out to already support this: `hap_to_cons_alle[hap]` is
an integer allele index, not a flag, and Phase 4 projects "any non-zero allele
index" as ALT, so a snarl with n alleles can be a single anchor where two reads
agree iff they carry the same allele. Multi-allelic snarls are therefore kept
whole, with observations carrying their allele index through Phase 3 unchanged.

| config | assessed | blocks | % chr20 | N50 | NGC50 | switchflips | Hamming |
|---|---|---|---|---|---|---|---|
| baseline | 47,895 | 366 | 70.2 | 217 | 198 | **16** | **19** |
| `--snarl-keep-whole` | 48,468 | **360** | **74.0** | **291** | **230** | 39 | 807 |

Coverage rises 70.2% -> 74.0% — the largest single gain anything has produced —
and NGC50 198 -> 230. But Hamming goes to 807.

That could have been an artifact of emitting multi-allelic records into a VCF
whose writer picks one ALT, so it was checked read-level against the assembly
truth BAM, which never touches our VCF representation:

| config | reads evaluated | discordant | hamming rate | phase sets |
|---|---|---|---|---|
| baseline | 175,843 | **248** | **0.001410** | 279 |
| `--snarl-keep-whole` | 185,034 | 6,322 | 0.034167 | 324 |

**25x worse, on 9,191 more reads.** Not an artifact — the approach genuinely
mis-assigns reads.

### 12.4 Why exact allele identity fails, and what would work

The premise "two reads agree iff they carry the same allele index" assumes reads
from one haplotype traverse a snarl identically. In the high-multiplicity repeat
snarls that dominate the missed sites, they do not: sequencing and alignment
noise scatter same-haplotype reads across several near-identical walks differing
by a repeat unit. Exact index equality then reads as disagreement, the consensus
never stabilises, and reads are confidently placed on the wrong haplotype — which
is exactly the 25x Hamming regression.

So the missing ingredient is **allele clustering**: group a snarl's walks by
sequence similarity into (ideally two) haplotype classes before using them as a
phasing signal, so a one-repeat-unit difference does not count as a different
haplotype. That is local haplotype resolution inside the snarl — closer to a POA
over the reads' walks than to a lookup — and is the genuine next step.

Both flags ship **off**. `--snarl-allele-phasing` is a correct representation fix
with a negligible effect; `--snarl-keep-whole` is a measured negative result that
localises the real problem precisely. Defaults are untouched and all six unit
test binaries pass.
