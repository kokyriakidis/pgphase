# How the phasing field evaluates, and how we should

Reference analysis of the HiPhase benchmark (the most recent and most directly
comparable), what it does and does not specify, and what we must add to be read
as comparable. Written after reproducing their metric on our data.

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
