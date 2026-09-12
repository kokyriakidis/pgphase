# Research plan: analyses, inputs, and what is still missing

Companion to `publication_plan.md`. That document argues the claims; this one is
the run manifest — every analysis, what it needs, what exists, and the order to
do them in.

Status legend: **[have]** ready now · **[part]** partially available ·
**[need]** missing, with the acquisition cost stated.

---

## 0. Asset inventory

### Tools

| tool | status | needed for |
|---|---|---|
| `vg` | **[have]** | giraffe, deconstruct, chunk |
| `pggaf` | **[have]** | GAF annotate/index |
| minimap2 v2.31, HipHap | **[have]** pinned submodules | truth BAMs |
| DeepVariant 1.10.0 | **[have]** docker image | all DV arms |
| samtools, docker | **[have]** | |
| **hap.py** | **[need]** `docker pull jmcdani20/hap.py:v0.3.12` | **every DV arm** |
| `pbmm2` | **[need]** conda; or use case-study BAM as-is | linear HiFi arm |
| whatshap / HiPhase / LongPhase | **[need]** `envs/*.yaml` exist, not installed | phaser baselines |
| `kmc` | **[need]** conda | giraffe haplotype sampling (Illumina) |
| **CIGAR refinement** | **[need]** — not in this repo | arms B and D |

hap.py is the single most blocking item: no DV arm can be scored without it.

### Graphs

| asset | status |
|---|---|
| `hprc-v2.1-mc-chm13-eval.gbz` + `.ri` + `.dist` | **[have]** 6.1 / 11.5 GB |
| `hprc-v2.1-mc-grch38-eval.full.gbz` | **[part]** 33.8 GB, **bare** |
| GRCh38 `.dist` / `.min` / `.hapl` | **[need]** rebuild: hours of CPU, large RAM, or copy from the cluster |

The pipeline uses **HPRC v2.1 throughout**, matching the CHM13 side, so the
v1.1 graph from the vg case study is *not* needed. One consequence: Illumina
results (2C) will not be numerically comparable to that case study's published
figures, which were produced on v1.1. Report 2C as an internal A/B (same graph,
same BAM, phasing on vs off) rather than as a delta against their table.

**This is the critical path for every GRCh38 graph arm.** Copying prebuilt
indexes from the eval cluster avoids the most expensive step in the plan.

### References and truth

| asset | status |
|---|---|
| CHM13v2.0 FASTA | **[have]** |
| GRCh38 no-alt FASTA | **[part]** downloading |
| GIAB HG002 v5.0q **CHM13** smvar VCF/BED | **[have]** |
| GIAB HG003 v4.2.1 **GRCh38** VCF/BED | **[have]** |
| GIAB HG002 v4.2.1 **GRCh38** VCF/BED | **[need]** trivial download |
| GIAB stratification BEDs v3.x | **[need]** trivial download; for stratified F1 |

### Reads and alignments

| asset | status |
|---|---|
| HG002 HiFi, CHM13 surjected BAM + GAF (chr1/12/18/20) | **[have]** |
| HG003 chr20 HiFi GRCh38 BAM (case study) | **[have]** 1.07 GB, 169,187 reads |
| HG003 GAF (graph alignment) | **[need]** requires GRCh38 giraffe indexes |
| HG002 GRCh38 BAM + GAF | **[need]** |
| Illumina (either sample) | **[need]** for the vg case study arm |
| ONT | **[need]** |

### Truth BAMs for phasing evaluation (external store)

`~/Downloads/pgphase-eval-data/truth/` — chr20, chr18, chr12 **[have]**;
chr1 **[part]** (~3 h to finish); chr19 **[need]**.

---

## 1. Phasing analyses (Contribution 1)

| # | analysis | inputs | status |
|---|---|---|---|
| P1 | graph vs BAM vs hybrid, chr20/18/12 | truth BAMs, GAFs | **[have]** |
| P2 | vs longcallD | longcallD phased BAM | **[have]** chr20 only; binary not on this box |
| P3 | vs whatshap / HiPhase / LongPhase | install from `envs/`, DV VCF as input | **[need]** install |
| P4 | gate ablations, `-q` sweep, frontier | — | **[have]** |
| P5 | **second sample (HG003)** | HG003 GAF + truth BAM | **[need]** — the weakest claim |
| P6 | chr1, chr19 | finish truth BAM / catalog | **[part]** |
| P7 | ONT | ONT reads + GAF | **[need]** |

P3 and P5 are what a reviewer will demand. P5 additionally unblocks the DV work,
so it is the highest-value missing item in the whole plan.

---

## 2. DeepVariant analyses (Contribution 2)

Design rule: **one variable per arm.** Each row is one DV run plus one hap.py run.

### 2A. Linear HiFi — HG003 chr20, GRCh38 (the headline)

Published baseline: SNP F1 **0.999067**, INDEL F1 **0.993679**.

| arm | BAM | phasing | needs | status |
|---|---|---|---|---|
| A1 vanilla | case-study BAM | DV internal | hap.py | **[have]** inputs |
| A2 refined | + refined CIGARs | DV internal | refinement step | **[need]** |
| A3 gHP | case-study BAM | **ours** | HG003 GAF | **[need]** |
| A4 refined + gHP | + refined | **ours** | both | **[need]** |

A1 is runnable now and reproduces the published number — worth doing first as a
harness check. The existing Slack result is A4 vs A1; A2 and A3 are what
attribute the gain between refinement and phasing.

### 2B. Graph-surjected HiFi — same sample, single mapping

| arm | BAM | phasing | status |
|---|---|---|---|
| B1 | giraffe-surjected | DV internal | **[need]** GRCh38 indexes |
| B2 | giraffe-surjected | **ours** | **[need]** |

Demonstrates one mapping suffices: giraffe emits the surjected BAM and the GAF
together, so read names match by construction. Comparing 2A against 2B also
answers whether pbmm2 or giraffe is the better substrate for DV.

### 2C. Illumina — the vg case study

Baseline uses **no phasing at all**, so ours is an addition, not a substitution,
and it shows the method is not HiFi-specific. Needs Illumina reads, `kmc`, and
the v2.1 GRCh38 graph + `.hapl`. **[need]**

Because we run v2.1 and the case study ran v1.1, this arm is an internal A/B
(phasing on vs off, same graph and BAM), not a delta against their published
numbers. Stating that plainly is better than an apples-to-oranges comparison.

### 2D. Single-mapping CHM13 ablation — **running now**

HG002 chr20, CHM13, GIAB v5.0q truth, giraffe-surjected BAM only. Two arms:
DV internal phasing vs our HP. First real answer on whether graph HP beats
DeepVariant's own phaser, with no refinement confound.

### 2E. Margin sweep against F1

Re-run 2A-A3 (or 2D-B2) at `--min-read-margin` 0/1/2/3. Margin 2 is the best
phasing but tags 82% of reads; margin 0 tags 95%. **F1 may peak at a
worse-phasing setting.** Cheap once one arm works — it is the same pipeline with
one flag changed, and it is the most interesting question in the paper.

### 2F. Stratified F1

hap.py with GIAB stratification BEDs: segdups, low mappability, homopolymers,
MHC. The mechanism predicts gains concentrate in segdups and low-mappability
regions; if they appear elsewhere, the causal story is wrong.

### 2G. ONT

Untested; expected to benefit more. Natural future work if time is short.

---

## 2.5 Comparability rules (apples-to-apples)

Every number in a table must differ from its neighbours in **exactly one**
dimension. The dimensions that silently break this:

| dimension | rule |
|---|---|
| reference | never mix CHM13 and GRCh38 rows |
| truth set | v4.2.1 (GRCh38) and v5.0q (CHM13) are different difficulties; deltas only within one |
| sample | HG002 and HG003 are not interchangeable |
| reads / BAM | same alignment for every arm in a table |
| DV model + version | PACBIO 1.10.0 throughout |
| **small model** | see below -- the subtle one |
| region | chr20 throughout |

### The small-model confound

`--disable_small_model` is **forced** on any arm that supplies external HP tags:
with `phase_reads=false` the small-model feature vector drops 106 -> 70 and the
shipped checkpoint raises `ValueError`. DeepVariant's published numbers were
produced with the small model **enabled**.

So "our HP vs the published baseline" confounds two changes at once. The fix is
three reference points, not two:

| arm | small model | phasing | role |
|---|---|---|---|
| A0 | **on** | DV internal | reproduces the published configuration |
| A1 | **off** | DV internal | the controlled baseline |
| B | **off** (forced) | **ours** | the treatment |

Report **B vs A1** as the result — one variable — and **A0** separately to show
the harness reproduces the published number and to state what disabling the
small model costs on its own. Quoting B against A0 would overstate or understate
the phasing effect by whatever the small model contributes.

---

## 3. Methodological findings that must be in the paper

Discovered while building this, and each would silently invalidate a naive
comparison:

1. **DeepVariant's small model is incompatible with external HP tags.** With
   `phase_reads=false` the small-model feature vector drops 106 -> 70 and the
   shipped checkpoint raises `ValueError`. **`--disable_small_model` must be set
   in both arms**, or the comparison confounds phasing with small-model on/off.
2. **`vg giraffe` emits the surjected BAM and the GAF from one mapping**, so no
   second alignment is needed and read names match by construction.
3. **Contig naming.** Surjected BAMs carry `GRCh38#0#chrN` / `CHM13#0#chrN`;
   references and truth VCFs use bare `chrN`. Strip the prefix by reheadering.
4. **CHM13 vs GRCh38 must not be mixed in one table.** Phasing work is CHM13;
   the DV case studies are GRCh38.

---

## 4. Critical path

```
hap.py image ──┬─> 2A-A1 (reproduce published baseline)      ~1 h
               └─> 2D    (CHM13 ablation, already running)

GRCh38 giraffe indexes ──> HG003 GAF ──┬─> 2A-A3, 2A-A4
   (copy from cluster, or rebuild)     ├─> 2B
                                       └─> P5 (HG003 phasing)

CIGAR refinement step ──> 2A-A2, 2A-A4
```

Three independent blockers: **hap.py** (minutes), **GRCh38 giraffe indexes**
(hours to rebuild, or a copy), **CIGAR refinement** (only you have it).

## 5. Suggested order

1. Pull hap.py; score 2D; run and score 2A-A1. *Proves the harness end to end.*
2. Obtain GRCh38 giraffe indexes. *Unblocks everything else.*
3. HG003 GAF -> 2A-A3 -> 2E margin sweep. *The headline claim.*
4. Add A2/A4 once the refinement step is available. *Completes the 2x2.*
5. P5 HG003 phasing benchmark, P3 phaser baselines. *Shores up Contribution 1.*
6. 2F stratification, then 2B / 2C / 2G as time allows.
