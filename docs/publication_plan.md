# Publication plan: pangenome-graph read phasing and its effect on variant calling

## Thesis

Two claims, each standing on its own, stronger together:

1. **Phasing reads in graph space beats phasing them in linear space.** A
   graph-native phaser, given only a pangenome alignment, produces
   substantially more accurate read haplotype assignments than linear
   (longcallD/pgphase-BAM), than the graph-augmented hybrid, and than
   DeepVariant's internal phaser.
2. **Those haplotype tags improve downstream small-variant calling.** Applying
   them to the BAM DeepVariant consumes improves F1 over the published
   baselines, on the tool authors' own case studies, without retraining.

The second claim is what makes the first matter to a reader who does not care
about phasing per se.

---

## Contribution 1 — graph-native read phasing

### C1.1 Method

Two gates, both cheap, each fixing a distinct and diagnosable failure:

- `--anchor-af-margin` — a site votes in k-means only when |AF - 0.5| <= F.
  Where the graph collapses paralogous loci, a site het on one copy and hom on
  the other sits near AF 0.25/0.75, clears the 0.20/0.80 depth filters, and then
  votes as a haplotype marker, so the pipeline separates *paralogs* rather than
  haplotypes -- confidently, which is why it is not caught by evidence gates.
- `--min-read-margin` — a read is committed only when its clean-SNP margin
  reaches a threshold. `init_assign_read_hap` otherwise commits on any non-zero
  score, and one-SNP-margin reads are 11% of reads but 69% of errors.

### C1.2 Headline benchmark (have)

chr20, identical truth BAM, identical harness, exclusions applied:

| pipeline | discordant | hamming | vs longcallD |
|---|---|---|---|
| longcallD (default) | 4,636 | 0.021662 | -- |
| pgphase BAM | 4,648 | 0.021638 | 1.0x |
| pgphase hybrid | 1,114 | 0.005282 | 4.2x |
| graph, margin 1 | 729 | 0.003670 | 6.4x |
| graph, margin 2 | 237 | 0.001353 | **19.6x** |

pgphase-BAM landing on longcallD's numbers is the control that shows the
harness is not flattering us.

### C1.3 Generalisation (partly have)

| | status |
|---|---|
| chr20, chr18, chr12 vs hybrid | **have** (4.7x / 5.4x / 13.2x at margin 2) |
| chr1 | truth BAM build was stopped; ~3 h to finish |
| chr19 | inputs ~70% generated; catalog + GAF remain |
| second sample (HG003, HG001) | **to do** -- the single biggest gap |
| ONT | **to do** |

One sample is the weakest point of C1. HG003 is already needed for C2, so it
comes almost free.

### C1.4 Ablations (have)

Each gate alone and together; `-q` sweep (the strongest single filter -- removing
it costs 3.6x more discordant reads for 3.1% more reads); the coverage/accuracy
frontier; and the operating-point table.

### C1.5 Limits, stated honestly (have)

A reviewer will ask where it fails, and we can answer precisely:

- **Coverage is heterozygosity-limited, not method-limited.** Of unphased reads,
  10.2% are observed at a median of **2** candidate sites -- and the BAM pipeline
  has exactly 2 there as well. Nothing is hiding. Hybrid's extra coverage is it
  committing anyway, at 7.4% error against its own 1.8%.
- **Runs of homozygosity.** chr20 44--45 Mb has 3 het SNPs against a chromosome
  median of 852/Mb; graph, BAM and hybrid independently agree (3/4/5).
- **Graph-blind satellite.** chr20 27.2--28.8 Mb has zero catalog sites over
  1.6 Mb because the pangenome carries no alternative haplotypes there (`AN`
  collapses ~450 -> 7 -> 3).
- **Alt-vs-alt hets are invisible.** 34,698 chr18 sites with `REF_COV = 0`, median
  alt depth 66, are scored AF = 1.0 and dropped. `--af-vs-site-depth` fixes the
  representation and *degrades* phasing, so it ships off -- an honest negative.

---

## Contribution 2 — effect on DeepVariant

The design principle throughout: **change exactly one thing per arm.**

### C2.1 Headline: PacBio model, HG003 chr20, GRCh38

Their case study, their data, their published numbers as the baseline.

- BAM: `GRCh38.m84039_241002_000337_s3.hifi_reads.bc2020.bam` (Revio SPRQ, 32x)
- Truth: GIAB HG003 v4.2.1 GRCh38; reference: GRCh38 no-alt analysis set
- Published: SNP F1 **0.999067**, INDEL F1 **0.993679**

Arms:

| arm | BAM | DV phasing | isolates |
|---|---|---|---|
| A vanilla | pbmm2 | internal | published baseline |
| B refined | pbmm2 + refined CIGARs | internal | refinement alone |
| C gHP | pbmm2 | **ours** (`phase_reads=false`) | phasing alone |
| D refined+gHP | pbmm2 + refined | **ours** | the combination |

The 2x2 matters: the existing Slack result reports D vs A, which cannot say
whether the gain is refinement, phasing, or their interaction. B and C separate
them and cost one DV run each.

### C2.2 Illumina: the vg case study

`hprc-v1.1-mc-grch38` + giraffe, WGS model. The case study uses **no phasing at
all**, so this is a clean addition rather than a substitution, and it shows the
method is not HiFi-specific. Arms: published baseline vs + our HP tags.

### C2.3 Single-mapping ablation (running now)

giraffe-surjected BAM only -- no pbmm2, no refinement -- HG002 chr20 in CHM13
against GIAB v5.0q CHM13v2.0 truth. Answers "does graph HP beat DeepVariant's
internal phaser on its own", with no refinement confound, and demonstrates that
**one mapping suffices**: giraffe emits the surjected BAM and the GAF together,
so the linear alignment DeepVariant consumes and the graph alignment we phase on
share read names by construction.

### C2.4 How much phasing accuracy can DeepVariant exploit?

The most interesting scientific question here, and one nobody has answered.
Margin 2 is our best phasing but tags only 82% of reads; margin 0 tags 95% and
still beats hybrid. Sweep margin 0/1/2/3 against **F1** rather than hamming.
Plausible outcome: F1 peaks at a *worse*-phasing, higher-coverage setting -- a
finding worth reporting either way.

### C2.5 Where phasing should matter most

Stratify hap.py by the GIAB stratifications: segmental duplications, low
mappability, homopolymers, MHC. The mechanism predicts gains concentrate in
segdups and low-mappability regions. If they appear elsewhere instead, our causal
story is wrong and we should know before a reviewer says so.

### C2.6 ONT

Untested, and the prior expectation is a larger benefit (noisier reads, weaker
internal phasing). High value, and the natural "future work" if time is short.

---

## Baselines to beat, per scenario

| scenario | baseline |
|---|---|
| phasing | longcallD, whatshap, HiPhase, LongPhase, pgphase-BAM, pgphase-hybrid |
| DV HiFi | published PacBio case-study numbers; DV internal phasing |
| DV Illumina | published vg case-study numbers |
| DV, phasing source | whatshap HP and hybrid HP, not just "no HP" |

The last row matters: "our HP beats no HP" is weak; "our HP beats whatshap HP and
hybrid HP on the same BAM" is the claim worth making.

---

## Figures

1. Phasing accuracy vs all phasers, three chromosomes (C1.2/C1.3).
2. Ablation: each gate's contribution, plus the coverage/accuracy frontier.
3. DV 2x2 on HG003 chr20 -- vanilla / refined / gHP / both (C2.1).
4. F1 vs phasing operating point (C2.4) -- the "how much can DV exploit" curve.
5. Stratified F1 by GIAB region class (C2.5).
6. Mechanism panel: the paralog-collapse locus, AF distributions, before/after.

---

## Risks

- **Margin 2 may lose on F1.** Fewer tagged reads could beat more accurate ones.
  Mitigated by C2.4 being a sweep, not a point.
- **Refinement may dominate phasing** in the D-vs-A gain. C2.1's B and C arms
  settle it; if refinement is the whole story, that changes the paper's emphasis.
- **One sample.** HG002-only phasing results are the most attackable claim.
  HG003 is needed for C2 anyway.
- **CHM13 vs GRCh38.** Phasing work is CHM13, DV case studies are GRCh38. Keep
  the two coordinate systems clearly separated in the text; do not mix them in a
  single table.
- **Gains may be small in absolute terms.** Going 0.99359 -> 0.99396 INDEL F1 is
  real but needs error bars -- report TP/FN/FP counts, and consider multiple
  chromosomes or replicates rather than chr20 alone.

---

## Status

| | |
|---|---|
| have | C1.1, C1.2, C1.4, C1.5; C1.3 for three chromosomes |
| running | C2.3 (HG002 chr20 CHM13, single mapping, 2 arms) |
| next | C2.1 -- fetch HG003 BAM + GRCh38 + v4.2.1 truth, run the 2x2 |
| then | C2.4 margin sweep, C2.5 stratification, C1.3 second sample |
| later | C2.2 Illumina, C2.6 ONT |

Blocking dependency: C2.1 needs the CIGAR-refinement step used in the existing
Slack result. Everything else here is reproducible from this repo plus public
downloads.
