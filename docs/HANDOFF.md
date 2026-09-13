# pgphase — project state, findings, and next steps

Self-contained handoff. Everything below was measured on HG002 chr20 unless
stated otherwise. Numbers that were later found wrong are marked **RETRACTED**
with the correction, because several of them circulated before being corrected.

---

## 1. What the project is

`pgphase` phases long reads using a **pangenome graph alignment (GAF) only** — no
linear-reference variant calling step. The research claim being built:

1. Phasing reads in graph space beats phasing them in linear space.
2. Those haplotype tags improve downstream DeepVariant calling.

Contribution 1 is what is currently being defended. Contribution 2 is deferred
(DeepVariant likely needs retraining on graph-phased reads to benefit properly).

Key pipelines (`./pgphase <subcommand>`):
- `collect-graph-variation` — the graph pipeline (GAF-driven). **This is the one
  under active work.**
- `collect-bam-variation`, `collect-hybrid-variation` — linear / hybrid baselines.
- `build-snarl-catalog` — builds the site catalog from a GBZ.

---

## 2. Current headline results (chr20, HG002, CHM13, GIAB v5.0q truth)

### Read-level, vs diploid-assembly truth BAM

| pipeline | discordant reads | hamming |
|---|---|---|
| longcallD (default) | 4,636 | 0.021662 |
| pgphase BAM | 4,648 | 0.021638 |
| pgphase hybrid | 1,114 | 0.005282 |
| **pgphase graph, margin 2** | **248** | **0.001410** |

### Variant-level, all phasers on the **same DeepVariant call set**

| method | assessed | blocks | % chr20 | N50 kb | NGC50 kb | switchflips | Hamming |
|---|---|---|---|---|---|---|---|
| whatshap 2.8 | 62,510 | 228 | 85.4 | 673 | 338 | 309 | 3,196 |
| HiPhase 1.6.0 | 62,265 | 194 | 87.2 | **856** | **644** | 172 | 2,279 |
| LongPhase 2.0.2 | 54,680 | 226 | 81.7 | 475 | 400 | 143 | 605 |
| graph -> DV calls (corrected) | 59,015 | 276 | — | 215 | 204 | **15** | **16** |
| hybrid -> DV calls (corrected) | 59,665 | 298 | — | 312 | 290 | 41 | 68 |
| **all-linear-site seed + conservative PS merge** | 59,447 | 247 | — | 401 | 374 | 39 | 79 |

The experimental linear-site path is still not a production command. **We win
accuracy, but still lose corrected contiguity:** its NGC50 is 374 kb against
LongPhase's 400 kb and HiPhase's 644 kb. That is the defensible claim today.

---

## 3. The critical methodological traps (read this before trusting any number)

Five separate wrong conclusions were produced by these; all are now fixed.

1. **Different call sets are not comparable.** Our native VCF has 47,895 het
   sites; competitors phase DeepVariant's ~62k. Restricting all tools to *our*
   sites showed competitors doing **better** than us (longcallD 3 switchflips vs
   our 16) — our apparent win was site selection. A size-matched *random* subset
   control is essential: random removal of 23% of sites removes ~20% of errors,
   removing *our* skipped 23% removes ~90%. **Always transfer our read HP tags
   onto the shared call set** (`scripts/phase_vcf_from_hp.py`) before comparing.
2. **HP tags are only meaningful within a phase set.** HP=1 in one block is
   unrelated to HP=1 in the next. Pooling reads by HP alone scrambles every site
   near a block boundary. `phase_vcf_from_hp.py` now groups by `(PS, HP)`.
3. **Phase blocks nest and overlap.** Gaps must be computed as the complement of
   *merged* blocks. Computing them as differences between consecutive sorted
   blocks produced a confident, wrong conclusion (see §4, RETRACTED).
4. **Position conventions differ.** The filtered-site dump's POS is the *snarl*
   position, not the variant position; catalog POS can sit tens of bp from the
   truth variant inside a repeat. Match with a ±25 bp window, never exactly.
5. **Input VCF phase must be cleared before transfer.** DeepVariant's input VCF
   already had 67,127 phased hets. The old transfer script left unsupported
   records unchanged, so old graph-to-DV results mixed inherited DeepVariant
   phase with pgphase evidence. The script now unphases every heterozygote before
   applying HP evidence. All pre-correction shared-call-set pgphase rows in the
   longer analysis are retracted.

---

## 4. The contiguity investigation (the main line of work)

**Question:** HiPhase emits 194 blocks where we emit 549 on the same call set.
Why, and can it be fixed?

### Eliminated causes
- **Chunk boundaries** — raising chunk size 20x (500 kb → 10 Mb) changes nothing
  (366 → 371 blocks). Blocks already span 2.5 default chunks.
- **Block-link threshold** — `--min-block-link-reads 1` gives +41% NGC50 for 13x
  worse Hamming. Wrong lever.
- **Adjacent-pair chain rule** — the break rule only linked consecutive het
  variants. Generalised to a window (`--block-link-window`, verified
  behaviour-preserving at 1). Moves 366 → 360 blocks. Not the cause.
- **Untagged reads as linking evidence** (`--link-by-alleles`) — lets reads link
  variants without a prior haplotype. Costs 13x Hamming for nothing.

### The actual reason
**84% of our block breaks (307/364) have no read spanning them at all**, median
gap 32.7 kb against a 15–20 kb HiFi read. No linking rule can repair a junction
no read crosses. Contiguity is limited by *which sites we phase*, not by how we
link them.

### RETRACTED: "the gaps are het-poor"
Computed without merging nested blocks. **Correction:** the 364 genuine gaps
total 19.6 Mb and contain **30,035 truth hets at 1,530/Mb — 1.8x the chromosome
average (852/Mb)**. We fail on the *most polymorphic* regions, not barren ones.

### HiPhase's bridges are correct, not reckless
Of our 466 junctions, a single HiPhase block spans 295, and its switch-error
rate inside them is **5.4%** — identical to its 5.0% baseline elsewhere.

---

## 5. Which sites we miss, and why

Truth hets inside our gaps that a competitor phases and we do not: **3,736**.
(Context: HiPhase phases only 9.4% of the 30,035 hets in our gaps, we phase 1.1%.
Everyone struggles there; they struggle 8.6x less.)

**What they are:** 76% SNVs, 91% in clean sequence (homopolymer run < 4), 1.8%
STR. Not an indel or repeat-resolution problem. All 3,736 are absent from our
VCF entirely — never emitted, not mis-genotyped.

**Why (±25 bp matching, deduplicated):**

| cause | n | share | median REF/ALT cov |
|---|---|---|---|
| `ref_only` — covered, every read matched the REF walk | 1,216 | 32.5% | 21 / 0 |
| no catalog record within 25 bp | 1,079 | 28.9% | — |
| `no_reads_in_chunk` — chunk-placement artifact | 722 | 19.3% | 0 / 0 |
| `high_af` — AF 1.0, graph REF allele not carried by sample | 428 | 11.5% | 0 / 22 |
| `low_depth` / `low_af` | 291 | 7.8% | — |

### RETRACTED: "the alternate haplotype's walk is absent from the catalog"
Instrumenting every traversal (`--verbose 1`, `[walk-match]`) shows **100.0% of
spanned traversals match an enumerated allele** (57,810,015 matched vs 11,817
unmatched). The matcher is not losing anything.

**What actually happens** (traced at `chr20:145067`): the catalog carries the
variant at 145061 — 6 bp away inside a TG microsatellite — as an LV1 snarl nested
in a parent with 6 alts. 51 reads span the child and *all 51 match allele 0*;
depth is 91. The 40 alt-carrying reads never span the child at all because they
traverse a **different allele of the parent**. The heterozygosity lives in the
parent's allele identity, one level up.

---

## 6. Multi-allelic snarls: the core structural finding

`build_graph_chunk` Phase 2 decomposes every snarl into independent ref-vs-alt_i
pairs "mirroring the BAM path". A read carrying alt_1 contributes **nothing** to
the alt_2 pair. For an alt_1/alt_2 heterozygote (no reference allele carried),
every pair reads AF = 1.0 and is filtered as `high_af`.

**Enrichment confirms this matters:**

| | catalog | missed sites |
|---|---|---|
| multi-allelic snarls (>1 ALT) | 8.2% | 33.6% |
| snarls with 8+ alts | 0.8% | **16.0%** |

8+ alt snarls are **20x enriched** among what we miss.

### Multi-allelic snarls are the *richest* phasing signal available
Measured against truth haplotype tags (`HO:Z:MAT/PAT`), sampling 2,400 snarls:

| alleles | median within-hap purity | purity ≥ 0.9 | MAT/PAT favour different alleles |
|---|---|---|---|
| 2 | 1.000 | 99.3% | **7.0%** |
| 4–7 | 0.972 | 61.3% | 33.2% |
| 8–15 | 0.885 | 48.8% | 52.7% |
| 16+ | 0.824 | 37.9% | **69.0%** |

Biallelic snarls are almost always homozygous in this sample (only 7% informative).
Multi-allelic snarls are informative 33–69% of the time. **This is where the
phasing signal is**, and the pipeline currently destroys it.

**Nesting level does NOT predict purity** (LV0 79.1%, LV1 80.5%, LV2 72.5% —
flat). So restricting to leaf-level snarls would not help; allele count and walk
length are the predictors.

---

## 7. Code changes made (all default OFF, defaults unchanged, 6/6 tests pass)

| flag | what it does | measured effect |
|---|---|---|
| `--block-link-window INT` [1] | break rule searches N preceding hets, not just the adjacent one | 366 → 360 blocks. Verified identical to old behaviour at 1. |
| `--link-by-alleles` | untagged reads link variants by allele pattern | 13x worse Hamming, no block gain |
| `--emit-nonanchor-hets` | emit/phase hets failing `--anchor-af-margin` (they still never anchor) | +3,500 sites, Hamming 19 → 168 |
| `--gaf-pad INT` [0] | widen per-chunk GAF read query | no effect (~0) |
| `--snarl-allele-phasing` | multi-allelic pairs contrast alt vs *other*, not alt vs ref | `high_af` 42,806 → 35,173; recovers 9 of 3,736 |
| `--snarl-keep-whole` | keep multi-allelic snarls as single n-allelic anchors | %chr20 70.2 → 74.0, NGC50 198 → 230, **but read discordance 248 → 6,322** |
| `--snarl-top2-frac F` [0.9] | require F of a snarl's reads on its top 2 alleles | see below |

**Bug fixes landed (affect default behaviour):**
- `scripts/phase_vcf_from_hp.py` — groups reads by `(PS, HP)` instead of `HP`.
  Pooling across phase sets was simply incorrect.
- Filtered-site dump now distinguishes `no_reads_in_chunk` (chunk artifact) from
  `ref_only` (real). Before: 1,829,985 "ref_only", of which 1,003,338 were
  phantoms — 860,440 of 969,625 zero-coverage site_ids also appeared elsewhere
  *with* coverage. This one relabel invalidated an entire earlier analysis.
- New diagnostics at `--verbose 1`: `[walk-match]`, `[rows]`, `[dedup]`,
  `[parent-gate]`. Plus `PGPHASE_DEBUG_SITE=<pos>` for per-read site tracing.

---

## 8. Whole-snarl phasing: regression mostly explained, still not default

`--snarl-keep-whole` should work — the k-means core is already allele-index
agnostic (`hap_to_cons_alle[hap]` is an integer; Phase 4 projects any non-zero
index as ALT). The previous catastrophic regression was partly real and partly
two implementation bugs:

1. After Phase 2 rebuilt `allele_counts` into candidate space, Phase 3 used that
   rewritten vector to decide whether the *original* source snarl was
   multi-allelic. For decomposed multi-allelic sites, reads on alt_j could fail to
   appear as allele 0 evidence against alt_i.
2. Whole n-allelic candidates were classified by collapsed non-ref AF
   (`alt_total / site_total`). An alt_1/alt_2 heterozygote with zero ref reads was
   therefore marked `CleanHom`, even though the whole-snarl representation had
   the two allele indices needed to phase it.

Both are fixed. Multi-allelic source identity is now preserved through the
candidate rewrite, n-allelic whole snarls classify hom/het by top allele fraction,
and n-allelic anchors get single weight rather than clean-SNP/indel double weight.

Fresh chr20 HG002 read-level runs after the fix:

| config | candidates | phase sets | reads evaluated | discordant | read-hamming |
|---|---:|---:|---:|---:|---:|
| baseline | 73,627 | 279 | 175,843 | **248** | **0.001410** |
| `--snarl-allele-phasing` | 77,570 | 282 | 176,196 | 518 | 0.002940 |
| `--snarl-keep-whole` (`top2 >= 0.90`) | 80,837 | 281 | 175,792 | 393 | 0.002236 |
| `--snarl-keep-whole --snarl-top2-frac 0.95` | 79,270 | 281 | 175,828 | 395 | 0.002247 |

So whole-snarl phasing is no longer a 25x accuracy failure; it is now near the
default, but still worse than baseline on read accuracy and not worth enabling
by default. The remaining likely issue is still allele clustering inside
high-multiplicity repeat snarls: exact allele-index equality is too brittle when
same-haplotype reads scatter across near-identical walks.

### Competitor-site gap diagnosis, refreshed with SITE_ID accounting

Added `--phase-sites-out FILE` to `collect-graph-variation`. It writes retained
graph-site candidates with source `SITE_ID`, allele counts, phase set, and hap
alleles. This matters because the main candidate TSV is normalized to variant
coordinates and does not preserve the original graph catalog ID; without this,
retained rows and filtered rows cannot be audited by the same key.

Also added `scripts/analyze_graph_gap_site_loss.py`, which classifies truth hets
outside pgphase merged phase blocks. With competitor VCFs supplied, it can target
truth sites near phased HiPhase/LongPhase/WhatsHap variants. The durable-ish
scratchpad competitor VCFs currently used are:

```
/tmp/claude-1000/-home-kokyriakidis-Downloads-pgphase/68d70dd6-387e-4e2f-885c-9efbd718ca83/scratchpad/phasers/chr20/
  hp.vcf.gz
  lp.vcf.gz
  ws.vcf.gz
```

For the exact-position competitor-site target (stricter than ±25 bp; closest to
the old 3,736-site handoff number), baseline now reports 3,197 truth hets phased
by at least one competitor and outside pgphase blocks:

| reason | baseline | `--snarl-allele-phasing` |
|---|---:|---:|
| no catalog record within ±25 bp | 1,012 | 1,012 |
| nearby catalog site retained, but unphased | 606 | 709 |
| `ref_only` | 753 | 756 |
| `no_reads_in_chunk` | 386 | 388 |
| `high_af` | 268 | 119 |
| `low_depth` | 143 | 144 |
| `low_af` | 17 | 22 |
| nearby catalog site already phased | 12 | 21 |

For the more forgiving ±25 bp competitor-site target, baseline has 5,595 misses:
1,915 no-catalog, 1,198 `ref_only`, 841 `no_reads_in_chunk`, 818 retained but
unphased, 499 `high_af`, 278 `low_depth`, 24 `low_af`, and 22 already phased
nearby. `--snarl-allele-phasing` cuts `high_af` to 236 but mostly converts those
sites into retained-but-unphased rows (1,028), so the next bottleneck is not site
discovery; it is making retained multi-allelic/repeat sites reliable anchors or
bridges.

Identifier-level catalog audit on the baseline diagnostic run:

- chr20 catalog records: 977,275, duplicate IDs: 0.
- Retained graph-site parent IDs: 70,581.
- Filtered graph-site parent IDs: 976,572.
- Only 3 catalog IDs had neither a retained nor filtered row.

So the pipeline is accounting for essentially all catalog records. The remaining
question is catalog completeness relative to the original GBZ/pantree graph:
about one third of competitor-site exact misses still have no catalog record
near the truth variant.

---

## Graph-locked native-BAM validation on chr12 and chr18

The direct chr20 recipe was tested on the available chr18 and chr12 graph GAFs,
surjected BAMs, catalogs, references, and diplinator truth. Native pgphase BAM
calls are much less selective than the chr20 DeepVariant source: chr18 yielded
12,007 exact-private GQ10 gap sites versus 958 on chr20. Direct joint phasing was
unsafe even after read-overlap filtering (565 sites, 3,141/226,453 discordant)
and after restricting to 114 clean balanced SNPs (3,149/226,080 discordant).
The errors were concentrated in two 3,000-plus-read blocks, so PS50 could not
remove them.

The selected fix is `scripts/merge_graph_hybrid_tags.py`. It treats joint
hybrid phasing as a proposal: graph-tagged reads retain their graph HP/PS
exactly, and graph-unphased reads are added only when at least ten shared reads
orient the hybrid block to one graph block with margin five, 90% purity, and
both haplotypes represented. This is truth-free and rejects the chr18 mixed
blocks.

| chromosome | graph reads/errors | graph-lock reads/errors | net reads/errors | Hamming |
|---|---:|---:|---:|---:|
| chr18 | 224,746 / 388 | 228,020 / 405 | +3,274 / +17 | 0.18% |
| chr12 | 410,699 / 257 | 415,431 / 286 | +4,732 / +29 | 0.07% |

For native BAM VCFs, extract with `--clean-snps-only`,
`--exclude-graph-positions`, `--min-vaf 0.30`, and `--max-vaf 0.70` plus the
existing BAM bridge options. This retained 111 sites on chr18 and 131 on chr12.
Direct joint output remains experimental; graph-lock output is the current
robust policy.
It adds reads to existing graph phase sets but deliberately does not merge
independent graph PS labels yet.

```bash
python3 scripts/merge_graph_hybrid_tags.py \
  --graph-bam graph.bam --hybrid-bam hybrid.bam --output locked.bam \
  --min-shared-reads 10 --min-vote-margin 5 --min-purity 0.90 \
  --require-both-haplotypes --threads 8
```

## 9. Ranked next steps

### Surjected-BAM private sites: measured opportunity

The proposed direction is sound: the graph alignment can be surjected to BAM,
and high-confidence heterozygous calls from that BAM can supply anchors that do
not exist in the snarl catalog. The existing `collect-hybrid-variation` command
already uses the surjected BAM, but it is **BAM-base + graph augmentation**. It
does not preserve the graph partition as the core while adding only private BAM
sites, and it rediscovers BAM candidates internally instead of accepting the
DeepVariant call set supplied to HiPhase/LongPhase/WhatsHap.

A fresh full-chr20 hybrid run produced 134,646 candidates and 64,184 phased
positions. Against the strict exact-position target of 3,197 truth hets phased
by at least one competitor but outside graph phase blocks:

- hybrid phases 808/3,197 at the exact site (25.3%);
- hybrid phase blocks cover 1,624/3,197 (50.8%);
- only 201/1,012 no-catalog sites are recovered exactly;
- all 808 exact recovered sites are emitted as phased hets, so there is no
  called-but-unphased pool; the other 2,389 sites never enter hybrid output.

This localizes the next missing layer to **site input**, not BAM availability or
post-call phasing. Competitors receive the DeepVariant VCF directly; pgphase has
no external linear-VCF seed path. The next prototype should therefore be a
graph-core mode that:

1. imports high-confidence heterozygous sites from a VCF made from the surjected
   BAM, initially the same DeepVariant VCF used in the competitor comparison;
2. keeps only sites absent from the graph catalog (with normalized allele-aware
   matching, not position-only matching);
3. obtains per-read alleles for those private sites from BAM CIGAR/base evidence;
4. locks existing graph-core haplotypes and phases private-site gap reads into
   local/disjoint phase sets first;
5. merges a private block into a graph block only with direct spanning-read
   orientation support on both haplotypes and a positive vote margin.

Do not feed all BAM sites into a fresh global k-means pass: prior hybrid tests
show that the missing reads concentrate in segmental duplications and have about
21% BAM phasing error. The additive/disjoint-PS design prevents those hard reads
from reorienting the accurate graph core. Evaluate recovery, Hamming, switches,
and NGC50 on the shared DeepVariant call set; candidate count alone is not a
success metric.

### Private-site and block-merge experiment (completed)

`scripts/augment_graph_catalog_with_linear_hets.py` injects biallelic linear
hets as synthetic catalog rows for an experimental hybrid run. By default it
now excludes exact REF/ALT alleles already present in the graph. The synthetic
walks are intentionally unsupported by GAF; their evidence comes only from the
surjected BAM and therefore passes the existing hybrid depth, AF, noise, and
phasing gates.

There are 22,084 exact private alleles among all 77,484 biallelic DeepVariant
hets, or only 4,206 after PASS/GQ20 filtering. Adding those 4,206 private sites
yielded 183 extra hybrid candidates but only one extra transferred phase call.
It improved hybrid N50/NGC50 from 312/290 to 314/296 kb with unchanged
switchflips/Hamming (41/68). Merging its weak PS overlaps reached 397/366 kb but
worsened Hamming to 379. Private-site absence is therefore real but not the
primary bottleneck under the current candidate classifier.

The larger stress test deliberately added all 77,484 biallelic sites, including
55,400 exact graph matches. It yielded 1,748 additional hybrid candidates and
253 additional transferred calls. Its larger gain means the useful missing
layer is **variant-first anchor usage at represented sites**, not simply catalog
completeness.

`scripts/phase_vcf_from_hp.py` now supports a parity graph over overlapping PS
calls. Conservative output calls use two reads per haplotype at 0.70 purity;
merge proposals may use one clean haplotype at 0.60 only when the other PS has a
full two-haplotype call. Requiring two agreeing sites, margin two, summed read
support two, and support four for PS-label jumps over 500 kb produced:

| configuration | assessed | blocks | N50 kb | NGC50 kb | switchflips | Hamming |
|---|---:|---:|---:|---:|---:|---:|
| corrected graph transfer | 59,015 | 276 | 215 | 204 | 15 | 16 |
| corrected hybrid transfer | 59,665 | 298 | 312 | 290 | 41 | 68 |
| all-site seed, no PS merge | 59,672 | 303 | 290 | 267 | 39 | 62 |
| **all-site seed + conservative PS merge** | **59,447** | **247** | **401** | **374** | **39** | **79** |

Aggressive single-haplotype merging reached raw NG50 595 kb, but Hamming rose
to 2,048 or more. The best coordinate-exclusion oracle reached N50/NGC50
425/393 kb at Hamming 232, still below LongPhase's corrected NGC50 and far below
HiPhase. Threshold tuning has therefore reached its useful limit on chr20.
Further contiguity requires genuinely new spanning evidence, likely longer
reads, trio/assembly information, or a variant-first joint phasing model; it
cannot safely be manufactured from unsupported PS adjacency.

### Graph + private-gap joint phasing (current direction)

The intended architecture is now implemented as an opt-in hybrid mode. First,
`scripts/extract_private_gap_sites.py` computes the complement of the finalized
graph phase blocks and emits only biallelic linear hets that are absent from the
graph catalog by exact POS/REF/ALT match. It clears inherited GT phasing and PS.
Then `collect-hybrid-variation --private-sites FILE` removes every BAM candidate
not in that VCF before profiles are built. BAM alleles and counts at graph-owned
candidates are explicitly cleared and repopulated only from GAF, so the joint
matrix is exactly **graph observations + whitelisted private BAM observations**.
BAM noisy-region MSA recall is disabled in this mode to prevent unlisted sites
from entering later.

HG002 chr20 has 364 graph gaps. Private PASS sites inside them:

| minimum GQ | offered sites | phased reads | read Hamming | read N50 kb | DV N50/NGC50 kb | DV switchflips | DV Hamming |
|---:|---:|---:|---:|---:|---:|---:|---:|
| graph only | 0 | 175,843 | 0.141% | 937 | 215/204 | 15 | 16 |
| 20 | 271 | 178,284 | 0.180% | 945 | 280/224 | 20 | 21 |
| **10** | **958** | **180,800** | **0.202%** | **952** | **300/280** | **25** | **27** |
| 0 | 6,607 | 181,401 | 0.212% | 965 | 303/280 | 32 | 36 |

GQ10 is the selected experimental point. Against graph-only it adds 4,957
truth-evaluable reads, raises read N50 by 15 kb, and raises corrected shared-call
NGC50 by 37%, while Hamming remains 27 versus LongPhase 605 and HiPhase 2,279.
GQ0 adds errors without improving NGC50. This is the first tested path that
fills only graph gaps, phases graph and private evidence together, and improves
both read and variant contiguity without sacrificing the graph accuracy claim.

Core audit at GQ10: 173,490 graph-only reads remain evaluable in the joint run.
Of these, private evidence fixes 25 previously discordant assignments and breaks
24 previously concordant assignments. The margin gate drops 2,380 marginal
graph reads and the joint run adds 7,333 reads, 7,184 of them concordant.

### All clean BAM sites: tested, not selected

`collect-hybrid-variation --graph-authoritative` now implements the broader
control requested after the gap-only experiment. It runs normal BAM candidate
discovery/classification first and retains every surviving BAM candidate. Every
candidate represented in the graph catalog is then made graph-owned: all BAM
profile alleles and count fields at that candidate are cleared and only GAF
observations repopulate it. BAM evidence survives only at non-graph sites, and
all graph plus non-graph candidates are phased jointly. `--private-sites`
implies the same authoritative ownership rule after applying its whitelist.

This broader mode did not beat selective gap filling on HG002 chr20:

| config | read margin | evaluated reads | read Hamming | read N50 kb | assessed | DV N50/NGC50 kb | switchflips | DV Hamming |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| **GQ10 gap-private** | **2** | **180,800** | **0.202%** | **952** | **59,148** | **300/280** | **25** | **27** |
| all clean BAM sites | 2 | 184,911 | 0.439% | 965 | 59,269 | 300/260 | 36 | 133 |
| all clean BAM sites | 3 | 169,685 | 0.431% | 1,014 | not run | not run | not run | not run |

All sites add about 4,100 evaluated reads at margin 2, but competing BAM-only
anchors inside coherent graph blocks more than double read error and reduce
corrected shared-call NGC50. Raising the read margin removes 15,226 reads
without recovering accuracy, so weak output assignments are not the primary
problem. Keep this mode for controlled experiments; the selected policy remains
private GQ10 sites only inside graph phase gaps.

### Fixed: first graph chunk loaded the full VCF

`load_sites_for_region` accepts 1-based inclusive VCF coordinates, but graph and
hybrid callers passed `region.beg - 1`. At the first chunk this became zero;
the tabix query formatter omitted coordinates when `beg == 0`, loading all
977,275 chr20 catalog records into that chunk instead of roughly 6,000. All
three VCF site-query callers now pass `region.beg`; GAF and FASTA queries retain
their separate 0-based conversion. The corrected GQ10 run reproduced 180,800
evaluated reads, 366 discordant, 952 kb read N50, and the exact prior shared-call
300/280 kb N50/NGC50 with 25 switchflips and 27 Hamming errors. The bug was a
large runtime/memory defect and a one-base chunk-overlap error, but did not
materially change this chr20 accuracy result because downstream active-region
gating discarded the misplaced first-chunk candidates.

### Regional diagnosis and supported bridges

Competitors are not generally better inside pgphase's six <60%-accuracy phase
sets. On the union of those exact intervals, pgphase has 68/1,369 discordant
reads (5.0%), versus HiPhase 644/2,552 (25.2%), LongPhase 943/3,068 (30.7%),
and WhatsHap 184/1,752 (10.5%). They phase more reads, but mostly with worse
accuracy. Two small intervals near 31.8 and 32.45 Mb are genuine competitor
wins and remain useful diagnostic targets.

The BAM-fallback hypothesis was tested with `--bam-authoritative-bed`. Broad
cenSat fallback worsened discordance 366 -> 383; terminal fallback was neutral;
an oracle two-interval fallback improved 366 -> 355 but could not be detected
from clipping, MAPQ, or GAF aligned fraction. Keep regional ownership opt-in.

Private sites do act as bridges. `extract_private_gap_sites.py --bam FILE`
retains only a supported read-overlap path from one graph block through private
sites to the next. It selected 203/958 sites across 112/364 gaps and improved
discordance 366 -> 318 and read N50 952 -> 971 kb, but lost about 2,900 evaluated
reads because one-sided/local gap islands were excluded. It is therefore a
useful high-precision mode, not the selected balanced mode.

The strongest chr20 operating point is all GQ10 private-gap sites plus
`--min-phase-set-reads 50`: 179,494 evaluated reads, 269 discordant (0.15%), no
phase set below 60% accuracy, and 991 kb read N50. Shared-DV evaluation retains
300/280 kb N50/NGC50 with 18 switchflips and Hamming 19. This threshold remains
experimental until reproduced on chr12/chr18 and another sample.

The two competitor-win intervals have now been root-caused. They are
paternal-only copy/cenSat alignments: 297/298 input reads in the broader first
window and 72/72 in the second are paternal by diplinator truth. The catalog
contains hundreds of correlated nested records under one parent snarl, and the
hybrid allele matrix shows graph-owned anchors repeatedly separating the same
paternal reads into a 5:4 copy partition. No private site causes the split.

Competitors do not recover two haplotypes. HiPhase makes the same 5:4 split on
the exact first-block reads. In the second interval it emits a one-sided block
with 20 paternal reads all tagged HP2 and no HP1; LongPhase abstains and
WhatsHap's containing block is 50% accurate. The apparent win is primarily
one-sided tagging or wider-block aggregation. `--min-phase-set-reads 50` safely
abstains on both pgphase blocks. A tested soft-mask mono-haplotype heuristic
fixed them but damaged a truly diploid centromeric block, so it was rejected.
`--phase-matrix-dump PREFIX` is available in hybrid mode for this audit.

```bash
python3 scripts/extract_private_gap_sites.py \
  --graph-phased-vcf graph.vcf.gz --graph-sites sites.vcf.gz \
  --linear-vcf deepvariant.vcf.gz --contig chr20 \
  --output private-gap.vcf --gaps-bed graph-gaps.bed

# High-precision alternative: emit only complete read-supported bridge chains.
python3 scripts/extract_private_gap_sites.py \
  --graph-phased-vcf graph.vcf.gz --graph-sites sites.vcf.gz \
  --linear-vcf deepvariant.vcf.gz --contig chr20 --bam surjected.bam \
  --min-bridge-reads 2 --output private-bridges.vcf

./pgphase collect-hybrid-variation \
  --ref ref.fa --bam surjected.bam --graph-sites sites.vcf.gz --gaf reads.gaf.gz \
  --private-sites private-gap.vcf --min-read-margin 2 --min-phase-set-reads 50 \
  --stitch-min-margin 0 --stitch-rule 0 -r chr20 -t 20 \
  -o candidates.tsv --phased-vcf-out phased.vcf -b phased.bam
```

1. **Allele clustering inside snarls.** Group a snarl's walks by sequence
   similarity so a one-repeat-unit difference is not a different haplotype. This
   is the principled version of the top2 gate — POA over reads' walks rather than
   exact index equality.
2. **Catalog completeness against GBZ/pantree.** The current run accounts for
   essentially every record in `chr20.sites.vcf.gz`, but 1,012/3,197 exact
   competitor-site misses have no catalog record within ±25 bp. Need compare
   `build-snarl-catalog` output to the source graph variation directly.
3. **Retained-but-unphased sites.** `--snarl-allele-phasing` converts many
   `high_af` misses into retained graph-site rows, but they still do not join
   blocks. Diagnose whether they fail because they are isolated, non-anchors, or
   conflicting anchors.
4. **Chunk placement by node, not VCF POS** (still visible as
   `no_reads_in_chunk`). Long snarl
   records are returned by tabix into chunks far from where their nodes sit.
   `--gaf-pad` does not fix this; site→chunk assignment must use node position.
5. **Decide the `--anchor-af-margin` operating point.** 0.12 → 0.50 buys +24%
   NGC50 for Hamming 19 → 49, still ~12x better than longcallD. Confirm on the
   shared call set before changing the default.
6. **Second sample (HG003) and more chromosomes.** Single-sample chr20 results
   are the most attackable part of the claim. Truth BAMs for chr18/chr12 exist.

---

## 10. Reproduction

```bash
# build
make all                       # also: make test_graph_phase etc., 6 test binaries

# the graph pipeline (current default operating point)
./pgphase collect-graph-variation \
  --ref /home/kokyriakidis/Downloads/pgbam-experiments/chm13v2.0.fa \
  --sites /home/kokyriakidis/Downloads/chr20.sites.vcf.gz \
  --gaf /home/kokyriakidis/Downloads/pgbam-experiments/HG002.chr20.annotated.coord.gaf.gz \
  --min-read-margin 2 --anchor-af-margin 0.12 \
  -o cands.tsv --phased-vcf-out out.vcf --phased-bam-out out.bam -t 20

# read-level accuracy (representation-independent — prefer this)
python3 scripts/evaluate_phase_accuracy.py out.bam \
  /home/kokyriakidis/Downloads/pgphase-eval-data/truth/chr20/diplinator_merged.bam \
  0 0 5 "" evaldir samtools "" "" "" ""

# variant-level, apples-to-apples: transfer HP onto the shared DV call set first
python3 scripts/phase_vcf_from_hp.py tagged.bam dv.vcf.gz out.vcf \
  --sample HG002 --region chr20 --support-cache support.tsv
whatshap compare --sample HG002 --tsv-pairwise cmp.tsv --switch-error-bed err.bed truth.vcf.gz out.vcf.gz
whatshap stats --block-list blocks.txt out.vcf.gz
# NGC50 = NG50 over (blocks BED minus switch-error BED), denominated on the chromosome
python3 scripts/compute_ngc50.py blocks.txt err.bed ngc.json --genome-size 66210255
```

**Key data paths** (this machine):
- reference `~/Downloads/pgbam-experiments/chm13v2.0.fa`
- catalog `~/Downloads/chr20.sites.vcf.gz` (977,275 records; INFO has LV/PS/AT, **no PA**)
- GAF `~/Downloads/pgbam-experiments/HG002.chr20.annotated.coord.gaf.gz`
  (pggaf format: 3 leading coord columns, so the **path is field 8**, MAPQ field 14, read name field 3)
- truth BAMs `~/Downloads/pgphase-eval-data/truth/{chr20,chr18,chr12}/diplinator_merged.bam` (`HO:Z:MAT/PAT`)
- phaser env: `micromamba run -n bench-phasers ...` (whatshap, bcftools, bedtools, HiPhase, LongPhase)

**Companion docs:** `docs/competitor_evaluation_analysis.md` (full detail, §5–§12),
`docs/publication_plan.md`, `docs/research_plan.md`, `docs/pantree_catalog_experiment.md`.
