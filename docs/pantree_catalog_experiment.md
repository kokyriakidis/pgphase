# pantree reference-tree catalog: chr20 evaluation runbook

## Outcome: not pursued. The catalog is not the constraint on chr20.

Step 1's pre-check was run against the existing chr20 artifacts and falsified the
hypothesis before any pantree run. Recorded here so the question stays closed.

**1. Anchor density does not separate good phase sets from bad ones.** Over the
284 evaluated phase sets in `chr20_bench/graph_eval/per_phase_set.tsv`, surviving
candidates per kb are 1.17 (perfect blocks), 1.12 (good), 0.92 (accuracy < 0.95),
and the fraction of catalog sites used is flat at 72--76% across all three. At
0.92/kb an 18 kb HiFi read still carries ~16 anchors.

**2. The blocks that looked site-poor are not catalog-poor.** PS 44779026 and
44820071 have 11,419 and 11,399 catalog sites in their spans (13.1 sites/kb) and
yield 121 and 118 candidates -- a ~99% loss rate. The worst-switching block on the
chromosome, PS 65995586 (29 switches), has the *highest* density measured:
21.2 catalog sites/kb, 1.89 candidates/kb, at adequate depth.

**3. The loss is heterozygosity, not filtering, coverage, or catalog gaps.** Run
with `--filtered-sites-out` over chr20:44.7--45.9 Mb (15,706 dropped, 147 kept):

| reason | count | share | median depth at site |
|---|---|---|---|
| `ref_only` | 13,882 | 88.4% | 66 |
| `high_af` | 1,679 | 10.7% | 64 |
| `low_af` | 117 | 0.7% | 60 |
| `low_depth` | 28 | 0.2% | 3 |

Depth is not the issue -- that window carries 5,003 alignments, *more* than the
well-phased control at chr20:29.0--30.2 Mb (4,267), with identical candidate DP
(median 56 in both). The sites are dropped because **they are homozygous in
HG002**: no read takes an alternate allele at 66x coverage.

**4. Root cause: runs of homozygosity.** Clean het SNPs per Mb across chr20
(median 852):

| window | all het candidates | clean het SNPs | vs. median |
|---|---|---|---|
| 27--28 Mb | 129 | 117 | 13.7% |
| **44--45 Mb** | 69 | **3** | **0.4%** |
| 45--46 Mb | 219 | 68 | 8.0% |

Every degraded phase set identified in the analysis falls inside these three
windows. There are **three heterozygous SNPs in the entire 44--45 Mb window**. No
site catalog, alignment method, or phasing model can phase a region where the two
haplotypes are identical. A richer catalog would add homozygous sites.

**Not to be confused with HipHap.** `github.com/jheinz27/hiphap` is the renamed
`diplinator` -- the Rust tool this repo's `evaluate_phase_accuracy.sh` uses to
assign each read to its best haplotype of a diploid assembly and emit the HapQ
score. It is part of the evaluation harness, unrelated to pantree or to variant
cataloging.

**Conclusion.** pantree would add candidate sites to the top of a funnel that
already discards 96% of what it has, in regions that fail for lack of
heterozygosity. Not pursued. The scripts and the shape test are kept: the shape
test is a useful regression check on `src/graph_sites.cpp` site validation, and
the comparison script works on any two catalogs.

**Caveat on `per_phase_set.tsv`.** `n_reads` counts reads assigned to that phase
set and evaluated against truth, *not* read depth, and `span_bp` is measured in
HG002 truth coordinates while `phase_set` is a CHM13 coordinate. Treating
`n_reads / span_bp` as coverage understates real depth by ~50x inside a run of
homozygosity, where few reads can be phased at all, and `PS + span_bp` is not a
valid CHM13 interval (2 of 284 blocks run past the contig end).

---

## Hypothesis

`build-snarl-catalog` wraps `vg deconstruct -a`, which defines variants by snarl
(superbubble) decomposition against the reference path. That representation cannot
emit a variant whose REF and ALT alleles are *both* absent from the linear
reference. pantree ([repo](https://github.com/oclb/pantree),
[Cell Genomics 2026](https://doi.org/10.1101/2025.08.04.668502)) defines a variant
instead as any graph edge outside a spanning reference tree, which does capture
those. The paper reports **+46.1% SNPs in segmental duplications** (2,313,929 vs
1,584,132) and 1,896,587 SNPs that `vcfwave` re-alignment also misses.

pgphase's residual phasing error is concentrated in exactly those regions --
CHECKPOINT.md records that switch errors are *haplotype ambiguity in segmental
duplications* and that beating the BAM pipeline on switches requires changing the
segdup model, not the thresholds. Paralog-distinguishing SNPs are frequently
non-reference-vs-non-reference differences, i.e. the class deconstruct drops.

**Hypothesis:** a pantree catalog supplies paralog-distinguishing anchors inside
segdups that the superbubble catalog structurally cannot, and those anchors reduce
switch/Hamming error there.

**What would falsify it.** Two things, and both are live:

1. CHECKPOINT.md's FN root-cause analysis found the chr20 misses are alignment and
   evidence limits -- indel representation (57/114) and het allele dropout (46/114)
   -- not site-catalog gaps. So **global F1 should not move**, and it moving is not
   the result being tested. Only segdup-stratified switch/Hamming counts are.
2. CHECKPOINT.md also records that in segdups the graph inherits the BAM's wrong
   paralog placement. If a read sits on the wrong copy, denser PSV anchors add
   confidently wrong evidence. **Switch errors getting worse is a real possible
   outcome**, and it is the outcome that settles the question against pantree.

Run step 1 before step 3. If step 1 shows no meaningful segdup-specific SNP gain,
stop -- the expensive phasing run has nothing to find.

## Prerequisites

```bash
# pantree (Python; not vendored -- research code, pinned by you, not by this repo)
git clone https://github.com/oclb/pantree.git && cd pantree && uv venv && uv sync

# Per-chromosome GFA from the same GBZ pgphase phases against
vg chunk --gbz --contig chr20 -x full.gbz -o /tmp/chunk_chr20
vg convert -f /tmp/chunk_chr20_graph_0_chr20.gbz > chr20.gfa

# Segdup annotations (already scripted)
./scripts/build_difficult_regions_bed.sh
```

**Cost.** The paper reports 36 GB RAM / 10 h wall for chr1; chr20 is roughly a
quarter of that. Compare against ~1 GB for the current per-chr `vg deconstruct`.
This is one-time preprocessing per graph version, like `build-snarl-catalog`, but
it does not fit the existing `scripts/build_catalog_allchr_v21.sbatch` resource
profile.

## Step 1 -- build the pantree catalog and diff it (cheap, decides the rest)

```bash
python3 scripts/pantree_to_pgphase_catalog.py chr20.gfa \
    -o chr20.pantree.sites.vcf.gz \
    --chrom 'CHM13#0#chr20' --ref-name CHM13 \
    --no-missingness \
    --stats chr20.pantree.stats.tsv --verbose
```

`pantree gfa2vcf` output is not directly usable: it carries no `AT` field, and `AT`
is the only thing `collect-graph-variation` matches reads against
(`src/graph_sites.cpp`, `src/graph_query.cpp`). The converter loads the graph once
through pantree's Python API and synthesizes, per variant edge, the reference-tree
path and the variant-edge path between the *same* pair of flanking handles -- the
boundary invariant `graph_site_validation_skip_reason()` enforces. POS/REF/ALT/VT/
AC/AN come from pantree's own writer internals, so the catalog stays comparable to
plain `pantree gfa2vcf` output.

Non-reference variants (pantree `REF='.'`, allele in `NR`) are written `REF=N` /
`ALT=*` with `PTNONREF`. They stay usable as phasing anchors -- the allele walks are
intact -- but `src/graph_collect.cpp` skips `*` alleles when deriving variant keys,
so they never reach the phased VCF. Use `--nonref drop` to exclude them entirely.

Then diff against the current catalog:

```bash
python3 scripts/compare_catalogs.py \
    --a chr20.sites.vcf.gz          --a-label deconstruct \
    --b chr20.pantree.sites.vcf.gz  --b-label pantree \
    --bed data/annotations/hg002v1.1.segdups.bed --bed-label segdup \
    --contig chr20 --contig-length 64444167 \
    --window 100000 \
    --out-prefix eval/pantree/chr20_catalog_cmp
```

POS/REF/ALT cannot be compared across the two catalogs -- deconstruct draws allele
boundaries at snarl endpoints, pantree at the reference-tree branch point, so the
same variation gets different coordinates and different flanking context. The
script therefore decomposes both catalogs into **canonical variant edges** (the
edges an alt walk traverses that its ref walk does not) and compares those sets.

**Decision rule.** Read `pantree_only`, class `SNP`, stratum `segdup`, as
anchors/Mb. A large gain concentrated there is the only result that justifies
going on. A gain that is mostly outside the BED, or mostly INS/DEL/REP (tandem
repeat expansion/contraction noise, which the paper notes dominates small indels),
does not support the hypothesis. `<prefix>.windows.tsv` breaks it down per 100 kb
so it can be lined up against the 37--38 Mb and 45--46 Mb windows already flagged
in CHECKPOINT.md.

## Step 2 -- verify the catalog shape is ingestible

```bash
./scripts/test_pantree_catalog_shape.sh
```

A pantree catalog differs structurally from a deconstruct one: every record is
biallelic (one variant edge = one binary site), there is no `LV`/`PS`/`PA` nesting
so no parent gating, and the VCF `ID` is the oriented variant edge and doubles as
the site key (`graph_site_key_str`, `src/graph_sites.hpp:145`).

This test reshapes the checked-in chr20 deconstruct catalog into exactly that form
(no pantree install needed), runs `collect-graph-variation` on both, and checks the
reshaped run still phases. On the `test_data/graph_chr20` fixture it yields 414
phased records vs 417 for the original -- **the contract holds with no C++
changes**. Run it after any change to `src/graph_sites.cpp` site validation.

## Step 3 -- phase and measure (only if step 1 passes)

Run the graph pipeline three ways over the same GAF: deconstruct catalog alone,
pantree catalog alone, and the two concatenated (sorted, re-indexed). The union is
the configuration that actually matters -- pantree is a *supplement* to the
superbubble catalog, not a replacement, because deconstruct's nested multiallelic
snarls carry structure a flat binary-edge catalog discards.

```bash
for cat in deconstruct pantree union; do
    ./pgphase collect-graph-variation \
        --ref ref.fa --sites chr20.${cat}.sites.vcf.gz \
        --gaf reads.coord.gaf.gz \
        --phased-vcf-out chr20.${cat}.phased.vcf \
        --phased-bam-out chr20.${cat}.phased.bam -t 16
done
```

Score with the existing harness, then stratify:

```bash
./scripts/evaluate_phase_accuracy.sh   # switch / Hamming vs truth
python3 scripts/phase_block_stats.py   # block NG50 / count
```

**Report switch and Hamming error inside the segdup BED separately from outside.**
A whole-chromosome average will dilute the effect being tested in both directions:
segdups are a small fraction of chr20, and any global F1 change is noise against
the FN root-cause finding above.

## Known limitations

- **Binary edges discard multiallelic structure.** pgphase's clean-het gating works
  on a site with `n` alleles and an allele fraction. Split into independent binary
  edges, a read that takes a third path through a complex bubble registers as "not
  alt" rather than as a distinct allele, which can bias AF and inflate false hets.
  This is the main reason to test the union rather than the replacement.
- **Non-reference POS collides.** pantree defines the position of a non-reference
  node as the highest reference position with a tree path to it, so an entire
  insertion branch collapses onto one coordinate. Many sites then share a POS.
  Tabix retrieval is unaffected; ordering assumptions downstream may not be.
- **vg-built graphs.** pantree's authors recommend Minigraph-Cactus graphs; on
  vg-built graphs the linear-reference walk may contain cycles, making POS
  ambiguous. HPRC MC graphs are fine. Pass `--no-missingness` regardless.
- **Private API dependency.** The converter imports pantree's `_VariantData` and
  `_VariantRecord` to keep POS/REF/ALT identical to `pantree gfa2vcf`. If pantree
  changes those, the converter errors with a message saying so rather than silently
  diverging. Pin the pantree commit you evaluated.
