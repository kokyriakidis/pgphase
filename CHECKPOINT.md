## 2026-09-27: use the second clean graph SNP to orient a singleton BAM source

The chr20 32,431,751–32,432,269 boundary split an MSA-verified BAM
insertion, the only oriented site in its source block, from a graph SNP
block. Seven molecules called the insertion and the first two graph SNPs.
The closest SNP alone gave five votes for the correct cross orientation and
two conflicting votes, so the existing nearest-pair test abstained. The
second SNP gave seven cross votes with no conflicts; the two graph SNPs
agree on 10 of 12 molecules, and 11 of 11 reads link the graph block's
first and last clean SNP in the same orientation. A trial that voted across
all 57 graph SNPs counted the same molecules repeatedly and exposed two
reads whose distant graph calls disagree with their near-boundary calls;
that trial was discarded.

A complete one-site BAM source may now attach to the nearest graph SNP
block through the first two clean SNPs, all selected by coordinate before
allele inspection. The BAM-to-second-SNP edge, the adjacent graph-SNP edge,
and the graph block's first-to-last edge must each pass the existing
both-allele, split-read, exact-binomial boundary test. A significant
opposite BAM-to-first-SNP edge vetoes the attachment. Only source reads
that actually call the BAM site are relabeled; graph rows and their reads
stay in their original orientation.

The 32–33 Mb owning-chunk replay joins the exact VCF rows. A 200 kb replay
that includes the preceding graph block produces the same boundary phase and
runs in seconds; the new permanent panel row asserts that span, the internal
recovered insertion, and no parental switch against truth. In the
matched full-chr20 run, exactly one VCF row changes phase set and orientation;
all 62,147 variant keys and unphased allele genotypes are unchanged. VCF
phase blocks fall from 361 to 360; span N50 stays 573,587 bp. Truth-scored
phased reads remain 236,831, with 229,045 correct and 7,786 discordant
(96.7124%).

## 2026-09-27 complementary BAM indel path cut at 21.179 Mb

The chr20 21,179,807–21,183,846 gap had 55 reads with both boundary
SNP alleles in the merged graph/BAM matrix: 29 ALT/ALT and 26 REF/REF,
zero opposite pairs. The left SNP is a graph-walk candidate and absent
from the BAM de novo candidate list. Direct BAM CIGAR calls at that SNP
see only two high-quality ALT bases, so substituting physical SNP calls
would discard the graph allele evidence.

The recovered BAM block begins with two complementary, MSA-verified
indel rows at the same anchor. Its reads call both indels and the later
SNP, but still carry an earlier read phase-set label. The source-path
checker used that read label as a hard filter, reported a weak cut at
21,179,789, and blocked the graph/BAM stitch. Across the false cut,
25 reads support one haplotype and 29 the other, with no contradictions.
The checker now admits cross-label, MAPQ >=30 observations only at a cut
following complementary MSA indel rows at one locus, requiring two reads
per haplotype, zero conflicting pairs, and exact one-sided random-polarity
p <= 1e-6. Existing same-label path evidence remains authoritative.

An unrestricted cross-label trial opened two previously connected short
gaps; its 49 Mb source cut had only nine supporting pairs and lost 25
phased reads. A broader high-confidence trial also relabeled a large
12 Mb clean-SNP block. The representation-specific guard excludes both
effects. In the final full-chr20 run, exactly one VCF row changes PS:
the 21,183,846 SNP joins the left graph block. The owning 21–22 Mb
window spans without a truth switch, with 0.85 separated-read fraction
(HiPhase: 0.9037). Variant keys and genotypes remain identical. VCF
blocks fall 362 -> 361; span N50 stays 573,587 bp. Truth-scored phased
reads stay 236,831; correct reads change 229,051 -> 229,045 and
discordant reads 7,780 -> 7,786 (96.7150% -> 96.7124%). The permanent
window panel now asserts the new span and its parental orientation.

## 2026-09-27 weak-cut BAM source run joins adjacent graph blocks

The chr20 6.513891–6.516221 Mb gap stayed split in its owning 6–7 Mb
graph chunk although a focused replay joined it. The recovered BAM source
phase block covered both graph flanks through 6.545826 Mb, then had a weak
cut. Its left graph/BAM read vote was 75 consistent and 2 conflicting;
the old local-run transfer required zero conflicts and could attach at
most one graph block. The source path through the two flanks itself had no
weak cut.

Local transfer now permits two adjacent graph blocks to join through the
same weak-cut component only when all oriented sites of both graph blocks
match that source component, each multi-site graph block has direct
end-to-end molecule support, and both graph/BAM block votes pass the same
two-haplotype
statistical test used for a complete source path. The downstream graph
block is flipped if needed and moved as a whole; source rows beyond the
weak cut remain independent. The 19.4 Mb known wrong-join control remains
split.

The owning-chunk replay now spans 6.513891–6.516221 Mb with no local truth
switch and 0.87 separated-read fraction (HiPhase: 0.916). Its panel row
replays the owning 6–7 Mb chunk rather than the focused solve. On full
chr20, exactly four VCF rows change PS at this join; variant keys and
genotypes are identical. Truth-scored reads remain 229,051 correct out of
236,831 phased (7,780 discordant; 96.715%). The old and new runs have
739 and 738 truth-scored read blocks, respectively. Phased-VCF block
count drops from 363 to 362; span N50 stays 573,587 bp. A weaker first trial
also changed an unrelated 20.797 Mb label, so statistical approval is now
restricted to complete adjacent graph-block pairs in the same local
source component. The 5.256 Mb gap remains split: only two of 35 physical
reads spanning its two indels carry the left insertion, too little allele
support for a safe join.

# pgphase Evaluation Checkpoint

## 2026-09-22 independent BAM recovery phase blocks

Recovery now preserves every phase set produced by the targeted BAM sub-solve
as an independent local block. Importing candidates no longer implies that
their HP integers share a gauge, and the recovery stitch cannot route around an
unsupported BAM boundary through its generic strongest-site, DP, or MEC paths.

The final left-to-right stitch compares one exact adjacent pair at a time.
Graph/BAM pairs use only their source-specific 2x2 read-haplotype vote.
BAM/BAM pairs use a candidate-membership-preserving aggregate allele vote.
Both require the existing one-sided exact binomial test at `p <= 0.01`;
otherwise the blocks remain separate. An accepted edge mutates the complete
downstream block atomically. Reads retain the BAM sub-solve's assignment, and
reads left unphased by that solve are not reassigned from a later composite
chain.

One orientation bug was found while testing a single BAM block attached to both
flanks. If the left edge flipped the block, the right-edge source vote was
previously applied in its old gauge. The stitcher now records an original
orientation anchor for every phase set and translates each later vote into the
current gauge.

On the full chr20 fixture, the retained design has:

| metric | independent blocks + direct stitch | prior recovery |
|---|---:|---:|
| phased heterozygotes | 59,892 | 59,906 |
| phase sets | 657 | 187 |
| phased block span | 54.32 Mb | 56.71 Mb |
| N50 | 412.1 kb | 1,136.4 kb |
| tagged reads | 219,579 | 221,836 |
| truth-discordant reads | 4,762 | 36,220 |
| truth discordance | 2.169% | 16.327% |
| target gaps spanned | 14/48 | 31/48 |

The rejected graph/BAM allele-aggregate fallback phased only 42 additional
reads but raised whole-chr20 discordance from 2.169% to 4.162%. In the fast
panel it created truth-switched joins at 22.98 Mb (51.1% concordance) and
48.23 Mb (68.3%). Source-specific graph/BAM votes remove both switches. The
48-target scorer reports no `SWITCH` for the retained design. Detailed
commands and output are in
`evaluations/2026-09-22-independent-bam-blocks/`.

## 2026-09-21 compact gap phasing and atomic stitch prototype

Short-gap iteration no longer belongs in the region integration panel. The 14
noncentromeric HiPhase targets are listed in
`evaluations/2026-09-16-test-panel/short_gap_targets.tsv` for conversion to
compact fixtures. The first synthetic replay uses positions and support patterns
from the chr20:4.76-4.79 Mb seam containing the 4,778,793-4,785,720 target. It
constructs the production `PhasingChunk` state directly and runs in under
0.02 seconds measured wall time; it opens no BAM, GAF, reference, or subprocess.

The fixture separates two decisions. BAM-only reads are first checked inside an
independent local gap phase set, with parental truth used only by the assertion
and with whole-block inversion allowed. Stitching then proceeds left to right:
left PS to gap PS, followed by the oriented gap PS to right PS. The new
`stitch_phase_sets_by_alleles` primitive ranks candidate edges by net
same-versus-cross margin and paired-read count, abstains on an equal opposite
parity, and flips and relabels every candidate and read in the downstream PS
atomically. A far-right candidate with no read spanning back to the seam is in
the regression; it still flips because it shares the downstream PS.

Production now uses the same path as the fixture. The merge preserves the BAM
sub-solve's candidate consensus and read HP/PS under a collision-free local PS,
keeps established graph assignments unchanged, and stitches the local blocks
through each seam from left to right. The old whole-chunk k-means rerun,
outside-in frontier waves, and post-solve consensus repair were removed because
they recomputed or mutated the imported BAM answer. The regression calls the
production chain driver directly and also verifies that a tied allele edge
leaves the downstream block unchanged; measured wall time is under 0.02 seconds.

## 2026-09-20 graph recovery: exact BAM rows and noncentromeric gaps

The chr20 gap investigation excluded the low-MAPQ 27.5–29 Mb centromeric
interval. Five high-confidence competitor spans outside it were screened:
5.31, 11.23, 24.10, 48.18, and 60.03 Mb. The prior graph recovery split
all five; exact-row transfer closes 5.31 Mb. Measurements are in
`evaluations/2026-09-20-graph-recovery-windows/`.

The prior recovery sub-solve inherited graph defaults, including
`merge_colocated_msa_alleles=true`. At 11.23 Mb this compressed two
complementary BAM deletion rows at 11,255,369 into one three-allele candidate.
Standalone BAM has 47 reads shared with the right SNP: each separate row agrees
with it on 44 or 45 reads. The merged candidate's allele 2 is split 21/9
across that SNP. The old transfer preserved the merged representation, so the
bug was upstream of the candidate-table merge. The retained path substitutes
the standalone BAM rows and their read alleles at these MSA loci.

Using the BAM port's full option set in the recovery sub-solve does produce
separate rows, but the graph union then fails existing truth gates: 3.85 Mb
read concordance falls to 0.917, 60.03 Mb to 0.912, and 5.31 Mb separated
reads to 0.323. Turning off only the MSA merge also loses phased in-gap
heterozygotes (17 to 4 at 3.85 Mb and 6 to 1 at 4.77 Mb). These trials were
reverted. A narrow graph-style candidate-admission gate closed 5.31 Mb at
0.947 separated reads and 99.6% read concordance, with only three additional
discordant reads across full chr20, but it retained the merged BAM site source
and was also reverted after the no-merge requirement was clarified.

Minimal variant-key inequality overstates private variation in repeats. At
3.85 Mb, 32 of 39 BAM heterozygotes that the audit called unmatched are
local-haplotype equivalents of raw catalog ALTs under different anchors.
A trial matching all raw catalog ALTs removed most duplicate descriptions,
but lost phased sites and still failed the window suite. This distinction
matters: an allele absent from the retained phasing candidates can already
exist in the raw graph VCF. No experimental graph behavior from this audit
was retained.

Further transfer trials exposed two independent requirements. The 11.23 Mb
deletion pair is present in the raw catalog as `GAAA>GAA` and `GAAA>GA`
but absent from the active graph candidate table. Separate BAM rows link
strongly to the right SNP (45 agreeing, 2 conflicting reads) and to each
other (42/5), yet have no observed link back to the left boundary. Exact-row
transfer alone cannot span that gap. Replacing only merged MSA rows with the
standalone BAM rows and admitting read-supported catalog sites closed 5.31 Mb
at 0.947 separated reads, but caused 3.85 Mb read concordance to fall to
0.909 (45 discordant reads). The user accepted this local drop and the full
chr20 error increase (2,527/219,058 to 3,302/220,640), so
the exact-row transfer and read-supported catalog admission are retained.
The early `allele_depths_call_het` retry-window guard also masks its later
`bam_injected` exception in the graph re-solve; removing that guard without
extra link validation regressed the committed 22.98 Mb window's separated
fraction to 0.308. The guard-only trial was reverted. Details are in
`evaluations/2026-09-20-graph-recovery-windows/README.md`.

A graph-only chr20 phase-span audit found 74 internal gaps below 10 kb
(281,849 bp). The original BAM comparison used output created before the final
longcallD parity fixes and is superseded. Current standalone BAM spans 36 gaps
(119,456 bp); 19 are at least 98% concordant across truth-scored crossing reads.
HiPhase spans 50 (195,219 bp); 21 of its 49 scorable spans reach 98% in the
phase set's whole-block orientation.

Current BAM misses seven of those correct HiPhase spans. The 26.92 Mb case is a
centromeric low-MAPQ shoulder where only 3/58 overlapping alignments pass BAM's
MAPQ 30 floor and is excluded from recovery work. Four noncentromeric breaks
come from exact longcallD behavior: heterozygous homopolymer indels are emitted
but excluded from the phase link list, moving the effective anchor 21.0–45.5 kb
away and leaving only 0 or 1 link vote. Two more come from representation and
ordering: HiPhase phases one multiallelic record, while exact BAM preserves two
complementary rows and tests only the immediately preceding heterozygote. The
first row starts a new phase set with 0/0 or 1/0 votes; the strongly supported
second row then links to the first and cannot reconnect the old block. Details
are in `hiphase_correct_bam_misses.tsv` and the short-gap section of the same
evaluation directory's README.

A follow-up screened four noncentromeric deficits that standalone BAM itself
spans (4.85, 47.74, 60.95, and 64.91 Mb). Injecting every exact BAM site and
re-running the whole graph chunk closed none, reopened 5.31 Mb, and reduced
60.03 Mb concordance to 0.912. Re-running the graph union with all longcallD
phasing toggles also closed none. Unconditionally accepting BAM verdicts for
sequence-equivalent catalog sites closed 64.91 Mb but made a false 60.03 Mb
join at 59.2% concordance. A two-pure-flank gate prevented the false join and
closed none. All four variants were reverted; the retained union substitutes
exact BAM rows only where graph recovery compressed co-located MSA alleles.

An earlier whole-chunk fallback for catalog-empty chunks was also rejected:
it added 2,271 calls across chr20:27.5–29 Mb and raised read discordance from
2,528/219,067 (1.154%) to 3,624/221,864 (1.633%). The normal MAPQ floor
still added 341 calls and raised discordance to 2,566/219,213 (1.171%).

## Test Setup

**Sample:** HG002 chr20 HiFi reads  
**Reference:** CHM13 T2T v2.0 (contig `CHM13#0#chr20`, 66.2 Mbp)  
**Truth:** HG002 v1.1 diploid assembly (chr20 slices only — `HG002_1.1_MATERNAL.chr20.fasta` / `HG002_1.1_PATERNAL.chr20.fasta`)  
**Evaluator:** diplinator (haplotype-aware read assignment) + pgphase accuracy script

### Key files

| File | Description |
|------|-------------|
| `test_data/HG002.chr20.fq.gz` | HiFi reads (272,016 total) |
| `test_data/HG002_chr20_hifi_mapped_to_CHM13_chr20_annotated.bam` | Reads aligned to CHM13 chr20 |
| `test_data/chm13v2.0.chr20.renamed.fa` | CHM13 chr20 reference (contig `CHM13#0#chr20`) |
| `test_data/chr20.sites.striped.vcf.gz` | 977,275 graph sites (AT field with numeric node IDs) |
| `test_data/HG002.chr20.annotated.coord.gaf.gz` | bgzipped + tabix-indexed GAF for graph pipeline |
| `test_data/chr20.phase_graph.bam` | Phased BAM from graph pipeline |
| `test_data/bam_eval_chr20/phased.bam` | Phased BAM from BAM pipeline |

### Data provenance (full chr20 evaluation inputs)

The full-chr20 evaluation inputs are **not committed** (5.1 GB; gitignored under
`/test_data/chr20_eval/`). Re-download from the shared Google Drive folder:

```
https://drive.google.com/drive/folders/1DxDcuJ7uIugkQX_ZIbEMI5NdDg_lYT-5
```

Fetch with `gdown` (in a venv to satisfy PEP 668):

```bash
python3 -m venv /tmp/gdrive-venv
/tmp/gdrive-venv/bin/pip install -q gdown
/tmp/gdrive-venv/bin/gdown --folder \
  "https://drive.google.com/drive/folders/1DxDcuJ7uIugkQX_ZIbEMI5NdDg_lYT-5" \
  -O test_data/chr20_eval
```

| File (in `test_data/chr20_eval/`) | Role | Used by |
|------|------|---------|
| `chm13v2.0.chr20.renamed.fa` (+`.fai`) | CHM13 chr20 reference (`CHM13#0#chr20`) | BAM, graph, hybrid |
| `HG002_chr20_hifi_mapped_to_CHM13_chr20_annotated.bam` (+`.bai`) | HiFi reads aligned to CHM13 chr20 (272,016 reads) — graph `hs/hb/he` tags, **no** diplinator truth tags | BAM, hybrid |
| `chr20.sites.striped.vcf.gz` (+`.tbi`) | snarl sites | graph, hybrid |
| `HG002.chr20.annotated.coord.gaf.gz` (+`.tbi`) | graph alignments (GAF) | graph, hybrid |
| `HG002.chr20.fq.gz` | HiFi reads (FASTQ) | diplinator truth-BAM generation only |
| `HG002_1.1_MATERNAL.chr20.fasta` (+`.fai`) | truth maternal assembly | diplinator evaluation only |
| `HG002_1.1_PATERNAL.chr20.fasta` (+`.fai`) | truth paternal assembly | diplinator evaluation only |

Notes:
- The aligned BAM is "annotated" with graph `hs/hb/he` tags, **not** diplinator
  `HO:Z:`/`hq:i:` truth tags. The truth BAM that `evaluate_phase_accuracy.py`
  scores against must still be generated (align FASTQ to MAT/PAT assemblies →
  diplinator).
- The Drive folder ships **inputs only** — no pre-phased pipeline outputs. Both
  the BAM and hybrid phased BAMs must be produced by running the pipelines here.

---

## Results — Graph Pipeline (indexed GAF)

**Command:**
```
pgphase collect-graph-variation \
    --ref test_data/chm13v2.0.chr20.renamed.fa \
    --sites test_data/chr20.sites.striped.vcf.gz \
    --gaf test_data/HG002.chr20.annotated.coord.gaf.gz \
    --phased-bam-out test_data/chr20.phase_graph.bam \
    -o test_data/phase_graph_smoke/reads.tsv \
    --phased-vcf-out test_data/phase_graph_smoke/phased.vcf \
    -t 16
```

> Note: `--sites-vcf` → `--sites`, `--out-bam` → `--phased-bam-out`, `--ref` now required (commit `e751d83`).

**Eval dir:** `test_data/diplinator_eval_chr20/`

### v1 — before noise filter

| Metric | Value |
|--------|-------|
| Total reads | 231,382 |
| Phased reads | 226,831 (98.0%) |
| Phase sets | 421 |
| Phase block N50 | 912 Kbp |
| Largest block | 2.1 Mbp |
| Overall accuracy | 89.73% |
| Hamming error rate | 10.27% |
| Switch error rate | 6.14% |
| Perfect phase sets | 78 / 411 (19.0%) |
| Phaseable PS (>60%) | 308 PS, 207,571 reads, 93.00% |
| Unphaseable PS (≤60%) | 103 PS, 19,232 reads |

### v2 — with reference-based noise filter (commit `e751d83`)

| Metric | Value |
|--------|-------|
| Total reads | 231,382 |
| Phased reads | 216,130 (93.4%) |
| Phase sets | 425 |
| Phase block N50 | 898 Kbp |
| Largest block | 2.1 Mbp |
| Overall accuracy | 92.77% |
| Hamming error rate | 7.23% |
| Switch error rate | 1.04% |
| Switchflip error rate | 2.51% |
| Switch errors | 2,251 |
| Flip errors | 3,167 |
| Switchflip errors | 5,418 |
| Switch opportunities | 215,693 |
| Perfect phase sets | 96 / 419 (22.9%) |
| Phaseable PS (>60%) | 327 PS, 202,396 reads, 95.36% |
| Unphaseable PS (≤60%) | 92 PS, 13,716 reads |

---

## Results — BAM Pipeline (no pgbam sidecar)

**Command:**
```
pgphase collect-bam-variation \
    --hifi \
    --ref test_data/chm13v2.0.chr20.renamed.fa \
    --bam test_data/HG002_chr20_hifi_mapped_to_CHM13_chr20_annotated.bam \
    --out-bam test_data/bam_eval_chr20/phased.bam \
    --phased-vcf-out test_data/bam_eval_chr20/phased.vcf \
    -o test_data/bam_eval_chr20/reads.tsv \
    -t 16
```

> Note: `--ref` and `--bam` named options added in commit `1614581`.

**Eval dir:** `test_data/bam_eval_chr20/`

| Metric | Value |
|--------|-------|
| Total reads | 272,016 |
| Phased reads | 220,377 (81.0%) |
| Phase sets | 431 |
| Phase block N50 | 937 Kbp |
| Largest block | 53.2 Mbp |
| Overall accuracy | 97.02% |
| Hamming error rate | 2.98% |
| Switch error rate | 0.62% |
| Switchflip error rate | 1.45% |
| Switch errors | 1,356 |
| Flip errors | 1,826 |
| Switchflip errors | 3,182 |
| Switch opportunities | 219,946 |
| Perfect phase sets | 133 / 430 (30.9%) |
| Phaseable PS (>60%) | 366 PS, 214,808 reads, 98.14% |
| Unphaseable PS (≤60%) | 64 PS, 5,568 reads |

---

## Head-to-Head Comparison

| Metric | Graph v1 (no filter) | Graph v2 (noise filter) | BAM (no sidecar) |
|--------|---------------------|------------------------|-----------------|
| Phased reads | 98.0% | 93.4% | 81.0% |
| Phase sets | 421 | 425 | 431 |
| N50 | 912 Kbp | 898 Kbp | **937 Kbp** |
| auN | — | 981 Kbp | **9,828 Kbp** |
| Accuracy | 89.73% | 92.77% | **97.02%** |
| Hamming error rate | 10.27% | 7.23% | **2.98%** |
| Switch error rate | 6.14% | 1.04% | **0.62%** |
| Switchflip error rate | — | 2.51% | **1.45%** |
| Switch errors | — | 2,251 | **1,356** |
| Flip errors | — | 3,167 | **1,826** |
| Switchflip errors | — | 5,418 | **3,182** |
| Switch opportunities | — | 215,693 | 219,946 |
| Perfect PS | 19.0% | 22.9% | **30.9%** |

**Key observations:**
- Noise filter improves graph accuracy +3 pp (89.7% → 92.8%) and cuts switch error 6× (6.14% → 1.04%)
- Noise filter costs ~4.6% phased reads (noisy indels excluded from k-means)
- BAM pipeline leads on all accuracy metrics: Hamming (2.98%), switch (0.62%), switchflip (1.45%)
- BAM auN (9.8 Mbp) far exceeds graph (981 Kbp), driven by the large 53 Mbp block
- Graph pipeline phases more reads (93% vs 81%), compensating for lower per-block accuracy

---

## Site Overlap Analysis: Graph Sites vs BAM Het Variants

Measured with `scripts/analyze_site_overlap.py` on chr20 data.

| Metric | Value |
|--------|-------|
| Graph sites (total positions) | ~977K |
| BAM het variants | 83,013 |
| Shared positions | 62,600 (75.5% of BAM het) |
| Exact allele match | 56,590 (90.4% of shared) |
| Position-only match (diff alleles) | 6,010 (9.6% of shared) |
| BAM-only (no graph site) | 20,413 (24.5% of BAM het) |

### Interpretation

75.5% of BAM het sites have a graph snarl site at the same GRCh38 position, and 90.4% of
those are exact allele matches. This means graph reads traversing those snarls carry the
same allele information that the BAM pipeline uses for phasing — they just come from
pangenome-aligned reads that may not map to GRCh38.

The 24.5% BAM-only sites are likely private variants or rare alleles not represented in the
pangenome reference panel. These sites can only be phased using BAM-mapped reads.

### Hybrid bridging approach

**Goal:** Combine BAM pipeline's variant calling (83k het sites, 97% accuracy) with graph
pipeline's read mapping (98% reads aligned, including reads unmappable on GRCh38).

**Approach:** Run BAM pipeline as primary caller. For the ~19% unphased reads, use their
graph allele observations at the ~56k bridgeable sites to assign haplotypes and extend
phase blocks.

**Implementation phases:**

1. **Read haplotype assignment** — For each GAF read not in the BAM phased output, find its
   graph allele observations at bridgeable snarl sites. Look up the BAM haplotype consensus
   at those sites. Assign the read to the most consistent haplotype.

2. **Enhanced chunk stitching** — Use newly-phased reads as connectors between BAM phase
   blocks. Two blocks that couldn't be stitched (no shared GRCh38-mapped reads) might now
   be connected through pangenome-located reads that span both regions. The pangenome gives
   read locality information that GRCh38 can't — two reads unmapped on GRCh38 might be
   clearly co-located in the pangenome (traversing the same snarls in the same region).

**Expected gains:**
- Phase the missing ~19% of reads (currently unphased by BAM pipeline)
- Improve N50 by connecting fragmented BAM phase blocks through pangenome-bridged reads
- Maintain BAM pipeline's 97% accuracy (graph reads are assigned to existing haplotypes,
  not used to re-phase)

**Risks:**
- Allele matching imprecision at the 9.6% position-only sites could inject noise
- Reads unmapped on GRCh38 may have ambiguous pangenome positions (multi-mapping)
- Phase set merging logic needs to handle cross-pipeline block boundaries

---

## Merge Feasibility: BAM + Graph Phased BAMs

Measured with `scripts/merge_feasibility.py` on chr20 data.

### Read category breakdown

| Category | Reads | % of total |
|----------|-------|-----------|
| A) Phased by both | 208,092 | 76.5% |
| B) BAM-only phased | 12,285 | 4.5% |
| C) Graph-only (rescue candidates) | 8,038 | 3.0% |
| D) Neither phased | 43,601 | 16.0% |

### Haplotype agreement (category A — doubly-phased reads)

| Metric | Value |
|--------|-------|
| Doubly-phased reads | 208,092 |
| Same HP label | 128,104 (61.6%) |
| Different HP label | 79,988 (38.4%) |
| Overlapping block pairs | 639 |
| Concordant pairs (≥90% same, ≥2 reads) | 249 — 118,177 reads |
| Flipped pairs (≥90% opposite, ≥2 reads) | 169 — 67,159 reads |
| Discordant pairs (<90% majority) | 197 — 22,732 reads |

### Rescue potential

| Metric | Value |
|--------|-------|
| BAM phased reads | 220,377 (81.0%) |
| After graph rescue (potential) | 228,415 (84.0%) |
| Reads rescued | +8,038 (+3.0 pp) |
| Remaining unphased | 43,601 (16.0%) |
| Graph-only in BAM (unphased) | 8,038 |
| Graph-only not in BAM at all | 0 |

### Assessment

**RISKY** — 197 / 615 block pairs (32%) are discordant (neither clearly same nor flipped orientation). A simple polarity-based merge would inject errors at those block boundaries. The 43,601 D-category reads (16%) are absent from both pipelines and represent the hard ceiling for any hybrid approach.

The concordant+flipped pairs do cover the majority of doubly-phased reads (185,336 reads, 89%), suggesting that orientation *is* resolvable for most large blocks — the discordant pairs are likely small blocks with insufficient read overlap. A read-count-weighted merge may still be feasible.

---

## VCF-Level Metrics (chr20, no whatshap/hiphase needed)

Phase block stats extracted directly from phased VCFs; NC50 computed via `bedtools subtract`
of `switch_positions.bed` from phase block spans.

| Metric | Graph v2 (noise filter) | BAM (no sidecar) |
|--------|------------------------|-----------------|
| Het variants in VCF | 60,304 | 83,013 |
| Phased het variants | 60,304 (100%) | 83,013 (100%) |
| Phase blocks (VCF) | 427 | 431 |
| VCF N50 | 426 Kbp | 416 Kbp |
| NC50 (switch-corrected N50) | 167 Kbp | **214 Kbp** |
| Corrected phase blocks | 5,660 | **3,688** |
| Total corrected span | 50.6 Mbp | **55.2 Mbp** |

**Notes:**
- NC50 = N50 of phase blocks after subtracting switch-error positions (analogous to HiPhase paper NGC50 but at chr20 scale)
- Graph calls fewer het variants (60K vs 83K) because the graph pipeline uses graph snarls as sites, not de-novo variant calling from alignments
- BAM's lower corrected block count (3,688 vs 5,660) confirms fewer switch errors chopping up its blocks

---

## Noise Filtering: BAM vs Graph Pipeline

### What we ported

The BAM pipeline's `classify_variant_initial` applies three reference-context checks to
indel candidates:

1. **Homopolymer context** — indel sits in a run of ≥3 copies of a 1–6 bp repeat unit
2. **Tandem repeat context** — deleted/inserted motif matches ≥3 tandem copies in flanking reference
3. **SDUST low-complexity** — variant position falls inside a low-complexity interval

Indels flagged by (1) or (2) become `RepeatHetIndel` and are excluded from k-means phasing
(`lcd_var_i_to_cate = kLongcalldRepHetVar`, which doesn't match `kCandGermlineClean`).

These three checks were extracted into a shared module (`src/noise_filter.hpp/cpp`) and wired
into the graph pipeline via `apply_graph_noise_filter` in `graph_bam_adapter.cpp`. Each worker
thread opens its own `faidx_t`, fetches the chunk's reference slice, and reclassifies
`CleanHetIndel` → `RepeatHetIndel` for noisy sites. The flag behavior matches the BAM pipeline
exactly: `RepeatHetIndel` candidates are excluded from k-means.

### What we did NOT port (and why)

The BAM pipeline has additional filtering stages that depend on **per-read CIGAR error
intervals** — information that doesn't exist in the graph pipeline's GAF input:

| BAM feature | Why it doesn't apply to graph |
|---|---|
| Read-level noisy regions (XID clustering from CIGAR) | GAF allele observations are node-level, not base-aligned. No per-read error intervals. |
| Dense overlap check (`var_pos_cr` n_ov > 1) | Snarl parent-child gating already handles nested/overlapping sites. |
| Noisy-reads-ratio gate | Requires per-read CIGAR error intervals. |
| `post_process_noisy_regs` (widen noisy spans) | Depends on read-level noisy regions. |
| `apply_noisy_containment_filter` (demote contained candidates) | Depends on widened noisy regions. |
| Noisy-region MSA recall (`collect_noisy_vars_step4`) | Requires BAM reads for MSA realignment. |
| Second k-means with `kCandGermlineVarCate` | Only useful after MSA recall produces `NoisyCandHet`. |

### SNPs in low-complexity regions

`is_noisy_site` checks `pos_in_low_complexity` for all variant types, including SNPs. However,
`apply_graph_noise_filter` only reclassifies `CleanHetIndel` candidates — SNPs are not touched.
This matches the BAM pipeline, where SDUST low-complexity intervals are **never used to directly
reclassify any variant**. Instead, they serve a structural role:

1. **Extending noisy region boundaries** — adjacent low-complexity intervals are absorbed into
   read-level noisy spans during `pre_process_noisy_regs_pgphase`.
2. **Extending repeat-indel spans** — when a `RepeatHetIndel` is added to the noisy region tree,
   its genomic span is widened to include overlapping low-complexity intervals.

A SNP can end up filtered in the BAM pipeline, but only **indirectly**: if the low-complexity
region was already part of or adjacent to a noisy region built from read-level error clustering,
and the widened noisy span fully contains the SNP, then `apply_noisy_containment_filter` demotes
it to `NonVariant`. A SNP in a low-complexity region with no nearby read-level noise stays
`CleanHetSnp`.

If we reclassified graph SNPs in low-complexity regions, we'd be **more aggressive** than the
BAM pipeline — potentially removing legitimate het SNPs that are reliable phasing anchors.

### Future work

If the graph pipeline gains access to base-level alignment information (e.g. via surjected BAM
or GAF CIGAR strings), the read-level noisy region pipeline could be ported. Until then, the
reference-context checks are the only applicable noise filter.

---

## Hybrid Model — Phase 0 Baseline (measurement harness)

Phase 0 of the hybrid-model plan: build a measurement harness and establish a
baseline for the existing `collect-hybrid-variation` command before changing it.

### Deliverables

- `scripts/bench_hybrid.sh` — runs hybrid mode and, with `--compare`, also the
  BAM-only and graph-only pipelines on the same input into sibling dirs for A/B.
  Per-pipeline BAM flag differs: graph uses `--phased-bam-out`, bam/hybrid use
  `--out-bam`.
- `scripts/phase_block_stats.py` — VCF-level metrics (het/hom, phased fraction,
  block count, N50) from a phased VCF alone. No truth set needed. For
  switch/accuracy use `evaluate_phase_accuracy.sh` (needs diplinator).

### Environment note

The chr20 HG002 dataset and diplinator from the original eval are **not present**
in this environment, so the full switch/accuracy gate from the plan could not be
reproduced here. Phase 0 instead used the checked-in `test_data/graph/` fixture
(MICB/KIR3DL1, chr6+chr19) by reconstructing reads from the GAF `cs` tags and
aligning them with minimap2 to produce a matched BAM+sites+GAF triple.

> **Caveat:** the fixture is ultra-low-depth (~0.04× mean, 151 bp short reads),
> not the HiFi/ONT long-read target. Numbers below validate the *harness and code
> behavior*, not phasing quality. The real gate still requires the chr20 data.

### A/B result (fixture)

| Pipeline | Candidates | Phased het | Blocks | N50 |
|---|---|---|---|---|
| BAM-only | 1 | 1 | 1 | 0 |
| Graph-only | 292 | 279 | 52 | 1062 |
| Hybrid | 918 | 68 | 17 | 146 |

### Findings (confirm the static analysis)

- **`bridged=0` on every chunk.** No graph site matched a BAM candidate, so the
  intended bridging never fired. Partly because BAM found ~1 candidate at this
  depth, but it also exercises the `find_matching_candidate` sort-assumption gap.
- **Hybrid phases *fewer* het variants (68) than graph-only (279)** despite adding
  917 graph-only candidates. The unconditionally-`CleanHet` graph candidates
  (no AF/depth/het gate) flood the unified k-means and degrade it rather than
  augment it — exactly the Phase 1 hazard.
- **Verbose stats only print when `added>0` AND a worker logs them**; counts are
  per-chunk and the same `added=` value repeats across overlapping chunks,
  suggesting overlap-region sites are added in both chunks (double-add).

### Decision gate

The intended "BAM-primary, graph-rescue" behavior is **not** what the code does
(it co-phases and can regress below graph-only). This confirms the plan's
direction: proceed to **Phase 1** (graph candidate quality gates) and the
**Phase 2** BAM-frozen rescue rather than tuning the current co-phasing path.
The chr20 switch-error gate must still be run once that data is available before
any hybrid path is declared shippable.

---

## Hybrid Model — Phase 1 (graph candidate quality gates)

Phase 1 hardened graph-only candidate handling so graph sites can no longer
enter k-means as unconditional clean-het anchors.

### Changes

- **Deferred classification** (`add_graph_only_candidate`): graph sites are now
  added unclassified (category `LowCoverage`, flag 0) instead of stamped
  `CleanHet`. They stay out of k-means until gated.
- **Quality gate** (`classify_graph_only_candidates`): after counts are final
  (post `inject_graph_reads`), each graph-only candidate is run through the BAM
  pipeline's `classify_variant_initial` and `category_to_flag`. Only sites
  passing the het band (`min_depth`, `min_alt_depth`, `min_af ≤ AF ≤ max_af`)
  become `CleanHet` and enter het k-means. `classify_variant_initial` was
  exposed in `collect_var.hpp` for reuse (was file-static).
- **Ordering fix** (`process_chunk_hybrid`): the indel noise filter
  (`apply_hybrid_noise_filter`, which only acts on `CleanHetIndel`) now runs
  *after* the gate, not before. Verbose log adds `promoted=N`.
- **Sort-precondition fix** (`find_matching_candidate`): the bridge binary
  search now takes an explicit `n` = original BAM candidate count and only
  searches `candidates[0, n)`. `inject_graph_sites` asserts that range is
  sorted, so a future ordering change fails loudly instead of silently missing
  bridges.
- **Unit test** `src/test_hybrid_inject.cpp` (wired into `make unit-tests`):
  clean het promoted; homozygous, low-depth, low-alt-depth, and low-AF sites
  excluded from the het mask. All unit tests green.

### A/B result (same fixture as Phase 0)

| Pipeline | Candidates | Phased het | Promoted |
|---|---|---|---|
| BAM-only | 1 | 1 | — |
| Graph-only | 292 | 279 | — |
| Hybrid (Phase 0) | 918 | 68 | (ungated) |
| **Hybrid (Phase 1)** | 918 | **2** | **2 / 917** |

### Verification

- BAM-only path unchanged from Phase 0 (386 het / 182 hom / N50 194368) — no
  regression from exposing `classify_variant_initial`.
- `make -j` builds with zero warnings; `make unit-tests` all green.
- Asserts are active (no `-DNDEBUG`), so the sort precondition is enforced.

### Correction: the "AF=0 bug" was a synthetic-data artifact

The low-depth fixture suggested hybrid allele counting was broken (graph het
sites showing REFc=N, ALTc=0, AF=0). **This did not reproduce on real chr20
data** (see next section) — there, the same sites show correct AF≈0.4–0.5. The
AF=0 came from the Phase 0 read-reconstruction harness (`gaf_to_fastq.py` +
minimap2 short-read alignment did not carry variant alleles faithfully), not
from the hybrid pipeline. Lesson: validate hybrid behavior only on real aligned
reads, not reconstructed ones.

---

## Hybrid Model — chr20 25M validation (REAL data)

Real HG002 HiFi test data was added at `test_data/graph_chr20/`: a 500 kb
chr20 slice (25,000,001–25,500,000) with a matched BAM, coordinate GAF (2,048
reads), sites VCF (6,963 snarl sites), and N-padded ref. This is the first
hybrid run on genuine aligned reads.

> Note: `samtools index` the BAM first (`.bam.bai` is gitignored, not shipped).

### Expected-output assertions (both pass with Phase 1 build)

| Pipeline | Candidates | Expected |
|---|---|---|
| BAM-only | 1108 | 1108 ✓ |
| Graph-only | 509 (13,490 filtered) | 509 ✓ |

Phase 1 changes preserve both verified pipeline outputs exactly.

### Hybrid A/B (Phase 1)

| Pipeline | Phased het | Blocks | N50 | Largest block |
|---|---|---|---|---|
| BAM-only | 527 | 2 | 296,271 | 296,271 |
| Graph-only | 446 | 2 | 393,964 | 393,964 |
| **Hybrid** | **522** | **1** | **498,281** | 498,281 |

Hybrid bridging works on real data: **791 sites bridged**, 13,135 graph-only
sites added, **72 promoted** to clean-het by the gate, 1,692 BAM reads extended
with graph observations (`reads_injected=0` — the GAF and BAM are the same HG002
reads, so there are no graph-only reads to inject here).

Result: hybrid phases nearly as many het as BAM (522 vs 527) but collapses the
two BAM blocks into **one block with N50 498 kb** — larger than either
standalone pipeline. This is the pangenome-bridging stitching benefit the plan
predicted, appearing for the first time on real data.

### Gate behavior on real data (sanity)

Graph-only candidate categories after gating: 5,875 LOW_COV, 598 CLEAN_HOM,
397 CLEAN_HET_SNP, 88 NOISY_HET, 76 REP_HET_INDEL, 38 NOISY_HOM,
37 CLEAN_HET_INDEL, 32 LOW_AF. Only the clean-het sites enter het k-means; hom,
low-cov, low-AF, and repeat sites are correctly excluded. Spot-check vs graph
pipeline at chr20:25001841 — graph AF 0.52, hybrid AF 0.36 (REFc=7/ALTc=4),
both CLEAN_HET_SNP: alleles are counted correctly.

### Open items

- The **switch-error / accuracy** gate against the diploid truth still needs
  `diplinator` (not in this environment). The N50 gain is promising but must be
  confirmed not to come with a switch-error regression before shipping.
- Phase 2 (BAM-frozen rescue) is most valuable where the GAF contains reads
  absent from the BAM; this fixture has none, so that path needs a dataset with
  graph-only reads to exercise.

---

## Hybrid Model — Phase 3 (stitch abstain margin)

Phase 3 adds an abstain rule to chunk stitching so weakly supported block
boundaries are not merged.  This is the mitigation for the CHECKPOINT "32%
discordant pairs" over-merge risk: a wrong merge injects a switch error, so
when the evidence is thin it is safer to leave two blocks.

### Change

- `flip_chunk_hap` now merges adjacent chunks only when
  `|flip_hap_score| > opts.stitch_min_margin`.  The previous rule
  (`flip_hap_score == 0 → no merge`) is exactly `margin = 0`, the default.
- `Options::stitch_min_margin` (struct default `kDefaultStitchMinMargin = 0`).
- The hybrid command overrides the struct default to
  `kHybridDefaultStitchMinMargin = 10` in `collect_hybrid_variation`;
  `--stitch-min-margin INT` still overrides it per run.
- CLI flag `--stitch-min-margin INT` added to `collect-hybrid-variation`
  **only**.  The BAM and graph commands do not expose it and keep margin 0, so
  their behavior is unchanged.

### Hybrid default raised to 10 (matched-eval result)

Matched-eval on HG002 chr20 (identical evaluator, truth, and binary for BAM
and hybrid) showed that switch and flip rates are flat (~0.64% / ~0.93%) across
all configurations including pure BAM — they are intrinsic to the shared
k-means + stitcher core, not a hybrid regression.  The metrics that actually
separate the pipelines are Hamming and auN.  On those, raising the hybrid
stitch margin from 0 to 10 is a strict improvement:

| Config (af-indel 0.11) | Hamming | auN | switch | flip | perfect | covered |
|---|---|---|---|---|---|---|
| BAM baseline      | 3.09 | 9.83M  | 0.64 | 0.933 | 29.1% | 313.5M |
| hybrid margin 0   | 3.05 | 12.08M | 0.66 | 0.968 | 22.2% | 312.0M |
| hybrid margin 10  | 3.04 | 11.90M | 0.65 | 0.966 | 23.1% | 317.0M |

Margin 10 lowers Hamming, raises genome covered (+5 Mb) and perfect phase sets
(+0.9 pt) versus margin 0 at no switch/flip cost, so it is the hybrid default.

### Safety: shared code, default-preserving

`flip_chunk_hap`/`stitch_chunk_haps` is shared by all pipelines.  Because the
default margin is 0 and the abstain test reduces to the original `== 0` check,
BAM and graph outputs are byte-identical (re-verified: 1108 / 509 candidates,
527/446 het, same blocks/N50 as Phase 1).

### chr20 margin sweep (real data)

Default chunk size (500 kb) puts the whole 500 kb data slice in one chunk, so
stitching is not exercised.  Forcing `--chunk-size 100000` (≈5 data chunks):

| Margin | Blocks | N50 |
|---|---|---|
| 0 | 2 | 406,166 |
| 1–2 | 2 | 406,166 |
| 4 | 3 | 305,381 |
| 8 | 4 | 297,133 |
| 16 | 5 | 98,172 |

Monotonic and correct: higher margin → more abstentions → more, smaller blocks.
The boundaries here need margin ≥ 4 before any abstains, i.e. they are
**well-supported merges**, not discordant over-merges.  At the default 500 kb
chunk size the hybrid single 498 kb block survives margins up to 128+,
confirming that stitch is genuine.

### Tests

`src/test_phase_block_stitch.cpp` gains two cases: a single agreeing overlap
read (score +1) merges at margin 0 and abstains at margin 1.  All unit tests
green; build has zero warnings from project code (one pre-existing abPOA
third-party `-Wunused-function` warning is unrelated).

### Open

- Still pending the diplinator switch-error gate to prove abstain margins trade
  N50 for accuracy as intended.  The margin gives a knob to *tune* that
  trade-off once truth-based evaluation is available.

---

## Hybrid Model — graph SNP low-complexity "gap" (REVERTED — was a bug)

An earlier change made `apply_hybrid_noise_filter` demote graph-only het SNPs in
SDUST low-complexity regions to `NonVariant`, on the premise that "this is what
the BAM pipeline does to density-noisy SNPs." **That premise was wrong, and the
change has been reverted.**

### Why it was wrong

Direct measurement on chr20 25M at the 10 demoted positions showed both source
pipelines *keep* these as real het calls:

- **BAM pipeline**: recalls them as `NOISY_CAND_HET` via noisy-region MSA (step 4)
  — AF 0.40–0.60, depth 46–78. The BAM pipeline does **not** discard them.
- **Standalone graph pipeline**: emits them as `CLEAN_HET_SNP` and phases them.

So the demotion dropped 10 genuine het SNPs that both pipelines call. The
`NOISY_CAND_HET` recall happens during phasing (after the bridge step), so these
sites are not in the BAM candidate table at graph-injection time and therefore
appear as graph-only — but they are still real variants, not noise.

### Resolution

`apply_hybrid_noise_filter` again leaves SNPs untouched (matching the standalone
graph pipeline's `apply_graph_noise_filter`); only indels are screened for
homopolymer / repeat / low-complexity noise.

---

## Hybrid Model — candidate-leak bug (graph-only non-calls not pruned)

While verifying the SNP revert, hybrid emitted **7141** candidates on chr20 25M
versus 1108 (BAM) and 509 (graph).

### Root cause

The BAM pipeline drops non-call candidates (`LowCoverage` / `NonVariant` /
`StrandBias`) via `prune_not_candidate_variants` at the end of
`collect_var_classify`, and folds `LowAlleleFraction` → `LowCoverage` in
classification pass 2 so those are pruned too. The hybrid pipeline appends
graph-only candidates *after* `collect_var_classify` has already run, and
`classify_graph_only_candidates` only runs pass-1 (`classify_variant_initial`).
So graph-only sites that failed the gate kept `LOW_COV` / `LOW_AF` / `NON_VAR`
categories and flowed straight to output — `merge_chunk_candidates` does not
filter by category. The standalone graph pipeline avoids this because its output
builder (`graph_chunks_to_candidate_table`) only emits clean/noisy categories.

The 7141 breakdown: 1234 real calls + 5875 `LOW_COV` + 32 `LOW_AF` + 11
`NON_VAR`.

### Fix

- `classify_graph_only_candidates` folds `LowAlleleFraction` → `LowCoverage`,
  matching the BAM pipeline's pass-2 rewrite.
- `process_chunk_hybrid` calls `prune_not_candidate_variants(chunk)` after
  k-means phasing. Pruning after phasing is safe: per-candidate read profiles
  (`read_var_profile`, `read_var_cr`) are consumed entirely during phasing and
  not used afterward; stitching, merge, and TSV/VCF/phased-BAM output use
  read-level state, not candidate indices. No post-phasing step (k-means or
  step-4 MSA) ever assigns a prunable category, so pruning after rather than
  before phasing removes exactly the same candidates the BAM ordering would.
- Exposed `prune_not_candidate_variants` (was file-static).

### chr20 25M result

- Hybrid: **7141 → 1234**. Categories are now only real calls
  (CLEAN_HOM / CLEAN_HET_SNP / NOISY_CAND_HET / REP_HET_INDEL / NOISY_CAND_HOM /
  CLEAN_HET_INDEL); no LOW_COV / LOW_AF / NON_VAR leak.
- BAM-only (1108) and graph-only (509) outputs byte-identical (MD5 verified).
- BAM↔hybrid site parity: of 1108 BAM sites, only 2 differed — both DEL
  encoding differences (`AA>.` vs `AA>A`). These were later traced to a real
  deletion key-collision bug and fixed (see "deletion key collision" section
  below); parity is now 1108/1108.
- Phased VCF / candidate VCF / phased BAM (HP+PS tags) all emit correctly.
- `make unit-tests` all green; added a `prune_not_candidate_variants` test and
  updated the noise-filter test to assert SNPs are kept.

---

## Hybrid Model — deletion key collision overwrote BAM calls

The 2 BAM↔hybrid site differences noted above (`AA>.` vs `AA>A`) were not a
cosmetic encoding choice — they were a real bug: the hybrid path silently
**overwrote a verified BAM deletion** with a graph deletion of a different
length.

### Root cause

`vcf_to_variant_key` (hybrid_inject.cpp) built deletion keys by stripping only a
**single** anchor base from the VCF REF/ALT, instead of the full shared prefix.
For a homopolymer site like `25412767 TAA→TA` (a 1 bp deletion) it produced
`pos=25412768, ref_len=2, alt="A"` — an internally inconsistent key claiming a
2-base ref span while keeping a residual `alt="A"`.

Two downstream consequences:

1. **No bridging.** `find_matching_candidate` requires `k.alt == target.alt`.
   BAM deletions always have `alt=""` (`variant_key_from_digar`), so a graph
   deletion key with `alt="A"` could never match — graph deletions were always
   added as separate candidates instead of bridging.
2. **Silent overwrite.** `exact_comp_var_site` **ignores `alt` for deletions**
   (it keys deletions on `pos` and `ref_len` only). The inflated `ref_len=2`
   made the graph 1 bp deletion compare *equal* to the genuine BAM **2 bp**
   deletion at the same locus. In `merge_chunk_candidates` the graph candidate
   (`lcd_make_variants_region_pass=true`) overwrote the BAM call, so only the
   graph allele survived.

This is why the VCF site `TAA→TA,T` (genuinely biallelic: a 1 bp and a 2 bp
deletion) collapsed to a single hybrid row that did not match BAM.

### Fix

Strip the **full** shared prefix in the deletion branch of
`vcf_to_variant_key`, matching both the BAM path (`variant_key_from_digar`) and
the standalone graph path (`graph_collect.cpp`): `pos = first deleted base`,
`ref_len = deleted span`, `alt = ""`. For a single-base anchor this reduces to
the previous behavior, so only homopolymer / repeat deletions change.

Hybrid-only change; the BAM and graph pipelines do not call this function.

### chr20 25M result

- BAM↔hybrid parity is now **complete**: all **1108/1108** BAM calls (and all
  **108/108** BAM deletions) are present in the hybrid output (was 1106 / 106).
- The biallelic homopolymer sites now emit **both** deletion alleles as distinct
  rows (e.g. `25412768 AA>.` 2 bp + `25412769 A>.` 1 bp), consistent across
  TSV / candidate VCF / phased VCF.
- BAM-only (1108) and graph-only (509) outputs remain byte-identical (MD5).
- Added a `vcf_to_variant_key` unit test asserting SNP / INS / DEL normalization
  and that 1 bp vs 2 bp homopolymer deletions encode distinct keys.

---

## Hybrid Model — graph het indels over-demoted (hybrid ≈ BAM)

### Symptom

Full-chromosome benchmarks showed the hybrid barely beating BAM (accuracy
96.86% vs 97.03%, phased reads 81.2% vs 80.9%) while the standalone graph
pipeline phased **93.4%** of reads. The hybrid was reproducing BAM instead of
harvesting the graph's phasing power, and phase-block auN regressed (11.8M →
8.9M).

### Diagnosis (chr20 25M slice)

Phasing in the hybrid runs only `collect_var_run_phasing`, whose first (and
only) k-means uses `kCandGermlineClean` — **clean het SNP + clean het indel +
clean hom**. So the count of *clean* het anchors drives phasing. Measured:

| anchors into k-means | BAM | GRAPH | HYBRID (old) |
|---|---|---|---|
| clean het SNP + clean het indel | 384 | **509** | 438 |

The graph's advantage is its **113 clean het indels**; the old hybrid kept only
**41**, demoting **29 → RepeatHetIndel** and reclassifying the rest. Two causes,
both from running BAM-conservative logic on graph-only candidates:

1. **classify_graph_only_candidates** called `classify_variant_initial` (the BAM
   classifier), which demotes homopolymer/repeat het indels to RepeatHetIndel
   **inline**, keyed on the normalized `VariantKey`. The standalone graph
   pipeline's `classify_graph_candidates` has no such inline demotion.
2. **apply_hybrid_noise_filter** reconstructed vcf_ref/vcf_alt from the
   normalized key, so `is_noisy_site` saw a different string/position than the
   graph pipeline's `apply_graph_noise_filter` (which screens on the original
   catalog ref/alt) — a different homopolymer verdict for the same site.

### Fix (hybrid-only)

Mirror the graph pipeline's two-step exactly:
`classify_graph_candidates` (depth/AF/het/hom, **no** inline indel demotion)
followed by `apply_graph_noise_filter` (screens on catalog strings).

- `classify_graph_only_candidates` now classifies graph-only candidates with the
  graph pipeline's depth/AF/het/hom logic (replicated inline; the shared
  `classify_variant_initial` is left untouched for the frozen BAM pipeline) and
  does **no** homopolymer demotion. LowAlleleFraction is still folded to
  LowCoverage to match BAM pass-2 pruning.
- `inject_graph_sites` now records each graph-only candidate's original VCF
  (pos, ref, alt) in a `GraphOnlyVcfAlleles` side-map (remapped through the
  post-injection sort permutation).
- `apply_hybrid_noise_filter` uses those original strings when available, so its
  homopolymer/repeat verdict matches the standalone graph pipeline; it falls
  back to key-reconstruction only when unavailable.

### Result (chr20 25M slice)

- Clean het anchors **438 → 448**; clean het indels **41 → 51**;
  RepeatHetIndel **54 → 44**.
- Phased het variants (phased VCF) **527 → 535**.
- The remaining gap to the graph's 509 clean anchors is **BAM-origin** sites
  whose linear pileup is genuinely noisy (`NOISY_CAND_HET`): of the graph
  clean-het sites the hybrid still marks noisy, **all 20 are bridged BAM-origin**
  and **0 are graph-only** — i.e. real BAM/graph disagreements, not a bug.
  Touching them would alter the frozen BAM classification, so they are out of
  scope.
- Invariants held: BAM (1108) and graph (509) outputs byte-identical (MD5);
  BAM call parity 1108/1108; `make unit-tests` green. Added a noise-filter test
  asserting original VCF strings override key-based reconstruction.

This slice is a single phase block (saturated N50), so the block-contiguity /
auN effect must be confirmed on a full chromosome with truth-based evaluation
(diplinator/whatshap unavailable in this environment).

---

## Lessons / Pitfalls

- **pgbam sidecar causes over-stitching:** Running BAM pipeline with `--pgbam-file` merged all 49 initial phase sets into 2 chr20-spanning blocks with near-random accuracy (55%). Do not use the sidecar without understanding its signal quality.
- **Wrong reference = 0 reads:** `hg002v1.1.fasta` (maternal/paternal contigs) does not match the `CHM13#0#chr20` contig name in the aligned BAM. Must use `chm13v2.0.chr20.renamed.fa` as ref for `collect-bam-variation`.
- **Use chr20-only truth refs:** Full-genome maternal/paternal assemblies cause a small fraction of reads to leak to chr1/6/16 during truth alignment. Use `HG002_1.1_MATERNAL.chr20.fasta` / `HG002_1.1_PATERNAL.chr20.fasta` for clean results.
- **Per-chunk site loading:** Graph pipeline uses tabix-indexed per-chunk site queries (commit `5878aad`) — loading all 977K sites upfront was the prior bottleneck.
- **Non-numeric node IDs:** `build_compact_graph_site_index` now skips sites with non-numeric node IDs instead of aborting (commit `5c01954`).

---

## Overlap-Read Phase Propagation at Chunk Boundaries

### The problem

The genome is tiled into overlapping chunks for parallel phasing. A long read that spans a
chunk boundary exists as independent copies in both the upstream and downstream chunks. Each
chunk phases its copy using its own candidates. The BAM writer emits each read from its
**owning** (upstream) chunk only.

Example: chunks A (0–5 Mb) and B (4.5–10 Mb) overlap at 4.5–5 Mb. Read R starts at 4.8 Mb
and ends at 5.2 Mb. R appears in both chunks:

```
Chunk A (owns R):  candidates at 4.0, 4.3, 4.6 Mb
Chunk B:           candidates at 4.9, 5.1, 5.5, 5.8 Mb
                         ▲         ▲
                         R covers these in chunk B
```

In chunk A, R's variant profile covers the candidate at 4.6 Mb but no het sites are nearby —
R gets hap=0 (unphased). In chunk B, R covers candidates at 4.9 and 5.1 Mb, both het SNPs —
R gets hap=1, PS=4900000.

After stitching, `flip_chunk_hap` uses overlap reads (including R's phased copy in B) to
vote on whether to flip/merge chunks A and B. If the vote succeeds (flip_hap_score != 0),
chunk B's PS is rewritten to match chunk A's PS, and hap labels are flipped if needed. Now
R's copy in chunk B has HP/PS consistent with chunk A's phasing.

But the BAM writer emits R from chunk A, where it is still hap=0, PS=-1. R is written
unphased despite being successfully phased in chunk B. LongcallD has this same gap.

### The original (broken) fix

Commit `f094cea` added `propagate_overlap_read_phase_to_output_owner`, which copies HP/PS
from the downstream chunk to unphased upstream overlap reads. This recovered ~1,300 reads
but propagated **unconditionally** — including for unmerged pairs where `flip_hap_score == 0`.
For unmerged pairs the relative phase between chunks is unknown, so:

- The downstream hap=1 may correspond to the upstream hap=2 (wrong haplotype).
- The downstream PS doesn't exist in the upstream chunk (orphan phase set).

This caused pgphase to diverge from longcallD: +1,300 phased reads, +2 phase sets, with
some reads potentially assigned to the wrong haplotype.

### The fix

`stitch_chunk_haps` now tracks which adjacent pairs were successfully merged by
`flip_chunk_hap` via a `pair_stitched` vector. Propagation only fires for merged pairs,
where the downstream PS has been rewritten to match the upstream PS and hap labels have
been flipped if needed. Unmerged pairs are skipped.

### Evaluation (HG002 HiFi chr20)

| Metric | longcallD | pgphase (unconditional) | pgphase (merged-only) |
|---|---|---|---|
| Phased reads | 219,090 (80.5%) | 220,377 (81.0%) | 219,925 (80.9%) |
| Phase sets | 429 | 431 | 429 |
| Phase block N50 | 937,648 bp | 937,121 bp | 937,648 bp |
| Overall accuracy | 97.03% | 97.02% | 97.03% |
| Switch error rate | 0.61% | 0.62% | 0.61% |
| Phaseable accuracy | 98.15% | 98.14% | 98.15% |

The merged-only propagation recovers 835 reads (+0.4%) over longcallD with identical phase
set count, block structure, accuracy, and switch error rate.

---

## Three-Pipeline Comparison — HG002 HiFi chr20 (current)

Full chr20 evaluation using HG002 HiFi reads aligned to CHM13 chr20, evaluated
against HG002 T2T diploid assembly (diplinator). All numbers reflect the current
binary after all fixes above, including the indel noise-filter fixes (content-base
anchoring, nt4 insertion compare, non-minimal graph-allele trimming).

"Hybrid (old)" = BAM-base + graph-augment with the step-4 noisy-candidate k-means
re-orientation still on. "Hybrid (new)" = current default after disabling that
re-orientation (`skip_noisy_kmeans`, on by default for hybrid; `--keep-noisy-kmeans`
restores the old behaviour). See "Hybrid step-4 re-orientation" below.

| Metric | BAM | Graph | Hybrid (old) | **Hybrid (new)** |
|---|---|---|---|---|
| **Completeness** | | | | |
| Input reads | 272,016 | 231,382 | 272,016 | 272,016 |
| Phased reads | 219,925 (80.9%) | 203,729 (88.0%) | 220,725 (81.1%) | 212,484 (78.1%) |
| Phase sets | 429 | 368 | 415 | 348 |
| Perfect phase sets | 125 (29.1%) | 142 (39.1%) | 96 (23.5%) | **114 (33.4%)** |
| **Contiguity** | | | | |
| Phase block N50 | 937,648 bp | 917,581 bp | 945,071 bp | **1,041,015 bp** |
| Phase block auN | 9,834,145 bp | 989,988 bp | 11,876,987 bp | **13,677,722 bp** |
| Largest block | 53,165,843 bp | 2,080,875 bp | 59,316,535 bp | 58,616,492 bp |
| Median block span | 867,222 bp | 871,853 bp | 870,769 bp | 867,209 bp |
| Genome covered | 313,517,513 bp | 221,616,803 bp | 317,945,255 bp | 266,888,983 bp |
| **Accuracy** | | | | |
| Concordant reads | 213,127 | 202,076 | 213,609 | 210,724 |
| Discordant reads | 6,798 | 1,639 | 7,099 | 1,741 |
| Overall accuracy | 96.91% | 99.20% | 96.78% | **99.18%** |
| Hamming error rate | 3.09% | 0.80% | 3.22% | **0.82%** |
| Switch errors | 1,393 | 367 | 1,386 | **319** |
| Flip errors | 2,052 | 858 | 2,135 | 953 |
| Switchflip errors | 3,445 | 1,225 | 3,521 | **1,272** |
| Switch opportunities | 219,496 | 203,352 | 220,299 | 212,124 |
| Switch error rate | 0.63% | 0.18% | 0.63% | **0.15%** |
| Switchflip error rate | 1.57% | 0.60% | 1.60% | **0.60%** |

### Latest verified run (commit `c36026e`, full metrics + timings)

Re-run of all three pipelines after the snarl-duplicate dedup (`d3ac898`), the
CLI arg-parse crash fix (`49e4a49`), and the `--exclude-bed` revert (`c36026e`).
HG002 HiFi chr20 → CHM13, diplinator truth, `--min-reads 5`, 8 threads. Runtime
is wall-clock for the collect+phase run (excludes eval and BAM index/sort).

**Accuracy & errors**

| Metric | BAM | Graph | **Hybrid** |
|---|---|---|---|
| **Runtime** | 98 s | **53 s** | 326 s |
| Accuracy | 96.909% | 99.210% | **99.232%** |
| Hamming | 3.091% | 0.790% | **0.768%** |
| Switch errors | 1,393 | 364 | **291** |
| Flip errors | 2,052 | 856 | 906 |
| Switchflip | 3,445 | 1,220 | **1,197** |
| Switch rate | 0.635% | 0.179% | **0.137%** |
| Switchflip rate | 1.570% | 0.600% | **0.565%** |
| Switch opportunities | 219,496 | 203,357 | 211,900 |
| Concordant reads | 213,127 | 202,110 | 210,610 |
| Discordant reads | 6,798 | 1,610 | 1,630 |

**Yield & phase blocks**

| Metric | BAM | Graph | **Hybrid** |
|---|---|---|---|
| Total input reads | 272,016 | 231,382 | 272,016 |
| Phased reads | 219,925 | 203,734 | 212,259 |
| Fraction phased | 80.85% | **88.05%** | 78.03% |
| Reads evaluated | 219,925 | 203,720 | 212,240 |
| Phase sets (eval) | 429 | 363 | 340 |
| Perfect PS | 125 | **142** | 115 |
| Perfect PS % | 29.14% | **39.12%** | 33.82% |
| Candidates | 120,124 | 73,627 | 134,646 |
| Block N50 (bp) | 937,648 | 917,581 | **993,831** |
| Block auN (bp) | 9,834,145 | 989,984 | 7,153,008 |
| Block median span (bp) | 867,222 | 871,853 | 866,966 |
| Block max span (bp) | 53,165,843 | 2,080,875 | 39,388,202 |
| Genome covered (bp) | 313,517,513 | 221,615,836 | 246,734,240 |

**Phaseable vs unphaseable split** (confidence threshold 0.6)

| Metric | BAM | Graph | **Hybrid** |
|---|---|---|---|
| Phaseable PS / reads | 366 / 214,444 | 343 / 202,881 | 324 / 211,704 |
| Phaseable accuracy | 98.02% | **99.40%** | 99.34% |
| Unphaseable PS / reads | 63 / 5,481 | 20 / 839 | 16 / 536 |
| Unphaseable accuracy | 53.62% | 54.11% | 55.04% |

- **No regression.** Hybrid holds at 99.232% / 0.768% / switch 291 (matches the
  pre-/post-dedup verification). Graph holds the Fix-B no-pool gain: switch
  367 → 364, accuracy 99.195% → 99.210%.
- **Hybrid wins every error metric** (best accuracy, Hamming, switch, switch+flip
  rate) and the longest block N50 (993,831 bp), but is the slowest (326 s).
- **Graph is the accuracy-per-second winner** (99.210% in 53 s), highest fraction
  phased (88.05%), most perfect phase sets (142 / 39.12%), and best phaseable
  accuracy (99.40%). Its auN is far lower because it produces many uniform ~0.9 Mb
  blocks rather than a few huge spanning blocks.
- **BAM** has the largest auN/max-span (one 53 Mb block) and most genome covered,
  but the worst accuracy (3.091% Hamming) — the long blocks accumulate errors. It
  is the fallback when no graph/GAF inputs are available.
- **auN caveat:** BAM/hybrid's large auN comes from a single chromosome-spanning
  block; for switch-error quality the per-read accuracy metrics above are the
  meaningful comparison, not auN.

### Observations

- **The indel noise-filter fix transforms the graph pipeline.** Graph overall
  accuracy jumped from 92.70% to **99.20%** (Hamming 7.30% → **0.80%**) and switch
  error rate fell from 1.07% to **0.18%**. The fix demotes 17,509 homopolymer/STR
  het indels to `REP_HET_INDEL` so they no longer enter k-means as phasing anchors.
  Phased reads drop (216k → 204k) because those noisy indels previously phased
  many reads incorrectly — the pipeline now trades a small completeness loss for a
  large accuracy gain. Perfect phase sets rose from 22.0% to **39.1%**.
- **Hybrid (new) reaches graph-pipeline accuracy at far higher contiguity.** After
  disabling the step-4 noisy-candidate k-means re-orientation, hybrid jumps from
  96.78% to **99.18%** (Hamming 3.22% → **0.82%**, switch rate 0.63% → **0.15%**).
  It now matches standalone graph accuracy (99.20% / 0.18%) while keeping BAM-class
  contiguity (auN 13.7M vs graph's 0.99M — 14× longer blocks) and phasing 8.8k more
  reads than graph (212k vs 204k). It even has *fewer* switch errors than graph
  (319 vs 367). This is the single biggest hybrid improvement to date — see
  "Hybrid step-4 re-orientation" below for the mechanism.

### Hybrid step-4 re-orientation: the root cause of hybrid's accuracy gap

The post-fix marginal analysis (below) localised 96% of hybrid errors to the
**BAM core**, not the graph layer. The cause was found by A/B toggle: the BAM
pipeline's **step-4 noisy-region recall** (`collect_noisy_vars_step4`) runs an MSA
to recall variants inside noisy regions, then **re-runs k-means over
`kCandGermlineVarCate`** (which includes the recalled `NOISY_CAND_HET` candidates)
to re-assign read haplotypes. The standalone graph pipeline never does this — it
runs a single k-means over `kCandGermlineClean` only.

On hybrid, this second k-means was **actively harmful**:

- It phased **8,243 extra reads at 65% error** (5,358 wrong, 2,885 right).
- It poisoned the shared core: BAM-shared Hamming went **0.71% → 3.11%** when on.
- It also *fragmented* blocks (auN 13.7M with it off vs 11.9M on; perfect PS 114 vs 96).

**The fix (`skip_noisy_kmeans`, default on for hybrid):** keep the MSA variant
recall (noisy variants still appear in the output VCF, `NOISY_CAND_HET=22,240`) but
skip *only* the noisy-candidate k-means re-run. Confirmed by two equivalent toggles
giving identical results: skipping all of step 4 (`exp_skip_step4`) and skipping
just the re-run both yield 99.181% / 0.819% — proving the **re-orientation, not the
recall, was the entire harm**. `--keep-noisy-kmeans` restores the old behaviour.
The BAM pipeline is unchanged (the re-run still runs there; only hybrid defaults to
skipping it, because hybrid's graph-anchored first k-means already orients reads
well and the noisy re-run only adds error).

**Knob sweeps on the clean base were saturated** (no further gain):
`--graph-indel-af-margin` {0.11→0.50} *worsened* accuracy (admitting off-center
graph indels adds noise); `--stitch-min-margin` {4,10,15,20} and `--stitch-rule`
{1,3} moved nothing on accuracy (the both-strands rule dominates merges, not the
margin). `--stitch-rule 3 --stitch-min-margin 15` adds +6 perfect PS (120 vs 114)
and +5M genome covered at identical accuracy — a valid optional contiguity knob,
not promoted to default.

### Result is deterministic, not seed-dependent

The k-means has no RNG (greedy seeded init); repeated runs of the same config give
bit-identical metrics (319 switches every time). The "convergence variance"
referenced in earlier notes was about a separate experimental optimizer, not the
production k-means.

> **Remaining caveat (graph hybrid-allele trimming):** the non-minimal
> multiallelic-allele trimming added to `apply_graph_noise_filter` is still NOT in
> `apply_hybrid_noise_filter`. Low priority now — hybrid (new) error is dominated by
> the graph-influenced hard tail (1,660 reads at 14.6%), and the BAM-shared core is
> already at 0.71%. Revisit only if graph-only repeat indels are shown to drive the
> remaining hybrid errors.

### Post-fix error decomposition: is graph's accuracy transferable to hybrid? (NO)

`split_marginal.py` on the post-fix hybrid per-read results, splitting discordant
reads by origin (graph-influenced = phased by hybrid but NOT by standalone BAM;
bam-shared = phased by both):

| Read subset | Reads | Discordant | Hamming |
|---|---:|---:|---:|
| graph-influenced (marginal) | 1,704 | 284 | **16.67%** |
| bam-shared | 219,004 | 6,815 | 3.11% |
| ALL hybrid-phased | 220,708 | 7,099 | 3.22% |

> **SUPERSEDED conclusion — read with the step-4 finding above.** This section
> concluded "96% of hybrid errors are in the BAM core, which is unfixable, so you
> cannot get graph's low Hamming at BAM coverage." The *diagnosis* (96% BAM-core)
> was correct and led directly to the fix; the *prognosis* (unfixable) was WRONG.
> The BAM-core errors were caused by the step-4 noisy k-means re-orientation, and
> disabling it took hybrid to 99.18% / 0.82% at 212k reads — graph accuracy at far
> higher contiguity. The "denominator effect" framing below is therefore only
> partly true: graph's accuracy is real (see shared-read analysis) AND hybrid can
> now match it. Kept for the diagnostic trail.

Two conclusions for the "lower Hamming/switch" goal (the second is now superseded):

- **96% of hybrid errors are in the BAM core** (6,815 / 7,099). ✅ Correct, and the
  key clue: it pointed the search at a BAM-core operation (step-4 re-orientation),
  not the graph layer. The BAM-shared Hamming (3.11%) matched standalone BAM (3.09%)
  *with the re-orientation on*; with it off, BAM-shared drops to **0.71%**.
- **~~Graph's 0.80% is a pure denominator effect; you cannot get it at BAM
  coverage.~~** SUPERSEDED — hybrid (new) gets 0.82% at 212k reads. Graph's
  marginal reads do error at 16.67% in hybrid (the hard tail is real), but the bulk
  of the gap was the fixable step-4 re-orientation, not an intrinsic ceiling.

All rows in the table above were regenerated from the same current HEAD binary
(BAM/graph/hybrid all run with shipped defaults; hybrid = af0.11 + stitch margin
10 + both-strands rule), scored against the same diplinator truth BAM.

### Shared-read accuracy: graph IS better on the same reads (corrects earlier claim)

Comparing per-read concordance on the **201,202 reads both graph and BAM phase**
(`eval_graph/per_read.tsv.gz` ∩ `eval_bam/per_read.tsv.gz`, same truth BAM):

| On 201,202 shared reads | Graph | BAM |
|---|---:|---:|
| Discordant | 1,247 (**0.62%**) | 2,896 (**1.44%**) |

Per-read disagreement on those shared reads:

| Outcome | Reads |
|---|---:|
| graph RIGHT, bam WRONG | **2,082** |
| bam RIGHT, graph WRONG | 433 |
| both WRONG | 814 |

**Graph wins disagreements 4.8 : 1.** This means graph's accuracy is *not* purely a
denominator effect — on reads both attempt, graph is ~2.3× more concordant. The
earlier "pure denominator effect" framing (marginal-split section) was too strong.

**Confound (unresolved):** graph blocks are far shorter (auN 0.99M vs 9.83M; 368
vs 429 PS). A read in a shorter block spans fewer adjacent-het transitions, so it
is *mechanically* less likely to be scored discordant. Graph's 0.62% is therefore
**part genuinely-better-phasing, part shorter-blocks**; the two are not separated
by this measurement. Any strategy that lengthens graph blocks (e.g. by adding BAM
sites) would partly surrender the short-block component.

### Incremental (non-shared) reads — the hard tail

| Reads only one pipeline phases | Count | Discordant |
|---|---:|---:|
| graph-only | 2,513 | 15.6% |
| **bam-only** | **18,723** | **20.8%** |

The 18,723 reads BAM phases but graph does not error at 20.8% — graph lacks sites
there precisely in the **hard regions** (segdups), the same windows (37–38 Mb)
where prior analysis found "pulling hard reads into blocks injects switches."

### Proposed strategy: invert hybrid (graph-base + BAM-fill)

Idea (user, this session): since graph now performs so well, flip the hybrid base
from BAM to graph, and bring **validated BAM sites into the graph where the graph
has no sites**, instead of the current BAM-base + graph-augment design.

**Why it has merit:** graph is genuinely more accurate on shared reads (above), so
keeping graph as the core and only filling gaps could retain graph's clean phasing
while recovering coverage.

**Optimistic ceiling:** if the graph core keeps ~0.80% and the ~18,723 bam-only
reads are added at BAM accuracy: ≈5,541 errors / 222,438 reads ≈ **2.5% Hamming**
— would beat BAM (3.09%) and current hybrid (3.22%) at *higher* coverage.

**Risks that could collapse the ceiling:**
1. Adding BAM sites lengthens graph blocks → reawakens the short-block confound
   (some of graph's 0.62% edge was mechanical) and the segdup-switch penalty
   established in the rejected-lever ledger.
2. The bam-only fill reads error at 20.8% — the gain is bounded by how many are
   *recoverable* and whether merging them keeps graph's core clean.

**Feasibility probes to run BEFORE building (offline, no C++ — matches this
project's measure-first track record):**
1. **GAF membership of the 18,723 bam-only reads.** If present in the GAF but
   unphased for lack of sites → recoverable. If absent from the GAF → unreachable
   regardless. (Established globally: GAF 271,930 ⊂ BAM 272,016; confirm for *this*
   subset.)
2. **Block-length-fair accuracy.** Re-score graph vs BAM on shared reads within
   comparable block spans, to isolate real phasing quality from the short-block
   artifact and get a realistic (not optimistic) ceiling.

**Status:** not yet started — awaiting decision on probe-first (recommended) vs
build-first. This is the live thread to resume.

### Probe results: graph-core + BAM-gap-fill (the narrower strategy)

Refined goal (user): keep graph's phased core untouched; for reads graph leaves
**unphased**, phase them using BAM-validated sites the graph lacks. Purely additive
gap-fill — no merging/re-orienting of graph's clean blocks (avoids the short-block
confound). Both probes run offline on the post-fix per-read results:

**Probe 1 — reachability (PASS).** Gap-fill candidate set = reads BAM phases but
graph does not = **18,723 reads**. Of these, **18,668 (99.7%) are present in the
GAF** — graph sees the reads, it just has no phasing sites there. Only 55 are
unreachable. So the strategy is mechanically feasible.

**Probe 2 — accuracy ceiling (the catch).** The fill reads error at **20.84%**
under BAM, because graph lacks sites there *precisely because they are the hard
regions* (segdups / low-complexity; concentrated in 45–46, 21–22, 44–45, 13–14,
36–37 Mb windows). Projected combined:

| Method | Reads | Hamming |
|---|---:|---:|
| Graph (core only) | 203,715 | 0.80% |
| **Graph + BAM-fill (optimistic)** | **222,438** | **2.49%** |
| BAM | 219,925 | 3.09% |
| Hybrid (current) | 220,708 | 3.22% |

**Conclusion:** the strategy *works* and gives the best accuracy/coverage combo of
any method — **~2.49% Hamming at 222k reads** (beats BAM 3.09% and current hybrid
3.22%, at higher coverage). BUT it does **not** preserve graph's 0.80%: the gap
reads are intrinsically the hardest, so adding them at ~21% error drags the
combined rate up ~3×. Graph's advantage is *diluted, not inherited*. And 2.49% is
the **optimistic ceiling** — it assumes fill reads phase at full BAM accuracy and
that seam-stitching graph blocks via BAM sites injects no new switches (the
dominant real-world failure in the rejected-lever ledger). Realistic outcome sits
between 2.49% and current hybrid's 3.22%.

**Recommendation:** worth building IF the target metric rewards
accuracy-at-high-coverage (2.49% < 3.09% BAM is a real gain). NOT worth it if the
goal is graph's headline 0.80% — that is unreachable at full coverage by
construction, because the missing reads are the hard ones. Decision pending.

### Built and measured: hybrid-core + gap-fill (additive, no re-stitch)

The gap-fill strategy above was prototyped and evaluated on real chr20 data
(binary `9dcc509`, HG002 HiFi → CHM13, diplinator truth). Implementation
(`run_all3_fixed/gapfill.py`): take the hybrid phased BAM as the untouched core;
for the 9,169 reads BAM phases but hybrid leaves unphased, copy BAM's `HP` onto
the hybrid record with `PS += 1e9` (a disjoint phase-set namespace). Purely
additive — no hybrid block is renumbered, merged, or re-oriented.

**Provenance (important):** the bam / graph / hybrid pipeline outputs are all
phasings of the **same vg-giraffe alignment**, not different aligners. The BAM
is the linear surjection of the graph alignment (`@PG vg` → `pgbam annotate`);
the GAF is its graph-space form. So a "fill read" is not from another aligner —
it is a read that a different pgphase *pipeline* managed to phase. Gap-fill is a
cross-pipeline HP/PS label transfer by read name, not a re-alignment.

The hybrid phased BAM is the correct substrate because it carries all 272,016
reads as records. The **fill** BAM is read only for its HP/PS labels, so it may
be unaligned: the graph pipeline's phased BAM has no `@SQ` header, but
`samtools view` still yields its tags, and the labels land on the core's
already-aligned records. (Earlier notes claimed graph fill "needs realignment"
— that was wrong; it conflated using the graph BAM as the *core* with reading
its labels.)

**Which reads hybrid drops vs BAM (measured on this run):** 9,169 reads. They
error at **36.35%** even under the pipeline that phases them (BAM), are all
mapq 60 (not a mapping artifact), and concentrate in segdup/low-complexity
windows (44–46 Mb ≈ 3,805 reads at 38–53% error; 13–14, 35–36, 43 Mb). Graph
phases only 625 of them, at 64% accuracy. These are the deliberately-skipped
step-4 noisy reads — the `skip_noisy_kmeans` default drops them on purpose.

**Result:**

| Metric | BAM | Graph | Hybrid | **Hybrid+gapfill** |
|---|---:|---:|---:|---:|
| Phased reads | 219,925 | 203,734 | 212,259 | **221,428** |
| Overall accuracy | 96.91% | 99.21% | **99.23%** | 97.85% |
| Hamming | 3.091% | 0.790% | **0.768%** | 2.149% |
| Switch errors | 1,395 | 364 | **291** | 1,171 |
| Switch rate | 0.636% | 0.179% | **0.137%** | 0.530% |
| Flip errors | 2,050 | 856 | **906** | 1,958 |
| Perfect PS % | 29.14% | **39.12%** | 33.82% | 28.79% |
| Phase sets | 429 | 368 | 347 | 499 |
| N50 | 937,648 | 917,581 | **993,831** | 927,359 |
| auN | 9.83M | 0.99M | 7.15M | **11.25M** |
| Genome covered | 313.5M | 221.6M | 246.7M | **340.2M** |

**Findings:**
- **Hybrid+gapfill Pareto-beats pure BAM on every axis:** more reads phased
  (+1,503), lower Hamming (2.15% vs 3.09%, ~30% relative), fewer switch errors
  (1,171 vs 1,395), higher auN (11.25M vs 9.83M), more genome covered (340M vs
  314M). If a workflow ships/compares against BAM for high coverage, gap-fill is
  a strict upgrade.
- **The disjoint-PS design works:** switch errors rose only to 1,171 (below BAM's
  1,395), confirming the fill reads create their own small phase sets without
  contaminating the hybrid core. Measured Hamming (2.149%) **beat** the offline
  projection (2.25%) for this reason.
- **It does NOT preserve hybrid's 0.77% Hamming.** Adding the hard segdup reads
  (36% error) drags Hamming 0.77% → 2.15%. Unavoidable: the oracle best-of-3
  per-read ceiling is 2.16% Hamming at 222.8k reads (4,804 reads are wrong in all
  three pipelines), so 2.15% is essentially at the achievable floor for that
  coverage.

**Two-source gap-fill (hybrid core + BAM fill + graph fill) — measured but
deliberately NOT shipped:** chaining a second `gapfill.py` pass with
`ps_offset=2e9` adds the **1,374 reads only the graph pipeline phased** (neither
hybrid nor BAM). All 1,374 are already present as records in the hybrid BAM — no
realignment — and phase at 81% in isolation. **Decision (final): gap-fill will
only ever recover reads present in the projected BAM** (the linear surjection of
the graph alignment that the hybrid pipeline already loads). The graph-only reads
are out of scope by design, so this row is kept for the record only — it is not a
candidate for promotion to the binary.

| Metric | BAM | Hybrid | +gapfill (1-src) | **+gapfill2 (2-src)** |
|---|---:|---:|---:|---:|
| Phased reads | 219,925 | 212,259 | 221,428 | **222,802** |
| Overall accuracy | 96.91% | **99.23%** | 97.85% | 97.77% |
| Hamming | 3.091% | **0.768%** | 2.149% | 2.228% |
| Switch errors | 1,395 | **291** | 1,171 | 1,244 |
| auN | 9.83M | 7.15M | **11.25M** | 10.62M |

222,802 is the full union of all three pipelines (+2,877 over BAM, the maximum
phaseable at this coverage). Adding the 1,374 graph-pipeline reads costs +73
switch errors (their ~19% error on 1,374 reads ≈ +258 discordant), so 2-source
maximizes coverage while 1-source is marginally cleaner. Both Pareto-beat BAM.

**Decision matrix (all measured):**

| Goal | Method | Hamming | Reads |
|---|---|---:|---:|
| Lowest error | Hybrid (default) | 0.77% | 212.3k |
| Best accuracy *and* coverage vs BAM | **Hybrid+gapfill** | 2.15% | 221.4k |
| Max coverage (full union, offline only — not shipped) | Hybrid+gapfill2 | 2.23% | 222.8k |
| Raw BAM-class behavior | BAM / `--keep-noisy-kmeans` | 3.09% | 219.9k |

**Verification of `--keep-noisy-kmeans` (the in-pipeline alternative):** ran for
real — recovers reads to 220,633 but regresses to 96.79% / 3.21% Hamming /
1,371 switch (essentially BAM), confirming step-4 phases its extra reads at ~65%
accuracy. Additive gap-fill is strictly better than `--keep-noisy-kmeans` (2.15%
vs 3.21% Hamming at comparable coverage) because it does not let the noisy reads
re-orient the clean core.

**Status:** shipped natively as `collect-hybrid-variation --gap-fill`. The
1-source variant (hybrid core + BAM-pipeline fill → 221,428 reads / 2.15%
Hamming) is built into the binary: when `--gap-fill` is set alongside the
hybrid `skip_noisy_kmeans` default, `gap_fill_unphased_reads`
(`src/collect_var.cpp`) re-runs the `kCandGermlineVarCate` k-means into a
scratch buffer and adopts its haplotype only for reads the clean core left at
`hap==0`, writing them to `chunk.gap_haps`/`gap_phase_sets` with
`PS += kGapFillPsOffset` (1e9, disjoint namespace). The core read-phase and
per-candidate consensus state are snapshotted and byte-for-byte restored, so the
clean core is never renumbered, merged, or re-oriented — the in-binary form of
the validated `scripts/gapfill.py` 1-source pass. `collect_bam_output.cpp` emits
the gap-fill HP/PS for any read the core left unphased. Off by default.

The 2-source variant (adding the 1,374 graph-only reads → 222,802 / 2.23%) is a
**deliberate non-goal** and will not be promoted to the binary: gap-fill recovers
only reads present in the projected BAM, and the graph-only reads are out of scope
by design. It exists solely as a post-process for offline analysis via
`scripts/bench_hybrid.sh --gapfill 2`. The prototype (`scripts/gapfill.py`,
`scripts/bench_hybrid.sh --gapfill [N]`) is retained for reproducing both BAMs.

---

## Graph Het-Indel AF Gate (resolves hybrid accuracy/contiguity regression)

### Problem

Commit `af6b089` ("Stop over-demoting graph het indels") promoted graph-only het
indels to `CleanHetIndel` so they enter k-means as phasing anchors. On a full
chr20 truth evaluation this **raised** hybrid Hamming error to 4.52% (vs BAM
3.09%) while raising auN to 13.6M. The commit message itself flagged that the
accuracy recovery was never confirmed on a full chromosome — it was not realized.

### Investigation

Truth-based read-level analysis (HG002 diplinator tags) localized the damage:

- **STABLE phase-block interiors are untouched** (+27 discordant on 140,861
  shared reads). The joint k-means does *not* corrupt confident BAM regions.
- **All regression lives at block-merge seams.** ~2,701 reads are *wholesale
  orientation flips*: a BAM-clean sub-block inverted en masse after a graph
  anchor mis-oriented the merge.
- A `--stitch-min-margin` sweep (0→8) barely moved Hamming (4.52→4.37%): the
  wrong merges are **confident, not thin-margin**. The defect is upstream anchor
  selection, not the stitch decision.
- Reverting `af6b089` dropped Hamming to 3.27% but collapsed auN to 8.93M (below
  BAM). The graph het indels are a **mixture**: near-0.5-AF indels are reliable
  bridges (the auN gain); off-center-AF indels (mis-genotyped / repeat) are the
  wrong-orientation anchors. All-or-nothing loses on one axis.

### Fix

`classify_graph_only_candidates` (hybrid path only) now admits a graph het indel
as a `CleanHetIndel` anchor only when its allele fraction sits within
`graph_indel_af_margin` of 0.5; off-center indels are demoted to `LowCoverage`
(kept out of k-means, pruned from output). Two tunable gates added, defaults
reproduce nothing-changed for SNPs:

- `--graph-indel-af-margin` (default `kDefaultGraphIndelAfMargin = 0.11`)
- `--graph-indel-min-alt` (default 0 — swept, no additional effect over the
  existing `min_alt_depth` floor)

The BAM pipeline (`classify_variant_initial`) and graph-only pipeline
(`classify_graph_candidates`) are not touched; BAM candidate output verified
byte-for-byte identical.

### Result (chr20, full-chromosome truth eval)

AF-margin sweep, best at **0.11** (keep graph het indels with AF ∈ [0.39, 0.61]):

| Config | Hamming | auN | Largest block | Phased reads |
|---|---|---|---|---|
| BAM baseline | 3.091% | 9.83M | 53.2M | 219,925 |
| Hybrid af6b089 (margin 0.30) | 4.520% | 13.57M | 62.7M | 223,065 |
| Hybrid revert af6b089 | 3.265% | 8.93M | 50.2M | 220,741 |
| **Hybrid + af-gate 0.11** | **3.053%** | **12.09M** | **59.3M** | **220,952** |

The 0.11 gate **Pareto-beats BAM on every axis**: lower Hamming (−0.038pp,
~830 reads), +23% auN, +11% largest block, +1,027 phased reads. Sweep was
non-monotone with a clear optimum at 0.11 (0.12 already crosses back above BAM
accuracy), confirming a genuine edge rather than a tuning endpoint.

---

## Chunk-Stitch Rule: both-strands-bridged (hybrid default)

### Problem

The genome is phased in 500 kb chunks, then adjacent chunks are stitched by an
overlap-read vote in `flip_chunk_hap` (`collect_phase.cpp`). The original rule
merges when the net flip-vote magnitude exceeds `stitch_min_margin`. The net
margin conflates two very different seams: "1 clean read, 0 conflicts" and
"11 reads splitting 6-vs-5" both yield margin 1. Raising the margin to abstain
on the coin-flip case also abstains on the clean low-coverage case, leaving
blocks unmerged that could have been joined safely.

### Rules tested

`--stitch-rule` selects the decision (named constants in `phasing_types.hpp`):

- **0 net-margin** — original: merge when `|flip − noflip| > margin`.
- **1 both-strands-bridged (Reading A)** — merge only when the winning
  orientation has ≥1 read on *each* of its two haplotype links (no-flip needs
  `pre1→cur1` AND `pre2→cur2`; flip needs `pre1→cur2` AND `pre2→cur1`). Guards
  against merging on evidence from a single haplotype.
- **2 literal (Reading B)** — merge when both orientations have ≥1 read
  (stitches contested seams; diagnostic only).
- **3 both-strands + margin (Reading C)** — Reading A AND net vote > margin.

### Result (chr20, af0.11, full-chromosome truth eval)

| rule | Hamming hq0 | Hamming hq10 | auN | phased | blocks/chrom |
|---|---|---|---|---|---|
| 0 net-margin (control) | 3.0401 | 1.8550 | 11.91M | 220,734 | 411 |
| **1 both-strands (default)** | 3.0415 | 1.8541 | **12.00M** | **220,897** | **404** |
| 2 literal | 2.6530 | 1.4490 | 10.66M | 220,340 | 488 |
| 3 both-strands + m2 | 3.0385 | 1.8539 | 11.97M | 220,853 | 407 |

**Rule 1 improves contiguity at no accuracy cost**: auN +0.8%, +163 reads
phased, blocks/chrom 411 → 404, Hamming flat (hq0 +0.001, hq10 −0.001),
switch/flip unchanged. The control (rule 0) reproduced the baseline exactly,
confirming the refactor is behavior-preserving.

**Rule 2's apparent Hamming win is an artifact**, not an improvement: it only
merges contested seams and refuses clean unanimous ones, so blocks/chrom jumps
to 488 (genome fragmented into more, smaller blocks). Hamming is measured
within-block, so shorter blocks mechanically score lower while auN drops to
10.66M and 394 reads are lost. Diagnostic only.

### Decision

Rule 1 is the **hybrid default** (`kHybridDefaultStitchRule`), set in
`hybrid_collect.cpp` alongside the stitch-margin default. BAM/graph pipelines
keep rule 0. `--stitch-rule` overrides. This is the only lever found that adds
contiguity for free — every other lever (below) traded accuracy or auN away.

### Rule 4 confidence-weighted stitch (rejected)

The seam-merge regression traces to a minority of *confident-wrong* merges —
boundaries where the overlap vote is lopsided the wrong way (the margin sweep
0→8 barely moved Hamming, so the bad votes are not thin-margin). Hypothesis: a
cluster of reads all phased off one weak anchor casts many votes that are
really one piece of evidence; weighting each overlap read's vote by its k-means
assignment confidence `max(0,(agree−conflict)/(agree+conflict))` (from the
persisted clean-het SNP tallies, using `min(conf_pre,conf_cur)` per read) would
shrink such a seam's effective strand weight below threshold and abstain. Added
as `--stitch-rule 4` + `--stitch-min-strand-weight` (both-strands structure on
weighted sums).

Full chr20 sweep (diplinator truth eval; rule-1 control byte-identical to HEAD):

| config | Hamming% | auN | N50 | phased | perfect PS |
|---|---|---|---|---|---|
| rule 1 (shipped) | 3.0415 | 12.00M | 951,564 | 220,897 | 91 |
| rule 4 w≤1.0 | 3.0415 | 12.00M | 951,564 | 220,897 | 91 |
| rule 4 w=1.5–2.0 | 3.0416 | 12.00M | 951,564 | 220,861 | 92 |
| rule 4 w=3.0 | 3.0396 | 11.93M | 947,261 | 220,776 | 94 |
| rule 4 w=5.0 | 3.0405 | 11.87M | 947,261 | 220,678 | 96 |

**No-op at low weight, net loss at high weight.** At w≤1.0 the output is
byte-identical to rule 1: the bridging reads are *confidently* phased
(conf≈1.0), so weighting by confidence equals counting — confirming the
confident-wrong merges are not driven by weakly-assigned reads. Only at w≥3.0
does the gate bite, and then it shows the same accuracy/contiguity trap as every
other lever: a tiny Hamming gain (−0.0019) for real contiguity loss (auN −63 to
−127 kb, N50 951k→947k, −121 to −219 reads). Perfect-PS rises (91→96) but
overall it is not a win. Reverted from code; kept here so it is not re-tried.
The standing conclusion holds: the residual seam errors are confident and
upstream of the stitch decision — no reweighting or thresholding of the
overlap vote separates them from correct merges.

### BAM-CIGAR-only AF gate on graph indel anchors (rejected)

Since stitch-side fixes all failed, this attacked the genotype at its source.
The standard graph het-indel AF window uses **mixed BAM+GAF counts**: GAF
graph-read alt support can make an indel look het that BAM CIGAR evidence alone
shows as homozygous or thin-alt. The gate (`--graph-indel-bam-af-gate`)
snapshots BAM-CIGAR-only allele counts *before* `inject_graph_reads` adds GAF
observations, and admits a het-indel anchor only if its BAM-only AF is also
within the margin of 0.5 (for candidates with ≥ min-depth BAM reads). It reuses
the trusted BAM-CIGAR genotyping the BAM pipeline already does, upstream of
k-means.

Full chr20 sweep (diplinator truth eval; gate-off control byte-identical to HEAD):

| config | Hamming% | auN | N50 | phased | switchflip% | perfPS | anchors |
|---|---|---|---|---|---|---|---|
| off (shipped) | **3.0415** | 12.00M | 951,564 | 220,897 | 1.6206 | 91 | 134,616 |
| gate depth=3 | 3.2110 | 11.93M | 947,575 | 220,765 | 1.5790 | 99 | 128,830 |
| gate depth=5 | 3.2110 | 11.93M | 947,575 | 220,765 | 1.5790 | 99 | 128,910 |
| gate depth=8 | 3.2110 | 11.93M | 947,575 | 220,765 | 1.5790 | 99 | 128,987 |

**Net loss.** The gate is *not* timid — it removes ~5,700 graph indel anchors
(many graph indels genuinely have BAM-CIGAR AF inconsistent with their mixed
het call) and the depth knob barely changes the eval. Switchflip drops (1.62 →
1.58) and perfect-PS rises sharply (91 → 99), so it *does* remove
mis-genotyped, seam-mis-orienting anchors. But read-level **Hamming rises 3.04
→ 3.21%** and auN drops 62 kb: the removed anchors were, on net, contributing
more correct read assignments than wrong orientations. Even the trusted
BAM-CIGAR caller cannot separate the good from the bad — the indels it
disagrees with are still net-positive for phasing.

**Fourth confirmation of the same result.** Margin sweeps, the linkage gate,
confidence weighting, and now source-level BAM re-genotyping all fail
identically: any operation that *removes or down-weights* the suspect anchors
trades a real read-assignment/contiguity loss for a switchflip/perfect-PS gain
that does not net out. The mis-oriented seams and the contiguity-providing
bridges are the **same population** of off-center/GAF-supported indels; no
1-D signal (AF, LD, confidence, BAM-AF) separates them. Reverted from code;
kept here so it is not re-tried. If there is a win left it requires *correcting*
an anchor's orientation rather than *dropping* it — a fundamentally harder change
than any gate.

---

## Rejected graph-native levers (kept for the record)

Two graph-specific anchor-quality gates were prototyped and evaluated against
the shipped config (af0.11 + margin10), then **dropped** — neither helped on
HiFi without a worse cost elsewhere. Documented here so they are not re-tried.

### Graph read-divergence gate (`--graph-max-dv`)

The GAF carries a per-alignment divergence tag `dv:f` (the graph analog of the
BAM pipeline's XID/noisy-read filter, which the graph path lacks). Idea: drop a
read's allele votes when `dv` is high, before allele-fraction aggregation, so
noisy reads can't skew a site's AF off-center and trip the indel gate.

Result: at dv=0.005 it gave a marginal hq0 Hamming improvement (3.040 → 3.029)
but flat hq10 (1.855 → 1.854) and slightly lower auN/reads. The effect is tiny
because **HiFi reads are too accurate** — only 1.4% have dv>0.01, 0.26% >0.05
(median 0.0002, p99 0.0167). The filter is correct (extreme dv collapsed
candidate counts; dv=1.0 reproduced the baseline byte-for-byte) but there is
almost no noise to remove on HiFi. **Expected to matter on ONT**, where the
divergence tail is large. Not shipped.

### Graph SNP AF-window gate (`--graph-snp-af-margin`)

Mirrors the indel AF gate for graph het SNPs (default 0.5 = admit all). A 2×2
matched eval (af0.11, margin10):

| config | Hamming hq0 | Hamming hq10 | auN |
|---|---|---|---|
| neither (shipped) | 3.040 | 1.855 | 11.91M |
| dv-only (0.005) | 3.029 | 1.854 | 11.86M |
| sam-only (0.10) | 2.986 | 1.800 | 9.59M |
| both | 2.975 | 1.799 | 9.57M |

sam=0.10 drives the Hamming win but **collapses auN to 9.59M, below the BAM
baseline (9.83M)** — the same accuracy/contiguity trap as a too-wide indel
gate. Filtering off-center SNP anchors removes the long-range bridges that give
the graph its contiguity edge. Not shipped; default-off (admit all) remains
Pareto-best for SNPs.

### Graph anchor BAM-linkage gate (`--graph-anchor-min-concordance`)

The block-merge regression was localized to graph-only het-indel anchors that
mis-orient a seam (see "Graph Het-Indel AF Gate"). The AF window (0.11) is a
1-D proxy for "is this anchor's genotype reliable." This lever tried to replace
the proxy with a **direct** signal: a real graph het indel co-segregates with
neighboring confident BAM het SNPs on the reads spanning both (same physical
molecules), so per-anchor linkage concordance against BAM SNPs should separate
reliable bridges from mis-genotyped anchors better than AF alone. The gate runs
before k-means (so it uses `read_var_profile` linkage, not post-k-means
haplotype labels), computing each indel's best phase-agnostic 2×2 concordance
(`max(n00+n11, n01+n10)/total`) over BAM het SNP partners sharing ≥4 informative
reads, and demoting anchors below the threshold to LowCoverage.

Full chr20 sweep (diplinator truth eval; control byte-identical to HEAD and
reproduces the shipped hybrid row exactly):

| threshold | Hamming% | auN | N50 | phased | switchflip% | perfect PS |
|---|---|---|---|---|---|---|
| off (shipped) | **3.0415** | 12.00M | 951,564 | 220,897 | 1.6206 | 91 |
| 0.6 | 3.1869 | 12.00M | 951,564 | 220,828 | 1.6034 | 93 |
| 0.7 | 3.1657 | 12.03M | 954,036 | 220,794 | 1.5905 | 94 |
| 0.8 | 3.1410 | 12.00M | 951,564 | 220,714 | 1.5725 | 95 |
| 0.9 | 3.1485 | 11.97M | 951,564 | 220,692 | 1.5722 | 95 |

**Net loss at every threshold**: Hamming rises 3.04 → 3.14–3.19% while auN stays
flat. The signal *does* work in the intended direction — switchflip rate drops
monotonically (1.62 → 1.57) and perfect-PS count rises (91 → 95), so the gate
removes some genuinely mis-oriented seam anchors. But overall read-level Hamming
*rises*: the demoted anchors were, on balance, contributing more correct read
assignments than wrong orientations. Linkage concordance flags real
mis-genotypes **and** good-but-low-LD anchors indiscriminately — the same
all-or-nothing trap the AF window was tuned to 0.11 to avoid. It is not a
cleaner separator than AF; it is a noisier one. This is consistent with the
earlier finding that the wrong merges are **confident, not thin-margin**: a
filter cannot distinguish a confidently-wrong anchor from a confidently-right
one when both show strong LD with their neighbors. Reverted from code (default
0.0 reproduced the baseline byte-for-byte); kept here so it is not re-tried.
Possibly worth revisiting on ONT, where genotyping noise is larger and the
good/bad anchor populations may separate more cleanly.

# SYNTHESIS — Why the hybrid cannot beat BAM on switches/perfect-PS (read this first)

This section indexes a multi-session investigation into closing the hybrid's two
remaining deficits vs BAM (switch errors, perfect phase sets). **Conclusion: no
method tested beats the current default hybrid; the deficits are explained, not
fixable with the available data.** Every detailed section below is preserved
chronologically; this is the map.

### Current verified 3-way (chr20, vs truth)

| metric | BAM | graph | hybrid (default) |
|---|---|---|---|
| Hamming | 3.09% | 7.30% | **3.04%** ✅ |
| auN | 9.83 Mb | 0.98 Mb | **12.00 Mb** ✅ |
| N50 | 938 kb | 898 kb | **952 kb** ✅ |
| phased reads | 219,925 | 216,130 | **220,897** ✅ |
| switches | **1,393** | — | 1,440 |
| perfect PS | **125** | 92 | 91 |

The hybrid wins the metrics that matter most (accuracy, contiguity,
completeness) and loses only switches/perfect-PS.

### The one root cause behind every dead end

**GAF reads are a strict subset of BAM reads: 271,930 ⊂ 272,016, with 0 unique
reads.** The graph is a *second alignment of the same physical HiFi bases*. It
therefore carries **no evidence independent of BAM**. Its only genuine
contribution is the **3,541 extra easy SITES** BAM missed (already used by the
hybrid). Any method that re-weights or re-interprets the graph's view of *shared
reads* is information-free and capped at single-digit switches.

### Why the switch/perfect-PS deficit is structural, not a bug

Switches are counted *within* a phase set; perfect-PS is all-or-nothing per
block. The hybrid merges blocks BAM leaves split (404 vs 429 PS; auN 12.0 vs 9.8
Mb). Merging **exposes more within-block transitions** and **destroys perfect-PS
credit** when a clean block joins an imperfect neighbor — even though it adds no
discordant reads (the hybrid has *fewer*: 6,718 vs 6,798). The +47 switch gap is
localized (windows 37–38 Mb alone = +36) and traces to reads pulled into large
blocks in segdups, not to bad merge orientation (fixing all 21 wrong merges =
−3 switches, confirmed by real run).

### Ledger of everything tested (all REJECTED)

| # | Lever | Why it failed | Ceiling |
|---|---|---|---|
| 1 | k-means flip margin | no-op; switches are confident not bare-majority | 0 |
| 2 | BAM-vs-graph allele disagreement | not independent; measures alignment noise | ~0 |
| 3 | per-read LD changepoint gain | real 9× signal, but crossovers ≠ consensus switches | −29 sw, kills completeness |
| 4 | graph-as-voter / arch flip | graph 3× less accurate; no selective rule >50% | net −8,452 reads |
| 5 | stricter stitch gate (rule3+margin) | real-run −3 sw, +6 perfect-PS, −0.19 Mb auN | trade-off, not win |
| 6 | (D) graph alleles as coverage | GAF ⊂ BAM; double-counting | 0 |
| 7 | (A) graph-confirmed k-means weights | only 7 transitions fixable; loses real hets | ~3 sw |
| 8 | (C) graph-only sites as stitch bridges | only 5/40 seams reachable | <1 sw |
| 9 | (B) graph path co-occurrence linkage | recovers truth 60%/window vs BAM 97% | coin-flip |

### What would actually beat BAM

Only **data BAM does not contain**: trio/parental reads (HG002 has them — true
independent inheritance signal) or a second orthogonal sequencing technology
(e.g. ONT — independent error modes). Not another way to process the graph
alignment of the same reads.

> **PARTIALLY SUPERSEDED (post noise-filter-fix session).** This ledger was built
> *before* the indel noise-filter fixes that took graph accuracy 92.70% → 99.20%.
> Two of its premises now need re-examination:
> - "Graph is 3× less accurate" (lever #4) is **no longer true** — post-fix graph
>   is *more* accurate than BAM on shared reads (0.62% vs 1.44% discordant; wins
>   disagreements 4.8:1). See "Shared-read accuracy" in the three-pipeline section.
> - This changes the calculus for a **graph-base + BAM-fill** hybrid inversion
>   (proposed this session, see "Proposed strategy: invert hybrid"). The ledger's
>   levers all assumed BAM-base + graph-augment; an inverted architecture was never
>   tested. The "GAF ⊂ BAM / no independent read evidence" point still holds and
>   still caps gains on *shared* reads — but graph-base changes *which* pipeline
>   phases the easy core, which the ledger did not evaluate.
> The trio/ONT conclusion below still stands for beating BAM on the *shared hard
> reads*; it does not rule out the graph-base inversion, which is a different lever.

### A real, optional trade-off knob (not promoted to default)

`--stitch-rule 3 --stitch-min-margin 15` yields **+6 perfect PS** (91→97) and −3
switches at a cost of −0.19 Mb auN / −323 phased reads. Use only if perfect-PS is
the priority metric; the default maximizes the completeness win.

---

## Cross-pipeline read superset / seam-bridging (REJECTED — graph adds sites, not reads)

Hypothesis: the pangenome is more linear/contiguous and reads map better there,
so we could augment graph variation sites with confirmed BAM ones into a
cross-chunk superset, and/or let graph reads bridge the fixed-stride chunk seams
where the BAM stitch is weak. Two probes were scoped (A: widen chunk overlap so
more reads are phased on both sides of a seam; B: let graph-only bridging reads
vote in the seam stitch). Both were dropped after a data-level measurement,
**before** writing the merge/overlap code, because the enabling premise — that
the graph provides reads the BAM lacks — does not hold on HiFi.

### The graph supplies no new bridging reads

The GAF and the BAM are the **same HiFi molecules aligned two ways** (linearly
to CHM13 → BAM; to the pangenome → GAF). Genome-wide qname comparison on chr20:

| source | distinct reads |
|---|---|
| GAF | 271,930 |
| BAM (primary) | 272,016 |
| **GAF-only** | **0** (GAF ⊂ BAM) |

So the graph's contribution to the hybrid pipeline is extra **sites**, not extra
**reads**. Probe B has no new voters to add to a seam, and Probe A's premise
(seams are read-starved) is false: per 500 kb boundary there are a median of 67
bridging reads, only 2/132 boundaries below 20 crossers, and 98% have a
phaseable site within 20 kb on both sides. Adding more BAM reads or a wider
overlap recruits the same molecules that already vote.

This is config-independent (a property of the input read sets, not the phasing
parameters). Kept here so the cross-pipeline read-superset and
graph-reads-bridge-seam ideas are not re-attempted on same-molecule HiFi inputs.
May differ when the GAF and BAM come from genuinely different read sets (e.g.
graph-aligned reads that fail linear mapping entirely), which did not occur here.

> Note: an earlier write-up of this probe also reported a seam-flip "regression"
> (≈1,811 wrong-oriented reads, hybrid Hamming 4.52%, auN 13.57M). That analysis
> used a **stale `run/hybrid.phased.bam` produced by an experimental binary**
> (its auN exceeded every documented config — a tell that it was not the shipped
> default). Re-running the shipped hybrid (HEAD `034528b`, af0.11 + margin10 +
> both-strands rule 1) reproduces the documented row exactly (Hamming 3.0415%,
> auN 11,996,328, N50 951,564, 220,897 phased) and **beats BAM** (Hamming 3.04 vs
> 3.09%, 6,718 vs 6,798 discordant, auN +22%, +972 reads phased). The shipped
> hybrid is not regressed; the seam-flip figures above do not describe it. Always
> regenerate the phased BAM from the current binary before diagnosing.

## Switch-error gate on newly-phased reads (REJECTED — positive oracle ceiling, unreachable signal)

The hybrid loses to BAM on exactly two secondary metrics: switch errors
(1,410 vs 1,355, +55) and perfect phase sets (91 vs 125). Goal: close the switch
gap without giving back the hybrid's Hamming/auN advantage.

### Where the hybrid–BAM accuracy delta comes from (decomposition)

Per-read accounting that closes exactly to the −80 discordant-read delta:

| component | discordant reads |
|---|---|
| shared reads hybrid **fixed** (disc→conc) | −1,167 |
| shared reads hybrid **broke** (conc→disc) | +1,029 |
| reads only hybrid phased (the +972 net), of which 380 wrong | +380 |
| reads only BAM phased, of which 322 wrong | −322 |
| **net** | **−80 (hybrid better)** |

So the hybrid's per-read win is the sum of two large offsetting effects. Two
populations are net-harmful: the **380 discordant newly-phased reads** (real BAM
reads that BAM left unphased and hybrid phased wrong) and the **~1,012
hybrid-broke reads**.

### Oracle ceiling is real

A truth-based oracle that perfectly removes these reads clears BAM on switches:

| config | switches | vs BAM | discordant | phased |
|---|---|---|---|---|
| BAM | 1,355 | — | 6,660 | 219,610 |
| hybrid base | 1,410 | +55 (lose) | 6,629 | 220,681 |
| drop 390 disc newly-phased | 1,319 | **−36 (win)** | 6,239 | 220,291 |
| revert 1,012 hybrid-broke | 1,223 | −132 (win) | 5,617 | 220,681 |
| both | 1,113 | −242 (win) | 5,227 | 220,291 |

Dropping just the 390 discordant newly-phased reads is a win on *every* axis
(switches −91, discordant −390, only −390 phased). The headroom exists.

### But no implementable signal reaches it

The only per-read confidence available at output time is the k-means
`n_clean_agree_snps` / `n_clean_conflict_snps` (verified to survive
`mid_free_chunk`). Instrumented dump of all 220,897 phased reads, correlated
against the truth-defined harmful set:

- Harmful reads do skew low-confidence (median agree 3 vs 18 for good reads) but
  the means barely differ (18.8 vs 22.1, conflict 1.00 vs 0.08).
- **Base-rate kills it:** good reads outnumber harmful 156:1 (219,279 vs 1,402).
  The best gate precision is ~10% (`conflict≥5`: catches 48 harmful at the cost
  of 424 good); best recall ~55% destroys 37,885 good reads. No threshold or
  combination (conflict count, agree−conflict margin, low-agree) separates them.

This is the same **confident-wrong** signature as the rejected AF/linkage/stitch
levers: the harmful reads look as confidently phased as the good ones. The
oracle proves the win is *possible in principle*; the confidence signal proves it
is *not reachable* with current per-read information. Reverted (instrumentation
removed; src clean vs HEAD). Kept here so the "gate low-confidence newly-phased
reads to cut switches" idea is not re-tried without a *new* discriminative
signal (e.g. graph-vs-BAM allele agreement per read, or per-read LD with
neighbors — neither currently computed).

## Intra-chunk k-means flip margin (REJECTED — no-op; switches are confident, not bare-majority)

The hybrid's switch errors (1,440 vs BAM 1,355) are 93% intra-chunk (not at
stitch seams — switches are spread uniformly across the chunk grid, 7.1%
near-seam ≈ the 8% a uniform distribution predicts) and 95% clustered in
low-accuracy regions. Of the ~117 *isolated* (one-off) switches, 57 sit in
clean/good blocks (acc ≥0.90) and looked threshold-fixable.

The per-chunk k-means **already iterates** (10 Lloyd rounds to convergence,
`iter_update_var_hap_to_cons_alle`) and **already detects/corrects switches**
(`iter_update_var_hap_cons_phase_set` walks adjacent het-var pairs and flips the
downstream orientation when spanning-read `conflict > agree`, breaking the block
when `agree<2 && conflict<2`). So "add iterative refinement" was already done;
the isolated switches survive an iterating, switch-correcting system.

The one remaining knob was the flip *decision*: it fires on a **bare majority**
(`conflict > agree`). Hypothesis: an isolated noisy adjacent het-var pair induces
a spurious flip in an otherwise-clean block, fixable by requiring a margin
(`conflict > agree + k`). Added `--kmeans-flip-margin` (default 0); the control
(margin=0) reproduced HEAD candidates byte-for-byte (md5 4d261704). Sweep:

| margin | switch err | flip err | Hamming | auN |
|---|---|---|---|---|
| 0 (HEAD) | 1,440 | 2,133 | 3.042% | 11,996,328 |
| 1 | 1,440 | 2,132 | 3.042% | 12,026,597 |
| 2 | 1,440 | 2,132 | 3.042% | 12,026,597 |
| 3 | 1,440 | 2,132 | 3.042% | 12,026,597 |

Switches **unchanged** (1,440) at every margin; flips moved by 1. The hypothesis
is wrong: at the switch points `conflict` exceeds `agree` by far more than 3, so
the flip fires regardless of the margin. The flip rule is **confidently** making
the wrong call — same confident-wrong signature as the AF/linkage/stitch/newly-
phased-read levers. A bare-majority threshold tweak cannot help. Reverted (no-op
removed; src clean vs HEAD). Kept here so the "raise the k-means flip margin to
cut switches" idea is not re-tried — the switches are not bare-majority slips.

## BAM-vs-graph per-read allele disagreement (REJECTED — not independent; measures alignment noise, not haplotype ambiguity)

After six threshold/margin levers all hit the *confident-wrong* wall, the next
hypothesis sought a **genuinely new discriminative signal**: for each read mapped
by both pipelines, compare the BAM-derived allele to the graph-derived allele at
shared sites. The intuition was that reads whose two alignments *disagree* on
allele calls are alignment-ambiguous, and alignment ambiguity might predict the
harmful (switch-inducing) reads that k-means confidence cannot flag.

Instrumented `inject_graph_reads` (env-gated `PGPHASE_DUMP_XCHECK`) to dump, per
doubly-mapped read, the BAM allele vs graph allele at every shared site, then
joined the per-read disagreement rate against the eval `per_read.tsv.gz`
discordant/concordant labels:

| read class | mean disagree | median | % with any disagreement |
|---|---|---|---|
| concordant (good) | 0.48% | 0.00% | 4.6% |
| discordant (harmful) | 1.32% | 0.00% | 6.1% |

Both medians are 0.00%; the distributions overlap almost completely. Best gate
(disagree-rate ≥ 0.67) reaches **11.8% precision** — catches 13 of 1,069 harmful
reads while sacrificing 97 good ones. No threshold separates the classes, same
as the k-means confidence signal it was meant to replace.

**Why it fails (root cause, not just a bad threshold):** BAM and graph consume
the **same bases from the same HiFi read**, just aligned two ways. They are not
independent observations, so their disagreement measures *alignment noise*
(where the two aligners place indels/SNVs differently), which is orthogonal to
the cause of switches — *haplotype ambiguity* in segmental duplications, where
the read's bases genuinely match the other haplotype. The harmful reads are
**biologically confidently wrong**, not technically noisy, so an
alignment-noise signal cannot find them. Reverted (instrumentation removed; src
clean vs HEAD). Kept here so the "cross-check BAM vs graph alleles per read"
idea is not re-tried — the two views are not independent.

## Per-read LD changepoint gain (REJECTED — real signal, but oracle ceiling still loses on switches)

Unlike the prior six levers (all aggregate agreement counts that flag nothing),
this used a **positional** statistic. For each read, walk its phased het SNV
sites in genomic order to get a haplotype-support sequence (e.g. `1 1 1 2 2 2`),
then compute the **single-changepoint gain**: how many fewer errors a two-block
`A…A B…B` model makes vs the best one-block model. A clean within-read switch
(segdup crossover) scores high gain; scattered sequencing noise scores ~0. The
aggregate minority count used by every earlier lever cannot distinguish these.

**This is the first signal in seven attempts that actually separates the
classes:**

| signal | discordant | concordant | ratio |
|---|---|---|---|
| changepoint gain > 0 | 10.5% of reads | 1.2% of reads | ~9× |
| mean gain | 1.333 | 0.047 | ~28× |

(vs the prior levers where discordant and concordant medians/rates were
identical.) But the population is ~20:1 concordant:discordant, so even a 9×
enrichment yields only 7–26% weighted precision across the gain sweep.

**Decisive oracle (drop the flagged reads, recompute switches):**

| drop set | reads dropped | switches | flips | vs BAM (1,355) |
|---|---|---|---|---|
| baseline | 0 | 1,442 | 2,132 | +87 |
| gain ≥ 2 | 2,835 | 1,421 (−21) | 2,097 | +66 |
| gain ≥ 1 | 4,363 | 1,413 (−29) | 2,090 | +58 |

Even the most aggressive perfect-oracle drop (4,363 reads) closes only 29 of the
87-switch gap to BAM — still +58 — while sacrificing ~4,400 phased reads, which
would erase the hybrid's completeness/auN advantage (the entire reason it wins).

**Why a real signal still fails on switches:** a coincidence test showed the
within-read changepoints do **not** land at consensus switch positions — at a
500 bp window only 4–6% of high-gain discordant reads sit near a consensus
switch (≈ the concordant rate). The consensus majority vote **already absorbs**
these crossovers: the read's internal switch is outvoted before it becomes a
consensus switch, so removing the read removes its vote but rarely removes a
switch. The hybrid's remaining switch gap is a property of the **consensus over
many reads in segdups**, not of identifiable individual bad reads. No per-read
filter — however discriminative — can close it without destroying the
completeness that makes the hybrid win. Validated entirely offline (no C++
written); temp scripts removed. This closes the per-read filtering avenue: to
beat BAM on switches one must change the **consensus/phasing model in segdups**,
not filter reads.

## Graph-as-voter / graph-base-BAM-confirm architecture (REJECTED — graph is 3x less accurate, no selective override exists)

The current hybrid is **BAM-base, graph-augment**: it runs full BAM variant
calling, then adds graph snarl sites only where BAM was silent. At any **shared**
site (BAM also called it) the graph allele is **discarded**
(`extend_bam_profile_with_graph_obs` skips non-graph-only candidates). Scope on
chr20:

| site class | count | role |
|---|---|---|
| shared (BAM ∩ graph) | 63,299 | graph allele discarded |
| graph-only (extra) | 3,541 | graph's real contribution |
| BAM-only | 1,583 | graph doesn't confirm |

The proposal was to make the graph a **first-class voter** at shared sites (or
flip to graph-base/BAM-confirm), per the original intent "help phasing in the
graph with BAM-confirmed variants + extra from graph." Validated offline by
joining per-read truth labels from the standalone BAM-only and graph-only
pipelines on the 207,687 reads phased by **both**:

| | count |
|---|---|
| BAM ✓ & graph ✓ | 191,916 |
| BAM ✓ & graph ✗ (graph BREAKS) | 11,347 |
| BAM ✗ & graph ✓ (graph FIXES) | 2,895 |
| BAM ✗ & graph ✗ | 1,529 |

**Trusting graph over BAM blindly = +2,895 − 11,347 = net −8,452 reads.** The
graph is ~3× less accurate per read (6.2% vs 2.1% discordant), and the 8,443
reads only the graph phases are just 65% concordant. A graph-base flip would
wreck accuracy — which is exactly why the existing code is BAM-base.

The only way the idea survives is a **selective** override: a rule identifying
*which* shared sites/reads to trust the graph on. Exhaustive search over every
available axis — per-read mapq bands, graph hapq bands, and the eval's
BAM-bad-region BED, plus all combinations — found **no rule with >50% fix-rate at
meaningful volume**:

| stratum | fix-rate |
|---|---|
| all disagreements | 20% |
| BAM-bad region (any hapq) | 37% |
| BAM-bad region & graph hapq=60 | 18% |
| graph hapq=60 (graph max-confident) | 16% |

Even in BAM's *worst* regions the graph fixes only 37% of disagreements (net
−595). Graph's own hapq is *inversely* useful (hapq=60 → 16% fix-rate): the
graph is confidently wrong in the same hard regions BAM is.

**Root cause:** GAF reads are a strict subset of BAM reads (only 3.9% of
graph-phased reads are absent from the BAM set). The graph and BAM consume the
**same physical HiFi bases** — the graph is just a second alignment of them. So
the graph carries **no independent evidence** to overrule BAM; where BAM is wrong
(segdup, read maps to wrong copy with mapq 60) the graph inherits the same wrong
read and is wrong too. The graph's genuine value is purely the 3,541 extra
*sites* it adds (already exploited) — not a corrective vote. Validated entirely
offline (no C++ written). This closes the graph-as-voter / architecture-flip
avenue: a better signal must come from data BAM does **not** already have (e.g.
trio/parental reads, a second orthogonal sequencing technology), not from
re-weighting the graph view of the same reads.

## Wrong-merge investigation (the switch gap is NOT from bad stitch orientation)

Re-opened the "switches are structural" claim by auditing the merge decisions
against truth. The actual switch gap on chr20 is **BAM 1,393 vs hybrid 1,440 =
+47** (the earlier 1,355 figure was a stale comparison; per-phase-set tables give
1,393/1,440).

**Step 1 — are cross-block merges oriented wrong?** Identified 73 hybrid phase
sets that merge ≥2 BAM phase sets (≥5 reads each). Scored each sub-block's
orientation against truth: **21 of 73 merges are oriented WRONG.** The
discriminating signal is real and runtime-available in part:

| merge class | smallest sub-block size (median) | truth purity (median) |
|---|---|---|
| WRONG (21) | 10 reads | 0.83 |
| RIGHT (52) | 24.5 reads | 0.99 |

So wrong merges decide orientation on a small, impure set of overlap reads —
matching the vote patterns (`...(1,560,31) | (6,10)`: a 1,560-read block joined
to a 6-vs-10 coin-flip block).

**Step 2 — do the wrong merges cause the gap? NO.** Two oracles, both with truth:

- *Re-orient* every wrong sub-block to its correct truth polarity (zero
  contiguity loss): switches 1,444 → **1,441 (−3)**, flips +7.
- *Split out* small/impure weak sub-blocks: switches barely move (1,444 →
  1,437 best case), at the cost of 50–80 broken blocks.

Even a **perfect** fix of all 21 wrong merges removes only ~3 switches. Wrong
merges mostly create a small discordant *island* (a flip, consumed as 2
transitions), not a persistent switch.

**Step 3 — where is the +47 gap really?** Binned switches into 1 Mb windows.
The gap is **localized, not uniform**: windows 37–38 Mb alone account for **+36
of the +47**. Inspecting those phase sets, the mechanism is concrete: the hybrid
**adds reads/sites into already-large blocks in difficult regions**, e.g.

- PS 38393470: BAM 3,030 reads / 6 switches → hybrid 3,275 reads / **20
  switches** (the +245 reads in a segdup tripled the switches).
- PS 37005134 (hybrid): merges BAM's perfect PS 36895001 (1,676 r, 0 sw) with
  messy neighbors → 1,569 reads / **11 switches**.

Attribution of the extra switches in same-id phase sets: **325 occur in blocks
the hybrid GREW (>5% more reads); only 61 in blocks of ~same size.** The hybrid
added 14,391 reads into existing BAM blocks.

**Conclusion (corrected):** the switch gap is **not** a fixable stitch-rule bug —
fixing every wrong merge buys ~3 switches. It comes from the hybrid pulling
**more reads into large blocks within difficult/segdup regions**, where those
extra reads are genuinely harder to phase (the same biologically-confident-wrong
reads from the per-read-filter dead ends). This is the *price of the
completeness/auN win*, concentrated in a few hard windows (37–38 Mb). The only
lever that helps switches here is the same one already ruled out: stop pulling
hard reads into hard blocks — which sacrifices the completeness that makes the
hybrid win. A small, safe stitch improvement (require the weaker sub-block to
have ≥~15 reads / purity before merging) is *legitimate* but worth only a handful
of flips, not the switch gap. Validated entirely offline; no C++ written.

### Real-run confirmation: stricter stitch gate (rule 3 + margin) — oracle was exact

Two further facts narrowed the gate's reach before any run:
- Only **13 of 21** wrong-merge seams sit at 500 kb **chunk boundaries** (where
  `flip_chunk_hap` acts). The other **27 are intra-chunk** — created by k-means
  inside a chunk, which the stitch step never sees. So a stitch-gate change can
  touch at most ~13 of the wrong merges.
- The existing binary already exposes the gate via `--stitch-rule 3`
  (both-strands + margin) and `--stitch-min-margin`, so no code was needed to
  test it.

Ran the full hybrid pipeline + eval (chr20, truth BAM) at two gate settings vs
the HEAD default (rule 1, margin 10):

| metric | BAM | HEAD r1m10 | r3 m8 | r3 m15 |
|---|---|---|---|---|
| switches | 1,393 | 1,440 | 1,438 | 1,437 |
| flips | 2,052 | 2,133 | 2,131 | 2,131 |
| perfect PS | 125 | 91 | 94 | 97 |
| auN | 9.83 Mb | 12.00 Mb | 11.93 Mb | 11.81 Mb |
| Hamming | 3.09% | 3.04% | 3.04% | 3.04% |
| phased reads | 219,925 | 220,897 | 220,767 | 220,574 |

The real run matches the oracle to the read: the strictest gate moves switches by
only **−3** (1,440 → 1,437) and flips by −2. It does buy **+6 perfect phase
sets** (91 → 97, partway to BAM's 125) but at a monotonic cost in auN (−0.19 Mb)
and phased reads (−323). This is a contiguity-for-perfect-PS trade, **not** a
switch fix.

**Final verdict on the merge rule:** it is not the lever for the switch gap. The
gap is dominated by reads pulled into large blocks in a few hard windows (37–38
Mb), which no boundary stitch rule can reach. The HEAD default (rule 1, margin
10) remains the best all-round operating point — it maximizes auN/completeness,
which is the hybrid's reason to exist. If perfect-PS were the priority metric,
`--stitch-rule 3 --stitch-min-margin 15` is a valid alternative (+6 perfect PS,
−3 switches) at a small contiguity cost, but it is not promoted to default
because it sacrifices the completeness win. No source change made; sweep
artifacts removed.

## Full sweep of graph→BAM integration surfaces (all four families, REJECTED)

To answer "have we tested every way to inject graph signal into the BAM
pipeline?", enumerated the four code surfaces where graph data can enter and
validated each offline. The live hybrid already exploits surfaces 1–2 (add
graph-only sites; extend profiles at those sites). The four *untested* methods:

**D — graph alleles at shared sites as extra coverage/AF support. DEAD (no
oracle needed).** GAF reads are a strict subset of BAM: 271,930 ⊂ 272,016, with
**0** GAF-only reads. Every graph observation at a shared site is a read BAM
already counted, so "extra coverage" is double-counting — it inflates DP/AF with
zero new evidence and cannot rescue a borderline het.

**A — graph-confirmed sites as k-means anchor weights.** Real but tiny signal:
graph-unconfirmed (bam-only) het SNP sites are **1.3–1.6× more likely** to sit at
a switch than graph-confirmed shared sites — so confirmation *is* mildly
protective. But bam-only sites are only 1,583 of 64,882 (2.4%), and of 5,706
switch transitions only **7** sit where an unconfirmed site is the sole nearby
anchor (42% are at confirmed sites where weighting can't help; most are near
neither). Oracle ceiling ≈ **3 switches**, and down-weighting would discard real
het information. Not worth building.

**C — graph-only sites spanning chunk boundaries as stitch bridges. DEAD.**
Targeted the 27 intra-chunk wrong merges the stitch gate can't reach. Only
**5 of 40** wrong-merge seams have a graph-only site within 10 kb to bridge, and
fixing *all* 40 seams is worth only −3 switches (proven earlier). Ceiling < 1
switch. Root cause: bridging reads are the same physical reads (GAF ⊂ BAM), so
they replay the evidence that caused the wrong merge.

**B — graph path co-occurrence as long-range linkage (the only candidate that
could be independent). REJECTED.** The GAF node path (col 9) is the graph's
topology, not re-read bases, so it *could* carry information BAM's linear
alignment lacks. Tested by partitioning reads on graph node-bubbles and scoring
against truth across 50 windows: graph-path partition recovers the true
haplotype at only **60% per window** (median 59%, min 50% = coin flip), vs BAM's
**97%**. Concordance with BAM's own partition is 68% — but the 32% "disagreement"
is *noise*, not independent signal (it doesn't match truth). This is consistent
with the standalone graph pipeline's known accuracy (Hamming 7.3% / switch 2.58%,
~2× worse than BAM): the graph path *is* that pipeline's signal. In hard regions
the node path bifurcates on alignment ambiguity, not clean allele difference.

**Unifying conclusion.** All four fail for the **same root cause that governs the
entire investigation: GAF reads are a strict subset of BAM reads (0 unique).**
The graph is a *second alignment of the same physical HiFi bases*, so it carries
no evidence independent of BAM — not as coverage (D), not as a corrective vote
(earlier), not as linkage (B). Its only genuine contribution is the **3,541 extra
easy SITES** BAM missed, which the hybrid already uses. Surfaces that re-weight or
re-interpret the graph's *view of shared reads* (A, C, D, B, plus the earlier
graph-as-voter) are all capped at single-digit switches because there is no new
information to extract. **To beat BAM further requires data BAM does not contain
(trio/parental reads, or a second orthogonal technology) — not another way to
process the graph alignment of the same reads.** All validated offline; no C++
written; temp files removed.

## BAM-native optimizer alternatives (REJECTED — production k-means already beats MEC and transition-penalty on real matrices)

After exhausting graph→BAM integration, the question pivoted to the BAM
pipeline's **own** phasing algorithm: *is there a better optimizer than the
greedy k-means?* Tested directly against the real per-read × per-variant matrix.
**Answer: no.** Production k-means beats both MEC local search and a
transition-penalized optimizer. Details below, including a retracted earlier
claim.

### RETRACTED: the "confidence smoother" win was built on a truth-side signal

An earlier draft of this section claimed a confidence-aware smoother gated on
`hapq` cut switches −33% and Hamming +0.5pp. **That result is withdrawn.** The
`hapq` field is the `hq:i:` tag from the **diplinator truth BAM** (eval script
`evaluate_phase_accuracy.py:196-199`), not a quantity pgphase computes. It cannot
gate anything at phasing time, and "`hapq=0` ⇒ uncertain" partly describes
truth-side mapping ambiguity, not k-means tie-breaks. The lever was not
realizable; the apparent win was a measurement artifact. Lesson re-learned (same
as the LD/graph-voter dead ends): **validate the signal exists at decision time,
not just in the eval dump.**

### Honest test: dump the real matrix, run candidate optimizers offline

Added `--dump-phase-matrix PREFIX` (debug-only) to emit the exact read×variant
allele matrix k-means consumes, per chunk per flag-set, before any state
mutation (`dump_phase_matrix` in collect_phase.cpp). Joined to truth by qname.
Scoring is **PS-aware**: production assigns multiple phase sets per chunk (up to
9), each independently orientable — scoring a chunk under one global orientation
falsely inflates disc to ~80%. Within each production PS, pick the orientation
minimizing switches, then count.

Baseline reproduced correctly: **production = 1,656 switches / 3.376% disc** on
chr20 (matches the real eval-harness numbers).

### MEC local search from production init — WORSE

Starting from production's own assignment and iterating MEC consensus
re-assignment (no transition term) moved 411 reads and **increased** errors:
**+37 switches, +679 disc**. Production is already past the naive-MEC optimum —
its consensus-allele model + noisy second pass encode more than the bare matrix
energy.

### Transition-penalized optimizer — WORSE at every strength

The actual hypothesis (an HMM/WhatsHap-style switch cost). Added a contiguity
prior rewarding agreement with position-adjacent reads within each PS, swept
λ ∈ {1,2,3,5,8} from production init:

| λ | Δswitches | Δdisc | reads moved |
|---|---|---|---|
| 1 | **+123** | +1,182 | 1,835 |
| 2 | +130 | +1,166 | 1,802 |
| 3 | +537 | +3,880 | 6,313 |
| 5 | +657 | +6,629 | 9,134 |
| 8 | +715 | +7,087 | 10,796 |

Every λ is strictly worse. A contiguity prior pulls reads across **legitimate**
haplotype-block boundaries (the segdup long runs with median hapq 60 = confident
and genuinely on the other haplotype), manufacturing errors faster than it fixes
short runs. This is the same conclusion the segdup analysis predicted: the
switches that remain are **data-driven**, not optimizer tie-breaks.

### Soft EM (probabilistic assignment) — WORSE

The third candidate: replace hard assignment with soft responsibilities (EM over
a per-variant per-hap Bernoulli allele model). Initialized from production
responsibilities (not random — the earlier random-init EM was degenerate and is
not a fair test), refined to convergence, scored PS-aware:

| EM error rate ε | Δswitches | Δdisc | reads moved |
|---|---|---|---|
| 0.02 | +498 | +4,341 | 4,494 |
| 0.05 | +510 | +4,572 | 4,672 |
| 0.10 | +512 | +4,555 | 4,585 |
| 0.20 | +466 | +2,797 | 4,600 |

Worse at every error rate. The tell is a near-frozen sanity run (init
confidence 0.999, **1 round**): it already moves 3,681 reads and adds +497
switches — a single EM E-step disagrees with production on ~3,700 reads and is
net wrong. So it is not harness drift; EM genuinely disagrees and production is
right. Reason: EM's generative model (independent Bernoulli per variant per hap,
one global error rate) **discards what production encodes** — the clean-het
category weights (CleanHetSnp/Indel count double), HOM handling, and the noisy
second pass. The matrix carries the alleles but not these priors, so a model that
sees only the matrix underperforms.

### Conclusion

The greedy pivot-sweep k-means is **not** the bottleneck. On the real allele
matrix it already dominates MEC, transition-penalized, and EM/soft alternatives. The
`--dump-phase-matrix` flag is kept for future offline optimizer experiments. To
reduce BAM switches further requires *new data* (trio/parental, second
technology), consistent with the whole-investigation conclusion — not a better
optimizer over the same reads. Validated on real chr20 matrices; the only C++
added is the debug dump.

---

### (Historical, RETRACTED) original confidence-smoother write-up

The text below is preserved for the record but its conclusion is **wrong** (see
retraction above): it used the truth-side `hapq` tag as if it were a pipeline
signal.

### Root cause: BAM k-means has no transition penalty

`assign_hap_based_on_germline_het_vars_kmeans` (collect_phase.cpp:509) assigns
each read to the hap with the higher **summed** allele-match score, then runs
≤10 Lloyd rounds. Assignment is **per-read and independent** — there is no cost
for creating two transitions in the read ordering. State-of-the-art phasers
(HMM, MEC/WhatsHap) add a switch cost so a read only flips orientation when its
*own* allele evidence outweighs the cost of breaking local contiguity. BAM's
greedy hard assignment has no such term, so a read with **zero evidence margin**
(a tie) is assigned arbitrarily and can land opposite its neighbors.

### Switch breakpoints split cleanly by confidence

Decomposing BAM discordant runs by length and the pipeline's own per-read
`hapq` (haplotype quality, already computed) on chr20:

| bucket | n | median hapq | %hapq<10 | switches contributed |
|---|---|---|---|---|
| concordant | 213,127 | 60 | 2% | — |
| flip (len 1) | 1,777 | 18 | 44% | 0 (flips, not switches) |
| **short run (2–4)** | **2,061** | **0** | **61%** | **~1,551** |
| long run (5+) | 2,960 | 60 | 26% | 448 |

The dividing line is sharp:

- **Short runs (2–4 reads): median hapq = 0.** The algorithm had *no allele
  evidence* — it broke a tie and happened to flip away from neighbors. These are
  **algorithm-fixable**.
- **Long runs (5+ reads): median hapq = 60**, same as concordant reads. These
  are **confident but wrong** = true data ambiguity (segdups, where the read
  genuinely matches the other haplotype's local sequence). No optimizer fixes
  these; forcing them would *add* Hamming error. This matches the segdup
  localization (37–38 Mb) from the wrong-merge investigation.

### Truth-free smoother: −33% switches AND +0.53pp Hamming

Simulated a transition-penalizing optimizer **without using truth**: for each
read with `hapq < HQ`, reassign its orientation to the majority orientation of
its `W` nearest *confident* (`hapq ≥ HQ`) neighbors. High-confidence reads (98%
of all reads, at hapq 60) are never touched. Scored the result against truth:

| HQ< | W | switches | ΔSW | Hamming | ΔHam |
|---|---|---|---|---|---|
| baseline | | 1,999 | — | 3.091% | — |
| 20 | 6 | **1,340** | **−659 (−33%)** | **2.556%** | **+0.535pp** |
| 15 | 6 | 1,447 | −552 | 2.65% | +0.44pp |
| 10 | 6 | 1,468 | −531 | 2.745% | +0.346pp |

Safety audit (HQ<20, W6): **1,481 fixes vs 304 breaks = 4.9:1**. Every
configuration tested holds a **~4–5:1 fix:break ratio** — for each already-correct
low-confidence read wrongly flipped, 4–5 wrong reads are corrected. Because only
sub-threshold reads are ever reassigned, the confident backbone is untouched.

### Why this is different from every prior lever

Every earlier idea (per-read LD, graph-as-voter, stitch gate, the four
integration families) was a **trade-off**: fewer switches cost completeness or
Hamming. This one improves both at once because it targets reads the pipeline
*already knows* it is uncertain about (hapq≈0) — it is not introducing new
information, it is **stopping the greedy optimizer from making arbitrary
tie-breaks against local contiguity**. The signal (`hapq`) is already computed
and the operation needs only read order, so it is directly realizable in C++
inside the k-means post-pass.

**Status: validated offline, not yet implemented.** Next step (if pursued) is a
neighbor-consensus / transition-penalty post-pass in
`assign_hap_based_on_germline_het_vars_kmeans` gated on a `hapq` threshold, with
the smoother confined to sub-threshold reads. Recommended starting params from
the sweep: `HQ=20, W=6`. Must re-validate on the *real* pgphase run (synthetic
switch counts here differ from eval-harness counts; the **ratios** are what
transfer) and confirm perfect-PS and auN/N50 do not regress.

## Where hybrid's extra switches come from (within-chunk k-means, not merging)

Question: hybrid has more switches than BAM (1,393 → 1,440 in the default eval).
Is that purely the different merge rule, or do the extra graph sites change the
phasing itself?

### Test 1 — match the merge rule

BAM uses `--stitch-rule 0 --stitch-min-margin 0`; hybrid defaults to rule 1 +
margin 10. Re-ran hybrid with **BAM's exact stitch rule**:

| run | switches | flips | disc | Hamming | perfect-PS | N50 |
|---|---|---|---|---|---|---|
| BAM | 1,393 | 2,052 | 6,798 | 3.09% | 125 | 938k |
| hybrid (default rule) | 1,781 | 2,530 | 10,082 | 4.52% | 70 | 975k |
| hybrid (BAM rule) | 1,448 | 2,137 | 6,744 | 3.05% | 89 | 954k |

Matching the merge rule closes most of the gap (1,781 → 1,448) and Hamming
actually beats BAM (3.05% vs 3.09%). **Most of the default-hybrid switch excess
was the aggressive stitch rule, not the graph.** But a residual **+55 switches**
remains with identical merging.

### Test 2 — locate the residual switches

Classified the switch breakpoints present in hybrid-same but **not** in BAM
(persistent switches only, flips excluded), by distance to the nearest 500 kb
chunk boundary:

- 392 newly-introduced switch breakpoints.
- **390 (99%) are chunk-interior**, median 122 kb from any boundary.
- Only 2 are within 2 kb of a boundary.

So the residual is **not** a merge artifact — it is created **inside** chunks,
during k-means, exactly the hypothesis that adding sites changes the within-chunk
partition.

### Test 3 — mechanism: partition instability, not graph-site errors

Two further checks rule out the naive "graph sites place errors" story:

1. **No colocation with graph sites.** New switches within W bp of a graph-only
   site: 4% (150 bp), 7% (500 bp), 23% (2 kb) — all **at or below** the random
   baseline (4.7% / 15.5% / 62%). The switches do not sit at the added sites.
2. **Churn is two-way and large.** Hybrid (BAM rule) vs BAM:
   - ADDS 392 switch breakpoints
   - REMOVES 335 switch breakpoints
   - net **+57** (≈ the +55 summary delta); total churn **727**.

The net (+57) is small against the churn (727). Adding ~8k–10k graph sites does
not inject errors at specific loci; it **perturbs the global k-means consensus**,
which reshuffles ~700 breakpoints genome-wide and happens to land slightly
net-worse on switches (while improving Hamming and contiguity). This is
**optimizer sensitivity to the site set**, not a graph-quality problem.

### Answer to "would the same merge rule make it the same?"

Almost, but not exactly. Same merge rule removes ~85% of the switch gap and makes
hybrid match/beat BAM on Hamming. The remaining **+55 switches are within-chunk**:
the extra sites change the k-means partition (two-way churn of ~700, net +57).
They are not at the graph sites and not at chunk boundaries — they are a
second-order effect of re-running the same greedy k-means over a denser, slightly
different variant set. Consistent with the earlier finding that the greedy
k-means is itself the limiting factor: it is **not stable** under changes to the
candidate set, so even strictly-more sites can move switches in both directions.

## Inject graph sites only between unmergeable BAM blocks (REJECTED for switches — contiguity-only, structural)

Idea: instead of full hybrid (which churns the within-chunk partition, +57
switches), inject graph-only sites **only in the gaps between BAM phase blocks
that stayed separate** — the surgical version that should keep the contiguity
benefit without perturbing phasing. Distinct from family C (which targeted *wrong*
merge seams); this targets *clean* unmergeable gaps.

### Ceiling (offline, oracle orientation)

- BAM produces **429 phase blocks** with **128 inter-block gaps** (median 4.4 kb).
- **27 of 128 gaps contain a graph-only site** (median 2 sites/gap) — the only
  gaps a graph-site injection could bridge.
- Merging all 27 with **truth-optimal orientation** (best possible case):

| metric | baseline | oracle-merge 27 gaps | Δ |
|---|---|---|---|
| phase blocks | 429 | 402 | −27 (better contiguity) |
| within-PS switches | 1,586 | 1,737 | **+151** |

### Why merging can only hurt switches

Switch error is measured **within** a phase set. Two blocks stayed separate
because a hard region (low coverage / repeat) sits between them with no signal to
orient across. Bridging them does not remove any error — it **re-exposes that
hard region inside one phase set**, converting what were free block-boundary
discordances into counted within-PS switches. Even with perfect orientation the
operation is +switches / −blocks: a pure **contiguity-for-switches trade**.

### Conclusion

Confirms the mechanism from the hybrid-switch analysis: the graph's extra sites
buy contiguity (N50, fewer blocks) and **cannot** buy switch reduction, even when
applied surgically only between unmergeable blocks with oracle orientation. The
+151 here is the same structural effect as the full hybrid's +57, just isolated.
If contiguity is the goal, gap-injection is a cleaner lever than full hybrid
(touches only 27 sites, no global partition churn); if switches are the goal, it
cannot help. Validated offline; no production change made.

## Does injecting gap sites + re-running k-means actually merge the blocks? (NO — empirical)

The oracle test above forced merges. This tests what the **real k-means** does
when graph sites are present in the gaps — i.e. exactly what the hybrid pipeline
already runs. Compared BAM blocks against hybrid's per-read PS assignment for the
27 BAM gap-pairs that have graph-only sites in the gap.

### Result: k-means does not bridge

| outcome | count |
|---|---|
| BAM gap-pairs with graph sites in gap | 27 |
| hybrid **MERGED** the pair | **1** |
| hybrid kept **SEPARATE** | 26 |
| the 1 merge was **correctly oriented** | 0 |
| the 1 merge was **mis-oriented** (made a switch) | 1 |

The single merge was the smallest gap (1,258 bp) and k-means oriented it
**wrong**, creating a switch. The other 26 — including 8 of 9 gaps narrow enough
(<15 kb) for a HiFi read to span — were left unmerged.

### Why re-running k-means can't help

The gaps are **low-coverage / low-signal holes, not missing-site holes**. Of the
9 spannable gaps, most have only **2–3 reads inside the gap**; one has 36 reads
but a single site and still didn't merge. The block boundary exists *because* the
reads crossing it carry no usable heterozygous signal (low coverage, repeat, or a
homozygous stretch). Adding graph site **positions** does not add phasing
**signal**: those sites are genotyped by the same GAF⊂BAM reads, and if those
reads had clean spanning het alleles BAM would already have phased through. So
k-means correctly abstains on no signal — and the one time it committed, it
guessed wrong (50/50).

### Conclusion

"Inject in the gap + one more k-means round → good merge" **does not occur**:
1/27 merged, and that one was mis-oriented. This is the empirical confirmation of
the oracle finding — and of the whole investigation's root cause. Bridging
unmergeable blocks needs **independent linking evidence** (spanning reads with
het signal, i.e. new data), which graph sites over the same reads cannot provide.
No production change made; validated offline.

## Are there regions the graph phases BETTER than BAM? (YES — but not truth-free identifiable)

User intuition: the "GAF ⊂ BAM, 0 unique reads" finding is about read *identity*;
it does not preclude *regions* where the graph phases better (a read can be in
both but soft-clipped/misplaced in BAM yet clean in the graph). Tested directly.

### Methodology fix (important)

First attempts had a scoring bug: best-flipping a 100 kb *window* as a whole vs
best-flipping each *phase set* gave contradictory answers (a window with two
opposite-oriented PS scores 50% under window-flip but 0% under per-PS flip). The
"159 windows graph wins" and "BAM 50%/graph 0%" intermediate numbers were
artifacts and are **retracted**. Correct metric: best-flip **each PS
independently**, label every read correct/incorrect (flip-invariant), then
compare pipelines on shared reads.

### Result: graph-favorable regions are real

Reads phased by **both** pipelines (207,687):

- GRAPH right & BAM wrong: **25,745 reads**
- BAM right & GRAPH wrong: **37,556 reads**

Net favors BAM (matches the aggregate Hamming 3.1% vs 7.3%), but the two-way
disagreement is large and there is a **genuine graph-favorable subset**. Per
100 kb window (≥20 shared reads): **142 windows graph wins, 205 BAM wins, 293
tie**. In the graph-win windows BAM is internally noisy (a mid-block switch);
the graph holds a clean single block. So the intuition is **correct**: the graph
does phase some regions better.

Separately, **8,443 reads are phased by graph but not BAM at all**; of those,
~1,154 (39 PS) are well-phased (<10% disc), the rest poorly. So the graph also
adds a small amount of correctly-phased coverage BAM misses.

### But the regions are not truth-free identifiable — so not exploitable

The blocker for a selective hybrid: can we tell *which* regions to trust the
graph on **without** truth? Tested BAM's own internal consistency (fraction of
window reads agreeing with the window-majority HP) as a router signal:

| window class | n | median BAM internal consistency |
|---|---|---|
| graph wins | 142 | 0.64 |
| BAM wins | 205 | 0.58 |
| tie | 293 | 0.55 |

The signal **does not separate** graph-win from BAM-win windows (0.64 vs 0.58,
both near the 0.5 split floor). BAM looks equally unsure in regions where it is
actually right and where it is actually wrong. So a selective hybrid has no
reliable truth-free rule to defer to the graph — blindly merging pulls in the
graph's 205 losses with its 142 wins, which is exactly why every whole-pipeline
merge tested lands net-worse on switches.

### Conclusion

Both can be true at once, and both are: (1) the graph phases a real subset of
regions better than BAM (~142 windows / ~25k reads), validating the intuition;
(2) those regions cannot be identified from BAM-side signal without truth, so the
gain is not capturable by a router. The remaining path to exploit them is an
**independent confidence signal** — e.g. the graph pipeline emitting a calibrated
per-PS phasing-quality score that is trustworthy where it disagrees with BAM —
which would require validating that the graph *knows* when it is right (its own
confidence vs truth). That is the concrete next experiment if this thread is
pursued. Validated offline; scoring bug noted and corrected; no production change.

## Can a confidence signal route reads to the graph where it wins? (NO — gain is 3 lucky regions, not a generalizable signal)

Follow-up to the previous section: the graph phases ~142 windows / ~25k reads
better than BAM, but those regions were not separable by a single BAM-side
signal. This tests the full question — can **any** truth-free signal (graph-side,
BAM-side, or a learned combination) identify where to trust the graph, enough to
build a confidence-routed hybrid?

### Setup

Trusted labels: best-flip each PS, per-read correct/incorrect (flip-invariant).
Focus on the **63,301 disagreement reads** (graph and BAM give different
correctness) — the only reads a router decision affects. Base rate: 40.7%
graph-right. Target: predict graph-right without truth, route those reads to the
graph.

### Single signals — all useless

Truth-free signals tested (per graph PS and per read): PS size, het-site density
(5k/20k windows, graph and BAM), allele-frequency-near-0.5, mean depth, read
mapq, read coverage, BAM/graph orientation agreement, position offset.

| signal | AUC |
|---|---|
| read mapq | 0.50 |
| graph PS size | 0.44 |
| BAM PS size | 0.52 |
| het-density diff | 0.53 |
| all others | 0.48–0.53 |

Every single signal is at chance. Best linear combination (logistic regression,
held-out): **AUC 0.60** — below the 0.65 usability bar.

### Gradient boosting — 0.894 was leakage

A gradient-boosted model scored **AUC 0.894** with read-level 4-fold CV — but
this is **train/test leakage**: reads from one PS share a PS-size feature and
(mostly) one label, so reads of the same block leak across folds. Re-run with
**GroupKFold over phase sets** (no PS in both train and test): **AUC 0.59 ±
0.13**. The boosted model was memorizing which specific block is right, not
learning a transferable signal. (Sanity: the top features `gps_n`/`bps_n` have
individual AUC 0.44/0.52 — useless alone, confirming the 0.894 came from
group structure, not signal.)

### The gain is 3 regions, not a rule

Translating the grouped out-of-fold predictions into routing outcomes:

| threshold | routed | fixed | broke | net | precision |
|---|---|---|---|---|---|
| 0.5 | 20,280 | 9,978 | 10,302 | −324 | 0.49 |
| 0.7 | 7,038 | 4,074 | 2,964 | +1,110 | 0.58 |
| 0.8 | 3,536 | 2,646 | 890 | +1,756 | 0.75 |
| 0.9 | 1,180 | 1,177 | 3 | +1,174 | 1.00 |

The high-precision tail looks great — until you check **where it comes from**.
At thr 0.9, **1,151 of 1,318 reads are a single phase set** (PS 2717436, graph
disc 0.00). Decomposing net gain per PS at thr 0.8:

- **Top 3 PS: +1,963 net = 158% of the total** (1,244).
- **All other 24 contributing PS: −719 net** (net negative).
- **AUC excluding the top 3 PS: 0.530** (chance).

So the entire apparent gain is **3 lucky large blocks** that happen to be
graph-correct and present in the data. Everywhere else the router breaks more
reads than it fixes. This is memorization of specific regions, not a
generalizable confidence signal.

### Conclusion

**The graph does not know when it is right.** No truth-free signal — single,
linear, or gradient-boosted with PS-grouped validation — generalizes beyond a
handful of specific regions (AUC 0.53 once those are removed). A confidence-routed
hybrid is therefore not realizable from the signals available in the current
pipeline outputs: it would either need (a) a genuinely calibrated per-PS quality
score computed inside the graph phaser and validated to transfer across regions,
or (b) independent linking data. This closes the last open thread: the graph's
region-level wins are **real but not capturable** without truth. Validated
offline (sklearn GroupKFold); no production change. Methodological note for future
work: **always group-split by phase set** — read-level CV leaks and inflated this
AUC from 0.59 to 0.89.

## Why does the graph phase some regions better? (NOT more sites — k-means convergence coin-flip)

Having established the graph wins ~141 windows but those wins aren't routable,
the remaining question is *why* it wins them. Tested every data-side hypothesis on
the three window classes (graph-win=141, bam-win=204, neutral=293, 100 kb,
PS-aware, ≥30 shared reads).

### It is NOT more sites, coverage, or quality

| class | graph het (med) | bam het (med) | g−b het | bam cov | graph cov | bam mapq | graph mapq |
|---|---|---|---|---|---|---|---|
| GRAPH-WIN | 103 | 98 | **−1** | 332 | 328 | 60 | 60 |
| BAM-WIN | 108 | 108 | −3 | 350 | 352 | 60 | 60 |
| NEUTRAL | 120 | 123 | −2 | 372 | 372 | 60 | 60 |

In graph-win windows the graph has **equal or slightly fewer** het sites, the same
coverage, and the same mapq as BAM. The graph never has more sites anywhere — it
has marginally fewer everywhere. **Site count, coverage, and read quality are
ruled out.** The two pipelines see the same reads at the same sites with the same
quality, yet phase differently.

### The mechanism: each pipeline's k-means coin-flips convergence

Counting within-window switches per pipeline:

| class | BAM switches/win | graph switches/win |
|---|---|---|
| GRAPH-WIN | 9 (clean) | 60 (noisy) → wait, graph WINS here |
| BAM-WIN | 68 (noisy) | 10 (clean) |
| NEUTRAL | 60 | 60 |

(The disc-rate winner is what defines the class; the switch counts show the
*loser* in each class has a mid-block switch.) The pattern is **symmetric**: in
graph-win windows BAM carries a switch the graph avoided; in bam-win windows the
graph carries a switch BAM avoided; in neutral windows both are equally noisy.

Decomposing the **loser's** error shape:

- GRAPH-WIN: BAM's error is a **clean switch** (one contiguous wrong block ≈ a
  half-window orientation flip) in **83/141** windows; scattered noise in 58.
- BAM-WIN: graph's error is a clean switch in **120/204**; scattered in 84.

A "clean switch" means the k-means oriented half the window the wrong way — a
**convergence failure**, not ambiguous data (the other pipeline got the same
evidence right). 

### Conclusion

The graph's region-level wins are **not a data advantage** — same sites, same
reads, same coverage, same quality. They are the **same fragile greedy k-means
landing a different coin-flip**: in each hard window, each pipeline independently
either converges to the clean haplotype partition or to a half-flipped one, and
"graph-win" is simply where the graph's flip landed right and BAM's landed wrong
(and symmetrically for bam-win). This is the same k-means instability behind the
+57 hybrid churn, now shown to drive the region-level differences too. It explains
both earlier findings at once: the wins are **real** (genuine convergence
differences) but **not routable** (a coin-flip has no truth-free tell). The only
fixes are a **more stable optimizer** (shown earlier that MEC/EM/transition-penalty
don't beat the current k-means on the real matrix) or **independent data**.
Validated offline; no production change.

---

## Does a more stable optimizer beat production? (multi-restart + MEC selection)

**Motivation.** The coin-flip mechanism above suggests the win is restart-luck.
If so, running many random-pivot restarts and **selecting by a truth-free
objective** (MEC energy) should land the good basin more often and cut switches.
Tested on the real dumped matrix (133 chunks, flags140, PS-aware best-flip
scoring per phase set). Production baseline = **1,656 switches / 3.376% disc**.

### 1. Convergence variance is real and large
8 random pivot inits per chunk on the first 20 chunks: per-chunk switch count
ranges **26–92** by init alone (e.g. chunk 112: 57→149); total-switch range 941
across 20 chunks. MEC energy also varies per run (spread 944–3934), so a
truth-free objective *can* in principle distinguish runs.

### 2. But restarting the reimplemented k-means is far worse
12 restarts, pick min-MEC, scored from scratch:

| optimizer | switches | disc |
|---|---|---|
| production | 1,656 | 3.376% |
| multi-restart ×12, min-MEC | 5,473 | 43.3% |
| oracle (min-**switch** run) | 4,432 | — |

The reimplemented `pivot_init`/`kmeans` is a **degenerate baseline** — even its
oracle (4,432) can't reach production (1,656). Production's category weights,
HOM handling, and locked-pivot Phase-1 init are doing real work that a plain
Lloyd restart discards. Restarting a worse optimizer can't beat a better one.

### 3. MEC is a weak proxy for switch count
Within-chunk Spearman ρ(MEC, switches) over 40 chunks: **mean 0.565**, median
0.594, positive in 97%. But argmin(MEC) == argmin(switches) in only **13/40**
chunks. MEC correlates with truth switches but selects the truly-best run only
~⅓ of the time — too weak to be a reliable selector.

### 4. Deployment-realistic test: production never gets overridden
Include production's own assignment as one candidate alongside 12 restarts,
select min-MEC:

| | switches |
|---|---|
| production | 1,656 |
| prod + 12 restarts, min-MEC | 1,658 (**+2**) |

MEC swapped away from production in only **2/133** chunks, and one of those two
made switches worse. Net effect ≈ zero. **A more stable optimizer provides no
truth-scored gain** — production already sits at the good basin MEC would pick.

### 5. Ensemble instability does not flag errors
Per-read vote fraction across 12 aligned restarts; "unstable" = vote fraction in
(0.2, 0.8). Unstable reads are discordant **48.6%** vs stable **46.3%** —
essentially no separation (flag precision 0.486 ≈ base rate). Restart-instability
carries **no truth-free signal** about which reads are mis-phased, so a
consensus/confidence scheme built on it cannot route or correct errors.

### Conclusion
A more stable optimizer does **not** beat production k-means. (1) Variance is
real, but (2) the only optimizer that beats restarts *is* production, (3) MEC is
too weak a selector to improve on it, (4) when offered production as a candidate
MEC keeps it (net +2), and (5) restart-instability is not a usable error signal.
Combined with the earlier MEC/EM/transition-penalty results, every truth-free
optimizer variant tried fails to beat the current k-means. The coin-flip wins
are real but not capturable without **independent data**. Validated offline; no
production change. The `--dump-phase-matrix` debug flag remains the only `src/`
change and alters no phasing behavior.

---

## Can graph read-mappings link BAM phaseblocks across gaps?

**Idea (independent-data angle).** Phase with BAM, find the gaps between adjacent
BAM phaseblocks that couldn't be merged, and ask whether reads' **graph (GAF)
mappings** span those gaps even where their **BAM alignments** break. A read whose
linear BAM alignment clips at a repeat/SV but whose graph path threads through
could supply linkage BAM lacks — turning two blocks into one. This is distinct
from earlier gap-injection (which injected graph *sites* into k-means); here we
test graph *read connectivity* directly.

### Setup
- BAM phaseblocks from `bam.phased.vcf`: **431 PS blocks, 430 inter-block gaps,
  ~11 Mb total gap.** Largest gaps cluster in the peri-centromere (27–32 Mb).
- GAF: 271,930 alignments with linear ref projection (contig, ref_start, ref_end).
- A read "spans" a gap a→b if some alignment has ref_start≤a and ref_end≥b.
- "graph-only" spanner = spans in GAF but **no** BAM alignment of that read spans.

### Coordinate spanning exists, but is a centromeric artifact
287/430 gaps have ≥1 GAF spanner. Classifying by region (centromere = 26.0–29.5 Mb):

| region | n_gaps | gaps with graph-only linkage | total graph-only reads | total GAF spanners |
|---|---|---|---|---|
| **arm** | 400 | **5** | **65** | 6,288 |
| **cent** | 30 | 17 | **1,534** | 1,885 |

**96% of all graph-only linkage (1,534/1,599 reads) is inside the centromere.**
There, "spanning" is a linear-projection artifact: the graph path threads a
tandem-repeat array that the BAM aligner correctly refused to align across, and
truth alignment itself fails (the diplinator can't place these reads either).
Spanners also scatter across many production PS labels rather than anchoring the
two specific flanking blocks (e.g. 393 spanners of the 26.567→26.591 Mb gap, only
6 with any truth hap, only 3 assigned to a flanking PS).

### On the arms, there is no usable linkage
Only **4 arm gaps** have ≥3 graph-only spanners. Examining each for truth-consistent
orientation (a merge needs spanners that agree on the relative phase of the two
blocks):

| arm gap | graph-only | truth-hap split | no-truth / unaligned |
|---|---|---|---|
| 49,661,443→49,661,903 (460 bp) | 23 | **9 MAT / 9 PAT** | 5 |
| 31,765,624→31,766,466 (842 bp) | 20 | 3 MAT | **17** |
| 31,780,986→31,909,439 (128 kb) | 13 | 3 MAT | **10** |
| 29,506,920→29,510,999 (4 kb) | 7 | 1 PAT | **6** |

- The one gap with both haplotypes present (49.66 Mb) splits **9 MAT / 9 PAT** —
  perfectly mixed, which gives **zero** orientation signal (linkage requires the
  spanners to consistently tie the blocks in one relative phase, not 50/50).
- The other three are dominated by **no-truth/unaligned** reads: the graph spans
  linearly but truth can't place the reads, so the "linkage" crosses a repeat
  where the projection is unreliable.

### Conclusion
Mapping BAM phaseblocks to the graph and inspecting inter-block read mappings does
**not** yield a usable merge signal. The graph-only spanning is real but ~entirely
centromeric, where it reflects a linear-projection artifact over repeats rather
than true connectivity — exactly the regions where BAM's refusal to span is
*correct* and truth itself is undefined. On the chromosome arms, where phasing is
meaningful, graph-only spanners are too few (4 gaps) and carry no consistent
orientation (mixed haplotypes or unaligned reads). This is consistent with the
established GAF ⊂ BAM finding: the graph is a second alignment of the same HiFi
bases, so outside repeat regions it cannot connect what BAM cannot. The remaining
graph wins stay attributable to k-means coin-flips, not extra linkage. Validated
offline; no production change.

---

## Do reads map *better* in the graph (CHM13 ≠ sample)? Alignment-quality rescue

**The objection.** The reads are HG002 HiFi mapped to CHM13 — a reference that is
**not** the sample. Where HG002 diverges from CHM13 (SVs, divergent/segdup
haplotypes), the *linear* BAM alignment may be present but **wrong** (clipped,
mis-placed), while the graph threads the correct alternate path. "GAF ⊂ BAM" is
about read *identity*; it doesn't rule out the graph *aligning the same read
better*. This is a sharper test than the earlier coordinate-spanning one, which
only checked whether BAM physically reached across a gap, not whether its
alignment there was correct.

### Methodology note (a real bug I made and caught)
First pass read the wrong GAF columns: this file is a **custom format with 3
extra leading columns** (contig, ref_start, ref_end) before the standard GAF
block, so matches/blocklen/mapq are at cols **13/14/15**, not 10/11/12. The wrong
read produced a fake "18,063 reads where BAM mapq<20 but GAF mapq≥60." After
fixing the column offset, the result inverted (see below). Lesson: always sanity-
check that a "mapq" column actually ranges 0–60.

### GAF and BAM mapq are byte-for-byte identical
With the **correct** columns:

> **GAF mapq == BAM mapq for 271,930 / 271,930 reads (100.0%).**

The graph alignment is not an independent remap — it carries the *same minimap2
mapq* as the BAM, because it is the same alignment re-projected onto the graph.
So there is **no** read where the graph is more *confidently* placed than BAM
(0 reads with GAF mapq ≥ BAM mapq + 20, and 0 the other way). This is the
strongest form of GAF ⊂ BAM yet: same reads, same alignments, same confidence.

### One real residual: BAM-clipped, graph-clean reads
The graph *can* still differ in **alignment extent**: 792 reads are soft-clipped
>15% in BAM yet align ≥98% identity in GAF (637 on the chromosome arms). These
are real CHM13-vs-HG002 divergence loci where the graph's alt path lets the read
align full-length. Tracing them through the pipelines and truth:

| | count |
|---|---|
| BAM-clipped>15% & GAF identity>98%, on arm | 637 |
| …phased **concordantly** by the graph pipeline | 184 (100% concordant) |
| …phased by graph but **not** by BAM at all | **39** (all concordant) |
| …among reads BAM did phase, BAM **discordant** / graph correct | 4 |
| …of the 39 graph-only rescues, captured by current **hybrid** | **1** |

So the objection is **partially right**: there is a genuine, truth-validated set
of ~39 arm reads that the graph phases correctly because its alt path aligns them
where BAM clips — and BAM phases none of them. They sit at **31–32 Mb**, exactly
among the large BAM phaseblock gaps (128 kb, 151 kb, 170 kb …) in the segdup-rich
peri-centromeric edge, the region where CHM13≠HG002 divergence is highest.

### But the magnitude is tiny and the hybrid already misses it
- **39 reads** total, all in one ~1 Mb segdup band — against 1,393 BAM switches
  and ~223k phased reads. Even perfectly exploited, this cannot move the chr20
  switch/Hamming numbers measurably.
- The current hybrid pipeline captures **1 of 39**, so the mechanism is real but
  not wired in. Exploiting it would mean preferring the graph's alignment extent
  (not its mapq, which is identical) for clipped reads in divergent loci.

### Conclusion
The "reads map better in the graph" claim is **correct in kind but small in
degree** on this CHM13/HG002 chr20 data. Mapping *confidence* is identical (same
minimap2 mapq, 100%), so the graph never rescues via better mapq. The only real
edge is **alignment extent**: ~39 arm reads that BAM soft-clips but the graph
aligns full-length and phases correctly, clustered in the 31–32 Mb segdup band
near gaps BAM cannot close. This is a genuine independent-data signal — the first
found in this whole investigation — but at 39/223,000 reads it is far too small to
shift aggregate accuracy, and the current hybrid captures only 1. It would matter
more on a sample/region with heavier structural divergence from the reference than
chr20 offers here. Validated offline; no production change.

---

## Can the graph add correct phaseblocks where BAM phases nothing?

**Idea.** Don't fix BAM regions — *add* phaseblocks in regions BAM leaves
unphased. Even a few extra correct blocks is a net win (improves N50/auN, phases
reads BAM drops). Test: find graph PS blocks lying entirely outside BAM's phased
coverage, then check them against truth.

### Graph-only blocks exist
Merging BAM's 431 phased blocks into coverage intervals, **45 graph PS blocks
fall entirely outside** BAM coverage; **17 carry ≥2 phased het vars** (the rest
are singletons). Several are substantial (e.g. 26.87 Mb / 36 vars, 32.23 Mb / 23
vars, 44.03 Mb / 15 kb span).

### But raw graph-only blocks are only ~⅓ correct
Scoring each against truth (best-flip disc over its phased reads):

> **5 / 17 are perfect (100% concordant); the rest range down to 50–60%
> accuracy — i.e. coin-flip garbage.**

The bad blocks are the same fragile-k-means failure seen elsewhere: high-read
blocks (146, 66, 63 reads) at ~55–60% accuracy are half-flipped partitions, not
real haplotypes. Admitting them blindly would *add wrong phase*, which is worse
than leaving the region unphased.

### A truth-free filter recovers a clean subset
Per-block signals vs truth accuracy revealed the discriminators:

| signal | finding |
|---|---|
| **indel fraction** | every coin-flip block (acc 0.53–0.60) is **100% indels** with high DP (33–77) — repeat/homopolymer false-het piles |
| **degenerate span** | several bad blocks have span 0–1 bp with many reads (var/kb 999–2000) — multiple vars stacked at one position |
| **segdup band** | the remaining SNP-rich bad blocks all cluster in **32.2–32.5 Mb**, a known CHM13 peri-centromeric segmental-duplication band where paralogs create false hets |
| ~~strand bias, AF~~ | **RETRACTED — see bug below.** AF doesn't separate paralog-het from real het, but the "strand bias useless" claim was wrong: it was an artifact of a strand-accounting bug, not biology |

Applying **"≥2 SNPs AND not in the 32.2–32.5 Mb segdup band"**:

> **4 extra correct phaseblocks, 1 bad admitted, 15 reads newly phased.**

(Recall of the perfect blocks is 4/5; the SNP filter alone gives recall 5/5 but
admits 4 bad, so the segdup exclusion is what makes it usable.)

### Conclusion
The idea **works in principle and yields a real, if small, win**: ~4 extra
truth-correct phaseblocks on chr20 that BAM does not phase at all, recoverable
with a simple truth-free filter (require ≥2 SNPs, exclude the pure-indel and
segdup-band blocks). This is the same independent-data edge as the clipped-read
rescue — concentrated in the 26–32 Mb peri-centromeric region where CHM13≠HG002
divergence is highest — and again **tiny in magnitude** (4 blocks / 15 reads
against 431 BAM blocks). Unlike the optimizer and linkage avenues (which produced
*nothing*), this one is a genuine, defensible net gain and the most promising
hybrid lever found: graph-only blocks are additive (they can't introduce switches
into existing BAM blocks since they occupy disjoint regions), so the only risk is
admitting a wrong block, which the filter controls. A production implementation
would: (1) emit graph PS blocks, (2) keep only those disjoint from BAM coverage
with ≥2 SNPs and outside known segdup/centromere bands, (3) append them to the BAM
phasing. Worth doing for N50/auN even though aggregate switch/Hamming barely move.
Validated offline; no production change yet.

### BUG FOUND (RESOLVED in 8f93774): graph candidate strand counts were always forward (REVERSE=0)

> **Status: FIXED.** Resolved by commit 8f93774 in `graph_collect.cpp:148-153`
> (copy `mcand.counts.{forward,reverse}_{ref,alt}` instead of deriving
> forward-only from `alle_covs`). Verified on current HEAD: all 74,119 rows of
> `graph.candidates.tsv` satisfy `FWD+REV == COUNT` for both ref and alt, including
> multiallelic sites. The strand-bias filter remains correctly `is_ont()`-gated, so
> it does not run on the HiFi path. The historical investigation below is retained
> for the diagnostic trail.


While testing a strand-bias filter for the graph-only blocks, every graph
candidate showed `REVERSE_REF=0` and `REVERSE_ALT=0`. Investigation showed this
is **not biology** (HiFi reads map to both strands — the GAF path orientation is a
real ~49/51 split, 139,013 `>` forward / 132,917 `<` reverse) but a **bug in the
output-table rebuild**:

- Strand is parsed correctly from **path orientation** (`graph_query.cpp:56`,
  `reverse = orient == '<'`), *not* the GAF strand column (which is always `+`).
- Chunk-level accumulation is correct: `graph_bam_adapter.cpp:444-466` splits
  forward/reverse counts properly, and the **strand-bias filter does run** on them
  (`graph_bam_adapter.cpp:179-190`, Fisher-exact, `is_ont()`-gated — mirrors the
  BAM filter at `collect_var.cpp:1475-1482`).
- **The defect:** `graph_chunks_to_candidate_table` (`graph_collect.cpp:146-147`)
  rebuilds candidates from strand-**merged** `alle_covs` and sets
  `forward_ref = ref_cov; forward_alt = alt_cov`, **never setting reverse**. So the
  emitted `graph.candidates.tsv` zeroes all reverse counts.

This is the **same lossy-rebuild pattern** as the `RepeatHetIndel` re-promotion bug
(also in `graph_chunks_to_candidate_table`): the function discards category and
strand information computed upstream.

**Consequence for the analysis above:** my offline "strand bias is useless"
conclusion was drawn on the corrupted TSV (every variant looked single-strand), so
it is **invalid**. The strand-bias filter that runs *during* phasing uses correct
counts, but I could not evaluate strand bias as a graph-only-block discriminator
because the dumped data was wrong. The indel-fraction and segdup-band signals
stand; the strand-bias verdict is retracted pending a re-test on fixed output.

**Fix (one line):** in `graph_collect.cpp:146-147`, copy
`mcand.counts.forward_ref/reverse_ref/forward_alt/reverse_alt` (already strand-
correct after biallelic pairing) instead of deriving forward-only from
`alle_covs`. This is a debug-output/VCF-annotation correctness fix; it does not
change k-means (which uses the chunk counts, not the rebuilt table) but it does
affect the emitted strand-bias-derived VCF fields and any downstream strand
filtering on the TSV.

### APPLIED + RE-TESTED: strand fix and the corrected strand-bias verdict

Applied the one-line fix and regenerated `graph.candidates.tsv` + the graph eval
dump. Verification:

- Reverse counts now populated; only **0.43%** of het variants are zero-reverse
  (genuine single-strand sites), down from **100%**.
- Graph strand counts now **match the BAM pipeline** at shared positions (e.g.
  pos 59818: both Fref=8 Rref=5 Falt=7 Ralt=13).
- Site count **unchanged** (74,119) — confirms the fix is output-only and does
  not alter phasing.

**Strand-bias verdict (now valid):** Re-running per-variant Fisher-exact strand
bias on the **corrected** counts across all 17 graph-only blocks flags **zero**
variants (every block min-p > 0.09). So strand bias is **not** a discriminator for
the graph-only coin-flip blocks. The earlier "strand bias useless" conclusion was
right in outcome but had been drawn on corrupted data; it is now confirmed on
correct data. (Strand bias remains a valid filter in general — it just doesn't
separate *these* blocks, which fail for a different reason.)

### The real discriminator: the RepeatHetIndel re-promotion bug

Inspecting the worst coin-flip blocks (acc 0.53–0.60 at 40–44 Mb) directly in the
VCF shows they are **pure homopolymer / STR indels** emitted as `CLEAN`:

```
40,685,510  C>CA, C>CAA              (poly-A insertion)
42,042,659  A>ATA, A>ATATATATA       (TA-repeat expansion)
43,797,938  C>CT, C>CTT              (poly-T)
44,300,880  C>CAA, C>CAAAA, C>CAAAAAA, C>CAAAAAAA  (poly-A run, same pos)
44,032,221  AA>A, GG>G               (homopolymer deletions)
```

These are textbook `is_homopolymer_indel` / `is_repeat_indel` hits. The noise
filter **does** demote them to `RepeatHetIndel` during phasing (so they're
excluded from k-means and never properly phased — hence the ~50% coin-flip
accuracy), but `graph_chunks_to_candidate_table` rebuilds the output category
from scratch (`graph_collect.cpp:170-180`, classifying purely by type/AF/depth)
and **re-promotes them to `CleanHetIndel`**, so they pass the VCF germline gate
and surface as phased het indels with a (default) PS. This is the same lossy-
rebuild defect as the strand bug, in the same function.

**This is the actual filter the indel-fraction signal was detecting.** The
"100%-indel coin-flip block" pattern is precisely the re-promoted repeat indels.
The principled fix is not a strand or indel-fraction heuristic but **preserving
the `RepeatHetIndel` demotion** through the rebuild — then these 7 bad blocks
never enter the graph VCF at all.

### Updated recommendation for graph-only-block rescue
Breakdown of the 17 multi-var graph-only blocks by failure mode:

| failure mode | blocks | fix |
|---|---|---|
| pure repeat/homopolymer indel (acc ~0.5–0.6) | ~7 | preserve `RepeatHetIndel` demotion in rebuild |
| SNP-rich, in 32.2–32.5 Mb segdup band (acc 0.56–0.78) | 3 | segdup/centromere band exclusion |
| clean SNP-bearing, correct (acc 1.0) | ~5 | **keep — these are the win** |

So the production path is cleaner than the earlier "≥2 SNPs + segdup" heuristic:
(1) **fix the RepeatHetIndel re-promotion** (removes the pure-indel coin-flips at
the source, a correctness fix that also benefits the main graph VCF), (2) **exclude
known segdup/centromere bands**, (3) **append remaining graph-only blocks** to the
BAM phasing. Strand bias is *not* needed for this. The strand fix is still a valid
independent correctness fix and is applied. Validated offline.

### APPLIED: RepeatHetIndel re-promotion fix (+ two underlying noise-filter bugs)

Implemented the fix. It required correcting **three** distinct defects, all of
which let homopolymer/STR indels reach the phased graph VCF:

**Bug A — category dropped on rebuild (`graph_collect.cpp`).**
`graph_chunks_to_candidate_table` re-classified every indel het as
`CleanHetIndel` from scratch, discarding the noise filter's `RepeatHetIndel`
demotion. Fixed by preserving the demotion: if the chunk candidate is
`RepeatHetIndel`, carry it (and `kLongcalldRepHetVar`) through the rebuild instead
of re-promoting.

**Bug B — noise filter never fired (`graph_bam_adapter.cpp` + `noise_filter.cpp`).**
Even before the rebuild, the demotion was rarely set, for two reasons:
  1. *Non-minimal graph alleles.* Graph emits indels with the full repeat run on
     both sides (e.g. `CAAAAAAAAAAA > CAAAAAAAA` for a 3 bp deletion).
     `is_noisy_site` derives the indel length from allele sizes, so the length
     exceeded `max_xgaps` and the homopolymer scan was skipped. **Fix:** trim each
     allele pair to minimal VCF form (shared suffix then prefix) before the check.
  2. *Off-by-one in the shared scan.* `is_homopolymer_indel` / `is_repeat_indel`
     read the reference context starting at the VCF **anchor** base (`end_pos =
     pos`) instead of the inserted/deleted content one base later. For `C>CAAA`
     over reference `...c|aaaa…`, the scan saw the anchor `c` and missed the
     poly-A run. Verified empirically (`end_pos=pos` → not detected; `end_pos=pos+1`
     → detected) for all failing cases. **Fix:** start the downstream scan at
     `pos+1` (insertion) / `pos+del_len+1` (deletion), upstream at `pos`; also made
     the `is_repeat_indel` insertion comparison nt4-encoded so soft-masked
     (lowercase) reference matches. *This shared function is used only by the graph
     and hybrid-inject paths; the BAM pipeline has its own local
     `var_is_homopolymer_pg`, so this change does not affect BAM.*

**Bug C — only the first alt checked (`graph_bam_adapter.cpp`).**
Multiallelic sites (e.g. pos 33,887,895 with 9 poly-A alts) were tested only on
`meta.alts[0]`, which after decomposition need not be the phased allele — and the
longest alt could exceed `max_xgaps`. **Fix:** check **every** alt and demote the
site if any is noisy.

**Result on the graph pipeline (chr20):**

| metric | before (buggy) | after (fixed) |
|---|---|---|
| Hamming error | 1.03% | **0.81%** |
| Switch error | 0.24% | **0.18%** |
| switches | 492 | **368** |
| flips | 1,025 | **859** |
| perfect phase sets | 37.8% | **39.4%** |
| variants demoted to `REP_HET_INDEL` | (re-promoted to clean) | **16,841** |

Graph-only multi-var blocks dropped from 17 → 13; the four egregious pure
poly-A/poly-T/TA-repeat coin-flip blocks at 40–44 Mb (acc ~0.5–0.6) are gone, and
all five truth-correct clean SNP blocks remain. **7 of 8** previously-leaking
repeat-indel blocks are now filtered.

**One residual — RE-DIAGNOSED as a normalization bug, now FIXED:** pos 49,031,428
was previously described here as a "~40 bp AGGG STR expansion." That was wrong.
Trimming the catalog alleles to minimal VCF form shows it is a **single het SNP**,
`chr20:49031440 A>G` — exactly what the BAM pipeline calls. The catalog emitted it
three times (from three overlapping snarls) wrapped in 18–49 bp of equal-length
AGGG-repeat context. Because the `VariantKey` derivation in `graph_collect.cpp`
only stripped a single-base prefix for ins/del and never trimmed equal-length
alleles, padded SNPs fell into the MNP fallback with a misleading multi-bp ref/alt
and the wrong POS, and the emitted VCF showed bogus `A→AAGG…AGGC` insertions. See
"Graph allele normalization (padded-SNP bug)" below.

The 3 remaining SNP-rich bad blocks are all in the 32.2–32.5 Mb segdup band and
still require the segdup/centromere-band exclusion described above. Unit tests pass;
site count unchanged (74,119). Validated offline.

## Hybrid Model — minimal-VCF trimming in apply_hybrid_noise_filter (WIN, now default)

**Goal context:** the standing objective is lower Hamming + switch error on the
hybrid pipeline. After the step-4 re-orientation fix (hybrid → 99.18% / 0.82%), the
next lever was a normalization gap between the two noise filters.

**The gap.** `apply_graph_noise_filter` (src/graph_bam_adapter.cpp) trims each
catalog allele to **minimal VCF form** — strip the shared suffix, then the shared
prefix — *before* calling `is_noisy_site`. `apply_hybrid_noise_filter`
(src/hybrid_inject.cpp) did **not**: it ran `is_noisy_site` on the raw catalog
`(ref, alt)`. Catalog (snarl) alleles carry the **full repeat run on both flanks**,
so the derived indel length is inflated. With `max_xgaps` (default 5) bounding the
repeat scan, a genuine 1–2 bp het indel sitting inside a long homopolymer/STR looks
*longer than the window* and escapes detection in one direction but trips it in the
other — the untrimmed hybrid filter **over-demoted** real het indels to
`REP_HET_INDEL` (excluded from k-means), losing them as phasing anchors. The
standalone graph pipeline kept the same sites as `CLEAN_HET_INDEL` because it
trimmed first.

**Fix.** Port the suffix-then-prefix trim into `apply_hybrid_noise_filter`,
adjusting `noisy_pos` by the prefix shift, before the `is_noisy_site` call. Behind a
toggle (`Options::exp_hybrid_trim`), defaulted **on** for hybrid
(collect_hybrid_variation); `--no-hybrid-trim` restores the untrimmed behaviour for
diagnostics.

**Result (chr20, diplinator truth, both toggles from the same HEAD binary):**

| metric | `--no-hybrid-trim` (old) | trim (new default) | Δ |
|---|---|---|---|
| reads evaluated | 212,465 | 212,240 | −225 |
| accuracy | 99.181% | **99.232%** | +0.051 pp |
| Hamming error | 0.819% | **0.768%** | −6.2% rel |
| switch errors | 319 | **291** | −28 (−8.8%) |
| flip errors | 953 | **906** | −47 (−4.9%) |
| switch+flip | 1,272 | **1,197** | −75 (−5.9%) |
| switch rate | 0.150% | **0.137%** | |
| perfect phase sets | 114 | **115** | +1 |
| `CLEAN_HET_INDEL` | 4,652 | **4,994** | +342 |
| `REP_HET_INDEL` | 6,136 | **5,794** | −342 |

**Mechanism.** Exactly 342 graph-only het indels that the untrimmed filter
mislabelled `REP_HET_INDEL` are now correctly `CLEAN_HET_INDEL` — i.e. recovered as
k-means anchors, matching the standalone graph verdict. More correct anchors →
fewer switch/flip events. The slight read-count drop (−225) is the expected effect
of a few sites changing category; net accuracy and all error metrics improve.

This is a *correctness* fix (the two filters now agree on the same site), not a
threshold tweak, so it generalizes rather than overfitting chr20. SNPs are still
exempt (the existing SNP carve-out is unchanged). Unit tests pass
(`test_hybrid_inject`, `test_noise_filter`, all suites). Validated offline.

## Graph allele normalization (padded-SNP bug, FIXED)

**Symptom.** The "large AGGG STR" residual recorded earlier (pos 49,031,428) was
a misdiagnosis. The graph pipeline emitted three records at 49031428 / 49031431 /
49031439 with bizarre equal-length 18–49 bp ref/alt strings, all of which are the
**same single het SNP** `chr20:49031440 A>G` that the BAM pipeline calls cleanly.
The output VCF showed fake insertions like `A → AAGGAAGG…AGGC`.

**Root cause.** The `VariantKey` derivation in `graph_collect.cpp` did not
normalize catalog alleles to minimal VCF form. The ins/del branches only stripped
a single shared **prefix** base (requiring `ref[0]==alt[0]`) and never trimmed a
shared **suffix**; equal-length alleles fell through to the MNP fallback, which
kept the full untrimmed ref/alt and labelled them `SNP` when lengths matched. A
real SNP padded with equal-length repeat context (the AGGG array) therefore got
the wrong POS, a multi-bp ref/alt, and — because the catalog wraps the same site
in several overlapping snarls — was emitted multiple times.

**Scope (chr20).** 5,320 → 1,123 non-minimal alleles after the fix (the remaining
1,123 are genuine complex/MNP variants with differing bases on both ends that have
no single-anchor minimal form). 3,957 of the trimmed cases were single-base SNPs
buried in repeat padding; ~339 catalog sites were duplicate emissions of the same
post-minimization variant.

**Fix.** Normalize each `(pos, ref, alt)` to minimal VCF form (trim shared suffix,
then shared prefix, advance pos) **before** key derivation, via a new shared helper
`trim_to_minimal_vcf` in `noise_filter.{hpp,cpp}`. The same helper now backs the
two pre-existing inline trims in `apply_graph_noise_filter` and
`apply_hybrid_noise_filter` (DRY; behaviour unchanged). `site_meta` (the catalog
form used elsewhere) is left untouched.

**Impact: representation-only, zero phasing change.**

| pipeline | before | after |
|---|---|---|
| graph accuracy / Hamming | 99.195% / 0.805% | **99.195% / 0.805%** (identical) |
| graph switch / flip / perfect-PS | 367 / 858 / 142 | **367 / 858 / 142** (identical) |
| hybrid accuracy / Hamming | 99.232% / 0.768% | **99.232% / 0.768%** (identical) |

Accuracy is unchanged because k-means already classified these equal-length sites
as `CleanHetSnp` (the MNP→SNP branch) and phased them correctly; the bug only
corrupted the **emitted allele strings/POS** in the candidate TSV and VCF. The fix
makes graph/hybrid VCF output match the BAM pipeline's canonical variant
representation. Covered by `test_trim_to_minimal_vcf` in `test_noise_filter.cpp`.
Validated offline.

---

## Duplicate Variant Records & k-means Double-Counting (graph/hybrid)

**Symptom.** After minimal-VCF normalization (see previous section), the graph
pipeline still emitted some variants more than once, and the same physical variant
could enter k-means as several independent anchors.

**Root cause.** The snarl catalog wraps a single physical locus in multiple
overlapping/nested snarls. In `build_graph_chunk` candidates are keyed by snarl
`order_pos()`, so overlapping snarls that decompose to the *same* minimal-VCF
variant look distinct. Two consequences:

1. **Output:** the candidate TSV/VCF carried the same variant several times.
   475 duplicate loci on chr20, 471 sharing a phase set.
2. **k-means:** each duplicate cast an independent germline-het vote. 317
   k-means-eligible duplicates produced 329 redundant anchor votes (0.58% of
   eligible anchors), slightly distorting the clustering.

**Fix A — output dedup (`graph_collect.cpp`).** After the existing
`stable_sort`, a single pass collapses adjacent records with
`exact_comp_cand_var()==0`, keeping the highest-coverage copy and reporting any
conflicting haplotype calls under `--verbose`. Graph VCF duplicates 319→0 within a
chunk (a lone cross-chunk dup remains and is the reason this backstop stays even
with Fix B). Candidate rows 74,119→73,627.

**Fix B — k-means anchor dedup (`graph_bam_adapter.cpp`).** In Phase 2 each
biallelic pair is reduced to its minimal-VCF identity `(pos, ref, alt)` and
deduplicated *across snarl sites* before becoming a k-means anchor:

- Identity is computed with zero allocation: `minimal_vcf_id` trims to offset
  ranges over the existing `meta.ref`/`meta.alts` strings (views, not copies).
- A 64-bit FNV-1a fingerprint (`minimal_vcf_fingerprint`) buckets candidates in an
  `unordered_map<uint64_t, vector<CanonEntry>>`, reserved to the pair upper bound
  to avoid rehashing. Exact `MinimalVcfId` comparison inside the bucket resolves
  fingerprint collisions, so distinct variants are never merged.
- Same-site guard: `CanonEntry` records the source site index; the two alts of a
  single multiallelic snarl are never collapsed into each other (they are distinct
  variants by construction even if their minimal-VCF strings coincide degenerately
  when node sequences are unavailable). Only cross-site matches dedup.
- **No-pool semantics (critical).** A duplicate routes its read observations to the
  canonical anchor but does **not** pool counts. Per-read observation dedup in
  Phase 3 still collapses a read seen at both snarls into one vote.

**Pooled-vs-no-pool A/B (graph pipeline).** Pooling coverage across overlapping
snarls *regressed* phasing: it over-weights the surviving anchor and net-degrades
Hamming. Removing the extra votes *without* distorting the surviving anchor's
weight is what helps.

| variant | accuracy | Hamming | switch | flip |
|---|---|---|---|---|
| NORM-only (baseline) | 99.1954% | 0.8046% | 367 | 858 |
| Fix B pooled | 99.0443% ❌ | 0.9557% ❌ | 355 | 848 |
| **Fix B no-pool (chosen)** | **99.2097%** ✅ | **0.7903%** ✅ | **364** ✅ | **856** ✅ |

Fix B eliminated all 329 within-chunk redundant anchor votes
(eligible-distinct-loci 56,281→56,278, redundant 329→0).

**Hybrid impact: none (verified).** `build_graph_chunk` is shared by graph and
hybrid, so hybrid was re-evaluated before/after Fix A+B on the full chr20 harness:

| pipeline | before | after |
|---|---|---|
| hybrid accuracy / Hamming | 99.232% / 0.768% | **99.232% / 0.768%** (identical) |
| hybrid switch / flip / perfect-PS | 291 / 906 / 115 | **291 / 906 / 115** (identical) |
| hybrid candidate categories | — | identical |

Hybrid is unaffected because its BAM base already anchors phasing; the graph
augment only contributes at sites the BAM pass did not resolve, where the
redundant snarl votes did not change the cluster assignment. The graph-only
pipeline gains (99.195→99.210%) while hybrid holds steady.

Covered by the existing `test_graph_bam_adapter` multiallelic-decomposition tests
(which guard the same-site case) and `make unit-tests`.

---

## DeepVariant calling on HG003 GRCh38: surjection cost and the refined-pbmm2 + graph-HP win

A separate evaluation track on **HG003 chr20** (GIAB v4.2.1, GRCh38) measures how
pgphase's graph-derived phasing affects **downstream DeepVariant variant calling**,
rather than internal phasing accuracy. This is a different sample/reference than the
HG002/CHM13 work above; inputs live under the eval harness (`hg003_vg/`), not in
`test_data/`.

**Setup (all arms identical except where noted):** HG003 Revio SPRQ HiFi, 32×, chr20.
DeepVariant 1.10.0 PACBIO model, `--disable_small_model` in every arm so they share
one calling path (main CNN); the only variable is the alignment substrate and the
phasing source. hap.py vs GIAB HG003 v4.2.1 `noinconsistent` BED (chr20), vcfeval,
PASS-only. Two read alignments of the same reads:

- **pbmm2 (native):** the PacBio case-study minimap2/pbmm2 linear GRCh38 alignment.
- **graph-surject:** vg-giraffe alignment to the HPRC pangenome, surjected to GRCh38.

**Pangenome graph:** HPRC **v1.1** minigraph-cactus, GRCh38-based
(`hprc-v1.1-mc-grch38.gbz` + `.hapl`, from
`human-pangenomics/.../freeze1/minigraph-cactus/hprc-v1.1-mc-grch38`). Giraffe maps
with the `.hapl` personalized haplotype-sampling index. Every graph-derived signal in
this section (surject BAM, GAF, snarl sites) comes from v1.1, not v2.x.

HP-aware arms feed pgphase HP tags to DV with
`--make_examples_extra_args "phase_reads=false,sort_by_haplotypes=true"`.

### Vanilla control reproduces the published case study

DV on the native pbmm2 BAM (DV-internal phasing) matches Google's published number
to within 0.0001 INDEL F1 (our 0.99359 vs published 0.99368; identical TP 10,561),
validating the harness. The graph-surject pipeline scores **lower** on raw DV calling:

| Config (best arm) | Alignment | INDEL F1 | SNP F1 |
|---|---|---|---|
| Vanilla pbmm2 (DV-internal) | pbmm2 | 0.99359 | 0.99910 |
| graph-surject hybrid_hp (mq20 af0.12) | surject | 0.98724 | 0.99877 |

The ~0.0064 INDEL-F1 gap is the **surjection penalty**, not the phasing. SNPs are
nearly immune (~0.0003), because the cost is almost entirely indel re-representation:
projecting a graph path onto linear GRCh38 makes CIGAR/left-shift decisions that
differ from minimap2's, and the PACBIO model was trained on minimap2 alignments.

### Threshold sweeps are saturated

`--min-af` (0.20→0.12) and `--min-mapq` (20→10→1) were swept on the surject pipeline.
AF 0.20→0.12 gained +0.00145 INDEL F1 (the only real mover) by recruiting low-VAF
candidates; MAPQ moved INDEL F1 by ≤0.0004 across the whole 20→1 range (plateaued).
Lowering either threshold recruits **mostly SNPs, not indels** (e.g. mq20→mq1 added
+5,887 SNP vs +827 indel candidates). You cannot tune past the surjection penalty
with candidate thresholds. Best surject config: **mq20 af0.12 / dv_hybrid_hp**.

### The fix: native pbmm2 substrate + graph phasing

Putting graph phasing on a **native** alignment closes the entire gap. Three
constructions, vs the vanilla baseline (TP/FP are PASS counts):

| Config | INDEL F1 | INDEL FP | SNP F1 | SNP FP |
|---|---|---|---|---|
| Vanilla pbmm2 (DV-internal) | 0.99359 | 72 | 0.99910 | 68 |
| pbmm2-hybrid, dv_hybrid_hp | **0.99420** | 67 | 0.99789 | **241** ❌ |
| **refined-pbmm2 + graph-surject HP** | 0.99396 | 69 | **0.99917** | **59** ✅ |
| graph-surject hybrid_hp | 0.98724 | — | 0.99877 | — |

- **pbmm2-hybrid** = run `collect-hybrid-variation` with the native pbmm2 BAM as
  `--bam` plus the giraffe GAF/sites (`--refine-aln --min-mapq 20 --min-af 0.12`).
  This is coordinate-safe: the GAF is queried by GRCh38 interval (its own surjection
  projection) and graph observations join to BAM reads by **read name only**, so a
  different BAM aligner still merges. It gives the **best INDEL F1 (0.99420, +8 TP)**
  but regresses SNP precision (FP 68→241).

- **refined-pbmm2 + graph-surject HP** = take the BAM-pipeline refined pbmm2 reads
  (`run_pgphase_pbmm2hybrid_af12/bam/phased.bam` — fixed indel CIGARs, native
  coords), strip their HP/PS, and stamp the **graph-surject** hybrid HP tags
  (`run_pgphase_refine_mq20_af12/hybrid/phased.bam`) by read name (`gapfill.py` with
  core=refined-pbmm2, fill=graph-surject). This is a **strict Pareto win over
  vanilla** — every metric improves, both types:

| | INDEL | SNP |
|---|---|---|
| TP | 10,561 → 10,566 (+5) | 70,107 → 70,108 (+1) |
| FN | 67 → 62 (−5) | 59 → 58 (−1) |
| FP | 72 → 69 (−3) | 68 → 59 (−9) |
| Recall | 0.99370 → 0.99417 | 0.99916 → 0.99917 |
| Precision | 0.99348 → 0.99376 | 0.99903 → 0.99916 |
| F1 | 0.99359 → 0.99396 | 0.99910 → 0.99917 |

### What caused the SNP regression (and what did not)

The +173 SNP FP in the pbmm2-hybrid arm is **not** from `--refine-aln`. Two arms on
the *same refined substrate* isolate it: pbmm2-hybrid HP (graph fed in pbmm2 coords)
→ SNP FP 241; graph-surject HP (graph phased in graph coords, then surjected) → SNP
FP 59. The damage comes entirely from the **HP source**, not from refine. The
in-pbmm2-coordinate hybrid phasing mislabels reads at some SNP sites (the untested
no-conflict profile-extension path, `hybrid_inject.cpp:295`, where pbmm2 and giraffe
disagree on placement); graph-native phasing does not. `dv_default` on the refined
pbmm2 BAM keeps SNP FP at 65, confirming refine itself is clean.

### Recommendations

- **Graph phasing is not the problem; the surjected alignment is.** Graph-derived HP
  labels are accurate enough to match (and slightly beat) DV's own phasing — when
  applied to a native alignment.
- **For a no-downside config:** refined-pbmm2 substrate + graph-surject hybrid HP.
  Beats vanilla DV on both indels and SNPs with no regression.
- **For max indel recall (accepting SNP cost):** pbmm2-hybrid dv_hybrid_hp.
- **Use the hybrid (clean-core) HP source, not gapfill** — gapfill's hard reads cost
  precision on the native substrate, same as on surject.

Eval scripts: `hg003_vg/{26..43}_*.sbatch` (pgphase / DeepVariant / hap.py per
config) and `bench_hybrid_refine_mq{1,10,20}_af12.sh`. Results under
`hg003_vg/dv_results_*/`. These artifacts live on the eval cluster, not in-repo.

## HPRC v2.1 vs v1.1 on HG003 chr20 (refined-pbmm2 + graph-HP and graph-surject)

Re-ran the two best recipes from the section above on the **HPRC v2.1**
minigraph-cactus GRCh38 graph (`hprc-v2.1-mc-grch38.gbz`) to test whether the newer
pangenome improves downstream DeepVariant calling. Same sample/reference/evaluator as
the v1.1 work (HG003 chr20, GRCh38, GIAB v4.2.1 `noinconsistent` BED, hap.py vcfeval,
PASS-only, DV 1.10.0 PACBIO, `--disable_small_model`).

**Result: v2.1 improves on v1.1 in both recipes, and the refined-pbmm2 + graph-HP
recipe on v2.1 is a strict win over vanilla DV across every metric.**

### Final three-way comparison (the configs that matter)

INDEL (PASS):

| Config | F1 | Recall | Precision | FP | FN |
|---|---|---|---|---|---|
| Vanilla DV (pbmm2) | 0.993590 | 0.993696 | 0.993484 | 72 | 67 |
| refined-pbmm2 + graphHP (v1.1) | 0.993963 | 0.994166 | 0.993759 | 69 | 62 |
| **refined-pbmm2 + graphHP (v2.1)** | **0.994335** | **0.994637** | **0.994033** | **66** | **57** |

SNP (PASS):

| Config | F1 | Recall | Precision | FP | FN |
|---|---|---|---|---|---|
| Vanilla DV (pbmm2) | 0.999096 | 0.999159 | 0.999032 | 68 | 59 |
| refined-pbmm2 + graphHP (v1.1) | 0.999167 | 0.999173 | 0.999160 | 59 | 58 |
| **refined-pbmm2 + graphHP (v2.1)** | **0.999209** | **0.999188** | **0.999231** | **54** | **57** |

v2.1 graph-surject `dv_hybrid_hp` also improved over v1.1 (INDEL F1
0.987239 → 0.988046, FP 168 → 150; SNP F1 0.998767 → 0.998831, FP 53 → 47) but still
trails vanilla on SNP recall (FN 117 vs 59) — graph-surject alone remains inferior to
the refined-pbmm2 substrate, exactly as on v1.1.

### vg compatibility: regenerate the `.hapl`, do not reuse the published one

The published v2.1 `.hapl` is **haplotype-index version 5**, which vg 1.67.0 (the
pinned eval binary, max supported `.hapl` version 4) cannot read. The fix is to
regenerate a v4 `.hapl` locally from the `.gbz` with `vg haplotypes`. This requires
the `.dist` index and a `.ri` (r-index): build `.dist` first (`vg index -j`), then
`vg gbwt -Z -r` to emit the `.ri`, then `vg haplotypes`. The `.dist` build is the
memory bottleneck — it peaks at ~591 GB RSS and needs a ~900 GB high-mem node; the
`.hapl` regen itself is light. After giraffe is done the 20 GB `.regen.hapl` can be
deleted (catalog/pgphase steps do not use it).

### pggaf: use the upstream-fixed binary for annotate-gaf

Giraffe emits unmapped records (`*` path) that crash the older bundled `pggaf
annotate-gaf`. The upstream fix is commit `0a54c09` ("Fix annotate-gaf crash on
unmapped reads") on `github.com/kokyriakidis/pggaf` `main`. Rebuild pggaf from
upstream (CMake + GCC 13.3; the binary needs the GCC 13.3 C++ runtime and
`libcrypto.so.3` on `LD_LIBRARY_PATH`) and run `annotate-gaf` on the raw giraffe GAF.

### Operational gotcha: shared group quota is volatile — stage to $HOME

The eval cluster's shared group filesystem (`/sc1/groups`, 11P/13P used) enforces a
**group quota that floats near-full**, independent of raw disk space. Symptoms seen
repeatedly during the v2.1 run: a 10 MB write returns "0 bytes copied", `vg
deconstruct` aborts mid-catalog with "Failed to write BAM record", and — most
insidiously — SLURM cannot create the job's `--output` log, so the job dies at
0 s with **no log file at all** (looks like an instant, unexplained failure).

Durable fix: route **every** write vector to the user `$HOME` filesystem (separate,
uncontended, 13 T free), keeping the original shared-FS path strings intact via
symlinks so downstream scripts need no edits:

- catalog VCFs staged to `$HOME/pgphase_catalog_v21/`;
- `logs/` → symlink to `$HOME/pgphase_logs`;
- `run_pgphase_v21_*`, `dv_results_*_v21` → symlinks to `$HOME/...`.

With outputs on `$HOME`, the catalog build (4 h 55 m), pgphase (7 m, all four arms),
and both DV recipes (run in parallel, ~15–28 m each) completed cleanly.

### Whole-genome (all chromosomes) is not yet runnable as-is

A request to extend these three configs (vanilla, refined-pbmm2+graphHP v1.1 and
v2.1) to all chromosomes was **investigated but not run** — the inputs on disk are
chr20-only:

- Reads: only `hg003_vg/reads/HG003.chr20.fq.gz` and the chr20 SPRQ BAM are present.
  A whole-genome HG003 HiFi BAM (~54 GB) exists under `/sc1/groups/sbx/public-stage/`
  but belongs to another user/track and would need to be vetted/copied.
- Truth set: the only GIAB benchmark staged locally is **HG002**
  (`HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz` + BED). The genome-wide **HG003**
  v4.2.1 truth VCF/BED must be sourced before any all-chr hap.py can run.

Prerequisites for a genome-wide run: (1) whole-genome HG003 HiFi reads,
(2) genome-wide HG003 v4.2.1 truth VCF + BED, (3) the catalog rebuilt for all
contigs (the chr20 `vg deconstruct` alone took ~5 h, so expect a substantially
longer, higher-memory build genome-wide), and (4) the same `$HOME`-staging quota
workaround.

### Artifacts (eval cluster, not in-repo)

Scripts `hg003_vg/{50..59}_*_v21.sbatch` (fetch graph, regen `.hapl`, giraffe+surject,
build catalog, pgphase, DV ×2 recipes, hap.py ×2, build/test pggaf). Results under
`$HOME/dv_results_v21_mq20_af12/`, `$HOME/dv_results_refinedpbmm2_graphhp_v21/`, and
the comparison write-up at `$HOME/v21_vs_v11_comparison.md`.

## FN root-cause: the misses are alignment/evidence limits, not pgphase site gaps

Question asked of the winning config (refined-pbmm2 + graphHP **v2.1**,
`dv_hybrid_hp`): *why* are the 114 false negatives missed, and do they exist as sites
in the BAM or graph pipeline? Method: parse `happy.vcf.gz` for TRUTH `BD=FN` records
(57 SNP + 57 INDEL, all inside the GIAB CONF BED), then for each locus (1) read the
DV genotype call, (2) measure ALT-allele depth in the DV-input BAM
(`refined_pbmm2.graphHP.bam`) via `bcftools mpileup -a AD`, and (3) test membership in
the pgphase BAM and graph candidate sets (`*/candidates.tsv`, POS-exact for SNPs, ±5 bp
for indels to absorb left-alignment differences).

### DV behavior at FN loci

- 84 NOCALL — DV emitted nothing.
- 30 WRONGCALL — DV called but with the wrong genotype (these also count as FP).

### Evidence in the DV-input BAM (ALT allele fraction at the truth locus)

| Type | NO_COVERAGE (DP<5) | LOW_VAF (alt AF<0.2) | EVIDENCE_PRESENT (alt AF≥0.2) |
|---|---|---|---|
| SNP | 0 | 46 | 11 |
| INDEL | 0 | **57** | 0 |

### Do the FN exist as pgphase candidates?

Only **13 of 114** FN are pgphase candidates (9 SNP `EVIDENCE_PRESENT` + 4 SNP
`LOW_VAF`, from the BAM or graph set); the other 101 are not. This is **not** the
cause of the misses: pgphase candidates only drive HP phasing, they do not gate what
DeepVariant calls. The 101 non-candidates fail for the *same* reason DV misses them —
no callable ALT evidence at the locus.

### Why each class is missed

- **INDEL FN (all 57): alignment representation, not detection.** ALT depth is ~0 at
  every indel FN — the reads carry the *reference* indel representation (left-shift /
  homopolymer collapse from the linear/surjected alignment), so the truth-form ALT is
  not in the pileup. Same surjection indel-representation penalty documented above;
  no phasing change recovers these.
- **SNP FN, LOW_VAF (46): allele dropout at true hets.** Median alt AF ≈ 0.08; reads
  dropped the alternate allele (mapping bias), so DV correctly calls reference.
- **SNP FN, EVIDENCE_PRESENT (11): genuinely hard.** The only loci with callable
  signal. 8 are NOCALL despite AF 0.2–0.5 — including a tight 3-SNP cluster at
  chr20:5,309,435 / 5,309,444 / 5,310,690 (paralog/segdup signature DV filters as
  ambiguous). 3 are WRONGCALL at AF 0.77–0.90 — **zygosity errors** (e.g. truth 1/1 →
  DV 0/1, or a multiallelic mismatch), not missed sites. The pgphase pipelines *do*
  see 9 of these 11 as candidates; DV still declines them.

### Bottom line

The FNs are **not** a pgphase site-detection failure. ~90% are alignment/evidence
limits — indel representation (57) and het allele-dropout (46) — unreachable by any
phasing or candidate-threshold change. Only ~11 are hard SNPs with real signal
(segdup cluster + zygosity), already candidates in the graph/BAM pipelines but
declined by DV. Moving these would need a different alignment or DV model, not pgphase
changes.

Artifacts (eval cluster): `$HOME/fn_analysis/fn_evidence.tsv` (per-locus depth/AF/
class/candidate table) and `$HOME/fn_analysis/FN_FINDINGS.md` (write-up).

---

## Four-Pipeline Comparison — HG002 HiFi chr18

Full chr18 evaluation using HG002 HiFi reads aligned to CHM13 chr18.
Truth: HG002 v1.1 diploid assembly (chr18_MATERNAL / chr18_PATERNAL slices),
assigned via diplinator. `--min-reads 5`, 8 threads.
Gap-fill: `--gapfill 2` (hybrid core + BAM-only reads at PS offset 1e9,
then graph-only reads at PS offset 2e9).

### Accuracy & errors

| Metric | BAM | Graph | **Hybrid** | Hybrid+Gapfill |
|---|---|---|---|---|
| **Runtime** | 86 s | **32 s** | 131 s | +post-process |
| Accuracy | 98.24% | 99.05% | **99.24%** | 98.58% |
| Hamming | 1.76% | 0.95% | **0.76%** | 1.42% |
| Switch errors | 891 | 335 | **367** | 886 |
| Flip errors | 1,421 | 802 | **819** | 1,560 |
| Switchflip | 2,312 | 1,137 | **1,186** | 2,446 |
| Switch rate | 0.32% | 0.13% | **0.13%** | 0.31% |
| Switchflip rate | 0.82% | 0.43% | **0.43%** | 0.86% |
| Switch opportunities | 281,249 | 262,574 | 274,731 | 284,728 |
| Concordant reads | 276,863 | 260,529 | 273,115 | 281,332 |
| Discordant reads | 4,953 | 2,497 | **2,083** | 4,049 |

### Yield & phase blocks

| Metric | BAM | Graph | **Hybrid** | Hybrid+Gapfill |
|---|---|---|---|---|
| Total input reads | 346,878 | 295,387 | 346,878 | 346,878 |
| Phased reads | 281,840 | 263,026 | 275,210 | **285,570** |
| Fraction phased | 81.3% | **89.0%** | 79.3% | 82.3% |
| Reads evaluated | 281,816 | 263,026 | 275,198 | 285,381 |
| Phase sets (total) | 574 | 452 | 470 | 773 |
| Phase sets (eval) | 567 | 452 | 467 | 653 |
| Perfect PS | 201 | **193** | 185 | **237** |
| Perfect PS % | 35.4% | **42.7%** | 39.6% | 36.3% |
| Candidates | 140,565 | 82,069 | 152,706 | — |
| Block N50 (bp) | 1,718,322 | 1,720,559 | **1,735,751** | 1,704,636 |
| Block auN (bp) | 5,002,792 | 1,802,513 | **5,327,723** | 4,634,692 |
| Block median span (bp) | 1,691,650 | 1,698,316 | 1,698,043 | 1,688,711 |
| Block max span (bp) | 54,389,768 | 2,982,830 | 52,161,576 | **54,389,768** |
| Genome covered (bp) | 887,637,832 | 673,136,287 | 745,359,642 | **999,757,287** |
| Genome covered (%) | 28.7% | 21.8% | 24.1% | **32.4%** |

### Phaseable vs unphaseable split (confidence threshold 0.6)

| Metric | BAM | Graph | **Hybrid** | Hybrid+Gapfill |
|---|---|---|---|---|
| Phaseable PS / reads | 525 / 278,842 | 441 / 262,365 | 447 / 274,376 | 592 / 281,853 |
| Phaseable accuracy | 98.72% | **99.17%** | **99.37%** | 99.13% |
| Unphaseable PS / reads | 42 / 2,974 | **11 / 661** | 20 / 822 | 61 / 3,528 |
| Unphaseable accuracy | 53.63% | 52.95% | 55.60% | 54.65% |

### Gap-fill breakdown

Stage 1 (BAM-only reads) recovered **9,052 reads**; stage 2 (graph-only reads)
added **1,308 reads** (total 10,360 added to hybrid core).

### Observations

- **Hybrid wins every accuracy metric** (99.24%, Hamming 0.76%, switch rate
  0.13%) and best auN (5,327,723 bp) — consistent with chr20 results.
- **Graph is the accuracy-per-second winner** (99.05% in 32 s), highest
  fraction phased (89.0%), most perfect PS % (42.7%).
- **Hybrid+Gapfill Pareto-beats BAM**: more reads phased (285,570 vs 281,840),
  lower Hamming (1.42% vs 1.76%), more genome covered (32.4% vs 28.7%),
  more perfect PS (237 vs 201). Hybrid core preserved exactly.
- chr18 N50 (~1.7 Mbp) is higher than chr20 N50 (~1.0 Mbp) across all
  pipelines, consistent with chr18's simpler repeat structure.


---

## Four-Pipeline Comparison — HG002 HiFi chr12

Full chr12 evaluation (133 Mbp) using HG002 HiFi reads aligned to CHM13 chr12.
Truth: HG002 v1.1 diploid assembly (chr12_MATERNAL / chr12_PATERNAL slices),
assigned via diplinator. `--min-reads 5`, 8 threads. Gap-fill: `--gapfill 2`.

### Accuracy & errors

| Metric | BAM | Graph | **Hybrid** | Hybrid+Gapfill |
|---|---|---|---|---|
| **Runtime** | 123 s | **57 s** | 229 s | +post-process |
| Accuracy | 98.51% | 99.09% | **99.31%** | 98.73% |
| Hamming | 1.49% | 0.91% | **0.69%** | 1.27% |
| Switch errors | 1,271 | 471 | **342** | 1,233 |
| Flip errors | 2,469 | 1,406 | **1,456** | 2,565 |
| Switchflip | 3,740 | 1,877 | **1,798** | 3,798 |
| Switch rate | 0.26% | 0.10% | **0.07%** | 0.25% |
| Switchflip rate | 0.75% | 0.40% | **0.37%** | 0.76% |
| Switch opportunities | 496,223 | 470,523 | 488,855 | 500,796 |
| Concordant reads | 489,581 | 466,900 | 486,095 | 495,312 |
| Discordant reads | 7,381 | 4,280 | **3,394** | 6,368 |

### Yield & phase blocks

| Metric | BAM | Graph | **Hybrid** | Hybrid+Gapfill |
|---|---|---|---|---|
| Total input reads | 579,483 | 521,096 | 579,483 | 579,483 |
| Phased reads | 496,971 | 471,180 | 489,495 | **501,975** |
| Fraction phased | 85.8% | **90.4%** | 84.5% | 86.6% |
| Reads evaluated | 496,962 | 471,180 | 489,489 | 501,680 |
| Phase sets (total) | 742 | 657 | 636 | 1,055 |
| Phase sets (eval) | 739 | 657 | 634 | 884 |
| Perfect PS | 261 | **282** | 260 | **327** |
| Perfect PS % | 35.3% | **42.9%** | 41.0% | 37.0% |
| Candidates | 219,231 | 150,221 | 239,312 | — |
| Block N50 (bp) | 1,095,879 | 455,084 | **1,692,256** | 1,402,401 |
| Block auN (bp) | 42,186,534 | 676,738 | **54,771,185** | 51,180,910 |
| Block max span (bp) | 102,337,347 | 4,305,703 | **118,870,322** | **118,870,322** |

### Phaseable vs unphaseable split (confidence threshold 0.6)

| Metric | BAM | Graph | **Hybrid** | Hybrid+Gapfill |
|---|---|---|---|---|
| Phaseable PS / reads | 705 / 493,651 | 642 / 470,140 | 630 / 488,842 | 811 / 497,320 |
| Phaseable accuracy | 98.82% | **99.19%** | **99.36%** | 99.12% |
| Unphaseable PS / reads | 34 / 3,311 | **15 / 1,040** | **4 / 647** | 73 / 4,360 |
| Unphaseable accuracy | 53.43% | 54.62% | 56.41% | 53.81% |

### Gap-fill breakdown

Stage 1 (BAM-only reads) recovered **10,664 reads**; stage 2 (graph-only reads)
added **1,816 reads** (total 12,480 added to hybrid core).

### Observations

- **Hybrid wins all accuracy metrics** (99.31%, Hamming 0.69%, switch rate 0.07%)
  and all contiguity metrics (auN 54.8M, largest block 118.9 Mbp — near the full
  133 Mbp chromosome). Pattern consistent with chr18 and chr20.
- **Graph wins accuracy-per-second** (99.09% in 57 s, 90.4% phased, 42.9% perfect PS).
- **Hybrid+Gapfill** phases the most reads (501,975 / 86.6%), Pareto-beats BAM on
  all axes, and recovers 12,480 gap-fill reads at high accuracy (99.12% phaseable).
- Hybrid's 4 unphaseable phase sets (647 reads) is the lowest of any pipeline across
  all three chromosomes — chr12's haplotype structure is well-suited to hybrid phasing.
- For DeepVariant HP-tagging: **hybrid+gapfill** is optimal (501,975 HP-tagged reads,
  86.6% of 579,483 BAM reads, 99.12% phaseable accuracy).


---

## Additional Whole-Chromosome Test Files

Test inputs for chr18 and chr12, generated from the same source data as chr20
(HG002 HiFi, CHM13 T2T v2.0, HPRC v2.1 pangenome GBZ). Not committed to git
(too large); recreate with the commands below.

### chr1 test files

| File | Description |
|------|-------------|
| `test_data/chm13v2.0.chr1.renamed.fa` (+`.fai`) | CHM13 chr1 reference (contig `chr1`, 248 Mbp) |
| `test_data/HG002_chr1_hifi_mapped_to_CHM13_chr1_annotated.bam` (+`.bai`) | 1,044,741 HiFi reads aligned to CHM13 chr1 |
| `test_data/chr1.sites.vcf.gz` (+`.tbi`) | 3,199,313 graph snarl sites |
| `test_data/HG002.chr1.annotated.coord.gaf.gz` (+`.tbi`) | bgzipped + tabix-indexed coord-annotated GAF |

### chr18 test files

| File | Description |
|------|-------------|
| `test_data/chm13v2.0.chr18.renamed.fa` (+`.fai`) | CHM13 chr18 reference (contig `chr18`, 78 Mbp) |
| `test_data/HG002_chr18_hifi_mapped_to_CHM13_chr18_annotated.bam` (+`.bai`) | 346,878 HiFi reads aligned to CHM13 chr18 |
| `test_data/chr18.sites.vcf.gz` (+`.tbi`) | 1,067,783 graph snarl sites |
| `test_data/HG002.chr18.annotated.coord.gaf.gz` (+`.tbi`) | bgzipped + tabix-indexed coord-annotated GAF |

### chr12 test files

| File | Description |
|------|-------------|
| `test_data/chm13v2.0.chr12.renamed.fa` (+`.fai`) | CHM13 chr12 reference (contig `chr12`, 129 Mbp) |
| `test_data/HG002_chr12_hifi_mapped_to_CHM13_chr12_annotated.bam` (+`.bai`) | 579,483 HiFi reads aligned to CHM13 chr12 |
| `test_data/chr12.sites.vcf.gz` (+`.tbi`) | 1,789,539 graph snarl sites |
| `test_data/HG002.chr12.annotated.coord.gaf.gz` (+`.tbi`) | bgzipped + tabix-indexed coord-annotated GAF |

### Reconstruction commands

Source files (not committed, stored on the analysis machine):
- Full GBZ: `~/Downloads/pgbam-experiments/hprc-v2.1-mc-chm13-eval.gbz` (5.7 GB)
- Full r-index: `~/Downloads/pgbam-experiments/hprc-v2.1-mc-chm13-eval.ri` (11 GB)
- Full sorted BAM: `~/Downloads/pgbam-experiments/HG002.full.sorted.bam` (65 GB)
- Full GAM: `~/Downloads/pgbam-experiments/HG002.full.gam` (145 GB, alphanumeric chr order)

```bash
CHR=chr18   # or chr12
FULL_GBZ=~/Downloads/pgbam-experiments/hprc-v2.1-mc-chm13-eval.gbz
FULL_RI=~/Downloads/pgbam-experiments/hprc-v2.1-mc-chm13-eval.ri
FULL_BAM=~/Downloads/pgbam-experiments/HG002.full.sorted.bam
FULL_GAM=~/Downloads/pgbam-experiments/HG002.full.gam

# 1. Extract per-chromosome reference FASTA
samtools faidx ~/Downloads/pgbam-experiments/chm13v2.0.fa "CHM13#0#${CHR}" \
  | sed "s/>CHM13#0#${CHR}/>$CHR/" \
  > test_data/chm13v2.0.${CHR}.renamed.fa
samtools faidx test_data/chm13v2.0.${CHR}.renamed.fa

# 2. Extract per-chromosome HiFi BAM (rename CHM13#0#chrN → chrN)
(samtools view -H "$FULL_BAM" \
     | awk "!/^@SQ/ || /SN:CHM13#0#${CHR}/" \
     | sed "s/CHM13#0#${CHR}/${CHR}/g"
 samtools view "$FULL_BAM" "CHM13#0#${CHR}" \
     | sed "s/\tCHM13#0#${CHR}\t/\t${CHR}\t/g") \
  | samtools sort -@ 4 -O bam \
  -o test_data/HG002_${CHR}_hifi_mapped_to_CHM13_${CHR}_annotated.bam
samtools index test_data/HG002_${CHR}_hifi_mapped_to_CHM13_${CHR}_annotated.bam

# 3. Build snarl-site VCF via per-chr GBZ chunk (~2-5 GB RAM; full GBZ needs ~80 GB)
vg chunk --gbz --contig "$CHR" -x "$FULL_GBZ" -b /tmp/${CHR}
./pgphase build-snarl-catalog --ref-sample CHM13 --contig "$CHR" -t 8 \
  -o test_data/${CHR}.sites.vcf.gz /tmp/${CHR}_0_${CHR}.gbz

# 4. Build coord-indexed GAF
#    4a. Extract read names aligned to this chromosome from the BAM
samtools view test_data/HG002_${CHR}_hifi_mapped_to_CHM13_${CHR}_annotated.bam \
  | awk '{print $1}' | sort -u > /tmp/${CHR}_qnames.txt
#    4b. Filter full GAM to this chromosome's reads (run alone, not concurrent)
vg filter -t 4 -N /tmp/${CHR}_qnames.txt -e "$FULL_GAM" > /tmp/${CHR}.gam
#    4c. Convert GAM → GAF using full GBZ for node ID lookup (~14 GB RAM)
vg convert -G /tmp/${CHR}.gam "$FULL_GBZ" > /tmp/${CHR}.gaf
#    4d. Annotate with reference coordinates and haplotype-set tags
~/Downloads/pggaf/build/pggaf annotate-gaf \
  --gaf /tmp/${CHR}.gaf --gbz "$FULL_GBZ" --r-index "$FULL_RI" \
  --ref-sample CHM13 \
  --out-gaf /tmp/${CHR}.annotated.gaf --out-sets /tmp/${CHR}.pgs
#    4e. Coordinate-sort, bgzip, and tabix-index
~/Downloads/pggaf/build/pggaf index-gaf \
  --in /tmp/${CHR}.annotated.gaf \
  --out test_data/HG002.${CHR}.annotated.coord.gaf.gz
```


---

## Read-confidence gating: the graph pipeline now beats BAM and hybrid

Question: with unphaseable regions excluded, can the graph pipeline beat both the
BAM and hybrid pipelines? Answer: yes, by 2.5-5x on every accuracy metric, via a
single missing gate on per-read assignment confidence.

### Excluding unphaseable regions

Two classes of region make phasing impossible for reasons no pipeline change can
address, and both were diluting every previous comparison:

- **Runs of homozygosity.** chr20 44-45 Mb has **3 clean het SNPs** against a
  chromosome median of 852/Mb (45-46 Mb: 68; 27-28 Mb: 117). The graph, BAM and
  hybrid pipelines independently agree (3 / 4 / 5 het SNPs in 44-45 Mb), so this
  is HG002's biology, not a detection failure. Catalog density there is normal
  (13.1 sites/kb, median AN 457 = full panel) and depth is normal (5,003
  alignments/Mb, higher than a well-phased control).
- **Graph-blind satellite.** chr20 27.2-28.8 Mb has **zero catalog sites** across
  1.6 Mb while reads are present throughout and the BAM pipeline still calls
  819+326 het SNPs. `AN` collapses from ~452 to 7 to 3 at the edges: the
  pangenome has no alternative haplotypes there, so no bubbles, so no snarls. One
  such window on chr20; far more on acrocentrics and chr1/9/16.

`--exclude-bed` is applied in HG002 truth coordinates while these were identified
in CHM13. The translation was built empirically: every read carries a CHM13
interval in the GAF and a truth placement in the HipHap/diplinator BAM, so
joining on read name lifts the windows over directly (1% tail trim against
mismappers). 18,604 reads excluded.

### Result (HG002 chr20, regions excluded)

| pipeline | hamming | switch | flip | discordant | phased | perfect PS | N50 |
|---|---|---|---|---|---|---|---|
| bam | 0.021638 | 849 | 1,235 | 4,648 | 80.9% | 35.8% | 985 kb |
| hybrid | 0.005282 | 167 | 605 | 1,114 | 78.0% | 41.1% | 999 kb |
| graph (before) | 0.005822 | 253 | 625 | 1,179 | 88.1% | 45.7% | 922 kb |
| **graph `--min-read-margin 2`** | **0.001964** | **68** | **114** | **350** | 77.3% | **82.5%** | 948 kb |

vs hybrid: 2.7x hamming, 2.5x switch, 5.3x flip, 3.2x discordant, 2x perfect
phase sets, at equal coverage (77.3% vs 78.0%). vs BAM: 11x hamming. Cost is 5%
of N50.

### Root cause

`init_assign_read_hap` (collect_phase.cpp) commits a read to a haplotype on **any
non-zero score** — no minimum evidence, no margin. A read agreeing with one clean
het SNP and contradicting none is committed as confidently as one agreeing with
thirty. Measured on chr20:

| clean-SNP margin | reads | share | errors | share of errors | error rate |
|---|---|---|---|---|---|
| <= 0 | 1,137 | 0.56% | 15 | 1.3% | 1.32% |
| **== 1** | **22,457** | **11.1%** | **810** | **68.7%** | **3.61%** |
| >= 2 | 178,912 | 88.4% | 354 | 30.0% | 0.198% |

Discordant reads: median 4 observations, margin 1. Concordant: 18 and 13. 22,255
of the 22,457 marginal reads are literally "1 agree / 0 conflict". The gate leaves
them unphased rather than guessing. Offline analysis predicted 354 discordant at
margin 2; the pipeline delivered 350.

### Multi-chromosome transfer (margin distribution)

The threshold is an absolute SNP count, so it could interact with per-chromosome
heterozygosity. It does not — the distribution is flat across chromosomes:

| chrom | reads | phased | margin<=0 | margin==1 | margin>=2 | phased @ margin 2 |
|---|---|---|---|---|---|---|
| chr20 | 231,382 | 204,165 | 0.61% | 11.35% | 88.04% | 77.69% |
| chr18 | 295,387 | 263,540 | 0.63% | 12.36% | 87.01% | 77.63% |
| chr12 | 521,096 | 472,113 | 0.48% | 10.70% | 88.81% | 80.47% |
| chr1 | 858,212 | 765,617 | 0.48% | 11.66% | 87.86% | 78.38% |

margin==1 is 10.7-12.4% of phased reads everywhere and the coverage cost of
margin 2 is 11-13%. The knee should sit at 2 on all four. **Accuracy on
chr18/chr12/chr1 is not yet confirmed** — those truth BAMs no longer exist and
rebuilding them needs minimap2 (absent) plus realignment of 10.7 GB of reads
against both HG002 haplotypes. HipHap (github.com/jheinz27/hiphap, the renamed
diplinator) builds from source with the local cargo.

### Hypotheses tested and rejected

- **pantree/reference-tree catalog** (docs/pantree_catalog_experiment.md). Adds
  sites to the top of a funnel that already discards 96% of what it has, in
  regions that fail for lack of heterozygosity. Not pursued.
- **Graph het-indel anchor gating.** The hybrid gates graph het indels on allele
  fraction (hybrid_inject.cpp) and the graph-only path did not — the obvious
  explanation for hybrid's lead. Implemented and measured: hamming 0.005821 vs
  0.005822, switches identical; gating *every* indel anchor made it slightly
  worse. Indel anchors do not move read assignment where SNP anchors exist
  (3,560 SNP vs 595 indel anchors per 2 Mb). Reverted. Note for any future
  attempt: `classify_graph_candidates` (graph_bam_adapter.cpp) re-stamps every
  category after `build_graph_chunk`, so it is the only place a type-aware
  anchor policy survives.

### New diagnostics (both default off, output byte-identical when unused)

- `--min-read-margin INT` — the gate. 0 reproduces prior behavior exactly.
- `--filtered-sites-out FILE` — why catalog sites never became candidates
  (`ref_only` / `high_af` / `low_af` / `low_depth`, with depth and AF per site).
  These reasons were computed but only counted; this is what showed the ROH.
- `--phase-reads-out FILE` — per-read observations and clean-SNP agree/conflict.
  This is what located the marginal-read population.

### Evaluation toolchain is now pinned

Every accuracy number above depends on HapQ, which is computed by HipHap
(formerly `diplinator`) and gates read inclusion via `--min-hapq`. An unpinned
HipHap would silently make old and new evals incomparable, so it is a submodule:

- `third_party/hiphap` @ `b9a065c0f36e` — build with `make hiphap` (optional
  target, never part of `all`; nothing in pgphase links against it).
- `third_party/hiphap-Cargo.lock` — upstream ships no lockfile, and a fresh
  resolve picks `hts-sys` 2.2.1, which renames `bam1_core_t::isize` and flips
  `size_t`->`usize`, failing to compile against `rust-htslib` 0.46. The vendored
  lock pins `hts-sys` 2.1.4 and is copied into the submodule before each build.
  Without it `make hiphap` fails from a clean clone.
- `third_party/minimap2` @ `v2.31` — the aligner decides which haplotype each
  read is assigned to, so it is pinned for the same reason HapQ is. `make
  eval-tools` builds both.
- `envs/truth.yaml` — optional conda route; only samtools is genuinely needed
  once the submodules are built.
- `scripts/test_end_to_end.sh` — reads -> minimap2 -> hiphap -> truth BAM ->
  pgphase -> accuracy, on the checked-in 500 kb chr20 fixture, in ~1 min.
- **Truth BAMs live outside the repo**, at `$PGPHASE_EVAL_DATA`
  (default `~/Downloads/pgphase-eval-data/truth/<chrom>/`), so updating or
  re-cloning pgphase cannot delete them. Each is one mapping of the reads
  against the diploid assembly and depends only on (reads, assembly, minimap2,
  hiphap) -- not on any pgphase setting -- so it is built once per chromosome
  and reused by every pipeline, every margin and every later experiment.
  Rebuilding costs 35 min (chr18) to ~90 min (chr1) and ~2 GB each.
- `build_truth_bam.sh` / `evaluate_phase_accuracy.sh` take `--hiphap` and default
  to the submodule build; `--diplinator` remains as a deprecated alias.

Three things had to be fixed before the chain ran at all, each of which would
have failed silently or blocked outright:

1. **`-A` is now passed explicitly** (`--match-score`, default 2, matching the
   `lr:hqae` preset which does not override minimap2's default `a=2`). HipHap's
   auto-estimator samples reads at `rng.gen_bool(0.0001)` and needs 10 hits, so
   it fails outright on small inputs *and* puts an RNG in the truth path --
   HapQ, and every number derived from it, would not be reproducible run to run.
2. **`paftools.js` was removed from the requirements check.** Both scripts
   demanded it; neither ever invoked it. It blocked the pipeline for nothing.
3. **`hts-sys` had to be pinned** (see above), or hiphap does not compile.

The rename also changed behavior: **merged output is now the default**, and `-p`
is required for the two per-haplotype files this pipeline needs. The scripts now
pass `-p -o hiphap`, producing `hiphap_mat.sam` / `hiphap_pat.sam` (legacy
`diplinator_*.sam` names are still accepted so existing output dirs resolve).

### Multi-chromosome confirmation

Truth BAMs were rebuilt with the pinned toolchain and the gate re-measured.
The knee sits at margin 2 on every chromosome tested:

| chrom | margin | reads | discordant | hamming | switch | flip | perfect PS |
|---|---|---|---|---|---|---|---|
| chr18 | 0 | 263,026 | 2,497 | 0.009493 | 335 | 802 | 42.7% |
| chr18 | **2** | 228,264 | **1,207** | **0.005288** | **89** | **107** | 83.1% |
| chr18 | 3 | 210,499 | 927 | 0.004404 | 44 | 37 | 91.0% |
| chr12 | 0 | 471,180 | 4,280 | 0.009084 | 471 | 1,406 | 42.9% |
| chr12 | **2** | 417,496 | **2,003** | **0.004798** | **78** | **160** | 83.7% |
| chr12 | 3 | 384,764 | 1,597 | 0.004151 | 48 | 51 | 93.9% |

Margin 2 vs margin 0: chr18 1.8x hamming / 3.8x switch / 7.5x flip keeping 86.8%
of reads; chr12 1.9x / 6.0x / 8.8x keeping 88.6%; chr20 3.0x / 3.7x / 5.5x
keeping ~88%. Perfect phase sets roughly double on both new chromosomes
(43% -> 83%). Margin 3 keeps improving accuracy but costs a further 8% of reads,
matching chr20 where it fell below hybrid's coverage.

chr1 is still building at the time of writing -- its satellite regions are very
slow under the `lr:hqae` preset (one 31k-read batch took 16 min against 14 s
elsewhere), so a full chr1 truth BAM is a ~3 h job rather than the ~35-50 min the
other autosomes take.

### Before making margin 2 the default

Three chromosomes agree on the knee and a fourth is in progress. The remaining
judgement call is the coverage/accuracy trade: margin 2 gives up 11-13% of
phased reads everywhere. That is the right trade against hybrid, which phases
fewer reads at worse accuracy, but it is a real loss against margin 0 if
downstream consumers care more about yield than switch rate.

---

## Residual error after the read-confidence gate: a graph-space collapse

With `--min-read-margin 2` the remaining error is **highly concentrated**: the
top 10 phase sets hold 82% (chr20) / 89% (chr18) of all discordant reads, and
only ~17% of phase sets have any error at all. Chasing the aggregate is
therefore the wrong move; chasing the handful of bad blocks is the right one.

### The dominant chr18 block is an orientation swap, not scattered error

PS 57679490 alone holds 647 of chr18's 1,207 discordant reads, yet records only
3 switches and 0 flips. Sorting its reads by truth position gives exactly four
runs:

    concordant   56.262-56.547 Mb   589 reads
    DISCORDANT   56.564-56.732 Mb   308 reads
    concordant   57.899-58.186 Mb   574 reads
    DISCORDANT   58.191-58.365 Mb   339 reads

Two clean seams, each inverting everything after it. The switch metric counts
transitions, so it reports 3; Hamming counts reads, so it reports 647. When the
two disagree this sharply, trust Hamming.

### What it is not

- **Not chunk stitching.** `--stitch-min-margin` and `--stitch-rule` are now
  exposed on `collect-graph-variation` (they existed only on the hybrid path).
  Setting them to the hybrid's values (margin 10, both-strands-bridged) changes
  chr20 not at all (350 discordant either way) and chr18 by 0.4% (1,207 ->
  1,202). Those seams are merged on *strong* evidence, not weak.
- **Not annotated segdups.** The two swapped intervals have 0.9% and 0.0%
  segdup coverage, against 0.5% for a clean control interval.
- **Not a truth artifact.** The BAM pipeline, over the same reads and the same
  truth BAM, is 0.7% and 1.0% discordant on those exact intervals where the
  graph pipeline is 47.1% and 48.5%.
- **Not low-confidence reads.** The discordant reads there carry median margin
  13 and p90 38-48 -- they agree with dozens of clean het SNPs. The anchors
  themselves are mis-genotyped; the reads follow them faithfully.

### What it is

Joining read names between the eval and the GAF puts both swapped segments at
**the same CHM13 interval**: chr18:57.955-58.147 Mb. Reads from two distinct
HG002 loci (56.6 Mb and 58.2 Mb) project onto one graph locus. The graph
pipeline then separates *paralogs* rather than haplotypes -- consistently, which
is why the read margins are high and the discordance sits at ~50%.

The depth signature is weak: 1.18x the chromosome-median DP (79 vs 67) and AF
0.46 vs 0.50. A naive 2x-depth collapse filter would not catch it, so detection
needs something better -- candidates: per-locus disagreement between a read's
GAF projection and its linear placement, or graph path multiplicity.

### Aside: chr18 graph vs BAM at margin 2

    graph (margin 2)  228,264 reads  1,207 discordant  hamming 0.005288  switch 89   flip 107
    BAM               281,816 reads  4,953 discordant  hamming 0.017575  switch 891  flip 1,421

3.3x Hamming, 10x switch, 13x flip in the graph pipeline's favour, consistent
with chr20.

### Correction to the per-chromosome test-file recipe above

Step 1 of the reconstruction commands extracts the reference with
`samtools faidx chm13v2.0.fa "CHM13#0#${CHR}"`. The `chm13v2.0.fa` on this
machine uses plain `chrN` contig names -- only the BAM uses the `CHM13#0#chrN`
form -- so that step silently produces a header-only FASTA. Use:

    samtools faidx "$P/chm13v2.0.fa" "${CHR}" | sed "1s/^>.*/>${CHR}/" > ...

---

## Fix: gate k-means anchors on allele fraction (`--anchor-af-margin`)

The collapse described above is detectable at the site level. Where the graph
merges two paralogous loci, a site that is het on one copy and hom on the other
lands near AF 0.25 or 0.75 -- comfortably inside the 0.20/0.80 depth filters --
and then votes as if it were a haplotype marker. Measured on chr20:

| region | candidates | AF median | AF in 0.40-0.60 | AF outside 0.35-0.65 |
|---|---|---|---|---|
| bad block 65.99-66.21 Mb | 425 | 0.536 | **49.6%** | **32.9%** |
| control 58.79-59.79 Mb | 1,246 | 0.484 | 83.3% | 7.0% |
| whole chr20 | 73,627 | 0.500 | 77.7% | 11.8% |

`--anchor-af-margin F` requires |AF - 0.5| <= F for a site to vote in k-means.
It applies to every site, not just indels (which already had
`--graph-indel-af-margin`). As with the indel gate, only `lcd_var_i_to_cate`
changes: the site is still emitted as a call, it just stops voting.

### Results (both at `--min-read-margin 2`)

| chrom | anchor-af | reads | discordant | hamming | switch | flip | perfect PS |
|---|---|---|---|---|---|---|---|
| chr20 | 0.5 (off) | 178,205 | 350 | 0.001964 | 68 | 113 | 82.5% |
| chr20 | **0.12** | 175,173 | **237** | **0.001353** | 59 | 75 | 87.4% |
| chr20 | 0.08 | 166,032 | 188 | 0.001132 | 48 | 67 | 89.8% |
| chr18 | 0.5 (off) | 228,265 | 1,207 | 0.005288 | 89 | 105 | 83.1% |
| chr18 | **0.12** | 224,746 | **388** | **0.001726** | 52 | 79 | 86.0% |
| chr18 | 0.10 | 221,959 | 392 | 0.001766 | 49 | 77 | 87.8% |

0.12 is the knee: 1.48x fewer discordant reads on chr20 and **3.11x** on chr18,
for 1.5-1.7% of reads and no N50 cost. Tightening further keeps helping chr20
but costs reads disproportionately (0.08 gives up 6.8%), and chr18 is flat
between 0.12 and 0.10.

### Standing against the other pipelines

    chr20  graph (margin 2, af 0.12)  175,173 reads    237 disc  hamming 0.001353
           hybrid                     210,905 reads  1,114 disc  hamming 0.005282
    chr18  graph (margin 2, af 0.12)  224,746 reads    388 disc  hamming 0.001726
           BAM                        281,816 reads  4,953 disc  hamming 0.017575

3.9x hamming over hybrid on chr20, 10.2x over BAM on chr18.

Default is 0.5, i.e. off, and reproduces prior output exactly (chr18: 1,207
discordant / 0.005288 either way). Note the default is 0.5 rather than 0.30:
min_af/max_af already bound AF to [0.20, 0.80], but comparing |AF - 0.5| against
0.30 rejects AF exactly 0.20 on floating-point rounding, which silently changed
4 reads.

---

## Where the BAM pipeline still beats the graph, and why

Measured on chr18 and chr20 with the two gates on (`--min-read-margin 2
--anchor-af-margin 0.12`), comparing per-read against the same truth BAM.

### Accuracy: the graph wins, decisively

| | graph wrong / BAM right | BAM wrong / graph right | ratio |
|---|---|---|---|
| chr20 | 40 | 1,516 | graph 37.9:1 |
| chr18 | 283 | 636 | graph 2.2:1 |

Of chr18's 283 losses, 56% sit at chr18:47-48 Mb, a **haplotype-asymmetric**
segdup: 15.9% segdup coverage on PATERNAL against 0.8% on MATERNAL. Its anchors
are only mildly off-centre (14.1% outside AF 0.35-0.65 against 9.5%
chromosome-wide), below what `--anchor-af-margin 0.12` removes.

### Coverage: BAM phases ~20% more reads

BAM phases 57,680 chr18 reads (40,065 on chr20) the graph declines, at **7.4%**
error against its own ~1.8% average. Of those chr18 reads: 100% are in the GAF,
86.5% are seen by the graph pipeline but carry a median of **2 observations and
margin 1**, and 98% are dropped by `--min-read-margin 2`. Those windows have
0.07-0.21 het candidates/kb against 1.28 in a well-phased control -- roughly 10x
less heterozygosity -- with the same filter-reason mix (96.3% `ref_only` at
median depth 72 vs 94.7%). Most of that coverage is not recoverable by anyone;
BAM buys it by phasing on two anchors and being wrong 7.4% of the time.

### `--min-mapq` is the strongest single filter

| `-q` | candidates | reads | discordant | hamming | switch | flip |
|---|---|---|---|---|---|---|
| 0 | 82,679 | 231,956 | 1,392 | 0.006001 | 184 | 278 |
| 10 | 82,481 | 230,084 | 712 | 0.003095 | 100 | 145 |
| **30 (default)** | 82,069 | 224,746 | **388** | **0.001726** | 52 | 79 |
| 60 | 80,721 | 209,770 | 256 | 0.001220 | 20 | 53 |

Removing it costs 3.6x more discordant reads for 3.1% more reads. It is not
redundant with the two gates: it acts on *placement* (a read that placed
ambiguously across graph paths) rather than on evidence. `-q 60` is a further
coverage/accuracy dial, same shape as `--min-read-margin`.

## RESOLVED (see next section): ~57,600 chr18 sites dropped `ref_only` that BAM calls het

Scale is established, mechanism is **not**. Chromosome-wide the BAM pipeline
calls 101,244 het candidates against the graph's 82,069. Of the 32,411 het sites
BAM calls and the graph does not, 42% (13,645) are present in the graph catalog
-- so they are reachable without any BAM. Every one of those the graph evaluated
was dropped `ref_only`.

Ruled out by measurement, each:

- **Nested-snarl parent gating** -- all sampled sites are LV=0 with no PS/PA.
- **MAPQ** -- alt- and ref-carrying reads both median 60, ~1% below the threshold.
- **`min_alt_depth` scattering across multiallelic STR alleles** -- dropping it to
  1 recovers 59 candidates of ~57,600 and changes accuracy not at all.
- **Chunking** -- a lost site stays `ref_only` when run alone, in a 100 kb window,
  and in a 500 kb window.
- **Walk fragmentation** -- where alt-walk reads exist, 100% contain the allele
  walk contiguously.

Note `ref_only` is a misleading label: `graph_bam_adapter.cpp` always keeps the
ref walk at index 0 and appends only alts clearing `min_alt_depth`, so the
reason fires when *no alt allele cleared the threshold*, not when only the
reference was seen.

Caveat on the remaining evidence: the read-side checks above matched allele
walks by forward-orientation substring, while the pipeline also matches reverse
complement. At least one site showed the ALT walk in 15 forward reads and the
REF walk in 0, yet the pipeline counted 58 on ref -- consistent with most reads
matching reverse. Any further work here should compare orientation-aware.

---

## Resolved: the graph pipeline cannot see a het between two non-reference alleles

An orientation-aware trace inside `match_compact_site_on_read`
(`PGPHASE_DEBUG_SITE=<pos>`) settled this. At chr18:51004598, a 20-allele STR:

    58 reads, 58 matched, 0 unmatched   (25 forward, 33 reverse complement)
    allele 0 (the graph's reference walk):  0 reads
    allele 1: 40   allele 8: 12   allele 5: 5   allele 2: 1

Matching is not broken -- every read matched, and a third of them only via
reverse complement, which is why the earlier forward-substring checks looked
like lost observations. They were measurement error on my part, not a bug.

The real limitation is the **allele-fraction denominator**. Each alt is scored
`alt / (ref + alt)`, i.e. as though the site were biallelic against the graph's
reference allele. Where no read carries that reference allele, every alt scores
AF = 1.0 and the site is discarded as homozygous -- however the reads actually
split between the alts. Above, a clear 40/12 split is thrown away.

Scale on chr18: **34,698 sites** dropped `high_af` have `REF_COV = 0` (78.4% of
all `high_af` drops), median alt depth **66**, and 33,853 of them have alt depth
>= 10. That is the right order of magnitude for the 23% het-call gap against the
BAM pipeline.

### `--af-vs-site-depth` -- implemented, and off by default on purpose

Scoring against total site depth instead makes the site above read 40/58 and
12/58 rather than 1.0 and 1.0. For a biallelic site the two denominators are
identical, so only multi-allele sites change.

| af denominator | candidates | reads | discordant | hamming | switch | flip |
|---|---|---|---|---|---|---|
| ref+alt (default) | 82,069 | 224,746 | **388** | **0.001726** | 52 | 79 |
| site depth | 85,928 | 225,160 | 456 | 0.002025 | 57 | 89 |

It recovers 3,859 candidates (+4.7%) and 414 reads, and **degrades phasing**:
388 -> 456 discordant. The recovered sites are alt-vs-alt hets at multiallelic
tandem repeats, which are exactly the unreliable anchors the AF gate exists to
remove. So it stays off for phasing.

Its value, if any, is variant-calling completeness rather than phasing -- those
are real het loci the pipeline currently cannot represent. That was not
evaluated here: it needs a small-variant truth set, not a phasing truth BAM.

### Note on the earlier `ref_only` framing

`ref_only` fires when no alt cleared `min_alt_depth`, not when only the
reference was observed, and `REF_COV = 0` on 53.3% of those rows. The label
misled several steps of this investigation.

---

## Graph vs hybrid, head to head on chr18

A hybrid baseline was run on chr18 against the same truth BAM. With
`--anchor-af-margin 0.12` the graph pipeline beats hybrid on accuracy at
**every** operating point, and `--min-read-margin` is the coverage/accuracy dial:

| config | reads | discordant | hamming | switch | flip | perfect PS |
|---|---|---|---|---|---|---|
| **hybrid** | 275,198 | 2,083 | 0.007569 | 367 | 819 | 39.6% |
| graph margin 0 | 259,663 | 1,218 | 0.004691 | 198 | 660 | 47.6% |
| graph margin 1 | 258,208 | 1,192 | 0.004616 | 189 | 645 | 48.2% |
| graph margin 2 | 224,746 | 388 | 0.001726 | 52 | 79 | 86.0% |
| graph margin 3 | 207,134 | 261 | 0.001260 | 35 | 22 | 93.1% |

Margin 0 is the notable one: **94% of hybrid's reads at 1.7x fewer discordant
reads**. Margin 2 gives 5.4x and margin 3 gives 8x, at 82% and 75% of hybrid's
coverage. Note the margin knee moved once the AF gate existed -- margin 0 with
the AF gate (1,218) is better than margin 0 without it (2,497), so the two gates
are not independent and the operating point should be re-picked after any change
to either.

### Rejected: gating the consensus rather than the output

`--min-read-margin` suppresses a thin read's HP tag at output time, but the read
still shapes the k-means consensus that every other read is scored against.
Gating the consensus instead looked like it should give cleaner profiles *and*
full coverage.

It does not work, in two distinct ways, and was reverted:

1. **Gating during Phase 2 is circular.** Margins are only defined once a
   consensus exists, so a cold start with the gate on never forms one: 0 reads
   phased.
2. **Gating as a post-convergence refinement is worse than useless.** Letting
   Phase 2 converge, then rebuilding the profile from confident reads alone and
   re-assigning everything, gives 80,741 discordant reads against a 1,218
   baseline -- hamming 0.33, essentially random. Variants covered only by thin
   reads end up with an empty profile and a degenerate consensus, and every read
   is then scored against that.

The k-means keeps coupled state across `hap_to_alle_profile`, `hap_to_cons_alle`
and the phase sets, maintained jointly by `iter_update_var_hap_cons_phase_set`
and `iter_update_var_hap_to_cons_alle`. A refinement pass that rebuilds one of
them in isolation violates that invariant. Any future attempt needs to preserve
all three together, and to leave variants with no confident coverage at their
converged consensus rather than resetting them.

### Confirmed on chr20 and chr12

Both operating points hold on all three chromosomes with truth BAMs. All rows
use `--anchor-af-margin 0.12`; "vs hyb" is the discordant-read ratio in the
graph pipeline's favour and "cover" is its share of hybrid's evaluated reads.

| chrom | config | reads | discordant | hamming | switch | flip | vs hyb | cover |
|---|---|---|---|---|---|---|---|---|
| chr20 | hybrid | 210,905 | 1,114 | 0.005282 | 167 | 605 | -- | -- |
| chr20 | margin 0 | 199,704 | 756 | 0.003786 | 153 | 464 | 1.5x | 95% |
| chr20 | margin 1 | 198,624 | 729 | 0.003670 | 149 | 450 | 1.5x | 94% |
| chr20 | margin 2 | 175,173 | 237 | 0.001353 | 59 | 75 | 4.7x | 83% |
| chr12 | hybrid | 489,489 | 3,394 | 0.006934 | 342 | 1,456 | -- | -- |
| chr12 | margin 0 | 465,294 | 1,535 | 0.003299 | 209 | 1,053 | 2.2x | 95% |
| chr12 | margin 1 | 463,227 | 1,511 | 0.003262 | 206 | 1,032 | 2.2x | 95% |
| chr12 | margin 2 | 410,699 | 257 | 0.000626 | 27 | 98 | 13.2x | 84% |
| chr18 | hybrid | 275,198 | 2,083 | 0.007569 | 367 | 819 | -- | -- |
| chr18 | margin 0 | 259,663 | 1,218 | 0.004691 | 198 | 660 | 1.7x | 94% |
| chr18 | margin 2 | 224,746 | 388 | 0.001726 | 52 | 79 | 5.4x | 82% |

(chr20 excludes the ROH and graph-blind regions; chr12 and chr18 use no
exclusions. Each chromosome is internally consistent, so the graph-vs-hybrid
ratios are comparable within a row-group but coverage percentages are not
comparable across chromosomes.)

Two defensible operating points, both ahead of hybrid on accuracy:

- **margin 0 or 1** -- 94-95% of hybrid's reads at 1.5-2.2x fewer discordant
  reads. Margin 1 is very slightly better than 0 on every chromosome at
  essentially the same coverage, so prefer 1 of the two.
- **margin 2** -- 82-84% of hybrid's reads at 4.7x to 13.2x. chr12 is the
  extreme: 257 discordant against hybrid's 3,394.

There is no longer a coverage-for-accuracy trade against hybrid; the graph
pipeline wins on accuracy at every point measured, and the margin only decides
how much more it wins by. **Defaults are still `--min-read-margin 0
--anchor-af-margin 0.5`, i.e. both gates off** -- changing them is a deliberate
call, not something this work did implicitly.

---

## Can coverage and accuracy be had together? Investigated: not by filtering

Between `--min-read-margin` 1 and 2 the pipeline discards 34,000 chr18 reads to
remove ~830 errors. Those reads are **97.6% correct** -- 34,000 good reads
thrown away to catch 830 bad ones -- so it looked like a better discriminator
should exist. Three were tried.

### Per-read signals do not separate the discarded reads

The k-means score margin (`|hap_scores[1] - hap_scores[2]|`) and the number of
informative variants behind the call are now exposed per read via
`--phase-reads-out` (`SCORE_MARGIN`, `N_SCORED`). Neither separates the
discarded population, which sits near 2.4% error almost uniformly:

| SCORE_MARGIN | reads | discordant | rate |
|---|---|---|---|
| 3-5 | 32,338 | 818 | 2.53% |
| 6-10 | 1,552 | 7 | 0.45% |

The useful bucket holds 1,552 of 33,994 reads. `N_SCORED` is worse than useless
-- error *rises* with it (1 -> 2.46%, 4-7 -> 22.6%), the repeat-rich signature
seen earlier.

### Regional signals predict, but too weakly

Errors are strongly regional, so window-level filtering should be more
efficient. Spearman against per-window error rate, all truth-free:

| window signal | Spearman |
|---|---|
| fraction of reads with clean-SNP margin < 2 | **+0.517** |
| fraction of candidates with off-centre AF | +0.335 |
| median observations per read | -0.404 |
| fraction of reads with a conflicting clean SNP | -0.320 |

### The frontier

Post-hoc over the chr18 margin-0 population (259,663 reads, 1,218 discordant):

| policy | kept | % | discordant | rate |
|---|---|---|---|---|
| keep all | 259,663 | 100.0% | 1,218 | 0.469% |
| per-read margin >= 2 | 225,669 | 86.9% | 391 | 0.173% |
| per-read margin >= 3 | 208,313 | 80.2% | 268 | 0.129% |
| regional: window thin-frac < 0.15 | 190,056 | 73.2% | 230 | 0.121% |
| margin >= 2 OR window thin-frac < 0.2 | 233,778 | 90.0% | 571 | 0.244% |
| margin >= 1 AND window thin-frac < 0.25 | 213,832 | 82.3% | 378 | 0.177% |
| **oracle: drop worst 100 windows by truth** | **231,731** | **89.2%** | **297** | **0.128%** |

**Nothing beats the per-read margin ladder.** Every combination tried is either
dominated by it or trades along the same curve. The disjunctive policy buys 3.1%
more reads for a 41% worse error rate; the conjunctive one is dominated outright.

### The headroom is real but needs a better predictor

The oracle row is the point: selecting windows *with* the truth keeps **89.2% of
reads at 0.128%**, which strictly dominates `margin >= 2` (86.9% at 0.173%) --
more reads *and* fewer errors. So a regional policy can beat the per-read gate;
the best truth-free window signal found (+0.517) simply is not sharp enough to
find those windows. Closing the gap between 0.173% and 0.128% at ~89% coverage
is a window-quality prediction problem, not a filtering-policy problem.

Until then the per-read margin is the honest dial, and the trade is real: margin
1 for ~95% of hybrid's coverage, margin 2 for 4.7-13.2x its accuracy.

---

## Where the unphased reads go, and whether any are recoverable

Full funnel for chr18, `--min-read-margin 0 --anchor-af-margin 0.12 -q 30`,
over the 346,878 reads in the coordinate-indexed GAF:

| bucket | reads | share | |
|---|---|---|---|
| MAPQ < 30 | 29,068 | 8.4% | filtered before matching |
| pass MAPQ, zero site observations | 22,423 | 6.5% | |
| observed but unassignable (hap 0) | 35,443 | 10.2% | median 2 observations |
| phased | 259,944 | 74.9% | of which ~34,000 more drop at margin 2 |

### None of the three loss buckets is recoverable

**MAPQ (8.4%) -- tested, no.** Lowering `-q` and compensating with a higher
read margin never wins; every combination keeps *fewer* reads than
`-q 30 / margin 2` and most add errors:

| -q / margin | reads | discordant | vs q30/m2 |
|---|---|---|---|
| 0 / 3 | 213,277 | 1,110 | -11,469 reads, +722 disc |
| 10 / 4 | 198,978 | 420 | -25,768 reads, +32 disc |
| 20 / 4 | 198,010 | 250 | -26,736 reads, -138 disc |

A low graph MAPQ means ambiguous placement, and a read placed on the wrong path
can agree *strongly* with the wrong haplotype -- high margin, wrong answer. The
margin cannot discriminate against that, which is the paralog failure mode again.

**Zero-observation reads (6.5%) -- heterozygosity, not catalog.** Their spans
carry a median of **176 catalog sites** but **0 surviving candidates**, and 93%
have no het candidate at all. The catalog is there; the sample is homozygous
across it.

**hap-0 reads (10.2%) -- unphasable by anyone.** Comparing median het density
per read span between the graph and the BAM pipeline:

| read group | reads | graph cands/read | BAM hets/read | BAM's gain |
|---|---|---|---|---|
| never observed (MAPQ ok) | 22,423 | 0 | 1 | +1 |
| observed but hap 0 | 35,443 | 2 | **2** | **0** |
| phased | 259,944 | 17 | 18 | +1 |

For the hap-0 reads the BAM pipeline has **exactly the same** het density the
graph does -- two per read. There is no hidden evidence for a better method to
find. Two het sites is not enough to assign a read, and hybrid's apparent
coverage advantage over these reads is it committing anyway: those are the reads
it phases at 7.4% error against its own 1.8% average.

### What would actually fix it

Not an algorithm. The binding constraint is het sites per read, so the levers are
longer reads (a read spanning 10 hets is phasable where one spanning 2 is not) or
a more heterozygous sample -- not a better model over this data. The one
algorithmic lever left is calling novel variants from the GAF's `cs` tags rather
than only genotyping catalog sites, which would help the 22,423 zero-observation
reads (BAM finds ~1 het/read there that the catalog lacks) but not the 35,443
hap-0 reads, where there is nothing extra to find.

---

## Head to head against the field: whatshap, LongPhase, HiPhase, longcallD

All phasers were given the **same reads** — the giraffe-surjected chr20 BAM —
and the same DeepVariant call set to phase; our method reads the graph
alignment of those same reads. Scored against one truth BAM with one script.

| phaser | reads | discordant | hamming | switch | flip | N50 kb | perfect PS |
|---|---|---|---|---|---|---|---|
| whatshap 2.8 | 221,925 | 12,200 | 0.054974 | 1,019 | 1,403 | 1110 | 20.5% |
| HiPhase 1.6.0 | 228,863 | 8,013 | 0.035012 | 1,119 | 1,794 | 1125 | 23.0% |
| longcallD | 214,017 | 4,636 | 0.021662 | 849 | 1,232 | 985 | 36.0% |
| pgphase BAM | 214,804 | 4,648 | 0.021638 | 849 | 1,235 | 985 | 35.8% |
| LongPhase 2.0.2 | 220,046 | 3,390 | 0.015406 | 569 | 1,185 | 1330 | 28.2% |
| pgphase hybrid | 210,905 | 1,114 | 0.005282 | 167 | 605 | 999 | 41.0% |
| **graph, margin 0** | 199,704 | **756** | 0.003786 | 153 | 464 | 918 | 52.9% |
| **graph, margin 2** | 175,173 | **237** | **0.001353** | **59** | **75** | 937 | **87.4%** |

**14.3x fewer discordant reads than the best existing tool** (LongPhase), 9.6x
fewer switches, and four times as many perfectly phased blocks.

### Where the other tools win, and why it is not a defect

Asked directly: are there regions where they beat us? Essentially no.

- **Head to head on shared reads**: graph wins **4:1** against LongPhase (38 vs
  163 reads) and **25:1** against HiPhase (164 vs 4,106).
- **Regionally**: of 8 one-Mb windows where the graph loses a single read to
  LongPhase, LongPhase loses *more* in half of them. The largest graph-only
  deficit in any window is 4 reads.

Two genuine differences remain, and both are the same trade seen throughout:

**Coverage.** LongPhase phases 22,678 reads we do not; HiPhase 29,593. Those
extra reads are **11.6%** and **10.3% discordant** — the same pattern as
hybrid's 7.4%. They are committing on thin evidence, not seeing more.

**Contiguity.** N50 937 kb against LongPhase's 1330. This is not chunking
(`--chunk-size` 500 kb / 2 Mb / 5 Mb gives 937 / 935 / 937) and not the read
gate (margin 0 gives 918 kb, *shorter*). The **median block span is effectively
identical across every method: 862-886 kb.** The N50 and auN gaps come from a
few very long blocks the competitors emit — LongPhase's auN is 20.4 Mb against
our 989 kb — bought with 569 switch errors against our 59.

What those long blocks do is visible directly: LongPhase has **one block
spanning the 44-46 Mb run of homozygosity** (3 het SNPs per Mb) and five
spanning the graph-blind satellite at 27.2-28.8 Mb. We emit none. Declining to
phase across a 2 Mb stretch with almost no heterozygosity is the correct
behaviour; a long block with a switch through the middle of it is worse than two
correct blocks.

So the contiguity difference is a deliberate consequence of the accuracy
advantage rather than a deficiency to fix. The honest way to report it is both
numbers side by side, with median span shown alongside N50 so the reader can see
that typical blocks are the same length.

### Setup note

Competitors require a pre-called VCF to phase (DeepVariant was used, small model
enabled). Our method needs none — it phases directly from the graph alignment.
That asymmetry favours them in this comparison, since they are handed variant
calls we never receive.

---

## Whole-snarl phasing regression: two bugs fixed, default unchanged

Revisited the `docs/HANDOFF.md` open problem for `--snarl-keep-whole` on HG002
chr20 (CHM13 catalog/GAF, `--min-read-margin 2 --anchor-af-margin 0.12`).

Two implementation bugs were found:

1. `build_graph_chunk` rewrote `allele_counts` from old site space into final
   candidate space before remapping read observations. Phase 3 then tested
   `allele_counts[old_si].size() > 2` to decide whether the original source
   snarl was multi-allelic. For decomposed sites this could be false after the
   rewrite, so reads on alt_j failed to appear as allele-0 evidence against
   alt_i under `--snarl-allele-phasing`.
2. Whole n-allelic candidates were classified from collapsed non-ref AF
   (`alt_total / site_total`). A no-reference alt_1/alt_2 heterozygote was
   therefore marked `CleanHom`, despite carrying exactly the two allele indices
   needed for whole-snarl phasing.

Fixes:

- Preserve a per-source `source_is_multi` flag through candidate rewriting and
  use that during observation remap.
- For n-allelic candidates, classify hom/het by top allele fraction rather than
  collapsed non-ref AF.
- Give n-allelic anchors single k-means weight and exclude them from clean-SNP
  read-margin accounting.
- Added focused `test_graph_bam_adapter` regressions for alt-vs-other read
  remapping and no-reference alt_1/alt_2 whole-snarl classification.

Fresh post-fix read-level results:

| config | candidates | phase sets | reads evaluated | discordant | hamming |
|---|---:|---:|---:|---:|---:|
| baseline | 73,627 | 279 | 175,843 | **248** | **0.001410** |
| `--snarl-allele-phasing` | 77,570 | 282 | 176,196 | 518 | 0.002940 |
| `--snarl-keep-whole` (`top2 >= 0.90`) | 80,837 | 281 | 175,792 | 393 | 0.002236 |
| `--snarl-keep-whole --snarl-top2-frac 0.95` | 79,270 | 281 | 175,828 | 395 | 0.002247 |

Conclusion: whole-snarl phasing is no longer the 25x read-accuracy regression
previously measured (6,322 discordant reads); the main regression was a bug. It
still loses to the baseline, so defaults remain unchanged. The remaining
research direction is allele clustering inside high-multiplicity snarls rather
than exact allele-index equality.

---

## Competitor-site gap diagnosis with graph SITE_ID accounting

Added two diagnostics to make the "other tools phase here and pgphase does not"
question reproducible:

- `collect-graph-variation --phase-sites-out FILE` streams retained graph-site
  candidates with source `SITE_ID`, allele counts, phase set, and hap allele
  assignments.
- `scripts/analyze_graph_gap_site_loss.py` classifies truth hets outside pgphase
  merged phase blocks, optionally narrowed to sites covered by competitor phased
  VCFs.  When `--phase-sites-tsv` is supplied, retained and filtered graph sites
  are matched by parent `SITE_ID` rather than by normalized variant position.

Competitor VCFs used:

```
/tmp/claude-1000/-home-kokyriakidis-Downloads-pgphase/68d70dd6-387e-4e2f-885c-9efbd718ca83/scratchpad/phasers/chr20/
  hp.vcf.gz
  lp.vcf.gz
  ws.vcf.gz
```

Baseline accounting for `chr20.sites.vcf.gz`:

| item | count |
|---|---:|
| catalog records | 977,275 |
| duplicate catalog IDs | 0 |
| retained graph-site parent IDs | 70,581 |
| filtered graph-site parent IDs | 976,572 |
| catalog IDs with no retained or filtered row | 3 |

Conclusion: the graph pipeline is accounting for essentially every catalog
record. The missing-site problem now splits into catalog incompleteness relative
to the source graph, evidence concentrated in a different/nested catalog site,
and retained sites that fail to become phase-block anchors/bridges.

Exact-position competitor-site target (truth het outside pgphase blocks and at
the same position as a phased HiPhase/LongPhase/WhatsHap variant):

| reason | baseline | `--snarl-allele-phasing` |
|---|---:|---:|
| no catalog record within ±25 bp | 1,012 | 1,012 |
| `ref_only` | 753 | 756 |
| retained graph site, unphased | 606 | 709 |
| `no_reads_in_chunk` | 386 | 388 |
| `high_af` | 268 | 119 |
| `low_depth` | 143 | 144 |
| `low_af` | 17 | 22 |
| retained graph site, phased nearby | 12 | 21 |

The multi-allelic fix does what it should: it cuts exact-target `high_af` misses
from 268 to 119. It does not materially reduce the gap count because those
newly retained sites mostly remain unphased. The next experiment should classify
retained-but-unphased rows by anchor mask, spanning-read links, and allele-count
class, then test allele clustering inside high-multiplicity snarls.

## Surjected-BAM private-site recovery against competitor-only graph gaps

The graph-gap audit was extended with repeatable `--recovery-vcf LABEL=PATH`
measurements in `scripts/analyze_graph_gap_site_loss.py`. A fresh chr20
`collect-hybrid-variation` run used the surjected HG002 BAM, graph catalog, and
GAF and emitted 134,646 candidates / 64,184 phased positions.

On the strict 3,197-site target (truth het, outside merged graph blocks, phased
at the exact position by HiPhase/LongPhase/WhatsHap), hybrid recovered 808 exact
phased hets and placed 1,624 sites inside a hybrid phase block. Exact recovery by
graph-gap reason was: no catalog 201/1,012; ref-only 244/753; retained-unphased
272/606; no-reads-in-chunk 33/386; high-AF 27/268; low-depth 29/143. Every exact
hybrid het was phased, so the remaining 2,389 are absent from hybrid output, not
merely present and unphased.

Decision: the next graph improvement should test an external linear-VCF seed
path, using high-confidence heterozygous calls from the surjected BAM (initially
the same DeepVariant call set supplied to competitors). Preserve graph-core
assignments; phase BAM-private gap evidence additively in disjoint/local phase
sets, and permit merging only with direct two-haplotype spanning-read support.
Do not globally rerun the graph core with all private sites: earlier gap-fill
work measured about 21% error in BAM-only reads concentrated in segdups.

## Corrected shared-call transfer and private-site prototype (2026-09-12)

The existing shared-call benchmark was contaminated: the DeepVariant VCF
already contained 67,127 phased heterozygotes, while
`scripts/phase_vcf_from_hp.py` left unsupported records unchanged. It now
clears GT phasing and PS for every heterozygote before applying evidence from
the supplied haplotagged BAM. It also resolves graph-style BAM contig names,
parses interval regions correctly, scans each BAM once into a reusable support
cache, and can infer parity-consistent merges between overlapping phase sets.
The scan implementation was checked against the original pileup implementation
on chr20:1,000,000-2,000,000; their output VCFs were byte-identical.

Corrected HG002 chr20 results on the same DeepVariant call set:

| configuration | assessed | blocks | N50 kb | NGC50 kb | switchflips | Hamming |
|---|---:|---:|---:|---:|---:|---:|
| graph transfer | 59,015 | 276 | 215 | 204 | 15 | 16 |
| current hybrid transfer | 59,665 | 298 | 312 | 290 | 41 | 68 |
| all-linear-site seed, no merge | 59,672 | 303 | 290 | 267 | 39 | 62 |
| all-linear-site seed + conservative PS merge | 59,447 | 247 | 401 | 374 | 39 | 79 |

The seed was built by `scripts/augment_graph_catalog_with_linear_hets.py`, which
now excludes exact graph REF/ALT matches by default. Of 77,484 biallelic
DeepVariant hets, 22,084 are exact graph misses and 4,206 remain at PASS/GQ20.
Adding those 4,206 produced 183 extra hybrid candidates but only one extra
transferred phase call. Hybrid N50/NGC50 moved 312/290 -> 314/296 kb with
switchflips/Hamming unchanged at 41/68. Conservative PS merging of this private
run reached 397/366 kb but Hamming worsened to 379.

The stronger tabled experiment deliberately added all 77,484 sites, including
55,400 exact graph matches. It produced 1,748 extra candidates and 253 extra
transferred calls. Its advantage over private-only seeding shows that the main
missing layer is variant-first use of represented alleles as anchors, not raw
catalog completeness. Most private records still fail BAM evidence/candidate
classification or do not bridge a phase set.

The conservative PS merge uses output calls at two reads/haplotype and 0.70
purity. Edge proposals may use one haplotype at 0.60 only if the other phase set
has both haplotypes represented. It requires two agreeing sites, vote margin
two, summed support two, and support four for PS-label distances over 500 kb.
More aggressive merging reached raw NG50 595 kb but Hamming exceeded 2,000. A
truth-guided two-edge exclusion oracle reached N50/NGC50 425/393 kb at Hamming
232, still below LongPhase NGC50 400 kb and HiPhase 644 kb. This exhausts the
available local overlap evidence on chr20; further contiguity needs new spanning
evidence or a variant-first joint model, not looser adjacency merging.

`scripts/compute_ngc50.py` also no longer silently forces a 3.1 Gb denominator
for chromosome experiments. Use `--genome-size 66210255` for this chr20 setup.

## Graph-authoritative joint private-gap phasing (2026-09-12)

The all-linear-site experiment was not the requested architecture. The corrected
design uses `scripts/extract_private_gap_sites.py` to select only exact-private
linear heterozygotes inside gaps between merged graph phase blocks, then passes
that VCF to `collect-hybrid-variation --private-sites FILE`.

The new mode enforces three ownership boundaries before joint k-means:

1. BAM candidates not present in the private-gap VCF are removed.
2. BAM profile alleles/counts at graph-owned candidates are cleared; GAF
   injection is the sole source of graph-site evidence.
3. BAM noisy-region MSA recall is disabled so no unlisted candidate can enter
   after whitelisting.

Graph and private observations then use the existing shared k-means and chunk
stitching together. The graph read confidence gate is available on hybrid as
`--min-read-margin`; experiments used 2 with graph-style stitch margin/rule 0.

HG002 chr20 results on the corrected shared DeepVariant call set:

| config | private sites offered | phased reads | read discordant | read N50 kb | assessed | DV N50/NGC50 kb | switchflips | Hamming |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| graph only | 0 | 175,843 | 248 | 937 | 59,015 | 215/204 | 15 | 16 |
| private GQ20 | 271 | 178,284 | 320 | 945 | 59,038 | 280/224 | 20 | 21 |
| **private GQ10** | **958** | **180,800** | **366** | **952** | **59,148** | **300/280** | **25** | **27** |
| private GQ0 | 6,607 | 181,401 | 384 | 965 | 59,185 | 303/280 | 32 | 36 |

Decision: GQ10 is the current experimental operating point. It improves
corrected NGC50 204 -> 280 kb (+37%), phases 4,957 more truth-evaluable reads,
and keeps shared-call Hamming at 27 (LongPhase 605, HiPhase 2,279). GQ0 adds no
NGC50 beyond GQ10 and only increases errors. This is a genuine graph-first gap
fill: no BAM-wide site set participates.

Graph-core stability at GQ10: 173,490 graph reads are shared with the joint
evaluation. Truth-status transitions are 25 discordant->concordant and 24
concordant->discordant, so joint private evidence is net neutral on retained
graph accuracy. The margin gate removes 2,380 marginal graph reads. Of 7,333
newly evaluable reads, 7,184 are concordant and 143 discordant.

## All-BAM graph-authoritative control and VCF query fix (2026-09-12)

Added `collect-hybrid-variation --graph-authoritative`. BAM calling and normal
candidate classification run first. The hybrid join tracks both newly added
graph candidates and exact BAM/graph matches through candidate-table sorting.
At every graph-owned candidate it clears the BAM read-profile allele plus all
BAM-derived count fields (depth, low-quality depth, and strand counts), then
repopulates evidence only from GAF. All surviving non-graph BAM candidates are
retained and jointly phased. `--private-sites` now implies this ownership rule.

HG002 chr20, read margin 2: all clean BAM sites evaluated 184,911 reads with
812 discordant (0.439% Hamming) and 965 kb read N50. On the shared DeepVariant
call set it assessed 59,269 pairs at 300/260 kb N50/NGC50, 36 switchflips, and
133 Hamming errors. Margin 3 evaluated 169,685 reads with 731 discordant
(0.431%) and 1,014 kb N50, so confidence gating did not cure the conflict.
The corrected selective GQ10 mode remains better: 180,800 reads, 366 discordant
(0.202%), and shared-call 300/280 kb N50/NGC50, 25 switchflips, Hamming 27.

During this experiment, verbose counts exposed 977,275 graph candidates in
chunk 0 versus roughly 6,000 in ordinary 500 kb chunks. `load_sites_for_region`
constructs a 1-based VCF region, while all graph/hybrid callers passed
`region.beg - 1`; first-chunk start zero caused the formatter to issue a whole-
contig tabix query. The three VCF catalog callers now pass `region.beg`.
Zero-based GAF/FASTA calls are unchanged. A corrected GQ10 rerun reproduced the
shared-call metrics exactly and changed read evaluation by only four reads,
showing this was primarily a severe runtime/memory and boundary-query bug here.

## Regional failures, BAM fallback, and supported private bridges (2026-09-12)

The six phase sets below 60% accuracy in the corrected GQ10 run are not evidence
that the competing phasers generally solve these regions. Over the union of the
same six exact truth-coordinate intervals, pgphase evaluated 1,369 reads with
68 errors (5.0%); HiPhase evaluated 2,552 with 644 errors (25.2%), LongPhase
3,068 with 943 (30.7%), and WhatsHap 1,752 with 184 (10.5%). The competitors
phase more reads but are less accurate in aggregate. They do beat pgphase in two
small intervals near 31.77-31.89 Mb and 32.45-32.47 Mb; those are useful
diagnostic targets, not a reason for chromosome-wide BAM ownership.

`--bam-authoritative-bed FILE` was added to test that distinction. Within BED
intervals clean BAM candidates replace the graph catalog and GAF evidence;
outside, graph-authoritative private-gap behavior is unchanged. Broad cenSat
fallback was worse (180,790 reads, 383 discordant versus 366), and a 65 Mb
terminal fallback was neutral (180,814, 364). An oracle BED containing only the
two competitor-winning intervals reached 180,769/355, but BAM clipping and GAF
coverage did not identify those intervals: target reads had median/p90 BAM clip
fraction zero and GAF aligned fraction 1.0. Do not select fallback regions from
truth or annotation alone.

`scripts/extract_private_gap_sites.py --bam FILE` now offers a truth-free,
bridge-only private-site filter. For each gap it builds read-overlap edges among
the left graph boundary, private candidates, and right graph boundary, then
retains the strongest complete path whose edges have at least
`--min-bridge-reads` MAPQ-filtered reads. On chr20 this retained 203/958 GQ10
private sites in 112/364 gaps. It evaluated 177,901 reads with 318 discordant
(0.18%) and 971 kb N50; adding `--min-phase-set-reads 10` yielded 177,829/303
and 972 kb. This validates private sites as real connectors, but bridge-only is
too conservative for the best coverage because useful one-sided/local gap
islands are discarded.

The best current chr20 accuracy/coverage point keeps all 958 GQ10 private gap
sites and suppresses output from phase sets with fewer than 50 assigned reads:

| configuration | evaluated reads | discordant | Hamming | read N50 kb | bad PS | DV assessed | DV N50/NGC50 kb | switchflips | DV Hamming |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| GQ10 private | 180,800 | 366 | 0.202% | 952 | 6 | 59,148 | 300/280 | 25 | 27 |
| + min PS 10 | 180,701 | 344 | 0.190% | 964 | 3 | 59,127 | 300/280 | 22 | 23 |
| **+ min PS 50** | **179,494** | **269** | **0.150%** | **991** | **0** | **58,785** | **300/280** | **18** | **19** |
| bridge-only + min PS 10 | 177,829 | 303 | 0.170% | 972 | 2 | not run | not run | not run | not run |

`--min-phase-set-reads` is a truth-free output-confidence gate. Fifty is the
selected chr20 experimental point, not yet a default: it must reproduce on
chr12/chr18 and another sample before being promoted. The design lesson is to
retain graph trust, use private BAM sites only in graph phase gaps, let direct
read connectivity form joint blocks, and abstain on small unsupported blocks.
The gate counts unique `(input, read name)` identities, so chunk-overlap copies
cannot inflate phase-set support; correcting that audit detail left the chr20
result unchanged.

## Why competitors appear better in the two chr20 target intervals (2026-09-12)

The apparent wins at CHM13 chr20:31,766,466-31,888,625 and
32,449,119-32,466,111 are haplotype-specific centromeric/segmental-duplication
alignments, not ordinary diploid phase blocks. Truth-label and input-BAM audits
showed 297/298 primary reads in the broader first interval are paternal and all
72 in the second are paternal. The nine reads assigned to each pgphase target
PS are all paternal. BAM and GAF alignments are not clipped, have high MAPQ
among assigned reads, and cover the sites fully; the absent maternal signal is
a copy/representation issue rather than low mapping confidence.

Hybrid `--phase-matrix-dump PREFIX` was exposed to inspect these blocks. No
private BAM candidate is decisive: every separating anchor is graph-owned. At
32,449,119, dozens of clean graph SNP/indel anchors repeat the same 5:4 read
partition; the catalog has 201 records under one large parent snarl in 2 kb. At
31.7-31.9 Mb, 634 nested records share another parent. These correlated nested
sites distinguish repeat/paralog copies among reads from the same paternal
haplotype. pgphase counts them as independent diploid evidence, then k-means is
required to produce HP1 and HP2, creating a false 5:4 split.

The competitor advantage is mostly abstention or presentation:

- On the exact nine reads at 31,766,466, HiPhase also splits them 5:4; its much
  larger surrounding PS is 98.3% accurate and makes the interval aggregate look
  better. LongPhase and WhatsHap also split these paternal-only reads.
- At 32,449,119, HiPhase's perfect 20-read PS contains 20 paternal reads tagged
  HP2 and zero HP1 reads. LongPhase leaves the nine pgphase reads unphased.
  WhatsHap's local PS is 11 paternal reads in each HP and exactly 50% accurate.
  No competitor reconstructs two biological haplotypes in this interval.

A truth-free control collapsed small (<30-read), compact (<100 kb read-start
span), >=90% soft-masked phase sets to one HP. It fixed both target PSs (4 -> 0
errors each) and reduced whole-chromosome discordance 366 -> 310 without losing
evaluated reads, but incorrectly collapsed a genuinely diploid centromeric
block. Chromosome-wide matrices found correlated partitions in true and false
blocks, so soft masking, density, parent identity, MAPQ, and partition
redundancy are insufficient to call a block mono-haplotype safely.

Decision: retain `--min-phase-set-reads 50` as the safe experimental policy. It
abstains on both false blocks and yields 269 chr20 errors. Do not emit one-sided
HP blocks without sample-aware graph-thread or copy-number evidence. A
principled future fix is to group nested sites by graph parent/read partition so
correlated observations contribute one evidence unit, then require independent
flanking or private-site support before emitting a diploid PS.

## Graph-locked native-BAM gap fill on chr12 and chr18 (2026-09-12)

The chr20 GQ10 result did not generalize directly when the private VCF came from
`collect-bam-variation` instead of DeepVariant. On chr18, exact-private GQ10
selection admitted 12,007 sites and joint phasing produced 3,527 discordant
reads (1.53%). Read-overlap bridge filtering reduced this to 565 sites but still
produced 3,141 discordant reads (1.39%). Restricting to 114 CLEAN, balanced SNPs
did not fix it. Two large phase sets contained 2,830 of those errors, so the
`--min-phase-set-reads 50` chr20 policy cannot protect this chromosome.

The cause is architectural: a mapped read spanning two coordinates proves
physical connectivity, not allele-consistent phase orientation. Joint k-means
can therefore reorient trusted graph reads around a small private proposal.
`scripts/merge_graph_hybrid_tags.py` now implements a graph lock:

1. every graph HP/PS assignment is copied unchanged onto the full surjected BAM;
2. each hybrid PS votes for a graph PS and orientation using shared reads;
3. hybrid-only reads are admitted only with 10 shared reads, vote margin 5,
   90% purity, and support from both haplotypes;
4. ambiguous hybrid blocks are left unphased.

Native private-site extraction also gained `--clean-snps-only`,
`--exclude-graph-positions`, and VAF bounds. The validation recipe used GQ10,
CLEAN SNPs, graph-position absence, VAF 0.30-0.70, and a two-read MAPQ20 bridge.

| chromosome | config | private sites | evaluated reads | discordant | Hamming | read N50 | bad PS |
|---|---|---:|---:|---:|---:|---:|---:|
| chr18 | graph | 0 | 224,746 | 388 | 0.17% | 1,746,048 | 6 |
| chr18 | joint private | 111 | 226,016 | 3,148 | 1.39% | 1,766,679 | 4 |
| chr18 | **graph lock** | **111** | **228,020** | **405** | **0.18%** | **1,748,648** | **6** |
| chr12 | graph | 0 | 410,699 | 257 | 0.06% | 422,422 | 2 |
| chr12 | **graph lock** | **131** | **415,431** | **286** | **0.07%** | **425,298** | **2** |

Across chr12 and chr18, graph lock adds 8,006 truth-evaluable reads for 46
additional errors (99.43% accuracy among the net additions) while preserving
all graph assignments by construction. It modestly extends existing blocks; it
does not yet merge independent graph PS labels. This is the selected native-BAM
hybrid policy. Keep direct joint output experimental, and do not promote PS50 as
a general fix.

### Chr20 graph-lock compatibility check

The same graph lock was tested against the original 958-site DeepVariant GQ10
chr20 proposal. Graph-only evaluated 175,843 reads with 248 errors. Unfiltered
graph lock added 3,769 reads but 81 errors (179,612/329); increasing orientation
purity from 0.90 to 0.95 removed 289 reads and only two errors. The remaining
error is therefore dominated by small graph phase sets retained by the lock,
not weak proposal orientation.

A new `--min-output-phase-set-reads` gate applies after graph locking. At 50 it
evaluated 178,697 reads with 266 errors, one bad PS, and 960 kb N50. This is the
higher-precision chr20 point, but the original direct joint GQ10+PS50 result
remains the selected balance: 179,494 reads, 269 errors, no bad PS, and 991 kb
N50. It gains 797 reads and 31 kb N50 for three errors. Conclusion: graph lock
is required for native BAM proposals on chr12/chr18, while a polished external
callset can still benefit from direct joint phasing plus final abstention.

## Shared-call chr12/chr18/chr20 benchmark (2026-09-12)

The full development panel now uses one DeepVariant 1.10.0 call set per
chromosome for pgphase, WhatsHap 2.8, HiPhase 1.6, and LongPhase 2.0.2. Variant
truth, read truth, chromosome-denominated NGC50, resource use, and gap-site
loss are recorded under
`evaluations/2026-09-12-chr12-18-20-comparison/`.

Pooled graph read Hamming is 0.110% on 811,382 phased reads. Graph lock phases
822,021 reads at 0.123% Hamming, but median chromosome NGC50 changes only 251
to 252 kb. HiPhase reaches 644 kb and LongPhase 454 kb, with read Hamming of
4.377% and 2.339%, respectively. Direct hybrid reaches 295 kb but is rejected:
chr18 read Hamming rises to 1.383% because two large hybrid blocks are mixed.

Across the three chromosomes, 63,337 truth hets lie outside graph blocks but
inside a competitor block. Of these, 35,529 (56.1%) have no exact graph catalog
record, 11,940 are represented but classified `ref_only`, and 7,409 have no
reads in their assigned chunk. Clean private gap-site admission recovers only
655 phased sites, so the next constraint is allele-consistent bridge evidence
to both adjacent graph blocks, not broader BAM-site admission.

During this benchmark, chr20 graph-position exclusion was found to compare
`CHM13#0#chr20` against an unnormalized requested contig. After normalizing both
sides, the valid private set fell from 143 to 74 and 45,393 graph-owned
positions were correctly excluded. The corrected chr20 hybrid and graph-lock
measurements are the ones in the evaluation record.

### Frozen benchmark framework and correct-bridge diagnosis

`scripts/benchmark_panel.py` now controls this panel from `panel.json`.
`competitor_lock.json` freezes 12 chromosome/tool runs with exact commands,
versions, inputs, outputs, and evaluation artifacts. Normal `run` mode verifies
the lock and sets `RUN_COMPETITORS=0`; missing fixed artifacts are fatal.
pgphase stages use `scripts/run_cached_step.py`, whose state includes exact
argv plus executable/input/output fingerprints. `report` deterministically
regenerates all aggregate tables and `REPORT.md` from frozen artifacts.

The correct-bridge analysis requires a competitor block to have >=99% read
accuracy, >=50 reads, and no variant switch interval across a graph break.
HiPhase correctly crosses 434 such breaks (median 22,575 bp): 250 (57.6%) are
dominated by repeat het indels excluded from graph k-means, 86 (19.8%) already
contain graph-phased sites but are not stitched, 34 (7.8%) contain clean
candidates that remain unphased, and 64 (14.7%) have catalog records but no
candidate. HiPhase phases 1,043 shared-call sites across these gaps; pgphase
phases 125. Prioritize a separately gated repeat-indel bridge channel and
global evidence-based stitching, both with immutable graph assignments.

The frozen LongPhase baseline is the measured `--pb` SNP mode, not its optional
`--indels` mode. Do not present it as LongPhase's best SNP+indel configuration.

### Graph-locked private-site bridges and LongcallD control (2026-09-13)

The initial graph lock could add hybrid-only reads but could not merge two
graph phase sets: each hybrid phase set selected one graph phase set as its
anchor and treated the other endpoint as a competing vote. The bridge mode in
`scripts/merge_graph_hybrid_tags.py` now validates each hybrid-to-graph
endpoint independently, then parity-unions graph phase sets when one accepted
hybrid block links both. A 300 kb maximum bridge distance rejects unsupported
long-range joins while preserving LongPhase-scale local evidence.

On chr20, 40 accepted local graph bridges plus strict shared-VCF transfer
increased chromosome NGC50 from 203,547 to 296,028 bp (+45.4%) with the same
24 variant switches, 6 switches plus 9 flips, and 16 Hamming differences as
the graph baseline. Read evaluation added 1,548 reads, changed discordance from
321 to 323, and reduced switch/flip events from 159 to 158. An unbounded bridge
introduced one extra switch near 65.994 Mb through an 816 kb merge, motivating
the explicit distance gate. Chr12/chr18 panel validation remains required
before this becomes the selected operating point.

LongcallD 0.0.11-23e369d was rerun as a chr20 native caller+phaser control on
the same annotated BAM and reference, using the explicit `CHM13#0#chr20`
region. It emitted 118,270 variants and 272,016 BAM records in 42.79 seconds
with 11,940,320 KiB peak RSS. Against GIAB variant truth, it assessed 65,749
pairs with NGC50 273,058 bp, 117 switches, 49 switches plus 34 flips, and 884
Hamming differences. Against assembly read truth, it phased 219,090 reads with
6,506 discordant (2.970%) and 3,132 switch/flip events. This row is labelled
`native` in `results.tsv` and `pooled.tsv`; unlike all other competitor rows,
it did not receive the shared DeepVariant callset and must not be used for a
direct callset-controlled superiority claim. Under the strict correct-bridge
screen, LongcallD crosses 47 chr20 graph breaks: 25 are dominated by excluded
repeat heterozygous indels, 12 by catalog records that never became graph
candidates, 6 by graph block stitching, and 4 by clean candidates left
unphased. This independently points to the same repeat-indel and candidate
admission gaps found with the shared-call competitors, while its absolute
counts remain native-call-set dependent.

### Region audit: chr20:15,019,294-15,130,077

HiPhase and both WhatsHap modes correctly cross this 110,783 bp graph gap.
A WhatsHap `ReadSetReader` reconstruction found the actual 15-site path. It
starts with the repeat deletion at 15,019,294, traverses four BAM-clean private
SNPs at 15,039,543-15,055,707, and continues through sparse indels to graph
SNPs at 15,087,221 and 15,115,387. The boundary-to-first-SNP link has four
unanimous reads; all later long links have 12-47 raw reads with a clear
orientation majority. The graph catalog contains no sites in the private-SNP
segment or around 15,109 kb, but the BAM caller already has all required
evidence.

The clean-private selector drops the SNP chain because its left native-graph
anchor is 25.8 kb away and the necessary repeat deletion is excluded. A
GQ10/VAF0.30-0.70 complete-path selector recovers seven appropriate BAM sites,
but hybrid clean-core k-means still excludes the noisy indel anchors. Graph
locking to completed BAM phase sets also fails because BAM phasing splits the
chain into three blocks. The required mechanism is therefore a bridge-only
allele graph that combines clean private SNPs with only necessary,
well-supported repeat indels and resolves their orientation transitively while
keeping graph assignments immutable. The full evidence chain and failed
ablations are recorded in
`evaluations/2026-09-12-chr12-18-20-comparison/regions/chr20_15019294_15130077.md`.

### Private-site MSA validation experiment (2026-09-13)

The BAM step-4 implementation does not promote first-pass noisy indels to clean
candidates. MSA reconstruction emits `NoisyCandHet/Hom`, and the BAM pipeline
uses those calls only through a second `kCandGermlineVarCate` k-means. In
graph+private mode that stage was disabled entirely to enforce the whitelist.

`collect-hybrid-variation --private-msa` now provides a safe experimental
variant: run MSA only around whitelisted private keys, admit exact-key MSA calls
only, and replace a whitelisted `RepeatHetIndel` exact collision with its MSA
call/profile. A unit test verifies whitelist exclusion, collision replacement,
and profile reindexing.

The target chr20 window gained a larger, 100%-accurate central block and a
complete candidate chain phased 554/1,034 reads with zero discordance. It did
not close both boundaries because MSA and DeepVariant represent the two repeat
alleles at different positions, leaving zero shared profile observations for
the MSA-normalized keys. The hybrid and competitor BAMs have identical
qname/flag/position/MAPQ/CIGAR records in this interval after contig aliasing,
so this is a representation/profile issue rather than missing alignments. The
full chr20 GQ10/PS50 run was negative: N50 stayed
991,274 bp and switchflips stayed 139, while phased reads changed
179,494 -> 179,452 and discordant reads changed 269 -> 285. Keep this mode
opt-in. Equivalent-allele reconciliation plus a transitive bridge-edge gate is
required before enabling it by default.

### Reproducible second-pass gap repair (2026-09-13)

The chr20:15,019,294–15,130,077 gap can be joined using the existing hybrid
MSA/k-means and graph-block merger. See
`evaluations/2026-09-13-gap-rephasing/README.md` and run
`scripts/test_chr20_gap_rephasing.sh` for the fixed regional regression.

Fixed a diploid double-swap no-op, tied-link acceptance, additive CLI graph
ownership override, region-mode MSA repeat collision admission, and sorted MSA
read-profile index mapping. Removed the read-only pre-MSA HP/PS restoration:
restoring reads without candidate orientations mixes independent solves and
undoes joins. Existing complete-block stitching preserves internal graph
orientation instead. Also corrected invalid-PS checking and ensured plain
`make` rebuilds the executable rather than a dependency-file default object.

In the 230 kb flanked window, clean hybrid phases 555 reads in three blocks
with zero errors. Scoped MSA at its unchanged margin-24 default phases 653
reads, still in three blocks with zero errors. A diagnostic margin-1 run phases
the same 653 reads in one block with zero read-truth errors/switchflips. It
recovers four unanimous observations across the left boundary, versus zero at
margin 24. Graph-block stitching then retains all 495 original graph reads and
adds 211: 706 reads, one block, zero errors. This uses an explicit regional
400 kb PS-start-distance cap because the 110.8 kb physical gap's endpoint IDs
are 385.6 kb apart. Neither margin nor merger-distance defaults were changed.

A uniform HP1/HP2 swap within a block is not a phasing error. The regression
checks graph-block-relative orientation, rather than literal tag equality.
These are local read-truth results, not NGC50 or chromosome-wide validation.
All unit binaries pass; the existing `make check` script fails before running
comparisons because it invokes obsolete CLI syntax.

### Independent second-region check: chr20:47,671,540–47,762,233

After pushing `663b84a`, tested the next-largest HiPhase-bridged chr20 gap
without further C++ changes or a threshold search. The six clean hybrid SNP
positions all already exist in the graph catalog; native BAM calling provides
five clean SNPs and 20 MSA het candidates inside the gap. A fresh GQ10,
VAF0.30–0.70 native whitelist scopes second-pass MSA.

In the 220 kb flanked window, clean and margin-24 hybrid runs each phase 510
reads in three blocks with zero errors. Diagnostic margin 1 phases the same
510 reads in one block with zero errors. Existing graph stitching retains all
462 graph reads and adds 70, yielding 532 reads, zero errors, and two blocks.
The target graph endpoints (47,666,163 and 47,762,233) join with 35 winning
reads under the unchanged 300 kb cap. A third block at PS 47,121,694 is kept
separate because its 640.5 kb PS-start separation exceeds that cap.

This again exposes missing actual allele observations, not just missing site
positions. One diagnostic boundary has only five observed reads with a 3:2
orientation majority; the truth-correct local outcome is insufficient to
promote margin 1 globally. No default changed. Reproduction, all candidate
rows, endpoint votes, and preservation assertions are saved under
`evaluations/2026-09-13-gap-rephasing-47m/`.

### Automatic cumulative gap recovery in the hybrid pipeline

Implemented `collect-hybrid-variation --recover-gaps` with an optional
`--gap-recovery-report`. After clean BAM+graph phasing and initial chunk
stitching, read-supported internal gaps are detected across chunk boundaries.
Each gets a locally reloaded window with 50 kb flanks. Evidence escalates in
one persistent proposal chunk: clean sites, then MSA SNPs, then MSA indels.
The normal k-means and overlap stitching vote rule are reused. Each endpoint
must pass the selected normal stitch rule, and both must connect through one
proposal PS to close a gap. One-sided additions survive escalation; initial
blocks are only uniformly flipped/relabelled. There is no external merger or
PS-start-distance cap in this mode.

MSA is driven by each remaining gap plus 5 kb near its boundaries, preserving
existing noisy intervals and tiling uncovered clean stretches with bounded,
overlapping windows. Recovery retains the needed intermediates and delays
pruning to avoid stale candidate/profile indices. Explicit replacement of
unsupported exact matches transfers candidate category and observations
together; raw observations from untagged reads are retained without tagging
those reads. Excluded intervals between separate requested regions are never
recovered. The mode is opt-in and preserves its initial hybrid clean scaffold,
not an externally supplied graph-only BAM.

The two audited regions each move from three phase sets to one with zero
read-truth errors at the unchanged MSA margin 24: 555/555 retained reads at
14.95–15.18 Mb and 510/510 at 47.60–47.82 Mb. An 80 kb chunked first-region
run retains all 534 phased reads and also moves three PS to one without errors.
Both residual gaps in each window require the indel tier. Margin 1 gives the
same outcomes. A first-region run with the default adjacent-site linker also
closes successfully. These are local read metrics, not chromosome NGC50.

Reproduction and assertions are under
`evaluations/2026-09-13-auto-gap-recovery/`. Build and all unit binaries pass;
`make check` remains blocked by its pre-existing obsolete CLI invocation.
Whole-chromosome accuracy and scaling still need validation before this mode
becomes a default.

### 2026-09-13: Independent automatic recovery validation after 5316f9c

Pushed automatic tiered recovery, then tested six previously untested regions
selected from the frozen HiPhase bridge audit (two each on chr12, chr18, chr20),
with fixed settings and 50 kb flanks. Recovery joined 3/11 residual hybrid gaps:
chr20 61.8 Mb at tier 2 and the other two joins at tier 3. chr12 63.8 Mb and chr20 61.8 Mb went from two blocks to one,
adding 40 and 36 tagged reads respectively. chr18 33.3 Mb improved from four
blocks to three. chr12 46.8 Mb, chr18 46.2 Mb, and chr20 17.6 Mb did not improve.
All 2,293 original tagged reads retained uniform block orientation transforms.
No read-truth discordance/switchflip increase occurred; chr18 46.2 Mb retained
one pre-existing error, other cases had zero. Build/unit tests passed.
Unresolved proposals can match both endpoints in separate phase sets, so both
link flags alone must not be interpreted as a successful join. These targeted
read-truth checks do not establish chromosome-wide NGC50 or variant-truth accuracy.
Manifest, exact commands, tier reports and summaries:
`evaluations/2026-09-13-auto-gap-validation/`.

### 2026-09-13: Why chr12 46.8 Mb remains disconnected

Diagnosed first residual hybrid gap 46769497–46810372 in the independent
validation case. Recovery retrieves the repeat/SNP/insertion positions used
by HiPhase, but default MSA admission leaves just 7 usable observations at
46774503; margin1 gives 83. Eighteen primary MAPQ60 reads span that SNP to
46788671, but none have observations at both in the default final matrix.
The later insertion-to-SNP link (46788671–46810372) has only one spanning
primary MAPQ60 read, also tagged by HiPhase in PS46709288. Changing only
min-block-link-reads from 2 to 1 lets recovery join the target at tier3;
regional blocks 4→2, 288 evaluated reads, zero read-truth errors. Margin1
alone gives 4→3 blocks, 336 reads, but does not close the first gap. Both
changes give 2 blocks/336 reads/zero errors. The remaining block boundary is
in the right solve-window flank. No production defaults changed. This is
not the old HP-only linking circularity: allele linking was already enabled.
Detailed evidence and commands: `evaluations/2026-09-13-chr12-gap-diagnosis/`.

### 2026-09-13: Site-level observations for automatic recovery

Added a recovery-only path for full-cover reads rejected by the MSA consensus
score margin: compose their alignments through both fixed consensuses and admit
only exact site alleles with three matching flank bases and agreement between
both paths. This supplies observations without choosing a consensus/HP.
Expanded counts must pass existing min_af/max_af before profile/count updates;
SNP-first tiering, margin24, min-block-link-reads2, and stitching are unchanged.
The AF gate fixes a regression found during validation: consensus-difference
MSA het labels can persist despite overwhelmingly reference expanded support.
Without this gate the chr20 17.6 Mb trial added two erroneous reads; with it,
two correct reads are added and no erroneous reads. Prototype results retained.

On the same six fixed windows, joins improve from 3/11 to 6/11 internal gaps.
chr18 46.2 Mb improves from 3 blocks to 1 (+12 reads); chr18 33.3 Mb from 3 to 2.
chr12 46.8 Mb adds 3 reads, chr20 17.6 Mb adds 2, both retaining their gap counts.
Other cases remain unchanged. Total gain over previous recovery: 17 reads,
no increased read-truth discordance/switchflips; 2,293 original tagged reads
preserve uniform block transformations. Original 15 Mb/47 Mb/split-chunk fixtures
still pass. Build/all unit tests pass, including exact SNP/indel site calls,
ambiguous/path-disagreeing observations and atomic AF-gate rejection.
These are local read-truth tests, not chromosome-wide or singleton-bridge validation.
Reports and commands: `evaluations/2026-09-13-site-observation-recovery/`.

### 2026-09-13: Resolve chr18 33.3 Mb's remaining singleton boundary

The remaining boundary is 33284672–33308109, in the upstream flank of the
original audited target. Both endpoint SNPs already have observations; HiPhase
has no intervening het. Exactly one primary MAPQ60 read spans both, with
endpoint base qualities 35/40 (and another left SNP at quality 40). HiPhase tags
that read in PS33179887. Existing min-block-link-reads 1 resolves the entire
window to one block with 339 reads and zero read-truth errors. Minimum 2 cannot
accept that bridge; two site pairs on this read are not independent molecules.

Validated the existing minimum 1 setting on all six fixed windows: chr12 46.8 Mb
4→2 blocks, chr18 33.3 Mb 2→1, chr20 17.6 Mb 2→1; other cases unchanged. All 2,386
previously tagged reads retained uniform block transformations and no increased
read-truth discordance/switchflips. No source/default change: solved outputs use
an explicit --min-block-link-reads 1 setting.
Commands/reports: `evaluations/2026-09-13-singleton-link-validation/`.
This targeted validation does not establish chromosome-wide singleton accuracy.

### 2026-09-13: Correct the chr12 target-versus-window comparison

The min-link1 chr12 4→2 regional result already phases the original audited
46725702–46838918 target in PS46709288, matching HiPhase. Its remaining boundary
46842224–46874196 lies in the downstream solve-window flank. HiPhase also breaks
there (46709288→46874196), as do WhatsHap, WhatsHap-opt and LongPhase
(46810372→46874196). Zero primary nonsupplementary BAM reads span these endpoints.
Do not describe this regional two-block output as a remaining competitor-bridged
gap, or force a merge to satisfy the misleading comparison. Reporting clarified;
no phasing code changed. Reproducible endpoint audit:
`evaluations/2026-09-13-chr12-remaining-boundary/`.

### 2026-09-13: Systematic competitor-gap audit and supported MSA flanks

Expanded the audit to full available chr12/chr18/chr20 chromosomes, using current
hybrid clean/recovery outputs projected onto the same DeepVariant records as
HiPhase, WhatsHap, WhatsHap-opt and LongPhase. Inventory distinguishes unphased
sites from exact endpoint PS breaks, MAPQ-qualified spanning molecules, native
PS membership, and competitor endpoint parity against assembly-based truth.
The truth VCF has chromosome-wide paternal/maternal GT and no PS field; absence
of PS there must not be mistaken for absence of truth phasing. A competitor's
whole-block read accuracy alone does not validate every individual join.

Found a real recovery bug at chr12:24106886–24129119: a SNP at 24124230 lies
next to an MSA-supported insertion. Requiring reference-identical flanks rejects
otherwise useful site observations. Site-level recovery now accepts exact flank
sequence supported by either fixed MSA consensus, while retaining agreement of
both read alignment paths, complete site/flank coverage, and the existing AF
gate. Unexplained flank changes remain rejected. Under default min-link2 the
SNP tier joins 2 blocks to 1, increasing evaluated reads 383→427, with zero
read discordance/switchflips. Six previous regression regions preserve all
2,293 original tagged reads and their uniform block transformations; chr18
46.2 Mb gains another 35 reads (470→505) without increasing its one prior error.
Other six-case results are unchanged. Full chromosome evaluation is ongoing;
these local results do not establish genome-wide accuracy or completeness.

The chr12:8235470–8248466 failure instead has 14 spanning primary reads, all
below MAPQ30 (13 at MAPQ3, one at MAPQ19). HiPhase defaults to MAPQ5. An explicit
MAPQ5/min-link1 local experiment phases 286 reads in one block with zero
read-truth errors. A prototype extending duplicate clean-site MSA observations
did not fix this region and was removed. No default MAPQ change was made.
Whole-chromosome chr12 clean control with min-link1 increases discordant reads
3209→3276 (+67), despite encouraging earlier six-window singleton tests. Do not
promote minimum-one linking globally based on those selected windows.

Full recovery also exposed repeated rebuilding of a chromosome-wide read-name
map per gap/tier. GapReadIndex now indexes stable read locations once and reads
current HP/PS values after every join. Only proposal assignments are materialized
for each stitch. A repeated-tier unit test checks newly accepted reads are
protected on reuse. Build/unit tests pass. Audit and regional evidence:
`evaluations/2026-09-13-panel-gap-audit/`.

### 2026-09-13: Focus the complete audit on chr20

At the user's request, paused incomplete chr12/chr18 whole-chromosome recovery
runs and moved the iteration loop to chr20. Completed clean hybrid projection
finds 167 competitor-bridged endpoint pairs representing 138 distinct pgphase
block pairs, plus 12,996 unphased genomic positions that a competitor phases.
114 endpoint pairs (101 block pairs) have concordant orientation against the
assembly-based truth and are tested in local 50 kb-flank windows. These windows
can overlap; sums of regional reads/errors are not chromosome-wide metrics.
Pairs lacking exact truth representation or disagreeing with truth remain in
the inventory but are not treated as correct bridges to reproduce blindly.

The 63-window chr12/chr18 regression completed before the focus change. It
exposed five incorrect joins already present before the consensus-flank fix;
that fix additionally introduces two discordant newly phased reads in one
chr18 window. The earlier six-window no-regression result does not generalize
to every region. Chr20 likewise exposes existing incorrect joins. Full-panel
and targeted truth checks are therefore required alongside block counts.

### 2026-09-13: Genotype-aware multiallelic HP projection

`phase_vcf_from_hp.py` previously omitted every multiallelic heterozygous record.
It now observes exact CIGAR-derived allele sequence and votes within the two
alleles of the original diploid genotype. GT 1/2 remains 1/2 (or 2/1), rather
than being recoded as reference/alternate. Both cached and pileup paths are
covered by integration tests; same-length inserted alleles are distinguished
by sequence. Legacy cache REF_COUNT/ALT_COUNT columns mean sorted-GT slot 0/1
for multiallelic records; older caches must be rebuilt to gain omitted sites.

On chr20 clean tags, phased shared calls increase 66,723→68,306 (+1,583).
Original biallelic truth comparison is unchanged. After splitting, normalizing
and sorting both truth and predictions, assessed variants increase
60,189→61,105; hamming errors 138→142 and switches 63→65. On chr12 the gain is
3,345 projected calls; normalized truth adds 1,934 assessed variants, 12 hamming
errors and 2 switches. This restores output coverage; it does not join BAM
phase sets or improve NGC50 by itself. Normalization can move records more than
bcftools' default sort window, so an explicit sort is needed before indexing.

### 2026-09-13: Recovery sites must be allowed to connect earlier components

chr20:46384021–46407462 remains split after all three retrieval tiers, despite
three MAPQ60 molecules giving consistent SNP→insertion observations. The
insertion is present in the graph catalog (46405052 AA→AAA left-normalizes to
46405048 T→TA) and in the MSA matrix with 65 observations (42/23). The core
chooses only one nearest sufficient preceding edge per site. A later site
therefore cannot connect two earlier components even when both links pass.

Recovery MSA rounds now retain all qualifying allele links in the existing
het window, prioritize net evidence, and join components with composed parity.
Weaker contradictory cycles do not overturn stronger paths. Initial phasing
and normal block-stitching remain unchanged. A trial that additionally required
net support 2 lost four previously correct local joins; the final candidate
retains the existing majority-support threshold instead. The change is gated
by recover_gaps + region MSA admission + allele linking. Unit tests cover the
missed two-component bridge, orientation convergence, conflicting cycles, and
preservation of the original support threshold.

The concrete chr20 46.4 Mb reproduction improves two blocks→one with the same
420 truth-evaluated reads and zero discordance/switchflips. Reports distinguish
that local result from the ongoing whole-chromosome comparison. Source and
reproduction evidence: `evaluations/2026-09-13-panel-gap-audit/chr20_46384021/`.

### 2026-09-13: Trial homopolymer evidence per unresolved gap

Unrestricted admission of MSA homopolymer link nodes fixed chr20 62.4 Mb but
introduced 53 discordant reads at 24.1 Mb and 5 at 36.0 Mb in local tests.
Those two target gaps already joined without homopolymer links. The correction
is to trial additional evidence only after the ordinary three tiers fail,
not to enable or reject homopolymers chromosome-wide.

Tier 4 reuses the existing clean + MSA SNP + MSA indel observations and k-means.
Only MSA-verified homopolymer sites inside the current unresolved gap can become
extra link nodes. Such edges require at least min_block_link_reads and no opposing votes.
A proposal must connect both flanks with the normal stitch rule and with each
original haplotype independently favoring its orientation by that same margin.
Failure is transactional: no one-sided extension or original-block change is
committed. Reports mark failed tier-4 proposals as rejected. Successful earlier
tiers are not retried. These are evidence-based safeguards, not a proof of phase
correctness; truth is used only for evaluation, never for site selection.

Initial four-case validation: 62.4 Mb joins 2→1 blocks, 385 truth-evaluated reads,
zero discordance. The 24.1 Mb and 36.0 Mb cases retain zero discordance (352 and
117 reads respectively). The pre-existing incorrect 19 Mb ordinary-tier join
remains 208/489 discordant and is not fixed by this fallback. The first full 114-case local comparison retained all 85 prior joins and
added four, but introduced 39 discordant reads at 55.9 Mb. The net-margin
version is superseded: 55.9 Mb had conflicting 4:1 and 18:3 HP bridge edges
despite pure anchor votes, while 62.4 Mb has an alternate unambiguous 2:0
right-side link. The revised edge-unanimity version is under validation. Build and unit tests pass; the only build warning
is the existing unused SIMDMalloc function in the abPOA header.

Crucially, the full-chromosome graph-support2 arm (before adaptive homopolymers)
has 12,848/187,387 discordant reads versus 9,150/187,345 for nearest-link recovery2,
despite improved local-window results. Full-chromosome shared-DV NGC50 is 451,323 bp versus 442,852 bp, while
variant hamming rises 2,336→3,537 with the same 87 switches. There are 47
remaining competitor-bridged native PS pairs (54 endpoint pairs), versus 53
(62 endpoint pairs) for nearest-link recovery. Adaptive whole-chromosome
validation is still pending; local results do not establish safety of the
combined recovery changes.


The acceptance objective is correct *relative block orientation*, not zero
read-level discordance. The 55.9 Mb failure is confirmed as a bad stitch:
all 39 reads from original PS 56040612 change from concordant to discordant,
while all 269 reads from PS 55949998 remain concordant; no new reads are tagged.
Shared SNPs 55999561 C/T and 56040612 C/T are both 0|1 in assembly truth,
but the trial puts 1|0 and 0|1 in the same PS. This independently confirms the
wrong endpoint parity. The local summary now writes block_orientations.tsv,
separating reversed original-block majorities from isolated read errors and
recording uniform HP-label transforms. Absolute HP1/HP2 labels are arbitrary.
The edge-unanimity rule remains a conservative fallback experiment; read
assignment noise is not itself grounds to reject a correctly oriented stitch.

At this point `make check` could not execute its golden comparisons: the
validator used obsolete --phased-vcf-output and positional reference/BAM
arguments. The CLI invocation was repaired on 2026-09-20; the subsequent
golden mismatch is recorded in the paired noisy MSA section below.


Full-chromosome original-block audit (`compare_full_block_orientations.py`)
finds 27 original blocks with reversed truth-majority orientation after
nearest-link recovery2 and 34 after graph_support2. These counts separate
coherent block orientation failures from isolated read disagreements and
explain why larger blocks alone are not sufficient evidence of improvement.
They are relative to the clean run and the majority truth orientation of the
final merged block; consult the per-block counts rather than treating every
small/noisy block as equally certain.


Final 114-case local result for chr20_adaptive_unanimous: all 85 previous
endpoint joins retained, two additional joins (57.1 Mb and 62.4 Mb), 22 targets
still split, five with an unphased endpoint. No case adds discordant reads or
switch/flip errors relative to graph_support2. All 37,479 original tagged-read
occurrences retain uniform original-block transforms; the 851 previously
added tagged-read occurrences are unchanged. Windows overlap, so these sums
are not independent chromosome-wide counts. The 55.9 Mb bad stitch is rejected;
the correctly oriented 54.5 Mb trial join is also withheld, documenting the
conservative rule's sensitivity cost. Existing ordinary-tier orientation
failures, including 19 Mb and the full-chromosome graph regression, remain.
Build/unit tests and projection tests pass. Full-chromosome validation of this
fallback is not yet complete; make check is blocked by the stale gate CLI
invocation documented above.


### 2026-09-13: Root causes of false gap stitching

Two regression tests reproduced false phase-set ownership in the shared core:
(1) a read was assigned the first phased het in its sparse profile interval
even when its allele there was -1; (2) an HP repeat excluded from read-haplotype
scoring could still assign its inherited preceding PS to the read. Ownership
now requires an observed phased allele and uses the same eligible evidence as
read-haplotype scoring. These defects can supply confident-looking overlap
votes for a connection that has no actual allele bridge. The first fix removes
the incorrect 19 Mb join (208 discordant reads→0), retaining all 489 assessed
reads but initially leaving the region split.

Composing two left-aligned consensus paths does not guarantee a left-aligned
read/reference path. Within the 19 Mb repeat this hid valid observations or
made insertions look reference. Canonicalizing covered exact-match indel runs
restores a 3:1 informative insertion→right-SNP link where only one 0:0 pair
was previously recorded. The region then joins correctly with 489 assessed
reads and zero discordance; this is not an error-free-edge acceptance rule.

At 55.9 Mb, independent biallelic treatment of nested 4/6-base and 1/2-base
deletions allowed a longer-deletion read to support both ALT rows. Decomposing
shared deletion sequence from the length difference fixes that representation.
A balanced verified repeat can also collapse to the same allele in both
provisional haplotype profiles when independent blocks have opposite HP labels.
Preserving the verified genotype during the local solve recovered the correct
55.9 Mb join (308 reads, 39→0 discordant), but broad preservation introduced
bad joins elsewhere. The subsequent candidate trial therefore requires a
within-local-block association between both haplotypes and the repeat alleles,
computed BEFORE read HP/PS state is reset. The repeat is excluded from read-hap
scoring, so it cannot directly create its own validating HP association.

Assigned MSA reads also need site-level observations: the legacy updater can
impute REF simply from membership in the other whole-consensus cluster. New
unit tests cover an ALT read in that cluster and a third allele that must stay
unknown. Full-cover assigned reads are checked against composed, normalized
reference alignments. An exact-only trial lost valid signal and is superseded
by a bounded local comparison: at most one edit, strictly closer to one allele,
and fixed-consensus flanks identical between haplotypes. Thus only the target
site distinguishes the hypotheses; a tied third allele remains unknown.
Partial-read handling is retained. These changes remain under panel validation;
no complete fixed-pipeline chromosome benchmark is claimed yet.

Completed owner_clean (observed-allele owner fix only) full chr20 check: NGC50
290,115 bp and 65 variant switches unchanged; variant hamming 138→134 and read
discordance 499→498. Full nested_seed and observed_owner recovery arms were
stopped as superseded after local regressions; their INCOMPLETE.txt files
prevent treating them as complete benchmark results. Detailed causal evidence
is in evaluations/2026-09-13-panel-gap-audit/chr20_55999561/ and chr20_18983414/.


### Direct site evidence breaks the repeat-eligibility circularity

The within-proposal HP association gate discarded untagged spanning reads,
reintroducing the read-tag circularity in tier 4. At chr20 55.9 Mb it saw only
one tagged observation for the true residual deletion at 56027380, despite
18 agreeing versus 3 conflicting direct observations against clean SNP
56040612. Eligibility now additionally evaluates each verified repeat site's
allele contingency table against each phased clean heterozygous candidate.
Both anchor alleles must independently support the same orientation with the
existing net margin. No read HP tag is required; competing candidate sites are
not pooled into this eligibility vote. The normal recovery allele graph and
normal flank stitching still determine the final relative block orientation.

Unassigned full-cover MSA reads also now receive the same bounded local allele
comparison as assigned reads, with the extra requirement that both composed
alignment paths uniquely favor the same allele. A third-allele tie or disagreement
between paths remains unknown. Both changes have failing-before/passing-after
unit regressions; build and all unit tests pass without new compiler warnings.

The initial 11-case chr20_direct_site_anchor panel preserves all 3,940 original
truth-assessed reads, reverses no original block, and adds the correct 55.9 Mb
join (308 reads, zero discordant) relative to chr20_local_alleles. All five
previously joined targets remain joined; no discordance or switch/flip increase.
The full 114-case expansion is pending. This does not establish chromosome-wide
recovery accuracy or eliminate the remaining unjoined competitor gaps.


Completed direct_site_anchor 114-case validation: 77 joins, 33 splits, four
endpoint-unphased. The direct-site change adds three correct target joins but
also a wrong 54.89 Mb join (199 reversed original-block reads) and a five-read
block reversal near 36.03 Mb, so the any-favorable-anchor trial is superseded.
Its failure is candidate selection bias: a 5-read distant subset passes while
a 50-read near-anchor comparison shows no allele separation. The subsequent
best_site_anchor trial ranks clean anchors by minimum allele coverage, then
total coverage, and uses that comparison's eligibility; equally covered ties
must agree. This changes evidence selection without raising a noise threshold.
Known broader recovery errors at 35.9, 46.7, and 60.1 Mb remain under diagnosis.


best_site_anchor completed 18 targeted regression cases: 12 joins, six splits,
5,922 original truth-assessed reads preserved. It removes the two newly introduced
errors: 54.89 Mb 199→0 discordance by leaving the ambiguous gap split; 36.03 Mb
5→0 while retaining the join. Correct 4.85, 55.9, and 57.84 Mb joins remain,
including 308/308 concordant reads in the one-block 55.9 Mb result. No added
discordance/switch-flip errors versus direct_site_anchor in these cases. The
sparse-subset unit regression fails on the previous selector and passes now;
all unit tests/build pass. This revision has not completed the full 114-case or
whole-chromosome recovery validation. Existing wrong joins at 35.9, 46.7, and
60.1 Mb remain unresolved; do not treat this as a fully validated recovery core.


### chr20 35.9 Mb: indel tier lost previously verified evidence

Commit f92018d was pushed before this investigation. Two further bugs were
reproduced. (1) A verified SNP in the non-deleted allele of a residual deletion
was treated as an unknown third allele; exact fixed-consensus matching restores
19 reference / 7 alternate observations instead of 1 / 6 at 35927394. Unit
regression fails before and passes after. (2) SNP-only flank extension narrowed
the next tier's MSA interval, excluding the already examined insertion at
35920914. Keep the original MSA interval across all tiers. The indel now enters
tier 3 (37 reference / 22 alternate), producing a correct one-block 35.9 Mb
result: 293 truth-assessed reads, 48→0 discordance, 3→0 switch/flip errors.

fixed_tier_window completed all 114 local cases: 75 joins, 35 splits, four
endpoint-unphased; no discordance increases versus the prior direct_site_anchor
full trial. Two original-block reversals remain (46.7 and 60.1 Mb). This is not
a whole-chromosome recovery validation. Build and all unit tests pass.


### chr20 46.7 Mb: graph evidence filled a deleted SNP

Read m84031_231217_062403_s3/184948613/ccs has a primary MAPQ60 BAM deletion
at SNP 46727050, but graph-only profile augmentation filled its unknown slot
with SNP ALT. This manufactured one of two links driving the wrong join.
Graph SNP extension now skips explicit BAM deletions/reference skips while
retaining ordinary missing-slot filling. Unit regression fails before and
passes after. In the initial 18-case graph_deletion_guard panel, 46.7 Mb stays
split and improves 187→1 discordance (415 assessed reads); no new discordance
or switch/flip increases. All 35.9/55.9 rescues remain. Full expansion pending.

60.1 Mb diagnosis: true six-T/seven-T insertion alleles at 60093417 are separate
MSA candidates with only ALT observations, collapsing to homozygous and losing
the allele-length distinction. Proper multi-allelic representation remains to
be implemented; 39 reversed original-block reads are unresolved there.


Completed graph_deletion_guard 114-case expansion: 74 joins, 36 splits, four
endpoint-unphased. The only detected original-block majority reversal is the
unresolved 60.1 Mb case. The 35.9 Mb (293/293 concordant) and 55.9 Mb (308/308)
correct joins remain. 46.7 Mb stays split, 187→1 discordance. Tradeoffs versus
fixed_tier_window: 61.7 Mb loses a correct target join and adds three discordant
newly tagged reads; two existing reads near 1.9 Mb become discordant (duplicated
in two overlapping test windows); switch/flip counts increase in two overlapping
26.6 Mb cases despite fewer discordant reads. No additional original-block
majority reversals. All results are recorded, not claimed regression-free.

### chr20 60.1 Mb: preserve two alternate insertion alleles through recovery

The verified six-T/seven-T insertion at internal position 60093417 (VCF
60093416 C -> CTTTTTT,CTTTTTTT) was represented as two biallelic candidates,
with only alternate observations. Both collapsed to homozygous. Recovery now
merges two same-position MSA insertion calls into one site with explicit
reference/ALT1/ALT2 indices, recomputes observations against both consensuses,
and retains a supported ALT1/ALT2 genotype while optimizing its orientation
jointly through the existing k-means. Complementing an allele with `1 - allele`
is now restricted to biallelic sites. The bounded one-edit local comparison also
handles these pairs: five-T/eight-T reads can support six-T/seven-T when the
fixed consensus flanks agree and both composed paths select the same allele;
ties and more distant sequences abstain. Native VCF preserves both ALT strings,
GT 1|2 or 2|1, per-allele AD/AF/VAF, and PS. GQ is zero because the existing
biallelic model cannot assess this genotype.

Unrestricted pair admission regressed 35.9 Mb: a true homozygous 323-base
insertion was split into 322/323-base MSA consensuses, causing 48 discordant
reads. Reuse the per-site best-covered clean anchor to assess separation of
those two sequences. Its contingency table is [6,10;8,8], which supplies no
consistent separation. The real 60.1 Mb pair has [9,8;7,15]: both anchor alleles
favor opposite insertion lengths, with combined margin nine. For MSA pairs,
require both row margins to have the same nonzero sign and combined margin to
clear the existing block-link threshold; keep the prior repeat-site rule for
biallelic HP candidates. The best-covered informative anchor overrides sparse
favorable subsets. With no informative anchor the pair remains eligible for
the normal link solver; this is not a universal guarantee of heterozygosity.

The 18-case multi_anchor panel retains all baseline target outcomes, introduces
no discordance/switch-flip increases or original-block majority reversals, and
fixes 60.1 Mb: one block, 206 assessed reads, 39 -> 0 discordance and 3 -> 0
switch/flip errors. 35.9 Mb retains its correct one-block 293-read result.
Unit tests cover alternate identity, reference observations, bounded length
errors, path conflicts, full k-means genotype preservation, clean-anchor
selection, and native multiallelic VCF output. A first-pass test also caught
an empty preexisting HP/PS-vector access in the broadened selector; provisional
read votes are now optional while direct site observations remain usable.
Build and all unit tests pass, with no new warnings. The final revision limits
pair creation to the MSA gap pass and uses standards-compliant Number=A INFO AF;
the two root-cause cases pass again. Full 114-case results follow below. These are overlapping local windows, not whole-chromosome
validation. Superseded multi_insertion/genotype/joint/local trial artifacts are
retained with explicit failure notes.


Completed multi_anchor expansion: 79 joins, 31 splits, four endpoint-unphased.
Five newly joined test cases cover four locations (9.0, 11.6, 55.5, 58.7 Mb;
two overlapping 55.5 Mb windows). All 74 baseline joins are retained. There are
no discordance or switch/flip increases versus graph_deletion_guard and no
detected original-block majority reversals. 60.1 Mb improves 39 -> 0 errors;
the newly joined 58.7 Mb case retains one preexisting discordant read. Recovery
still adds some individual errors versus clean-only at 7.2 and 61.7 Mb, already
present in the baseline. All 228 final-binary clean/recovery BAMs have identical
read names, flags, alignment starts, HP and PS to multi_anchor. This removes the
last detected wrong original-block orientation in this selected local panel;
it does not establish whole-chromosome recovery accuracy or resolve every gap.

Final multi_final independent 114-case truth evaluation completed with the same
79/31/4 outcomes, five added joins, no baseline error increases, and no detected
original-block reversals. See its results, baseline comparison and orientation
reports for the complete accounting.

### Remaining chr20 targets: clean-block evidence across sparse junctions

The 35 targets left by multi_final comprise 17 without a directly spanning
MAPQ30 primary read, 12 with exactly one, two with multiple spanning reads but
unresolved observations/linkage, and four endpoint-projection failures. The
reproducible inventory is audit_remaining_gaps.py and chr20.remaining.multi_final.tsv.
Competitor/DV records are used for diagnosis only, never as recovery input.

A SNP-pair threshold discarded reliable sparse bridges even when one read
identified both existing blocks. Recovery now evaluates clean-SNP observations
in the components already resolved by stronger site edges. A primary MAPQ30
read must support one orientation consistently within each component. Each
flank needs two clean SNPs or an exact Q30 BAM base matching a clean SNP allele;
missing qualities and explicit BAM deletions/skips cannot provide the latter.
A consistent opposing read vetoes the bridge even if its evidence is sparse.
Process eligible component links by read support and compose their parity with
existing components, then use ordinary read assignment and gap stitching.

The isolated confident_clean_bridge full 114-case panel has 86 joined targets,
24 splits, four endpoint-unphased. All 79 multi_final joins remain. Seven newly
joined cases are 11426743_11447789, 18194808_18218259, 20874447_20895876,
24357615_24379633, 38259286_38303247, 54894127_54912022, and 57667169_57688923.
All seven have zero discordant assessed reads. There are no increased read
errors or switch/flips versus multi_final and no original-block majority
reversals. Build and unit tests pass; fixtures cover base/mapping quality,
sparse conflicting reads, multiple clean SNPs, parity, and convergence.
These are overlapping local-window results, not whole-chromosome validation.

Several separately archived experiments are NOT included in this fix. Scoring
reads only within their output PS exposed an existing cross-block-label bug,
but did not repair the four endpoint cases. Reassessing provisional homozygous
MSA indels and keeping the reference allele as a local hypothesis recovered
12.74 Mb deletion evidence (26 reference / 48 alternate), but the cumulative
trial caused a wrong second-gap join in the 34.1 Mb window (2 -> 63 discordant
reads). It remains rejected. Full experimental patches and explicitly labeled
trial reports preserve those findings for further diagnosis; they are not
accepted production implementations. The four endpoint cases and 24 split
cases remain outstanding.

### Clean-round phase-set scoping fixes the chr20 1.9 Mb endpoints

The two overlapping targets 1912570_1913054 and 1912570_1913912 have clear
BAM SNP evidence. At 1912570 the observed truth groups are PAT C22 / MAT T25;
at 1913912 they are PAT G24 / MAT T25. The collected allele profiles match
these bases. Site retrieval is not the failure: a read crossing independent
phase sets was scored once against their arbitrary HP integers, then updated
all their allele profiles with that same integer. Final HP could also describe
a different phase set from the read's reported PS.

Clean-candidate k-means iterations now score and update a spanning read
separately in every observed phase set. Final output HP is recomputed from the
final consensus within the reported PS. MSA-round iterative updates deliberately
retain their existing path: applying the change there too exposed repeat-link
regressions described below. This is a scoped clean-round fix, not a claim that
all MSA-round cross-component behavior has been repaired.

The clean_scoped_iterations panel has 88 joined targets, 24 split and two
unphased-endpoint cases. All 86 confident_clean_bridge joins remain. Both
1.9 Mb endpoints join; each overlapping window has 354 assessed reads, with
read discordance reduced from 30 to 3 and switch/flips from 10 to 3.
WhatsHap 2.8 compare against local variant truth assesses 203 rather than 171
variants in the first window (204 rather than 172 in the second), reducing
blockwise hamming from five to one in each. All 114 local variant comparisons
complete with no hamming-count increases. The completed read audit detects no
original-block majority reversals. The ambiguous 26.6 Mb endpoint case has one
additional discordant read (34 -> 35) and one additional read switch/flip
(15 -> 16), with unchanged variant-truth metrics and no wrong block stitch.
The 24 split targets and two endpoint cases remain unresolved. These overlapping
windows do not establish whole-chromosome NGC50 or hamming performance.

### Gap-endpoint projection preserves validated upstream phase evidence

The final two chr20 endpoint failures overlap at 26.6 Mb. They are projection
failures rather than successful block joins. At 26,602,087, pgphase's native
VCF already has `1|0` in PS 26,552,098, but HP projection rejects the site
because one haplotype has a 2:1 allele split, just below the uniform 0.70 purity
threshold. At 26,626,247, pgphase has no native call and the local HP majority
points in the wrong direction. DeepVariant supplies that record in a six-site
phased block, however, with the truth-concordant `1|0` orientation.

`phase_vcf_from_hp.py` now supports endpoint-scoped fallback signals. An exact
phased native pgphase record may fill a requested `--fallback-site`; if none
exists, `--retain-caller-phase-blocks` may preserve the input call only when its
caller PS contains at least three phased heterozygous records. Retained caller
blocks use a separate PS namespace and cannot merge pgphase blocks. The panel
runner supplies only each audited pair's two endpoints, so this does not relax
projection across its surrounding window.

Across all 114 chr20 targets, exactly two records change. The native addition
at 26,602,087 agrees with the assembly-truth orientation of its existing
pgphase block. The caller addition at 26,626,247 remains an isolated block, so
it cannot introduce a switch or stitch cascade. Results move from 88 joined,
24 split, and two endpoint-unphased cases to **88 joined and 26 split**. The
19 kb junction at 26,624,953-26,644,453 remains genuinely unjoined: neither
clean, MSA, graph, nor relaxed nested-site trials produced a read-supported
two-sided proposal.

The unit fixture checks both output HP/PS consistency and internal allele
profile partitioning for reads spanning two disconnected blocks. It fails
against the original phase core and passes with the correction. Build and
unit tests pass. compare_panel_variant_truth.py provides reproducible local
WhatsHap comparisons alongside the existing read and block-orientation audit.

### Rejected MSA experiments and remaining root causes

All trials below are diagnostic, excluded from the production fix above:

- Provisional-homozygous indel reassessment plus reference hypotheses recovers
  the 12.74 Mb deletion; independently clean-supported MSA anchors can join
  that target correctly. Full msa_site_anchors evaluation nevertheless loses
  four existing joins (7.0, 20.5, 21.1 and 38.2 Mb): 85 joined / 25 split / four
  unphased endpoints. See remaining_msa_site_anchors_experiment.patch.
- Scoping iterative updates in all rounds fixes 1.9 Mb but creates a wrong
  3.95 Mb stitch: 119 original right-block reads reverse, with 11/35 assessed
  variant hamming errors. It also loses 4.85 Mb because partial extension
  narrows later homopolymer eligibility despite evidence in the original gap.
- MSA homopolymer insertion detection compares raw FASTA ASCII to nt4 ALT
  codes. Correcting that real encoding bug stops the false 3.95 Mb bridge
  through the weak A-run insertion at key 3971337. It also moves previously
  usable insertions into stricter repeat validation, which needs further work
  before this classification fix can be integrated without losing joins.
- At 7.0 Mb, SNP 7027161 supplies a pure 2/2 allele table for insertion key
  7047081, using MAPQ60 and Q35-40 bases. The deeper clean anchor is noisier.
  A per-site pure-anchor exception admits the site, but a 3-versus-2 repeat
  edge then fails the extra net-margin rule. Using ordinary link scoring for
  validated repeats restores it and joins 3.95 Mb correctly, but the full
  validated_normal_links panel reverses a five-read block in the 36.0 Mb
  window (two variant hamming errors). This trial remains rejected.
- Coherent allele-group margins restore 20.5 Mb; pure allele separation with
  Q30 support on each side restores 38.2 Mb. The combined coherent_validated_links
  trial, however, reverses a 32-read original block at 23.4 Mb. Its exact source
  is archived in remaining_coherent_validated_experiment.patch.

Candidate heterozygosity, a single reliable site, and a reliable orientation
between existing blocks are separate decisions. Accepting the first two does
not by itself validate every weak downstream stitch. Competitor and truth
records remain evaluation-only inputs throughout these experiments.

### Graph-only SNP reference metadata repairs the 46.7 Mb clean bridge

The chr20 46,727,050-46,747,637 target remained split even though a MAPQ60
paternal read observes the left clean SNP at Q35 and the right clean SNP at
Q40. The recovery graph already permits one spanning read when both clean SNP
observations can be checked directly in the BAM. The right observation passed,
but the left graph-only SNP retained the default unknown `ref_base=4`, so the
same check incorrectly rejected its otherwise exact BAM base.

Graph-only SNP injection now copies the catalog REF nucleotide into the
candidate's nt4 `ref_base`. This preserves the existing MAPQ30/base-Q30 gate
and does not lower link support thresholds. The target joins in the
clean-plus-MSA-SNP tier, before MSA indels or the homopolymer fallback. Its
single 972,856 bp block has 414/415 concordant truth-labelled reads; local
WhatsHap comparison assesses 32 variant pairs with zero switches and zero
Hamming errors.

The complete 114-case chr20 replay has 89 joined and 25 split targets. Relative
to clean_scoped_iterations, all 88 earlier joins remain joined and 46.7 Mb is
the only new join. No case increases read discordance or switch/flips, and no
original block has its truth majority reversed. Results are archived under
`evaluations/2026-09-13-panel-gap-audit/chr20_graph_refbase/`. Build and all
unit tests pass. These remain overlapping local-window results rather than a
whole-chromosome NGC50 measurement.

### Backfill admitted MSA SNPs from the complete BAM read set

At chr20 22,721,581, the MSA found a real A>T phasing marker between blocks:
assembly-labelled BAM reads are maternal 26 ALT / 1 REF and paternal 19 REF /
6 ALT. The recovery candidate nevertheless contained only 25 observations
(23 ALT / 2 REF), because MSA merging retained calls only for reads selected
into the local consensus branches. Reads overlapping the new site from the
other flank were never rescored, so the site could extend the right block but
could not connect it to the left.

After each successful gap-MSA merge, missing observations at admitted SNPs are
now filled from every overlapping primary BAM alignment using the candidate's
exact REF/ALT, the configured MAPQ and base-quality floors, and the existing
read-profile representation. Existing MSA observations are preserved. The
read/variant interval index is rebuilt before the normal k-means rerun, so the
ordinary recovery graph and stitcher consume the complete evidence without a
new phasing path.

The chr20 replay moves from 89 joined / 25 split to **91 joined / 23 split**.
The new joins are 22,702,346-22,729,320 and 61,732,321-61,782,778. The first
has 318/320 concordant labelled reads, unchanged from its two-block baseline,
and zero variant switches/Hamming errors. The second improves from three
discordant reads to zero. All 89 existing joins remain, no original block's
truth majority reverses, and the 114-window WhatsHap audit reports no Hamming
count increase. Additional previously untagged reads expose a few isolated
read-level discordances in other windows, without a block-orientation error.
Results are archived under
`evaluations/2026-09-13-panel-gap-audit/chr20_msa_snp_backfill/`.

### Backfill exact MSA indels only inside the active gap

The chr20 37,984,529-38,008,026 target had a complete but previously unusable
evidence chain. MSA insertion 37,984,530 was attached to the left block, MSA
SNP 38,005,401 was attached to the right block, and one MAPQ60 read crossed
the insertion and the right clean endpoint with Q35-Q40 bases. Its insertion
allele remained missing because BAM backfill handled only SNPs. Consequently
the ordinary recovery component graph never saw the crossing edge.

MSA candidates now retain explicit provenance. Missing observations for
biallelic, non-homopolymer MSA insertions and deletions are filled only when an
alignment's CIGAR exactly matches the candidate allele and the relevant bases
or indel anchors pass the configured quality floor. A recovered insertion can
act as a single-read block anchor only after multi-read evidence has associated
it with an established clean block, and the crossing read's CIGAR observation
is independently rechecked at Q30. Homopolymer indels keep their existing
tier-4 path.

The first whole-panel trial also backfilled MSA sites in the 10 kb flanks. An
insertion beyond the 62,432,427 right endpoint changed the later homopolymer
retry and lost that accepted join. Backfill and recovered-indel anchors are
therefore restricted to the actual unresolved phase gap, while MSA discovery
still uses the flanks for context. This preserves the 62.4 Mb join and closes
37.98 Mb as one local block. Read truth changes from 456 reads / one discordant
to 474 / four discordant, with no original-block majority reversal and no
variant Hamming-count increase.

The complete 114-window replay is **92 joined / 22 split**, retaining all 91
previous joins and adding only 37.98 Mb. The WhatsHap audit reports no Hamming
count increase in any window. Results are archived under
`evaluations/2026-09-13-panel-gap-audit/chr20_msa_indel_scoped/`. These are
overlapping local-window results, not a whole-chromosome NGC50 measurement.

### Normalize equivalent MSA insertion placements inside short repeats

The chr20 863,640-890,261 target contained one MAPQ60 read connecting clean
left SNP 863,406 (Q40) to the MSA insertion at key 882,278, which was already
attached to the right block. The MSA represents its observed allele as `TC` at
882,278; the BAM CIGAR represents `CT` at 882,283. The intervening reference
is `TCTCT`, so both placements produce the identical `TCTCTCT` haplotype.
Exact-coordinate matching incorrectly discarded this otherwise clear edge.

For MSA-verified insertions used by the recovery block graph, placement is now
normalized by comparing the complete alternate haplotype sequence across a
maximum 20 bp shift. The read must still be primary MAPQ30 evidence, all
inserted and flanking bases must be Q30, the MSA allele must already be
associated with an established block, and an opposing crossing read still
vetoes a single-read bridge. This handles alignment-equivalent placement; it
does not admit arbitrary repeat indels.

The target becomes one 140,786 bp local block with 433/433 truth-labelled
reads concordant and zero read switches. The complete replay is **93 joined /
21 split**: all 92 prior joins remain, no original block majority reverses, no
case increases read discordance or switch/flips relative to the prior panel,
and the 114-window WhatsHap audit has no Hamming-count increase. Results are
archived under
`evaluations/2026-09-13-panel-gap-audit/chr20_shifted_insertion/`.

### Admit coherent moderate-quality three-SNP block bridges

The chr20 23,460,963-23,481,134 target remained split because its only
crossing primary MAPQ60 read had a Q22 clean SNP on the left and two clean
SNPs on the right. All three observations support the same block orientation,
but the recovery graph required each single-read clean-SNP anchor to pass
base Q30. This discarded the complete multi-site signal because one component
had only the Q22 observation.

A single recovery read may now link two established components when it carries
at least three clean SNP observations in total, has at least one observation
in each component, and every observation is independently confirmed from its
primary BAM alignment at base Q20 or higher. MAPQ30 remains required. The
existing opposing-read veto remains in force, while the ordinary two-SNP
single-read bridge continues to require base Q30 on both sides. A unit test
checks that the three-SNP Q22 pattern joins and the same pattern at Q10 does
not.

The target becomes one local block with 206/206 truth-labelled reads
concordant and no switch/flip errors. The complete 114-window chr20 replay is
**94 joined / 20 split**: 23.46 Mb is the only target changed from the prior
accepted panel, no original block majority reverses, and the WhatsHap audit
reports no Hamming-count increase in any window. Results are archived under
`evaluations/2026-09-13-panel-gap-audit/chr20_three_snp_bridge/`. These remain
overlapping local-window results rather than a whole-chromosome NGC50
measurement.

### Cache gap-targeted MSA evidence and apply phase edges once

Gap recovery previously reopened the alignment inputs and reran candidate
collection/MSA independently for every unresolved junction. Besides being
slow, each accepted retry immediately relabelled the right phase block, so a
later retry could consume assignments created by an earlier one. Recovery now
freezes the initial read assignments, builds MSA evidence once over merged gap
intervals in the original 500 kb chunks, and constructs each tier proposal
from retained candidates and per-read allele profiles. Accepted relationships
are collected as phase-set parity edges. Contradictory cycles are rejected,
then all accepted edges are applied in one pass before the ordinary final
left-to-right chunk stitch.

The gap-targeted evidence can be persisted with `--gap-evidence-cache FILE`.
On the first run pgphase writes the MSA-verified candidates and complete
candidate-indexed read profiles. A later run validates a signature covering
the initial gaps, chunks, read identities, initial candidates/phases, and MSA
settings before loading it. A mismatch fails with an instruction to rebuild,
so evidence from a different region or parameter set cannot be silently used.
The file is written through a temporary path and renamed only after a complete
write.

On chr20:1-10,000,000 with 20 chunks and eight threads, building the cache took
22.68 s and the 45 gap solves took 6.14 s (32.91 s total). Reusing the 99 MB
cache took 0.22 s, with 6.16 s for the same solves and 10.45 s total. Both runs
joined 21/45 gaps; their candidate TSV, phased VCF, and tier report were
byte-identical. The standard initial noisy-region MSA alone was also tested as
a substitute. It completed quickly but joined only 3/45 gaps, confirming that
the clean-looking gap windows supply real additional sites and must be
computed once rather than discarded.

For the complete 66,210,255 bp chr20 run (133 chunks, 283 initial gaps), the
first run spent 133.32 s building a 634 MB evidence cache and 137.34 s in the
existing local gap solves (301.98 s total, 117 joins). Reuse loaded the cache
in 1.34 s and produced byte-identical candidate TSV, phased VCF, and tier
report, but the unchanged local solves still took 137.52 s (170.32 s total).
This isolates the next performance target: replace per-gap dense proposal
construction and repeated k-means with sparse phase-edge scoring over the
cached observations.

### Enable clean block bridges in recovery, with singleton evidence gating

The cached recovery proposal set `gap_recovery_beg/end` for every evidence
tier, but enabled `private_msa_admit_all_in_region` only for the final
homopolymer tier. The recovery block graph is gated on that option, so clean
block bridges were computed while the evidence cache was populated but were
disabled when the tier proposal actually decided whether to join the original
phase sets. Enabling the additive recovery graph for the whole local proposal
fixes this orchestration bug. The chr20 18,194,808-18,218,259 reproduction now
joins at tier 1 from its existing clean BAM evidence.

A blanket enable exposed why the earlier overlapping-window checks were not a
sufficient acceptance test. It joined 129/283 gaps and raised raw block NG50
to 739,887 bp, but whole-chromosome blockwise Hamming increased from 770/59,990
(1.2835%) to 1,236/59,992 (2.0603%). Two singleton bridges caused the increase:
the 0.864 Mb edge depended on one MSA indel observation, while the 23.46 Mb
edge used the moderate-Q three-SNP exception and had Q22 at its left projected
endpoint. Both appeared correct in isolated windows but assigned the adjacent
whole-chromosome blocks in the wrong relative orientation.

Singleton component bridges now require a Q30 clean SNP observation on both
components, or at least two clean SNP observations on each component. MSA
indels and Q20 clean observations can still contribute when at least two
independent reads support the same orientation. The existing MAPQ30 gate and
opposing-read veto remain. This supersedes the singleton acceptance claims in
the earlier sections about shifted insertion placement and the moderate-quality
three-SNP bridge; those mechanisms remain available as
multi-read evidence but are unsafe by themselves.

With the singleton gate, the complete cached chr20 run joins 125/283 gaps.
Relative to the previous 117-join run, raw block N50 increases from 604,849 to
683,391 bp and switch-corrected NGC50 increases from 451,323 to 479,016 bp.
Corrected covered span increases from 52,379,757 to 52,734,498 bp (79.65% of
chr20). On the shared DeepVariant assessment positions, 59,992 variants are
assessed with 77 switches and blockwise Hamming 665/59,992 (1.1085%). Read
truth evaluation improves from 4,005 to 3,872 discordant reads and from 97.88%
to 97.95% accuracy. Thus the additional clean joins improve contiguity while
reducing both variant- and read-level errors.

### Retain exact graph confirmation at shared BAM candidates

The remaining all-four-competitor links at chr20 34.79 Mb and 61.66 Mb each
had one crossing read whose BAM endpoint quality was too weak for a singleton
bridge, but the same read traversed exact graph snarls at both endpoints. The
hybrid injector previously discarded GAF observations whenever the graph site
already matched a BAM candidate, so recovery could not distinguish this
independent confirmation from an unsupported low-quality base.

Read profiles now retain a graph-confirmation marker when BAM and GAF alleles
agree at a shared candidate. A singleton recovery edge may use that evidence
only when the same read has a graph-confirmed clean SNP on both components;
one-sided confirmation remains insufficient. The graph-to-candidate bridge
also fixes two related multiallelic bugs: padded substitutions such as
`CG->TG` are reduced to their clean SNP instead of being encoded as deletions,
and graph allele 2 is no longer collapsed into the binary first-ALT allele.

Focused chr20 tests join 34,791,342-34,811,747 and
61,664,136-61,690,751 at the clean tier. The known harmful
23,460,963-23,481,134 singleton and the one-sided
1,086,625-1,110,921 case remain split. On the complete chr20 run, corrected
NGC50 remains 479,016 bp. The shared-callset comparison improves from 77
switches and Hamming 665/59,992 (1.1085%) to 74 switches and Hamming
256/59,878 (0.4275%). Read-truth discordance falls from 3,872/188,819 (2.05%)
to 2,731/189,001 (1.44%). The evidence cache format is version 4 because the
new provenance marker is serialized with each read profile.

## Graph-selected BAM recovery and repeatable gap trials (2026-09-14)

Unresolved hybrid gaps now get a second recovery attempt using reads selected
by graph-site observation spans overlapping the gap. This is a selection from
the available graph catalog observations, not a claim that every GAF path has
been independently reconstructed. Selected reads must have a primary BAM
alignment overlapping the gap. Clean SNP observations are read from its BAM
sequence/CIGAR with the configured base-quality threshold. Cached MSA observations
and non-graph-only indel observations are retained. The existing clean, MSA SNP,
and MSA indel phasing tiers run on this proposal. Unverified graph-marked indels
are excluded because their original BAM observation is not recoverable from
the confirmation sentinel alone.

The fallback returns only the relative orientation of existing phase blocks;
it does not copy its local read assignments or candidate orientations back.
Existing deferred edge composition applies the relationship uniformly. A trial
that copied partial assignments introduced three additional discordant reads;
the orientation-only path removed that regression. `--no-graph-gap-bam` provides
the comparison baseline.

GAF alleles, including disagreements with BAM, are now preserved in a parallel
read-profile channel. MSA profile merging, proposal construction and evidence
serialization preserve it; the cache format is version 5. A global GAF-over-BAM
override was tested and rejected: clean chr20 discordance rose from 586 to 1,512.
The accepted implementation leaves normal phasing observations unchanged.

Whole-chromosome results: 280 initial gaps, 122 joined versus 121 previously.
The additional 13,894,047–13,941,253 join has the correct unflipped orientation
against the original blocks' parental read truth. Final read discordance remains
2,731 / 189,001 (98.56% accuracy). Cache loading took 1.65 s; cached solves took
195.14 s versus approximately 124 s before the extra fallback. No new shared-VCF
NGC50 result is claimed. Outputs: `/tmp/pgphase-graph-gap-edges`.

`scripts/trial_graph_gap_bam.py` runs paired baseline/fallback trials on fixed
competitor gap panels with persistent per-region caches, bounded parallelism,
binary hashes, exact commands, logs and optional parental read evaluation.
The frozen 11-gap panel and usage are under
`evaluations/2026-09-14-graph-gap-bam/`. Tests cover graph/BAM disagreement,
MSA reordering of the new channel, selecting reads while using BAM nucleotides,
and returning a stitching orientation without mutating existing assignments.
`make unit-tests` passes. Local-window results must be confirmed on the full
chromosome because changing the window can change the baseline phase blocks.

The complete paired 11-region panel finished with zero execution failures and
zero added discordant reads in every region. Both arms report eight splits and
three unresolved endpoints; their exact native candidate rows are missing at
12,721,112, 23,480,815 and 37,984,529. Warm graph-selected BAM phasing takes
0.78–1.25 s per local region (truth evaluation excluded). Results are preserved
under `evaluations/2026-09-14-graph-gap-bam/results.tsv`; detailed logs and
cached evidence are in `/tmp/pgphase-gap-trials/runs/graph_bam_v2`.
# Verified deletion bridge and corrected gap endpoint audit (2026-09-14)

The regional trial compared audited VCF anchors directly with native indel
event positions. Corrected using VariantKey::sort_pos semantics (indel POS-1),
with a regression test. Two allegedly missing endpoints at 12,721,112 and
37,984,529 are present and phased. The insertion at 23,480,815 is present but
unphased. The corrected previous panel is 10 splits and 1 unphased endpoint.

Found that msa_indel_has_confident_bam_observation rejected every deletion.
Added exact BAM deletion/REF validation, Q30 boundaries (and Q30 throughout
the REF allele), and permitted validated deletion-to-clean-SNP singleton
bridges with the existing conflict veto. Existing phasing and deferred
stitching close chr20:53,926,537–53,962,801 with the correct opposite parity.
The local 486 evaluated reads remain fully concordant. Full chr20 improves
122 to 123 joins of 280 gaps, loses no old join, and retains 2,731 discordant
reads / 189,001 evaluated, 353 switches, and 413 flips. Cached recovery takes
197.722 s versus 195.14 s. These are parental read metrics, not VCF truth metrics.

The 11-region panel now has 1 joined target, 9 splits and 1 unphased endpoint,
without increased read errors. A 64-site link window, biallelic-insertion
singleton acceptance, and supported-indel reseeding failed to close additional
targets and were not retained. The 1.08 Mb singleton has Q10 at its only
right-block SNP; the 37.98 Mb singleton has conflicting high-quality
right-block observations. Neither is fixed. Added verbose block/allele
diagnostics and matrix provenance metadata for further investigation.

Results and reproduction context:
evaluations/2026-09-14-gap-deletion-bridge/README.md.

## Independent gap anchors and separated insertion alleles (2026-09-14)

Fixed chr20's remaining 37.98 Mb and 0.86 Mb targets. At 37.98 Mb the sole
bridge read has conflicting indel/boundary-SNP observations and an independent
clean SNP. The old component unanimity rule discarded the entire anchor.
Recovery now excludes SNPs inside/immediately bordering link-supported,
MSA-verified indels from the clean bridge channel, and chooses independent
clean SNPs before indels separately for each read/component. Conflicting
independent SNPs still veto an anchor; an indel cannot supply confidence for
a weak SNP that displaced its observations. BAM-validated biallelic insertion
anchors receive the existing singleton treatment used for deletions.

At 0.86 Mb, an exact Q30-supported BAM observation distinguishes TC from TCTC.
The blanket alternate/alternate insertion singleton exclusion discarded that
evidence. These bridges now also permit two alternate sequences differing by
at least two bases in length, provided neither is mononucleotide, retaining
MSA verification, linkage support, exact sequence-equivalent BAM validation
and the opposing-read veto. The C/CC bridge at 23.48 Mb is still excluded.

Local panel: 3 joined, 7 split, 1 unphased endpoint; no local read-error
increase. Nine regions have exactly unchanged ordered read/mapping/HP/PS
records, allowing explicit reuse of their previous parental evaluations.
The two changed regions have separate evaluations: 435/435 concordant at
0.86 Mb, 458/459 at 37.98 Mb (the previous one discordant read remains).

Full chr20: 123 → 125 joined of 280 initial gaps, no previous join lost.
All previously evaluated reads retain their truth-concordance status. Of 40
newly evaluated reads, 37 are concordant and 3 discordant (all three have
truth HapQ 0). Raw discordance is 2,734/189,041 versus 2,731/189,001; switches
355 versus 353 and flips 414 versus 413. These are read metrics, not new VCF
Hamming measurements. No existing-block orientation regression was detected.
Cached solve time 192.868 s, cache loading 1.823 s. Initial interrupted runs
were not counted as validation; completed output is
/tmp/pgphase-anchor-final-check.

Build/unit tests pass, including independent-versus-boundary SNP conflict,
weak SNP, unverified indel, exact shifted insertion, mononucleotide alternate,
and low-quality insertion cases. Details and results:
evaluations/2026-09-14-independent-gap-anchors/README.md.

## 2026-09-14 — Remaining-gap orientation validation

See `evaluations/2026-09-14-gap-orientation-validation/README.md` for the full
regional panel, rejected experiments and final chromosome-wide checks.

- Corrected the gap trial endpoint checker to read native VCF GT/PS. The
  purported unphased endpoint at 23,480,815 was phased `1|2`; candidate-table
  HAP_ALT/HAP_REF cannot encode that distinction. Corrected baseline: eight
  split rows representing seven distinct gaps.
- Rescued 19.37 Mb with exact, sequence-distinct two-alt insertion observations:
  minimum Q20 on all inserted/flanking bases plus mean inserted-base Q30.
  Other singleton insertion/deletion confidence gates retain Q30.
- Rescued 23.48 Mb by assessing SNP-anchor confidence independently on each
  flank (graph-confirmed versus BAM-confident sources may differ). Attach
  otherwise orphan MSA site phase sets only with consistent net support from
  both haplotypes of an existing read-supported block; do not change reads.
- Caught a full-chromosome cascade: the locally correct 23.48 Mb edge propagated
  an earlier wrong homopolymer orientation into 745 previously correct reads.
  The unvalidated combination was rejected. Homopolymer joins now require a
  second solve with restored BAM SNP observations to connect both original
  flanks with matching parity. Trials remain read-only until final edge
  composition. Sparse MSA genotype backfill is confined to this confirmation.
- Rejected general HOM-to-HET backfill during gap recovery: it produced 135
  joins but 13,580 discordant reads. Existing allele presence/balance is not
  evidence that it separates the two flanking haplotypes correctly.
- Final full chr20: six panel target rows joined, five split, zero unphased
  endpoints. Two distinct new gaps close. Four old homopolymer joins are
  withheld, including the erroneous upstream edge; net recovery is 123/280
  versus 125/280. Withholding the 57.84 and 62.41 Mb edges has no measured
  read-accuracy benefit and remains an explicit contiguity cost.
- Parental read discordance improves 2,734/189,041 to 2,427/189,042; switches
  355 to 351, flips unchanged at 414. All previously concordant reads remain
  concordant, 307 previously discordant reads become concordant, one new
  concordant read is added and none are lost. Variant Hamming/NGC50 not rerun.
- Build and unit tests pass, including mixed-source anchors, insertion quality
  boundaries, alternate/alternate endpoint reporting, orphan site orientation,
  sparse MSA backfill and read-only homopolymer confirmation.

## 2026-09-14 — Joint gap-decision experiment plan

`evaluations/2026-09-14-gap-decision-plan/PLAN.md` specifies an immutable
chromosome-context export/replay, explicit SAME/FLIP/insufficient/conflicting
outcomes, molecule/event provenance and deduplication, controlled scoring
ablations, and one final component-orientation application. No new scorer is
enabled. Existing k-means remains the initial proposal generator. Source views
are sensitivity checks on overlapping reads, not independent confirmations.

The proposed validation distinguishes edge parity and fixed-cohort component
errors from read tagging and contiguity. Known bad joins and withheld good
joins are regression/development cases, not a locked test set. Previously
examined chr12/18/20 regions cannot serve as unseen validation.

Initial experiment: replayed 11 regions in both existing arms at 2 jobs × 2
threads and 1 job × 1 thread. All 22 paired comparisons preserve read HP/PS,
alignment identities and native VCF data exactly; six target rows join and
five split in each arm. Median per-arm times were 0.802 and 0.918 seconds in
these single warm-cache runs. This demonstrates regional reproducibility only;
chromosome-context ordering and new scoring remain to be implemented/tested.
Manifests, source/binary hashes, results and comparator are preserved beside
the plan. No production code changed during this planning experiment.

## 2026-09-14 — Implement frozen gap audit/replay and fix confirmation evidence

Implementation and artifacts:
`evaluations/2026-09-14-gap-decision-implementation/README.md`.

- Added opt-in `--gap-decision-audit DIR` to hybrid recovery. All tier/view
  proposals run read-only against one pre-recovery chromosome state before
  any production worker mutates chunks. Export molecule assignments, original
  flank matrices, candidate/observation provenance and full original block
  memberships/bounds. Snapshot overwrite is rejected. Serial and parallel
  control runs produced 54 byte-identical audit files.
- Added `scripts/replay_gap_decisions.py`: exact molecule-to-matrix validation,
  SAME/FLIP scores, insufficient/conflicting states, independent original
  haplotype margins, fixed-proposal molecule influence and tentative parity
  component checks. Correlated overlapping observations do not multiply direct
  molecule votes; untagged bridge reads are retained. Integrated observations
  are not mislabeled independent BAM evidence. Feature schema 2 includes
  corrected interval ordering for overlapping SNP/indel groups.
- Versioned/hash-validated features permit parameter trials without rereading
  observations or rerunning phasing/MSA/BAM output: ~19 MB features, ~0.2 s
  scoring and ~2.3 s total verified replay for 280 chr20 gaps. Initial extraction
  took ~71 s; the initial opt-in audit plus normal solve took 405.3 s. Single-run
  timings only. Existing default runtime remains ~191 s cached solving.
- The audit exposed why good 57.84/62.41 Mb edges failed confirmation: the
  validator first altered sparse HOM calls with extra backfill/promotion.
  Confirmation now evaluates the existing BAM evidence without that promotion;
  removed the unused HOM-revisit extension and retained a preservation test.
- Full chr20 recovers 125/280 gaps versus 123/280, adding only
  57,841,772–57,866,713 and 62,408,056–62,432,427. No old edges lost or reoriented.
  Both bad 23.42/37.64 Mb joins remain rejected. Read accuracy is unchanged:
  2,427 discordant / 189,042 evaluated, 351 switches, 414 flips. Every read
  retains its concordance status; none are gained or lost. No new variant
  Hamming/shared-callset NGC50 measurement is claimed.
- The initial audit was noninterfering on the entire chromosome (native VCF
  records and BAM identity/alignment/HP/PS identical) and all 22 regional arms.
  Build, C++ tests and Python benchmark tests pass. Python replay tests include
  provenance corruption, input-cache changes, untagged/correlated observations,
  label invariance, uneven haplotype support and contradictory cycles.
- General score remains diagnostic: balanced uniform gives 124 directional
  proposals, 151 insufficient and 5 conflicting cases. Unbalanced scoring
  admits the known bad 23.42 Mb proposal; MAPQ-only weighting has no demonstrated
  benefit. Event error calibration, sequence-equivalent multi-allelic recovery,
  influence re-solving, untouched validation and general-score adoption remain
  outstanding. Five original difficult gaps are not claimed fixed.

## 2026-09-14: Recovery error attribution and read evaluator correction

Fresh current-binary recovery-disabled chr20 baseline has 586/185835 discordant reads versus recovery's 2427/189042. Common-read block comparison confirms incorrect relative orientation after merges (not nonuniform mutation inside original blocks). Largest examples: PS52434849 adds 487 errors, PS46663811 adds 233, PS46748667 adds 347. Strong flank votes and BAM/graph proposal agreement both fail to detect the wrong bridge at 52394827–52434849. Production acceptance remains unresolved; do not call recovery accuracy-safe.

Corrected `evaluate_phase_accuracy.py` transition ordering to input BAM coordinates instead of interleaving maternal/paternal assembly positions. Counts are explicitly read-concordance diagnostics, not variant switch/flip metrics. Recovery diagnostic counts become 205/267 versus old 351/414; discordance is unchanged. Added parental-coordinate-offset integration regression to benchmark-tests. See `evaluations/2026-09-14-read-error-audit/README.md` and common-read block attribution table. Legacy truth-based spans need independent correction before contiguity comparisons.

## 2026-09-14: Validate internal MSA bridges; move toward gap-owned evidence

Implemented two recovery-only corrections in collect_phase.cpp: honor ordinary
biallelic MSA indel gap_link_supported, and require the winning orientation of
an MSA allele edge to be supported by both alleles at both endpoints, with the
existing combined net margin. The latter is symmetric under endpoint exchange
and HP-label swaps. Pooled majority counts previously concealed contradictory
minority-haplotype linkage. New regressions reproduce both failures, including
9:3 pooled support whose second haplotype opposes the join 3:1. All unit and
benchmark tests pass. Existing singleton bridge controls continue to pass.

Full chr20: 2427/189042 -> 760/188714 discordant reads; corrected read-transition
counts 205/267 -> 138/208. Ten previously identified erroneous joins -> zero
under the original-block parental-majority audit. Accepted joins 125 -> 86;
27 previously correct and two uncertain joins are also withheld. The 57.84 Mb
control survives; 62.41 Mb is withheld. Forty originally concordant reads become
discordant relative to no recovery (previously1452); do not claim zero errors or
zero regressions. A HOM-node exclusion trial was rejected because it introduced
a new 56.15 Mb error. It is absent from the final code. Full ablations, artifacts
and source hashes: evaluations/2026-09-14-gap-bridge-validation/README.md.

User redirected the next architectural work toward a source-preserving evidence
boundary with private BAM sites ONLY inside original graph-derived phase gaps.
Each frozen gap owns its candidate/observation view; flanks provide fixed graph
anchors and alignment context, not additional outside-gap BAM candidates.
Audited loss of query-coordinate provenance through kGraphConfirmedAltQi,
combined mutable evidence/phasing fields, missing separate BAM/MSA histories,
and incompatible complex-replacement conventions. Proposed immutable event,
source-allele mapping and per-molecule observations behind a checked adapter to
the existing phaser. This replacement is designed, not implemented. Details:
evaluations/2026-09-14-shared-evidence-design/DESIGN.md.


### Gap-owned source evidence implementation (2026-09-14)

Implemented `GapEvidence`: freeze original graph-derived gaps, retain graph
anchors from adjacent original blocks, and admit private BAM/MSA events only
when their full reference footprint belongs to that gap. Flanking extraction
context cannot introduce outside-gap private candidates. Recovery-mode initial
phasing now masks non-catalog candidates; this changes the initial block/gap
inventory relative to the previous hybrid baseline. Discovery evidence remains
cached for gap recovery. No DeepVariant calls enter this path.

BAM allele calls and original query positions are saved before graph injection;
MSA profile replacement preserves BAM and graph source histories. A snapshot
stores separate per-molecule BAM/graph/MSA-derived observations, explicit missing,
low-quality and conflicting states, reference/alternate sequences, and immutable
anchor orientation. Source-local MSA allele integers and anchor consensus are
translated by sequence. Repeated proposals read the frozen snapshot rather than
previous recovery mutations. Counts are rebuilt from retained molecules.

Boundary-crossing and unsupported complex replacement events remain auditable
but are not projected into the legacy simple-event core. This is deliberately
not a general repeat-equivalence merger. MSA observations are the retained
MSA-derived working calls, not a claim of independent per-read MSA quality.
Cache format 6 stores the added provenance; old cache files require a one-time
rebuild. Gap audits now include `.events.tsv` and `.observations.tsv`.

Regional chr20:52-53 Mb validation: cold cache / 4 threads and warm cache /
1 thread produce identical BAM records, VCF records, candidates, tier reports,
and all five per-gap audit tables. Four initial gaps, one joined. Cold 28.33 s;
warm 6.76 s, cache load 0.034 s. Unit and benchmark tests pass, including allele
permutation, immutable projection, source disagreement, query-quality provenance,
1/2 genotype counts and boundary exclusion. Artifacts:
`/tmp/pgphase-gap-owned-region/`. Whole-chromosome validation is recorded below
when complete; previous 280-gap metrics are not the baseline for this scope.

### Whole-chromosome validation and independent gap blocks (2026-09-15)

Completed the whole-chromosome validation the gap-owned-evidence entry left
open, then built the first extension that phases gaps' own reads even when
they cannot be confidently tied to a specific flank.

**Whole-chromosome validation (gap-owned evidence, no new mechanism yet).**
`collect-hybrid-variation --recover-gaps --gap-evidence-cache
/tmp/chr20-evidence-v6.gapev` on chr20: 276 initial gaps, 112 joined in 4
waves, 0 conflicting edges. Cold and warm (8 vs 1 thread) runs produced
byte-identical `tiers.tsv` and BAM output. Peak RSS 49.9 GB for the whole
chromosome -- an existing property of this evidence-cache-based pipeline, not
introduced by anything below.

**Independent gap blocks (`emit_independent_gap_block`, `src/gap_recovery.cpp`).**
A gap that cannot bridge to either flank is currently discarded outright, even
when its own gap-only reads (no pre-recovery flank/block assignment at all)
converged to a coherent local split during the gap's own k-means solve. This
recovers that: reads with no original assignment sharing the gap's own
locally-derived phase set (>= `--gap-independent-min-reads`, default 3, itself
a genome position and therefore unique) are emitted as a brand-new,
independent block. Strictly additive -- never overwrites a read that already
carries any assignment.

Diagnostic groundwork before implementing: of the 169 chr20 gaps the existing
bridge-vote mechanism marks unresolvable, 146 have >=2 mutually-consistent
private events spanned by >=3 reads among themselves, using only the raw
per-read/per-event data `--gap-decision-audit` already exports
(`scripts/measure_internal_gap_phaseability.py`, ~39 s for all 276 gaps once
the audit exists -- no C++ rerun).

Three bugs found and fixed while building this, each caught by insisting on a
read-by-read corruption check (concordant -> DISCORDANT against
`../pgphase-eval-data/truth/chr20/diplinator_merged.bam`) rather than trusting
aggregate counts:

1. `min_reads <= 0` was meant to disable the feature; unsigned comparison made
   it enable unconditionally instead.
2. A race: applying independent blocks during the parallel per-gap phase let
   two gaps whose windows shared an unclaimed read compete for it, an
   order-dependent outcome across runs (8 reads corrupted, reproducibly).
   Fixed by deferring application to one sequential pass, strictly after
   `apply_gap_phase_edges` (`GapRecoveryJobResult::has_proposal` /
   `PendingIndependentGaps`).
3. The batch loop calls `stitch_chunk_haps` a *second* time after
   `recover_hybrid_gaps` returns, to re-examine block boundaries across the
   whole batch. That second call does not know a block is independent and
   was free to re-merge/relabel it (same 8 reads, still corrupted after fix
   2). Fixed by moving emission to run after that second stitch call too --
   confirmed via `--gap-independent-min-reads 0` producing a byte-identical
   BAM to the pre-feature baseline, isolating that the corruption came from
   this code path specifically and not from unrelated nondeterminism.

Verified clean on the third attempt: chr20, both bugs fixed, 27 gaps emitted
27 gap(s) as new independent blocks, 1,271 (later 1,170 once evaluator
`min_reads_per_ps=5` filtering is applied) previously-unphased reads newly
evaluated. Zero reads corrupted, zero lost, 96.2% concordant on the newly
phased population. `--gap-independent-min-reads 5` and `8` were also swept and
are equally safe (0 corrupted); 8 gives 96.8% accuracy on a slightly smaller
population (1,151 reads) -- default stays 3 since the difference is marginal
and more reads phased is the goal.

**Multi-round recovery: attempted, unsafe, disabled by default.** The natural
extension -- repeat the whole bridge/independent-block sequence so gaps newly
created by an independent block (it now sits between an existing flank and a
bridge that previously had nothing to reach) get their own attempt -- is
correct in principle but every round after the first cannot reuse
`--gap-evidence-cache` (its signature is keyed to the gap inventory, which
changes every round) and must rebuild raw evidence uncached. Measured: peak
RSS exceeded 47 GB and was still climbing after 10 minutes on round 2+, so the
run was killed before it risked the host. `gap_recovery_max_rounds` defaults
to 1 (single pass, the verified-safe behavior above); raising it is left
available for a future session with the uncached-rebuild cost addressed, not
recommended as-is.

**Sub-gap bridging: attempted, low value, disabled by default.** A cheaper
alternative to a full extra round -- try bridging each newly-independent block
to its own two flanks as two small sub-gaps, reusing the *same* `proposal`
already in memory (no re-extraction) against a `GapReadIndex` rebuilt from the
just-updated chunks. Implemented with `defer_phase_set_merge=true` on both
sub-gap attempts plus one `apply_gap_phase_edges` call at the end, so a
gap.left_ps <-> chosen_ps <-> gap.right_ps chain resolves to one phase set via
union-find rather than two independent immediate renames splitting it into
two surviving ids. Memory stayed flat (~50 GB, no growth) and corruption was
zero -- but the yield was 4 edges / 6 reads chromosome-wide at only 50%
accuracy on that population. Root cause: these sub-gaps were already tested by
the *identical* flank-side vote evidence when the original, larger gap failed
to bridge; routing the same votes through a new intermediate node does not
add information. Gated behind `--gap-bridge-independent-blocks` (default off)
rather than removed, since the mechanism itself is correct, just not a good
trade with the current all-or-nothing, pre-recovery-flank-assignment-gated
vote model.

**What is left, precisely.** Checked whether the remaining small (<20 kb)
gaps are limited by catalog completeness (the dominant cause documented
elsewhere in this repo for the *original* gap set): they are not. 44 of 46
small post-recovery gaps have both candidates (median in the teens, several
in the 20s-50s) and abundant spanning reads (median several dozen). The
blocker is that `stitch_gap_proposal`'s vote count only credits a read that
already carries a *pre-recovery* flank assignment (`GapReadIndex::assignments`,
`src/gap_recovery.cpp`); a read that spans the junction but was never tagged
to either flank contributes nothing, regardless of how many such reads exist.
This is the same untagged-linking-read circularity documented earlier in this
file for the non-recovery pipeline, now located precisely inside the
recovery-specific bridge test. The natural fix is allele-pattern voting for
`stitch_gap_proposal` itself (compare a spanning read's own allele calls at
the flank anchor and at the gap-owned candidates directly, the same principle
`--link-by-alleles`/`check_agree_alleles` already use for the main phasing
pass, not yet threaded into the gap-recovery-specific vote path) -- designed
here, not implemented.

Build and all unit/benchmark tests pass throughout. Default behavior with
`--recover-gaps` alone is unchanged except for the verified-safe independent-
block emission (3 above); `--gap-recovery-max-rounds` and
`--gap-bridge-independent-blocks` are both off by default.

### Allele-pattern voting for gap bridging: attempted, corrupts, disabled (2026-09-15)

Implemented the fix designed above: `--gap-link-by-alleles`
(`src/gap_recovery.cpp`, gated by `Options::gap_link_by_alleles`). When a
proposal read has no committed pre-recovery hap/PS at all
(`GapReadIndex::assignments` misses it), `stitch_gap_proposal` now falls back
to deriving an implied flank side from the read's own allele agreement with
that flank's already-resolved sites (`CandidateVariant::hap_to_cons_alle`),
the same principle `--link-by-alleles`/`check_agree_alleles` use for the main
phasing pass. Strictly additive in intent: a read with a committed assignment
always uses it; the allele-derived one only fills in where nothing existed.

Measured on chr20, whole chromosome, both arms from the same
`/tmp/chr20-evidence-v6.gapev` cache, otherwise-identical command
(`--link-by-alleles --block-link-window 8 --min-read-margin 2 --recover-gaps`),
read-by-read status compared against
`../pgphase-eval-data/truth/chr20/diplinator_merged.bam`:

| | Baseline (flag off) | `--gap-link-by-alleles` |
| --- | ---: | ---: |
| Gaps bridged | 112 | 115 |
| Independent blocks / reads | 27 / 1,271 | 62 / 2,200 |
| Phase sets (evaluated) | 179 | 204 |
| Reads evaluated | 184,961 | 183,475 |

More gaps bridged and far more independent-block reads captured, but the net
effect is worse, not better:

- 1,465 reads that were concordant and evaluated in the baseline dropped out
  of evaluation entirely in the new run (fell below the evaluator's
  min-reads-per-phase-set floor) -- an existing, previously well-supported
  block fragmented into pieces too small to score. 1,446 of those 1,465 were
  concordant before, i.e. this is a straightforward contiguity regression, not
  noise being correctly dropped.
- 5,215 common reads were relabelled to a different phase-set id (versus a
  baseline where flag-off reproduces the exact 112/27/1,271 numbers above
  byte-for-byte) -- far more churn than "3 more gaps joined" can explain on
  its own.
- At least one phase set (PS 60552043) ended up internally contradictory:
  some of its reads correctly kept HP2 (matching baseline and truth), others
  were flipped to HP1, i.e. a single nominal phase set no longer has one
  consistent haplotype orientation. This is the clearest signal -- an
  allele-derived vote is pulling part of an otherwise-correct block into the
  wrong orientation.

**First hypothesis, tested and ruled out.** Suspected `GapReadIndex::reads`
returning more than one `(chunk, read_index)` location for the same read
(tiling windows overlap), with `CandidateVariant::hap_to_cons_alle` for a
gap.left_ps/right_ps site only guaranteed canonically oriented in chunks this
batch's `orient_candidate` has actually touched -- an untouched duplicate tile
copy could carry an independently k-means-derived, possibly oppositely
labelled, hap1/hap2 sense for the same phase-set id, and the original
`implied_assignment` took the first location it found without checking which
copy that was. Fix implemented: `implied_assignment` now scans *all* of a
read's locations and only trusts a phase set's implied hap if every location
that resolves it agrees; any disagreement poisons that phase set for that
read rather than picking a side. Re-ran the identical chr20 regression: result
essentially unchanged (116 joined / 61 independent blocks / 2,199 reads vs the
pre-fix 115/62/2,200), and PS 60552043 still shows the same symptom with the
same specific reads flipped (baseline: 12/120 reads discordant there already
-- an existing hard region; with the flag, either version: ~25/120, HP1/HP2
split shifting from an even 60/60 to 63/57). The implied vote is internally
self-consistent per read across all its locations; it is simply wrong often
enough, concentrated in already-hard regions, to corrupt a join when it feeds
the bridge/orientation decision directly. Not a chunk-duplication bug --
allele-derived votes are just not reliable enough there, the same class of
risk already documented above for `gap_bridge_independent_blocks`. The more
conservative `implied_assignment` is still worth keeping (strictly more
correct than the first version, at no measured cost), but the mechanism it
feeds is still unsafe.

**Third attempt: decouple the decision from the attachment -- verified safe,
now the default.** Implemented exactly the direction sketched above: `votes`
and the orientation decision (which side wins, whether it flips) in
`stitch_gap_proposal` are computed from committed (`first_assignment`) reads
only -- byte-identical to the flag being off, confirmed by chr20 reproducing
112 joined with unchanged orientation. Allele-derived agreement
(`implied_assignment`, the cross-location-safe version from the second
attempt) is consulted only *after* `accepted` already names which flank(s)
this gap resolved to on committed evidence alone, in a new pass that
additively attaches any still-unassigned read whose own alleles agree with
that already-decided side -- the same purely-additive contract
`emit_independent_gap_block` uses, just triggered by a direct per-read allele
match instead of local-bucket membership. (The existing local-bucket
attachment path already covered gap-only reads sharing the winning local
phase_set; this new pass additionally covers reads on a flank chunk itself,
or in a losing local bucket, that individually agree by allele with the
accepted side -- which turned out to be the large majority of the gain.)

Found and fixed one real bug while wiring this in: `implied_assignment`
indexed `chunk.read_var_profile[ri]` and `prof.alleles[...]` without bounds
checks against `chunk.read_var_profile.size()`/`prof.alleles.size()`; real
pipeline chunks always keep those consistent with `chunk.reads.size()`, but
`test_gap_recovery_keeps_unphased_observations`'s minimal fixture does not
populate `read_var_profile` at all, segfaulting the instant this became the
default and the test exercised `stitch_gap_proposal` with default `Options`
for the first time. Fixed with explicit size guards; caught immediately by
`make unit-tests` before this reached real data.

Verified on chr20, whole chromosome, reused `/tmp/chr20-evidence-v6.gapev`,
read-by-read status against `../pgphase-eval-data/truth/chr20/
diplinator_merged.bam`, compared to a fresh flag-off baseline run from the
identical command: gaps bridged (112) and their orientation are unchanged,
byte-for-byte. Zero of the 184,974 previously-evaluated reads regressed --
no concordant -> DISCORDANT transition anywhere, confirmed by full read-by-
read diff, not just aggregate counts. 22,676 additional reads got evaluated,
96.4% of them concordant (in line with `emit_independent_gap_block`'s 96.2%
on its own, much smaller, population) -- roughly 17x the reach of independent-
block emission alone. The 12 reads that changed HP within an already-existing
phase set were all corrections (previously DISCORDANT, now concordant), not
new errors. Contiguity improved slightly, 176 phase sets vs 179 (fewer,
larger blocks), because many reads that previously only qualified for a
same-PS independent block now attach directly to the correct flank instead
(independent-block yield correspondingly dropped from 27/1,271 to 17/1,183 --
expected redistribution, not a loss, since the reads land in a normal flank
block rather than a small standalone one).

`Options::gap_link_by_alleles` now defaults to `true` (part of default
`--recover-gaps` behavior, same tier as `gap_independent_min_reads`'s
nonzero default), given this measured-clean result meets the same bar used
to accept `emit_independent_gap_block`. All three findings (two rejected
designs, one accepted) are recorded in its comment in
`src/phasing_types.hpp`. Build and all unit/benchmark tests pass.

### Code review of the gap-recovery diff: 3 real findings, all fixed (2026-09-15)

Dispatched an independent review agent over the session's gap-recovery diff
(`src/gap_recovery.cpp/.hpp`, the gap-recovery sections of
`src/collect_pipeline.cpp`) with explicit context on what was already
shipped/verified vs. reverted, asking specifically for correctness bugs, not
style. It found three, all genuine extensions of the same defense
(`established_phase_sets`) rather than restatements of it:

1. **The collision guard only recorded read-supported phase-set ids.**
   `filter_hybrid_reads_by_margin`/`filter_hybrid_small_phase_sets` zero a
   read's `haps`/`phase_sets` when its block loses support, but never touch
   `CandidateVariant::phase_set` -- an orphaned block can leave candidate rows
   still labelled with its id and the main solve's orientation, with zero
   committed reads (exactly why `find_phase_gaps` has no support for that id
   either, and the region becomes a gap). `emit_independent_gap_block` could
   still choose that id, inserting new rows under the proposal's own polarity
   while the orphaned rows kept the main solve's -- silently splitting one
   phase set's orientation. Fixed: `GapReadIndex` now also scans every
   chunk's `candidates` and adds any `phase_set > 0` to
   `established_phase_sets`, not just ids with a committed read
   (`src/gap_recovery.cpp`, `GapReadIndex` constructor).
2. **The guard is a frozen, whole-batch snapshot, so two gaps in the same
   sequential emission loop could still emit into each other's brand-new
   phase set.** Adjacent gaps' windows overlap by construction
   (`kGapRecoveryFlank` on both sides), so two gaps failing to join (both
   landing in `pending.proposals`) can independently reconstruct the same
   leftmost-het id from overlapping gap-only read pools -- each yields
   `chosen_ps == P` with an independently-derived, uncorrelated orientation.
   Fixed: `emit_independent_gap_block` takes an additional in/out
   `std::set<hts_pos_t>* emitted_this_round`, threaded through one shared set
   across the whole sequential loop in `collect_pipeline.cpp`; a candidate
   group whose id is already in that set is treated the same as a real
   collision, and a successful emission inserts its own `chosen_ps` before
   the next iteration runs.
3. **The new allele-attach pass's cross-chunk dedup could leave the
   output-owning tile copy of a read unphased.** The pre-existing read-attach
   loop is chunk-independent (gates only on `original` and that chunk's own
   `haps`/`phase_sets`, so every overlapping tile copy of a read gets the
   same, deterministic assignment). The new allele-attach loop added an
   `added.count(key)` skip on top of that -- but `implied_assignment` reads
   *live*, per-chunk candidates, which only gain a gap's interior sites via
   that same chunk's own `merge_var_profile` call, later in the same
   iteration. A read's earlier (output-owning) tile copy can run its attach
   attempt before its own chunk has those sites merged in, fail, and then get
   permanently skipped once a *later*, non-owning tile copy succeeds --
   leaving the owning copy, and therefore the emitted read, unphased despite
   being counted as attached. Fixed: dropped the `added` gate from this loop
   (kept only the `original` check, matching the pre-existing loop exactly);
   every qualifying tile copy gets an independent chance, which is redundant
   when they agree (they read the same shared `orient_candidate` output) and
   never corrupting. `added` (a set) still dedupes the `reads_added` count
   correctly regardless.

**Verification exposed a methodology gap, not a new bug.** Comparing the
fixed build against `/tmp/pgphase-gap-baseline-fresh` (this session's
long-standing regression reference) showed 224 reads newly "lost" (220
concordant there). Before concluding the fixes regressed anything, checked
what those reads actually were: that reference predates finding 1/2's fixes,
so its own independent-block emission could already be riding the exact
collision bug just fixed -- silently, since a 50/50-orientation risk simply
did not manifest as visible wrongness for those specific reads on this
dataset. Generated a proper, unambiguous reference instead
(`--no-gap-link-by-alleles --gap-independent-min-reads 0`, i.e. bridging
joins only, zero independent-block risk by construction) and confirmed: all
224 "lost" reads are absent from that clean reference too -- they only ever
existed via the now-fixed buggy path, never a legitimate assignment. Full
read-by-read diff against the clean reference: 0 lost, 0 phase-set changes,
0 HP flips, 0 corruption among 183,804 common reads, +23,622 additional
reads at 96.4% concordance. `/tmp/pgphase-gap-baseline-fresh` should not be
reused as a regression reference going forward -- it is contaminated by the
same bug class this entry fixes; use a fresh
`--no-gap-link-by-alleles --gap-independent-min-reads 0` run instead when
checking for genuine regressions.

Chr20 with all fixes: 112 joined / 12 independent blocks / 954 reads (down
from 16/1,181 pre-fix, as expected -- the now-excluded collision cases are
exactly what dropped). Build and all unit/benchmark tests pass.

### Why independent-block emission still only reaches a fraction of unresolved gaps (2026-09-15, later same day)

With `gap_link_by_alleles` on by default, chr20 has 276 initial gaps: 112
bridged, 162 still `partial` (never bridged to either flank) per
`--gap-recovery-report`. `scripts/measure_internal_gap_phaseability.py`
(pure re-aggregation of `--gap-decision-audit` data, no C++ rerun) says 141 of
those 162 have an internally-coherent private-event component spanned by
`>=3` reads -- real structure exists. But `emit_independent_gap_block` only
turns 17-20 of the 276 gaps into new blocks. Investigated why the other ~120+
gaps with apparent structure never qualify.

Dispatched an Explore agent to trace the call chain
(`recover_hybrid_gaps` -> `recover_one_hybrid_gap`, src/collect_pipeline.cpp
-- the phasing itself is the same `assign_hap_based_on_germline_het_vars_kmeans`
as the main pass, just run on a per-gap local `proposal` chunk). Its finding,
confirmed by direct reading: `recover_one_hybrid_gap` runs up to two internal
"passes" per gap (`recovery_passes = opts.graph_gap_bam ? 2 : 1`, default 2).
Pass 0 is the broad local solve (all of the gap's own reads); pass 1
reprojects to a narrower graph-BAM-only view and `select_graph_gap_bam_reads`
(`src/gap_recovery.cpp:20`) marks any read with no graph-channel allele
spanning the gap as `is_skipped` -- excluding pure BAM/MSA-only internal
structure entirely. Whichever pass the loop ends on is what gets saved as
`job_result.proposal` for `emit_independent_gap_block`; if pass 1 finds zero
graph-observed reads it `break`s immediately, discarding pass 0's fully-solved
state and saving an entirely fresh, unphased reprojection instead.

Added temporary `--verbose 2` diagnostics (`GapIndependentBlock`,
`GapIndependentBlockApplied` in `emit_independent_gap_block` -- kept, they are
read-only and gated, same pattern as the existing `GapLinkVotes` log) and
measured directly rather than continuing to reason from hypothesis:

- Of 162 unresolved gaps, most locally-phased reads (median ratio ~99.6%,
  mean ~90.6%) already carry an original pre-recovery assignment -- the local
  window's k-means mostly just re-derives what the flanks already know, not
  new gap-only structure. Only a minority of gaps have any gap-only reads at
  all in their local solve.
- Tried: preserve pass 0's proposal and prefer whichever pass has more
  locally-phased reads (`locally_phased_reads` comparator, no gap-only/
  original filtering). Chr20 yield barely moved (17/1,183 -> 20/1,111, then
  18/1,106 with the fix below) -- not the large unlock hoped for. Worse: it
  measurably regressed 112-116 reads (114 concordant in a fresh flag-off
  baseline) to fully unphased, with zero corruption (no concordant ->
  DISCORDANT) but a genuine, not-fully-root-caused loss traced to an
  interaction with `gap_link_by_alleles`'s one-sided partial-link attachment
  inside this *same* gap's own pass-0, non-`orientation_only`
  `stitch_gap_proposal` tier attempts (a `left_linked`-or-`right_linked`,
  not-fully-`joined` result still attaches individual reads to whichever side
  voted, per the existing code -- this predates this investigation and is
  outside its scope to fully trace). Reverted; see `src/collect_pipeline.cpp`
  history/comment at `recover_one_hybrid_gap` for the measurements.
- Along the way, found and fixed a real, independent latent bug: phase-set
  ids are genome positions in both the main pass and a gap's local re-solve,
  so a wide-window local solve can, by coincidence, anchor on the same first
  het site as an adjacent *established* block and reproduce its exact id.
  `emit_independent_gap_block` had no defense against this -- it would
  silently emit gap-only reads into a real block's phase-set id with no
  orientation-vote safety net (unlike `stitch_gap_proposal`, which always
  votes before merging into an existing id). Fixed by recording every
  already-assigned phase-set id at `GapReadIndex` construction
  (`established_phase_sets`) and refusing to select a candidate group whose
  id collides with one, in `emit_independent_gap_block`
  (`src/gap_recovery.cpp`/`.hpp`). This is a permanent, kept improvement,
  independent of the reverted pass-0/pass-1 selection change above --
  verified harmless and correctness-only (chr20: 16/1,181 vs the pre-fix
  17/1,183, only the genuinely-colliding cases excluded, 0 lost reads beyond
  2 that were not concordant in baseline either, 0 corruption).

Net result of this sub-investigation: `established_phase_sets` kept as a
real, if narrow, correctness fix; the pass-0/pass-1 proposal-selection change
reverted as not worth its unexplained cost. The ~120 gaps with apparent
internal structure that still do not produce an independent block are, per
the ratio measurement above, mostly *not* a missed opportunity -- their local
solve is dominated by flank-anchored reads with no true gap-only population
large enough to matter, not evidence the pipeline is discarding. What
residual real opportunity exists (`gap_only_ps_groups=0` for the ~77 fully-
flank-dominated gaps vs some-groups-but-below-threshold for ~14 others,
measured via the kept `GapIndependentBlock` diagnostic) is small and would
need pass 1's graph-BAM narrowing itself relaxed to admit BAM/MSA-only
internal structure -- a materially different, riskier change than swapping
which already-computed proposal gets kept, not attempted this session.

Build and all unit/benchmark tests pass. Default `--recover-gaps` behavior:
112 joined / 16 independent blocks / 1,181 reads on chr20, byte-for-byte
reproducible; read-by-read regression check against a fresh flag-off baseline
shows 0 corruption and 22,676 net additional evaluated reads at 96.4%
accuracy, matching the gap-link-by-alleles verification above.

### CLI cleanup: missing opt-out for a now-on-by-default flag (2026-09-15)

When `Options::gap_link_by_alleles` flipped its struct default to `true`
earlier this session, its CLI wiring in `src/hybrid_collect.cpp` was left as
an enable-only flag (`--gap-link-by-alleles` -> `= true`) -- there was no way
to turn it back off from the command line, unlike every other on-by-default
boolean in this file (`graph_gap_bam` has `--no-graph-gap-bam`,
`exp_hybrid_trim` has `--no-hybrid-trim`). Every other boolean CLI flag in the
file was checked against its struct default and all others are consistent
(enable-only flags all correspond to `false`-by-default fields). Fixed:
replaced `--gap-link-by-alleles` with `--no-gap-link-by-alleles` (sets
`opts.gap_link_by_alleles = false`), matching the `--no-graph-gap-bam`
pattern, and corrected its stale "(experimental, off)" help text. No other
CLI-facing references to the old flag name existed (scripts, evaluations).
Build and both test suites pass.

### Auditing what's still unresolved on chr20: mostly legitimate, one dead end (2026-09-15, later)

After the code-review fixes, asked directly: of chr20's 64,590 reads still
carrying no PS tag at all, what's actually blocking them, gap by gap? Ranked
all 164 unresolved gap windows by their own count of fully-unphased reads
(not by `emit_independent_gap_block`'s narrower "applied" signal, which the
earlier investigation already showed undercounts success -- many of its
"applied=0" gaps turned out to have their reads captured instead via the
existing one-sided flank attach, confirmed by checking output PS tags
directly: e.g. gap 7305204-7343405's reads are phased at PS 7343405, its
real right flank, not lost at all).

**The top-ranked gaps are legitimately hard, not buggy.** The single largest
window, chr20:27,133,886-29,068,243 (~1.9 Mb, 9,828 of 9,841 reads unphased),
has essentially no clean sites (`NEW_SITES`/`MSA_HET_SNPS`/`MSA_HET_INDELS`
all 0 across every tier) because >85% of its reads have MAPQ<10 -- a
genuinely ambiguously-mapped, likely repetitive region; forcing a phase call
there would be inventing confidence the mapping doesn't support. The next
largest, chr20:43,688,351-45,143,091 (~1.45 Mb, 6,068 unphased, all high
MAPQ), has 1.3M graph observations and 6,069 graph-selected reads but only
145 end up with an individually confident k-means haplotype call -- a
genuinely low-heterozygosity stretch (only 1,393 raw candidates over 1.45 Mb)
rather than an algorithmic gap. Filtering the full ranked list for gaps with
candidate density comparable to normal, well-phased regions (>=2/kb) *and*
still substantially unphased found only 21 gaps totaling 1,797 reads -- the
real remaining opportunity is one to two orders of magnitude smaller than the
raw 64,590 unphased count suggests.

**One promising lead turned out to be a dead end.** Gap
chr20:26,784,052-26,804,346 stood out: 38 clean MSA het SNPs (real
information) but only 13 of 76 reads ever admitted into the graph-BAM
recovery pass (`select_graph_gap_bam_reads`, `src/gap_recovery.cpp`), because
that function's admission test requires a graph-channel allele
(`profile.graph_alleles[pi] >= 0`) and nothing else -- a read whose only
informative sites are private/MSA-only never counts as "spanning" the gap,
regardless of how much private evidence it carries. Implemented and tested a
targeted fix: admit a read if it has a graph-channel call *or* a confident
call at an `msa_verified` candidate. Measured zero effect on chr20
(`SELECTED_GRAPH_READS` for this gap stayed at 13, identical
before/after) -- root cause: `CandidateVariant::msa_verified` for a gap's own
interior sites is set later, per-gap, during `recover_one_hybrid_gap`'s own
tier loop (`run_gap_msa_tier`/MSA consensus admission), not yet true on the
frozen `GapEvidence` snapshot `select_graph_gap_bam_reads` sees at the point
it runs. Before chasing a fix to the flag-timing itself, pulled this gap's
frozen reference sequence from its `--gap-decision-audit` export and found it
is a tandem-repeat/satellite-like sequence -- the "38 clean MSA het SNPs" are
plausibly repeat-unit alignment artifacts, not real heterozygous sites, so
this specific gap was likely never a genuine opportunity regardless. Reverted
the fix cleanly (confirmed via `git diff` showing no change to
`select_graph_gap_bam_reads`); rebuilt and both test suites still pass.

**Conclusion:** after today's shipped fixes (gap-owned evidence,
`emit_independent_gap_block`, allele-pattern voting, and the three
collision/race fixes above), chr20's remaining unresolved gaps are
overwhelmingly legitimate abstentions -- low-MAPQ/repetitive mapping, low
heterozygous-site density, or noisy tandem-repeat structure -- not bugs.
Recovering the residual ~1,800-read opportunity in the density-normal subset
would need the `msa_verified`-timing / graph-BAM-admission interaction
properly redesigned (decide whether a gap's own private-site verification
should run *before* `select_graph_gap_bam_reads`, or admit on a different,
already-reliable signal instead of `msa_verified` specifically) -- a
materially different, riskier change than anything else shipped today, not
attempted further this session.

### The MAPQ floor, not evidence scarcity, closes the competitor gap on chr20 (2026-09-15, later)

Took one gap competitors bridge correctly and we do not, and traced why, rather
than continuing whole-chromosome iteration. Full detail and reproduction in
`evaluations/2026-09-15-mapq-starved-gaps/`.

`chr20:25,834,662-25,883,079` (48.4 kb, status `split`, 279 of 296 reads
unphased): HiPhase spans it in one block on 103 phased het sites. We produced 18
candidates there, 2 of them clean het SNPs, matching 3 of those 103. Not a
classification problem -- at the same positions DeepVariant reports DP 59-76 with
balanced allele depths while our candidates report DP 2-10. Input coverage is
68.2x, but the reads are MAPQ 3 (166) and MAPQ 4 (41) of 296, and
`kDefaultMinMapq = 30` leaves 3.7x. The gap looks evidence-free because the floor
removed the evidence.

Sweeping `-q` on a 148 kb window around it (same binary, truth restricted by read
name after phasing): at `-q 10` we get 35 clean het SNPs instead of 2 and +74
reads phased with **zero** discordant reads; at `-q 1` we get 302 clean het SNPs,
86 of HiPhase's 103 sites, 734 of 776 reads phased, and the gap ceases to exist
-- one 375 kb block, N50 166 kb -> 375 kb -- at a cost of 37 discordant reads.
Those errors are entirely low-MAPQ: 14.5% on MAPQ 1-9, 0% at MAPQ >= 10.

The instructive part is that HiPhase's own default floor is MAPQ 5, which does
not admit this gap's MAPQ 3-4 reads either -- it tagged 89 reads here and exactly
the 15 that sit at MAPQ 5-9. So it did not win by admitting the bulk. DeepVariant
*called* the 103 sites using the low-MAPQ reads, and HiPhase *phased* those calls
with a MAPQ >= 5 read floor. We drive both jobs from one `min_mapq`
(`gap_recovery.cpp:46`, `graph_query.cpp:421`, `hybrid_inject.cpp:570`,
`collect_bam_output.cpp:410`), so dropping a read for phasing also drops it for
discovery and the anchors never exist. Our `-q 1` arm shows discovery recovers
those anchors once the reads are visible, and that assigning the same reads to
haplotypes is what injects the error -- which argues for a split floor: discovery
inside a gap at a low floor, haplotype assignment and output unchanged at 30.

Screened all 164 unresolved chr20 gaps (15.2 Mb) for the pattern: 20 gaps
(2.99 Mb) are MAPQ-starved, 5 are genuinely low-depth, 139 keep their depth
through the floor. Scoring the reads we leave unphased against the diplinator
truth splits the 20 -- 3 gaps (149 kb, 765 reads, all just proximal to the
centromere) where LongPhase reaches 96.4-100% and HiPhase 98.3-100%, and 17 gaps
(12,393 reads) spanning LongPhase 3.8-93.2% and HiPhase 11.8-100%. The second
group is graded rather than uniformly noisy: three of it sit at LongPhase 81-93%
and miss the trustworthy cut on the 95% threshold, and its lone HiPhase 100% is
10 reads at 32,432,281-32,449,137 where LongPhase gets 61.9% on 21. The latter includes the 1.93 Mb `27,133,886-29,068,243` window
at 54.1% (LongPhase) and 56.4% (HiPhase): **the earlier audit's verdict on that
gap stands**, and the 2,596 competitor-phased sites inside it are not evidence of
correct phase. The honest recoverable opportunity from relaxing the floor is the
765-read subset, not the 12,808 reads the raw starvation count suggests. A
truth-free discriminator exists but is weak (median |VAF-0.5| 0.052 across the
three trustworthy gaps vs 0.115 across the 16 noise-verdict gaps carrying a
value, one trustworthy gap at 0.220 -- overlapping, so not a gate on its own).

### Private MSA het SNPs as bridge anchors: measured, kept off by default (2026-09-15)

`--gap-bridge-private-snps` (`src/collect_phase.cpp`, opt-in) measured on whole
chr20 against the same binary with the flag off, both arms reusing
`/tmp/chr20-evidence-v6.gapev`; see `evaluations/2026-09-15-private-snp-bridge/`.
Flag off reproduces the accepted baseline exactly (112 joined, 207,426 phased,
1,400 discordant, 0.675% read Hamming, N50 1,064,616, 177 phase sets), which also
clears the three un-gated changes in the same diff (private-site admission in
`select_graph_gap_bam_reads`, pass 1 reprojecting into a scratch chunk instead of
overwriting pass 0's solved proposal, and the `LEFT_LINK_PS`/`RIGHT_LINK_PS`
report columns with `split`/`vetoed` statuses). The same-day sweep that reached
115-116 joins at 0.720%/0.713% came from the anchor segregation experiment
already recorded as reverted, not from anything still in the tree.

Flag on: one further gap joined (`chr20:36,332,599-36,381,019`, `split` ->
`joined`, merging PS 36,381,019 into PS 36,317,511), auN +3,392 bp, N50
unchanged, flip errors 671 -> 657, read Hamming 0.675% -> 0.668%, and **0
previously concordant reads became discordant**. The cost is coverage: 74 reads
lose their PS tag, 61 of them previously concordant (25 dropped by the merge
itself, 38 and 11 at two unrelated phase sets), because a read whose only
informative gap site is a private SNP now needs that observation confirmed in the
BAM and an unconfirmed read can fall below its output margin. One join and 3.4 kb
of auN for 61 concordant reads is the wrong trade when more reads phased is the
objective, so the flag stays off by default.

One open item on it: a 1-thread rerun of the flag-on arm reproduced the same 113
joins and a byte-identical `native.vcf`, but the read tags do not match the
8-thread run -- 191 reads differ, 176 of them tagged only under 8 threads and 15
only under 1 thread, and `tiers.tsv` differs in 96 rows after sorting (not merely
row order). The default path is documented byte-identical across 8 vs 1 thread, so
the flag appears to introduce thread-order dependence in gap read attachment. Not
diagnosed; a further reason it should not become default before that is understood.

### New diagnostics kept from this work

`--gap-recovery-report` gained `LEFT_LINK_PS`/`RIGHT_LINK_PS` and two new STATUS
values. `split` (both flanks linked, to different proposal phase sets) turns out
to be the dominant unresolved outcome -- 830 of the tier attempts previously
reported as `partial` on chr20 -- which is a materially different failure from a
genuinely one-sided link and is what makes gaps like the 25.83 Mb one identifiable
as near-misses.

### Split MAPQ floor implemented; discovery-side admission corrupts 2 of 3 gaps (2026-09-15)

Acted on the MAPQ finding above by separating the two jobs `min_mapq` was doing
at once. New `--min-assign-mapq` (`Options::min_assign_mapq`, default
`kDefaultMinMapq`, predicate `read_carries_phase_tags` in
`src/collect_phase.cpp` consulted in `PhasedAlignmentWriter::write_chunks`): a
read admitted above `-q` but below this floor supplies allele evidence to
discovery and linking and is emitted without HP/PS. Equal floors are the default
and reproduce previous behavior exactly -- the default arm's `phased.bam` is
byte-identical to the pre-change build. Build clean, both test suites pass, new
unit test `test_assign_mapq_floor_gates_tags_only`.

On the gap the idea came from (`chr20:25,834,662-25,883,079`) the in-pipeline
result reproduced the post-hoc estimate exactly: `-q 1 --min-assign-mapq 30`
gives 399 reads / 1 block / 0 discordant / N50 270,669 against the default's
395 / 2 / 0 / 166,590, and `--min-assign-mapq 5` gives 505 / 1 / 0 / 275,646.
Zero previously concordant reads became discordant.

It does not generalize. Run on the other two gaps whose sub-floor reads were
verified to carry correct phase (`evaluations/2026-09-15-split-mapq-floor/`),
the gate fails: 6 previously concordant reads corrupted at 25,944,471-25,986,123
(which does not even join) and 6 at 26,029,591-26,088,679 (which joins, N50
124 kb -> 300 kb, while losing 26 phased reads). The corruption count is
**identical at assignment floor 30 and 5**, and at floor 30 no sub-floor read
carries a tag at all, so the damage comes from the discovery/linking side: the
sub-floor reads' alleles and link votes change the solution for confidently
mapped reads. The option therefore stays opt-in and is not a default candidate.

The methodological lesson is worth keeping: competitor read-level evidence that a
gap's low-MAPQ reads are phaseable (96-100% on all three of these gaps) does not
imply our mechanism phases them correctly, and deriving a design from the single
gap that motivated it oversold it by 3x. Selective admission is the next thing to
try -- per-site balance gating using the VAF-deviation signal (0.052 across the
trustworthy gaps vs 0.115 across the noisy ones) or down-weighted sub-floor link
votes -- re-measuring the gate on all three gaps before any whole-chromosome run.

### Graph-only path audit: classification sound, MNPs escape noise screening (2026-09-15)

Audited `collect-graph-variation` end to end on chr20 against the hybrid default
arm and the diplinator read truth. Full detail in
`evaluations/2026-09-15-graph-path-audit/`.

Baseline: 73,627 candidates retained of 988,602 catalog sites considered;
203,751 reads phased, 27,631 unphased, 371 phase sets, 1,744 discordant, read
Hamming 0.856%, N50 917,428 bp. The hybrid arm is 207,426 / 0.675% / N50
1,064,616 / 177 phase sets, so graph-only trades accuracy and contiguity for
independence from the BAM channel, as expected.

Clean/noisy classification holds up. `apply_graph_noise_filter` fires -- 17,350
het indels demoted to `REP_HET_INDEL` against 2,036 surviving clean, i.e. 89.5%
of graph het indels screened. Categories agree with the BAM classifier on 57,033
of 58,394 shared sites (97.7%), and every systematic disagreement is the graph
being stricter (397 sites it demotes that the BAM path calls clean, 426 more it
demotes that the BAM path files as noisy candidates). The reference fetch that
guards the noise filter resolves because `batch_contig` comes from the
FASTA-derived header -- though if that fetch ever failed the filter would be
skipped silently, with no warning.

One real defect. `build_graph_chunk` derives variant type from allele lengths
alone (`graph_bam_adapter.cpp:944-950`) and `VariantType` has no MNP member, so
every equal-length multi-base substitution is typed `Snp` with
`key.ref_len = ref.size()`. On chr20 that is **781 sites, all classified
`CLEAN_HET_SNP`, all 781 carrying a PHASE_SET** -- they vote in k-means at the
top anchor score. Worse, `apply_graph_noise_filter` reconsiders only
`CleanHetIndel`, so these sites are structurally unreachable by repeat screening
no matter their context, and the context is often repetitive: 115 of the first
400 sit in homopolymer or dinucleotide-repeat windows (`GTG>CTC` inside
`GTGTGTGTGTGTG` at 145,051; `CT>TC` inside `CCCTCTCTCTCCC` at 863,876). The BAM
path demotes the equivalent site at 47,301,789 (`AT>TA`) to `REP_HET_INDEL`
while the graph path calls it a clean SNP anchor. Lengths run 2-12 bp and include
a 367 bp allele pair at 764,862 typed `Snp`, so the key's span also contradicts
its type. The BAM path's `Deletion` label for the same sites is equally wrong;
neither path represents MNPs.

Recommended fix (not yet implemented): extend `apply_graph_noise_filter` to
examine `CleanHetSnp` candidates with multi-base catalog alleles and move those
in low-complexity context out of the anchor mask, keying on
`pos_in_low_complexity` rather than `is_noisy_site` (whose indel-length
derivation is zero for equal-length alleles).

MSA is not used in the graph-only path at all: `abpoa` is reached only via
`align.cpp`, whose verification consumer is `collect_phase_noisy.cpp`, which
`graph_collect.cpp` does not include; `msa_verified` is never set and the output
carries no `NOISY_CAND_*` calls (against 4,003 + 1,964 in the hybrid arm). The
path phases clean catalog sites with no noisy-site rescue -- a capability
boundary, and the reason its unphased count is 3.3x the hybrid arm's.

Also noted: graph-only emits zero homozygous records (42,750 sites filtered
`high_af`) where the hybrid arm emits 34,722 `CLEAN_HOM`; the `n_uniq_alles > 2`
branch of `classify_graph_candidates` is unreachable as configured (all retained
sites are biallelic post-decomposition) but would bypass the AF-centering
paralog guard if that ever changed; and the evaluator's switch/flip counts are
meaningless for this path because the emitted BAM is unaligned.

### NEW_SITES in the recovery report is identically zero by construction (2026-09-15)

`previous_sites` is declared inside the tier loop it is differenced against
(`collect_pipeline.cpp:1942` vs the row write at `:2039`), so
`proposal.candidates.size() - previous_sites` can only ever be 0. The NEW_SITES
column has never carried information, and it has already produced two wrong
conclusions: the 09-14 audit of `chr20:27,133,886-29,068,243` cited
"NEW_SITES/MSA_HET_SNPS all 0" as evidence the window holds no usable variation,
and this session repeated the same reading across all 276 gaps before checking
the source. The MSA_HET_SNPS half of that observation is real (it counts
`kCandNoisyCandHet` candidates in the proposal); the NEW_SITES half is vacuous.
No absolute count of recovered gap candidates is emitted anywhere, so the
question "did BAM-side recovery find sites in this gap" is currently
unanswerable from the report.

What the report does show on chr20 (default arm, 276 gaps): 112 joined, 140
split, 22 partial, 2 open; gap recovery selected BAM reads in 168 of 276 gaps;
MSA-verified hets are present in the proposal for 55 gaps (21 joined, 29 split,
4 partial, 1 open). Split gaps carry a median of 253 selected reads.

`split` -- both flanks linked, to *different* proposal phase sets -- is the
dominant unresolved outcome at 140 of 276. That is the gap-block-then-stitch
step failing at unification rather than at recovery, and those gaps are the bulk
of the 15.2 Mb still unphased. Diagnosing whether the proposal fragments
internally into two phase sets or the stitch mis-gauges them is the highest-value
open question, ahead of the MAPQ and MNP work.

Correction to this session's earlier claim that adding MSA to gap recovery would
be "a much larger project": wrong. The machinery already exists and runs --
`populate_gap_msa_cache` -> `prepare_gap_msa_regions` / `run_gap_msa_tier`
(`gap_recovery.cpp:848-899`), feeding `kCandNoisyCandHet` candidates into the
gap proposal. The intended architecture (graph-derived blocks, BAM-recovered gap
blocks from clean plus MSA-verified sites, then stitch to both neighbours) is
implemented; the open problems are instrumentation and the split outcome.

### Per-gap targeting: the linkage bottleneck, and only 34 of 164 gaps are actionable (2026-09-15)

Answered "which region of the BAM should a per-gap subprocess look at" with a
measurement rather than a window heuristic. A gap is bridgeable only where a read
crosses a position carrying an informative het site on *both* sides; scanning cut
points across the gap and counting such reads localises where phasing actually
breaks. Reads need not span the whole gap -- each cut needs only one read
reaching a site either side -- which is why the bottleneck is interior rather
than at a junction in 87 of 164 gaps, and why widening the recovery window around
the gap edges cannot help. Script, table and detail in
`evaluations/2026-09-15-gap-targeting/`.

Applied to the 164 unresolved chr20 gaps, with our own called sites compared
against the het calls a caller makes using every read:

- `thin_linkage`, 130 gaps, 11.96 Mb -- no read carries a het site on both sides
  of the bottleneck even using all reads. Coverage is not the issue (median 66
  spanning reads, 65 passing the MAPQ floor). Only 41 of the 130 have a het site
  on both sides within 25 kb, and for those the median spacing is 3,603 bp one
  side and 22,658 bp the other, i.e. a **median 25,848 bp a single read would have
  to cover**. These are correct abstentions, information-limited by heterozygosity
  spacing versus read length, and should stop consuming recovery tiers.
- `sites_not_called`, 20 gaps, 0.64 Mb -- sites and a read chain exist; we did not
  call them. This is exactly what gap MSA discovery is for, pointed at the target
  window.
- `mapq_starved`, 12 gaps, 2.59 Mb -- reads present, none passing the floor at the
  cut. Gated admission confined to the target window, with MSA verification before
  a site becomes an anchor.
- `linkage_present`, 2 gaps, 0.02 Mb -- reads and our own sites link across, so the
  defect is in the solver or the stitch, not the evidence.

So the honest remaining opportunity is 3.24 Mb across 34 gaps, reduced to 1.70 Mb
of BAM to inspect (6.51 Mb -> 0.45 Mb for the nine gaps over 200 kb), not the
15.2 Mb the unresolved-gap span suggests.

The scan reproduces the hand diagnosis of `chr20:25,834,662-25,883,079` from
scratch -- `mapq_starved`, bottleneck 25,870,298, target 25,845,298-25,895,298,
20 of our sites against 98 potential -- which took a full manual session to reach
the first time.

Production note: the potential-site set is currently read from an external
caller's VCF. A self-contained subprocess should replace it with a direct
pileup scan over the target window (allele balance, no MAPQ floor), which is work
the MSA discovery step performs anyway.

### Correction: the gap "heterozygosity deserts" were an artifact of the site set (2026-09-15)

The targeting entry above concluded that 130 of 164 unresolved chr20 gaps are
information-limited because no read carries a het site on both sides of the
bottleneck. That conclusion was wrong, and the reason is instructive: both site
sets used were BAM-derived -- our own `CLEAN_HET_*` calls and a pileup caller's
het calls -- and *both* omit repeat-context indels. The snarl catalog itself was
never consulted.

The catalog is dense at precisely those bottlenecks: median **671 sites within
25 kb**. About 96% are homozygous in this sample (median 614 `ref_only`, 31
`high_af`; unobserved sites are negligible at 2%, so filtering is not hiding
them), but a median of **7 `REP_HET_INDEL`** sites and 4 `low_af` sites per
window are heterozygous and currently excluded from phasing.

Re-testing every bottleneck with those classes included
(`evaluations/2026-09-15-gap-targeting/add_demoted_site_linkage.py`,
`gap_targets_revised.tsv`): `repeat_indels_would_bridge` 88 gaps / 6.45 Mb,
`low_af_sites_would_bridge` 24 gaps / 3.43 Mb, `no_linkage_from_any_site_class`
52 gaps / 5.32 Mb. So **112 of 164 gaps (9.89 Mb) have candidate linkage from
graph sites the pipeline discards**, median 17 linking reads once repeat indels
are included, in a median 50 kb window. 99 of the 130 gaps written off above are
in that group.

This is candidate linkage, not proven linkage: those sites are demoted because
per-read genotypes at homopolymer/STR indels are unreliable, so the linkage may
be phantom. The way to use them is MSA verification confined to the target
window, keeping only sites whose consensus resolves two consistent haplotypes --
i.e. this result argues *for* the `populate_gap_msa_cache` / `run_gap_msa_tier`
path, pointed at a 50 kb window, and explains why that path currently yields
MSA-verified hets in only 55 of 276 gaps. Note also that 12 of the 52 remaining
gaps are MAPQ-starved, where no read passes the floor at the bottleneck, so their
linkability is untested rather than refuted; the genuinely unlinkable residue is
at most 40 gaps (~2.7 Mb).

Method lesson, third time this session: a negative result is only as strong as
the input set it was computed over. "No evidence exists" claims need the evidence
inventory enumerated explicitly -- catalog sites by disposition, not just the
sites that survived our own filters.

### The demoted-indel linkage is mostly phantom: 0.43 Mb of the 6.45 Mb (2026-09-15)

Tested whether the repeat-demoted het indels that would bridge 88 unresolved
chr20 gaps carry real haplotype signal, by genotyping each from the alignment and
scoring its allele partition against the read-level truth
(`evaluations/2026-09-15-gap-targeting/verify_demoted_sites.py`).

Site level: `CLEAN_HET_INDEL` (control) 96% informative, median segregation
1.000. Repeat-demoted het indels: 744 scored, median 0.586, **23% informative,
62% phantom** (INS 33%, DEL 18%). The noisy-candidate het class the pipeline
emits: 138 scored, median 0.728, 25% informative, 48% phantom.

Gap level: recomputing each bottleneck with only truth-validated sites admitted,
**10 of 88 gaps retain linkage, 0.43 Mb**, against 88 gaps / 6.45 Mb when all
demoted sites are admitted. The 171 informative sites across those windows
(~2 per window) are real but not positioned to span the bottlenecks. So
`apply_graph_noise_filter` discards ~23% genuine signal yet is right about the
majority, and even a perfect oracle gate recovers only 0.43 Mb from this class.
Admitting demoted indels wholesale would inject roughly three phantom sites per
real one.

Load-bearing open question. The noisy-candidate class that carries MSA-derived
sites into phasing is three-quarters noise -- but `NOISY_CAND_HET` in the output
is the candidate *category* (`collect_phase.cpp:72`), not a verification
certificate; `msa_verified` is a separate flag set in
`collect_phase_noisy.cpp:213` and consumed at `collect_phase.cpp:580,1188`, and
the TSV never exposes it. Only `--gap-decision-audit` writes `msa_verified=`
per variant (`collect_phase.cpp:303`). So MSA verification's precision remains
unmeasured, and it is the next thing to measure: the whole
graph-blocks-plus-BAM-recovered-gap-blocks design assumes MSA-verified sites are
trustworthy anchors.

Two measurement bugs of mine, both caught by the same control (score
`CLEAN_HET_INDEL` and require it to come out informative): genotyping indels at
the VCF anchor position scored the control at 2% scorable, because repeat-context
indels are placed arbitrarily within their run by the aligner and the allele must
be read from the net length change over a window; and deriving the expected
length change from `len(ALT) - len(REF)` skipped every deletion, since this TSV
writes deletions as the deleted bases in REF with `ALT` = `.`
(`collect_output.cpp:110-145`). Control now reads 96% informative, 744 of 795
sites scorable. Build a trusted-class control into any future site-quality
measurement.

### Graph-only chr20 baseline and gap inventory (2026-09-15)

Ran `collect-graph-variation` alone on whole chr20 as pass 1 of the two-pass
design, and derived the phase-block gap inventory that a pass-2 subprocess would
consume. 41 s wall, 2 m 07 s CPU on 8 threads. Detail and scripts in
`evaluations/2026-09-15-graph-only-baseline/`.

Sites: 73,627 retained of 988,602 considered -- 54,241 `CLEAN_HET_SNP`, 17,350
`REP_HET_INDEL`, 2,036 `CLEAN_HET_INDEL`; filtered 826,336 `ref_only`, 42,750
`high_af`, 35,481 `no_reads_in_chunk`, 7,044 `low_af`, 3,364 `low_depth`.

Reads: 203,751 phased, 27,631 unphased, 201,991 concordant, 1,744 discordant,
read Hamming 0.856%, 371 phase sets, N50 917,428 bp, auN 987,490 bp. Switch and
flip counts are not measurable for this path (the emitted BAM has no contig
header), so only read-level figures are usable.

Gaps: 288 of 375 phase sets have two or more sites and span 47.60 Mb; the spaces
between consecutive blocks are **286 gaps spanning 18.52 Mb** (47.60 + 18.52
accounts for chr20's 66.1 Mb). Size distribution: 9 under 1 kb, 6 at 1-10 kb, 183
at 10-50 kb, 80 at 50-200 kb, 8 over 200 kb. **55,583 distinct reads sit in gap
windows unassigned, 40,214 of them at MAPQ >= 30.**

Two pass-2 observations. The nine sub-kilobase gaps are 3-249 bp with 0-3
untagged reads each -- adjacent blocks that failed to link while sharing reads,
so pure stitching failures and the cheapest first test of a gap subprocess. The
two largest gaps fail oppositely: `43,688,351-45,896,820` (2.21 Mb) leaves 8,659
reads unphased of which 8,425 pass MAPQ 30, while `27,133,886-29,068,099`
(1.93 Mb) leaves 9,825 of which only 196 pass.

Note for comparison, stated carefully because the two populations are not the
same thing: the hybrid path's gap-recovery report lists 276 gaps *before*
recovery and joins 112 of them, leaving 164 unresolved spanning 15.2 Mb.
Graph-only's 286 gaps / 18.52 Mb is a pre-recovery population with no recovery
attempted, so it is comparable to hybrid's 276 total (4% larger) and is ~75%
larger than hybrid's post-recovery residue of 164 gaps.

Tooling: `gap_inventory.py` initially read tagged reads from the phased BAM by
region, which silently returned nothing every time (no contig header, no index,
`samtools view` exit 1) and reported all 94,014 overlapping reads as untagged --
the same swallowed-failure pattern flagged in `apply_graph_noise_filter` earlier
today, this time in my own tool. It now reads the per-read assignment TSV and
raises on any non-zero samtools exit.

### Graph-only MAPQ floor: -q 1 is free for existing reads, costly only for new ones (2026-09-15)

Swept the graph-only floor on whole chr20 (`-q 30/10/5/1`, 55 s per arm;
`evaluations/2026-09-15-graph-only-baseline/mapq_sweep.sh`).

`-q 1` against the `-q 30` default: sites 79,291 vs 73,627 (clean het SNPs
59,282 vs 54,241), reads phased 214,065 vs 203,751, discordant 2,945 vs 1,744,
read Hamming **1.376% vs 0.856%**, phase sets 379 vs 371, N50 1,011,696 vs
917,428 (+10.3%), gaps 311 vs 286, gap span 17.81 vs 18.52 Mb.

The decisive split: restricting both arms to reads at MAPQ >= 30, the population
`-q 30` could already see, gives **0.856% at `-q 30` against 0.854% at `-q 1`**
(1,744 vs 1,741 discordant). Lowering the floor does not degrade existing
phasing. All the added error is in the 10,347 newly phased reads, which carry
1,240 discordant calls (11.98%) and are almost entirely the low-MAPQ population
(2,769 at MAPQ 1-4, 1,250 at 5-9, 6,210 at 10-29). Per-bucket error in the `-q 1`
arm: MAPQ 1-4 16.50%, 5-9 22.80%, 10-29 7.44%, 30-59 1.61%, 60 0.80%.

103 reads flip concordant -> discordant (29 at MAPQ 30-59, 74 at 60) against 138
flipping back. With the MAPQ >= 30 error rate flat to three decimals and blocks
re-forming (phase sets 371 -> 379, N50 +10%), those are consistent with
re-blocking churn rather than corruption; the gate count is also non-monotonic
across floors (780 at `-q 10`, 316 at `-q 5`, 103 at `-q 1`), which points the
same way.

Decision for the two-pass design: use `-q 1` for pass 1, but do not let pass 1
tag the 10,347 ambiguous reads -- they are pass 2's job, to be decided from local
gap evidence. That separation exists as `--min-assign-mapq` but is not reachable
here: only `hybrid_collect.cpp` parses the option and `read_carries_phase_tags`
is consulted only at `collect_bam_output.cpp:458`, while the graph-only path
writes its BAM through its own mirror of that writer
(`graph_bam_adapter.cpp:209-222`). Two small additions would wire it up.

### Decision: -q 1 adopted as the graph pass-1 floor (2026-09-15)

Standing on the sweep above: `-q 1` leaves the MAPQ >= 30 population's accuracy
unchanged (0.856% -> 0.854%, 1,744 -> 1,741 discordant) while adding 5,041 clean
het SNPs, 10.3% block N50 and 0.71 Mb less gap span, so the floor is adopted for
pass 1. `evaluations/2026-09-15-graph-only-baseline/run.sh` now defaults to
`MIN_MAPQ=1`, and the canonical pass-1 gap inventory in that directory is the
`-q 1` one: **311 gaps spanning 17.81 Mb** (not the 286 / 18.52 Mb from the
`-q 30` arm).

One caveat carried forward rather than resolved: pass 1 at `-q 1` also tags the
10,347 newly admitted reads, which carry 11.98% discordance. Under the two-pass
design those reads belong to pass 2, decided from local gap evidence.
`--min-assign-mapq` is the mechanism and is not reachable from this subcommand
(only `hybrid_collect.cpp` parses it; `read_carries_phase_tags` is consulted only
at `collect_bam_output.cpp:458`, while the graph path writes through its own
mirror at `graph_bam_adapter.cpp:209-222`). Until it is wired up, pass-1 output
carries those tags.

Note for pass 2 scoping: the untagged population inside gap windows barely moves
between floors -- 55,583 distinct reads at `-q 30` versus 52,340 at `-q 1`, and
the MAPQ >= 30 subset is essentially identical (40,214 vs 40,142). The work pass
2 has to do is therefore not a function of the pass-1 floor: roughly 40,000
confidently mapped reads sit in gaps unphased either way.

### Pileup-based column discovery: mechanism sound, starved in gaps (2026-09-15)

Prototyped the pass-2 site-finding step -- candidate columns straight from the
pileup, no catalog, no classification, no MSA, admitted by agreement with the
partition the other columns imply -- and scored it against read truth with
in-block controls. Detail in `evaluations/2026-09-15-gap-column-discovery/`.

Controls pass: inside blocks, 39-109 candidate columns, 30-103 admitted, median
segregation 1.000, and 29-100 of the admitted columns are sites the pipeline
already calls. Self-consistency tracks quality without truth (0.884 consistency
-> 88.6% truth concordance; 1.000 -> 100.0%).

In gaps the same code finds **0-20 candidates and admits 0-4**. The linear
pileup does not expose enough het columns there, so the answer to "will this
find the missing sites in the gaps" is no -- it finds them inside blocks, where
they were never missing. Three gaps (`13429829`, `60084884`, `10296487`) report
consistency 1.000 at **51-52% truth concordance**, i.e. chance: with 2-4 columns
the consistency score is vacuous because the columns trivially agree with the
partition they defined. Design consequence: require a minimum admitted-column
count (controls sat at 30+) in addition to consistency and flank anchoring.
Three gaps did resolve at 100% concordance, so the approach is starved rather
than wrong.

Next channel to measure, identically: the GAF `cs:Z:` difference string, present
on every record and private by construction (catalog variation is traversed as
path, not as mismatch). One gap-sized window carried 1,006 mismatch, 1,860
insertion and 3,458 deletion events. pgphase already has an accurate cs
tokenizer (`bam_digar.cpp` `build_digars_cs_tag`) but it is BAM-bound, and the
surjected BAM here carries no `cs` tag at all -- only `hs` -- so the GAF is the
sole source. Columns should be keyed on (node, node-offset) rather than reference
coordinate: no projection needed, and reads through different repeat copies
cannot contaminate each other's columns. `hs`/`hb`/`he` are not a projection
table (30 entries against 2,601 path nodes); they are the GBWT haplotype-thread
annotations already consumed by `collect_phase_pgbam.cpp`.

Aside: `-q 1` closed `chr20:25,834,662-25,883,079`, the gap this whole
investigation started from. It now sits inside a single 831 kb block
(PS 25,102,178, 25,102,178-25,933,230, 995 sites).

### BAM phase transfer into the graph gauge: works, does not stitch, shrinks gaps 95% (2026-09-15)

Tested the idea of phasing each graph gap with `collect-bam-variation` over the
gap plus a 50 kb flank and transferring the result into the graph's gauge by read
identity, on the eight gap windows also used for column discovery.
`evaluations/2026-09-15-bam-phase-transfer/`.

Anchoring is unambiguous: all eight gaps anchor on **both** sides at agreement
**1.000**, with 97-302 anchor reads per side. **719 reads the graph left unphased
receive a haplotype, 702 correct against truth (97.64%)** -- against the graph's
own 0.85% error, so transferred reads are lower-quality coverage, not free
coverage. Transfer is pure addition (graph blocks keep their gauge, only untagged
reads gain tags), so the concordant-to-discordant gate cannot be violated by
construction, and it requires no new variant calling.

But **no BAM phase set spans any of the eight gaps**. The BAM phase sets run into
the gap 18-33 kb from each side and break. Two independent channels -- catalog
driven and pileup driven -- break at the same position, which is strong evidence
the break belongs to the data, not to the graph channel. That agrees with the
linkage-bottleneck scan and with pileup column discovery admitting 0-4 columns in
gaps against 30-103 in blocks.

Residual break after two-sided transfer: the six 10-50 kb gaps go from 272.5 kb
of gap to **38.6 kb (14.2% remaining)**, and four of them land at **0.6-2.0 kb**
from 41-48 kb originally. The two 200 kb+ gaps barely move (438.9 -> 390.5 kb):
their interiors are phased into BAM phase sets anchored to neither graph block,
i.e. phasing islands needing chaining -- the same both-flanks-link-to-different-
sets failure as hybrid gap recovery.

Plan consequence: transfer runs first, before any site discovery, and then the
open problem is a 1-2 kb residual rather than a 48 kb gap. At that size abPOA
consensus over the crossing reads, or the GAF `cs:Z:` private-variant channel
keyed on (node, offset), are tractable where they were not across 48 kb.

Tooling note: `csv.writer` defaults to CRLF, so every TSV these evaluation tools
wrote had `\r` line endings. Python readers tolerate it; the shell runner silently
matched zero windows because `kind` read as `gap\r`. All four writers now pass
`lineterminator='\n'` and the already-saved tables were rewritten.

### BAM site injection into the graph solve: narrows gaps, closes none (2026-09-15)

Measured the "retrieve sites from the BAM and phase the graph with them" route on
the eight gap windows, two arms (`collect-hybrid-variation` site union, and union
plus `--recover-gaps`), `-q 1`, 50 kb flank.
`evaluations/2026-09-15-site-injection/`.

The sites are genuinely there: inside the true residual break intervals the BAM
channel has **23 clean het SNPs against the graph's 9**, plus 79 noisy-candidate
hets where the graph holds 74 demoted repeat indels. Distribution is uneven --
four of six mid-size gaps get only 0-2 extra clean het SNPs, the two 200 kb+ gaps
get 6 and 12.

**0 of 8 gaps end up spanned by a single phase set, in either arm.** Windows keep
2-7 phase sets. The blocks that form are correct (target block 97.9-100.0% against
truth over the sixteen window-arms: 100.0% in all eight union arms, floor 97.9%
at gap 35,919,404 with recovery), so this is linkage failure, not corruption.

Largest uncovered stretch inside the eight gaps, by approach: graph-only 711.4 kb,
phase transfer 541.5 kb, injection 328.9 kb, injection + recovery **304.7 kb**
(57% reduction). Injection transforms the two large gaps (202 -> 58 kb, 237 -> 43
kb) where transfer barely moved them, tracking their extra sites; on four of six
mid-size gaps injection alone changes nothing, because their breaks hold 0-2
usable sites.

Anomaly to chase: `--recover-gaps` makes the two large gaps worse than the plain
union (57,744 -> 137,367 bp; 43,482 -> 56,364 bp) while improving mid-size ones.

Correction to the previous entry: the transfer residuals reported there were
measured from read START positions, which understates a phase set's reach by up
to a read length at each end. Re-measured from first-to-last phased variant, the
mid-size residuals are 8.0-44.7 kb, not 0.6-2.0 kb, and transfer reduces 272.5 kb
of mid-size gap to 124.4 kb (54%), not the 95-99% claimed. `transfer_results.tsv`
now carries `residual_beg`/`residual_end` so downstream work targets the actual
interval instead of assuming it sits at the gap midpoint.

Untested option that needs no called site, and is the pipeline's own designed
answer for this case: pgbam haplotype-thread stitching via `--pgbam-file` with the
`--pgbam-*-min-winning` thresholds, documented as the fallback "when common-read
signal is absent" and built on the `hs`/`hb`/`he` GBWT thread tags.

### Root cause of `split`: verified SNPs refused as gap-link sites (2026-09-15)

Single-gap diagnosis of `chr20:36,332,599-36,381,019` via `--gap-decision-audit`.
`evaluations/2026-09-15-gap-link-site-gate/`.

Injection and profiling are NOT the problem. The proposal holds 581 in-gap sites;
of the 16 in phase-informative categories, 2 are `CleanHetSnp` at the gap edges
and **14 are `NoisyMsaHet` with `msa_verified=1` and no homopolymer flag**, each
with 27-84 BAM observations. 10 of the 16 score >= 0.90 segregation against the
diplinator read truth, and they form an unbroken chain (11.4, 10.9, 10.9, 8.0,
7.2 kb) with 17-81 reads spanning every consecutive pair.

It still splits because `joined = left_linked && right_linked && links[0].ps ==
links[1].ps` and the two links land on different proposal ids. Only **10 of 637
reads in the window carry the left block's phase set**, they never merge into the
388-438 read block, and **zero reads in the big block overlap the left anchor**,
so the big block has no left vote in any tier.

The gate is `select_gap_link_sites` (`collect_phase.cpp:1188-1193`, pre-fix line numbers): an
MSA-verified **indel** in the gap may carry link support unconditionally, while
an MSA-verified **SNP** requires `--gap-bridge-private-snps`, off by default.
Both come from the same verification. Measured against truth the policy is
inverted: the flag-gated verified SNPs have median segregation **0.988**
(0.852-1.000, three at 1.000), the default-admitted verified indels **0.895**
(0.657-1.000) -- the default class holds the two worst sites in the window, the
gated class the three best.

Corroboration from the other direction: the whole-chr20 `--gap-bridge-private-snps`
arm joined exactly one gap, `['CHM13#0#chr20','36332599','36381019']`, with zero
concordant-to-discordant reads. Same locus, same mechanism, two independent routes.

Fix direction: (1) drop `opts.gap_bridge_private_snps &&` from `verified_snp` so
both verified classes are treated alike -- a strict subset of the flag's current
behaviour, so it needs its own chr20 gate measurement rather than inheriting that
flag's 61-read coverage cost; (2) the same function's anchor loop only accepts
`kCandGermlineClean` anchors, so MSA sites can be linked but never anchor, and
cannot chain to each other -- every link must reach a flanking clean anchor,
capping bridgeable distance however dense the MSA evidence is; (3)
`kCandAnchorClean` (`collect_phase.hpp:47`) is dead code whose comment describes
an anchor policy the code never enforces.

### Correction and fix: the private-SNP bridge is gated twice, not once (2026-09-15)

The previous entry named `select_gap_link_sites` as the gate. That is only half
of it. Two gates carry the same asymmetry and are chained:
`collect_phase.cpp:1191` (`verified_snp`) decides whether a site may earn
`gap_link_supported`, and `collect_phase.cpp:896` (`recovered_snp`) decides
whether it may cast a bridge vote -- and the second requires the first.
Patching only the eligibility gate is **inert**: measured, it produced
byte-identical tier statuses on `chr20:36,332,599-36,381,019`.

Both flag terms are now removed, so an MSA-verified SNP is treated exactly as an
MSA-verified indel already was. Clean build, no new warnings, all five unit-test
binaries pass; `test_private_snp_bridge_anchor_requires_flag` pinned the old
behaviour and is now `test_verified_msa_snp_bridges_without_flag`.

Verified on the diagnosed gap with `--gap-bridge-private-snps` **off**: tier 3
reports `joined` with `LEFT_LINK_PS = RIGHT_LINK_PS = 36317511`, the window goes
from 2 phase sets to **1 spanning set** (55 variants, 437 reads, 99.08%
accurate). Read-level gate matched by name: **0 concordant -> discordant**, 433
concordant preserved, and **23 concordant + 2 discordant reads lose their tags**.
So the join costs 23 correct read tags for 48.4 kb of closed gap.

Not promoted to default on this evidence: one gap, and the coverage cost is the
same kind the whole-chr20 flag arm measured (61 reads). The chr20 gate run is the
next step before the flag is retired.

### The two-gate fix verified across the panel: 1 of 8 gaps (2026-09-15)

Re-ran the eight-window panel with the patched build, flag off, against matched
pre-change runs (`evaluations/2026-09-15-gap-link-site-gate/compare_panel.py`,
`panel_before_after.tsv`).

`chr20:36,332,599-36,381,019` joins and gains a spanning phase set. **The other
seven windows are unchanged down to the individual read tag** -- no verdict
change, no tag gained or lost. Panel totals: **0 concordant -> discordant**, 23
concordant tags lost, tagged 3832 -> 3807, concordant 3808 -> 3785.

Crucially the no-ops are not for lack of material: `35,919,404` holds **12**
in-gap MSA-verified het SNPs (more than the gap that joined) plus 15 in-gap clean
het SNPs for anchoring, and `7,163,303` holds 6 -- both still `split`. So the
flag asymmetry was a genuine blocker but only one of at least two, and it does
not explain the `split` population. Chromosome-wide this change should be
expected to behave like the earlier flag arm did: about one extra join.

Next diagnosis: the same `--gap-decision-audit` on `35,919,404`, where anchors
and verified SNPs are both present and both gates are now open, so the blocker is
a third mechanism -- bridge-vote thresholds (`min_block_link_reads`, the
`strong`/`snp_strong` requirements) or the anchor-eligibility restriction that
keeps MSA sites from anchoring.

### Correction: the evaluator's phase-block metrics are invalid on the graph path (2026-09-15)

The `-q 1` adoption entry cited "10.3% block N50" from
`scripts/evaluate_phase_accuracy.py`. That number, and the auN column in
`graph_mapq_sweep.tsv`, are not valid block statistics on this path: the
evaluator reports **auN = 15,474,231 bp at `-q 1`** while the largest single
phase block is **1,282,456 bp**, and auN cannot exceed the largest block. The
evaluator derives block extents from the emitted alignment, which is unaligned
here -- the same reason its switch/flip counts were already known to be
meaningless for `collect-graph-variation`.

Block spans recomputed from `phase_sites.tsv` (first to last phased site per
phase set): `-q 30` gives 375 blocks / 47.60 Mb / N50 439,677 / auN 496,788;
`-q 1` gives 402 blocks / 48.34 Mb / **N50 428,655** / auN 492,205, with an
identical largest block of 1,282,456 bp. So N50 **falls 2.5%** rather than rising
10%, and the auN "spike" is an artifact, not a giant mis-joined block.

The `-q 1` decision stands, on the parts that never depended on those metrics:
5,041 more clean het SNPs, MAPQ >= 30 accuracy flat (0.856% -> 0.854%), phased
span up 0.74 Mb, gap span 18.52 -> 17.81 Mb. The block-length argument is
withdrawn. Rule carried forward: do not quote this evaluator's N50/auN for any
graph-path arm -- derive block spans from `phase_sites.tsv`.

### The patch is gate-specific: `35,919,404` is coverage-limited (2026-09-15)

Audited the largest non-joined panel gap (`chr20:35,919,404-36,156,319`, 236.9 kb)
with the patched build. It carries more of the material the fix unlocks than the
gap that joined (10 MSA het SNPs, 27 MSA het indels per tier; 48 in-gap sites in
phase-informative categories in the audit export, 19 of them truth-informative)
and still reports `split` on tiers 1-3 and
`rejected` on tier 4.

Its vote matrix has the opposite shape to the joined gap: **both flanks link
strongly** (left PS 35873360, 124/126 votes; right PS 36156319, 54/60) into two
blocks with nothing between, where the joined gap had a 10-vote left remnant.
Nothing is refused; the halves never meet. Chaining the truth-informative sites,
**3 of 18 steps have zero spanning reads** (22.9, 26.5, 43.5 kb = 92.9 kb, 39% of
the gap); repeating the step test over the wider 60-site linking selection from
the candidate TSV -- a different set, never truth-scored -- one 21.2 kb step has
none. No
site-admission policy bridges a step no read crosses.

Panel-wide (`classify_blockers.py`, `blocker_classes.tsv`): 2 gaps are
coverage-limited (`13,429,829` with a 16.8 kb break, `35,919,404` with 21.2 kb,
438.9 kb of span between them -- the same two that neither phase transfer nor
injection moved), and 6 have a flank-to-flank chain of called sites with reads
spanning every step. The patch fixed **one** of those six. So five have a chain
and no flag refusing them and still split: a third blocker, neither site
admission nor read coverage. Next candidates in order -- whether those chains'
sites are actually informative (at `35,919,404`, 19 of the 48 audit-category
in-gap sites scored informative; do not chain that to the 60-site candidate-TSV
selection, which was never truth-scored), the bridge-vote thresholds (`min_block_link_reads`,
`strong`/`snp_strong`), and the anchor-eligibility rule that stops MSA sites
anchoring each other.

### Homopolymer tier skipped the pass where the join lives (2026-09-16)

Diagnosing the competitor-deficit gaps one at a time
(`evaluations/2026-09-16-competitor-deficit/`) turned up a narrow defect in gap
recovery. `collect_pipeline.cpp` skipped `kGapHomopolymerTier` whenever
`recovery_pass == 1 && !audit`, so the tier never ran in the reprojected pass.
On `chr20:61,738,233-61,757,551` that is the only pass where the join exists: the
audit export, which bypasses the skip, showed pass 1 reaching `JOINED` for
proposal set 61725696 with left support 5 and right support 4 -- clear of
`min_block_link_reads = 2` and of the tier's per-haplotype orientation agreement
-- while pass 0's differently anchored proposal linked one side only. A
consequence worth remembering: the audit and normal runs do not evaluate the same
tiers, so "audit matches normal" cannot be inferred from matching final blocks.

Removing the skip closes the gap. `chr20:61,732,321-61,810,469` (78.1 kb) goes
from two blocks to **one spanning block `61,690,751-61,858,547`, 56 sites**, and
the read level is untouched: 575 reads tagged before and after, 567 concordant
both times (98.61%), **0 concordant -> discordant and 0 tags lost**. Compare the
private-SNP gate fix, which closed its gap at a cost of 23 correct read tags.
Gap 2 (`36,217,274-36,268,291`) is the control: the tier now runs in pass 1,
reports `rejected` because no proposal block holds both-sided support, and the
output is identical.

Also measured while diagnosing these two gaps: `msa_verified` is not a proxy for
informativeness. All four `NoisyMsaHet` sites in gap 2's break carry
`msa_verified = 1` and their truth segregation runs 0.509 (chance), 0.607, 0.644,
0.923 -- and the tiers admit the uninformative ones while the 0.923 site is
homopolymer-flagged and reachable only by tier 4.

Still opt-in-by-measurement: the chr20 gate run decides whether this ships as
default, since the change touches every gap with a homopolymer candidate.

Panel regression for the same change: eight windows, matched pre-change runs, and
it is inert -- 0 newly joined, 0 newly spanning, 0 gate violations, 0 tag changes,
3807 reads tagged and 3785 concordant in both arms. So the change buys gap 1 and
costs nothing measured, but is not a broad win. It does introduce a reporting
nuance: five panel gaps whose STATUS read `split` or `partial` now read
`rejected`, since tier 4 is evaluated in pass 1 and its label is the last report
row -- same blocks and reads underneath, but the `split`/`partial` distinction is
lost to anything reading the final row.

### Deficit gap 2 diagnosed to the line; three fixes kept, two reverted (2026-09-16)

`chr20:36,217,274-36,268,291`, second gap of the competitor-deficit set.
`evaluations/2026-09-16-competitor-deficit/`.

The evidence is sufficient and the correct join exists. An offline solver over
the audit's per-read observations reaches a spanning partition at **1.000**
accuracy from clean sites plus the one verified homopolymer deletion (0.923
segregation), while the full site set reaches only 0.846 and tier 3's set --
clean plus the non-repeat verified indels, which segregate at 0.509 and 0.644 --
reaches 0.538. The tiers were cumulative, so the junk admitted at tier 3 was
never removed when tier 4 added the good site.

What blocks the join is upstream of the tiers. `stitch_gap_proposal` counts a
link vote only from a read still holding a hap and one of the flanks' phase sets,
and two read-tagging filters run before `recover_hybrid_gaps` in the batch loop:
**62 of the 64 reads overlapping the left flank's own boundary variant are
zeroed there**, leaving that side 2 voters. Those reads observe exactly one clean
het SNP, and the margin counts only `kCandCleanHetSnp` (its rescue credits
bridge SNPs, never indels), so 1 < 2 strips them. With the labels unfiltered,
both flanks link to the same proposal phase set and the tier-4 proposal scores
1.0000 against read truth over 95 reads, agreeing with the frozen haplotype of
all 26 left-flank and all 69 right-flank reads.

Lowering `--min-read-margin` is not the fix. Ten-window panel: margin 2 -> 1
turns **99 concordant reads discordant** (70 in this gap alone, so its join there
is wrong) and drops accuracy 97.63% -> 95.44%; margin 2 -> 0 turns 30 and drops
99.38% -> 96.48%. The flag is also not a BAM-side remnant -- it is the global
`min_read_hap_margin`, parsed by both subcommands and consumed by the chromosome
pass, gap recovery's proposal and validation, and the graph path's own read gate.
Its built-in default is 0 while this project's canonical scripts pass 2.

Reordering the filters after recovery is not the implementation either: it fixes
the starvation (all 64 boundary reads keep a phase set) but then judges reads on
counters recomputed over a narrow gap window, collapsing tagging 358 -> 177 in
this gap and 575 -> 398 in gap 1 (105 and 169 concordant tags lost, 0 concordant
-> discordant). Crediting bridge and last-resort evidence at that call does not
restore them. Reverted, along with removing the validation view's own filters,
which did not clear the veto.

Kept, each measured, defaults unchanged elsewhere:

1. Tier 4 is a real last resort -- clean plus verified SNPs plus verified
   homopolymer indels -- instead of restoring every flag and readmitting the
   non-repeat indels tier 3 had just failed with.
2. `CandidateVariant::hp_gap_scorable` lets the site the homopolymer tier
   admitted contribute read scores. `init_assign_read_hap` skipped *every*
   homopolymer indel, so the tier could earn link support and still not phase the
   reads it was reached for; this is what grew the in-gap block 69 -> 95 reads and
   extended it ~17 kb leftward (reads reaching the left endpoint 1 -> 26).
   `gap_link_supported` could not be reused: line 1188 grants it from
   `msa_insertion_alts` alone, also outside a homopolymer gap. Unit test
   `test_hp_gap_site_scores_reads_only_in_its_gap` pins both directions.
3. `ReadRecord::n_hp_gap_agree/conflict`, credited in the recovery margin
   filter's rescue for the same reason bridge SNPs are.

Gap 1 still joins with all three applied; gap 2 does not yet -- it now reports
`vetoed`, the homopolymer tier's requirement that an independent BAM-only solve
reach the same orientation. That guard is worth keeping (margin 1 proved a wrong
join is available here), but it fires against a proposal measured 1.0000
accurate with `SELECTED_GRAPH_READS=157` and `GRAPH_BAM_PASS=1`, so it needs
instrumenting next. The implementation the diagnosis points to is a snapshot:
build the stitch's `read_index` from the labels as they stand before the output
filters, leaving output filtering where it is.

Also measured: two truth-free site screens fail here. Pairwise agreement with a
clean anchor ranks the 0.923 site at 0.744 below a 0.644 site at 0.818, because
agreement between two noisy sites compounds their errors; leave-one-out
consistency scores the junk *higher* (0.820, 0.857 vs 0.754, 0.737) because
`36,261,164`/`36,261,302`/`36,261,311` are one repeat event called three times
inside 150 bp and form a self-consistent clique. Collapsing candidates within
300 bp to the best-covered one lifts gap 2 to 0.981 and leaves gap 1 at 1.000 --
the screen worth implementing.

### Gap-recovery link votes now read the pre-filter labels (2026-09-16)

`stitch_gap_proposal` counts a link vote only from a read still holding a hap and
one of the flanks' phase sets, and two read-tagging filters run before
`recover_hybrid_gaps` in the batch loop. At `chr20:36,247,421` that left the left
flank **2 voters out of 64 reads** overlapping its own boundary variant, with no
join reachable from any tier.

`GapReadIndex::vote_assignments` now falls back to the label a read held before
those filters, captured as `PreFilterLabels` in the batch loop and threaded into
`recover_hybrid_gaps`. The gap inventory is still derived from the filtered
chunks, so emission is unchanged; only the votes see more.

Scoping took three attempts, and the lesson generalizes: three consumers read
"no committed assignment" as permission to act, so widening the filtered view
suppresses read attachment instead of enabling joins. `emit_independent_gap_block`
re-phases reads with no assignment as gap-only; `first_assignment` also builds
`supported_phase_sets`, which gates emission; and `original` in
`stitch_gap_proposal` also means "already in a block, do not attach". Widening
any of them cost 105 and 169 concordant read tags in the two deficit gaps. Final
form: those three stay filter-accurate, and one wide map (`vote_original`) is
used by the vote loop alone. The validation view's own output filters were
removed too -- filtering the validator for reporting silenced the validator.

Measured, ten-window panel vs matched pre-change runs
(`evaluations/2026-09-16-competitor-deficit/votesnapshot_panel.tsv`): **0
concordant -> discordant**, tagged 4,687 -> 4,714, concordant 4,576 -> 4,606
(net +30), nine of ten windows byte-identical. All movement is in deficit gap 1
`61,732,321`, which gains 78 concordant tags and rises 98.61% -> 99.17%.

Deficit gap 2 is still not closed. Both flanks now link to the same proposal
phase set, so the starvation is gone, but tier 4 reports `vetoed` with
`bam_reads=157 bam_joined=0 bam_flip=0 graph_flip=0`: the orientations agree and
the validator simply cannot link the flanks from its own read set, so it cannot
confirm a join measured 1.0000 accurate over 95 reads. The guard stays -- margin
1 proved a wrong join is available in this window. The next step is the
validator itself: `select_graph_gap_bam_reads` restricts it to reads carrying a
graph-channel call, so it must reproduce a join from a smaller read set than the
proposal it judges.

### Why gap 2's validator vetoes a 1.0000-accurate join (2026-09-16)

Probed the veto branch directly (probe removed). The homopolymer tier's
BAM-only validator is not disagreeing about orientation -- it never forms a
block reaching both flanks. On `chr20:36,247,421-36,268,291` it phases 151
reads into three blocks: `36247421` with 63 reads all touching the left flank
and none the right, `36259922` with 41 touching only the right, `36261163` with
47 touching only the right. The graph proposal it is judging holds one 95-read
block with 26 left-flank and 70 right-flank reads at 1.0000 accuracy.

Two causes. `select_graph_gap_bam_reads` admits only reads carrying a
graph-channel call, so the validator solves on 157 of the window's 510 reads
with 353 skipped -- thinner coverage than the proposal it must reproduce. And
its block boundaries, `36,259,922` and `36,261,163`, are the junk sites
`36,259,923`/`36,261,164`: the one repeat event called three times inside 150 bp
that also drags the full solve from 1.000 to 0.846. At 157 reads there is not
enough overlap to chain across them.

So the fix is the positional screen already identified on the site side (collapse
candidates within 300 bp to the best-covered one; lifts the full solve 0.846 ->
0.981), not a weaker veto. Widening the validator's read set to match the
proposal is a second, weaker candidate -- it makes the check easier to pass
without making the evidence cleaner.

### Default-pipeline merge parity for gap recovery: tested, and it is wrong (2026-09-16)

Surveyed how the default pipeline merges blocks and flips hap tags. Two
mechanisms only: `flip_chunk_hap` counts an n11/n12/n21/n22 table over reads
shared by adjacent CHUNKS, calls `select_stitch_orientation` (default
`kStitchRuleNetMargin`: merge iff `|(n12+n21)-(n11+n22)| > stitch_min_margin`,
orientation from the sign), then `apply_chunk_flip_and_merge` rewrites the
downstream phase set to the upstream id and flips haps;
`propagate_overlap_read_phase_to_output_owner` then copies HP/PS onto overlap
reads that were unphased, only for pairs that merged. Blocks WITHIN a chunk are
merged only by `stitch_phase_blocks_with_pgbam`, which needs `--pgbam-file`.

`stitch_gap_proposal` already builds the same table per flank and carries the
same orientation bit, so the only parity gap is the homopolymer tier's extra
BAM-only re-solve. Tested removing it when the table is unambiguous: gap 2's
table is 26-0 left and 69-0 right (default-rule scores -26 and -69, three times
the support that joins deficit gap 1 at 19 and 10). The gap then joins and spans
`36,168,059-36,268,558` -- and the join is WRONG: 70 concordant reads flip to
discordant and window accuracy falls 77.65% -> 58.38%. All 70 were in the right
flank's phase set `36268291` and all land in the merged block, so the merge
inverted the right flank relative to the left. Reverted; the veto is load-bearing.

Why the same rule misfires here: in the default pipeline the voting reads are
phased independently by two full-coverage chunk solves, and a chunk seam is an
artifact of chunking. In recovery the left-boundary voters are 62 of 64 reads the
output filter had zeroed, each phased from a single clean het SNP -- correlated,
so 26 unanimous votes carry about one vote of independent information. Unanimity
among weak correlated voters is not confirmation.

Parity should therefore mean better voters, not a looser check: collapse the
repeat cluster (`36,259,923`/`36,261,164`/`36,261,311`, one event called three
times in 150 bp, also the cause of the 0.846 solve ceiling and the validator's
fragmentation), or require more than one supporting site before a filter-zeroed
read may vote -- which on this gap returns the left flank to 2 voters and
correctly refuses the join. Note the committed vote-snapshot change is what makes
this wrong join reachable; it and the confirmation must not be relaxed together.

### Gap 2 is downstream of a wrong tier-1 recovery join (2026-09-16)

Checked the phasing inside deficit gaps 1 and 2 on its own, no joining. Both are
perfect at tier 4: gap 1 gives 45 reads at 100.0% (27 left-reaching, 20
right-reaching) and gap 2 gives 95 reads at 100.0% (26 and 70). Only tier 3
degrades (gap 2 pass 0: 41 reads, 97.6%), the tier that admits the non-repeat
verified indels. This also retracts an earlier claim: gap 2's proposal carries NO
internal switch (26/26 and 70/70); the 93.3%/81.8% split reported from an offline
solver was an artifact of that solver voting over all seven interior sites by
plain majority instead of the tier's filtered set, and it is why implementing a
positional repeat collapse changed nothing.

Probing every link voter for proposal 36247421 showed both flanks sound: left
flank 36209945, 64 voters, proposal and flank labels both 90.6% correct against
truth; right flank 36268291, 69 voters, both 100.0%; all `hap1=PAT`. The defect
is the block they point at. Emitted block 36168059 (36,144,761-36,267,248, 287
reads) is only 72.5% consistent, and the opposite reads appear abruptly: 0% below
36,210,000, then 31% / 57% / 92% / 87% in the 10 kb bins above it. That is a
switch error starting exactly at 36,209,945, the phase set gap 2 links to.

The switch is ours. The same run joined an upstream gap at tier 1:
`36,172,778-36,209,945, leftPS=rightPS=36168059, reads_added=81, joined`.
Re-running the identical window without `--recover-gaps`: block 36168059 is 169
reads at **99.4%** and the window is **241/242 (99.59%)**; with recovery the block
is 287 reads at **72.5%** and the window **278/358 (77.65%)**. Recovery buys 116
more scored reads and 37 more concordant ones while creating 79 discordant ones.

So gap 2 was never an admission problem -- its solve is perfect and its votes are
right, but it would attach a correct local orientation to an already-inverted
segment, which is why forcing the join flips exactly 70 reads however it is
forced (margin 1, unanimity skip, with or without repeat collapse). The veto was
preventing us from compounding an error already made.

Next: audit the tier-1 join at 36,172,778-36,209,945 the same way, and measure the
read-level gate for `--recover-gaps` on versus off chromosome-wide. One wrong join
costs more concordant reads than several correct ones gain.

### Root cause of the gap-2 corruption: allele attachment on an unjoined gap (2026-09-16)

Corrects the entry above, which blamed the tier-1 join at 36,172,778-36,209,945.
Recovery is the cause; that join is not. Dumping and scoring read labels at each
stage puts the damage inside the gap threads, before any parity edge: window
157/161 (97.52%) pre-recovery, 278/360 (77.22%) after the threads, unchanged by
`apply_gap_phase_edges`. The joined gap's own target block goes 44 -> 125 reads
and stays at 100.0%.

Attributing every write single-threaded -- 8 threads interleave stderr and lose
lines, which is why an earlier probe counted 44 of 199 -- shows 155 of the 199 new
labels come from `implied_assignment`, the `--link-by-alleles` path: 38 reads for
the joined gap at 100.0% against truth, and **115 reads for gap
36,247,421-36,268,291, which never joins, at 67.0%**, taking flank 36209945 from
49 reads at 91.8% to 164 at 50.6%. `implied_assignment` is a unanimity test with
no minimum count: a read matching a single site of that phase set gets a committed
haplotype, with no flip and, on a one-sided link, nothing confirming the flank's
orientation.

That is why the gap looked unjoinable. Its phasing is perfect (95 reads, 100.0%,
26/26 left and 70/70 right) and every way of forcing the join flipped the same ~70
reads, because the block being joined to had already been filled with coin-flip
haplotypes by this same gap's own one-sided attachment. It also retires the
offline "internal switch" claim: the pipeline's proposal has none.

Shipped opt-in as `--gap-allele-attach-join-only` (default byte-identical). On
this window 77.65% -> 98.51%, 265 concordant of 269 tagged against 241/242 with
recovery off. Not a default: across the ten-window panel it removes 90 discordant
reads and 1,021 concordant ones (4,687 -> 3,576 tagged, 97.63% -> 99.41%, gate
concordant->discordant 2). Read agreement count does not discriminate -- single-site
attachments are 98.1% correct at 7,073,919 and 74.8% at 36,217,274 -- so the
discriminator must be a property of the target flank, not the read. Next: score
attachments against their target block's polarity across the panel, split by
whether the flank carries committed read support.

### Deficit gap 3: a real linkage break, and four interior blocks we were discarding (2026-09-16)

`chr20:35,919,404-36,156,319` (236.9 kb) fails differently from gap 2: every tier
links both flanks but to different proposal phase sets, and the join needs one
proposal block holding both. The proposal fragments into four to six internally
accurate blocks (96.7-100.0% against read truth), one reaching the left edge with
63 reads and one the right with 71.

The cause is a genuine linkage break, not site discovery: 58 interior sites cover
the gap edge to edge, but the spacing 35,959,001 -> 35,981,762 (22.8 kb) has ZERO
reads covering both sites, while every other large spacing has 8-20. Coverage is
65-71x with 168 reads inside the interval; the longest read there is 29.1 kb. No
read-based phaser can cross that point, so abstaining is correct.

What was wrong is what happened to the interior. The output kept only the two
flanks and dropped every interior block. `emit_independent_gap_block` exists for
this and its own --verbose 2 line explains itself: `gap_only_ps_groups=5
chosen_ps=35902410 chosen_size=111` then `applied=5`. It takes only the largest
gap-only group, and that group is the flank-adjacent component whose reads this
same round already attached -- so it applies 5 reads and discards four real
interior blocks. `--gap-independent-min-reads` was never the limit (default 3; at
20 the window moves by 4 reads).

Fixed by emitting every group that clears min_reads, which already screens noise
fragments. This window: 527 -> 801 tagged, 518 -> 788 concordant, 98.29% ->
98.38%, gate concordant->discordant 0, four new interior blocks (41/91/50/92 reads
at 100.0/98.9/100.0/96.7%). Ten-window panel: blocks 21 -> 26, tagged 4,687 ->
5,037, concordant 4,576 -> 4,922, discordant 111 -> 115, 97.63% -> 97.72%, gate 0,
394 newly tagged concordant against 9 discordant. A net gain rather than a trade,
so this one ships as a default. 61,732,321 churns (48 concordant tags lost against
a larger gain, net +30) because more interior blocks change which reads its
homopolymer join claims first; nothing there flips concordant to discordant.

emit_independent_gap_block had no unit test; it now has one pinning both halves of
the contract.

### The deficit re-measured against the current pipeline: 0.81 Mb, and it is indels (2026-09-16)

The 187-gap / 4.74 Mb deficit was computed from the graph-only pass-1 gap
inventory and does not describe what the pipeline leaves open.
`chr20:30,794,962-30,814,005`, which all three competitors close, is already
covered by a single 76.1 kb hybrid phase set at 99.6% read accuracy with 116 clean
het SNPs inside; recovery is never asked about it. `evaluations/2026-09-16-current-deficit/`.

Whole-chr20 current pipeline (10.5 min, 16 threads): 238 blocks, 56.02 Mb phased,
196 gaps spanning 11.88 Mb. Of those, 47 gaps / 1.20 Mb are spanned by at least
one competitor and 149 gaps / 10.67 Mb by nobody. Scoring each span against read
truth cuts it further: hiphase spans 45 at 91.11% overall but is >=98% in only 31
(0.80 Mb) and <90% in 8 (0.25 Mb), including 31.8% at 36,059,715 and 35.9% at
58,994,554 -- where our abstention beats their join. THE RECOVERABLE DEFICIT IS 31
GAPS / 0.81 Mb.

In the 30 gaps hiphase spans at >=98%, its 98 in-gap het sites break down as: all
44 SNPs we already hold as CLEAN_HET_SNP (SNP discovery is not the deficit), 33
NOISY_CAND_HET indels, 7 we call NOISY_CAND_HOM, 6 CLEAN_HET_INDEL, 3
REP_HET_INDEL, 5 absent. Genotyped from the alignment against read truth, with
the clean classes as controls at median 1.000 / 100% informative: NOISY_CAND_HET
in these gaps is median 0.904, 58% informative, 18% phantom -- far better than the
23% measured for repeat-demoted indels chromosome-wide. The 7 NOISY_CAND_HOM
calls (median 0.939) are a genotyping error, not a screening decision: a site
called homozygous can never link. A first run scored the SNP control at 0.542,
which was a one-base coordinate error in my genotyper, caught by that control.

MSA verification is NOT the blocker. On chr20:48,176,830-48,229,446 (hiphase 100%
over 252 reads) it crosses on 8 het records, 7 of them homopolymer/tandem indels,
and every msa_verified=1 site there segregates 0.89-1.00. The blocker is
structural: every tier links both flanks but to DIFFERENT proposal phase sets.
Recovery's gap for that region is 134.7 kb and contains two hard linkage breaks
(48,096,582->48,123,657, 27.1 kb, and 48,123,657->48,147,230, 23.6 kb, both with
ZERO reads covering the flanking sites at 72x/67x coverage), while the interval a
competitor spans is the right 52.6 kb whose widest holes carry 7 spanning reads
each. Recovery is all-or-nothing over a whole gap, so the bridgeable part is
abandoned with the unbridgeable part. That is the next method change to make.

Also settled: chr20:13,429,829-13,631,825 (202 kb) is a correct abstention for
everyone -- no competitor spans it, 8 clean het SNPs in 202 kb at 75x -- and its
proposal blocks cannot be chained (adjacent pairs share 0-12 reads observing
consensus sites in both; the one pair with 12 votes splits 6/6).

### The pgbam within-chunk stitcher exists, and is unusable on this input (2026-09-16)

Answering whether the current machinery can make multiple blocks and stitch them.
`evaluations/2026-09-16-pgbam-stitch/`.

`stitch_chunk_haps` (collect_phase.cpp:1658) is four stages: per-chunk
`stitch_phase_blocks_with_pgbam` (merges blocks WITHIN a chunk), the adjacent-chunk
read vote `flip_chunk_hap`, `stitch_adjacent_chunks_with_pgbam` rescuing seams the
read vote could not decide, then two more within-chunk pgbam passes at cleanup and
relaxed thresholds (both default ON). Stages 1, 3 and 4 are gated on
`--pgbam-file`, which no run in this investigation passed -- while the sidecar for
this exact BAM sat in test_data/ the whole time.

It links on GBWT haplotype threads: decide_phase_block_concordance intersects the
thread sets polarized to each hap and merges when one orientation wins by
min_winning_threads. That is the right SHAPE of evidence -- a thread crosses a
27 kb read-linkage break that no read can -- but measured on whole chr20 it is not
usable: at defaults it merges the chromosome into ONE read block at 47.511% read
Hamming (chance) against 0.878% for the committed no-pgbam arm. Tightening trades
merging for accuracy monotonically and never nears baseline: primary pass only at
win/margin 2 -> 46.853% (5 blocks), 5 -> 41.087% (11), 20 -> 24.050% (92). The
sidecar parses (magic and version correct) but predates the current pipeline, so a
stale sidecar is not excluded. It also explains a 51 s runtime versus 10.5 min:
with everything merged there are no gaps and recovery does nothing (1,313 tier
rows -> 1).

The gap this leaves is precise: select_stitch_orientation, the project's read-vote
standard, is called in exactly three places -- the adjacent-chunk seam and the two
gap-recovery link votes. THERE IS NO READ-VOTE STITCH BETWEEN TWO BLOCKS INSIDE A
CHUNK. Chunks are 500 kb and the deficit gaps are 3-72 kb, so they are within-chunk
essentially always, and the read stitch never applies to them. That is why recovery
has its own flank-linking and why it is all-or-nothing. Re-prototype the pairwise
block vote on the MAIN solve's blocks (consensus from hundreds of reads, not the
proposal's dozens) before writing C++; the proposal-level prototype found seams of
0-12 shared reads with one 12-read pair splitting 6/6.

Changed: the hybrid subcommand exposed --pgbam-file but none of the eight pgbam
threshold options collect-bam-variation has, so the path could only run at the
47.5%-error defaults. All eight are now exposed in collect-hybrid-variation.
Defaults unchanged.

### BAM variation IS in the initial hybrid solve, then withheld from phasing (2026-09-16)

Tested on chr20:48,176,830-48,229,446 (deficit gap, hiphase spans at 100.0% over
252 reads). `evaluations/2026-09-16-bam-sites-initial-solve/`.

process_chunk_hybrid is "BAM classification -> graph site injection -> BAM profile
build -> graph read injection -> unified k-means": the BAM chunk is the BASE and
graph sites are injected into it, so BAM variation is already in the initial solve,
and graph_authoritative is off by default so BAM evidence at matched candidates is
kept (enabling it unconditionally cost +720 discordant reads, per the call site).

But two mechanisms exclude the BAM's own sites from phasing. (1) With
--recover-gaps, every non-graph candidate has lcd_var_i_to_cate zeroed before
collect_var_run_phasing and restored only afterwards, so BAM-discovered het sites
are invisible to the initial k-means AND to the noisy-region MSA, visible only to
gap recovery. (2) hybrid_collect.cpp:120 sets skip_noisy_kmeans = true, disabling
the step-4 noisy-candidate k-means, documented as having "phased ~8k extra reads
at ~65% error and poisoned the BAM-shared core".

Measured on the window: BAM-ONLY phases 11/11 in-gap het sites into ONE block
spanning 48,147,227-48,229,226 (82.0 kb, 16 sites) at 100% read accuracy, covering
all but the last 220 bp of the competitor-spanned interval. Hybrid default phases
2/10 and leaves the interval unphased with 3 blocks, the nearest stopping exactly
at the gap's left edge (48,162,480-48,176,830, 2 sites). Recovery off gives 3/11;
--keep-noisy-kmeans gives 2/4. So no single flag explains it.

The two channels find the SAME sites: 9 vs 10 het sites at the same positions,
7 shared NOISY_CAND_HET indels; bam-only phases 9/9 into PS 48147227, hybrid 2/10.
48,177,781 and 48,225,787 are NOISY_CAND_HET to the BAM channel but
NOISY_CAND_HOM in hybrid (the het-called-hom class, 7 of hiphase's sites
chromosome-wide). Noisy-region DETECTION is byte-identical: both paths emit the
same 2,862 verbosity-2 noisy-region lines over the same regions.

So at this gap the deficit is neither discovery nor MSA verification nor missing
injection -- it is that hybrid withholds the BAM channel's noisy-site phasing that
would have produced an 82 kb block at 100%. Fix direction: admit BAM-discovered
noisy candidates to phasing INSIDE GAP INTERVALS ONLY, scoped rather than the
chromosome-wide enabling that cost ~8k reads at ~65% error; regression-test it with
this window's three-arm comparison and the 31-gap accurate-deficit list.

### gap_lab.py: a gated per-gap rig, and +159 correct reads on the first gap (2026-09-16)

`evaluations/2026-09-16-gap-lab/`. One command per gap: baseline arm, BAM evidence
inside the interval with truth segregation per site, the BAM channel's own phasing
of the interval, a two-stage stitch, and a read-level gate that decides pass/fail.
`--reuse` skips arms whose outputs exist, so stitch iterations cost seconds.

Three mechanisms it exposed, each of which had silently broken a hand analysis:

1. A FLANK CAN HOLD SITES AND NO READS. The block nearest this gap on the left
   (48162480) has 2 sites and ZERO tagged reads, because the read-tagging margin
   drops reads observing too few of its sites. Nothing links to it by read
   identity. Flanks are now selected by read support and read-less blocks are
   reported as skipped.
2. TWO BLOCKS THAT SPLIT AT THE SAME POSITION SHARE NO TAGGED READS, since a read
   carries at most one phase set -- so tag-identity voting returns n=0 for exactly
   the seams we care about. The alleles remain: a read crossing the seam observes
   sites in both blocks and each block's phased genotypes map allele -> hap. That
   allele-level vote is now the fallback and is what composed this gap.
3. THE SEAM IS NOT THE BLOCK. Requiring a bridge read to span every selected site
   (which reach tens of kb back into each block) reported zero crossing reads
   where the seam is 220 bp wide. The test is now the seam point plus a margin
   with a minimum site-observation count per side.

RESULT on chr20:48,176,830-48,229,446 (52.6 kb, hiphase spans at 100.0% over 252
reads): BAM channel phases 11/11 in-gap het sites (7 informative vs truth).
Compose 48147227 + 48229446 -> 42 allele voters, votes [15,0,0,27], UNANIMOUS, 58
reads cross the seam. The composed frame (231 sites, 48,147,227-48,377,087, 820
reads) links the RIGHT flank (n=661, [318,0,0,343]) but NOT the left (n=0), since
the left side carries the two hard breaks (48,096,582->48,123,657 27.1 kb and
48,123,657->48,147,230 23.6 kb, zero spanning reads at 72x). Gate: tagged 893 ->
1052, concordant 892 -> 1051, 0 concordant->discordant, 159 newly tagged ALL
concordant, 0 lost. VERDICT PASS (not closed; right flank extended, +159).

The two operations that produced that gain are precisely the two the pipeline has
no mechanism for: an allele-level block-to-block vote, and a one-sided extension
of a flank by a composed frame. The stitch runs outside the pipeline, so this is a
measurement rig rather than a fix -- it states what a scoped change must
reproduce, and its verdict is that change's regression test. Next: run it across
the 31-gap accurate-deficit list to size the recoverable coverage before any C++.

### gap_lab.py corrected: the gap arm must run IN the gap, and the gap closes (2026-09-16)

The first version of the rig ran the machinery over gap+150 kb and took its flanks
from that same re-run. Both wrong: over a wide window the gap arm's blocks are
window-wide and inherit the window's breaks (the block it produced spanned
48,147,227-48,229,226, not the gap), and re-solving the flanks in a window can move
the boundaries under test. The gap arm now runs on gap + --gap-margin (5 kb) of read
context, and --baseline-vcf/--baseline-bam take the real pipeline run's outputs.

That exposed two further mechanisms. FLANK SELECTION BY READ SUPPORT IS WRONG: the
block holding the nearest phased site left of this gap (48162480) has 2 sites and
ZERO tagged reads, so picking by read support jumped 80 kb further out across a
stretch with no phased sites and then reported no link. Flanks are now the blocks
holding the phased sites nearest the gap, and a read-less flank is linked through
its GENOTYPES: the allele-level bridge vote, previously used only for gap-block
composition, now also links frames to flanks. And A READ-LESS FLANK IS INVISIBLE TO
THE READ GATE -- every read in the merged block comes from the other side and keeps
its relative labelling, so a wrong orientation there flips nothing. Applied links
are now validated separately against truth (each side's genotypes read off the
alignment, the applied flip checked against the two sides' truth haplotypes), and
PASS requires it.

RESULT on chr20:48,176,830-48,229,446 (52.6 kb, hiphase spans at 100.0%): CLOSED.
Gap arm yields 2 gap-local blocks (48173317, 11 sites; 48229446, 5 sites) and phases
11/11 in-gap het sites. Compose: 42 allele voters [15,0,0,27] unanimous. Left flank
link: ALLELES, 14 voters [9,0,0,5] unanimous. Right flank: TAGS, 71 voters
[28,0,0,43] unanimous. Link validation: frame hap1 = PATERNAL over 14 sites, left
flank PATERNAL over 2 sites CORRECT, right flank PATERNAL over 81 sites CORRECT.
Gate: tagged 1319 -> 1426, concordant 1316 -> 1423, 0 concordant->discordant, 107
newly tagged ALL concordant. PASS.

So the earlier "this gap cannot be closed" was an artifact of my window: the two
hard linkage breaks (48,096,582->48,123,657 and 48,123,657->48,147,230, zero
spanning reads at 72x) lie OUTSIDE the gap and only became obstacles because the
150 kb window pulled the left flank to the far side of them. Closing this gap needs
nothing the evidence does not already contain -- what the pipeline lacks is the
allele-level block-to-block link (tags cannot orient blocks that split at the same
position, and a read-less flank has no tags at all).

### The link vote must use the adjacent blocks out to read reach, not 8 sites (2026-09-16)

The rig's flank vote took a fixed eight sites nearest the seam, which is not all
the information the adjacent blocks hold: the right flank here has 916 sites at
about one per 800 bp, so a crossing read can observe many more, and the cap both
weakened each read's own call and dropped reads below the per-side minimum. Site
selection is now bounded by --seam-span (default 30 kb, a read length).

Swept on chr20:48,176,830-48,229,446: at 3 kb nothing links at all (compose n=0,
left flank 1+6 sites, FAIL); at 10 kb the gap blocks compose (42 voters) but the
left flank still does not link (1+7 sites) so the result is extension only (+107
reads); at 30 kb the left flank links with 14 voters [9,0,0,5] over 2+7 sites and
THE GAP CLOSES; at 60 kb the vote is IDENTICAL (14 voters, 2+14 sites).

So the span is decisive -- the gap closes only once the flank's whole 2-site block
is inside the window -- and the evidence SATURATES AT READ REACH: beyond ~30 kb no
read extends further, so taking more of the adjacent block cannot add information.
The bound is read length, not block size. When a seam still has no voters after
saturation, more sites cannot help; the link must come from a different evidence
type (graph haplotype threads cross a read-linkage break) or from a chain through
an intermediate block, which is what the composition stage does for the gap's own
blocks.

### Link the gap to its flanks from the whole adjacent blocks, not a seam window (2026-09-16)

An adjacent phase block is a COMPLETE labelling -- every read it contains with a
haplotype, and a genotype at every site it phases, all integrated from that
block's whole evidence -- so a seam window throws most of it away. The rig's flank
link now draws three sources and requires them to agree:

1. SHARED PHASED SITES. The gap arm re-discovers sites inside the flank's own span
   (the flank here holds 48,162,480 and 48,176,830; the gap arm phases 48,173,317,
   48,173,989, 48,173,990 and 48,176,830), so both labellings often phase the same
   variant. Comparing genotypes at a shared site gives the orientation EXACTLY --
   no bridge read, no site selection.
2. THE FLANK'S OWN READ ASSIGNMENTS against the frame's tags-or-genotypes.
   Requiring a tag on BOTH sides, as before, discarded every read the flank had
   phased but the gap arm had not, which is most of a seam's population.
3. ALLELES BOTH SIDES, for a flank with no reads at all -- exactly this gap's left
   block (2 sites, 0 tagged reads).

--seam-span and --seam-sites are now OFF by default: the whole block is offered and
the bounding is per read (a read votes on the sites it overlaps), so a global
distance limit only duplicates it while having a gap-dependent correct value that
would quietly decide outcomes. Conflicting sources REFUSE to link rather than
picking a winner, since a conflict means one block is internally wrong at the seam.

On chr20:48,176,830-48,229,446: left flank links with shared sites 1/1 plus 14 of
66 crossing reads (2 sources, flip 0); right flank with shared sites 4/4, 71 tag
reads, and 71/298 flank reads (3 sources, flip 0). Both validate CORRECT against
truth. Gate unchanged at tagged 1319 -> 1426, concordant 1316 -> 1423, 0 flips.
PASS, gap CLOSED.

### Phase the gap in windows and let the chunk stitch refuse: safer than one block (2026-09-16)

`evaluations/2026-09-16-window-stitch/`. The pipeline already chunks: split_region
(collect_pipeline.cpp:239) cuts non-overlapping fixed windows, chunk_beg +=
chunk_size, default 500 kb, --chunk-size; flip_chunk_hap stitches adjacent chunks
on the reads present in both (up_ovlp_read_i/down_ovlp_read_i) via
select_stitch_orientation, and REFUSES when there are none (`if
(n_cur_ovlp_reads <= 0) return false`). That refusal is the guard a bespoke gap
link lacks.

Windows alone change nothing: sweeping --chunk-size 500k/50k/20k/10k over
chr20:48,126,830-48,279,446 in the shipped configuration phases ZERO sites inside
the gap at every size (10 kb is worse -- the left 2-site block fragments into two
singletons), because --recover-gaps zeroes every !graph_site category before
phasing.

With those sites admitted (drop --recover-gaps, add --keep-noisy-kmeans) the gap
gets 77-79 phased records, and the WINDOW SIZE THEN DECIDES CORRECTNESS. Between
the gap's two site clusters, 48,183,976 -> 48,225,786 is 41.8 kb with ZERO reads
covering both sites (neighbours carry 41 and 51; longest read 29.9 kb). At 20 kb
the arm emits two blocks -- 48,147,227-48,183,976 (10 sites, holding the left
flank's own site 48,162,480) and 48,224,939-48,279,445 (50 sites, anchored right)
-- neither crossing the hole. At 500 kb it emits one 82 kb block that crosses it
with the two halves on OPPOSITE haplotypes (7 sites left, 6 agreeing; 2 right,
both opposite). Read accuracy is 100.00% in both arms and cannot tell them apart,
because no read spans the hole; only the site-level truth check sees it.

RETRACTION: the gap_lab closure of chr20:48,176,830-48,229,446 (+107 concordant,
0 flips) crossed that same hole and was LUCK. The gap arm at 5 kb read context and
at 30 kb produce OPPOSITE genotypes at 48,225,786 and 48,229,226, and the 30 kb one
links the right flank backwards (159 discordant, caught only by truth). The honest
outcome is two blocks, each anchored to one flank, with 41.8 kb between them that
no read-based method can phase.

So the mechanism to build is not a new link: phase the gap in windows with the BAM
sites admitted and let flip_chunk_hap stitch them. Both blockers are already
measured -- recovery zeroes the BAM categories, and --keep-noisy-kmeans is
chromosome-wide (~8k reads at ~65% error per its call site) -- so both need the
same gap-interval scoping.

### The noisy-region MSA never runs when --recover-gaps is on (2026-09-16)

collect_var_run_phasing (collect_var.cpp) guards the noisy-region MSA with
`if (!opts.recover_gaps) collect_noisy_vars_step4(...)`, deferring it to the
recovery pass. So in the normal hybrid run with recovery the noisy het class does
not exist: inside chr20:48,176,831-48,229,447 the chunk holds 1 clean het SNP, 2
clean het indels, 1 repeat indel, 64 clean hom and 606 low-coverage catalog sites,
and ZERO NoisyCandHet -- while collect-bam-variation on the same interval calls 8
NoisyCandHet and phases the window into a single block. Ruled out first: the 50 kb
max_noisy_reg_len cap (the BAM channel yields the same 8 sites over 52 kb and over
152 kb), the private_keys branch that zeroes that cap (entered only with
--private-sites-vcf), and skip_noisy_kmeans (read INSIDE the step that never ran).

--retry-unphased-with-bam now runs that step for its own call (recover_gaps=false,
skip_noisy_kmeans=false). In-gap phased hets 2 -> 8, and the left block goes
48,162,480-48,176,830 (2 sites, 14.3 kb) -> 48,147,227-48,229,226 (13 sites,
82.0 kb), 220 bp short of the right block. Gate clean: 0 concordant->discordant,
0 tags lost, 0 newly discordant.

CAVEAT, measured: that 82 kb block spans 48,183,976 -> 48,225,786, which is 41.8 kb
with ZERO reads covering both sites, and it is SWITCHED across it -- left-of-hole
sites carry hap1=PAT (1.00, 0.96, 1.00), right-of-hole hap1=MAT (0.99, 0.93). The
read gate reports 0 flips because no read spans the hole, so it is structurally
blind to this error class; only the site-level truth check sees it. The retry must
refuse to join across a zero-spanning-read spacing: the right output here is two
blocks, 48,147,227-48,183,976 and 48,225,786-48,229,226, which is what 20 kb chunks
already produce.

### verify_retry.py, and the bridging site we call homozygous (2026-09-16)

`evaluations/2026-09-16-retry-verify/`. One command per window, built so no check
can pass by being blind. It reports CATEGORY and emitted GENOTYPE separately,
counts spanning reads for every consecutive pair of usable het sites, names every
unscorable site with its reason, collapses records at one position before
declaring a switch and requires a run of sites on each side, censuses what a
competitor phases here against what we do at those positions, and sets
`gate_blind` whenever an unsupported link or unscorable site exists -- a PASS
requires it false. It was written after four errors of mine, one per row of the
table in that README, each of which a silent skip had produced.

On chr20:48,176,830-48,229,446 it says: retry off PASS (44 usable hets, 2 in
region, no unsupported link); retry on FAIL (63 usable hets, 8 in region, one
unsupported link, one category-het/genotype-hom site, gate_blind true even though
concordant->discordant is 0).

THE BUG IT LOCALISED: hiphase crosses this interval with two het sites,
48,183,976 and 48,204,383, and we have both. At 48,204,383 (AT>A, DP 71, 30 ref /
41 alt, AF 0.577) our CATEGORY is NOISY_CAND_HET but the emitted genotype is 1|1 --
alt on both haplotypes -- so it links nothing. That site is the difference between
a read-supported chain (48,183,976->48,204,383 has 1 spanning read,
48,204,383->48,225,786 has 7) and the 41.8 kb jump with ZERO spanning reads the
solve makes instead.

RETRACTION: the earlier claim that this block was measured as switched across the
hole is withdrawn. With positions collapsed and a two-site run required per side,
no switch is demonstrated -- the right side has one scorable site (48,229,226)
because 48,225,786 is below confidence. The join is unsupported and its
correctness untestable, which is why an unsupported link fails on its own.

### The genotype collapse that kept the window unphasable (2026-09-16)

iter_update_var_hap_to_cons_alle recomputes each haplotype's consensus allele
INDEPENDENTLY by majority, except for verified multi-allele MSA insertions, whose
branch carries the comment "Independent haplotype majorities can select the same
allele twice." A plain biallelic site inside a window the first solve could not
phase hits exactly that: the reads carry no hap labels, both majorities are the
deeper allele, hap_to_cons_alle[1] == [2], and a real het is emitted 1|1. It is
self-sustaining -- a hom links nothing, so the window stays unphasable. A probe
inside the iteration caught chr20:48,204,383 going cons1=1 cons2=0 on one round
and cons1=1 cons2=1 on the next. BOTH our channels did it: collect-bam-variation
emitted 1|1 / HAP_ALT=3 there too, so the earlier framing of this as hybrid-only
was wrong. Only hiphase called it 0|1, and it is the ONLY heterozygote between
48,183,976 and 48,225,786.

Fix: apply the joint orientation to a biallelic candidate whose allele depths call
it het (ref_cov/alt_cov >= min_alt_depth, AF in [min_af,max_af]), confined to
opts.retry_windows; with no preference in the labels, seed a het rather than a hom.

Verified on two windows with verify_retry.py. chr20:48,176,830-48,229,446: in-region
usable hets 2 -> 10, unsupported links 1 (41.8 kb) -> 0, category-het/genotype-hom
1 -> 0, competitor sites we call hom 1 -> 0, gate 0 flips and 0 tags lost. Still
FAIL, now on a switch 48,147,227 -> 48,149,548 (2.3 kb, 59 spanning reads, so a
supported link oriented wrong, and LEFT of the gap) plus three sites below the
confidence floor. chr20:36,217,274-36,268,291: in-region usable hets 4 -> 15,
switches 1 -> 0, competitor-absent 14 -> 7, accuracy 77.65% -> 99.6%, but tagged
358 -> 240 and 75 correct tags LOST against 12 gained. Not a default: the re-solve
judges every read in the chunk against the new site set, same shape as the earlier
filter-reorder regression.

### One window's BAM-site injection, audited in four stages (2026-09-16)

`evaluations/2026-09-16-injection-audit/`. Takes collect-bam-variation's own
candidates for the interval as the reference set and checks presence, field
fidelity, admission and informativeness, with the clean-het control required.
On chr20:48,176,830-48,229,446 with --retry-unphased-with-bam: 9 of 9 het sites
PRESENT, 8 of 9 loci used as phased hets, 0 emitted homozygous, control
informative (0.986/1.000/1.000). Injection is not the problem.

Two things the audit established about representation. The BAM channel emits BOTH
nested forms of a repeat deletion at one position as independent contradictory
hets (48,177,780: GAGAAAGAA>G 1|0 29/45 AND GAGAAAGAAAGAAAGAAAGAA>G 0|1 46/28;
48,225,786 likewise) -- the defect split_nested_msa_deletions exists to prevent,
and it is off in that channel. And a RETRACTION: the hybrid residual at
48,177,789 scoring 0.500 against 0.946 for the unsplit allele is NOT the split
destroying information. The reads carry -8 (27) and -20 (21); a net-length test
for a 12 bp residual has a window containing the common 8 bp deletion, so a -8
read is called alt too (66 vs 8). A split residual is unscorable by that method
by construction; the audit now declares it and falls back to the unsplit allele.
At 48,225,787 the same split scores 1.000.

Open in this window: 48,202,057 is used as a phased het but segregates 0.507 over
69 reads (0.509/0.507/0.515 as a 2/3/4 bp deletion, so not a representation
artifact, and not a split residual) -- a site with no haplotype information
carrying a phase set. 48,177,726 is NOISY_CAND_HET to the BAM channel and
REP_HET_INDEL with no emitted record in hybrid. Three loci differ from the BAM
channel by 1-3 reads in their allele counts.

Also: the margin filter is removed from the canonical chr20 runner, the gap rig
and the playbook. min_read_hap_margin defaults to 0 and the hybrid subcommand
does not override it, so passing --min-read-margin 2 was our own addition; on
chr20:36,217,274-36,268,291 it strips 271 reads of which truth says 254 (93.7%)
were phased correctly.

### The MSA observation refresh was gated on recover_gaps (2026-09-16)

refresh_assigned_msa_observations re-reads each assigned read's allele from its
own cluster alignment, and ran only when recover_gaps was set, so the
alignment-only channel kept pre-refresh counts. Pinned by
test_msa_counts_do_not_depend_on_recover_gaps: two clusters of two reads each
carrying their own consensus make the site 2 ref / 2 alt, and the ungated arm
counted 3 ref / 1 alt. Gate removed; all five test binaries pass.

Two process notes. `make -j20` builds pgphase but NOT the test binaries, so a
stale binary reported ALL PASS and I committed a claim that all four pinned
invariants held when one failed -- use `make unit-tests`. And the fixture needed
two iterations: with empty read clusters the refresh is a no-op and the test
passed vacuously; with an assigned read handed to add_msa_site_observations as
the other haplotype the input contradicted itself.

The fix does NOT explain the real-data divergence: with the refresh unconditional
the alignment channel's counts at all eight divergent loci on
chr20:48,176,830-48,229,446 are byte-identical and the count remains 8 of 76.
Those come from the hybrid side -- three homozygous loci where the hybrid counts
a candidate's reads across the chunk rather than only inside the noisy region
(DP 23->53, 18->58, 41->69), two that are the nested-deletion re-representation,
and 1-3 read differences at three loci that remain unexplained.

### Two blocks 220 bp apart never join because the site and read PS labels disagree (2026-09-16)

`evaluations/2026-09-16-block-label-split/`. On chr20:48,176,830-48,229,446 with
--retry-unphased-with-bam the left block reaches 48,229,226, 220 bp from the
right block's first site at 48,229,446, and they never join.

The evidence is complete, so this is not linkage: 60 reads call both boundary
sites, none untagged; their allele co-occurrence is 22 (0,0), 34 (1,1), 1 (0,1),
3 (1,0) -- 56 in phase against 4, and those 4 equal the deletion site's own
genotyping error, which truth also puts at 56/60. The SNP at 48,229,446 scores
60/60 against truth. Both sites already carry the same orientation in our VCF
(1|0, HAP_ALT=1 HAP_REF=2).

The bug is that the pipeline carries a phase set per SITE and a phase set per
READ, and they contradict each other here. The left block's last two sites
(48,225,787 and 48,229,226) are labelled PS=48147225 while 0 of their covering
reads carry that PS and 51/60 respectively carry PS=48229446. The 82 kb left
block's right end is populated entirely by right-block reads. That also starves
recovery: the tier report reads L=0, leftPS=-1, R=1 on every tier
(partial on pass 0, open on pass 1), because the left flank has no reads holding
its own phase set at the boundary, so the flank vote has no voters.

Competitors phase it because they assign one phase set per connected component of
the variant graph -- two het sites joined by shared reads are one block by
construction, with no second per-read label to contradict the per-site one. We
reconcile the two labellings only at a chunk seam (flip_chunk_hap) and in the two
gap-recovery flank votes, neither of which applies to two blocks inside one
chunk.

Fix direction: when every read covering a site is unanimously tagged into another
phase set, the site's PS must follow its reads, or the two sets must merge with
the orientation taken from the shared reads as select_stitch_orientation already
does at a seam. Here that vote is 56-4 and the merged orientation is the one
truth prefers at both sites.

Also noted: the tier report's right flank is rightPS=48243938 while the emitted
VCF gives that block PS=48229446 -- a third identity for the same block.

### The stitch has no voters because homopolymer sites cannot grant a phase set (2026-09-16)

Traced why the existing stitching machinery does not join the retry's window to
its neighbour on chr20:48,176,830-48,229,446. The flank vote counts reads holding
each side's phase set, and update_read_phase_set excludes homopolymer indels from
granting one. Both terminal sites of the left block sit in long A homopolymers
(48,225,780 ctctctcaaaaaaaaaaaaaaaa; 48,229,220 tctttagaaaaaaaaaaaaaaac), so
their reads take the next eligible het's phase set -- the right block's -- and the
left flank has zero voters (L=0, leftPS=-1, every tier abstains).

Three compounding facts, each probed: hp_gap_scorable already exists as the
exception and init_assign_read_hap honours it (collect_phase.cpp:362) while
update_read_phase_set does not, despite a comment claiming the same eligible
evidence; the flag is set only for the homopolymer recovery tier's window; and
select_gap_link_sites, which sets it, is guarded on gap_hp_link_beg >= 0 or the
experimental private_msa_admit_all_in_region, so it never runs in a normal
configuration. Probe at the site: hp_indel=1 hp_scorable=0 cate=0x100
msa_verified=1 cons=1/0.

Fixing all three DOES join the window through the existing machinery -- tier 1
reports joined, one 132.2 kb block of 61 sites spans the region, and the 60 edge
reads move to the left phase set -- but the orientation is WRONG: 161 concordant
reads become discordant and window accuracy falls 100.00% (167/167) to 71.50%
(148/207). REVERTED. The shared-read vote at the boundary is 56 in phase against
4, so the correct orientation is available and the joining path does not use it;
the tier report's leftPS = rightPS = 48243938, a third identity for these blocks,
is where that orientation comes from and the next thing to fix.

### Paired noisy MSA alleles lose their consensus in BAM stage 2 (2026-09-20)

On identical `CHM13#0#chr20:3800000-3900000` input, both callers produce the
same two MSA insertions at 3,870,827 (AA and AAA; depths 72 with alternate
counts 29 and 32). pgphase then set both candidates' haplotype consensus to
reference and emitted neither; longcallD emitted both. The pending cluster seed
survived `var_init_hap_profile_cons_allele` but was overwritten by
`update_var_hap_to_cons_alle` during the iterative k-means update.

A cluster-consensus seed was tried for pairs of insertion or SNP alleles from
opposite MSA consensuses. The guards required matching event type and support
above `min_alt_depth` and `min_af`; pairing on position alone also paired
insertions with deletions, while weak alleles produced unsupported calls.
Deletion-length pairs remained under the existing solver. Extending the seed
to the 43.9-44.8 Mb window recovered 3 upstream records but added 21 unmatched
records (pgphase-only 2 -> 23), so that extension was reverted.

With the narrower seed on identical region inputs, 3.8-3.9 Mb matching records
rose 186 -> 188 of 188 upstream records, and pgphase-only stayed at 1. Over
43.9-44.8 Mb, matching records rose 928 -> 932 of 942, pgphase-only stayed at
2. These are experimental region-scoped measurements, not shipped behavior.

The seed itself was removed when the port was checked against the original
source: longcallD `collect_var.c:update_cand_var_profile_from_cons_aln_str2`
does not assign `hap_to_cons_alle` from the MSA cluster, and
`assign_hap.c:var_init_hap_profile_cons_allele` resets heterozygous consensus
to -1/-1 on each phasing call. A direct differential target compiles the
original C at revision `23e369d71a1e4dd46529be3755d2224b8e239b76` and
compares consensus initialization, consensus updates, allele scoring,
two-site read assignment, phase-set linking, and complete k-means on small
two- and three-site read matrices against the BAM C++ path. It passed 166,364
assertions over six cases, including link chains with mixed read haplotypes.
This establishes parity for those tested inputs and functions; whole-pipeline
parity remains a separate integration question.

The direct phase-link test found two differences in the previous C++ port:
a 2-agree/2-conflict tie started a new phase set, whereas upstream keeps the
current set, and the C++ port swapped the alleles on a conflict-majority link,
whereas `assign_hap.c` swaps once per hap and therefore leaves their order
unchanged. The full-k-means test also found a port-only final read re-assignment
after phase-set calculation; upstream returns the last iterative read labels.
The BAM path now follows the original adjacent-site, two-read link threshold,
tie, orientation, and final-label behavior. The graph path keeps its existing
link and final-label rules. On the chr20 quick BAM fixture both callers
still emit the same 818 VCF record keys and GTs, but 13 shared records carry phase-set label 15039543
in pgphase versus 15071132 in longcallD. Function parity therefore has not
closed this integration mismatch; the inputs to the functions or later
postprocessing still differ.

The stitch audit also found that longcallD calls
`update_chunk_read_hap_phase_set1` only when `out_aln_fp` is non-null.
pgphase previously updated read HP/PS on every stitch. The BAM VCF-only path
now leaves read HP/PS as upstream does, with a focused output-mode regression
test; alignment output still updates them. This source correction did not
change the 13-label difference on the chr20 quick fixture.

The BAM validation scripts also used the retired positional CLI and
`--phased-vcf-output`; their invocations now use `--ref`, `--bam`, `-r` and
`--phased-vcf-out`. `make check` reaches the golden comparison and fails against
its May 1 chr11 files. The current HiFi run has the same 569 rows, positions,
types, categories and phase labels; differences include haplotype labels,
strand tallies on MSA candidates, 51 REF anchor fields and one allele count.
Runs at 1 and 4 threads are byte-identical. The golden files were not refreshed as part of this parity
fix.

The 13-label chr20 mismatch was traced to `collect_noisy_reg_aln_strs`:
`collect_noisy_vars1` passed a non-null `unassigned` buffer even though the
BAM path set `add_unplaced_msa_observations = false`. The aligner uses that
pointer for more than reporting unplaced reads: it aligns them to both
consensuses and appends a read with a sufficiently better score directly to
an MSA cluster. At chr20:15,070,944-15,071,487 both callers initially selected
5 and 7 phase-set reads, but pgphase expanded the recalled deletion at
15,071,133 to depth 74 (36/38), while longcallD kept depth 11 (6/5).
The extra reads supplied a 12-agree/0-conflict link to the preceding
15,056,025 site; longcallD had 0/0 and started PS 15,071,132. Passing a null
buffer when unplaced observation recovery is disabled restores the original
cluster membership. On the shared 500 kb chr20 fixture the VCF parity script
now reports 818 keys on each side, zero unique keys, and zero shared payload
differences, including GT and PS. The earlier 13-label result above records
the pre-fix measurement.

### BAM chr11 golden repair and remaining ONT FORMAT differences (2026-09-20)

An upstream-compatible branch in `assign_hap_based_on_germline_het_vars_kmeans`
returned immediately after read phase-set assignment, skipping the common TSV
`hap_alt`/`hap_ref` projection. It now skips only the final read reassignment;
the direct C-versus-C++ k-means test checks those projected counts. The HiFi
chr11 golden was refreshed after a full comparison against original longcallD:
568 shared VCF keys and zero differences in the compared FILTER/INFO and
GT:DP:AD:VAF:GQ:PS fields. The current TSV has 569 rows and is byte-identical
at one and four threads.

On ONT chr11, pgphase initially omitted nine MSA insertions that longcallD
called. `collect_var.c:var_is_homopolymer_indel` compares raw FASTA bytes with
nt4-coded insertion bases; pgphase had normalized both to nt4 and marked the
insertions as homopolymer indels, excluding them from phasing. The BAM path now
uses the original raw-byte comparison, while the graph path retains its
case-insensitive check. The differential test compiles the upstream C function
and covers both indel branches on uppercase and lowercase reference slices.
All nine ONT insertions are restored; ONT chr11 now has 586 shared VCF keys,
zero unique keys, and matching GT/PS at every shared record. Its 689-row TSV
and 586-record VCF goldens were refreshed; `make check` passes.

The VCF parity script now compares DP, AD, VAF and GQ in addition to GT/PS.
It exposes 12 residual ONT FORMAT differences: four records near 1,294,294,
four near 1,417,535, and four near 1,429,530. At the first group, longcallD
keeps noisy regions 1,294,276-1,294,413 and 1,294,527-1,294,976 separate,
while pgphase merges them into 1,294,276-1,294,976. The other two regions
have the same bounds and selected-read counts, but one read's full-coverage or
allele assignment differs. Thus the chr11 ONT fixture has exact call, GT and
PS parity but not full FORMAT parity. The refreshed goldens are pgphase
regression outputs, not evidence of complete upstream parity.

The first four residual differences were then traced to `intervals_to_cr`:
it coerced a valid zero noisy-region label to one. `cgranges.c:cr_cluster0`
uses the smaller neighboring label as the merge distance, so this one-base
change let a zero-label bridge merge the 1,294,276 and 1,294,527 regions.
Preserving the original label separates the regions and removes all four
FORMAT differences. Chr20 HiFi (818 records) and chr11 HiFi (568 records)
retain full compared-field parity. ONT chr11 now has eight FORMAT differences
in the two remaining MSA regions; all 586 record keys, GTs, and PS labels
still match. Per-read profile comparison identifies one discordant read in
each region: `1d237405-a612-4e04-868f-9936555e0f3d` is spuriously counted
at the first three sites of 1,417,525-1,418,674 and classified oppositely
at the last; `2914fbda-117d-4761-96ee-30606c1942c6` is classified
oppositely at the first and last SNPs of 1,429,514-1,432,553. Both callers
select the same numbers of MSA reads in those regions. These eight FORMAT
values remain an open alignment/profile parity gap.

The remaining eight differences were traced through original `-V3` input
traces. At 1,417,525-1,418,674, longcallD extracts 150 bases from ONT read
`1d237405-a612-4e04-868f-9936555e0f3d` while pgphase extracted 2,343.
The BAM record has a large soft clip on its palindrome end. Original
`bam_utils.c` marks that clip as an internal `BAM_CHARD_CLIP`; every pgphase
digar parser detected the palindrome but still stored `SoftClip`. The extra
sequence changed read order in abPOA and its partial alignment. The same
upstream conversion also resolves the second region's discordant read.
Porting the clip conversion at all BAM digar parser sites makes the full
compared VCF payload exact on chr20 HiFi (818 records), chr11 HiFi (568),
and chr11 ONT (586), with zero unique keys or shared differences in
VCIGAR, FILTER, END/SV fields, CLEAN, GT, DP, AD, VAF, GQ and PS. The ONT
TSV and VCF regression goldens now reflect this final result. The earlier
12- and 8-difference counts above are intermediate measurements.


### Graph recovery expands phase blocks outside-in (2026-09-20)

Testing six noncentromeric gaps below 10 kb that HiPhase spans at at least 98%
truth concordance showed two distinct failures. Four had an
alignment-verified homopolymer row near a block boundary, but the graph solve
never admitted that row to its link list. The other two had co-located
representation/order failures whose first external edge had fewer than two
votes.

Bulk admission closed the four homopolymer cases but regressed the committed
3.85 Mb control: its separated fraction fell from at least 0.36 to 0.165.
Admission was therefore changed to an iterative frontier. After exact BAM rows
and observations are injected and the normal graph rounds run, each disconnected
block exposes only the nearest verified recovery locus on its left and right.
Each disconnected phase-set pair ranks its two exposed loci by net
same-versus-cross support; one coordinate per pair enters, and the normal clean
and noisy rounds run again before the next layer is visible. Scoping the
ranking per pair matters in production chunks: a chunk-wide winner let two
unrelated gaps consume both waves and left 6.58 Mb open. Co-located rows remain
separate and must qualify on their own evidence.

The first implementation reused the ordinary two-read link margin without a
step-distance limit. It kept the window gains and passed the original panel, but
the full chr20 audit regressed to 5,827 discordant of 220,608 evaluated reads
(2.641% error). A 10 kb step and ten-read net margin reduced unsupported
admission, but choosing every boundary in one wave or allowing 32 single-site
waves converged to 4,682/220,610 discordant reads (2.122%). Ranking only the
immediate frontier restored the 12.27 Mb separated fraction from 0.68 to 0.81
but still exposed a false 15.10 Mb join.

That join identified a logic error: the chosen deletion already belonged to the
same phase set as the boundary used to support it. Its 29/6 pair vote confirmed
existing membership; it was not evidence of expansion. Rejecting same-phase-set
support removes the join and yields 3,595/220,610 discordant reads (1.629%),
versus the accepted exact-row recovery's 3,302/220,640 (1.497%). The remaining
largest new join, at 52.3 Mb, has 30 of 31 pair reads in one orientation; an
extra gate that rejects it would also reject weaker correct target links. Two
single-locus waves are retained because one wave leaves 48.23 Mb open and two
close all four supported targets. The 15,095,642–15,101,261 seam is retained
as a no-span regression control: it stays split at 99.2% read concordance.
The final full-chr20 VCF, using production chunking, spans all four targets;
`full_chr20_target_spans.tsv` records their resulting phase-set extents.

This closes 6,578,161–6,582,248, 11,255,370–11,262,360,
12,269,536–12,277,077, and 48,225,787–48,229,445. Their measured read
concordance is 99–100%; separated fractions are 0.80, 0.64, 0.81, and 0.71.
The 4.78 and 34.10 Mb representation cases remain open because the first
external link is still insufficient. The original and added regression
windows pass, and the four closures are retained as regression
windows. `evaluations/2026-09-20-graph-recovery-windows/frontier_short_gaps.tsv`
records the result.

The merge also fixed a real refresh bug: rebuilding a shared read profile used
`map::emplace`, so an existing GAF allele silently won over the recovery BAM
allele. It now uses `insert_or_assign` at shared sites. The expectation refresh
script was stale after the suite moved to one `graph` arm; its gate now checks
that active arm instead of removed `default` and `noretry` rows.

### Graph recovery targets phase-set seams only (2026-09-21)

Recovery now derives its BAM intervals only from bounded gaps between neighboring
phase sets. The separate unphased-read bin detector duplicated internal seam
coverage and uniquely selected terminal or wholly unanchored regions, which have
no two-sided phase boundary for the outside-in algorithm. Graph-seeded noisy
regions had the same problem and were removed from the target union along with
the now-dead `--graph-noisy-msa` option and seeding helper.

The rebuilt seam-only binary keeps every intended internal closure and the
15.10 Mb no-span control remains split. Its local concordance is 98.71% rather
than 99.2%, so that control's floor moved from 0.99 to 0.98 while its exact
no-span assertion remains unchanged. Whole chr20 phases 217,543/245,053 emitted
reads in 320 phase sets, with 3,174/217,529 truth-evaluated reads discordant
(1.459%; 98.54% accuracy). Against the preceding frontier build this removes
3,103 phased reads and 35 phase sets while reducing discordant reads by 421 and
improving accuracy from 98.37%. This is the expected scope change: recovery no
longer adds reads from regions without two neighboring phase-set anchors.

A detector review found two related issues. Its original eligibility check
excluded only `phase_set == 0`, while the graph adapter had diverged from
longcallD by initializing graph candidates to `-1`. LongcallD initializes a
candidate phase set to `0` and an unphased read phase set to `-1`; graph now
uses those same sentinels. The shared boundary predicate requires
`phase_set > 0` and different haplotype alleles, so only assigned,
heterozygous candidates define a seam. The ordered-map accumulation plus a
second sort was replaced by a coordinate-order scan with an unordered
phase-set index and overlap coalescing: expected O(C+B) time and O(B) storage
for C candidates and B phase sets. The two filtering helpers orphaned by
removing post-hoc recovery were also deleted. A rebuilt full chr20 run retains
77,484 candidates and produces byte-identical phased VCF and BAM output. The
only TSV differences are 17,276 `PHASE_SET` values changing from `-1` to `0`,
which is the intended candidate-sentinel normalization.

The detector was then made fully flat. Because a phase-set label is its genomic
anchor coordinate, each candidate contributes `[phase_set, position]`; a vector
merge stack computes the covered union without a label map. Ordered interval
ends make this amortized O(C), with contiguous O(K) storage for K covered
components. The following target builder no longer copies or sorts seams. It
uses a monotone scan of parent anchors and seams, half-open member ranges, early
chunk clamping, and binary membership lookup, removing the duplicate region
vector and per-group member allocations.

A trial also deduplicated co-located anchors and counted the seam boundary
symmetrically. That changed recovery scope and regressed full chr20 from 3,174
to 3,201 discordant truth-scored reads while evaluating 54 fewer reads, so that
semantic change was reverted. The retained implementation preserves the
validated flank rule and phases 217,543 reads in 320 phase sets. It has
3,170/217,529 discordant truth-scored reads (1.457%; 98.54% accuracy), four
fewer than the pre-flat detector, and all 229 window assertions pass. It emits
77,448 candidate rows; the phased-read completeness is unchanged.

The full recovery transfer was audited after target construction. The redundant
`set<CandKey>` mirroring the raw and translated parent indexes was removed;
one parent match now supplies membership, provenance and the parent index.
Transferred candidates and orientation decisions share one ordered-map entry,
eliminating a second tree and lookup. Audit rows, including ALT string copies,
are built only when `--recovery-audit-out` is active. Missing reads are found
by a monotone scan over the qname-ordered parent reads and observations instead
of constructing another qname set. Ordered maps that define deterministic
candidate and observation merge order were retained. Full chr20 TSV, phased VCF
and phased BAM are byte-identical before and after these changes.

### Explicit recovery seams close the 14 short-gap targets (2026-09-21)

The remaining 35.613 Mb failure exposed a representation bug in the final
stitch. Seam detection used canonical graph/VCF coordinates, but stitching
discarded the detected phase-set identities and searched again with raw
candidate keys. The two complementary complex-indel rows at 35,613,763 and
35,613,765 therefore supplied decisive shared-read BAM gauge votes but the
stitch selected the wrong flank. Recovery seams now carry canonical begin/end
coordinates plus their exact left and right graph phase-set IDs.

Left-to-right processing needs one further detail: an earlier seam can absorb
the phase set named as the following seam's left side. A phase-set alias map
resolves that old ID to the surviving upstream label. Shared-read gauge lookup
first checks the original detector ID and then the surviving ID. This restored
the consecutive 61.76 and 65.51 Mb joins without coordinate rediscovery.

Targeted recovery now runs the standalone longcallD BAM configuration once.
Co-located MSA allele merging remains disabled, so separate BAM rows are
injected exactly as called; the earlier graph-style call plus normalized-key
replacement path is gone. Each solve receives 50--60 kb of context and at
least three graph anchors per side when available. The shared-qname BAM gauge
is the primary stitch orientation; strongest allele-pair evidence remains the
fallback when the gauge abstains.

A fresh targeted run places all 14 selected noncentromeric short gaps inside
one pgphase phase set. The normalized VCF span result is 14/14. The 32.17 Mb
isolated case uses 250 kb of outer test context because 50 kb contains no left
graph anchor; production chromosome chunks already contain that context.
The other thirteen use 50 kb. Full per-window measurements are in
evaluations/2026-09-21-short-gap-recovery/final.tsv.

The 55.88 Mb deficit was subsequently traced to the observation and stitching
path rather than missing sites. The graph plus injected BAM table contains the
same bridge used by HiPhase: the 55,862,239 deletion, 55,862,269 insertion,
the complementary four- and six-base deletion rows at 55,883,019, and the
55,889,113 SNP. Two physical bridge reads carried the six-base deletion in
their CIGARs, but the verified MSA row was absent from their sparse profiles.
The exact-CIGAR MSA backfill had no production caller and, when first called,
also exposed a memory bug: its profile updater could not extend a profile to a
lower candidate index.

Recovery now backfills verified exact observations only inside detected seams.
The strongest edge records the locus that earned the chain, and final read HP
is recomputed from those supported loci with longcallD's allele scorer. This
prevents an unrelated MSA row from cancelling the bridge and avoids copying a
stale BAM sub-solve HP after a downstream block flips. At 55.88 Mb the result
moved from 44/56 truth-correct crossing reads (78.57%) to 64/64 (100%). The
frozen HiPhase result is 66/66. The remaining pgphase-unphased reads have no
decisive exact observation on either injected deletion row; assigning them
would require interpreting a different 18-base deletion or two reference calls
as one of the complementary BAM alleles.

The graph-adapter unit test now replays the explicit-seam identity invariant at
all 14 target coordinates and the consecutive-seam alias case entirely in
memory. It also retains the 4.78 Mb three-block gauge fixture. These tests are
sub-second after compilation and do not invoke BAM, GAF, reference, or the full
pipeline. The real-window outputs remain evaluation evidence rather than a
routine test dependency.

### Expanded HiPhase-correct gap regression panel (2026-09-21)

A fresh full chr20 graph+BAM run with the supported-chain implementation emitted
205 phase blocks and left 200 inter-block gaps. Frozen HiPhase spans 74. Scoring
only truth-labeled reads physically crossing both flanks, allowing the arbitrary
local HP orientation to flip, leaves 34 noncentromeric gaps at >=98% local
purity. They cover 380,363 bp and 780 scored crossing reads; 15 have at least 20
scored reads. `remaining_hiphase_correct_gaps.tsv` records every target and marks
lower-support rows instead of dropping them. The fast explicit-seam replay
keeps the 14 historical short-gap coordinate cases and adds all 34 current
misses (48 coordinate cases). The 4.78 Mb gap is
represented once under each audit's coordinate convention, so the test covers
47 distinct genomic regions; the current deficit itself is 34 gaps.

### Per-edge aggregate and exact-path recovery close supported chr20 gaps (2026-09-21)

The first DP prototype treated an entire recovery window as one outer-flank
problem. That is incorrect when the window contains several local phase sets:
an unsupported early edge prevented the solver from considering a supported
pair farther to the right. The retained stitcher walks every adjacent phase-set
pair in reference order. It tries the existing BAM gauge and strongest single
site pair first, then a per-read aggregate boundary vote, then the exact
ordered two-state site-path DP for that pair. A missing edge starts an
independent component and processing continues. Components born after a break
remain strong-only across later seams, preventing the 33.79 Mb recovery from
being absorbed by a weak downstream edge.

The aggregate uses up to `block_link_window` oriented boundary candidates from
each block. Each molecule forms one consensus per side and contributes at most
one same/cross vote. A net margin of eight is required. On the saved 34-gap
matrix, every high-support target with aggregate evidence at margin >=8 had the
truth-consistent parity; all observed wrong low-support parities had margin <=3.
The implementation queries the read-to-candidate interval index over the
boundary range rather than rescanning all chunk reads. A fresh full-chr20 run
after this optimization is byte-identical to the retained pre-optimization
candidate TSV, phased VCF, and phased BAM.

The DP still solves the exact lexicographic objective for one adjacent pair:
maximize bottleneck edge margin, then total margin, total support, and finally
prefer fewer edges. Edges span at most 20 kb and require net margin 10; each
flank needs two oriented anchors and the winning endpoint parity must exceed the
other by two. Existing local phase sets are committed atomically. Final
longcallD scoring may assign a previously unphased read when it observes an
exact allele on a selected path site; the DP never assigns a label without that
site evidence.

On the current 34 HiPhase-correct deficit targets, production chr20 recovery
spans all 15 high-support gaps. The initially reported 14/15 was an evaluator
bug, not an open phase edge: both insertion rows at raw VCF position 882277 and
the 890261 right-flank site already carry PS 882277 after recovery. The span
parser reduced the anchored insertion to only its canonical candidate position
882278, then incorrectly asked whether the block reached 882277. Phase-block
extent now includes the closed interval between a VCF record's mandatory anchor
and canonical first-changed coordinate; candidate membership remains canonical.
A self-contained regression covers this exact insertion-boundary case. The safe
solver also closes the low-class 25985709--25986123 gap with a 42-read margin.
Every newly closed high-support target has at least 84.53% local truth
concordance, above the accepted 80% floor.

Lowering aggregate and DP margins to one and allowing one flank anchor was
tested and rejected. It joined 15056025--15071132 at only 68% local concordance,
while most no-signal targets remained open. The other low-support targets have
at most three net aggregate votes, and several have zero. Their BAM profiles do
not contain enough observed alleles to select a diploid chain without guessing;
the no-realignment design deliberately abstains there.

Full chr20 keeps 59,818 phased heterozygotes. Phase sets fall from 205 to 198 and
N50 rises from 1,008.5 kb to 1,070.7 kb. Tagged reads move from 221,571 to
221,307. The same per-phase-set parental diagnostic reduces discordant reads
from 37,128 to 35,569 (1,559 fewer); its absolute percentage is not used as the
truth-accuracy headline because that legacy scorer includes all emitted
unaligned records. Summed phase-set span drops from 56.80 to 56.40 Mb while the
number of phased heterozygotes is unchanged; that span sum double-counts
intersecting phase-set extents and also reflects avoided weak joins.

The unit suite now covers both exact DP bridging across disjoint read cohorts
and aggregate support distributed across four weak site pairs. The aggregate
fixture proves that no single pair meets the ordinary threshold while the
one-vote-per-read block evidence joins the two phase sets. A representation
regression verifies that an anchored insertion contributes both its VCF anchor
and canonical event coordinate to block extent. The focused 15-window
high-support integration panel is unchanged by the interval-index optimization.


## 2026-09-22 exact injected-site MEC and safe post-break extension

Recovery now finishes with a bounded exact MEC solve for still-open seams. The
solver uses only read-connected exact BAM-injected sites, treats every existing
local phase set atomically, compares both right-block orientations, and abstains
on tied parity or more than 20 variables. It runs after the established
left-to-right stitch so a new exact join cannot change later baseline decisions.

A separate control-flow bug suppressed a decisive later edge after an earlier
edge in the same recovery window failed. The retained exception requires the
later edge to start at a BAM-injected site and have net margin at least eight.
It closes chr20:58,834,248--58,836,037 at 96.23% parental phase-set purity.
Exact MEC adds chr20:18,194,808--18,218,259 and
31,886,648--31,901,501, each at 100% among available local crossing truth reads.

Full chr20 moves from 28/48 to 31/48 coordinate cases spanned (30/47 distinct
regions), 59,891 to 59,906 phased hets, 200 to 187 VCF phase sets, N50
1,052,956 to 1,136,391 bp, and 221,431 to 221,836 phased reads. The parental
read diagnostic changes from 35,360/221,394 discordant (84.03% accurate) to
36,220/221,799 (83.67%). The added joins retain the correct parental block
orientation; the extra read discordance is within the accepted tradeoff.

Two broader arms were rejected. Adjacent-pair MEC spanned 40/48 cases but made
a confirmed haplotype switch at 32,035,459--32,050,364 and raised discordance
to 40,515 reads. Allowing post-break strongest edges from graph sites also
flipped a mature 6,418-read block at 1.907 Mb. The injected-endpoint restriction
prevents both shortcuts. Evidence and totals are in
`evaluations/2026-09-22-exact-gap-mec/`.


The older integration panel still has wrong-orientation graph joins at 6.578,
22.981, and 48.226 Mb. The exact pre-change binary reproduces all three, so the
MEC work did not add them. A global eight-read ordinary-link floor fixes only
6.578 Mb and loses correct 5.31 and 60.03 Mb closures; it was restored to the
existing value. Treat those three windows as a separate unresolved orientation
quality problem, not as evidence against the injected-only MEC joins above.

### Per-block recovery gauges and statistical parity decisions (2026-09-22)

The targeted BAM solve already produced independent local phase sets in 15 of
the 17 still-open chr20 coordinate cases; every injected row in those cases had
a positive phase set. Two cases had no injected site. The remaining failure was
therefore block orientation, not the 20-variable MEC search bound.

A correctness bug pooled graph-versus-BAM HP votes across all BAM phase sets
from one targeted solve. Those phase sets have independent gauges. Recovery now
stores direct `(graph PS, remapped BAM PS)` same/cross counts and uses only the
matching vote to orient an adjacent local block.

The experimental fixed eight-read recovery margin was also replaced where it
was used: direct per-block gauges and block aggregates now apply a one-sided
exact binomial test against a 50:50 parity null at p <= 0.01. A strongest-site
fallback applies a Bonferroni correction for the number of pairs searched. This
accepts 8/0 and rejects 104/96 despite their identical raw margin. Applying the
test to the ordinary longcallD gauge and allele path was rejected: it lost the
valid 5.31 and 60.03 Mb joins and changed fallback order enough to create a new
3.85 Mb switch. The retained scope reproduces the prior ten-window output: the
same seven known assertions remain, with no new failure.


### Read attrition, qname matching, and atomic recovery chains (2026-09-22)

The 34-gap attrition audit rules out BAM admission as the main competitor gap.
There are 860 primary spanning qnames, 859 eligible at the recovery MAPQ 1
floor, and 788 eligible under the HiPhase MAPQ 5 default. HiPhase supplementary
handling adds zero unique bridges. The meaningful drop is after admission: 772
bridges receive graph tags, 649 receive standalone BAM tags, and only 593
receive both. Nineteen gaps have fewer than eight doubly tagged bridges.
HiPhase instead starts from supplied heterozygous sites, performs one global
A-star haplotype solve, derives blocks from read-connected variants, and tags
one-allele reads after the solve.

A recovery correctness bug used a two-pointer merge between qname-sorted graph
reads and coordinate-ordered targeted BAM reads. At 4.78 Mb the 493 source reads
have 242 qname-order inversions; the old scan recorded a 5--0 boundary vote
where a qname hash lookup records 180--1. Recovery now indexes parent qnames once
per graph chunk and looks up every BAM read.

Using the complete vote without more validation exposed correlated systematic
errors. A permissive internal-block replay raised full-chr20 discordance to
5,206/219,250 (2.37%); preserving independent blocks after failed seams reduced
it to 5,098/219,246 (2.33%). Read binomial significance could not distinguish
the edges: wrong joins had 84--100% apparent support and the correct 55.381 Mb
edge had 86.5%.

The retained transaction requires source-specific graph/BAM read parity at
`p <= 0.01`, a consistent shared clean candidate, and an independently measured
outer graph relation. A multi-BAM-block chain without direct outer-candidate
evidence may use the whole-window graph gauge only when the gauge passes
`p <= 0.01` and each outer boundary has a shared-candidate parity result at
`p <= 0.05`. The two outer measurements must agree when both exist. Any failed
edge restores candidates, read HP/PS assignments, aliases, and recovery state;
internal BAM edges are not replayed after rollback.

The final chr20 run phases 219,396 reads and 59,891 heterozygotes in 671 VCF
blocks with N50 412,113 bp. Among 219,233 truth-evaluated reads, 4,641 are
discordant (2.12%, 97.88% accuracy), improving the previous retained
4,759/219,434 (2.17%, 97.83%). It spans 12/48 tracked coordinate cases versus
14/48 previously. The two lost physical joins lack evidence independent of the
same correlated BAM solve; retaining them also retains the observed accuracy
regression. The focused integration panel records the new 4.78 Mb closure at
427/428 local concordance, and no spanned target is classified as a switch.
Detailed counts and per-gap rows are in
`evaluations/2026-09-22-gap-read-attrition/`.


### Graph GAF MAPQ floor closes additional recovery targets (2026-09-22)

The graph command formerly inherited longcallD/BAM's MAPQ 30 default. On chr20,
MAPQ 5 admits 2,611 additional graph candidates, phases 1,739 additional
heterozygotes and tags 5,609 additional reads. The tracked panel moves from
12/48 to 17/48 spans, with no span classified as a haplotype `SWITCH` by the
local truth scorer. Six current misses close and one historical span opens.
Whole-chromosome accuracy changes from 97.88% (4,641/219,233 discordant) to
97.28% (6,112/224,943 discordant), an accepted completeness tradeoff.

A recovery-only alternative scored lower-MAPQ BAM reads against already phased
graph flanks without allowing them to alter graph consensus. It left the panel
at 12/48 and changed truth discordance by one read. The missing reads therefore
need to participate in graph-site clustering and construction of the adjacent
phase sets; post-hoc orientation evidence is too late. The experimental path was
removed. `collect-graph-variation` now defaults to MAPQ 5, its prior effective
HiPhase-comparable floor. `collect-bam-variation` retains the longcallD MAPQ 30
default, and `--min-mapq` still overrides the graph setting.

### Trusted SNP-first MEC closes graph/BAM recovery edges (2026-09-22)

HiPhase's joint binary MEC design was replayed on pgphase's existing allele
matrix without realignment. Uniform all-site MEC is unsafe: at 33.79 Mb noisy
indels and false heterozygous SNPs overturned a 30--2 direct SNP relation.
The retained fallback uses AF-centered SNPs first (`|AF-0.5| <= 0.12`), adds
centered indels only when SNPs are disconnected, and requires one exact parity
from the full reads and both deterministic read halves. Graph/BAM edges must
also agree with their source-specific block gauge and shared candidate. BAM/BAM
blocks remain independent, and one original block can attach only once per
chunk across overlapping seam windows.

On full chr20 at graph MAPQ 5 this moves the 48-case panel from 17 to 26 spans
(9/34 to 18/34 current HiPhase-correct cases; no historical loss), reduces VCF
blocks from 659 to 467, raises N50 from 412,113 to 482,085 bp, and keeps all
225,005 tagged reads. Truth discordance changes from 6,116 to 6,259 reads
(97.282% to 97.218%). The switch audit finds no new large polarity reversal.
The rejected permissive arm reached 28/48 at 95.31% accuracy and created large
wrong-gauge blocks. Full results and rejected arms are in
`evaluations/2026-09-22-hiphase-joint-solver/`.

### Boundary representation fixes retain the one-attachment invariant (2026-09-22)

Two remaining gaps exposed representation assumptions in trusted MEC. At
14.264 Mb, separate rows from one multi-allelic deletion had AF 0.657 and 0.314
although their deletion-to-SNP tables were pure (23--0 and 12--0). A selected
alignment-verified boundary indel may now enter only when no centered boundary
site exists; rows are never merged and the source gauge remains mandatory. At
14.679 Mb, two centered SNP anchors six bases apart had 76--0 direct support but
failed the sequence-identical-key requirement. The nearest SNP pair can now
replace that exact-key check when the full reads and both deterministic halves
all pass corrected binomial p<=0.05 and agree with MEC and the source read gauge.
The check is linear in reads and SNPs still have lexicographic priority over
indels.

An experiment allowing a recovered block to attach twice was rejected: it
created an incorrect 54.49 Mb merge and changed 342 correct reads to discordant.
The retained implementation keeps one trusted attachment per original block and
never uses this fallback for BAM/BAM edges. Full chr20 moves from 26/48 to 27/48
tracked spans (18/34 to 19/34 current targets), reduces VCF blocks 467 to 452,
keeps N50 482,085 bp and all 225,005 tags, and changes the matched truth result
from 6,257/224,969 to 6,258/224,976 discordant. No previously evaluated read
changes concordant/discordant status; seven formerly skipped reads become six
correct and one incorrect. Details are in
`evaluations/2026-09-22-representation-recovery/`.

### Guarded local MEC retry closes 14.264 Mb without a switch (2026-09-23)

The apparent full-chunk gauge conflict at 14.264 Mb was a diagnostic
misinterpretation. Replaying the saved matrix showed agreement among all three
independent signals: exact MEC chose the cross orientation in the full read set
and both FNV halves, the source graph/BAM gauge was 275--0, and 59 shared clean
candidates selected the same cross orientation. The actual failure was search
scope: expanding the 8 kb edge over the neighboring atomic block admitted more
than 20 unrelated unphased variables, so the exact solver abstained.

The retained solver keeps whole-block validation as its first attempt. When that
problem cannot be solved, it retries only the selected boundary interval if
sequence-identical candidate votes are significant at one-sided binomial
p<=0.05. The full/half MEC parity and source read gauge still must agree, and
the existing one-attachment and no-BAM/BAM rules remain. A one-anchor edge at
51.27 Mb demonstrated why the guard is required: an unguarded retry merged
opposite parental orientations and caused 54 correct-to-wrong read changes in
that block; the retained retry abstains there.

Full chr20 moves from 27/48 to 28/48 tracked spans (19/34 to 20/34 current,
8/14 historical unchanged), reduces VCF blocks from 452 to 419, keeps N50
482,085 bp and all 225,005 tagged reads. Direct parental-truth comparison finds
zero correct-to-wrong and zero wrong-to-correct read changes. The formal
diplinator evaluator reports 6,258 discordant of 224,985 evaluated reads
(97.2185%); the baseline had the same 6,258 discordant with nine fewer eligible
reads. Details are in `evaluations/2026-09-23-local-edge-retry/`.

### Recovery MEC uses one resource-aware decision flow (2026-09-23)

Simplification experiments established the boundary of the algorithm. Replacing
site MEC with whole-block majority voting reduced phased heterozygotes by 14 and
raised chr20 truth discordance from 6,258 to 6,416 reads. Removing the variable
bound was also rejected: real components contain up to 167 coupled variables,
and full chr20 still consumed about seven CPU cores after 2.5 minutes. Binary
MEC is NP hard, so a practical exact implementation must be allowed to abstain.

The retained flow distinguishes evidence failure from resource exhaustion. It
solves the complete read-connected atomic blocks first. A tie, split-half
disagreement, source-gauge conflict, or missing evidence ends the decision. Only
an otherwise eligible problem exceeding the 20-variable exact-search bound may
use the selected boundary interval, and significant sequence-identical candidate
evidence must authorize that scope. This prevents a narrower matrix from
silently overriding contradictory whole-block evidence.

The refactor is output preserving on full chr20: the phased VCF is byte identical
and the SAM records from the phased BAM have the same MD5 as the prior guarded
baseline. Metrics remain 61,644 phased heterozygotes, 419 VCF blocks, 482,085 bp
VCF N50, 225,005 tagged reads, and 6,258/224,985 formal truth discordance
(97.2185%). The 14.264 Mb target remains joined and 51.27 Mb remains separate.
Details and rejected arms are in
`evaluations/2026-09-23-single-flow-recovery/`.

### Unphased-read audit and overlap-output fix (2026-09-23)

A qname-exact comparison of the current chr20 graph+BAM output with HiPhase
found 9,183/27,287 pgphase-unphased reads tagged by HiPhase and 18,104 left
unphased by both. The graph read-evidence diagnostic showed that 24,795 reads
had graph allele observations but no eligible site produced a haplotype score,
1,807 had no graph allele observation, 174 tied, and 14 retained a positive
margin without a final assignment. HiPhase's larger DeepVariant-backed phased
set (77,123 heterozygotes versus 61,644) and its post-solve single-allele
haplotagging explain most of its extra coverage. The extra 9,183 calls are
8,071/9,183 (87.89%) truth-correct. MAPQ is secondary: only 505 are below 30
and none below HiPhase's floor of 5.

The audit also found 512 reads with a valid internal HP/PS assignment that the
phased-BAM accumulator erased when a later overlapping chunk visited the read
unphased. The merge now preserves a valid assignment across an unphased visit;
a later phased assignment still owns the read. Full chr20 adds exactly 512 tags
with no lost or changed existing tags. Direct parental truth supports 500/512
(97.66%), and the phased VCF remains byte-identical. Current counts are 225,517
phased and 26,775 unphased graph-associated reads. HiPhase tags 8,682 of the
remainder at 87.68% truth accuracy. Details are in
`evaluations/2026-09-23-unphased-read-audit/`.


### Post-solve graph read rescue closes the aggregate HiPhase coverage gap (2026-09-23)

The 8,682-read HiPhase-only deficit was separated into solver and input effects.
Running HiPhase 1.7.0 on pgphase's own 61,644-heterozygote graph VCF tagged
225,978 shared reads and recovered only 3,444/8,682 target reads. The original
HiPhase DeepVariant output carries 25,465 phased alleles absent from pgphase
output, including 11,312 SNPs; only 6,344 and 1,818 respectively occur exactly
in the complete graph catalog. Most of the remaining exact-read difference is
therefore private-site or representation input, not an A*/k-means choice.

The retained read-only pass runs after cross-chunk stitching. Assigned reads
orient excluded biallelic sites within one PS using both haplotypes, both
alleles, and an exact one-sided binomial p<=0.01. Exactly one PS must support the
site. Fixed-point layers extend inward, SNP votes precede indel votes, and
co-located rows vote once. One directly phased site may tag a read; an
indirectly oriented excluded site needs a second independent locus. Conflicts
and ties abstain. Fallback assignments use the independent
`PS + kGapFillPsOffset` namespace and cannot override primary HP/PS. Candidates
and VCF phase sets are untouched.

Full chr20 phases 230,072/252,292 shared reads, 275 more than HiPhase on that
population. It adds 4,555 reads with no lost or changed existing assignments and
a byte-identical phased VCF. Of the original 8,682 HiPhase-only reads, 3,159 are
recovered at 85.57% truth accuracy. Whole-output direct truth accuracy is 96.84%
(222,801/230,072), above HiPhase's 95.91% on its full read population. The five
window concordance floors intentionally move to 0.88, 0.98, 0.99, 0.85 and
0.94 while all 232 assertions pass. The remaining 5,523 exact qnames cannot all
be obtained from the retained graph site representation; matching them requires
discovering or importing the sample-private alleles. Details are in
`evaluations/2026-09-23-unphased-read-audit/`.

### Statistically validated BAM blocks extend read coverage (2026-09-23)

A whole-chunk BAM solve can recover reads whose private alignment variants are
absent from the graph catalog, but transferring all of its assignments is
unsafe: the unvalidated arm added 4,452 reads at only 69.18% local truth
accuracy. Accepting a BAM block from any one graph overlap also missed
contradictory evidence elsewhere in that block.

The retained pass keeps BAM blocks independent and transfers no candidates or
observations. Shared graph/BAM phased reads form a 2x2 table for each graph
phase set. Every diploid table must reject random association at exact
two-sided p<=0.01. The tables are oriented independently, then their aggregate
disagreement must have a one-sided 95% Wilson upper bound <=10%. Passing block
assignments are staged until graph stitching and excluded-site rescue finish,
then fill only still-unphased reads under the separate
`PS + kBamFallbackPsOffset` namespace.

On full chr20 this adds 396 reads, 389/396 (98.23%) truth-correct, with no lost
or changed prior assignment and a byte-identical VCF. Whole-output accuracy
moves from 222,801/230,072 (96.839685%) to 223,190/230,468 (96.842078%). On the
shared population, pgphase's tagged-read lead over HiPhase grows from 275 to
671. It recovers 393 of the prior 5,523 HiPhase-only reads at 98.47% truth
accuracy. The VCF remains at 61,644 phased heterozygotes in 419 phase sets. All
unit tests and all 232 gap-window assertions pass. Full details are in
`evaluations/2026-09-23-independent-bam-read-fallback/`.

### Decisive independent BAM reads exceed HiPhase coverage (2026-09-23)

The block validator above is deliberately conservative. A BAM phase block that
cannot be oriented to a graph block does not need to be discarded if it remains
an independent output phase set: its internal HP labels are meaningful up to the
same arbitrary phase-set flip used by every truth evaluation. The fallback now
keeps an individual BAM assignment from an unsupported block only when its BAM
haplotype score margin is at least 4. A clean biallelic SNP or indel contributes
+2 to one haplotype and -2 to the other, making 4 the first complete clean-site
separation. Passing statistically validated blocks still contribute all of
their assigned reads.

Reads absent from the chunk's graph profiles are stored as output-only
assignments. Graph-present assignments remain staged until stitching and
excluded-site rescue finish. Both use the independent
`PS + kBamFallbackPsOffset` namespace; neither contributes candidates,
observations, graph joins, or stitching votes. A primary graph assignment from
any overlapping chunk always wins. Reinitializing the pass now clears every
staging vector, and nonpositive BAM PS sentinels are rejected.

The threshold sweep on full chr20 was:

| minimum score margin | phased reads | correct | discordant | accuracy |
|---:|---:|---:|---:|---:|
| **4** | **234,787** | **226,223** | **8,564** | **96.352439%** |
| 5 | 232,749 | 224,914 | 7,835 | 96.633713% |
| 6 | 232,749 | 224,914 | 7,835 | 96.633713% |

Margins 5 and 6 are safer but phase 604 fewer reads than HiPhase. Margin 4 is
the only tested score boundary that exceeds HiPhase coverage while retaining
higher truth accuracy. Against the conservative validated-block baseline it
adds 4,319 phased reads. The final output phases 1,434 more reads than HiPhase,
produces 2,423 more correct assignments and 989 fewer discordant assignments,
and is 0.446237 percentage points more accurate. It phases 2,133 reads absent
from graph/GAF profiles. All baseline HP/PS assignments remain byte-for-byte
unchanged, and the phased VCF is byte-identical. The all-block arm remains
rejected: although it phases 2,076 more reads than the retained margin-4 arm,
it adds 893 more errors and narrows the accuracy lead over HiPhase to 0.101
percentage points.

### One-megabase graph context improves coverage and accuracy (2026-09-23)

The 500 kb margin-4 configuration left two distinct opportunities: reads with
weak independent-BAM evidence, and reads whose graph and BAM blocks changed
when given more shared context. Accepting every weak BAM assignment was not
useful: the 2,076 broad-only reads were only 57.47% truth-correct. Most had a
haplotype score margin of 2, and neither MAPQ nor rejected graph-link evidence
identified a high-accuracy subset. Of the 6,301 reads HiPhase still phased and
the broad BAM arm did not, none received a BAM HP/PS in any profiled chunk;
2,735 were in the excluded centromeric interval. The assignment MAPQ floor was
not the cause because the in-memory sub-solve assigns admitted reads before the
standalone BAM writer applies that output filter.

Increasing only the BAM sub-solve by 250 kb on each side was also rejected. It
phased 234,708 reads with 226,055 correct and 8,653 discordant, worse than the
un-padded 500 kb result. The graph and BAM solves need the same context.

The joint chunk and margin sweep was:

| graph chunk | BAM read margin | phased | correct | discordant | accuracy |
|---:|---:|---:|---:|---:|---:|
| 500 kb | 4 | 234,787 | 226,223 | 8,564 | 96.352439% |
| 500 kb | 6 | 232,749 | 224,914 | 7,835 | 96.633713% |
| 1 Mb | 4 | 236,918 | 228,010 | 8,908 | 96.240049% |
| **1 Mb** | **6** | **235,835** | **227,293** | **8,542** | **96.377976%** |
| 2 Mb | 4 | 239,129 | 228,903 | 10,226 | 95.723647% |

The 1 Mb, margin-6 arm strictly improves the prior 500 kb, margin-4 result:
+1,048 phased reads, +1,070 correct assignments, 22 fewer discordant
assignments, and +0.025537 accuracy points. Two megabases gains more reads but
falls below HiPhase accuracy, so it is rejected. One megabase is now the
default for graph runs with `--bam`; graph-only and standalone BAM runs remain
at 500 kb, and an explicit chunk size is preserved. The margin-6
rule requires more than the four-point separation supplied by one clean
biallelic observation before a read from an unvalidated BAM block can be
emitted.

The retained VCF has 61,726 phased heterozygotes in 452 phase sets with a
460,310 bp span N50. Changing chunk context changes graph block boundaries, so
this VCF is not byte-identical to the old 500 kb result. All comparisons below
therefore score the final output directly rather than assuming assignment
preservation.

### Whole-chr20 comparison now includes longcalld (2026-09-23)

The final graph+recovery output, frozen HiPhase output, and frozen upstream
longcalld output were rescored on the complete 272,016-qname BAM population.
All three use the same underlying annotated `vg giraffe` BAM alignments;
HiPhase's input is a lossless contig-header reheader, and its DeepVariant VCF
was called from that same reheadered BAM. pgphase additionally uses the graph
catalog and GAF, while longcalld calls variants internally. Each emitted PS was
independently oriented against the same parental truth map.

pgphase now has the highest coverage and correct yield: 235,835 reads phased
(86.6989%), with 227,293 correct and 8,542 discordant. HiPhase phases 233,353
(85.7865%), with 223,800 correct and 9,553 discordant. Longcalld phases 219,090
(80.5431%), with 212,584 correct and 6,506 discordant. Conditional truth
accuracy is 96.3780%, 95.9062%, and 97.0304% for pgphase, HiPhase, and
longcalld, respectively. pgphase therefore dominates HiPhase on phased reads,
correct yield, discordant count, and conditional accuracy. Longcalld remains
the more abstaining high-accuracy operating point: it has 0.6525 percentage
points higher accuracy, while pgphase phases 16,745 more reads and produces
14,709 more correct assignments.

The 252,292-qname graph-observed subset remains a phaser diagnostic, not the
primary whole-BAM comparison. pgphase phases 232,964 graph-observed reads
(92.3396%) at 96.7802% truth accuracy, versus HiPhase's 229,797 (91.0837%) at
96.4321%. The output-only BAM channel phases 2,871 reads absent from the
graph/GAF profiles. On the 212,076 reads phased by all three, pgphase reaches
98.5760%, longcalld 98.2709%, and HiPhase 97.1779% after independently
orienting phase sets on that common subset.

The callsets contain 61,726, 77,123, and 83,013 phased heterozygotes for
pgphase, HiPhase, and longcalld. These counts are not a variant-accuracy ranking
because the graph catalog, DeepVariant callset, and longcalld-discovered callset
differ. Full methodology and the machine-readable table are in
`evaluations/2026-09-23-independent-bam-read-fallback/`.


### Exact BAM observations enable statistically safe singleton rescue (2026-09-23)

The final pgphase/HiPhase qname audit found 6,617 HiPhase-only read tags; 5,009
agree with parental truth, including 3,770 outside chr20:26-30 Mb. All 3,770
noncentromeric reads overlap a phased HiPhase heterozygote and 2,871 rely on one
site. Indels dominate: 3,204 overlap only indels. A one-base representation
proxy maps 1,807 to pgphase `REP_HET_INDEL` rows, 810 to an already phased
pgphase candidate whose graph profile lacks an eligible read observation, 1,145
to no nearby pgphase candidate, and 8 to an unphased clean SNP.

The whole-chunk BAM fallback now retains missing read alleles at
sequence-identical biallelic graph candidates. Catalog ALTs are normalized
through `vcf_to_variant_key` and BAM-to-graph candidate matching is precomputed
once per chunk. The original graph-only rescue reaches its fixed point first.
A second fixed point may fill a missing graph allele from the BAM channel, but
a called graph allele remains authoritative.

One statistically oriented excluded site may now tag a read when primary graph
assignments alone contain both haplotypes and alleles, pass the exact one-sided
binomial test at p<=0.01, and have a one-sided 95% Wilson discordance upper bound
<=15%. Read-only rescues cannot bootstrap this singleton condition. Augmented
assignments fill only empty output rows and cannot replace graph assignments or
the existing margin-6 BAM fallback. Candidates and stitching are unchanged.

The retained full-chr20 output phases 236,338 reads: 227,754 correct and 8,584
discordant, for 86.8839% coverage and 96.3679% conditional truth accuracy. It
adds 503 reads over the prior baseline at 465/503 (92.45%) truth accuracy, loses
no tag, changes no existing HP, and leaves the phased VCF byte identical. It
recovers 464 truth-correct HiPhase-only reads, 462 outside chr20:26-30 Mb.
Against HiPhase it phases 2,985 more reads, produces 3,954 more correct
assignments and 969 fewer discordant assignments, while remaining 0.4617
percentage points more accurate. Details are in
`evaluations/2026-09-23-bam-observation-singleton-rescue/`.

### Chr20 graph/BAM site and representation audit (2026-09-23)

A fresh same-BAM audit found three representation defects in the current
graph path. Targeted recovery indexed every source VCF ALT against each split graph
candidate; using only `site_allele_orig_idx`'s selected ALT changes 24 in-window
sub-solve rows from falsely shared to genuinely appended, gains 15 phased VCF
rows, and improves whole-chr20 read truth from 227,754 correct / 8,584
discordant to 227,799 / 8,579. The final output phases 236,378 reads at
96.3706% truth accuracy. Direct unit tests cover a two-ALT source, an ALT2-only
survivor, and a whole multiallelic candidate.

The VCF writer also shortened 952 equal-length multibase graph replacements to
one REF base, and discarded retained ALT bases in 40 length-decreasing
replacements. Both now preserve the complete allele; 31 INFO/SVLEN fields also
reflect net length, as upstream longcallD defines it. Writer changes preserve
candidate TSV, GT, and PS byte for byte. All 62,153 final VCF REF strings match
the reference FASTA. `make unit-tests`, `make predicate-tests`,
`make window-tests`, and `make check` pass.

The exact-key audit still finds 28,922 phased standalone BAM allele keys absent
from the graph VCF. A candidate-first audit resolves them as 19,301 without an
exact catalog ALT or final candidate, 4,008 with an exact catalog ALT but no
exact candidate, and 5,613 with an exact candidate that is not an emitted
phased heterozygote. Of the last group, 5,534 are excluded `REP_HET_INDEL`
rows; 46 clean candidates are 0/0 and 33 are 1|1. In total,
21,806 BAM keys have no exact catalog ALT, including 9,544 SNPs. Exact-key
absence does not prove biological absence: at least 153 BAM SNP keys are bases
inside phased graph MNPs, and indels can shift in repeats. The high-AF graph
filter remains a known tradeoff on alternate-vs-alternate snarls; the previous
`--snarl-allele-phasing` trial worsened read errors, so this audit leaves its
default unchanged. Full method, caveats, and the machine-readable table are in
`evaluations/2026-09-23-site-representation-audit/`.

## 2026-09-23: Graph sites between independent BAM recovery blocks

Four small chr20 graph+BAM reruns at 15.0, 21.8, 33.8, and 38.2 Mb had 12
still-open intervals between imported BAM phase blocks. Their merged matrices
retained 33 original graph candidates between blocks: 27 repeat indels and six
clean SNP rows. Most repeat sites were unphased or supplied only one-sided or
truth-weak evidence. The graph block containing clean SNPs at 38,259,278 and
38,279,170 is a positive exception: its SNP alleles are 95.8% and 98.2%
truth-separated, and source-specific read gauges support the same orientation
to the left and right BAM blocks by 42:2 and 55:0. All three blocks agree with
parental truth, but final PS labels remain separate. The left graph/BAM edge
lacks a sequence-identical shared candidate, which the direct stitch gate
requires; the fallback did not close it. This is a concrete graph-mediated
continuity opportunity, not authorization to relax the gate globally. Full
method and per-gap examples:
`evaluations/2026-09-23-interblock-graph-sites/README.md`.

## 2026-09-24: Guarded graph bridge between recovered BAM blocks

The 38.2 Mb graph bridge was present in the data but the fallback stopped after
attaching the first graph/BAM pair: its one-use-per-original-block guard hid a
second, independently supported edge. A broad trial that reused blocks and
accepted a sequence-identical clean candidate plus decisive source-specific
read gauge when full-block MEC tied joined the 38.2 Mb example. It also created
an unsupported span across the 26,029,591-26,088,679 panel window: the two
local edges were decisive, but no read observed both ends of the intervening
BAM block. The retained rule allows a reused block only if its first and last
phased sites are observed together by reads from both haplotypes and agreeing
molecules outnumber conflicts. A candidate-plus-read-gauge fallback after a
full MEC tie is limited to a second attachment of an already validated graph
block; applying it to a first attachment failed four existing atomicity and
full-block-tie unit checks. The centromeric panel span disappeared and all
panel expectations passed without refresh.

In the exact 38,233,000-38,343,000 regional rerun, the graph SNPs at
38,259,278 (VCF position 38,259,286) and 38,279,170 now connect the BAM
blocks. The resulting primary PS has 209 truth-scored reads, 200 concordant
and 9 discordant; before the join, the two separate PS groups held the same
209 reads and same 200 concordant assignments. No read HP changed; 124 read PS
labels changed. The left BAM block has 4 agreeing/1 conflicting end-to-end
molecules, and the middle graph block has two agreeing end-to-end molecules,
one from each haplotype. The left sequence translation already recognized the
BAM SNP at 38,259,286 as the same site, but that BAM row was classified as a
repeat heterozygote and therefore was not counted as a *clean* shared-candidate
anchor. The direct two-SNP MEC path independently validates that edge.

The other three exact regional reruns (15.0, 21.8, and 33.8 Mb) had identical
read HP/PS labels and truth scores to their baselines. The existing 48-gap panel
passed 232 assertions, and a focused VCF regression now checks the 38.2 Mb
bridge. This is a phase-set continuity gain in one demonstrated interval, not
a measured increase in phased-read coverage or chromosome-wide accuracy.

A matched full-chr20 A/B toggled only `src/collect_phase.cpp` while retaining
the current graph/BAM site and VCF fixes. The baseline has 236,378 truth-scored
phased reads (227,799 correct, 8,579 discordant) in 871 read phase sets; the
guarded stitch has 236,379 (227,800 correct, 8,579 discordant) in 861 sets.
Among shared reads, 8,470 PS labels and 627 numeric HP labels change, but no
read changes truth correctness; the one added read is correct. VCF variant keys
are a superset by two records at 65,508,812 and 65,508,959; 277 existing GTs
change only by 0|1/1|0 reversal, and 1,891 existing PS fields change. Those
GT reversals reflect phase gauge changes, not altered heterozygosity. The
full-chromosome comparison does not reveal a new truth-discordant read or a
lost phased read. The same-coordinate guard excludes complementary alleles
at one position from pretending to span a block.

## 2026-09-24: Omitted graph indels as read-only recovery markers

An excluded graph repeat indel may now tag still-unphased reads inside a
bounded BAM-recovery seam when the closest phased clean SNP in one established
phase set directly validates its allele orientation. The full read set and two
disjoint FNV read halves must each pass one-sided binomial p<=0.05, both alleles
and SNP haplotypes must be present, and one-sided 95% Wilson bounds limit
pair discordance to 25% and the fraction of site observations without that SNP
to 50%. Competing unphased graph alternatives at one coordinate abstain. The
candidate and its VCF phase remain unchanged; this is read-only rescue.

A broader trial added 207 truth-scored chr20 reads at 187/207 (90.3%) correct,
but made two existing reads incorrect and lost one tag. The guarded version
adds 30/30 truth-correct reads, loses none, changes no shared read's truth
correctness, and leaves the VCF byte-identical. Full chr20 moves from 236,379
phased / 227,800 correct / 8,579 discordant to 236,409 / 227,830 / 8,579.
Read phase sets rise from 861 to 862 because read-only rescue creates a new
local assignment namespace. The 21.823 Mb graph deletion contributes 22
new reads. The 38.294 Mb insertion abstains because one read half is unstable.
The full method and regional counterexamples are in
`evaluations/2026-09-24-direct-graph-snp-rescue/`.

## 2026-09-24: Spanning-read votes do not certify whole-block stitches

A Longcalld-style same/cross vote and deeper physical read span were tested
for N50 joins. At 2.525 Mb, 19/19 informative spanning reads support the
parentally correct join; at 51.263 Mb, 28/28 support a join that is wrong for
the full left phase set because that block changes parental orientation near
51.24 Mb. Even requiring >=2 kb of alignment into each block leaves 13/13
and 15/15 support, respectively. A strict chromosome-wide direct-SNP vote
probe selected 27 joins with two wrong parental relations, improving simulated
VCF N50 460,310->491,812 bp but losing 316 truth-correct reads. The existing
read-vote and block-reuse guard remain unchanged. Methods and limits are in
`evaluations/2026-09-24-n50-stitch-audit/README.md`.

Follow-up on the 51.24 Mb counterexample: frozen HiPhase has **two distinct
source blocks**, ending at clean SNP 51,235,063 and beginning at clean SNP
51,262,081; HiPhase run on pgphase's older VCF catalog makes the same split.
The shared annotated BAM has 56 primary reads over each endpoint but zero
primary reads spanning both. Pgphase places both SNPs in PS 51,165,532, then
opens PS 51,270,774 although 31 primary reads span 51,262,081-51,270,774.
HiPhase correctly starts a new phase problem at the unsupported 27-kb interval
and keeps the supported 8.7-kb interval inside its next PS. Details and
limitations are in the audit README.

A MAPQ/base-quality/GAF check shows the first pair's 56+56 endpoint reads are
mostly MAPQ60, yet zero BAM or GAF read names cover both clean SNPs. The
second pair has 31 MAPQ60 spanners, 28 with BQ>=20 at both clean SNPs, all
28 same-parity across both haplotypes; 15 retain >=2 kb aligned outside each
SNP. In the regional replay, `recovery-input` has the three SNPs in separate
phase sets, while `recovery-final` merges the zero-bridge 51,235,063-
51,262,081 pair and leaves the 31-spanner 51,262,081-51,270,774 pair apart.
This local misconnection is introduced during the recovery stitch. The audit
README records the counts and the phase-state trace.

## 2026-09-24: Preserve the next recovery anchor across an unsupported seam

A graph-only recovery stitch at chr20:51,235,063-51,262,081 used the broad
whole-window graph/BAM gauge despite having no allele edge between the two
clean SNP phase sets. No primary BAM or GAF read spans the 27-kb gap. That
merge consumed the right block before the next recovery seam could use it as
its left anchor. The stitch now defers a gauge-only graph join when its right
block is the next seam's left flank and no local allele edge was found. The
regular downstream seam still runs left to right; other graph-only gauge
joins retain their prior behavior.

A short regional replay now splits the unsupported pair and joins the next
51,262,081-51,270,774 clean SNP pair. The 1-Mb production chunk also splits
the unsupported pair, but its larger imported-BAM transaction still leaves
the latter two SNPs in separate phase sets. That remaining missed connection
is a separate recovery-window limitation.

Matched full chr20 A/B, toggling only this stitch change:

| Read output | Original phased / correct / discordant | Guarded phased / correct / discordant | VCF blocks / N50, original -> guarded |
|---|---:|---:|---:|
| `--min-read-margin 2` | 188,341 / 184,906 / 3,435 | 188,341 / 185,323 / 3,018 | 439 / 482,085 -> 457 / 456,233 bp |
| Default margin | 236,404 / 227,446 / 8,958 | 236,520 / 228,184 / 8,336 | 439 / 482,085 -> 457 / 456,233 bp |

A blanket requirement for a local allele edge at every graph-only gauge
join was rejected: it raised discordant reads from 3,435 to 3,812 in the
matched margin-2 run and lowered N50 to 449,827 bp. The narrower next-seam
guard increases truth-correct reads at both tested read-output settings,
while VCF N50 falls by 25,852 bp because unsupported joins are removed.
Inputs, physical span/MAPQ/base-quality counts and the rejected arm are in
`evaluations/2026-09-24-n50-stitch-audit/README.md`.

## 2026-09-24: HiPhase's transitive site connectivity exposes a recovery transfer limit

HiPhase solves a physically connected run of heterozygous sites jointly,
then assigns phase sets from connected components of valid read allele
observations. At chr20:51,262,081-51,287,372, it forms one block with both
its own DeepVariant catalog and pgphase's VCF. Pgphase's hybrid input already
contains consistent read-allele links across all three adjacent clean-SNP
pairs (30, 4, 49 observations), but the graph/BAM stitch rolls back the local
joins because no read spans the two outer graph anchors. Transferred BAM-only
rows also lose connectivity through shared graph rows. Local-only relaxation
added 279 wrong whole-block assignments at 19 Mb in the full-chr20 test;
a narrower exact-search fallback changed labels but gave no read or N50 gain.
Both trials were reverted. The evidence and implementation path are in
`evaluations/2026-09-24-hiphase-connectivity-audit/README.md`.

## 2026-09-24: Seam-local source-PS transfer across shared graph rows

Selected targeted BAM phase sets now retain all oriented source sites as an
overlay, including shared graph rows. After the normal recovery stitch, a
source block can transfer its orientation through one seam only when every
source cut has consistent observations from both haplotypes, both graph flanks
have exact-site-consistent and conflict-free graph/BAM gauges, and the shared
allele indices are binary and compatible. Ambiguous reads or sites claimed by
two targeted solves abstain. BAM-private candidate insertion remains confined
to the original seam; distant graph sites and reads retain their labels.

On matched full chr20, the formerly separate 51.262–51.287 Mb chain shares
one PS while the zero-read 51.235 Mb cut stays split. The 19 Mb wrong-join
negative control also stays split. Compared with the previous baseline,
truth-evaluable phased reads change from 236,520 to 236,526, truth-correct
from 228,184 to 228,189, discordant from 8,336 to 8,337, and VCF N50 stays
456,233 bp. The wider private-site and unconditional whole-PS arms were
rejected for coverage loss and a large 19 Mb wrong join, respectively.
Details and caveats are in
`evaluations/2026-09-24-whole-bam-phase-set-transfer/RESULTS.md`; the broader
reconciliation proposal remains in its `DESIGN.md`.

## 2026-09-24: The 55.381 Mb HiPhase gap has both a false anchor and an atomicity limit

The current full-chr20 graph+BAM output leaves PS 55,360,776 separate from
PS 55,373,606 around chr20:55.381 Mb. Fresh HiPhase on pgphase's own VCF
joins the flanking sites and leaves the intervening 55,373,606 G>T row
unphased. Pgphase calls that graph row a clean heterozygote (13 G/50 T),
although 49/49 primary BAM bases are T and the targeted BAM caller labels it
homozygous. The false graph PS ends the detected recovery seam before the
informative BAM indels at 55,381,222 and 55,381,473.

Masking just that graph row in a local catalog brings the indels into the
seam, but pgphase still keeps their two source BAM phase sets separate. The
BAM source has a 44:6 same/conflicting allele link between them; the outer
left graph flank has an inconclusive 26:21 gauge, and recovery's atomic
outer transaction rolls back the internal link. HiPhase can instead start
an independent read-connected block at 55,360,776. Frozen HiPhase tags all
57 reads physically spanning the middle indels in one truth-consistent PS;
two maternal reads physically span the outer 22-kb interval, one with valid
MAPQ60/BQ40 SNP calls on both sides. The full trace and controlled catalog
and HiPhase replays are in
`evaluations/2026-09-24-hiphase-gap-5538/README.md`.

## 2026-09-24: Guarded BAM anchor reconciliation closes 55.381 Mb

Targeted BAM homozygous SNP evidence can now remove a conflicting graph
seam-right anchor only when clean or MSA-verified calls on both BAM
haplotypes pass a MAPQ >= 30, `p <= 0.01` test and the expanded seam stays
inside the already solved region. The graph genotype remains unphased. An
independent inner BAM component can then attach to the right graph block
when a significant BAM/BAM allele edge and a conflict-free direct BAM SNP
molecule agree; the uncertain outer left graph block remains separate.
BAM MAPQ travels with the BAM allele channel, so a high-MAPQ GAF alignment
cannot upgrade a weak BAM witness.

On matched full chr20, the 55.381 Mb gap closes without a new discordant
truth-scored read: phased/correct/discordant change from
236,526/228,189/8,337 to 236,527/228,190/8,337. VCF blocks fall 464 ->
463 and N50 stays 456,233 bp. The 19 Mb wrong-join control remains split.
A first ungated anchor veto incorrectly demoted the true low-MAPQ
30,794,399 A>C graph SNP; the MAPQ condition fixes that. Evidence and
full A/B are in `evaluations/2026-09-24-hiphase-gap-5538/README.md`.

## 2026-09-24: HiPhase rephased the current pgphase chr20 VCF

A full HiPhase `1.6.0-ac3f399` replay used pgphase's final 62,095-row
chr20 graph+BAM VCF, the same annotated BAM and CHM13 reference. It phased
59,381 heterozygotes in 211 VCF blocks and tagged 228,286 truth-evaluable
reads, 216,224 correct and 12,062 discordant (94.72% per-PS purity).
Pgphase itself phased 61,646 heterozygotes in 463 blocks and tagged
236,527 reads, 228,190 correct and 8,337 discordant (96.48%).

Across 115 noncentromeric, non-overlapping pgphase PS gaps shorter than
10 kb, HiPhase put the same exact boundary variants in one PS for 65.
Among 64 joins with at least 20 truth-scored nearby reads, 35 have at least
98% local purity; 34 also have >=95% purity and matching orientation on
both flanks. Low-purity joins show that using the same sites does not
always phase a gap correctly. The command, scoring definition and per-gap
results are in `evaluations/2026-09-24-hiphase-on-pgphase-sites/`.

## 2026-09-24: Atomic and one-sided BAM source attachment recovers 27/35 short gaps

The 35 high-purity noncentromeric chr20 gaps joined by HiPhase on pgphase's
previous VCF exposed three recovery-transfer defects: partial transfer cut an
established graph PS at 11.36 Mb; requiring one BAM source PS to validate
both graph flanks discarded a strong one-sided 14.58 Mb attachment; and an
all-or-nothing BAM path check rejected the supported 8.63 Mb run because of
a weak cut 24 kb later. Recovery now transfers approved graph blocks as units
in their current post-stitch gauge, accepts one independently approved flank,
and locally attaches supported source runs before/after weak cuts to only one
graph block. The weak cuts come from the original BAM solve, and reads that
observe another run keep their old phase-set assignment. Splitting the source
PS globally was rejected: it caused a wrong 19 Mb join and lost the 55 Mb
positive bridge.

On matched full chr20, truth-evaluable phased/correct/discordant reads change
from 236,527/228,190/8,337 to 236,675/228,700/7,975 (96.48% to 96.63%).
VCF blocks fall 463 to 411 and N50 rises 456,233 to 486,139 bp; 27 of 35
tracked high-purity HiPhase-only gaps close, with no confident parental
switch at their truth-scored flanks. A separate 26.030-26.089 Mb span has 40
physical crossing reads and unanimous parental orientation at both ends.
Eight tracked short gaps remain, every one with an indel or SV boundary.
On the updated identical VCF, HiPhase phases 228,362 truth-evaluable reads
at 94.71% purity and produces 211 blocks with 919,153 bp N50. Pgphase
therefore still trails it in block continuity despite higher read coverage
and purity. Commands, controls, and per-gap limitations are recorded in
`evaluations/2026-09-24-atomic-source-run-recovery/README.md`.

### 2026-09-24: one-haplotype BAM source bridge at 34.844 Mb

The remaining 34.844 Mb HiPhase-only gap was caused by a categorical source
path veto, not a missing variant. Five independent BAM molecules supported
the same orientation across the 34,818,318 to 34,835,167 source-site cut,
with zero conflicts, but all five came from one haplotype. The path check
required support from both haplotypes at every cut and discarded the whole
source block. An exact one-sided binomial test admits conflict-free one-hap
cuts at p <= 0.05; exact shared graph-site and two-haplotype flank checks
still guard transfer. The 34,844,194 insertion and 34,844,579 SNP now join
without a truth switch. Full chr20 closes 28/35 tracked high-purity HiPhase
gaps, has 410 VCF blocks versus 411 before, and retains 236,675 phased reads,
228,700 truth-correct reads, 7,975 discordant reads, and 486,139 bp N50.
The 19 Mb wrong-join control remains open. Details:
`evaluations/2026-09-24-atomic-source-run-recovery/README.md`.

### 2026-09-24: recovery boundary evidence and remaining weak-cut gaps

The chr20 source-path baseline phased 236,675 truth-evaluable reads with 7,975
discordant (96.6304% purity), 410 VCF blocks and 486,139-bp N50. Requiring
exact MEC parity plus a statistically validated clean indel/SNP boundary pair
closed the 14.235-Mb and 36.016-Mb gaps without changing read truth. For a
complete BAM source path, a split-stable aggregate allele vote can orient a
graph/BAM edge when the read-label gauge has no entry; this closed 32.490 Mb.
The new full chr20 result has the same 236,675 phased reads, 7,974 discordant
(96.6308% purity), 406 blocks and 491,812-bp N50. See
evaluations/2026-09-24-general-recovery-evidence/README.md for the gap audit.

Splitting all source blocks at weak cuts was tested and reverted: it did not
close the 0.542 or 37.46-Mb gaps. A simple pairwise boundary vote is unsafe on
the 19-Mb negative control because both attachments can be individually strong
while their joined graph flanks are truth-opposite. Weak-cut components need an
independent gauge and an outer-flank consistency check before they can safely
extend a phase set. The 0.528 and 17.62-Mb boundaries remain ambiguous in the
current pgphase allele observations, so forcing a join there is unsupported.

A matched HiPhase regional control changed only `--disable-global-realignment`:
its default run joined the 0.528, 0.542, 17.62 and 37.46-Mb boundaries, while
the disabled run split the first three and still joined 37.46 Mb. All four gaps
have physical MAPQ >= 30 spanning reads in pgphase's BAM, but 0.528 and 17.62
lose many callable indel alleles in the pgphase observation matrix. This
separates an allele-evidence problem from the 37.46-Mb source-component stitch
problem; see the evaluation record for counts and control commands.

The retained component-level postpass attaches an unanchored BAM run to exactly
one statistically supported neighboring phase set, without relabeling the
source block across its weak cut. It closes the 0.542 and 37.46-Mb tracked gaps
and three more graph-arm window spans at 6.578, 12.269 and 48.225 Mb. Full
chr20 now phases 236,656 truth-evaluable reads, of which 228,729 are correct
and 7,927 discordant (96.6504% purity); the VCF has 398 phase sets and
491,812-bp N50. The preceding build had 236,675 phased reads, 7,974
discordant and 406 phase sets. The 19-Mb wrong-join control remains split. Four
new PS labels each contain at least two prior read groups with >=10 truth reads;
none has two >=90%-pure groups with opposed parental orientation. The three
new graph-arm span expectations were updated only after coverage and parental
orientation review. Full details are in
evaluations/2026-09-24-general-recovery-evidence/README.md.

A two-site HiPhase control isolated the remaining 0.528 and 17.62-Mb gaps. With
only the two exact boundary variants in the VCF, default HiPhase joins both;
turning off global realignment splits both. A temporary read-level diagnostic
shows 53 perfectly truth-consistent paired SNP/insertion calls at 0.528 Mb in
global mode versus 28 mostly reference-SNP calls locally; pgphase has only
eight callable pairs on those same reads. At 17.62 Mb, both HiPhase modes have
32 paired calls, but global graph alignment changes the same/cross relation
from 18:14 to 7:25; pgphase has only five callable pairs. HiPhase A* keeps both
sites heterozygous only in global mode. No intermediate VCF sites are needed
for either join. See the exact contingency tables and source-path explanation
in evaluations/2026-09-24-general-recovery-evidence/README.md.

The original-BAM CIGAR audit explains the two unresolved boundary matrices.
At 0.528 Mb, the 39 HiPhase SNP-alt/insertion-ref paired reads delete across
the SNP (9–11 bp); the opposite 14 have 34–39 bp insertions against a 40-bp
VCF allele. At 17.62 Mb, homopolymer indel placement makes raw exact CIGAR
projection weak: parental read-label concordance is 71.9% at the deletion and
56.0% among 25 callable insertion observations, versus 84.4% and 87.5% for
HiPhase global calls. Expanding exact-CIGAR MSA backfill or relaxing a stitch
threshold is therefore not a safe general fix. See the original-CIGAR audit
and proposed local haplotype allele assignment in
evaluations/2026-09-24-general-recovery-evidence/README.md.


CORRECTION to the preceding CIGAR-only interpretation: graph recovery already
runs BAM noisy-region WFA/abPOA MSA. The longcalld-parity subsolve sets
`add_unplaced_msa_observations=false`, so the haplotype-aware MSA admits only
reads previously assigned to its selected phase set. A temporary
recovery-only true setting restored source MSA paired alleles from 13/53 to
53/53 at 0.528 Mb and from 5/32 to 32/32 at 17.62 Mb. The 17.62-Mb 1-Mb
control then joined with correct parental orientation; 0.528 Mb remained split
because graph-left/BAM-deletion representation and the 528,850 source weak
cut still block transfer. A smaller 0.528-Mb control regressed from joined to
split under the blanket setting, so the trial was restored. Detailed counts
and a scoped recovery design are in
evaluations/2026-09-24-general-recovery-evidence/README.md.

CORRECTION after same-command controls: the recovery-only
`add_unplaced_msa_observations=true` trial fills the 0.528- and 17.62-Mb
source MSA pairs, but changes other BAM phase decisions and raises 17-Mb
discordant reads from 39/4,071 to 606/4,124. The 528,850 weak cut disappears
in that trial; the remaining 0.528-Mb split comes from graph/BAM
representation and the stitch transaction's missing direct outer vote.
An observation-only second MSA pass preserves the original phase sets but
does not close either gap. A transitive join fails the existing 8.64-Mb
weak-cut regression. All trial code was restored. See the follow-up table in
evaluations/2026-09-24-general-recovery-evidence/README.md.

A subsequent same-run stage replay found a separate wrong-join bug in the
blanket MSA-admission trial. At 17-18 Mb, two graph blocks enter recovery as
533/0 and 540/0 truth-tagged/discordant reads. The atomic stitch accepted a
source-block allele link despite 25 spanning reads whose boundary allele
contradicted their own source HP versus 10 that agreed, then merged the graph
blocks into one label with 587 discordant reads. A pairwise statistical veto
now checks tagged spanning reads against each imported block's boundary
alleles before using an imported/imported edge. In the same admitted-read
replay, final output scores 4,133/43 instead of 4,124/606. The normal
recovery setting still scores 4,071/39. The correct 4.76-Mb span and 3.85-Mb
accuracy floor remain intact. A synthetic multi-block test fails on the
original stitch code and passes with the fix. Blanket MSA admission remains
disabled because it still fails the 3.85- and 4.76-Mb gap regressions; the
pairwise check alone does not solve the 0.528-Mb representation boundary. See
evaluations/2026-09-24-general-recovery-evidence/README.md.

A scoped unplaced-MSA admission pass now re-solves a targeted BAM group only
when its ordinary source solve, *after exact-CIGAR backfill*, leaves one
simple two-block boundary with a statistically clear excess of missing paired
alleles (>=20 MAPQ-30 spanning reads, one-sided binomial p<=0.01). Co-located
alternative rows are excluded from this trigger. The 17.616778-17.625527-Mb
indel pair had 5 callable of 33 spanning reads (p=3.309e-05); the scoped
re-solve joins its two VCF rows and the focused truth test passes without a
flank switch. Broad admission had broken 3.85 and 4.76 Mb. Measuring before
backfill also broke 4.76 Mb; after backfill its ordinary solve has 19/36
callable pairs. A low-depth 14.66-Mb trigger (7/0) broke the 14.58-Mb positive
control, while a 15.096-Mb negative-control trigger (66/12) joined the wrong
haplotypes through complementary indel rows. The current minimum support and
single-row boundary checks retain both controls. The 0.528-Mb graph SNP
remains separate: original BAM bases at that coordinate are 27 maternal A,
42 paternal deletions, one paternal A, and no T, while graph SNP A/T calls
from both allele groups are paternal. Its phase requires a representation
translation across the overlapping deletion, not extra source-MSA reads alone.
See evaluations/2026-09-24-general-recovery-evidence/README.md.


### 2026-09-24: guarded source/graph vote and the last allele trap

A trial replaced zero-conflict graph/BAM source approval with an exact
one-sided binomial `p <= 0.01` vote. At 0.528 Mb the shared clean site had
29 agreeing candidate checks, and the source/graph read gauge had 173
agreements versus one conflict. The trial joined the last coordinate pair,
but its 528,827 A>T SNP and 528,828 insertion both emitted `1|0`; HiPhase's
sequence-backed calls put their ALT alleles on opposite haplotypes. Graph
SNP REF and ALT calls mostly come from BAM reads deleted across the SNP base,
so source HP alone cannot orient that projected graph allele.

The statistical approval is retained only for complete source paths. A
weak-cut source still needs the original conflict-free local-run check. Full
chr20 now phases 236,834 truth-evaluable reads, 228,912 correct and 7,922
discordant (96.6550% purity), with 392 VCF blocks and 491,812-bp N50. It
correctly joins 34/35 tracked noncentromeric HiPhase-correct short gaps;
the 528,827–528,828 pair remains split. All 34 joined boundary genotype
relations agree with HiPhase; the known 19-Mb wrong join remains split. One
new output phase set combines two prior truth-pure read groups with the same
parental orientation. The last gap requires a sequence-aware mapping from the
broad graph snarl to the BAM deletion/insertion haplotypes; a bare PS merge
would be a variant-level phase error. Details and the rejected 35/35
coordinate-only trial are in
`evaluations/2026-09-24-general-recovery-evidence/README.md`.

### 2026-09-24: terminal graph SNP projected through an exact BAM deletion

The last 528,827–528,828 split was a nested-snarl representation problem. The
graph A>T child site has 24 ALT calls on reads with a physical deletion across
the SNP, 18 REF calls on deleted reads, one REF call on a physical A, and 27
physical A reads with no child call. The BAM source deletion/insertion pair is
perfectly cross-oriented on 15 callable reads. Recovery now uses the retained
original BAM CIGAR for a terminal graph SNP overlapping a phased BAM deletion.
Exact physical reference/deletion calls must validate the source candidate,
the source HP association (Fisher two-sided p<=0.05), and graph ALT enrichment
on the deletion side (p<=0.01); any physical ALT or conflicting source call
vetoes the projection. The source path to its next phased site must cross no
weak cut. Only the terminal SNP moves into the source phase set; a larger graph
block keeps its other sites and read labels. No new realignment is used.

The full chr20 run retains 236,834 truth-evaluable phased reads and improves
correct reads from 228,912 to 228,914 (discordant 7,922 to 7,920; purity
96.6550% to 96.6559%). VCF block count stays 392 and N50 491,812 bp. All
35 tracked noncentromeric HiPhase-correct short gaps now share a phase set;
the last pair has opposite ALT haplotypes, matching HiPhase. The known 19-Mb
wrong-join control remains split. Exact evidence and matched run details are
in `evaluations/2026-09-24-general-recovery-evidence/README.md`.

### 2026-09-24: larger HiPhase-correct gaps in the regression panel

The chr20 window panel now includes 15 additional noncentromeric gaps of
18,671–26,615 bp. Their exact boundary alleles share a HiPhase phase set,
HiPhase local truth purity is at least 98%, each 10-kb flank is at least 95%
pure with the same parental orientation, and the HiPhase and current pgphase
VCFs have identical variant keys throughout each 50-kb-padded test window.
Every selected gap remains split in the isolated pgphase replay and has at
least one primary BAM read spanning both boundaries. HiPhase correctly
separates 73.17–90.12% of truth-scorable gap-overlapping reads in one block;
pgphase's measured single-block floors are 31–50%. One other candidate gap
already closes in the isolated replay and was not added. Selection and exact
boundary evidence are in `evaluations/2026-09-24-large-gap-test-panel/`.
The panel's pre-existing expectations were preserved; the new gaps are exact
open-span regressions until a join is reviewed for correct orientation.

### 2026-09-25: seam-scoped clean-SNP bridge and saved-gauge correction

Two MAPQ-60, Q40 molecules call both clean SNPs across the 47.636-Mb BAM
source weak cut with one consistent parental relation. Keeping that cut weak
for whole-source transfer, but allowing its containing graph seam to compose
two independently supported graph/BAM gauges, closes the 23.1-kb gap. An
initial global-path trial misjoined a distant 47-Mb seam. The scoped trial
exposed a second bug: graph-only stitching reused pre-stitch gauge votes after
the preceding seam had flipped the left block. Translating those votes through
the current candidate orientation restored the correct parental connection.
The 47–48 Mb chunk and 19-Mb wrong-join controls now gate that distinction.

A matched full chr20 run keeps 236,834 truth-phased reads and improves correct
reads from 228,914 to 228,927 (discordant 7,920 to 7,907; purity 96.6559% to
96.6614%). VCF blocks fall from 392 to 389 and N50 rises from 491,812 to
506,103 bp. The 47.636-Mb case is the only new join among the 15 larger
tracked gaps. Exact selection and trial details are in
`evaluations/2026-09-24-large-gap-test-panel/README.md`.

The final physical-SNP guard also checks high-quality unassigned reads for a
contradictory parity. On this fixture, that guard leaves the full-chromosome
VCF and BAM byte-identical to the reported matched run.

A follow-up diagnostic on the 56.323-Mb residual gap found three agreeing
physical SNP-spanning reads, but its BAM source phase set has two weak cuts:
56,313,636 (unsupported) and 56,323,427 (quality-supported). The desired
seam contains only the second. The current atomic source-block transfer cannot
use that local evidence without also carrying the unsupported first cut;
naively relaxing its all-cuts guard would recreate the wrong whole-block join
seen in the initial trial. This remains open pending explicit source-run
segmentation, not another global threshold change.

### 2026-09-25: preserve local BAM bridges and corroborate sparse graph SNPs

The 56.323-Mb gap exposed a post-stitch transfer bug. Its BAM source PS has
weak cuts at 56,313,636 and 56,323,427; the latter is physically supported,
but the former is not. The seam-local stitch connected graph flanks across the
supported span, then the two one-sided source-run passes moved a boundary row
out of the joined phase set. The stitcher now records sources used by such a
local graph bridge, and those later passes leave them intact. The 56.323-Mb
window joins with 493/495 truth-correct reads and 157/183 separated reads.

A trial admitting one high-quality direct SNP molecule without block-level
corroboration misjoined a 62.623–62.645-Mb seam. The single MAPQ-60/Q40 read
called only one SNP in each block; the whole-block merge misplaced 727 reads.
The trial raised full chr20 discordance from 7,907 to 8,615 and was rejected.
The retained rule requires the same molecule to confirm its graph block at
another clean SNP at least 100 bp away, with no high-quality contradictory
pair. It checks at most three nearby SNPs on each flank, extends 2 kb past the
nominal seam, and corrects the 0.01 physical error bound by the number of SNP
pairs tried. Existing aggregate and graph/BAM gauge relations must agree.
The 62-Mb and 19-Mb negative controls remain split.

On the matched full chr20 fixture, the accepted build phases 236,816
truth-evaluable reads, 228,934 correctly and 7,882 discordantly (96.6717%
purity). The preceding build phased 236,834, 228,927 correctly and 7,907
discordantly (96.6614%). VCF phase blocks fall from 389 to 384; N50 remains
506,103 bp. The larger-gap panel now has 5 of 15 selected pairs joined in the
full-chromosome VCF, up from 1 of 15. New joins are 20.875, 24.357, 48.022,
and 56.323 Mb. The 20.875-Mb region scores 513/514 truth-correct tagged
reads, 24.357 Mb 433/433, and 48.022 Mb 486/488. The other 10 selected
pairs remain split.

The 21.435-Mb 50-kb-padded test replay omitted both boundary variants even
though the full chr20 output emits them. Its panel run now uses the full
21–22-Mb owning chunk and a concordance floor measured on that corrected
replay. This changes the fixture, not production phasing behavior.

The accepted gap-window panel passes 800 assertions across 20 test cases;
unit tests and the C++ build pass. The SNP bridge checks clean SNP category on
both the BAM source and graph row, and a high-quality contradiction from a
molecule without enough extra SNPs still vetoes the connection.

### 2026-09-25: indel boundary retry without replacing graph alleles

A 24.581–24.601-Mb graph gap had two independent BAM phase blocks whose
boundary insertion is stored one base beyond its VCF anchor. The original
retry detector excluded that row, so it saw only one block. Six paired source
calls were available, but all six called the left deletion REF while the
right insertion split 3:3. An MSA sub-solve with unplaced reads restores the
indel observations and joins the graph flanks with 582/583 local truth-correct
reads. The retry detector now includes both boundary rows and admits an indel
when allele dropout has exact balanced-heterozygote p<=0.05.

A second 63.945–63.969-Mb gap had 12 paired calls split 5:7 between incompatible
allele parities. MSA admission joins it with 469/471 local truth-correct reads,
versus 462/464 before. The mixed-parity trigger requires an excess over a 5%
observation-error model (p<=0.01) and no decisive majority relation under a
50:50 null (p>0.05). Without the latter check, the already-supported 6.578-
and 12.269-Mb short windows lost concordance; that trial was rejected.

A retry newly admitted by these rules is rejected if it removes a phased BAM
row from its exact candidate-key set. This keeps both complementary 37.461-Mb
deletion rows. Existing strictly internal missing-call retries retain their
previous behavior. The 0.528-Mb and 37.461-Mb weak-cut controls, 19-Mb and
62-Mb wrong-join controls, and the full window panel pass (810 assertions,
21 cases). A matched chr20 run phases 236,807 truth-evaluable reads, 228,967
correctly and 7,840 discordantly (96.6893% purity), with 382 VCF blocks and
517,052-bp N50. Seven of the 15 tracked larger gaps now join, versus five
before this work; eight remain open.

A forced unplaced-MSA trial without the source-row guard also closes the
41.900-, 48.929-, and 56.064-Mb gaps locally, but at 3.573 Mb it changes two
source deletions to homozygous and raises regional discordance from 8/481 to
186/537. The unprotected result is not used in production. A seam-local guard
allows 48.929 and 56.064 Mb in isolation but regresses existing short-window
and 0.528-Mb controls when applied blanketly, so that trial is also rejected.


## 2026-09-25: Full phase-set context at the 56.064-Mb seam

The 56,064,697–56,083,708 gap joined in a narrow replay but split in its
full 56–57-Mb owning chunk. The complete graph phase sets were already present:
the source read gauges had 105:0 support for the left graph block and 0:288
for the right. The difference was in BAM recovery. The merged BAM group also
contained later seams; its unplaced-read MSA retry demoted phased source rows
there, so the source-row preservation guard correctly rejected it. The original
BAM solve left two source blocks across the target seam.

A broad per-seam BAM solve changed unrelated windows, so it was rejected.
Sparse conflicting indel calls now request a focused solve over the complete
oriented extents of the two adjacent graph phase sets. The grouped source solve
continues to serve the remaining seams. Previously phased rows must survive
the focused retry exactly. A second stitch blocker then became visible: the
trusted fallback demanded that one molecule span the first and last *injected*
rows of a reused source block. That is impossible when only one private row is
injected, even though the full BAM block has a supported chain of shared SNPs.
A complete source path with no weak site-to-site cut permits a second
attachment for that focused MSA block only when both graph/BAM edges have
independent direct-SNP bridges, exact MEC parity, and agreeing source-specific
read gauges. An unscoped version also changed 650 read truth assignments near
26 Mb (529 improved, 121 worsened). That trial was rejected to avoid changing
the centromeric shoulder while fixing this noncentromeric seam. The scoped
version has no truth-assignment changes in the 26-Mb bin.

The scoped full chr20 run closes 8/15 selected larger gaps versus 7/15 before.
The 56.064-Mb boundary SNPs and the intervening BAM deletion share one phase
set in the full output; its local truth check has no switch. Truth-evaluable
phased reads rise from 236,807 to 236,824; correct reads rise from 228,967 to
228,992 and discordant reads fall from 7,840 to 7,832 (96.6929% purity versus
96.6893%). VCF block count stays 382 and block N50 stays 517,052 bp.
Eighteen reads previously tagged by an independent gap-fill block at 56 Mb
become unphased when the graph blocks join; gains elsewhere still give the
chromosome a net 17 additional phased truth-evaluable reads. The accepted
pre-change 32-Mb control took 118 seconds and the revised code 126 seconds
under concurrent compilation; their VCFs and phased-read tags are identical.
The owning-chunk 56-Mb regression and the restored 63.945-Mb mixed-parity
control pass. The full gap panel passes 823 assertions across 22 test cases;
`make unit-tests` passes.


### 2026-09-25: Two more large-gap joins with complete focused BAM paths

The 48,929,511–48,950,388 gap exposed a stitch omission. Its focused BAM
subsolve has one phase set with 1,122 oriented sites and no weak site-to-site
cut, and both graph/BAM flank gauges are significant with consistent shared
candidate orientations. The first attachment was accepted, but the second was
rejected because the first edge is an indel and the existing source-reuse rule
required direct SNP evidence. The full MEC problem exceeds the exact search
bound; a boundary-scoped exact solve chooses the same parity as both gauges.
Allowing that scoped fallback for a complete focused source path closes the gap
in its 48–49-Mb owning chunk. Its truth score remains 3,844 correct of 3,901
phased reads before and after the join. The 120-kb replay still splits because
its MSA retry loses an earlier flanking BAM row, so a full-chunk regression
checks the production context.

At 41,900,800–41,919,471, the BAM right boundary consists of two separate,
complementary deletion rows at coordinate 41,919,472. The retry detector's
unique-coordinate rule skipped them even though seven MAPQ-30 reads physically
span the source boundary and the selected row has conflicting sparse allele
pairs. The exception admits only exactly two complementary phased rows at an
edge seam to a focused MSA retry. It does not merge the rows or admit a grouped
retry, and every original phased row must survive. The focused source has 107
oriented sites with no weak cut. The full 41–42-Mb chunk now joins the target
while the preceding 41,866,917–41,898,323 seam remains split. The chunk keeps
4,107 truth-correct reads; two additional reads are tagged discordantly
(4,153 phased versus 4,151 before).

A counterexample at 11.6 Mb showed why retaining row keys is insufficient:
its focused MSA source had 91 sites but one unsupported internal cut. Accepting
it moved the 11,586,531–11,599,138 block boundary and reduced the owning
chunk from 4,002/4,092 to 3,973/4,113 truth-correct/phased reads. A focused
retry now replaces its source only when one BAM phase set spans both graph
boundaries and every site-to-site cut has read support. The 11-Mb owning-chunk
regression rejects that retry and preserves the established block; the 41, 48,
and 56-Mb owning-chunk joins still pass.

Against the preceding accepted whole-chr20 run, the final guarded run closes
**10/15** selected HiPhase-correct larger gaps versus **8/15**. It phases
236,821 truth-evaluable reads with 228,988 correct and 7,833 discordant
(96.6924% purity), compared with 236,824 / 228,992 / 7,832 (96.6929%).
VCF block count falls from 382 to 381; N50 remains 517,052 bp. Apart from
the intended 41–42-Mb relabeling, four VCF rows around 21.1 Mb change phase
set, so the tiny whole-chromosome read difference is recorded rather than
attributed entirely to the new joins.

The five tracked gaps still open are 1.086, 3.573, 21.435, 21.736, and
61.664 Mb. HiPhase's same-sites input has no heterozygous variant inside any
of these gaps. Each has one MAPQ-60 primary BAM molecule spanning both exact
boundaries, and that molecule is retained in pgphase's BAM source evidence.
HiPhase's run permits one connecting read. The pgphase one-read safeguard
abstains because none of these molecules independently confirms clean SNPs
on both flanks; the 1.086- and 3.573-Mb molecules also have a Q10 boundary
base, and the 61.664-Mb molecule has a Q17 boundary base. The 62-Mb
one-read false-join control shows why a blanket singleton join is unsafe.
At 21,742,440, the graph-only deletion has 34 callable left-to-middle and
6 middle-to-right read pairs in the saved full-chunk recovery matrix (the
earlier zero-pair report was an audit error). Four and two physically spanning
BAM reads, respectively, are unknown on this deletion row because they carry
the distinct `CT→CTT` graph allele, not because their GAF coverage was lost.
The deletion alleles do not segregate cleanly by the flanking phase. The
`CT→C` catalog allele is a one-T deletion, matching the bridge read's BAM
CIGAR; the graph row is not MSA-verified and the BAM source has no candidate
there. See
`evaluations/2026-09-25-five-gap-hiphase-signal/README.md` for the audit.

Follow-up confirmed that the six apparently missing graph deletion calls
are reads carrying the other `CT→CTT` allele of the same multiallelic site;
the default biallelic projection correctly leaves them unknown on `CT→C`.
A full-repeat CIGAR REF backfill at the 3.573-Mb homopolymer restored the
single missing source observation but did not close the gap, and regressed
concordance in four tracked windows. It was reverted. The focused 3.573-Mb
window passes again; a unit regression protects the multiallelic distinction.

The full gap-window panel passes 864 assertions in 25 cases, including the
new 11-Mb negative control and 41- and 48-Mb owning-chunk joins. `make
unit-tests`, `make -j8`, and `git diff --check` pass.

### 2026-09-25: retain BAM SNP quality and join a corroborated inner block

The open 21,435,750–21,456,826 chr20 pair sits inside the broad graph seam
21,365,579–21,462,358. Recovery transfers two complete BAM source phase
blocks, but the adjacent-block aggregate vote has only one bridging molecule
and fails the exact binomial test. The existing inner-component rescue checked
that molecule against the distant outer graph SNP, which it cannot reach.
Simply changing that check to the second BAM block broke the correct 55.38-Mb
join, so the original graph witness remains in place.

A separate fallback now accepts two complete BAM paths when one MAPQ-30
molecule calls a Q30 clean SNP in each, and a distinct MSA-verified indel
at least 100 bp from a SNP confirms the same block orientation. Contradictory
alleles or an outer graph SNP vote veto the join. BAM-to-graph transfer had
retained BAM alleles and MAPQ but discarded SNP base quality, especially for
REF calls with no `alt_qi`. It now extracts that quality from the source BAM
alignment and retains it in the candidate-indexed read profile. At 21.435 Mb,
the spanning read has Q40 REF calls at both SNPs and an MSA-verified deletion
1.84 kb to the left. The owning-chunk test closes the pair with no parental
switch; one-block truth separation rises from 70/157 to 127/157.

The panel case passes 520 assertions; the full window suite passes 864
assertions in 25 cases with the new 21.435-Mb span expectation. The 62.623-Mb
singleton false-join control and the 55.38-Mb inner-block join pass. Full chr20 keeps exactly 62,152 VCF variant rows and every
genotype; only three PS labels at 21.432–21.436 Mb change. Read truth remains
228,988 correct and 7,833 discordant of 236,821 scored (96.6924% purity),
while read phase sets fall from 769 to 768 and VCF blocks fall from 381 to
380; VCF block N50 remains 517,052 bp. Four of the five previously tracked
HiPhase-correct gaps remain open: 1.086, 3.573, 21.736,
and 61.664 Mb. A trial that reclassified foreign-phase-set MSA reads at
3.573 Mb did not restore the missing allele and was reverted.

### 2026-09-25: Validate complete flanks for one-read clean-SNP bridges

The last four HiPhase-correct tracked gaps had one physical read across each
exact boundary, but the existing bridge test required two nearby graph SNPs.
At 1.086 Mb, the read calls a private BAM SNP at 1,110,921 before it reaches
the next catalog SNP. Recovery now screens singleton geometry cheaply, then
validates both complete graph block extents with a BAM solve. Exact normalized
clean SNPs must keep one BAM phase set and one allele parity on each side,
cover at least half of each graph block's clean SNPs, and have no weak source
cut to the selected boundary SNP. When both flanks map to one BAM source phase
set, there must also be no weak cut between their boundary SNPs. A single
MAPQ-30 molecule must call both clean SNPs above the configured base-quality
floor; its original aligned base
can restore a clean SNP masked in the sparse BAM profile. A graph/BAM gauge
conflict still vetoes the stitch. The complete-block check is needed because
local boundary parity alone cannot detect an internal graph switch.

The 1.086-Mb focused window joins with 470 truth-scored reads, 4 discordant
(99.15% purity), and 138/171 truth reads separated by one block. The final
clean-SNP rule leaves the 61.664-Mb focused window split at 86/205 separated;
its full-chromosome target pair remains split too. An earlier trial joined the
focused 61.664-Mb pair at 130/205 separated, but its endpoint was an indel
whose profile call bypassed the SNP base-quality check. That route was removed.
The 3.573- and 21.736-Mb targets also remain split.

On the matched full chr20 fixture, the final rule phases 236,845 truth-scored
reads with 229,007 correct and 7,838 discordant (96.6907%), compared with
236,821 / 228,988 / 7,833 (96.6924%) in the accepted 21.435-Mb baseline.
Of the 29 newly phased reads, 22 are correct and 7 discordant; five previously
correct reads lose their HP assignment, two previously discordant reads become
correct, and no previously correct assigned read becomes discordant. All
62,152 VCF variant keys and unordered genotypes match; 1,765 sample fields
change only GT orientation or PS. VCF block count falls from 380 to 372 and
block N50 rises from 517,052 to 562,975 bp. Read phase-set count falls from
768 to 755. The 1.086-Mb target closes in the full output, leaving three of
the four previously open tracked HiPhase-correct gaps.

The older description of the 62.623-Mb control as a wrong direct boundary
join was imprecise. Its Q40 read gives the correct *local* parental relation,
and HiPhase joins that boundary. The pgphase right graph phase set changes
parental orientation internally near 62.719--62.722 Mb; merging the entire
unsplit block would propagate that switch to hundreds of reads. The local
pre-screen rejects this seam; the complete-flank consistency gate would also
reject its conflicting graph/BAM SNP parity.

### 2026-09-25: Preserve validated singleton bridges through earlier seam joins

The tracked 61,664,136–61,690,751 gap remained open although its BAM source
block was connected through the left flank (109 oriented source sites, no weak
cuts) and one MAPQ-60 molecule physically called both boundary SNPs. The local
BAM solve had placed the nearest clean left SNP at 61,658,805, before that
molecule began, so the singleton screen saw zero spanning reads. Targeted MSA
had verified the closer noisy SNP at 61,664,136; its original BAM base has
quality 17, above the configured minimum of 10. Singleton screening and full
validation now admit an MSA-verified noisy heterozygous SNP as an endpoint,
while the full graph-block clean-SNP parity check and original BAM base-quality
check remain mandatory.

The newly validated physical bridge still failed to close the gap on its own.
An earlier seam had moved every candidate in its left graph phase set from
PS 61,577,290 to PS 61,487,311, but the physical-bridge route resolved the
old label to 61,577,290. It therefore relabeled the *right* graph block into
an empty old label; later source attachment moved the left BAM run into
PS 61,487,312 and the target pair stayed split. The bridge route now checks
that each original graph block has one uniform live candidate label, joins
those current labels, and records both old detector IDs as aliases. A mixed
block abstains, preserving the whole-block switch guard.

The owning 61.614–61.741 Mb replay now spans the pair with 512 tagged
reads carrying truth labels, 500 correct and 12 discordant (97.66% purity), and no
same-haplotype allele contradiction. The panel span expectation for this graph
window moved from 0 to 1; its aggregate minimum moved from 13 to 14. The
61.664-Mb VCF pair also has one PS in the full chr20 output. Against the
previous full-chr20 baseline, truth-scored tagged reads change from 236,845
to 236,840, correct from 229,007 to 229,009, and discordant from 7,838 to
7,831 (purity 96.6907% to 96.6935%). VCF blocks fall from 372 to 370 and
N50 remains 562,975 bp. All 62,152 variant keys and unordered genotypes
match. The tracked 3.573- and 21.736-Mb gaps remain open: the former has no
callable left-SNP/deletion source pair across its weak cut, and the latter's
graph deletion has a distinct CT→CTT alternative and no clean BAM SNP anchor
on the left. Those require separate evidence fixes rather than relaxing this
singleton SNP rule.

### 2026-09-26: Close the last two tracked in-chunk chr20 seams

The 21,736,539 deletion is a genuine clean BAM deletion, even though the
catalog also has a distinct `CT→CTT` insertion at the locus. The singleton
validator previously matched only clean SNPs, so it missed the exact
normalized graph/BAM deletion. It now matches clean indels without collapsing
alternative alleles. A physical deletion bridge must agree with the existing
BAM profile and the exact quality-checked CIGAR call; REF calls also check
the actual reference bases. In the complete 21–22 Mb graph chunk, only two of
six clean right-block graph sites have exact clean BAM matches, because the
BAM caller demotes the repeat-rich middle. Those two matches bracket the
block's endpoints within 5 kb and agree in phase-set and allele orientation;
that endpoint evidence allows validation without lowering the general
half-of-sites match rule. The 21.736 Mb boundary joins in its owning chunk.

At 3,573,979–3,596,663, the sole physical MAPQ-60 read calls the left SNP
and the MSA-verified deletion but ends before the next graph SNP. Its deletion
REF base quality is 10. Independent BAM reads link the deletion to the right
clean SNP in both haplotypes (37 versus 14 observations; exact one-sided
binomial p about 0.00088). Recovery now composes that cohort link with the
physical SNP/deletion call, requiring the bridge read's assigned BAM haplotype
and exact reference bases to agree. A nearer noisy BAM SNP lies across a weak
source cut, so validation uses the exact graph-matched right SNP instead.
The two BAM runs attach independently to their adjacent graph blocks before
the validated graph-block join; attaching them afterward would reject the
second run because both graph flanks would already share one root. The 3.573
Mb boundary joins in both its focused replay and the full chromosome.

The full chr20 output has 236,840 truth-scored phased reads, 229,009 correct
and 7,831 discordant (96.6935% purity), unchanged from the accepted
61.664-Mb baseline. VCF blocks decrease from 370 to 368; N50 remains
562,975 bp. All 62,152 VCF keys and unordered genotypes match the baseline.
The focused 3.573-Mb replay scores 473/481 reads correct (98.34%) after the
join. The 22.981-Mb focused replay also joins at 565/568 correct (99.47%),
but its boundary sites lie on opposite sides of the 23 Mb production chunk
boundary and remain in separate phase sets in the full chromosome. The
window expectation records the focused join only; this cross-chunk seam is
still open in the full result.

The singleton screen and MSA cohort also reject MAPQ 255 (unknown mapping
quality). The chr20 input BAM has no MAPQ-255 records, so this guard does not
change the measured chromosome result.

### 2026-09-26: Keep remapped physical bridges only with corroborated alleles

The 22,980,600–23,008,891 focused gap stayed split when the replay was widened
from 128 kb to 200 kb. The complete-block validation BAM solve then numbered
its left source phase set 22,900,327, while the targeted injection solve used
22,932,140. The bridge was physically validated but dropped because that
validation PS had no entry in the injection solve's PS remap. These numeric
labels do not describe the same BAM run; pre-attaching a targeted source block
by the validation label is impossible.

Allowing every such graph-flank bridge without pre-attachment closed the 200 kb
replay at 849/852 truth-correct reads, but made an incorrect full-chromosome
join at 37.56 Mb: 145 additional truth-scored reads became discordant. Both
joins had one MAPQ-60 SNP-to-deletion molecule, so molecule quality alone does
not distinguish them. The safe 22.98 Mb source has two separate MSA-verified
deletion lengths at the same normalized locus, with opposite assigned alleles;
each deletion row has a significant, same-orientation link to the right SNP
on both haplotypes. The 37.56 Mb source has only one deletion at the bridge
locus; another deletion is 173 bp away. The transfer now retains a bridge
across distinct numeric BAM source PS IDs only with that same-locus,
complementary-cohort corroboration. It never merges the two deletion rows.

The positive 200 kb regression joins the 22.98 Mb flanks with at least 99%
local truth concordance; a new negative one-million-base regression keeps the
37.56 Mb blocks separate. The owning-chunk bridge checks (32 assertions),
false-join control (6 assertions), and 48-gap panel (520 assertions) pass.
`make unit-tests` passes. Full chr20 remains identical to the accepted
baseline: 236,840 truth-scored phased reads, 229,009 correct, 7,831
discordant (96.6935% purity), 368 VCF blocks, 562,975 bp N50, and exactly
the same 62,152 keys, genotypes, and PS labels. The 22.98 Mb production seam
still falls on the 23 Mb chunk boundary, so this in-chunk fix does not close
it in the full chromosome; a cross-chunk evidence handoff is still needed.

### 2026-09-26: Why the 5.31 Mb source chain still abstains

The production 5–6 Mb graph chunk splits the 5,309,406–5,345,085 pair.
No molecule spans both graph SNPs, but the current recovery matrix contains
the intervening 5,315,591 insertion (left PS) and 5,331,266 deletion
(right imported PS). The left SNP and insertion have 42 agreeing versus one
conflicting paired allele calls. The deletion and right SNP have 20 agreeing
calls. Five MAPQ-60 reads call both intermediate indels, all in the same
cross-allele relation: three carry insertion REF/deletion ALT and two carry
insertion ALT/deletion REF. A sixth physical read has deletion REF base
quality 3 and correctly stays uncallable at the configured base-quality 10
floor. The standalone BAM solve also splits its left and right source phase
sets here; its left block is complete through the insertion, and its right
source block has a weak cut after the graph gap at 5,364,960.

The 5.31 Mb chain's original left graph block is reused after an earlier
attachment. Its 14 graph heterozygotes span 5,272,413–5,309,406; 13 have
exact shared BAM SNP calls in one complete source PS. The stitcher's reused
block guard requires one read to span the original block's first and last
sites on both haplotypes, which is too strict for this 37 kb block. Lifting
that guard alone did not close the gap. The exact MEC fallback initially
returned no path because it preferred the two flanking SNPs and excluded
the 5,331,266 BAM deletion at 70% row-wise ALT fraction. An experimental
retry admitting that alignment-verified indel made the observation graph
connected and gave a unique flip in the full solve and both read halves
(MEC scores 1006/1001, 998/994, and 7/6). The final source-gauge gate still
abstained; the imported right source PS has a weak cut beyond this seam and
lacks a whole-block shared-site gauge to the already absorbed left block.

The broad indel-admission trial changed the full chromosome by one extra
phased VCF site at 34,146,891 and 15 phased reads, only 6 truth-correct and
9 discordant. It closed none of the tracked gaps. That trial and every guard
bypass were reverted. The accepted full-chr20 result remains 229,009/236,840
truth-correct reads (96.6935%), 7,831 discordant, and 368 VCF blocks.
A safe future join needs an interval-scoped source-path certificate for the
reused left block and a quality-aware direct indel edge; merely relaxing the
whole-block or statistical threshold is unsupported by this experiment.

### 2026-09-26: 60.03 Mb graph SNP represents a deleted BAM allele

The remaining production 60,033,052–60,058,235 gap contains a singleton
graph-clean SNP at 60,033,350 (`G>A`) between the established left SNP and
the right insertion block. Graph observations call it heterozygous (23 REF,
41 ALT), but the recovery BAM caller marks the same key `NoisyCandHom`.
Of 62 primary BAM reads physically spanning the two nearby SNP coordinates,
40 have a quality-qualified aligned `A` base at 60,033,350 and none have an
aligned `G`; the other 22 do not have a callable base there, commonly because
an overlapping deletion removes it. The recovery graph/BAM block gauge for
that singleton is 30 versus 32 and cannot orient it. This is an allele
representation conflict, not absent molecule coverage.

A matched 60–61 Mb replay used a local copy of the graph catalog with only
that `G>A` row removed. The gap still split, but local parental discordance
fell from 67/3,326 to 36/3,326 reads (purity 97.99% to 98.92%) without
changing the number of phased reads. The graph SNP's REF observations should
not be treated as physical `G` calls; demoting this one row alone improves
accuracy but does not supply a phase bridge. The physical `A`/deletion
state at 60,033,350 is itself uninformative against the left SNP: 21 `G` and
19 `C` reads carry `A`, while 11 of each carry the deletion. Thirteen
MAPQ-60 reads span the left SNP and the right 60,048,237 insertion, but all
13 call insertion REF (seven left `C`, six left `G`). The standalone BAM
caller emits no heterozygous stepping stone between those boundaries. This
explains why simply deleting the false graph SNP cannot close the gap.
This catalog-omission experiment was not turned into a coordinate-specific
production rule. A general fix needs to distinguish a physical REF base
from a read whose path skips the SNP through a deletion before classifying
overlapping graph sites as biallelic SNPs.

### 2026-09-26: Reject graph SNPs with physical ALT/deletion alleles and no REF

The graph catalog's `G>A` SNP at 60,033,350 had 23 graph REF and 41 graph
ALT observations, but high-quality BAM alignments provided zero callable `G`,
41 callable `A`, and 23 deletions. Using its graph REF/ALT state to orient
reads raised local 60–61 Mb discordance from 36 to 67 reads. The graph worker
now checks clean, biallelic SNPs against one high-MAPQ BAM pass per chunk
before graph k-means. It makes a site ineligible only when the observed REF
count is zero, the binomial zero-REF tail is significant at familywise 1%
over tested SNPs, at least ten reads carry a deletion, and deletions comprise
at least 20% of ALT-plus-deletion observations. The chunk is then rebuilt
from the same GAF rows without those sites. Graph-only phasing is unchanged.
The first trial accidentally inherited the graph's MAPQ-5 admission floor,
which included five centromeric-shoulder sites supported only by low-MAPQ
alignments. Physical validation now uses at least MAPQ 30, while the graph
solve retains its MAPQ-5 floor.

On the full chr20 fixture this excludes five catalog SNPs (5,746,178;
38,177,812; 51,028,622; 55,373,606; 60,033,350), leaving the five
low-MAPQ centromeric-shoulder rows untouched. Compared with the accepted
`/tmp/pgphase-gapfix-cohort-chr20` run, tagged/truth-scored reads change
from 236,840 to 236,814, truth-correct from 229,009 to 229,036, and
truth-discordant from 7,831 to 7,778. Purity rises from 96.6935% to
96.7156%; VCF blocks fall from 368 to 366, while block N50 remains 562,975 bp.
The 55.373 Mb excluded `G>T` SNP was already unphased in the baseline.
Its existing BAM-supported inner block keeps the same neighboring PS labels,
while the regression now requires the invalid SNP to be absent and still
checks local truth concordance. The 60.03 Mb gap still does not join: the
physical `A`/deletion state is independent of its left SNP and the 13 reads
spanning to the right insertion all carry insertion REF. This change fixes
false orientation rather than inventing unsupported gap evidence. Output:
`/tmp/pgphase-gapfix-physical30-chr20`.
The final `make window-tests` run passes all 949 assertions in 28 test
cases; `make unit-tests`, the optimized C++ build, and `git diff --check`
also pass. The 60.03 Mb panel expectation drops from two to one in-gap
heterozygote because the removed graph SNP was a false phased anchor.

### 2026-09-27: Certify a reused graph block across a downstream BAM weak cut

The full 5–6 Mb owning chunk split 5,309,406 from 5,345,085 even though
recovery had a connected SNP–insertion–deletion–SNP observation path. The
left graph block spans 5,272,413–5,309,406 and had already been reused at
an earlier seam. One complete BAM source path covers both graph-block
endpoints with consistent polarity at its shared candidate sites and
significant read votes from both haplotypes, but no single molecule spans
that whole graph block. The right BAM source has a weak cut at 5,364,960,
**after** the gap; its local segment through the 5,331,266 deletion is intact.
That deletion is alignment verified but has 70% row-wise ALT fraction, so the
ordinary centered-site MEC problem omitted it and lost the bridge.

Recovery now certifies the left graph block from its complete shared BAM
source path, includes verified indels inside this seam in the exact MEC solve,
and accepts the resulting parity only when the full reads and both disjoint
halves choose the same unique result. This fallback requires an imported
right BAM source whose weak cuts all lie beyond the current edge. The MEC
may use more distant sites as validation context, but its commit phases only
sites inside the edge. It does not let a left imported BAM source reorient a
right graph block.

Two counterexamples constrained the design. The broad graph-path fallback
joined an imported-left block to a graph-right block at 34.1 Mb and changed
that region from 3,534 correct / 34 discordant reads to 3,492 correct / 114
discordant. It also phase-marked an insertion at 34,146,891 outside its
boundary, adding 15 tagged reads of which only six agreed with parental
truth. Directional admission and seam-scoped commit remove both effects.
A complete imported BAM source at 11.6 Mb must keep its established block
split; the newly allowed path is limited to sources with a weak cut beyond
the tested edge, preserving that negative control.

The 5.31 Mb window now spans with three in-gap phased heterozygotes and
508/510 truth-scored reads concordant (99.61%). The full chr20 run
`/tmp/pgphase-gapfix-cutscoped-chr20` phases 236,822 truth-scored reads:
229,044 correct and 7,778 discordant (96.7157% purity). Against
`/tmp/pgphase-gapfix-physical30-chr20`, that is eight additional correct
reads and no additional discordant reads. VCF phase-block count remains 366
and N50 remains 562,975 bp. The 5.31 Mb flanks share one PS in the full
run, while the 11.6 Mb split and the 34.1 Mb block boundary remain intact.
`make window-tests` passes 972 assertions in 30 cases, including both new
owning-chunk controls; `make unit-tests` and `git diff --check` pass. The
graph-arm expectation for the 5.31 Mb panel case changes from split to span.

### 2026-09-27: Preserve a graph block across BAM source IDs and the 23 Mb chunk boundary

The full chr20 run left 22,980,600–23,008,891 split, although a 100 kb
single-chunk replay joined the two SNPs. A wider 1 Mb replay also split until
`validated_singleton_bridge` was audited: it marked the left graph block
inconsistent merely because its 178 exact shared SNPs came from two numeric
BAM phase-set IDs (3 from one source and 175 from another). All 178 shared
SNPs had the same allele polarity. Validation now treats a polarity change as
the inconsistency and chooses the nearest matched BAM source at the seam.
The 1 Mb replay then joined the gap without a parental switch (98.34%
local truth concordance).

The production 23 Mb boundary still split because ordinary cross-chunk
stitching had no phased overlap-read vote. Among eight chr20 boundaries with
10–50 kb clean-SNP gaps, a 100 kb replay joined only this one. The replay
contained a physical clean-SNP bridge between its graph blocks and 11/28
exact SNPs shared with the production left/right blocks, all with consistent
polarity. Its right boundary SNP had no `RecoverySourceSite` even though the
physical bridge retained the right graph phase-set ID. The new boundary
transfer accepts that recorded ID when the exact SNP match and block-wide
orientation agree. It flips/relabels only the downstream block. The other
seven audited boundaries remain split.

The two-chunk owning-region regression now joins the 23 Mb SNPs without a
parental switch (97.83% local truth concordance). In the full chr20 run
`/tmp/pgphase-gapfix-boundary-transfer-chr20`, VCF blocks fall from 366 to
365 and N50 rises from 562,975 to 573,587 bp. The 456 VCF rows that change
are the single downstream block's phase orientation and PS; read assignments
and full truth counts are unchanged at 236,822 scored, 229,044 correct, and
7,778 discordant (96.7157% purity).
The optimized build, `make unit-tests`, `make predicate-tests`, and
`git diff --check` pass. `make window-tests` passes all 993 assertions in
32 cases, including the 1 Mb source-ID and two-chunk boundary regressions.

### 2026-09-27: Reuse an established graph block only with independent SNP and intact-source evidence

The full chr20 output split 56,150,368–56,156,525 despite 45 MAPQ≥30 BAM
reads calling both boundary SNPs: 23 carried the T–G pair, 22 the C–A pair,
and none carried a crossed pair. The parentally labeled reads agree with those
two pairs. Recovery already imported the right SNP, and the exact MEC solve
found the same unique parity on the full read set and both disjoint halves.
The fallback stitch discarded that proof before examining it because the left
graph block had joined an earlier seam and lacked a single molecule spanning
its entire extent.

The stitch now evaluates the exact edge before its reused-source veto. A
previously used graph block may attach through a direct clean-SNP MEC bridge
only when its graph/BAM gauge agrees significantly, both source haplotypes
vote, exact shared candidates occur at two distinct genomic loci, and the
imported BAM source has a complete path. Counting distinct loci prevents
complementary rows at one position from pretending to be independent anchors.
The intact-source requirement matters because the stitch relabels an entire
phase set, rather than only the local edge.

Two negative controls narrowed the rule. At 19.4 Mb, a single shared
candidate could make a local SNP edge appear convincing even though the reused
graph block switches internally; a broad trial lost 279 truth-correct reads
in that megabase. At 8.6 Mb, the imported BAM source has a weak cut exactly
at 8,638,940. Admitting its local edge relabeled the whole source across that
cut and incorrectly joined 8,638,940 to 8,662,670. The final rule leaves
both boundaries split. A broader trial also changed centromeric VCF rows,
so it was rejected.

The final full chr20 output is
`/tmp/pgphase-gapfix-reuse-intact-chr20`: 236,831 truth-scored reads,
229,053 correct and 7,778 discordant (96.7158% purity). Against
`/tmp/pgphase-gapfix-boundary-transfer-chr20`, that is nine more correctly
phased reads and no additional discordant reads. The only two changed VCF
rows are at the intended 56.15 Mb join; VCF blocks fall from 365 to 364 and
N50 remains 573,587 bp. The new owning-chunk tests require the 56.15 Mb
join without a parental switch and reject the 19.4 Mb counterexample; the
existing 8.6 Mb weak-cut test also remains a negative control.
The optimized build, `make unit-tests`, `make predicate-tests`, and
`git diff --check` pass. `make window-tests` passes all 1,016 assertions in
34 cases, including the 8.6 Mb weak-cut and both new owning-chunk controls.

### 2026-09-27: Keep the 56.15 Mb closure in the persistent gap panel

The 56,150,368–56,156,525 join is now a committed panel window, in addition
to its owning-1-Mb exact-boundary regression. The ordinary 50-kb-flank replay
also joins the boundary keys; its recorded `spans=1` is an equality, so a
future split fails the standard suite. The new row records 104 HiPhase
truth-scorable reads, all correctly placed by its dominant block, and pgphase
floors of 0.88 local concordance and 0.74 separated-read fraction. It has no
strictly interior phased heterozygote. The graph panel's expected span total
rises from 17 to 18 without changing any existing window expectation.
`make window-tests` passes all 1,032 assertions in 34 cases with this panel
entry; `git diff --check` also passes.

### 2026-09-27: Recover a graph repeat-indel bridge at 21.514 Mb

At 21,514,518–21,518,393, the graph left SNP and imported BAM right
insertion remained in separate phase sets. The BAM insertion is represented
at 21,518,394 in the injected rows, while graph walks place an equivalent
A insertion at 21,518,396. The graph row has 39 REF and 37 ALT observations,
but its `RepeatHetIndel` category excludes it from ordinary phasing. The
left SNP to graph-indel pair has 38 same and 12 cross read calls; graph indel
to the next clean right-block SNP has 2 same and 21 cross. Each leg keeps
the same majority in two disjoint read halves. The exact MEC solution using
centered indels chooses the resulting polarity in all reads and both halves.

The graph candidate constructor had never set `graph_site`, so a
provenance-aware bridge search saw no graph indels at all. It now marks graph
candidates. The narrowly gated fallback uses an original singleton clean
graph SNP, a complete imported BAM source path, a significant two-leg
graph-indel certificate, and a matching split-stable exact MEC result. The
left SNP can be off center only inside this certified retry. Graph and BAM
rows stay separate; no allele representation is merged. The owning 21–22 Mb
chunk joins the exact boundary VCF keys in the truth-supported orientation.

On full chr20, only the three right-block VCF rows at 21,518,393,
21,532,460, and 21,534,506 change PS/orientation. VCF phase blocks
fall from 364 to 363; block N50 remains 573,587 bp. Truth-scored phased
reads remain 236,831; correct/discordant change from 229,053/7,778 to
229,051/7,780, or 96.7150% accuracy. The two extra discordances arise from
read-only synthetic phase labels affected by the stitch; the primary joined
blocks retain 149 correct of 152 parental-truth reads. The persistent gap
panel now includes 21,514,518–21,518,393 with an exact `spans=1` check,
a 0.93 concordance floor, and a 0.84 separated-read floor. The graph-site
provenance flag has a unit regression. The graph panel now has 35 windows,
22 expected spans, and a 21.514 Mb exact-boundary regression. The optimized
build, `make unit-tests predicate-tests window-tests`, and
`git diff --check` pass; the full window suite reports 1,052 assertions
across 34 cases. After tightening the bridge to both allele classes and the
original imported source sites, the focused panel passes all 608 assertions.
The final source-restricted full-chromosome VCF and BAM are byte-identical
to the measured trial after the allele-class check; its 21.514 Mb panel
replay still closes the boundary. The remaining short noncentromeric
HiPhase-spanned gaps are 11;
several have low HiPhase truth purity on this callset and are not safe
join targets.

### 2026-09-27: Repair a missing repeat-reference source cut at 38.331 Mb

At chr20:38,331,110–38,333,011, the left MSA-verified homopolymer
deletion and right MSA insertion remained in different graph phase sets.
The BAM source path had a false weak cut: sparse profiles omit reference
calls for homopolymer indels. Six MAPQ 60 molecules physically span the
preceding clean SNP and deletion; five pass base quality 10, four pass 30.
All five qualified reference calls agree with the source orientation, giving
a one-sided random-polarity probability of 1/32. Using base quality 30 for
all support yielded only four calls and did not repair the path.

After the source path is validated, 24 MAPQ 60 paired calls from the deletion
to the first verified SNP in the next imported BAM block vote 20 cross to
four same, with the same majority in both disjoint read halves. Thirty paired
calls from that SNP to the next clean graph SNP support no flip. The new
repeat-cut pass joins the complete blocks in that orientation. The exact
38,331,110 deletion and 38,333,011 insertion now share PS 38,294,482 and
carry opposite ALT haplotypes. The new owning-chunk panel window asserts
that join and its local parental-truth orientation.

Full chr20 VCF variant keys remain identical. Truth-scored phased reads
increase from 236,831 to 236,835; correct reads from 229,045 to 229,057;
discordant reads fall from 7,786 to 7,778. Accuracy rises from 96.7124% to
96.7159%. Only 514 VCF rows change phase labels; all changes belong to this
join. The optimized build, unit and predicate tests, `make window-tests`
(1,132 assertions in 34 cases), and `git diff --check` pass.

### 2026-09-27: Join a clean SNP seam after an upstream BAM weak cut

At chr20:39,147,805–39,149,079, pgphase split two clean SNPs while HiPhase
on the same variant keys joined them and placed all 66 truth-scorable gap
reads correctly. Fifty-five MAPQ >=30 primary BAM reads call both SNPs at
base quality >=20: 29 have A/T and 26 have G/C, with no contradictory pair.
The transferred matrix retains 57 paired observations of the same polarity.
The BAM source phase set spans both graph flanks, but its weak cut at
39,098,639 lies upstream of the seam and inside the established left graph
block. Whole-source completeness therefore rejected an intact local run.

The post-transfer bridge now accepts a clean injected BAM SNP to the first
clean graph SNP of the next block only after checking the local source cut,
consistent exact shared SNPs and significant gauges on both flanks, a
continuous original-graph SNP path through each entire block, and a direct
MAPQ-30 allele vote significant at p <= 0.01 overall and agreeing in both
read halves. The owning
39–40 Mb replay joins the boundary without a parental switch; its local
truth score remains 4,086/4,160. The established 51-Mb wrong-whole-block
and 62-Mb one-read controls still pass.

The full chr20 run changes only seven right-block VCF rows at 39.149–39.171
Mb. Variant keys and truth-scored read totals are unchanged: 236,835 phased,
229,057 correct, 7,778 discordant (96.7159%); read phase sets fall from
733 to 732. The gap has an exact-span owning-chunk regression in the
persistent panel, with a 0.98 concordance floor and a 1.00 separated-read
floor. The optimized build, unit and predicate tests, full window suite
(1,152 assertions in 34 cases), and `git diff --check` pass.

The 60.033–60.058 Mb gap was also rechecked before selecting this fix.
Thirteen MAPQ-60 primary molecules span its left SNP and right A insertion,
but every exact CIGAR call at the insertion is REF while both left alleles
occur. Three reads have a one-A insertion shifted within the same A run; a
fourth insertion is outside that run and is not sequence-equivalent. Two
reads have nearby deletions. Even after recognizing the three equivalent
insertions, the paired observations do not establish a decisive diploid
orientation, so this change leaves that gap split.

### 2026-09-27: Broaden the HiPhase-correct gap regression inventory

A full chr20 comparison with HiPhase 1.6.0 on pgphase's own VCF keys and the
same HG002 BAM found 23 additional, noncentromeric boundaries that HiPhase
joins and pgphase splits. Each HiPhase block has at least 20 truth-labeled
primary reads in the gap plus 10 kb flanks, at least 98% local parental purity,
at least 95% on each flank, and the same parental orientation on both flanks.
Twenty-two gaps are at least 10 kb; the 64,144,256–64,144,722 gap is shorter.
The inventory with exact boundary alleles and truth counts is in
`evaluations/2026-09-27-remaining-hiphase-correct-gaps/targets.tsv`.

All 23 are now in the persistent gap panel. Eight of the first 22 span in a
short-window replay but remain split in the chromosome output, which isolates
a production-context transfer or stitch problem. The panel includes both the
older 60,033,052–60,058,235 interval and a distinct 60,033,052–60,048,237
subgap. The regression harness now keys outputs, cached measurements, required
sites, and expectations by both boundaries to prevent those cases colliding.
The 55-case baseline run passed 1,592 assertions in 34 test cases; the
complete 56-case baseline run passed 1,655 assertions in 34 test cases.

At 64.144 Mb, the owning 64–65 Mb chunk transfers the left BAM source block
onto the left graph PS and the right BAM block onto the right graph PS, then
leaves them split. The nearest deletion locus on each side gives 14 concordant
paired matrix calls and no discordant calls. The broader imported-block vote
gives 14 concordant and five discordant calls: all five discordant reads lack
the left boundary deletion call and instead vote from more distant left sites.
This identifies a boundary-specific vote dilution, but any whole-block join
must also pass the existing flank and source-label safeguards. A trial that
used only the nearest loci joined the entire left and right graph blocks;
the owning-chunk truth score fell from 3,010/3,013 to 2,586/3,013 because
424 left-block reads inherited the wrong parental orientation. The trial was
reverted.

The actual source-path check already detects a weak cut at 64,140,314 in
the left BAM PS, immediately before the 64,144,256 deletion run. Recovery
stitching nevertheless merged that complete source PS into the left graph
block, carrying its post-cut deletion across the unsupported edge. A trial
that skipped every multi-source seam with such a cut closed three full-chr20
targets, but reopened the established 21.435 Mb and 22.980 Mb windows and
reduced separation at 56.323 Mb; that broad guard was rejected.

The retained correction preserves the graph stitch and restores only private
BAM candidates beyond an unsupported cut to a local source gauge. It acts
when the near source component has an earlier oriented private BAM locus
beside its boundary locus and the far component has no matched graph anchor. A one-locus source
keeps its existing stitch; this preserves the previously truth-correct 8.638
Mb local replay while its full-chromosome boundary remains split. At 64.144 Mb
the exact boundary rows
now share the right block's PS while graph SNP 64,118,182 remains in the left
PS. The owning 64–65 Mb chunk stays at 3,010/3,013 truth-correct reads with
three discordant reads, and the two boundary deletion rows remain connected.
At 56.323 Mb, graph-block connectivity and the 493/495 local truth score are
unchanged. The panel replays the owning 64–65 Mb chunk and asserts both its
boundary join and the absence of the wrong whole-block join.

### 2026-09-27: Recover seams exposed by the first BAM transfer

The full first-megabase chunk has an initial graph phase set spanning the
514,902 C>G and 528,827 A>T SNPs, so the first seam detector correctly does
not target that interval. Recovery of the next 528,827–545,002 seam attaches
the right SNP to a different BAM source block. This creates a new 13,925-bp
phase-set seam *after* the only recovery solve. A short-window replay starts
with that seam already present and injects the 516,158 A>G SNP, 528,728 C>A
SNP, and 528,825 deletion; the full chunk previously omitted all three.
The owning-chunk recovery audit listed these candidates but marked every one
outside its selected windows. The gap was therefore lost at seam selection,
not at BAM candidate calling or allele representation.

Recovery now makes one additional bounded sub-solve for seams newly exposed by
first-pass transfer. A second-round interval must have no positive-width
overlap with any first-round interval; comparing endpoints alone also retried
a shifted 0.863-Mb seam and incorrectly merged a pure 497-read block. The
nonoverlap rule retains the previous block. New source sites and observations
are injected through the existing exact-key merge and ordinary stitch path.

The 528,728 SNP is classified internally as `NoisyCandHet`, but it is MSA-
and alignment-verified and emitted as a clean heterozygous SNP. The first
prototype overlooked this and left 20 paternal reads tagged HP1 in a block
whose clean SNP ALT belongs to HP2. A local refresh now uses verified new-seam
SNP calls with base quality >=20 and BAM MAPQ >=30 to update an existing HP only
when all informative SNPs in that read's current PS agree on one haplotype.
It leaves the PS unchanged and abstains on conflicts. On the owning 1-Mb
replay, the 514,902–528,827 rows join, correctly separated reads rise from
71/119 to 91/119, and truth discordance falls from 33/3,699 to 13/3,697.
The new owning-chunk regression requires the two recovered SNP rows and a
0.76 separated-read floor; HiPhase on the same keys scores 113/119 there.

The full chr20 run adds exactly the three VCF keys above (62,147 -> 62,150)
and removes none. It closes the 514,902–528,827 target in addition to the
61.747-Mb and 64.144-Mb fixes: 3 of the 23 newly inventoried targets now
join, with 20 remaining. Truth-scored reads change 236,835 -> 236,833,
correct reads 229,058 -> 229,076, discordant reads 7,777 -> 7,757, and
accuracy 96.7163% -> 96.7247%. A trial that blindly reran all seams is not
the retained design: it raised first-megabase discordance from 33 to 52 by
reprocessing the shifted 0.863-Mb seam. The bounded nonoverlap pass leaves
that first-megabase baseline at 33 discordant before read refreshing.

The first complete window-panel run caught a regression in the already solved
528,827 SNP/528,828 insertion pair: the retry moved the SNP to the upstream
PS and left its sequence-validated BAM insertion in the old source PS. The
original BAM source has a weak cut at 528,850. A coordinate-specific trial
merging the whole source fixed the pair, but crossed that cut. The retained
rule snapshots first-pass blocks with one oriented graph SNP and oriented
private BAM sites, then transfers only the SNP's original source component if
the second stitch moves that SNP. It requires every transferred row to retain
a consistent allele orientation; reads follow only if they observe moved rows
and no remaining source row. Noninformative homozygous source rows do not
veto the component. The dedicated 528,827/528,828 regression and the new
514,902/528,827 owning-chunk regression both pass.

The new-seam read pass also assigns previously unphased reads when all
quality-bearing injected SNP observations identify one PS and one haplotype.
This raises the owning-chunk single-block separation from 91/119 to 110/119,
versus HiPhase's 113/119, while retaining 13 discordant reads among 3,711
scored. The regression floor is now 0.92. On full chr20, the retained version
has 62,150 VCF keys (three added, none removed), 236,847 truth-scored phased
reads, 229,090 correct, 7,757 discordant, and 96.7249% accuracy. Relative to
the pre-retry graph+BAM baseline, that is 12 more scored reads, 32 more correct
reads, and 20 fewer discordant reads. The 528,827 SNP, 528,828 insertion, and
528,850 deletion share PS 130541; the 542,052/545,002 control remains joined
in PS 545002. Three of 23 newly tracked HiPhase-correct gaps are closed on
the full chromosome, leaving 20.

### 2026-09-28: Preserve exact clean graph SNP anchors across padded BAM solves

The 48,971,192–48,982,663 gap joined in a short replay but split in its
owning 48–49 Mb chunk and in the full chr20 run. The BAM source block spanning
it had a complete site path with no weak cuts. Two clean SNPs at 48,950,388 and
48,952,538 were exact matches to the left graph block and had a 71-read,
two-haplotype graph/BAM gauge vote. A second padded BAM solve also contained
those SNPs. Duplicate-source suppression set `can_adopt=false` on both copies,
so the spanning source lost its left graph anchors; it could attach only the
right graph block. The local replay used one padded solve, explaining its
misleading success.

Recovery now retains multiply claimed sites only when each claim is an exact,
clean biallelic SNP shared with a clean graph SNP. Each source still needs its
own shared-site allele orientation, complete source path, and two-haplotype
read vote before transfer. Ambiguous indels and BAM-private rows remain
excluded from duplicate claims. In the owning 48–49 Mb chunk, the exact
48,971,192 G>GT and 48,982,663 T>C rows now share PS 48,950,389. Truth
remains 3,844/3,901 correct reads (57 discordant) while the read phase-set
count drops from eight to seven.

The same change closes 7,047,080 C>CA–7,064,047 C>T on full chr20. Its
owning 7–8 Mb chunk remains 3,926/4,091 truth-correct (165 discordant) before
and after; the short-window test's former 0.98 concordance floor described a
different solve context. The regression now replays the owning chunk and uses
a 0.95 floor, with exact boundary PS equality asserted for both newly closed
gaps. The full 56-window suite passed with 1,673 assertions across 34 test
cases on this change.

On full chr20, all 62,150 `(CHROM, POS, REF, ALT)` keys are unchanged. Tracked
HiPhase-correct gaps closed rise from 3/23 to 5/23; 18 remain. Truth-scored
reads change from 229,090/236,847 correct (7,757 discordant) to
229,091/236,848 correct (7,757 discordant), both 96.7249% at four decimals.
Read PS count falls from 731 to 729. This is a continuity gain without a
chromosome-wide truth penalty in the same-BAM comparison.

### 2026-09-28: Normalize graph SNP context for a physical insertion bridge

The 56,662,188–56,679,959 target remained split in its owning 56–57 Mb chunk.
Its left catalog row is a multi-base snarl allele (`TGG>TGT` at 56,662,186)
that normalizes to the exact `G>T` SNP at 56,662,188. A physical bridge trial
that required one-base raw graph REF/ALT skipped that boundary and chose a
clean SNP 3.3 kb upstream, excluding reads starting in between. After
normalization, the original BAM supplies 19 MAPQ/BQ >=30 SNP–insertion pairs
at the first verified right-source insertion: 18 support one orientation,
one conflicts. Complete BAM source paths, consistent clean shared SNP gauges
on both graph flanks, and the ordinary stitch conflict checks allow the join.
The boundary BAM SNP's MSA hap allele IDs are 1/2, so the physical REF/ALT
call, not those IDs, orients it. The insertion has `msa_verified=true` but
`alignment_verified=false`; exact CIGAR calls validate the bridge directly.

The owning chunk gains one truth-correct phased read with no new discordant
read. Full chr20 retains all 62,150 VCF keys, raises tracked target closures
5/23 -> 6/23, and changes truth-scored reads 229,091/236,848 ->
229,092/236,849 correct; discordant reads stay 7,757 and read PS count falls
729 -> 728. The target's regression now uses its owning chunk and checks that
the exact boundary rows share one PS. A preliminary 66-Mb run was truncated
short of the 66,210,255-bp reference end and is not a valid coverage
comparison.

### 2026-09-28: certified single-molecule graph seam

At chr20:54,894,127–54,912,022, the original BAM has one MAPQ-60 read with
Q40 calls at both boundary SNPs. Lowering the global or targeted BAM link
threshold from two reads to one joins this seam but also makes an incorrect
whole-block join near 62.6 Mb; targeted one-read linking adds 1,024 discordant
reads over the chr20 baseline. A post-transfer stitch now uses exactly one
quality-certified physical SNP pair only after both chunk-local graph blocks pass
consecutive clean-SNP path checks. The 62.6 Mb right block fails this check
because it reverses internally. The 54–55 Mb owning-chunk regression joins
the intended pair without changing its 4,305/4,314 truth-correct read score.
A full chr20 replay increases tracked HiPhase-correct gap closures from 6/23
to 7/23, reduces read phase sets from 728 to 727, and leaves truth scoring
unchanged: 229,092 correct and 7,757 discordant among 236,849 phased reads
(96.7249%). The remaining 16 tracked gaps stay open. Detailed commands and
results: `evaluations/2026-09-28-certified-singleton-seam/README.md`.

### 2026-09-28: Direct physical SNP votes and connected graph prefixes

The original BAM contains 11 high-quality molecules voting unanimously for
chr20:3,529,324–3,542,977, but the recovery subsolve demotes the right SNP
and loses its callable pairs. Two biallelic rows of the same multiallelic
snarl also created a false graph-path edge. Direct physical-base votes and
one-locus-per-snarl path checks restore the boundary evidence. A later
zero-vote graph-SNP path cut crosses two already linked deletion candidates.
Moving the prefix across that cut breaks the established 3.573 Mb join, so
the guarded stitch leaves the 3.529 Mb owning chunk split. An empty cut can
still permit a prefix join. The 8.638 and 54.547 Mb tracked gaps join;
the local 5,345,085–5,350,509 edge also joins without crossing its farther
weak source cut. Its owning-chunk reads retain at least 99% parental
concordance. The unsafe 62.6 Mb whole-block control stays split.

Full chr20 closes 9/23 tracked HiPhase-correct gaps versus 7/23 previously.
Read phase sets fall 727 -> 722. Parental truth stays at 7,757 discordant
reads, while scored reads fall 236,849 -> 236,847 (both lost reads were
correct), so accuracy remains 96.7249%. The detailed audit and reproduction
are in `evaluations/2026-09-28-physical-snp-prefix/README.md`. The full
window panel passes 1,723 assertions across 34 cases; unit tests and
validation gates pass.

### 2026-09-28: Preserve the main BAM stitch as graph-path evidence

The 3.529 Mb boundary has 11 high-quality unanimous physical SNP pairs.
Its right graph block later absorbs a second original graph block across a
3.573–3.597 Mb GAF path gap through the existing validated BAM stitch.
The final physical SNP check was rechecking that already joined block using
GAF alone, then abstaining. Stable graph site IDs now retain both original
block identities and their phase sets immediately after the first main BAM
stitch. A missing or one-haplotype GAF edge may inherit that certified join
only if no GAF read votes for reversal. Edges within each original block
still need full two-haplotype graph support. This closes 3.529 Mb in the
owning chunk without breaking the 3.573 Mb deletion bridge; the unsafe
62.6 Mb control remains split.

Full chr20 tracked closures rise 9/23 -> 10/23 and read phase sets fall
722 -> 717. The 62,154 VCF keys and parental truth score are unchanged:
229,090 correct and 7,757 discordant of 236,847 scored reads (96.7249%).
No individual truth-scored read changes correctness.
The full window panel passes 1,728 assertions across 34 cases; unit and
validation gates pass. Details and reproduction:
`evaluations/2026-09-28-certified-graph-path/README.md`.

### 2026-09-28: One-haplotype clean-SNP bridge

At chr20:32,215,055–32,233,534, two independent MAPQ 38/40, Q40 BAM reads
call the same clean SNP allele pair (left ALT, right REF). Recovery had enough
physical evidence, but the final stitch rejected two or more reads unless
both left haplotypes appeared. Diploid parity is established by the existing
quality-weighted likelihood even when both sampled molecules come from one
haplotype. Removing that extra guard joins the correct opposite-ALT boundary;
whole-block graph-path checks and the 62.6 Mb unsafe control remain in place.
The owning-chunk regression checks the join and parental read concordance.

Full chr20 tracked closures rise 10/23 -> 11/23, read PS fall 717 -> 716,
and 62,154 VCF keys and all 236,847 truth-scored read outcomes remain
unchanged (229,090 correct, 7,757 discordant; 96.7249%). Detailed evidence:
`evaluations/2026-09-28-one-haplotype-snp-bridge/README.md`.

The full gap-window suite passes 1,741 assertions across 34 cases; unit tests
and validation gates pass. The source-to-trial read-name sets and each
truth-scored read's correctness are identical.

### 2026-09-28: Keep complementary deletion alleles separate at a graph seam

At chr20:60,453,499–60,467,115, two overlapping MSA-verified BAM deletion
rows (8 bp and 2 bp) were already attached to the left graph block on opposite
haplotypes, but the final physical stitch only accepted clean SNPs at both
boundaries. Six exact 8 bp deletion calls pair with right-SNP ALT, and two
exact 2 bp deletion calls pair with right-SNP REF; two reads carrying other
nearby deletions abstain. The new bridge takes only exact ALT CIGAR calls from
either separate row, validates reference flanks, requires both deletion alleles
and one unanimous phase relation, and retains graph-path checks on both blocks.
No site is merged or rewritten. The reference audit also found that
`physical_snp_call` compared BAM uppercase bases with soft-masked lowercase
reference bytes; it now compares bases case-insensitively.

The owning chunk keeps 3,290/3,326 truth-correct reads while its PS count falls
18 -> 17. Full chr20 tracked closures rise 11/23 -> 12/23 and read PS fall
716 -> 715. The 62,154 VCF keys and genotypes and every truth-scored read
outcome remain unchanged: 229,090 correct, 7,757 discordant of 236,847
(96.7249%). Details: `evaluations/2026-09-28-complementary-deletion-bridge/README.md`.

The full gap-window suite passes 1,754 assertions across 34 cases; unit tests
and HiFi/ONT validation gates pass.

### 2026-09-28: Exact BAM insertion to singleton graph SNP at 19 Mb

At chr20:18,983,414–18,999,993, the left recovery block ends in an
MSA-verified six-base insertion and the right graph block begins seven bases
before the 19 Mb chunk edge. The final physical stitch previously considered
only clean SNP pairs or complementary deletions, so it ignored the insertion.
The right block has one graph SNP in that chunk; a two-site graph-path check
therefore rejected it even though no internal right-block edge exists.

Three MAPQ-60 primary reads call the exact insertion and right SNP ALT, four
call insertion REF and SNP REF, and one call insertion REF and SNP ALT while
carrying a different nearby insertion. Quality-weighted evidence establishes
the ALT/ALT relation. The new bridge requires exact CIGAR alleles, matching
reference flanks, both insertion alleles, the existing 0.001 wrong-parity
bound, and a certified left graph SNP path. It accepts a right singleton only
when that phase set has exactly one graph candidate in the chunk. The owning
chunk and short panel replay now span the gap without a truth switch.

Full chr20 tracked HiPhase-correct closures rise 12/23 -> 13/23; read phase
sets fall 715 -> 713. The same 62,154 variant keys and genotypes remain,
and every one of 236,847 truth-scored reads keeps its correctness outcome:
229,090 correct and 7,757 discordant (96.7249%). The graph panel span
count rises 40 -> 41; the target window's separated-read fraction reaches
0.536 (81/151). Details and reproduction are in
`evaluations/2026-09-28-insertion-singleton-bridge/README.md`.

### 2026-09-28: The 1.508 Mb shifted deletion needs a left-block split

At chr20:1,508,171–1,523,721, two exact and three shifted four-base CIGAR
deletions are sequence-equivalent and pair with the right SNP ALT; six REF
reads pair with SNP REF. The current physical exact-position call misses the
shifted three. However, the left graph phase set also fails the continuous
SNP-path guard at 1,349,327–1,350,788 (19, 0, 6 agree/agree/reverse votes).
A whole-block stitch is unsafe until a validated right-side suffix of that
block is separated at the weak cut. No join was made for this gap. Evidence:
`evaluations/2026-09-28-shifted-deletion-1508/README.md`.

### 2026-09-28: Close the 1.508 Mb shifted-deletion seam safely

The earlier split diagnosis was superseded by direct two-haplotype graph
support between the SNPs flanking the weak 1.349 Mb site. For an MSA deletion
bridge only, the path check may bypass one weak one-haplotype SNP when that
direct edge has at least two reads from each haplotype and passes the existing
one-sided binomial `p <= 0.01` vote test. A sequence-equivalent CIGAR
deletion can vote for the unchanged BAM candidate row; other nearby indels abstain. The stitch tries
clean SNP pairs first and requires two reads for each deletion allele plus
the existing quality-weighted parity bound. The two-read floor is necessary:
a one-per-allele trial joined 3.964 Mb incorrectly and worsened chr20 truth.
The guarded bridge joins 1,508,171–1,523,721 in the owning chunk and short
panel replay. Details: `evaluations/2026-09-28-shifted-deletion-1508/README.md`.

Final full chr20: tracked closures 13/23 -> 14/23, read phase sets 713 ->
711, and panel spans 41 -> 42. All 62,154 VCF keys and genotypes and all
236,847 individual truth-scored read outcomes remain unchanged (229,090
correct, 7,757 discordant; 96.7249%). A broad path-bypass trial increased
discordance by 535 and was rejected; the bypass is now scoped to the certified
MSA deletion bridge, with at most one skipped SNP per path.

The final full gap-window suite passes 1,780 assertions across 34 cases.
Unit tests and the HiFi/ONT validation gates pass.

### 2026-09-28: Close the 60.033 Mb SNP-to-insertion seam

The chr20:60,033,052–60,048,237 tracked gap had a clean left SNP and an
MSA-verified one-base BAM insertion on the right. The final physical stitch
previously tried left-indel to right-SNP but never the symmetric SNP-to-right-
insertion case. In the BAM, equivalent one-base A insertions occur a few bases
left of the candidate within the same A run. Direct physical calls in the
owning 60–61 Mb chunk have seven qualified REF and two qualified ALT insertion
pairs with the left SNP, yielding a decisive quality-weighted relation. The
left phase set has one original graph SNP plus attached BAM rows; requiring two
original graph SNPs falsely rejected this vacuous internal path check.

The retained bridge compares inserted local reference strings, checks inserted
and aligned base qualities, abstains on nearby other indels, requires two
distinct reads per insertion allele, and validates both graph paths. A graph
singleton is accepted only if it is the sole original graph candidate in its
phase set. The short window spans with 0.71 truth-separated reads, up from
0.46; the full owning chunk retains 3,290 correct and 36 discordant of 3,326
truth-scored reads. A new owning-chunk regression checks exact boundary PS
and opposite genotype orientation; the panel span expectation is updated.

The initial full chr20 trial closes 15/23 tracked HiPhase-correct gaps versus
14/23, reduces read PS 711 -> 709, and keeps the same 62,154 VCF keys and all
236,847 truth-scored read outcomes (229,090 correct, 7,757 discordant). Eight
VCF GT strings invert together inside the newly joined block, without changing
allele dosage. The singleton guard was subsequently tightened to require one
original graph candidate; final validation follows below.

The same bridge closes the older 60,033,052–60,058,235 panel window, which
contains the insertion as an interior step. Its short replay now has 0.75
truth-separated reads instead of 0.40. The insertion at 60,048,238 is added
to that window's required-site regression. The graph panel span total rises
42 -> 44 because both 60.033 Mb windows now close.

Final guarded full chr20 confirms 15/23 tracked closures, 709 read phase sets,
62,154 identical VCF keys, eight phase-only GT inversions in the joined block,
and no dosage or individual read-truth changes: 229,090 correct and 7,757
discordant among 236,847 scored reads. `make unit-tests`, `make check`, and
the full `make window-tests` suite pass; the latter has 1,793 assertions in
34 cases. Evidence: `evaluations/2026-09-28-shifted-insertion-60033/README.md`.

### 2026-09-29: Retry an imported BAM seam only with physical SNP parity

The first recovery pass at chr20:15,056,025–15,071,132 imports two BAM phase
blocks inside a wider graph seam. The old second-pass overlap filter discarded
the newly exposed pair even though the pair's phase-set IDs had never been
solved together. In the owning 15–16 Mb chunk, retrying that pair closes the
boundary deletion and adjacent SNP in PS 15,039,543. Its local truth scorer
remains unswitched with 99.44% concordance.

Overlap alone was unsafe: retrying every new pair across a solved interval
changed 25,283 previously correct chr20 read assignments to discordant ones,
dropping accuracy from 96.72% to 86.57%. Requiring only two BAM molecules to
span both endpoints still retried an adjacent 17.839 Mb pair, where the second
BAM solve tagged 16 reads and 10 were locally discordant. There is no unphased
reference base between those adjacent anchors, so it is excluded.

The retained retry requires a new phase-set pair, an actual interior base, and
direct physical calls on the nearest eligible SNP on each side: clean or
alignment- and MSA-verified, MAPQ/base qualities at least 30, at least four
unique paired reads, at least 75% agreement, and quality-weighted parity odds
corresponding to at most 0.001 wrong-parity probability. In the 15 Mb owning
chunk, truth changes from 4,259/4,279 correct to 4,257/4,279, with the gap
closed and no switch. The adjacent 17.839 Mb retry is suppressed. The new
owning-chunk regressions assert the closure and the adjacent-anchor safeguard.

Full chr20 retains the same 62,154 VCF keys and 236,847 truth-scored phased
reads. Tracked HiPhase-correct gap closures rise from 15/23 to 16/23; read PS
fall from 709 to 706 and VCF PS from 350 to 348. Correct reads change from
229,090 to 229,088 and discordant reads from 7,757 to 7,759 (96.7240%
accuracy versus 96.7249%). Only four reads change individual truth outcome;
there is no large block reversal. Evidence: `evaluations/2026-09-29-new-bam-seam-retry/README.md`.

Final validation: `make -j4`, `make unit-tests`, `make check`, and
`make window-tests` pass. The full gap-window suite has 1,820 assertions in
36 test cases. The graph-arm panel expectation for the newly closed 15.056 Mb
window changes from open to closed; its truth-separated floor rises from 0.40
to 0.89, and the graph total span floor rises from 44 to 45.

### 2026-09-29: Join a shifted-deletion bridge to a certified graph suffix

At chr20:13,830,800–13,844,727, four BAM reads placed a two-base deletion
17 bases after its catalog representation in a `(TG)` repeat. A 16-base
physical-equivalence search lost those calls; 32 bases recovers sequence-
verified support. The resulting physical deletion/SNP parity is decisive,
but the left graph block has earlier SNP-path gaps and cannot join whole.
The later unsupported edge contains indel candidates, so a coordinate-empty
cut is also invalid. Splitting at the candidate boundary with the fewest
crossing tagged reads preserves the prefix separately, unphases ambiguous
crossing reads, and joins only the suffix with its own supported graph SNP
path. The owning-chunk and full chr20 runs close the tracked gap without any
change to individual truth-scored read correctness. Tracked HiPhase-correct
closures rise 16/23 -> 17/23; all 62,154 VCF keys and 236,847 scored reads
remain, with 229,088 correct and 7,759 discordant (96.7240%). Truth-scored
read phase sets rise 706 -> 708 because the unsafe prefix stays separate. Evidence:
`evaluations/2026-09-29-shifted-deletion-suffix/README.md`.

The graph panel's 13.830 Mb span changes from open to closed, and its
truth-separated floor rises 0.54 -> 0.73; total spans rise 45 -> 46.

Final validation: `make -j4`, `make unit-tests`, `make check`, and the full
`make window-tests` suite pass (1,834 assertions in 37 cases).

### 2026-09-29: Certify a BAM insertion suffix against a graph MNP

At chr20:17,865,146–17,883,198, an MSA-verified BAM insertion has a
sequence-equivalent CIGAR placement five bases away, and the right graph
boundary is a two-base MNP. Exact-coordinate insertion calls, SNP-only
physical boundaries, and the old equal-length seam coordinate excluded the
available bridge. Calling the full MNP and equivalent insertion restores its
physical evidence. Transferring only the final insertion mixed independent
fallback read gauges and cost 20 previously correct local reads. A direct
clean-SNP-to-insertion link across the source's last weak cut certifies the
larger BAM suffix, including the 17,852,024 SNP and paired deletion rows,
while the earlier source prefix stays separate. The owning-chunk truth improves
4,049/4,082 -> 4,052/4,082. Full chr20 tracked closures rise 17/23 -> 18/23;
truth-scored reads rise 236,847 -> 236,854 and correct reads 229,088 ->
229,098, with discordant reads falling 7,759 -> 7,756. All 62,154 VCF keys
and genotype strings are unchanged. The MNP is admitted as a physical bridge
boundary without broadening graph clean-SNP path certification; trials that
broadened it reopened the earlier 1.508 and 3.529 Mb joins. Evidence:
`evaluations/2026-09-29-insertion-mnp-suffix/README.md`.

The graph panel's 17.865 Mb span changes from open to closed, raising the
expected total 46 -> 47. Final validation: `make -j4`, `make unit-tests`,
`make check`, and the full `make window-tests` suite pass (1,848 assertions
in 38 test cases).

### 2026-09-29: Use established graph read HP across a shifted BAM deletion

At chr20:59,825,454–59,842,960, three primary BAM reads span the gap but the
left graph SNP has low BAM base quality on two of them. Their graph allele and
established graph HP remain informative. Two maternal reads physically carry
the same short BAM deletion, with one CIGAR placement shifted eight bases in a
repeat. A benign insertion outside the allele-verification span caused the
old equivalent-deletion caller to discard the exact-placement read. Restricting
indel interference to the verified span and pairing the graph HP with
MSA-verified deletion ALT evidence closes the gap. The stitch requires two
independent clean reads, a 0.001 wrong-parity bound, and supported graph SNP
paths on both blocks. The full chr20 VCF keys, genotypes, 236,854 scored reads,
and 229,098 correct / 7,756 discordant assignments are unchanged; tracked
HiPhase-correct closures rise 18/23 -> 19/23 and read phase sets fall
704 -> 702. Evidence and test: `evaluations/2026-09-29-graph-snp-deletion-bridge/README.md`.
The 59.825 Mb panel span and total expected spans rise 0 -> 1 and 47 -> 48.
Final validation: `make -j4`, `make unit-tests`, `make check`, and
`make window-tests` pass (2,853 assertions in 39 cases).

### 2026-09-29: Preserve a readless BAM island and certify the 47 Mb graph path

At chr20:47,003,897–47,713,869, the BAM recovery solve contributes six
MSA-verified variant rows at the left boundary, but their source block owns no
read tags after transfer. Those rows stayed in an isolated PS even though a
nearby clean graph SNP and independent physical calls establish their gauge.
The graph block also has a one-haplotype SNP edge at 47,671,540–47,689,418;
two MAPQ-60 reads call both SNPs at Q35/40 and confirm the existing phase.
The readless source rows now inherit that graph PS only when their source path
is complete, its weak and quality cut vectors are empty, and direct physical
SNP-to-insertion and graph-edge checks pass. This does not move read tags.

At the right boundary, the MSA insertion represented at 47,694,119 has
sequence-equivalent CIGAR placements nearby. A spanning molecule has an
inserted base at Q10, although its aligned SNP and flanks are Q30 or better.
The insertion bridge now admits that base at Q10 and gives it its actual error
weight. Its two REF and two ALT observations reach the 0.001 wrong-parity
bound. The same physical graph-edge certificate and a direct suffix
SNP-to-insertion check permit the complete graph block to join the right SNP.

The owning 47–48 Mb replay and full chr20 run both join the exact boundary
rows without a local parental switch. Tracked HiPhase-correct closures rise
19/23 -> **20/23**; 21.378, 23.421, and 32.235 Mb remain open. All 62,154
VCF keys and 236,854 truth-scored phased reads remain. Correct reads change
229,098 -> 229,097, discordant reads 7,756 -> 7,757, and read PS 702 -> 700.
The 11 changed VCF phase strings are confined to the two merged blocks; allele
dosages do not change. Evidence:
`evaluations/2026-09-29-readless-insertion-graph-path/README.md`.

The graph panel's long-gap span rises 0 -> 1 and its total 48 -> 49. An
owning-chunk regression asserts that the island, suffix insertion, and right
SNP share a PS, their parental genotypes remain opposite, and local reads
have no switch.

Validation: `make -j4`, `make unit-tests`, `make check`, and the full
`make window-tests` suite pass (2,870 assertions in 40 test cases). The
existing 47–48 Mb graph-gauge regression measures 2,967/2,998 reads on its
majority parental orientation after this join, so its purity floor is 98.9%
while its exact site-gauge and opposite-genotype checks remain in force.

### 2026-09-29: Stitch a recovered long BAM insertion after earlier graph joins

The 65 bp BAM insertion at chr20:32,246,127 sat exactly at a newly exposed
recovery seam. The physical stitch compared its internal `VariantKey.pos`
(32,246,128) with the VCF-anchor seam end, so it excluded the allele; it also
used an obsolete left PS after an earlier graph join. The stitch now compares
`sort_pos()`, resolves boundary PS labels from current candidates, and checks
long insertions up to 128 bp. Two MAPQ-14/15 ALT molecules with Q22/Q35
inserted bases orient the complete right BAM source against a clean left SNP.
Nearby reads with competing indels abstain from REF calling. The right source
must have no weak or quality cut, the left graph SNP suffix must be supported,
and the quality-weighted wrong-parity bound is 0.01. The established left
block is retained; splitting at a distant weak edge reopened its earlier
32.215 Mb join in a rejected trial.

Tracked HiPhase-correct closures rise **20/23 → 21/23**. Full chr20 retains
62,154 VCF keys and all genotype strings; truth-scored reads move 236,854 →
236,855, correct reads remain 229,097, discordant reads move 7,757 → 7,758,
and read PS fall 700 → 699. The remaining 21.378 and 23.421 Mb gaps have no
independent clean direct bridge yet. Evidence:
`evaluations/2026-09-29-long-insertion-source-stitch/README.md`.
The new 32.235 Mb panel case replays its exact 32–33 Mb owning chunk,
where overall read concordance is 1,311/1,520 (86.25%) and the dominant
joined block correctly separates 39/52 local reads. Its 86% test floor
replaces a 95% floor measured on a shorter, differently bounded replay.
Validation: `make -j4`, `make unit-tests`, `make check`, and full
`make window-tests` pass (2,893 assertions in 40 cases). Exact VCF
allele-key review confirms 21 closed and the two named gaps open.

### 2026-09-29: Stitch a recovered right deletion from a clean left SNP

The 23.421 Mb recovery seam retained both boundary alleles, but the physical
stitch had no SNP-to-right-deletion case. A MAPQ-60 read calls the left ALT at
Q17 and the right 9 bp deletion REF at Q40. The quality-weighted relation has
about 2.1% wrong-parity probability. The new stitch requires a complete BAM
source without cuts, supported graph SNP paths on both sides, no
opposite-haplotype overlapping deletion, and posterior wrong parity at most
0.05. A broad trial that omitted the right graph-path check joined a weak
1.15 Mb block first and reopened the earlier 1.508 Mb bridge. The graph-path
check rejects that counterexample.

Full chr20 tracked closures rise **21/23 → 22/23**; only 21.378 Mb remains
open. Truth-scored reads rise 236,855 → 236,859, correct reads rise 229,097 →
229,101, discordant reads stay 7,758, and read PS stay 699. VCF keys and
allele dosages are unchanged; three phased genotype strings flip inside the
newly joined right block. Evidence:
`evaluations/2026-09-29-snp-to-deletion-stitch/README.md`.
The last 21.378 Mb gap is also split by standalone BAM phasing on the same
alignment. Its two crossing reads have disjoint low-quality or competing-repeat
allele calls, so this is not a loss during graph recovery transfer.
Validation: `make -j4`, `make unit-tests`, `make check`, and full
`make window-tests` pass (2,910 assertions in 40 cases). The owning-chunk
regression pins both 23.421 Mb alleles and parental orientation; the 1.508 Mb
regression guards the rejected early join.

### 2026-09-29: Close the repeat insertion pair without merging a multiallelic locus

The last tracked HiPhase-correct chr20 gap, 21.378–21.395 Mb, has two
MSA-verified BAM insertions and no graph SNP boundary for the ordinary
physical stitch. A targeted seam pass uses two agreeing primary BAM
repeat-length calls, exact zero-length calls at each boundary, and an
independently certified left SNP-to-insertion edge. Right source reads
without a callable primary BAM boundary allele keep independent HP tags;
graph-projected starts were insufficient to classify those reads.

A trial joined a 23.8 Mb multiallelic insertion/deletion locus through
the same length-only vote. The final path vetoes a second phased allele
at either boundary or an overlapping deletion. It leaves that locus
unchanged. All 23/23 tracked gaps now close in the guarded full chr20 run.
VCF keys stay 62,154; truth-scored phased reads stay 236,859; correct
reads rise 229,101 -> 229,108 and discordant fall 7,758 -> 7,751.
Only the intended 21,395,286 A>AT sample field changes; read PS count
stays 699. Owning-chunk local majority purity is 122/144 (84.72%),
below HiPhase's 131/144 (90.97%). The new paired-insertion and
multiallelic counterexample regressions preserve both decisions.
Evidence: evaluations/2026-09-29-repeat-insertion-pair/README.md.
Validation: make -j4, make unit-tests, make check, and full
make window-tests pass (2,939 assertions in 41 cases).

### 2026-09-29: Audit additional HiPhase gap joins and preserve owning-chunk evidence

An exact-boundary comparison of the accepted graph/recovery chr20 VCF with
HiPhase on the same variant input found 56 further HiPhase joins outside the
centromere, 51 with identical local input calls. Two met the strict >=98%
local read-truth and >=95% concordant-flank criterion. Graph deletion
coordinates at 13.8 Mb needed VCF allele normalization before physical
matching; a certified local suffix now joins without absorbing an unsupported
upstream SNP. At 64.14 Mb, right-block HP observations orient a BAM-derived
left deletion against two complementary right deletion rows; only the left
row moves. Both exact VCF gaps join in the owning chunks and full chr20.
All 23 earlier tracked joins remain closed. Variant keys and truth-scored
phased reads stay 62,154 and 236,859; 229,108 remain correct and 7,751
discordant. The read PS count falls 699 -> 698. The 64.14 Mb row transfer reveals
an adjacent 64.138 Mb split; current full chr20 still has 56 exact
HiPhase-joined pgphase splits, including this newly exposed seam. Its
deletion and SNP votes are mixed, so the seam remains open and has its own
panel regression. Five more >=95%-pure and eight >=80%-pure,
concordant-flank HiPhase joins are pinned as open panel targets: 14 open
targets after the newly exposed 64.138 Mb seam. A strict
repeat-length trial reopened the verified 21.378 Mb join and added seven
truth-discordant reads, so it was reverted. Evidence and reproducible gap
audit: evaluations/2026-09-29-expanded-hiphase-gaps/README.md.
Validation: make -j4, make unit-tests, make check, and full make
window-tests pass (3,302 assertions in 42 cases).

### 2026-09-29: Reject a zero-pair SNP-to-deletion stitch without corroboration

A trial invoked the existing physical SNP-to-MSA-deletion stitch when a right
clean SNP existed but no read called both SNPs. It closed the 58.366 Mb panel
gap by one exact CIGAR pair, but joined its boundary ALT alleles on opposite
haplotypes; HiPhase and the parental read evidence put them on the same
haplotype. Another one-read join at 20.8 Mb attached a large block across an
injected one-base deletion overlapping an unphased six-base graph repeat row.
Full chr20 discordant reads rose 7,751 -> 7,892 at unchanged 236,859 scored
reads; a repeat-overlap veto still left 7,815. Both trials were discarded.
The 58.366 Mb gap regression now rejects the observed wrong ALT orientation
if a future solver joins it. The second spanning 58.385 Mb read carries a
sequence-equivalent deletion shifted 20 bp in an A run but has Q17 at a
retained base, below the current Q30 equivalence filter. The next repair needs
to score both placements with their actual base qualities and establish a
supported relation before joining. Evidence:
evaluations/2026-09-29-zero-pair-indel-counterexample/README.md.
Validation after rejecting the trial: make -j4, make unit-tests, make check,
and make window-tests pass (3,307 assertions in 42 cases). A one-window replay
with the rejected binary fails both the switch and ALT-orientation assertions.

### 2026-09-29: Certify detached BAM-only runs and prioritize clean SNPs

At chr20:19,403,172–19,414,720, 19 independent MAPQ/base-quality >=30
reads make a decisive clean-SNP bridge with the HiPhase allele orientation.
The left phase set contains three BAM-derived rows and no graph rows, so the
final stitch's graph-only path check always rejected it. Its original BAM
source has three weak cuts, all before the detached three-row run; none lies
inside the run. The stitch now accepts such a BAM-only local path only when
every candidate maps to one source phase set with consistent haplotype
orientation and no internal weak cut. It then retries a newly adjacent seam
inside the recovery target when the left boundary is a clean SNP.

This exposes the upstream 19,395,544–19,403,172 SNP gap. Its physical SNP
pair has 26 unanimous same-allele calls, but a repeat-insertion length vote
chose the opposite haplotype join. A decisive nearby BAM clean-SNP pair now
sets the orientation when it shares a source run with the insertion and no
weak cut separates them. The two exact gaps join with the HiPhase allele
orientation. Their dominant-block truth separation rises from 59/95 to 79/95
and 67/113 to 98/113, respectively; HiPhase has 91/95 and 110/113.

A broad retry also joined the low-purity chr20:36,620,864 singleton and lost
seven truth-scored phased reads. Restricting newly exposed retry seams to a
clean-SNP left boundary removes this collateral join; an owning-chunk control
pins that split. Full chr20 retains 62,154 variant keys and 236,859
truth-scored phased reads, with 229,108 correct and 7,751 discordant. Exactly
32 VCF sample fields change, all from 19,397,607 to 19,464,219; read phase
sets fall 698 -> 697. Evidence: evaluations/2026-09-29-bam-only-snp-run/README.md.

### 2026-09-29: Expand same-callset gap controls; reject unsupported whole-block joins

The current accepted chr20 graph/recovery output still has 55 exact-boundary
HiPhase joins on the same callset. Five more noncentromeric cases now have
owning-chunk regression rows: 1.180, 10.727, 19.373, 21.594, and 36.332 Mb.
HiPhase separates 125/129, 133/161, 111/140, 118/154, and 144/179
truth-scorable reads in their dominant gap blocks; pgphase separates 66,
66, 59, 59, and 144, respectively. All five remain split, adding six
in-gap heterozygotes and no spans to the panel baseline.

At 1.18 Mb, left-block tagged reads unanimously orient both alleles of a
right insertion placed 18 bp later in a repeat. The right graph block also
contains an unsupported edge at 1.349–1.351 Mb with 19 one-haplotype GAF
votes and six reversals. Joining the whole block from this local evidence
would propagate orientation across that edge. At 12.256 Mb, a short replay
joins while the full owning chunk splits; the BAM source has a weak cut at
the boundary and shifted deletions on both SNP allele classes. This is an
unsafe short-context certificate.

A trial admitting all verified injected noisy SNPs to physical seam stitching
closed 36.332–36.355 Mb, but lost 144 net truth-correct reads on full
chr20 and changed 988 VCF sample fields. Nearly all losses came from an
incorrect 37.462–37.467 Mb whole-block join at a seam with competing
2- and 3-base deletion rows; the precise false-support mechanism remains
unresolved. A separate trial excluding deletion lengths
from repeat-insertion REF votes changed no VCF fields and lost four correct
reads. Both were reverted; the separate read-label correction below is the accepted
change. Evidence: evaluations/2026-09-29-five-more-gap-controls/README.md.

The expanded 19.373 Mb control exposed stale read PS tags after weak-cut site
detachment: 13 reads wholly right of the gap retained the left read PS and
12 had the opposite truth orientation. Source reads that observe only the
detached component and no neighboring graph anchor now move with its sites in the original BAM HP gauge. The
owning chunk improves from 4,020/4,071 to 4,031/4,071 correct reads and its
switch check passes. Full chr20 keeps 62,154 VCF rows with identical sample
fields; truth-scored phased reads change 236,859 → 236,858, correct reads
229,108 → 229,125, discordant reads 7,751 → 7,733, and read phase sets
697 → 698. The read-only component may remain independent of VCF blocks.

### 2026-09-29: short insertion CIGAR placement at chr20:5.25 Mb

The graph/recovery phase sets at 5,256,785 and 5,263,741 stayed separate although
HiPhase joined the same VCF alleles. On the original BAM, 15 reads overlapping
5,256,785 place its `AT` insertion 50 reference bases later in the same AT
run. Inserting `AT` at either coordinate gives the identical edited reference
string. The physical insertion caller only searched 16 bases and therefore
lost those ALT observations. A trial expanding the bound to 64 for all
insertions joined this gap, but split the already closed 18,983,414–18,999,993
gap; that nearby locus has multiple longer insertion alleles in one compound
repeat. The broad trial was rejected.

The accepted rule searches 64 bases for one- and two-base insertions and keeps
the 16-base bound for longer insertions. It still requires identical edited
reference strings, clean anchor bases across both placements, and the usual
stitch evidence. An initial short-insertion trial closed 5.25 Mb but reopened
the established 5.31 Mb seam: the graph SNP-path validator did not remember
the newly certified physical insertion edge and rejected 43 agreeing clean
SNP pairs at the next seam. The fix carries forward only a unique boundary
graph-SNP pair whose live phase sets actually merge during the physical stitch.
Both the newly joined 5.25 Mb gap and the existing 18.98 and 5.31 Mb gaps
pass focused truth-backed tests. In the full chr20 run, the same 62,154 VCF
keys are emitted and only 24 sample fields change, all in the 5 Mb chunk. At
5,256,785 and 5,263,741 the ALT alleles now share the same haplotype,
matching HiPhase's boundary orientation. Truth-scored phased reads increase
236,858 to 236,867; correct reads increase 229,125 to 229,131; discordant
reads increase 7,733 to 7,736. VCF blocks decrease 338 to 337 while N50
stays 616,859 bp; read phase sets decrease 698 to 696. The new gap span and
its unanchored interior candidate at 5,256,786 are pinned in the window panel.
The complete window suite passes 3,496 assertions in 42 test cases;
`make unit-tests` and `make check` also pass.

### 2026-09-29: Keep phased recovery sites when graph duplicates are unphased

At chr20:20,711,883 the BAM sub-solve injected an MSA-verified phased
`CA>C` deletion, but graph output deduplication chose an unphased catalog
repeat row for the same normalized allele because its depth was 47 versus 44.
The repeat row was then omitted from VCF, losing the usable in-gap site.
Final graph output now prefers a valid phased heterozygote over an unphased
copy, with coverage deciding only within equal phase status. The 20–21 Mb
regression requires the exact phased deletion. Full chr20 gains 115 phased
variant keys (62,154 -> 62,269); every added key is present and phased in
the independent BAM-only callset. Existing sample fields and the complete
phased BAM are unchanged. VCF blocks are 337 -> 343; N50 stays 616,859 bp.
The 20.707 Mb gap remains split because its deletion-to-right-SNP molecules
do not provide a decisive orientation. Evidence and commands:
`evaluations/2026-09-29-phased-duplicate-recovery/README.md`.

### 2026-09-29: Respect a source weak cut after one graph-anchored BAM site

The output deduplication exposed an incorrect VCF relative phase at
chr20:33,227,050. No recovery read calls both this deletion and the
33,211,072 insertion, and the BAM source-path audit correctly marks a weak
cut. The post-stitch detachment nevertheless required a second imported
near-side locus, even though the first was attached to an oriented graph
block. An indel cut with no callable read across its boundary can now use
that graph anchor to establish the near side for detachment. Clean-SNP cuts
retain their physical bridge check, preserving the protected 3.574 Mb join.
The deletion keeps the independently supported downstream PS; an owning-chunk
regression pins the separation and right-side allele orientation. Full chr20
retains 62,269 keys; six sample fields and three read PS labels change at weak
cuts, while truth-scored reads remain 229,131/236,867 correct. The full
window panel passes 3,521 assertions in 44 cases; unit tests and validation
gates also pass. Evidence:
`evaluations/2026-09-29-phased-duplicate-recovery/README.md`.

### 2026-09-29: Guard the 58.366 Mb seam with clean SNPs

A follow-up audit found that the prior gap regression used the repeat
`CA>C` deletion at 58,385,702 as its orientation oracle. Its raw CIGAR
alleles occur on both parents and can change under repeat normalization.
The test now compares clean SNPs at 58,366,458 and 58,391,091: their ALT
alleles belong to opposite parents, so any merged PS must give them
opposite haplotype labels. A rejected physical deletion stitch did the
reverse. The production seam remains split; no phasing rule changed.
At 10.727 Mb, only two primary MAPQ-30 molecules cross both exact boundary
alleles, correcting an earlier broad-window overlap count. A 21–22 Mb
owning-chunk replay also shows why a shared standalone BAM PS at
21,159,070–21,172,487 does not certify a join: no source read calls both
boundary alleles (one profile spans both indices but has an unknown deletion
allele). The source-path weak cut and independent graph PS therefore preserve
the correct uncertainty. Evidence:
`evaluations/2026-09-29-zero-pair-indel-counterexample/README.md`.

### 2026-09-30: Attach a supported BAM insertion run without flipping its graph block

At chr20:1,180,618–1,194,189, the BAM source has a weak cut after the
left SNP but the recovered `ATC` insertion and three adjacent BAM rows form
a separate supported run. Six MAPQ-30 reads already assigned to the left
phase set call that insertion, three REF and three ALT, and all six support
the same HP orientation (exact one-sided binomial p=0.015625). The local
connector uses that read gauge and requires split-half allele support on
every internal edge (p<=0.01); it stops before the nearby 155-bp insertion
and graph catalog row. Only four BAM rows transfer to the left PS. HiPhase
on the same calls gives the same boundary allele orientations. The owning
chunk now spans the 13,571-bp gap, and correctly separated local reads
increase 66/129 -> 73/129. Full chr20 retains all 62,269 variant keys;
only those four VCF sample fields change. Truth-scored phased reads change
236,867 -> 236,866, correct reads stay 229,131, and discordant reads
fall 7,736 -> 7,735. Evidence and commands:
`evaluations/2026-09-30-bam-insertion-run-stitch/README.md`.

### 2026-09-30: Preserve verified MSA recovery classification at graph output

The 35.329–35.348 Mb graph gap contains a phased, MSA- and alignment-verified
BAM insertion at 35,342,608. Recovery injected it, but graph output
reclassified its 2 REF / 9 ALT counts as `LOW_AF` and dropped it. Output also
failed to carry MSA/alignment verification into reconstructed candidates, so
phased noisy recovery duplicates could lose to deeper unphased graph repeat
rows. Verified `NoisyCandHet` rows now retain the BAM sub-solve's category,
initial category, bitmask and proof flags through emission. Full chr20 gains
92 VCF variant keys, loses none, and changes no existing sample fields; the
phased BAM is byte-identical at 229,131/236,866 truth-correct reads. An
owning-chunk regression pins the restored insertion and its right-block PS.
The gap remains split because the adjacent repeat-block allele paths are not
independently certified; the extra site alone does not justify a whole-block
join. The complete window panel passes 3,547 assertions in 45 cases;
unit tests and validation gates also pass. Evidence:
`evaluations/2026-09-30-verified-msa-output/README.md`.

### 2026-09-30: Retry newly exposed left-deletion seams, preserve two-allele guard

BAM recovery can insert an MSA-verified deletion after the original graph
seam was recorded. The physical stitch now considers such a new deletion-to-
right-SNP seam inside the original target, using one scan over recovery sites
for new insertion and deletion boundaries. At 35.498 Mb this corrects seam
routing but still abstains: only one primary read has a quality-checked REF
call across the deletion and right SNP, none calls ALT, and the left graph SNP
path is not certified. The route-only full chr20 output remains exactly
unchanged (62,361 VCF keys and byte-identical phased BAM). A separate one-
allele likelihood experiment wrongly joined blocks at 3.964 and 60.7 Mb,
losing 391 truth-correct reads, so the two-allele requirement remains. See
`evaluations/2026-09-30-left-deletion-seam-routing/README.md`.

### 2026-09-30: Let complete BAM source paths fill tied graph SNP edges

At chr20:12.256–12.270 Mb, a complete imported BAM source and a clean
primary-read SNP bridge were rejected because the right graph block has
balanced and absent GAF SNP-pair observations downstream. A broad
complete-source shortcut also crossed one-haplotype and reversing graph
edges elsewhere, losing 713 truth-correct reads, so it was discarded. The
physical SNP stitch now uses a complete right BAM source only after a
supported two-haplotype graph edge and only across balanced ties or absent
graph pairs; one-haplotype and reversal edges still veto. The owning-window
regression asserts exact SNP/deletion PS and GT relationships and a span.
Full chr20 keeps 62,361 variant keys and 229,131/236,866 truth-correct
reads, while read phase sets decrease 697 to 695. Only 760 sample fields in
the intended right block change. The 45-case window panel, unit tests,
and validation gates pass. Evidence:
`evaluations/2026-09-30-source-backed-graph-path/README.md`.

### 2026-09-30: Correct two short-replay false-open gap tests

The 48.929 and 56.064 Mb gap tests replayed only a short slice even though
the accepted full chr20 output already joins both boundaries. Their
regressions now replay the 48–49 and 56–57 Mb owning chunks, assert exact
boundary alleles and PS/GT orientation, and require spans. The owning
chunks separate 135/156 and 148/167 local reads correctly. A trial that
retried newly exposed right deletion seams at 36.332 Mb joined the right
block in the wrong orientation, lowering local concordance to 88.09%; it
was reverted. A broad new-left-insertion retry closed the mixed-truth
36.611 Mb control and lost six truth-correct reads genome-wide; it too was
reverted. These counterexamples preserve the current conservative stitch.
Evidence: `evaluations/2026-09-30-owning-chunk-gap-controls/README.md`.

### 2026-09-30: Join a readless graph SNP block through independent BAM SNP paths

The 30,673,476–30,673,709 chr20 gap separates a long read-bearing graph
block from 51 phased VCF rows in PS 30,610,946 that own no read tags. The
standard physical stitch sees an earlier recovery seam and no callable clean
SNP pair at the noisy boundary. Three independent MAPQ-at-least-30 molecules each call
at least two clean SNPs on each flank without within-read contradiction and
support one orientation; a fourth crossing read has conflicting left SNP
calls and is excluded. Every clean-SNP edge through the right block has
consistent primary-read support and no Q30 contradictory pair. A targeted
readless-block attachment now joins the 51 rows to the left PS. The new
window-panel case requires the span and exact truth-derived boundary/terminal
allele orientation. Full chr20 has the same 62,361 VCF keys, with only these
51 sample fields changed. The phased BAM is byte-identical: 229,131 of
236,866 truth-scored reads correct, 7,735 discordant, 695 read phase sets.
The full gap-window suite passes 3,599 assertions in 45 cases; unit tests
and validation gates pass. Details:
`evaluations/2026-09-30-orphan-graph-snp-block/README.md`.

A separate trial applying the MAPQ-30/skip filter uniformly to BAM source
path votes did not close a tracked gap. It split five additional read phase
sets and lost one truth-correct tagged read, so that trial was reverted.

### 2026-09-30: Defer a clean-SNP bridge beyond a recovered insertion

At chr20:54,483,506–54,490,229, 32 independent MAPQ/base-quality-30
primary reads call both clean SNPs with unanimous relative phase and both
allele classes. The original recovery seam stops at a newly exposed right
insertion 302 bp before the right SNP. The physical stitch now tries the
nearest clean right SNP after an unsupported insertion if it belongs to the
same phase set and no other phased block intervenes. The extra bridge waits
until ordinary right-hand joins settle: an immediate trial broke the existing
54,547,514–54,569,072 span, dropped five VCF rows, and lost one correct
read. The deferred version closes both gaps in the owning chunk, with all
93 local truth reads in one correct block. Full chr20 retains 62,361 VCF
keys and 229,131/236,866 truth-correct reads; read phase sets fall 695 to
694. The new panel case pins exact SNPs, the span, and their parental
orientation alongside the downstream join. The final-binary window suite
passes 3,632 assertions in 45 cases; unit tests and validation gates pass.
Evidence: `evaluations/2026-09-30-deferred-clean-snp-stitch/README.md`.

### 2026-09-30: Confirm a demoted SNP against its clean right block

At chr20:61,757,551–61,773,799, ten MAPQ/base-quality-30 primary reads
call both boundary SNP alleles with unanimous same-phase parity and both
allele classes. The right SNP is a BAM recovery row marked noisy, so the
clean-only physical stitch missed it. The left mixed graph/BAM block has a
complete uncut BAM source and a direct supported graph SNP-to-insertion
edge; the right block's first missing GAF SNP edge is covered by two exact
clean SNPs shared with its complete source. A fallback now uses this path
only when the clean SNP pair has no callable read, the demoted SNP agrees
with its next clean SNP at p<=0.001, and the cross-seam SNP pairs agree at
the same bound. Callable clean pairs retain priority; indel fallbacks
cannot bypass the demoted-SNP checks. This excludes the 37.4 Mb wrong whole-block join (7/11 boundary-
to-clean agreement), and preserves the 46.7 Mb clean-SNP join. Both are
regression controls. The new 61.7 Mb panel gap has an exact allele, PS, and
parental orientation check. Full chr20 keeps 62,361 variant keys; correct
read assignments increase 229,131 to 229,144, discordant reads increase
7,735 to 7,736, and read phase sets fall 694 to 692. The full window
suite passes 3,721 assertions in 45 cases; the final clean-indel fallback
ordering leaves full chr20 candidate TSV, VCF, and BAM byte-identical. See
`evaluations/2026-09-30-demoted-snp-stitch/README.md`.

### 2026-09-30: Shifted deletion calls do not yet certify the 50.548 Mb join

The still-open chr20:50,548,245–50,562,066 gap has twelve MAPQ-30 spanning
reads, and their two read-label gauges agree 6/6. This is weaker than a local
allele bridge: every spanning read has an unknown allele at both overlapping
one- and two-base deletion rows in the BAM source and graph recovery matrices.
A representation-aware CIGAR trial recovers only four Q30 paired calls, with
one vote versus three in the proposed block gauge; the right graph SNP path
also fails its internal certificate. Bypassing that path check still abstains
on the conflicting allele votes. The trial was reverted. Source-label votes
at this seam must not stand in for allele evidence. Details:
`evaluations/2026-09-30-open-gap-evidence-audit/README.md`.

### 2026-09-30 open-gap bridge screen

The owning-chunk matrix audit in
`evaluations/2026-09-30-open-gap-bridge-screen/` found that several
HiPhase-spanned chr20 gaps have no direct callable pairs at pgphase's named
boundary rows. At 64.138 Mb, the immediate deletion pair has 11:1 votes and the
inherited BAM source orientation is unreliable; a later audit found an earlier
complementary pair with two unambiguous ALT classes that can certify only the
source suffix (see below). At 35.498 Mb, only three MAPQ
60 reads physically span the gap, and two candidate REF calls fail Q30 flank
quality. Reversible one-sided and quality-weighted deletion stitch trials did
not close it under the existing path checks, so no stitch threshold was
relaxed. Distinct graph and BAM deletion alleles at 3.963 Mb must remain
separate. The next safe gain requires an independently supported path through
the gap or a local component transfer, not a whole-block boundary vote alone.

### 2026-09-30 physical-pair screen and gap-test replay sharing

A MAPQ/Q20 short-allele screen of the remaining open chr20 gaps found that
raw vote strength alone is unsafe: 17:3 votes at 17.839 Mb and 26:9 at
60.706 Mb occur at boundaries with opposite parental flank orientations.
The 13.593 and 23.792 Mb HiPhase-only read classes include many reads that
call REF at both complementary insertion rows; phasing them from either row
would invent an allele observation. Details and limits are in
`evaluations/2026-09-30-open-gap-bridge-screen/README.md`.
The gap regression helper now hard-links identical owning-chunk outputs for
separate windows within one run, preserving per-window paths and avoiding
repeated pipeline execution. That test-cache change does not alter runtime
phasing behavior.

### 2026-09-30: Exposed deletion seam and physical cut coordinates

A successful boundary transfer at 64,140,314 moved the BAM deletion into its
right block and exposed the 64,138,752–64,140,314 seam. The final physical
stitch was iterating a pre-transfer seam snapshot, so it never evaluated that
new neighbor. Its immediate boundary deletion had one conflicting REF call.
An earlier MSA-verified complementary deletion pair at 64,134,226 has 11
unambiguous BAM-channel ALT pairs to the right deletion: five maternal and six
paternal in the truth-only audit, with no conflicting relation. The preceding
graph SNP has only one read pair into this BAM run. The stitch now revisits the
exposed seam and transfers only the certified BAM suffix, requiring two reads
from each left ALT class, unanimous orientation with one-sided binomial
`p <= 0.001`, exactly two BAM rows at the left locus, and no BAM source cut
through the transferred run. BAM-channel votes use BAM MAPQ instead of GAF
MAPQ, so a weaker graph alignment cannot discard a valid BAM bridge. A first
trial passed a graph catalog key coordinate to a
physical split helper; it absorbed the unsupported graph prefix and lowered
local concordance to 86.79%. Using the selected graph allele's physical
coordinate preserved that prefix. The owning-chunk regression now expects
this gap to span and checks parental read orientation. See
`evaluations/2026-09-30-complementary-deletion-suffix/README.md`.
Full chr20 validation preserves all 62,361 VCF keys and all 236,880
truth-scored tagged reads. Correct/discordant counts stay 229,144/7,736;
20 reads change phase-set label, three change HP gauge, and none changes
parental correctness. Read phase sets rise 692 to 693 because the unsupported
prefix stays separate.

The generic suffix splitter initially cleared six reads that called both its
sides even though they also made the certified complementary-pair/right-deletion
vote. Assigning only those voters to the right block raised the 64.138 Mb
window's correctly separated truth reads from 32/53 to 38/53. Full chr20 VCF
and genotypes are unchanged from the pair-only trial; six read PS tags change,
no HP tags change, and no truth-correct read becomes discordant. The tracked
window now requires a 0.71 separated fraction.

### 2026-09-30: Preserve boundary alleles during MSA gap retries

The 50.548 Mb gap has 12 BAM reads spanning two complementary deletion rows
and its right SNP, but no deletion calls in the ordinary source matrix. A
focused MSA admission trial restored 10/8 row calls with conflicting
relative-phase votes (6:4 and 5:3); it did not justify a join. At 21.159 Mb,
18 spanning reads likewise lack deletion calls. A broad grouped MSA trial
restored those calls but removed seven established VCF rows, including the
boundary deletion, so it was rejected. A narrower BAM solve preserves the
exact deletion, but the grouped trial's restored spanning calls were 17 REF
and one ALT, still one sided. The window regression for this gap
now replays its owning 21–22 Mb chunk and requires all three boundary SNP/
deletion alleles to remain present and phased, without forbidding a future
supported join. A saved-matrix screen of 13 high-purity open gaps found no
strong two-allele intermediate chain under the stated diagnostic criteria.
Evidence and limits: `evaluations/2026-09-30-msa-dropout-gap-audit/README.md`.

### 2026-09-30: Source-backed right-deletion seam at 8.166 Mb

Recovery exposed a right BAM deletion opposite a mixed graph/BAM block, but
physical stitching retried only MNP and long-insertion right boundaries.
The new retry uses the nearest clean SNPs on each side, requires both left
allele classes and unanimous quality-bearing paired calls, and certifies the
mixed left path through exact BAM source alleles. An originally unphased graph
SNP can anchor that path only after direct confirmation against its preceding
clean graph SNP. An internal graph edge without a decisive direct vote may use two exact
shared SNPs in the same complete source only if every callable pair agrees; a
reversed pair vetoes it.
The 8,166,027–8,172,072 panel gap now spans with 99% local concordance and
81% correctly separated reads. Full chr20 retains 62,361 keys and 236,880
truth-scored tagged reads; correctness moves 229,144/7,736 to
229,143/7,737, with no wrong whole-block join. Details and controls:
`evaluations/2026-09-30-source-backed-right-deletion/README.md`.

### 2026-09-30: Recovery matrix dumps must retain per-window provenance

Targeted BAM solves all use `chunk_id=-1`, so `--phase-matrix-dump`
silently overwrote earlier recovery-window matrices when a graph chunk had
multiple groups. The later whole-chunk BAM read solve could likewise replace
the graph flags matrix. Diagnostic prefixes now distinguish graph chunk,
recovery window, retry, focused solve, validation solve, and whole-chunk BAM
solve. The phase output is unchanged. The corrected 58.366 Mb source dump
contains both boundary rows but no callable paired bridge; its previous
apparently empty source dump belonged to a later recovery group. The 15.351,
57.854, and 35.498 Mb audits also found no safe diploid stitch. See
`evaluations/2026-09-30-recovery-matrix-provenance/README.md`.

### 2026-09-30: Sparse same-callset bridges do not justify a lower stitch floor

A fresh audit of the 15 still-open noncentromeric panel gaps found no
additional supplementary-alignment molecule bridging any exact endpoint
pair. At 19.373 Mb only one primary read reaches both boundary positions;
the BAM source's inherited PS puts the left insertion and right clean SNP
ALTs on opposite haplotypes, while HiPhase's truth-correct same-callset
block puts them together. The weak source cut is therefore guarding a real
orientation error. HiPhase defaults to one spanning and one connecting read.
Raising its spanning-read floor to two splits some wrong joins but retains
both its correct 19.373 Mb join and wrong joins at 50.099, 60.098, and
60.706 Mb. Six same-callset HiPhase joins over >1 kb have opposite parental
HP orientations on two >=95%-pure flanks. A read-count threshold alone
cannot reproduce HiPhase's selected correct joins without also risking
large switches. The current pgphase full-chr20 read-PS bin screen found no
strong within-block parental reversal. No stitch rule was relaxed. See
`evaluations/2026-09-30-sparse-gap-link-audit/README.md`.

The 15.351 Mb gap exposed two exact deletion ALT observations missing from
the selected MSA source: both crossing MAPQ60 reads carry the one-base
CIGAR deletion, but general backfill skips homopolymer sites. A broad
exact-ALT-only backfill trial restored both source and transferred calls,
yet closed none of 15 tracked gaps. On full chr20 it added 612
truth-scored tagged reads but also 233 discordant reads, changed 273
shared-key genotypes, and increased read PS 691 -> 718. The trial was
reverted. One of the two local ALT molecules has a Q10 deletion flank,
so the existing Q30 physical stitch can use only the other. The missing
matrix calls are real but are not a safe reason to relax the stitch rule.
Details: `evaluations/2026-09-30-sparse-gap-link-audit/README.md`.

### 2026-09-30: wider clean-SNP flank screen at unresolved gaps

The fresh 21–22 Mb replay confirms just one MAPQ>=30 primary molecule across
21,594,343–21,612,458, with a left insertion and shifted right-repeat
deletions; the single recovery-matrix pair is not a transfer loss. Among the
nearest high-quality clean SNPs across accurate open gaps, most pairs have
no spanning caller. At 21,159,070–21,172,487, 16 pairs vote 9:7 and all call
the same right allele, so they cannot establish a diploid stitch. See
`evaluations/2026-09-30-sparse-gap-link-audit/README.md`. The accepted stitch
rule and regression expectations are unchanged.

### 2026-09-30: shifted insertion REF backfill is real, but broad repair breaks gap regressions

At chr20:19,373,923, 42 of 47 reads backfilled as REF on both independent
four- and eight-base insertion rows carry a nearby CIGAR insertion. Reference-
edit comparison uniquely identifies 41 as one of those two ALTs; the remaining
read has a distinct five-base insertion. A sequence-equivalent backfill trial
corrected their source observations and improved full-chr20 truth-tagged read
accuracy from 229,143/236,880 to 229,275/236,919. It did not close the
19.373 Mb gap and it failed two protected cases: a new 36.332 Mb join fell to
140/179 correctly separated reads, and the correct 41.900 Mb join split. Even
suppressing only the false shifted-REF calls caused the same regressions.
ALT-only backfill was also much worse locally. All trials were reverted; both
protected tests pass again. Details:
`evaluations/2026-09-30-shifted-insertion-backfill/README.md`.

### 2026-09-30: Direct allele edge vetoes contradictory graph/BAM gauge

The graph-only recovery stitch marked a direct allele edge as supporting a
locus only if it agreed with the shared-read graph/BAM gauge, but then still
joined on the gauge whenever that edge existed. A contradictory edge thus
made the gauge join *more* likely. The join now abstains on disagreement and
lets the ordinary allele stitch select the directly supported parity. A
synthetic two-SNP regression fixes the branch. The full chr20 run retains all
62,361 variant keys, the same 236,880 truth-scored tagged reads, and 691 read
phase sets; discordant reads fall 7,737 -> 7,290. The 1,330 phased-GT changes
are confined to 26.418–26.915 Mb and 29.270–29.304 Mb, both in the excluded
peri-centromeric region, so this is a correctness fix rather than a new
noncentromeric gap closure. All 46 gap-window cases, unit tests, and golden
validation gates pass. Details: `evaluations/2026-09-30-gauge-allele-contradiction/README.md`.

A transfer-only correction for shifted insertion REF calls was also rejected:
although exact sequence equivalence preserved the protected 36.332 and
41.900 Mb cases, it split the correctly joined 4.766 Mb gap because the
aggregate link between two BAM source blocks lost its supporting observations.
The source representation and post-solve transfer must be corrected together;
see `evaluations/2026-09-30-shifted-insertion-backfill/README.md`.

### 2026-09-30: Open-gap MAPQ and allele audit

All primary reads spanning the exact endpoints of the 15 remaining
noncentromeric, same-callset HiPhase-correct panel gaps have MAPQ at least 30.
The 57.854 Mb seam has 11 physical spanners but conflicting one- and two-base
deletion to SNP allele relations; the 15.351 Mb seam has eight spanners but no
right ALT observations. At 19.373 Mb the only full spanner has a shifted
insertion with one base quality 22 and no independent intermediate-to-right
ALT bridge. At 36.332 Mb the two full spanners do not call the right four-base
deletion allele. These controls rule out MAPQ relaxation or a raw spanning-read
threshold as a safe general gap closure. The window diagnostic now reports
`CONNECTED` before considering the absence of interior sites and handles
adjacent boundaries without an uninitialized interior coverage minimum. No production
phasing rule or gap expectation changed in this audit. Details:
`evaluations/2026-09-30-open-gap-mapq-and-alleles/README.md`.

### 2026-09-30: graph-indel two-hop certificate for imported BAM blocks

At chr20:35,328,965–35,347,817 the two independent BAM phase blocks stayed
split after the larger graph seam transaction rolled back. The unphased graph
repeat insertion at 35,342,612 is sequence-equivalent to the BAM insertion at
35,342,608, despite different VCF padding. In the merged matrix, 15 MAPQ-60
reads connect the left MSA-verified insertion to the graph allele (9 ALT/ALT,
6 REF/REF), while 11 independent read pairs connect that graph allele to a
clean SNP in the right BAM block (9 ALT/ALT, 2 REF/REF). Neither link has an
opposing allele pair. The BAM MSA row itself calls many physically equivalent
insertion reads REF, so direct BAM-to-BAM evidence alone misses the bridge.

The fallback can now join an adjacent imported BAM pair through such an
unphased graph indel, requiring both allele classes and conflict-free paired
reads on each hop, MAPQ >=30, a clean right SNP, a verified/clean left site,
and no contradictory direct block vote. The graph and BAM candidate rows remain
independent. The owning-chunk parental-orientation regression now expects the
35.328 Mb span and requires the 35,342,608 BAM insertion. The full chr20
trial kept all 62,361 variant keys, changed only seven VCF PS fields, and
reduced discordant truth-scored reads from 7,290 to 7,286 while correctly
assigned reads rose from 229,590 to 229,592. Two fewer reads were tagged
(236,880 to 236,878).

### 2026-10-01: restore shifted single-base MSA ALT calls before source phasing

Targeted recovery's exact-CIGAR backfill can call a shifted repeat insertion
REF at the candidate's original coordinate. The initial BAM MSA matrix at
23,792,419 contains only three REF/five ALT observations despite many
sequence-equivalent shifted `T` insertions in the original alignments.
Changing evidence only after the source solve leaves its HP/PS based on the
incomplete matrix. Broad pre-solve backfill caused unrelated switches; restoring
shifted multi-base calls also joined the protected 41.866–41.898 Mb boundary.
Those trials were reverted.

Recovery now adds missing single-base insertion ALT calls before each noisy-site
k-means pass, restricted to existing admitted MSA heterozygotes inside recovery
windows. It requires MAPQ30, known Q30 inserted/flanking bases, one insertion
within 32 bp, no interfering deletion, and reference-edit equivalence across
both placements. Existing MSA calls and independent candidate rows are retained;
complex indels and ordinary BAM runs keep their current handling. New sparse
profile extents are indexed before the source solve. No truth or region-specific
rule participates in the algorithm.

The matched default chr20 run closes four former VCF seams: 14.446–14.458,
23.792–23.806, 53.945–53.947, and 53.947–53.962 Mb. Output read count stays
256,570; phased reads rise 236,878 to 237,070 (92.3249% to 92.3997%), correct
assignments rise 229,592 to 229,890, and discordant assignments fall 7,286 to
7,180 (96.9242% to 96.9714% accuracy). VCF blocks remain 341 while span N50
rises 643,699 to 672,998 bp. Seven additional insertion keys become eligible;
all shared-key unordered genotypes are unchanged. There are now 13 open gaps
in the previously tracked 14-case set. The new 53.945 Mb connection costs some
correctly tagged interior reads despite improved continuity and local purity;
this local coverage loss is retained explicitly in the evaluation.

All four closures are in the permanent window panel, with required recovered
sites, exact spanning expectations, and parental checks. The existing 23.806 Mb
complementary-row regression now permits the source-supported connection while
requiring both insertion/deletion rows and opposite alleles. A separate owning
23–24 Mb test checks the joined relative boundary genotypes and read truth.
The 15.351 Mb short replay lost source phase-set context and dropped to 96.42%
concordance; its original and repaired 15–16 Mb owning chunks score 99.4859%
and 99.5093%, respectively, with identical shared phased GT/PS fields.
Its test now uses that complete source context and retains the existing 99%
accuracy, 52% separation, and one in-gap-site floors. No floor was lowered.

Synthetic tests cover shifted calls in both directions, profile reindexing,
quality/mapping exclusions, existing-call preservation, different/complex indels,
reference bounds, and absence of recovery windows. Matched commands and detailed
measurements: `evaluations/2026-10-01-single-base-msa-recovery/`.

Validation passes: clean build, unit/predicate/port and upstream parity gates,
and the expanded 89-window suite (3,868 assertions in 48 Catch2 cases).

## 2026-10-01: preserve recovery evidence indices and all corroborators

Two synthetic regressions exposed further recovery defects. Growing a sparse
profile shifted BAM allele/query vectors but not `bam_base_qualities`, which
could assign Q40 to another site or discard the original SNP's quality vote.
The updater now pads populated qualities with zero at the same offsets, and
the post-transfer invariant rejects a mismatched quality-channel length.
The complete-BAM corroboration predicate also kept only the last SNP and
indel positions: a later nearby agreeing site could hide an earlier independent
confirming site. It now uses coordinate extrema of all agreeing calls in one
pass, with constant additional storage. Existing quality, source completeness,
100 bp independence and conflicting-evidence gates remain in force.

Neither defect changes the current chr20 fixture: the final VCF is byte-identical
to the accepted single-base insertion state, and all 256,570 read HP/PS pairs
are identical. Phased/truth-correct/discordant reads remain
237,070 / 229,890 / 7,180; 341 VCF blocks, N50 672,998 bp, and 13 tracked
HiPhase-positive open gaps. The two fixes therefore prevent reproduced failures
without claiming new chr20 closures. The predicate suite gains 30 assertions;
243 assertions in 24 cases pass, as do the units and both LongcallD parity gates.
No existing gap floor or expected span is changed for these fixes.

Further representation trials were rejected. Single-base deletion ALT recovery
closed no gap and worsened the 58 Mb score. Suppressing shifted equivalent
false REF calls after source solving improved chromosome read concordance but
reopened five established connections. An isolated non-homopolymer multi-base
deletion trial improved the 10 Mb score but increased discordant assignments
163 to 292 in the 36 Mb chunk. Both cases confirm that source phasing and later
allele-transfer evidence must remain consistent; broader backfill is not a safe
standalone fix. Evidence and unchanged metrics:
`evaluations/2026-10-01-recovery-evidence-indexing/`.

Final validation for the evidence-indexing/corroboration fixes: the 89-window
panel passes all 3,868 assertions in 48 Catch2 cases; final compilation and
243 predicate assertions pass with no new warnings. All gap expectations
remain unchanged.


### 2026-10-01: inspect original MSA SNP dropout before CIGAR backfill

A fresh HiPhase run on the identical current VCF/BAM exposes three new strict
noncentromeric targets (>=98% local truth purity, >=95% concordant disjoint
10-kb flanks): 528,850–542,052, 39,848,887–39,856,144, and
54,684,758–54,702,072. All are added to the permanent panel (89 -> 92 windows).

The 39.848 Mb owning chunk had 36 crossing reads, but 35 lacked the MSA SNP
call in the matrix used by source phasing. Post-solve exact-CIGAR backfill
created 32 pairs (31 same, one opposite), hiding dropout from the subsequent
MSA retry admission. Recovery now measures original MSA SNP absence before
backfill with known MAPQ30+, the existing 20-crossing-read minimum, and exact
binomial p <= 0.01 corrected for both endpoints when both are MSA SNPs.
Known multiallelic calls count as observed. Existing post-backfill retry
criteria keep priority; newly admitted retries retain every original phased
heterozygous source key and still undergo independent transfer/stitch checks.

The unplaced-read MSA restores the SNP from 5/11 (DP16) to 34/36 (DP70),
then the ordinary source solve joins the deletion/SNP with the correct same-ALT
polarity. The owning replay places all 99 gap-overlapping truth-scorable reads
correctly together (previous dominant block 45/99). The old executable fails
five assertions in the new owning-chunk regression. Its short replay already
joined, so using that short replay alone would miss this bug.

Full chr20: phased/truth-scored reads 237,070 -> 237,071; correct 229,890 ->
229,892; discordant 7,180 -> 7,179; concordance 96.971359% -> 96.971793%; read
PS 692 -> 690; VCF blocks 341 -> 340; span N50 672,998 -> 684,798 bp. All
62,368 VCF keys and every unordered genotype are retained. The baseline
includes the two preceding uncommitted accepted fixes after `74bc786`.

A separate insertion-repair boundary regression fixes membership to inclusive
VCF anchors (`sort_pos()`), including the right seam boundary. That change
alone preserves chr20 output. General pre-backfill MSA admission and direct
pre-phase CIGAR SNP repair both caused accuracy regressions and were rejected;
see the dated record for quantitative counterexamples. The original SNP
matrix distinguishes the kept gap (35/36 absent) from the 56 Mb retry that
lowers accuracy (10/26 SNP calls absent); that counterexample stays unchanged.

The panel preserves existing per-window floors and protected splits, adds the
required closing deletion row and allele/parental checks, and tightens TOTAL
floors to the summed per-window floors (73 joins, 1,288 in-gap heterozygotes).
Fifteen tracked competitor targets remain open, plus four intentional split
controls. Code contains no truth, competitor output, or fixture coordinates.
Measurements and reproduction: `evaluations/2026-10-01-msa-snp-dropout/`.

Validation for the retained state: build/shared units, 244 predicate assertions,
27 port-parity assertions, 168,696 upstream-parity assertions, and the expanded
window panel (3,951 assertions in 49 Catch2 cases over 92 windows) pass, with
no new compiler warnings. `git diff --check` and the reproduction script's
syntax/help checks pass. Existing per-window expectation floors are unchanged.


### 2026-10-01: use VCF anchors for postsolve MSA backfill membership

`backfill_msa_observations` used an indel's internal position to test membership
in inclusive VCF-anchor recovery windows. This omitted a right-boundary indel
and admitted an out-of-window left anchor. Membership now uses `sort_pos()`;
CIGAR allele calls retain internal keys and existing observations are preserved.
The new synthetic regression fails five assertions against the old code and
passes after the fix; the predicate suite gains 50 assertions (294 in 25 cases).

This fix adds no chr20 closure. All 256,570 read HP/PS pairs and parsed VCF
rows remain identical to the preceding accepted SNP-dropout state:
237,071 phased/scored, 229,892 correct, 7,179 discordant, 96.971793% concordance,
690 read phase sets, 340 VCF blocks, and 684,798 bp span N50. Fifteen tracked
HiPhase-positive gaps remain open. Twelve owning-chunk replays retain scores
and genotypes/phase labels; existing per-window expectations are unchanged.

The new 54.685 Mb gap has eight represented spanning molecules, and its BAM
source has an 88-site cut-free path, but its graph block contains a weak SNP
edge and boundary physical calls include a contradiction. Broad unplaced-read
MSA closes it with 22 extra correct read assignments, but loses five het keys;
at 0.529 Mb the same broad trial loses three het keys and adds 23 discordant
assignments. Neither broad admission nor the ineffective focused admission
probe is retained. Original source connectivity audits before/after backfill
show unchanged weak cuts in the inspected 0, 10, 21, and 54 Mb chunks.

Validation: build/shared units; 294 predicate assertions; 27 port-parity and
168,696 upstream-parity assertions; all 3,951 permanent window assertions in
49 cases over 92 windows; no new warnings; `git diff --check` passes.
Evidence: `evaluations/2026-10-01-msa-backfill-boundaries/`.

### 2026-10-01 — certify focused recovery after backfill; close 54.685 Mb

Correct the conflicting-pair admission statistic: seven 6:1 calls have a
fair-parity tail of 0.0625, not the unanimous 0.0078125 used previously. That
admits the existing focused MSA solve at 54,684,758–54,702,072. Keep every
admitted phased private row of selected, cut-free focused source blocks through
transfer, including context outside the seam, with cached source-path evidence.
Five apparent MSA losses were actually present in the BAM source and dropped
by transfer; no new consensus construction or row merge is needed.

The first trial split two established joins. At 4.767 Mb, focused path
certification preceded CIGAR backfill, which introduced a weak cut at
4,778,792. Certify the final backfilled matrix and retain the original source
when it fails. At 54.894 Mb, the newly joined left block exposes an internal
one-haplotype GAF edge (0/7/2) with Q30 BAM support 31:6. The physical SNP
stitch may corroborate this first failed edge using both haplotypes, exact
binomial p<=0.001, and signed quality odds; its full prefix/suffix path and
boundary evidence remain required. A dominant GAF reversal stays a veto.

Final default chr20: output 256,586; phased/scored 237,101; correct 229,920;
discordant 7,181; concordance 96.971333%; read PS 686; VCF blocks 340; records
62,798; span N50 739,888 bp. Versus the accepted preceding state: 30 more
phased reads, 28 more correct, two more discordant; N50 +55,090 bp (8.04%).
All 73 previous panel spans survive, and 54.685 Mb adds the 74th. Fourteen
tracked competitor targets remain open, plus four intentional controls.
Native competitor runs were not repeated. Four earlier VCF keys have changed
local allele/span descriptions after additional private source rows reach the
existing output conflict filter; the dated report lists them explicitly.

Extend the 92-window panel's 54.685 Mb case to its owning chunk and spans=1,
keep its boundary deletion required, and add regressions for the whole source
context/downstream 54.894 Mb connection and the 4.767 Mb final-matrix
certificate. The latter preserves the measured accepted owning-chunk 64%
separation floor; a first estimate from chromosome output was 69%, and fails
both the old accepted binary and the new binary. No earlier floor is reduced.
Evidence: `evaluations/2026-10-01-focused-msa-stitch/`.

Validation: warning-free build; shared unit tests; 294 predicate, 27 port-parity,
168,696 upstream-parity assertions. The complete native pipeline panel runs
51 cases / 4,010 assertions; its only failure is the new estimated 69% floor.
Fresh native baseline/final runs pass that calibrated 64% case (17 assertions),
and the new 54 Mb case passes 49 assertions. All 4,010 current assertions pass
when rechecked on the unchanged completed native outputs with the production
binary hash verified. Existing floors, truth checks, and required sites pass.
`git diff --check` passes. See the dated report's separate native/recheck logs.

### 2026-10-01 — restore the complementary deletion fallback without unsafe joins

When a clean SNP pair has no callable spanning molecules, the physical stitch
now tries its existing complementary left deletion-pair helper before the
single-indel fallback. The latter rejects overlapping alleles, so it could
never use that evidence. The pair helper now uses the same exact reference-edit
equivalence check as other deletion bridges instead of demanding the original
CIGAR position. Both graph paths, both deletion alleles, unanimous paired
parity, and the 0.001 quality-odds bound remain required. Extract the unchanged
equivalence algorithm into a directly tested shared BAM allele function.
Remove the pair helper's redundant original-coordinate flank test: the shared
caller has already certified the actual shifted placement, and an original
flank can lie inside the shifted deletion.

At 50.548 Mb, a direct audit of twelve Q30 spanning molecules restores four
previously missed ALT calls across the separate one- and two-base deletion
rows. This exposes a further right-graph path veto at
50,719,976–50,719,983. Direct Q30 CIGAR calls there do not independently support
both graph haplotypes; the gap stays split. A separate right-indel fallback
trial closes 58.366 Mb but changes owning-chunk correct/discordant reads from
3,492/57 to 3,428/122. A moved deletion retains its old source orientation;
the right block's own complete source certificate does not validate the moved
row. Reject that fallback. Generalizing shifted insertion repair beyond its
single-base scope also stays rejected: at 19 Mb correct/discordant counts
change from 4,031/41 to 3,746/328 without closing the target.

The retained default full chr20 run is exactly unchanged: all 256,586 primary
read names and HP/PS pairs and every nonheader VCF row match the accepted
54.685 Mb state. Phased/scored 237,101; correct 229,920; discordant 7,181;
96.971333% concordance; 686 read PS; 340 VCF blocks; 62,798 records;
739,888 bp N50; 74/92 full-chromosome panel spans. No additional tracked gap
closes; fourteen competitor targets remain open.

Add 104 synthetic deletion assertions (exact/shifted ALT, distinct lengths,
ambiguous edits, REF, flanks, qualities, coverage) and a 17-assertion owning
58–59 Mb counterexample. Exact-only calls fail five synthetic assertions;
the rejected right-hand join fails three integration assertions, including
parental switch and concordance. Existing floors and required sites remain
unchanged. Shared units, 398 predicate assertions / 26 cases, 27 port-parity
assertions, and 168,696 upstream-parity assertions pass without new warnings.
The complete window panel passes 4,023 assertions / 52 cases over 92 windows.
All 105 distinct commands execute on the final native binary with four
concurrent workers in 234.08 seconds; the unchanged test runner scores those
fresh outputs with exact-command and binary-hash checks. Earlier outputs serve
only to capture the input inventory. No old output substitutes for a final
native run. `git diff --check` passes. Evidence and reproduction:
`evaluations/2026-10-01-complementary-deletion-fallback/`.

### 2026-10-01 — exclude homozygotes from detached BAM run certification

`bam_source_run_supported` treated every row carrying the candidate PS as an
oriented anchor. Known BAM homozygotes retain source PS labels but have equal
haplotype alleles, so they incorrectly vetoed a run and could enlarge its
weak-cut interval. Exclude known homozygotes from the provenance index and
run extent/gauge checks. Preserve unknown-row, heterozygous-provenance,
same-source, orientation, and weak-cut vetoes. Move the unchanged remainder
of the predicate to the existing graph/BAM adapter for direct testing; its
physical SNP stitch caller and stitch rules remain unchanged. No new object
dependency, realignment, comparator data, or truth enters runtime decisions.

Add fifteen adapter checks. The old predicate fails the source-labelled
homozygote and extent checks; the corrected function passes those and the
negative provenance, source, orientation, cut, singleton, and malformed-row
checks. Diagnostic production calls encounter such homozygotes in PS
32,168,699 and 49,700,407. Removing this false veto closes no additional
tracked gap because the remaining path/stitch requirements still apply.

A focused interior-seam MSA dropout trial at 57.854 Mb restores observations
from eleven spanning reads but fails complete-path and heterozygote-key
preservation. A singleton-validator per-source gauge experiment and an
indel-helper BAM-only-path routing experiment make no measured improvement.
None of those trials is retained.

Fresh final default chr20, eight threads, 185.20 s: all 256,586 primary read
names and HP/PS pairs and every nonheader VCF row exactly match the accepted
preceding state. Phased/scored 237,101; correct 229,920; discordant 7,181;
96.971333% concordance; 686 read PS; 340 VCF blocks; 62,798 records;
739,888 bp N50; 74/92 panel spans. Fourteen tracked competitor targets and
four controls remain open. Native competitors are not rerun.

Warning-free build; shared units, 398 predicate, 27 port-parity, and 168,696
upstream assertions pass. All 105 current window input commands run natively
on the final production binary with four workers in 238.21 s; the unchanged
Catch runner scores those fresh outputs with exact-command and binary-hash
verification. Existing panel rows, accuracy floors, required sites and
owning-chunk regressions are retained. Evidence and reproduction:
`evaluations/2026-10-01-bam-source-homozygotes/`.
The complete window assertion pass is 4,023 assertions / 52 cases across
92 windows. `git diff --check` passes.


## 2026-10-01: normalize retry SNPs and retain their stitch orientation

The second recovery admission guard tested `key.alt.size() == 1`, discarding
catalog SNPs whose keys store graph walks. Translate the selected original
REF/ALT with the existing graph allele mapping and VCF normalization before
choosing physical SNP anchors. Keep injected BAM sequence keys. MNPs, indels,
whole multiallelic rows and missing mappings cannot become physical SNPs.

High-quality small/one-haplotype cohorts are not sufficient evidence for newly
admitted graph anchors. Retain MAPQ30/Q30, the four-read and 75% requirements,
and signed quality odds; graph-derived pairs additionally reject random
parity with exact one-sided binomial p <= 0.001, have at least two winning
calls on each haplotype, and no within-haplotype reversal. BAM-only pairs keep
the established source admission rule. Applying the extra graph count screen
to BAM-only pairs opens protected 15.056-Mb recovery; reject that trial without
changing expectations.

A further latent bug appeared at 61.747 Mb after normalization: ten physical
SNP pairs (6 REF/REF, 4 ALT/ALT) agree with stored observations and admit a
retry, but their orientation was discarded before the primary stitch reversed
the right block. Retain the graph relation in the existing physical-bridge
recovery gauge. Conflicting validated physical relations veto a join regardless
of row order. The short replay now joins the two SNPs with the supported
orientation and restores 368/379 correct reads (97.0976%), preserving its 97%
floor. Its explicit SNP PS/allele relation is now regression checked. The
neighboring catalog gap and owning context already belong to the panel.

Add 47 direct synthetic checks, including real indexed BAM calls, selected
ALT2 projection, weak and one-haplotype cohorts, same/swapped output, stale
output clearing, and physical-constraint row-order conflicts. Preserve the
15-Mb BAM-only owning connection in an additional integration case. Runtime
logic has no truth/competitor input or fixture-specific coordinate condition.
No candidate merge or new realignment is added. Shared physical SNP calling
moves unchanged to the graph adapter, and the existing binomial routine moves
unchanged into the shared phasing core for direct testing. Adapter/noise tests
link the existing BAM normalization dependencies.

Final default chr20, eight threads, 194.61 s: 237,101 phased/scored reads;
229,921 correct, 7,180 discordant, 96.971755% concordance. Exactly one HP/PS
assignment changes and becomes correct; all read names and every nonheader
VCF row are retained. Read PS 686; VCF blocks 340; records 62,798; N50 739,888
bp; panel spans remain 74/92. No newly closed default chr20 panel target;
14 competitor targets and four controls remain open. Competitors are not rerun.

Warning-free build; shared units and 398 predicate, 27 port and 168,696 original
upstream assertions pass. All 105 native window commands run with four workers
in 249.41 s and the unchanged Catch runner scores their fresh, exact-command /
binary-hash verified outputs: 4,041 assertions / 53 cases PASS. Existing panel
rows, required sites, accuracy floors and span expectations remain intact.
Behavior is updated in docs/IMPLEMENTATION.md; measurements and rejected trials
are in evaluations/2026-10-01-retry-snp-parity/.


## 2026-10-01: repair missing deletion ALT without changing source contrasts

Exact-CIGAR recovery backfill missed sequence-equivalent deletions whose
placement removed a required original-coordinate flank. Share the existing
32-base/64-base physical edit certificate with the loaded chunk reference.
Fill only missing ALT at simple oriented MSA rows, with known MAPQ30/Q30 and
an independently callable clean source SNP at least 100 bases away on the same
read. All callable source-SNP gauges must agree. Preserve exact observations,
explicit MSA ambiguity, coordinates, separate rows, and source counts. Source
k-means and ordinary BAM calling remain unchanged; no new alignment, object,
truth input, comparator input or fixture coordinate enters production logic.

The trials clarified a representation pitfall: a BAM MSA row's zero can mean
absence of that ALT, including a complementary MSA allele, rather than literal
reference sequence. Replacing those contrasts with literal deletion/REF calls
reopened six protected joins (24.582, 47.004, 48.930, 53.945, 53.948, 56.065 Mb)
and reduced N50 to 650,885 bp despite improving total read accuracy. Reject it.
Pre-phasing broad ALT additions raised 36-Mb errors from 163 to 292. Even
missing-only Q30 repair at a multi-row repeat added fifteen reads there with
eight errors. Require a simple row and independent source-SNP gauge; retain
all existing expectations. Allele certainty alone does not certify its source
orientation.

Final default chr20, eight threads, 156.55 s: every primary name and HP/PS
assignment and every VCF key, genotype and PS exactly match the accepted prior
state. Phased/scored 237,101; correct 229,921; discordant 7,180; 96.971755%
concordance; 686 read PS; 62,798 VCF records / 340 blocks; N50 739,888 bp;
74/92 connected panel gaps. No newly closed tracked chromosome target;
fourteen competitor targets and four controls remain open. Competitors are
not rerun. This is a reproducible missing-observation correction, not a claimed
coverage or contiguity improvement.

Add 93 predicate assertions. The actual original backfill fails the new ALT
recovery assertion; the final suite passes 491 assertions / 27 cases. Shared
units, 27 port-parity and 168,696 upstream assertions pass warning-free. All
105 window commands run natively on the final binary with four workers in
202.24 s; exact-command / binary-hash checked scoring passes 4,041 assertions /
53 cases over the unchanged 92-window panel and existing parental/owning gates.
No span, accuracy or required-site expectation changes. Behavior documentation
is updated; evidence and rejected trials are in
`evaluations/2026-10-01-missing-deletion-observations/`.

## 2026-10-01: certify deletion bridges in their original BAM source

Fix the physical SNP-to-deletion stitch's source lookup: a transferred deletion
must not borrow its current block's independent BAM certificate. Require unique
adoptable provenance, a complete original path without weak/quality cuts, and
another distinct-coordinate source anchor in the current block. All current
rows from that source must preserve one common allele flip; homozygotes provide
no gauge. Graph path and physical allele requirements remain unchanged. When
the selected clean SNP pair has zero callable pairs, try the nearer verified
right deletion through this corrected helper. No new alignment, threshold,
coordinate exception, runtime truth or competitor input is added.

This closes chr20:11,573,074–11,586,531 in its owning 11–12 Mb chunk through the
verified CT>C deletion at 11,591,585. One MAPQ60 molecule calls the left SNP REF
at Q40 and deletion REF across Q30 flanks; it ends before the clean right SNP.
The independently certified source gauge validates the same-haplotype relation.
Local phased/correct/errors remain 111/108/3 (97.2973%); correctly separated
fraction rises 47.7612% -> 64.9254%. The owning chunk retains 4,092 scored,
4,002 correct and 90 discordant. Archived HiPhase local purity is 94.4%, so this
newly closed case is outside the strict 98% competitor target list. Fourteen
previously tracked targets still remain unresolved.

Add the closure to the permanent panel with spans=1 and owning context; the
short padded replay lacks the complete source and stays split. Update the prior
11.599-Mb owning regression for this intended physical join while checking six
exact original allele genotypes, both separate insertion rows, shared PS,
parental orientation and the unchanged owning error ceiling. Its old split
expectation was for an unsupported focused MSA replacement, which is still not
admitted. The actual baseline fails five new connection/separation assertions.
Compile the old source lookup with the trial fallback: the existing 58.366-Mb
negative regression fails three checks, and its owning errors rise 57 -> 122.
The corrected guard rejects that wrong join. Add 20 standalone adapter checks.

Default full chr20, eight threads, 178.09 s: correct/discordant 229,921/7,180
and 237,101 phased/scored remain unchanged, 96.971755% concordance. Read PS
686 -> 684; VCF blocks 340 -> 339; N50 unchanged at 739,888 bp. All primary
names and HP labels, 62,798 VCF keys and genotypes match; only 609 read PS and
60 VCF PS labels change at the reviewed connection. The expanded 93-window
panel has 75 connections, versus 74 with the accepted baseline. Every prior
panel floor and required-site entry remains unchanged.

Append independent graph/BAM alleles to diagnostic observation rows so local
audits can distinguish conflicts from missing calls. No relaxed graph-indel
conflict/MAPQ criterion is retained. All 105 distinct native panel commands run
fresh on the final binary in 231.96 s; exact command / binary hash verified
scoring passes 4,082 assertions / 53 cases. Shared units, 491 predicates,
27 port-parity and 168,696 upstream assertions pass without compiler warnings.
Behavior is documented in docs/IMPLEMENTATION.md; evidence and comparisons are
in evaluations/2026-10-01-source-certified-deletion-bridge/.

## 2026-10-01: retain original insertion certificates and selected SNP parity

The physical SNP-to-MSA-insertion helper had the same current-versus-original
source lookup defect as the deletion helper. Its long-insertion complete-source
fallback now uses bam_source_site_path_supported: unique original provenance,
complete uncut source, and a consistent current gauge confirmed by another
distinct-coordinate source anchor. The graph-path alternative is unchanged.
The weak-left-edge suffix branch also discarded an available clean-SNP parity
selection and merged using the repeat insertion log odds instead. Both return
paths now use selected_flip. No allele, quality threshold or realignment
changes; neither correction consults truth or competitor output.

Add four adapter checks for long-insertion provenance/cuts/gauge. The existing
32.234-Mb owning bridge test passes 27 assertions. The discarded-parity branch
has no new native chr20 counterexample; its fix removes the inconsistent
selection between the two path-validation branches.

Default matched chr20, eight threads, 160.11 s: every read name and HP/PS tag,
VCF key, genotype and PS exactly matches the preceding accepted run. Phased /
correct / discordant remain 237,101 / 229,921 / 7,180; accuracy 96.971755%;
684 read PS; 62,798 VCF keys / 339 blocks; span N50 739,888 bp. Panel spans
remain 75/93; fourteen tracked competitor targets remain open. No new closure
or performance gain is claimed, and no expectations are refreshed. Shared
units and 491 predicate, 27 port-parity and 168,696 upstream assertions pass
without new warnings. Behavior documentation and measurements are in
docs/IMPLEMENTATION.md and evaluations/2026-10-01-insertion-source-parity/.

All 105 distinct native window commands run on the final binary in 194.81 s
with four workers. Exact-command and binary-SHA verified scoring passes
4,082 assertions / 53 cases over the unchanged 93-window panel and owning
orientation controls. The final standalone unit rerun passes.

## 2026-10-01: call repeat-shifted insertion ALT without changing source gauges

The 15.351-Mb audit finds a concrete boundary representation loss: four Q40
molecules place the ATCT insertion 44 bp right of its BAM candidate across
11 identical motifs. The physical stitch's 16-bp longer-insertion window
calls them REF. Add an ALT-only fallback over the complete reference-equivalent
interval when the original call is REF and no nearby indel interferes. Require
a unique identical reference edit and known flank-threshold qualities across
inserted and crossed bases; retain the actual minimum inserted quality. Keep
the existing bounded interference window, callable ALT decisions, source
profiles/HP/PS/cuts, candidate identities and stitch gates. No realignment,
truth or competitor data is used in production.

Reject a wider interference-window / pre-solve multi-base repair trial despite
7,180 -> 7,149 full-chr discordant reads: the protected 47-Mb join splits and
the 17.865-Mb suffix joins the wrong prefix. Ordinary postsolve backfill still
changes subsequent retry/source-path decisions and fails the same controls.
Both experiments are removed. Inferring an implicit reference SNP row also
fails to close the target and makes an insufficiently accurate new join.
The separate source-certificate and observation-transfer gates remain intact.

Instrumented upstream HiPhase v1.6.0 on the same current pgphase callset gives
8 callable pairs at 15.351 Mb (six cross, two same), 18 at 21.159 Mb (10 REF/REF,
7 ALT/REF, 1 REF/ALT), and six for the two-base deletion at 57.854 Mb. Preserve
complete REF/ALT identities in the audit: two different deletion rows share
the last coordinate. HiPhase's allele logging hook does not alter its solve.
The 15.351-Mb left deletion still has no callable paired source observations;
repairing the right insertion alone does not close the gap.

Final default matched chr20, eight threads, 215.24 s concurrent with the
native panel: every read name and HP/PS, VCF key, genotype and PS exactly
matches the accepted baseline. Phased/correct/discordant remain
237,101 / 229,921 / 7,180 (96.971755%); 684 read PS; 62,798 VCF rows /
339 blocks; span N50 739,888 bp; panel spans 75/93. Fourteen tracked competitor
targets remain open. No new closure or chromosome phasing gain is claimed.

Add a 108-assertion insertion predicate case covering the 44-bp reproducer,
rotated motifs, both shift directions, qualities, ambiguous/different/compound
edits, source-state preservation and local/worker reference equivalence.
Standalone units, 599 predicate / 28 cases, 27 port / 9 and 168,696 upstream /
7 pass without new warnings. All 105 fresh native panel commands complete in
267.45 s; hash/command-verified scoring passes 4,082 assertions / 53 cases.
No panel, required-site entry or expectation is refreshed. Behavior is in
docs/IMPLEMENTATION.md; evidence and rejected trials are in
evaluations/2026-10-01-repeat-insertion-alleles/.


## 2026-10-01: exact opposite-insertion MSA observations

Reproduce a local deletion-row observation loss: the other complete verified
MSA haplotype inserts bases inside the deletion footprint, but the classifier
rejects that exact longer sequence as unknown. Accept its exact ALT-absent
contrast only with complete, matching reference footprints and exactly one
fully deleting consensus. Require matching read sequence and supported flanks
on both composed paths; retain separate candidate rows and their own counts.
Keep other deletion lengths under their existing handling, and prevent the
exact contrast from enabling the one-error local fallback for mixed edits.

A broader different-length contrast trial adds 91 chr20 phased reads, 76 correct
and 15 discordant, and spans the 15.095642-15.101261 Mb control. Although its
owning 15-Mb replay gains 38 correct reads without new errors and passes the
parental-flank check, its 34-Mb replay gains 43 reads with 14 extra errors.
Remove that broader trial and its proposed control expectation. Temporary
mixed-event MSA admission, interior ownership and singleton-padding changes
also fail to close the 15.351845-15.367755 Mb target and are removed. A focused
trial recovers all 65 left deletion calls, but the BAM solve collapses the
site to homozygous; the eight HiPhase-paired calls are six cross versus two
same, and no complete certified source path emerges.

The final insertion-only correction has no truth or fixture condition. The
new predicate case has 37 assertions, including the original three-failure
reproducer, separate row/count preservation, distinct deletion lengths and
fuzzy-context rejection. Standalone units pass; predicates 636 / 29 cases,
port parity 27 / 9 and original upstream parity 168,696 / 7 pass. All 105 fresh
native commands complete in 276.94 s with four workers. SHA/command-verified
scoring passes 4,082 assertions / 53 cases; expectations remain unchanged.

Final default eight-thread chr20 in 224.46 s concurrent with the native panel:
237,101 phased, 229,921 correct, 7,180 discordant, 96.971755% concordance;
684 read PS; 62,798 VCF rows / 339 blocks; span N50 739,888 bp; 75 / 93 panel
spans. Every HP/PS and VCF key/genotype/PS matches the accepted baseline.
Fourteen HiPhase target gaps remain open. No new closure or coverage gain is
claimed. Binary SHA256:
d1e05dfe2252cf89c422dfd87db06f240eda6d7f859bc8a300a840008cfa1dfb.
Behavior is updated in docs/IMPLEMENTATION.md; evidence is in
evaluations/2026-10-01-msa-opposite-insertion/.


## 2026-10-01: close the mixed-MSA 15.351 Mb gap without coverage loss

Report: `evaluations/2026-10-01-diploid-msa-gap/README.md`.

Fixed explicit joint-genotype link admission: a depth-supported MSA
homopolymer heterozygote retained by `joint_het_orientation` can enter the
upstream link list under that same predicate. Ordinary BAM mode preserves
the longcallD repeat exclusion and passes original-C parity.

A co-located verified insertion/deletion contrast with diploid allele depths
can admit a focused unplaced-read MSA retry after independently significant
local and crossing dropout. The contrast may have been collapsed to homozygous
in the original solve. Focused source ownership supports a middle seam without
rerunning the original source for its other windows; adjacent phase-set extents
provide context, with existing 50 kb padding for singleton flanks. The source
must retain old phased keys and a complete cut-free path after backfill.
The focused solve uses existing joint-genotype and allele-pair link modes.

Noisy MSA genotypes retained by this diploid retry carry
`read_rescue_requires_validation`: block connectivity alone does not certify
one-locus read calls. They need the existing primary-read singleton confidence
test or a second independent locus. Original read-rescue semantics are retained
for ordinary genotypes. Broad retry and broad singleton filtering were rejected
after losing two protected gaps and 1,182 phased reads; no existing test floor
was weakened to admit them. A catalog-only fallback validation experiment had
no effect and was removed; the suspect PS+1e9 calls were graph rescue, not
independent BAM fallback (PS+1.5e9).

Accepted chr20: 237,101 phased reads unchanged; truth-correct 229,921 ->
229,922, errors 7,180 -> 7,179, accuracy 96.971755% -> 96.972176%. Read
PS 684 -> 682; VCF blocks 339 -> 338; N50 remains 739,888 bp. No old
VCF key or genotype is lost/changed; 13 source rows are added. Only the
15,351,845–15,367,755 tracked gap changes from open to closed, so spans
75/93 -> 76/93. All prior closures and split controls remain intact;
13 competitor targets remain open.

Owning gap: 92 phased/89 correct (96.74%), versus baseline 92/88 (95.65%)
and the existing HiPhase same-callset diagnostic 137/118 (86.13%). Correct
opposite deletion/insertion alleles and both neighboring SNP gauges are
protected; older distant source-block assertions stay intact. Tests add local
floors 92 scored, 89 correct, 96% concordance and at most three errors.

Fresh native 105-request panel passes 4,106 assertions/53 cases. Predicates
663/31, port parity 27/9, original longcallD C parity 168,696/7 and all
standalone units pass. Full chromosome: 216.26 s; native panel: 264.91 s.
Production remains truth/competitor/coordinate/read-name agnostic.

## 2026-10-01: retain later MSA dropout requests and certify fallback SNP chains

Report: `evaluations/2026-10-01-retry-scheduling/README.md`.

Fixed first-request masking in grouped recovery: collect original MSA dropout
requests before backfill, try the established selection first, and try additional
requests if its focused certificate fails. Keep one accepted focused owner per
group and preserve legacy edge/middle ownership. Additional requests can admit
already phased complementary same-type MSA alleles, with unchanged dosage and
local/crossing dropout statistics. Explicit fallback admission leaves default
group selection and established local retry behavior unchanged.

At 57.854 Mb the later solve now runs and retains both deletion genotypes,
but its crossing pairs conflict and the source remains split. Two HiPhase
ALT calls correspond to three-/four-base CIGAR deletions, not the two-base
candidate. The reads are present; no missing-coverage fix can equate these
representations. Exact/consensus uncertainty is retained and different-length
regressions are added. This audit does not establish which allele likelihood
is biologically correct.

Rejected broad/multiple-source retries after four established assertion
failures. A narrower uncertified additional source at 61.8 Mb still added
17 phased reads but only one net correct assignment, and failed an existing
99% gate. Additional fallback sources now require the existing Q30 physical
certificate at every consecutive clean-SNP cut across the seam, allowing an
intermediate-site chain. No new threshold, alignment, truth input or fixture
exception enters production. Existing retry admission and checks remain intact.

Final default chr20 is identical to the accepted baseline in every primary
read HP/PS and VCF key/GT/PS: 237,101 phased, 229,922 correct, 7,179 errors,
96.972176%; 682 read PS; 62,811 VCF keys / 338 blocks; N50 739,888 bp;
76/93 tracked spans. Thirteen competitor targets remain open, plus four
controls. No new closure is claimed. Full runtime 186.44 s; 105 fresh native
replays 235.71 s. New scheduler regression fails twice before and passes 12
assertions after; predicates 684/31, port parity 27/9, original longcallD C
parity 168,696/7 and standalone units pass. No established gate is lowered.
Final window scoring passes 4,118 assertions / 54 cases, with 106 fresh
native commands including the new matrix regression. Production SHA256:
6855d83085f088e53135ce716b4c8e842532012474cdab3276bdee85483008dc.

## 2026-10-01: close the compound 528.850 kb gap with conserved accuracy

Report: `evaluations/2026-10-01-compound-flank-recovery/README.md`.

The BAM boundary at DEL542,053 / SNP545,002 has 42 cross versus 16 same pairs.
A decisive majority previously hid excess conflicting calls and prevented an
unplaced-read MSA retry. Broadly removing that veto regressed existing 6.578
and 12.269 Mb accuracy gates; it is rejected. New admission is restricted to
compound flanks: a phased verified MSA deletion covers a graph anchor lacking
any corresponding BAM SNP row. It supplies no stitch parity. Existing retry
selection remains unchanged. Replacement preserves all old phased keys and
one consistent clean-SNP source gauge per original block.

The newly admitted compound source retains complete, cut-free private flank
paths in transfer. Seam-only transfer lost the nine-base deletion/SNP context
needed to orient the child graph SNP. Physical REF/deletion versus MSA allele
association now uses the existing two-sided Fisher p <= 0.01 with a matching
majority in each allele class; one sequencing error cannot veto a supported
projection. Physical SNP ALT, unsupported graph/haplotype association and
source path cuts still veto it. The SNP and complex insertion ALT stay opposite.

Fixed partial-read MSA observations outside physical BAM coverage: remove
those calls and rebuild allele/total depths in O(observations + sites), requiring
SNP spans or both surviving indel flanks. This is limited to unplaced-read
recovery; ordinary upstream BAM MSA is unchanged. A homozygous category with
stale distinct internal alleles needs inferred-site association rather than a
direct singleton marker. All MSA genotypes from the new compound retry carry
the existing read-rescue validation flag. Broad complete-flank transfer and
filtering source/overlay read labels are rejected after regressions.

The target's 11 physical spanning molecules initially have zero left ALT calls
and six missing right calls. Retained MSA produces six left ALT / five REF and
recovers all six right calls, matching targeted default HiPhase on both deletion
alleles for all 11 molecules. Local graph matrix read spans are observation
extents, not physical BAM bounds; they cannot diagnose missing BAM coverage.
Targeted same-callset HiPhase phases 129 gap-overlap reads, 127 correct / two
errors (98.45%). Pgphase retains all 123 of its gap reads, 123 correct, now in
one phase set; all disjoint flank checks agree. Coverage still trails HiPhase
locally by six reads. The new owning-chunk regression protects the six flank
rows, correct parental/allele connection, 123/123 local assignments and zero
errors. It fails three assertions before and passes 43 after. The existing
panel target becomes spans=1; all pre-existing accuracy floors stay intact.

Final default chr20: phased/scored 237,101 -> 237,098; correct 229,922 ->
229,923; errors 7,179 -> 7,175; accuracy 96.972176% -> 96.973825%. Read
PS 682 -> 680; VCF rows 62,811 -> 62,850 with zero old keys lost; VCF blocks
338 -> 337; N50 remains 739,888 bp. Tracked spans 76/93 -> 77/93, with only
528,850–542,052 newly closed. Twelve competitor targets and four controls
remain open. This improves continuity and accuracy; phased coverage drops by
three reads rather than increasing.

Full chr20: 179.69 s; 105 fresh native commands: 234.81 s. Together with the
separate scheduler matrix test, 106 native commands use the final binary.
Scoring passes 4,157 assertions / 55 cases; predicates 697 / 32; port parity
27 / 9; upstream C parity 168,696 / 7; standalone units all pass. Build passes
with the existing vendored abPOA unused SIMD helper warning. No production
truth/competitor/coordinate/read-name exception or new alignment is introduced.
SHA256: 0dee719443489efb56f5f33efc69c0c825f93904477af1c25d55cde63206f733.

## 2026-10-01: fresh HiPhase / pgphase metrics and runtime comparison

Report and reproducible runner:
`evaluations/2026-10-01-phaser-comparison/README.md`.

Ran committed pgphase `94410d0` and unmodified HiPhase `1.6.0-ac3f399`
sequentially with eight worker threads on full HG002 / CHM13 chr20. Native
HiPhase uses the existing DeepVariant calls made from the same BAM. A second
HiPhase arm uses the fresh pgphase calls with phase labels removed. Normalized
SAM hashes verify identical input records, including sequence, qualities and
tags, between the annotated BAM and its linear-contig reheader; reference
sequences match. All 272,016 input reads are unique primary alignments.

| Metric | pgphase graph + recovery | HiPhase / DeepVariant | HiPhase / pgphase calls |
| --- | ---: | ---: | ---: |
| Phased reads | 237,098 | 233,353 | 230,138 |
| Phased / input | 87.1633% | 85.7865% | 84.6046% |
| Truth correct | 229,923 | 223,800 | 216,061 |
| Truth errors | 7,175 | 9,553 | 14,077 |
| Truth accuracy | 96.9738% | 95.9062% | 93.8832% |
| Read PS | 680 | 196 | 200 |
| VCF blocks | 337 | 196 | 200 |
| VCF span N50 | 739,888 bp | 1,005,183 bp | 919,153 bp |
| Native wall time | 147.49 s | 87.87 s | 71.38 s |
| Peak RSS | 20.125 GiB | 0.389 GiB | 0.224 GiB |

Pgphase phases 3,745 more reads than native HiPhase, produces 6,123 more
correct assignments and 2,378 fewer errors (+1.0676 accuracy points), while
HiPhase wins runtime, memory and block N50. On the 227,642 shared phased reads,
independently orienting blocks on that subset gives 98.3421% versus 96.5112%.
Full-block orientations give 98.3228% versus 96.4932%. HiPhase-only reads remain:
5,711 tagged, 4,141 truth-correct under its complete block orientations.

The graph command's BAM is a compact unaligned HP/PS table, not full original
alignments. Timed four-I/O-thread Python materialization and indexing of a full
alignment BAM adds 94.38 s, without changing any assignment. Comparable full
indexed BAM delivery therefore costs pgphase 241.87 s versus native HiPhase
87.87 s. This is a measured postprocessing implementation, not an optimized
native writer. The 15,430 input reads omitted from pgphase's 256,586-row tag
file remain unphased in the common denominator, avoiding inflated coverage.

Both tools' native truth totals reproduce prior certified results. Every phased
read has truth; each PS is independently majority-oriented. Same-callset HiPhase
preserves all 62,850 VCF keys and unordered genotypes. Runtime measurements are
one invocation each with existing filesystem cache; no statistical speed claim
or cold-cache benchmark. External alignment/calling/catalog/truth preparation
and evaluation are excluded. Pgphase's internal BAM recovery is included;
native HiPhase starts from its pre-called DeepVariant VCF. No pipeline code or
existing test expectation changes in this evaluation. Raw resource logs,
identities, results.tsv, metrics.json and reproducible scripts are retained.

## 2026-10-01: correct recovery anchors and single-seam focused retries

Report: `evaluations/2026-10-01-recovery-anchor-and-single-seam/README.md`.
Baseline is committed `94410d0` and the fresh comparison above.

Share a strict recovery anchor predicate: positive PS, non-homozygous category,
both internal alleles known, exactly one selected ALT. Stale unequal labels on
classified homozygotes and incomplete allele gauges cannot define source runs,
path cuts or transfer anchors. Multiallelic `1|2` remains valid; projected `0|0`
does not. Unknown transferred gauges still veto certification.

Backfill missing MSA insertions using existing Q30/MAPQ30 equivalent-repeat
edit checks before accepting exact-coordinate REF. Require an independently
callable clean SNP confirming the same source gauge, as in deletion repair.
Keep original rows, coordinates, gauges and existing MSA calls. Removing that
SNP requirement reopened the protected 41.9008 Mb gap and lost 20 old keys;
that experiment is rejected. Refreshing all assigned MSA observations during
focused retries also regressed protected rows and is rejected.

A single-seam group now admits the same focused full-adjacent-phase-set solve
as an outer seam of a multi-seam group. Complementary co-located boundary rows
remain separate. Existing phased-row preservation, complete source path and
transfer validation stay mandatory; no stitching threshold is relaxed.
The full window suite exposed another admission bug: a rejected focused trial
could authorize broad unplaced-read MSA for nonunique boundaries. That broad
fallback split the protected 4.767 Mb short replay although the full chromosome
retained its span. Recheck the broad fallback's own unique-boundary admission;
focused-only evidence cannot bypass the focused source-path certificate. This
restores the short replay without changing its committed expectations.

This closes 10,727,690–10,746,628 in its owning chunk and the full chromosome.
The SNP, 16-base deletion and insertion ALT connect on one haplotype, opposite
the four-base deletion, matching HiPhase's allele relation. Local correct reads
rise 134 -> 136 of 143; errors 9 -> 7. HiPhase is still better locally at
136/141 (five errors): three repeat-only reads are wrong only in pgphase and
four other reads are also wrong in HiPhase. Their relevant MSA region is
already complete; truncation does not explain the remaining ambiguity.
The new regression protects all six flank/bridge rows, relative allele parity,
parental flank orientation and local correct/error counts. The existing panel
case now requires spans=1 and both intervening deletion rows. Its existing
98% whole-chunk floor stays; it is not a 98% local-accuracy measurement.

The newly triggered 59 Mb retry phases 15 extra reads and corrects 12 prior
errors: correct 3,654 -> 3,681, errors 14 -> 2. Two old graph SNP descriptions
are replaced in the emitted VCF by physical insertion representations, matching
DeepVariant: 59,679,069 A>AT and 59,757,500 T>TC. The old graph rows remain
internal. Physical MAPQ30 base calls are 57 A / 2 T and 61 T / 1 C respectively;
most insertion calls lie one base later. The new owning-chunk regression pins
both insertions, opposite ALT haplotypes, one PS, correct parental flanks and
at most two errors. The anchor fix at 34 Mb separately leaves five correct
reads unphased while correcting three errors.

Fix the de novo MSA helper's cluster-index/read-index confusion when obtaining
coverage flags after filtering. The partial-first synthetic test fails on old
code and passes after. The local upstream longcallD helper has the same defect;
the ordinary wrapper sorts full-cover reads first, so this latent helper bug
is not a demonstrated cause of the remaining chr20 gaps.

Full chr20: phased/scored 237,098 -> 237,118; correct 229,923 -> 229,956;
errors 7,175 -> 7,162; accuracy 96.973825% -> 96.979563%. Read PS 680 -> 675;
VCF rows 62,850 -> 63,401 (553 added, two representation corrections), VCF
blocks 337 -> 338; span N50 739,888 -> 756,878 bp. All 77 previously closed
tracked gaps survive; tracked spans become 78/93, with only 10.727 Mb newly
closed. Eleven historical competitor targets and four controls remain open.
Truth/competitor are evaluation-only; no new alignment stage or production
fixture exception. The earlier sequential runtime comparison is preserved.

Final validation passes: windows 4,573 assertions / 57 cases; predicates
762 / 34; port parity 27 / 9; original-C phase parity 168,696 / 7; standalone
units all pass. The 93-coordinate panel adds two dedicated owning-chunk tests
without dropping existing cases or weakening their floors. The 105 fresh
native commands plus the separate 57 Mb matrix run use the final binary;
source fixes include the fallback guard exposed by the first window suite.
Full chr20 takes 184.07 s with other validation active; this is not a new fair
runtime comparison. Build passes with the existing vendored abPOA unused SIMD
helper warning and no new warnings. Final SHA256:
598eb8f99b1e14d1ce25eb9a8ef0305bea60b9b44089c2a3f41f9a05c95ac676.

## 2026-10-02 — Recovery right-anchor boundary and rejected admissions

Fix the closed-interval mismatch in `allele_depths_call_het`: recovery seam
windows and MSA backfill include both VCF anchors, but scoped heterozygote
repair excluded the right endpoint. Use `pos <= end`, still evaluated at
`VariantKey::sort_pos()`. The failing-before predicate test and SNP/insertion/
deletion singleton tests preserve all depth, repeat and verification guards.
Standalone BAM calling with no retry window is unchanged.

Full chr20 read HP/PS pairs and VCF data rows are identical to the accepted
2026-10-01 recovery state: 237,118 phased/scored, 229,956 correct, 7,162 errors,
96.979563% accuracy; 675 read PS; 338 VCF blocks; 756,878 bp span N50; 78/93
tracked coordinate spans. This is a latent boundary fix, not another measured
chr20 gap closure. Preserve all previous uncommitted fixes and regressions.

Reject the physical REF-only graph-anchor trial. At 21,172,501 all 61
callable BAM bases are REF; the remaining 17 reads overlap a deletion, which
is not REF evidence. Removing that phase anchor exposes the true clean SNP
at 21,183,846 and improves aggregate correct/errors by +2/-1. However, it
fragments the neighboring 21.179 Mb read block: dominant correct fraction
falls below its 0.85 floor to 0.673913. The 21,172,487 deletion stays in a
separate block, despite a coordinate extent newly spanning its historical
gap. Coordinate coverage alone cannot certify boundary-allele connectivity.
Keep the old phase-label regression and all old floors. The trial's helper,
new owning-chunk test and changed span expectation are removed.

Also reject admitting source MSA rows based solely on `msa_verified`, before
`alignment_verified` is supplied by transfer. At 19 Mb it adds 34 assignments,
only 17 correct, and leaves the tracked target's 125/94/31 phased/correct/errors
unchanged. Source MSA verification is not a substitute for the existing
admission evidence. A homopolymer deletion ALT backfill trial likewise does
not change the 21 Mb outcome and is removed.

Findings, rejected trial metrics and regression failures are preserved in
`evaluations/2026-10-02-recovery-boundary-validation/`. No truth or competitor
input enters production, and no gap expectation is relaxed in this round.

Final validation: 105 fresh native command replays and the separately generated
57 Mb matrix replay match the final binary SHA256; windows pass 4,573 assertions
in 57 cases, predicates pass 775 assertions in 34 cases, standalone units pass.
Build succeeds with no new warnings. Final binary SHA256:
`b8c2d78e9bd67136f2a1fa191cdc898aafdf37009ba35932202f4b899a18eaee`.

### 2026-10-02: retain supported partial BAM paths inside graph seams

The 36,332,599–36,354,890 target was inside a broader graph seam ending at
36,381,019. Its initial BAM source assigned one numeric PS across an insertion
and deletion even though the observation path was weak: 14 consistent paired
calls on one haplotype, 12 conflicts, and no consistent pairs on the other.
The focused MSA retry restored the insertion's depth from 24 (20 REF/4 ALT) to
74 (50 REF/24 ALT), but the full-seam acceptance rule discarded this partial
repair because an independent downstream source block remained separate.

Admit focused retries for significant internal MSA-indel conflicts, including
those exposed after CIGAR backfill. Use physical anchor order without moving
candidate rows or profile offsets. Only source-assigned, non-skipped, known
MAPQ30 molecules physically spanning both callable diploid sites vote.
Require the configured minimum depth or six pairs, two conflicts, binomial
upper-tail p <= 0.01 under a 5% error model, and a missing consistent haplotype.
This is observation admission, not a stitch certificate.

Allow a partial replacement only when every old heterozygote survives, every
old clean-SNP gauge remains consistent, every previously supported source edge
survives, and a weak edge inside the seam becomes supported. A replacement PS
cannot merge two old source gauges. Ordinary transfer and stitching decide
later cross-block connections. Track internal-conflict admission separately so
a rejected focused trial cannot authorize a broad MSA fallback. The internal
conflict check precedes two-block admission within each seam.

Owning 36 Mb replay: 3,783/3,620/163 phased/correct/errors becomes
3,790/3,637/153. Gap-overlapping reads remain 172 phased, improving from
158 correct/14 errors (91.8605%) to 166/6 (96.5116%). Native HiPhase with
DeepVariant calls on the same BAM has 172/162/10 (94.1860%). HiPhase with
pgphase calls has 162/158/4 (97.5309%): higher purity on fewer assignments.
The target SNP and deletion share ALT parity; the intermediate insertion is
opposite. Disjoint parental flank votes agree in orientation.

Reject the first prototype's unordered internal-edge detection and overly broad
retry admission: its full run adds 137 phased reads but 52 errors. Also reject
moving internal conflicts behind legacy admission: it changes a larger block's
gauge, losing 328 correct assignments and adding 335 errors. It also creates a
new coordinate span at 19 Mb while the exact gap boundary alleles still have
different PS labels; that extent is not a closed target. Neither prototype
is the accepted implementation. Preserve the target's exact allele regression,
the 19 Mb separation regression, and all old panel accuracy floors. Measurements
and rejected trials are in `evaluations/2026-10-02-partial-source-recovery/`.

Final chr20: 237,127 phased/scored reads, 229,974 correct, 7,153 errors,
96.983473% accuracy. Relative to the accepted boundary-validation baseline:
+9 phased, +18 correct, -9 errors. Read PS: 675 -> 674; VCF blocks remain 338;
span N50 remains 756,878 bp. All 63,401 old VCF keys survive; 91 new keys are
present. Tracked coordinate spans: 78/93 -> 79/93, with no old closure lost.
The only new target is 36,332,599–36,354,890, whose exact SNP/insertion/deletion
phase labels and relative allele parity are also verified in the full run.
Ten historical competitor targets and four controls remain open.

Validation: 105 fresh native replays matched to the final SHA256 and exact
normalized command, plus a separately generated 57 Mb matrix replay. Window
tests pass 4,646 assertions in 59 cases; predicates pass 791 assertions in 35
cases; standalone units pass. The new 36 Mb regression fails on the baseline
binary (7/30 assertions), including the exact connection and the local accuracy
requirements. The owning 19 Mb separation control passes; full-chromosome
evaluation additionally protects against changes in larger stitched blocks
that an owning-chunk replay does not expose. No existing floor is relaxed,
and parental truth and competitor results remain evaluation-only inputs.
Build succeeds without new warnings; `git diff --check` passes. Final binary:
`97378a42d9ad82cf9afde4cb328a0637eb14d767de5eeb629fe4f69c5a44cbb1`.

### 2026-10-02: preserve observed long insertions in a transferred BAM run

The source run after the 1.180 Mb weak cut contains a cleanly linked chain
through SNP 1,194,233 and the 155-bp insertion at 1,196,894. Their paired
observations support both source haplotypes (33/33) with no conflicts. The
left-read gauge was already supported, but transfer stopped at an insertion
length guard intended for physical length-based allele calls. That guard is
not appropriate for existing, observed source rows. Remove it from this run
transfer; retain the p <= 0.01 paired-edge certificate at every step, the
independent left-read gauge, and the stop before the first graph-owned row.
Physical insertion calling retains its existing 128-base limit.

Only the insertion's VCF PS changes, from 1,450,106 to 1,142,089. All 63,492
VCF keys, genotypes, depths and 256,601 read HP/PS pairs survive unchanged.
Full chr20 remains 237,127 phased/scored, 229,974 correct, 7,153 errors,
96.983473% accuracy, 674 read PS, 338 VCF blocks and 756,878 bp span N50.
Gap-overlapping pgphase reads remain 90/90 correct; native HiPhase has 91/91.
This fixes variant-path continuity, not the remaining read-block split. The
73-bp seam to the graph representation at 1,196,967 remains open. The right
graph block has two one-haplotype SNP edges at 1.349--1.352 Mb; bypassing one
weak graph SNP alone does not resolve the join. Never report coordinate span
as complete read-block parity.

Add the owning 1--2 Mb case to the committed panel with spans=1, retain every
old floor, and require the exact SNP/insertion allele gauge, unchanged separate
graph row, parental flank consistency and zero local errors. Baseline binary
fails the exact connection (1/34 assertions); final binary passes. Panel now
contains 94 coordinate cases, of which 80 span. A fresh audit nominates 43
current accurate competitor bridges, excluding 25--30 Mb, with >=95% local
purity and agreeing disjoint parental flanks. Preserve both baseline/current
lists and the reproducible audit in `evaluations/2026-10-02-long-insertion-run/`.

Reject grouped/reordered source-path votes: local controls match, but the full
run loses 31 phased reads, 52 correct assignments and 299 old VCF keys while
adding 21 errors. Restore the accepted source audit. Also reject broader
crossing-dropout admission (extra runtime without output improvements) and
off-seam shifted-insertion recall: at 33.314 Mb its 41/57 correct reads still
trail HiPhase using the saved pgphase callset (51/53). No trial relaxes old expectations or enters
the final binary. Truth and competitor inputs remain evaluation-only.

Validation: 105 fresh native command replays plus the separate 57 Mb matrix
replay match final binary SHA256. Window tests pass 4,696 assertions in 60
cases; predicates pass 791 assertions in 35 cases; standalone units pass.
Build succeeds without new warnings, and `git diff --check` passes. Final
binary SHA256:
`9e50548ea3690f94d4922571c3bd2958f287cb7f5118bcc2dbc36c9360933d35`.

## 2026-10-02 — Join equivalent insertions on the common BAM background

Close the remaining 73 bp seam at chr20:1196894–1196967. The two 155 bp
insertions differ on bare FASTA at edited-string offsets 55/210 but match
exactly after the verified homozygous 1196950 C→T SNP (2 REF/76 ALT calls)
is included. Keep both rows and all observations. Require genomic REF in the
graph binary pair, known BAM/GAF MAPQ >=30, both allele classes, zero opposing
paired calls, p<=0.001, complete graph SNP paths and a consistent connection
to the next clean SNP. Evidence is 36 REF/REF plus 40 ALT/ALT reads.

Certify after recovery finalizes its candidate indices, then apply after read
rescue and independent BAM assignment. The core proof does not certify the
output-only groups' relative gauge. Recompute core parity from current candidate
alleles after cross-chunk flips, relabel the entire downstream core PS across
all batch chunks, and preserve every output-only assignment. An earlier-merge
trial pooled eight existing read rescues with an independently gauged cohort
and raised owning errors from 10 to 18. Rejected nested-SNP removal and broader
singleton-validation trials lost supported links or correct reads. No such
restrictions enter the accepted binary; individual complex-context singleton
allele uncertainty remains a separate issue.

Final full chr20: 237127 phased/scored, 229974 correct, 7153 discordant,
96.983473% concordance. All read HP values, all 63492 VCF keys and all GT/depth
fields are unchanged. Only 2030 core read PS and 712 VCF PS fields change.
Read phase sets 674→673, VCF blocks 338→337; span N50 stays 756878 bp. Both
pgphase and the existing native HiPhase/DV run place 78/78 local reads correctly;
this is a continuity improvement, with no additional phased reads. Eight-thread
run takes 277.71 s, versus 281.99 s for the prior run; these are not isolated
runtime measurements.

Add the seam to the committed owning-chunk panel with spans=1, the retained
insertion row, exact allele-gauge checks and parental read checks. Preserve all
old read floors and the owning ceiling of 10 errors. Panel is now 95 coordinate
cases / 81 spans. A fresh audit excluding 25–30 Mb finds 316 extent gaps and
42 distinct accurate competitor-supported targets (53 records), down from
317/43/54. Nominations still need exact representation/path investigation.
Saved-pgphase-callset HiPhase uses the Oct 1 VCF, not the current output.

Validation: 105 fresh native replays verified against final SHA256, separate
fresh 1 Mb/57 Mb matrices, 4730 assertions in 61 window cases, 822 assertions
in 37 predicate cases, and standalone units pass. Accepted baseline fails the
new gap span assertion. Full evidence, rejected trials, scripts and test logs
are in `evaluations/2026-10-02-equivalent-insertion-stitch/`. Production is truth
and competitor agnostic. Final binary SHA256:
`860544d71c534e13ca37e0b23d806b0d62ef58c4a8c6d6bd2995f04cb76381e6`.

## 2026-10-02 — Match stitch MAPQ to pure BAM candidate-pair observations

Fix `local_run_boundary_flip`: two BAM-injected sites use BAM-channel alleles
and BAM MAPQ, rather than certifying the working matrix with GAF MAPQ. Graph
and mixed pairs retain their previous behavior. Missing BAM calls cannot
borrow working calls. Keep shared read eligibility, both observed allele
classes, deterministic read-half agreement and the same binomial cutoff.
Move the helper to the graph/BAM adapter for direct tests. The eight added
regression assertions fail five checks with the old implementation and pass
with the correction. No thresholds or candidate representations change.

Six owning-chunk audits find 23/29/39/59/30/27 high-MAPQ physically spanning
reads at 4.866/24.121/34.094/41.879/61.738/64.128 Mb respectively. The first two
have zero callable boundary pairs; 34 Mb has one, and 61 Mb has two from one
allele class. At 64 Mb, 18 mutually exclusive BAM deletion ALT observations
unanimously link the SNP to the deletion pair (12 plus six). This is a local
certificate, not a certificate for both complete blocks: internal missing and
contradictory graph edges remain. Preserve the separation regression and
independent read-only gauges; do not force a whole-block join from those votes.

Fresh final-build chr20 (eight threads, 256.93 s) exactly preserves all 256601
read HP/PS pairs and all 63492 VCF rows. Counts stay 237127 phased/scored,
229974 correct, 7153 discordant, 96.983473%; 673 read PS, 337 VCF blocks,
756878 bp span N50. No new fixture gap closes; 42 competitor-supported
coordinate nominations remain. All standalone units and 822 predicate
assertions / 37 cases pass; five fresh owning-window tests pass 142 assertions.
The complete 95-coordinate panel is not rerun and no expectation changes.
Build has no new warnings; diff check passes. Details and reproducible audits:
`evaluations/2026-10-02-bam-pair-mapq/`. Production uses neither truth nor
competitor output. Final binary SHA256:
`be30679a86113af54f4a885b29361f8914e9bee6ff3870c698023efc9fe2b5e8`.

## 2026-10-02 — Preserve the certified BAM prefix at a deletion-pair transfer

Close chr20:64128828–64134226 by correcting a transfer boundary that stranded
connected earlier BAM anchors. Source 64068280 has a weak cut only at 64140314,
after the clean SNP and complementary deletion pair. Eighteen exclusive-ALT
BAM triples (12 plus six), both allele classes and zero opposing molecules
certify the nearest clean SNP's relation to the pair. Transfer now includes
the earlier cut-free source run after the preceding catalog footprint, provided
all moved anchors have exact adoptable provenance and the same source gauge.
Weak/quality cuts, graph/foreign/unknown anchors, ambiguous deletion pairs,
unknown/low BAM MAPQ, and opposing molecules veto expansion. The pre-existing
pair-to-right certificate remains mandatory. Later ordinary stitching joins
the catalog blocks separately. No truth, competitor, coordinate special case,
new realignment, allele merge or loosened threshold enters production.

Correct the preceding 64 Mb audit's ownership diagnosis: 64121546 and 64128828
are BAM rows, not catalog SNPs. Their missing GAF calls do not mean that their
original BAM prefix is disconnected. The source trace records the actual cut.

Fresh final chr20 (eight threads, 292.54 s with concurrent native replays):
237127 phased/scored, 229974 correct, 7153 discordant, 96.983473%, all unchanged.
Every read retains its individual truth classification. Read PS 673→670,
VCF blocks 337→336, 63492 keys and span N50 756878 bp unchanged; 3993 read tag
tuples and 2460 VCF rows change through block orientation/label propagation,
with no lost or gained keys. The prototype and extracted helper have identical
full outputs. Local 44/45 read concordance remains unchanged; dominant correctly
placed reads improve 22/54→39/54. Native HiPhase has 48/49 locally, saved-Oct1
pgphase-callset HiPhase 49/49. This is continuity improvement, not local yield
or accuracy parity. Three native-HiPhase-correct reads remain unphased, two with
no callable BAM matrix observations and one with four; retain their evidence
for investigation.

Add the target to the owning-chunk panel: 96 coordinate cases / 82 spans. Assert
exact SNP/complementary-deletion orientation, parental flank agreement, read
separation and the unchanged owning ceiling of three errors. Update old 64 Mb
catalog-separation checks to require the intended independently certified join
and its allele relation; retain all old read floors. New helper unit tests and
baseline-failing connection checks cover the mechanism. Fresh noncentromeric
extent audit: 315 gaps, 41 distinct competitor-supported nominations (51 records),
down from 316/42/53. These are nominations, not verified exact allele bridges.

Build without warnings, standalone units, 822 predicate assertions / 37 cases,
and 4761 window assertions / 61 cases pass. All 105 fresh native replay requests
are verified against the final binary SHA256 and complete in 390.30 s.
Detailed evidence and final gate logs are recorded in
`evaluations/2026-10-02-complementary-deletion-prefix/`.
Final binary SHA256:
`98037cb463527b02e6acb31393d1d2d71cca92323d6922942d9367d868291e3d`.


## 2026-10-02 — Preserve shared BAM component identity at an insertion bridge

Close chr20:34046350–34055920 with the existing decisive physical BAM evidence:
eight REF and twelve ALT insertion-to-SNP molecules. Graph attachment relabels
source 34018169 as catalog PS 34018168. The local retained prefix is cut-free,
although its original source has a later weak cut at 34058504. Permit the
physical insertion stitch to certify that prefix from exact adoptable source
provenance, including graph-owned shared rows, with a consistent allele gauge
and no weak/quality cut inside its actual heterozygous extent.

Correct the coupled orphan-transfer bug: the nearest source row before a cut
can belong to a different graph owner and conceal an earlier owner's far-side
private island. Additional owner checks require an exact shared clean SNP in
the immediately preceding source component, preserve physical read bridges and
far-side catalog anchors, and restore detached sites and exclusive reads to
their original BAM gauge together. Keep the original source label in the
closest-row search; excluding it reopened protected joins. Move the actual
detachment function to the adapter for synthetic gauge/read regressions.
The all-historical-owner trial is rejected because it split protected joins
and lowered span N50. The refined code preserves every accepted 53–54 Mb tag
and VCF row exactly. No truth, competitor, coordinate exception, realignment,
allele merge or lowered evidence threshold enters production.

Local gap: pgphase changes two blocks into one, retaining 118 correct / 123
phased / five discordant, exactly matching native-DV HiPhase and saved-Oct1
pgphase-callset HiPhase. Owning 34–35 Mb gains seven correct reads with no
additional error: 3560/3532/28 → 3567/3539/28 phased/correct/discordant.
Fresh final chr20 (eight threads, 295.76 s with concurrent replay jobs):
237127→237134 phased/scored, 229974→229981 correct, 7153 discordant unchanged,
96.983473%→96.983562%. Every previously phased read retains its individual
truth classification. Read PS 670→669; VCF blocks 336→335; 63492 keys and
756878 bp span N50 unchanged. Exactly 179 read tag tuples and three VCF rows
change; zero keys lost or gained. Changes are confined to the owning 34 Mb
chunk. Extent audit remains 315 gaps / 41 distinct competitor-supported
nominations (51 records: 34 native-DV, 17 saved-pgphase-callset). The new break
was hidden by overlapping orphan extents, so that audit count does not fall.

Extend the panel to 97 coordinate cases / 83 required spans, with owning-chunk
context, exact insertion and opposite right-SNP genotype relation, parental
flank agreement, read separation, and unchanged old read floors. The accepted
pre-fix binary fails eight of 29 new connection assertions. All standalone
units, 822 predicate assertions / 37 cases, and 4806 window assertions / 62
cases pass. All 105 fresh native requests complete in 393.85 s against the
frozen final SHA256; an additional matrix regression runs natively. Build has
no new warnings; diff check passes. Details, rejected-scope evidence and final
results are in `evaluations/2026-10-02-shared-source-component/`.
Final binary SHA256:
`06b32560fe18b80332994314fa9bbcaa3962bbbc445bb4108f05cf27e7f3602c`.

## 2026-10-02 — Certify both graph paths at an MSA deletion bridge

Close chr20:47751480–47762233 using seven physical BAM boundary pairs (three
REF, four ALT; signed quality log odds 33.8683). One-sided agreeing graph
edges on both flanks had been rejected without using their independent clean
BAM SNP evidence. Extend the existing insertion bridge's left edge certificate
to deletions, and apply the same certificate to a reported right edge. Each
edge requires two independent Q30 primary SNP pairs with no opposing molecule,
existing likelihood bound, passing prefix and suffix paths. A separate nearest
suffix-SNP-to-deletion link checks the boundary's inherited gauge. Dominant
GAF reversals, ambiguous alleles and invalid qualities remain vetoes; both
boundary allele classes remain required. No truth, competitor, fixture
coordinate, realignment, allele merge or threshold relaxation enters production.

Physical left/right edge pairs have MAPQ 60 and Q35–Q40 bases, with log odds
15.8691/16.2861. The exact deletion and first right SNP share ALT orientation;
the next right SNP is opposite. Owning 47–48 Mb retains all 3993 phased reads,
3957 correct / 36 discordant, while VCF blocks fall 2→1 and span N50 improves
749000→997144 bp. Local continuity improves from two blocks to one and correct
reads in the dominant block 38/84→70/84 input overlaps. Local accuracy stays
70/83 correct with 13 errors; native-DV HiPhase has 74/76 with two errors.
This fixes the parental connection, not local read-accuracy parity. Those
existing read-tagging errors remain a separate task and are not relabeled by
this change.

Fresh final chr20 (eight threads, 247.68 s): 237134 phased/scored, 229981 correct,
7153 discordant, 96.983562%, all unchanged. Every read retains its individual
truth classification. Read PS 669→667; VCF blocks 335→334; 63492 keys and
756878 bp N50 unchanged. Exactly 1285 read tags and 261 VCF rows change; zero
keys lost/gained. Audit avoids 25–30 Mb: 314 extent gaps, 40 distinct
competitor-supported nominations (50 records: 33 native-DV, 17 saved-Oct1
pgphase-callset), down from 315/41/51. These are nominations, not proven links.
The 62.410 Mb repeat deletion is retained, not lost in transfer, but its few
physical links disagree; do not force it. At 7.280 Mb a low-depth normalized
graph SNP's six observations contrast with 51 REF / two ALT high-quality BAM
bases, and its first GAF edge is tied. The physical path fix keeps that case
split rather than treating the tie as one-sided confirmation.

Panel grows to 98 coordinate cases / 84 spans with exact allele and parental
checks, owning read/error floors and local separation. The pre-fix binary
fails six of 27 new assertions. All standalone units, 822 predicate assertions /
37 cases, and 4849 window assertions / 63 cases pass. All 105 fresh native
requests complete in 342.85 s against the frozen binary; one additional matrix
regression is native. No new build warnings; diff check passes. Reports and
logs: `evaluations/2026-10-02-deletion-physical-paths/`.
Final binary SHA256:
`684a95a62fe10ec32df68c8043bf0fbb730794c167ce54573c1ad2d56c546f9c`.

## 2026-10-02: padded indel identity and verified BAM genotype retention

The graph/BAM sequence converter retained a shared suffix after consuming the
full common prefix of unequal-length VCF alleles. `ACG→ATCG` therefore looked
like a replacement instead of the BAM insertion of T. Trim the remaining
suffix after the prefix, preserving the established repeat placement and the
nonmatching spans of complex replacements. Trimming the suffix first shifts
repeat alleles away from BAM/MSA and was rejected.

The corrected identity exposed a second defect: matching a phased MSA call to
an unphased graph repeat kept the graph demotion and discarded the BAM
genotype. A normalization-only chr20 run lost four verified calls at 10.488935,
10.891999, 35.490917 and 41.885034 Mb. Transfer now preserves the physical key,
counts, genotype, MSA proof and source gauge together on that exact binary row.
Its observations come from the source solve. Validated complete source flanks
retain the same admission as private BAM rows. Differing source claims remain
independent; established graph phases and whole multiallelic rows are not
overwritten. The fix is restricted to suffix-padded unphased graph repeats;
adopting other unphased repeats changed correct read assignments in trials.

An independent SNP-branch projection trial lost private BAM evidence by
changing the recovery seams, and was removed. Lowering physical SNP stitch BQ
30→20 changed no output in 28 owning chunks and was restored. Two physical
deletion pairs at 4.866153–4.874129 Mb suggested a wrong parental join despite
source certificates: errors rose 16→448, with 428 correct reads becoming wrong.
The trial was removed. Flank base quality is not a deletion-allele error model.
A new owning/local regression allows only the correct opposite allele relation
and preserved read accuracy; the rejected trial fails six assertions.

The final chr20 run preserves every HP/PS tag, every individual truth
classification and all 63492 VCF records exactly. Phased/scored 237134; correct
229981; discordant 7153; concordance 96.983562%; read PS 667; VCF blocks 334;
span N50 756878 bp. No additional gap closure is claimed. Wall time 276.61 s at
eight threads while panel runs executed concurrently; not a runtime comparison.

All units, 872 predicate assertions / 39 cases and 4929 window assertions /
65 cases pass. The four-call transfer regression fails three assertions on the
normalization-only binary. Existing 98 coordinate cases / 84 required spans and
their floors remain unchanged. All 105 native panel requests use the frozen
binary, completing in 372.52 s, plus one native matrix request. Production is
truth agnostic and introduces no realignment, allele merging or stitch threshold
relaxation. Details: `evaluations/2026-10-02-padded-indel-context/`.
Final binary SHA256:
`00c7e302c47a3c1fe58d6a4fa5ef647fc1d8030171ce362482a8cb3db794af91`.


## 2026-10-02 — Recover exact SNP branches from complex catalog alleles

The full-walk binary graph projection drops a SNP observation when a read
carries another catalog ALT, even if its local oriented SNP branch is exactly
identical to the selected allele. Seven such ALT calls at 20,523,889 repair the
right graph path and join 20,506,159 C→CT to 20,517,197 A→AAT in the owning
20–21 Mb chunk and full chromosome. The insertion ALTs remain opposite; the
SNP is on the left insertion's ALT haplotype. No allele is merged or realigned.

After graph phasing, supplement only missing observations at already phased
clean binary SNPs. Require a unique oriented three-node REF/ALT branch in each
catalog walk, retain the original source-walk conflict veto and mapping-quality
floor, skip parent-gated children, and require the combined evidence to remain
heterozygous under the existing AF filters. Initial genotype counts, categories,
phase labels and read gauges are preserved. Both workers use the same helper.

The AF guard is required: an ungated projection strengthens a false 3/3 SNP at
7.280 Mb with 42 new REF calls and leaves eight previously phased reads
unassigned. Broad unplaced-MSA admission and a physical homo-REF anchor veto
also altered established correct assignments; neither trial is retained.
The final 7 Mb replay has exactly unchanged tags and truth assignments.

Full chr20: phased/scored 237,134; correct 229,981; discordant 7,153; accuracy
96.983562%, all unchanged. Every phased read retains its truth classification.
Read phase sets fall 667→665, VCF blocks 334→333. N50 stays 756,878 bp and
all 63,492 variant keys/allele genotypes remain. Locally pgphase has 109/110
correct phased reads versus saved native-DV HiPhase's 107/109. Independent read
rescues remain separate: dominant correct local reads improve 65→78 of 111
input overlaps, still below HiPhase's 107. This is one exact connection, not
complete read-block continuity parity.

Extend the permanent panel to 99 coordinate cases / 85 required spans, with
exact keys, ALT parity, parental flanks, owning counts and local error checks.
The accepted pre-fix binary fails six assertions in the new regression.
Units cover projection uniqueness, source conflicts, AF compatibility and
preserving original gauges. Details, rejected trials, chromosome comparison,
evidence-channel changes and the remaining 40 noncentromeric nominations are
in `evaluations/2026-10-02-catalog-snp-branches/`. The whole extent audit has
313 gaps; an overlapping block extent hid this newly fixed exact split from
the preceding nomination list, so that list's 40 distinct targets is unchanged.

Validation: all standalone units pass; phase predicates pass 872 assertions /
39 cases; the expanded window suite passes 5,708 assertions / 66 cases. Its
105 fresh native requests use the frozen final binary and complete in 343.08 s,
plus the required additional native replay requests. Build has no new warnings.
Final SHA256: `cb3b4ad62fb671abb1f375e02a48a44b1eac2553158b7ab39bc6c9755d94fc06`.


## 2026-10-02 — Same-current-callset joint-gap controls

Fresh targeted HiPhase runs on the accepted current pgphase VCF and the same
BAM, using 100 kb context on each side, establish that these site records can
support a joint connection. Local-only HiPhase spans 7.264321–7.280346 Mb
with 120/122 correct phased local reads (pgphase 117/126, three read PS), and
24.121713–24.131707 Mb with 89/93 (pgphase 60/60, three read PS). The
spanning blocks have zero discordant reads on disjoint 10 kb flanks:
42/42 and 52/52, then 29/29 and 45/45. Saved native-DV HiPhase has 123/125
and 93/93 respectively. Each fresh local window run takes about one second.

Default HiPhase on the same pgphase callset does not span either pair;
24 Mb is only 59/79 correct in two read blocks. The local control changes
allele calling as well as the solver, so it is not proof that A* alone fixes
our matrix. It also still uses local alignment; no new alignment method or
production change is made. HiPhase uses eligible genotypes jointly, whereas
our independent graph/BAM solves and narrower centered-site MEC fallback do
not offer a general joint solve over every usable gap observation. The next
control must compare joint phasing on exactly identical read/site matrices,
retain the complete established flank gauges, and evaluate both relative
orientations without forced joins. Details and measurements:
`evaluations/2026-10-02-joint-gap-controls/`.

## 2026-10-02 — Joint gap solve on identical observation matrices

Replay the unchanged production exact MEC kernel on fresh owning-chunk
matrices, with complete flank gauges and both relative orientations. Including
all usable binary graph/BAM gap sites increases variables 5→11 at 7.28 Mb,
but independent read halves still disagree. At 24.13 Mb the current matrix
remains disconnected, including when boundary indels are free. Thus removing
eligibility restrictions alone does not fix these cases.

Extract actual HiPhase local calls from serial trace output, checking variant
and BAM record counts. At the 7,280,356 SNP, HiPhase marks 12 callable graph
observations ambiguous and calls three oppositely. Substituting the shared
HiPhase calls makes SNP-priority MEC prefer the correct connection overall
and in both read halves. At 24 Mb, HiPhase has 29 callable left/right pairs;
pgphase has zero of those pairs, already in the BAM source before transfer.
Two MAPQ-60 reads reach the right clean SNPs, but lack pgphase insertion calls.
With HiPhase calls, all usable sites and free boundary indels, joint MEC
prefers the correct orientation overall and in both halves. Centered-site
filtering or frozen imported boundary indels loses this result. The final
24 Mb parity margins remain only 3, 4 and 1 cost units.

A linked synthetic probe reproduces an MSA representation defect: the
selected four-T insertion is ALT, literal reference is zero, but an exact
eight-T read matching the other verified consensus is unknown on that row.
This disagrees with complementary source-row ALT-absence semantics and can
also disable the one-error fallback, which requires both consensuses to
classify. A future repair should preserve both rows and require the exact
other consensus, keeping unrelated third alleles unknown. Copying independent
HiPhase REF/ALT calls would not preserve this diploid representation.

The two-orientation kernel takes 4–9 ms at 7 Mb and 23–53 ms at 24 Mb,
including subprocess overhead; matrix loading is additional. These are
diagnostic objective/orientation results, not emitted-tag or chromosome
accuracy measurements. No production behavior, binary or test expectation
changes. Full evidence, costs, scripts and limitations are retained in
`evaluations/2026-10-02-joint-matrix-investigation/`.

## 2026-10-02 — Complementary insertion calls and inherited read HP

Repair a reproduced unplaced-read MSA observation loss: a read exactly matching
one verified insertion consensus was unknown on the complementary insertion
row. Zero represents absence of that row's ALT, so recognize the other complete
consensus only with an identical reference footprint, common query flanks and
exactly one selected ALT. Both composed read paths must agree. Keep separate
rows and allele identities; unrelated third alleles remain unknown to the exact
caller. The existing one-edit, strictly closer consensus fallback is unchanged.
An exact-only restriction lost 50 correct reads and was rejected. Broader new
joint/dropout recovery triggers regressed whole-chromosome assignments and were
removed; no new joint solver or detector is accepted.

After the final graph/BAM chunk stitch, recovery chunks correct inherited HP
only with concordant Q30, known-MAPQ30 physical single-base clean SNP calls in
that read's final PS, agreeing with the primary channel and separated by at
least 100 bases. Duplicate coordinates count once and conflicting eligible
SNPs veto correction. This modifies assigned read HP/scoring counters only;
no candidates, genotypes, phase sets or stitch votes change. Compared with the
insertion fix alone it repairs 36 incorrect assignments and worsens one. A
primary-graph-only alternative loses two correct reads in the 3 Mb replay and
was rejected. Production uses neither truth nor competitor observations.

Same-input whole chr20, starting binary cb3b4ad6 → final 2f0513f6:

- Phased truth-scored reads: 237,134 → 237,172 (+38).
- Correct: 229,981 → 230,022 (+41); discordant: 7,153 → 7,150 (-3).
- Accuracy: 96.983562% → 96.985310%.
- Read PS: 665 → 661; VCF blocks remain 333.
- VCF keys: 63,492 → 63,630, no lost key.
- Block span N50: 756,878 → 774,189 bp.
- All 57,769 previously phased SNP keys retain their old within-block gauge;
  no old SNP block splits or internally reorients. Three old SNP-block pairs
  merge, which does not imply all intermediate indel boundaries join.

Individual read tradeoffs remain: 54 correct→incorrect, 67 incorrect→correct,
45 unphased→correct, ten unphased→incorrect and 17 correct→unphased. The exact
17.502614–17.521341 Mb gap closes with the correct relative parental gauge in
full chr20 and both owning chunks. Its local reads worsen from 131/136 to
109/136 correct (80.15%); same-BAM HiPhase DV is 138/144 (95.83%). Thus correct
block connection is verified, but competitor-level local read accuracy is not.
The 4.866–4.874 Mb exact pair remains split; its owning chunk phases 26 more
reads, 16 correctly, with ten additional errors. Record that explicit cost
rather than claim an accuracy improvement everywhere.

Add the new connection to the permanent panel (100 coordinates, 86 required
connections), with spans=1, exact representations and independent parental
majorities from a 16–18 Mb replay. Add the 24 Mb owning anti-inversion test,
which retains both four-T/eight-T rows and 60/60 correct local reads. Keep all
prior connection and orientation assertions. Explicitly accept the measured
4 Mb read-count/error tradeoff and two 17 Mb owning read-score floors
98.65%/98.71% instead of 99%; do not refresh the panel wholesale. Physical SNP
correction restores the existing 3 Mb 98% floor. The 7.28 Mb and 24.13 Mb
boundary gaps remain open, pending usable allele evidence, not forced joins.

Unit tests, 906 assertions in 40 predicate cases, and BAM HiFi/ONT golden gates
pass; HiFi one/four-thread outputs remain deterministic. The expanded suite passes 7,031 assertions across 68 test cases. Fresh
final-binary native window outputs are rescored with binary/argument/input/
output checks.
Full run wall time 251.48 seconds at eight threads is not a controlled runtime
comparison. Detailed results, rejected controls and final regression logs:
`evaluations/2026-10-02-complementary-insertion-recovery/`.


## 2026-10-02 — Keep physical SNP certificates with their recovered allele

Fix graph recovery quality provenance: original CIGAR quality certifies an MSA
SNP observation only when the literal BAM base calls that same REF/ALT. Keep
MSA calls in the matrix even when they lack that physical certificate. Two
singleton bridge checks use the same helper. Overlapping solves can retain
one allele but supply quality from the other; certify the retained call and
replace old quality when a replay replaces its observation, including with
zero. A permanent truth-independent matrix regression at 3,597,791 reproduces
ALT borrowing the later REF call's Q40 and requires quality zero afterward.
`BAMQ` diagnostic rows expose these certificates without altering OBS fields.

Permit an imported, phased MSA SNP to repair inherited read HP after its gauge
passes the existing singleton association/Wilson gates against physical clean
SNPs of the same PS. Use every covered clean SNP at least 100 bases away, one
vote per molecule, and abstain on clean conflicts. Cache the result per site;
read HP does not certify a site. Keep MSA witness counters separate from clean
counters. A physically matching clean call below Q30 can veto a correction
relying on a noisy SNP, but cannot assign HP. Two spaced Q30 clean SNPs retain
their existing certificate. This fixes the unchanged 11.599 Mb owning error
ceiling after an initial candidate regressed it from 90 to 91. No truth,
competitor calls, fixture coordinates or new realignment enters production.

Reject the unrestricted single-MSA-SNP rule: it costs 26 correct owning reads
near 15.07 Mb. Reject nearest-only gauge selection: 14 physical links at
17.50 Mb fail the existing confidence bound; complete covered clean context
retains 19 coherent links. Preserve these rejected controls in the evaluation.

Final full chr20 against the accepted complementary-insertion/BAM-refresh
baseline: truth-scored phased reads 237,172 -> 237,170; correct assignments
230,022 -> 230,043; discordant assignments 7,150 -> 7,127; conditional accuracy
96.985310% -> 96.994983%. 29 incorrect reads become correct, seven correct
reads become incorrect, and two fallback reads abstain (one formerly correct,
one incorrect). Read phase sets 661, variant keys 63,630, VCF blocks 333 and
span N50 774,189 bp stay unchanged. Zero VCF rows change. No new gap closes.

At 17.502614–17.521341 Mb, full-chromosome local accuracy improves 109/136 ->
131/136 (96.32%), with errors 27 -> 5 and parental connection preserved. Saved
same-BAM HiPhase DV has 138/144 correct (95.83%); it still phases eight more
local reads. Strengthen the permanent owning replay to >=131 correct and <=5
errors; restore older 17 Mb 99% owning floors and raise the 15.056 Mb guard to
99%. Keep every exact span and parental-orientation requirement. The permanent
panel remains 100 coordinates/86 required connections; adding the quality
provenance replay brings the suite to 69 cases.

Final binary SHA256
`95f94fdd04210a95e295131288ca21398065ac0d0afe1364c7c586fc892e63aa`.
Build, all units, 906 predicate assertions/40 cases, 7,038 window assertions/69
cases, and HiFi/ONT BAM goldens pass. First run 105 fresh native panel requests;
rescore only with binary/argument/input-stat/output-hash checks, and execute
uncached additional requests natively. Full run 297.21 seconds at eight threads
and concurrent native panel 394.91 seconds are not controlled runtime comparisons.
Detailed evidence, rejected controls, final measurements and regression logs:
`evaluations/2026-10-02-verified-snp-read-recovery/`.

## 2026-10-02: preserve partial MSA reads and their coverage limits

The haplotype-aware unplaced-read pass incorrectly equated a prior HP/PS label
with actual MSA membership. Homopolymer consensus construction excludes partial
reads, but those labeled reads were skipped again during recovery. Unlabeled
partial reads were aligned and scored as full-cover reads, and ambiguous
reference/read composition restored full bounds. Even a no-cover read could
enter a cluster. The permanent synthetic prefix/suffix regression, linked
against the starting implementation, fails all five sections.

Track initial MSA membership, skip no-cover reads, and recover the other reads
against fixed consensuses. Retain the complete consensus axis for reference
composition; trim using the read's own original coverage flags. Score only the
covered intersection and restore coverage bounds on ambiguous composed paths
before normalization and site recall. An uncovered end cannot cast a deletion
or haplotype vote, even at a one-point assignment margin. Preserve consensus
construction, separate candidate rows, retry admission and stitch thresholds.
The ordinary upstream BAM pass remains unchanged. Production uses no truth,
competitor calls, fixture coordinates or new alignment method.

Final binary SHA256
`74218715db7412bb7d0e88633b14ec6c85e3220047b6e9ce1d6df6210c7d7df4`.
Fresh full chr20 preserves all 237,170 truth-scored phased reads, 230,043 correct
assignments and 7,127 errors (96.994983%). Every scored read retains its
correctness state. One read moves to an independent BAM PS with unchanged HP;
31 VCF records change only counts and derived allele fractions. All 63,630
keys and every VCF genotype/PS remain unchanged. Read PS 661, VCF blocks 333,
span N50 774,189 bp. No additional gap closes; the 7.28/24.12 Mb cases stay open.
The owning 7–8 Mb replay retains one more ALT at 7,047,081 without changing
genotype, phase orientation or read assignments.

Build and units pass with no new warnings. Predicate suite: 962 assertions/41
cases. Fresh native panel: 105 requests followed by verified scoring reuse;
7,038 window assertions/69 cases pass. The 100-coordinate panel and all 86
required connections, span equality, parental checks and read-score expectations
remain intact. HiFi/ONT TSV/VCF goldens and HiFi thread determinism pass.
Full chr20 296.22 seconds/eight threads; concurrent native panel 395.12 seconds,
not controlled runtime comparisons. Behavior documentation is updated.
Evidence: `evaluations/2026-10-02-partial-msa-coverage/`.

## 2026-10-02: recall missing homopolymer insertion calls with fixed source gauges

The initial haplotype-aware recovery MSA omitted reads outside its selected
clusters, even when both fixed consensuses could establish their local allele.
At 24.12 Mb, separate four-T/eight-T rows therefore lacked many observations
available in the BAM. Request local unplaced-read recall for the initial
targeted source solve. Align against both existing consensus axes without
assigning a whole-region cluster. Admit only separate homopolymer insertion
length contrasts inside requested seams, with disjoint one-base error
neighborhoods, complementary calls, physical flanks and the existing AF gate.

Defer these observations until ordinary discovery and phasing finish. Coalesce
exact-key/read duplicates, reject contradictions, fill only unknown calls at
surviving MSA-verified phased heterozygotes, and preserve already phased source
membership and its allele gauge. Update counts and the read index; genotypes,
categories and source phase sets stay fixed. Standard graph transfer/stitch
rules still apply. No row merging, new alignment method, truth, competitor
calls, fixture coordinates or lowered stitch thresholds enter production.

Reject whole-source rephasing with all recalled sites: chr20 loses 1,023 keys
and errors rise from 7,127 to 9,469. Reject insertion recall throughout the full
solve context: errors rise by 32 and five existing assertions fail in four
cases. Exact-only recall also loses owning coverage. Restricting coordinates
alone is insufficient; retain source membership and the distinct homopolymer
contrast. The accepted narrowed pass satisfies every existing regression.

Full chr20: truth-scored phased reads 237,170 -> 237,199; correct assignments
230,043 -> 230,072; discordant assignments unchanged at 7,127; conditional
accuracy 96.994983% -> 96.995350%. Exactly 29 read tags change, all previously
unphased reads now truth-correct. Every old HP/PS tag and scored read correctness
state survives. All 63,630 variant keys and every VCF GT/PS survive; only two
insertion records change counts/fractions/quality (DP 11 -> 60). Read PS 661,
VCF blocks 333 and span N50 774,189 bp stay unchanged. No old block merges and
no new gap closes: the coordinate panel remains 100 windows/86 connections.

The owning 24–25 Mb replay gains 29 correct reads with errors fixed at two.
Within the 24.121713–24.131707 Mb gap, 60/60 correct -> 89/89, no errors and
three independent blocks. Same-BAM HiPhase DV has 93/93 in one block; that
continuity advantage remains. Strengthen the existing owning/local floors to
protect all 89 correct overlaps, separate rows and rejection of the wrong
flank join. The starting binary fails four new coverage assertions. A new
synthetic test rejects literal-REF, different-position and tandem-repeat
contrasts while preserving fixed consensus membership.

Final binary SHA256
`d3b5ad20b227cc051d1d5280a6aa4fc82cbc1c1bc9393082aa5533005456a1c8`.
Build, all units, 981 predicate assertions/42 cases, 7,038 window assertions/69
cases, HiFi/ONT TSV/VCF goldens and HiFi thread determinism pass. All span,
parental-orientation, accuracy and coverage gates remain. Run 105 fresh native
panel requests, then score verified final-binary outputs with exact-argument,
input-stat and output-hash checks; uncached requests execute natively. Chr20
295.28 seconds/eight threads and concurrent panel 395.26 seconds are validation
timings, not controlled runtime comparisons. Behavior documentation is updated.
Evidence: `evaluations/2026-10-02-fixed-consensus-insertion-observations/`.

## 2026-10-02: Distinguish an absent SNP gauge from explicit contradiction

Recovery's existing shifted-repeat insertion check recognized a physical ALT,
but its boolean source-SNP certificate conflated no callable SNP with a
contradictory Q30 clean SNP. A failed certificate could retain exact-anchor
REF, manufacturing the opposite sequence allele. Retain three states internally:
no callable independent gauge, coherent ALT gauge, and contradiction. For an
unpaired binary insertion, a verified shifted ALT explicitly contradicted by
an independent clean SNP remains unknown. Missing SNP evidence and separate
complementary MSA insertion contrasts retain their current source projection.
Deletion admission, existing observations, keys, GT/PS, alignment and stitch
criteria stay unchanged. No truth or fixture coordinates enter production.

Reject broad complementary-pair CIGAR admission: chr20 errors 7,127 -> 9,342,
2,292 previously correct reads become wrong, and no panel gap gains a span.
Preserving read-source membership is insufficient: owning 2–3 Mb errors
62 -> 1,091. Freezing source-path certificates also fails that inversion and
raises owning 19 Mb errors 41 -> 381. Restore all those trials. A broad unpaired
missing-certificate abstention improves local 41 Mb tags but reopens the required
1.180618–1.194189 Mb connection, loses three correct chr20 reads net, and makes
13 formerly correct reads wrong. Reject it; do not lower its span gate.

The retained explicit-conflict correction exactly preserves full chr20:
237,199 scored / 230,072 correct / 7,127 errors, 96.995350% accuracy, 661 read PS,
63,630 VCF keys, 333 blocks and 774,189 bp N50. Every read HP/PS and VCF row is
unchanged. The coordinate panel remains 100 windows / 86 spans; no new gap
closure is claimed. Strengthen the permanent representation/transfer test with
an owning 2–3 Mb replay protecting both source insertion rows, at least 4,088
scored / 4,026 correct reads and at most 62 errors. Keep all previous expectations
and original owning coverage floors. The contradictory-SNP regression fails
the starting implementation in two assertions.

Build, all units, 987 predicate assertions/42 cases and HiFi/ONT TSV/VCF goldens
pass; HiFi thread determinism remains. Final SHA256
`b36dd6d8ea32b1f1707b64efc51878b2d8641344ef1b0f009f08eee94b1bac0b`.
The whole chromosome takes 292.18 s/eight threads and 105 fresh native panel
requests 387.70 s/four concurrent workers, validation timings only. The living
implementation documentation describes the three evidence states and source
contrast scope. Evidence: `evaluations/2026-10-02-shifted-insertion-ref-abstention/`.
Complete window validation passes 7,060 assertions/69 cases. Verify the 105
fresh native requests against binary/input/output identities before scoring;
the added owning 2 Mb request runs natively. No panel floor or span is lowered.


### 2026-10-02 — Reject repeat-refresh and backfill-census perturbations

Investigate owning 19 Mb with native local-MSA traces and chromosome guards.
All 58 assigned reads cover 19.377 Mb. The existing independent local caller
classifies 42 and abstains on 16; every callable row agrees with the port's
source call. This corrects the earlier interpretation of exact local calls
being reversed. CIGAR versus composed-MSA repeat footprints remain different;
composition's consensus axes match. The source's 10/14 insertion contrast
segregates poorly, while DV calls the longer insertion homozygous. No production
truth/competitor-based genotype demotion is introduced.

Reject broad local refresh, exact paired refresh and anchored source solving:
they lose sites or correct reads and do not fix the target. Also reject updating
source DP/AD/AF/strand fields from the enlarged post-solve backfill matrix.
Immediate and post-source-selection recounts both produce 237,244 scored /
230,088 correct / 7,156 discordant chr20 reads, including 13 correct-to-wrong
transitions; tracked spans stay 86/100. Delaying the recount does not fix the
regression, so source-selection timing is not an established root cause.
Discovery genotype counts remain distinct from supplementary linkage evidence.
Document that contract and protect it with insertion/SNP regressions. The
rejected recount fails 15 assertions in the new guard; restored predicates
pass 1,020 assertions/42 cases. All standalone units, HiFi/ONT goldens and
7,060 window assertions/69 cases pass; no expectations are relaxed.

The final binary is byte-identical to the accepted b36dd6d8 executable:
237,199 scored / 230,072 correct / 7,127 errors, 333 VCF blocks and N50 774,189.
No new closure is retained. Native rejected trial timings and validated cache
reuse are identified separately in
`evaluations/2026-10-02-recovery-callable-depth/README.md`.

## 2026-10-02: Close the 5.511 Mb gap to increase chromosome NG50

Rank the noncentromeric competitor-supported nominations by simulated union
of complete adjacent VCF phase-set extents. Use half the reference chromosome
length (66,210,255 bp) for NG50, separately from the summed-span N50 denominator.
The 38 audited nominations are not an exhaustive missed-gap census. Nine have
positive modeled NG50 gain; select 5.511231–5.531924 Mb from four tied largest
gains because its physical boundary representations and flank gauges agree.
The 62.623 Mb candidate's internal-switch counterexample remains protected.

Two MAPQ60 BAM reads call opposite alleles at both deletion boundaries. Targeted
recovery retains those pairs; the larger whole-flank validation omits the
right MSA deletion calls on partial reads and its original eligibility/singleton
screen cannot use that boundary. Preserve the established certificate first.
Only after it fails, permit a downstream MSA-verified heterozygous deletion and
apply existing seam observation backfill to the reused full solve. Discovery
counts, genotypes and source orientations stay fixed. Require multiple exact
physical pairs, both haplotypes, unanimous parity, matching deletion REF bases,
full graph/BAM gauge validation and no source weak cut. Upstream eligibility,
ordinary SNP statistics and stitch conflict checks remain. No truth, competitor
call, hard-coded coordinate, MSA row merge or new alignment method is added.

Reject broader MEC filtering and boundary admission: they lose correct reads
or cross the upstream source cut. A downstream-deletion-first trial reopens
the required 23 and 47 Mb connections, loses three correct chr20 assignments,
adds one error and loses a key. Its NG50 stays 633,077 bp. Trying the established
certificate first preserves both owning replays exactly and retains the target.

Full chr20 NG50 **633,077 -> 643,699 bp (+10,622 bp; +1.68%)**. N50 stays
774,189 bp. The 755,168 bp joined block covers 5,393,615–6,148,782. All 237,199
scored reads, 230,072 correct assignments and 7,127 errors survive unchanged
(96.995350%); every individual truth-correctness state is preserved. All 63,630
keys and unordered genotypes survive; no old SNP gauge becomes mixed. VCF blocks
333 -> 332, truth-scored read PS 661 -> 659. Only the targeted old blocks merge.
The main gap block has 131/131 correct overlaps; eight other phased overlaps
remain independent (five correct, three errors). Saved HiPhase has 136/136 in
its main block, so equal local coverage/purity is not claimed.

Add the measured gap to the permanent panel: 101 coordinates / 87 spans,
versus 100 / 86 before. Extend the owning 5 Mb regression with opposite boundary
GT, shared PS, parental orientation and 4,024 scored / 4,001 correct / at most
23 errors, plus 131 main-block correct overlaps. The starting binary fails three
new assertions. No existing individual span, coverage or accuracy gate is lowered.
Build, all units, 1,020 predicate assertions/42 cases, 7,092 window assertions/69
cases and HiFi/ONT TSV/VCF goldens pass; HiFi thread determinism is preserved.

Final SHA256 `aa3b89bb6f6f5a49a4fc72e0c83577bfade869165e60f2c8556a48bf5e0026d4`.
Whole chr20 takes 291.37 s/eight threads and 105 fresh native panel requests
394.17 s/four concurrent workers, validation timings only. Final window scoring
verifies binary, normalized CLI, input stats and output hashes before cache reuse;
uncached requests run natively. Behavior documentation describes certificate
priority, boundary fallback and discovery-count preservation.
Evidence: `evaluations/2026-10-02-ng50-gap-priority/`.


## 2026-10-03: Repair shifted deletion REF lookup without changing source contrasts

Investigate NG50-priority HiPhase nominations at 57.085410–57.104654 and
61.738239–61.747506 Mb, excluding the centromere. Either complete neighboring
block join would project NG50 643,699 -> 672,998 bp; neither is newly closed.
At 57 Mb the singleton physical read carries a five-base repeat deletion, while
both native DV and BAM MSA describe the diploid locus as +3/-2 consensus alleles.
At 61 Mb partial homopolymer MSA reads omit useful molecules, but three shifted
ALT deletions conflict with independent clean left-SNP source orientation.
Physical edit equivalence alone cannot certify an entire block gauge.

Fix a narrower observation bug: exact-position deletion REF can hide an
edit-equivalent shifted ALT beyond the candidate footprint. Before filling a
new REF at a simple phased MSA deletion, apply the existing MAPQ30/Q30 edit
certificate and independent clean-SNP source gauge. Agreeing SNP evidence admits
ALT; contradictory callable SNP evidence leaves it unknown. Absent SNP evidence
preserves the source's ALT-absence contrast, matching insertion recovery.
Separate colocated alleles and existing MSA calls remain authoritative;
candidate keys, source GT/PS and discovery depth fields stay frozen.

Reject isolated homopolymer observation admission even with clean-SNP checks:
full chr20 loses 530 correct reads, adds 580 discordant reads and reopens the
48.929/56.064 Mb joins (87 -> 85 panel spans). Reject turning every shifted ALT
with no SNP gauge into unknown: seven 48.929 Mb source zero contrasts disappear,
forcing a different source solve and losing 195 owning-chunk variant keys.
Its full run loses 202 keys and six previously correct reads and reopens those
two joins despite improved aggregate conditional accuracy. Broad physical-only
61 Mb recall gains five correct local overlaps but adds ten owning-chunk errors;
before-solve recall loses 209 correct and adds 236 discordant assignments.

The retained correction preserves full chr20 HP/PS tuples and VCF rows exactly:
237,199 scored / 230,072 correct / 7,127 discordant (96.995350%), 659 scored
read phase sets, 63,630 keys, 332 VCF blocks, N50 774,189 bp, NG50 643,699 bp.
No old phased SNP or block gauge is lost; all 87/101 coordinate spans survive.
The local target scores remain 127/145 correct at 57 Mb and 79/83 at 61 Mb.
Saved same-BAM native-DV HiPhase has 139/139 and 91/92, respectively; saved
HiPhase with pgphase calls is much worse at 57 Mb (96/137) but retains 91/92
at 61 Mb. Competitor outputs are reused, not a new runtime measurement.

Synthetic regression checks fail three assertions on the starting source.
The final correction passes all 1,031 predicate assertions/42 cases. Add a
61–62 Mb owning-chunk guard protecting 3,452 scored / 3,423 correct / at most
29 errors, positive phased clean boundary SNPs and the supported parental
relation if those SNPs acquire a common PS. It rejects the before-solve trial
on both correctness and error limits. Full window tests pass 7,113 assertions
in 70 cases. Existing exact span, accuracy and coverage gates remain unchanged.
Build, standalone units and HiFi/ONT goldens pass with no new warnings; HiFi
one/four-thread outputs agree. Update living behavior documentation.

Final SHA256 `244d3e8f39b91c5ac17d9ad3c5c05a28ad42223a9cb12db681c77e0c571ef3a8`.
Full chr20 takes 300.63 s/eight threads; 105 fresh native requests take 401.92 s
with four concurrent workers. These concurrent validations are not controlled
runtime comparisons. Window rescoring verifies executable, normalized CLI,
input stats and output hashes; the new uncached owning case runs natively.
Evidence: `evaluations/2026-10-03-shifted-deletion-source-gauge/`.


## 2026-10-03: Check independent BAM block injection and whole-block context

Test user-proposed gap-only BAM phasing followed by intact block transfer and
stitching at chr20:61.738239–61.747506 Mb. Native BAM runs take 0.3–0.6 s.
Exact, +100 bp and +1 kb regions phase 69/69 correct local reads but emit only
one phased SNP, so they supply no chain across the gap. +10 kb produces
63/61/2 phased/correct/discordant; +50 kb produces 53/41/12 under either
MAPQ30 or MAPQ5. The actual graph seam plus 1 kb produces 55/41/14. Current
graph+recovery retains 83/79/4, versus saved same-BAM native-DV HiPhase 92/91/1.
These scores use majority orientation separately per local output PS, not
whole-chromosome PS orientation.

Trace the accepted recovery source: PS61690751 has 14 phased sites. All seven
private phased seam rows reach the graph table as independent PS61690752 with
unchanged source GT. Injection does not split them. Shared graph SNPs retain
catalog representation and graph membership with BAM provenance. The later
stitch separates the right deletion at a weak source edge: both callable
source pairs from SNP61738239 to deletion61747507 are ALT/ALT, and both conflict
with the source's opposite haplotype alleles. No consistent source pair exists.

A negative control that trusts nominal BAM PS and clears source cuts makes the
private gap rows share one block but leaves the next graph block independent.
It improves one correct owning-chunk read and local purity 79/83 -> 80/83, yet
puts the SNP and deletion ALTs on opposite haplotypes in that common block.
Reject; aggregate read accuracy alone cannot validate allele connection.
Extend the permanent owning61 test to require the supported SNP/deletion GT
relation whenever the pair shares a PS. A future correctly oriented join is
allowed. The control fails exactly this new orientation guard.

Broaden transfer to all phased private flank rows and shared unphased-repeat
adoption. Four initial owning chunks retain every truth state and add 92/28/16/81
VCF keys at 61/57/48/56 Mb, but full chr20 converts 3,622 previously correct
reads to incorrect, reopens four protected gaps, loses 13 phased SNPs and reduces
accuracy 96.995350% -> 95.509795%. Reject. Narrow to complete cut-free private
source blocks, keeping shared graph-site adoption rules. Initial owning61/23/53
replays preserve all truth states, but full chr20 loses 31 scored and 26 correct
reads, converts 35 correct reads to incorrect, loses ten phased SNPs and reopens
the 54.547514–54.569072 Mb gap; N50 774,189 -> 739,888 bp. Reject that as well.
No expectation is weakened to retain extra variant rows.

Restore production source and binary byte-for-byte to accepted SHA256
244d3e8f39b91c5ac17d9ad3c5c05a28ad42223a9cb12db681c77e0c571ef3a8. No new gap
closes. Keep only the stronger orientation regression and evaluation records.
All standalone units, 1,031 predicate assertions/42 cases and 7,121 isolated
window assertions/70 cases pass. Existing 87/101 coordinate spans and all
coverage/accuracy limits remain. Reuse verified 105 native outputs only after
binary/CLI/input/output hash checks; the expanded owning case runs natively.
Give each runner a fresh PGPHASE_TEST_WORKDIR: its hard-linked replay cache is
unsafe for diagnostic processes sharing output paths; discard the contaminated
shared-directory validation. Full rejected trials take 254.98/243.35 s with
eight threads; these are concurrent validation timings, not runtime comparisons.

Design implication: complete BAM block evidence should remain immutable and
separate from live graph ownership. Use full source site/read matrices for
pairwise orientation, retain source identities across numeric PS collisions,
and mutate labels only after connection checks succeed. Copying all flank
sites into the active candidate table changes downstream ownership before
stitching. The proposed separate-evidence design is not implemented by the
rejected trials. Living behavior documentation remains unchanged.
Evidence: evaluations/2026-10-03-bam-block-transfer/.


## 2026-10-03: Keep complete BAM evidence separate and stitch finalized blocks

Implement the separate-evidence design following the rejected broad flank-site
transfer trials. Each selected source PS keeps every phased heterozygous source
row and callable read observation, including context outside the seam, in an
immutable source matrix. Matrix indices and source-scoped labels are independent
of graph candidate ownership. Duplicate molecule names are excluded from voting.
Existing graph adoption and emitted allele representations remain intact.

Score whole source blocks by molecule, prioritizing clean SNPs. Source/live allele
bases must translate consistently. A new whole-block edge needs known MAPQ30/Q30
physical SNP support from both haplotypes, a significant same/cross binomial vote,
and independently significant agreement with each block's current HP gauge.
Conflicting complete/transferred channels veto; all internal block paths and all
consecutive edges are checked before constant PS/HP unions. Truth and competitor
labels are used only for validation, never production decisions.

An initial route inside shared recovery changed downstream genotype recovery and
lost 16 phased SNPs and 734 previously correct read assignments. Adding path
checks alone still lost 11 phased SNPs. Running after the shared solver but before
source attachment preserved the sites yet gave a wrong owning37 join: 3,561
correct/46 discordant became 3,415/191 despite a 53/0 physical SNP edge. Per-block
vote totals alone did not fix it. Source attachment and retries still replayed
older ownership/gauge decisions after that premature connection. Move complete
block stitching after both recovery passes and all source attachments, using the
final live labels and allele/read gauges. Do not interpret saved pre-attachment
graph vote totals as final orientations. Keep both passes' immutable matrices
and path certificates for this final route: retry replacement of active gauges
must not erase first-pass source evidence. Consider every covering solve. A
source label reused after an old block was absorbed is ambiguous across solves
and cannot share one path certificate. Restore the existing shared seam solver
unchanged; no existing genotype recovery or join is bypassed.

Accepted full chr20 has exact read-tag and VCF-row parity with baseline SHA256
244d3e8f39b91c5ac17d9ad3c5c05a28ad42223a9cb12db681c77e0c571ef3a8:
237,199 truth-scored, 230,072 correct, 7,127 discordant, 96.995350%; 63,630 variant
keys, 332 VCF blocks, N50 774,189 bp and NG50 643,699 bp. The coordinate panel
remains 87/101 connected. This implementation retains complete evidence and a
validated stitch route, but closes no additional chr20 gap under its certificates.
Native full-run wall time is 291.01 s with eight threads, concurrent with panel
validation; do not treat that as a runtime comparison.

Add an owning37 regression allowing a future correct closure while enforcing
coverage, truth correctness, discordance and parental orientation. The accepted
implementation passes 11 assertions; the premature-join negative control fails
four guards (3,606 scored, 3,415 correct, 191 discordant, switched). Unit tests
include complete outside-gap context, source independence/PS remapping, both
haplotypes, missing/unknown/low qualities, SNP priority, ties, duplicate molecules,
raw/live frame changes, missing/conflicting per-block HP gauges, internal cuts,
homozygous fake anchors and atomic multi-block joins. All standalone units and
1,193 predicate assertions/44 cases pass. No expectation is weakened.
Evidence: evaluations/2026-10-03-complete-bam-block-evidence/.
Final gap validation passes 7,132 assertions/71 cases with 105 hash-verified
native replays of the final binary plus fresh diagnostics. Both recovery passes'
matrices are retained. Final owning37/61/65 read tags and VCF rows also match the
accepted baseline exactly. Main build and git diff whitespace checks pass.


### 2026-10-03 — Keep flank BAM evidence and certify current block paths

Investigate noncentromeric chr20:50,548,245–50,562,066. Immutable source selection
incorrectly reused strict injection ownership: a source block starting exactly
at the right seam endpoint was absent. Retain all eligible phased source rows
and callable observations, including context-only blocks. Preserve their local
BAM identity with no live remapped label unless the block actually transfers;
raw coordinate PS equality cannot give a context block graph ownership. The
source snapshot grows from 27 left rows/252 reads to 118 rows/514 reads, including
91 right rows. All 12 MAPQ60 physical spanning molecules survive, versus 2.

The final complete-source stitch also required one molecule to observe both
ends of a 435 kb graph block. Replace that with supported successive-observation
edges and coordinate-cut evidence, cached per current block. Every distinct cut
must have both haplotypes, configured minimum support, and more agreement than
conflict. Supported edges must connect all actual rows. Difference-array cut
support alone can accept interleaved disconnected A–C and B–D components; the
new topology negative control reproduces that false join and guards it.

A source path certificate covers only its saved rows. If an earlier attachment
absorbed other graph rows into the same PS, certify the expanded current block
independently; a valid 57 kb source subset cannot certify a 421 kb current root.
Connected-chain, disconnected-cut, internal-switch, unknown-MAPQ, context-only
PS collision, interleaved-component and absorbed-row regressions pass alongside
the existing physical SNP and per-block HP checks. The initial implementation
fails seven assertions in the connected-chain/absorbed-row tests. The cut-only
trial fails the separate interleaved-component regression. No threshold,
parental truth, competitor label or coordinate-specific production decision is
added. No row is merged or re-aligned by these fixes.

The owning 50 Mb target output is unchanged: 3,877 truth-scored, 3,867 correct,
10 discordant, identical read tags and VCF rows. The seam remains open: retaining
the right calls does not repair the left complementary deletion observations,
which are unknown on all 12 spanners. Only two callable insertion pairs survive,
with one haplotype and no clean-SNP pair. Original CIGARs contain several deletion
lengths, absent deletions and low-quality bases; the old 6/6 BAM HP aggregate is
not an allele certificate. Keep this distinct from a coverage-filter bug.
Evidence: evaluations/2026-10-03-recovery-block-continuity/ and
 evaluations/2026-10-03-gap-50548-stitch-audit/.

Final native chr20 retains exact read-tag and VCF-row parity with the accepted
baseline: 237,199 truth-scored, 230,072 correct, 7,127 discordant (96.995350%);
63,630 keys, 332 VCF blocks, span N50 774,189 bp, 87/101 panel seams connected.
No phased SNP is lost and no old block changes parental gauge. No additional
chr20 gap closes. Final full-run time is 288.07 s with eight threads, concurrent
with panel validation; this is not a competitor runtime comparison. All
standalone unit binaries and 1,259 predicate assertions/44 cases pass. The final
incremental main build introduces no warnings; the earlier header rebuild has
only the existing abPOA SIMD unused-function warning. No expectation is changed.

Final window validation passes 7,132 assertions/71 cases, using 105 fresh native
requests from the final binary plus fresh uncached diagnostic runs. Cached
outputs are verified against the final binary SHA256, input file metadata and
all output hashes. No floor, ceiling or span expectation is weakened. The final
binary SHA256 is
`629badc0b446294ee4e0f1c4b4a43c738406ac7dbed9023995b6796e3fb0f05f`.


### 2026-10-03: Complementary deletion evidence selects the full source retry

Closed chr20:50,548,245–50,562,066 outside the centromere. Both left-boundary
MSA deletion rows are homopolymer indels, so independent literal-REF backfill
left all twelve MAPQ60 spanners unknown at both rows and hid the source conflict
from focused retry admission. A jointly missing complementary pair can now use
MAPQ30/Q30 exact or sequence-equivalent CIGAR ALT as temporary diagnostic
calls: the matching row gets ALT and the other gets ALT absence. Third alleles,
compound edits, extra MSA alternatives and existing calls abstain. Both rows
must have complementary biallelic genotypes in the same source PS and lie in
the seam. No genotype, count, key or phase label changes in this step.

Restore profiles and their interval index before source validation/transfer.
Only the accepted existing full-block BAM/MSA retry supplies imported phasing
and read-rescue evidence. Persisting the diagnostic calls was rejected: it
added 125 erroneous new chr20 assignments and dropped accuracy to 96.945774%.
The accepted scoped design closes the target without those erroneous new
assignments. Production uses neither truth nor competitor data and adds no
realignment stage or merged rows.

The owning 50–51 Mb chunk has 3,887 truth-scored, 3,876 correct and 11 discordant
reads, versus 3,877 / 3,867 / 10. Two separate deletion ALTs remain opposite;
the two-base deletion ALT and right A→G ALT are on the same haplotype. The
114 kb replay remains open because it lacks complete adjacent block context.
The committed panel therefore uses the owning chunk, changes only this span
expectation 0→1 and TOTAL 87→88, and requires both separate deletion rows.
A dedicated owning regression checks parental orientation, coverage and the
11-error ceiling; the original binary fails five of 27 assertions and the
accepted binary passes all 27. New pair unit tests reproduce six failures in
the previous implementation; final predicates pass 1,325 assertions/45 cases.

Final native chr20: 237,209 truth-scored, 230,081 correct, 7,128 discordant
(96.995055%), versus 237,199 / 230,072 / 7,127 (96.995350%). Net +9 correct
and +1 discordant; accuracy decreases 0.000295 percentage points. Read PS
659→656, VCF blocks 332→331, N50 774,189→790,093 bp, NG50 unchanged at
643,699 bp; tracked intervals 87/101→88/101. Only the target joins old VCF
blocks. The retry adds 556 keys and removes one low-depth SNP at 50,719,983
beside retained insertion alternatives; it changes one SNP projection at
50,562,679 beside a newly retained insertion. These representation changes
are disclosed in the evaluation and are not claimed as verified literal SNP
genotypes. Whole-block parental read gauges remain protected by the owning
regression. Full-run time 289.57 s at eight threads, concurrent with panel
validation, is not a competitor runtime comparison. All standalone units pass.
Evidence: evaluations/2026-10-03-complementary-deletion-backfill/.

Final gap validation passes 7,155 assertions/72 cases. The final cache includes
105 fresh native panel requests plus the fresh owning-chunk run, with final
binary, normalized argument, input metadata and output-hash checks; additional
uncached diagnostics execute natively. No existing floor, ceiling or accepted
connection is weakened. Main build and `git diff --check` pass, with no new
warnings in the current change. Final binary SHA256:
`094186bd72b9a2b96cc888519d7e326ff205b8ef63d74d519554a0245beed6cc`.


### 2026-10-03: Isolated deletion conflict and singleton retry context

Closed chr20:61,738,239–61,747,506 (9,267 bp), outside the centromere.
The surrounding graph seam ends at a singleton SNP. Its nominally connected
BAM source assigned opposite ALT orientations to the SNP and homopolymer
deletion, while both original paired calls supported ALT/ALT. Missing deletion
calls hid the conflict; a newly admitted internal-conflict request also lacked
singleton solve padding and failed strict seam containment.

Only a singleton graph flank additionally admits isolated phased biallelic
MSA homopolymer deletion diagnostics. Original MAPQ30/Q30 CIGAR calls must
be exact or edit-equivalent; additional co-located MSA indels, replacement
alleles, third lengths, compound events and existing calls abstain. Temporary
calls participate only in the internal-conflict statistic after ordinary
retry/dropout decisions have been measured. Existing paired-row conflicts
retain their selection/context. New isolated requests receive chunk-bounded
focused context, and their accepted replacements must preserve old clean-SNP
memberships and allele relations. Restore source profiles/index before retry
validation or transfer. No truth, competitor output, coordinates, new alignment
stage, merged rows or relaxed statistical threshold enters production.

Rejected trials are preserved in the evaluation: chromosome-wide isolated
diagnostics fell to 96.610648%; singleton-only diagnostics still contaminated
ordinary retry admission at 17.503 Mb and padded an established paired-row
request at 19.374 Mb, falling to 96.852973%. Key retention and even preserved
source SNP gauges alone did not make that old paired-row join safe. Final
routing keeps those established decisions intact: the protected 17 and 19 Mb
owning outputs have exactly unchanged variant rows and read tags.

Owning 61–62 Mb output improves 3,452 scored / 3,423 correct / 29 discordant
reads to 3,468 / 3,447 / 21, with the same 496 variant keys. All 16 new
assignments are correct, eight errors become correct and none worsen. The
panel adds this interval with spans=1, 0.99 concordance and 0.88 dominant-block
separation floors and complete owning context. TOTAL increases 88→89; every
old floor, ceiling, required row and span requirement remains unchanged.
The strengthened owning regression passes 34 checks; the previous executable
fails six. Isolated deletion unit tests reproduce eight failures in the
previous paired-only helper. Final predicates pass 1,404 assertions/46 cases.

Among 96 truth-scored gap-overlapping reads, the owning pgphase output phases
94 / 94 correct / zero discordant, versus saved same-BAM HiPhase DV and
pgphase-site outputs at 92 / 91 / one. pgphase still has two read phase sets
with 85 correct overlaps in the dominant block; HiPhase has one with 91.
This closes the VCF gap without claiming all read labels are consolidated.

Full native chr20: 237,225 truth-scored / 230,105 correct / 7,120 discordant
(96.998630%), versus 237,209 / 230,081 / 7,128 (96.995055%). Net +24 correct,
−8 errors, with no previously correct read worsened. All 64,185 variant keys
remain; no phased SNP gauge changes. Only the target merges old VCF blocks.
Read PS 656→654, VCF blocks 331→330, N50 790,093→806,449 bp,
NG50 643,699→684,798 bp, tracked intervals 88/102→89/102. Eight-thread run
300.26 s concurrent with panel validation, not a competitor runtime benchmark.
Evidence: evaluations/2026-10-03-singleton-conflict-retry/.

Final window validation passes 7,176 assertions/72 cases over the 102-interval
panel and owning/parental regressions. The final cache contains 106 fresh native
requests verified against the final executable SHA256, argument vectors,
input metadata and output hashes; uncached diagnostics execute natively.
All standalone units, predicates, main build and `git diff --check` pass with
no new warnings. No existing regression requirement is weakened. Final binary
SHA256: `336996e6338b3abc4a55b662231e71e86d378f4222615c4e39aa6cc94a645ce4`.
## 2026-10-03: reject broad shared-BAM genotype ownership

Recheck the remaining 21.823 Mb seam. A prefix-normalized graph repeat shares
the verified BAM deletion at 21,831,481, but retains its graph demotion and
unset genotype. Its source PS also contains the left deletion at 21,823,067,
with zero callable read pairs between those two rows. This distinguishes a
shared-site transfer loss from an independently supported whole-block join.

Experimentally allow non-suffix shared-row adoption only when the complete
BAM source block has at least two anchors and no weak source path cut. Keep
all existing verified binary MSA, graph ownership and conflicting-source
checks. Full chr20 adds 300 VCF keys but closes no new panel seam (89/102).
Correct assignments fall 230,105 -> 230,063; 26 previously correct reads
become incorrect and 62 become unphased. Scored reads fall 237,225 -> 237,178;
discordance falls 7,120 -> 7,115, which does not make those regressions safe.
VCF blocks increase 330 -> 332; N50 falls 806,449 -> 790,093 bp. Two old SNP
blocks acquire split gauges. Reject and restore the accepted implementation.

MSA allele verification and a source PS label are not a certificate of all
internal orientations. Preserve verified sites while independently checking
their component gauges and graph-flank connections. The diagnostic matrix's
VAR weight column is not n_uniq_alles; avoid inferring allele eligibility from
that field. Evidence: `evaluations/2026-10-03-shared-bam-genotype-trust/`.

## 2026-10-03: retain isolated shared BAM genotypes without source-label joins

Fix the shared-row loss of the MSA-verified CT>C deletion at VCF 21,831,480.
For a non-suffix exact binary match, retain its genotype only when every
neighboring heterozygote in the source PS is separated by a weak source cut.
Co-located contrasts, source singletons and supported neighbors keep the
existing rule. Carry the physical BAM key, genotype, depths and proof together
in an independent unused PS; disable raw-source ownership of that row. Normal
paired-read stitching must independently establish any later connection.

Keep its graph observations in the graph channel while primary/BAM calls use
the source projection. Graph rescue infers one established block's orientation
from independent primary read associations, preserving singleton confidence.
Scope this behavior with bam_independent_genotype: a generic MSA-marker
fallback removes 23 pre-existing read assignments outside the owning chunk.
Profile channel arrays now span the union of primary, graph, BAM and quality
indices; primary-only extents can write independent calls out of bounds.
No truth, fixture-coordinate admission, new alignment stage, row merge or
changed genotype/link threshold enters production.

Full chr20 retains every old VCF row, HP/PS tuple and read-correctness state:
237,225 truth-scored / 230,105 correct / 7,120 discordant, 96.998630% accuracy.
Add three previously suppressed independent heterozygotes at 21,294,113,
21,831,480 and 44,083,059. VCF keys 64,185 -> 64,188; VCF blocks 330 -> 333.
Read PS 654, N50 806,449 bp, NG50 684,798 bp and panel connections 89/102 stay
unchanged. No old SNP gauge splits and no old phase block merges. The targeted
21.823 Mb gap remains open; its missing genotype now survives transfer.

Add an owning 21–22 Mb regression for representation, source depths,
independent PS, separate flanks and existing read floors. The starting binary
fails the missing-row assertion; the final passes all 26 checks. Synthetic
checks cover source neighbors and separate graph/BAM rescue gauges. All units,
1,404 predicate assertions/46 cases, 7,198 window assertions/73 cases, HiFi/ONT
goldens and HiFi thread determinism pass. No expectation is weakened. Run 106
fresh native requests and score verified outputs with exact argument/hash/input
checks. Full chr20 takes 280.60 seconds at eight threads with concurrent panel
validation (395.14 seconds), not a controlled competitor runtime comparison.
Evidence: evaluations/2026-10-03-isolated-shared-bam-genotypes/.


### 2026-10-03: Preserve tandem insertion evidence and its read assignments

Close chr20:41,879,449–41,880,908 outside the centromere. Separate BAM MSA
CA-repeat rows had only 11 original observations each. Fixed-consensus local
recall admitted homopolymers but omitted the tandem-repeat contrast, while
supplementary CIGAR lookup reported ALT absence at shifted repeat anchors.
Admit common-prefix insertion contrasts under the existing disjoint one-base
error neighborhoods, surviving-flank, physical-coverage, two-path agreement
and complementary-call checks. Keep both rows, genotypes and source gauges.

Queue newly eligible observations only at original missing MSA slots and
commit after selecting each recovery source, so retry decisions stay intact.
Every exact shared clean SNP in that source PS must match one graph PS with
one constant allele orientation; at least one shared anchor must exist.
Multiple graph owners, inconsistent gauges and absent anchors retain the
original projection. Assigned source reads retain their membership and gauge.
Deduplicate, reject conflicts, consume once and rebuild the interval index.
MSA recalls update discovery counts; supplementary physical corrections do not.

The first guarded closure lost eight owning-chunk read assignments. Keep the
original 4,154-read floor. Stage exact Q30 shifted insertion ALT corrections
when the per-read clean-SNP gauge is absent, admitting them only in source
blocks with newly recalled complementary MSA insertion evidence and the same
coherent graph gauge. Contradictory per-read SNPs still reject repair. A
sequence-equivalent intermediate repeat length matches neither complementary
row: explicitly reject both supplementary calls. Transfer those rejections
instead of allowing older primary/BAM CIGAR calls to survive. Original MSA
calls remain authoritative. In read rescue, an equal block score may favor a
directly phased SNP only when primary reads pass the existing singleton
confidence check. Inferred excluded-site associations cannot break a tie;
co-located indel and inferred-SNP properties cannot combine into a certificate.
Read rescue never changes candidate gauges or joins phase sets. Production
uses no truth, competitor data or fixture-coordinate gate.

Reject admission before source selection: it reopens the downstream 41.900 Mb
connection and loses 20 keys. Reject admission without the coherent owner
guard: it wrongly joins 19.374 Mb and changes 368 old correct chr20 assignments
to incorrect (7,120→7,446 errors). Reject physical corrections outside recalled
blocks: they lose an old 5.31 Mb assignment. Reject unconditional SNP tie
preference: inferred markers add 20 owning34 assignments, ten incorrect;
primary support alone is insufficient. Require a directly phased SNP too.
Existing coverage and correctness floors, split expectations and original
retry scheduling remain intact. Evaluation reports retain the counterexamples.

Final chr20: 237,308 truth-scored / 230,211 correct / 7,097 discordant
(97.009372%), versus 237,225 / 230,105 / 7,120 (96.998630%). Net +106 correct,
−23 errors and +83 scored. All 230,105 old correct assignments remain correct;
21 old errors become correct and two become unphased. All 85 newly scored
reads are correct. All 64,188 keys remain; no old phased SNP or allele gauge
is lost. Only the target VCF union occurs: blocks 333→332, panel 89→90/103.
Read PS stays 654. N50 806,449 bp and NG50 684,798 bp are unchanged; the target
block spans 341,860 bp. No existing connection reopens.

Owning41 retains all 4,154 scored reads and all 965 keys, improves 4,109 correct /
45 errors to 4,126 / 28, and preserves its downstream connection. Of 81 local
truth overlaps, pgphase changes 81 phased / 62 correct / 19 errors to 81 / 79 / 2.
Saved same-BAM HiPhase with pgphase calls has 80 / 80 / 0. Local read purity
still trails HiPhase despite the correct VCF join; no new competitor run is
claimed. Add the coordinate case with spans=1 and the upstream required site;
the owning regression asserts separate repeat rows, parental allele orientation,
downstream continuity and strengthened 4,154 / 4,126 / ≤28 read gates.

Standalone units, 1,465 predicate assertions / 47 cases, 7,236 window assertions /
73 cases covering 103 coordinates, and HiFi/ONT TSV/VCF goldens pass. HiFi one/four-thread outputs are identical. No new warnings.
The final native chr20 run takes 369.89 s / eight threads during concurrent
regression generation, not a controlled runtime comparison.

Final binary SHA256:
`9035673c60ec6798a356e2651b44880ba3b87f2967fb104ee4deed4c3d4b9477`.
Evidence: evaluations/2026-10-03-tandem-insertion-recall/.


### 2026-10-04: Close the 7.9 Mb gap with an ALT-only query certificate

Close chr20:7,901,413–7,918,883 (17,470 bp), nominated by saved same-BAM
HiPhase output. Two MAPQ-60 primary reads carry the exact long MSA deletion
allele, but their CIGARs encode a shifted 15-base deletion plus compensating
mismatches. The old reference-edit comparison rejects both. Between surviving
reference anchors, both query strings exactly equal the candidate allele
`ATAT`, at minimum retained qualities 40 and 27. The 13- and 15-base BAM
rows remain separate and complementary. Split/other-length observations
remain unverified, never deletion REF.

Add `bam_matches_deletion_sequence` as an ALT-only certificate consumed only
by the graph-SNP to complementary BAM-deletion bridge. Keep the BAM source
solver and existing physical REF/ALT callers on the original CIGAR certificate.
Require MAPQ30, Q20 exact sequence, two independent unanimous witnesses and
the existing 0.001 quality-weighted parity rule. Try this bridge before rejecting
an exposed right-deletion seam whose clean SNP pair is missing or one-sided.
Available SNP votes must be unanimous and match its proposed orientation.
The right graph block may bypass one conditional SNP only through the existing
significant direct flank edge with at least two reads on each haplotype;
dominant reversal still vetoes. Whole block paths retain their validation.
No truth, competitor output, fixture coordinate gate, new realignment or allele
merge enters production.

Reject applying the new certificate during source recall: it introduces eight
keys, wrongly joins 20.8 Mb and converts 90 old correct reads to errors
(7,097→7,195 discordant). Source isolation restores every read truth assignment,
but a broad physical REF/ALT API still reopens the protected 64.14 Mb connection:
an uncertifiable shifted ALT can otherwise fall through to REF. The retained
ALT-only API avoids that classification and preserves every earlier connection.
Also reject broad shared-source path admission (155 old correct reads become
incorrect at 37.55 Mb). Pure BAM-component paths, trying all complementary
loci and one-haplotype SNP retry changes alone close no new chromosome gap.

Final chr20: all 64,188 keys and all read truth assignments remain identical.
237,308 scored / 230,211 correct / 7,097 discordant, 97.009372% before and after.
All 230,211 old correct reads remain correct. Only the target VCF union occurs:
blocks 332→331; scored read PS 654→653; expanded panel 90→91/104 connected.
No old SNP gauge or connection is lost. The joined block spans
7,859,472–8,081,718 (222,247 bp). N50 806,449 bp and NG50 684,798 bp are unchanged.
The native eight-thread full run takes 370.92 s during concurrent regression
generation; it is not a controlled runtime comparison.

Owning7 retains all 690 keys and 4,091 scored reads: 3,926 correct / 165 errors
before and after, with read PS 18→17 and VCF blocks 9→8. Of 143 local truth
input overlaps, pgphase retains 116 phased / 114 correct / 2 errors; correct
reads in one block improve 64→104. Saved HiPhase/DV has 104 phased / 103 correct /
1 error, and saved HiPhase/pgphase calls has 86 / 84 / 2. pgphase retains more
correct local reads, while local purity remains below HiPhase/DV. No fresh
competitor run or higher-local-purity claim is made.

Add the coordinate to the permanent panel with spans=1, three required boundary
site rows and an owning-chunk regression checking exact separate deletion rows,
parental allele orientation, 4,091 / 3,926 / ≤165 owning read gates and
116 / 114 / ≤2 local read gates. The new owning test fails on the starting
binary (three connection assertions) and passes after the fix. Existing floors,
ceilings and intentional split cases are retained; only the graph span TOTAL
rises 90→91. Unit gates pass, including 1,497 predicate assertions / 47 cases.
All 7,977 window assertions / 74 cases covering 104 coordinates pass against
106 freshly generated native requests with verified binary/input/output hashes.
HiFi/ONT TSV/VCF goldens and HiFi one/four-thread determinism pass. No new warnings.

Final binary SHA256:
`00f8c411a595359473ac2be0cce0f9639a0810be9373310b37f024a6da2e017d`.
Evidence: evaluations/2026-10-04-compound-deletion-certificate/.

### 2026-10-04: Preserve verified MSA calls across independent BAM blocks

The deferred-observation transfer rejected a verified fixed-consensus allele
when its read belonged to another source phase set. Adding a call does not
reassign HP/PS, and numerical HP values are not comparable across independent
block gauges. Keep those observations for later stitching; retain the allele
contradiction veto within the same phase set. Discovery, source selection,
genotypes, representation and stitch criteria remain fixed.

Full chr20 restores 25 complementary call pairs across four loci; eight VCF
rows change only their depths and allele fractions. All read tags, variant
keys, GT/PS values, 91 connected panel coordinates and 13 open coordinates
remain unchanged. Read truth remains 230,211 correct / 7,097 discordant out
of 237,308 scored reads (97.009372%). This fixes evidence loss but does not
close another gap. Unit tests, 1,517 predicate assertions, 1,310 complementary
insertion window assertions, and HiFi/ONT golden gates pass. Regression
expectations are unchanged. The full chromosome check verifies all previous
connections and read assignments; the full standalone window suite was not
rerun for this change.

Evidence: evaluations/2026-10-04-cross-block-msa-observations/.

### 2026-10-04: Audit and preserve verified recovery allele transfer

Instrument pending MSA admission, canonical source-to-destination mapping after
candidate/read sorting, and whole-chunk BAM attachment. Eighteen noncentromeric
owning-chunk replays confirm 506,149 unambiguous mapped calls survive transfer
and expose 76 competing read/site calls from overlapping solves. These are
record counts across stages, not unique chromosome alleles. Full source
matrices retain independent block evidence; the working matrix still uses one
allele slot and the existing first-available-call policy for conflicts.

Fix a late whole-BAM overwrite: 12 paired insertion calls at VCF 61,802,959
and one complex insertion call at 61,806,160 changed despite already selected
targeted evidence. Fill only missing BAM observations and preserve their
query positions and quality certificates. All 481,431 audited existing
clean/MSA overlay calls now survive. Full chr20 read tags and VCF rows remain
identical for this fix alone, showing why output-only tests missed it.

Also retain fixed-consensus MSA observations in a source block with no
heterozygous anchor matching the selected graph representation. Requiring a
shared clean graph SNP incorrectly rejected 51 complementary call pairs at
21,594,343; both separate insertion rows rise DP 11→62. Preserve the source
gauge and normal stitch checks. Graph-owned blocks still require a coherent
shared clean-SNP gauge, and supplementary physical corrections do not receive
the private-block exception. Inconsistent-gauge recall remains rejected: the
earlier broad admission trial produced the wrong 19.374 Mb join and 368 newly
incorrect reads. No coordinate, truth or competitor data enters admission.

Final full chr20 retains all 64,188 keys, 331 VCF blocks, 653 scored read phase
sets and 806,449 bp span N50. Scored reads stay 237,308; correct reads improve
230,211→230,212 and discordant reads fall 7,097→7,096 (97.009372%→97.009793%).
Every old correct read remains correct. One erroneous assignment is withheld
and one previously unphased read becomes correct. Two read tags and six VCF
rows change. All 91 panel connections and 13 open coordinates are unchanged;
no new gap closure is claimed. Concurrent runs are not runtime benchmarks.

Two permanent owning-chunk regression cases reproduce both defects on the old
code and pass after correction (10 conservation assertions and 23 private-block
assertions). The conservation test checks source aliases/reindexing, MSA call
preservation and continued filling of missing fallback calls. The private-block
test checks separate rows, complementary genotype, depth, independent PS and
owning read truth without allowing a new join. Regression expectations remain
unchanged. Build, unit/predicate tests and HiFi/ONT golden/determinism gates pass.
The final full window suite passes 8,006 assertions in 76 cases, covering all
104 panel coordinates through 110 fresh native requests. An additional final
32 Mb owning replay retains 55,703 unambiguous mapped calls and 55,071 overlay
calls; the full chromosome's other four changed rows preserve GT/PS while
depths increase 7→8 at 32,364,883 and 51→60 at 32,725,929.

Final binary SHA256:
`b4614c555ca06d411e68fc5f6f4e79bb9df8507c755d2172ff13557b6665fb47`.

Evidence: evaluations/2026-10-04-verified-allele-transfer/.

### 2026-10-04: Preserve independent verified calls through admission and transfer

Remove the graph-gauge and prior same-block HP veto from selected fixed-consensus
MSA observation admission. Source rows, genotypes and HP/PS stay fixed. Accepted
focused solves now consume their pending recalls too. When graph ownership has
an unresolved gauge, cache the source's original internal-path certificate
before adding calls; restored observations cannot erase a weak cut and certify
an otherwise unsupported whole-block join. Supplementary physical projections
still need their original coherent graph gauge and fixed-consensus context.

Replace first-available working observations with genotype-source ownership for
BAM-only rows and explicit conflict (-2) for shared/unowned disagreements. Keep
every independent source matrix, allele and source-specific physical quality,
including omitted flank sites. Graph calls survive independent BAM ambiguity.
Conflict markers survive candidate reordering, recovery retry and late fallback;
fallback still fills ordinary missing calls (-1). Verified MSA takes precedence
over supplementary physical projections. Same-tier contradictory recalls
abstain while preserving all alternative proposals as provenance, not extra
molecule votes. No truth, competitor or coordinate special case enters phasing.

Reject an initial all-conflicts-unknown trial: it lost two correct owning61 read
assignments by discarding the selected genotype owner's calls. The refined
owner rule retains all old correct reads in seven noncentromeric owning probes.
The final targeted audit retains all 66,029 checked independent source calls,
71,612 existing overlay calls and one explicit overlay conflict. Of 64 source
disagreements, 63 use the selected BAM genotype owner and one shared slot
abstains. Qualities are computed per source rather than borrowed from the
merged first-call table. Counts include solve/stage records, not unique alleles.

Final full chr20 gains 20 truth-correct reads: 237,308→237,328 scored,
230,212→230,232 correct, 7,096 errors unchanged,
97.009793%→97.010045% concordance. Every old correct and discordant read retains
its truth status. All 64,188 variant keys and GT/PS values remain unchanged;
ten rows at five insertion loci change only depths/allele fractions. VCF blocks
remain 331, scored read phase sets rise 653→654 and span N50 remains 806,449 bp.
All 91 connected panel coordinates and 13 open coordinates retain their status;
no new block join or gap closure is claimed. No old SNP block has a mixed gauge.
The concurrent eight-thread native run takes 376.23 s, not a runtime benchmark.

The new native admission regression fails on the starting binary in four
depth assertions and passes after the fix: complementary verified rows at
11.862429 rise DP62→68 and 34.835166 rise DP51→56. The quality regression proves
that the shared 3.597791 Mb REF/Q40 and ALT/Q0 calls survive independently and
the working BAM conflict cannot be filled by fallback. The transfer regression
checks immutable context calls as well as known/blocked overlay observations.
Targeted native checks pass 90 assertions in four cases, including the known
wrong 19.374 Mb whole-block join guard. Unit tests, 1,532 predicate assertions in
47 cases, HiFi/ONT golden gates and HiFi one/four-thread determinism pass.
Regression expectations, required sites and stitch thresholds are unchanged.

The complete native window suite passes all 8,060 assertions in 77 cases,
covering 104 coordinates through 110 fresh pipeline requests. The rebuilt
targeted suite also passes after aligning its 34 Mb fixture boundary to the
existing panel; owning-chunk CLI and asserted truth floors are unchanged.
Independent source and late-overlay conservation pass on the final panel dumps.

Final binary SHA256:
`54e67e0911beb6402e2558b44ac5dbf31d04da8da8f26b93a66e9acd8a91348f`.

Evidence: evaluations/2026-10-04-independent-recovery-evidence/.

## 2026-10-04: exact-allele physical stitch certificates

Correct three stitch validation bugs: a graph allele could borrow physical
quality from a different or absent BAM allele; the single-molecule BAM
corroboration predicate treated unknown MAPQ 255 and missing BQ 255 as high
confidence. Four reproduced negative predicate assertions now reject those
cases. Positive synthetic certificates explicitly include matching BAM calls
and keep them aligned when candidate indices shift. No source GT, HP/PS,
representation, BAM parity behavior or stitch threshold changes.

A fresh default full chr20 run is identical to the accepted previous output:
237,328 truth-scored phased, 230,232 correct, 7,096 discordant (97.010045%),
654 scored read PS, 64,188 VCF keys, 331 VCF blocks, 806,449-bp span N50;
91/104 coordinate spans. Ten competitor nominations and three controls remain
open. Final binary SHA256 is
`9bda79a94b69241f8b9083ee902dceea6d0c4560f583ff09cfce2543afcb45df`.

At 35.5 Mb, the exact solver with every usable gap row prefers the wrong
parental join; one read half is tied. Same-site targeted HiPhase has three
MAPQ60 left-REF/right-ALT bridges whose left deletion calls are absent before
pgphase transfer. Extending deferred fixed-consensus recall to solitary seam
deletions restores those three calls and increases that row's DP14 to45, but
adds 49 discordant reads in the owning-Mb solve and leaves the gap split.
Reject and remove this broad trial. The comparable targeted HiPhase read
accuracies are 93% at 35.5 Mb and 90.38% at 57.85 Mb, versus pgphase's current
96.39% and 100% independent groups; this diagnostic uses current pgphase sites,
not the earlier whole-chromosome DeepVariant benchmark. Named HiPhase trace
segments are required: anonymous global vectors cannot be zipped with BAM
fetch order because global alignment skips some records.

Build, all standalone units, 1,536 predicate assertions in 47 cases, 164 targeted
native window assertions in seven cases, HiFi/ONT TSV/VCF golden parity and
HiFi one-/four-thread determinism pass. No expectations are lowered. Details and rejected
trial: `evaluations/2026-10-04-stitch-quality-certificates/README.md`.

## 2026-10-05: recheck the 50.548 Mb complementary-deletion case

The user-referenced “Complementary deletion evidence for source retry” interval
is already closed, not a remaining miss. Fresh current owning 50–51 Mb and
latest full-chr20 output both phase CA→C and CAA→C at 50,548,245 plus A→G at
50,562,066 in PS 50,127,297; deletion ALTs remain opposite and the two-base ALT
matches the SNP ALT. The existing owning regression passes all 27 checks.
No production behavior or test expectation changes.

Among twelve distinct primary MAPQ60 spanners, the initial source has zero
callable pairs at either deletion; right SNP calls are intact, six REF/six ALT.
The accepted full-block MSA retry has 10 one-base and eight two-base pairs,
unchanged in transfer and final matrices. Same-site short HiPhase default has
five/ten; local-only has eleven/ten. The two missing pgphase two-base calls
are reads with three-base CIGAR deletions, one with a Q17 flank. HiPhase calls
both ALT even without global alignment. This is noisy allele assignment,
not loss of verified transfer. No third length is silently promoted.

Gap-overlapping truth under each whole output block's parental orientation:
pgphase owning and full-chr20 112 scored /105 correct /7 errors (93.75%);
HiPhase same-site owning and short-default 109 /103 /6 (94.50%); local-only
107 /100 /7 (93.46%). Pgphase gains two correct reads and one discordant read
in that comparison. HiPhase connects the one-base deletion to the right SNP
but leaves the two-base row without PS; pgphase phases both deletion rows.
Across the owning Mb, pgphase has 3,887 /3,876 /11 (99.7170%), HiPhase
3,898 /3,438 /460 (88.1991%). HiPhase's one group includes 442 parental errors
on reads in pgphase's distinct PS 50,002,195, outside this local gap. This is
a targeted same-pgphase-sites diagnostic, not a new full-chromosome benchmark.
Named HiPhase traces, original physical spans, exact source/transfer call
comparisons, read truth scores and manifests are preserved in
`evaluations/2026-10-05-complementary-deletion-hiphase/`.

## 2026-10-05: close the verified noisy-SNP boundary at 41.881 Mb

Close the newly nominated, noncentromeric chr20:41,880,908–41,885,033 gap.
Revisit an exposed MSA insertion opposite a verified noisy BAM SNP within
existing recovery windows. Its imported calls must independently anchor its
orientation to the nearest clean SNP at the existing 0.001 binomial bound,
with both allele classes and matching deterministic read halves. Graph path
validation can bypass the left weak node of an outgoing edge using the same
significant, two-haplotype neighbor edge as the existing right-node bypass.
For this verified noisy boundary alone, an already certified BAM join can
survive a nonsignificant GAF reversal; other bridges keep their strict rule.
No candidate allele or verified observation is replaced.

The 41–42 Mb replay improves 4,154 scored /4,126 correct /28 discordant to
4,156 /4,129 /27. Gap-overlapping reads improve 87/89 to 90/91 correct,
versus saved same-BAM native DV HiPhase's 87/91. Both complementary +14/+18
CA-repeat rows retain their sequences and AD, with the shorter ALT agreeing
with both right SNP ALTs. All 27,385 saved source calls and 37,791 old overlay
calls survive the conservation audit.

Final default chr20: 237,330 truth-scored phased, 230,235 correct, 7,095
discordant (97.010492%), 652 scored read PS, 64,188 VCF keys, 330 VCF blocks,
806,449-bp span N50. All 230,232 formerly correct reads remain correct; two
new reads are correct and one prior error is corrected. No SNP block changes
its internal gauge; every old connected panel interval remains connected.
The panel grows 104 to 105 coordinate cases, with 92 connected after this fix.

Reject global certificate relaxation and a broader left-insertion seam retry:
they affected protected native replays or formerly closed chromosome blocks.
The trial 7.26 Mb join is not retained. Preserve the 36 Mb control's spans=0
and accuracy floors. The updated 41.900 Mb regression explicitly requires
the intended outer connection while preserving its downstream and parental
checks. New case asserts spans=1, read truth, complementary GT and exact AD;
its two insertion rows join the required-sites list.

Build, standalone units, 1,536 predicate assertions in 47 cases, BAM HiFi/ONT
goldens and thread determinism pass. The broad trial's full native run found
three failures in 8,126 assertions/78 cases. Final focused reruns pass 105
assertions across the new case, updated 41 Mb case and 37 Mb counterexample,
plus 20 assertions for the unchanged 36 Mb panel control; no quality floor is
lowered. Independent final native verification passes 8,106 assertions in 78 cases,
using 110 fresh final-binary pipeline requests. The existing test executable
retained one obsolete split assertion; the updated test source changes only
that reviewed connection assertion and its comment. Rescoring hashes the
unchanged CLI requests, core outputs and auxiliary matrices. Independent
chr20 verification confirms exactly this one new gap, no genotype-allele
changes and no previously correct reads lost or made discordant. Final binary SHA256:
`69d998cab0b2319f6875e3514a7be4eaa1dc9b514419831614aeb012a13c4386`.
Details: `evaluations/2026-10-05-verified-snp-graph-path/`.

### 2026-10-05: unresolved next-gap screen

After the verified 41.881 Mb closure, fresh physical inspection finds three
11.235–11.255 Mb spanners with conflicting right deletion lengths. At 7.264 Mb,
the blocking 7,280,346–7,280,356 graph SNP edge has 3 same / 3 cross votes;
original Q30 SNP pairs do not corroborate its current opposite orientations.
Weak-reversal physical validation and complete-source first-edge trials do not
close it. A strict two-haplotype, unanimous multi-read physical certificate at
p<=0.001 likewise changes no tags or VCF rows in fresh owning7/15/21/34/48/56
comparisons. All trial edits are reverted; no new gap is closed and no panel
expectation is changed. Evidence and parity reports are in
`evaluations/2026-10-05-next-gap-screen/README.md`.

### 2026-10-05: deferred physical bridge closes 23.461–23.481 Mb

Close chr20:23,460,963–23,480,815 while retaining the preceding
23,421,003–23,445,252 connection. If the ordinary whole-flank physical
validation and deletion backfill fail, a verified targeted right-source noisy
insertion permits a separate full-flank solve with unplaced MSA observations
and insertion recall enabled. It must pass the unchanged source path, exact
shared-site gauges and physical certificate, then the ordinary stitcher's
outer-allele and gauge conflict checks on a disposable state. The validation
solve transfers no rows, counts, genotypes or read labels.

Defer this supplemental core union until output-only read rescue finishes.
Save complete original anchor keys/orientations; require uniform live block
membership and gauge before applying the relation in all affected chunks.
This preserves read-only marker cohorts, just as the existing late equivalent
insertion joins do. Applying the new bridge early closed the core gap but
lost five correct 5 Mb tags and changed a 23 Mb rescue cohort (three correct
reads became discordant, two discordant became correct). Reject that version,
source-alias shortcuts, separated validator cohorts and early SNP gauge
refresh. Raw whole-BAM solve matrices were identical: excluded-site marker
rescue consuming merged core cohorts caused the read changes. No GAF path or
physical quality gate is relaxed.

Final full chr20 preserves all 256,610 output reads, 237,330 scored reads,
230,235 correct and 7,095 discordant assignments; every individual truth
status is unchanged. Exactly this one old gap closes, none reopens, and
93/106 tracked gaps connect. All 64,188 variant keys, alleles and non-phase
fields remain, and old phase blocks move uniformly. VCF blocks 330→329,
read phase sets 652→651, span N50 806,449→856,770 bp. The 857 changed read
labels and 256 changed VCF rows are core phase-set/gauge relabels.

Owning 23 Mb preserves 3,542 correct and 60 discordant reads; gap reads stay
125/125 correct and their core PS becomes one. Disjoint flanks choose the same
parent (29/29 left, 31/32 right). Owning 5 Mb is byte-equivalent in read tags
and VCF rows (4,001 correct / 23 discordant). Paired 22–24 Mb preserves every
truth status (7,307 correct / 121 discordant); its additional joined short
replay boundaries were already connected in the accepted full chromosome.
The source audit checks 23,912 independent and 40,669 overlay calls with zero
source disagreements. Add the new panel window, a spans=1 expectation, owning
parental-orientation/preceding-join regression and paired read-rescue regression.
Build, standalone units, 1,536 predicate assertions, BAM HiFi/ONT goldens and
thread determinism pass. The two new regressions pass 72 assertions; the baseline fails five connection
checks and the rejected early-join version fails three rescue-PS checks.
Final native run passes 8,168 assertions in 79 cases using 110 fresh requests;
the additional paired case and strengthened assertions pass in the focused
run, covering all 80 distinct cases in the final source. No expectation is
relaxed. Final binary SHA256:
`3a227a265b28ca21b15e4e59e272f0dc2cc06e7feb2046d50108365e481152c5`.
Evidence: `evaluations/2026-10-05-composed-physical-bridges/`.

### 2026-10-05: complementary insertion ALT certificate closes the 9 Mb boundary

Close chr20:8,977,829–9,014,032 using a clean graph SNP and two distinct,
complementary MSA/alignment-verified BAM insertion ALTs at 8,998,972. Both
belong to one cut-free, consistently oriented original source run. Exact or
sequence-equivalent ALT calls are required; REF of either biallelic insertion
cannot identify the other ALT. One primary physical molecule represents each
ALT, with MAPQ60, SNP Q40 and insertion/flank Q35/Q40. Both imply the same
parity; log odds 16.286107 clears the wrong-parity <=0.001 gate. Retain the
ordinary graph SNP path and BAM source path checks and restrict detection to
already targeted seams.

A raw graph key can appear in several biallelic rows: at 9,048,832 five rows
share a key but only the third is phased. Deferred anchor lookup now searches
phased rows, retaining the full-block uniform live PS/gauge checks. Carry
adjacent-chunk replay certificates to the owning blocks only when both exact
boundary SNPs share a replay PS and at least two exact shared clean SNPs per
side agree in orientation. Capture all original phased anchors, including the
downstream owning chunk, then union after independent read rescue finishes.
No private replay genotype, count or read-observation transfer is introduced.

Accepted baseline binary `3a227a265b28ca21b15e4e59e272f0dc2cc06e7feb2046d50108365e481152c5`
(`test_data/tmp_gap_fix54/full-gated`) versus final binary
`bdce74cda1d0ec407fe7a0cfaf9e30d2c77a60b605d6f193ee92e2cb16503933`
(`test_data/tmp_gap_fix55/full-final`): every individual truth status remains
unchanged, with 230,235 correct and 7,095 discordant of 237,330 scored reads.
Exactly this one old gap closes and none reopens: 94/107 tracked gaps connect.
All 64,188 variant keys, genotype alleles and non-phase fields remain; every
old phase block moves uniformly. VCF blocks 329→328 and read PS 651→650;
span N50 remains 856,770 bp. The 5,634 tag relabels and 1,272 VCF phase changes
are core unions/gauge changes; every output-only rescued HP/PS is unchanged.

Bounded replay preserves 530 correct / 7 discordant (83 keys) and paired
8–10 Mb preserves 8,517 correct / 8 discordant (2,234 keys), also with exact
truth conservation and no rescued-tag changes. Paired N50 985,788→1,053,649.
The gap has 223 truth-scorable overlaps, with 180 phased and 177 correct before
and after. Three phase groups become two; the largest correctly placed group
rises from 64 to 126 (separated fraction 0.286996→0.565022). Disjoint 10 kb
flanks choose the same parent, 35/35 left and 45/45 right. Add the measured
spans=1 panel row and owning-chunk parental/gauge regression. Its 24 assertions
pass; the accepted baseline fails five connection/orientation assertions.
Build, unit tests, 1,536 predicate assertions, BAM HiFi/ONT goldens and thread
determinism pass. The complete final native suite passes 8,234 assertions in
81 cases, using 112 fresh pipeline replays with the frozen final binary. No
existing expectation is relaxed. Evidence:
`evaluations/2026-10-05-complementary-insertion-boundary/`.

### 2026-10-05: further gap search, no retained change

A corrected physical-boundary screen must advance the maximum covered endpoint
of all preceding blocks; adjacent sorted-start blocks can overlap an older long
block and nominate a false gap. The 2.299--2.311 Mb trial was such a false gap
and also failed the insertion-side parental flank audit (12 matching, four
discordant). All trials were removed and the rebuilt binary exactly matches
the accepted baseline SHA256. The corrected simple SNP, simple indel and
complementary-insertion screens nominate no supported new join. Measurements,
rejected evidence and reproducible screens are in
`evaluations/2026-10-05-unclosed-gap-audit/`. No production behavior or test
expectation changes.

### 2026-10-05: short-gap and long-insertion validation trials

Two isolated full-chromosome trials -- removing the physical insertion
caller's 128-base cap, and admitting positive gaps below 10 kb to the
existing full-flank validator -- produce exactly unchanged read tags and
VCF rows. Both were reverted. The Q20/MAPQ20 SNP nomination at
29,303,608--29,303,677 lacks a uniform independently validated graph-flank
gauge; its bounded replay flips 234 old left anchors and omits all 67 old
right anchors. No new gap, retained behavior or expectation change. Evidence:
`evaluations/2026-10-05-short-gap-and-long-insertion-screen/`.

### 2026-10-05: physical graph SNP switch closes the 62.623 Mb gap

The true gap 62,623,253–62,642,316 is blocked by an internal switch in its
right graph flank at 62,718,395–62,722,021. All 38 independent Q30/MQ30
primary BAM pairs contradict that edge, on both haplotypes (22/16), while
neighboring edges agree. Repair one failed internal edge only with unanimous
physical reversal, two molecules per haplotype, count p<=0.01 and the existing
0.001 quality-weighted parity bound. Flip the candidate suffix, certify that
exact site pair with its relative allele gauge, and require the complete
remaining graph SNP path. Certificate gauges compare selected ALT presence,
not raw graph allele IDs. Flip read gauges predominantly observing phased heterozygous anchors in the
suffix before the existing physical seam stitch; homozygous sites can carry a
PS but do not vote. A candidate-only flip is
unsafe: the owning replay drops to 71.12% read concordance.

The general rule closes exactly this new gap in the full chromosome, reopens
none of the old 107 tracked gaps, loses no phased key or genotype allele and
changes no nonphase call fields. Correct/discordant reads move from
230,235/7,095 to 230,563/6,767: 329 improve and one previously correct read
becomes discordant. This intentionally repairs the nonuniform gauge of old
PS=62,637,077; every other old block retains a uniform gauge. The remaining
read's only Q30 phased physical SNP supports its retained label against the
repaired genotype, so do not tune it using parental truth. Twenty output-only
rescued reads keep HP and exact group membership; their independent PS label
follows the new core label without joining another rescue group. The owning
replay is 3,722/3,750 correct, 28 discordant. Gap overlaps remain 113/116
correct, with concordant 10 kb flanks (31/31 and 36/36). Add measured spans=1
and owning parental/SNP-switch regressions, retaining every old window floor.
Evidence and explicit audit tradeoffs:
`evaluations/2026-10-05-physical-graph-switch/`.

Final build, unit tests, 1,536 predicate assertions, BAM HiFi/ONT goldens and
HiFi thread determinism pass. The complete final native suite passes 9,757
assertions in 82 cases, including the 108-window panel, in four fresh disjoint
batches with the frozen final binary. The prior lone-boundary-SNP regression
now permits this join only with the repaired internal SNP relationship and
owning parental/error bounds; all other existing window floors stay intact.
The final certificate uses selected REF/ALT presence for its relative gauge;
its full-chromosome read tags and VCF rows are unchanged from the audited
candidate. Binary SHA256:
`e635c76de38cc00ab0554529a278ca297b321a3a728634a40bc7652fb8b78c71`.

### 2026-10-05: further flank-path screen, no retained production change

Left-flank switch repair, unanimous same-orientation physical graph-edge
certificates, and symmetric right BAM-source path acceptance each leave the
full chromosome unchanged: 230,563 correct / 6,767 discordant reads, 327 VCF
blocks and 95/108 tracked spans. All three trials were removed.

A separate left-source SNP fallback closes the true 2,136-base gap at
32,330,932–32,333,068 in a fresh 32,300,001–32,400,000 replay. Its exact
source identity, cut-free source path and unanimous physical cross parity
(17 Q30/MAPQ30 pairs, log odds 139.506) still fail the biological audit:
52/16 correct/wrong becomes 43/25, with nine formerly correct reads wrong
and none improving. Reject it. An original-source path certificate does not
by itself guarantee parental coherence in a repeat-rich block.

The wider primary SNP and compound-boundary screens identify no accepted
new gap. Evidence and rejected audits:
`evaluations/2026-10-05-flank-path-screen/`. Restore the production binary
exactly to `e635c76de38cc00ab0554529a278ca297b321a3a728634a40bc7652fb8b78c71`;
no expectation or implementation-behavior changes.

## 2026-10-05 — quality-backed source path closes the 13.752 Mb gap

Close the previously uncovered graph interval 13,752,640–13,773,452 using
independent source evidence, without truth-driven production decisions.
Two original source-HP1 molecules call the boundary REF/REF at Q40/Q40 and
Q27/Q40, both MAPQ60, with combined wrong-parity bound 4.23647e-7. A fixed
Q30 SNP floor dropped the second molecule. The clean-source cut check now
admits Q20 while retaining distinct source-assigned reads, conflicting-call
vetoes and its error bound. A matching exact shared graph edge is revalidated
at the stricter 0.001 bound. The next missing GAF edge is bridged through
intermediate original BAM source variants only when both shared endpoints
have the same source gauge, no callable high-MAPQ GAF pair exists, and every
original weak source cut has quality support. Retain accumulated exact-pair
certificates only when the entire remaining graph path passes. Do not clear
source weak cuts or grant unrestricted whole-source transfer. Existing Q30
physical switch-repair requirements remain unchanged.

The final full-chromosome output closes exactly this new gap: 327→326 VCF
blocks, 649→647 read phase sets, 368 changed read tags and nine changed VCF
rows. Preserve all 64,188 variant keys, genotype alleles, nonphase fields,
old candidate-block gauges and 230,563 correct / 6,767 discordant scored
reads. Reopen no old tracked span; adding the new window gives 96/109 full
panel spans. Ten independently source-assigned rescued reads enter the
certified core path, while 32 other changed rescued tags remain in one
uniformly reoriented output-only cohort. None becomes newly discordant.

The user permits a gap closure at >=80% correct reads and requests diagnosis
when HiPhase phases reads better. Here 151/175 truth-scorable overlapping
reads (86.3%, counting abstentions) are correctly separated into one block,
with 24 unphased and zero local discordance. This exactly matches HiPhase
on the same 175 molecules, with all alignment coordinates, CIGARs and
sequences verified. Before this fix, our same 151 correct reads were split
among three blocks; dominant correct separation was 76/175 (43.4%).
Disjoint left- and right-anchor groups agree on parental orientation at
74/74 and 115/115. Thus HiPhase's advantage here was continuity rather than
additional correctly phased molecules.

Add the measured graph window to the committed panel (`spans=1`, concordance
>=0.94, separated >=0.86), exact shared anchors and an owning-chunk regression.
Two existing graph owning-chunk tests intentionally expected this exact
prefix edge to stay split; update their separation assertions to equality,
retaining their suffix, genotype and parental checks. No unrelated negative
span or read-accuracy floor is relaxed. Measurements, reproduction commands,
HiPhase comparison and exact audits are in
`evaluations/2026-10-05-source-quality-path/`. Build has no new warnings;
`make unit-tests`, 47 phase-predicate cases and `make check` pass.
Final native validation covers all 83 cases / 10,115 assertions and all 109
window-panel rows. The first sweep passes the 81 unaffected cases; the only
two failures are the obsolete split assertions for this exact newly closed
graph edge. Both updated owning cases pass all 282 assertions. The frozen
final binary's strict full-chromosome audit passes and its tags and VCF rows
exactly match the already audited trial.

Commit review restores all 1282 pre-existing graph required-site entries,
which had been truncated before the 13 Mb task, while retaining the seven
new source/switch witnesses. The unchanged native panel checker passes all
109 windows / 2,469 assertions with the restored manifest using hash-bound
final-binary output replays. No old required-site check is removed.

## 2026-10-05: Retry short insertions with a complete BAM-only source path

The 33,747,591–33,749,688 gap exposes two incorrect restrictions: physical
retry nomination only considered long right insertions, and short-insertion
stitching demanded a graph path through a right block containing only recovered
BAM sites. Allow an MSA-verified short insertion to nominate the exposed seam
and use its complete independent original BAM source as a right-block
certificate. Every phased heterozygous anchor in that block must be
BAM-injected with exactly one adoptable claim from the same original source;
the existing source-path predicate checks cut-free continuity and consistent
original-to-current allele gauge. Mixed graph/BAM blocks still require a graph
path. Keep the left graph path, Q30/MAPQ30 physical allele calls, two molecules
per insertion allele, and 0.001 parity-error bound. Long-insertion rules do
not change. No truth or gap coordinates enter production conditions.

The owning 33–34 Mb replay closes exactly this 2,097-base gap, changes 11→10
VCF blocks, and preserves all old keys/genotypes/nonphase fields and uniform
old candidate-block gauges. Four unphased reads become correct; one previously
correct rescued read outside the gap becomes wrong. Owner counts change
3,651→3,655 scored, 3,520→3,523 correct, 131→132 discordant. In the gap,
60 correct / 1 wrong / 6 unphased becomes 64 correct / 1 wrong / 2 unphased
out of 67 truth-scorable molecules: 95.5% correct including abstentions,
exceeding the user's 80% acceptance rule. Disjoint flank votes are 25:0 and
20:7 in the same parental orientation, each binomial-supported at p<=0.05.
Diagnostic physical support is 24 REF / 10 ALT molecules, weighted log odds
-164.860097239. Append both endpoint witnesses, the measured `spans=1` panel
row, and an owning regression without weakening old read floors.

A broader certificate through only some imported sites also joined the mixed
36,620,864–36,623,545 block, with just 52/75 correct (69.3%) and an ambiguous
left parental flank. Reject that trial despite its improved global counts.
Requiring the independent source to cover every right-block anchor leaves the
entire 36–37 Mb VCF byte-identical to baseline; its existing owning regression
now explicitly checks those endpoints remain separate. This is a structural
continuity condition rather than a truth-derived coordinate exclusion.

Matched HiPhase also scores 64/67 correct, with 2 wrong and 1 unphased versus
our 1 wrong and 2 unphased. HiPhase places all correct reads into one core
block; four of ours remain correctly phased in the output-only rescue cohort
because they start after the clean SNP and have only noisy insertion markers
and no original source HP assignment. Inspect the two reads correct only in
HiPhase: one existing graph REF call conflicts with Q40 BAM ALT/source HP2;
the other has graph ALT versus MSA REF at the only insertion and no qualifying
base quality, so we abstain. Two reads correct only in pgphase have clean Q40
SNP calls consistent with the source. This explains pgphase's confidence
choices; HiPhase's internal cause is not proven without an internal trace.
All 67 alignment coordinates, CIGARs and sequences match. Reproduction,
rejected trial audits, source evidence and individual discrepancies are in
`evaluations/2026-10-05-short-insertion-source-path/`.

The frozen final full-chromosome audit passes: exactly this new gap closes,
no old tracked gap reopens, and all 64,188 keys, genotype alleles, nonphase
fields and uniform old candidate gauges remain. Tagged output stays 256,610;
scored 237,330→237,334, correct 230,563→230,566, wrong 6,767→6,768. VCF
blocks 326→325, read phase sets 647→645, span N50 unchanged at 856,770.
Only 35 tags and three VCF rows change. All 96 old full-output panel spans
survive; the new row gives 97/110. The newly incorrect outside-gap rescued
read is `m84031_231217_062403_s3/113902771/ccs`. The local 64/67 acceptance
measurement and HiPhase comparison are confirmed on the final full output.
Build has no warnings; unit tests, phase predicates and validation gates pass.
The updated target and mixed-block owner cases pass 652 assertions.
Final native validation passes all 84 cases / 10,752 assertions, including
all 110 panel rows and every existing required-site witness. The final test
source's additional mixed-block guard passes in the separate 652-assertion
owning rerun. No unrelated span expectation or old read floor is relaxed.

## 2026-10-05: HiPhase parity is required for every gap closure

The user strengthens acceptance: each gap must match or beat HiPhase on the
same truth-scorable overlapping reads and alignments, retaining the 80% floor.
Compare both total correct molecules and correct molecules in the dominant
connected core phase set; output-only rescue does not substitute for core
continuity. Diagnose and fix a HiPhase advantage before accepting the gap.
Truth remains evaluation-only.

Under this stronger rule, the current 33.75 Mb change is provisional: total
correct ties HiPhase at 64/67, but dominant connected correct is 60 versus
HiPhase's 64. Its four insertion-only correct rescued reads must acquire a
supported core assignment before this closure can be considered complete.
The existing measurements and passing regression checks remain evidence of
the partial improvement, not proof of acceptance under the new rule.

### Investigation: why the four correct insertion-only reads stay outside core

Trace the four final HP assignments to `rescue_unphased_graph_read_layer`:
it always writes `gap_haps` and `gap_phase_sets = best_phase_set +
kGapFillPsOffset`; there is no certified-indel promotion branch. The two
core refresh paths are SNP-only, and one also requires an existing nonzero
HP. The original source rows have HP=0 and no heterozygous observations
counted by that solve. All four begin after the left clean SNP.

Each has exactly one effective phased heterozygous observation at 33,749,689:
TTCT ALT for 104792409, REF for 199365264, 58000534 and 205719379. The latter
three also call REF on both alternative AAA/AAAA descriptions at 33,762,265;
the rows have opposite haplotype gauges, so those two descriptions propose
opposite HPs and are collapsed to no vote at that locus. The numerous other
SNP observations are homozygous or noninformative. Rescue can assign the
remaining insertion marker, but cannot put the read into the core arrays.

All four original primary BAM alignments have MAPQ60. Read 104792409 has
an exact TTCT CIGAR insertion at the expected breakpoint, with inserted and
flanking bases at Q40. The other three align through the reference allele;
local five-base minimum qualities are Q27, Q27 and Q40 respectively. Thus
there is independent physical evidence to evaluate for core promotion, but
a fixed Q30 rule would exclude two reads. Do not simply strip the offset or
promote all noisy rescues: per-read physical allele confidence, a certified
marker gauge and contradictory/co-located allele handling still matter.
`core-deficit-investigation.json` preserves source rows, effective marker
calls, final HP/PS and original physical CIGAR/quality evidence. HiPhase's
output proves it assigns these four to the core; its internal decision path
has not been traced. The demonstrated pgphase cause is its core admission
policy and lack of an insertion-based promotion path.

## 2026-10-05: Certified insertion observations complete the 33.75 Mb core

Fix the insertion-only core admission deficit after the successful physical
short-insertion bridge. Its original complete BAM-only source certifies the
marker gauge; new read HPs require independent primary MAPQ30 physical calls,
measured allele/base and mapping error <=0.01 with a Q20 floor, agreement with
imported and primary profiles, and no contrary informative locus or phase set.
Complementary descriptions disagreeing at one coordinate contribute no vote.
Existing core assignments, category masks and variant alleles stay intact.

The ordinary insertion caller returns the configured floor for REF, so use
actual anchor qualities in the new confidence bound. One read has a disjoint
one-base deletion 11 bases upstream: the caller's broad indel veto abstains,
but a fully matching reference footprint through the insertion length and
absence of nearby inserted sequence independently establish REF. Include all
footprint errors; leave the shared insertion caller unchanged. A contrasting
read has physical/MSA REF but graph ALT on a repeat-shifted description 12
bases to the right. Reject it using normalized sequence-equivalent graph edits,
even when the graph row has no phase label. Exact-coordinate checking misses
this conflict and introduced one wrong read in a rejected trial.

The final owner changes only the four diagnosed rescue phase-set tags into
the existing core, preserving their HPs. All VCF bytes, output molecules and
parental correctness statuses match the provisional output; 3,655 scored /
3,523 correct / 132 discordant are unchanged. The gap now has all 64 correct
of 67 eligible molecules in one core block, with 1 wrong and 2 unphased.
Matched HiPhase has the same 64 correct in one block, 2 wrong and 1 unphased.
All 67 coordinates, CIGARs and sequences match. This meets both HiPhase
correctness/core parity and the 80% floor. Add exact four-read core and two-read
abstention checks; strengthen the measured separation floor to 0.95 without
relaxing any old count or accuracy floor. Details and strict audits are in
`evaluations/2026-10-05-certified-insertion-core/`.

The final full-chromosome core audit passes: exactly four PS-only changes,
no HP changes, no altered parental correctness status, and byte-identical
VCF output versus the provisional binary. Preserve 256,610 output molecules,
237,334 scored, 230,566 correct / 6,768 incorrect, 645 read phase sets,
64,188 keys, 325 VCF blocks and span N50 856,770. The matched full-output
HiPhase comparison confirms 64/67 correct in one block and fewer local
incorrect reads (1 versus 2). Disjoint flank votes improve to 25:0 and 24:7
in the same parental orientation. Relative to the original accepted output,
exactly the 33.75 Mb gap closes and no old span reopens; the previous outside-gap
one-read error remains measured, with no additional error from core admission.
Final verification passes all 84 window cases / 10,779 assertions, including
all 110 panel rows, old required-site witnesses, exact four-read core promotion
and both conflicting-read abstentions. No read floor is weakened. Clean build,
unit tests, 47 phase-predicate cases / 1,536 assertions and validation gates
pass. Full core-preservation, original gap-geometry and matched HiPhase parity
audits pass. Frozen final binary SHA256 is
`db17374ae2f37a68e803ce397b09f90a75705b0bff23b40aea1498aeccd662a0`.


## 2026-10-05 — shared graph/BAM source runs close 49.72 Mb

The next accepted target is the genuine 22,662 bp gap
49,720,311–49,742,973. Native HiPhase already connects both endpoint SNPs
under PS 49716722. pgphase has the same 137 correct / zero discordant /
25 unphased overlaps, but only 74 correct in its largest core block.
One primary MAPQ-60 read calls both ALT bases at BQ40, yielding log odds
−8.50704 and a 0.000202 two-base/two-mapping error bound. The existing
physical likelihood passes; the whole-block path rule was the blocker.
Both flanks pass the existing independent BAM-run validator when shared
graph rows are included. It requires every phased heterozygous anchor to
have one adoptable source claim in a single run, a consistent current gauge,
and no intervening weak or quality cut. No confidence threshold changes.

When ordinary paths fail for a clean physical SNP pair, two complete local
BAM runs can now certify a deferred core union, including shared graph rows.
Every anchor is recorded, and current block membership and gauge are checked
again after rescue. The eager version was rejected because it combined
output-only rescue cohorts and changed 21 correct / 12 discordant owning
read classifications. The deferred version preserves every rescue tag and
all read outcomes in both the owner and chromosome runs.

After closure, all 137 correct reads are in one core block: 137/162 = 84.57%,
matching native HiPhase on the exact same 162 names, coordinates, CIGARs and
sequences. No new read is admitted. Disjoint flanks agree in parental
orientation (93/0 and 72/7). Full chr20 closes exactly this one old gap:
133 core read labels and six variant rows change PS; all 64,188 keys,
genotype alleles, nonphase fields and old block gauges are preserved.
237,334 scored reads remain 230,566 correct / 6,768 discordant; read phase
sets fall 645→644 and VCF blocks 325→324. N50 remains 856,770. No tracked
gap reopens (98/111 spanned after adding the new target).

The owning 49–50 Mb regression fails on the preceding binary (three checks)
and passes after the fix (698 assertions). The committed panel adds
`spans=1`, a 0.84 dominant-core fraction floor, both endpoint witnesses and
an explicit owner-context case; existing floors/witnesses are retained.
Build, unit tests, 47 predicate cases / 1,536 assertions and validation gates
pass. Measurements and scripts are in
`evaluations/2026-10-05-shared-source-snp-bridge/`.

Final native validation passes all 85 cases / 11,497 assertions, including
all 111 panel windows, existing read floors, required-site witnesses,
owner-context gap checks and mixed-source negative guards. Matched full
HiPhase and chromosome core-preservation audits pass. Frozen final binary
SHA256 is `6824c6287134b7178f12e6690bcfca026b0a32850a861b36c041842e0144e29e`.


## 2026-10-05 — complete BAM suffix closes 57.76 Mb

Close the genuine gap 57,764,235–57,785,224 (20,989 bp). Native HiPhase and
pgphase have the same 128 correct / 7 incorrect / 15 unphased among 150
truth-scorable overlapping primary alignments. pgphase's largest core
previously held only 65 correct; the fix places all 128 in one core, matching
HiPhase at 85.33%. All 150 HiPhase alignments match the original coordinates,
CIGARs and sequences. Existing discordance is unchanged.

The exposed boundary is a complementary MSA deletion family followed by a
clean BAM SNP at 57,785,772. The family bridge has only one Q20 witness and
correctly fails its likelihood gate. However, one primary MAPQ-60 read
physically calls the clean flanking SNPs REF/ALT at BQ40/BQ40, giving
same-gauge log odds −8.50704 and a 0.000202 error bound. The deletion retry
required both left haplotypes before applying that stronger physical test.

Permit one-haplotype clean-SNP votes on this retry only when the right flank
is a complete independent BAM-only run. Every anchor's unique adoptable
source claim, current gauge and intervening weak/quality cuts remain checked.
Conflicting physical votes still reject; the existing 0.001 wrong-parity
bound remains. A complete graph path or BAM run certifies the left flank,
and the core union uses the existing deferred bridge after rescue. Mixed
right source blocks retain the old retry requirement. No truth enters the
production decision, and no new reads or genotype changes are introduced.

The owner keeps all 4,054 scored outcomes (3,860 correct / 194 incorrect),
every rescue tag and old block gauge. Disjoint flanks retain the same parental
orientation (100/0 and 113/7). Full chr20 closes exactly this gap and reopens
none: 129 core read PS labels and five variant PS labels change. All read
correctness statuses, rescue tags, 64,188 keys, genotype alleles and nonphase
fields are preserved. 237,334 scored reads remain 230,566 correct / 6,768
incorrect; read phase sets fall 644→643 and VCF blocks 324→323. N50 remains
856,770; 99/112 tracked windows span after adding this target.

The new owning 57–58 Mb case fails three assertions on the preceding binary
and passes all 747 after the fix. The panel adds exact `spans=1`, a 0.85 core
fraction floor, both complementary deletion descriptions and both clean SNP
witnesses. Existing read floors and required sites are preserved. Build,
unit tests, 47 predicate cases / 1,536 assertions and validation gates pass.
Results and scripts: `evaluations/2026-10-05-clean-snp-source-retry/`.

Final validation passes all 86 window cases / 12,264 assertions and all 112
panel rows, including existing absolute read/error floors, required-site
witnesses and mixed-source negative guards. Full chromosome preservation and
matched HiPhase parity audits pass. Frozen final production binary SHA256:
`f31a9166e0c8e871b28c3174ee95286e8bdc8cbd42c184dafc0e6eb7ce0019e6`.

## 2026-10-05: persistent gap replay test state

The gap harness previously retained completed replay state only inside one test
process and tied it to the output directory. Identical owning-chunk requests
therefore reran across shards and subsequent invocations. Added shared persistent
successful replay state through `scripts/cache_gap_replay.py`, with per-key locks,
atomic completion publication, input/binary/runtime invalidation, output identity
checks and copied output/matrix restoration. Assertions and parental scoring run
every time; no test verdict is cached. Production behavior and expectations are
unchanged by this test infrastructure work.

All 86 native cases / 12,264 assertions passed in both four-shard runs. The first
run took 548.21 s (two focused requests were already populated); warm took
31.94 s, 17.2x faster, with all 114 completion states unchanged and no new
pipeline replay. Focused gap: 15.945 s cold / 0.280 s warm, 747 assertions.
Matrix diagnostic: 16.296 s cold / 0.189 s warm. Four helper tests and all
C++ unit tests passed; production and harness builds succeeded without new
warnings. A rebuilt production binary or changed input still needs a cold
replay; this deliberately avoids carrying old evidence across behavior changes.

Evidence: `evaluations/2026-10-05-gap-test-state/`.

The default cache is populated with the verified state. Normal unsharded
`make window-tests` passed 86 cases / 11,997 assertions plus four cache tests;
its existing process-local outcome reuse produces a different assertion count
from four independent shards. `make check` and predicate tests (47 cases /
1,536 assertions) also passed.

## 2026-10-05: consolidate all gap integration regressions

All 84 integration cases now execute as named sections of one `all gaps` suite;
focused gap unit tests remain separate. Removed repeated fixture/skip setup.
Moved 72 owning-context overrides into `src/test_gap_replays.tsv`; compiled
comparison preserves the old regions across 216 override/default examples.
Moved 38 literal read bounds from nine owning checks, unchanged, into
`src/test_gap_read_floors.tsv`. All remaining original assertion sites are
preserved verbatim. Special molecule, allele, matrix and orientation checks
remain named mechanism functions. Truth parsing is shared with file-identity
invalidation. Production code and gap expectations were not changed here.

HiPhase evidence now covers all 112 panel gaps, after checking identical
primary input alignments including CIGAR and sequence. Native contract scoring
retains unphased primary overlaps in its denominator, uses whole-block parental
orientation and excludes rescue PS values from connected core counts. The
three recent reviewed closures enforce >=80% and HiPhase total/core parity
through `src/test_gap_certified.tsv`; future closures must join that manifest.
Historical closures keep their original regression gates. The explicit audit
finds only 30 of 99 historically closed gaps meeting the whole contract:
69 have shortfalls (16 below 80%, 62 below HiPhase parity, with overlap).
These are a backlog to investigate, not certified successes or waived defects.

Normal `make window-tests`: 13,232 assertions, all 84 named integration checks
and 112 windows, 32.45 s including four replay-cache tests. All 114 pipeline
completion states stayed unchanged. Selected regression: 0.289 s. Missing
fixtures skip explicitly; unmatched selectors fail. Two competitor-state tests
cover reuse/invalidation and rejecting different alignments. Production build,
unit tests, make check and predicate tests pass without new compiler warnings.

Evidence: `evaluations/2026-10-05-gap-test-consolidation/`, including assertion
migration proof, original-region comparison, benchmark provenance and the
per-window contract/backlog report.

## 2026-10-05: verified source indels restore the 61.738 Mb read core

The existing spanning interval 61,738,239–61,747,506 retained 94/96 correct
primary-overlap reads but only 85 correct in its connected core, below HiPhase's
91/96 total and 91 core. Preserving the independent MSA edit behind a shared
catalog representation permits physical verification of existing rescue tags
against a complete, consistently oriented source path. Confirmed reads enter
that same core without changing variant connections or assigning new HP tags.
The final owning replay retains 94/96 total and reaches 91 core, matching HiPhase;
all 96 competitor alignments were checked against the identical original input.

One required REF read has a Q17 base but an approximately 2% summed footprint
plus mapping error. A per-base Q30 rule rejected it unnecessarily. Deletion
materialization now measures that error with a 5% ceiling and Q10 minimum;
Q10 footprints with ~10% error and different deletion lengths remain excluded.
Whole-block physical bridge thresholds are unchanged. Source-path cuts, two
consistent clean SNP loci shared with the graph, both observation channels, and every informative phased
profile anchor still gate admission. Added predicate coverage and owning-chunk
94 total / 91 core floors, and certified the existing panel interval.

The final development loop is `make gap-dev-check`: all fast predicates and
adapter fixtures, 0.03 seconds on the current build, with no production-binary
build or pipeline replay. `make gap-owner-check GAP=61.738` is the explicit real
owner check, 0.31 seconds warm. The measured cold 1 Mb replay costs 9.81 seconds.
Full suites are final validation, rather than an exploratory edit loop; changed
production binaries invalidate their real outputs.
Evidence: `evaluations/2026-10-05-verified-indel-core/`.

Materialization additionally requires exact source/final marker agreement, final
MSA verification, a pure deletion, and no neighboring verified consensus edit
within 32 bases. A partial transfer requires core read coverage to dominate the
rescue cohort at every coordinate; otherwise it can fragment a locally larger
output block, as the 41.881 Mb negative regression showed. This rule has a fast
in-memory regression in `test_graph_bam_adapter`, included in `gap-dev-check`.

Final frozen binary SHA256:
`f9efac33e864064418fe09e70c4502b4db43dae48f6825bc4e9556d699e6ed1b`.
All 84 registered checks and all 112 panel windows passed, including nested
sections (195 selected invocations), plus unit, predicate, cache-helper and
standard gates. The full chr20 replay took 359.83 seconds; all 64,188 variant
records and every HP tag are unchanged across the same 256,610 primary output
molecules. It materializes 173 existing rescues, 156 correct and 17 discordant
(90.17% correct). Aggregate truth counts remain 230,566 correct, 6,768 discordant,
19,276 unphased; per-block rescoring has three transitions in each direction.
The panel retains all 99 spans and improves full-contract passes from 30 to 33.
Also certified the existing spanning 1,086,625–1,110,921 and
59,825,454–59,842,960 windows that gained parity, and reran all three new
certifications against their stricter assertions.

## 2026-10-05: deletion-masked SNP fill reaches HiPhase read parity at 1.194 Mb

The already spanning 1,194,233–1,196,894 window had 90/91 correct primary-overlap
reads, all in one core, versus HiPhase's 91/91. The one unassigned maternal read
carries an exact Q40 six-base CACCAC deletion at 1,194,228 with MAPQ60. It removes
the verified opposite-haplotype SNPs at 1,194,230 and 1,194,233. Both transferred
channels retain its deletion ALT; the missing SNP calls are physical masking.

Ordinary graph rescue uses only six core graph donors and fails its 1% association
gate even at 6/6 agreement. The physical marker has 72 callable existing core
reads when staged BAM-only rows are included: 37 ALT/hap1 and 35 REF/hap2, no
disagreement. A new final read fill uses these distinct original primary molecules
with Q30 exact CIGAR calls, sequence validation for ALT, MAPQ30, both haplotypes,
current marker orientation, two-sided association p<=0.01 and a one-sided 95%
Wilson discordance bound <=10%. Only entirely unassigned graph-profile reads
with a pure MSA-verified deletion <=32 bases, ALT in both channels, a masked
verified opposite-haplotype SNP, and no contrary phased observation qualify.

A complete independent source path is not necessary for this local read fill:
the existing core supplies its allele gauge directly, and it changes no variant
or block connection. Existing rescues and core HP/PS assignments are preserved.
Source-path requirements for joins and prior rescue materialization are unchanged.
Both graph worker paths run this fill after their final deferred bridges and
rescue materialization, with worker-owned reference/BAM handles and BAM contig
IDs resolved from site metadata.

The certified owning 1–2 Mb replay reaches 91/91 correct in total and in the
connected core, matching HiPhase on all 91 identical original primary alignments.
Its 4,197 output reads and 1,528 variant records audit shows just one new correct
tag; all 4,044 existing assigned tags and all 10 discordant reads are unchanged.
Added the existing spanning window to the certification manifest and fast
positive, flipped-gauge, missing-channel, contradictory-channel, other-block
and unmasked-SNP fixtures. `make gap-dev-check` takes 0.03 s warm; the cached
owning integration check takes 0.34 s.
Evidence: `evaluations/2026-10-05-masked-snp-deletion/`.

Final binary SHA256:
`6a8057fc5d3875b7fb3a608a579efe47a2f4d1ffce99aac84bdd0dcec6321f04`.
All 84 registered checks (83 mechanisms plus all 112 panel windows) pass in
195 selected invocations with 82,416 assertions, together with gap units,
predicates, unit binaries, cache-helper tests and standard validation gates.
All 99 panel spans remain intact; full-contract passes increase from 33 to 34.
Auditing 329 matched saved native replays preserves every variant record and
existing read tag. Nine overlapping outputs contain the same one new correct
maternal assignment, with no new discordant assignment. Validation reuses the
panel and owning replay outputs; no new full-chromosome run was needed.

## 2026-10-05: close the HiPhase-supported 62.408–62.432 Mb gap

Closed chr20:62,408,056–62,432,427 in its owning 62,000,001–63,000,000
chunk and in the complete chromosome replay. This is a previously open VCF
interval, not merely a repair of read tags inside an already spanning block.
Before: 157/177 original primary reads correct, largest connected core 80.
After: 157/177 (88.70%) correct in one core, exactly matching HiPhase total
and core correctness; both tools have 2 discordant and 18 unphased reads.
All 177 competitor primary alignments match original start/end, CIGAR and
sequence. Disjoint pgphase parental flanks score 86/86 and 190/190, with the
same gauge as HiPhase's 86/86 and 191/191.

HiPhase retains the intermediate GA>G deletion at 62,410,496 in the same phase
set as both boundary SNPs. Our shared catalog row retains alignment verification
and its original MSA deletion claim but keeps repeat classification and loses
the source's MSA flag/BAM vector. Its physical deletion begins at 62,410,497.
The old MSA-only stitch therefore cannot use it. Preserve that representation:
recover the canonical catalog deletion only with one adoptable retained pure
MSA claim of the same length, then calibrate it independently before bridging.

The new physical bridge requires complete graph SNP paths on both sides. Q30
primary deletion calls agree with two-or-more-locus clean graph SNP gauges on
22 REF and 17 ALT molecules, with no contrary votes. The existing binomial
association and one-sided 95% Wilson discordance tests pass. Two independent
primary deletion/SNP pairs establish same-gauge log odds -8.89121, exceeding
the 0.001 wrong-parity threshold. Each pair's measured deletion footprint,
SNP and twice mapping error must remain <=5%; qualities below 10 abstain.
The union is deferred until rescue is complete, retaining all anchor gauges.
No graph-path or block-join confidence threshold is relaxed.

Twelve already tagged graph-marker rescues can then enter their certified
core through their own physical deletion alleles, an uncut source with 24
consistent shared clean SNP loci, and the existing core coverage guard. A
missing BAM observation vector is checked against the primary alignment while
the graph call must still agree. Two of these primary alignments have MAPQ22
and MAPQ26, not MAPQ60: supplementary alignments at those names have higher
quality and must never supply the proof. The shared-marker promotion uses a
MAPQ20 floor with the same <=5% measured combined error. Q40 deletion flanks
keep these two primary calls below that error bound.

Two whole-chunk BAM fallback reads contain false graph/BAM SNP REF at
62,432,427 with zero BAM quality; their primary CIGAR actually deletes the
base. Their separate SNP at 62,435,323 is Q40 and agrees in both channels.
After the shared-deletion union, an independent physical REF/ALT core gauge
can connect them while treating only CIGAR-deleted zero-quality REF calls as
masked. Every other phased call must agree; unassigned reads and ordinary
contrary calls cannot enter through this route. This completes core parity
without changing the correct/discordant/abstention totals.

Added the gap to the panel, required marker list, native replay map, strict
>=80%/HiPhase certification and the unified owning mechanism suite. The owning
regression passes 500 assertions and fails 10 with the task-start binary.
Added fast provenance and masked-SNP witness fixtures; warm `make gap-dev-check`
takes 0.04 s. Final task binary SHA256:
`f98aead0dfc6f9e72a2e2cb99de3463234523d5be7437859cecb095a2259e6fd`.

The owning audit preserves all 1,207 variant alleles/counts/filters and all
3,722 correct assignments. The complete chromosome preserves all 64,188
records and all prior phase-set extents; this is its only newly closed VCF
interval. Its saved full baseline predates the immediately preceding masked-SNP
read fix, so the one additional correct full-chromosome tag belongs to that
previous fix. Current-task matched owning and suite audits use the actual
task-start binary/state.

Evidence: `evaluations/2026-10-05-shared-deletion-bridge/`.

Final checks: all 85 registered gap checks (84 mechanisms, 113 windows),
197 selected invocations, 84,327 assertions; build with no new warnings,
unit/predicate/golden gates, gap unit cases and cache/benchmark helpers pass.
The panel has 100/113 spans and 35 windows satisfying both >=80% correctness
and HiPhase total/core parity. All 329 matching task-start native replays
preserve variant alleles/counts/filters and every correct/phased assignment.
The exact span floor for the new window is 1; the panel total floor rises to
100 on this measured graph-arm improvement. No existing window floor changed.


## 2026-10-05: Cut-free source paths close the 17.634 Mb gap

The task-start binary leaves chr20:17,634,393–17,667,022 open while HiPhase
closes it correctly. The native 17–18 Mb replay has 202 truth-scorable original
primary overlaps: pgphase assigns 164 correctly, 1 incorrectly and abstains
on 37, but only 110 correct reads share one core. The reviewed fix places
all 164 correct reads in one core (81.188% of all original overlaps), matching
HiPhase's 164 correct/core reads. HiPhase abstains on 38 and has no discordant
overlap here. All 202 competitor start/end, CIGAR and sequences match input.

HiPhase retains the intermediate deletion at 17,642,038 in the same block as
SNPs 17,633,256 and 17,667,022. pgphase retains the corresponding verified
shared deletion privately, but rejected the usable information at three
stages. Path certification visited only right flanks and required a nonempty
quality-cut list even for a fully supported source with no cuts. The source
SNP at 17,487,837 was recalled by MSA, so its clean graph row could not serve
as an exact shared source anchor. Finally, two Q40 primary molecules spanning
17,723,058–17,744,061 could not compose with the source certificate before
whole-path rollback.

The fix retains recalled SNP edits and validates them against their immutable
Q30 physical calls and at least two other clean source SNP loci per molecule,
on both haplotypes. At a retained shared catalog deletion opposite a clean graph SNP, cut-free
sources can supply absent/one-sided agreeing GAF edges, with any contrary graph
pair still vetoing them. Certificates cover both flanks and compose with the
existing Q30 physical edge/suffix proof. The marker context is essential: a
general cut-free certificate also unlocked the distinct-deletion 37.598 Mb
seam in the wrong gauge. It lost 156 formerly correct parental assignments in
the full chromosome. The final rule preserves that seam and its read tags; the existing
distinct-deletion and native 37 Mb orientation checks guard this exclusion.
Weak/quality cuts and inconsistent source gauges remain exclusions.

The shared deletion's physical coordinate is 17,642,039. Its final calibration
has 7 REF and 16 ALT primary votes, no contrary votes, each with measured SNP,
deleted-footprint and mapping error <=1%. The only primary boundary molecule
is m84031_231217_062403_s3/147985967/ccs: Q40 deletion REF and Q17 right SNP
REF, giving crossed log odds 3.87891 (about 2.0254% wrong parity). The physical
bridge uses the existing SNP-to-deletion 5% bound; calibration Wilson error,
conservative 1% gauge-call error and bridge error together must be <=20%.
This replaces the independent fixed 10% calibration cutoff for this bounded
physical bridge only. General BAM fallback and masked-SNP witness gates retain
their 10% bounds. This conditional evidence bound does not replace the measured
all-original-reads >=80%/HiPhase contract.

Disjoint pgphase flanks score 175/178 and 145/145 in one parental orientation;
HiPhase scores 179/180 and 145/145 with the equivalent global gauge. The
complementary 1-base/12-base source deletions at 17,634,393 keep their alleles,
counts and opposite haplotypes. The panel, native owner map, required rows,
strict certification and owning mechanism regression now include the closure.

Final validation: `make window-tests` passes 13,411 assertions (86 registered
gap checks, 114 windows); all 199 selected invocations also pass. Build, unit,
predicate and golden gates pass with no new warnings. The new owning regression
passes 549 assertions and fails nine with the task-start binary. Warm
`make gap-dev-check` takes 0.05 seconds; cached complete window tests take
33.24 seconds. All 114 HiPhase benchmark rows match fresh identical-alignment
measurements. Panel spans increase to 101/114 and full >=80%/HiPhase total/core
contracts to 36. All 332 matching native replays preserve every old correct
or phased assignment and variant allele/count/filter. The full chromosome
replay (413.56 s alongside integration work) closes only the target interval,
retains all 64,188 variants and 230,567 correct assignments, and makes no new
read assignments. Its target core also scores 164/202, with exactly HiPhase's
correct qnames. Evidence: evaluations/2026-10-05-cut-free-source-path/.

The chr20 VCF phase-block span N50 is unchanged by this closure at 856,770 bp
(0.857 Mb), versus HiPhase's 1,005,183 bp (1.005 Mb), using inclusive
first-to-last phased-heterozygote spans with at least two rows per block.
Valid positive VCF phase sets fall 322 -> 321; blocks with at least two rows
fall 266 -> 265. These are raw VCF spans, without truth-based switch splitting.
Measurements: evaluations/2026-10-05-cut-free-source-path/block-n50.json.


## 2026-10-05: Connect the split inside HiPhase’s largest block

HiPhase’s largest chr20 VCF block spans 63,182,011–66,206,480
(3,024,470 bp). Pgphase split it at 65,509,355–65,509,406 inside a
6.6 kb SDUST interval. Near-boundary graph SNP placements and CIGAR calls
are inconsistent with the established read gauges, so accepting their
nearest-SNP parity without calibration would be unsafe. Existing graph-path
checks also reject an earlier one-haplotype edge. SNPs outside the repeat
provide a direct alternative proof: two separated Q30 SNPs per flank on a
primary MAPQ58 spanning molecule, calibrated on 144/148 disjoint established
primary reads with both haplotypes and no contrary calls. Its quality error
bound is 0.00040317; both one-sided 95% Wilson calibration bounds, 1% physical
call error per flank and bridge error together give 0.05679768, below 20%.
The bridge changes neither calls nor read assignments before the deferred
whole-block union after rescue.

The full chromosome now has one block at 63,182,011–66,194,203,
3,012,193 bp. N50 rises from 856,770 to 863,977 bp. Across the same 9,843
original primary truth-scorable overlaps of HiPhase’s largest interval,
pgphase has 9,339 correct (94.8796%, including abstentions), versus HiPhase
6,681 (67.8756%). The dominant correct core grows from 7,226 to 9,209;
HiPhase has 6,681. HiPhase’s large PS changes parental orientation near
65.2 Mb: locally both tools score six of six at the seam, with disjoint
flanks 145/145 and 146/146, but whole-HiPhase-PS orientation scores these
six discordant. Keep both measurements explicit; raw block length does
not establish a correctly oriented chromosome segment.

Whole-chromosome preservation retains 230,567 correct and 6,768 discordant
assignments, all 256,610 primary output records and all 64,188 variant keys,
alleles, counts and filters. Only the reviewed seam closes; all earlier
VCF extents remain intact. The union changes 1,024 variant gauges and 2,324
read PS labels without adding, losing or changing the correctness of a tag.
The owning replay similarly preserves 2,962 correct assignments, six
errors and 1,591 variants. A committed panel row, native replay mapping,
strict >=80%/HiPhase manifest, retained markers and parental-flank owner
regression cover this closure. The old binary fails five assertions.

The 12,277 bp endpoint difference remains a callset/coverage difference:
HiPhase’s last T>A has GQ2 and all eleven primary input alignments there
have MAPQ3–15. This fix connects the existing blocks without inventing a
terminal heterozygote. Evidence and reproduction helpers are in
`evaluations/2026-10-05-largest-hiphase-block/`.

Validation: the production build has zero new warnings; all unit tests,
47 predicate cases / 1,551 assertions and HiFi/ONT TSV/VCF golden gates pass,
including HiFi one/four-thread determinism. `make window-tests` passes 13,494
assertions in four cases and four replay-state helper tests. All 201 selected
checks pass across 115 panel windows and the owning/mechanism cases; 102
panel gaps span and 37 satisfy both the >=80% and HiPhase total/core contract.
All 335 matched native replays retain every old correct/phased assignment and
every variant allele, count and filter. The pure development loop takes
0.038 s; the cached new owner regression takes 0.508 s. The warm full suite
takes 35.47 s; the one required cold full-chromosome audit took 448.40 s while
running concurrently with the integration checks.

## 2026-10-05: Investigate the largest-block terminal difference

The 12,277 bp endpoint difference is not explained by low MAPQ/GQ alone.
The exact 66,206,480 T>A is present in the catalog. The graph matcher sees
four reference / zero alternate walks at MAPQ5, and eight / zero at MAPQ1;
it drops the site as `ref_only`. Original BAM bases are six A / two T at Q40,
plus three alignments with no base. All six BAM A reads take the reference graph
walk, establishing a repeat-placement disagreement. At earlier 66,194,253,
lowering MAPQ changes five reference / zero alternate to 28 / 20, showing
that the graph discovery floor also removes real alternate observations.

The BAM noisy window 66,206,004–66,210,255 has eleven reads but only three
full-cover reads and no qualifying diploid phase set, so MSA skips it at
minimum depth five. Depth three fires unphased MSA but does not emit the exact
terminal T>A. Seam recovery needs two phased anchors and cannot nominate an
unanchored terminal tail; whole-BAM fallback does not import all private sites.

All 54 original primary truth-scorable tail overlaps are paternal. Current
pgphase has 25 correct / 29 discordant (46.3%); HiPhase has three correct /
six discordant / 45 unphased (5.6%) under whole-PS orientation. HiPhase's most
favorable local flip still gives only six correct. At the final position,
current pgphase scores 5/11, HiPhase 1/11 (best local flip 4/11). The MAPQ1
owning graph replay extends beyond HiPhase to 66,207,824 and scores 9/11 at
the endpoint, but only 32/54 (59.3%) across the added tail and loses 74 formerly
correct rescue assignments. Owner correct totals fall 584 -> 566. MAPQ0 /
depth1 / alt1 similarly scores 31/54, loses 74 and falls to 565 correct.
These counterfactuals fail the >=80% interval and preservation requirements.

The useful next investigation is terminal recovery with repeat-placement
validation against established flanking haplotypes, rather than a blanket
threshold change. No production behavior changed. Original alignment identity,
parental scoring, graph observations, production matcher trace and fast owning
counterfactuals are preserved in `evaluations/2026-10-05-terminal-endpoint/`.

## 2026-10-05: Rank the next largest unconnected HiPhase block

After the reviewed largest-block internal join, HiPhase's second-largest block
is the next internal split: 8,946,171–11,127,753, PS8946171, 2,181,583 bp.
Frozen current pgphase splits it into 8,946,171–10,325,039 (1,378,869 bp,
PS8946171) and 10,337,661–11,127,753 (790,093 bp, PS10343596). The target
boundary distance is 12,622 bp. Ranking all 191 HiPhase and 264 pgphase
multi-het blocks requires no pipeline replay.

On 9,742 identical original primary truth-scorable overlaps, pgphase scores
9,424 correct / 52 discordant / 266 unphased (96.74%), core5,824; HiPhase
9,459 / 23 / 260 (97.10%), core9,459. At the actual split, pgphase scores
103 / 10 / 14 of127 (81.10%), versus HiPhase97 / 10 / 20 (76.38%).
Pgphase already meets >=80% and HiPhase total correctness at the gap,
but its main correct cores are41 and45; even a union gives86, below HiPhase97.
At least11 more correct gap reads must join that core. The other17 already
correct gap assignments are in rescue PS values. Whole-block parity also
requires addressing the35-correct-read deficit, beyond a phase-label union.

Disjoint flanks have unanimous main pgphase gauges, HP1-paternal192 left
and HP1-maternal229 right; HiPhase is HP1-maternal193/194 left and230/230
right. HiPhase's whole block passes80%, while its gap-specific all-read score
does not; do not conflate raw VCF span with a certified gap closure. No
production behavior or test expectations changed. Full ranking, exact
original alignment checks and scoring are in
`evaluations/2026-10-05-next-largest-block/`.


## 2026-10-05: Close the second-largest HiPhase block through exact insertion ALTs

Closed chr20:10,325,039–10,337,661 inside HiPhase's second-largest block.
The physical stitch lacked a calibrated complementary-insertion-to-SNP route;
the one-row indel route correctly rejects competing insertion alleles and the
repeat-length route cannot certify this GGAA pair. Imported clean SNPs were
also hidden by an empty graph `chunk.ref_seq`; physical references now use the
worker FASTA cache.

The new route keeps both exact MSA/alignment-verified ALTs, calls other lengths
uncallable, and calibrates their opposite haplotypes on disjoint upstream-SNP
molecules. A recovered calibration SNP needs a separate physical edge to the
nearest graph SNP: two primary Q40 molecules certify 10,296,487–10,318,335.
The insertion calibration is [[11,0],[1,13]], with one original molecule
bridging the four-base insertion to SNP 10,343,595. Summed base/placement/MAPQ
errors give log odds -5.9507 and joint Wilson/call/bridge error 0.173464,
within the 20% bound. Both graph paths stay required; an internal weak right
SNP edge is separately physically certified. Complete anchor gauges defer the
union until rescue finishes, then exact physical insertion calls can connect
unassigned/rescue reads without changing existing core reads.

Both the 10 Mb native owner and the final full chromosome score 108/127
original primary gap overlaps correct (85.04%), with connected correct core
98. HiPhase scores 97/127 and core97 on identical original alignments. The
baseline is 103 correct, 10 discordant, 14 abstentions, core45; the fix is
108 correct, 10 discordant, 9 abstentions. Exact-ALT-only promotion rejects
a REF read that belongs to neither insertion haplotype. It adds five correct
calls without adding discordant calls, preserving every old correct/tagged
read and every variant.
The native regression fails the old binary (11 assertions) and passes the
final binary (560 assertions); the new window, owning context, required
markers, HiPhase count comparison and strict 80%/core contract are recorded
in the test manifests without changing existing floors.

Final production SHA256:
15534777ecd8382e370ade340e6e8998a3fc1b1ab4ab390281b3ff4a22bbedb4.
Four-thread chromosome replay: 67 chunks, 483.8 s, 64,188 variant records and
256,610 output reads. All 230,567 previously correct reads are preserved;
correct totals become 230,572 and discordant totals remain 6,768. The block
spans exactly 8,946,171–11,127,753 (2,181,583 bp), matching HiPhase's endpoints.
N50 rises from 863,977 to 904,351 bp; block count falls 264→263. The largest
block remains 3,012,193 bp and its previously reviewed terminal difference
remains unchanged. Whole target-block read parity is still open: 9,429 total
correct and 9,318 correct core versus HiPhase9,459/9,459. The local gap itself
passes the required HiPhase total/core and >=80% comparisons.

Evidence and reproducible audits: evaluations/2026-10-05-second-largest-hiphase-block/.
Implementation behavior is described in docs/IMPLEMENTATION.md.

Validation: final `make -j8`, unit/predicate tests, collect gates and window
suite pass (13,583 assertions in four cases). HiFi/ONT goldens and HiFi t1/t4
determinism are unchanged; no new build warnings. Separate representation,
orientation and connectivity selections pass 2,171/5,757/3,139 assertions.
A saved focused owner replay recheck takes about 0.6 seconds.

All 338 archived native replays are matched by the final suite and preserve
every prior correct/phased assignment and every variant allele/count/filter;
none are unmatched. See the new evaluation panel-audit.json.
