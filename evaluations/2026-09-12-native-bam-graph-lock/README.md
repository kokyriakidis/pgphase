# Native-BAM Graph-Lock Evaluation

Date: 2026-09-12

Sample: HG002 HiFi

Reference: CHM13v2.0

Chromosomes: chr12, chr18, and chr20

Evaluated base revision: `c1f103f` plus the graph-lock and extraction changes
committed with this record. Environment: Python 3.12.3, pysam 0.24.0, samtools
1.19.2, and htslib 1.21.

## Question

Can clean heterozygous sites called from the surjected BAM fill gaps in graph
phasing without allowing BAM evidence to damage trusted graph assignments?

## Inputs

All input paths are relative to the repository except truth BAMs:

| chromosome | reference | graph catalog | GAF | surjected BAM | truth BAM |
|---|---|---|---|---|---|
| chr18 | `test_data/chm13v2.0.chr18.renamed.fa` | `test_data/chr18.sites.vcf.gz` | `test_data/HG002.chr18.annotated.coord.gaf.gz` | `test_data/HG002_chr18_hifi_mapped_to_CHM13_chr18_annotated.bam` | `~/Downloads/pgphase-eval-data/truth/chr18/diplinator_merged.bam` |
| chr12 | `test_data/chm13v2.0.chr12.renamed.fa` | `test_data/chr12.sites.vcf.gz` | `test_data/HG002.chr12.annotated.coord.gaf.gz` | `test_data/HG002_chr12_hifi_mapped_to_CHM13_chr12_annotated.bam` | `~/Downloads/pgphase-eval-data/truth/chr12/diplinator_merged.bam` |

The exact native-BAM pipeline is executable as `commands.sh`; the two failed
chr18 controls are in `controls_chr18.sh`. The original DeepVariant-private
chr20 pipeline and graph-lock controls are in `commands_chr20.sh`. Outputs are written to
`/tmp/pgphase_native_bam_graph_lock` by default; set `OUT_ROOT` to retain large
intermediates elsewhere.

## Configurations

- `graph`: graph-only baseline with read margin 2 and anchor AF margin 0.12.
- `private_gq10`: chr18 control using every native GQ10 private gap candidate.
- `bridge_gq10`: chr18 control retaining a coordinate-spanning read path.
- `clean_snp_joint`: direct joint phasing with graph-absent CLEAN SNPs, GQ10,
  VAF 0.30-0.70, and a two-read MAPQ20 bridge path.
- `graph_lock`: graph assignments copied unchanged; only graph-unphased hybrid
  reads from blocks with 10 shared reads, margin 5, purity 0.90, and both
  haplotypes represented are added.
- `graph_lock_final50`: the same lock followed by a 50-read gate on the final
  graph-labelled output phase sets.

## Results

| chromosome | configuration | private sites | evaluated | discordant | Hamming | N50 bp | bad PS |
|---|---|---:|---:|---:|---:|---:|---:|
| chr18 | graph | 0 | 224,746 | 388 | 0.17% | 1,746,048 | 6 |
| chr18 | private GQ10 | 12,007 | 230,819 | 3,527 | 1.53% | 1,730,887 | 17 |
| chr18 | bridge GQ10 | 565 | 226,453 | 3,141 | 1.39% | 1,755,742 | 4 |
| chr18 | clean SNP joint | 111 | 226,016 | 3,148 | 1.39% | 1,766,679 | 4 |
| chr18 | **graph lock** | **111** | **228,020** | **405** | **0.18%** | **1,748,648** | **6** |
| chr12 | graph | 0 | 410,699 | 257 | 0.06% | 422,422 | 2 |
| chr12 | clean SNP joint | 131 | 413,505 | 405 | 0.10% | 542,044 | 5 |
| chr12 | **graph lock** | **131** | **415,431** | **286** | **0.07%** | **425,298** | **2** |
| chr20 | graph | 0 | 175,843 | 248 | 0.14% | 937,061 | 8 |
| chr20 | joint GQ10 PS0 | 958 | 180,804 | 366 | 0.20% | 951,797 | 10 |
| chr20 | **joint GQ10 PS50** | **958** | **179,494** | **269** | **0.15%** | **991,274** | **0** |
| chr20 | graph lock PS0 | 958 | 179,612 | 329 | 0.18% | 943,722 | 8 |
| chr20 | graph lock, PS50 proposal | 958 | 179,563 | 326 | 0.18% | 943,722 | 8 |
| chr20 | graph lock, purity 0.95 | 958 | 179,323 | 327 | 0.18% | 943,722 | 8 |
| chr20 | graph lock + final PS50 | 958 | 178,697 | 266 | 0.15% | 960,028 | 1 |

Graph lock added 8,006 truth-evaluable reads across chr12 and chr18 for 46
additional errors: 99.43% accuracy among net additions. It preserved all graph
HP/PS assignments by construction.

## Updates And Decisions

1. Native pgphase BAM calls are not interchangeable with the polished
   DeepVariant callset used for the successful chr20 experiment. GQ10 alone
   admitted 12,007 chr18 sites.
2. Coordinate-spanning read support is not phase support. It reduced the site
   set but did not prevent two large chr18 mixed blocks.
3. PS50 cannot repair this failure because the dominant incorrect phase sets
   contain more than 3,000 reads.
4. Private extraction now supports CLEAN-SNP, graph-position-absence, and VAF
   gates. These improve proposal quality but do not make direct joint output
   safe on chr18.
5. Graph lock is the selected native-BAM policy. Direct joint output remains an
   experimental proposal source.
6. Graph lock extends reads assigned to existing graph PS labels. It does not
   yet merge independent graph blocks; such merging needs allele-consistent
   support on both flanks.
7. Chr20 confirms graph-lock safety but not universal dominance. With the
   polished DeepVariant private set, direct joint GQ10+PS50 remains the selected
   balance: it has 797 more reads and 31 kb higher N50 than locked final-PS50,
   for three additional errors. Locked final-PS50 is the higher-precision point.
8. Raising graph-lock orientation purity from 0.90 to 0.95 removes 289 reads
   but only two errors. The residual errors are primarily inherited small graph
   phase sets, not weakly oriented additions.

## Stored Artifacts

Each configuration directory contains the evaluator's `summary.txt`,
`summary.json`, `worst_phase_sets.tsv`, and `bad_phase_regions.bed`. The compact
cross-configuration ledger is `results.tsv`. Large alignments and per-read truth
tables are intentionally excluded from git.

The stored chr20 private input has SHA-256
`e67e119d371a54e7edcd8a99d75ce5d2073228e4189aad43e37e8c9504e0a3ff`.
