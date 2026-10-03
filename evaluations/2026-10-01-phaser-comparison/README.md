# Fresh pgphase / HiPhase chr20 comparison

Date: 2026-10-01 (America/Los_Angeles).

## Setup

- HG002 HiFi, complete CHM13 chr20, **272,016 input reads**.
- pgphase commit `94410d04ba8acda80b56feb9a3ec562249c2328d`, default
  graph + BAM recovery, existing graph catalog and GAF.
- Unmodified HiPhase `1.6.0-ac3f399`, default realignment,
  `--ignore-read-groups`. Native arm uses the DeepVariant calls made from the
  same annotated BAM; diagnostic arm uses this fresh pgphase run's exact
  VCF keys and unordered genotypes with GT phase separators and PS removed.
- Eight worker threads per tool; HiPhase default four I/O threads.
- Runs execute sequentially on an Intel Core i9-7900X (10 cores / 20 logical
  CPUs). One fresh invocation per arm, existing filesystem cache retained.
  These are single-run measurements, not cold-cache or repeated estimates.
- Normalized SAM SHA256 agrees between the original and linear-reheadered
  BAM, including sequences, qualities, flags, CIGARs and tags; only contig
  prefixes are normalized for comparison. The reference sequences also agree.
  All input records are primary and have unique query names.
- All timings include each executable's output work. Alignment, catalog
  construction, DeepVariant calling, truth preparation, diagnostic VCF
  preparation, and evaluation are outside the measured tool invocations.
  pgphase's timed command includes its internal BAM calling/recovery; HiPhase
  starts from an already called VCF. No DeepVariant runtime is inferred here.

Outputs are `/tmp/pgphase-hiphase-comparison-2026-10-01/`.
`identity.json` records binary hashes, input provenance, hardware and run order.
The four resource files retain the actual commands and GNU time measurements.
`run.sh` reproduces all runs into a new output directory and refuses reuse.

## Whole-input read accuracy and continuity

| Metric | pgphase graph + recovery | HiPhase / DeepVariant | HiPhase / pgphase calls |
| --- | ---: | ---: | ---: |
| Input reads | 272,016 | 272,016 | 272,016 |
| Phased reads | **237,098** | 233,353 | 230,138 |
| Phased percentage | **87.1633%** | 85.7865% | 84.6046% |
| Unphased reads | **34,918** | 38,663 | 41,878 |
| Truth-correct reads | **229,923** | 223,800 | 216,061 |
| Truth errors | **7,175** | 9,553 | 14,077 |
| Accuracy among phased reads | **96.9738%** | 95.9062% | 93.8832% |
| Correct reads / all input reads | **84.5255%** | 82.2746% | 79.4295% |
| Read phase sets | 680 | 196 | 200 |
| Phased VCF heterozygotes | 62,443 | 77,123 | 60,035 |
| VCF phase blocks | 337 | 196 | 200 |
| VCF block span N50 | 739,888 bp | **1,005,183 bp** | 919,153 bp |
| Largest VCF block | 1,723,925 bp | **3,024,470 bp** | 2,788,940 bp |

The native comparison adds **3,745 phased reads**, **6,123 correct assignments**
and removes **2,378 errors** in pgphase. Conditional accuracy improves by
1.0676 percentage points. HiPhase retains the block-continuity advantage.

The diagnostic arm preserves all **62,850** input VCF keys and unordered
genotypes, including the homozygous records. Its 200 blocks leave 2,408 input
heterozygotes unphased. Phased-site counts in the native DeepVariant arm are
not a variant-accuracy ranking because its input callset differs.

## Runtime and memory

| Measurement | pgphase | HiPhase / DeepVariant | HiPhase / pgphase calls |
| --- | ---: | ---: | ---: |
| Native command wall time | 147.49 s | 87.87 s | 71.38 s |
| User CPU time | 923.03 s | 719.06 s | 652.63 s |
| System CPU time | 87.90 s | 41.38 s | 42.47 s |
| Peak resident memory | 20.125 GiB | 0.389 GiB | 0.224 GiB |
| Additional full-BAM write + index | 94.38 s | Included | Included |
| Total with full indexed BAM | **241.87 s** | **87.87 s** | **71.38 s** |

The graph command emits an **unaligned tag-only BAM**, containing read names
and HP/PS, not original alignment sequences. HiPhase writes and indexes a full
alignment BAM. `materialize_tags.py` therefore copies the pgphase assignments
onto every original BAM record and indexes the result, timed separately with
four read/write I/O threads. All assignments are verified unchanged. This
Python postprocessing is a measured way to deliver comparable outputs, not
an optimized native implementation or isolated solver timing.

Compared with native HiPhase / DeepVariant, pgphase's command is **1.68 times
slower**, or **2.75 times** including the measured full-BAM step. Its peak
resident memory is **51.69 times greater**. Accuracy and read coverage lead;
runtime, memory and block N50 do not.

HiPhase uses its default automatic local fallback after broad global
realignment failure: five blocks in the DeepVariant run and eight in the
same-callset run. Neither diagnostic logging nor modified competitor code is
used. All timed commands exit successfully.

## Scoring definitions and overlap

`score.py` joins output assignments to the complete original BAM population by
query name. It requires HP 1/2 and positive PS. The pgphase tag-only BAM has
256,586 rows; its 15,430 omitted input reads remain in the denominator and
count as unphased. Materializing the full BAM does not add phase assignments.
The common primary-only population is identical to the whole input population.

Truth is the existing `test_data/derived/chr20_truth_hap.tsv`. Within each
emitted `(chromosome, PS)`, choose the majority parental orientation and count
the minority assignments as errors. Every phased read has an available truth
label. Accuracy divides correct assignments by phased/scored reads; coverage
divides phased reads by all 272,016 input reads. Truth is used only for
evaluation. This metric does not estimate switches or flips, and the competing
tools' differing block sizes remain visible in the continuity metrics.

N50 uses the span from first to last phased heterozygote in each VCF PS,
inclusive, weighted by these block spans; it is not NG50 or a truth-corrected
continuity measure. Read PS count includes independent read-only BAM blocks.

For pgphase versus native HiPhase:

- Both phase **227,642** reads. Independently orienting each tool's blocks on
  this common subset yields **98.3421%** for pgphase and **96.5112%** for HiPhase.
- Using the full-population block orientations instead gives 98.3228% and
  96.4932%; the accuracy advantage persists under this calibration.
- pgphase alone phases **9,456** reads, **6,099** correct under its full-block
  orientation. HiPhase alone phases **5,711**, **4,141** correct.
- The shared-read accuracy gain accounts for another **4,165** correct
  assignments; the total correct-yield gain is 6,123.

The fresh pgphase totals independently reproduce the certified pre-comparison
237,098 / 229,923 / 7,175 counts. Native HiPhase reproduces the established
233,353 / 223,800 / 9,553 counts. `metrics.json` includes both cohorts, overlap,
genotype preservation and BAM/reference identity checks; `results.tsv` is the
compact table. No production pipeline behavior or regression floor changes in
this evaluation.
