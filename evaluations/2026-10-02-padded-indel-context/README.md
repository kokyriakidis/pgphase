# Padded indel identity and BAM genotype retention

## Correctness fixes

`vcf_to_variant_key` consumed the complete common prefix of unequal-length
alleles but kept their remaining common suffix. For example, `ACG→ATCG`
became a replacement of `CG` with `TCG`, rather than insertion of `T`.
This prevented an exact match to the physical BAM allele. The conversion now
trims the remaining suffix after consuming the prefix. Consuming the prefix
first preserves the established BAM/MSA placement inside repeats. Complex
replacements keep their nonmatching reference and ALT spans.

Correcting identity exposed a transfer defect: matching an independently phased
BAM MSA call to an unphased graph repeat retained the graph demotion and lost
the BAM genotype. The normalization-only chromosome run lost four VCF calls,
although every read assignment remained identical. The affected selected
catalog alleles are:

| Raw catalog POS | Selected REF→ALT | Retained BAM VCF allele |
|---|---|---|
| 10488935 | GCTTTTTTTT→GTT | 10488935:GCTTTTTT→G |
| 10891998 | CCAC→CCC | 10891999:CA→C |
| 35490915 | TATATTTTTTTTTTTTTTT→TATTTTTTTTTTTTTTTT | 35490917:TA→T |
| 41885027 | ATATATATATTTTTTT→ATATATATTTTTTTT | 41885034:TA→T |

Transfer now retains the verified BAM physical key, genotype, counts and MSA
proof together on the same binary allele row, in an unused source-scoped PS.
The parallel catalog metadata still identifies that selected allele. Its
observations come from the source solve; unphased graph repeat calls cannot
silently become observations of an independently recalled MSA consensus.
Candidates in validated complete source flanks use the same admission rule as
private BAM rows. Conflicting source claims remain independent.

The change applies to suffix-padded, unphased binary graph repeats with an exact
MSA allele match. It does not overwrite established graph phases or interpret
whole multiallelic genotypes as binary calls. Production uses no truth,
competitor output or fixture coordinates. It introduces no realignment,
allele merging or relaxed stitch thresholds.

## Experiments rejected

- Projecting other catalog walks through a unique single-base SNP branch
  increased SNP calls, but changed recovery seams and lost private BAM evidence.
  Both initial and supplemental projection trials were removed.
- Lowering physical SNP bridge base quality from 30 to 20 changed no output
  across 28 owning chunks. The threshold remains 30.
- Two deletion-only physical ALT/ALT pairs at 4.866153–4.874129 Mb suggested
  a flip with an apparently strong base/mapping likelihood. Both source paths
  passed their certificates, but the proposed join reversed the parental
  connection. Owning-chunk errors rose 16→448; 428 previously correct reads
  became wrong. Local errors rose 0→13. Native-DV HiPhase has 50/51 correct
  there. The trial was removed. Base quality of deletion flanks does not measure
  the deletion representation's allele error probability.
- Adopting all matching unphased repeats changed previously correct read
  assignments at 33 and 48 Mb. Those trials were removed. The retained fix
  addresses the padded-context identity defect and preserves existing handling
  for other catalog rows.

## Regression checks

- Physical-key tests cover pure insertions/deletions, complex replacements,
  minimal alleles, SNPs and unchanged repeat placement: 872 assertions / 39 cases.
- Adapter units check genotype/key/proof transfer and reject already phased,
  unverified, homozygous, multiallelic and nonrepeat ownership changes.
- A new owning-chunk test requires all four exact BAM alleles to remain phased,
  with their verified MSA category, no duplicate key and unchanged read accuracy
  floors. The normalization-only binary fails three assertions, one per chunk.
- A new 4.866 Mb regression permits a future join only with the correct opposite
  boundary allele relation and preserved owning/local read accuracy. The rejected
  deletion connector fails six assertions; the accepted baseline passes.
- Existing span expectations and the 98-coordinate / 84-required-span panel
  are retained. No extra gap closure is claimed by this change.

All standalone units and 4929 window assertions / 65 cases pass. All 105 fresh
native panel requests completed in 372.52 s against the frozen final binary;
one additional matrix regression ran natively. Native outputs are reused only
for identical commands from this same binary, with SHA verification before and
after each request. `validation.json` records the counts and additional command.

Owning-chunk comparisons and rejected-trial logs are recorded alongside this
report. `measure_output_parity.py` compares every read's tags and parental
classification, all emitted VCF keys, block counts and span N50. Truth is used
only for evaluation. `compare_local_gap.py` limits scoring to read names that
physically overlap the nominated gap in the input BAM.

## Final chromosome check

The corrected binary is byte-for-byte equivalent to the preceding accepted
output at the read-tag and VCF-record level. All four verified MSA calls remain
present. There are no new joins or individual read-classification changes.

| Metric | Before | After |
|---|---:|---:|
| Output read names | 256601 | 256601 |
| Phased/scored reads | 237134 | 237134 |
| Truth-correct reads | 229981 | 229981 |
| Discordant reads | 7153 | 7153 |
| Read concordance | 96.983562% | 96.983562% |
| Read phase sets | 667 | 667 |
| VCF records / keys | 63492 | 63492 |
| VCF phase blocks | 334 | 334 |
| Span N50 | 756878 bp | 756878 bp |

Full chromosome wall time was 276.61 s with eight threads. Other native panel
runs were executing concurrently, so this is a correctness run, not a runtime
comparison. Input paths and commands match the previous chr20 graph+recovery
evaluation. The frozen binary SHA256 is
`00c7e302c47a3c1fe58d6a4fa5ef647fc1d8030171ce362482a8cb3db794af91`.
Exact results: `chr20-comparison.json`. The intermediate loss is retained
separately in `suffix-only-chr20-comparison.json`.
