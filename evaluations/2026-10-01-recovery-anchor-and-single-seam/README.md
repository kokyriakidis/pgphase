# Recovery anchors, shifted insertions and single-seam retries

Date: 2026-10-01. Baseline: `94410d04ba8acda80b56feb9a3ec562249c2328d`.
Inputs and competitor runs are the identical-BAM chr20 comparison in
`../2026-10-01-phaser-comparison/`. Truth is used only for evaluation and tests.
Production contains no truth, competitor, coordinate or read-name exception.

## Correctness changes

1. **One recovery anchor predicate.** A positive PS and unequal internal
   haplotype labels were insufficient: classified homozygotes could retain
   stale unequal labels, and an unknown allele could look like REF. Recovery
   now requires two known alleles, exactly one selected ALT, and a
   non-homozygous category. Multiallelic graph `1|2` remains eligible; a row
   projecting to `0|0` does not. Source selection, observation-path cuts and
   transfer certificates share this predicate. Unknown transferred gauges
   continue to veto certification.
2. **Shifted insertion backfill.** A missing MSA call could receive exact-
   coordinate REF even when the BAM carries an equivalent insertion elsewhere
   in the repeat. Backfill can certify that unchanged edit at known MAPQ30/Q30
   and require an independently callable clean SNP in the same source gauge.
   Existing MSA calls, separate complementary rows, candidate keys and gauges
   are retained. The SNP test is shared with the existing deletion repair.
3. **Focused retry for a single seam.** Complementary co-located BAM boundary
   rows were eligible for focused retries only when another padded seam
   touched their group. A single seam now gets the same full adjacent graph
   phase-set context and unchanged source-path, row-preservation and transfer
   validation. No stitch threshold is relaxed.
4. **Fallback admission scope.** A failed focused trial with complementary
   boundaries could incorrectly authorize broad unplaced-read MSA. That broad
   solve lacks the focused path certificate and split the protected 4.767 Mb
   short replay despite passing full-chromosome span checks. It must now
   independently satisfy the original unique-boundary admission rule. This
   retains the old 4.767 Mb span and accuracy floors.
5. **De novo MSA coverage index.** Per-read alignment bounds used the cluster
   number to index original read coverage flags. They now use the clustered
   read's input index after filtering. A synthetic partial-first input
   reproduces the defect. This error also exists in the checked local
   longcallD `align.c` helper. The ordinary wrapper sorts full-cover reads
   first, so this is an API correctness fix, not a demonstrated cause of the
   remaining chr20 gaps.

## Full chr20 result

| Metric | Baseline | Fixed |
| --- | ---: | ---: |
| Input primary reads | 272,016 | 272,016 |
| Phased / truth-scored reads | 237,098 | 237,118 |
| Truth-correct reads | 229,923 | 229,956 |
| Truth errors | 7,175 | 7,162 |
| Accuracy among phased reads | 96.973825% | 96.979563% |
| Read phase sets | 680 | 675 |
| VCF rows | 62,850 | 63,401 |
| VCF phase blocks | 337 | 338 |
| VCF block span N50 | 739,888 bp | 756,878 bp |
| Tracked spans | 77/93 | 78/93 |

All previously closed tracked gaps remain closed; only
**10,727,690–10,746,628** changes from open to closed. Eleven historical
competitor targets and four controls remain open. Native HiPhase's accuracy
is 95.9062%, with 233,353 phased reads and 1,005,183 bp N50 in the preserved
comparison. Its continuity advantage remains. These recovery runs overlap
other validation work and are not a new runtime comparison.

The compact tag BAM grows from 256,586 to 256,601 rows. Missing input reads
remain unphased in the common 272,016-read denominator. Truth accuracy chooses
one majority orientation per read PS; it is not a switch-error metric.

### 10.727 Mb: correct connection, remaining read ambiguity

The original BAM source has a left SNP and a right block with complementary
4- and 16-base deletion rows plus a four-base insertion. The boundary has
sparse conflicting paired observations. A grouped retry requires unique
boundaries; the focused retry keeps both deletion rows and uses the complete
adjacent graph blocks to establish the source path.

The accepted source spans both graph flanks without a weak cut. The final
left SNP, 16-base deletion and insertion ALT share one haplotype; the 4-base
deletion ALT is opposite. Both graph SNP flanks retain their relative parity,
and disjoint parental flanks agree. This matches HiPhase's allele connection.

| Gap-overlap assignments | Baseline pgphase | Fixed pgphase | Native HiPhase |
| --- | ---: | ---: | ---: |
| Phased reads | 143 | 143 | 141 |
| Correct | 134 | 136 | 136 |
| Errors | 9 | 7 | 5 |
| Accuracy | 93.71% | 95.10% | 96.45% |

The owning 10–11 Mb chunk gains ten phased reads, eight correct and two
discordant. It retains every previous VCF key and adds 197 private-context
rows. The new regression pins all six flank/bridge rows, complementary
deletion parity, no parental flank switch and the local 136-correct / at-most-
seven-error result. The existing panel case now requires a positive span and
both intervening deletion rows; its existing 98% whole-chunk floor stays.

**Local accuracy has not matched HiPhase.** Three repeat-only reads receive
the wrong MSA haplotype in pgphase but the correct one in HiPhase. They reach
neither clean outer SNP, and carry compound repeat indels. The existing MSA
already covers the relevant variant footprint; truncation is not the cause.
Four other wrong reads are also wrong in native HiPhase. Do not treat raw
indel vote agreement as proof of a correct per-read assignment.

### 59 Mb: graph child SNPs represent insertions

The additional focused solve reduces owning-chunk errors from **14 to 2**,
phases 15 additional reads, and raises correct assignments from 3,654 to 3,681.
Two old graph SNP keys disappear, while 356
verified private rows enter the VCF. This is a supported representation
correction, not candidate loss in transfer:

| Position | Old graph row | Physical MAPQ30 base calls | New MSA / DeepVariant row |
| --- | --- | --- | --- |
| 59,679,069 | `A→T` | 57 A, 2 T | `A→AT` |
| 59,757,500 | `T→C` | 61 T, 1 C | `T→TC` |

Most insertion calls lie immediately after the reference base (25 inserted
T and 21 inserted C calls respectively). The child graph rows retain their
observations internally. The existing conflicting-allele output filter keeps
the insertion representation. A new owning-chunk regression requires both
insertions, opposite ALT haplotypes in one PS, no parental flank switch and
at most two errors across the chunk.

The shared-anchor fix separately changes the 34 Mb chunk: five correctly
phased reads become unphased and three wrong assignments become correct.
The whole-chromosome gain includes that small coverage cost.

## Rejected experiments

- Certifying shifted insertions without independent SNP source-gauge evidence
  reopens the protected 41.9008 Mb gap and removes 20 old keys. Edit equivalence
  alone cannot establish the source haplotype. This trial is not retained.
- Refreshing every assigned MSA observation during focused retries also loses
  protected 41 Mb rows. Existing assigned MSA observations remain authoritative.
- At 37.458 Mb, fresh HiPhase does not phase both historical boundary keys:
  native DeepVariant calls the left deletion homozygous REF and leaves the
  multiallelic right deletion unphased. Same-callset HiPhase skips the left
  deletion and the two-base right deletion. Coordinate read-block coverage
  does not establish that it solved pgphase's complete allele chain. The
  focused pgphase trial has conflicting cuts and remains rejected.

## Reproduction

Build `make -j8`. With the fixture inputs and truth map prepared:

```bash
./pgphase collect-graph-variation \
  --ref test_data/chm13v2.0.chr20.renamed.fa \
  --bam test_data/HG002_chr20_hifi_mapped_to_CHM13_chr20_annotated.bam \
  --sites test_data/chr20.sites.striped.vcf.gz \
  --gaf test_data/HG002.chr20.annotated.coord.gaf.gz \
  -t 8 -r 'CHM13#0#chr20' -o /tmp/recovery-candidates.tsv \
  --phased-vcf-out /tmp/recovery-phased.vcf \
  --phased-bam-out /tmp/recovery-phased.bam
make unit-tests predicate-tests parity-tests upstream-parity-tests
make window-tests
```

Full native outputs and fresh window replays are under
`/tmp/pgphase-gap-next17/verified{,-native}/`. `metrics.json` and
`identities.json` record the before/after counts and exact tracked span delta.
Validation logs and the final binary identity are retained alongside this
report. Native window outputs were computed with the final binary; scoring
reuses only matching argument lists and the same SHA256.

Final validation: 4,573 assertions / 57 window cases; 762 / 34 predicate
cases; 27 / 9 port-parity cases; 168,696 / 7 original-C phase-parity cases;
standalone units all pass. The panel retains 93 coordinate cases; two dedicated
owning-chunk regressions were added, giving 57 Catch2 cases rather than 57
coordinate gaps. There are 105 fresh native panel commands plus the separate
57 Mb matrix regression. The final full run takes 184.07 seconds while panel
validation is also active; this timing must not replace the sequential runtime
comparison. The build has the existing vendored abPOA unused SIMD helper
warning and no new warnings.

Final binary SHA256:
`598eb8f99b1e14d1ce25eb9a8ef0305bea60b9b44089c2a3f41f9a05c95ac676`.
