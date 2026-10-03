# Preserve the certified BAM prefix of a complementary deletion bridge

## Bug and correction

The owning 64–65 Mb recovery solve preserves a BAM source beginning at
64,068,280. Its sole weak cut is at 64,140,314, downstream of the clean SNP
at 64,128,828 and complementary deletion pair at 64,134,226. The existing
rightward deletion-pair transfer starts at the pair, leaving earlier connected
BAM anchors behind. It creates an artificial phase boundary inside the audited
source run.

The correction extends that transfer to the first BAM anchor after the complete
physical footprint of the preceding catalog alleles. It requires the nearest
clean BAM SNP to agree unanimously with mutually exclusive deletion ALT calls,
both SNP allele classes, distinct known-MAPQ-30 BAM reads, and a one-sided
unlinked binomial tail no greater than 0.001. Every transferred anchor must
have adoptable provenance in the same original source and orientation; the
interval must contain neither a weak nor a quality cut. Double-REF/double-ALT
calls abstain. The existing independent pair-to-right certificate is still
required. No candidate representation changes, new realignment, truth input,
competitor input, coordinate exceptions, or threshold relaxation enter production.

The local prefix transfer leaves catalog anchors independent. Subsequent
ordinary stitching closes the larger blocks with its own evidence checks.
Independent output-only read rescue gauges retain their separate handling.

`source-trace.log` records pre-fix source provenance/cuts from an exploratory
instrumented build; the instrumentation is absent from the final code. The
preceding BAM-pair-MAPQ evaluation incorrectly described the 64.128 Mb SNP
and the nearby 64.121 Mb SNP as graph-owned when examining missing GAF edges.
They are BAM-owned. This investigation corrects that diagnosis: missing GAF
observations at those BAM rows do not invalidate their cut-free BAM source.
The implementation still validates graph-block connections separately.

## Evidence and results

At **64,128,828–64,134,226**, 27 high-MAPQ primary BAM reads physically span the
boundaries. Eighteen exclusive-ALT triples are usable: 12 SNP-REF/first-deletion
ALT and six SNP-ALT/second-deletion ALT. All agree; their fair-parity tail is
2^-18. Seven double-REF triples and two missing deletion pairs abstain.

| Measurement | Before | After |
|---|---:|---:|
| Local phased / correct / discordant reads | 45 / 44 / 1 | 45 / 44 / 1 |
| Local dominant block's correct reads / truth-scorable overlaps | 22 / 54 | 39 / 54 |
| chr20 phased/scored reads | 237,127 | 237,127 |
| chr20 correct / discordant | 229,974 / 7,153 | 229,974 / 7,153 |
| chr20 read concordance | 96.983473% | 96.983473% |
| Read phase sets | 673 | 670 |
| VCF phase blocks | 337 | 336 |
| VCF keys | 63,492 | 63,492 |
| VCF span N50 | 756,878 bp | 756,878 bp |

Every scored read retains its individual correct/discordant classification,
not merely the aggregate totals. There are 3,993 changed read HP/PS pairs and
2,460 changed VCF rows from legitimate block gauge/label changes, with zero
lost or gained keys. The extracted helper reproduces the initial prototype
outputs exactly in both the owning chunk and the whole chromosome.

Native HiPhase/DV phases 49 local reads, with 48 correct and one discordant;
saved-pgphase-callset HiPhase phases all 49 correctly. The latter still uses the
Oct 1 callset. Thus this fixes continuity, **not local read-yield or accuracy
parity**. Three reads correctly phased by native HiPhase remain unphased by
pgphase: one has four callable BAM matrix observations, the other two have
none. All three have MAPQ 60 and ambiguous/missing calls at the deletion pair.
`remaining-local-reads.json` preserves that next investigation. The fourth
additional native-HiPhase assignment is locally discordant.

The panel grows from 95 to **96 coordinate cases**, with **82 spanning cases**
(previously 81). The new target uses the owning 64–65 Mb replay. Exact assertions
require the SNP ALT to match the two-base-deletion ALT and oppose the one-base
ALT. Parental flank agreement, read separation, the existing owning ceiling of
three errors, earlier deletion joins, and both complementary rows are protected.
Two old separation assertions are deliberately replaced with exact joined
orientations because the independently validated catalog join is the intended
improvement. No previous read-quality floor is lowered.

A fresh extent/parental-flank audit excludes 25–30 Mb. It finds **315 extent
gaps and 41 distinct competitor-supported nominations (51 records)**, down from
316/42/53. Thirty-four are supported by native HiPhase, seventeen by the saved
pgphase callset, with overlap. These nominations still require representation
and exact allele-path analysis; they are not all solved.

Final binary SHA256:
`98037cb463527b02e6acb31393d1d2d71cca92323d6922942d9367d868291e3d`.
The fresh eight-thread chr20 run takes 292.54 s while native replays share the
machine; this is not an isolated runtime comparison. All 105 fresh native
replay requests use this SHA256 and complete successfully in 390.30 s.

## Reproduce

```bash
make -j8
make unit-tests predicate-tests
PGPHASE_TEST_WORKDIR=/tmp/pgphase-prefix-tests make window-tests
```

The baseline binary fails the new owning-chunk exact-connection regression.
The unit tests exercise the actual prefix helper, including source/cut/gauge
mismatches, conflicting or ambiguous calls, unknown/low MAPQ, duplicate names,
intervening catalog/third-allele rows, and a missing nearest SNP call. The
committed logs record the final window and predicate gates.

All evaluation scripts expose `--help`. Example (with the fixture Python
interpreter that provides pysam):

```bash
python measure_output_parity.py \
  --before /tmp/pgphase-gap-next22/final-full \
  --after /tmp/pgphase-gap-next23/final-full --output results.json
python measure_prefix_transfer.py \
  --before /tmp/pgphase-gap-next22/final-full \
  --after /tmp/pgphase-gap-next23/final-full \
  --matrix /tmp/pgphase-gap-next23/final-64/matrix.chunk0.recovery-final.tsv \
  --output prefix-evidence.json
python audit_competitor_gaps.py \
  --pg-vcf /tmp/pgphase-gap-next23/final-full/phased.vcf \
  --pg-bam /tmp/pgphase-gap-next23/final-full/phased.bam \
  --competitor-root /tmp/pgphase-hiphase-comparison-2026-10-01 \
  --output current-competitor-bridges.json
```
