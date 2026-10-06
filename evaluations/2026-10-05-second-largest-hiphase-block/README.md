# Second-largest HiPhase block: complementary insertion bridge

Target: chr20:10,325,039–10,337,661, the internal split in HiPhase's
8,946,171–11,127,753 block (2,181,583 bp). Baseline binary SHA256:
`477aaecb8adadf8f6f20bc380676de81f15d51f8eacbdb787c2cf5a8bf18140e`.
Final binary SHA256:
`15534777ecd8382e370ade340e6e8998a3fc1b1ab4ab390281b3ff4a22bbedb4`.

## Cause and correction

HiPhase represents the left locus as a multiallelic heterozygote:
`T>TGGAA,TGGAAGGAA`, GT `2|1`, GQ52, DP65, AD6,34,21. Pgphase keeps
separate complementary insertion rows. Its physical stitching had no
calibrated complementary-insertion-to-SNP route; the single-indel route
rejects competing alleles, and the repeat-length route accepts only simpler
motifs. No original molecule calls both bounding clean SNPs.

A separate bug hid recovered physical SNPs behind `chunk.ref_seq`, which
graph chunks do not populate. Their reference bases now come from the
worker's FASTA cache.

The new route uses exact ALT observations, with other lengths uncallable.
It certifies recovered SNP 10,318,335 against the preceding graph SNP
10,296,487 on two independent primary Q40 molecules. Twenty-five disjoint
calibration molecules then give the complementary insertion gauge
`[[11,0],[1,13]]`. One primary molecule physically connects the four-base
insertion to right SNP 10,343,595. Its summed inserted-base, placement,
SNP-base and mapping error gives log odds -5.9507. The one-sided Wilson
calibration bound is 0.160866; adding the 1% calibration-call allowance and
bridge error gives a joint error bound of 0.173464, below 20%.

Both graph paths remain required. The right block's weak edge at
10,843,933–10,859,317 is separately certified from primary BAM SNP calls;
the remaining SNP suffix must pass its existing path check. Whole source
membership does not stand in for these physical edges. Anchor gauges are
retained until read rescue finishes, then verified before the deferred union.
Exact, independently checked insertion observations can promote reads into
that union, while existing core assignments are preserved.

## Native owner certification

Replay: `CHM13#0#chr20:10000001-11000000`, four threads. All 127 original
primary truth-scorable overlapping reads remain in the denominator.

| Measure | Baseline pgphase | Fixed pgphase | HiPhase |
|---|---:|---:|---:|
| Correct reads | 103 | 108 | 97 |
| Discordant reads | 10 | 10 | 10 |
| Unphased reads | 14 | 9 | 20 |
| Correct / all original reads | 81.10% | 85.04% | 76.38% |
| Dominant connected correct core | 45 | 98 | 97 |
| Spans the gap | no | yes | yes |

HiPhase's whole target block is high confidence overall (97.10% original-read
correctness); its local gap accuracy is below 80%. Pgphase must independently
pass 80% and both HiPhase count comparisons, which it does.

The owner audit preserves every existing correct or phased read and all
1,075 VCF records, including alleles, counts, filters and non-gauge fields.
It adds five correct assignments without adding a discordant assignment.
Exact-ALT-only promotion rejects an otherwise incorrectly assigned REF read,
which belongs to neither allele of this complementary locus. Three
previously correct rescue reads enter the connected core. Left and right
parental orientations are checked on disjoint original-primary flanks.
See `owner-contract.tsv`, `owner-preservation.json`, and `owner-results.json`.

The native-owner regression fails on the baseline binary (11 failed
assertions) and passes on the final binary (560 assertions). The gap is added
to the window panel, owning replay manifest, required-marker manifest, HiPhase
counts and strict 80%/core certification manifest. Existing floors are kept.

## Reproduction

Full chromosome outputs use the normal four-thread graph command, without a
matrix dump, in `test_data/tmp_gap_fix71/frozen_final/0/`. Run `audit_block.py`
with the bench-phasers Python environment to compare those outputs with the
frozen previous chromosome and indexed HiPhase BAM. It checks original BAM
geometry, CIGAR and sequence identity, parental orientations, local and
whole-block accuracy, preservation, block spans and N50.

Run the focused regression with:

```sh
LD_LIBRARY_PATH=/usr/lib/x86_64-linux-gnu make gap-owner-check GAP="second largest"
```

## Full chromosome result

The final production binary processes 67 owning chunks with four workers in
483.8 seconds. Pgphase now has one PS8946171 block spanning exactly
8,946,171–11,127,753 (2,181,583 bp), matching HiPhase's endpoints. It retains
all 2,220 original phased heterozygous rows in the two joined pgphase blocks.

| Measure | Previous pgphase | Fixed pgphase | HiPhase |
|---|---:|---:|---:|
| Chr20 block count | 264 | 263 | 191 |
| Chr20 N50 bp | 863,977 | 904,351 | 1,005,183 |
| Largest block bp | 3,012,193 | 3,012,193 | 3,024,470 |
| Correct original reads over the target block | 9,424 | 9,429 | 9,459 |
| Target block connected correct core | 5,824 | 9,318 | 9,459 |

The reviewed gap itself passes both HiPhase comparisons. Read accuracy across
the whole 2.18 Mb block still trails HiPhase by 30 correct reads; this closure
does not certify parity for every read elsewhere in that block. There are
261 abstentions in pgphase and 260 in HiPhase; pgphase has 52 discordant
reads versus 23 in HiPhase across that full interval.

The full chromosome audit checks 9,742 identical original HiPhase alignments
against the original BAM, including start/end, CIGAR and sequence hashes. It
preserves all 230,567 previously correct pgphase assignments, all previously
tagged reads, and all 64,188 variant records and non-gauge fields. Correct
assignments become 230,572; discordant assignments remain 6,768.
Both disjoint flank cohorts retain their correctness and consistent parental
orientation. The full chromosome and native owner have the same 108/127
correct total and 98 connected correct core in the gap. See `results.json`.

## Validation

`make -j8`, `make unit-tests predicate-tests check`, and `make window-tests`
pass with the final binary; the complete window suite passes 13,583 assertions
in four cases. HiFi/ONT goldens and HiFi one/four-thread determinism remain
unchanged. No new build warnings are introduced. Separate representation,
orientation and connectivity selections pass 2,171, 5,757 and 3,139 assertions.
The saved owner replay allows a focused repeat in about 0.6 seconds.

All 338 archived native replays have matching final outputs and preserve every
previous correct/phased read and every variant allele, count and filter.
No prior replay is unmatched. See `panel-audit.json`; `audit_panel.py` repeats
the comparison against a completed window-suite output directory.
