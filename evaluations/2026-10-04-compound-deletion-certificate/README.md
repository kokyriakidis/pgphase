# Compound CIGAR deletion certificate

## Target and signal

Target: chr20:7,901,413–7,918,883 (17,470 bp), outside the excluded
centromere. Saved same-BAM HiPhase results nominate a correctly oriented join.
Truth and competitor outputs are evaluation inputs only.

The BAM MSA retains two separate deletion rows at VCF anchor 7,918,883:
`ATATTTTTTTTTTT>A` and `ATATTTTTTTTTTTTT>A`. Their genotypes are complementary.
Two primary MAPQ-60 reads encode the longer allele as a 15-base CIGAR deletion
at 7,918,882 with compensating mismatches. Deleting the CIGAR reference segment
alone gives a different sequence from deleting the candidate segment. Comparing
the complete query between the surviving reference anchors instead gives the
exact candidate allele `ATAT` on both reads. Their minimum retained base
qualities are 40 and 27. See `target_signal.json`.

The maternal read with split deletions and the paternal read with a 16-base
deletion do not certify that allele. No different deletion length supplies a
REF call, and no candidate rows are merged or realigned.

## Retained design

- `bam_matches_deletion_sequence` is an ALT-only certificate. Its false result
  means unverified, never REF. The BAM source solver and existing physical
  REF/ALT callers retain their original CIGAR edit certificate.
- Only the graph-SNP to complementary BAM-deletion bridge consumes this
  certificate. Reads must match exactly one separate ALT row, pass MAPQ30 and
  the Q20 sequence floor, and agree unanimously across at least two independent
  molecules. The existing quality-weighted 0.001 parity rule remains.
- A transfer-created right-deletion seam tries this bridge before rejecting
  its missing or one-haplotype clean-SNP pair. Existing SNP votes veto an
  inconsistent proposed orientation. Clean SNPs remain the preferred bridge.
- One conditional graph SNP in the right block may be bypassed only through
  the existing significant direct flank edge, with both haplotypes represented
  by at least two agreeing graph reads. A dominant reversal remains a veto.

Production uses neither truth nor competitor output and has no fixture
coordinate gate. Both complete block paths are checked before joining.

## Local measurements

| Output | Gap-overlapping phased | Correct | Errors | Correct in one block |
|---|---:|---:|---:|---:|
| Accepted pgphase baseline | 116 | 114 | 2 | 64 |
| Correctly joined pgphase | 116 | 114 | 2 | 104 |
| Saved HiPhase / DeepVariant | 104 | 103 | 1 | 103 |
| Saved HiPhase / pgphase calls | 86 | 84 | 2 | 84 |

There are 143 truth-scorable input overlaps. The fix changes continuity,
with no new read errors. pgphase retains more correct reads than either saved
competitor arm; its local purity remains below HiPhase / DeepVariant. This is
not a fresh competitor run or a claim of higher local purity.

In owning 7–8 Mb, retain all 690 variant keys and 4,091 scored reads:
3,926 correct and 165 discordant before and after. Read PS count falls 18→17,
VCF blocks 9→8. The owning regression checks separate exact deletion rows,
SNP/deletion parental orientation, local read counts and errors, and owning
read counts and errors. It fails on the starting binary (three orientation
assertions) and passes after the fix. The panel adds the coordinate case with
`spans=1`, plus explicit boundary-site retrieval requirements.

## Rejected experiments

Broad shared-source path admission made a wrong whole-block join at 37.55 Mb
and converted 155 old correct assignments to errors. Trying all complementary
deletions, pure BAM component certificates, and one-haplotype SNP retry changes
alone did not close a new chromosome gap. None is retained.

Applying the query certificate during BAM source recall changed source
classification, introduced eight variant keys, and made a wrong 20.8 Mb join:
230,211→230,156 correct and 7,097→7,195 errors; 90 old correct reads became
incorrect. Source isolation restored all read truth assignments, but a broad
physical REF/ALT API still reopened the protected 64.14 Mb connection. A failed
query certificate can otherwise fall through to a REF call. Reject both trials
and expose the new check as ALT-only to the paired bridge. Preserve their
reports in `full-compound-parity.json`, `full-compound-blocks.json`, and the
`rejected-stitch-only-*` reports.

## Reproduction

Build with `make -j8`, then run `make unit-tests predicate-tests`, `make check`,
and `make window-tests`. The complete fixture inputs and derived truth map are
required. The new owning test is:

```
./test_gap_windows 'shifted compound CIGAR deletions close the owning 7.9 Mb gap'
```

For repeated scoring, `replay_cached_panel.py` checks the final executable
SHA256, input size/mtime, every exact CLI argument except its output directory,
and every cached output hash. Unknown requests run natively. The cache contains
fresh runs of the final binary, not baseline output.

## Final chr20 verification

The only new VCF union is the nominated gap. All earlier connections and
all 64,188 variant keys remain. VCF blocks fall 332→331; scored read phase
sets fall 654→653. All 237,308 scored reads retain their truth assignments:
230,211 correct / 7,097 discordant, 97.009372% before and after. All 230,211
old correct assignments remain correct. No phased SNP or existing SNP gauge
is lost. The expanded panel spans 91/104 coordinates, versus 90/104 before.
The remaining 13 open coordinates include three intentional split controls.

The joined block spans 7,859,472–8,081,718 (222,247 bp). VCF N50 remains
806,449 bp; NG50 remains 684,798 bp. See `full-parity.json`,
`block-transfer.json`, and `metrics.json`. The native eight-thread chr20 run
takes 370.92 seconds during concurrent regression generation; this is not a
controlled runtime comparison.

Final binary SHA256:
`00f8c411a595359473ac2be0cce0f9639a0810be9373310b37f024a6da2e017d`.

Validation: standalone units, 1,497 predicate assertions / 47 cases,
7,977 window assertions / 74 cases covering 104 coordinates, HiFi/ONT
TSV/VCF goldens and HiFi one/four-thread determinism pass. The cache contains
106 freshly generated native requests. Existing expectations are retained;
only the new coordinate and the graph span TOTAL 90→91 are added.
