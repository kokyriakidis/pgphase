# Further gap search: short validation windows and long insertion alleles

No new gap closes and no production change is retained. The accepted baseline
is `test_data/tmp_gap_fix55/full-final`, binary SHA256
`bdce74cda1d0ec407fe7a0cfaf9e30d2c77a60b605d6f193ee92e2cb16503933`.
All pre-existing source changes are preserved.

Two isolated experiments were evaluated across the complete chromosome:

1. Remove the 128-base allele-length cap from
   `physical_equivalent_insertion_call`. Retain all sequence-equivalence,
   base-quality, nearby-indel interference, physical bridge, and source-path
   checks. The 29 and 31 Mb owning replays and full chromosome are unchanged.
2. Admit positive gaps shorter than 10 kb to the existing full-flank physical
   validator (`kMinSingletonGap = 1`). Retain its singleton/deletion geometry,
   shared-site gauges, source cuts and physical-call checks. The full
   chromosome has no newly spanning old gap and exactly unchanged read tags
   and VCF rows. The 37.459--37.462 Mb control remains unconnected.

Both full comparisons retain 256,610 output reads, 237,330 truth-scored reads,
230,235 correct / 7,095 discordant assignments, 650 read PS labels, 64,188 VCF
keys, 328 VCF blocks and span N50 856,770 bp. The saved parity JSONs record
zero tag changes, zero VCF-row changes, and zero lost or gained keys. The
snapshotted trial binary hashes and exact VCF hashes are in `manifest.json`.
No threshold or expected span was retained from either experiment.

Additional diagnoses:

- The 29,303,608--29,303,677 boundary gives 16 REF/ALT and three ALT/REF pairs
  at MAPQ >=20 and base quality >=20, with no opposite pair. Its downstream
  SNP is an MSA/noisy BAM candidate. This is a nomination, not a certificate
  of the full graph flank. An independent 29,114,041--29,355,992 replay gives
  261 matching and 234 flipped old left-block genotypes, with six left anchors
  missing. All 67 anchors of the old right BAM component are absent. It cannot
  justify a uniform union of the two original blocks. The bounded replay's
  aggregate read concordance is 69.47%, also unsuitable as a trusted replay.
- There are no secondary or supplementary records in the annotated input
  BAM (272,016 primary alignments); split-alignment linkage is not an available
  alternative for these gaps.
- Replays padded to 100 kb on either side of the 34, 38 and 56 Mb chunk
  boundaries still do not span their nominated gaps.
- The 24.122--24.132 Mb insertion/deletion boundary has several insertion
  lengths beyond its two candidate ALTs, and contradictory deletion events.
  The 37.990--38.005 Mb private deletion pair likewise contains multiple
  physical lengths. Neither is a clean diploid ALT certificate.

The Python screens use the corrected maximum-covered-endpoint gap frontier
and pysam. Run from the repository root with the accepted baseline and test
inputs available; output goes under `test_data/tmp_gap_fix57/`. The full runs
use the standard reference, annotated BAM, striped sites VCF and coordinate
GAF with `collect-graph-variation -t 8`, no region restriction, phased VCF and
BAM output. Trial logs and intermediate matrices remain under
`test_data/tmp_gap_fix57/`. No panel expectation is changed.
