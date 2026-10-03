# Complete recovery evidence and block continuity

Target: chr20:50,548,245–50,562,066, outside the centromere.

## Bugs and fixes

1. Complete-source evidence was selected using injection ownership. A BAM
   block beginning at the right gap boundary contributed no strictly interior
   candidate and disappeared from the immutable matrix. Retain all eligible
   source heterozygotes and callable observations. Context-only blocks preserve
   their original local PS but receive no live imported label, preventing a
   numeric collision from adopting graph ownership.
2. The graph continuity check required one molecule to call both ends of a
   hundreds-of-kilobases block. Certify its overlapping read-observation edges
   and all coordinate cuts instead. Require both haplotypes, the existing
   minimum support and more agreement than conflict. Supported edges must
   connect every row: interleaved A–C and B–D components cannot certify one
   haplotype gauge. No read tags substitute for allele observations.
3. A BAM source certificate was reused after other graph rows had been
   absorbed into the same live PS. Restrict it to blocks actually covered by
   the saved source rows; certify an expanded current block independently.

## Target evidence

The immutable matrix now contains both source blocks: 27 eligible left rows
and 91 right rows. All 12 MAPQ60 physical spanning reads survive, versus 2
before. `snapshot-audit.json` records the measured retention.

The owning 50 Mb output preserves all read tags and VCF rows, with 3,877
truth-scored reads, 3,867 correct and 10 discordant. The gap remains open.
The restored right calls do not create missing left-boundary allele calls:
all 12 source rows were unknown at both complementary deletion candidates.
The two callable left insertion profiles supply only one haplotype and no
clean-SNP pair. The CIGAR audit also contains repeat-length disagreements,
absent deletions and low-quality bases. Retaining coverage does not certify
an allele or a whole-block orientation. Do not force the old 6/6 read-label
aggregate into a stitch.

## Verification

The original implementation fails the connected-overlapping-chain test and
accepts unsupported rows covered only by an unrelated source certificate.
The initial cut-only repair passes those tests but incorrectly accepts
interleaved disconnected components; the added topology regression reproduces
that defect before the connectivity guard. The accepted implementation must
pass all three negative controls as well as the existing physical SNP,
per-block HP, source-frame and atomic-stitch checks.

The production algorithm contains no truth labels, competitor labels or
coordinate-specific decisions. Parental truth is used only to score outputs.

Final full chromosome (`full-parity.json` and `block-audit.json`): exact read-tag
and VCF-row parity, 237,199 truth-scored, 230,072 correct, 7,127 discordant
(96.995350%), 63,630 keys and 332 VCF blocks. Span N50 is 774,189 bp and the
coordinate panel remains 87/101 connected. No phased SNP is lost, no old block
changes gauge and no additional gap closes. Final wall time is 288.07 s at
8 threads, concurrent with panel validation; no runtime comparison is implied.
All standalone units and 1,259 predicate assertions/44 cases pass.

Final window validation passes 7,132 assertions/71 cases, using 105 fresh native
requests from the final binary plus fresh uncached diagnostic runs. Cached
outputs are verified against the final binary SHA256, input file metadata and
all output hashes. No floor, ceiling or span expectation is weakened. The final
binary SHA256 is
`629badc0b446294ee4e0f1c4b4a43c738406ac7dbed9023995b6796e3fb0f05f`.
