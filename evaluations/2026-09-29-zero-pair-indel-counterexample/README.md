# Zero-pair SNP fallback: why a one-read deletion stitch is unsafe

The physical stitch prefers a clean SNP pair. A trial also called its existing
SNP-to-MSA-deletion helper when a right SNP was present but no primary read
called both SNPs. This changed 313 VCF sample fields across three loci and
joined chr20:58,366,458–58,385,702 in the owning chunk. The same BAM and VCF
were used by pgphase and HiPhase. No production stitch change was kept.

| Full chr20 measure | Accepted baseline | Zero-pair trial |
|---|---:|---:|
| Variant keys | 62,154 | 62,154 |
| Truth-scored phased reads | 236,859 | 236,859 |
| Truth-correct | 229,108 | 228,967 |
| Discordant | 7,751 | 7,892 |
| Read phase sets | 698 | 694 |

The new stitch found one Q40, MAPQ-60 physical pair at each of three
SNP-to-deletion seams: chr20:11,573,074–11,591,586,
20,782,736–20,801,420, and 58,366,458–58,385,703. All had the same
quality-weighted log odds (absolute value 6.80953). The 11.57 Mb join changed
no truth score; the 20.8 Mb join incorrectly merged a large block, and the
58.366 Mb join put the two boundary ALT alleles on opposite haplotypes.
HiPhase puts those two ALTs on the same haplotype. The owning-chunk gap panel
reported 81% local read concordance and 40% one-block correctly separated
reads for the trial, compared with HiPhase's 82% correctly separated reads.
A span check alone would have accepted the wrong allele relation.

At 20.8 Mb, the injected one-base deletion overlaps an unphased six-base
graph repeat deletion. A trial vetoed overlapping repeat rows and blocked the
join in a 20–21 Mb replay. A full chr20 replay required comparing within-chunk
coordinates without relying on matching candidate `tid`; this blocked the
20.8 Mb join but still left 7,815 discordant reads because of the wrong
58.366 Mb join. The accepted baseline remains better.

At 58.385 Mb, two MAPQ-60 primary reads span the left SNP and deletion locus.
One calls the SNP REF and an exact one-base deletion at 58,385,703. The other
calls the SNP ALT and places the same one-base deletion 20 bp downstream in an
A homopolymer. Reference-edit comparison says the deletions are sequence
equivalent. Its first retained A has base quality 17, so pgphase's Q30 shifted-
deletion filter abstains; the exact read alone then gives the wrong parity.
MAPQ 5 adds no additional spanning pairs. The present stitch does not try
this zero-pair fallback, which is the safe result. The gap-window regression
now rejects a joined VCF block that puts these exact ALT alleles on opposite
haplotypes.

A correct future join needs independent allele support or a quality-aware
model for both equivalent deletion placements, with a chromosome-scale truth
check. All trial binaries and outputs are under `/tmp/fallback-full/`,
`/tmp/repeat-veto-full/`, and `/tmp/repeat-veto-no-tid-full/` while available.

Validation of the accepted tree: `make -j4`, `make unit-tests`, `make check`,
and `make window-tests` pass (3,307 assertions in 42 cases). A one-window
replay using the rejected trial binary fails both the parental switch check
and the new exact boundary-orientation assertion.

Additional source trace: the two MAPQ-60 reads spanning the 58.366 Mb SNP
and 58.385 Mb deletion end their BAM sub-solve allele profiles at candidate
index 282, immediately before the deletion at index 283, despite physically
spanning it. The exact-CIGAR caller classifies the SNP-REF read as deletion
ALT and the SNP-ALT read as deletion REF. The second read's deletion is shifted
20 bp within the A run, so its exact-CIGAR REF call is a representation
artifact. Raw CIGAR backfill cannot safely repair this gap; a local MSA or
sequence-equivalence call with explicit base-quality handling is needed.

## Clean-SNP recheck (2026-09-29)

The repeat deletion is a poor orientation control. Among MAPQ-30 primary
reads, exact deletion CIGAR calls occur on both parents (10 maternal, seven
paternal), and REF calls likewise occur on both (23 maternal, 21 paternal).
The earlier regression compared the left SNP with this deletion and therefore
encoded an uncertain repeat-allele relation. The next clean SNP at
58,391,091 is the stable right-block control: the left 58,366,458 ALT occurs
on 26 paternal and zero maternal reads, while the right SNP ALT occurs on
32 maternal and one paternal read in the parental allele audit. HiPhase puts
their ALT alleles on opposite haplotypes (`0|1` left, `1|0` right). A rejected
trial oriented both ALT alleles onto the same haplotype. The window regression
now checks the two clean SNPs whenever a future stitch assigns them one PS.
This check uses only the expected allele relationship at runtime; parental
truth is used to establish the test expectation, never for phasing.

At 10,727,690–10,746,628, the broad recovery-window overlap count also
overstated direct bridge coverage. Only two distinct MAPQ-30 primary reads
span the exact left SNP and right insertion coordinates, and neither has a
local insertion CIGAR call. The BAM source PS split is therefore not a
missing 35-read allele bridge.

A 21–22 Mb owning-chunk replay under `/tmp/pgphase-21159-audit/` shows another
shared-source-PS trap. Standalone BAM assigns both exact 21,159,070 `C>T` and
21,172,487 `CA>C` rows PS 21,148,338. In its read-observation matrix, only
one read profile covers both candidate indices, and the deletion allele is
unknown on that read; there are zero callable boundary pairs. Recovery keeps
the left SNP in PS 21,148,338 and places the deletion with the right graph
block in PS 21,172,501. The apparent BAM PS continuity is an inherited label
across a weak source edge, not an observed allele path.
