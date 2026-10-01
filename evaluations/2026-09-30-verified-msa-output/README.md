# Preserve verified BAM MSA sites in graph output

At chr20:35,342,608, targeted recovery appends the BAM-derived
`G>GAGATAGAT` insertion inside the still-open 35,328,965–35,347,817 gap.
Its sub-solve calls `NoisyCandHet` after MSA and alignment verification,
with 2 REF and 9 ALT reads and `LOW_COV` as its initial category. The
inserted candidate remains phased with the right-side deletion. Final graph
output used to reapply its allele-fraction discovery threshold, classify it
`LOW_AF`, and omit it. This is an emission bug: it erases a source candidate
that has already passed the BAM sub-solve's classification. In an
independent primary-read CIGAR audit requiring MAPQ at least 30 and high-quality
inserted bases, 30 ALT calls have paternal truth and 26 REF calls have
maternal truth; two paternal reads call REF and seven reads have another
nearby indel. The physical evidence supports the recovered heterozygote
despite its 2/9 sub-solve depth balance.

Preserving the verified source category exposed a second bug. Output builds a
new `CandidateVariant` for each ALT but did not copy `msa_verified` and
`alignment_verified`. Duplicate selection therefore judged valid phased
`NoisyCandHet` rows unverified and chose deeper unphased graph repeat rows at
three other positions in the 35–36 Mb chunk. The output reconstruction now
keeps those proof flags and the `bam_injected` origin. The owning chunk gains
two exact phased keys (35,149,175 and 35,342,608), loses none, and retains all
previously phased rows. Its phased BAM is byte-identical.

| Full chr20 | Before | After |
|---|---:|---:|
| VCF variant keys | 62,269 | 62,361 |
| Added / removed keys | — | 92 / 0 |
| Changed existing sample fields | — | 0 |
| Truth-scored phased reads | 236,866 | 236,866 |
| Truth-correct / discordant reads | 229,131 / 7,735 | 229,131 / 7,735 |

The full phased BAM is byte-identical. The accepted input is
`/tmp/pgphase-run-bridge-final-full/`; the measured output is
`/tmp/pgphase-verified-msa-output-full/`. The final code also preserves the
BAM source's `INIT_CAT` (rather than replacing it with `NoisyCandHet`); that
metadata-only refinement was verified in the owning-chunk replay after the
full measurement.

The boundary itself remains split. A fresh 35–36 Mb owning-chunk replay
shows the new insertion in the right phase set and the left insertion in a
separate phase set. The final physical stitch currently revisits a newly
exposed insertion seam only for insertions over 64 bases. The eight-base
seam is therefore omitted, but simply broadening that scan would be unsafe:
the read-observation edges inside each adjacent repeat block fail the
independent two-half path check. Primary-read CIGARs put seven ALT calls on
left HP1 and five REF calls on left HP2, with one discordant REF call on HP1;
this orients the insertion locally, but does not certify both full blocks.
The exact variant-presence and PS relation are pinned by a separate owning-
chunk regression. The existing narrow-window accuracy floor and split-span
expectation remain intact.

Validation: `make -j4`, `make unit-tests`, `make check`, the focused
owning-chunk regression, and the complete `make window-tests` panel pass
(3,547 assertions in 45 cases).
