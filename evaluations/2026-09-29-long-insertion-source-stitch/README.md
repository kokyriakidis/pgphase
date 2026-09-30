# Chr20 long-insertion seam after recovery

The graph seam at chr20:32,234,664–32,246,127 remained open after recovery.
The right boundary is a 65 bp BAM insertion. Its internal candidate position is
32,246,128 but its VCF anchor and seam boundary are 32,246,127. The physical
stitch selected candidates by internal position, dropping this exact-boundary
insertion. Recovery also exposed a new seam after an earlier stitch had merged
the left graph block: using the seam's old phase-set label then missed its
current left SNP.

Two independent primary reads (MAPQ 14 and 15) carry the exact insertion ALT
and call the left graph SNP REF at Q40. Their inserted base minima are Q22 and
Q35. Two nearby apparent insertion REF reads carry other indels in the 65 bp
repeat neighborhood, so their REF calls are ambiguous. With those calls
excluded, the two ALT reads give log odds -5.03383 for the phase relation,
above the 0.01 wrong-parity threshold (log odds 4.59512). The right BAM source
path is complete and has no weak or quality cuts. The left graph block has an
earlier weak edge at 32,207,031–32,233,534, but its SNP suffix from
32,233,534 to the seam is supported. Its earlier 32,215,055–32,233,534 join
was already established.

The implementation now uses the insertion VCF anchor for seam membership,
resolves current phase sets at the boundary, and admits long insertions up to
128 bp to the physical CIGAR check. For insertions over 64 bp, a competing
indel within one insertion length makes an apparent REF call abstain. The
long-insertion bridge requires two ALT molecules, actual inserted-base
qualities, MAPQ at least 5, a 0.01 wrong-parity bound, and a complete right
BAM source. It verifies the left boundary-side graph SNP suffix and orients
the right source against the already established left block.

A trial that split the left block at its distant weak edge closed the new gap
but reopened the previously closed 32,215,055–32,233,534 gap, leaving the
tracked total at 20/23. Retaining the established left block closes the new
gap and keeps the earlier join, giving **21/23** tracked HiPhase-correct
closures. The remaining open targets are chr20:21,377,985–21,395,286 and
chr20:23,421,003–23,445,252.

Full chr20, before → after:

| Measure | Before | After |
|---|---:|---:|
| VCF variant keys | 62,154 | 62,154 |
| VCF genotype strings changed | — | 0 |
| Truth-scored phased reads | 236,854 | 236,855 |
| Truth-correct reads | 229,097 | 229,097 |
| Discordant reads | 7,757 | 7,758 |
| Read phase sets | 700 | 699 |
| Tracked correct gap closures | 20/23 | 21/23 |

Only three VCF sample fields change, all phase-set labels. The cleaned final
replay produced byte-identical phased VCF and BAM files to the successful trial. One additional read is tagged discordantly; three
existing read truth outcomes also change (two correct to wrong, one wrong to
correct). The net truth-correct count is unchanged. The target and earlier
32 Mb joins are asserted together in the owning-chunk regression.

The panel now replays the complete 32–33 Mb owning chunk for this window. In
that wider replay its truth-scored BAM has 1,311/1,520 correct reads (86.25%)
across all phase sets; the dominant block correctly separates 39/52 local
reads (75.0%), compared with HiPhase's 38/52 on the tracked target. The
previous 95% concordance floor came from a shorter padded replay and is not
comparable to the owning-chunk measurement; its new floor is 86%.

## Remaining boundary evidence

At 23.421 Mb, one primary MAPQ-60 read spans the left SNP and right 9 bp
deletion. It calls left ALT at Q17 and the full right REF at Q40; that relation
agrees with the parental blocks, but its estimated SNP error alone is about
2%. HiPhase phases the two boundaries together while leaving the interior
23,437,698 deletion unphased. This is one direct boundary observation, not an
independently corroborated bridge.

At 21.378 Mb, two paternal MAPQ-60 reads span the boundary insertions. One
has a distinct 2 bp insertion near the left 6 bp repeat insertion and Q3 at
the right REF flank; the other has Q10/Q22 near the left REF flank and Q40 at
the right. Neither supplies clean paired REF calls under the current physical
allele check. Counting the distinct 2 bp insertion as the 6 bp site's REF
allele would conflate repeat alleles. These gaps remain split pending a
representation-aware, independently supported link.
