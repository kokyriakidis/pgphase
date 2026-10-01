# Readless graph SNP block at chr20:30.673 Mb

The baseline full chr20 VCF has PS 29,440,378 at `30,673,476:G>C` and
PS 30,610,946 at `30,673,709:T>A` through `30,676,069:T>C`. The right
PS contains 51 phased VCF rows and no tagged reads. HiPhase on the same
pgphase callset places both boundary variants in one phase set, with 37/38
truth-scored nearby reads correct.

The BAM recovery audit covers 30,660,651–30,791,916, but the physical stitch
selects an older seam and has no clean SNP exactly at the noisy left
boundary. Three independent MAPQ-at-least-30 primary reads span the 233-bp gap and
call at least two clean SNPs per side with unanimous within-read haplotype
votes (4/5, 5/12 and 4/8 left/right SNP calls). They agree on the block
orientation. A fourth crossing read calls 2 left SNPs on one haplotype and
4 on the other, so it cannot vote. All 16 consecutive clean-SNP edges of
the 2.3-kb right block have at least one consistent Q30 physical read and no
contradictory pair.

A readless-right-block connector uses those physical calls and a 0.001
combined wrong-parity bound, then transfers the whole right PS with its
relative orientation. The owning 30–31 Mb chunk closes the gap. Its phased
BAM stays at 3,559/3,573 correct truth-scored reads; the VCF boundary and
terminal allele relationship matches parental read evidence. The committed
window panel adds `30673476-30673709`, asserts `spans=1`, and checks those
exact allele relationships.

Full chr20 changes exactly 51 VCF sample fields, all in
30,673,709–30,676,069. Variant keys remain 62,361. The phased BAM SHA256
is unchanged (`764d0eb80ee5b042556c688e49296e5824dda99f001e60f1f94320721ad9dcc1`):
236,866 truth-scored reads, 229,131 correct, 7,735 discordant, 695 read
phase sets. This improves VCF block continuity without changing read
coverage or read accuracy.

A scan of the currently open HiPhase-correct gaps found only this seam with
four crossing reads that each call at least two clean SNPs on both sides.
At 21.159–21.172 Mb, nine versus seven one-SNP read votes disagree and no
read calls two clean SNPs on each side; at 64.138–64.140 Mb, the clean-SNP
flanks have only one-site support. The orphan-block rule therefore does not
extend to those weaker cases in the chromosome trial.

A separate uniform MAPQ-30 source-path filter trial split five additional
read phase sets, lost one truth-correct read, and closed no tracked gap. It
was reverted. The 36.343–36.354 Mb BAM source cut was also audited: all 26
paired reads were MAPQ 60, with 14 supporting and 12 contradicting the
insertion-to-deletion path. That cut is not a low-MAPQ artifact.

Verification: `make -j8 pgphase`, `make -j8 unit-tests`, `make check`,
and `make window-tests` all pass. The full window suite reports 3,599
assertions in 45 test cases; the focused new panel row passes 29 assertions.
