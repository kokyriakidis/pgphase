# Confirm a demoted BAM SNP before joining graph recovery blocks

The chr20 graph/BAM output split 61,757,551 G>A from 61,773,799 A>G,
although the same-callset HiPhase run put them in one phase set. The owning
61–62 Mb replay showed why: the latter SNP is a phased BAM recovery row
classified `NOISY_CAND_HET`, so the physical stitch skipped it and tried the
next clean graph SNP at 61,782,778. No MAPQ/base-quality-30 primary read
calls both clean boundary SNPs. Ten independent reads do call the left clean
SNP and the demoted right SNP: six REF/REF and four ALT/ALT, with no
contradictory pair. Both blocks use the same parental orientation.

The left block contains one clean graph SNP, a recovery-adopted graph
insertion, and a complete BAM source. The usual graph-only path checker
requires multiple clean graph SNPs, while the BAM-only checker rejects a
mixed graph/BAM block. The direct graph SNP-to-insertion edge supports their
current orientation, and the source has no weak or quality cut. On the
right, the first clean graph SNP edge has no GAF pair, but both endpoint
SNPs are exact shared rows of one complete, uncut BAM source with consistent
allele orientation. The new fallback uses these two existing path proofs
only after its physical SNP pairs pass the independent read checks.

A broad trial wrongly joined a 37.4 Mb block. Its noisy boundary SNP
agrees with the next clean SNP on only 7 of 11 high-quality paired reads;
four contradict it. That trial lowered full chr20 truth-correct reads from
229,131 to 229,001. The final fallback requires the demoted SNP to agree
with the next clean SNP at a one-sided random-parity probability <=0.001
and to include both allele classes. At 61.7 Mb, 25 of 26 paired reads agree;
only one contradicts. The established indel fallbacks still try the original clean SNP before
the demoted SNP is considered; none may use the demoted allele to bypass
its own confirmation.

The first guarded trial also displaced an existing correct 46.727–46.749 Mb
clean-SNP join. One high-quality primary read calls both clean alleles, so
the final implementation keeps that clean pair ahead of the nearer noisy
candidate. The 37.4 Mb blocks stay split, the 46.7 Mb blocks stay joined,
and the 61.7 Mb blocks join in their owning chunks. All three are in the
window regression panel with exact allele and phase-set checks.

The 61.7 Mb window has 111 truth-scorable overlapping reads. HiPhase tags
110 and places 109 correctly in one block. The old owning pgphase replay
tagged 95, placed 91 correctly across several blocks, and its largest block
placed 49 correctly. The final replay tags 109, places 104 correctly, and
its largest block places 76 correctly. The target VCF SNPs and the already
closed 61,747,506–61,757,551 gap share one phase set without reversing their
parental allele relationship.

| Full chr20 measure | Previous | Final |
|---|---:|---:|
| VCF variant keys | 62,361 | 62,361 |
| Truth-scored tagged reads | 236,866 | 236,880 |
| Truth-correct / discordant reads | 229,131 / 7,735 | 229,144 / 7,736 |
| Read phase sets | 694 | 692 |
| Phased-heterozygote VCF phase sets | 344 | 343 |

Exactly 555 VCF sample fields change, from 61.770 to 62.351 Mb; the variant
keys and all phased genotype strings are unchanged; only PS labels change. The old 49/49-correct and
2,182/2,209-correct read blocks combine into one 2,231/2,258-correct
block, exactly the sum of their correct and discordant counts. The other
14 newly tagged reads contain 13 correct and one discordant assignment.
Runtime logic uses no parental truth.

Reproduction inputs are the shared chr20 fixture. Previous outputs:
`/tmp/pgphase-next-snp-final-full/`. Final outputs:
`/tmp/pgphase-gap61757-ordered-full/`. The final ordering adjustment
produces byte-identical candidate TSV, phased VCF, and phased BAM to the
scored `/tmp/pgphase-gap61757-cleanpriority-full/` run. The wrong
whole-block trial is
`/tmp/pgphase-gap61757-full/`. Focused owning-chunk outputs are under
`/tmp/pgphase-gap61757-cleanpriority-{37,46,61}/`.

Validation: `make -j8 pgphase`, `make unit-tests`, and `make check` pass.
The complete window suite passes 3,721 assertions in 45 test cases; after
the final fallback-order adjustment, the 37.4, 46.7, and 61.7 Mb focused
windows pass 29, 29, and 31 assertions. `git diff --check` is clean.
