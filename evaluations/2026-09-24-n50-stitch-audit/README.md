# Chr20 phase-block stitch evidence (2026-09-24)

## Question and data

Can a Longcalld-style count of reads spanning neighboring blocks increase the
chr20 graph+BAM phase-block N50 without losing the read-truth accuracy lead?
The baseline is `/tmp/pgphase-direct-snp-rescue/tight_chr20/` against
`test_data/derived/chr20_truth_hap.tsv`: 236,409 phased reads, 227,830
truth-correct, 8,579 discordant; 445 VCF phase blocks, N50 460,310 bp.

Longcalld's `iter_update_var_hap_cons_phase_set` starts a new within-chunk
phase set only if both agreement and conflict counts are below two. Its
`flip_variant_hap` joins overlapping chunks on any non-tied net read vote.
Pgphase already counts same/cross read links in `RecoveryBlockGaugeVote`,
checks statistical parity in the graph/BAM stitch, and uses a non-tied
Longcalld-style vote for ordinary chunk boundaries. The graph/BAM fallback
also guards against reusing an unvalidated middle block.

## Direct-spanning-read probe

For each adjacent VCF phase-set pair with at most a 30-kb gap, inspect up to
eight clean SNPs from the last 20 kb of the left block and the first 20 kb of
the right. A primary read votes only if it observes informative alleles in
both blocks at base quality >=20. Its physical span into each block is the
intersection of its reference alignment interval with that block's VCF span.
This is a diagnostic on the source BAM, not a replacement for the production
candidate/observation model.

| Boundary | Informative reads | Reads extending >=2 kb into each block | Truth result |
|---|---:|---:|---|
| 2.525 Mb, 173-bp gap | 19/19 same parity, both hap links | 13/13 same parity | Consistent with parental orientation; likely missed join |
| 51.263-51.271 Mb, 8.692-kb gap | 28/28 same parity, both hap links | 15/15 same parity | **Wrong as a whole-block join**: left block switches internally near 51.24 Mb |

At 51.17 Mb, the left phase set's reads before 51.23 Mb overwhelmingly use
one parental orientation, while reads from 51.24 Mb onward use the opposite.
Reads crossing its right boundary validate only the local end. Requiring a
longer physical read span does not catch that switch: an additional 3-kb
minimum still leaves 8/8 apparent supporting reads at this bad boundary.

## Why HiPhase is correct at 51.24 Mb

The frozen DeepVariant/HiPhase `blocks.native.tsv` records separate source
blocks 180 (50,127,297-51,235,063) and 181
(51,262,081-52,267,263). Its VCF ends one PS at the clean SNP 51,235,063
and starts the next at the clean SNP 51,262,081. The same split occurs when
HiPhase is run on pgphase's own older VCF catalog. In the annotated source
BAM, 56 primary reads cover each endpoint, but **zero primary reads cover
both**; the farthest primary read from the left endpoint ends at 51,255,234.
HiPhase's block generator requires a mapping or qualifying supplemental
alignment to bridge consecutive included variants (default one spanning read).
The two source block indexes in the frozen output show it broke the problem
before phasing; no unsupported phase relation was asserted across this 27-kb
interval. Its subsequent clean SNP pair 51,262,081-51,270,774 has 31 primary
reads physically spanning both sites, and both variants share PS 51,262,081.

In contrast, pgphase outputs clean SNPs 51,235,063 and 51,262,081 in the same
PS 51,165,532, then opens a new PS at 51,270,774. It has therefore already
asserted an unsupported long-range relation *inside* its left block. The
28/28 local boundary vote is consistent with HiPhase's second block, but
merging pgphase's entire left block would carry its earlier, opposite
orientation along. Among tagged reads fetched over 51.05-51.50 Mb, the two HiPhase PSs are
824/833 and 1024/1033 truth-correct, respectively, after independently
orienting each PS against parental truth. The observed zero-spanner interval explains this
specific split; the broader rule remains a physical/read-allele connectivity
check, not a truth-based decision.

## MAPQ, clean-site, and graph-alignment audit

The endpoints 51,235,063, 51,262,081, and 51,270,774 are all pgphase
`CLEAN_HET_SNP` calls. There is no other heterozygous pgphase VCF call
between the first two. In the actual pre-stitch hybrid `recovery-input`
matrix, the intervening rows are three homozygous candidates and one
unphased deletion, so none supplies a phased heterozygous step across the
27,017-base gap. The `flags5004` diagnostic is emitted by the BAM subsolve
with the same dump prefix and does not describe the graph-only state.

| SNP pair | Primary BAM coverage at endpoints | Primary reads covering both | MAPQ of spanning reads | Clean SNP observations, BQ >= 20 | Alignment reach beyond both SNPs |
|---|---:|---:|---|---:|---|
| 51,235,063 -> 51,262,081 | 56 / 56 | 0 | none | 53 / 55 individually; 0 joint | farthest left-site read ends 51,255,234 (6,847 bases short of right SNP) |
| 51,262,081 -> 51,270,774 | 56 / 63 | 31 | all 60 | 28 joint, all 28 same-parity (11 and 17 from opposite haplotypes) | left flank 214 / 4,055 / 14,279 bases; right flank 164 / 3,085 / 16,091 bases (min / median / max) |

MAPQ at the first SNP is 54 reads at 60, one at 52, one at 16; at the
second SNP it is 54 at 60, one at 59, one at 52. The second pair's 31 physical
spanners all have MAPQ 60. Of those 31, 28 observe the expected allele at
both clean SNPs with base quality >=20; the other three have a low-quality or
missing/other call. The minimum base quality across the two SNPs in the 28
informative reads is 22 (median 40). Requiring MAPQ >=30 or MAPQ >=60 retains
all 28 informative connections. Requiring >=2 kb aligned *outside each SNP*
retains 15. The coordinate-indexed GAF has exactly the same read-name bridge
counts, 0 and 31, with MAPQ 60 for all 31 records over the second pair.

The regional replay's pre-stitch `recovery-input` matrix shows three separate
phase sets at these three SNPs. Its `recovery-final` matrix merges the first
pair across the read-disconnected interval and leaves 51,270,774 separate.
Thus the unsupported relation in this output is introduced by the recovery
stitch, and it misses the stronger clean-SNP connection to the right. An
instrumented replay found the precise branch: the graph-only gauge path called
`merge_phase_sets_in_place` even when `strongest_phase_set_edge` found no local
link. Its broad graph/BAM gauge votes were decisive because reads on either
side independently matched the BAM solve, although no molecule connected the
two graph blocks. The retained fix defers that gauge-only join when the right
block is needed as the left flank of the next recovery seam. In the short
51,185,063-51,312,081 replay, the original code placed all three SNPs in one
PS; the fixed code leaves 51,235,063 apart and puts 51,262,081 with
51,270,774. The production 51-52 Mb chunk also splits the unsupported first
pair; its larger second seam still does not join the latter two SNPs.

A broad guard requiring an allele edge for every graph-only gauge join was
rejected. On matched full chr20 at `--min-read-margin 2`, it kept 188,341
phased reads but increased truth-discordant reads from 3,435 to 3,812 and
reduced VCF N50 from 482,085 to 449,827 bp. The retained next-seam guard
keeps the same 188,341 phased reads, **reduces** discordant reads to 3,018,
and gives N50 456,233 bp (457 VCF blocks versus 439). These are matched
runs with the same input and command; the change in read accuracy comes from
changed phase-set/read assignments, not a change in tagged-read count.
At the default read-output margin, the exact matched A/B changes 236,404
phased / 227,446 correct / 8,958 discordant to **236,520 phased /
228,184 correct / 8,336 discordant**. VCF N50 is 482,085 -> 456,233 bp.
The older saved 236,409-read baseline used a different command state and is
kept as historical context, not used for this matched attribution.

A chromosome-wide diagnostic applied a stricter direct-SNP rule (both hap
links, two-sided binomial p<=0.01 overall, p<=0.05 in each of two read-name
folds, and conflict <=25% of supporting reads). It selected 27 adjacent joins;
two disagree with the independently established parental relation. Simulated
phase-set merges changed N50 from 460,310 to 491,812 bp but reduced
truth-correct reads from 227,830 to 227,514 (316 lost). This simulation
assumes the original per-read HP assignments remain fixed and reorients each
resulting merged PS for its best parental agreement. It is a diagnostic upper
bound on the quality of simply relabeling phase sets, not a full pipeline run.

## Conclusion

A same-versus-flipped count is useful to orient a *local boundary*, and
physical overlap length can filter reads that barely touch an edge. Neither
proves the phase of an entire middle block when no molecule spans its ends.
Before relaxing the existing reuse guard, a joining algorithm must validate
an unbroken internal phase path or split the middle block at a weak/internal
switch. The direct-SNP count relaxation was not retained; the narrow
next-seam guard above is a separate correctness fix.
