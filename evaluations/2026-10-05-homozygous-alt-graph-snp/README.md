# Withdraw a physically homozygous graph SNP that splits a supported block

The closed seam is **36,614,185–36,623,545**, a **9,360 bp** boundary distance
inside HiPhase's eighth-largest chr20 block, **36,137,653–37,393,302**
(1,255,650 bp). The baseline is fix74, SHA256
`9c55c8ac6b67240284980560f308b104fee1dc776b10fd0b89988a10aef251c9`.
The fixed optimized binary is
`a81c2c28df9da273bb5fca1de3ee6b715a96d84489019d2971061b57e1635b41`.

## Selection

`../2026-10-05-next-largest-block/` contains the ranked spans and original-read
seam screening. The first block has an unsupported terminal difference; the
second, third and fourth are covered. Every internal seam in the fifth block
and both seams in the sixth fail the HiPhase >=80% all-original-read rule,
even with the best local orientation. The seventh block's 2,709 bp terminal
difference qualifies, but a normalization trial did not reach HiPhase's
40-correct-read core and lost six previously correct assignments. That trial
was discarded. The eighth block is the next qualifying internal connection
with a demonstrated fix. Its earlier **36,268,558–36,286,778** seam remains
open; this change does not claim to close the entire HiPhase block.

## Cause and fix

The catalog's **36,620,864 G>A** has graph-walk counts REF=11, ALT=44. Those
counts produced a heterozygous graph anchor and a separate phase set. Original
primary MAPQ30 reads instead supply **54 Q20-or-better A bases, zero G bases,
zero deletion calls and no third base**. HiPhase handles this as homozygous ALT.
The existing physical REF-absence validator required a substantial deletion
allele, so it never withdrew the false heterozygote when there was no deletion.
`audit_physical.py` reproduces the certificate without parental truth. Even
using all 977,275 catalog rows as a conservative family size, the adjusted
zero-REF probability is **5.424966e-11**, below the 1% gate.

Zero-deletion homozygous ALT now passes that validator. Any callable REF or
third base still vetoes it, and the existing mixed ALT/deletion rule retains
its 10-call/20% deletion requirement. A catalog flag survives reconstruction.
The initial graph solve establishes the existing gauges; physical
reclassification then retires the homozygous SNP before supplementation and
recovery. Phase sets with surviving heterozygotes keep their gauges and read
assignments. A singleton containing only the withdrawn anchor loses its
unsupported labels and must be reassigned through other evidence.

Reclassifying before the initial solve changed marginal calls in the already
closed 52 Mb owner (four fewer correct reads). The final order preserves all
its original read assignments and accuracy floors. This is a production
ordering rule, with no truth input or coordinate exception. Graph-only runs
have no physical flag and retain their behavior.

The homozygous catalog call is retained, including **AD=11,44 / DP=55**.
Biallelic emission now preserves its validated classification instead of
reclassifying from graph counts: **CLEAN_HOM, GT=1/1, no PS**. These AD counts
remain graph-walk counts; they are not substituted with the physical base
counts used to validate the genotype.

## Native owner and permanent regressions

The denominator is every original primary truth-scorable overlap, including
abstentions. Rescue PS>=1e9 is excluded from the connected correct core.

| 9,360 bp seam | Before pgphase | Fixed pgphase | HiPhase |
|---|---:|---:|---:|
| Original scorable reads | 102 | 102 | 102 |
| Correct | 75 | 93 | 91 |
| Discordant | 26 | 7 | 6 |
| Unphased | 1 | 2 | 5 |
| Largest correct connected core | 31 | 92 | 91 |
| Correct / all original reads | 73.53% | 91.18% | 89.22% |

The overlapping original **36,611,593–36,620,864** control also closes:
**84 correct /5 discordant /2 unphased**, with **83 correct core reads** out
of 91, versus HiPhase's 82 correct/core. Both windows are certified against
the all-read >=80% rule and HiPhase total/core floors. The control's span
expectation intentionally changes 0→1; its measured core floor increases.
Other old floors remain unchanged.

The native 36,000,001–37,000,000 owner retains all **4,013 output reads**.
Correct reads increase **3,637→3,660**, discordant reads decrease **153→134**,
and abstentions decrease **223→219**. Transitions include 20 wrong→correct,
6 abstaining→correct, 2 wrong→abstaining and **3 correct→wrong**. These three
old labels belonged to the invalid singleton SNP phase set; two reads have
ambiguous short repeat insertions and the third carries a high-quality
opposite-parent SNP observation. No truth-based override retains their old
labels. The change improves aggregate accuracy but does not preserve every
previous correct label.

Changed recovery context also removes one earlier MSA call,
**36,608,713 C>CAA** (38 reads, 13/25 graph/recovery counts). There are
**614→613** native VCF records, so this change does not claim complete variant
preservation. The original catalog homozygous SNP survives with its counts.

One/four-thread candidate TSVs, VCFs and read assignments match exactly.
The owning regression asserts the homozygous category/genotype/counts and
absence of PS, correct/core HiPhase parity, spanning, disjoint parental
orientation and relative gauges of the two true SNPs and insertion anchors.
A unit regression checks the zero-REF gate, unchanged neighboring gauges and
retirement of singleton read labels. The old binary fails **12 of 532** assertions in the new owner test; the
fixed binary passes all 532 in **0.38 s** using its verified replay state.

## Full chromosome and validation

The final 67-chunk run creates **36,286,778–37,397,733**, length
**1,110,956 bp**, with 870 phased heterozygotes. Every previous block extent
remains covered. Blocks decrease **260→259**, N50 increases
**904,351→934,592 bp**, and the largest block remains **3,012,193 bp**.
The seam and overlapping control reproduce their native read/core results.
All 5,118 HiPhase target-block alignments match the original geometry, CIGAR
and sequence. Independent 50 kb flanks in the joined pgphase core agree on
parental orientation **208:0** and **195:0**; neither cohort overlaps the
seam or the other flank.

Across the entire 5,118-read HiPhase block, pgphase improves
**4,810→4,833 correct**, with **83 wrong /202 unphased**, and its correct core
improves **2,930→4,252**. HiPhase has **4,906 correct /67 wrong /145 unphased**
and core 4,906. Thus full-block read parity remains short by **73 correct
reads /654 core reads**, in addition to the earlier open seam and leading
endpoint difference. The earlier seam still has **92 correct /1 wrong
/16 unphased**, core 68 of 109, versus HiPhase's 100 correct/core.

All **256,610 output read names** remain. Chromosome totals improve from
**230,612 correct /6,732 wrong /19,266 unphased** to
**230,643 correct /6,705 wrong /19,262 unphased**. The full transition audit
reports **22 old correct→wrong**, **47 wrong→correct**, **6 unphased→correct**
and **2 wrong→unphased**. All losses originated in the two withdrawn singleton phase sets: **3** at
36,620,864 and **19** at **45,859,664 A>G**. The latter also has decisive
physical evidence: **72 high-quality G bases and no A/deletion/third base**.
Its 45,835,000–45,875,000 original-read cohort improves **92→100 correct**,
with **66→58 wrong and 80 unchanged abstentions** out of 238. It remains an
open, inaccurate region; this change does not certify a closure there.
`audit_side_effects.py` reproduces the losses and physical evidence. The net improvement is
31 correct reads, not preservation of every old assignment.
There are **64,484→64,483 VCF records**: in addition to the native missing
C>CAA call, A>G at 64,561,498 is removed and C>T at 29,309,711 appears.
These changes are reported explicitly in `results.json`.

The optimized build has zero warnings. All unit tests, **1,551 predicate
assertions in 47 cases**, HiFi/ONT TSV/VCF golden gates and HiFi thread
determinism pass. The complete window suite passes **16,031 assertions in
four cases**, and all four replay-cache helper tests pass. Current-binary
cache fingerprints and results are audited for all **222 prior native output
labels /117 independent requests**. Existing accuracy floors are preserved;
the two intended span changes have strict all-original-read HiPhase checks.
The newly added expectation stores the measured connected-core fraction
92/102, rather than the total-correct fraction 93/102.

## Reproduction

```bash
make -j8
make unit-tests predicate-tests check
make gap-owner-check GAP='homozygous graph SNP'
make window-tests

/home/kokyriakidis/micromamba/envs/bench-phasers/bin/python \
  evaluations/2026-10-05-homozygous-alt-graph-snp/audit_physical.py
/home/kokyriakidis/micromamba/envs/bench-phasers/bin/python \
  evaluations/2026-10-05-homozygous-alt-graph-snp/audit_owner.py

evaluations/2026-10-05-homozygous-alt-graph-snp/replay.sh
/home/kokyriakidis/micromamba/envs/bench-phasers/bin/python \
  evaluations/2026-10-05-homozygous-alt-graph-snp/audit_block.py
/home/kokyriakidis/micromamba/envs/bench-phasers/bin/python \
  evaluations/2026-10-05-homozygous-alt-graph-snp/audit_panel.py
```
