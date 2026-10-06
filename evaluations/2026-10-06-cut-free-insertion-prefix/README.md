# Cut-free insertion prefix at 36.268 Mb

The 18,220 bp internal gap **36,268,558–36,286,778** belongs to the eighth-largest HiPhase block, **36,137,653–37,393,302** (1,255,650 bp). The original alignments contain a usable short-insertion bridge. pgphase incorrectly demands a complete original BAM source even though its only weak cut, **36,343,992**, lies beyond the exact shared graph SNP that anchors this insertion to the existing right block.

## Selection

[Selection](selection.json) preserves the fresh ranking and all larger uncovered seams. Larger split blocks at ranks five and six have HiPhase accuracy below 80% at their internal seams. The seventh-ranked block has a separate 2,709 bp terminal extension; this repair targets the largest qualifying **internal split**, not that endpoint. The eighth-ranked block retains a separate 2,430 bp leading interval with poor HiPhase accuracy.

## Repair

- Certify only the original BAM prefix from the recovered `36,286,778 T>TG` marker through its first exact shared clean graph SNP, `36,317,511 G>C`. Every intervening anchor in the live run needs unique original-source provenance and the same allele gauge. Any internal weak or quality cut vetoes the prefix. A complete right graph path separately certifies the rest of the established block.
- Accept one high-quality REF molecule and two high-quality ALT molecules with the existing MAPQ/base-quality 30 and wrong-parity bound of 0.001. Their weighted log odds are −10.3839; the required magnitude is 6.9068. The already joined left block retains its supported prefix and independently verifies its boundary-side SNP suffix after the weak GAF edge `36,247,421–36,268,291`.
- Retain the physical marker and both block gauges until recovery and rescue finish. Revalidate the final union before promoting reads from original primary alignments. Missing imported calls are not contradictory calls.
- A one-base homopolymer marker can use the net length of one physical CIGAR edit inside the reference run when repeat-base quality or one-base length slippage prevents an exact call. All aligned and inserted repeat bases must match the motif; a competing edit or a length more than one base from the nearest allele vetoes it. Matching high-quality anchors outside both sides of the run certify mapping/base confidence. Other informative loci and equivalent graph descriptions still veto disagreement. Variant descriptions, counts and observations are preserved.

Truth is used only by evaluations and tests.

## Native owner

The optimized native **36–37 Mb** replay scores **100 correct, 1 discordant and 8 unphased of 109 original primary truth-scorable overlaps**: **91.7431%**, matching HiPhase's **100/1/8**. All **100** correct reads belong to one connected core; before the fix only **68** did. The before-fix total was **92** correct.

Across all **4,013** primary owner reads, correct assignments rise **3,660 → 3,668**, discordant stay **134**, and unphased fall **219 → 211**. No previously correct or phased assignment is lost. All **613** owner VCF records retain their alleles, counts and filters. One- and four-thread candidates/VCF bytes and BAM assignments are identical. [Owner audit](owner-results.json).

The eight remaining unphased gap reads lack a usable heterozygous marker; HiPhase also leaves those eight unphased.

## Reproduction

```bash
make -j8
make unit-tests predicate-tests check
make gap-owner-check GAP=insertion-prefix
make window-tests
bash evaluations/2026-10-06-cut-free-insertion-prefix/replay.sh
/home/kokyriakidis/micromamba/envs/bench-phasers/bin/python evaluations/2026-10-06-cut-free-insertion-prefix/audit_owner.py
/home/kokyriakidis/micromamba/envs/bench-phasers/bin/python evaluations/2026-10-06-cut-free-insertion-prefix/audit_block.py
```

The new committed panel row, HiPhase measurement, certified closure, exact span expectation, marker checks and native-owner parental-orientation regression all cover this seam. The owner regression passes **589 assertions**; the original binary fails **11** of them. With saved replay state the owner check completes in **0.73 seconds**. The in-memory predicates pass **1,551 assertions in 47 cases**; all unit tests and HiFi/ONT golden/determinism checks pass. Optimized compilation introduces no warnings.

## Whole-chromosome result

The joined block spans **36,140,083–37,397,733**, **1,257,651 bp** with **900** phased heterozygotes. The older 1,380 bp orphan inside its left portion remains separately labelled; it does not split the spanning main phase set. All older block extents remain covered. Block count falls **259 → 258**; **N50 increases 934,592 → 944,186 bp**. Largest block stays **3,012,193 bp**.

Across all **256,610** primary output names, correct assignments rise **230,643 → 230,651**, discordant remain **6,705**, and unphased fall **19,262 → 19,254**. No previously correct or phased assignment is lost. All **64,483** VCF records preserve alleles, counts and filters; changes are consistent phase gauges. [Full preservation audit](full-preservation.json), [HiPhase comparison, geometry/sequence/quality/MAPQ identity and parental flanks](results.json).

The repaired gap matches HiPhase. The entire HiPhase block still has **65 fewer total correct reads and 162 fewer connected correct core reads**, and its **2,430 bp leading extension remains open**. This measurement does not certify that remaining interval or claim full-block read parity.

The full gap suite passes **16,115 assertions in four cases**, plus all four replay-cache helper tests. Existing expectations were preserved; the new window adds a measured exact span and stricter read/core floors. Completed native states remain reusable for subsequent development.

[Panel audit](panel-audit.json) verifies the final binary for **223 native output labels / 117 independent requests**. Only **five** labels change, all sharing the repaired owning context. Every earlier correct assignment, variant key and regression floor is preserved; no other full-chromosome block seam closes as a side effect.
