# gap_lab.py -- one command per gap: evidence, phasing, stitch, gate

Every gap diagnosis before this was assembled by hand, so each experiment's
verdict depended on which checks were remembered. `gap_lab.py` runs the whole
loop for one gap and prints a pass/fail:

```sh
python3 evaluations/2026-09-16-gap-lab/gap_lab.py \
  --gap-left 48176830 --gap-right 48229446 --flank 150000 \
  --truth-map /tmp/truth_hap.tsv --seam-window 500 --max-new-discordant 5 \
  --work /tmp/gaplab/g48176830 --out verdict.json
```

Stages:

1. **Baseline** -- the pipeline as shipped over the gap plus a flank. Its blocks
   define the gauge, and any change has to beat it.
2. **Evidence** -- every het candidate the BAM channel finds inside the interval,
   with its category and, against read truth, how well its allele partition
   segregates. A site that does not segregate is not evidence, whatever its
   category says.
3. **Gap phasing** -- the BAM channel over the same window. This is the existing
   machinery, not a reimplementation, so whatever it produces is what a scoped
   change inside the pipeline would have to reproduce.
4. **Stitch**, in two stages: compose the gap channel's own blocks with each
   other, then link each composed frame to the flanks. Composing first matters --
   asking one block to reach both flanks is the all-or-nothing test recovery
   already fails.
5. **Gate** -- read-level scoring before and after. A stitch that turns
   concordant reads discordant is a regression however many reads it adds.

`--reuse` skips arms whose outputs exist, so iterating on the stitch costs
seconds rather than a pipeline run.

## Three things it found in its first run, all of which were wrong in my hand analyses

**A flank can hold sites and no reads.** The block nearest the gap on the left,
`48162480`, has 2 sites and **zero tagged reads** -- the read-tagging margin drops
reads that observe too few of its sites. Nothing can link to it by read identity,
so flanks are now chosen by read support and read-less blocks are reported as
skipped.

**Two blocks that split at the same position share no tagged reads.** A read
carries at most one phase set, so tag-identity voting returns n=0 between blocks
that break together -- which is the state of every gap examined here. The alleles
are still present: a read crossing the seam observes sites in both blocks, and
each block's own phased genotypes say which haplotype each allele belongs to.
That vote is now the fallback, and it is what composes this gap.

**The seam is not the block.** The first version of that vote required a read to
span every selected site, which reaches tens of kb back into each block while the
seam here is 220 bp wide; it reported zero crossing reads. The test is now the
seam point plus a margin, with a minimum number of site observations per side.

## Result on chr20:48,176,830-48,229,446 (52.6 kb)

The deficit gap hiphase spans at 100.0% over 252 reads.

| stage | result |
|---|---|
| baseline | 3 blocks, 893 tagged reads; does not span the gap |
| flanks (read-supported) | left `48043584` (28 sites), right `48229446` (182 sites); `48162480` skipped, 0 reads |
| BAM evidence in interval | 11 het sites, **11 phased by the BAM channel**, 7 informative against truth |
| compose | `48147227` + `48229446`: 42 allele voters, votes **[15, 0, 0, 27]** -- unanimous, 58 reads cross the seam -> composed, flip 0 |
| link | composed frame (231 sites, 48,147,227-48,377,087, 820 reads) links the **right** flank at n=661, votes [318, 0, 0, 343]; **no link to the left flank** (n=0) |
| gate | tagged 893 -> 1052, concordant 892 -> **1051**, accuracy 99.89% -> 99.90%, **0** concordant->discordant, **159 newly tagged, all concordant**, 0 lost |
| verdict | **PASS** -- not closed, right flank extended, +159 concordant reads |

So the gap does not close, and for a reason already established: the left side
carries two hard linkage breaks (`48,096,582->48,123,657`, 27.1 kb, and
`48,123,657->48,147,230`, 23.6 kb) with zero reads covering the flanking sites at
72x coverage. No read-based method crosses those.

But the composed frame covers the whole competitor-spanned interval and adds
**159 correctly phased reads with no read flipped**, which the shipped pipeline
leaves unphased. The two operations that produced it are exactly the two the
pipeline has no mechanism for: an allele-level block-to-block vote, and a
one-sided extension of a flank by a composed frame.

## What to use it for

The stitch here runs outside the pipeline, so this is a measurement rig, not a
fix: it says what a scoped change would have to reproduce, and its verdict is the
regression test for one. Next candidates, in order of evidence:

- Run it across the 31-gap accurate-deficit list
  (`evaluations/2026-09-16-current-deficit/deficit_scored.tsv`) to size the total
  recoverable coverage before touching C++.
- Then the scoped pipeline change this and the previous experiment point at:
  admit BAM-discovered noisy candidates to phasing inside gap intervals only, and
  give recovery an allele-level block-to-block link so a frame that reaches one
  flank still extends it.
