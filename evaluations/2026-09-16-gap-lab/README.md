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

## Correction: the first version did not run the machinery in the gap

The first version ran the phasing machinery over the gap plus a 150 kb flank and
took its flanks from that same re-run. Both are wrong, and they were my errors:

- **The machinery must run on the gap interval.** Over gap+150 kb its blocks are
  window-wide and inherit the window's own breaks -- the block it produced spanned
  48,147,227-48,229,226, which is not the gap's phasing at all. The gap arm now
  runs on the gap plus `--gap-margin` (5 kb) of read context, so its blocks are
  gap-local, the same shape as a recovery proposal.
- **The flanks are the frozen chromosome-wide blocks.** Re-solving them in a
  window can move the very boundaries under test. `--baseline-vcf` /
  `--baseline-bam` now take the real pipeline run's outputs and the windowed
  re-run is only a fallback.

Fixing those exposed two more, and the gap closed once all four were fixed.

## Four mechanisms it found, each of which had silently broken a hand analysis

**A flank can hold sites and no reads -- and it is still the flank.** The block
nearest this gap on the left, `48162480`, holds 2 sites and **zero tagged reads**:
the read-tagging margin drops reads that observe too few of its sites. Choosing
flanks by read support instead picked a block 80 kb further out, across a stretch
holding no phased sites at all, and then reported no link. Flanks are now the
blocks holding the phased sites nearest the gap, and a read-less flank is linked
through its own genotypes rather than its tags.

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

**A read-less flank is invisible to the read-level gate.** Every read in the
merged block comes from the other side and keeps its own relative labelling, so a
wrong orientation across that link cannot flip anything and the gate passes
regardless. The applied links are therefore validated separately: each side's
genotypes are read off the alignment, compared with the read truth, and the
applied flip is checked against the two sides' truth haplotypes. `--pass`
requires that check as well as a clean gate.

## What the link uses from the adjacent blocks

The first version took eight sites nearest the seam. That was wrong twice over,
and the corrections were the same insight applied twice: an adjacent phase block
is not a few sites near its edge, it is a complete labelling -- every read it
contains with a haplotype, and a genotype at every site it phases -- and all of
that is already integrated from the whole block's evidence.

So the flank link now draws on three sources, in order of directness, and they
are required to agree:

1. **Shared phased sites.** The gap arm re-discovers sites inside the flank's own
   span -- here the flank's two sites are 48,162,480 and 48,176,830 while the gap
   arm phases 48,173,317, 48,173,989, 48,173,990 and 48,176,830 -- so both
   labellings frequently phase the same variant. Comparing the two genotypes at a
   shared site gives the relative orientation exactly, with no bridge read and no
   site selection at all.
2. **The flank's own read assignments** against the frame's tags-or-genotypes. A
   phase set already carries every read it contains with a haplotype; there is
   nothing to re-derive on that side. The earlier version required a tag on
   *both* sides, which discarded every read the flank had phased but the gap arm
   had not -- most of the population at a seam.
3. **Alleles on both sides**, for a flank that carries no reads at all. The block
   adjacent to this gap is exactly that case: 2 sites, zero tagged reads.

`--seam-span` and `--seam-sites` are now off by default: the whole adjacent block
is offered and the bounding is per read, since a read votes on the sites it
overlaps and silently ignores the rest. A global distance limit only duplicates
that, and its correct value is gap-dependent -- a site-sparse flank needs a wider
one -- so a fixed default would quietly decide outcomes. They remain as runtime
knobs. The earlier sweep that motivated them is still instructive: capped at eight
sites the left flank offered one site and never linked, at a 10 kb limit the gap
managed only an extension, and only with its whole 2-site block in view did it
close.

**Correction to that reading.** The eight-site arm cannot be the reason the left
flank offered one site: the cap was `sorted(left_sites)[-n_sites:]` with
`n_sites = 8`, a pure count cap with no distance term, and this flank has only
**2** sites in total. Two is fewer than eight, so the cap never bound, and
whatever reduced the flank to a single site in that arm was not it. The cause was
not established, and the conclusion the sweep is quoted for -- that a fixed
default would quietly decide outcomes -- rests on the 10 kb arm, where the limit
does bind.

Sources that fire must not contradict each other. A shared-site comparison and a
read vote are independent evidence about one orientation, so a conflict means one
of the two blocks is internally wrong near the seam, and the rig refuses to link
rather than picking a winner.

## Result on chr20:48,176,830-48,229,446 (52.6 kb): CLOSED

The deficit gap hiphase spans at 100.0% over 252 reads.

| stage | result |
|---|---|
| gauge | whole-chr20 pipeline run: 238 blocks, 1,319 tagged reads in the gauge window |
| flanks | left `48162480` (2 sites, **0 tagged reads**), right `48229446` (916 sites, 859 reads) |
| gap arm | machinery on `48,171,830-48,234,446` -> 2 gap-local blocks; **11/11** in-gap het sites phased |
| compose | 42 allele voters, **[15, 0, 0, 27]** -- unanimous, 58 reads cross the seam -> flip 0 |
| link left | shared sites **1/1 agree**; alleles **14 of 66** crossing reads -> 2 sources agreeing, flip 0 |
| link right | shared sites **4/4 agree**; tags **71 reads**; flank tags x gap alleles **71/298** -> 3 sources agreeing, flip 0 |

The two- and three-source rows above are from `gap_lab.py` **after** the
shared-site vote was added. An earlier background run of the same window,
dispatched before that code existed, printed a single `alleles` source on the
left (`n=14 votes=[9, 0, 0, 5] [alleles (66 cross, 2+14 sites)] -> link,
flip=0`); it is the same verdict reached on one source rather than two, and it
should not be cited as multi-source agreement.
| link validation | frame hap1 carries PATERNAL over 14 sites; left flank PATERNAL over 2 sites (**CORRECT**), right flank PATERNAL over 81 sites (**CORRECT**) |
| gate | tagged 1,319 -> 1,426, concordant 1,316 -> **1,423**, accuracy 99.77% -> 99.79%, **0** concordant->discordant, **107 newly tagged, all concordant**, 0 lost |
| verdict | **PASS -- gap CLOSED, both links validated, +107 concordant reads** |

The earlier conclusion that this gap cannot be closed was an artifact of the
first version's window: the two hard linkage breaks
(`48,096,582->48,123,657` and `48,123,657->48,147,230`, zero spanning reads at
72x) lie **outside** the gap, and only became obstacles because the 150 kb window
pulled the left flank to the far side of them.

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
