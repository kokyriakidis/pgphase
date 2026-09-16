# Why two blocks 220 bp apart never join: the site and read phase-set labels disagree

On `chr20:48,176,830-48,229,446` with `--retry-unphased-with-bam`, the retry now
recovers the bridging site and the left block reaches `48,229,226` -- **220 bp**
from the right block's first site at `48,229,446`. The two never join, and the
reason is not linkage.

## The evidence is complete

| | |
|---|---|
| reads spanning both boundary sites | **60**, none untagged |
| their allele co-occurrence | 22 `(0,0)`, 34 `(1,1)`, 1 `(0,1)`, 3 `(1,0)` -> **56 in phase, 4 against** |
| `48,229,446` (CLEAN_HET_SNP) vs read truth | **60/60, 1.000** |
| `48,229,227` (NOISY_CAND_HET) vs read truth | **56/60, 0.933** |
| orientation of both sites in our own VCF | identical: `1\|0`, `HAP_ALT=1 HAP_REF=2` |

The four exceptions are not a phase inconsistency: they equal the deletion site's
own genotyping error, which truth puts at the same 56/60. So the shared reads
give the relative orientation directly, and it is the orientation truth prefers
at both sites.

## The bug: two labellings that can disagree, and nothing reconciles them

The pipeline carries a phase set per **site** (`CandidateVariant::phase_set`, the
`PS` in the emitted record) and a phase set per **read** (the `PS` tag). At this
boundary they contradict each other:

| site | its `PS` | reads tagged LEFT | reads tagged RIGHT | untagged |
|---|---|---:|---:|---:|
| 48,183,976 | 48147225 | 75 | 0 | 0 |
| 48,202,056 | 48147225 | 7 | 0 | 62 |
| 48,204,383 | 48147225 | 1 | 1 | 69 |
| **48,225,787** | 48147225 | **0** | **51** | 24 |
| **48,229,226** | 48147225 | **0** | **60** | 0 |
| 48,229,446 | 48229446 | 0 | 60 | 0 |

The left block's last two sites are labelled with the left phase set while
**every read covering them is tagged into the right one**. So the 82 kb "left
block" has a right end populated entirely by right-block reads.

That is also exactly why recovery cannot close it. The tier report for the
boundary reads:

```
gap 48229226-48229446  tier=1..3 pass=0  L=0 R=1  leftPS=-1  rightPS=48243938  status=partial
gap 48229226-48229446  tier=1..3 pass=1  L=0 R=0  leftPS=-1  rightPS=-1        status=open
```

`L=0, leftPS=-1`: the left flank has no reads holding its own phase set at the
boundary, so the flank vote has no voters and every tier abstains. The split is a
labelling artifact that then starves the mechanism meant to repair it.

## Why the competitors phase it

A caller like HiPhase assigns one phase set per connected component of its
variant graph: two heterozygous sites joined by shared reads are in the same
block by construction, and there is no second, independent per-read phase-set
label that can contradict the per-site one. With 60 reads spanning 220 bp and
both sites heterozygous, these two sites are trivially one component. HiPhase
spans the whole region in a single block and, measured against the same read
truth, its span here is correct.

Our pipeline maintains both labellings and reconciles them only in two places --
`flip_chunk_hap` at a **chunk seam**, and the two flank votes inside gap
recovery. Neither applies to two blocks inside one chunk, so they can diverge
and stay diverged.

## The fix this points at

Reconcile the labels rather than inventing a new link: when every read covering a
site is unanimously tagged into a different phase set, the site's phase set must
follow its reads -- equivalently, merge the two phase sets, taking the
orientation from the shared reads exactly as `select_stitch_orientation` already
does at a chunk seam. Here that vote is 56 against 4, and the merged orientation
is the one truth prefers at both sites, so the merge is safe on this window.

Worth noting separately: the tier report's right flank is `rightPS=48243938`,
not the `48229446` the emitted VCF carries for that block -- a third phase-set
identity for the same block, and a sign the reconciliation problem is not
confined to this boundary.

## Why the existing stitching machinery does not fire, traced to the line

The machinery is the right one and it is already wired: the retry re-phases the
window, and gap recovery's flank vote is exactly the "stitch to the previous
window" step. That vote counts **reads holding each side's phase set** -- and
`update_read_phase_set` (`collect_phase.cpp`) refuses to let a homopolymer indel
grant a read its phase set:

```cpp
// Use the same eligible evidence as init_assign_read_hap. An
// excluded repeat can inherit a preceding PS without a link.
if (var.is_homopolymer_indel || var.lcd_var_i_to_cate == kCandNoisyCandHom ||
    (!var.msa_insertion_alts.empty() && !var.gap_link_supported)) continue;
```

Both of the left block's terminal sites sit in long A homopolymers
(`48,225,780: ctctctcaaaaaaaaaaaaaaaa`, `48,229,220: tctttagaaaaaaaaaaaaaaac`), so
every read over them skips them and takes the phase set of the next eligible het
-- the right block. Hence `L=0, leftPS=-1` and every tier abstains. The stitch
has no voters because its input was excluded upstream, not because it is missing.

Three things compound, and a probe confirms each:

1. The exception for exactly this case exists -- `CandidateVariant::hp_gap_scorable`
   -- and `init_assign_read_hap` honours it (`collect_phase.cpp:362`), so such a
   site can assign a read's **haplotype**. `update_read_phase_set` omits it,
   though its comment claims the same eligible evidence, so the same site can
   never grant a **phase set**.
2. The flag is set only for the homopolymer recovery tier's window
   (`in_hp_gap`), never for a retry window.
3. `select_gap_link_sites`, which sets it, is itself guarded on
   `opts.gap_hp_link_beg >= 0 || (opts.recover_gaps && opts.private_msa_admit_all_in_region)`
   -- neither holds in a normal run, so the function never executes. A probe at
   the site confirms the rest is eligible: `hp_indel=1 hp_scorable=0 cate=0x100
   msa_verified=1 cons=1/0`.

## Enabling it joins the window -- with the wrong orientation

Fixing all three (parity in `update_read_phase_set`, the flag for retry windows,
and running the selection when a retry window exists) produces exactly the
intended behaviour:

| | before | after |
|---|---|---|
| blocks over the region | 2 (82.0 kb + 50.0 kb) | **1 (48,147,225-48,279,445, 132.2 kb, 61 sites)** |
| spans the region | no | **YES** |
| reads at `48,229,227` | `PS=48229446` (60) | **`PS=48147225` (60)** |
| tier 1 pass 0 | `L=0 R=1 leftPS=-1 status=partial` | **`L=1 R=1 leftPS=rightPS=48243938 status=joined`** |

And it is **wrong**. The read-level gate:

```
tagged 393 -> 433   concordant 393 -> 232   conc->disc 161   lost 3   new 3c/40d
```

**161 concordant reads become discordant**, and window accuracy falls from
**100.00%** (167/167) to **71.50%** (148/207). Reverted; the reverted build
reproduces two blocks and 167/167.

So the missing piece is not only that the site cannot carry a phase set. The
orientation the joining path then chooses disagrees with the reads: the shared-read
vote at the boundary is **56 in phase against 4** (and those 4 equal the deletion
site's own error rate), yet the merge inverts one side. The tier report naming
`leftPS = rightPS = 48243938` -- a third phase-set identity for these blocks --
is where to look next: the orientation is taken from that path, not from the 60
reads that span both sites.

## Why it mis-joined: a label merge with no orientation flip

Diagnosed with the built-in diagnostics rather than new probes. `--verbose 2`
prints the flank votes:

```
GapLinkVotes  48229226 48229446  proposal_ps=48243938  side=0 (left)   0,4,1,1    -> straight 1, flipped 5
GapLinkVotes  48229226 48229446  proposal_ps=48243938  side=1 (right)  0,42,39,0  -> straight 0, flipped 81
```

The left flank is oriented by **6 voters, 5 against 1**, the right by 81
unanimous -- but that is not the cause. Raising `--min-block-link-reads` to 6, 8
or 12 still reports tier 1 `joined` and leaves read consistency at 54-56%, so the
flank vote is not what corrupts the window.

The read tags say what does:

```
390 reads tagged in both arms:   HP changed for 1,   PS changed for 230
left  population (161 reads):  PS 48147225 -> 48147225   HP changed for 0
right population (173 reads):  PS 48229446 -> 48147225   HP changed for 0
after the merge:  left  hap1 = PATERNAL (160 against 1)
                  right hap1 = MATERNAL (173 against 0)
```

The right population is **relabelled into the left block's phase set with no flip
applied**. Both populations keep the haplotype labels they were given under their
own orientation, so the merged block carries two contradictory hap-to-parent
mappings: 233 reads right, 160 wrong, which is the 53.6% measured (before the
change each block was 100.0% on its own).

That is exactly the hazard the exclusion's own comment names -- *"An excluded
repeat can inherit a preceding PS without a link."* Letting a homopolymer indel
grant a phase set lets a read inherit the neighbouring block's label **without
any evidence that the two sides share an orientation**, and nothing downstream
supplies the missing flip.

### What the fix has to be

Two steps, and only the second was attempted:

1. **Compute the relative orientation from the shared reads and apply it.** At a
   chunk seam this is `select_stitch_orientation` followed by
   `apply_chunk_flip_and_merge`, which rewrites the downstream phase set to the
   upstream id *and flips the hap labels as it does so*. Within a chunk nothing
   performs this step.
2. Relabel the phase sets as one block.

The evidence for step 1 is present and unambiguous -- 60 reads span the two
boundary sites, 56 in phase against 4 -- so a within-chunk merge that reuses
`select_stitch_orientation` on those reads and then `apply_chunk_flip_and_merge`
is the shape this window needs. Relabelling without it produces a spanning block
at chance accuracy.

## Root cause: the left block already contains an internal switch

Two hypotheses were tested and both are **wrong**, so they are recorded as such:

- *The site chain crossed a spacing with no spanning reads.* It did not. Every
  consecutive pair in the left block is read-supported: 18.1 kb carries 7 reads,
  2.3 kb carries 58, 21.4 kb carries 7, 3.4 kb carries 51, and the 220 bp
  boundary carries 60. The earlier "41.8 kb with zero spanning reads" was the
  **pre-fix** site set, before `48,204,383` was recovered.
- *Site and read labels disagree within a block.* They agree: HP against the
  allele a read carries is 100.0% at `48,162,480`, `48,183,977`, `48,229,446` and
  `48,232,579`, and 93.3% at `48,229,227` (its own error rate).

Per-site orientation against read truth gives the answer:

| site | our GT | consistency | hap1 carries |
|---|---|---:|---|
| 48,162,480 | `1\|0` | 1.000 | PATERNAL |
| 48,176,830 | `1\|0` | 0.986 | PATERNAL |
| 48,183,976 | `0\|1` | 1.000 | PATERNAL |
| 48,177,788 | `0\|1` | **0.500** | no information |
| 48,202,056 | `0\|1` | **0.507** | no information |
| 48,204,383 | `0\|1` | 0.958 | PATERNAL |
| **48,229,226** | `1\|0` | 0.933 | **MATERNAL** |
| 48,229,446 | `1\|0` | 1.000 | MATERNAL |
| 48,230,918 | `1\|0` | 0.983 | MATERNAL |

**The left block's body is PATERNAL-on-hap1 and its terminal site is
MATERNAL-on-hap1** -- in phase with the right block, not with its own body. That
switch is in the shipped build, before any merge is attempted.

Everything else follows from it:

- Per-block read consistency reads 100.0% because the mis-oriented terminal
  sites have **no reads of their own** -- every read covering them belongs to the
  right block -- so the error is invisible in read space.
- The boundary vote reports "straight" (`n11=34 n12=3 n21=1 n22=22`, the same
  56-against-4 measured by hand) because it compares the *mis-oriented* site with
  the right block, i.e. the right block against itself.
- Merging therefore inverts the left block's body: 233 reads right, 160 wrong,
  the 59.29% measured, where each block alone had been 100.0%.
- `verify_retry.py` reports `switches: 0` because its detector needs a run of
  scorable sites on each side, and `48,202,056` is unscorable while `48,225,787`
  is absent from the VCF. This is exactly the `gate_blind: True` it flagged, and
  the reason that flag exists.

Note where the switch sits: between `48,204,383` and `48,229,226`, the only
intervening evidence is `48,202,056` and `48,177,788`, both carrying **no
haplotype information** (0.507 and 0.500). So the single open finding from the
injection audit -- a site used as a phased het at chance segregation -- is not a
side note. It is in the chain exactly where the orientation flips.

### Where the fix belongs, and what is not yet known

The parity that orients a site into a block is applied in the emit loop of the
phase-set assignment (`collect_phase.cpp`, `parity[hi]` swapping
`hap_to_cons_alle[1]/[2]`), and the union-find that can set it merges components
on a bare majority, `flip[child] = ... ^ (edge.conflict > edge.agree)`, with no
minimum count or margin. But that path logs `GapCleanBlockLink` under
`--verbose 2` and logged **zero edges** for this window, so it is not what
oriented `48,229,226` here. The path that did is not yet identified and is the
next thing to pin -- not to guess at.

A merge is not safe on this window until that is fixed: the orientation the
merge would need is the one the block body carries, and the block's own terminal
site contradicts it.
