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
