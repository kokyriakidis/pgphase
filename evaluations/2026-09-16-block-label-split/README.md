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
