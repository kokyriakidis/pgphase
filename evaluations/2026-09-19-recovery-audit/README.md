# What the recovery finds, and where it is lost

A structure that records every candidate the recovery sub-solve finds inside a
recovery window, and what the merge did with it: `RecoveredCandidate`
(`collect_pipeline.hpp`), written by `--recovery-audit-out FILE`.

## The window

`chr20:55,313,902-55,358,363`. Hiphase spans it with ONE block at **100.0% over
66 reads**; we tag 46 reads in one block and place **20 of them wrong**. Both
parents are present, so unlike chr20:45,391,776 this window is scorable and the
competitor proves the signal is there.

The signal is one site: hiphase's `55,336,493 GTGTG>G`, **separation 0.971**
against read truth (alt 0 MAT / 15 PAT, ref 19 MAT / 1 PAT). Our nine records
in the span top out at 0.638.

## What each arm holds there

| arm | at 55,336,46x |
|---|---|
| alignment | `INS C>GT` DP 48 (32/16) **and** `DEL GTGT>.` DP 29 (3/26), both `NOISY_CAND_HET`, both phased into PS 55,331,014 |
| graph | `INS C>GT` DP 47 only, `REP_HET_INDEL`, PS -1 |

The deletion is MSA-reconstructed, which a graph chunk cannot do (no digars),
and the insertion is demoted by the reference-context rule. So the graph arm
has neither.

## Where the recovery loses them

With `--graph-noisy-msa` the sub-solve re-reads the BAM and finds both. The
audit over this 90 kb region:

| | count |
|---|---:|
| candidates found in recovery windows | 168 |
| already known to the parent | 39 |
| rejected as outside a window | 74 |
| appended to the chunk | 90 |
| **with usable metadata** | **90** |

Both target candidates are appended with correct anchored forms (`META_REF=C`
for the insertion, `CGTGT` for the deletion) and correct phase sets, and all 90
appended candidates round-trip through `vcf_to_variant_key` (90/90 ok). So
neither translation nor injection is where they die.

They die in the parent chunk's re-solve. At the writer:

```
WR2 idx=131 pos=55336460 type=1 cate=6 h=[0,0] ref_cov=32 alt_cov=17
WR2 idx=132 pos=55336460 type=2 cate=6 h=[1,1] ref_cov=3  alt_cov=26
```

Each haplotype takes its majority allele independently, both majorities land on
the same side, and the writer skips a site whose two consensus alleles are
equal (`graph_collect.cpp:218`). A 32/17 heterozygote is emitted as nothing.

This is precisely the collapse `allele_depths_call_het` guards against, and the
guard is scoped by `retry_windows` -- which the sub-solve now has and the
parent did not.

## Two fixes measured, one kept

**Carrying the recovery windows into the parent re-solve** emits the insertion
(records 54 -> 56 on the window). Chromosome-wide on the default path it is a
regression: read hamming 1.153% -> 2.067% (corrected for unscorable blocks,
0.967% -> 1.883%) with blocks 339 -> 498. Kept behind `--graph-noisy-msa` only;
the default path is byte-identical over whole chr20.

**Pinning the sub-solve's consensus** for injected sites -- carrying the het
call rather than re-deriving it -- is worse: emitted records on the window fall
54 -> 46 and the window collapses to a single site. The reason is a gauge
mismatch: `hap_to_cons_alle` is expressed in the SUB-SOLVE's haplotype labels,
while the parent resets reads and re-labels them, so a pinned consensus asserts
an arbitrary orientation as fact. Allele depths and per-read alleles are
orientation-free and are carried; the consensus cannot be.

What remains for the deletion is representation, not admission: our biallelic
record splits 3 ref / 26 alt because most reads in a GT tract carry some other
length, so every depth-based het test rejects it, while hiphase's record at the
same locus splits 15/19 and is informative.

Unit 3/3, predicate 151/151, window 125/125.

## Do we need a redesign? One specific part, and the evidence says which

The BAM sites ARE usable in the graph gap region, and nothing rejects them:
discovery inside a recovery window is BAM-only by construction, 90 of 90
appended candidates carry usable metadata, and all 90 round-trip. The problem
is downstream of injection, and every patch applied to it has been measured:

| what was tried | default path | flagged arm |
|---|---|---|
| baseline | **1.153%** (corrected 0.967%), 339 blocks | -- |
| carry recovery windows into the parent re-solve | 2.067% (1.883%), 498 blocks | 4.205% |
| key the escape on `bam_injected` provenance instead | **byte-identical** | 4.843%, 501 blocks |
| pin the sub-solve's `hap_to_cons_alle` | -- | window records 56 -> 46 |
| admit repeat indels by link purity (earlier) | 1.435% / 1.886% / 2.038% | -- |
| stage 2 with a chunk-wide round 2 (earlier) | -- | 4.197% |

Six different admission mechanisms, one shared failure: each lets the injected
sites participate in a chunk-wide re-solve whose read labels are re-derived
from scratch, and the solve gets worse rather than better.

Provenance keying is kept because it is the narrowest key and leaves the
default byte-identical; the window list is removed.

### What the evidence actually points at

The merge keeps the sub-solve's candidates and per-read alleles and THROWS AWAY
its phase sets, then asks the parent to re-derive phase for those sites. That
is the one structural choice all six failures share:

- Pinning failed for a specific, diagnosable reason -- `hap_to_cons_alle` is
  expressed in the SUB-SOLVE's haplotype labels while the parent resets reads
  to its own gauge -- and nothing in the merge supplies the map between the two
  gauges.
- That map is exactly what the cross-chunk stitch already computes:
  `select_stitch_orientation` votes with reads two blocks share and flips one.

So the redesign that follows from the measurements is bounded: treat a
recovery sub-solve's result as a BLOCK TO STITCH rather than as sites to
re-solve. Keep its phase sets, orient them against the parent block by the
shared-read vote, relabel, and do not re-run the parent k-means over them. That
supplies the gauge pinning lacked, and it stops the injected sites perturbing
every read's label, which is the mechanism behind all six regressions above.

Not implemented. It should be measured on chr20:55,313,902-55,358,363, where
hiphase spans one block at 100.0% over 66 reads and we place 20 of 46 wrong,
against the default's 1.153% / 0.967% corrected.

## The stitch design, implemented and measured

`--stitch-recovered` (off by default) does what the previous section proposed:

1. **Orient** each recovery sub-solve against the parent on the reads they
   share -- the vote `select_stitch_orientation` runs between chunks. On
   `chr20:55,290,000-55,380,000` it is decisive: **same 157, cross 14, no flip**
   over 171 shared reads, so the sub-solve's gauge already matched the parent's.
2. **Import as a block**: the parent keeps the read labels it already has and
   only incorporates the imported sites (`anchored`), instead of re-solving
   from scratch.

### Whole chr20

| arm | tagged | read blk | disc | hamming | corrected | VCF blk |
|---|---:|---:|---:|---:|---:|---:|
| default | 219,055 | 326 | 2,526 | **1.153%** | **0.967%** | 339 |
| `--graph-noisy-msa` | 225,718 | 369 | 10,931 | 4.843% | 4.691% | 501 |
| `--graph-noisy-msa --stitch-recovered` | 225,645 | 336 | 9,451 | **4.188%** | 4.043% | **350** |

The stitch import is the best of every recovery-import variant tried: it cuts
the flagged arm's misplaced reads by 1,480 and takes VCF blocks from 501 back
to 350, almost the default's 339 -- so importing rather than re-solving does
repair the fragmentation the re-solve caused. It does not make the arm usable:
4.188% against the default's 1.153%. The flag stays off and the default path is
byte-identical over whole chr20.

### What the measurements refuted

- **Gauge mismatch was NOT the reason pinning failed.** The vote says the two
  gauges already agree (157 vs 14, no flip), and carrying the sub-solve's
  per-site consensus still cost 8 emitted records on the window (54 -> 46). A
  consensus is defined against the read set it was derived from; as a fixed
  constraint in the parent it contradicts reads the sub-solve never saw. Only
  the orientation is applied now, never the consensus.
- **Anchored incorporation alone is inert** on the window (54 records, same
  block structure as the flag alone); its value is chromosome-wide, in the
  block count.

The target window `chr20:55,313,902-55,358,363` is still not closed: the
insertion at 55,336,460 is not emitted by any arm that does not also carry the
recovery windows into the parent re-solve, and that costs 1.153% -> 2.067% on
the default path. What blocks it is representation, not import -- our biallelic
deletion record splits 3 ref / 26 alt in a GT tract where hiphase's record at
the same locus splits 15/19.

## Is the representation wrong in the BAM pipeline, or only on import?

**Only on import, and not in the import itself.** Traced end to end at
`chr20:55,336,460`.

**The alignment pipeline's representation is correct and complete.** It emits
the locus as a complementary pair:

```
55336460  C     > CGT    GT=0|1  PS=55331014  AD=32,16
55336460  CGTGT > C      GT=1|0  PS=55331014  AD=3,26
```

hap2 carries the 2 bp insertion, hap1 the 4 bp deletion -- the right
description of a GT tract where both haplotypes differ from the reference.

**The import preserves it.** Probing the parent chunk immediately after the
merge:

```
PROF pos=55336460 type=1 inj=1 cons=[0,1]
PROF pos=55336460 type=2 inj=1 cons=[1,0]
```

Exactly the alignment's pair, in the parent's own candidate table, with
metadata and phase sets. So neither translation, nor anchoring, nor injection
damages the record.

**Two steps AFTER the merge destroy it.**

1. The second k-means round (`kCandGermlineVarCate`, unanchored) recomputes the
   imported sites -- they are `NoisyCandHet`, so this round owns them -- and
   each haplotype takes its majority allele independently: `(0,1)` and `(1,0)`
   become `(0,0)` and `(1,1)`. The writer skips a merged site whose two
   consensus alleles are equal, so the whole locus disappears.
2. `drop_superseded_colocated_records` demotes the deletion as a duplicate
   description of a colocated merged record.

**Anchoring that second round fixes the window and does not generalise.** With
it, the insertion is emitted, joins the existing block `PS 55,331,014`, and the
target window becomes ONE 9-site block instead of two. Chromosome-wide it costs
read hamming 4.188% -> 6.160% and blocks 350 -> 568, so it is not kept; the
round stays unanchored and the loss is recorded rather than traded away.

Exempting an injected record with a decided heterozygous consensus from the
colocated-record demotion is kept (it is the correct rule -- such a record
describes the other haplotype, not a duplicate call) but is inert here: the
deletion is already gone by then.

| arm | tagged | read blk | disc | hamming | corrected | VCF blk |
|---|---:|---:|---:|---:|---:|---:|
| default | 219,055 | 326 | 2,526 | **1.153%** | **0.967%** | 339 |
| `--graph-noisy-msa` | 225,718 | 369 | 10,931 | 4.843% | 4.691% | 501 |
| `+ --stitch-recovered` | 225,645 | 336 | 9,451 | **4.188%** | 4.043% | **350** |
| `+ anchored second round` | 225,789 | 450 | 13,909 | 6.160% | 5.888% | 568 |

Default byte-identical over whole chr20. Unit 3/3, predicate 151/151,
window 125/125.
