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
