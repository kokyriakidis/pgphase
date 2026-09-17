# One panel window end to end: chr20:24,105,188-24,142,287

37.1 kb. Hiphase spans it at 100.00% over 192 scored reads. It is the window that
gains **nothing** in every arm tried, which makes it the cleanest signal in the
panel.

## What is inside the gap

19 candidates, identical in the hybrid and the alignment channel:

| class | n | phase information |
|---|---:|---|
| `CLEAN_HOM` (0 reference / all alt) | 15 | none by construction |
| `CLEAN_HET_SNP` | 2 | the gap's own left and right boundaries |
| `CLEAN_HOM` | -- | -- |
| `NOISY_CAND_HET` | 2 | the only interior evidence |

Provenance: **13 of the 19 are catalog-matched, 6 are alignment-only**. Of the
four heterozygotes, **three are alignment-only and one -- the right boundary
`24,142,287` -- is catalog-matched**.

Hiphase's phased heterozygotes inside the gap are exactly three:
`24,105,188`, `24,121,713` (`C>CTTTT,CTTTTTTT,CTTTTTTTT` at `3|1`) and
`24,142,287`. So it has no more interior sites than we do -- one.

## Is the interior site injected correctly, with the right representation?

Yes, and this is the multiallelic work paying off. `24,121,714` is now a single
record carrying both alleles, `ALT = TTTT,TTTTTTTT`, DP 61 -- the same locus that
a few commits ago was two biallelic records admitting 13 of 74 reads.

Scored against read truth, it is the best site in the window:

| allele | reads called | parental split | purity |
|---|---:|---|---:|
| +4 (`TTTT`) | 13 | **13 paternal, 0 maternal** | **1.000** |
| +8 (`TTTTTTTT`) | 16 | **16 maternal, 0 paternal** | **1.000** |

29 reads carry a callable allele, and the two alleles separate the haplotypes
perfectly. Hiphase's `3|1` picks the same two lengths. So the representation, the
depths and the per-read assignments are all correct here.

## Why the window still does not phase

Two causes, and only the first is ours to remove.

**1. The hybrid never lets the interior sites into the solve.** Both
`NOISY_CAND_HET` records sit at `PHASE_SET = 0`:

| site | alignment channel | hybrid |
|---|---|---|
| `24,121,714` (multiallelic, informative) | `PS = 24,103,779` | **`PS = 0`** |
| `24,131,708` | `PS = 24,103,779` | **`PS = 0`** |

`hybrid_collect.cpp` sets `skip_noisy_kmeans` unconditionally, so the hybrid's
left block ends at the gap's left boundary (`24,103,779-24,105,188`, 2 sites)
while the alignment channel's reaches `24,121,713` (3 sites). The chain then
steps 37.1 kb from boundary to boundary and reports `agree = 0, conflict = 0`.

**2. With them in, the right half still cannot be crossed.** The alignment
channel chains the left half correctly and then stops:

| link | spacing | reads spanning | reads with a callable allele at **both** | verdict |
|---|---:|---:|---:|---|
| `24,105,188` -> `24,121,714` | 16.5 kb | 11 | **6** | determined, and correct (4 reads paternal/allele 1, 2 maternal/allele 2) |
| `24,121,714` -> `24,142,287` | 20.6 kb | **2** | **0** | no link |
| `24,121,714` -> `24,131,708` | 10.0 kb | 29 | 8 | link exists |
| `24,131,708` -> `24,142,287` | 10.6 kb | 29 | 24 | link exists but carries no orientation |

The zero is not an artifact of a strict allele call -- it holds at every flank
width from 25 bp down to 0, because only **two reads span the 20.6 kb at all**.

And the site that *is* linked to the right boundary by 24 reads is a phantom.
`24,131,708` is a 1 bp deletion at 39/29, allele fraction 0.427 -- het depths by
every threshold -- but its net-length distribution across the tract does not
separate the haplotypes:

```
  +0   17 maternal / 15 paternal
  -1   14 maternal /  8 paternal    <- the emitted allele, purity 0.636
  -2    6 maternal /  2 paternal
  +1    4 maternal /  2 paternal
```

Its pairing with the right boundary is 6 / 3 / 7 / 8 across the four allele
combinations: no orientation. It is never emitted in either arm
(`HAP_ALT = HAP_REF = 0`), which is the correct outcome.

## What this window says

Injection, representation, depth and per-read assignment are all correct at the
one interior site that matters. The window is unphased because of the hybrid's
`skip_noisy_kmeans` exclusion, and beneath that because its informative pair is
20.6 kb apart with two spanning reads and no read yielding an allele at both,
while the only bridging site is uninformative.

That makes it the wrong window to chase for a phasing fix: after removing the
exclusion there is no read-based evidence left to use. Hiphase spans it from the
same three sites, so its join across an interval two reads span is not something
our read-level machinery can reproduce by using more sites -- which is consistent
with the site-level check that its five spans of this kind were all correct.
