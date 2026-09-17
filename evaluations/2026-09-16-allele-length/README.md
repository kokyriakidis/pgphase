# Why the remaining panel windows do not span: exact-allele matching in repeats

With `--keep-noisy-kmeans` and the MSA depth fix, two of the six panel windows
span. This is why the other four do not, worked on `chr20:24,105,188-24,142,287`,
the window that gained **nothing** in that arm (3 in-gap hets, 0 reads tagged) --
the cleanest signal in the panel.

## It is not a site deficit, and not coverage

We phase the same three sites hiphase does inside the gap:

| position | hiphase | us |
|---|---|---|
| 24,105,188 `A>G` | `0\|1` | `0\|1` |
| 24,121,713 | `CTTTT,CTTTTT` -> **`3\|1`** | `C>CTTTT` -> `0\|1` |
| 24,142,287 `A>G` | `1\|0` | `0\|1` |

Mean depth across the gap is **70.9** at MAPQ >= 1 and **70.3** at MAPQ >= 30, so
neither coverage nor mapping quality is involved. The chain is:

```
24,103,779   agree=0   conflict=0
24,105,188   agree=53  conflict=0
24,121,714   agree=6   conflict=0     <- 6 usable reads at 70x
24,142,287   agree=0   conflict=0     <- the break, 20.6 kb
24,142,446   agree=76  conflict=0
```

## The measurement: the haplotypes differ by net length, and we match exactly

Net length change over a +/- 25 bp window at `24,121,713`, across the 74 reads
that fully cover it, split by parental truth:

| net | reads | MAT | PAT |
|---:|---:|---:|---:|
| +8 | 16 | 16 | 0 |
| +7 | 10 | 10 | 0 |
| +6 | 8 | 8 | 0 |
| +9 | 4 | 4 | 0 |
| +10 | 1 | 1 | 0 |
| **+5** | 7 | 2 | 5 |
| +4 | 13 | 0 | 13 |
| +3 | 11 | 0 | 11 |
| +2 | 3 | 0 | 3 |
| +1 | 1 | 0 | 1 |
| **0 (reference)** | **0** | 0 | 0 |

Two things follow. **No read carries the reference allele** -- both haplotypes
have an insertion here, so the locus is not ref-versus-alt at all, and a
biallelic record must misrepresent it. And the haplotypes **separate by net
length**: maternal reads sit at +6 and above, paternal at +4 and below, with only
the seven reads at +5 mixed. A net-length assignment separates 67 of 74 reads.

Our record is `C>CTTTT`, exactly +4. Matching that allele exactly admits **13 of
74** reads, and the chain scored **6**. So the site is not information-poor: the
scoring discards 82% of the reads that do carry the signal, and the 20.6 kb link
is left with no evidence at all.

## Same defect, third window

This is the signature already measured at `48,225,786`, where the seven reads
spanning the failing link separated perfectly by net length -- paternal -7, -6,
-6, -5 against maternal -1, 0, 0 -- while our records there were a 1 bp and a
4 bp deletion, so no paternal read matched either exactly and the link scored one
usable read of seven. With 74 reads instead of 7, the mechanism is now
unambiguous.

It also explains the records carrying `AF = 1.000` with **zero reference reads**
(`39,846,791`, two such records with opposite genotypes; `48,225,787`): each
biallelic record only ever sees the reads whose indel matches its own allele, so
its allele fraction is computed over a subset that by construction contains no
reference read.

One root cause therefore accounts for the starved link votes at full coverage
and for the impossible allele fractions, and for the unspanned windows measured
here.

**Not for all of them, though.** An earlier version of this sentence claimed it
accounts for *every* remaining unspanned window. The measurement behind it
establishes the net-length mechanism on the windows examined in this document;
`55,843,827` was separately diagnosed with a different defect -- of the eight
sites the retry admits there, all carry `msa_verified = 1` and only one
segregates at or above 0.90 against read truth, which is a site-quality problem
rather than an allele-length-matching one. Two mechanisms, not one, and the
allele-length fix should not be expected to close that window.

## The fix, and why it is not bolted on here

Assign a read to the nearer haplotype by **net length across the variant's repeat
tract**, rather than requiring an exact match at the variant key, in
`get_var_allele_i_from_cons_aln_str`. Emitting the locus as multiallelic, as
hiphase does, is the complementary half on the output side.

This changes how every indel site is called, not only the ones in gaps, so it
gets its own pass with the panel and the chromosome-wide hamming baseline as
gates -- the same discipline that caught the reclassification attempt being
destructive. The two candidate ambiguities to settle there: which reads the +5
class should go to when a tract's two haplotypes are adjacent in length, and
whether the tract bounds come from the variant's own repeat annotation or from a
fixed flank.
