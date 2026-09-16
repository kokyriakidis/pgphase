# The deficit on stock defaults, and one window worked end to end

Every earlier deficit set was derived from a non-stock configuration. This one is
derived from the pipeline as it ships: `collect-hybrid-variation` with no extra
flags, whole chr20.

| | gaps | span |
|---|---:|---:|
| our gaps | 963 | -- |
| a competitor spans it | 744 (hiphase), 685 (longphase), 327 (whatshap) | -- |
| **a competitor spans it at >= 99% over >= 100 scored reads** | **104** | **2.62 Mb** |

`find_current_deficit.py` and `score_deficit.py` (from the sibling directory) do
the derivation and score every competitor span against the diplinator read truth,
so a window only qualifies when the competitor is demonstrably right. Several
competitor spans are worse than chance -- hiphase is under 90% in 459 of the
gaps it spans -- and those are not opportunities.

## Target: chr20:55,843,827-55,889,113

45.3 kb, 252 reads, hiphase spans it at **100.00%** over 220 scored reads. Our
stock run leaves it open as two blocks that end and begin exactly on the gap
bounds (`55,824,242-55,843,827`, 6 sites; `55,889,113-55,939,036`, 71 sites).

**We already hold every site hiphase crosses on.** Its five in-gap hets are
`55,843,827 G>A`, `55,862,239 TAATGGC>T`, `55,862,269 T>CAGTAAATTAATT`,
`55,883,019 AATATAT>A,AAT` and `55,889,113 A>T`; our candidate table has all of
them with textbook heterozygous depths -- `55,862,240 AATGGC>.` at 31/29
(AF 0.483), `55,862,270 T>CAGTAAATTAATTATC` at 31/29, `55,883,020 ATAT>.` at
17/18 -- each classified `NOISY_CAND_HET` and each left at `PHASE_SET 0`.

The reason is one line: `hybrid_collect.cpp:128` sets
`opts.skip_noisy_kmeans = true` unconditionally, so that class never enters
phasing in the hybrid subcommand.

## Bug found and fixed: `--retry-unphased-with-bam` was a no-op

The scoped mechanism for exactly this case did nothing. It was gated on
`!discovery_flags.empty()`, and `discovery_flags` is only populated under
`recover_gaps`; worse, the re-solve was gated on `readmitted > 0`, which counts
*restored* categories. Without recovery nothing had been zeroed, so `readmitted`
stayed 0 and the re-solve -- the part that admits the noisy class -- never ran.
The option therefore appeared to work and changed nothing.

Fixed: the retry runs on its own, the category restore stays recovery-specific,
and the re-solve triggers whenever a window failed. The detector then finds
exactly the right interval, `55,843,827-55,889,114 (45.3 kb)`, holding 614
candidates of which 600 are catalog sites.

| the window | tagged | blocks | in-gap hets | concordance |
|---|---:|---:|---:|---|
| hiphase (target) | 650 | 2 | 5 | 99.85% |
| stock | 413 | 2 | 2 | 99.76% |
| stock + retry | **501** | 2 | **10** | **97.80%** |

The gap shrinks from 45.3 kb to 10.8 kb and 79 newly tagged reads are correct,
against 9 wrong and one previously-correct read flipped.

## Why it cannot be the default yet: MSA verification does not discriminate

Chromosome-wide the trade is bad:

| whole chr20 | tagged | blocks | read hamming |
|---|---:|---:|---:|
| stock | 212,320 | 452 | **0.559%** |
| stock + retry | 218,989 | 521 | **2.723%** |

+6,669 reads for roughly five times the error. The cause is the admission gate,
and this window measures it for the first time. Every one of the eight interior
sites the retry admits carries `msa_verified = 1` (category `0x100`), yet scored
against read truth:

| site | msa_verified | homopolymer | segregation |
|---|---:|---:|---:|
| 55,861,267 `CT>C` | 1 | 1 | 0.825 |
| 55,862,239 `TAATGGC>T` | 1 | 0 | **0.689** |
| 55,862,269 `T>TCAGTAAATT` | 1 | 0 | **1.000** |
| 55,867,093 `AT>A` | 1 | 1 | **0.516** |
| 55,871,836 `C>CA` | 1 | 0 | 0.721 |
| 55,882,616 `CT>C` | 1 | 1 | **0.559** |
| 55,883,019 `AATAT>A` (two records) | 1 | 0 | **0.507** |
| 55,883,023 `TAT>T` | 1 | 0 | **0.507** |

Two of eight are informative and four are phantoms. So `msa_verified` passes
phantoms as readily as real sites, and neither allele fraction nor the
homopolymer flag separates them -- `55,862,239` and `55,862,269` share
AF 0.483 and segregate 0.689 and 1.000.

The consequence is visible in the chain votes: `55,867,094` votes 24 against 20
and `55,871,837` votes 11 against 7 -- near-even edges of the kind already shown
to invert parity -- and the residual break at `55,882,617` reports `agree=0
conflict=0` despite **27 reads spanning** `55,871,837 -> 55,882,617`, while
hiphase links across 20.8 kb on just 2 spanning reads.

Next: the two structures the data does separate on are positional. Three of the
phantoms (`55,883,019` twice and `55,883,023`) are one repeat event called three
times inside 4 bp, and the homopolymer-flagged trio scores 0.825/0.516/0.559. A
cluster-collapse plus homopolymer screen would leave `55,862,239`, `55,862,269`
and `55,871,836` -- keeping the informative site. That is the next thing to
measure, and it has to be measured chromosome-wide, since the 0.559% baseline is
what any default must not damage.
