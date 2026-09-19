# What the pipeline is for: reads correctly separated

Truth-het counting was the wrong headline. The question is whether a region
gets phased, whether its reads end up correctly separated, and if not what the
competitor has that we do not. This record measures all three.

## The metric

`separated()` = of the reads overlapping a window that the parental truth map
can score, the fraction that ONE phase set places on the right parent.
Contiguity and correctness in one number: purity alone hides fragmentation --
two immaculate half-blocks separate nobody -- and block count alone hides
switches. It is now asserted per window and per arm as `min_separated`.

## The deficit is fragmentation, not error

Scored over every gap in the in-chunk arm's chr20 output against hiphase's
phased BAM, both against the read truth map:

| | windows |
|---|---:|
| hiphase separates >= 10 points more of the reads | **123** |
| we separate >= 10 points more | 3 |

The pattern is the same in almost every one:

| window | ours | hiphase |
|---|---|---|
| 55,999,194 | 100.0% purity, 2 blocks, 63 of 246 reads -> **25.6%** | 99.5%, 1 block, 214 -> **86.6%** |
| 23,460,963 | 100.0%, 2 blocks, 64 of 227 -> 28.2% | 99.5%, 1 block, 199 -> 87.2% |
| 12,717,796 | 100.0%, 2 blocks, 74 of 219 -> 33.8% | 100.0%, 1 block, 201 -> 91.8% |
| 5,309,406 | 98.6%, 2 blocks, 74 of 172 -> 42.4% | 99.4%, 1 block, 169 -> 97.7% |

Where we phase, we are essentially never wrong. We split. That is why the
earlier truth-chain selection found 11 targets and this finds 123: the chain
test asks whether a join COULD exist, this asks what the output achieves.

The guard window makes the recovery's own case on the same metric:
26,029,591 separates **0.99** of its reads under `--in-chunk-recovery` against
**0.03** under the post-hoc default.

## What the competitor has that we do not

For the two panel windows with the most headroom, every het site hiphase
phases inside the gap, traced through our own channels:

| site | in catalog | filtered as | our candidate |
|---|---|---|---|
| 4,766,934 | no | -- | **CLEAN_HET_INDEL** |
| 4,785,719 | yes | `low_depth` | -- |
| 4,786,025 | yes | `high_af` | REP_HET_INDEL |
| 4,791,668 | yes | `high_af` | -- |
| 60,048,247 | yes | `ref_only` | -- |
| 60,056,028 | yes | `ref_only` | REP_HET_INDEL |
| 60,058,235 | yes | not filtered | **CLEAN_HET_INDEL** |

Three classes, in priority order:

1. **We hold the site, classified CLEAN_HET_INDEL, and it still does not
   produce a phased record inside the gap** (4,766,934 and 60,058,235). We
   emit 0 and 1 records respectively across those gaps. This is the strongest
   lead: the evidence is present and admitted-looking, and something between
   the candidate table and the output drops it.
2. **The catalog holds the site and the graph channel's own counts reject it**
   -- `high_af`, `ref_only`, `low_depth` -- while the alignment sees a
   heterozygote there. Recovery is supposed to rescue exactly this and does
   not.
3. **`REP_HET_INDEL`**, the repeat class, which the solve excludes by
   construction.

None of the seven is a linking failure. The classification the window gate
prints for these windows -- NO SITES -- is correct, and the sites are missing
for three different reasons that need three different fixes.

## Next

Class 1 on 4,766,934: we hold a CLEAN_HET_INDEL at a position hiphase phases,
inside a gap where we emit no records at all. Follow that candidate from the
table to the writer and find where it is lost.

## The shipped default is the bad configuration

Found while A/B-ing the emission fix. Every in-chunk measurement in this
project was taken with `--no-anchored-stage2`; the flag defaults ON, so the
configuration a user gets is not the one that was measured. Whole chr20, same
binary, only that flag differing:

| arm | tagged | read blocks | discordant | read hamming | VCF blocks |
|---|---:|---:|---:|---:|---:|
| `--in-chunk-recovery` (anchored stage 2 ON, the default) | 219,247 | 422 | 9,497 | **4.332%** | 491 |
| `--in-chunk-recovery --no-anchored-stage2` | 219,059 | 323 | 2,543 | **1.161%** | 324 |

Anchoring stage 2 costs 6,954 extra misplaced reads and 99 extra read blocks
here -- it is worse on both axes, contiguity and correctness. This also means
the window panel, whose `inchunk` arm passes only `--in-chunk-recovery`, has
been scoring the 4.332% configuration: its separated() floors of 0.23-0.53
against hiphase's 0.86-0.98 are measured on the arm nobody should run.

Not yet fixed, because which way to fix it is a real decision: turn the default
off for this path, or find why anchoring hurts when recovery merges sites. The
anchored round was adopted on the evidence that it helps the graph arm WITHOUT
recovery, and that evidence has not been re-examined against these numbers.
