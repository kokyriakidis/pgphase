# Are BAM variations in the initial phasing? Yes -- and then they are withheld

Tested on one window: `chr20:48,176,830-48,229,446` (52.6 kb), a deficit gap
hiphase spans at 100.0% over 252 reads.

## The hybrid solve already carries BAM variation

`process_chunk_hybrid` (`collect_pipeline.cpp:778`) is documented as "BAM
classification -> graph site injection -> BAM profile build -> graph read
injection -> unified k-means". The BAM chunk is the **base** and graph sites are
injected into it, then one k-means runs over the union. Injection is not the
missing step, and it runs in the other direction from the obvious guess.

`graph_authoritative` is off by default, so BAM read evidence at graph-matched
candidates is kept rather than replaced (`backfill_graph_candidate_counts`); the
comment at the call site records that enabling it unconditionally cost +720
discordant reads chromosome-wide.

## Two places the BAM's own sites are then excluded from phasing

1. **With `--recover-gaps`, every non-graph candidate has its category zeroed
   before the solve** and restored only afterwards:

   ```cpp
   discovery_flags.push_back(candidate.lcd_var_i_to_cate);
   if (!candidate.graph_site) candidate.lcd_var_i_to_cate = 0;
   ```

   So a BAM-discovered het site is invisible to the initial k-means and to the
   noisy-region MSA, and visible only to gap recovery.

2. **Hybrid disables the step-4 noisy-candidate k-means by default**
   (`hybrid_collect.cpp:120`, `opts.skip_noisy_kmeans = true`), the stage that
   orients and phases `NOISY_CAND_HET` sites. The comment records why: it "phased
   ~8k extra reads at ~65% error and poisoned the BAM-shared core".

## Measured on the window

| arm | blocks | spans target | tagged | read accuracy | in-gap het phased |
|---|---:|---|---:|---:|---:|
| hybrid, default (recovery on) | 3 | no | 464 | 99.78% | **2 / 10** |
| hybrid, `--recover-gaps` off | 2 | no | 205 | 100.00% | **3 / 11** |
| hybrid + `--keep-noisy-kmeans` | 2 | no | 205 | 100.00% | 2 / 4 |
| **BAM-only channel** | 2 | no | 392 | **100.00%** | **11 / 11** |

The BAM-only channel phases **every** het site in the interval into a single
block spanning **48,147,227-48,229,226 (82.0 kb, 16 sites) at 100% read
accuracy** -- covering all but the last 220 bp of the interval a competitor
spans. Hybrid leaves that interval unphased and fragments into three blocks, the
nearest of which stops exactly at the gap's left edge (`48,162,480-48,176,830`,
2 sites).

Neither switch explains it on its own: turning recovery off (which stops the
category zeroing) moves 2/10 to 3/11, and restoring the noisy k-means moves
nothing. So the suppression is not a single flag.

## The sites are the same; only their fate differs

Comparing candidate tables position by position in the interval:

| | het sites | phased |
|---|---:|---:|
| BAM-only | 9 | **9** (all into PS 48147227) |
| hybrid | 10 | **2** |

Same positions, same categories, 7 shared `NOISY_CAND_HET` indels. Two sites
differ in genotype: `48,177,781` and `48,225,787` are `NOISY_CAND_HET` to the BAM
channel and `NOISY_CAND_HOM` in hybrid -- het called homozygous, the same error
class found at 7 of hiphase's sites chromosome-wide, and a site called homozygous
can never link.

And noisy-region **detection** is identical: run at verbosity 2, both paths emit
the same 2,862 noisy-region diagnostic lines over the same regions. Both find the
evidence; one phases it.

## What this means for the deficit

At this gap the deficit is not discovery, not MSA verification (every
`msa_verified = 1` site here segregates 0.89-1.00 against read truth), and not
missing BAM injection. It is that the hybrid solve withholds the BAM channel's own
noisy-site phasing, which on this window would have produced an 82 kb block at
100% accuracy.

The fix direction this points at is scoped rather than global: let BAM-discovered
noisy candidates participate in phasing **inside gap intervals only**, which is
where the BAM channel demonstrably gets them right, instead of the chromosome-wide
enabling that previously cost ~8k reads at ~65% error. That is a narrow change
with an existing gate to test against -- the 31-gap accurate-deficit list and the
read-level concordance gate.

Not attempted here: the per-gap version of that change. The measurement above is
what justifies building it, and the same three-arm comparison on this window is
its regression test.
