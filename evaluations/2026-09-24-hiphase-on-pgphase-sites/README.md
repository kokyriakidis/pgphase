# HiPhase on pgphase's current chr20 sites

Date: 2026-09-24. HiPhase `1.6.0-ac3f399` (the version in the frozen
competitor lock) was run on the unmodified annotated HG002 chr20 BAM and
CHM13 reference, with pgphase's final full-chr20 graph+BAM VCF as its *only*
variant input. The VCF was bgzip-compressed and tabix-indexed without changing
records. HiPhase used `--threads 8 --ignore-read-groups` and otherwise default
phasing thresholds. Inputs and outputs are under
`/tmp/hiphase-on-pgphase-final-chr20/`; the pgphase baseline is
`/tmp/pgphase-try5538-final/`.

## Whole-chromosome result

| Metric | pgphase graph+BAM | HiPhase on the same pgphase VCF |
|---|---:|---:|
| VCF heterozygotes phased | 61,646 | 59,381 |
| VCF phase sets | 463 | 211 |
| Truth-evaluable tagged reads | 236,527 | 228,286 |
| Truth-correct tagged reads | 228,190 | 216,224 |
| Truth-discordant tagged reads | 8,337 | 12,062 |
| Per-PS read truth purity | 96.48% | 94.72% |

Each PS is oriented independently against parental read truth before counting
correct and discordant reads. Both tools use the same BAM and variant records;
HiPhase extracts and phases alleles with its own method, including its
realignment path. This experiment does not isolate the phase optimizer alone.

## Current short pgphase gaps

Group pgphase's phased heterozygous VCF records by PS and take each PS's
minimum and maximum variant positions. Sort these intervals and keep each
consecutive, non-overlapping pair whose boundary sites are less than 10 kb
apart. Exclude pairs with either endpoint in the 26.0-29.5 Mb centromeric
interval. This gives 115 current short gaps. HiPhase closes a gap when the
*exact same* two boundary variants `(POS, REF, ALT)` are phased into one
HiPhase PS. It closes 65/115.

For each closed gap, count truth-labeled primary reads carrying that HiPhase PS
in the interval from 10 kb before the left boundary through 10 kb after the
right boundary. Local purity is the larger parental-orientation count divided
by all such reads. Of the 65 joins, 64 have at least 20 truth-scored local
reads; 35/64 have at least 98% local purity. Of those 35, 34 also have at
least 95% purity on each 10-kb flank and the same majority orientation on
both flanks. Some joins are much weaker: 17,839,393-17,839,398 has 57.7%
local purity (137 reads), and 49,764,957-49,768,952 has 63.3% (139 reads).
The per-gap evidence is in `short_gap_replay.tsv`.

The previously investigated 51.262-51.287 Mb and 55.381 Mb chains are also
joined by HiPhase on pgphase's VCF; pgphase's current implementation has
closed both. HiPhase keeps the unsupported 51.235-51.262 Mb boundary split.
Thus shared sites suffice for *some* missed joins, but neither universal
closure nor universal correctness follows from the two example windows.
