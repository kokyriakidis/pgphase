# Private MSA het SNPs as gap-bridge anchors (`--gap-bridge-private-snps`)

Measures the opt-in flag added to `src/collect_phase.cpp` on whole chr20, against
the same binary with the flag off. Both arms reuse `/tmp/chr20-evidence-v6.gapev`,
so the only difference between them is the flag.

```sh
./evaluations/2026-09-15-private-snp-bridge/run.sh
python3 evaluations/2026-09-15-private-snp-bridge/compare.py \
  --baseline /tmp/pgphase-private-snp-bridge/off \
  --run /tmp/pgphase-private-snp-bridge/on \
  --output evaluations/2026-09-15-private-snp-bridge/off_vs_on.json
```

## Flag-off reproduces the established baseline

The working tree also carries three changes that are *not* behind the flag
(`select_graph_gap_bam_reads` admitting private sites, pass 1 reprojecting into a
scratch chunk instead of overwriting pass 0's solved proposal, and the new
`LEFT_LINK_PS`/`RIGHT_LINK_PS` report columns with `split`/`vetoed` statuses).
The flag-off arm is identical to the accepted baseline on every headline metric
— 276 initial gaps, 112 joined, 207,426 phased reads, 1,400 discordant,
0.675% read Hamming, N50 1,064,616 bp, 177 phase sets — so those three changes
cost nothing at the default operating point. An earlier same-day sweep that
reached 115-116 joins at 0.720%/0.713% Hamming came from the anchor
segregation experiment recorded as reverted in `src/collect_phase.cpp`, not
from anything still in the tree.

## What the flag does

| metric | off | on | delta |
|---|---:|---:|---:|
| gaps joined (of 276) | 112 | 113 | +1 |
| phased reads | 207,426 | 207,373 | -53 |
| reads evaluated | 207,417 | 207,364 | -53 |
| discordant reads | 1,400 | 1,386 | -14 |
| read Hamming | 0.675% | 0.668% | -0.007 pp |
| switch errors | 249 | 249 | 0 |
| flip errors | 671 | 657 | -14 |
| block N50 | 1,064,616 | 1,064,616 | 0 |
| block auN | 1,234,653 | 1,238,045 | +3,392 |
| phase sets | 177 | 176 | -1 |

Read-by-read against the fixed baseline cohort:

- **0 previously concordant reads became discordant** — the regression gate holds.
- 3 discordant reads became concordant; 21 reads newly phased (19 concordant).
- 74 reads lost their PS tag entirely (61 of them previously concordant).

## Where the single join and the read loss come from

The newly joined gap is `chr20:36,332,599-36,381,019`, which moves from `split`
(both flanks linked, but to different proposal phase sets) to `joined`. Its right
flank PS 36,381,019 (3,978 reads) merges into PS 36,317,511 (197 reads), giving
4,150 of the 4,175 reads — 25 reads drop out of the merge. The remaining 49 lost
reads sit at two unrelated phase sets (38,236,063: 38 reads; 1,142,089: 11), which
is the read-level side of the same flag: a read whose only informative gap site is
a private SNP now needs that observation confirmed directly in the BAM
(`clean_snp_has_bam_observation`), and an unconfirmed read contributes nothing and
can fall below its output margin. The flag therefore both admits new bridge
anchors and suppresses unconfirmed private-SNP read evidence; the error rate
improves because the suppressed population was enriched for discordant reads.

## Verdict

Safe but a poor coverage trade: one extra gap join and 3.4 kb of auN, paid for
with 61 concordant reads losing their phase tag. Since "more reads phased is the
goal" decided the `--gap-independent-min-reads` default, the flag stays **off by
default** and remains available for per-region diagnostic use. `--gap-bridge-private-snps`
requires `--recover-gaps` and errors out otherwise.

## Files

- `run.sh` — both arms plus read-level evaluation.
- `compare.py` — per-read transitions, coverage change and gap-status movement.
- `off_vs_on.json` — the comparison above.
