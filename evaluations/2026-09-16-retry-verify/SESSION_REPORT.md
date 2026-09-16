# Session report: gap phasing in pgphase, 2026-09-16

16 commits, `281ec0d..ccf74a4`, pushed to `origin/main`. 99 files, 17,363
insertions; 818 of those in `src/` across 9 files, the rest evaluations. All five
unit-test binaries pass; every source change is off by default or scoped to a
window.

## What shipped in `src/`

| change | commit | default |
|---|---|---|
| MSA-verified SNPs may bridge recovery gaps (two chained gates, not one) | `814dbe7` | on |
| The homopolymer tier also runs in the reprojected pass | `814dbe7` | on |
| Gap-recovery link votes read pre-filter read labels | `306a58c` | on |
| Allele attachment requires a joined gap | `84cef4c` | **off** (`--gap-allele-attach-join-only`) |
| Every qualifying independent gap block is emitted, not just the largest | `a302c67` | on |
| pgbam stitch thresholds exposed in the hybrid subcommand | `fdae47f` | unchanged |
| Retry a window the solve could not phase or connect, BAM sites admitted there | `7018da4`, `0a9640a` | **off** (`--retry-unphased-with-bam`) |

## The three findings that matter

**1. The noisy-region MSA never runs when `--recover-gaps` is on.**
`collect_var_run_phasing` guards it with `if (!opts.recover_gaps)`, deferring it
to the recovery pass, so the noisy het class does not exist in the normal solve.
Inside `chr20:48,176,831-48,229,447` the chunk held 1 clean het SNP, 2 clean het
indels, 1 repeat indel, 64 clean hom and 606 low-coverage catalog sites and
**zero** `NoisyCandHet`, while `collect-bam-variation` on the same interval calls
8 and phases the window into one block. Three other explanations were checked and
ruled out first: the 50 kb `max_noisy_reg_len` cap (the BAM channel yields the
same 8 sites over 52 kb and 152 kb), the `private_keys` branch that zeroes that
cap (entered only with `--private-sites-vcf`), and `skip_noisy_kmeans` (read
*inside* the step that never ran). Running the step in the retry took in-gap
phased hets from 2 to 8 and grew the left block from 14.3 kb to 82.0 kb.

**2. The site the competitor bridges with, we call homozygous.** HiPhase crosses
that interval using two het sites. We have both. At `48,204,383` (`AT>A`, DP 71,
30 ref / 41 alt, AF 0.577) our category is `NOISY_CAND_HET` and the emitted
genotype is `1|1` — alt on both haplotypes — so it links nothing. That one site
is the difference between a read-supported chain (1 spanning read, then 7) and
the 41.8 kb jump with zero spanning reads the solve makes instead. This is the
next thing to fix; no linkage policy matters until the genotype is right.

**3. Admitting BAM sites chromosome-wide is a bad trade, so the window
confinement is load-bearing.** On chr20 it doubles read Hamming error, 0.878% to
1.837% (1,831 to 3,415 discordant reads), fragments blocks 238 to 408, and spans
14 of 196 gaps. Confining admission to windows the solve actually failed is what
separates the gain from that cost.

## Deficit, re-measured against the pipeline

The earlier target list came from the graph-only pass-1 inventory, which does not
describe what the pipeline leaves open. Against a fresh whole-chr20 run: **47 of
196 gaps (1.20 Mb)** are spanned by a competitor, 149 gaps (10.67 Mb) by nobody.
Scoring those spans against read truth — some are worse than chance, where our
abstention is the better answer — leaves **31 gaps, 0.81 Mb** where a competitor
spans at ≥98% over ≥50 reads. That is the real target set. In those gaps we
already hold **every** substitution the competitor uses; the shortfall is
entirely indels, dominated by the noisy-candidate class.

## Retractions

Three claims of mine were withdrawn after measurement, each in the committed
record as well as here:

- **"Zero phased sites inside the gap at every chunk size."** Never measured; the
  sweep printed no such column. It is 2 phased hets against 64 homozygous.
- **"The gap closes, +107 concordant reads, no flips."** The frame spanned a
  41.8 kb interval no read crosses; the same arm at a different read context
  produced the opposite orientation. Luck, not a closure.
- **"Confirmed switched across the hole."** With positions collapsed and a run of
  two sites required per side, no switch is demonstrated — the right side has one
  scorable site. The join is unsupported and untestable, which is weaker than
  switched.

## Tooling, so these do not recur

`evaluations/2026-09-16-retry-verify/verify_retry.py` is the regression test for
any gap-phasing change: one verdict per window, built so no check passes by being
blind. CATEGORY and emitted GENOTYPE are separate columns; spanning reads are
counted for every consecutive site pair and a zero fails on its own; unscorable
sites are listed with a reason and fail; a switch needs a run of sites per side
after collapsing one position's records; and `gate_blind` is set whenever an
unsupported link or unscorable site exists, so a clean read gate cannot carry a
PASS. It reproduced all four of the errors above on first run.

The four lessons are also in the `pgphase-gap-diagnosis` skill, so the next
session starts with them.

## State and next steps

Working tree: `.devcontainer/devcontainer.json` (pre-existing, untouched), the
gap-lab tool and its verdicts, and stray `.md` copies at the repo root made for
artifact export.

1. Why `48,204,383` is genotyped `1|1` at AF 0.577 — the bridging site.
2. Make the retry refuse a link across a spacing with no spanning reads, so the
   output is two honest blocks rather than one unverifiable one.
3. Re-verify with `verify_retry.py`, then run the 31-gap target list.
