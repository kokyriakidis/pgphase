# verify_retry.py: a window verdict no check can pass by being blind

Every wrong conclusion in the retry work came from a check that skipped
something silently instead of failing. This tool exists so that cannot happen
again: one command per window, and a `FAIL` whenever the evidence for a claim is
absent rather than negative.

```sh
python3 evaluations/2026-09-16-retry-verify/verify_retry.py \
  --vcf /tmp/retry-on/native.vcf \
  --bam test_data/HG002_chr20_hifi_mapped_to_CHM13_chr20_annotated.bam \
  --truth-map /tmp/truth_hap.tsv \
  --region-left 48176830 --region-right 48229446 \
  --phased-bam /tmp/retry-on/phased.bam \
  --baseline-phased-bam /tmp/retry-off/phased.bam \
  --competitor-vcf "$CD/hiphase/phased.vcf.gz" \
  --output verdict.json
```

## The four errors it makes impossible

| error made | how it was made | what the tool does instead |
|---|---|---|
| "the site is phased and usable" | read `PHASE_SET` from the candidate TSV, where a candidate carries a phase set even when the solve put the alt allele on both haplotypes | reports CATEGORY and emitted GENOTYPE separately; a category-het / genotype-hom site is a named finding |
| "the hole is an artifact of the old site set" | saw a site inside the interval and assumed it bridged, without counting reads | counts spanning reads for every consecutive pair of usable het sites; a link with zero is `unsupported_link` and fails |
| "the gate is clean, so the join is fine" | the read gate cannot see a switch nothing spans | sets `gate_blind` whenever an unsupported link or unscorable site exists in the window, and a PASS requires it false |
| "confirmed switched across the hole" | compared single sites, including one below confidence and one position against itself | collapses records at a position, then declares a switch only where a run of `--switch-run` sites agrees on each side |

Unscorable sites are counted, listed with the reason (`too_few_scored_reads`,
`tied_orientation`, `below_min_confidence`) and fail the verdict. Dropping one
silently is what hid the site this window turns on.

## What it says about chr20:48,176,830-48,229,446

| | retry off | retry on |
|---|---:|---:|
| usable het sites in window (in region) | 44 (2) | **63 (8)** |
| links / unsupported links | 42 / 0 | 61 / **1** |
| unscorable sites | 0 | 2 |
| category-het / genotype-hom | 0 | **1** |
| switches (run of 2 each side) | 0 | **0** |
| competitor sites we do not use | 34 absent | **21 absent, 1 called hom** |
| gate: concordant -> discordant | 0 | 0 |
| `gate_blind` | false | **true** |
| verdict | PASS | **FAIL** |

The retry recovers real evidence -- 13 of the competitor's sites move from
absent to phased, and in-region usable hets go 2 -> 8. It fails on two specific
things, both actionable:

**1. The site the competitor bridges with, we call homozygous.** HiPhase crosses
this interval using two het sites, `48,183,976` and `48,204,383`. We now have
both. At `48,204,383` (`AT>A`, DP 71, 30 ref / 41 alt, AF 0.577) our category is
`NOISY_CAND_HET` but the emitted genotype is `1|1` -- alt on both haplotypes --
so it links nothing. That single site is the difference between a read-supported
chain (`48,183,976 -> 48,204,383` has 1 spanning read, `48,204,383 ->
48,225,786` has 7) and the 41.8 kb jump with **zero** spanning reads that the
solve actually makes.

**2. The join it makes is unsupported, and its correctness is untestable.**
`48,183,976 -> 48,225,786` is 41.8 kb with no read covering both ends. An
earlier note here claimed that block was measured as switched; that is withdrawn
-- with positions collapsed and a two-site run required, no switch is
demonstrated. The right-hand side of the hole has one scorable site
(`48,229,226`); `48,225,786` is below confidence. So the orientation there is
neither confirmed nor refuted, which is exactly why an unsupported link fails on
its own rather than waiting for a truth check that cannot settle it.

## Next

Two changes, in this order, each re-verified with this tool:

1. Find why `48,204,383` is genotyped `1|1` at AF 0.577. It is the bridging site,
   and no linkage policy matters until its genotype is right.
2. Make the retry refuse a link across a spacing with no spanning reads, so the
   output is two blocks rather than one unverifiable one.
