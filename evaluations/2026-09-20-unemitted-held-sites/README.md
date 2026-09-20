# The 254 records we hold and do not emit: two divergences removed, three causes eliminated

After the BAM path became a faithful port (`d712c02`), 370 of upstream's phased
records are still missing from ours out of 118,270. Bucketed on a fresh run:

| bucket | records |
|---|---:|
| we hold a candidate and emit nothing | **254** |
| no candidate at all | 85 |
| we emit the position with a different allele | 31 |

The 254 are 207 `NOISY_CAND_HET`, 46 `CLEAN_HET_SNP`, 1 `CLEAN_HET_INDEL`, all
phased, and **251 of them carry `HAP_ALT/HAP_REF = 0/0`** -- both haplotype
consensus alleles are reference, so no ALT is written and the record is dropped
for carrying no alternate. Our depths are not the problem: over the 254 our DP
equals upstream's at 118 sites and exceeds it at 135.

## The mechanism, measured

At `chr20:3,870,827` (INS, `NOISY_CAND_HET`, DP 72, 43/29, phased) the two
haplotype profiles are `[25,17]` and `[31,27]` -- the alternate is a minority on
BOTH haplotypes, so the argmax returns reference twice. The totals (42 and 58
against DP 72) give the arithmetic away: 28 reads are unlabelled, and an
unlabelled read is counted toward both haplotype profiles (ours and upstream
alike, `assign_hap.c:299-302`), so it dilutes both.

Comparing labelling against upstream at that locus, over the 75 reads spanning
it: we label 42 (56%), all 42 truth-consistent; upstream labels 69 (92%), 55 of
69 truth-consistent. Over 20 sampled sites from the 254 the pattern holds and
is not a one-off -- ours 221 of 1,086 spanning reads labelled (20%) at 65.6%
truth-consistency, upstream 488 (45%) at 60.2%. Upstream labels about twice as
many reads at these loci and is somewhat less accurate per label, and that
density is what lets its profiles resolve. Per-site it goes both ways
(`37,336,382`: ours 87%, upstream 15%; `44,704,972`: ours 0%, upstream 100%).

## Two divergences found and removed, both inert here

Both are real differences from upstream, both are now off for the BAM path, and
both left whole-chr20 output **byte-identical** -- so neither is the cause:

| divergence | upstream | measured effect |
|---|---|---|
| `read_to_cons_allele_score` returned 0 for any candidate with `msa_insertion_alts` unless `gap_link_supported`, silencing the whole noisy-het class | no such gate; every candidate in the mask votes (`assign_hap.c:127-147`) | inert: `msa_insertion_alts` is only populated by the refresh path, which the faithful configuration turns off |
| the complement inference ran only at sites with exactly two allele slots | unconditional (`assign_hap.c:141-142`) | inert: the candidates at these loci are biallelic, so the inference already applied |

They are kept removed because the directive is a faithful port, and recorded as
inert so nobody re-measures them expecting a gain.

## A third cause eliminated by instrumenting the scorer

Tracing every score a read receives before it goes unlabelled, at
`chr20:3.86-3.88 Mb`: every variant in those reads' span in the clean round is
homozygous (`cons = [1,1]`, score `+1` on both haplotypes, 3,644 of 3,650 traced
observations). A homozygous site cannot distinguish haplotypes, so it scores
without counting toward `n_vars_used`, and the read is correctly left unlabelled
-- upstream does the same (`assign_hap.c:173`). The clean round's abstention is
right; the question is entirely about the noisy round.

## What is left

The region's only het candidates are the 17 `NOISY_CAND_HET` sites, so
everything depends on the noisy round, and there the state is self-reinforcing:
unlabelled reads dilute both profiles, the diluted argmax returns reference
twice, and a site whose consensus is `[0,0]` then cannot label any read. Whether
upstream escapes this by its pivot choice, by its sweep order, or by seeding
consensus values earlier is the next measurement, and it needs upstream's own
per-iteration state rather than more reading of ours -- the instrumented copy at
`/tmp/lcd-instr` is the tool for it.

Not yet attempted, and worth stating so it is not assumed done: no change in
this record improves the 254.
