# Site-level evidence for automatic gap recovery

Automatic recovery now lets reads rejected by the whole-consensus MSA score
margin contribute individual allele observations. This applies only inside
`--recover-gaps`; the existing score margin (24), minimum block-link reads (2),
k-means, and stitching thresholds are unchanged.

For each rejected full-cover read, compose its alignments through both fixed
MSA consensuses back to reference coordinates. At each MSA heterozygous site,
require an exact reference or alternate allele, three matching reference/read
bases on each flank, and agreement between both alignment paths. Unexpected
alleles, incomplete coverage, ambiguous bases, noisy flanks, and path disagreement
abstain. No consensus cluster or HP assignment is manufactured. SNP-only tier 2
still excludes added indel observations; tier 3 can use them.

Before adding observations, the expanded reference/alternate counts must pass
the existing min_af/max_af heterozygosity limits. Counts and read profiles are
updated together only for accepted additions. Original consensus-derived
observations are retained. Ordinary phasing then determines whether the new
observations provide useful connectivity.

## Regression caught during implementation

The first working prototype lacked the expanded allele-frequency gate. In the
chr20 17.6 Mb case it added 40 phased reads, including two discordant reads.
The new reads encountered MSA SNPs labelled heterozygous from differences
between consensuses, despite overwhelmingly reference support after expansion.
This exposed why agreement between two alignment paths alone is insufficient.
The final gate uses existing genotype thresholds and is applied to all regions;
it is not a region-specific exception. Prototype summaries are preserved in
`before_af_gate/`. The final chr20 check adds two reads and has zero read-truth
errors. A unit test verifies that rejected additions mutate neither counts nor
profiles.

An earlier integration check also caught that `make_ref_read_aln_str` leaves
alignment coverage bounds unset. The new full-cover-read path sets its own
bounds before allele calling; it does not change the existing composition
function's behavior for other callers.

## Reproduction

Run `python3 evaluations/2026-09-13-site-observation-recovery/run.py`, then
`python3 evaluations/2026-09-13-site-observation-recovery/summarize.py`.
The manifest selects exactly the same six regions and settings as
`../2026-09-13-auto-gap-validation/`; comparison.tsv compares against that frozen
recovery baseline, not merely against the clean pass. Data requirements are
unchanged. Heavy BAMs and input matrices go to `/tmp/pgphase-site-observation-recovery`.
The final binary fingerprint and per-arm commands are retained here.

Also run the original two chr20 cases and split-chunk fixture with:
`OUT=/tmp/pgphase-site-observation-regression bash evaluations/2026-09-13-auto-gap-recovery/run.sh`.
The fixture verifier checks block preservation, tier ordering, and read truth.
Build and all unit tests passed, including SNP/indel allele calls, path
inconsistency, coverage/flank rejection, third alleles, and the AF admission gate.
The only rebuild warning came from the existing abPOA SIMDMalloc header.

This is targeted read-truth validation. Both alignment paths use the same read
and are a consistency check, not independent observations. No singleton bridge
fallback or chromosome-wide accuracy/NGC50 claim is introduced.

## Final results versus the previous recovery implementation

| Region | Previous → new blocks | Additional phased reads |
|---|---:|---:|
| chr12 46.8 Mb | 4 → 4 | 3 |
| chr12 63.8 Mb | 1 → 1 | 0 |
| chr18 46.2 Mb | 3 → 1 | 12 |
| chr18 33.3 Mb | 3 → 2 | 0 |
| chr20 61.8 Mb | 1 → 1 | 0 |
| chr20 17.6 Mb | 2 → 2 | 2 |

Recovery now joins 6/11 internal gaps across these windows, versus 3/11
previously. It phases 17 additional reads over the previous recovery version
(93 over the clean initial pass). No region increases its discordant-read or
switchflip count. The one existing error in chr18 46.2 Mb remains. All 2,293
initially tagged reads retain one uniform PS/parity transformation per initial
block. SNP evidence remains cumulative, tier-2 indels remain excluded, and
successful joins stop further escalation. Original 15 Mb, 47 Mb, and split-chunk
fixtures all still reach one block at both tested MSA margins with zero errors.

At chr12:46774503, the final input matrix has 10 usable observations (6 ref,
4 alt), up from 7 (4 ref, 3 alt). Its target gap remains open under the unchanged
two-read rule; this change does not manufacture evidence for its singleton link.
