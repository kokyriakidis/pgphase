# Next target: HiPhase's third-largest chr20 block

The latest frozen pgphase chromosome output still splits HiPhase's
**24,103,779–25,944,471** block, PS24103779, into three blocks.
HiPhase's span is **1,840,693 bp**, its third-largest, with 2,648 phased
heterozygote rows. Its second-largest block is now covered exactly. The
largest block retains the separately investigated 12,277 bp terminal
difference; it has no remaining internal split.

| Pgphase PS | First–last phased heterozygote | Span | Rows |
|---|---|---:|---:|
| 24103779 | 24,103,779–24,121,713 | 17,935 bp | 4 |
| 24517362 | 24,131,707–25,855,631 | 1,723,925 bp | 2,102 |
| 25805896 | 25,855,633–25,964,482 | 108,850 bp | 229 |

Two connections are missing: **24,121,713–24,131,707** (9,994 bp boundary
distance) and **25,855,631–25,855,633** (2 bp). The third pgphase block extends
beyond HiPhase's right endpoint. Closing only one connection would leave the
HiPhase span split.

## Original-read comparison

Eligibility comes from the original primary BAM alignments, with parental
truth available. Unphased reads remain in every denominator, including reads
absent from a tool's output. HiPhase's start/end, CIGAR and sequence SHA256
match all **7,591** original alignments overlapping the block. Correctness
uses each output phase set's whole-cohort parental orientation. Rescue phase
sets >=1e9 do not count as connected core.

| Interval | Tool | Correct | Discordant | Unphased | Correct / all | Largest correct core |
|---|---|---:|---:|---:|---:|---:|
| Whole HiPhase block | Pgphase | 7,365 | 48 | 178 | 97.02% | 6,705 |
| Whole HiPhase block | HiPhase | 7,150 | 1 | 440 | 94.19% | 7,150 |
| First seam, 115 reads | Pgphase | 89 | 0 | 26 | 77.39% | 29 |
| First seam, 115 reads | HiPhase | 93 | 0 | 22 | 80.87% | 93 |
| Second seam, 91 reads | Pgphase | 53 | 38 | 0 | 58.24% | 37 |
| Second seam, 91 reads | HiPhase | 20 | 0 | 71 | 21.98% | 20 |

**The first seam is the next qualifying gap:** HiPhase exceeds 80% and
connects all 93 correct reads. Pgphase needs at least 93 total correct and 93
correct connected-core reads. A union of its existing cores would connect
only 40 correct reads; 49 other correct reads currently have rescue labels.
It also needs at least four additional correct assignments to match HiPhase.

The first seam already has an owning-chunk regression,
`gap_complementary_insertion_recall_cannot_invert_its_snp_flanks`. Keep its
allele, read and orientation checks. The chromosome audit finds 115 eligible
original overlaps; the regression's historical minimum of 89 scored/correct
reads is not the denominator for the new 80% comparison.

The second seam is a separate accuracy problem. HiPhase's VCF phase set spans
it, but its local read assignments do not meet 80%. A raw VCF block connection
must not be reported as a read-certified closure. Pgphase's 38 discordant
reads require investigation before accepting a join there.

Disjoint 50 kb flank cohorts establish pgphase's gauges independently of
each seam's overlapping reads. First seam: the left core has 54 unanimous
hap1-maternal votes, the right 210 unanimous hap1-paternal votes. Second
seam: the middle core has 151 unanimous paternal votes on the left, and the
right core has 191 maternal / one paternal vote on the right. The complete
flank counts, including other phase sets and abstentions, are in next-block.json.

## Boundary calls and reproduction

At the first seam pgphase retains complementary 4T/8T insertions at
24,121,713 and a heterozygous CT>C deletion at 24,131,707. HiPhase phases the
4T/8T insertion alleles in a multiallelic record (GT3|1, GQ3) but calls the
deletion homozygous reference. These observations identify a candidate for
the calibrated insertion approach; they do not yet certify a physical bridge.
At the second seam, pgphase has phased T>C and T>G SNPs. HiPhase's nearby
repeat SNP/indel records are NoCall. Exact records are in boundary-calls.json.

The production binary is unchanged, SHA256
`15534777ecd8382e370ade340e6e8998a3fc1b1ab4ab390281b3ff4a22bbedb4`.
Inputs are pgphase's `test_data/tmp_gap_fix71/frozen_final/0`, HiPhase's
`test_data/tmp_gap_fix48/competitor/hiphase_dv`, the original annotated BAM,
and `test_data/derived/chr20_truth_hap.tsv`. Rank all 191 HiPhase blocks against
all 263 pgphase blocks and reproduce the original-read audit with:

```bash
/home/kokyriakidis/micromamba/envs/bench-phasers/bin/python \
  evaluations/2026-10-05-third-largest-block-target/find_next_block.py
```

Rankings are in ranked-blocks.json and the selected block, both seams,
alignment checks and parental metrics are in next-block.json. This task
identifies and audits the next target; production phasing and test expectations
are unchanged. No pipeline rerun or rebuild was needed.
