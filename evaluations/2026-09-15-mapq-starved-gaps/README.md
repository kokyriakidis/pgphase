# Why competitors phase chr20 gaps we leave open: the MAPQ floor

Diagnosis of one gap the competitors bridge correctly and we do not, followed by
a chromosome-wide screen for the same pattern. All of it reads frozen outputs and
the surjected BAM; only the MAPQ sweep reruns pgphase, on a 148 kb window.

## The gap

`chr20:25,834,662-25,883,079` (48.4 kb), our status `split` — both flanks link,
to different local phase sets — with 279 of 296 reads carrying no PS tag.
HiPhase spans it inside one block (PS 24,103,779) using 103 phased het sites
(90 SNP, 13 indel).

We called 18 candidates in that window, 2 of them `CLEAN_HET_SNP`, and matched
only 3 of HiPhase's 103 sites. The cause is not classification and not variant
discovery: at the same positions DeepVariant reports DP 59-76 with balanced
allele depths (26,33 / 40,35 / ...), while our candidates report DP 2-10.

Input coverage in the gap is 68.2x. The read population is
**MAPQ 3 (166 reads) and MAPQ 4 (41 reads)** out of 296 — see
`gap_mapq_histogram.tsv`. Our default floor is `kDefaultMinMapq = 30`
(`src/phasing_types.hpp:31`), applied at gap read selection
(`src/gap_recovery.cpp:46`), GAF parsing (`src/graph_query.cpp:421`), injection
(`src/hybrid_inject.cpp:570`) and BAM output (`src/collect_bam_output.cpp:410`).
That floor leaves 3.7x of the 68.2x, which is why the gap looks evidence-free.

## The MAPQ floor sweep (`mapq_floor_sweep.tsv`)

Same 148 kb window, same binary, only `-q` changed. Read truth is the
diplinator BAM, restricted by read name after phasing.

| -q | gap depth | clean het SNPs | of HiPhase's 103 | phased reads | discordant | window Hamming | blocks | gap joined |
|---:|---:|---:|---:|---:|---:|---:|---:|:--|
| 30 | 3.7x | 2 | 3 | 395 | 0 | 0.00% | 2 | no |
| 20 | 4.2x | 3 | 4 | 414 | 0 | 0.00% | 2 | no |
| 10 | 13.4x | 35 | 26 | 469 | 0 | 0.00% | 2 | no |
| 5 | 16.5x | 34 | 26 | 494 | 0 | 0.00% | 2 | no |
| 1 | 68.2x | 302 | 86 | 734 | 37 | 5.04% | 1 | **yes** |

At `-q 1` the gap stops existing: one 375 kb block instead of two (N50 166 kb ->
375 kb), 734 of 776 reads phased, and 86 of HiPhase's 103 sites recovered, 85 as
clean het SNPs. The cost is 37 discordant reads, and they are entirely a
low-MAPQ phenomenon — 14.5% error on the MAPQ 1-9 population, **0% at MAPQ >= 10**.

## HiPhase runs at MAPQ 5, which does not admit these reads

This gap's reads sit at MAPQ 3-4, below HiPhase's own default floor of 5. HiPhase
tagged 89 reads in the window and exactly 15 of them are MAPQ 1-9 — the 15 reads
at MAPQ 5-9. So HiPhase did not phase this gap by admitting the MAPQ 3-4 bulk.

The asymmetry is where the signal is: DeepVariant **called** the 103 sites using
the low-MAPQ reads (DP 59-76), and HiPhase then **phased** those calls using only
MAPQ >= 5 reads. We apply one floor to both jobs, so discarding the reads for
phasing also discards them for discovery, and the anchors never exist. Our own
`-q 1` arm confirms discovery works once the reads are present (302 clean het
SNPs), and that assigning those same reads to haplotypes is what injects the
error.

## How general is it (`unresolved_gaps.tsv`, `lowmapq_verification.tsv`)

Of 164 unresolved chr20 gaps spanning 15.2 Mb: 20 gaps (2.99 Mb) are MAPQ-starved
(input depth >= 20x, >= 50% of it below MAPQ 30), 5 have genuinely low input
depth, and 139 keep their depth through the floor and are limited by something
else.

Scoring the reads we leave unphased against truth splits the 20 sharply:

- **3 gaps, 149 kb, 765 reads** — competitors reach 96-100% (`25,834,662-25,883,079`,
  `25,944,471-25,986,123`, `26,029,591-26,088,679`, all just proximal to the
  centromere). Real, recoverable.
- **17 gaps, 12,393 reads** — LongPhase spans 3.8-93.2% and HiPhase 11.8-100.0%
  on the reads we drop, so this group is not uniformly noisy: three of the 17
  sit at LongPhase 93.2%, 81.7% and 81.2% (26,407,782-26,417,736;
  29,496,882-29,510,999; 32,521,867-32,907,753) and fail the trustworthy test on
  the 95% threshold rather than on read count, and the HiPhase 100.0% belongs to
  32,432,281-32,449,137 on only 10 reads, where LongPhase manages 61.9% on 21.
  The separation from the trustworthy group is therefore clean at the extremes
  and graded in between. The group also includes the 1.93 Mb
  window `27,133,886-29,068,243`, where LongPhase manages 54.1% and HiPhase
  56.4% on the reads we drop. The earlier CHECKPOINT audit was right to call
  that one a legitimate abstention; the 2,596 competitor-phased sites inside it
  are not evidence of correct phase.

A truth-free discriminator is visible but weak: median |VAF-0.5| at competitor
het calls is 0.052 across the three trustworthy gaps against 0.115 across the 16
noise-verdict gaps that have a value, with one trustworthy gap at 0.220 —
overlapping ranges, so it cannot gate on its own. NM tags are absent from the surjected
BAM, so per-read divergence was not available as a second axis.

## Implication

The actionable change is a split floor: keep the MAPQ 30 requirement for
haplotype assignment and output, and let candidate discovery inside a gap see
reads down to a lower floor. `min_mapq` is currently one knob driving both
(four call sites above). The `-q 10` row shows the conservative half is free in
this window — 35 clean het SNPs instead of 2, +74 reads phased, zero discordant
— while joining the gap needs the MAPQ 3-4 reads as anchors without letting them
own a haplotype call.

## Files

- `screen.py` / `unresolved_gaps.tsv` — depth-vs-floor screen over all unresolved gaps.
- `verify_lowmapq.py` / `lowmapq_verification.tsv` — truth scoring of the dropped reads per gap.
- `mapq_floor_sweep.tsv`, `gap_mapq_histogram.tsv` — the single-gap sweep and its read population.

## The right amount of signal: an asymmetric floor (`asymmetric_floor.tsv`)

Chasing HiPhase's site list is the wrong target -- it phases DeepVariant calls on
a linear BAM. (Though the MAPQ values are not the difference: for all 296 reads in
this window our surjected BAM and the linear BAM agree exactly, 207 reads at MAPQ
3-4 in both, so the mapping ambiguity is real and shared, and HiPhase's 89 tagged
reads are exactly its MAPQ >= 5 population.) The question is how much of our own
signal to admit, and for which job.

Tested by post-hoc filtering the `-q 1` output so that only reads at or above an
assignment floor keep HP/PS -- an estimate of a split floor without touching the
caller:

| configuration | phased reads | blocks | discordant | window Hamming | N50 |
|---|---:|---:|---:|---:|---:|
| single floor 30 (current default) | 395 | 2 | 0 | 0.00% | 166,590 |
| single floor 10 | 469 | 2 | 0 | 0.00% | 170,529 |
| single floor 5 (HiPhase's floor) | 494 | 2 | 0 | 0.00% | 176,361 |
| single floor 1 | 734 | 1 | 37 | 5.04% | 375,181 |
| discover 1 / assign 5 | 505 | **1** | **0** | 0.00% | 275,646 |
| discover 1 / assign 10 | 479 | **1** | **0** | 0.00% | 274,608 |
| discover 1 / assign 30 | 399 | **1** | **0** | 0.00% | 270,669 |

Letting MAPQ 3-4 reads inform candidate discovery and block linking, while
withholding haplotype tags from them, joins the gap with zero discordant reads and
lifts window N50 from 167 kb to 271-276 kb. At `assign 30` the contiguity gain
arrives without a single low-MAPQ read receiving a tag; relaxing the assignment
floor to 5 adds 110 more phased reads, still at zero error.

Caveat: this is a post-hoc tag filter, so the low-MAPQ reads did participate in
the k-means clustering and in the flank-link votes that produced the single block.
A real implementation has to decide deliberately whether they vote on orientation
or only supply candidate sites, and it needs validating on the other two
trustworthy gaps and then whole-chr20 before any default moves.
