# Closing the seam inside HiPhase’s largest chr20 block

HiPhase’s largest phased-heterozygote VCF block is chr20:63,182,011–66,206,480,
PS 63182011, 3,024,470 bp and 5,820 rows. The task-start pgphase output splits
this region into 63,182,011–65,509,355 and 65,509,406–66,194,203. The target is
therefore the 51 bp seam at 65,509,355–65,509,406.

The seam lies inside SDUST interval 65,507,467–65,514,083. Near-boundary graph
SNP placements disagree with physical calls and with the established read
haplotypes. A nearest-SNP bridge alone would not validate those gauges. The
existing whole-graph path check also fails at an earlier one-haplotype edge,
65,152,522–65,168,109. The new route uses physical SNPs outside the repeat,
calibrated on independent reads assigned to each existing phase set. It keeps
all SNP calls, candidate counts, alleles, filters and existing read assignments.

The native 65,000,001–66,000,000 replay supplies 144 calibration molecules on
the left (75/69 haplotypes) and 148 on the right (68/80), with no contrary
calibration. A separate primary MAPQ58 molecule spans the repeat and calls at
least two separated Q30 SNPs on each flank. Its conservative union error is
0.00040317; the joint calibration and molecule error bound is 0.05679768. Both calibration Wilson upper bounds, one percent call error per
flank, and bridge error together stay below 20%. The union is deferred until
read rescue finishes and applies to complete, consistently oriented phase sets.

The seam has six original primary alignments; every one is truth-scorable.
Both pgphase and HiPhase correctly phase all six when independently oriented
on the disjoint local flanks. Each tool has 145/145 correct left-flank and
146/146 correct right-flank reads. The full HiPhase PS changes parental
orientation near 65.2 Mb, and is mostly hap1-maternal before that switch but
hap1-paternal at the seam. Its whole-PS orientation therefore calls these six
reads discordant. `hiphase.json` and the benchmark manifest use that consistent
whole-PS scoring method; `orientation.json` preserves the additional local
comparison. The local comparison remains six of six, rather than pretending
that a global phase switch means HiPhase cannot bridge this seam.

Native preservation: 2,962 correct and six discordant assignments remain;
all 1,591 variant records retain their alleles, counts and filters. The union
changes 578 variant PS labels and 1,518 read PS labels, with no new assignment,
lost tag or lost correct read. The new owner regression fails five assertions
with the task-start binary. The pure calibration fixtures cover same and
reversed gauges, contrary molecules, missing bridges, weak qualities, missing
diploid calibration, reversed calibration, tiny cohorts and excessive error.

The HiPhase terminal marker at 66,206,480 is a T>A call with GQ2. All eleven
original primary alignments there have MAPQ3–15. Pgphase’s last emitted
heterozygote is 66,194,203. The follow-up [endpoint investigation](../2026-10-05-terminal-endpoint/README.md)
traces the actual cause: reference-only GAF observations, low-MAPQ alternate
observations at earlier tail sites, and skipped BAM realignment at the
chromosome end. Low GQ/MAPQ alone does not establish that extension is
unsupported. A MAPQ1 replay extends beyond HiPhase, but gets only 32/54
tail reads correct and loses 74 previously correct owner assignments.

Full-chromosome and panel preservation results follow in the JSON artifacts.

The frozen full chromosome output has one PS 63896021 across
63,182,011–66,194,203: 3,012,193 bp and 4,939 phased heterozygote rows. This
covers 99.594% of HiPhase’s raw largest-block span and retains the existing
endpoint difference of 12,277 bp. Chr20 N50 improves from 856,770 to 863,977 bp;
HiPhase N50 is 1,005,183 bp. Spans are inclusive first-to-last phased
heterozygote positions; singletons do not count as multi-site blocks.

Across all 9,843 original primary truth-scorable alignments overlapping
HiPhase’s largest interval, pgphase has 9,339 correct, 404 discordant and 100
unphased (94.8796% correct, abstentions included). Its dominant connected
correct core increases from 7,226 to 9,209. HiPhase has 6,681 correct, 2,915
discordant and 247 unphased (67.8756%), with a dominant correct core of 6,681.
All 9,843 HiPhase alignments match the input geometry, CIGAR and sequence.
Each tool’s phase sets are oriented independently on their whole output read
cohort; local orientation is reported separately rather than replacing this
consistent whole-block measurement.

Whole-chromosome preservation retains all 256,610 primary outputs,
230,567 correct assignments, 6,768 discordant assignments and 64,188 variant
records. The union changes 1,024 variant gauges and 2,324 read PS labels.
There is no new assignment, lost tag or lost correct read. All previous VCF
extents remain contained in the final blocks. The only new closure is
65,509,355–65,509,406 (`full-closures.json`).

The fast development commands are `make gap-dev-check` and
`make gap-owner-check GAP="largest HiPhase"`. The owner regression reruns its
546 assertions against saved native output; its corrected task-start replay
fails five assertions. The baseline and final production SHA256 values are
in `full-closures.json`.

Validation: the production build has zero new warnings; all unit tests,
47 predicate cases / 1,551 assertions and HiFi/ONT TSV/VCF golden gates pass,
including HiFi one/four-thread determinism. `make window-tests` passes 13,494
assertions in four cases and four replay-state helper tests. All 201 selected
checks pass across 115 panel windows and the owning/mechanism cases; 102
panel gaps span and 37 satisfy both the >=80% and HiPhase total/core contract.
All 335 matched native replays retain every old correct/phased assignment and
every variant allele, count and filter. The pure development loop takes
0.038 s; the cached new owner regression takes 0.508 s. The warm full suite
takes 35.47 s; the one required cold full-chromosome audit took 448.40 s while
running concurrently with the integration checks.
