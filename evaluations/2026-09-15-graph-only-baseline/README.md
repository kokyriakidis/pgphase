# Graph-only chr20 baseline and its phase-block gap inventory

Pass 1 of the two-pass design: `collect-graph-variation` alone, no BAM channel.
`run.sh` reproduces the run and the truth evaluation; `gap_inventory.py` derives
the gaps a follow-up subprocess would have to close.

41 s wall, 2 m 07 s CPU on 8 threads.

## Sites

73,627 candidates retained of 988,602 catalog sites considered:

| category | sites |
|---|---:|
| `CLEAN_HET_SNP` | 54,241 |
| `REP_HET_INDEL` | 17,350 |
| `CLEAN_HET_INDEL` | 2,036 |

914,975 filtered: `ref_only` 826,336, `high_af` 42,750, `no_reads_in_chunk`
35,481, `low_af` 7,044, `low_depth` 3,364. So 83.6% of the catalog is homozygous
reference in this sample and a further 4.3% homozygous alt (87.9% together); the
phasing-relevant catalog is 7.4% of it.

## Reads and accuracy against the diplinator truth

| metric | value |
|---|---:|
| reads phased | 203,751 |
| reads unphased | 27,631 |
| concordant | 201,991 |
| discordant | 1,744 |
| read Hamming | 0.856% |
| phase sets | 371 |
| block N50 | 917,428 bp |
| block auN | 987,490 bp |

`switch_errors` and `flip_errors` both report 0 and are **not measurable** for
this path: the emitted BAM carries no contig header, so the evaluator cannot
compute them. Only the read-level figures above are meaningful.

## Gap inventory

288 of the 375 phase sets carry two or more sites; those blocks span 47.60 Mb,
and the spaces between consecutive blocks give **286 gaps spanning 18.52 Mb**
(47.60 + 18.52 accounts for chr20's 66.1 Mb).

| gap size | gaps | span | untagged reads |
|---|---:|---:|---:|
| < 1 kb | 9 | 0.00 Mb | 11 |
| 1-10 kb | 6 | 0.03 Mb | 339 |
| 10-50 kb | 183 | 5.60 Mb | 12,478 |
| 50-200 kb | 80 | 6.50 Mb | 17,324 |
| > 200 kb | 8 | 6.39 Mb | 25,795 |

**55,583 distinct reads sit in gap windows without a haplotype assignment,
40,214 of them at MAPQ >= 30.** (Summing the per-gap column gives 55,947, since a
read can overlap two gaps when the block between them is shorter than a read.)

Two features worth noting for pass 2. The nine sub-kilobase gaps are 3-249 bp
wide with essentially every read already tagged (0-3 untagged each) -- adjacent
blocks that failed to link despite sharing reads, i.e. pure stitching failures
and the cheapest possible first test. And the two largest gaps fail for opposite
reasons: `chr20:43,688,351-45,896,820` (2.21 Mb) has 8,659 untagged reads of
which 8,425 pass MAPQ 30 -- confidently mapped reads left unphased -- while
`chr20:27,133,886-29,068,099` (1.93 Mb) has 9,825 untagged of which only 196
pass, the pericentromeric mapping-ambiguity case.

## Tooling note

`gap_inventory.py` first read tagged reads from the phased BAM by region, which
silently returned nothing for every gap (no contig header, no index, `samtools
view` exit 1) and reported all 94,014 overlapping reads as untagged. It now reads
`--phase-reads` (the per-read assignment TSV, authoritative) and raises on any
non-zero `samtools` exit.

## Comparing with the hybrid path

These two gap populations are not the same measurement. The hybrid
gap-recovery report lists 276 gaps *before* recovery and joins 112, leaving
**164 unresolved spanning 15.2 Mb**. The 286 gaps / 18.52 Mb above are
pre-recovery with nothing attempted, so the like-for-like comparison is 286
against hybrid's 276 total (4% larger, as expected without the BAM channel's
sites); against hybrid's post-recovery residue of 164 gaps it is ~75% larger.

## MAPQ floor sweep (`mapq_sweep.sh`, `compare_arms.py`)

Whole chr20, graph-only, at four floors. 55 s per arm including evaluation.

| arm | sites | phased | discordant | read Hamming | phase sets | N50 | gaps | gap span | gate |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| `-q 30` (default) | 73,627 | 203,751 | 1,744 | 0.856% | 371 | 917,428 | 286 | 18.52 Mb | — |
| `-q 10` | 76,374 | 209,963 | 2,645 | 1.260% | 381 | 920,605 | 309 | 18.06 Mb | 780 |
| `-q 5` | 77,012 | 211,255 | 2,688 | 1.272% | 378 | 920,605 | 301 | 17.96 Mb | 316 |
| `-q 1` | 79,291 | 214,065 | 2,945 | 1.376% | 379 | 1,011,696 | 311 | 17.81 Mb | 103 |

"gate" counts reads concordant at `-q 30` that become discordant in the arm.
Clean het SNPs rise from 54,241 to 59,282 between `-q 30` and `-q 1`; repeat het
indels barely move (17,350 to 17,830).

### The extra error is entirely in the newly admitted reads

Splitting the `-q 1` arm by each read's input mapping quality:

| input MAPQ | phased | discordant | error |
|---|---:|---:|---:|
| 1-4 | 2,769 | 457 | 16.50% |
| 5-9 | 1,250 | 285 | 22.80% |
| 10-29 | 6,210 | 462 | 7.44% |
| 30-59 | 13,986 | 225 | 1.61% |
| 60 | 189,850 | 1,516 | 0.80% |

Restricting both arms to reads at MAPQ >= 30 -- the population `-q 30` could
already see -- gives **0.856% at `-q 30` against 0.854% at `-q 1`** (1,744 vs
1,741 discordant on 203,751 vs 203,836 reads). So lowering the floor does not
degrade the phasing that already existed. The 10,347 reads newly phased at
`-q 1` carry 1,240 discordant calls between them (11.98%), and they are almost
all the low-MAPQ population (2,769 at MAPQ 1-4, 1,250 at 5-9, 6,210 at 10-29).

The 103 gate violations sit at MAPQ 30-59 (29) and 60 (74), against 138 reads
corrected in the other direction. With the MAPQ >= 30 error rate flat to three
decimal places, those flips are consistent with blocks re-forming -- phase sets
go 371 to 379 -- rather than with systematic corruption. The
gate count being non-monotonic across floors (780 at `-q 10`, 316 at `-q 5`, 103
at `-q 1`) points the same way.

### Consequence for the two-pass design

`-q 1` is the right floor for pass 1 on the evidence above: 5,041 more clean het
SNPs, 0.74 Mb more phased span, 0.71 Mb less gap span, and no cost to the reads
already phased. (An earlier version of this sentence claimed "10% more contiguity"
from the evaluator's N50; that metric is invalid on this path and the claim is
withdrawn -- see the correction section below, where site-derived N50 falls 2.5%.) What it should not do is *tag* the 10,347 ambiguous reads at
11.98% error -- under the two-pass design those reads are pass 2's job, to be
decided from local gap evidence rather than from a chromosome-wide floor.

That separation already exists as `--min-assign-mapq` (added earlier today) but
is **not reachable from this subcommand**: only `hybrid_collect.cpp` parses the
option, and `read_carries_phase_tags` is consulted only in
`collect_bam_output.cpp:458`. The graph-only path writes its phased BAM through
its own mirror of that writer (`graph_bam_adapter.cpp:209-222`), which never
consults it. Wiring it up is two small additions: expose the option on the graph
subcommand, and consult the predicate in the graph BAM writer.

## Correction: the evaluator's block-size metrics are invalid for this path

The `N50` column in the sweep table and the auN column in `graph_mapq_sweep.tsv`
come from `scripts/evaluate_phase_accuracy.py`, which derives block extents from
the emitted alignment. That alignment is unaligned on this path, and the numbers
are not trustworthy: the evaluator reports **auN = 15,474,231 bp at `-q 1`**
while the largest single phase block is **1,282,456 bp**, and auN cannot exceed
the largest block. Its N50 for the same arm (1,011,696 bp) is likewise more than
double the site-derived value.

Recomputing block spans directly from `phase_sites.tsv` (first to last phased
site per phase set):

| arm | blocks | phased span | N50 | auN | largest block |
|---|---:|---:|---:|---:|---:|
| `-q 30` | 375 | 47.60 Mb | 439,677 | 496,788 | 1,282,456 |
| `-q 1` | 402 | 48.34 Mb | **428,655** | 492,205 | 1,282,456 |

So block structure is essentially flat across the floor change -- N50 **falls
2.5%** and auN falls 0.9%, with an identical largest block -- and the earlier
"N50 rises 10%" reading was an artifact of the evaluator, not a real gain. It is
withdrawn.

What survives unchanged as the case for `-q 1`, because none of it depends on the
evaluator's block metrics: 5,041 more clean het SNPs, the MAPQ >= 30 population's
accuracy flat at 0.856% -> 0.854%, phased span up 0.74 Mb, and gap span down from
18.52 Mb to 17.81 Mb. The decision stands on those; the block-length argument
does not.
