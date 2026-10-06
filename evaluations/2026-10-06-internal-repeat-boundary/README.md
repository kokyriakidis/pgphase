# Recover the eleventh HiPhase block's internal repeat chain

The selected gap is **12,717,796–12,740,002 (22,206 bp)** inside HiPhase's
**11,796,969–12,955,838 (1,158,870 bp)** block, ranked eleventh by length.
The ranking uses the preceding leading-source-boundary fix. All larger
remaining seams fail HiPhase's 80% original-read correctness condition;
`selection.json` records those measurements.

pgphase finds the **12,721,112 A>AT** repeat insertion but excludes it from
clean phasing. The graph deletion carries **26 REF/49 ALT observations**;
its equivalent MSA **12,735,894 TA>T** call collapses to HOM with only eight
ALT observations. Consequently neither can connect the SNP desert. The
short recalled insertion at 12,740,002 also confuses the local rescue gauge.
HiPhase's DeepVariant input has both true repeat heterozygotes (GQ 35 and 4)
and calls the downstream insertion homozygous (GQ 6).

The repair reconstructs this insertion/deletion chain using original primary
alignments, without truth. It requires a pure bounded homopolymer run, matched
Q20 external flanks, at most one physical edit and one base of length slippage.
Current core reads independently calibrate the left insertion, with Q20 edits
and agreeing graph calls. Seven distinct molecules cross both repeat markers:
five REF/REF and two ALT/ALT, with no contrary pair. Separate right molecules
observe the deletion plus at least two consistent, separated Q20 clean SNPs;
five agree and one disagrees, exceeding the 80% agreement rule. The existing
left calibration uses its association test and a 20% Wilson/call-error bound;
the right side uses observed agreement and quality likelihood, not a 95%
accuracy claim. The predicate tests cover contrary bridges, single-class links,
weak calibration, weak quality, sparse right evidence and agreement below 80%.

The collapsed deletion becomes a phased heterozygous nonanchor with physical
REF/ALT observations and consistent depth/strand counts. Existing core read
assignments are preserved. Physically validated unassigned/rescued molecules
enter the joined core. The union propagates across already stitched chunks.
The two-chunk regression caught stale repeat PS labels: unlike clean anchors,
excluded repeat rows retain their old labels. This pass discovers seams from
surviving non-repeat anchors and infers the marker gauge from current reads.

| Tool | Correct | Discordant | Unphased | Original scorable | Correct connected core |
|---|---:|---:|---:|---:|---:|
| Before | 112 | 6 | 60 | 178 | 74 |
| Fixed | **164** | **0** | **14** | **178** | **164** |
| HiPhase | 164 | 0 | 14 | 178 | 164 |

Abstentions remain in the denominator; output-only rescue phase sets do not
count toward the connected core. Correctness is **92.13%**, with total and
connected-core parity against HiPhase on identical original alignments.
The committed panel, exact `spans=1` expectation, measured HiPhase benchmark,
required deletion row, certified manifest and owning/two-chunk regression
record this closure. The saved pre-fix binary fails the new regression.

The native owner changes from 4,147 correct/21 discordant/108 unphased to
4,199/15/62 over 4,276 primary reads. One- and four-thread candidate tables,
VCF records and BAM tags are identical. The two-chunk regression preserves the
preceding leading-gap fix and requires its deletion, the left SNP and the
right SNP to remain in one phase set. Its cached 613 assertions take **0.61 s**.

The full chromosome has **230,738 correct**, **6,676 discordant** and
**19,196 unphased** among 256,610 primary reads: 46 abstentions and six
previously discordant reads become correct. No previously correct or phased
assignment is lost; read-status changes are confined to the selected gap.
All previous variant alleles, counts and filters are preserved, and only
12,735,894 TA>T is added (`full-preservation.json`). There are 64,480 VCF rows.

Every previous block extent is retained. Block count falls **257→256**;
**N50 rises 944,186→944,265 bp**; the largest block remains **3,012,193 bp**.
The joined block spans **11,792,550–12,954,878 (1,162,329 bp)**. Across the
whole selected HiPhase block, pgphase has **5,056 correct/core 5,030**, compared
with HiPhase **5,049/core 5,049** over 5,177 original scorable reads. This gap
has full parity, while the whole block remains 19 connected-correct reads
short and has a separate **12,954,878–12,955,838 (960 bp)** terminal gap. No
additional whole HiPhase span becomes covered (`results.json`).

The optimized build, units, all 1,551 predicate assertions, new bridge fixtures,
HiFi/ONT goldens and determinism, independent HiPhase benchmark helper tests,
and cache helper tests pass. The full suite passes **16,878 assertions in four cases**, plus four cache
helper tests. All **227 previous native output labels (119 independent
requests)** are audited under the final binary and runtime; only the two
selected owner/continuation labels change, with no lost correct assignments
(`panel-audit.json`). All existing required-site rows and read floors are
preserved.
