# Close the eleventh HiPhase block’s leading source boundary

The selected gap is **11,796,979–11,813,446 (16,467 bp)** inside HiPhase’s
**11,796,969–12,955,838 (1,158,870 bp)** block, ranked eleventh by length.
The fresh ranking uses the preceding terminal-source-extension output.
Larger uncovered boundaries fail HiPhase’s 80% original-overlap rule, even
with their best local orientation; `selection.json` records that screen.

## Cause and repair

The catalog at 11,813,617 describes several repeat alternatives. Its selected
minimal **11,813,622 T>C** SNP has 5 ALT/11 graph observations, but original
MAPQ30 primary BAM bases contain **61 T, 8 C, no deletion or third base**.
The physical ALT fraction is 11.6%, below the 20% heterozygote threshold.
Decomposition appends an ALT suffix to the displayed site ID; physical
validation looked up that ID rather than the original catalog key, and its
single-candidate gate also skipped this selected SNP/deletion pair. The false
SNP boundary hid the BAM-recalled **T>TAC insertion** and downstream
**11,813,668 T>C** SNP.

HiPhase’s DeepVariant input already contains the correct T>TAC insertion
and the downstream T>C SNP (GQ 14 and 35 respectively). pgphase’s BAM MSA
also recalls them; its false catalog boundary prevented their use.
`physical-evidence.json` records the original CIGAR base counts and HiPhase rows.

Validation now uses the original catalog key. A padded mixed SNP/indel site
with one selected clean SNP and exactly one other selected ALT is retired after the
initial solve when the physical ALT fraction is below `--min-af`, REF is
observed, no deletion/third base is observed, and the corrected one-sided
balance test is <=0.01. Surviving gauges and observations are preserved.
The REF-absence test family retains its previous eligibility and correction.

Broader trials also applied this binary rule to single-candidate sites and
complex decompositions with additional selected branches. Full-chromosome
and panel checks caught lost tagged reads and a reopened existing gap near 7 Mb.
That trial was rejected; the final rule requires exactly two selected ALT rows in the mixed decomposition
that caused this bug. Existing regression floors were not loosened.

The leading BAM block contains complementary overlapping **14 bp and 4 bp
deletions**. Raw CIGAR placement can shorten those calls to 12 bp and 2 bp.
The final bridge compares both expected deletion haplotypes over one reference
window, using sequence edit distance and net deletion length together. Only
period-one/two repeat differences qualify; ties, reference-length calls and
conflicting sequence/length evidence abstain.

Four independent original molecules agree on the block connection. Separate
left-only and right-only cohorts calibrate both existing read gauges. The
existing calibrated repeat/SNP predicate enforces its 20% joint uncertainty
bound, with MAPQ/BQ20 calls and <=0.001 joint wrong-parity bound. The right BAM
SNP retains a cut-free source prefix to an exact shared graph SNP. This connects
an independently calibrated BAM-only block without requiring a graph path
inside it. Existing rescued cohorts join the connected core; physical deletion
calls correct conflicting assignments within that core. Truth is evaluation only.

The full chromosome audit exposed a second defect: applying this late union
only in its owning chunk split the pre-existing 904,351 bp continuation at
the next chunk boundary. `before-global-union.json` records the failing
block extents. The union now propagates its already-global phase-set IDs
and gauge through every chunk. The owning regression additionally replays
11–13 Mb and requires the left deletion, right boundary and 12,717,796 SNP
to remain in one phase set.

## Native results

The denominator includes every original primary truth-scorable gap overlap,
including unphased reads. Rescue phase sets do not count toward connected core.

| Tool | Correct | Discordant | Unphased | Scorable | Dominant correct core |
|---|---:|---:|---:|---:|---:|
| Before | 106 | 12 | 18 | 136 | 61 |
| Fixed | **126** | **2** | **8** | **136** | **126** |
| HiPhase | 126 | 2 | 8 | 136 | 126 |

The owner has 4,589 primary reads: **4,004→4,024 correct**, **90→80 discordant**,
**495→485 unphased**. Twenty-one reads become correct. One previously correct
read becomes discordant inside this gap, and one unphased read becomes discordant.
The final two discordant gap reads are the same two discordant reads as HiPhase:
`m84031_231217_062403_s3/74975318/ccs` and
`m84031_231217_034919_s2/20453286/ccs`. Both physically carry the opposite
SNP base to their parental truth; production does not use truth to override them.
Every previously tagged read remains tagged, and unrelated native allele
counts/filters are unchanged. One- and four-thread candidates, VCF records and
BAM assignments match (`owner-results.json`).

Disjoint cohorts that observe the left and right markers independently verify
parental orientation. The new panel row, certified manifest, exact `spans=1`
expectation and native 11–12 Mb plus 11–13 Mb owning regression enforce total correctness
and connected-core parity. The saved old binary fails 12 assertions in the new
regression. `make gap-benchmark` independently regenerated HiPhase’s new
136/126/126 row on identical original alignments and preserved all prior rows.

## Full chromosome and preservation

The joined block spans **11,792,550–12,717,796 (925,247 bp)**. Every previous
block extent is retained. Block count falls **258→257**, **N50 stays 944,186 bp**,
and the largest block stays **3,012,193 bp**.

Across all 256,610 primary reads, correct **230,666→230,686**, discordant
**6,692→6,682**, unphased **19,252→19,242**. Twenty-one reads become correct;
the one previously correct read that becomes discordant is inside this gap
and is discordant in HiPhase too. Every previously phased read stays phased.
All remaining old variant allele counts and filters are unchanged
(`full-preservation.json`). The global union changes 3,960 read tags, mostly
phase-set names; read status changes are confined to the repaired gap.

The same physical rule withdraws seven other under-supported catalog SNPs:
2,282,198 T>A; 42,290,009 A>G; 47,003,854 A>G; 50,431,256 A>T;
62,940,433 A>G; 63,734,075 C>T; and 63,793,665 C>T. Original primary MAPQ30
base counts are recorded in `physical-evidence.json`. Recovery also exposes
47,003,826 G>GAGAA and 50,880,007 T>TA. Those changes preserve the existing
read assignments and block extents. The full VCF has 64,479 records.

Over the entire 1,158,870 bp HiPhase block, pgphase has **5,004 correct/core
4,003** versus HiPhase **5,049 correct/core 5,049**, among 5,177 original
truth-scorable primary overlaps. The remaining total/core deficits are 45 and
1,046 reads, associated with the separate boundaries below. No additional
whole HiPhase block becomes covered by this repair.

The optimized build has zero new warnings. Units, all 1,551 predicate assertions,
HiFi/ONT golden outputs and determinism pass. The full suite passes **16,794
assertions in four cases**, plus four cache-helper tests. The owning single-
and two-chunk regression runs **1,131 assertions in 0.81 s** from saved states.
`panel-audit.json` verifies all **225 previous labels / 118 current requests**
under the final binary; no additional previously correct read is lost.

Final optimized binary SHA256:
`8c58ce2ef080d20d37aa7593c6abb259522236095c08507f171e99475ec56f98`.

## Remaining boundaries

This closes the leading 16,467 bp seam. Two separate boundaries in the same
HiPhase block remain open: **12,717,796–12,740,002 (22,206 bp)** and the
**12,954,878–12,955,838 (960 bp)** endpoint. `selection.json` records their
read evidence. The complete HiPhase block remains split at those separate boundaries.

## Reproduction

```sh
make -j8
make gap-dev-check
make gap-owner-check GAP=repeat-snp-validation
make unit-tests predicate-tests check
LD_LIBRARY_PATH=/usr/lib/x86_64-linux-gnu make window-tests
bash evaluations/2026-10-06-leading-source-block/replay.sh
/home/kokyriakidis/micromamba/envs/bench-phasers/bin/python evaluations/2026-10-06-leading-source-block/audit_block.py
/home/kokyriakidis/micromamba/envs/bench-phasers/bin/python evaluations/2026-10-06-leading-source-block/audit_physical.py
/home/kokyriakidis/micromamba/envs/bench-phasers/bin/python evaluations/2026-10-06-leading-source-block/audit_owner.py
/home/kokyriakidis/micromamba/envs/bench-phasers/bin/python evaluations/2026-10-06-leading-source-block/audit_preservation.py
```
