# Long repeat anchors close the next qualifying HiPhase block

The next qualifying uncovered HiPhase block is **chr20:14,719,378–15,335,938
(616,561 bp)**, ranked 36th. Every larger remaining seam has HiPhase below
80% correct on all original truth-scorable primary overlapping reads, including
abstentions. `ranked-blocks.json` and `screened-seams.json` preserve the task-start
ranking and the identical-read screen.

## Cause and repair

The block has two production seams: **15,023,123–15,039,543 (16,420 bp)** and
**15,100,456–15,101,262 (806 bp)**. Both tools already have the relevant alleles.

At the first seam, a noisy recalled one-base deletion occurs on both parents.
The reliable MSA source is **15,019,255 TACACAC>T**, before that boundary.
Its physical six-base deletion is left aligned at 15,019,256, while the original
CIGAR places it at 15,019,295: a 39-base equivalent shift beyond the old
32-base placement search. The boundary-only nomination missed the earlier
source, and a local alignment window failed to cover the full AC run.

Recovery now searches verified source deletions of 4–32 bases within 25 kb of
the downstream SNP. It normalizes only a physical copy through at most 64
equivalent positions, preserving the original candidate. Full repeat REF/ALT
sequences and external 16-base Q20 anchors identify the allele, with at most
two bases of net-length slippage. Independent original MAPQ30 molecules ending
before the downstream SNP calibrate the source against an upstream clean SNP;
separate molecules bridge both classes to the downstream SNP. Existing
calibrated-indel confidence bounds, unanimous parity, and wrong-parity odds
at most 0.001 are required. An existing agreeing source rescue can retain Q10
anchors with at most 20% call error; it never calibrates a join or gains a new
haplotype from that weaker evidence.

The second seam uses **15,101,262 A>ATG / A>ATGTGTGTG**, a two-base versus
eight-base tandem insertion. Recovery accepted only one versus two motif copies,
and required SNP gauges within 10 kb. The downstream SNP is 14,125 bases away.
The TG reference run is 53 bases long and shifts insertion placement to its far
end, outside the old local flank window.

Recovery now accepts larger whole-copy two-base motifs and gauges within 20 kb.
Only classes separated by more than one motif copy extend alignment across the
full bounded reference run; the established two/four-base path retains its
local window. Complete sequence and signed nearest length must agree, with
bounded negative slippage. Clean phased SNPs still veto assignments. Contrary
noisy recalled rows cannot erase a certified separated-class call.

The overlapping **15,095,642–15,101,261** control also joins. Reaching HiPhase's
read count required applying the certified source flip before retaining a
separated-tandem rescue, and physically reconstructing the shifted verified
one-base insertion **15,109,300 C>CA**. Two independent diploid cohorts confirm
that existing insertion's gauge against the downstream SNP. Each cohort must
pass p<=0.01 association and both agree; pooled 95% discordance bound plus call
error is at most 20%. The stricter terminal insertion predicate is unchanged.
Each assigned read needs its own Q20 physical call with at most 1% error.
Physical SNPs, certified complementary tandem calls and clean SNP profiles veto
contrary corrections. The latter veto preserves the correctly assigned short
repeat molecule when its one-base insertion call is noisy.

`audit_physical.py` reconstructs certificates from original sequence, quality,
CIGAR and mapping quality without parental truth. The deletion calibration is
**[[22,0],[0,27]]** with two diploid agreeing bridges and log odds **−16.209563**.
The insertion calibration is **[[8,0],[0,6]]**, with twelve agreeing reverse
bridges and log odds **85.255646**. The retained one-base insertion cohorts are
**[[1,8],[8,0]]** and **[[0,8],[8,0]]**; their joint error bound is **13.4998%**.
`physical-certificates.json` records every molecule and its role.

## Results

| Original overlapping primary reads | Before correct/core | Fixed correct/core | HiPhase correct/core |
|---|---:|---:|---:|
| 16,420 bp seam; 126 scorable | 99/58 | **101/101** | 101/101 |
| 806 bp seam; 73 scorable | 44/18 | **70/70** | 68/68 |
| Overlapping control; 101 scorable | 69/38 | **95/95** | 95/95 |

Fixed correctness is **80.16%, 95.89%, and 94.06%** respectively, with
unphased reads in each denominator. The repaired seams and the newly connected
control meet or exceed HiPhase in total-correct and connected-core reads.

The full production core now extends **14,719,378–15,569,349 (849,972 bp)**,
covering the complete HiPhase target. Across the entire target, 2,813 original
reads score: pgphase has **2,711 correct, eight discordant, 94 unphased, and
2,694 correct in the connected core**; HiPhase has **2,728 correct, eight
discordant, 77 unphased, all 2,728 correct in its core**. Thus whole-block read
parity remains **17 total-correct / 34 core-correct reads short**, outside the
per-seam closure certificates. Structural closure does not establish whole-block
read parity.

The full chromosome retains every previous block extent and has **251 blocks**,
**N50 955,496 bp**, and largest block **3,012,193 bp**. All 82,287 variant keys,
alleles, counts and categories are preserved. Original-read correctness changes
**230,783→230,811**, discordant **6,646→6,644**, unphased **19,181→19,155**.
All previous correct and phased labels are preserved; all status changes occur
in the affected neighborhood. One previously unphased read is discordant after
assignment, alongside 25 newly correct and three corrected reads.

The 15–16 Mb owner has **4,300 correct, 14 discordant, 230 unphased** versus
4,272/16/256 before; one/four-thread outputs are identical. The 14–16 Mb
continuation has **7,429/32/1,342** versus 7,401/34/1,368. Both preserve all
previous correct and phased labels and variant evidence.

## Verification

The committed panel contains both production seams; the existing overlapping
control now has an intentionally changed exact `spans=1` expectation. All three
are certified against the measured HiPhase total-correct and connected-core
counts, with the 80% denominator including unphased reads. The control retains
its original short replay. The two production seams retain the full 15–16 Mb
owner, and the target endpoints retain the 14–16 Mb continuation. The native
fixture checks parental orientation on independent sides, exact denominators,
source alleles/counts, and the affected rescue/correction/veto reads.

The original failing binary reproduced both production seams in the new
regression. Unit fixtures test actual source cohorts, opposite gauges, missing
allele classes, sparse association, and confidence-bound rejection. Earlier
window floors stay intact. Full-chromosome and every previous native replay
are audited for lost correct/phased labels and changed candidate evidence.
Measurements and validation totals follow in `results.json`,
`owner-results.json`, `continuation-results.json`, `full-preservation.json`,
`panel-audit.json` and `test-runtime.json`.

Reproduction uses `replay.sh`, then the `audit_*.py` scripts with the local
bench-phasers Python and pysam. `make unit-tests predicate-tests check` and
`make window-tests` run the validation gates.

Final validation: `make -j10` builds without errors or warnings; all unit tests
pass, phasing predicates pass **1,551 assertions / 47 cases**, and HiFi/ONT TSV
and phased VCF goldens match. The final gap suite passes **17,534 assertions /
four cases**, including replay-cache self-tests. The new fixture passes **751
assertions**; the task-start binary fails **nine assertions**, reproducing the
bug. The earlier 0.865 Mb fixture retains its original floors and passes **682
assertions**. The HiPhase comparator self-tests pass and all **131 panel windows**
are measured on identical alignments. Cached final runtimes are **1.12 seconds
for the focused regression** and **40.33 seconds for the full suite**; cold
pipeline execution is separate from those cached times.

The previous replay audit passes **238 output labels / 126 independent
requests**, all checked against the final production binary. Four labels
change: the affected short controls and their owning/stitched contexts. No previous correct or phased label, block extent, variant key,
or candidate evidence is lost (`panel-audit.json`).
