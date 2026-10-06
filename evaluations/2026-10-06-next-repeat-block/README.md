# Compound insertion prefix closes the next qualifying HiPhase block

The next qualifying open gap is **chr20:45,876,157–45,896,820 (20,663 bp)**,
inside HiPhase's **45,866,904–46,636,707 (769,804 bp)** block, ranked 30th.
The fresh ranking uses the preceding terminal insertion/deletion repair.
Every larger remaining seam fails HiPhase's 80% correctness rule on all
original primary truth-scorable overlapping alignments. `selection.json`,
`ranked-blocks.json` and `screened-seams.json` preserve the screen.

## Cause and repair

Both tools have the compound repeat alleles at **45,866,904 G>GACAGACAGACACACACAC /
G>GACAGACAGACACACACACAC**, the intermediate **45,883,706 C>CA** marker and the
**45,896,820 A>G** SNP. HiPhase's DeepVariant input marks the intervening
**45,876,157 G>GT** insertion homozygous, while pgphase's MSA source calls it
heterozygous. Original primary MAPQ30 molecules carry that one-base T insertion
on both parental haplotypes. The compound insertion lengths also overlap
between haplotypes, so that source's read assignments cannot replace a
physical bridge certificate.

The late repeat bridge only recognizes complete period-one/two insertion
motifs and requires a clean SNP at the other boundary. The two compound
insertions share a nonrepeat prefix followed by a differing AC copy, and the
reliable intermediate marker is a graph homopolymer insertion. A second defect
charges the low-quality shared inserted prefix against phasing confidence,
even though it is present in both alleles. That rejects the original short-allele
bridge with prefix Q10 bases despite its Q40 flanks and REF marker call.

The repair verifies both complete insertion sequences in 16-base flanks,
requires sequence distance and net length to agree, limits compound length
slippage to two bases, and rejects absence of the compound insertion. Only
distinguishing inserted suffix bases contribute class error; shared prefix
bases still enter the complete sequence comparison. Unknown qualities abstain.
Original primary nonduplicate MAPQ30 molecules need Q20 anchors, Q10
distinguishing inserted bases and at most 1% physical base/mapping error per
marker.

Separate nonbridging cohorts calibrate the left source and right marker gauges.
The left uses the existing calibrated indel predicate; the right uses the
existing diploid repeat-insertion predicate. Both predicates must independently
accept the same parity. Both graph-marker allele classes need physical bridges,
and every bridge must agree. After the globally propagated union, each original
Q20 marker call with at most 1% total error can assign or correct its own read;
phased clean SNP conflicts veto an assignment. Candidate alleles, counts,
categories and variant keys are preserved. Parental truth is evaluation only.

`audit_physical.py` independently reconstructs the production evidence without
parental truth: left matrix **[[27,0],[3,14]]**, right **[[15,1],[1,4]]**, and four
agreeing bridges, **one REF / three ALT**. The log odds are **28.629897**;
left joint error bound **16.8816%** and right diploid joint bound **14.2523%**
both pass their existing 20% predicates. `physical-evidence.json` lists the
molecules, original calls, qualities and error products. No zero-length compound
call enters calibration.

## Native results

| Original gap reads, including abstentions | Before | Fixed | HiPhase |
|---|---:|---:|---:|
| Primary truth-scorable | 141 | 141 | 141 |
| Correct | 111 (78.7%) | **129 (91.5%)** | 128 (90.8%) |
| Discordant | 20 | 8 | 4 |
| Unphased | 10 | 4 | 9 |
| Correct in one connected core | 72 | **129** | 128 |

The repair exceeds HiPhase's total-correct and connected-core counts. It retains
four more discordant reads while abstaining on five fewer reads. Two already
correct right-source rescues with original inserted-base Q22/Q27 enter the core
through their own matching physical A insertion; Q30 calibration remains
unchanged.

The 45–46 Mb owner retains 4,050 primary names: **1,559→1,577 correct**,
**797→785 discordant**, **1,694→1,688 unphased**. All 12 corrected and six newly
phased reads overlap the gap. One- and four-thread candidates, VCF records and
BAM assignments match. All existing variant evidence and all previously correct
or phased assignments are preserved (`owner-results.json`).

The 45–47 Mb continuation retains 8,033 primary names: **5,405→5,423 correct**,
**814→802 discordant**, **1,814→1,808 unphased**, with the compound pair,
boundary SNP and 46,636,707 SNP in one phase set. A short padded continuation
lacks the original BAM phase-set context, so both owner and full continuation
regions are recorded in `src/test_gap_replays.tsv`. Original bridge molecules
and a disjoint downstream SNP cohort verify the same parental orientation.
The noisier compound-only cohort is reported separately in the full audit.

The committed panel, certified closure manifest, identical-alignment HiPhase
measurement, exact span expectation, source/SNP keys, native read bounds and
owning/continuation regression all cover this gap. The old binary fails
**13 assertions** of the new regression; the repair passes **677 assertions**.

## Reproduction

```bash
export LD_LIBRARY_PATH=/usr/lib/x86_64-linux-gnu
make -j8
make unit-tests gap-dev-check check
make gap-owner-check GAP=45.876
make window-tests
bash evaluations/2026-10-06-next-repeat-block/replay.sh
/home/kokyriakidis/micromamba/envs/bench-phasers/bin/python evaluations/2026-10-06-next-repeat-block/audit_physical.py
/home/kokyriakidis/micromamba/envs/bench-phasers/bin/python evaluations/2026-10-06-next-repeat-block/audit_owner.py
/home/kokyriakidis/micromamba/envs/bench-phasers/bin/python evaluations/2026-10-06-next-repeat-block/audit_continuation.py
/home/kokyriakidis/micromamba/envs/bench-phasers/bin/python evaluations/2026-10-06-next-repeat-block/audit_block.py
/home/kokyriakidis/micromamba/envs/bench-phasers/bin/python evaluations/2026-10-06-next-repeat-block/audit_preservation.py
```

## Full chromosome and preservation

The joined core spans **45,854,428–46,636,707 (782,280 bp)** with 575 phased
heterozygotes, covering the complete 769,804 bp HiPhase target. Whole chr20
block count falls **254→253**. **N50 remains 955,496 bp** and the largest block
remains **3,012,193 bp**.

Among 256,610 primary names, correct assignments rise **230,765→230,783**,
discordant fall **6,658→6,646**, and abstentions fall **19,187→19,181**. Every
previously correct and phased assignment, variant key, allele count/filter and
block extent is preserved. All 18 status improvements overlap this gap; outside
it only a consistent joined phase-set gauge changes. `results.json` and
`full-preservation.json` record original-alignment identity and every transition.

Across all 3,129 original truth-scorable reads overlapping the whole HiPhase
block, pgphase has **3,035 correct / 26 discordant / 68 unphased**, versus
HiPhase's **3,042 / 19 / 68**. The separate residual is **7 correct reads /
41 connected-core correct reads**. The noisier compound-only cohort remains
**36 correct / 15 discordant**; completing the block geometry and satisfying the
repaired gap's parity contract does not claim that whole-block deficit is fixed.

All **17,348 gap assertions in four cases**, four cache-helper tests, unit
checks, **1,551 predicate assertions**, and HiFi/ONT golden/determinism gates
pass. HiPhase independently measures **129 committed windows** on identical
original alignments. The preservation audit checks all **236 prior native
labels / 124 independent requests** with the final binary; none change
(`panel-audit.json`). Prior committed fixture rows remain unchanged apart from
additive TOTAL counts (`committed-fixture-preservation.json`).

The focused saved-state regression takes **0.60 seconds**; the complete suite
and cache helpers take **39.65 seconds** (`test-runtime.json`). Final production
SHA256: `0d64a8ab810cd3e6dbfddda6471abe57dbd5b95b416bcccf92b99e72d40e0688`.
