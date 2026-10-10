# Complete complex tandem calls at the 19.373 Mb overlapping source

HiPhase's rank-18 block is chr20:18,442,592–19,464,219 (1,021,628 bp).
The baseline is the previous accepted source-provenance repair, binary
`10ee7e7a82564a4877e99783a339c9b4b999b1b0c5c5f78fc6975946fa5f8712`,
full output `test_data/tmp_gap_fix87/final/0`.
The fresh audit checks adjacent phased-site PS transitions, including overlaps,
in all of the 50 largest HiPhase blocks. This is the largest remaining
connection with HiPhase correctness >=80% where pgphase trails its total or
connected-core correct count. Larger blocks have no qualifying deficit at
their remaining transitions (`ranked-blocks.json`, `screened-seams.json`).

| Interval | Original reads | Before correct / core | After correct / core | HiPhase correct / core |
|---|---:|---:|---:|---:|
| 19,365,840–19,373,922 | 89 | 66 / 60 | 87 / 86 | 86 / 86 |
| 19,373,922–19,377,345 | 67 | 39 / 33 | 60 / 58 | 54 / 54 |

Correct fractions include unphased original primary reads: 97.75% and 89.55%.
The first interval retains one discordant and one unphased read; the second
retains seven discordant reads. Total correctness meets or exceeds HiPhase,
and the dominant connected core matches or exceeds HiPhase at each interval.
`audit_block.py` compares the exact original read names, alignment geometry,
flags, mapping quality, sequence and base qualities with HiPhase's BAM.

## Defect and repair

The useful MSA insertion contrast at 19,373,923 consists of one/two `TTCC`
copies (VCF 19,373,922). The physical caller's one-/two-base motif limit
excludes this four-base repeat and the preceding five-base `TTTAT` repeat.
The source's compound 10/14-base insertion at 19,377,346 segregates poorly:
complete sequence controls find the long allele in both parents. HiPhase
calls its longer insertion homozygous, while using the informative repeats.
Pgphase's artificial compound contrast consequently misorients reads even
though its main phase set already spans the separate `TTCC` orphan.

A new bounded physical caller handles pure three- to eight-base motifs and
reference/ALT or complementary ALT contrasts. It includes the entire reference
repeat (up to 256 bases in either direction) and external 16-base flanks.
Nearest net-length class and unique complete-sequence distance must agree,
with no more than two bases of slippage and one additional shared sequence
edit. Compound edits and mixed substitutions/deletions abstain.

MAPQ30 SNP-bearing molecules ending before the orphan independently calibrate
the upstream repeat: 12/1 and 0/17 core-haplotype counts. Separate MAPQ30
crossing molecules support both classes with 12 and 10 agreeing bridges.
The existing calibrated-repeat predicate rejects random gauge association,
contrary bridge parity and excessive joint error. Each read then supplies its
own bounded repeat call; contradictory clean phased SNPs or physical repeat
calls veto assignment. No truth or competitor calls enter production.

An initial hard MAPQ30 assignment floor withheld two already correct short
rescue reads at MAPQ21/17, leaving the first core at 84 correct. Calibration
keeps MAPQ30, while assignment checks each molecule's combined mapping and
external-base error against 20%. Those two complete, high-quality repeat
calls now enter the existing core. They cannot nominate a block join.

The native 19–20 Mb owner corrects 21 discordant labels, preserving all 4,033
old correct reads and all phased labels. Correct/discordant/unphased counts
become 4,054/20/154 from 4,033/41/154. Twenty-six HP/PS labels change, including
five old correct rescue labels. Candidate TSVs and VCF records remain
unchanged, and one-/four-thread candidates, VCF records and BAM assignments
are identical. `owner-results.json` and `physical-proof.json` record the checks.

The whole HiPhase target remains split at 19,377,345–19,395,544. HiPhase gets
102/131 correct there (77.86%); the new stage changes no variant phase sets
or block-link certificates. The existing 19 Mb unsupported-join regression
remains active. The nonsegregating compound variant calls themselves remain
unchanged. This repair accepts the two measured read connections, not the
entire HiPhase block.

Across the whole target, correct reads improve 4,373→4,394, versus HiPhase
4,389. Connected-core correctness improves 3,949→3,975, still 414 below
HiPhase's 4,389. Full chr20 gains 21 correct reads (230,901→230,922), with
discordant reads falling 6,642→6,621 and 19,069 unphased unchanged.
Every previously correct/phased read, all 64,483 complete VCF records
(including GT/PS), and every old block extent are preserved. No read-status
change occurs outside the target neighborhood. N50 remains 955,496 bp,
with 250 blocks and a largest block of 3,012,193 bp.

## Regression and reproduction

Both intervals enter the committed panel, certified manifest, identical-input
HiPhase reference and measured expectations. Their shared native owner retains
phase-set context. The regression asserts measured total/core floors, both
independent parental flanks, repaired molecules, and the two short rescues.
The starting binary fails 18 assertions; the repair passes 710. The warm
regression takes 0.76 seconds. Fifteen new pure-predicate checks exercise
motif eligibility, bounded slippage, sequence/length agreement and abstention.

The first complete suite exposed a separate historical score bug: `separated`
selected the phase set with the most tagged reads and then counted its correct
reads. Correct rescue promotions made that selection change, reporting 58
despite the old core's 59 correct reads remaining intact. The scorer now
selects the greatest correct-read count. A synthetic promotion counterexample
fails one of 16 assertions before this correction and passes all 16 after it.
The existing 0.42 floor and all other expectations remain unchanged; the
original-primary 80% and HiPhase certification logic is unchanged.

The final complete gap suite passes 17,841 assertions in five cases, and all
four replay-cache tests pass. All units, 1,551 phasing-predicate assertions,
HiFi/ONT goldens and thread-determinism checks pass. `validation.json` records
the final executable identity and the full native-panel preservation audit.
All 246 native output labels / 128 current independent replay requests preserve
every previously correct/phased read and complete VCF record. Eight output
labels change only read tags in the affected owner region (`panel-audit.json`).

```bash
LD_LIBRARY_PATH=/usr/lib/x86_64-linux-gnu make gap-dev-check
LD_LIBRARY_PATH=/usr/lib/x86_64-linux-gnu make gap-owner-check GAP=19.373
LD_LIBRARY_PATH=/usr/lib/x86_64-linux-gnu make window-tests
bash evaluations/2026-10-06-next-overlapping-block/replay.sh
python evaluations/2026-10-06-next-overlapping-block/audit_block.py
python evaluations/2026-10-06-next-overlapping-block/audit_preservation.py
python evaluations/2026-10-06-next-overlapping-block/audit_owner.py
python evaluations/2026-10-06-next-overlapping-block/audit_panel.py
```

Run from the repository root with pysam available for the Python audits.
The original input and HiPhase output are unchanged. Full chromosome metrics,
whole-target deficits and validation are recorded in the JSON reports here.
