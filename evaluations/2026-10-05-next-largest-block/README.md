# Next largest HiPhase block still split in pgphase

The next target after the reviewed largest-block internal closure is HiPhase's
second-largest chr20 block, **8,946,171–11,127,753**, PS8946171:
**2,181,583 bp**, with 2,510 phased heterozygote records.

Current frozen pgphase output has two blocks with these exact outer endpoints:

| Phase set | First–last phased heterozygote | Span | Het records |
|---|---|---:|---:|
| 8946171 | 8,946,171–10,325,039 | 1,378,869 bp | 1,320 |
| 10343596 | 10,337,661–11,127,753 | 790,093 bp | 900 |

The missing connection is **10,325,039–10,337,661**, a boundary distance of
**12,622 bp**. The full chromosome output, rather than a shortened replay,
identifies this split. Neither boundary pair currently appears in the committed
window panel. This task identifies the target; it does not change the solver
or claim a closure.

## Read correctness and closure requirements

The whole interval has 9,742 original primary truth-scorable alignments.
HiPhase's output was checked against every corresponding original start/end,
CIGAR and sequence SHA256. All 9,742 match. Each tool's phase-set orientation
is measured on its whole output cohort; scoring includes abstentions.
Read eligibility uses original alignments, never pgphase's synthetic output
BAM. Rescue phase sets >=1e9 are excluded from connected-core counts.

| Whole block | Correct | Discordant | Unphased | Correct / all | Largest correct core |
|---|---:|---:|---:|---:|---:|
| Pgphase | 9,424 | 52 | 266 | 96.74% | 5,824 |
| HiPhase | 9,459 | 23 | 260 | 97.10% | 9,459 |

At the split, the same 127 original primary truth-scorable reads score:

| Gap | Correct | Discordant | Unphased | Correct / all | Largest correct core |
|---|---:|---:|---:|---:|---:|
| Pgphase | 103 | 10 | 14 | 81.10% | 45 |
| HiPhase | 97 | 10 | 20 | 76.38% | 97 |

Pgphase already meets 80% total correctness at this gap and exceeds HiPhase's
correct-read count. It fails the connected-core comparison. Its two main
cores have 41 and 45 correct gap reads; a pure phase-set union would connect
86, still below HiPhase's 97. A certified closure needs to connect at least
11 additional correct reads, while preserving existing correct assignments.
The 17 other currently correct gap reads are in rescue phase sets. To match
HiPhase across the whole 2.18 Mb interval also requires addressing the
35-read total correctness deficit; a label union alone cannot do that.

Disjoint +/-50 kb flank cohorts independently confirm the gauges. On the
left, pgphase's main PS has 192 unanimous orientation votes for hap1-paternal;
on the right its main PS has 229 unanimous votes for hap1-maternal. HiPhase's
same PS has 193 maternal / one paternal orientation vote on the left and 230
maternal on the right. A bridge must account for pgphase's opposite gauges.
These cohorts exclude all 127 gap reads. All-read flank correctness is 193/224
on the left and 230/230 on the right for both tools; lower left completeness
is retained in the comparison rather than omitted.

HiPhase's full block is well above 80%, but its gap-specific all-read score is
76.38%. Therefore the raw block rank and the read certificate are distinct
facts: this is the largest remaining internal connection by VCF span, not a
claim that HiPhase's own gap assignments meet the user's 80% criterion.

## Reproduction

`find_next_block.py` ranks all 191 multi-heterozygote HiPhase blocks against all
264 corresponding pgphase blocks, then audits the next internal split. The
already-investigated largest block has only a terminal endpoint difference;
it is explicitly excluded from selecting the *next* target. Rankings for all
blocks are in `ranked-blocks.json`; the target and read evidence are in
`next-block.json`.

```bash
/home/kokyriakidis/micromamba/envs/bench-phasers/bin/python \
  evaluations/2026-10-05-next-largest-block/find_next_block.py
```

Inputs are the frozen `tmp_gap_fix69/frozen_final/0` output, HiPhase's
`tmp_gap_fix48/competitor/hiphase_dv` output, the original annotated BAM and
the derived parental truth map. Production SHA256 remains
`477aaecb8adadf8f6f20bc380676de81f15d51f8eacbdb787c2cf5a8bf18140e`.
No pipeline replay, rebuild or full test-panel run was necessary. Rank,
alignment identity, parental scoring and connected-core checks passed;
Python compilation and diff whitespace checks passed.
