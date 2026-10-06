# Next target: HiPhase's fourth-largest chr20 block

The latest frozen pgphase output still splits HiPhase's
**52,696,410–54,419,033** block, PS52696410, into two blocks.
HiPhase's span is **1,722,624 bp**, its fourth-largest, with 2854 phased
heterozygous records. The second- and third-largest HiPhase spans are now
covered by single pgphase blocks. The largest retains its separately reviewed
12,277 bp terminal difference; it has no internal split. The complete ranking
reports that endpoint difference rather than treating it as closed.

| Pgphase PS | First–last phased heterozygote | Span | Rows |
|---|---|---:|---:|
| 52696410 | 52,696,410–52,696,940 | 531 bp | 2 |
| 52711825 | 52,711,825–54,419,033 | 1,707,209 bp | 2125 |

The missing connection is **52,696,940–52,711,825**, a **14,885 bp** boundary
distance. The two outer endpoints already match HiPhase's target exactly.

## Original-read comparison

Eligibility comes from original primary BAM alignments with parental truth.
Unphased reads, including names absent from a tool's output, remain in the
denominator. HiPhase's geometry, CIGAR and sequence SHA256 match all **7030**
original alignments overlapping the block. Each phase set's full chromosome
read cohort establishes its parental orientation. Rescue PS >=1e9 does not
count toward a connected core. No original overlaps lack parental truth.

| Interval | Tool | Correct | Discordant | Unphased | Correct / all | Largest correct core |
|---|---|---:|---:|---:|---:|---:|
| Whole HiPhase block, 7030 reads | Pgphase | 6935 | 18 | 77 | 98.65% | 6799 |
| Whole HiPhase block, 7030 reads | HiPhase | 6963 | 10 | 57 | 99.05% | 6963 |
| Seam, 136 reads | Pgphase | 119 | 8 | 9 | 87.50% | 47 |
| Seam, 136 reads | HiPhase | 117 | 0 | 19 | 86.03% | 117 |
| Reads physically spanning both seam endpoints, 9 reads | Pgphase | 8 | 1 | 0 | 88.89% | 4 |
| Reads physically spanning both seam endpoints, 9 reads | HiPhase | 8 | 0 | 1 | 88.89% | 8 |

**This is a qualifying next gap:** HiPhase exceeds 80% correctness and connects
117 correct reads. Pgphase already exceeds its local total-correct count but
needs at least **117 correct connected-core reads** to accept the closure.
Joining the existing left/right local cores would connect only **59** correct
reads (47+12); **60** additional correct reads have rescue labels (8+52).
At least 58 of these need reliable core attachment if no additional reads are
corrected. The eight current discordant reads also require investigation.

Whole-block parity additionally requires at least **28** more correct original
reads and raising the connected correct core from 6799 to at least 6963.
A local seam closure alone must not be described as full read parity with
HiPhase without measuring that remaining deficit.

Only two original truth-scorable primary reads end within the disjoint 50 kb
left flank: pgphase has one correct and one discordant there, while HiPhase
has two correct. The current left flank cannot establish a reliable parental
orientation independently. The right disjoint flank has 218 correct pgphase
core reads; HiPhase connects 222. A future fix must validate orientation with
adequate independent physical evidence rather than claiming the pooled
read score proves the join. Original alignments physically spanning both
endpoints are present; this is not a gap without molecule coverage.

## Boundary calls and next investigation

Pgphase phases GGA>G at 52,696,410 and G>A at 52,696,940 in the small left PS.
HiPhase phases the same alleles in its full target PS, with the opposite
numerical HP gauge. Pgphase's next block starts with an A insertion,
G>GA at 52,711,825 (NOISY_CAND_HET; depth 12, REF7/ALT5). HiPhase has no record
for that exact insertion. Both tools phase the TTGTG>T deletion at 52,715,881,
but pgphase depth is 23 (REF13/ALT10), compared with HiPhase depth 79
(REF35/ALT25). The nearby HiPhase indel records include homozygotes, RefCalls
and a NoCall, rather than additional phased markers through the seam.

These differences identify what to inspect next: recovery source observations,
the short insertion's independent phase path, the shared deletion's retained
calls and physical bridging molecules. They do not yet establish a production
bug or certify a safe union. Exact nearby records are in boundary-calls.json.
There is no existing panel row for this seam; a fix must add a native-owner
regression and measured span/read/parental checks.

## Reproduction and provenance

Production is unchanged at commit fb064fe; binary SHA256:
`7e50015b328beeed83607d502a1e6d6a928d2fea8dbbfb77328c3e5610e76be7`.
Inputs are the completed 67-chunk output in
`test_data/tmp_gap_fix73/frozen_final/0`, HiPhase in
`test_data/tmp_gap_fix48/competitor/hiphase_dv`, the original annotated BAM,
and `test_data/derived/chr20_truth_hap.tsv`.

```bash
/home/kokyriakidis/micromamba/envs/bench-phasers/bin/python \
  evaluations/2026-10-05-fourth-largest-block-target/find_next_block.py
```

Ranked-blocks.json compares all 191 HiPhase blocks with all 261 pgphase blocks.
Next-block.json records the selected target, original-read scores, disjoint
flanks, nine spanning molecules, alignment checks and input fingerprints.
Validation.json records the assertions and acceptance targets. This task
identifies the next target; no pipeline rebuild/rerun or production/test
expectation change is needed.
