# Restore exact SNP branches from alternative catalog alleles

## Defect

The default binary graph projection retains observations of the catalog REF
walk and the selected complete ALT walk. A different ALT walk is missing in
that binary profile, even when it traverses the identical SNP branch and differs
only elsewhere in the snarl. This can make an internal graph path appear
one-sided after BAM recovery, blocking an otherwise supported connection.

In the owning 20–21 Mb chunk, the selected SNP at **20,523,889 A→T** is
catalog ALT 1 of the snarl at 20,523,884. Its oriented local branch is
`>115761253>115761255>115761256`. Other catalog walks preserve this branch
while changing preceding substitutions or deletions. Seven reads have a callable
ALT in the BAM channel but no graph call at this SNP. Recovering their exact
graph branch repairs the path proof and joins the seam at
**20,506,159 C→CT / 20,517,197 A→AAT**. The two insertion ALTs remain on
opposite haplotypes. The right SNP at 20,523,889 is on the same ALT haplotype
as the left insertion. The later SNP at 20,540,928 is opposite.

Another nine REF calls are restored at a distinct catalog SNP near 20,323,422;
they do not create another join. `supplemental-evidence.json` records all 16
changed observation-channel tuples in the owning recovery matrix.

## Correction

After graph k-means and before BAM recovery, supplement only already phased
binary clean SNP candidates:

1. Normalize the selected VCF allele and require a single-base substitution.
2. Require REF and the selected ALT walks to differ at exactly one internal
   node, with all other oriented nodes identical.
3. Locate the oriented three-node REF or ALT branch in each catalog walk.
   Exactly one occurrence supplies an allele; absent, repeated, opposing or
   differently oriented branches abstain.
4. Apply the existing mapping-quality floor, skipped-read rule and source-walk
   conflict veto. Deduplicate each read/candidate before adding observations.
   Conflicting full walks remain unknown even if they project to the same
   local branch. Conditional child sites cannot bypass parent gating.
5. Require the original genotype counts plus the proposed distinct-read calls
   to remain within the original minimum and maximum allele-frequency limits.
6. Fill missing observations only and rebuild the profile interval index once.

Candidate counts, categories, genotypes, initial PS labels and read HP gauges
remain those of the original full-walk solve. Reported genotype depths therefore
continue to describe that solve, while recovery also sees these supplementary
exact branch observations. Existing recovery/stitching decides whether to join.
This is not a new alignment, a new variant representation, an ALT-to-REF
shortcut, a search-bound change or a relaxed stitch threshold. The production
function takes no truth, competitor calls, coordinates or fixture identifiers.
Both graph worker paths use it; graph-only runs retain their existing behavior.

## Results

| Measurement | Accepted before | Accepted after |
|---|---:|---:|
| Owning 20–21 Mb phased / correct / discordant reads | 3,876 / 3,840 / 36 | unchanged |
| Owning VCF blocks | 7 | 6 |
| Local overlapping phased / correct / discordant reads | 110 / 109 / 1 | unchanged |
| Local read phase sets | 4 | 2 |
| Correct local reads in the dominant block / input overlaps | 65 / 111 | 78 / 111 |
| chr20 truth-scored phased reads | 237,134 | unchanged |
| chr20 correct / discordant reads | 229,981 / 7,153 | unchanged |
| chr20 read concordance | 96.983562% | unchanged |
| chr20 read phase sets | 667 | 665 |
| chr20 VCF blocks | 334 | 333 |
| chr20 VCF keys | 63,492 | unchanged |
| chr20 span N50 | 756,878 bp | unchanged |

Every phased read retains its individual truth classification. There are 272
changed HP/PS tuples and 11 changed VCF rows from phase-label/gauge propagation,
with no lost or gained variant keys and unchanged allele genotypes.

Saved native-DV HiPhase phases 109 of the 111 input gap overlaps, with 107
correct and two discordant (**98.17%**). pgphase phases 110 with 109 correct
and one discordant (**99.09%**). HiPhase still places more of these reads in
one block (107 versus 78 correct); the independent BAM read rescue blocks are
preserved. This fixes the exact core connection without claiming read-block
continuity parity. The competitor is the saved October 1 run using the same
BAM, not a new runtime benchmark or a fresh run on pgphase's changed callset.

The full run takes 254.45 s during concurrent validation, not a controlled
runtime comparison. The chromosome output is retained at
`test_data/tmp_gap_next27/accepted/`. Binary SHA256:

`cb3b4ad62fb671abb1f375e02a48a44b1eac2553158b7ab39bc6c9755d94fc06`

The extent audit falls from 314 to 313 gaps. **40 distinct noncentromeric
competitor-supported nominations remain**, with 50 records. The new seam was
not among the preceding 40 nominations because an overlapping phase-set extent
obscures the exact split. These nominations do not establish an exact allele
bridge; they must still be investigated. The new exact connection is now a
permanent panel and owning-chunk regression despite that extent-audit limitation.

## Rejected experiments

- Admitting unplaced MSA reads across all targeted solves changed existing allele
  consensuses. At 24 Mb it lost the +8 insertion description and moved 12
  truth-correct reads to wrong assignments. The flag is restored to its accepted
  value; standalone BAM behavior is unchanged.
- A homozygous physical REF anchor veto changed established blocks at 7 and
  33 Mb, including 36 correct-to-wrong transitions in the 7 Mb replay. It is
  absent from production.
- Exact branch supplementation without the genotype-compatibility check joined
  the 7.264–7.280 Mb nomination but left seven correct reads and one wrong read
  unphased and exchanged one correct/wrong assignment pair. The weak graph SNP
  had only 3 REF / 3 ALT calls; the projection added 42 REF calls, reducing its
  inferred ALT frequency to 6.25%. The final AF guard rejects that projection.
  The accepted 7 Mb replay preserves every output tag and truth classification.
- At 24.121 Mb, many CIGAR insertions are +3/+7 T bases near the existing +4/+8
  MSA alleles. Only one molecule links an insertion class to clean right-side
  SNPs in the physical audit. No new repeat class or unsupported whole-block
  join is added.

## Regression and reproduction

The new exact-key, ALT-parity, owning-read and parental-flank regression has
763 assertions; the accepted pre-fix binary fails six of them, including the
span and dominant-block improvement. Existing expectation floors and ceilings
are preserved. The panel grows to **99 coordinate cases / 85 required spans**.

Standalone units cover exact REF/ALT branches, unknown and reversed nodes,
repeated branches, duplicate/conflicting source walks, mapping quality,
idempotence, preserving known calls and gauges, the AF guard, unphased rows and
conditional children. All standalone unit binaries pass. Phase predicates pass
872 assertions / 39 cases. Fresh native panel requests use the frozen final
binary and are recorded in `native-validation.json`. Window tests pass
**5,708 assertions / 66 cases** across the expanded panel.

```bash
make -j8
make unit-tests
./test_phase_predicates
make window-tests
```

The new owning-chunk regression can be run alone with
`./test_gap_windows '[snp-branch]'`. Its replay includes the full 20–21 Mb chunk;
a short replay is not used as a substitute for the established phase-set gauges.

Comparison and nomination scripts are reused from the sibling
`2026-10-02-padded-indel-context` and `2026-10-02-deletion-physical-paths`
evaluations; their `--help` lists explicit input/output arguments. Truth is used
only in evaluation and regression checks.
