# Joint gap phasing on unchanged observations

Investigate whether using all usable graph/BAM gap sites fixes the remaining
splits. Use accepted binary SHA256
`cb3b4ad62fb671abb1f375e02a48a44b1eac2553158b7ab39bc6c9755d94fc06`.
Fresh owning-chunk snapshots are under
`test_data/tmp_gap_next27/joint-control-{7,24}/`; diagnostic work is under
`test_data/tmp_gap_next29/`. No production behavior, binary, or test expectations
change in this investigation.

## Controlled solve

Extract the existing `solve_recovery_mec` kernel unchanged, then evaluate both
relative orientations of the complete established flank blocks. Free gap sites
have binary observations with both alleles present; genotyped homozygotes and
true multiallelic observations are excluded. Include graph repeat indels and
off-center sites in the all-site arm. Variables extend from the first left
anchor to the first right SNP; all oriented sites of both flank PS remain
available, including far flank context. This is an evaluation model, not a
literal invocation of the production eligibility/stitch wrapper.

Compare centered extras (`|AF−0.5| ≤0.12`) against all such extras. Compare
uniform weights against strict SNP priority: one SNP mismatch costs more than
all possible indel mismatches. Evaluate the full evidence and the two existing
FNV read halves separately. Parental truth is used only afterwards to identify
the correct relative gauge. The kernel remains exact for this bounded replay;
these fixtures need at most 11 variables, below the 20-variable replay guard.

| Control | 7.264321–7.280346 Mb | 24.121713–24.131707 Mb |
|---|---|---|
| Current observations, centered extras | Correct full optimum; read halves disagree | Disconnected |
| Current observations, all extras | Correct full optimum; read halves disagree | Disconnected |
| HiPhase local calls at shared sites, all extras; original flank rows fixed | SNP priority gives correct full and both-half optima | Full optimum correct; halves disagree |
| HiPhase local calls, all extras; boundary indels free, complete SNP flank gauges fixed | Same as preceding row (fixed flank rows were SNPs) | Correct full and both-half optima |

At 7 Mb, admitting all extras increases variables 5→11 but does not fix the
current matrix's contradictory SNP evidence. With SNP priority, current
all-site costs for same/flip are 5,510/3,427 overall, 4,219/1,302 in half 0,
and 1,291/2,125 in half 1. The correct connection is flip. Substituting the
actual HiPhase local shared-site calls gives 5,109/3,022, 3,809/2,138, and
1,299/884: all three prefer the correct flip. Uniform weights still disagree.

At 24 Mb, the original matrix remains disconnected even with all extras or
free boundary indels. With HiPhase calls and free boundary indels, centered
extras prefer the wrong connection overall. Including all extras admits the
second insertion and deletion, plus two graph repeat rows: five variables.
Costs are 123/120, 59/55, 64/63 with uniform weights, and
15,085/15,082, 5,840/5,836, 9,245/9,244 with SNP priority. The correct flip
wins in full and both halves. Margins remain small; this does not establish a
chromosome-wide safe stitch rule or emitted read-tag accuracy.

The two-orientation kernel calls take about 4–9 ms at 7 Mb and 23–53 ms at
24 Mb, including subprocess overhead. Matrix loading is additional. Full
pipeline replays take 11.04 and 8.34 seconds; neither is a chromosome benchmark.

## Where the observations differ

Recover genuine HiPhase local allele vectors from serial `-vvv` trace output
without modifying HiPhase. Match complete vectors to indexed BAM fetch order
using its exact flag/MAPQ filter; assert both variant counts and record counts.
These runs use the same accepted current pgphase VCF and original BAM as the
preceding joint-call controls. `observation-comparison.json` records each
shared-site discrepancy and bridge read.

At 7 Mb, the 7,280,356 SNP has 12 callable graph observations that HiPhase
marks ambiguous and three opposite allele calls. Calls at the left 7,264,321
SNP all agree. Thus the difference is boundary allele assignment, not absence
of a site record. Changing this observation subset makes the existing exact
kernel agree across read halves under SNP priority.

At 24 Mb, HiPhase has 29 reads with callable pairs between the left SNP/
insertion sites and right deletion/SNP sites; pgphase has zero of those pairs.
Two MAPQ-60 reads reach the clean right SNPs as well as the left insertion:

- `m84031_231217_062403_s3/7996640/ccs`: CIGAR inserts seven T bases, minimum
  insertion BQ 27. HiPhase calls the eight-T row ALT; pgphase has no insertion
  call. Three right SNP calls agree between tools. The deletion call differs.
- `m84031_231217_062403_s3/132189840/ccs`: CIGAR inserts one T base, BQ 10.
  HiPhase calls both insertion rows REF; pgphase has no insertion calls.
  Three right SNP calls agree. This is noisy allele evidence, not an exact
  match to either four- or eight-T haplotype.

The loss precedes transfer: the BAM recovery-source matrix already has zero
paired insertion/deletion calls. Injection preserves its callable counts:
four-T row 7/14, eight-T row 10/15, deletion row 19/10. Source PS labels alone
therefore do not prove a heterozygote-only read chain.

## Reproducible MSA representation defect

`msa_site_event_allele` recognizes literal REF or the selected insertion ALT.
For a four-T/eight-T diploid contrast, an exact read of the other consensus
becomes unknown on the selected row, although this source row uses zero as
absence of that ALT. The one-error fallback also requires both consensuses to
have valid, distinct row alleles, so this inconsistency can disable fallback.
`probe_msa_contrast.cpp` calls the actual linked public helper on complete
synthetic alignments and prints:

```
selected=1 other_verified_consensus=-1 literal_reference=0
```

This proves one information-loss mechanism, not that every missing call in
the real window has this cause. HiPhase's separate four-T binary row also
calls many eight-T reads ALT because it compares four-T against literal REF;
copying those binary calls directly would violate pgphase's complementary
source-row semantics.

## Design implication

Use all usable gap observations jointly, but first preserve their actual
allele identities. Repair complementary insertion recognition against the
existing two MSA consensuses: exact other-consensus absence can project to
zero, while an unrelated third allele stays unknown. Validate noisy calls
against the diploid contrast, rather than independently comparing each ALT
only with reference. Keep the two rows and their original representations.

Protect established clean SNP flank gauges; uncertain imported boundary
indels need to be gap variables in a joint solve. An AF close to 0.5 should
not be a prerequisite for every intermediate site. Confidence must come from
the joint read evidence and ambiguity handling. Changing eligibility alone
does not recover the missing calls or correct boundary calls.

The HiPhase-call overrides are diagnostic and change the observation matrix;
only the original-channel arms are same-matrix solver comparisons. These
replays do not reproduce HiPhase's 80/10 local SNP/indel weights, global
realignment, or A* queue handling. No new alignment method is introduced into
pgphase, and truth never enters the kernel or candidate selection.

Reproduction helpers and machine-readable costs are retained alongside this
report. Scripts intentionally contain evaluation-only fixture indices and
paths. `build_replay_solver.py` generates the saved `mec_replay.cpp`; compile
with `g++ -O3 -std=c++17 -Wall -Wextra`. Run `replay_joint_matrix.py` with
chunk 7 / PS 7236076,7280339 or chunk 24 / PS 24103779,24517362; add
`--hi-local-calls` and/or `--free-flank-indels` for the diagnostic controls.
