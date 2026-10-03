# Reject false literal REF calls for shifted, unpaired MSA insertions

## Defect and correction

The exact-anchor lookup can return REF even when the existing reference-edit
verifier recognizes the selected insertion at a rotated repeat placement.
The independent physical SNP check previously returned the same `false` for
no callable SNP and an explicit contradictory SNP. In both cases, backfill
retained the preceding REF projection. A contradictory source orientation
cannot establish that opposite sequence allele.

Retain three states internally: no callable independent gauge, coherent ALT
gauge, and a contradiction. At the existing known MAPQ30/Q30, edit-equivalence
and spaced clean-SNP checks, an explicit contradiction leaves an unpaired
binary insertion unknown. Coherent evidence still admits ALT. Missing SNP
evidence and separate complementary MSA rows retain their existing source
projection; zero on a complementary row can mean ALT absence rather than
literal REF. The deletion caller retains its original behavior.

No candidate key, clustering algorithm or stitch threshold changes. There is
no new alignment, row merging, parental truth, competitor input or fixture
coordinate in production. The synthetic regression covers absence versus
contradiction and preserves source labels. Its contradictory-SNP assertions
fail the starting implementation. This correction does not claim to solve all
uncertified insertion projections; broadening it failed continuity checks.

## Owning-chunk results and guards

Fresh native replays use the same BAM, reference, catalog, GAF and four threads.
Truth is used only afterwards, orienting HP independently for each output PS.
The final narrowed implementation exactly preserves the read tags and VCF rows
of the 1–2, 19–20 and 41–42 Mb owning replays. The 1.18 Mb connection survives.
Retain the original 41 Mb coverage/correct/error floors and ceilings.

Add the 2–3 Mb owning replay to the permanent representation/transfer test. It
protects both separate insertion rows and the 4,088 scored / 4,026 correct /
62-error limits. The rejected expanded-call trial produced 1,091 errors there,
so this guard detects its wrong whole-block connection. Small coordinate
replays cannot substitute for the owning gauge in this counterexample.

The rejected missing-certificate version looked promising locally: at 41 Mb
it corrected seven erroneous tags and abstained on two former errors without
losing a correct read or changing VCF rows. At 19 Mb one former error abstained.
But full chr20 lost three correct assignments net, made 13 formerly correct
reads wrong, and reopened 1.180618–1.194189 Mb. Reject it and restore the original
coverage expectations rather than weakening that existing span gate.

## Rejected experiments

1. Jointly fill missing complementary insertion calls from exact or rotated
   Q30 CIGAR edits, preserving original calls and source GT/PS. At 19 Mb this
   adds 25 read observations per row (DP 7 -> 32), but leaves the target gap and
   local tags unchanged. Full chr20 errors rise from 7,127 to 9,342, with 2,292
   formerly correct reads becoming wrong. VCF keys stay fixed and none of the
   100 panel coordinates gains a span. Reject this implementation.
2. Restrict new paired calls to the established source membership and gauge.
   The 2–3 Mb owning replay still changes 1,062 formerly correct reads to wrong;
   its errors rise 62 -> 1,091. This restriction is insufficient.
3. Freeze source-path certificates before admitting the paired calls. This
   fails to prevent the 2 Mb inversion and increases owning 19 Mb errors from
   41 to 381. The source-certificate timing hypothesis is not a sufficient fix.
4. Apply the literal-REF abstention broadly to complementary insertion rows.
   Owning 41 Mb loses 20 keys, gains one key, and loses six formerly correct
   reads through changed retry selection. Reject this broader interpretation;
   preserve the separate-row source semantics.
5. Reject an unpaired insertion whenever a source certificate is missing,
   without distinguishing absence from contradiction. Full accuracy rises
   96.995350% -> 97.000628%, but correct reads fall 230,072 -> 230,069 and panel
   spans fall 86 -> 85. No keys disappear, but a previously fixed gap reopens.
   Local accuracy alone is insufficient; reject this version.
6. Remove the deferred homopolymer recall's source-allele gate at 24 Mb. It
   restores a physical seven-T/eight-T bridge and yields the correct whole-MEC
   optimum, but one independent read half ties. Owning tags stay unchanged and
   the gap stays open. Restore the original admission gate; do not force a join.

Preserved JSON files record the complete counts and correctness transitions.
Trials are under `test_data/tmp_gap_next34/`; they are not accepted outputs.

## Remaining representation issue

At 19.377 Mb, the separate MSA tandem-repeat insertion observations discriminate
parental reads poorly even though neighboring clean SNP gauges are coherent.
Several high-quality reads carrying the longer CIGAR allele have contrary
primary MSA calls. These existing primary calls are deliberately preserved;
restoring missing alleles alone does not repair their source clustering gauge.
HiPhase DV has 111 correct / three erroneous local overlaps in one block versus
pgphase's starting 94 correct / 31 erroneous overlaps in three blocks. A proper
joint representation/orientation repair is still needed before those flanks
can be joined safely. This investigation does not establish that read coverage
alone or one nominal BAM source PS certifies that connection.

## Final chromosome validation

Starting binary SHA256:
`d3b5ad20b227cc051d1d5280a6aa4fc82cbc1c1bc9393082aa5533005456a1c8`.
Final binary SHA256:
`b36dd6d8ea32b1f1707b64efc51878b2d8641344ef1b0f009f08eee94b1bac0b`.
Outputs are `test_data/tmp_gap_next27/run-length-final/` and
`test_data/tmp_gap_next34/full-contradictory/`.

| Metric | Before | After |
|---|---:|---:|
| Truth-scored phased reads | 237,199 | 237,199 |
| Correct assignments | 230,072 | 230,072 |
| Discordant assignments | 7,127 | 7,127 |
| Conditional read accuracy | 96.995350% | 96.995350% |
| Read phase sets | 661 | 661 |
| VCF keys | 63,630 | 63,630 |
| VCF phase blocks | 333 | 333 |
| Span N50 | 774,189 bp | 774,189 bp |

Every HP/PS tag and every VCF row is unchanged. No candidate key disappears,
no old SNP gauge changes and no old block merges. The coordinate panel remains
86 spanned windows out of 100; fourteen remain open, including four controls
and ten historical competitor targets. The newer 7.28/24.12 Mb owning boundary
cases outside that coordinate panel also remain split. No new closure or chr20
accuracy improvement is claimed for the retained fix.

Build, units, 987 predicate assertions/42 cases and HiFi/ONT TSV/VCF goldens
pass, including HiFi one/four-thread determinism. The starting code fails two
contradictory-SNP assertions. No new build warning or dependency is introduced.
Original panel expectations and owning coverage floors are preserved. The new
2 Mb owning guard protects against the rejected expansion's wrong join.

The native chromosome run takes 292.18 seconds/eight threads, and 105 fresh
native panel requests take 387.70 seconds with four concurrent workers. These
are concurrent validation timings, not a controlled runtime comparison.
`native-cache.json` records the final binary, exact normalized CLI arguments,
input stats and native output hashes. `replay_cached_panel.py` verifies all of
those before scoring a cached output; unmatched requests execute natively.
The added 2–3 Mb owning regression is such a new native request.

The complete window suite passes **7,060 assertions in 69 cases**, including
the new owning replay and every existing span, orientation, accuracy and
coverage gate. No panel expectation is refreshed or lowered.
