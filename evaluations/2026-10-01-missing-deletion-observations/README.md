# Missing deletion observations: preserve source contrasts and gauges

## Retained correction

Recovery backfill called MSA deletions only at their exact BAM CIGAR position.
When the alignment placed an equivalent deletion on a neighboring repeat base,
the candidate's required flank could be deleted and the call stayed missing.
The same edit already had a verified sequence-equivalence caller in physical
stitching. Reuse that certificate against the source chunk's loaded reference.
The shared implementation retains the existing 32-base shift and 64-base
length limits; worker-reference physical stitching is unchanged.

The retained repair fills only a missing ALT. It preserves every callable exact
contrast, existing MSA decisions, candidate coordinates, separate allele rows,
and source counts. Newly recovered evidence requires known MAPQ >=30 and base
qualities >=30, no different MSA indel at the same internal position, an
oriented source deletion, and an independently callable clean SNP at least
100 bases away on the same molecule in the same source PS. Every callable
clean SNP checked in that source gauge must agree; a reversal vetoes repair.
The stored query index names a surviving quality-checked flank. The function
is private to the existing backfill implementation. No new object dependency,
realignment, runtime truth/competitor input, or coordinate exception is added.
Ordinary BAM calling and source k-means remain unchanged.

## Representation finding and rejected experiments

**An MSA row's zero need not mean literal reference sequence.** Separate
complementary BAM rows use ALT-absence contrasts: another MSA allele can be
zero at this row. A physical reference-versus-deletion certificate used to
replace those zeros changes their meaning and breaks established source gauges.
The earlier suspicion that every different indel called zero was a false REF
was therefore too broad. Preserve those source contrasts, as the user requested.

| Trial | Observed failure | Disposition |
|---|---|---|
| Add verified shifted deletion ALTs before source k-means | Owning 36-Mb discordant reads 163 -> 292; no tracked target closes | Rejected |
| Replace all small-deletion backfill with literal edit calls | Chr20 correct +102 and discordant -53, but six protected joins open and N50 falls 739,888 -> 650,885 bp | Rejected |
| Keep exact ALT but reinterpret repeat-shifted REF/competing edits | Same six protected joins open; source rows/labels change | Rejected |
| Preserve callable contrasts and restore missing ALT without a source gauge check | At 36 Mb, Q30 repair adds 15 scored reads, eight discordant; multi-row repeat cohorts remain unreliable | Rejected |
| Preserve exact contrasts; simple row, Q30/MAPQ30, independent source-SNP gauge | All protected chromosome joins, VCF rows and HP/PS assignments retained | Retained |

The six reopened cases in the broad trials were 24.582, 47.004, 48.930,
53.945, 53.948 and 56.065 Mb. Their accepted owning-chunk outputs were also
replayed with the original binary to verify that these were actual regressions,
not differences in replay context. Expectations were never weakened.
The fourteen unresolved competitor targets were inspected outside the
centromeric exclusion. No additional tracked chromosome gap closes.

## Default chr20 comparison

Accepted baseline: `/tmp/pgphase-gap-next8/parity-full`.
Retained native run: `/tmp/pgphase-gap-next9/final-full`, eight threads,
156.55 seconds, identical annotated BAM, reference, catalog, GAF and defaults.
Competitors were not rerun.

| Metric | Baseline | Retained |
|---|---:|---:|
| Primary read names | 256,586 | 256,586 |
| Phased/truth-scored reads | 237,101 | 237,101 |
| Truth-correct / discordant reads | 229,921 / 7,180 | 229,921 / 7,180 |
| Per-read-PS concordance | 96.971755% | 96.971755% |
| Read PS | 686 | 686 |
| VCF records / blocks | 62,798 / 340 | 62,798 / 340 |
| VCF span N50 | 739,888 bp | 739,888 bp |
| Connected committed panel gaps | 74/92 | 74/92 |

All primary read names and HP/PS assignments match exactly. Every VCF key,
genotype and PS matches; no key is lost or added. This fixes a reproducible
observation loss, without claiming another chr20 closure or coverage gain.
Fourteen competitor targets and four regression controls remain open.

## Regressions and verification

- Add **93 predicate assertions** for left/right shifted ALT recovery,
  different-length abstention, exact ALT and zero preservation, complementary
  MSA rows, missing/low MAPQ and base quality, source-SNP absence/reversal/quality,
  and existing MSA decisions. Link the final regression against the actual
  pre-change caller: its missing-ALT assertion fails. The corrected predicate
  suite passes **491 assertions / 27 cases**.
- Shared unit tests, **27 port-parity** and **168,696 original upstream-parity**
  assertions pass; warning-free production build and `git diff --check` pass.
- Every one of the **105 distinct native window commands** runs on the final
  binary with four workers in **202.24 s**. The unchanged Catch runner scores
  those fresh outputs with exact-command and production SHA256 verification:
  **4,041 assertions / 53 cases PASS**, covering all 92 committed windows and
  their existing owning-chunk and parental-orientation checks.
- No panel row, required-site expectation, accuracy floor or ceiling changes.
  No new window is claimed closed, so no fabricated span expectation is added.
  Living behavior is updated in `docs/IMPLEMENTATION.md`; metrics, identities,
  rejected-trial scores and validation logs are kept here.

Reproduce the fast regression with:

```bash
make predicate-tests
./test_phase_predicates '[deletion-backfill]'
```

Run all integration gates with `make window-tests`. The native command inventory
is `native-inputs.jsonl`; its output paths are scratch destinations and can be
replaced when repeating the commands. Earlier outputs were used only to capture
the command inventory, never substituted for final native runs.
