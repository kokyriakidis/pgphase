# Complementary deletion fallback: preserve evidence and reject unsafe joins

## Retained fixes

The physical stitch preferred a clean SNP pair, but when that pair had no
callable spanning reads it tried only a single left MSA indel. That helper
deliberately rejects complementary overlapping deletion rows. The existing
complementary-pair helper is now tried first in that zero-pair fallback. Its
requirements are unchanged: continuous graph SNP paths on both sides, a
verified ALT of exactly one deletion row, observations of both deletion alleles,
unanimous parity, and the existing quality odds bound of 0.001. A callable
clean SNP pair still takes priority.

The complementary helper previously required an exact CIGAR position although
other physical deletion bridges already recognized equivalent placements. It
now uses the same reference-edit equivalence check within 32 bases for
deletions up to 64 bases. The candidate coordinate and length remain unchanged.
A different deletion length, another edit within the verified span, a mismatch,
or insufficient base quality causes abstention. Separate rows are not merged,
and ambiguous REF calls at a complementary locus never cast votes.
The pair helper also stops rechecking flanks at the original candidate
placement: one of those bases can lie inside a correctly shifted deletion.
The shared caller already checks the actual placement and all surviving
reference bases across both coordinates at Q30.

The existing equivalent-deletion implementation is extracted to
`bam_equivalent_deletion_allele` in `collect_phase_noisy.cpp`, with the graph
wrapper retaining its worker-local reference cache. This makes the allele
check directly testable. No new alignment, filtering threshold, truth input,
competitor input, or fixture coordinate enters production code.

## Findings and rejected trials

- **50,548,245–50,562,066:** twelve primary BAM molecules span the boundary
  sites. Some one- and two-base deletions are placed two bases right of the
  MSA rows. Directly auditing those twelve reads at Q30 changes the one-base
  row from zero ALT / six REF / six unknown to one ALT / five REF / six unknown;
  the two-base row changes from zero ALT / one REF / eleven unknown to three
  ALT / one REF / eight unknown. Four ALT observations were lost by exact
  position matching; they remain separate allele rows. The counts are in
  `boundary-call-audit.tsv`, with the calling tool in `audit_boundary_calls.cpp`.
  Restoring the missing complementary fallback exposes a further
  path veto at graph SNPs 50,719,976–50,719,983 inside the right block. The
  strict physical certificate cannot resolve that edge: Q30 CIGAR pairs are
  45 `A/T` and one `A/A`, so they do not independently support both haplotypes
  of those graph SNPs. This target remains split; raw CIGAR homozygosity is
  insufficient to reinterpret the graph/MSA alleles in this complex region.
- **58,366,458–58,385,702:** a trial added the right-hand SNP-to-indel fallback
  when the farther clean SNP pair had zero calls. It joined the target, but
  owning-chunk truth counts changed from 3,492/3,549 correct to 3,428/3,550:
  64 fewer correct assignments and 65 more discordant. The deletion moved
  from source PS 58,305,921 to live PS 58,391,092 while retaining its source
  allele orientation. Certifying the right block's own source does not certify
  that moved deletion's relation to the right SNP path. Exact one-base
  deletion calls paired with the right SNP occur on both SNP alleles (7/7).
  The right-hand fallback is **not retained**.
- Extending missing-only shifted-insertion repair beyond its existing
  single-base scope did not close any of the five inspected owning chunks.
  At 19 Mb it changed 4,031 correct / 41 discordant to 3,746 / 328 despite
  Q30 exact edit equivalence. Candidate allele certainty alone does not validate
  the source block's inherited gauge. This trial is **not retained**.

## Default chromosome results

Baseline: `/tmp/pgphase-gap-next5/final-full`; retained production run:
`/tmp/pgphase-gap-next6/verified-full`. Same annotated BAM, graph catalog, GAF,
reference, chr20, defaults, eight threads.

All **256,586** primary read names and HP/PS pairs are identical. Every
nonheader VCF row is identical. `identity.json` records the exact comparison;
`metrics.json` records the independent truth scoring.

| Metric | Baseline and retained code |
|---|---:|
| Phased/truth-scored reads | 237,101 |
| Truth-correct reads | 229,920 |
| Discordant reads | 7,181 |
| Per-read-PS concordance | 96.971333% |
| Read phase sets | 686 |
| VCF records / blocks | 62,798 / 340 |
| VCF span N50 | 739,888 bp |
| Connected committed panel gaps in full chr20 | 74/92 |

**No additional tracked gap closes. Fourteen competitor targets remain open.**
Native competitor runs were not repeated. These fixes restore eligible
evidence without weakening the path and allele guards that currently abstain.

## Regressions

- The new synthetic deletion case has **104 assertions** across ten sections:
  exact ALT, right- and left-shifted ALT, distinct deletion lengths, clean REF,
  interfering and unrelated edits, flank mismatches, unknown/low quality, and
  incomplete read coverage. Replacing the call with the old exact-only caller
  fails five assertions. The full predicate suite has 398 assertions / 26 cases.
- The new owning 58–59 Mb case has **17 assertions**. It checks exact boundary
  keys, keeps the left SNP separate, verifies the moved deletion remains in
  the right PS, and checks parental orientation, concordance, and separation.
  The rejected right-hand fallback fails three assertions: the forbidden join,
  parental switch, and concordance. The accepted implementation passes.
- Shared unit tests, 27 port-parity assertions, and 168,696 upstream-parity
  assertions pass. No previous window floor or required-site expectation is
  lowered. The panel remains 92 windows, with this additional owning-chunk
  negative regression.
- The complete final panel passes **4,023 assertions in 52 Catch2 cases**.
  All 105 distinct native pipeline commands were executed on the final binary
  with four concurrent workers in **234.08 seconds**. The unchanged Catch2
  runner then scored those freshly generated outputs through a wrapper that
  verifies the complete command arguments and production binary SHA256.
  This is a fresh native replay, followed by cached assertion checking; earlier
  outputs were used only to capture the command inventory, never as final
  test results. `native-replays.log`, `panel.log`, and `validation.json` record
  the result. Assertion totals include fewer repeated input-index requirements
  when an earlier case has already loaded the same window's read spans.

Reproduce the new fast checks with:

```bash
make predicate-tests
./test_gap_windows '[moved-deletion]'
```

The complete window suite is `make window-tests`. Logs accompany this report.
