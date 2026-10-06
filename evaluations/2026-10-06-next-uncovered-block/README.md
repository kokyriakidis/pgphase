# Calibrated tandem insertions close the twenty-third-largest HiPhase block

The next qualifying open seam is **chr20:865,572–882,277 (16,705 bp)**,
inside HiPhase's **130,540–1,117,887 (987,348 bp)** block. Larger remaining
open seams fail HiPhase's 80% correctness rule on all original overlapping
primary truth-scorable reads. `selection.json` records the ranking and screening.

Before the fix, pgphase splits this span into 130,540–865,572 and
882,277–1,117,887. Its gap has 50/71 correct reads, 20 abstentions and one
wrong assignment; only 28 correct reads share one connected core. HiPhase
has 58/71 correct, seven abstentions and six wrong, all 58 correct in one core.
The original BAM alignments, CIGAR, sequence, qualities and flags are verified
identical for the competitor comparison. Abstentions remain in the denominator.

## Cause and repair

The right boundary already contains two complementary MSA- and alignment-
verified BAM alleles, **882,277 A>ATC / A>ATCTC**, with original counts
17 REF / 11 ALT and 13 REF / 15 ALT (DP28 each). The physical CIGAR insertions
are rotated CT/CTCT at 882,283, equivalent to TC/TCTC at 882,278. Recovery's
physical complementary-repeat caller supported homopolymers only, and the
existing short-insertion path required additional SNP-spanning molecules.
The clean left SNP is **863,406 C>A**; the first clean right witness is the
whole substitution **890,261 TT>AC**. This SNP desert has only one physical
SNP/insertion bridge. Catalog padding and replacement encoding also obscure
the minimal SNP and whole-substitution keys.

The new serial pass recognizes complementary two-/four-base repeat classes,
uses minimal SNP keys and whole physical substitution calls, and independently
calibrates the insertion gauge on downstream molecules that cannot call the
left SNP. Original primary, nonduplicate MAPQ30 molecules require Q10 inserted
bases, Q20 physical flank anchors and SNP/substitution bases, known qualities,
and summed base/mapping error <=1%. A 16-base flank alignment compares both
complete alleles; negative net lengths, tied distances and inconsistent length
classes abstain. Zero net insertion may enter the calibrated shorter class,
without creating a REF variant genotype. Both classes must pass the existing
p<=0.01 association check. Bridges must unanimously agree with odds >=999:1.
The calibration Wilson upper bound plus a 1% call allowance and bridge
posterior must jointly stay <=20%.

Independent CIGAR replay (`audit_physical.py`) reconstructs **23 calibration
molecules**, matrix `[[10,0],[1,12]]`, and one bridge. Its log odds are
**7.238005**, wrong-parity posterior **0.0718%**, and joint error bound
**18.4014%**. No parental truth enters production.

After the certified union, unassigned reads need their own physical class,
MAPQ30 and <=20% call error to join. Existing source rescues can retain their
assignment at MAPQ20 if the unique physical class agrees; a contrary rescue
requires a new physical call within the 20% error limit. An excluded graph
repeat at 883,357, with no alignment verification, formerly vetoed several
valid physical assignments. Such a row cannot veto this independent physical
certificate; clean and verified conflicting observations still veto it.
Unsupported reads retain their original tags. Candidate counts, alleles and
categories are unchanged.

## Results and preservation

| Full-chromosome gap result | Before | Fixed | HiPhase |
|---|---:|---:|---:|
| Original primary scorable reads | 71 | 71 | 71 |
| Correct | 50 | **59 (83.1%)** | 58 (81.7%) |
| Unphased | 20 | 11 | 7 |
| Discordant | 1 | 1 | 6 |
| Correct in one connected core | 28 | **58** | 58 |

The fixed geometry covers **130,540–1,117,887 exactly**, preserving every
prior block extent. Both disjoint marker cohorts agree on parental orientation.
The native owner and stitched 0–2 Mb continuation retain the SNP, insertion
and terminal variant gauge. Single-thread and four-thread owner candidates,
VCF records and BAM assignments are identical.

Full chr20 has **255 blocks**, down from 256. **N50 rises from 944,265 to
955,496 bp**; the largest block remains 3,012,193 bp. Among 256,610 primary
reads, correctness changes from 230,738 to 230,749, discordance from 6,676 to
6,674 and unphased from 19,196 to 19,187. Nine unassigned reads become correct
and one physically contradicted wrong rescue is corrected. All previously
phased assignments remain phased. Variant keys and common evidence are identical.

Removing validated reads from an old rescue PS changes its small residual
majority. One unchanged maternal read becomes scored discordant, and two
unchanged paternal reads become scored correct. This is explicitly reported
in the preservation audit: those three HP/PS assignments do not change, and
no previously correct read whose tags change becomes wrong. The audit uses
the committed C++ scorer's deterministic tie rule. An earlier Python audit
used encounter order for ties and reported 60/71 in the native owner; the
consistent score is **59/71**. No existing regression bound was relaxed.

The seam satisfies HiPhase parity. Across the entire 987 kb span, pgphase
still trails HiPhase by **15 correct reads / 38 connected-core correct reads**;
closing the geometry does not claim that separate whole-block deficit is fixed.

## Verification

The committed panel includes the new gap with measured exact `spans=1`,
independently measured HiPhase total/core counts, the 80% contract and parental
orientation. The new owner regression preserves the original insertion counts
and checks the following chunk's endpoint gauge. The old binary fails it.
All prior required-site, read-floor and HiPhase rows are preserved.

Build and all unit tests pass with no new warnings. The complete gap suite
passes **17,054 assertions in four test cases**, plus four cache-helper tests,
in **38.35 seconds** from saved states. All 1,551 in-memory predicate assertions
and HiFi/ONT golden output and thread-determinism gates pass. HiPhase measurements cover 126 panel
windows on identical original alignments. The native regression passes 637 assertions from
saved states in 0.58 seconds. All 230 prior output labels (120 independent
requests) are audited under the final binary; four labels including the owner
change, with no lost phased assignment or variant evidence. Full-suite and native-label audit results
are recorded in `test-runtime.json` and `panel-audit.json`.

Evidence: `physical-evidence.json`, `owner-results.json`, `continuation-results.json`,
`results.json`, `full-preservation.json`, `panel-audit.json` and `selection.json`.
Run `bash replay.sh [output-directory]` for the full chromosome;
`make gap-owner-check GAP=tandem-insertion` for the focused saved-state check.
Use `LD_LIBRARY_PATH=/usr/lib/x86_64-linux-gnu` for pipeline/tests and a Python
with pysam for the audit scripts.
