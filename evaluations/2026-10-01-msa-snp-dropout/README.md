# Original MSA SNP dropout hidden by CIGAR backfill (2026-10-01)

## Cause and retained fixes

The owning 39–40 Mb chunk kept the 39,848,887 deletion and 39,856,144 SNP
in different BAM source phase sets. Both boundary sites were present.
Thirty-six primary MAPQ30+ molecules physically crossed them, but the original
MSA SNP matrix lacked a call on 35 of those molecules. The SNP's original
counts were 5 REF / 11 ALT (DP 16).

Recovery ran exact-CIGAR backfill before deciding whether MSA needed unplaced
reads. That backfill supplied 32 callable boundary pairs (31 same-allele,
one opposite), concealing the original dropout. These observations arrived
**after** the source HP/PS solve; they did not repair that solution.

Recovery now inspects the original source matrix for MSA SNP dropout before
backfill. It uses the existing 20-crossing-read minimum, known MAPQ30+, and an
exact upper binomial tail at p <= 0.01 for a majority of absent calls. The
probability is corrected for the two endpoint tests when both are MSA SNPs.
Known multiallelic alleles count as observations. This triggers a retry, never
a phase connection. Existing post-backfill admission decisions retain priority.
A retry requested only by original SNP dropout must preserve every phased
heterozygous source key; source transfer and stitch validation remain in force.
No truth, HiPhase output, or fixture coordinate enters production logic.

The admitted unplaced-read MSA restores the SNP to 34 REF / 36 ALT (DP 70),
and the ordinary source solve connects it with the deletion at the correct
same-ALT polarity. All 99 truth-scorable gap-overlapping reads are correctly
placed together in the owning replay; the prior dominant block placed 45/99.
The full chromosome also connects the pair.

A separate coordinate defect in the single-base insertion repair is fixed:
recovery windows use inclusive VCF anchors, while an insertion's internal key
is one base to the right. Membership now uses `sort_pos()`, so an insertion
anchored exactly at the right boundary is eligible. Its synthetic regression
uses a one-coordinate window and checks exclusion on either side. This change
alone leaves chr20 output unchanged.

## Matched chromosome result

The baseline includes the preceding uncommitted single-base insertion and
recovery-evidence indexing fixes after commit `74bc786`. These measurements
must not be attributed to a comparison against that commit alone.

| Measure | Before | After |
| --- | ---: | ---: |
| Output reads | 256,570 | 256,570 |
| Phased, truth-scored reads | 237,070 | 237,071 |
| Truth-correct assignments | 229,890 | 229,892 |
| Discordant assignments | 7,180 | 7,179 |
| Read concordance | 96.971359% | 96.971793% |
| Read phase sets | 692 | 690 |
| VCF keys | 62,368 | 62,368 |
| VCF blocks | 341 | 340 |
| Span N50 | 672,998 bp | 684,798 bp |

Every VCF key and unordered genotype is retained. A final run with corrected
two-endpoint testing and multiallelic call handling has identical HP/PS tags
and VCF fields to the accepted experiment. The input has no MAPQ255 records;
excluding unknown MAPQ from the new admission check cannot change this fixture.
The owning 56–57 Mb accuracy counterexample is identical to baseline.

## Refreshed competitor targets and permanent tests

HiPhase 1.6.0 (`ac3f399`) was rerun on the **current pgphase VCF**, identical BAM
and reference, with eight threads, default realignment and
`--ignore-read-groups`. This is the matched-callset diagnostic used to find
missing connections. The native DeepVariant benchmark was not rerun.

The established audit finds 62 exact-key joins across pgphase VCF seams outside
26–29.5 Mb. Three previously unlisted cases satisfy >=98% local read purity,
>=95% purity on both disjoint 10-kb flanks, at least five reads per flank,
consistent parental orientation, and identical input variant keys:

| Gap | Crossing molecules | HiPhase local correctness | Status |
| --- | ---: | ---: | --- |
| 528,850–542,052 | 11 | 203/205 | Open |
| 39,848,887–39,856,144 | 36 | 181/181 | Closed by this fix |
| 54,684,758–54,702,072 | 8 | 221/225 | Open |

All three are now in the committed panel, **89 to 92 windows**. The solved
case runs its owning chunk: its short replay already spanned the gap before
this fix. It has exact `spans=1`, unchanged accuracy floors, 100% gap separation,
a required deletion row, and a dedicated source-depth / DEL–SNP polarity /
parental-orientation regression. The old binary fails five assertions in that
regression, including separation 45/99 rather than 99/99 and the DP16 SNP.
The two open cases retain measured split expectations. Prior window floors
are not lowered. TOTAL floors are tightened to the sum of all per-window
floors: 73 joins and 1,288 in-gap heterozygotes. There are 15 tracked competitor
targets still open, plus four intentional split controls.

`current_hiphase_audit.tsv` contains the complete audit; `new_targets.tsv`
contains the three additions. Competitor separated fractions in the panel use
gap-overlapping, truth-scorable reads, so they differ from audit local purity.

## Rejected trials

- Moving all admission before backfill lost 12 heterozygous keys and raised
  chromosome discordance from 7,180 to 7,382. Requiring every retry to retain
  source rows reopened established joins and reduced N50 to 643,699 bp.
- Preserving legacy admission while adding general original-matrix admission
  retained keys and closed this gap, but raised discordance to 7,201 and failed
  the 56.064 / 56.662 Mb panel accuracy floors. Checking whether any source
  path spanned the group did not eliminate that counterexample.
- Direct Q30 SNP backfill before source phasing closed this gap but raised
  discordance to 8,207, changed some genotypes, and reopened an established
  61.758 Mb join. Requiring five clean reference-matching bases on each side
  still gave 8,107 discordant reads. These changes and their proposed unit tests
  were discarded. High base quality cannot validate every noisy-region allele
  placement, and adding calls before k-means can affect other MSA decisions.
- The retained SNP-dropout admission distinguishes the two cases without truth:
  35/36 SNP calls absent at the new closure, versus 10/26 at the 56 Mb
  counterexample. Ambiguous indel calls keep their existing recovery path.

## Reproduction and validation

`commands.sh BEFORE_BINARY AFTER_BINARY OUTPUT_DIRECTORY` repeats the matched
full chr20 runs. Derive parental truth with `scripts/make_truth_hap_map.sh`.
Scoring orients HP independently per read phase set and excludes secondary
and supplementary alignments. Graph output BAM is unaligned; primary read names
and truth are used for chromosome scoring, with original BAM coordinates for
local separation and flank orientation. No switch/flip estimator is used.

Temporary binaries, matrices, comparisons and logs are in
`/tmp/pgphase-gap-next3/`; the final chromosome run is `final-full/`.
The established audit is reproduced by
`evaluations/2026-09-29-expanded-hiphase-gaps/audit.py` with its HiPhase paths
pointed at the fresh same-callset output. Before and after metrics are saved
in `summary.tsv` and `metrics.json`.

Build, shared units, predicate and upstream/port parity gates pass with no new
warnings. The predicate suite has 244 assertions in 24 cases. The expanded
window panel passes 3,951 assertions in 49 Catch2 cases across 92 windows,
including all previously established joins and protected splits. `git diff
--check` passes.
