# Close the first seam in HiPhase's third-largest block

Target: chr20:24,121,713–24,131,707, a 9,994 bp boundary distance inside
HiPhase's 24,103,779–25,944,471 block (1,840,693 bp). Baseline production
SHA256 is `15534777ecd8382e370ade340e6e8998a3fc1b1ab4ab390281b3ff4a22bbedb4`;
final SHA256 is `81ad9a79ed25b71494b8e7cf4288ff31f56c98f133bd908af4c5cf924ef1775a`.

## Cause and fix

The left core contains one original graph SNP, C>T at 24,103,779, a physically
confirmed recovered A>G SNP at 24,105,188, and complementary MSA insertion
rows C>CTTTT / C>CTTTTTTTT at 24,121,713. The graph-path validator rejected
the singleton even though it has no internal graph edge. The recovered SNP
still needs independent physical links to that graph anchor; 50 distinct
primary Q30 SNP pairs confirm its existing gauge with no contrary pair.

The original primary molecules connecting the insertion locus to the right
A>G SNP at 24,142,287 have **7T and 1T** insertions. Neither is an exact 8T/4T
ALT. Their sequence evidence nonetheless separates the two non-reference
classes. HiPhase assigns both to the correct block, whereas an exact-only
binary call cannot use them and must not interpret the other length as REF.

The new fallback is confined to complementary verified homopolymer insertion
alleles. It compares CIGAR-bounded query sequence to both edited reference
strings in a 16 bp flank window using edlib, requires a unique nearer allele,
and permits at most one shared sequence edit beyond length slippage. A
positive net insertion is required; zero length and sequence ties abstain.
Where the length class itself is unique, the sequence decision must agree.
Distinct upstream SNP-bearing molecules calibrate the two classes, separately
from molecules overlapping the right bridge SNP. Both graph paths retain
independent certificates. The right path's weak 24,402,574–24,422,148 edge is
checked from primary physical SNP calls and its suffix must still pass.

The measured calibration matrix is **[[4,0],[0,5]]**, with association p=0.00390625.
Both class-specific bridge witnesses unanimously require a gauge reversal.
Their augmented one-sided 95% Wilson bounds are 0.4163759 and 0.4614185;
the coherent wrong-parity bound is their product, **0.1921235 <=0.20**. Each
includes the conservative 1% gauge error and its class's worst measured
insert/base/placement/mapping call error (0.0028973 and 0.100302). A contrary
molecule, a missing class, a tiny calibration or an excessive joint bound
rejects the union. This certificate requires **both** independent class
witnesses; it does not relax the existing sparse exact-ALT bridge's rules.
`audit_bridge.py` reproduces these measurements without loading parental truth.

A second defect prevented read promotion after a valid union: the final
marker check used `exact_comp_var_site`, whose nucleotide comparator can
collapse distinct graph-walk allele keys at the same position. It could select
an excluded row before the actual phased anchor. It now matches complete
keys and phased anchors, identically to deferred bridge application.

The completed union can materialize physically agreeing repeat-class reads.
New assignments require MAPQ30, inserted bases >=Q10, query window endpoints
>=Q20, summed call error <=20%, and no contrary phased profile locus.
Existing MSA rescues can enter only their own original core and retain their
inherited HP; sequence evidence must agree. Their accepted low-quality bases
can remain, with unknown qualities still rejected. Both ordinary phase-tag
storage and the separate `gap_haps` / `gap_phase_sets` storage are checked.
Existing core assignments, candidate genotypes, allele sequences, counts and
filters are preserved. Parental truth is evaluation evidence only.

## Owning-chunk result and permanent regression

The native 24,000,001–25,000,000 owning replay retains all 4,061 output reads
and all 1,234 variant records. Correct reads increase **3,922→3,926**, discordant
reads remain **2**, and every prior correct/tagged assignment is retained.
The four new correct assignments are precisely the reads HiPhase tags that
pgphase previously missed. All 49 correct rescue reads overlapping the seam
now enter the connected core.

| Original primary seam reads | Before pgphase | Fixed pgphase | HiPhase |
|---|---:|---:|---:|
| Scorable, abstentions included | 115 | 115 | 115 |
| Correct | 89 | 93 | 93 |
| Discordant | 0 | 0 | 0 |
| Unphased | 26 | 22 | 22 |
| Largest connected correct core | 29 | 93 | 93 |
| Correct / all | 77.39% | 80.87% | 80.87% |

The closure is in the panel, strict certification and HiPhase manifests,
required-marker list and read floors. The existing
`gap_complementary_insertion_recall_cannot_invert_its_snp_flanks` regression
now requires the connection, 93 total/core correct reads, preserved insertion
parity, and disjoint parental flank agreement. No existing expectation is
weakened. The old binary fails 12 assertions; the final focused regression
passes 1,794 assertions. The in-memory adapter tests cover length ties,
non-reference classes, sparse calibration, contradictory witnesses, reversed
gauges and the 20% joint bound.

The owning pipeline runs in about 9.8 seconds; its repeated cached check takes
about 0.38 seconds while rerunning every assertion. Native owner audits reuse
the already geometry-verified HiPhase measurements from the target ranking;
the chromosome audit independently rechecks all original competitor alignments.

## Full chromosome result

The frozen 67-chunk replay confirms the same **93/115 correct reads and
93-read connected core**, with no discordant seam reads. Pgphase now joins
24,103,779–25,855,631 into a **1,751,853 bp** block with 2,106 phased
heterozygous records. Block count decreases **263→262**; N50 remains
**904,351 bp**.

All 256,610 output reads and 64,188 variant records are retained. Correct
reads increase **230,572→230,576**, discordant reads remain **6,768**, and
no previously correct or phased assignment is lost. Alleles, counts and
filters are preserved. All 7,591 HiPhase alignments in the target block
match the original BAM geometry, CIGAR and sequence.

The complete HiPhase block is still split at its separate 2 bp seam,
25,855,631–25,855,633. Its pgphase score is unchanged at 53 correct and
38 discordant among 91 original overlaps. Across the whole HiPhase block,
pgphase now has 7,369 correct, 48 discordant and 174 unphased reads; its
largest connected correct core is 6,821, versus HiPhase's 7,150. This fix
certifies the first seam only.

## Validation

The final production binary SHA256 is
`81ad9a79ed25b71494b8e7cf4288ff31f56c98f133bd908af4c5cf924ef1775a`.
The build, unit tests, 1,551 phase-predicate assertions, HiFi/ONT TSV and VCF
goldens, and HiFi one/four-thread determinism checks pass without new warnings.
The complete window suite passes **13,630 assertions in four cases**; the
cache helper's four tests also pass. The independent physical certificate,
owner and chromosome audits pass.

All **219 saved native outputs**, covering 115 independent replay requests,
have matching final outputs. Every previously correct read and every variant
allele, count and filter is preserved, with no increase in discordant reads.
Three output labels change: the 24 Mb owning replay gains four correct reads,
and two overlapping 10 Mb panel replays withdraw the same previously
discordant read, `m84031_231217_034919_s2/112199276/ccs`. That read remains
in the output and now abstains. Each of those replays retains 4,364 correct
reads while discordant reads decrease 46→45. This withdrawal is recorded
explicitly in `panel-audit.json`; phased-tag preservation is false for those
two labels, while correct-read preservation remains true. Existing regression
floors and ceilings are preserved.

The unusually expensive 31.95–33.05 Mb replay took about 315.2 seconds on the
baseline and 303.6 seconds on the fixed binary. The broad check needs fresh
states after a binary change; four concurrent cache warmers completed the
remaining independent requests while the suite retained all its assertions.
Repeated owning checks use the completed state in 0.38 seconds.

## Reproduction

```bash
make -j8
make gap-dev-check
make gap-owner-check GAP='complementary insertion recall'
make unit-tests predicate-tests check
make window-tests

/home/kokyriakidis/micromamba/envs/bench-phasers/bin/python \
  evaluations/2026-10-05-complementary-repeat-insertion/audit_bridge.py
/home/kokyriakidis/micromamba/envs/bench-phasers/bin/python \
  evaluations/2026-10-05-complementary-repeat-insertion/audit_owner.py

evaluations/2026-10-05-complementary-repeat-insertion/replay.sh
/home/kokyriakidis/micromamba/envs/bench-phasers/bin/python \
  evaluations/2026-10-05-complementary-repeat-insertion/audit_block.py
/home/kokyriakidis/micromamba/envs/bench-phasers/bin/python \
  evaluations/2026-10-05-complementary-repeat-insertion/audit_panel.py
```

The full replay uses four threads and the same original BAM, reference,
striped site catalog and indexed GAF as the baseline. Its output goes to
`test_data/tmp_gap_fix72/frozen_final/0`; the baseline is
`test_data/tmp_gap_fix71/frozen_final/0`. Owner outputs are in
`test_data/tmp_gap_fix72/baseline/24` and `final_checked/24`.

Follow-up: the separate 25,855,631–25,855,633 seam described above is now addressed
by the [physical graph-SNP contradiction fix](../2026-10-05-repeat-copy-seam/README.md).
The measurements here remain the baseline before that follow-up.
