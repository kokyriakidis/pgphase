# Close the two-base repeat-copy seam

Target: chr20:25,855,631–25,855,633, the remaining seam inside HiPhase's
third-largest block, 24,103,779–25,944,471 (1,840,693 bp).
The earlier 24,121,713–24,131,707 closure remains independently certified.

## Cause and production change

This is a phase-set boundary between adjacent variants, not a two-base interval
that lacks spanning reads. There are 91 original truth-scorable primary overlaps.
The left terminal graph T>C at 25,855,631 is a false SNP inside a homopolymer
catalog allele: four MAPQ30-or-better graph ALT molecules delete the SNP in their
original CIGAR, with Q30 REF bases on each side. A physical REF witness validates
the graph REF class. No credible physical ALT supports the SNP.

Removing that anchor alone exposes another false graph SNP, 25,853,217 G>A.
Both sites cover exactly the same 20 graph molecules, but their allele
partitions differ. Four independent graph ALT molecules physically carry Q30 G;
the mapping-plus-base error product is 0.000031971816. Five high-quality witnesses
validate its graph REF class. All 20 molecules are maternal according to the
held-out truth diagnostic. Neither production eligibility nor the physical
certificate uses parental truth.

These copy-specific anchors restrict BAM recovery to a region where REF calls
from different repeat copies dominate the true shared haplotype SNPs. The fix
sets both physically contradicted graph candidates to the existing nonanchor
bitmask and repeats the initial graph solve before recovery. Original counts,
categories and observations remain available. The existing BAM recovery then
uses the valid SNPs to connect the blocks and correct the read assignments.

Exclusion requires a terminal homopolymer SNP, at least two independent deleted
ALT molecules with MAPQ30 and Q30 matching flanks, a physical REF-class witness,
and a reference-only partner on the identical graph molecule set. The partner
requires at least two Q30 contrary molecules, a MAPQ30 REF-class witness, and
mapping-plus-base error product <= 0.001. Known positive lower MAPQ can support
this partner only with individual error < 0.5. Any Q30 physical ALT within the
cohort, or MAPQ30/Q30 ALT outside it, vetoes the certificate. Nonprimary,
duplicate, QC-failed and unknown-quality evidence abstains.

A global `--anchor-af-margin 0.2 --emit-nonanchor-hets` experiment also closes the
seam, but the shipped change uses the physical certificate above and retains
the existing global thresholds. Demoting only the deleted terminal SNP fails:
the reference-only companion becomes a replacement false anchor.

## Native owner and regression

The native 25,000,001–26,000,000 replay improves 3979 correct /74 discordant /
43 unphased to 4015 /38 /43 among 4096 truth-scorable output reads. No old correct
or phased assignment is lost. At the seam, 53 correct /38 discordant becomes
88 correct /3 discordant /0 unphased among the same 91 original overlaps;
all 88 correct reads enter one connected core (96.70%). HiPhase has 20 correct /
0 discordant /71 unphased, with core 20, on these identical original alignments.
The remaining three wrong assignments are explicitly retained in the score.

All 1154 retained VCF records preserve alleles, counts, filters and genotypes.
The false terminal T>C is withdrawn. Recovery adds 297 calls within
25,844,279–25,855,615, yielding 1451 records. Disjoint parental flanks establish
the same orientation. One- and four-thread candidate TSVs, variant evidence
and read assignments are identical.

The gap is committed in the panel, native-owner override, HiPhase and strict
certification manifests, required marker list and read floors. The regression
passes 2013 assertions on the fixed binary and fails 12 on the baseline. Its
cached rerun takes 0.96 seconds. Existing floors and ceilings are retained;
the graph TOTAL span count increases 103→104 for this added gap.

In-memory tests cover missing/ALT-vetoed evidence, inadequate independent
molecules, absent REF-class witnesses and excessive mapping/base error.

## Full chromosome and HiPhase comparison

The frozen 67-chunk, four-thread chr20 replay joins the entire target into
PS24103779: 24,103,779–25,964,482, **1,860,704 bp**, with 2631 phased
heterozygous records. This extends 20,011 bp past HiPhase's target endpoint.
Block count falls 262→261; **N50 remains 904,351 bp**.

On all 7591 original truth-scorable primary overlaps of HiPhase's target,
pgphase improves 7369 correct /48 discordant /174 unphased to 7405 /12 /174.
Its connected correct core improves 6821→7380, exceeding HiPhase's 7150.
All 7591 competitor alignments match the original geometry, CIGAR and sequence.
Pgphase has 97.55% correct across all original overlaps, versus 94.19% for
HiPhase. HiPhase has fewer wrong assignments (one), while abstaining on 440
reads; this comparison counts abstentions in the denominator. Disjoint
parental flanks agree on the joined phase-set orientation.

The previous insertion seam retains 93 correct /0 discordant /22 unphased
and core 93 among 115 original overlaps. Chromosome-wide, all 256610 output
reads remain; correct 230576→230612, discordant 6768→6732, unphased 19266
unchanged. No old correct or phased assignment is lost. All retained variant
alleles, counts and filters remain; 64188→64484 emitted rows reflects 297
new recovery calls and withdrawal of the one physically contradicted SNP.

Final binary SHA256:
`7e50015b328beeed83607d502a1e6d6a928d2fea8dbbfb77328c3e5610e76be7`.

Final build and unit tests pass without new warnings. The phasing suite passes
1551 assertions in 47 cases; HiFi/ONT TSV and VCF goldens and HiFi one/four-thread
determinism pass. The full window suite passes 15161 assertions in 4 cases,
including four cache-helper tests. All 220 saved task-start native outputs
(115 independent requests) preserve correct and phased assignments, with no
increase in discordant reads and retained variant evidence unchanged. Four
concurrent cache warmers use the existing locks and fingerprint rules; every
suite assertion runs. See `validation.json` and `panel-audit.json`.

## Reproduction

Use the bench-phasers Python environment with pysam for the audits:

```bash
bash evaluations/2026-10-05-repeat-copy-seam/replay.sh
/home/kokyriakidis/micromamba/envs/bench-phasers/bin/python evaluations/2026-10-05-repeat-copy-seam/audit_physical.py
/home/kokyriakidis/micromamba/envs/bench-phasers/bin/python evaluations/2026-10-05-repeat-copy-seam/audit_block.py
LD_LIBRARY_PATH=/usr/lib/x86_64-linux-gnu make gap-owner-check GAP=graph-snp-contradiction
LD_LIBRARY_PATH=/usr/lib/x86_64-linux-gnu make unit-tests predicate-tests check window-tests
/home/kokyriakidis/micromamba/envs/bench-phasers/bin/python evaluations/2026-10-05-repeat-copy-seam/audit_panel.py
```

Baseline executable: `test_data/tmp_gap_fix73/pgphase-before`.
Baseline full chromosome: `test_data/tmp_gap_fix72/frozen_final/0`.
Final executable: `test_data/tmp_gap_fix73/pgphase-final`.
Final full chromosome: `test_data/tmp_gap_fix73/frozen_final/0`.
The physical audit uses the baseline owner matrix to preserve the original
pre-exclusion graph cohorts. Panel preservation uses the task-start snapshot
`/tmp/pgphase72-final` and requires every audited replay to match the final
production binary fingerprint.
