# Representation repair: step 2

## Scope

Add a shared, reference-validated identity for repeat-shifted alleles and compound
replacements. Integrate it first into the whole-chunk BAM observation overlay.
The starting executable is the verified step-1 binary
`088b5c56e17d04828af0ccbfede26b9ca3ff65476a62f796066af7a6ac63eb9d`;
its source snapshot and executable remain in `test_data/tmp_representation_step2/`.
The production executable for this evaluation is
`5c47ca57f979344cce86d167b18225b20d9ffc7accc47a99dd602c7236694ee7`.

`normalize_candidate_identity` in `src/allele_identity.*` validates REF against
reference bases, compares uppercase ACGT, removes common context, and left-aligns
an edit by repeatedly extending both alleles on the left when one is empty, then trimming the
shared suffix. This rotates tandem motifs rather than comparing insertion length
or proximity. The identity retains the complete changed ALT and reference length,
including compound replacements. Unknown reference context, symbolic/ambiguous
alleles and unchanged alleles supply no new alias. The reference callback uses
the existing worker-local reference cache; there is no shared faidx access or
arbitrary flank cutoff.

The overlay retains its exact raw and selected-sequence indexes and adds a
collision-aware normalized index. A normalized graph identity must be unique;
an exact match cannot bypass normalized ambiguity or identify a different row.
A newly normalized match also requires a unique BAM source description. Graph
aliases require an actual REF/selected-ALT pair, with original allele zero as
REF. No ALT1/ALT2 pair is relabeled as REF/ALT by this stage.

Candidate physical coordinates are kept separately from their canonical alias.
Each transferred shifted observation must cover both original edits and their
outside flanks, and the read must already have an independent graph call at
that candidate. A normalized identity alone cannot certify a missing source
REF class; such calls await complete-contrast/confidence handling. The overlay fills the BAM channel, retaining an independent
conflicting graph call. It does not alter candidate genotypes/counts or import
source HP/PS as a graph gauge. Missing normalization context preserves the
established exact-matching rule. Seam transfer remains exact: enabling normalized
seam anchors is coupled to the next step's molecule deduplication and explicit
source-gauge handling.

## Fast verification

The standalone identity binary exercises 3,328 complete-haplotype fixtures per
run. Each edit is independently applied to the complete reference haplotype;
equivalent resulting sequences must receive the same identity, the canonical
edit must reproduce that sequence, and normalization must be idempotent. Focused
checks cover homopolymer shifts, rotated tandem motifs, padded SNPs, compound
replacements, different equal-length alternatives, reference mismatches,
symbolic alleles, unknown left context, contig-start anchoring, index conflicts,
permanent normalized ambiguity and short-read coverage.

One hundred executions average **2.9 ms**, including process startup. Mutations
that disable repeat shifting, ignore canonical conflicts, ignore read coverage
or bypass the independent-graph-call requirement fail 772, two, one and one
checks respectively. See `fast-checks.json` and
`mutation-checks.json`. Unit tests and the existing 47 phasing predicate cases
(1,551 assertions) pass. The build adds no warnings. HiFi/ONT TSV and VCF golden
checks and HiFi t1/t4 determinism pass (`gates.log`). All four replay-cache tests
pass (`cache-tests.log`).

`check_normalization.py` independently unanchors `bcftools norm` results and
compares them with the production helper for all **64,483** frozen chromosome
records. All identities agree, including all 504 realigned records; the batch
helper takes about 0.28 s. The old-record tag maps reordered normalized records
to their source rather than assuming stable output order. The original FASTA
has lowercase sequence. Uppercasing both alleles before comparison reveals
**19 duplicate pairs**, including nine missed by the previous case-sensitive
normalization audit. The original ten remain. Positions and identities are
preserved in `normalization-checks.json`. All 19 pairs agree on GT/PS; eight
pairs differ in AD, reinforcing that their overlapping molecule depths must not
be summed. This stage recognizes those identities;
it does not yet merge their candidate rows or sum their depths.

## Owning matrix verification

The 5.31 Mb owning regression passes 351 assertions. An additional full
5,000,001–6,000,000 owning replay dumps the independent graph/BAM channels and
is compared with the frozen same-region investigation matrix. The production
matcher finds 173 additional normalized aliases, transfers 9,870 BAM calls,
and rejects 1,446 unsupported observations: three lack coverage and 1,443
previously covered calls lack an independent graph observation. All 1,825
candidate states and all 4,256 incoming read states remain identical. Graph
observations and SNP quality certificates remain identical; zero quality
padding is treated as absent, rather than as a new certificate. All 9,870
changes add a BAM observation to a previously absent BAM slot; no effective
working allele changes, no known BAM call is replaced, and no existing BAM call
is removed. The 697 added graph/BAM disagreements remain separate from the working graph call.
`audit_matrix.py` and `matrix-checks.json` preserve this check. The owning
output also preserves all 4,256 primary HP/PS pairs and their parental status:
4,001 correct, 23 discordant and 232 unphased. Candidate TSV and VCF records
match (`owner-output-checks.json`).

## Remaining representation work

Merge equivalent descriptions with one observation per original molecule and
retain disagreements/provenance; then preserve complete sample ALT1/ALT2
contrasts and MSA consensus windows. A normalized alias proves equivalent
alleles, not independent read confidence or source orientation. These stages
must precede broader normalized seam transfer and the joint graph/BAM solve.
No new gap closure is claimed by the normalization stage.

## Rejected admission and corrected checks

The initial integration passed the >=80% criterion on its new output rescues
(29 correct, four discordant), but failed two stronger existing owning error
ceilings. No expectations were updated. The normalizer proves equivalent
alleles; it does not certify a missing source reference class or independently
phase its reads. Requiring an independent graph call at the matched candidate
fixes both owning regressions with the original bounds. Their final runs pass
1,007 and 333 assertions (`owner-short-final.log`, `owner-deletion-final.log`).
A fast guard regression and a mutation that removes it preserve this lesson.
The rejected executable, full results, competitor comparison and failing tests
are saved separately under `provisional/` and the evaluation work directory.

## Final panel and chromosome verification

All **107 named checks pass**, using the original expectations and read floors.
Cold parallel verification takes 705.4 s; the standard `make window-tests` run
then passes **18,138 assertions in five cases** using the completed replay
states. `window-checks.json`, `standard-windows.log` and `panel-comparison.json`
preserve the results. All **141 native panel measurements** are unchanged.
The HiPhase measurement cache validates unchanged input/indices, truth, panel,
helper and competitor identities and output checksum before reuse;
`hiphase-verified.tsv` is byte-identical to the committed measurements.

The full chromosome preserves all **256,612 primary HP/PS assignments** and
parental classifications: 230,965 correct, 6,620 discordant and 19,027 unphased.
All 82,287 candidate rows are byte-identical and all 64,483 VCF data records
match. The 250 variant blocks retain **N50 955,496 bp** and largest block
**3,012,193 bp**. No rescue tags or connected-core assignments change.
`full-comparison.json` records the comparison against the frozen step-1 output.
The full and native HiPhase comparisons for all 141 windows also match the
preceding investigation exactly (`panel-audit.json`). This preserves existing
closures and deficits: 59 of the 118 windows where HiPhase reaches 80% pass the
full span/80%/total/core contract, and 59 remain deficient. Normalized identity
recognition does not claim to close them.

## Reproduction

```bash
make -j$(nproc)
make unit-tests gap-dev-check
python3 evaluations/2026-10-06-representation-step2/check_normalization.py
make gap-owner-check GAP="a recovered short insertion retains"
make gap-owner-check GAP="independent BAM SNP pairs certify"
make check
python3 evaluations/2026-10-06-representation-step1/check_windows.py \
  --binary ./pgphase --out test_data/tmp_representation_step2/windows-final
make window-tests
bash evaluations/2026-10-06-shared-insertion-source/replay.sh \
  test_data/tmp_representation_step2/full-current
/home/kokyriakidis/micromamba/envs/bench-phasers/bin/python \
  evaluations/2026-10-06-representation-step2/audit_full.py
/home/kokyriakidis/micromamba/envs/bench-phasers/bin/python \
  evaluations/2026-10-06-representation-step2/audit_panel.py
```

The owning-matrix audit requires the `--phase-matrix-dump` replay of
CHM13#0#chr20:5,000,001–6,000,000, saved in
`test_data/tmp_representation_step2/matrix5/`; run `audit_matrix.py` after it.
The frozen before output and matrix are retained as evaluation inputs.
