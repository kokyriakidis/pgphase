# Verified SNP evidence for recovered read assignments

## Three reproduced defects

1. Recovery attached the quality of the original CIGAR base to an MSA-derived
   allele without checking that they described the same base. A Q40 REF or
   third-base observation could therefore be advertised as a physical ALT
   certificate. The 17 Mb source audit finds 128 Q30 disagreements in the
   investigated interval; these include homozygous recalled rows and are not
   128 proven wrong phase assignments. The allele remains available as MSA
   evidence, but its original-CIGAR quality is not a physical certificate.
2. Overlapping recovery solves can disagree. At 3,597,791, the first solve
   retains ALT for `m84031_231217_034919_s2/16389926/ccs`, while another source
   calls REF. Its original BAM base is REF at Q40, MAPQ60. Separate observation
   and quality maps retained the first ALT but attached the later REF quality.
   Certify the retained allele, and replace an old quality whenever a replay
   replaces the observation, including with zero. The permanent matrix replay
   reproduces `40 == 0` before the transfer fix and requires no certificate for
   the retained ALT afterward. MSA alleles themselves remain unchanged.
3. Read HP inherited from indel clustering could ignore an admitted, phased
   MSA SNP even when its original BAM base agreed at high quality. Near
   17.50 Mb, most locally incorrect reads carry a Q35/Q40 observation at
   17,499,524. The old post-stitch correction admitted only clean SNPs and
   required two loci, excluding this verified noisy-region SNP.

`bam_snp_observation_quality` now uses the existing physical base caller and
retains quality only for the matching binary REF/ALT observation. Ambiguous
nucleotides, third bases, uncovered/deleted positions, missing quality and
contradictory MSA calls return zero. Both physical singleton-bridge checks and
source quality transfer use this certificate. Neither candidate rows nor the
original MSA allele observations are replaced. Standalone BAM behavior and
its golden outputs remain unchanged.

An imported, source-phased single-base MSA SNP can correct read HP only after
physical links validate its allele gauge against clean phased SNPs of the same
PS. Every covered clean SNP at least 100 bases away contributes; a molecule
with conflicting clean calls abstains and each molecule supplies one vote.
Both allele and haplotype classes must occur. Reuse singleton rescue's exact
one-sided association test (p <= 0.01) and one-sided 95% Wilson discordance
bound <= 15%. A cached result is reused for that site throughout the chunk.
The physical/primary observations must agree at Q30 and known MAPQ30. Any
eligible contradictory SNP vetoes the individual correction. A physically
matching clean call below Q30 can veto an HP override relying on a noisy SNP,
but cannot assign HP by itself. The existing spaced Q30 clean-SNP certificate
is sufficient. Ordinary clean SNPs retain the spaced two-locus rule. MSA witnesses are kept separate from
clean-SNP counters; no genotypes, PS labels, variants or joins are changed.

No truth, competitor calls, fixture coordinates or new alignment method enters
production. Read HP itself is not site-certification evidence. The certificate
uses candidate gauges and paired physical allele observations.

## Rejected controls

- The quality correction alone changes three previously correct tags to
  incorrect ones, preserving all rows, blocks and phased counts. These MSA
  observations can be useful, but the old CIGAR-quality claim was invalid;
  do not mistake removal of that claim for an accuracy gain by itself.
- Allowing every verified MSA SNP without a site gauge check gives a net
  26 correct reads, but reverses 41 previously correct reads. One site near
  15.07 Mb costs 26 correct reads. That unrestricted rule is rejected.
- Using only the nearest clean SNP excludes five usable bridge molecules at
  17.50 Mb. Its 14 conflict-free paired calls cannot meet the existing 15%
  Wilson upper bound, so it loses the intended correction and regresses four
  net correct reads chromosome-wide. Do not relax the statistic: use all
  available covered clean SNPs, yielding 19 coherent physical links in the
  targeted audit. The nearest-only rule is rejected.
- The complete-gauge candidate initially adds one error to the 11.599 Mb
  owning replay (91 against its unchanged ceiling of 90). The affected read
  has a contrary physical clean SNP at Q10, plus a Q35 noisy SNP. Ignoring
  the clean call permits an unnecessary HP reversal. Keep a weak, physically
  matching clean call as a veto against a correction relying on a noisy SNP;
  it cannot assign HP. Two spaced Q30 clean SNPs keep their existing positive
  certificate. The unchanged 90-error ceiling passes with this veto. The
  rejected candidate and failing gate remain recorded.

Controls and actual per-read transitions are preserved in the JSON files.
The owning 15 Mb probe using 15–16 Mb has different input bounds from the
permanent replay's 14.95–16 Mb context; their aggregate read counts must not be
compared as equal-input results. The final panel uses the original exact CLI
argument vectors.

## Local result and regression design

With complete clean-SNP context, the 17.502614–17.521341 Mb connection remains
intact with the same parental orientation. Local truth-scored overlaps are
unchanged at 136; correct assignments increase 109→131 and errors fall 27→5.
Conditional read accuracy is 96.32%. Saved same-BAM HiPhase DV has 138/144
correct (95.83%), but still phases eight more overlapping reads. Accuracy
parity is achieved there; coverage parity is not.

Strengthen the permanent 16–18 Mb regression from 109 correct/27 errors to at
least 131 correct/at most five errors. Its dominant-block floor rises 0.70→0.85;
the exact span, representations and parental majority checks are retained.
Restore the two older 17 Mb owning-score floors to 99% and strengthen the
15.056 Mb owning-score floor to 99%, guarding the ungated counterexample.
No span expectation is weakened; the panel still has 100 coordinates and 86
required connections.

Adapter units exercise exact physical-quality provenance; verified MSA
correction; conflict with clean SNPs; unsupported/reversed gauges; one-sided
and single-molecule evidence; retaining farther anchors when the nearest call
has no quality; weak clean-SNP contradiction; retaining the spaced clean-SNP
certificate; and preserving unassigned reads. The single-SNP correction
reproduction fails before the change (`msa-snp-before.txt`). Existing predicate
and BAM golden gates are retained. The overlapping-source certificate test is
truth independent: it checks the retained allele and its quality, rather than
accepting a haplotype because it happens to match a parental label.

## Final chromosome and panel validation

Final binary SHA256:
`95f94fdd04210a95e295131288ca21398065ac0d0afe1364c7c586fc892e63aa`.
Accepted baseline: `2f0513f66ddcbc2a06f0a1b37d3674b9acb4f6df6c735540c7c8df8af0a3a03e`.

| Measurement | Accepted baseline | Final |
|---|---:|---:|
| Truth-scored phased reads | 237,172 | 237,170 |
| Truth-correct read assignments | 230,022 | 230,043 |
| Discordant assignments | 7,150 | 7,127 |
| Conditional read accuracy | 96.985310% | 96.994983% |
| Read phase sets | 661 | 661 |
| VCF variant keys | 63,630 | 63,630 |
| VCF phase blocks | 333 | 333 |
| VCF span N50, bases | 774,189 | 774,189 |

The final full chromosome has 38 changed HP/PS tags and zero changed VCF rows.
29 previously incorrect reads become correct; seven previously correct reads
become incorrect. Two independent BAM fallback reads become unphased, one
previously correct and one incorrect. The changed SNP HP no longer validates
their old fallback connection, so leave them unassigned rather than force a
join. Net result: 21 more correct assignments and 23 fewer errors. The final
atomic quality transfer leaves the full chromosome HP/PS and VCF results equal
to the complete-gauge candidate with the clean-call veto, while repairing the
observable quality provenance invariant.

Final local 17.50 Mb result is independently rescored from the full chromosome:
131/136 correct, five errors, one block; HiPhase has 138/144 correct, six errors,
one block. Coverage is not yet equal. No new gap closes in this change; the
7.264321–7.280346 and 24.121713–24.131707 Mb cases remain open.

Validation on the final binary:

- Build succeeds with no new warnings.
- `make unit-tests`: all pass.
- `make predicate-tests`: 906 assertions in 40 cases pass.
- `make check`: unchanged HiFi TSV/VCF goldens and thread determinism;
  unchanged ONT TSV/VCF goldens.
- `make window-tests`: 7,038 assertions in 69 cases pass. The added certificate
  replay fails before the atomic transfer fix and passes after. The permanent
  panel retains 100 gap coordinates and 86 required connections. Span equality,
  parental orientation and the existing 11.599 Mb error ceiling remain intact.

The panel first executes 105 fresh native requests at the exact recorded input
bounds, arguments and thread counts. Rescoring uses binary, argument, input-stat
and output-hash checks; additional requests absent from the cache run natively,
including the new matrix regression. Accepted-baseline outputs are not reused.
The full final run takes 297.21 seconds at eight threads; the concurrent native
panel takes 394.91 seconds with four workers. These are validation timings,
not a controlled runtime comparison. The cache records native file paths under
`test_data/tmp_gap_next31/atomic-native/`; full output is under
`test_data/tmp_gap_next27/atomic-quality-snp-gauge/`.
