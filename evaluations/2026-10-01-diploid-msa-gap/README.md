# Focused diploid MSA recovery and singleton read evidence

## Reproduced defects

The chr20:15,351,845–15,367,755 target remained split although the existing
unplaced-read MSA could provide a source path. Its original MSA boundary matrix
lost all eight crossing deletion calls despite deep local coverage. Admission
required 20 crossing reads and did not recognize this local/dropout combination.
The seam also lies inside a merged recovery group; focused retry ownership
previously allowed only an edge seam. Its right graph block is a singleton,
which left the focused solve without right-side BAM context.

A compact replay of the saved observation matrix reproduced a second defect:
`joint_het_orientation` retained the noisy heterozygous deletion genotype, but
`upstream_assign_hap` still excluded that homopolymer genotype from the link
list. With that explicit option honored and the existing `link_by_alleles`
mode enabled in the focused solve, the two boundaries share one source PS in
the opposite allele orientation supported by physical read pairs. Ordinary BAM
mode retains its upstream exclusion and passes upstream parity tests.

A third defect became visible after joining: read rescue treated every phased
MSA row as an unconditional singleton. Most newly rescued gap reads observed
only the noisy complementary deletion/insertion locus. A genotype retained for
block connectivity does not certify a reliable one-locus read assignment.

## Correction

- Require a co-located verified opposite-type MSA indel row with biallelic
  heterozygous allele depths. The contrast may have a collapsed genotype;
  requiring it to be already phased would conceal the original genotype bug.
- Admit a focused retry for that mixed MSA boundary only when local and crossing
  cohorts independently show significant majority dropout at known MAPQ30.
  Local depth is at least 20; configured minimum depth applies to both cohorts.
  Exact binomial upper tails use p <= 0.01, with two-endpoint correction for
  crossing dropout. Admission supplies no orientation certificate.
- Allow a middle seam to own the focused source. The original grouped source
  still owns the other seams and excludes the focused interval from transfer.
  Use complete adjacent graph phase-set extents and existing 50 kb padding for
  singleton flanks, all clamped to the owning chunk.
- For this original-indel-dropout retry, use the existing joint genotype and
  allele-pair link modes. Retain every previously phased source key and require
  a complete source path without weak cuts after backfill. Separate MSA rows
  remain separate.
- A noisy MSA indel retained by the focused diploid retry records
  `read_rescue_requires_validation` and needs
  the existing primary-read singleton test or a second independent locus.
  Singleton support requires both haplotypes/alleles, exact association p <=
  0.01, matching orientation, and a 95% Wilson discordance bound <= 15%.
  Rescue-derived labels cannot establish singleton confidence.

No new realignment, truth input, competitor input, coordinate branch or read
name rule is introduced into production phasing.

## Owning-chunk observations

Native 15–16 Mb replay, same BAM/GAF/reference/catalog and parental truth:

| Run | Gap truth reads | Phased | Correct | Concordance |
|---|---:|---:|---:|---:|
| Accepted baseline | 151 | 92 | 88 | 95.65% |
| Focused retry without singleton correction (rejected) | 151 | 121 | 100 | 82.64% |
| Broad retry with general singleton correction (rejected) | 151 | 98 | 89 | 90.82% |
| Narrowed mixed-MSA retry with provenance-specific singleton check | 151 | 92 | 89 | 96.74% |
| HiPhase same-callset diagnostic | 151 | 137 | 118 | 86.13% |

In the final owning-chunk output, the deletion `CT>C` at 15,351,845 has
GT `0|1`; the `A>AATCT` insertion at 15,367,755 has GT `1|0`. The actual
adjacent left SNP `A>T` at 15,329,499 and right SNP `C>T` at 15,380,509
also have GT `1|0`. All four share PS 15,518,376. The distant older source
block remains separately labeled 15,039,543; its original four-site connection
is preserved. Numeric HP orientation is arbitrary, so the regression checks
these relative allele relations rather than a fixed `0|1` label. Local disjoint left
and right flanks retain the same parental orientation; read concordance is
100% and 99.48%, respectively. The owning chunk retains 4,280 scored reads,
with 4,260 correct and 20 discordant versus 4,259 correct and 21 discordant
in the baseline.

An experiment validating independent BAM fallback against original catalog-only
labels made no difference and was removed. The suspect `PS + 1e9` assignments
were graph read rescues, not independent BAM fallback (`PS + 1.5e9`). The
whole-chunk BAM matrix was unchanged.

## Regressions

- The genotype/link unit reproducer fails before the link-list correction;
  explicit joint mode now connects, ordinary BAM mode stays split, depth
  failure stays split, and crossed allele pairs preserve opposite orientation.
- The noisy MSA singleton reproducer fails before the rescue correction. Clean
  diploid support still permits singleton rescue; noisy repeat/nonrepeat rows
  remain unassigned without changing their candidate phase state.
- Admission tests cover local/crossing dropout, configured depth, missing
  profiles, multiallelic calls, MAPQ255/low MAPQ/skipped reads and either endpoint.
- The existing owning-chunk gap regression gains exact boundary GT/PS checks,
  a disjoint-flank parental switch check, and local floors of 92 scored/89
  correct/96% concordance with at most three errors. Its spans expectation is
  changed deliberately from zero to one. Earlier source connection assertions
  remain intact.

The broad trial is rejected: chromosome phased count fell 237,101 -> 235,919
and two protected chromosome gaps (5.31 and 56.064 Mb) reopened. It also joined
the 15.095 Mb control and failed six assertions in three protected test cases.
Its new local 98-read floor was provisional and never accepted. The narrowed
trial preserves the accepted baseline's 92 local phased reads and improves
correct assignments from 88 to 89. All preexisting test floors remain intact;
the new local check is set from that measured baseline coverage and improved
accuracy. Ordinary retry admission, singleton flank geometry and read-rescue
semantics remain on their established paths.

## Final chromosome and validation

| Metric | Accepted baseline | Final |
|---|---:|---:|
| Primary output reads | 256,586 | 256,586 |
| Phased/scored reads | 237,101 | 237,101 |
| Truth-correct reads | 229,921 | 229,922 |
| Discordant reads | 7,180 | 7,179 |
| Read concordance | 96.971755% | 96.972176% |
| Read phase sets | 684 | 682 |
| VCF blocks | 339 | 338 |
| VCF variant keys | 62,798 | 62,811 |
| VCF span N50 | 739,888 bp | 739,888 bp |
| Tracked spans | 75/93 | 76/93 |

Only one read HP changes, from incorrect to correct. Read names and the set of
phased read names remain identical. No original VCF key is lost and no original
VCF genotype changes; 13 source rows are added. The only tracked span change
is 15,351,845–15,367,755. All 75 previously closed gaps and all four split
controls are retained. Thirteen of the fourteen previously open competitor
targets remain open. Phased-read percentage is unchanged at 92.406055%.

- Final binary SHA256: `588e94e129a673c42d7c6b42a8803cda037ac4e85e19413676aef8517ef2d73b`.
- Full chr20 native run: 216.26 s with eight workers.
- Fresh native panel: 105 unique requests, four concurrent workers, 264.91 s.
  Catch2 scores those exact outputs through an argument-and-SHA-checked cache;
  no pipeline result comes from an earlier binary.
- Final window gate: **4,106 assertions / 53 cases pass**.
- Predicate gate: **663 assertions / 31 cases pass**.
- All standalone unit tests pass.
- Port parity: **27 assertions / nine cases pass**.
- Original longcallD C parity: **168,696 assertions / seven cases pass**.
- Build has no new warnings; `git diff --check` passes.

Durable evidence: `chr20-metrics.json`, `chr20-changes.json`,
`local-read-summary.json`, `local-reads.tsv`, `owning-boundaries.vcf`, and
`validation.json`. Local read accuracy is conditional on reads phased by each
tool; HiPhase tags more local reads (137 versus 92). The comparison is to the
existing same-callset targeted diagnostic, not a new chromosome-wide HiPhase
run. No N50 improvement or complete competitor coverage parity is claimed.
