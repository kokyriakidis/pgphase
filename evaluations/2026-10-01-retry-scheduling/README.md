# Preserve later recovery requests without accepting unsafe source joins

## Retained changes

Recovery previously asked a grouped BAM matrix for its first retry request.
An earlier unsuccessful focused attempt could suppress a later, independently
supported original-MSA dropout request. Recovery now records those original
requests before CIGAR backfill, tries the established selection first, and
tries the additional requests when that selection fails certification.
Each attempted seam has a distinct matrix dump prefix.

One accepted focused source owns its seam; the original matrix serves the
remaining seams. Established group admission, local retry mode, edge ownership
and grouped retry behavior are preserved. Simultaneously importing overlapping
focused contexts is not retained.

An additional dropout request may use distinct, already phased complementary
MSA alleles of the same indel type. They must share a positive source PS and
pass the existing diploid depth/fraction and independent local/crossing dropout
checks. This is explicit fallback admission; it does not change the established
first choice or turn a different deletion length into REF/ALT evidence.

An additional source must pass the existing key-preservation and complete
observation-path checks, plus an independent physical clean-SNP chain check.
From the last clean SNP at or before the seam to the first at or beyond its
right end, every consecutive SNP cut must pass the existing Q30 certificate:
two independent source-assigned molecules, known MAPQ30/Q30, no physical
contradiction, and combined wrong-parity bound at most 0.01. Intermediate SNPs
permit a chain without one molecule spanning the entire gap. No new alignment,
allele merge, threshold, truth input, competitor input or fixture condition is
introduced into production.

## 57,854,341–57,866,713: missing evidence and the attempted solve

The owning 57–58 Mb chunk groups this seam with earlier seams. Before the fix,
only the earlier failed focused trial was written. The corrected scheduler
also reaches this complementary one-/two-base deletion pair.

The new trial retains both distinct heterozygous deletion rows and recovers
observations, but its right SNP remains in another source PS. The one-base
row has nine callable crossing pairs: five cross versus four same. The
two-base row has four: three same versus one cross. These are conflicting
relations, not an independently certified bridge. The trial is rejected.
`target-trial-rows.tsv` preserves the separate keys and genotypes.

Two reads that HiPhase calls two-base ALT remain missing in pgphase's MSA
matrix. The original BAM contains a three-base deletion on
`m84031_231217_062403_s3/131269058/ccs` and a four-base deletion on
`m84031_231217_034919_s2/10358603/ccs`, at MAPQ60. Both production exact and
reference-equivalent deletion callers abstain for the one-/two-base targets.
The read is present; the physical allele length differs. HiPhase assigns it
to a nearby allele, while pgphase's current exact/consensus representation
leaves it unknown. This audit does not establish which likelihood assignment
is biologically correct. `boundary-call-audit.tsv` and its C++ caller preserve
the evidence; synthetic tests protect the different lengths from being
silently converted to the target ALT or REF.

## Rejected experiments

- Trying all raw and post-backfill seams and accepting multiple focused
  sources added 77 phased reads, 44 correct and 33 discordant. Concordance
  fell from 96.972176% to 96.959246%. It joined a tracked split control and
  failed four existing assertions across three cases, including the previously
  closed 15.351 Mb gap's local accuracy gate. This version is removed.
- Preserving the first choice and accepting an uncertified later source added
  17 phased reads with only one net additional correct assignment and 16
  additional discordant. It failed the established 99% owning-window gate
  for 61,757,551–61,773,799. Its new source brought seven extra rows near
  61.818–61.823 Mb and merged surrounding blocks. Retaining old keys and an
  MSA observation path was insufficient. The physical SNP-chain certificate
  rejects this source and restores that gate.

No existing expectation, required-site gate, accuracy floor or discordance
ceiling is weakened. No new tracked closure is claimed from these trials.

## Accepted default chromosome result

Same annotated BAM, graph catalog, GAF, reference and defaults, eight threads.
Baseline: `/tmp/pgphase-gap-next14/final`; retained run:
`/tmp/pgphase-gap-next15/certified`. Runtime: 186.44 seconds concurrent with
native panel replays. Production SHA256:
`6855d83085f088e53135ce716b4c8e842532012474cdab3276bdee85483008dc`.

Every primary read name and HP/PS pair is identical to baseline. Every VCF
key, genotype and PS is identical. No tracked gap changes state.

| Metric | Baseline and retained code |
|---|---:|
| Phased / truth-scored reads | 237,101 |
| Truth-correct / discordant reads | 229,922 / 7,179 |
| Per-read-PS concordance | 96.972176% |
| Read phase sets | 682 |
| VCF keys / blocks | 62,811 / 338 |
| VCF span N50 | 739,888 bp |
| Full-chr20 tracked spans | 76/93 |

**Thirteen competitor targets and four split controls remain open.**
The scheduler now attempts the missing recovery; it does not claim that the
conflicting observations establish a correct join. Competitors were not rerun.

## Regression verification

The new owning 57–58 Mb scheduling case fails two assertions on the accepted
pre-fix binary: only one trial is present and the later pair is absent. The
retained implementation passes all 12 assertions, preserving separate
complementary genotypes and rejecting the conflicting join.

Predicate tests pass 684 assertions in 31 cases, including explicit
complementary fallback admission, unchanged default admission, different-PS,
wrong-polarity, homozygous/identical-key rejection, and three-/four-base
physical deletion checks. Shared standalone units pass. Port parity passes
27 assertions in nine cases; original longcallD C parity passes 168,696 in
seven. The pre-existing 61 Mb accuracy regression rejects the uncertified
fallback and passes with the physical chain certificate.

All 105 established native replay commands completed on the retained SHA
with four workers in 235.71 seconds. The new matrix regression is run natively
while the remaining Catch2 assertions use those newly produced outputs after
verifying exact arguments and SHA. Reproduction:

```bash
make predicate-tests
LD_LIBRARY_PATH=/usr/lib/x86_64-linux-gnu ./test_gap_windows '[retry-scheduling]'
make window-tests
```

The committed panel retains all 93 windows and every established expectation.
The scheduling case extends the owning-chunk regression coverage; this change
closes no new gap requiring a new positive panel expectation.

Final full window suite: **4,118 assertions in 54 cases pass**. There are
106 fresh native commands in total, including the new matrix regression
(16.28 s). `validation.json` and `panel.log` record the final gates.
