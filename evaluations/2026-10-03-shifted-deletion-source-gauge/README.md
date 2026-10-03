# Shifted-deletion recovery and source-gauge preservation

## Scope

Investigate two noncentromeric HiPhase-supported nominations with the largest
modeled NG50 gain: chr20:57,085,410–57,104,654 and
61,738,239–61,747,506. Joining either complete neighboring VCF block pair
would raise span NG50 from 643,699 to 672,998 bp. Ranking does not certify a
haplotype connection. Neither target is closed by the retained correction.

## Defect and correction

An exact-position CIGAR lookup can call REF when the same deletion lies beyond
that row's footprint. The existing missing-deletion recovery returned this
REF immediately, without checking reference-edit equivalence. At a simple
phased MSA deletion, examine the existing Q30/MAPQ30 reference-edit certificate
before filling an exact-position REF. An independent clean SNP at least 100 bp
away must agree with the source block's ALT gauge to admit the shifted ALT.
If a callable SNP contradicts that gauge, leave the new observation unknown.

A missing SNP gauge must preserve the source's ALT-absence contrast, as existing
insertion recovery already does. At 48.929 Mb, seven reads carry the equivalent
deletion but no callable independent clean source SNP. Removing their zero
contrasts displaces the accepted source solve, loses 195 owning-chunk variant
keys and reopens a protected gap. Preserve them. Separate colocated MSA
alternatives, existing observations, candidate keys, discovery counts,
genotypes and phase labels retain their source semantics. Homopolymer admission
stays unchanged. No realignment or merge is added. Truth and competitor calls
are used only for evaluation and never by production.

Synthetic regressions reproduce the exact REF/shifted ALT conflict, an agreeing
SNP, a contradictory SNP, an absent SNP gauge and separate colocated MSA alleles.
The starting implementation fails three assertions; the correction passes all
1,031 predicate assertions across 42 cases. Add a 61–62 Mb owning-chunk test
that protects 3,452 scored / 3,423 correct / at most 29 discordant assignments.
A future join is permitted, but clean boundary SNPs must keep their supported
parental relation when they acquire a shared PS. Existing span, coverage and
accuracy expectations are unchanged.

## Why the targets remain open

At 57.085 Mb, only one MAPQ60 primary read physically spans the boundary sites.
Its CIGAR deletes five bases from an A run. The recovery source represents two
MSA consensus alleles, +3 A and −2 A. Native DV has the same +3/−2 representation:
this is not proof of a missing candidate row. The singleton read's −5 event
needs a reliable consensus-allele assignment, rather than promotion into a new
whole-block gauge. The current exact physical certificate declines it.

At 61.738 Mb, the isolated source deletion is homopolymeric and its partial-read
MSA observations omit useful spanning molecules. Shifted two-base CIGAR
placements can be edit-equivalent to this candidate. But three spanning reads
carry the deletion and a clean left SNP whose source haplotype contradicts the
source deletion gauge. Simply adding those molecules under the frozen source
orientation is unsafe. The fix checks this conflict rather than forcing a join.

`local-57.json` and `local-61.json` compare gap overlaps using majority truth
orientation separately per local PS. They use the saved October 1 same-BAM
HiPhase outputs, not a new competitor run. These local metrics do not substitute
for chromosome scoring, which derives each PS orientation using all scored
reads. Whole-block continuity and whole-chunk error limits remain separate gates.

## Rejected trials

- Admit missing observations on isolated phased homopolymer deletions after
  the source solve, with independent Q30/MAPQ30 clean-SNP agreement for REF or
  ALT. Local target outputs are unchanged, but the full chromosome goes from
  230,072 to 229,542 correct and 7,127 to 7,707 discordant reads. Two protected
  panel joins reopen (48.929 and 56.064 Mb); spans fall 87 to 85. Reject.
- Remove that admission but turn every edit-equivalent ALT without a callable
  source SNP into unknown. Full chromosome loses 202 variant keys and six
  previously correct assignments, reduces scored reads by nine, and reopens
  the same two joins. Reject despite aggregate conditional accuracy increasing.
- Admit isolated homopolymer physical identities without independent source
  SNP agreement after solving. The local 61 Mb gap gains five correct reads,
  but its owning chunk adds ten errors and an extra VCF block. Reject.
- Add those physical identities before the noisy k-means solve. The 61 Mb
  owning chunk loses 209 correct assignments and adds 236 errors. Reject.

The JSON artifacts record these measured failures. Interrupted or abandoned
runs are not counted as successful full evaluations. No expectation is lowered
to accept a trial.

## Validation

Full chromosome output is unchanged: 237,199 scored / 230,072 correct /
7,127 discordant reads (96.995350%), 659 truth-scored read phase sets,
63,630 VCF variant keys, 332 VCF phase blocks, N50 774,189 bp and NG50
643,699 bp. No HP/PS tuple or VCF sample row changes, no key disappears,
no old phased SNP is lost and no source block acquires a mixed SNP gauge.
All 87/101 coordinate-panel spans survive. The current two targets stay open.

Build, all standalone unit tests and 1,031 predicate assertions/42 cases pass.
HiFi/ONT TSV and VCF goldens pass, with HiFi one/four-thread determinism.
The fresh full chr20 run takes 300.63 s at eight threads while other validation
jobs execute; this is not a controlled runtime comparison. All 105 fresh native
panel requests complete in 401.92 s with four concurrent workers. Rescoring
verifies the final executable hash, normalized CLI, input stats and every cached
output hash. Requests absent from the manifest execute natively. The complete
window suite passes **7,113 assertions / 70 cases**, including the new owning
61 Mb guard. The coordinate panel remains 101 cases / 87 connected spans;
no closure expectation or accuracy/coverage limit is weakened.

The new 61 Mb owning-chunk guard rejects the before-solve physical-identity
trial on two assertions: correct reads fall to 3,214 and discordant reads rise
to 265. `rejected-owning61-regression.log` retains this counterexample.

Starting executable SHA256:
`aa3b89bb6f6f5a49a4fc72e0c83577bfade869165e60f2c8556a48bf5e0026d4`.
Starting chromosome output: `test_data/tmp_gap_next36/full-preserved/`.
Final executable SHA256:
`244d3e8f39b91c5ac17d9ad3c5c05a28ad42223a9cb12db681c77e0c571ef3a8`.
Final output: `test_data/tmp_gap_next37/full-gauge-preserved/`.

Reproduce the full run from the repository root with default recovery:

```bash
./pgphase collect-graph-variation \
  --ref test_data/chm13v2.0.chr20.renamed.fa \
  --bam test_data/HG002_chr20_hifi_mapped_to_CHM13_chr20_annotated.bam \
  --sites test_data/chr20.sites.striped.vcf.gz \
  --gaf test_data/HG002.chr20.annotated.coord.gaf.gz \
  -t 8 -r 'CHM13#0#chr20' -o /tmp/shifted-deletion/candidates.tsv \
  --phased-vcf-out /tmp/shifted-deletion/phased.vcf \
  --phased-bam-out /tmp/shifted-deletion/phased.bam
make unit-tests predicate-tests window-tests
```

`chromosome-parity.json` and `block-transfer.json` report full-chromosome truth,
allele-gauge and continuity checks. `ranking.json` independently verifies NG50.
