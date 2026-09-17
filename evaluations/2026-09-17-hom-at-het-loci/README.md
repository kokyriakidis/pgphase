# Three in-gap loci the competitor calls heterozygous and we call homozygous

From the in-gap classification: of the 18 heterozygotes hiphase phases inside the
panel gaps, three land on a record of ours categorised `NOISY_CAND_HOM`. They
have **three different causes**, measured against read truth.

| locus | truth composition | our record | cause |
|---|---|---|---|
| `55,883,019` | MAT -6 (26 reads), PAT -4 (27), **no reference** | leftover 4 bp deletion, 1/36, AF 0.973 | a duplicate of a locus already merged correctly |
| `5,339,363` | MAT +38 (21), PAT +34 (21), **no reference** | -5 deletion, 0/40, AF 1 | the emitted alleles are not the locus' alleles |
| `12,735,894` | MAT -1 (36), PAT **reference** (22) | -1 deletion, **DP 8** of 75, 0/8, AF 1 | 22 reference reads receive no observation |

## Fixed: the duplicate at 55,883,019

The merge had already produced the right record -- `AATATAT -> AAT,A` at `1|2`,
the same REF and the same two alleles hiphase emits as `A,AAT` at `2|1`, with
depths 29/33 against its 31/29. But a **third** co-located deletion survived
alongside it, and it was emitted too:

```
before   55,883,019  REF=AATAT    ALT=A        1|1:37:1,36:0.973     <- spurious
         55,883,019  REF=AATATAT  ALT=AAT,A    1|2:62:0,29,33        <- correct
after    55,883,019  REF=AATATAT  ALT=AAT,A    1|2:62:0,29,33
```

Two things kept it alive. The merge required **both** co-located records to be
`NoisyCandHet`, but a homozygous verdict on one of a co-located pair is the
*symptom* -- that record scores the other haplotype's reads against its own
allele, which is how it reached 36 alt against 1 reference at allele fraction
0.973 -- so the precondition refused exactly the loci that need merging. And the
merge runs per noisy region, so a record produced by a different region or an
earlier pass is out of its reach entirely.

The fix is therefore in two parts: one heterozygous verdict between the pair is
enough to merge, and a chunk-level pass after `collect_noisy_vars_step4` demotes
any deletion record at the same position whose length is already one of a merged
record's alleles, so the prune drops it. The allele is not lost -- it is allele 1
or 2 of the record that remains.

Measured: stock defaults **byte-identical** on the panel (0 concordant to
discordant, 0 records lost or gained), `--retry-unphased-with-bam` 4 of 6 spanned
at 99.49% with 0 concordant to discordant.

## Located, not fixed: the other two

**`5,339,363` -- the emitted alleles are not the locus' alleles.** Truth is +38
maternal against +34 paternal over 65 reads with nothing at reference; we emit a
5 bp deletion and, 6 bp away, a 14 bp insertion. Both records therefore have no
reference read and an allele fraction of 1, and the deletion is classified
homozygous. This is the allele-set completeness defect the retrieval audit
measured across the panel (15 of 24 het indel loci hide a mode carried by 5 or
more reads); no representation fix reaches it, because the allele that should be
there was never proposed.

**`12,735,894` -- the reference reads are not counted.** Here a reference
haplotype does exist: 22 of 75 reads sit at reference length against 36 at -1.
We record DP 8 with 0 reference, so the allele fraction is 1 and the site
classifies homozygous. This is the zero-reference counting defect: a read that
covers the window but carries a different event, or the reference, receives no
observation rather than a reference one. A fix was attempted and reverted -- the
per-cluster pass assigns reference to any read whose cluster consensus lacks the
variant, which is itself false at a multiallelic locus, and a unit test asserts
that a third allele must be recorded as unknown rather than as an invented
reference vote.

## Fixed: the verdict that blocked its own correction, at 12,735,894

`add_msa_site_observations` exists to compose the observations the MSA clustering
never made, raising a site's depth toward its true coverage. It accepted only
`NoisyCandHet` records -- and a record whose **missing reference observations are
exactly what made it read homozygous** is `NoisyCandHom`, so it was skipped and
its counts could never be corrected. A probe at the locus printed nothing at all,
which is how the circularity showed up.

Truth there is 22 of 75 reads at reference length against 36 carrying the
deletion -- a reference haplotype genuinely exists, so this is not the
multiallelic case. The record reported:

```
before   DP=8    0 ref / 8 alt    AF=1.000    NOISY_CAND_HOM
after    DP=61   22 ref / 39 alt  AF=0.639    NOISY_CAND_HOM
```

Letting the pass reach homozygous records corrects the counts. Across the six
panel windows it corrects **9 records**, one of them from `DP 7` to `DP 72`
(`48,210,485`, ref/alt 1/6 -> 19/53).

| arm | spanned | concordance | c -> d |
|---|---:|---:|---:|
| stock defaults | 0 of 6 | 99.68% | **0** (byte-identical) |
| `--retry-unphased-with-bam` | **4 of 6** | **99.49%** | **0** |

### The promotion is a separate question, and it fails

Re-deriving the het/hom call from the corrected counts does finish the job at
this locus -- the record becomes `NOISY_CAND_HET` and is emitted
`12,735,894 TA>T 0|1`, matching hiphase exactly and matching truth. But as a
blanket rule it admits enough additional noisy heterozygotes to take the panel
from 0 concordant-to-discordant reads to **380**, at 83.87% concordance with the
retry, and it loses 3 concordant tags. Reverted.

That is the same wall every broad admission of this class has hit today, and it
is consistent with the census: the noisy class in these intervals contains both
the informative sites and the bridges, and no threshold tested so far separates
them. So the counts ship corrected -- which is what every allele-fraction-gated
consumer downstream reads, including the graph-only classifier that re-derives
the category -- and the verdict does not change until the bridge problem has a
gate.
