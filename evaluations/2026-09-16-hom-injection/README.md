# Injecting homozygous graph sites, and what actually happens at one

## What injection adds

Measured over the six panel windows, against `collect-bam-variation` run
standalone so the only difference is the catalog claim. Injection adds 58
candidates the alignment channel does not have:

| category of the added candidate | n | can it phase? |
|---|---:|---|
| `CLEAN_HOM` | **40** | no |
| `REP_HET_INDEL` | 14 | no, demoted by construction |
| `CLEAN_HET_INDEL` | 4 | yes |

So **54 of the 58 additions carry no phase information**, and 69% of them are
homozygous. Only `CleanHet{Snp,Indel}` enter het k-means.

Checked against read truth, the homozygous verdicts are mostly right: the
ALT-carrying reads split between the parents roughly evenly (`5,290,557` MAT 36 /
PAT 45, `5,327,179` MAT 26 / PAT 41, `24,083,214` MAT 29 / PAT 32), which is what
a homozygous ALT site looks like. One is wrong -- `12,795,753` at AF 0.803, just
past `max_af = 0.80`, has 25 ALT reads split MAT 2 / PAT 23, i.e. a real het.

## Dropping them: tried, and it removes the wrong record

`classify_graph_only_candidates` assigns `CleanHom` when AF exceeds `max_af`, and
`CleanHom` survives the prune pass, so these records reach the output. Folding
that verdict into `LowCoverage` instead (so the prune pass removes it) is
behaviourally inert for phasing -- the panel is byte-identical, 0 gate violations,
same blocks, same 2,794 tags at 99.68% -- and removes 5-12 candidates and 10-30
VCF records per window.

But two of the removed records were emitted `0|1`, and one of them is real:
`55,896,395 CTTTTTTTTTTTTTTTT>C`, 28 ALT reads **all paternal**, purity 1.000.
Reverted.

## Why a homozygous candidate is emitted as a het

This is the mechanism, and it is worth more than the patch was. The locus
`55,896,39x` is a T homopolymer where **neither haplotype is reference**:
maternal reads carry -3/-4 (32 reads), paternal -16/-17 (32 reads).

**The alignment channel alone represents this correctly**, as two alleles on
opposite haplotypes:

| record | DP | ref/alt | AF | category | emitted |
|---|---:|---:|---:|---|---|
| `55,896,395 CTTT>C` (-3) | 58 | 30/28 | 0.483 | `NOISY_CAND_HET` | `1\|0`, PS 55,883,019 |
| `55,896,395 CTTTTTTTTTTTTTTTT>C` (-16) | 30 | 0/30 | 1.000 | `NOISY_CAND_HET` | `0\|1`, PS 55,883,019 |

**The hybrid keeps only one of them, and mislabels it:**

| record | DP | ref/alt | AF | category | emitted |
|---|---:|---:|---:|---|---|
| `-3` | 58 | 30/28 | 0.483 | `NOISY_CAND_HET` | **not emitted, PS 0** |
| `-16` | 37 | **2/35** | **0.946** | **`CLEAN_HOM`** | **`0\|1`**, PS 55,889,113 |

Two separate things go wrong, and they answer the question directly.

1. **The `-16` record's counts change under injection**, from 30 reads at 0/30 to
   37 at 2/35. Its AF moves 1.000 -> 0.946, which is still above `max_af`, so it
   is classified `CleanHom`. It was never het *by allele fraction* in either arm --
   AF 1.000 in the alignment channel is equally hom-looking. What made it a usable
   het in the alignment channel was its **category** (`NOISY_CAND_HET`), not its
   AF, and injection is what replaces that category with `CleanHom`.
2. **The emitter does not consult the category.** It writes the phased genotype
   from the haplotype assignment, so a `CLEAN_HOM` candidate that acquired a phase
   set is emitted `0|1` regardless. That is why "hom in the BAM" and "emitted
   `0|1`" are not contradictory: they are two different fields, decided in two
   different places, and here they disagree.

So the record is not a het that injection created. It is a **hom-classified record
carrying a het genotype** -- and it is the only record left at a locus the
alignment channel had described correctly with two records. Pruning it therefore
looks like losing a het, because the locus' other allele was already gone.

## Where the fix belongs

Not in the prune, and not in the classifier's AF window. The locus has two
alleles, both non-reference, and any single biallelic record must call one of them
"reference" and land outside a het AF window. The fix is the multiallelic
representation already identified: emit the locus once with both alleles, so the
maternal -3 and paternal -16 are the two ALTs rather than two competing records of
which injection keeps the one that looks homozygous.

Only after that does "do not inject hom" become safe to enforce, because then a
`CleanHom` verdict will mean the locus really is homozygous rather than meaning
one allele of a het locus was measured against the other.

## The `0 / 30` records: a het that no allele-fraction test could have called

`0 ref / 30 alt` at AF 1.000 is not a heterozygous call, and it was not made by
one. Two different mechanisms decide the category and the counts, and they
disagree.

**The category comes from comparing the two haplotype consensuses, not from an
allele fraction.** `update_cand_var_profile_from_cons_aln_str2`
(`collect_phase_noisy.cpp:586`) walks the variant lists of the two MSA cluster
consensuses in parallel: a variant present in one consensus and absent from the
other is stamped `NoisyCandHet` outright, because it is what distinguishes the
haplotypes. That is sound reasoning and strictly better than an AF test -- the
`-16` deletion at `55,896,39x` really is on one haplotype only.

**The counts are then accumulated only over reads matching that record's own
allele.** A read carrying the other haplotype's `-3` deletion matches neither the
`-16` allele nor the reference across that window, so
`update_cand_var_profile_from_cons_aln_str2` records **no observation** for it
rather than a reference observation, and it is excluded from `DP` entirely. So
`DP` becomes the alt count, `REF_COUNT` is 0, and `AF` is identically 1.000.

Chromosome-wide this is not a corner case:

| category | records | with zero ref reads | AF exactly 1.000 |
|---|---:|---:|---:|
| `CLEAN_HET_SNP` | 59,787 | **0** | -- |
| `CLEAN_HET_INDEL` | 2,543 | **0** | -- |
| `REP_HET_INDEL` | 1,986 | **0** | -- |
| **`NOISY_CAND_HET`** | 25,518 | **1,769 (6.9%)** | **all 1,769** |

Only the MSA-constructed class is affected, which is consistent with the
mechanism: the clean classifiers derive the category from AF, so an AF of 1.000
would have made them `CleanHom` and they can never be in this state.

**The category is usually right; the counts are what is wrong.** Sampling 40 of
the 1,769 and scoring each against read truth: 31 are scorable at >= 10 ALT reads
and **22 of those 31 (71%) are real heterozygotes** at purity >= 0.90 -- including
`24,121,714` (DP 6 against 74 covering reads, 13 ALT reads, purity 1.000), the
locus the allele-length work was already stuck on. Across the sample the **median
`DP` is 0.393 of true coverage**: these records see under 40% of the reads that
cover them.

### Why this is the same defect as the `CLEAN_HOM` promotion

It closes the loop on the previous section. Injection did not turn a homozygote
into a het. `classify_graph_only_candidates` **re-derives** the category from the
allele fraction, and the allele fraction is 0.946 because of this counting defect,
so a consensus-derived het verdict is overwritten with an AF-derived `CleanHom`
one. Whichever code path touches the record last wins, and the record then carries
a `CleanHom` category with a `0|1` genotype.

Every AF-gated consumer is exposed the same way: the graph-only classifier, the
noise filter, and any threshold on `allele_fraction` sees a homozygote at 1,769
sites that are mostly real hets.

### The fix, and it is smaller than the multiallelic one

A read that **covers** the record's window but carries a different event should be
counted -- as reference in a biallelic record, or as the other ALT in a
multiallelic one -- not dropped. Counting it as reference at `55,896,395` gives
roughly 30 alt / 68 covered, AF ~0.44, which agrees with the consensus-derived het
verdict instead of contradicting it. That alone makes `AF` usable for
MSA-constructed sites and stops the `CleanHom` overwrite, and it is separable from
(and a prerequisite for) multiallelic emission. It changes `DP`/`AF` on 1,769
chr20 records, so it needs the panel and the 0.559% chromosome-wide gate.
