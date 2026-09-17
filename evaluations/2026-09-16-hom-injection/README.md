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
