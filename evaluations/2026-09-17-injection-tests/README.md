# Testing that the alignment channel's sites are injected correctly

Three separate test cases, because they fail for different reasons and one
"injection works" assertion would hide which:

| test case | tag | asks |
|---|---|---|
| the alignment channel's sites all reach the hybrid | `[completeness]` | is every site the alignment channel finds present in the hybrid's table? |
| shared sites keep their alleles | `[representation]` | is each locus described by the same alleles, with no second description added? |
| read counts and alleles per site are correct | `[counts]` | are the per-site counts internally consistent and within the reads available? |

`src/test_bam_site_injection.cpp`, Catch2, `make window-tests`. Both channels are
run once per window and cached, so a test case costs assertions rather than
pipeline runs: 3 cases, 78 assertions, ~12 s over the six panel windows.

## What the counts test compares against, and what it does not

It does **not** compare the hybrid's counts against the alignment channel's.
That looks like the obvious test and it is wrong: the two channels legitimately
differ by a read or two, and neither is the reference. Measured against the
reads at offset 0 (the convention the candidate tables use):

| site | alignment channel | hybrid | reads in the BAM |
|---|---|---|---|
| `12,754,263 G>A` | 32 / 18 | **32 / 16** | **32 / 16** |
| `12,797,914 C>T` | **35 / 24** | 35 / 22 | **35 / 24** |

At one site the hybrid matches the reads and the alignment channel is inflated;
at the other it is the other way round. An assertion that the hybrid's counts may
only grow -- which is what this test first contained -- fails on both and
establishes nothing. Cross-channel drift is now reported, not asserted.

What is asserted is what the reads can settle:

- **internal consistency** -- `DP == REF_COUNT + ALT_COUNT`, `AF ==
  ALT_COUNT / (REF_COUNT + ALT_COUNT)`, and each strand tally summing to the
  count it accompanies;
- **depth within coverage** -- `DP + LOW_QUAL_COUNT` may not exceed the number
  of reads **overlapping** the site's reference span. Overlap, not coverage at
  the anchor base: a read touching part of a 24 bp deletion can be observed
  there without covering the anchor, which is why an anchor-based bound was
  exceeded by one read at `55,905,752` by a record that is not over-counted.

## Three defects the tests found on their first run

Recorded in `src/test_bam_site_injection_allow.tsv` so the tests gate on "no new
ones" rather than being red. Removing a row is how a fix gets locked in.

### 1. The catalog's multi-base substitutions are injected as spurious deletions

At `chr20:39,895,615` the alignment channel calls **two adjacent SNPs**, which is
the correct description:

```
39895615 G>A  DP=61 31/30 CLEAN_HET_SNP
39895616 C>A  DP=59 31/28 CLEAN_HET_SNP
```

The catalog claims the same two bases as one record, `REF=GC ALT=AA`. Injected,
it is typed `DEL`, left-shifted a base, and **emitted as a heterozygous 2 bp
deletion beside the two SNPs**:

```
hybrid VCF:  39895614 CGC>C  0|1:60:31,29   <- spurious
             39895615 G>A    0|1:61:31,30
             39895616 C>A    0|1:59:31,28
```

Same reads, described twice, and the extra record is not merely redundant: it is
a deletion call at a locus that carries no deletion. `39,896,922` is the same
case (`TG>CT` against SNPs `T>C` and `G>T`). This is the class the
redundant-site drop already handles for indels the alignment channel called at
the same position, but that check keys on TYPE, and an equal-length substitution
claim typed `DEL` does not match a `SNP`.

### 2. MSA-derived candidates carry no strand tally

29-45 candidates per window have `REF_COUNT` or `ALT_COUNT` above zero with
their strand fields at `0+0` -- for example `48,147,225 SNP C>G` at 44 / 7 with
no strand at all. The same defect was found and fixed on the injection path,
where four accumulation sites were changed to mirror each count increment; the
MSA path was not touched. It matters for the same reason: the ONT strand-bias
screen computes an expected value from `forward_alt + reverse_alt` and declines
to test when it is not positive, so these candidates are silently exempt from a
screen every initially-discovered candidate faces.

### 3. One record's depth exceeds the reads that overlap it

`55,905,752`, a 53 bp deletion, reports `DP = 57` where **56** reads overlap its
span. Three independent counts agree on 56: this test's htslib overlap count, a
Python check over the same BAM, and `samtools view -c` over the span at any
MAPQ and primary-only. The mechanism was not established; one record in six
windows, and the bound is tight enough that it would have gone unnoticed
without this check.
