# Testing that the alignment channel's sites are injected correctly

Three separate test cases, because they fail for different reasons and one
"injection works" assertion would hide which:

| test case | tag | asks |
|---|---|---|
| the alignment channel's sites all reach the hybrid | `[completeness]` | is every site the alignment channel finds present in the hybrid's table? |
| shared sites keep their alleles | `[representation]` | is each locus described by the same alleles, with no second description added, and no multiallelic record genotyped homozygous? |
| read counts and alleles per site are correct | `[counts]` | are the per-site counts internally consistent and within the reads available? |
| REF matches the reference sequence | `[alleles]` | does each candidate's REF agree with the FASTA under that type's convention, and does every emitted record obey the VCF convention? |
| emitted records agree with the candidate table | `[consistency]` | are the emitted depths the candidate's depths, and is every candidate inside the region asked for? |
| a wider window finds the same interior sites | `[stability]` | does retrieval of an interior site depend on how much flank was requested? |
| an injected claim's reference bases are consumed whole | `[claims]` | does a record carrying a claim's ALT consume the claim's whole REF? |
| a site the alignment channel had arrives unchanged | `[fidelity]` | are clean-class sites carried byte-identically, and is nothing demoted? |

`src/test_bam_site_injection.cpp`, Catch2, `make window-tests`. Both channels are
run once per window and cached, so a test case costs assertions rather than
pipeline runs: 8 cases, 127 assertions, ~35 s over the six panel windows (the
stability case adds one widened run per window).

## Representation is more than keeping the alleles

Merging co-located records into one multiallelic record is the **default** in
the alignment channel -- `merge_msa_insertion_alleles` and
`merge_msa_colocated_deletions` are each called unconditionally -- and it is the
correct representation, not something to split back into biallelic records. The
27 multiallelic records the channel emits across the six windows are mostly
`1|2` or `2|1` with **zero reference depth**:

```
55,883,019  AATATAT > AAT,A    GT=2|1  AD=0,29,33
55,896,395  CTTTT…  > CTTT…,C  GT=1|2  AD=0,35,34
```

Split those into two biallelic records and each allele is measured against a
reference no read carries: the allele fraction runs to 1 and both classify
homozygous. That is the hom-at-het defect the merge exists to prevent, and
`55,883,019` is one of the loci it was diagnosed on.

Two assertions follow from that, added after the first version of these tests
checked only that the alleles survive into the candidate table.

### A multiallelic record with reads on both alternates is not homozygous

Checked in **both** channels. The alignment channel is where these records are
built, so a wrong genotype there is the origin; the hybrid is where it would be
carried. Looking only at the hybrid hides the defect completely whenever the
hybrid does not emit the record -- which is currently every one of them.

Three are caught and listed in the allowance:

| locus | REF>ALT | GT | AD |
|---|---|---|---|
| 55,795,217 | `CAAAA…>CA,C` | `1\|1` | 0, **40**, **28** |
| 55,815,775 | `CAAAA…>CAA,C` | `1\|1` | 0, **18**, **10** |
| 5,379,662 | `CTTTT…>CT,C` | `1\|1` | 0, **45**, **31** |

These are the loci the joint two-haplotype orientation branch does not reach.

### How many multiallelic loci the hybrid actually emits

**0 of 27.** Surviving into the candidate table is not the same as reaching the
output: the hybrid's solve excludes the noisy class (`skip_noisy_kmeans`, a
hybrid-only override), so a merged `NOISY_CAND_HET` record gets no phase set and
is never written. Recorded as a per-window floor rather than asserted equal,
because admitting that class chunk-wide is a measured bad trade -- but recorded,
so a change in either direction is visible.

## Threads: a floor of four, not a default

Every pipeline run these tests make uses at least four threads.
`PGPHASE_TEST_THREADS` can raise it; a lower value is clamped up rather than
honoured, so a test run cannot be made accidentally serial. Each invocation is
written to `cmd.txt` beside its logs, so the thread count actually used is
checkable rather than assumed.

Four rather than more, because that is where the scaling stops paying. Measured
on `chr20:30,000,000-35,000,000`, varying only `-t`:

| threads | real | user | sys | avg parallelism `(user+sys)/real` | speedup `real_1/real_n` |
|---:|---:|---:|---:|---:|---:|
| 1 | 161.6 s | 149.8 s | 12.3 s | 1.00 | 1.00x |
| 4 | 46.7 s | 155.5 s | 16.2 s | 3.68 | 3.46x |
| 12 | 43.6 s | 163.7 s | 18.4 s | 4.18 | 3.71x |
| 20 | 44.0 s | 164.1 s | 18.6 s | 4.15 | 3.67x |

An earlier version of this table omitted the `sys` column, which made the
parallelism figures impossible to reconstruct from the row they sit in and look
inflated by 8-12%: they are `(user+sys)/real`, not `user/real`. Both definitions
are given above. The useful number for how long a run takes is the speedup
column, which tops out at **3.7x**; the parallelism column includes system time
and reads ~4.2.

Total CPU work is flat and the VCF is byte-identical at every thread count, so
this is wall time only. The ceiling is structural: `collect_pipeline.cpp:618`
caps workers at the number of chunks in a batch and joins between batches, and
chunk costs are uneven.

At panel-window size the suite is short enough that the thread count barely
shows -- the full run is 12.4 s wall against 12.3 s user, so the floor is a
guarantee about how the tests run rather than a speed-up here.

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

## Fixes

### 1. Multi-base substitution claims are no longer injected (fixed)

`augment_chunk_with_graph_sites` now drops a claim whose REF and ALT are the
same length and longer than one base when **every** differing offset already
carries the matching SNP in the alignment channel. The existing redundant-site
drop could not see these: it keys on `TYPE`, and `vcf_to_variant_key` types an
equal-length substitution as `DEL`, which never matches a `SNP`.

Dropped only when every offset is already called. A partially covered claim
still names a base the alignment channel did not call, and discarding it would
lose that base, so those keep the previous behaviour.

Measured: `39,895,614 CGC>C` and `39,896,921 TTG>T` are gone from the emitted
VCF and the two SNP pairs they sat beside are unchanged, down to their depths
and genotypes. Panel 0 concordant->discordant.

### 2. MSA-derived candidates now carry strand tallies (fixed)

`update_variant_depth_fields` derives `ref_cov`/`alt_cov` from `alle_covs` and
leaves the four strand fields alone, so every candidate the MSA path built or
re-counted carried coverage with no strand. `derive_msa_candidate_strand_counts`
fills them in one sweep from the read profiles and `ReadRecord::reverse`, called
once after `collect_noisy_vars_step4` -- derived, not accumulated at each site
that touches a count, which is the contract the graph-only backfill already
follows.

A record is filled only when the derivation reproduces its **own** `ref_cov` and
`alt_cov`. The counts pass through several stages after the profiles are
written, and apportioning a strand split across counts it does not match would
be inventing the tally rather than measuring it.

Measured across the six panel windows: **211 of 211** strandless candidates
filled, none skipped, so the derivation agreed with every record's own counts.
The allowance ceilings are now 0.

### 3. The over-counted record is not a double count (diagnosed, not fixed)

A probe on both per-read accumulation sites shows `55,905,752` receives **57
increments from 57 distinct reads**, all through the with-phase-set consensus
path. No read is counted twice, so the "accumulate onto existing counts"
hypothesis is wrong -- a guard against re-counting a read that already has an
observation was written, measured inert, and reverted.

What remains is narrower and still a defect: 57 distinct reads are credited with
an allele at a site only **56** reads overlap, at any MAPQ and primary-only, by
three independent counts. One credited read cannot see the site. `full_cover` is
computed against the consensus alignment string using a running
`delta_ref_alt` offset rather than against the record's reference span, so the
window a read is tested against is not the span the record names. Establishing
that is a change to the MSA coordinate mapping, which is why it is recorded
rather than guessed at. The allowance row stays, with this mechanism.


## What the four retrieval cases found

**`[stability]` and `[consistency]` are clean.** A window widened by 20 kb on
each side finds exactly the same candidates in the shared interior, in all six
windows -- so retrieval of an interior site does not depend on how much flank
was asked for. Emitted depths equal the candidate's depths everywhere, and no
candidate falls outside the requested region.

**`[alleles]` found a convention, then a defect.** The candidate table does not
use one coordinate convention for all types:

| type | POS is | REF holds |
|---|---|---|
| SNP, DEL | the site | the reference string starting at POS |
| INS | one past the anchor | the ANCHOR base, at POS-1, with ALT holding only the inserted bases |

186 SNP and 21 DEL candidates match the reference at POS; 30 of 36 INS match at
POS-1, the other 6 ambiguous because the anchor sits in a homopolymer. The first
version of this test asserted the VCF convention for insertions and failed all
30 while saying nothing true -- the assertion was wrong, not the pipeline. The
emitted VCF is checked separately and obeys the VCF convention in every record.

Two candidates break even that: `12,680,257` and `55,905,390` hold the base at
their own POS rather than the anchor.

**`[claims]`: a claim whose REF runs past the anchor is emitted with only the
anchor consumed.** Those two candidates are the candidate-table symptom of it.

| locus | catalog claims | we emit | left unconsumed |
|---|---|---|---|
| 12,680,256 | `AT > AAATAAAATAAAATA` | `A > AAATAAAATAAAATA` | `T` |
| 55,905,389 | `CT > CCG` | `C > CCG` | `T` |

A claim whose REF extends past the anchor says those extra bases are
**replaced**. Emitting the claim's ALT while consuming only the anchor leaves
them in place, so the record asserts `AAATAAAATAAAATA-T` where the claim said
`AAATAAAATAAAATA`. The alleles look carried and the sequence is wrong.

The check is deliberately narrow: only a record carrying a claim's ALT
**verbatim** while consuming fewer reference bases is flagged. Of 298 claims in
these windows with REF past the anchor and an insertion-like ALT, 11 have a
record at the same position and only these 2 match that signature -- the rest
are different events the alignment channel called there, which is not this
defect. A first, looser version of this measurement flagged 10 of 11 by
comparing any record at the same position, which conflated the two.


## Are the alignment channel's sites injected "as they are"?

Measured over the 1,517 candidate keys the two channels share:

| the site's class in the alignment channel | shared | counts changed | promoted | demoted |
|---|---:|---:|---:|---:|
| `CLEAN_HET_SNP` | 461 | **0** | - | **0** |
| `CLEAN_HOM` | 727 | **0** | - | **0** |
| `CLEAN_HET_INDEL` | 13 | **0** | - | **0** |
| `NOISY_CAND_HET` | 173 | 16 | 37 | 0 |
| `NOISY_CAND_HOM` | 143 | 64 | 68 | 0 |

**A clean-class site is carried byte-identically** -- all nine count fields, in
1,201 of 1,201 -- and nothing is demoted anywhere. Every one of the 105 category
changes is a promotion to a clean class, which is what honouring a graph claim
on an existing record is supposed to do.

**Noisy-class sites are re-measured, and that is a repair.** 80 of 316 have
their counts recomputed, mostly upward. At the four largest changes the reads
settle which value is right:

| locus | alignment channel | hybrid | reads overlapping |
|---|---|---|---:|
| 12,732,381 | DP 11 (0/11) | DP 77 (0/77) | 78 |
| 12,735,985 | DP 8 (1/7) | DP 72 (0/72) | 75 |
| 5,325,616 | DP 7 (0/7) | DP 70 (0/70) | 72 |
| 24,123,429 | DP 7 (0/7) | DP 69 (0/69) | 71 |

The alignment channel shows the depth-starvation signature -- a count a tenth of
the coverage -- and the hybrid lands within a few reads of it. The hybrid has
more sites in the region, so the noisy-region MSA has more to work with and
assigns an allele to more reads. Asserting identity here would pin the starved
value as correct, so `[fidelity]` does not: it requires identity for the clean
classes, forbids demotion everywhere, and reports noisy drift. The bound that
does apply to those sites -- depth within the reads that overlap them -- is
asserted in `[counts]`.
