# Does longcallD phase a multiallelic site, or split it?

Read from the source and measured on upstream's own chr20 output, which is in
the eval data (`##source=longcallD version=0.0.11-23e369d`).

## It splits

Whole chr20, phased records only:

| tool | phased records | ALT carrying a comma | `GT 1\|2` or `2\|1` | positions with 2+ rows |
|---|---:|---:|---:|---:|
| **longcallD 0.0.11** | 118,270 | **0** | **0** | **2,401** |
| pgphase alignment arm (2026-09-12 snapshot) | 118,270 | 0 | 0 | 2,401 |
| pgphase graph arm (snapshot) | 50,954 | 0 | 0 | 72 |
| pgphase hybrid (snapshot) | 93,311 | 0 | 0 | 233 |
| **hiphase** | 77,123 | **2,063** | **1,718** | 0 |

Upstream never emits a multiallelic record. Where both haplotypes are
non-reference it writes two biallelic rows at the same position -- 2,401 of them
on chr20, e.g.

```
34435952  REF=AAT    ALT=A     GT=0|1
34435952  REF=AATAT  ALT=A     GT=1|0
29696172  REF=G      ALT=GAA   GT=1|0
29696172  REF=G      ALT=GAAA  GT=0|1
```

hiphase is the tool that emits one record with two ALTs.

The code agrees with the measurement. `collect_var.c:1523-1559` allocates room
for two ALTs on a het and loops `hap = 1..2`, appending one per non-reference
haplotype with `GT[hap-1] = ++n_alt_allele`, and `vcf_utils.c:156-159`
comma-joins them, so the FORM exists. But the alt sequence comes from a single
candidate (`alt_len = cand_vars[cand_i].alt_len`) -- one ALT per candidate --
so a locus whose haplotypes carry different alleles is two candidates and
therefore two records. The comma path is reachable only when both ALTs are the
same string, which is why the count is zero in practice.

## Which makes our alignment arm the divergence, not the graph arm

The 2026-09-12 snapshot of our alignment arm matches upstream exactly: same
118,270 phased records, same 2,401 split positions, zero multiallelic. The
current binary does not:

```
55267796  REF=C      ALT=CTT,CTTTTT  GT=2|1
55272344  REF=AACAC  ALT=AAC,A       GT=2|1
55273343  REF=C      ALT=CT,CTT      GT=2|1
```

That is `merge_msa_insertion_alleles` / `merge_msa_colocated_deletions`, added
deliberately in this project: split into two biallelic rows, each allele is
measured against a reference no read carries, allele fraction runs to 1 and
both halves classify homozygous.

So the ordering is upstream splits, our alignment arm was changed to merge, and
the graph arm still splits. The graph arm is behind our own alignment arm, not
behind upstream.

## Current state, all figures from runs of the current binary

The table above mixes the Sep 12 eval-data snapshots with upstream; this one is
the binary as it stands, whole chr20, phased records only. The graph-arm
snapshot figure (72 two-row positions) and the current figure (400) are
DIFFERENT RUNS eleven days and many commits apart, and an earlier revision of
this document put them in one paragraph without saying so.

| run | phased | comma-ALT | `GT 1\|2` | positions with 2+ rows |
|---|---:|---:|---:|---:|
| longcallD 0.0.11 (upstream) | 118,270 | 0 | 0 | 2,401 (4,824 records) |
| alignment arm, Sep 12 snapshot | 118,270 | 0 | 0 | 2,401 |
| **alignment arm, current** | 116,179 | **1,664** | **1,657** | **806** (1,612 records) |
| graph arm, Sep 12 snapshot | 50,954 | 0 | 0 | 72 (145 records) |
| **graph arm, current default** | 62,352 | **0** | **0** | **400** (800 records) |
| **graph arm, current routed** | 69,906 | **0** | **0** | **653** (1,306 records) |
| hiphase | 77,123 | 2,063 | 1,718 | 0 |

Two things this makes precise.

**The divergence from upstream is 1,664 records.** That is the size of the
merge in the arm that has it: upstream writes 2,401 two-row positions on chr20,
we write 1,664 merged records and still 806 two-row positions.

**The merge is incomplete even where it exists.** 806 positions in our own
alignment arm still carry two biallelic rows. Either those are loci the merge
should cover and does not, or they are a class that must not merge (a SNP and
an insertion at one base are not two alleles of one event). That is a separate
question from propagating the merge to the graph arm, and it is answerable from
these 806 positions.

The previous session note framed this as "the graph arm never emits a
multiallelic record", implying it alone was wrong. That was measured on our
arms only; upstream behaves the same way, and the intended behaviour is our
merge, not upstream's split.

## Bringing the alignment arm to parity: what it takes, measured

The alignment arm is a port, so the merge is a divergence to reconcile. It is
now behind `merge_colocated_msa_alleles` (`--merge-colocated-msa-alleles`,
default ON -- see below), which makes the comparison runnable.

**Turning the merge off moves most of the way to parity.** Whole chr20 against
upstream's own output, matching on `(POS, REF, ALT)`:

| alignment arm | records | identical to upstream | upstream-only | ours-only | merged positions |
|---|---:|---:|---:|---:|---:|
| merge ON | 116,179 | 114,118 | 4,152 | 2,061 | 1,532 |
| **merge OFF** | **117,696** | **117,167** | **1,103** | **529** | **0** |
| upstream | 118,270 | -- | -- | -- | 0 |

Comma-ALT records go 1,664 -> 0 and two-row positions 806 -> 2,335 against
upstream's 2,401. Identical records reach 99.1%. It also *improves* the graph
arm: read hamming 1.153% -> 1.095%, misplaced reads 2,526 -> 2,401.

**But it is not parity, because our split is incoherent where upstream's is
not.** Counting positions at which ONE haplotype claims two different ALT
alleles -- the rule the window suite asserts, ALT-side only, REF-side claims
ignored:

| chr20 | positions putting two ALTs on one haplotype |
|---|---:|
| upstream longcallD | **0** |
| ours, merge OFF | **412** |
| ours, merge ON | 8 |
| graph arm, merge OFF | 0 |

Upstream writes two rows per locus that are complementary by construction: one
haplotype carries the ALT, the other the reference. Ours, with the merge off,
writes 412 positions where a haplotype is assigned two different alleles at
once, and the gap-window suite fails on exactly those.

So the parity task is NOT "stop merging". It is "make the split complementary
the way upstream's is", and the 412 positions are the work. Until then the
merged form is the one to ship, which is why the option defaults ON.

The residual 1,103 upstream records we lack with the merge off decompose as 916
at positions we do not emit at all, 148 indels where the position is present
but the allele differs, and 24 SNPs likewise -- a separate parity bucket from
the representation question.

## Two coherence fixes, and the real blocker to parity

### The alignment writer never ran the conflict resolution

`drop_conflicting_haplotype_alleles` was wired into the graph writer only
(`graph_collect.cpp`), so the alignment arm emitted positions where one
haplotype carries two different ALT alleles: 8 with the merge on, 412 with it
off. Upstream scores 0, because it writes one ALT per candidate and derives the
genotype from that candidate's own haplotype alleles, so its rows are
complementary by construction. Now wired into both writers.

### Dropping is the wrong resolution when both alleles have support

At `882,277` the two records are `A>ATC` (11 alternate reads) and `A>ATCTC`
(15), both called `1|0` against 2 reference reads. Both alleles are real, so
the locus is heterozygous with two non-reference alleles and the answer is one
allele per haplotype -- what upstream emits. Dropping the contained form costs
798 records on chr20 and takes two-row positions to 1,542 against upstream's
2,401.

`make_colocated_alleles_complementary` moves the less supported allele to the
other haplotype instead, and only when both carry at least `min_alt_depth`
alternate reads; anything weaker still falls to the drop.

Measured with the merge ON, whole chr20:

| | alignment arm | graph arm |
|---|---|---|
| positions putting two ALTs on one haplotype | **8 -> 0** | 0 -> 0 |
| phased records | 116,179 -> 116,170 | 62,352 -> 62,458 |
| read hamming | -- | **1.153% -> 1.153%** |
| tagged / misplaced | -- | 219,055 / 2,526 unchanged |

Read placement does not move; the fixes change which records are emitted and
on which haplotype, not how reads are assigned. Window suite 125/125.

### What still blocks parity: the split's genotyping, not its coherence

With the merge off the gap-window suite fails on `in_gap_hets`, not on
conflicts:

```
window 3,852,321  in_gap_hets  4 >= 17   FAILED
window 4,766,928  in_gap_hets  1 >=  6   FAILED
window 5,309,406  in_gap_hets  1 >=  3   FAILED
```

Split into two biallelic rows, each allele is measured against a reference no
read carries -- allele fraction runs to 1 and both halves classify homozygous,
so they never reach the output as heterozygotes at all. That is the defect the
merge was built to avoid, and it is upstream of every representation fix: the
coherence passes above cannot help a record that was never called het.

So parity requires genotyping a split allele against the LOCUS -- the other
allele and the reference together -- rather than against the reference alone,
which is what upstream's cluster-consensus comparison does implicitly. Until
that is done the merge stays on by default, and the option exists so the
comparison can be re-run in one command.
