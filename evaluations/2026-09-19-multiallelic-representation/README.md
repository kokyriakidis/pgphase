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
the graph arm still splits -- 0 multiallelic records on chr20 against 400
positions carrying two het rows. The graph arm is behind our own alignment arm,
not behind upstream.

The previous session note framed this as "the graph arm never emits a
multiallelic record", implying it alone was wrong. That was measured on our
arms only; upstream behaves the same way, and the intended behaviour is our
merge, not upstream's split.
