# Phased BAM recovery rows hidden by graph repeat duplicates

The exact chr20:20,711,883 `CA>C` deletion occurs inside the open
20,707,556–20,731,688 graph gap. The targeted BAM solve phases it and the
recovery audit reports `APPENDED=1`, `META_BUILT=1`, and verified alignment.
The graph catalog has an unphased repeat row for the same normalized allele.
Before the fix, final output deduplication retained that 47-read graph row over
the 44-read phased BAM row solely by total coverage, then the VCF filter omitted
the repeat row. The deletion disappeared even though recovery had injected it.

Final graph output now prefers a valid phased heterozygote over an unphased
copy of the same normalized variant. Coverage still decides between copies
with the same phase status. The owning 20–21 Mb chunk emits the deletion once,
phased in the left block. A targeted window regression asserts this exact row.
The right SNP at 20,731,688 stays in another phase set: only four MAPQ-30
primary reads physically span the deletion and right SNP, and their paired
alleles do not establish a decisive parity. This change repairs site loss
without forcing an unsupported stitch.

On complete chr20, variant keys increase **62,154 → 62,269**. All 115 added
keys are present as phased exact alleles in the independent BAM-only output;
no key is removed and no sample field changes on the 62,154 shared keys.
Phased VCF blocks change **337 → 343**, while span N50 remains **616,859 bp**.
At this output-only stage the phased BAM was byte-identical to the accepted
baseline.

Inputs and outputs: accepted baseline
`/tmp/pgphase-physical-cert-final-full/`, this run
`/tmp/pgphase-dedup-phased-full/`, BAM-only comparison
`/tmp/pgphase-full-bam-gap-audit/`. The retained fix uses no parental truth in
runtime decisions.

## Source weak-cut transfer follow-up

The initial deduplication exposed a wrong relative phase at
chr20:33,227,050 `TA>T`: the local BAM sub-solve assigned the deletion and
chr20:33,211,072 `T>TA` to one source PS with opposite ALT haplotypes, while
the BAM-only whole-chromosome solve put both ALT alleles together. In the
isolated recovery matrix, **no read calls both sites**. The source-path audit
correctly marks a weak cut between them, but transfer had kept the deletion
in the left graph PS because it required *two imported near-side sites* before
detaching the far side. The near-side imported insertion already had a phased
graph anchor, so that extra requirement discarded valid boundary information.

Transfer now accepts an oriented graph anchor as near-side support only
when the cut follows an indel and no recovery read calls both sides. Cuts
following a clean SNP retain the later physical SNP-to-deletion bridge check;
this preserves the independently certified 3.574-Mb join. The 33.227-Mb
deletion remains phased with the supported downstream block, in a different
PS from the 33.211-Mb insertion. A new owning-chunk regression checks this
separation and the downstream allele orientation. Runtime decisions use
read/site evidence only; parental truth was used to audit the result.

Complete chr20 still has **62,269** variant keys, including all 115 gained
from deduplication. Against the dedup-only run, six existing VCF sample
fields change near weak cuts, and three reads change PS labels. Truth-scored
phased reads remain **236,867**, with **229,131 correct** and **7,736
discordant**. Phased VCF blocks change **343 → 344**. Of the
115 added sites, 100 have a nearby shared BAM-only anchor in the same two
blocks and all 100 have consistent relative phase; the remaining 15 have no
such shared anchor. The whole-chromosome output is
`/tmp/pgphase-weak-cut-conservative-full/`.

Validation: `make -j4`, `make unit-tests`, `make check`, and the complete
`make window-tests` panel pass. The window panel has **3,521 assertions in
44 cases**, including the added duplicate-site and weak-cut regressions and
the protected 3.574-Mb physical bridge.
