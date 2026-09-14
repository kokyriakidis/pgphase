# Why the chr12 46.8 Mb recovery failed

Investigated the frozen HiPhase-bridged gap chr12:46725702–46838918,
using the same solve window chr12:46675702–46888918 and binary as the
six-region validation. The first residual hybrid break is
46769497–46810372. This is an evidence-admission and minimum-link-support
failure, not a failure to invoke the recovery tiers or the final stitcher.

## Concrete evidence

HiPhase phases these heterozygous variants together in PS 46709288:

| VCF position | Alleles | HiPhase GT |
|---|---|---|
| 46769497 | G / C | 0\|1 |
| 46774475 | G / GATATATATAT / GATATATATATATATAT | 1\|2 |
| 46774503 | C / T | 1\|0 |
| 46788671 | T / TAA / TAAA | 2\|1 |
| 46810372 | A / G | 0\|1 |

Recovery retrieves the repeat alleles, SNP, and later insertion. Its internal
insertion coordinates are one greater than the VCF anchor positions above.
In the final default gap-0 input matrix, only 7 reads have usable observations
at the first repeat and nearby SNP. The SNP has 4 reference / 3 alternate
observations. With only `--private-msa-margin 1` changed, this rises to 83 reads
(45 reference / 38 alternate at the SNP). The default threshold is 24.
`src/align.cpp` admits extra reads against the two MSA consensuses only when
the score difference meets this threshold. Thus retrieving a site does not
mean the spanning reads have usable allele observations there.

There are 18 primary, nonsupplementary, MAPQ60 BAM reads spanning the first SNP
and the insertion at 46788671, but zero shared usable observations between
those sites in the default final matrix. This admission loss is real, although
other recalled sites can provide alternate local connections.

The decisive later link is sparse in the BAM itself. Exactly one primary,
nonsupplementary read spans 46788671–46810372 (21,701 bp):
`m84031_231217_062403_s3/22155749/ccs`, alignment [46787973,46811749), MAPQ60.
It has observations at the corresponding recovered insertion and next clean
SNP. HiPhase tags that same read HP1, PS46709288. Its block includes both variants.
This supports a single-read bridge explanation; we did not instrument HiPhase's
internal solver to establish its exact chosen edges.

The pgphase link rule in `iter_update_var_hap_cons_phase_set` requires at least
`min_block_link_reads` agreeing or conflicting observations (default 2).
The test already uses `--link-by-alleles --block-link-window 8`, so the older
HP-tag-only circularity is not the explanation for this run.

## Controlled interventions

| Change from default recovery | Regional blocks | Evaluated reads | Discordant / switchflip |
|---|---:|---:|---:|
| None | 4 | 288 | 0 / 0 |
| MSA margin 24 → 1 only | 3 | 336 | 0 / 0 |
| Minimum block-link reads 2 → 1 only | 2 | 288 | 0 / 0 |
| Both | 2 | 336 | 0 / 0 |

Margin 1 alone extends the left block and closes a later gap, but leaves the
first target gap open. Minimum link support 1 closes the first gap at tier 3.
The original sites at 46769497, 46810372, 46818641 and 46838918 all receive
PS46709288 in that run, spanning the audited target. A separate gap farther
right in the solve-window flank remains, hence two regional blocks.

This isolates the two-read requirement as a binding restriction in this case.
It does not establish that lowering it globally is accurate: a singleton edge
has no independent replication. No production defaults or algorithm were changed.
A narrowly controlled singleton-bridge policy would need broader truth validation.
The earlier MSA observation loss merits separate evaluation of local allele
calling, rather than treating whole-window consensus separation as the only
way a read can supply an allele observation.

## Reproduction and limits

Each arm contains its pgphase command, tier report, and read-truth summary.
Heavy outputs and matrices are under `/tmp/pgphase-gap-46m`. Evaluate using
`scripts/evaluate_phase_accuracy.py` with the same truth subset from the
six-region validation (`/tmp/pgphase-auto-gap-validation/chr12_46725702/truth.bam`).
Graph-only frozen candidates in the first residual interval include a clean
SNP at 46774503 (DP17, AF0.705882, unphased) and a repeat indel there (excluded).
HiPhase VCF/BAM are the frozen chr12 comparison outputs. These diagnostic
threshold changes are local causal probes, not a new accuracy benchmark.
