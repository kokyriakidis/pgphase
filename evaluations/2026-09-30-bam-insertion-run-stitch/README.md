# A BAM insertion run attaches at chr20:1.18 Mb

The accepted baseline `/tmp/pgphase-weak-cut-conservative-full/` leaves
chr20:1,180,618 `T>A` and 1,194,189 `T>TATC` in separate VCF phase sets.
HiPhase, using the same variant keys, joins them with both ALT alleles on the
same haplotype. The independent BAM solve gives the boundary sites one source
PS, but marks a weak cut after the SNP. Recovery was preserving that cut
while attaching the downstream BAM rows to the large right graph block.

The recovery-input matrix (`/tmp/pgphase1180-matrix/`) shows that the
right insertion at candidate 237 and the next deletion and two SNPs are
linked by split-half supported allele edges. The next 155-bp insertion is near another graph insertion and could
represent the same event, so it is left in the right block. Six MAPQ-30
reads already assigned to the established left PS have a callable insertion
allele: three REF, three ALT. All six select the same orientation, with
one-sided binomial p=0.015625. The new stitch moves only the four short-run
BAM rows after checking that the left graph path is supported; it does not
move the right graph block or any of its read labels. A premature qname
deduplication in the trial initially suppressed valid later allele profiles;
deduplication now follows a callable allele.

| Full chr20 | Baseline | Local-run stitch |
|---|---:|---:|
| Variant keys | 62,269 | 62,269 |
| Truth-scored phased reads | 236,867 | 236,866 |
| Truth-correct reads | 229,131 | 229,131 |
| Discordant reads | 7,736 | 7,735 |
| Read phase sets | 697 | 697 |

Only four VCF sample fields change: the insertion at 1,194,189, the deletion
at 1,194,227, and SNPs at 1,194,230 and 1,194,233 move from PS 1,450,106
to PS 1,142,089 with unchanged genotypes. The graph insertion at 1,196,967
stays in PS 1,450,106. The owning-window test reports 73/129 correctly
separated reads, up from 66/129; HiPhase has 125/129 in that window.

The final output is `/tmp/pgphase-run-bridge-final-full/`. Reproduce its read score
with `python3 /tmp/score_pgphase_truth.py` while that temporary audit script
and the full outputs remain present. The window regression pins the
exact boundary keys, their relative GT, and the retained graph cut. Runtime
uses read alleles and labels; parental truth only scores the result.
