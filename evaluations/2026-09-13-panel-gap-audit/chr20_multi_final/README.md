# Final MSA insertion-pair recovery validation

The final binary includes multiallelic insertion identity, bounded local allele
observations, joint genotype orientation, and best-clean-anchor eligibility.
It also restricts pair creation to the MSA gap pass and emits per-ALT INFO AF
with Number=A. Build and all unit tests pass. Native VCF at 60093416 correctly
reports C -> CTTTTTT,CTTTTTTT, GT 2|1, AD 0,22,31, and the stitched phase set.

All 228 clean/recovery BAMs across the 114 local cases were compared with the
multi_anchor frozen binary: read names, flags, starts, HP and PS are identical.
The independent 114-case read-truth evaluation completed: 79 joins, 31 splits,
four endpoint-unphased. All 74 baseline joins remain, with five additional
joined cases across four locations (9.0, 11.6, 55.5 and 58.7 Mb). No discordance
or switch/flip increases versus baseline and no original-block majority
reversals are detected. 60.1 Mb improves 39 -> 0 discordance among 206 reads.
Existing individual recovery errors at 7.2 and 61.7 Mb remain. Results are in
results.tsv and summary.json. Comparison baseline is chr20_graph_deletion_guard.
These are overlapping local windows, not whole-chromosome accuracy/NGC50 results.
