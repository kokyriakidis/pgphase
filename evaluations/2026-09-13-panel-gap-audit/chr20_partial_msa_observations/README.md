# Remaining-gap diagnostic trial

Rejected cumulative indel/read-scoring trial. Full 114-case panel: 78 joins, 32 splits, four endpoint-unphased. It introduces a wrong original-block orientation near 34.1 Mb (2 to 63 discordant reads), loses correct 12.26 and 57.78 Mb joins, and retains the unresolved endpoint cases. Do not use this trial as a validated recovery implementation. Its failures motivated isolating the clean bridge and preserving the original HP genotype-update behavior.

The comparison baseline is chr20_multi_final. These are overlapping local
windows; they do not establish whole-chromosome accuracy or NGC50.
