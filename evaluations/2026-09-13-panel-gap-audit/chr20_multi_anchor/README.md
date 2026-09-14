# Per-site validation of MSA insertion pairs

Completed 114 overlapping local chr20 cases against chr20_graph_deletion_guard:
79 joined targets, 31 split, four endpoint-unphased. Five newly joined cases
cover four locations: 9.0, 11.6, 55.5 (two overlapping windows), and 58.7 Mb.
No target joins are lost, no discordance or switch/flip increases versus the
committed baseline, and no original-block majority reversals are detected.
60.1 Mb retains its join and improves 39 -> 0 discordance, 3 -> 0 switch/flips
among 206 assessed reads. 35.9 Mb retains its correct 293-read one-block result.

The fix preserves two alternate insertion alleles and their per-read evidence,
optimizes their genotype orientation jointly, and uses the best-covered clean
anchor to reject MSA length artifacts lacking allele separation. No coordinates
or truth/competitor labels enter the acceptance rule.

Recovery still adds some individual discordance relative to clean-only phasing
at 7.2 and 61.7 Mb, already present in the committed baseline. The 58.7 Mb newly
joined case retains one preexisting discordant read. These are local read-truth
comparisons, not full-chromosome hamming or NGC50 validation. The frozen binary
predates the final output AF correction and restriction of pair creation to the
MSA gap pass; see chr20_multi_final for final-binary confirmation.
