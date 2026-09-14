# Graph/BAM SNP conflict validation

Graph profile extension abstains when an existing BAM read explicitly deletes
or skips a SNP position. Includes the fixed tier window and overlapping-SNP
allele corrections. Build and all unit tests pass; the deletion/skip regression
fails against the previous implementation.

Completed 114 local cases: 74 joined targets, 36 split, four endpoint-unphased.
Only the 60.1 Mb case retains a detected reversal of an original block.
At 46.7 Mb, the wrong join is prevented: 415 assessed reads, 187→1 discordance,
4→1 switch/flip errors, target split. Correct 35.9 Mb (293 reads, zero errors)
and 55.9 Mb (308 reads, zero errors) joins remain.

Tradeoffs versus fixed_tier_window: the previously joined 61.7 Mb target now
stays split and adds three discordant newly tagged reads; two reads near 1.9 Mb
become discordant (the same reads appear in two overlapping test windows).
Switch/flip counts increase in two overlapping 26.6 Mb cases despite fewer
discordant reads. No additional original-block majority reversal is detected.
These limitations are retained in results.tsv, baseline_comparison.tsv, and
block_orientations.tsv. This is not a whole-chromosome accuracy claim and the
60.1 Mb multi-allelic insertion problem remains unresolved.
