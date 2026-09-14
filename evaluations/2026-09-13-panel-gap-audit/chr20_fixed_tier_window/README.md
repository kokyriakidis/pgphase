# Preserve the gap MSA window across evidence tiers

The MSA search interval is fixed when the local recovery window is created.
SNP-only extensions can still move stitching boundaries, but cannot drop indels
before the later tier gets to use them. Includes the exact-consensus correction
for a verified SNP within a deletion's non-deleted allele.

Completed 114 cases: 75 joined targets, 35 split, four endpoint-unphased.
Compared with the prior full direct_site_anchor trial, no region increases read
discordance. The previously wrong 35.9 Mb join is corrected (293 assessed reads,
48→0 discordance). Two original-block orientation failures remain, at 46.7 and
60.1 Mb. These are local overlapping windows, not whole-chromosome metrics.
See chr20_35919404/README.md for the causal evidence.
