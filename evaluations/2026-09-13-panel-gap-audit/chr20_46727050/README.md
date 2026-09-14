# chr20 46.7 Mb: graph SNP vote over an explicit BAM deletion

The incorrect recovery join was supported by two left-to-right reads:

- m84031_231217_062403_s3/184948613/ccs: left SNP 46727050 was encoded as
  ALT although its primary MAPQ60 BAM alignment has a one-base deletion there.
  The local reference AAG aligns to A-G; there is no nucleotide observed at the
  SNP. The graph-only site filled the BAM unknown with a traversal-derived ALT.
- m84031_231217_034919_s2/78578596/ccs: left SNP ALT is observed, but its
  MSA insertion observation near 46747599 conflicts with the right SNPs.

Graph augmentation now abstains at a SNP when an existing BAM read explicitly
deletes or skips its reference position. Other missing graph-site observations
can still be filled. The unit regression verifies deletion and reference-skip
abstention, normal-site filling, and no inflated allele counts. It fails before
and passes after the change.

Initial 18-case graph_deletion_guard validation: the 46.7 Mb region stays split,
retains 415 truth-assessed reads, and improves 187→1 discordance and 4→1
switch/flip errors. This removes the amplified wrong-block orientation; it does
not yet reproduce competitor contiguity. No region adds discordance or switch/
flip errors relative to fixed_tier_window. The 35.9 and 55.9 Mb rescues remain.
The expanded 114-case validation is complete; see chr20_graph_deletion_guard/README.md
for the remaining 60.1 Mb block reversal and individual-read/contiguity tradeoffs.
