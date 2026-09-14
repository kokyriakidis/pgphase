# Clean-block evidence for remaining gaps

Completed all 114 local chr20 cases: 86 joined targets, 24 split, four
endpoint-unphased. All 79 baseline joins remain. Seven additional correct
joins occur at 11.4, 18.2, 20.9, 24.4, 38.3, 54.9, and 57.7 Mb. No increased
read discordance or switch/flip counts versus multi_final; no original-block
majority reversals. Existing individual recovery errors at 7.2 and 61.7 Mb
remain. Original read assignments within each baseline block are preserved.

After ordinary supported site edges establish components, evaluate each
primary MAPQ30 read against their resolved clean-SNP orientations. A flank
needs at least two consistent clean SNPs, or an exact Q30 BAM base matching
an observed clean SNP. Missing qualities and BAM deletions/skips cannot
supply the single-site alternative. Any consistent opposing bridge read
vetoes the proposed direction, including a read with fewer supporting sites.
Process component bridges by supporting read count, composing parity with
existing stronger edges. Normal read phasing and gap stitching still finish
the proposal; no competitor/truth calls enter the acceptance rule.

The unit fixtures check multi-site support, high/low base and mapping quality,
opposing sparse reads, orientation composition and convergence. This isolated
change excludes the indel and phase-scoped-read experiments. Results are from
overlapping local windows, not whole-chromosome NGC50 or hamming validation.
