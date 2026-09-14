# Direct clean-site anchors for gap-local repeat recovery

A verified repeat can qualify through direct read allele co-occurrence with a
clean phased heterozygous site even when those reads lack HP tags. Each allele
of the clean site must favor the same association by the existing net margin.
Candidates are assessed separately; normal graph phasing and flank stitching
remain responsible for the final join. Truth and competitor calls are used only
for evaluation, not site selection or orientation.

Initial 11 regression cases: six target joins, five splits; no original block
majority reversal and no added discordance or switch/flip errors relative to
chr20_local_alleles. Newly rescued 55.9 Mb: 308 truth-assessed reads, zero errors.
Full 114-case expansion pending; this is not a whole-chromosome benchmark.

## Completed broad validation: superseded experimental trial

All 114 cases completed: 77 target joins, 33 splits, 4 endpoint-unphased.
Compared with chr20_unassigned_local, this trial adds correct joins at 4.85,
55.9, and 57.84 Mb, but adds a wrong join at 54.89 Mb (199 discordant original
block reads) and reverses a five-read block near 36.03 Mb. Existing wrong joins
at 35.9, 46.7, and 60.1 Mb also remain. It is not a validated production solution.

At 54.89 Mb, a 0/2 versus 3/0 distant-anchor subset qualified a repeat despite
the much fuller 18/8 versus 16/8 comparison at the nearby clean SNP. The latter
shows essentially identical allele ratios on both haplotypes. A subsequent
best_site_anchor trial uses the best-covered clean comparison (maximize the
less-covered anchor allele, then total coverage) rather than any favorable
comparison, with equally covered ties requiring consistent eligibility.
