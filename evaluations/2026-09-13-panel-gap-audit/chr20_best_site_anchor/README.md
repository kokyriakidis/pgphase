# Best-covered direct site evidence

This trial evaluates verified repeats against individual clean heterozygous
anchors using read allele observations, including untagged reads. It chooses
the clean comparison with greatest minimum allele coverage, then total coverage;
equally covered comparisons must agree on eligibility. This prevents a favorable
small subset from overriding a fuller comparison with no haplotype separation.
The ordinary k-means, recovery allele graph, and flank stitching remain in use.
Truth and competitor data are evaluation-only.

Completed 18 regression cases: 12 target joins, 6 splits; all 5,922 original
truth-assessed reads preserved. Compared with direct_site_anchor, 54.89 Mb
changes from a wrong join (199 discordant reads) to a split (zero); 36.03 Mb
stays joined and improves from 5 discordant reads to zero. No case adds discordance
or switch/flip errors. Correct new joins at 4.85, 55.9, and 57.84 Mb remain.
The 55.9 Mb case has one block, 308 assessed reads, and zero discordance.

Unresolved wrong joins remain at 35.9, 46.7, and 60.1 Mb. The 114-case panel was
completed for the superseded direct_site_anchor trial; this final selection
revision has only the 18-case validation. No corrected whole-chromosome recovery
benchmark or production-safety claim is made. Build and all C++ unit tests pass;
the sparse-subset regression fails against the previous selector and passes here.
