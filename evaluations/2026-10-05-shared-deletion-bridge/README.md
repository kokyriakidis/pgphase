# Close chr20:62,408,056–62,432,427 through a retained shared deletion

The interval is open in the task-start binary and closed correctly by HiPhase.
Both programs initially assign 157 of the same 177 truth-scorable original
primary overlaps correctly. pgphase splits those reads across two cores and
independent rescue labels; its largest connected core contains only 80 correct
reads. The fix joins the variant blocks and places all 157 correct reads in one
core: **88.70% of all original primary overlaps**, including abstentions.
Both tools retain 2 discordant reads and 18 abstentions. Competitor start/end,
CIGAR and sequence match all 177 original input alignments.

| Measurement | Before | After | HiPhase |
|---|---:|---:|---:|
| VCF spans interval | No | Yes | Yes |
| Correct original reads | 157 | 157 | 157 |
| Correct reads in one core | 80 | 157 | 157 |
| Discordant original reads | 2 | 2 | 2 |
| Unphased original reads | 18 | 18 | 18 |

The graph catalog deletion `GA>G` at 62,410,496 retains the verified independent
MSA deletion but keeps repeat classification. Its physical edit begins at
62,410,497. Requiring the graph row's own MSA flag or the original source key
therefore misses valid evidence. HiPhase retains this intermediate deletion
in the same phase set as both boundary SNPs.

The new bridge preserves the catalog representation and verifies its physical
gauge: 22 REF and 17 ALT Q30 primary molecules agree with their independent
clean graph SNP gauges, with no contrary votes. Both graph SNP paths are
complete. Two primary molecules call deletion REF and the right SNP REF,
producing same-gauge log odds −8.89121, above the existing 0.001 confidence
requirement. No quality or path threshold is relaxed for block joining.

The union remains deferred until rescue finishes. Twelve already rescued
marker reads then pass their own physical deletion checks, the complete
independent source, consistent shared SNPs and existing coverage guard. Two
have primary MAPQ22/26; their measured combined error remains below 5%, with
Q40 deletion flanks. Shared deletion rescues therefore use a MAPQ20 floor
plus the same measured error bound. Their graph observations agree with the
physical alleles even though the retained graph representation has no BAM
observation vector.

Two existing whole-chunk BAM fallbacks falsely report graph SNP REF at
62,432,427 although their primary CIGAR deletes that base. A separate Q40
SNP at 62,435,323 agrees in both graph and BAM profiles and independently
matches the final core's physical REF/ALT gauge. The fix connects those tagged
reads only after the certified shared-deletion union, with physical deletion
checks at every ignored zero-quality REF and no other contrary phased call.

The owning replay is 62,000,001–63,000,000. Disjoint parental flanks score
86/86 and 190/190 for pgphase, with the same orientation as HiPhase's 86/86
and 191/191. The owning audit preserves all 1,207 variant alleles, counts and
filters, all 3,722 correct assignments and all prior phased assignments.
The changed phase-set labels and global gauge swaps are audited explicitly.

The new panel row, required marker list, native owner, strict certification
manifest and owning mechanism regression protect actual closure, >=80%
correctness, HiPhase total/core parity and parental orientation. The new
regression passes 500 assertions and fails 10 assertions with the task-start
binary. Fast in-memory fixtures include shared-deletion provenance and masked
REF/witness failures; `make gap-dev-check` takes 0.04 seconds warm.

Task-start binary SHA256:
`6a8057fc5d3875b7fb3a608a579efe47a2f4d1ffce99aac84bdd0dcec6321f04`.
Final binary SHA256:
`f98aead0dfc6f9e72a2e2cb99de3463234523d5be7437859cecb095a2259e6fd`.

Measurements and reproducible audits are in `hiphase.json`, `orientation.json`
and `owner-audit.json`. Final suite and full-chromosome audit results are
recorded separately in this directory.

Final validation passes all **85 registered gap checks**, including 84 mechanisms
and all 113 panel windows, across 197 selected invocations and 84,327 assertions.
All standard build/unit/predicate/golden gates, three gap unit cases, four
cache-helper tests and two benchmark-helper tests pass. The tightened panel
totals check also passes 1,018 assertions. The panel has 100/113 measured spans
and 35 windows meeting both >=80% correctness and HiPhase total/core parity.
All 329 matching task-start native replays preserve variant alleles/counts/
filters and every correct or phased assignment; recorded changes are reviewed
PS unions rather than discarded tags.

The full chr20 replay takes 351.70 seconds and confirms the same 157/177
connected-core result. `full-closures.json` verifies this is the only new VCF
closure, with no previous extent lost. `full-audit.json` compares the saved
full chromosome from before the immediately preceding masked-SNP fix: all
64,188 variants and all previous correct/phased reads survive. Its one added
correct tag is the previously documented masked-SNP read, not an additional
claim for this task. `panel-audit.json` uses the actual immediately preceding
suite state and isolates the current changes.
