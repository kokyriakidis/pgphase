# Gap-local homopolymer bridges with unambiguous edge support

This successor to `../chr20_adaptive_hp/` changes only the extra homopolymer
edges: they must meet the existing minimum read count and have no opposing
allele votes. Ordinary clean/MSA links retain their existing majority rule.
Only unresolved gaps reach this fallback after all three ordinary tiers;
only verified homopolymer sites within the remaining gap are eligible.
Normal k-means and stitching are reused. Both original flanks must anchor the
same proposal block, and failed trials leave original assignments untouched.

This is a conservative evidence-selection trial, not a zero-read-error policy.
Correct relative block orientation is the evaluation target. Individual read
noise and global HP1/HP2 label swaps are distinct from reversing one joined
block. `block_orientations.tsv` records this distinction against read truth.

In the prior 55.9 Mb trial, both bridge edges had conflicting votes (4:1 and
18:3). All 39 right-block reads became discordant while all 269 left-block
reads stayed concordant. No new reads were tagged. Shared SNPs 55999561 and
56040612 are both 0|1 in assembly truth, but the trial emitted 1|0 and 0|1 in
one PS: an incorrect stitch, independently of read-level accuracy estimates.
The revised trial rejects this join. The 62.4 Mb gap has an alternative 2:0
right-side link and remains joined; 57.1 Mb also remains joined. The previously
correct 54.5 Mb fallback join is conservatively withheld, which is a real
sensitivity cost and not evidence that its prior orientation was incorrect.

Reproduce with a fresh label after changing source:

```sh
python3 evaluations/2026-09-13-panel-gap-audit/run_link_panel.py --chromosomes chr20 --label chr20_adaptive_unanimous --all-concordant --workers 8
python3 evaluations/2026-09-13-panel-gap-audit/summarize_link_panel.py --label chr20_adaptive_unanimous --baseline-label chr20_graph_support2
```

The binary is frozen in `/tmp/pgphase-chr20_adaptive_unanimous/pgphase` and its
hash is recorded here. The 114 local cases include overlapping windows and
are not independent chromosome-wide observations. Full-chromosome validation
of this fallback has not been run. The underlying graph-support2 arm has a
separate full-chromosome hamming regression documented in `../chr20.metrics.tsv`.


Final 114-case local result for chr20_adaptive_unanimous: all 85 previous
endpoint joins retained, two additional joins (57.1 Mb and 62.4 Mb), 22 targets
still split, five with an unphased endpoint. No case adds discordant reads or
switch/flip errors relative to graph_support2. All 37,479 original tagged-read
occurrences retain uniform original-block transforms; the 851 previously
added tagged-read occurrences are unchanged. Windows overlap, so these sums
are not independent chromosome-wide counts. The 55.9 Mb bad stitch is rejected;
the correctly oriented 54.5 Mb trial join is also withheld, documenting the
conservative rule's sensitivity cost. Existing ordinary-tier orientation
failures, including 19 Mb and the full-chromosome graph regression, remain.
Build/unit tests and projection tests pass. Full-chromosome validation of this
fallback is not yet complete; make check is blocked by the stale gate CLI
invocation documented above.
