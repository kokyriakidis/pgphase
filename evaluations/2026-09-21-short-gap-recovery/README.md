# Short-gap recovery validation (2026-09-21)

`final.tsv` records the targeted chr20 validation after the explicit-seam and
consecutive-seam gauge fixes. The panel is the 14 noncentromeric gaps in
`../2026-09-16-test-panel/short_gap_targets.tsv` that HiPhase spans at at least
98% local crossing-read truth purity.

Each row was produced by `collect-graph-variation` with BAM recovery enabled.
The outer region was 50 kb on each side, except the 32.17 Mb target, which used
250 kb because a 50 kb isolated region contains no established left graph
anchor. Candidate positions were normalized by removing the common REF/ALT prefix,
matching `vcf_to_variant_key()`. Phase-block extent includes the interval from
the raw VCF anchor through that canonical position, so an anchored indel covers
its breakpoint under either representation.

Truth scoring joined HP/PS from the emitted unaligned phased BAM to the primary
input alignment by qname, retained reads physically covering both normalized
gap boundaries, and compared HP with
`test_data/derived/chr20_truth_hap.tsv`. Haplotype numbers were allowed one
orientation flip per local phase set.

All 14 gaps are inside one pgphase phase set after recovery. The later
supported-chain correction changed the 55.88 Mb row from 44/56 (78.57%) to
64/64 (100%). The frozen HiPhase result there is 66/66. The union already held
the HiPhase bridge sites; two physical bridge reads were missing exact deletion
observations in the sparse matrix, and stitching copied stale HP labels from
the pre-bridge BAM sub-solve. Exact-CIGAR observation backfill plus final read
scoring on the locus that earned the strongest chain edge fixes both defects.
The remaining pgphase-unphased crossing reads have no decisive exact allele on
the two separate injected deletion rows.

The fast unit suite does not rerun these BAM/GAF regions. It replays the common
explicit-seam identity invariant at all 14 coordinates, plus the consecutive
seam alias case, in memory. The real-window results in this directory are the
integration evidence.

## Expanded current-deficit panel

`remaining_hiphase_correct_gaps.tsv` expands the fast replay beyond the original
14 short-gap regressions. It was derived from a fresh full chr20 graph+BAM run
of the current binary (`pgphase` SHA-256
`7070008f8aa9626e1bd9dd3e7077056fd3e290b8aacbc12c3e7058f1a179a280`).
The emitted VCF SHA-256 was
`7dc3624b2e4575cdfddd7a3f7fc3b4c4e66f110b16b891368d64a7bdf6c7f3a4`;
the locked HiPhase VCF SHA-256 was
`a5ac0e009867d3b6786b3e50f41c4f0e16d221921a7ea9ccaf5bff75183dac68`.

The current run has 200 gaps between emitted phase blocks. HiPhase spans 74 of
them. After excluding the centromeric interval and requiring at least 98% local
truth purity among reads physically crossing both flanks, 34 gaps remain. They
cover 380,363 bp and 780 truth-scored crossing reads. Fifteen have at least 20
scored reads; the other 19 remain in the target set with `support_class=low` so
the test does not silently discard a valid but sparsely observed chain.

The in-memory unit replay now covers 48 coordinate cases: the original 14
historical regressions plus all 34 current misses. The 4.78 Mb biological gap
appears in both lists under their respective coordinate-normalization conventions,
so these are 47 distinct genomic regions. The current deficit itself is 34
gaps. The replay exercises explicit seam identity and final left-to-right
stitching without rerunning BAM/GAF discovery. The TSV is the integration
inventory and retains the measured support for each current miss.


## Per-edge aggregate and exact path fallback

`dp_path_results.tsv` compares the exact pre-DP full-chr20 result with the
retained per-edge recovery. Recovery evaluates every adjacent phase-set pair,
so an unsupported edge no longer hides a supported pair later in the same
window. Distributed boundary evidence is combined as one vote per molecule;
when that abstains, the exact ordered two-state DP searches injected BAM sites
between that pair.

Production chr20 closes all 15 high-support current deficits and one low-class
gap with 42 net aggregate votes. All newly closed high-support gaps have at
least 84.53% local truth concordance. The initially reported miss at 882277 was
a span-evaluator bug: the recovered VCF already assigns the insertion rows at
raw position 882277 and the right-flank site at 890261 to PS 882277. Block
extent now includes both an indel's VCF anchor and canonical event coordinate.
Weak thresholds were tested and rejected after producing a 68% local join.

The full comparison used identical 500 kb chunks and 16 workers. The final
candidate TSV, VCF, and BAM are byte-identical before and after replacing the
aggregate's all-read scan with the existing candidate-interval index. Hashes
and chromosome-level counts are recorded in `dp_path_results.tsv`.
