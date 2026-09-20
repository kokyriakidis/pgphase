# Scoping `joint_het_orientation` to the alignment path: fixes two known sites, breaks ten more

## The hypothesis

`allele_depths_call_het` (`collect_phase.cpp:744`) already contains the fix
for a biallelic candidate whose own allele depths clearly call it
heterozygous but whose independent per-haplotype majority vote collapses it
to a false homozygote -- the mechanism behind chr20:48,204,383 (documented in
`b44fd96`, `01b4732`). It is live only inside `opts.retry_windows` or when
`opts.joint_het_orientation` is set, and nothing in the shipped
`collect-bam-variation` CLI ever sets either, so the fix is dead in the
default alignment path.

A second, related failure was found independently this session:
chr20:576,038 (`A>AT`, a 1bp insertion into a 12bp A-run, DP 60, 46 ref / 14
alt) is a homopolymer indel, excluded from `init_assign_read_hap`'s read
scoring entirely. Its `hap_to_alle_profile` still accumulates real votes
via `update_var_hap_profile` (which has no homopolymer guard) --
hap1=[33 ref, 1 alt], hap2=[13 ref, 13 alt] -- but hap2's vote is an exact
tie, and the argmax tie-break prefers reference. Both haplotypes land on
REF, so `hap_alt=hap_ref=0` and the site is silently dropped from the
phased VCF entirely, despite 13 of its 14 alt reads sitting in one
haplotype's bucket. `allele_depths_call_het` would also rescue this case:
its predicate (ref_cov/alt_cov >= min_alt_depth, AF in [min_af,max_af]) is
satisfied and does not itself check `is_homopolymer_indel`.

Losing candidates like these upstream of a phase set removes exactly the
kind of linking evidence that keeps the greedy, iterative per-read k-means
solve anchored to one consistent orientation. On whole chr20 (default
build, `ac59de2`), truth-scored against parental origin (HG002 trio, read
truth derived per `scripts/make_truth_hap_map.sh`), 4 of 295 scored phase
sets showed a genuine internal orientation switch -- most cleanly
PS=130540 (807 clean-SNP-anchored genotypes, switches once between
656,411 and 658,439, ~2kb apart, confirmed against longcallD which does
not switch there).

The fix tested: set `opts.joint_het_orientation = true` only inside
`collect_bam_variation()` (`collect_pipeline.cpp`), after CLI parsing,
leaving the `Options` struct default (`false`) untouched -- so the graph
arm and every unit/predicate test that constructs its own `Options`
directly is unaffected.

## The measurement

Unit tests, predicate tests, and window tests all pass unchanged (125/125
window assertions, including the graph-arm window that regresses when the
struct default itself is flipped instead of scoping the change to the
CLI). Locally, both target sites are fixed exactly as predicted:
48,204,383 -> `0|1` with its own phase set; 576,038 -> `1|0`, phased into
PS=130540; the PS=130540 switch itself resolves (656,411 and 658,439 both
read MATERNAL-on-hap1 after the fix).

Chromosome-wide the picture reverses. Re-running the same truth scan
(50,821 clean-SNP-anchored genotypes, same purity/depth thresholds) finds
**14 of 285 scored phase sets switching, not 4** -- `switches_baseline_ac59de2.tsv`
and `switches_with_joint_het_orientation.tsv` are the two full per-site
scans. PS=130540 leaves the list, as do the other three originally-found
switches (one of them, PS=26549599, is the same thin-evidence site as
before -- 4 scored SNPs -- and stays a switch either way). Ten new phase
sets switch that did not before: 3883778, 6992739, 8106314, 19502530,
20617179, 20782613, 21817199, 23007537, 37413502, 55824242 (plus
57104654, 61405579, 61690751 -- 13 new in total against 3 fixed). Several
are lopsided in a way the earlier four were not -- PS=19502530 is 937
PATERNAL against a single MATERNAL outlier at its first scored site,
PS=23007537 is 439 against 2 -- consistent with the joint-orientation
branch now also mis-resolving isolated noisy sites it previously left
alone (correctly unresolved, in effect) rather than only rescuing the
collapsed ones it was aimed at.

Net: 3 fixed, 13 broken. Reverted. `git diff` is empty on `src/`.

## Why the broad predicate is not safe as a default

`allele_depths_call_het`'s AF/depth gate is necessary but not sufficient to
tell a genuinely heterozygous site from a noisy one whose independent
per-haplotype vote merely landed on a plausible-looking split by chance --
exactly the same class of failure it is meant to fix, just in the opposite
direction. This matches the shape of every other attempt at widening this
admission in the project history (`6cea093`'s noisy-MSA branch gate,
`01b4732`'s own retry-window scoping): the fix is real for the sites it
targets and not free elsewhere, and "elsewhere" is not visible to the
window-test panel because that panel is a curated dozen-ish loci, not a
chromosome-wide scan. The window-test panel passing is necessary but not
sufficient evidence for a change like this one.

## Status

Not shipped. The two root causes (48,204,383-style independent-majority
collapse; 576,038-style homopolymer-excluded-from-scoring-but-not-from-voting
tie-break) are both confirmed and reproducible, and `allele_depths_call_het`
is a real fix for both in isolation, but applying it wherever its predicate
holds is a net-negative change chromosome-wide. Any real fix needs either a
narrower trigger (only apply the joint-orientation resolution when it does
not also break an already-consistent neighbor -- i.e. verify against the
existing link chain before committing, not just against the site's own
depths) or per-site scoping (something closer to `retry_windows`, but
populated automatically rather than left empty by default). Neither is
built.

Co-authored-by: Claude Sonnet 5 <noreply@anthropic.com>
