# A real, root-caused bug in `refresh_assigned_msa_observations`, fixed without the accuracy cost

## The bug

`refresh_assigned_msa_observations` (`collect_phase_noisy.cpp`, called from
`make_vars_from_msa_cons_aln`'s two-cluster branch) exists to correct
cross-cluster count contamination in `NoisyCandHet` candidates (its own
comment: "with two clusters of two reads each carrying their own consensus
the site is 2 ref / 2 alt, and without the refresh it was counted 3 ref /
1 alt"). It does this by unconditionally resetting every such candidate's
`alle_covs` to all-zero, then rebuilding both the per-read profile entries
and the aggregate counts from a from-scratch per-read reclassification,
`call_local_msa_allele`.

That classifier requires a *consistent, valid* allele call from both
haplotype consensuses before it will attribute any read at all
(`alleles[0] < 0 || alleles[1] < 0 || ... return -1`). That holds when both
haplotypes carry the same kind of event -- the ordinary case -- but not
when a homopolymer is unstable enough that the two consensuses disagree on
event TYPE, not just allele: one haplotype an insertion, the other a
deletion, relative to reference, at the same run. When that happens, every
read fails classification, and the unconditional reset plus
`update_read_var_profile_with_allele(vi, -1, ...)` (an unconditional
overwrite, not a conditional one) both zeroes the aggregate counts AND
destroys the correct per-read entries `update_cand_var_profile_from_cons_
aln_str2` had already established moments earlier -- discarding real data,
not just failing to improve on it.

This is not a porting question. `refresh_assigned_msa_observations` and
`call_local_msa_allele` (flank-slicing, edit-distance-based local allele
calling) have no counterpart anywhere in longcallD's source -- this is
entirely this project's own added mechanism, and the bug lives entirely
inside it.

## Measured, one locus

chr20:411,654 (`T>TTT`, a 1bp insertion into an 18bp T-homopolymer;
sibling haplotype consensus at the same run is a 1bp deletion, not an
insertion): before any fix, `total_cov=0`, `hap_alt=hap_ref=0`, and the
site is silently absent from the phased VCF. longcallD calls it
`0|1:51:35,16:0.314:60:130540` -- DP 51 matches the raw BAM overlap count
exactly (51 reads, all comfortably above the MAPQ floor).

## Three attempts

**Fix 1 (aggregate-only):** save `alle_covs`/`total_cov` before the reset,
restore them if the rebuild leaves the variant at zero. Verified at
411,654: DP/AD/AF now match longcallD exactly. But `hap_alt`/`hap_ref`
stayed `0/0` and the site still never reached the VCF -- the per-read
profile entries (which the actual haplotype-resolution step reads, not
the aggregate counts) were still being clobbered with -1.

**Fix 2 (per-read skip on failure):** in the per-read reclassification
loop, skip the call to `update_read_var_profile_with_allele` entirely
when `call_local_msa_allele` returns -1, instead of writing -1 over the
existing entry. Fixed 411,654 exactly (`0|1:51:35,16:0.313725:60:130540`,
matching longcallD including GT and PS) and, chromosome-wide, raised
record identity from 99.1% to 99.58%. But it also regressed read-phasing
accuracy from 99.28% to 98.45% (discordance more than doubling, 1,570 ->
3,399 reads) -- reverted, and the reason matters: applying the fallback
*per read* means that for a variant where the refresh partially succeeds
(most reads classify fine, a handful don't), those handful keep their
old, potentially cross-contaminated values. That undoes the refresh's own
purpose for every variant with even one stray unclassifiable read --
which is most of them -- to fix the rare fully-unclassifiable one.

**Fix 3 (shipped -- per-candidate, not per-read):** a pre-pass determines,
per candidate, whether `call_local_msa_allele` returns a valid call for
*any* read at all. Only candidates with zero valid calls anywhere are
skipped entirely (left exactly as `update_cand_var_profile_from_cons_aln_
str2` established them, matching Fix 1+2's effect for the fully-broken
case). Every other candidate goes through the refresh exactly as before,
per-read fallback included -- so the contamination-fixing behavior this
function exists for is fully intact everywhere it was already working.

## Whole chr20, fix 3 on top of both already-shipped fixes

Record-level parity against longcallD's own output:

| | our records | identical | ours-only | upstream-only |
|---|---:|---:|---:|---:|
| shipped (anchored_stage2 + merge off) | 117,896 | 116,868 (99.1%) | 1,028 | 1,402 |
| + fix 2 (per-read, reverted) | 118,099 | 117,596 (99.58%) | 503 | 674 |
| **+ fix 3 (per-candidate, shipped)** | 118,038 | **117,009 (99.13%)** | 1,029 | 1,261 |

Smaller than fix 2's gain -- most of fix 2's record-parity win came from
partially-classifiable candidates being fully rescued, which is exactly
the part that cost accuracy. Fix 3 only recovers the genuinely-broken
(zero-valid-calls) subset, which is real but smaller.

Read accuracy against the diplinator truth BAM:

| | accuracy | discordance |
|---|---:|---:|
| shipped (anchored_stage2 + merge off) | 99.28% | 0.72% (1,570/217,638) |
| + fix 2 (per-read, reverted) | 98.45% | 1.55% (3,399/219,589) |
| **+ fix 3 (per-candidate, shipped)** | **99.27%** | **0.73% (1,595/217,732)** |

No meaningful change from the pre-fix baseline (25 more discordant reads
out of 94 more evaluated -- within noise). The regression fix 2 caused is
gone.

## Status: shipped

`collect_phase_noisy.cpp` (`refresh_assigned_msa_observations`), not
scoped to a CLI path -- this is a pure bug fix, not a representation or
algorithm choice, so unlike `anchored_stage2` and `merge_colocated_msa_
alleles` there is no graph-arm tension to scope around. `make window-tests`
(125/125, graph-arm windows included) confirms that directly. Unit and
predicate suites also pass unchanged.

Co-authored-by: Claude Sonnet 5 <noreply@anthropic.com>
