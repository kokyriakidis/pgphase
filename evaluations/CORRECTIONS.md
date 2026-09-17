# Corrections

Claims made in this project's session narration that were wrong or overstated,
with the authoritative value. They are recorded here because the documents they
belong to were either already correct or do not contain the claim: a reader
comparing a summary against these writeups should know which statements did not
survive checking.

Each entry names the quantity, the claim, and what the underlying tool output
actually says.

## The 9.89 Mb opportunity was not measured whole

**Claimed:** "Most of that 9.89 Mb is phantom", as the headline of a result.

**Actual:** 9.89 Mb is the sum of **two** revised-cause buckets --
`repeat_indels_would_bridge` (88 gaps, 6.45 Mb) and `low_af_sites_would_bridge`
(24 gaps, 3.43 Mb), which is the 112 gaps the figure is quoted over; a third
bucket, `no_linkage` (52 gaps, 5.32 Mb), is separate and not part of it. The
verification run tested only the 88-gap, 6.45 Mb repeat-indel subset; the
24-gap low-AF subset was never touched. `gap_targeting.md` scopes its own headline correctly to
"0.43 Mb of the 6.45 Mb"; the narration did not.

## The proposal's substantive block was not 388-438 reads in every tier

**Claimed:** "In every tier the proposal ends up with two phase sets: a 20-read
remnant and the substantive block of 388-438 reads."

**Actual:** the per-tier breakdown gives 224 reads in tier 1 (PS 36381019), 438
in tier 2 and 388 in tier 3. Tier 1 falls outside the quoted range. The
conclusion it supported -- no tier joins, and the left block never merges into
the dominant right block -- is unaffected.

## "7 reads spanning the hole" was a count of sites, not reads

**Claimed:** "The 7 reads I measured spanning that hole carry allele evidence
but no tags."

**Actual:** the only 7 established at that point was
`verdict['gap_het_informative'] = 7`, the number of the interval's 11 het
**sites** that segregate at or above 0.90 against truth. No output reported
seven reads spanning anything. Sites and reads were conflated in one sentence.

## The insertion's net length was +13, not +14

**Claimed:** a table gave the record at `5,339,369` (`REF=C`,
`ALT=TGTGTGTGTGTGTG`) a net length of `+14`.

**Actual:** under the convention that same table used for its two deletions --
net = `len(ALT) - len(REF)` with `.` read as length 0, giving `-4` for `REF=ATAT`
and `-5` for `REF=CACGC` -- the insertion is `14 - 1 = +13`. The `+34 on 21
reads` reported for this locus in `evaluations/2026-09-16-best-chain/README.md`
is a different quantity (the net length truth assigns the haplotype) and is not
affected.
