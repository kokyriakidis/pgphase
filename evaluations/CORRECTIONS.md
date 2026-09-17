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

## The three multiallelic loci do not share one truth verdict

**Claimed**, in commit 80094d3's message: "Truth cannot validate the orientation
at these three loci: they are 1-base-different deletions inside 14-17 bp
homopolymers and read net lengths do not segregate (purity 0.536 to 0.615)."

**Actual**, per-allele, from the measurement that sentence summarised:

| locus | ALT1 | ALT2 |
|---|---|---|
| 55,795,217 | −14, MAT 13 / PAT 15, purity **0.536** | −15, MAT 9 / PAT 13, purity **0.591** |
| 55,815,775 | −15, MAT 2 / PAT 0, purity **1.000** | −17, MAT 7 / PAT 1, purity **0.875** |
| 5,379,662 | −13, MAT 16 / PAT 23, purity **0.590** | −14, MAT 16 / PAT 10, purity **0.615** |

The quoted range covers 55,795,217 and 5,379,662 only. At 55,815,775 the net
lengths **do** segregate — 1.000 and 0.875 — so "read net lengths do not
segregate" is false there and the range is wrong.

That locus is unscorable for a different reason, which the commit message should
have given instead: only 2 and 8 reads carry the two alleles with a truth label,
and **both alleles map to MATERNAL**, which cannot be right for a heterozygote.
So its orientation is unvalidated because the truth-scored support is too thin
to place either allele, not because the alleles fail to separate.

The conclusion the sentence supported is unchanged: at none of the three loci
does truth establish which haplotype carries which allele, so what the fix
removes is a genotype contradicted by its own depths — `1|1` while 28 to 31
reads carry allele 2 — rather than a demonstrably wrong phase.

## Two panel figures were stated without their convention

**Claimed**, in the claim-over-repeat-screen record and commit 1b89d06: "in-gap
phased hets 12 -> 14" and "39 newly concordant reads", under "Panel, stock
defaults".

**Actual:** both numbers are real but measured two different ways, and neither
was labelled.

`score_panel.py` counts hets with the gap bounds **inclusive**, so each window's
two flank anchors count and the panel total starts at 12. `test_gap_windows`
counts **strictly** inside, `pos > gap_left && pos < gap_right`, and reads 0.
The delta is **+2** either way -- `5,315,591` and `12,721,112`, one per window.
Quoting 12 -> 14 next to a harness that reports 0 -> 2 for the same-sounding
metric invited exactly the mismatch it got.

"39 newly concordant" is a set difference: 39 reads are newly tagged and
concordant, while 1 read that was concordant lost its tag, so concordant reads
net **+38** (2,785 -> 2,823) against +37 tagged. Discordant reads fall 9 -> 8.

## "The hybrid still emits nothing there" was arm-specific

**Claimed**, in the injection-tests record and in commit b3bd99a's message, about
`chr20:55,919,945` after the duplicate fix: "The hybrid still emits nothing
there -- that is the `emitted_multi` baseline".

**Actual:** true only where the noisy class is excluded from the solve. With the
class admitted the same record is phased into PS 55,815,793 and emitted as
`55919944 C>CA,CAA GT=1|2`, carrying both alleles. The candidate is present in
both cases with identical counts; only admission differs. Stated without the
qualification the claim is false for a supported configuration -- and since the
re-solve is now the default, it is false for the default one.

The `emitted_multi` allowance rows remain correct as written: they record what
the arm with the class excluded emits.

## RETRACTED: the "three defects the retry-by-default change exposed" were my test's false positives

Commit dcc18a0 recorded three defects that making the re-solve the default was
said to expose. **All three were artifacts of rules in the injection test suite,
not pipeline defects.** Each was checked against the reads afterwards and the
pipeline's output is correct in every case.

**1. `48,243,089 A>TT` is a correct complex record, not a malformed indel.**
The test required a length-changing ALT to begin with REF. That is a convention
for simple indels, not a VCF requirement, and this locus is not a simple indel:
reads at the anchor split **31 A / 31 T**, a genuine heterozygous substitution.
So `alt_ref_base = T` is right and the record says what it should -- the alt
haplotype carries T at the anchor plus an inserted T. The rule flagged a correct
record.

**2. The depths at that locus do not disagree.** Two candidates sit there, and
each is emitted at its own anchor with its own counts:

| candidate | emitted |
|---|---|
| `POS=48243089 INS ALT=TT` DP 64, 39/25 (injected claim) | `48243088 AA>ATT AD=39,25` |
| `POS=48243090 INS ALT=T` DP 57, 25/32 (alignment) | `48243089 A>TT AD=25,32` |

The test matched an emitted record to a candidate by equal POS. An insertion
emits at `POS - 1`, so it paired the record with the neighbouring candidate and
reported its depths as a mismatch.

**3. `55,846,004` is a different event from the claim, not an under-consumed
one.** The catalog claims `CT>CTTT`, net **+2**. Our record `C>CTTT` is net
**+3**. Reads: **+3 on 47 of 74 (63.5%)**, +2 on 17 (23.0%). The majority event
is the one we emit; the claim describes a minority allele. The test's signature
-- a record carrying the claim's ALT while consuming fewer reference bases --
conflates two events that happen to share an ALT string, which is easy in a
homopolymer.

The lesson is about the checks, not the pipeline: two of the three rules encoded
a convention as a requirement, and the third assumed a coordinate convention the
emitter does not use. The suite did find real defects earlier -- the spurious
deletions from substitution claims, the missing strand tallies, the insertion
REF -- so this is not an argument that it was worthless, only that these three
entries were wrong and are withdrawn.
