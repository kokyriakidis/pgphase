# Close chr20:17,634,393–17,667,022 with a calibrated shared deletion

The interval is open in the task-start binary and closed by HiPhase. The fix
joins the existing variant blocks and puts **164/202 original primary reads
(81.19%) in one correct connected core**, including unphased reads in the
denominator. This matches HiPhase's total and connected-core correct counts;
the correct qname sets are identical. All 202 competitor start/end positions,
CIGARs and sequences match the original input alignments.

| Measurement | Before | After | HiPhase |
|---|---:|---:|---:|
| VCF spans the interval | No | Yes | Yes |
| Correct original reads | 164 | 164 | 164 |
| Correct reads in one core | 110 | 164 | 164 |
| Discordant original reads | 1 | 1 | 0 |
| Unphased original reads | 37 | 37 | 38 |

The fix changes phase-set connectivity, not read assignments or variant calls.
pgphase's existing discordant read remains; HiPhase abstains on that read.
The native owning replay is **17,000,001–18,000,000**. It retains all 999 variant
alleles, counts and filters, all 4,047 correct assignments and every previously
phased assignment. It changes 123 variant gauges and 592 read tags through a
consistent whole-block gauge swap, with no new read assignments.

## Why HiPhase could connect the interval

HiPhase retains a deletion at 17,642,038 and the flanking SNPs in one phase set.
pgphase already has that physical evidence, but its shared graph row retains
repeat classification rather than the source MSA flag. The earlier shared
marker bridge also demanded two clean SNPs per calibrating molecule, although
only one clean SNP lies close enough to this deletion. Separately, GAF omits
intermediate source variants, and a recalled SNP at 17,487,837 lacked retained
MSA provenance. These exclusions made both graph paths appear unsupported.

The fix retains source SNP MSA edits. A recalled SNP must match its canonical
catalog SNP and immutable Q30 physical calls, with two distinct independent
clean source SNP loci per molecule, both haplotypes and unanimous fair-parity
probability <=0.001. A cut-free source may then supply absent or one-sided
agreeing graph edges. Contrary graph calls, weak source cuts and quality cuts
still veto this route. A source prefix can compose with the existing Q30
physical edge proof, with at least two primary molecules, wrong-parity bound
<=0.001 and a complete suffix. All certificates roll back if the whole graph
path still fails.

This additional path route applies to both flanks of a retained shared catalog
deletion of at most 32 bases opposite a clean graph SNP. That marker context
is required: enabling it for all cut-free sources incorrectly unlocked the
distinct-deletion 37.598 Mb seam and lost 156 formerly correct assignments.
The final rule leaves that seam unchanged. The existing distinct-deletion and 37 Mb orientation regressions already
replay the exact native 37,000,001–38,000,000 owner.

The shared deletion's physical edit begins at 17,642,039. Seven REF and sixteen
ALT primary molecules calibrate it against a clean graph SNP, without contrary
votes. Each calibrating call has base qualities >=20 and a summed SNP,
deleted-footprint and twice-mapping error <=1%. Both haplotypes and alleles
must be represented, with two-sided association p<=0.01. The cross-seam
molecule `m84031_231217_062403_s3/147985967/ccs` calls deletion REF at Q40 and
the right SNP REF at Q17 (MAPQ60). Its crossed log odds are 3.87891, or about
2.0254% wrong parity. Physical bridge calls use the existing SNP-to-deletion
5% error bound, and all qualifying pairs must agree. Calibration's one-sided
95% Wilson discordance bound, conservative 1% physical gauge error and bridge
error together must be <=20%. This replaces the separate fixed 10% calibration
cutoff for this bounded shared-marker bridge. General BAM fallback and
masked-SNP witness gates retain their 10% bounds.

The conditional runtime error bound does not replace the measured requirement
for >=80% correct among all truth-scorable original primary overlaps, nor the
HiPhase total/core comparison. The deferred union preserves rescue assignments
and the complementary one-base/twelve-base source deletions at 17,634,393.
No private catalog deletion is forced into the emitted VCF.

Disjoint pgphase flanks score 175/178 and 145/145 in one parental orientation.
HiPhase scores 179/180 and 145/145 with an equivalent global gauge. The panel,
native replay map, required markers, strict certification manifest and owning
mechanism regression all include the new gap. Read floors also preserve its
4,047 correct reads and absolute ceiling of 29 discordant reads in the native
owner, with >=99% concordance. Twenty-four fast
fixtures cover recalled SNP provenance and shared-deletion calibration.

Task-start binary SHA256:
`f98aead0dfc6f9e72a2e2cb99de3463234523d5be7437859cecb095a2259e6fd`.
Final binary SHA256:
`96e95ef60a81898b58989268e4eecf7f758766376cfd000121caf9d09d5f8b70`.

See `hiphase.json`, `orientation.json` and `owner-audit.json` for measured
read sets, parental checks and the native preservation audit. Final full-panel
and chromosome validation is recorded in `validation.json`, `panel-audit.json`,
`full-audit.json`, `full-closures.json` and `gap-contract.tsv`.

Final validation passes all **86 registered gap checks**, including 85 mechanisms
and all 114 panel windows. The 199 selected invocations pass, and the official
`make window-tests` passes 13,411 assertions in four test cases (three unit
cases and the unified gap suite) in 33.24 seconds with native states cached.
The strengthened new regression passes 549 assertions and fails nine with the
task-start binary. Build, all unit tests, 47 predicate cases (1,551 assertions),
golden validation gates, four cache-helper tests and two benchmark-helper
tests pass. Warm `make gap-dev-check` takes 0.05 seconds. All 114 committed
HiPhase benchmark rows match a fresh measurement on identical alignments.

The panel now has 101/114 measured spans and 36 windows passing both >=80%
all-original-read correctness and HiPhase total/core parity. All **332 matching
native replays** preserve variant alleles/counts/filters and every previous
correct or phased assignment. The complete chr20 replay takes 413.56 seconds
with other integration work running concurrently. It confirms the same
164/202 connected-core result and identical correct qnames to HiPhase.
`full-closures.json` verifies that the target is the only new VCF closure and
every previous extent survives. `full-audit.json` preserves all 64,188 variant
records and all 230,567 correct assignments among 256,610 primary outputs;
there are no lost phased assignments or new read assignments.

Reproduce the focused regression with
`make gap-owner-check GAP="17.634"`; use `make gap-dev-check` during edits and
`make window-tests` for the complete panel. The comparison and audit scripts
in this directory use the recorded frozen outputs; `compare_hiphase.py --full`
checks the chromosome output against the same competitor.

Chr20 VCF phase-block N50 remains **856,770 bp** before and after the closure;
HiPhase is **1,005,183 bp** under the same definition. The valid VCF phase-set
count falls 322 to 321. These are inclusive first-to-last phased-heterozygote
spans for blocks with at least two rows, without truth-based switch splitting.
See `block-n50.json` for input paths and complete measured statistics.
