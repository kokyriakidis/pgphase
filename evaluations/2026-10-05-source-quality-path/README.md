# Quality-backed source path closes 13,752,640–13,773,452

The accepted graph output leaves this 20,812-base interval uncovered by any
older block, although an independently phased BAM source crosses it. The
owning 13–14 Mb replay and the full chr20 trial close exactly this new gap.
No parental labels enter production decisions.

## Evidence and behavior

The initial failed graph edge is the exact shared clean-SNP pair
13,752,640 C>T / 13,773,452 G>C. Two distinct primary MAPQ-60 molecules
call REF/REF and agree with original source HP1. Their base-quality pairs
are 40/40 and 27/40. A fixed Q30 floor discarded the second molecule;
allowing Q20 while retaining the original error-bound test exposes the
independent support. The product of conservative per-molecule mapping and
base-error bounds is 4.23647e-7. No other qualifying canonical pair conflicts.
`bridge-evidence.json` records both reads and original source assignments;
`inspect_bridge.py` reproduces this check from the original BAM and source
dump.

That edge alone is insufficient: the graph path next fails at
13,773,452–13,804,846. There is no callable high-MAPQ GAF pair there, but
intermediate source variants form an independently supported BAM chain. Both
endpoints are exact adoptable shared clean SNPs from the same original source
and agree on its orientation. Every weak cut in that original source has
independent quality support. The missing-edge certificate cannot override
an agreeing, one-haplotype or reversing callable graph pair. The first edge
requires a stricter physical error bound of 0.001. All new certificates are
rolled back unless the complete remaining graph path passes. Source weak
cuts remain recorded; no unrestricted whole-source transfer is enabled.

The existing physical internal-switch repair retains its Q30 requirement,
two molecules per haplotype, count and likelihood bounds, candidate/read
suffix flips and complete remaining-path gate. The new source certificates
agree with the existing orientation and do not perform suffix flips.

## Measurements

The owning replay changes 11 VCF blocks to 10. All 2,835 correct / 179
previously discordant scored reads retain their status. The gap contains
175 truth-scorable reads; 151 are phased correctly before and after, while
the dominant correctly separated count rises from 76 to 151 (43.4% to
86.3%). Disjoint groups observing only the left or only the right endpoint
have 74/74 and 115/115 consistent parental votes in the same orientation.

The full trial preserves all 256,610 tagged output reads, 237,330 scored
reads, 230,563 correct / 6,767 discordant reads, 64,188 variant keys, genotype
alleles, nonphase VCF fields and uniform old candidate-block gauges. It
changes 327 VCF blocks to 326 and 649 read phase sets to 647; span N50 stays
856,770. Exactly 368 read tags and nine VCF rows change. No old tracked span
reopens. All 95/108 previously spanned panel windows remain spanned. Adding the
newly closed gap makes that 96/109; it was not an old panel window.

Ten formerly output-only reads enter the now-certified core block. All ten
carry the independent original BAM source HP: the two physical bridge reads
and eight right-source reads. The other 32 changed output-only tags retain
one rescue cohort, uniformly reoriented into the joined block's gauge. This
is an explicitly measured change in core membership; it does not change any
read's parental correctness. `evaluate.py` rejects lost calls, mixed old
candidate gauges, new discordance, unrelated rescue moves, opposite rescued
HP gauges, extra new gaps and reopened old spans.

The Q20-only replay and a trial that certifies the first edge without the
source chain have identical owning outputs to baseline. A full trial of the
latter also changes no calls (`full-source-loop-parity.json`). Thus lowering
the quality floor or skipping the complete-path check alone is not the fix.

The committed panel adds the new gap with measured `spans=1`, concordance
floor 0.94 and separated floor 0.86. Its replay uses the owning chunk; the
separate owning regression checks absolute read correctness, discordance,
parental orientation and the identical REF/ALT gauge of the shared anchors.
Inputs, binaries and output paths are recorded in `manifest.json`.

Two prior owning-chunk regressions explicitly expected this earlier edge to
remain split. Their separation assertions are updated to equality for this
measured graph-arm improvement; their existing suffix, genotype and parental
checks remain. No old read-accuracy floor or unrelated negative span changes.

The user permits closing a gap when at least 80% of its reads are correct.
The new gap exceeds this even counting abstentions: 151/175 (86.3%) of all
truth-scorable overlapping reads are correctly separated into one block,
and all 151 phased reads are correct. Existing stronger measured regression
floors remain in place.

## Matched HiPhase comparison

`compare_hiphase.py` compares the same 175 truth-scorable original molecules
against the existing DeepVariant/HiPhase output, verifying all 175 alignment
coordinates, CIGARs and sequences match after contig renaming. HiPhase and the
fixed graph arm correctly phase the same 151 reads in one block, leave the
same 24 unphased and have zero local discordance. There is no HiPhase-only
correct read in this gap. Before the fix, our 151 correct reads were spread
across three phase sets. HiPhase's advantage was continuity, now recovered
by retaining the Q27 physical bridge under the error-bound test and the
supported intermediate BAM chain when GAF has no direct pair. The panel
records the measured `hiphase_dv` separated fraction 0.8629.

`hiphase-variants.json` records its eight phased heterozygotes from
13,751,618 through 13,804,846. All eight are already present in our original
BAM source with the same relative genotypes, including the two intermediate
SNPs at 13,787,628 and 13,792,302 and the deletion markers. Our graph-path
validator considers clean graph SNPs, so it demanded a direct pair across
13,773,452–13,804,846 instead of using these intermediate source observations.
This is an evidence-consumption gap, not a need to obtain additional HiPhase
input variants. The source-chain certificate uses the existing independent
source path under the explicit missing-GAF and complete-path checks.

## Final validation

The frozen final binary passes the full-chromosome strict audit and exactly
matches the trial's read tags and VCF rows (`full-final-audit.json`,
`final-trial-parity.json`). No new compiler warning is introduced.
`make unit-tests`, all 47 phase-predicate cases / 1,536 assertions and
`make check` pass, including HiFi/ONT golden outputs and thread determinism.
All 109 window-panel rows and all 83 native regression cases are validated,
covering 10,115 assertions with the updated expected join. The first sweep
passes all 81 unaffected cases and identifies only the two old split
assertions for this exact newly closed edge. Rerunning those two updated
owning cases passes all 282 assertions. `native-validation.json` records the
coverage and the deliberately changed graph-arm assertions.

## Commit review: preserve the original required-site checks

Commit review detected a pre-existing truncation of the required-site file.
Restore all 1282 original graph required-site entries and retain the seven
new source/switch witnesses. The existing native window checker passes all
109 windows / 2,469 assertions against the restored file, rescoring the
hash-bound immutable final-binary outputs. No original required-site entry is
removed. The replay manifest and helper are in `test_data/tmp_gap_fix60/`, and
`native-restored-required.log` records the result.
