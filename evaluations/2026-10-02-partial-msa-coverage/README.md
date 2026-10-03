# Coverage of unplaced reads in haplotype-aware MSA

## Reproduced bugs

The optional recovery pass skipped any read carrying the selected PS and HP,
even when homopolymer consensus construction had excluded that partial read.
Those labels did not establish actual MSA membership. The pass also aligned
all other reads as full-cover reads, scored uncovered ends as deletions, and
gave ambiguous composed paths full reference/read bounds. Even a read with
neither end covered could be assigned to a consensus.

The permanent synthetic regression uses two fixed homopolymer consensuses
and a prefix or suffix fitting both equally. With a one-point assignment
margin, an uncovered repeat length must still provide no haplotype vote.
The final test linked against the starting alignment implementation fails
five sections: labeled prefixes/suffixes disappear, unlabeled prefixes/suffixes
are placed in a cluster, and a no-cover read is placed too. After repair,
partial reads remain available as ambiguous observations with the proper
bounds; no-cover reads abstain. The first default-margin reproduction also
exposes the incorrect full bounds on ambiguous reads.

## Repair

Record which input reads actually receive initial MSA alignments. When the
caller requests unplaced-read observations, use that membership rather than
HP/PS to decide which reads remain. Skip reads with neither end covered.

Reuse the existing full-axis WFA alignment against the two fixed consensuses,
then trim according to each read's original coverage flags. The partial WFA
entry point crops the consensus axis; using it directly would misplace the
later reference composition. Restrict scoring to the covered intersection,
so an absent end cannot act as deletion evidence. Restore coverage bounds
after composing ambiguous reference/read paths, before normalization or site
recall. Full-cover scoring is unchanged.

No rows are merged, no consensus construction or stitching threshold is
changed, and no new alignment method is added. This pass is requested by
targeted graph recovery. Ordinary BAM operation does not request it. Production
uses no truth, competitor calls or fixture coordinates.

## Chromosome comparison

Starting binary SHA256:
`95f94fdd04210a95e295131288ca21398065ac0d0afe1364c7c586fc892e63aa`.
Final binary SHA256:
`74218715db7412bb7d0e88633b14ec6c85e3220047b6e9ce1d6df6210c7d7df4`.
Use the same reference, graph catalog, GAF, BAM, default options and evaluation
truth map as the accepted baseline. Full native outputs are
`test_data/tmp_gap_next27/atomic-quality-snp-gauge/` and
`test_data/tmp_gap_next27/partial-coverage/`.

| Metric | Before | After |
|---|---:|---:|
| Truth-scored phased reads | 237,170 | 237,170 |
| Correct assignments | 230,043 | 230,043 |
| Discordant assignments | 7,127 | 7,127 |
| Conditional read accuracy | 96.994983% | 96.994983% |
| Read phase sets | 661 | 661 |
| VCF variant keys | 63,630 | 63,630 |
| VCF phase blocks | 333 | 333 |
| VCF span N50 | 774,189 bp | 774,189 bp |

Every truth-scored read retains its correctness state. One read's PS changes
from a graph label to an independent BAM label, with HP unchanged and correct
in both blocks. No VCF GT or PS changes, and no variant key is lost or added.
The 31 changed VCF records contain changed counts and derived allele fractions.
For example, the owning 7–8 Mb replay recovers an additional insertion ALT
observation at 7,047,081: DP 80 -> 81 and ALT 32 -> 33, retaining its genotype
and PS. This demonstrates retained evidence, not an additional gap closure.
The 7.264321–7.280346 and 24.121713–24.131707 Mb gaps remain open.

## Verification

Build succeeds with no new warnings; the existing unused-function warning in
abPOA's SIMD header remains. Unit tests pass. The expanded predicate suite
passes 962 assertions in 41 cases. HiFi and ONT TSV/VCF goldens pass, including
HiFi one/four-thread determinism. Starting-code failures are preserved in
`regression-before.txt`; aggregate comparison and individual changes are in
`chr20-parity.json` and `changed-records.json`.

The regression panel first runs 105 fresh native requests with exact recorded
arguments, bounds and thread counts. `native-cache.json` records final-binary,
input-stat and output-file hashes. `replay_cached_panel.py` verifies those
identities before scoring reuse; requests absent from the cache run natively.
Starting-version outputs are never reused. Full chr20 takes 296.22 seconds
at eight threads, and the concurrent panel takes 395.12 seconds with four
workers. These are validation timings, not a controlled runtime comparison.
The final window suite passes 7,038 assertions in 69 cases, preserving the
100-coordinate panel and all 86 required connections. Existing span equality,
parental-orientation checks and read-score expectations are unchanged.
