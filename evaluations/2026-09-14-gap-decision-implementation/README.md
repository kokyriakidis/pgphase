# Frozen gap evidence and offline decision experiments

## Implemented

`collect-hybrid-variation --recover-gaps --gap-decision-audit DIR` runs a
read-only audit before any production recovery mutates chromosome chunks.
It exports all applicable clean/MSA/indel/homopolymer tiers and both configured
observation views, without stopping at the first join. Production acceptance
still uses the existing rules. Audit workers share frozen inputs, have private
proposal state and write separate gap files. Existing snapshots are not
silently overwritten. A per-contig manifest records the cache signature and
relevant stitching options; the trial runner also records the binary hash.

The per-contig `.members.tsv` freezes every original tagged molecule, including
those outside gap windows; `.blocks.tsv` records original phased-site bounds
and haplotype read counts. Each gap additionally exports:

- `.evidence.tsv`: schema, reference slice, candidate allele representations,
  frozen original read HP/PS, input/molecule identity, mapping flags/quality,
  integrated and graph allele observations, query indices and available BQ.
- `.reads.tsv`: every molecule's proposed HP/PS/skipped state for each view/tier.
- `.votes.tsv`: both original-flank 2x2 vote matrices for every proposal block,
  together with the existing rule's trial-level joined/flip result.

The integrated observation channel can contain graph injection; it is **not**
an independent BAM-only channel. Query index -2 denotes graph confirmation;
BQ -1 is unknown, not low error probability. View 1 restores BAM SNP calls on
graph-selected reads; its MSA observations can still share the same molecules
as view 0. The exports are alternatives, not independent replicates.

`scripts/replay_gap_decisions.py` verifies that molecule assignments reproduce
every exported vote matrix exactly. It reports SAME/FLIP costs, explicit
insufficient/conflicting evidence, per-haplotype margins, molecule influence,
and consistency of tentative connected components. Missing one source is not
counted as opposite-orientation evidence; repeated tiers never multiply votes.
One molecule losing all support is distinguished from reversing orientation.
Contradictory components have no proposed orientation assignments.

Direct-observation diagnostics additionally include untagged bridge reads.
Overlapping event descriptions share one local event; conflicts stay visible.
This overlap grouping is not sequence-equivalence normalization. SNP type 8
matches the native BAM_CDIFF representation, covered by a regression test.

A versioned feature cache separates expensive observation extraction from
score tuning. It retains molecule contributions independently of thresholds
and checks SHA-256 hashes of all frozen source files before reuse. Changing
weighting/margins needs no phasing, MSA, BAM output or observation parsing.
Changing input bytes, targets or feature schema invalidates reuse.

## Experiment already produced a safe production fix

The all-tier chromosome audit explains why 57.84 and 62.41 Mb were withheld.
Their unmodified BAM-view proposals support the same relative orientation on
both original haplotypes. The old confirmation step first backfilled/promoted
sparse homozygous MSA calls, changing the proposal rather than merely checking
its orientation. That promotion is removed, including its now-unused helper
extension. The ordinary MSA observation backfill still leaves HOM calls alone.

The bad 23.42 Mb proposal remains rejected: its BAM view has only one left
haplotype supported by a single read, despite much larger right-side counts.
The bad 37.64 Mb proposal lacks a BAM connection. The existing independent
per-haplotype support requirement and same-parity confirmation remain active.

Full chr20 output `/tmp/pgphase-confirm-original`:

| Metric | Previous accepted | Confirmation fix |
| --- | ---: | ---: |
| Recovered initial gaps | 123/280 | 125/280 |
| Evaluated reads | 189,042 | 189,042 |
| Discordant reads | 2,427 | 2,427 |
| Switches | 351 | 351 |
| Flips | 414 | 414 |

Only 57,841,772–57,866,713 and 62,408,056–62,432,427 are newly joined. No old
edge is lost or changes parity. Every read retains its truth-concordance
status; none are newly tagged or lost. Both known bad edges stay rejected.
No new variant Hamming or shared-callset NGC50 measurement is claimed.

## Validation and measured speed

The initial audit reproduced the prior full chr20 native VCF records and BAM
read identity/alignment/HP/PS records exactly. The same noninterference check
passed all 22 regional arm comparisons (44 BAM/VCF comparisons). C++ tests
verify exported vote matrices and read-only behavior. Python tests cover
missing-versus-conflicting evidence, label gauges, duplicate/corrupted data,
overlap grouping, untagged bridges, per-haplotype imbalance, parity cycles,
cache reuse across thresholds and rejection of changed inputs.

The first chromosome feature extraction took about 70 seconds and produced a
19 MB feature cache. Cached scoring took about 0.2 seconds; loading and hash
verification brought a complete replay to about 2.3 seconds. These are single
runs. The first export plus normal production solve took 405.3 seconds versus
191.2 seconds for production alone; audit is opt-in and amortized across trials.

## Usage

Add this to the usual cached chr20 hybrid command:

```sh
--recover-gaps --gap-evidence-cache /tmp/chr20-evidence-v5.gapev \
--gap-decision-audit /tmp/chr20-decision-audit
```

Then prepare/reuse scoring features:

```sh
python3 scripts/replay_gap_decisions.py \
  --audit /tmp/chr20-decision-audit \
  --feature-cache /tmp/chr20-decision-features.json \
  --output /tmp/chr20-gap-decisions.json
```

Subsequent runs can change `--min-margin`, `--min-haplotype-margin` or the
experimental `--weighting mapq` without rebuilding the evidence. Run regression
tests with `make unit-tests` and `make benchmark-tests`. The regional runner
accepts `--decision-audit` to preserve all exports alongside its usual outputs.

## Not yet validated or enabled

The general score is diagnostic, not a calibrated posterior or an automatic
replacement for production acceptance. Large flank margins can conceal a
single fragile bridge. Current leave-one-molecule checks operate on fixed
proposal assignments; they do not re-solve k-means or rebuild MSA without that
molecule and therefore cannot establish independence from its contribution to
proposal construction. Direct observations only orient against already phased
sites; they do not solve every untagged-read chain through private gap sites.

Complete sequence-equivalent multiallelic-event reconstruction, calibration of
repeat/alignment errors, event-removal re-solving, untouched-region validation
and production adoption of a general joint score remain outstanding. Alternative nonadjacent block edges are not generated yet. The five previously unresolved
original gaps are not claimed fixed. Competitor calls and parental truth remain
diagnostic/evaluation inputs only; they never enter the phaser or replay score.


Final audit concurrency checks compared 54 snapshot files across four control
regions and both existing arms: all files, including full original block
membership tables, were byte-identical with serial and parallel execution.
Reusing an audit directory returned an actionable error and left every frozen
file hash unchanged. The 62.41 Mb small-region rerun does not reproduce its
full-chromosome join, reinforcing why the chromosome-context export is needed.

`replay_ablations.json` preserves three diagnostic arms on the same 280-gap
snapshot. Balanced uniform scoring reports 124 directional candidates,
151 insufficient cases and 5 conflicts. Dropping the per-haplotype requirement
reports 126 directional candidates and includes the known wrong 23.42 Mb
proposal. MAPQ weighting with the balanced requirement reports 123 directional
candidates and does not establish an accuracy improvement. These are proposal
classifications, not accepted production edges or newly resolved gaps.
