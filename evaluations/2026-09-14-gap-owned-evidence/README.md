# Gap-owned evidence validation

The recovery-enabled hybrid pipeline starts with catalog sites, freezes its
graph phase gaps, and projects only gap-owned private BAM/MSA events plus graph
anchors into the existing phaser. Candidate discovery still uses the surjected
BAM; DeepVariant is not an input. Raw discovery/MSA work is cached once.

Implementation: `src/gap_evidence.hpp/cpp`, source channels in
`ReadVariantProfile`, cache format 6, and frozen snapshots in
`recover_hybrid_gaps`. Unsupported complex replacements and events crossing a
gap boundary are audited but excluded from legacy projection. This does not
implement a general equivalence resolver for overlapping or shifted events.

The two audit tables added per gap are:

- `.events.tsv`: event interval, role, reference/alternate dictionary.
- `.observations.tsv`: molecule, separate BAM/graph/MSA-derived calls, status,
  original BAM query position and base quality when available.

## Regional reproducibility

Artifacts: `/tmp/pgphase-gap-owned-region/{cold,warm,baseline,comparison}`.
Exact commands are in each run's `command.json`.

HG002 chr20:52,000,001–53,000,000 against CHM13; 500 kb chunks. Cold cache
with four threads and warm cache with one thread produce identical BAM records,
VCF records, candidates, tier reports and all per-gap audit tables. Four gaps,
one correctly oriented join, five original blocks becoming four. Cold elapsed
28.33 s (cache build 24.26 s); warm elapsed 6.76 s (cache load 0.034 s).

Parental read truth: `../pgphase-eval-data/truth/chr20/diplinator_merged.bam`.
The frozen original membership gives 2,518 phased/evaluated reads and one
discordant read; recovery gives 2,629 and one. One original member is no longer
evaluated after recovery; do not interpret the net gain as exact membership
preservation. Original blocks have uniform orientation transformations on their
common evaluated reads. The other three gaps remain split.

## Whole chromosome procedure

Cold audited run: `/tmp/pgphase-gap-owned-chr20`.
Warm final-code run: `/tmp/pgphase-gap-owned-chr20-warm`.
Cache: `/tmp/chr20-evidence-v6.gapev`.

The new initial inventory is 276 gaps; previous 280-gap results used a different
initial hybrid core and are not a controlled baseline. Reconstruct the actual
pre-recovery read tags from the frozen membership audit:

```bash
python3 evaluations/2026-09-14-gap-owned-evidence/validate.py \
  /tmp/pgphase-gap-owned-chr20 /tmp/pgphase-gap-owned-chr20-baseline
```

This checks every private event's footprint against its owning original gap and
requires every original phased molecule to exist in the output. Evaluate both
baseline and recovery using `scripts/evaluate_phase_accuracy.py`, then use
`evaluations/2026-09-14-gap-bridge-validation/compare.py` with the frozen audit
to assess accepted block parity. Its switch/flip counts are read-concordance
transitions, not variant switch errors. No truth information enters recovery.

## Checks

`make -j8 pgphase unit-tests` and `make benchmark-tests` pass. Tests cover source
disagreements, original BAM coordinates/quality, MSA replacement provenance,
allele-dictionary permutations for observations and anchor consensus, 1/2
genotype counts, immutable projection and gap ownership.

`make check` stops before a comparison: the existing gate script passes the
unsupported `--phased-vcf-output` option. The current CLI accepts
`--phased-vcf-out`. This unrelated gate incompatibility was not changed here.
