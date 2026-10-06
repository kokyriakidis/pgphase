# Verified source indels restore connected-core parity

Target: chr20:61,738,239–61,747,506. This interval already had a spanning
variant phase set, but its read core fell short of HiPhase. Its existing panel
entry and exact span/orientation assertions are retained; the certification
manifest now requires the full 80% / total / core contract.

| Output | Correct / all scorable primary overlaps | Connected-core correct |
|---|---:|---:|
| Before | 94/96 (97.92%) | 85 |
| After | 94/96 (97.92%) | 91 |
| HiPhase | 91/96 (94.79%) | 91 |

HiPhase has one discordant and four unphased reads; pgphase has no discordant
and two unphased reads. The comparison checks all 96 HiPhase primary input
alignments for identical start, end, CIGAR, and sequence. Orientation is measured
independently against parental truth over the entire replay / competitor output,
and rescue-offset phase sets do not count as a connected core.

The original independently verified BAM edit could lose its physical identity
when a shared catalog row retained the graph description. Preserve that edit in
source provenance and require the final marker to retain the same edit and MSA
verification before checking already tagged rescues after physical bridges.
Materialization requires a complete source path, two consistently oriented clean
SNP loci shared with the graph, matching graph/BAM observations, and no contrary
phased profile anchor. Pure, isolated edits are required. A partial transfer is
vetoed if the rescue cohort has greater local coverage than its core anywhere
along the cohort. The 41.881 Mb regression demonstrated that otherwise a larger
output cohort can fragment, reducing contiguity. This coverage rule has an
in-memory regression. The transfer creates no HP assignments or variant joins.

The decisive reference read, `m84031_231217_034919_s2/46208283/ccs`, has mapping
quality 60 and one Q17 base in the verified deletion footprint. A blanket Q30
cutoff rejected it even though its summed base plus mapping error is about 2%,
below the 5% per-call bound. Q10 footprints with about 10% error and differing
repeat deletion lengths remain rejected. The fast predicate regression tests
both measured-error cases; the owning-chunk regression retains 94 total and
91 core correct as explicit floors.

## Development cost

Use `LD_LIBRARY_PATH=/usr/lib/x86_64-linux-gnu make gap-dev-check` for the
edit loop. It runs all 1,551 predicate assertions and the adapter fixtures in
0.03 seconds on the current build, with no production-binary build or alignment
replay. `PREDICATE="[deletion]"` selects the Catch2 predicates when useful.
After the logic passes, run `make gap-owner-check GAP=61.738`: 0.31 seconds with
saved state. The measured cold 1 Mb owning replay took 9.81 seconds. Production
binary changes invalidate real pipeline outputs; run broad validation once the
patch is stable. An older binary's output cannot test a production change.

## Final validation

All 84 registered checks (83 mechanisms and the 112-window panel) passed in
195 selected invocations, retaining every nested section. Unit, predicate,
cache-helper, and standard validation gates passed. All 99 spanning panel
intervals remain spanning; 33 now meet the full 80% / total / core contract,
compared with 30 before. Besides the target, 1,086,625–1,110,921 and
59,825,454–59,842,960 gained full-contract parity. All three now enter the
certification manifest and pass its stricter assertions.

The final full chromosome replay took 359.83 seconds. All 64,188 variant
records, including genotype and phase-set fields, are identical; all 256,610
primary output molecules remain present and every HP tag is unchanged. The
transfer moves 173 existing rescues into connected cores: 156 correct and 17
discordant (90.17% correct). Whole-output truth counts remain 230,566 correct,
6,768 discordant, and 19,276 unphased. Per-block parental rescoring has three
correct-to-discordant and three discordant-to-correct transitions; totals remain
unchanged. Existing regression floors and discordance ceilings were retained.

`hiphase.json` records the matched-input comparison. `compare_hiphase.py` and
`audit_full.py` reproduce the audits using the retained test-data artifacts.
`validation.json`, `full-audit.json`, and `gap-contract.tsv` record the frozen
binary, broad checks, and per-window measurements.
