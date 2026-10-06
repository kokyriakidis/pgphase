# Unified gap integration suite

All 84 former integration cases now run as named sections of `all gaps`, sharing
one fixture and selector. Three focused gap unit cases remain independent.
The committed panel still contains all 112 gaps. Special owning-chunk, allele,
matrix, physical molecule and orientation regressions retain their assertions.

72 owning-context overrides moved from nested conditionals into
`src/test_gap_replays.tsv`. Compiled comparison against the old expressions
confirmed identical regions for 216 override and neighboring/default examples.
38 original literal read bounds from nine owning regressions moved unchanged to
`src/test_gap_read_floors.tsv`. All other 1,163 original assertion sites remain
verbatim. No production behavior or window expectations changed.

HiPhase was measured for all 112 gaps after verifying identical primary input
alignments (position, end, CIGAR and sequence) for every truth-scorable overlap.
Measurements include unphased reads. The preparation helper caches measurements
against its source and input/competitor/truth/panel identities plus output content.
A reused measurement does not rerun any test assertions.

Native scoring additionally records correct primary overlaps using whole-block
parental orientation and excludes rescue phase sets from dominant core counts.
The three reviewed closures in `src/test_gap_certified.tsv` enforce >=80% primary
correctness, total correct >= HiPhase, and dominant core correct >= HiPhase.
Historical gaps retain their original regression gates: 30/99 closed gaps meet
the strict contract and 69 have a shortfall. Those gaps are a visible backlog,
not newly certified successes; see historical-shortfalls.tsv.

`make window-tests` passed 13,232 assertions in four Catch2 cases, including all
84 named integration checks and all 112 panel gaps, in 32.45 seconds. All 114
pipeline completion states stayed unchanged. One selected regression took
0.289 seconds. Unmatched filters fail, and missing fixtures produce an explicit
skip while focused unit cases still run. Four replay cache tests, two competitor
state tests, production build, C++ unit tests, make check and 47 predicate cases
all passed. Compiler output contains no new warnings.

Run everything: `make window-tests`.
Select a mechanism: `PGPHASE_GAP_FILTER=clean-snp-source-retry make window-tests`.
Refresh HiPhase evidence: `make gap-benchmark HIPHASE_BAM=... BENCH_PYTHON=...`,
where BENCH_PYTHON has pysam installed. Review any changed competitor evidence.
`PGPHASE_GAP_BENCHMARK` and `PGPHASE_GAP_CERTIFIED` override the respective
manifests. Fixture-changing runs need benchmark evidence for those fixtures.

The per-run gap-contract.tsv reports every panel window. Legacy tags are selector
metadata in the unified runner; use PGPHASE_GAP_FILTER for former case names or
tags rather than the removed standalone integration case names.
