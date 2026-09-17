# A regression test framework for the chr20 gap windows

## Why these are integration tests

Every defect this session found in the gap windows was invisible to the unit
tests, and most were invisible to read-level accuracy too:

| defect | what a unit test saw | what read accuracy saw |
|---|---|---|
| phantom site bridging the k-means, inverting a 13-site flank | nothing | nothing -- 100.00% |
| a homozygous verdict blocking the correction that would overturn it | nothing | nothing |
| a duplicate record beside a correctly merged multiallelic locus | nothing | nothing |
| a join across a 41.8 kb interval no read crosses | nothing | nothing -- 100.00% |

They are properties of the pipeline's output on real data, so the tests that
assert them run the pipeline. One window costs about 1.2 s, which puts the whole
panel at roughly 25 s: a test run, not a benchmark script.

## What is asserted

`src/test_gap_windows.cpp`, built against the Catch2 v2.13.10 single header
vendored in `third_party/catch2/` next to the other vendored dependencies. One
section per arm and window over the committed panel, and a second test case for
the panel totals.

| quantity | assertion | why |
|---|---|---|
| `spans` | **equality, both directions** | a span appearing where none is expected is a coin flip across an interval no read crosses. `chr20:48,176,830-48,229,446` was once reported CLOSED at 100% read accuracy with its halves on opposite haplotypes. |
| in-gap hets | floor | phased heterozygotes strictly inside the gap -- the quantity the retry moves. The default phases only the boundary sites, which is why its floor is 0. |
| tagged | floor | reads carrying `HP`, i.e. coverage gained. |
| concordance | floor | scored per phase set, since `HP` labels are only meaningful relative to their own set. |
| discordant | ceiling | the absolute count as well as the rate: a rate can be held up by coverage while reads go wrong. |

Per-window floors are generous on their own, so the totals case asserts the
panel-level numbers the work is actually steered by.

## Truth is derived, not committed

The truth BAM is aligned to the parental assemblies, so a read's haplotype is
the contig it aligns to -- there is no tag to read and no CHM13-coordinate
region query. `scripts/make_truth_hap_map.sh` does one pass (~12 s) and writes
`test_data/derived/`, applying the same rule `scripts/evaluate_phase_accuracy.py`
applies: `-F 0x904`, haplotype from the `_MATERNAL`/`_PATERNAL` suffix with the
`HO:Z:` fallback, the `hq:i:` floor, first qualifying record per read wins. The
rule is shared on purpose so a window test and a whole-chromosome evaluation
cannot disagree about truth.

The MAPQ floor matters and defaults to 0: a floor of 1 drops 804 records whose
placement on the parental assemblies is ambiguous, and excluding them here while
the chromosome-wide evaluation keeps them would make the two disagree about 0.3%
of reads. The generated map is byte-identical to the one every measurement this
session used -- 272,016 reads, zero disagreements.

## Expectations are emitted by the code that asserts them

`src/test_gap_windows_expect.tsv` is written by the test binary itself
(`PGPHASE_EMIT_EXPECTATIONS`), so the numbers cannot drift from the measurement;
a second implementation of the scoring could disagree with the first and neither
would be obviously wrong. Current state:

| arm | spanned | in-gap hets | tagged | concordance | discordant |
|---|---:|---:|---:|---:|---:|
| default | 0 of 6 | 0 | 2,794 | 99.67% | 9 |
| `--retry-unphased-with-bam` | **4 of 6** | 23 | 3,149 | 99.49% | 16 |

Refreshing the file is how a regression gets committed, so
`scripts/refresh_gap_window_expectations.sh` exists to make that deliberate and
the commit is expected to say which arm moved and why.

## Verified by breaking it, not only by passing

| injected fault | result |
|---|---|
| a window expected to span stops spanning | fails on that window, `false == true` |
| a window spans where none is expected | fails on that window, `true == false` |
| tagged floor raised by 20 | fails naming the window, `543 >= 563` |
| truth map absent | skip naming the file and `scripts/make_truth_hap_map.sh` |
| expectations absent | skip naming the file |

Two defects in the framework surfaced that way rather than from inspection:

- A concordance floor printed with rounding-to-nearest sat **above** its own
  measurement (0.9976 for a measured 0.99757869), so a freshly generated
  expectations file failed against the run that produced it. Floors are now
  rounded down.
- A Catch2 `INFO` stays active for the remainder of a section, so every later
  failure in that section wrongly reported the expectation row as missing. That
  path uses `FAIL` directly now.

Both are the kind of defect a framework hides if it is only ever run green.

## Running it

```
./scripts/make_truth_hap_map.sh     # once, ~12 s
make window-tests                   # ~25 s
./test_gap_windows "[totals]"       # tag filters work as usual
```

`make unit-tests` is unchanged at 4/4, and no pipeline source is touched by this
change.
