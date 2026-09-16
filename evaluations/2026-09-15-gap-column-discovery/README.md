# Can column discovery find usable sites in the gaps? Not from the linear pileup

Prototype of the pass-2 site-finding step: candidate columns taken straight from
the pileup with no catalog lookup, no classification and no MSA, then admitted by
how well each agrees with the read partition the other columns imply. Scored
against the read-level truth, with in-block windows as controls -- the same
trusted-class control discipline that caught two genotyping bugs earlier.

```sh
python3 evaluations/2026-09-15-gap-column-discovery/discover_columns.py \
  --bam test_data/HG002_chr20_hifi_mapped_to_CHM13_chr20_annotated.bam \
  --windows windows.tsv --truth-map truth_hap.tsv \
  --called-sites /tmp/graph-only-mapq/q1/variants.tsv \
  --site-output discovered_sites.tsv --window-output discovery_windows.tsv
```

| window | kind | candidates | admitted | consistency | reads | truth concordance | median segregation | already called |
|---|---|---:|---:|---:|---:|---:|---:|---:|
| `inblock_1499101` | control | 109 | 103 | 1.000 | 246 | 100.0% | 1.000 | 100 |
| `inblock_514902` | control | 49 | 30 | 0.976 | 123 | 91.1% | 1.000 | 29 |
| `inblock_865232` | control | 39 | 39 | 0.884 | 193 | 88.6% | 1.000 | 36 |
| `gap_35919404` | gap | 20 | 3 | 1.000 | 130 | 71.5% | 0.750 | 0 |
| `gap_36332599` | gap | 9 | 4 | 0.986 | 80 | 100.0% | 0.994 | 0 |
| `gap_13429829` | gap | 7 | 4 | 1.000 | 216 | 51.4% | 1.000 | 2 |
| `gap_60084884` | gap | 2 | 2 | 1.000 | 126 | 52.4% | 0.993 | 1 |
| `gap_10296487` | gap | 2 | 2 | 1.000 | 138 | 52.2% | 0.993 | 1 |
| `gap_16673759` | gap | 1 | 1 | 1.000 | 53 | 100.0% | 1.000 | 1 |
| `gap_7163303` | gap | 1 | 1 | 1.000 | 66 | 100.0% | 1.000 | 1 |
| `gap_7073919` | gap | 0 | 0 | — | 0 | — | — | 0 |

## The mechanism works; the linear pileup does not feed it

The controls validate the method: inside existing blocks discovery finds 39-109
candidate columns, admits 30-103, and those columns have median segregation
1.000 against truth -- with 29-100 of them being sites the pipeline already
calls, so discovery rediscovers the known-good set rather than inventing one.
Self-consistency also tracks quality without truth: the control at consistency
0.884 came out at 88.6% truth concordance, the one at 1.000 at 100.0%.

In the gaps the same code finds **0-20 candidate columns and admits 0-4**,
against 30-103 in blocks. The linear pileup simply does not expose enough
heterozygous columns there.

## A failure mode the design has to guard against

With only 2-4 admitted columns, self-consistency reads 1.000 while the partition
is meaningless: `gap_13429829`, `gap_60084884` and `gap_10296487` all report
consistency 1.000 and land at **51.4%, 52.4% and 52.2% truth concordance** --
chance. A handful of columns trivially agree with the partition they themselves
defined. So consistency is necessary but not sufficient: a minimum admitted-column
count is required before a partition may be trusted (the controls sat at 30+),
alongside the flank-anchoring requirement.

Three gaps did resolve cleanly (`gap_36332599`, `gap_16673759`, `gap_7163303` at
100% concordance), so the approach is not dead -- it is starved.

## Consequence

Answering the question that prompted this: **no, pileup-based discovery does not
recover the missing sites in the gaps.** It recovers them inside blocks, where
they were never missing. The remaining candidate is the graph channel: every GAF
record carries a `cs:Z:` difference string relative to the traversed path, so its
mismatches are private variation by construction, and one gap-sized window
carried 1,006 mismatch, 1,860 insertion and 3,458 deletion events. That has to be
measured the same way -- same windows, same truth scoring -- before it is
believed.

## Aside: `-q 1` closed the gap this investigation started from

`chr20:25,834,662-25,883,079`, the gap diagnosed by hand over a full session, is
no longer a gap. At `-q 1` it sits inside a single 831 kb block (PS 25,102,178,
spanning 25,102,178-25,933,230 with 995 sites), which is why it is absent from
the window list above.
