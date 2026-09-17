# The window under active work

The panel is `chr20:5,309,406-5,345,085` -- one window, the one whose gap the
re-solve closes. The five other windows and the panel-wide injection suite were
removed to keep each iteration to a few seconds against a single target; both
are in git history and restorable:

```
git checkout dcc18a0 -- evaluations/2026-09-16-test-panel/panel.tsv
git checkout dcc18a0 -- src/test_bam_site_injection.cpp src/test_bam_site_injection_allow.tsv
```

`make window-tests` now runs in about 4 seconds (two arms over one window, 38
assertions).

## What it asserts

| arm | flags | spans | in-gap phased hets |
|---|---|---:|---:|
| `default` | none -- the re-solve is on | **1** | 2 |
| `noretry` | `--no-retry-unphased-with-bam` | 0 | 1 |

Plus, in code rather than in the expectations file so they cannot drift: no
clean het inside the gap may be left without a phase set, no block may switch
across the gap (each end placed independently on at least five scored reads at
at least 90% agreement), and read concordance stays above a loose 0.95 floor.

`src/test_gap_windows_required.tsv` names the two sites the closure rests on,
each asserted retrieved and used:

| site | class | reached by |
|---|---|---|
| `5,315,591 C>CT` | `CLEAN_HET_INDEL` | the catalog claim, once it survives the repeat screen |
| `5,331,265 TAGAC>T` | `NOISY_CAND_HET` | the re-solve admitting the noisy class -- no catalog record covers it |

Losing either fails the suite even if the gap still closes by some other route.
