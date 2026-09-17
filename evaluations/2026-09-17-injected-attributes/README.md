# Do injected candidates carry all their attributes?

The injection verifier checks the candidate table as a whole. It never checked
the fields of the candidates injection **adds**, so those were tested directly:
the 58 graph-only candidates injection adds across the six panel windows, against
the invariants a candidate's own fields must satisfy.

| invariant | failures |
|---|---:|
| `DP == REF_COUNT + ALT_COUNT` | 0 / 58 |
| `AF == ALT / (REF + ALT)` | 0 / 58 |
| `DP > 0` | 0 / 58 |
| strand REF fields filled when `REF_COUNT > 0` | **39 / 58** |
| strand ALT fields filled when `ALT_COUNT > 0` | **58 / 58** |
| `LOW_QUAL_COUNT` filled at `DP > 20` | **56 / 58** |

Depth, counts and allele fraction are exactly consistent. **The strand tallies
and the low-quality depth are not populated at all.**

## Why that is not only a reporting gap

`classify_graph_only_candidates` applies the ONT strand-bias screen to graph-only
candidates, and reads exactly those fields:

```cpp
} else if (opts.is_ont() &&
           [&] {
               const int fa = c.forward_alt;
               const int ra = c.reverse_alt;
               const int expected = (fa + ra) / 2;
               if (expected <= 0) return false;
               return fisher_exact_two_tail(fa, ra, expected, expected) <
                      opts.strand_bias_pval;
           }()) {
    cat = VariantCategory::StrandBias;
```

With both tallies at zero, `expected` is zero and the lambda returns false, so
**the screen can never fire for an injected candidate** while every
alignment-derived candidate is subject to it. The screen is ONT-only, so the HiFi
defaults measured throughout this work are unaffected -- but the filter is wired
up and unable to work.

## The fix is not one function, which is why it is not made here

Filling the tallies in `backfill_graph_candidate_counts` from
`chunk.reads[prof.read_id].reverse`, by the same rule the BAM pipeline uses at
`collect_var.cpp:483`, was implemented and measured: ALT tallies went from 58
empty to 25, and **29 records then had strand sums that disagreed with their own
alt counts**. Alt coverage for these candidates accumulates at several
independent sites -- `hybrid_inject.cpp` increments `alt_cov` at three places
besides the backfill, and `graph_bam_adapter.cpp:127` assigns
`forward_ref = ref_cov; forward_alt = alt_cov`, putting every read on the forward
strand as a placeholder rather than measuring it.

Partial tallies are worse than absent ones: at zero the Fisher test no-ops, while
at partial values it runs on numbers that are wrong, and can return a
`StrandBias` verdict a correct count would not support. So the change was
reverted. A coherent fix has to populate the strand fields at every site that
contributes coverage to a graph-only candidate, including replacing the
all-forward placeholder in the graph adapter, and be validated in `--ont` mode
where the screen actually runs.

## Fixed: every site that accumulates coverage now carries the strand

The earlier attempt patched only `backfill_graph_candidate_counts` and produced
strand tallies that disagreed with their own counts, because coverage for an
injected candidate accumulates at more than one place. All of them are now
consistent, and the design rule is that **each strand tally mirrors the count it
accompanies**, so `forward + reverse == count` holds by construction wherever
coverage is added rather than being re-derived afterwards:

| site | what it counts | strand source |
|---|---|---|
| `extend_bam_profile_with_graph_obs` | graph observations applied to an existing BAM read | `chunk.reads[read_i].reverse` -- the read's own alignment |
| graph-read injection | a read present in the graph and not in the BAM | `GraphReadAllele::reverse`, the traversal orientation |
| `backfill_graph_candidate_counts` | profile entries at graph-only candidates | `chunk.reads[prof.read_id].reverse` |

Two further defects fell out of the same reading. The local `ReadObs` struct
dropped the `reverse` field `GraphReadAllele` provides, so the strand was
discarded before any of these sites could use it; and an injected `ReadRecord`
never set `reverse` at all, leaving every graph-only read on the forward strand.
The backfill also discarded low-quality observations entirely (`allele < 0`
covers both -1 uninformative and -2 low quality), so `low_qual_cov` stayed at
zero; -2 now counts as depth without counting as an allele vote.

### Measured over the 58 injected candidates in the six panel windows

| invariant | before | after |
|---|---:|---:|
| `DP == REF_COUNT + ALT_COUNT` | 0 fail | **0 fail** |
| `AF == ALT / (REF + ALT)` | 0 fail | **0 fail** |
| ALT strand tally present | **58 fail** | **0 fail** |
| REF strand tally present | **39 fail** | **0 fail** |
| ALT strand sums to `ALT_COUNT` | -- | **0 fail** |
| REF strand sums to `REF_COUNT` | -- | **0 fail** |

The tallies are measurements rather than a placeholder: across those candidates
the ALT observations split **1,324 forward against 1,509 reverse**.

Four records report a DP one greater than the number of reads covering the site
with 25 bp of flank either side; each matches the unflanked covering count
exactly and is unchanged from before the fix, so that is the check's flank
requirement, not double counting.

`--ont` now runs with real inputs on these candidates: 118 carry a non-zero ALT
strand tally in that window, and none is flagged `STRAND_BIAS`, which is the
expected outcome for HiFi reads -- the screen can now reach a verdict instead of
declining to test.

Panel, both arms: stock defaults **byte-identical** (0 concordant-to-discordant,
0 records lost or gained), `--retry-unphased-with-bam` 4 of 6 spanned at 99.49%
with 0 concordant-to-discordant. Suite 5/5. The verifier carries a seventh check,
`attributes`, so this stays gated: **0 across all six windows**.

### Still a placeholder, out of this path

`graph_bam_adapter.cpp:127` sets `forward_ref = ref_cov` and
`forward_alt = alt_cov`, putting every read on the forward strand. That is
`add_graph_candidate`, called only from the standalone graph chunk builder and
not from hybrid injection, so it is untouched here and recorded as the remaining
instance of the same pattern.
