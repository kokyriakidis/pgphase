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
