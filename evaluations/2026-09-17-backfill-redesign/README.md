# Backfill, rewritten: derive the counts, do not accumulate them

## What was wrong with the shape, not just the fields

Coverage for a graph-only candidate was accumulated **incrementally at every site
that touched a read profile**, and each of those sites re-implemented the count
update:

| site | what it counted |
|---|---|
| `backfill_graph_candidate_counts` | BAM profile slots, swept before Phase B |
| `extend_bam_profile_with_graph_obs` | graph observations applied to a doubly-mapped read |
| the graph-read injection loop | observations on a read present only in the graph |

Three consequences followed from the shape itself, and this session hit all
three:

- **A field added to the counts has to be added in three places.** That is why
  the strand tallies were missing at all of them while depth and allele counts
  were maintained, and why patching one site produced tallies that disagreed
  with their own counts.
- **Whether a read is counted twice is a question answered by inspection**, and
  it has to be re-answered whenever a site is added, because a slot written by
  one site is visible to the sweep in another. It was answered wrongly: see
  below.
- **The sweep was not idempotent.** It accumulated onto whatever was already
  there, so running it twice doubled every field.

## The rewrite

Profile mutation writes alleles. **Nothing else writes counts.** The counts are
derived from the final profile state by a single sweep that zeroes first, so they
are a pure function of the profiles -- idempotent, and double counting is not
expressible, because a profile holds one slot per candidate and no other path
adds coverage.

One primitive is the only place an observation becomes coverage:

```cpp
static void accumulate_observation(VariantCounts& counts, int allele, bool reverse) {
    if (allele == kProfileAlleleLowQual) { ++counts.low_qual_cov; return; }
    if (allele < 0) return;                       // uninformative
    ++counts.total_cov;
    if (allele == 0) { ++counts.ref_cov; reverse ? ++counts.reverse_ref : ++counts.forward_ref; }
    else             { ++counts.alt_cov; reverse ? ++counts.reverse_alt : ++counts.forward_alt; }
}
```

The ordering had to change with it. The sweep ran at `collect_pipeline.cpp:988`,
**before** Phase B injected graph reads at 990 -- which is precisely why the
injection sites needed counts of their own. It now runs after Phase B and before
`classify_graph_only_candidates`, where the existing comment already claimed
counts were final.

## What the rewrite found: a real double count

The only field that moves is `LOW_QUAL_COUNT`, from non-zero to zero on 8 of the
58 injected candidates -- and it is a double count being removed, not information
being lost. `extend_bam_profile_with_graph_obs` overwrites a slot whose previous
allele is negative, and `-2` (low quality) is negative: so the pre-Phase-B sweep
counted that read as low-quality depth, Phase B replaced its slot with the
graph's allele, and the injection site counted the same read again as an allele
vote.

Measured against the number of reads covering each site:

| locus | covering reads | old `DP` | old `DP + LOW_QUAL` | new `DP` | new `LOW_QUAL` |
|---|---:|---:|---:|---:|---:|
| 5,290,557 `SNP` | 81 | 80 | **86** | 80 | 0 |
| 12,795,753 `DEL` | 68 | 66 | **70** | 66 | 0 |
| 12,790,908 `INS` | 68 | 68 | **71** | 68 | 0 |

`depth_with_low_quality = total_cov + low_qual_cov` is what
`classify_variant_initial` tests against `min_depth`, so those candidates were
being admitted on a depth that exceeded their own coverage.

## Evaluation

| check | result |
|---|---|
| unit suite | 5/5, with a new test for the contract |
| panel, stock defaults | byte-identical: 0 concordant-to-discordant, 0 records lost or gained |
| panel, `--retry-unphased-with-bam` | 4 of 6 spanned, 99.49%, 0 concordant-to-discordant |
| injection verifier | 1,575 candidates: dropped 0, duplicated 0, missing 0, attributes 0 |
| `DP == REF_COUNT + ALT_COUNT` | 0 / 58 fail |
| `AF == ALT / (REF + ALT)` | 0 / 58 fail |
| strand tallies sum to their counts | 0 / 58 fail |
| strand tally present when count > 0 | 0 / 58 fail |

The new unit test asserts the contract rather than a value: the sweep derives
depth, allele, strand and low-quality counts from three profiles carrying `0`,
`1` and `-2`; sweeping twice gives the same counts; and the strand tallies sum to
their counts which sum to depth.

## Cleaning the stale parts of the hybrid path

Done on evidence rather than intuition -- the compiler and the call graph decide
what is dead, and every cited reference was checked rather than assumed.

**Removed or corrected:**

| item | why it was stale |
|---|---|
| `hybrid_inject.hpp` contract for `backfill_graph_candidate_counts` | said *"Call AFTER collect_var_build_profiles, BEFORE inject_graph_reads"* -- the opposite of what the code now requires, and the most dangerous kind of stale comment, since it would lead the next change straight back into the bug |
| its description "scans existing BAM profiles and accumulates the missing counts" | it no longer accumulates and no longer reads BAM profiles alone; it zeroes and derives from every profile |
| `collect_pipeline.cpp:982` comment | described the backfill but was left sitting above `clear_bam_evidence_at_graph_candidates` when the call moved |
| `redundant` counter in `augment_chunk_with_graph_sites` | set but never read (`-Wunused-but-set-variable`); the behaviour it counted is covered by the verifier's candidate totals |
| two doc comments saying support is "accumulated" from BAM and graph reads | the mechanism is now a derivation |

**Checked and left alone, with the evidence:**

- `-Wall -Wextra -Wunused` on both hybrid translation units: **0 warnings** after the
  counter went, so there is no other dead code or unused parameter there.
- Every function declared in `hybrid_inject.hpp` has callers.
  `clear_bam_evidence_at_graph_candidates` looks unreachable but is not -- it is
  the `graph_authoritative` branch, and that option does have a CLI setter
  (`hybrid_collect.cpp:370`).
- `graph_only_pre_sort`, `all_graph_pre_sort` and `pre_sort_vcf_alleles` are all
  read; so is every field of the local `ReadObs`, including the `reverse` field
  added with the strand fix.
- All three `file.cpp:NNN` references cited in the hybrid sources still point at
  the code they claim: the two allele-fraction expressions and the homopolymer
  flag assignment.

Behaviour-neutral, as a cleanup should be: suite 5/5, and both arms emit
**identical records** to the pre-cleanup build -- 0 lost and 0 gained on stock
defaults and on `--retry-unphased-with-bam`.
