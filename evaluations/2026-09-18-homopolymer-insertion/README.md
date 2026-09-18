# No insertion was ever flagged as a homopolymer indel

`var_is_homopolymer_indel` (`collect_phase_noisy.cpp:163`) compared the
reference as a **raw FASTA byte** against an **nt4-coded** alt base:

```cpp
const uint8_t ins_base0 = base_to_nt4(alt[0]);        // 0..3
...
if ((uint8_t)(unsigned char)chunk.ref_seq[idx] != ins_base0) return false;   // 'a' == 97
```

`'A'` is 65 and `'a'` is 97; neither can equal 0..3, so the insertion branch
returned false at `i == 0` for every input. The deletion branch compares raw
bytes on both sides and worked. The result: the pipeline flagged homopolymer
**deletions** and silently never flagged homopolymer **insertions**.

longcallD carries the same mismatch -- `collect_var.c:1730` tests
`chunk->ref_seq`, a `char*` from `faidx_fetch_seq`, against an nt4 abPOA
consensus base -- so the old comment claiming parity ("same as LCD") was
accurate and this fix is a knowing divergence from upstream. Upstream's other
detector, `var_is_homopolymer` (`collect_var.c:306`), is a proper STR test and
is already ported correctly here as `var_is_homopolymer_pg`
(`collect_var.cpp:1169`); it is used for CLASSIFICATION, not for this flag.

Soft-masking compounded it: these loci read `gtctcaaaaaaaaa` and
`ttgatttttttttt` in CHM13, so a raw-byte comparison would have failed on case
even with matching encodings. Both branches now compare through `base_to_nt4`,
which is case-insensitive, and both gained bounds guards.

## Why it matters: the flag gates four consumers

| consumer | site | effect of a missed flag |
|---|---|---|
| `init_assign_read_hap` | `collect_phase.cpp:377` | the site scores reads it should not |
| `iter_update_var_hap_cons_phase_set` | `:551` | it enters the link list |
| `update_read_phase_set` | `:868` | it confers phase-set eligibility |
| `select_init_var` | `:198` | it can be chosen as the sweep pivot |

## The seam it fixes

`chr20:55,843,827-55,889,113`, where hiphase spans at 99.83% and we inverted a
flank. The chain crosses a 20.8 kb step through three interior noisy calls that
hiphase declines (`./.`, `0/1` unphased, `./.`). Two segregate at chance against
read truth -- 0.525 at `55,871,837` and 0.571 at `55,882,617` -- and the first
is an insertion in an 8 bp A-run, so it was flagged 0 and scored reads.

| arm | blocks | spans | in-gap hets | reads | accuracy | flanks |
|---|---:|---|---:|---:|---:|---|
| before | 2 | yes | 6 | 409 | 56.48% | L=PAT 0.99 R=MAT 1.00 **switched** |
| after | 2 | yes | 5 | 377 | **99.73%** | L=PAT 0.99 R=PAT 1.00 consistent |

## Whole chr20, `collect-bam-variation` default

| arm | tagged | blocks | phased hets | discordant | read hamming |
|---|---:|---:|---:|---:|---:|
| base | 222,255 | 358 | 81,050 | 8,146 | 3.665% |
| fixed | 216,955 | 365 | 80,926 | 1,427 | **0.658%** |

-6,719 discordant reads, a 5.6x error reduction, for 2.4% fewer tagged reads
and 7 more blocks. Both binaries were built from the same tree with only this
file differing, so the comparison isolates the fix.

Scope of the leak, measured in the 105 kb seam window: 6 of 12 noisy MSA
insertions sit in a reference run of >= 5 and were all flagged 0 before, against
12 of 19 deletions correctly flagged. So it was systematic, not one site.

Unit tests 4/4, window tests 66/66, both unchanged.

## The test that would have caught it

`make predicate-tests` -> `src/test_phase_predicates.cpp`, Catch2 against the
vendored v2.13.10 header, 65 assertions in 7 test cases, ~2 s. It needs no BAM
and no reference file: each case builds a `PhasingChunk` carrying only
`ref_beg`/`ref_end`/`ref_seq`, which is the entire input this predicate reads.

Coverage:

| function | cases |
|---|---|
| `base_to_nt4` | both cases of ACGT, `U`, and that `N`/`n`/`-` land outside 0..3 so callers can reject them |
| `var_is_homopolymer_indel`, insertions | single and multi-base into a run; mixed alt rejected; alt not matching the run rejected; run shorter than 5 rejected; no run rejected; empty alt rejected |
| `var_is_homopolymer_indel`, deletions | single and multi-base inside a run; span crossing a base boundary rejected; soft-masked run accepted |
| `var_is_homopolymer_indel`, guards | SNP, position before the slice, context running off the end, ambiguous base in context, empty slice |
| `select_stitch_orientation` | all four rules: net-margin (clear/tie/margin/null opts), both-strands (one empty link, both qualifying, equal abstain), literal, both-strands-margin |

Two cases are regressions for the bug above: a lowercase run, and the literal
chr20:55,871,832 context `gtctcaaaaaaaaa` with the insertion at 55,871,837.

Fault injection: restoring the raw-byte comparison fails **6 assertions in 1
test case** -- lines 63, 66, 67, 93, 94 and 103, i.e. every insertion case
including both regressions -- while the other 59 still pass. Restored by file
copy rather than `git checkout`, and re-run through `make predicate-tests` so
the binary actually relinks: `make -j20` alone rebuilds the object but leaves a
stale test binary, which made a first injection attempt report a false pass.

Exposing the two predicates required moving them out of the anonymous namespace
in `collect_phase_noisy.cpp` to namespace scope, declared in
`collect_phase_noisy.hpp`. They have no other callers outside that file.
