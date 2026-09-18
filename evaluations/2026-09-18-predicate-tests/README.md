# Unit tests for the admission predicates, and a second case bug

`make predicate-tests` grew from 65 assertions to **105 in 10 test cases**,
covering three more functions. Two needed exposing: `var_is_homopolymer_pg` and
`var_is_repeat_region_pg` lost `static` and are declared in `collect_var.hpp`;
`allele_depths_call_het` likewise in `collect_phase.hpp`.

| function | what the cases pin |
|---|---|
| `var_is_homopolymer_pg` | deletion inside a long run; soft-masked reference; a unit-2 STR; non-repeat context rejected; indel longer than `xid` not judged; empty slice |
| `var_is_repeat_region_pg` | deletion whose motif repeats three times downstream; insertion consistent with the tract; different motif rejected; non-repeat rejected; over `xid` not judged; three copies running past the slice; position before the slice; a deletion straddling a soft-mask boundary; a run of N rejected |
| `allele_depths_call_het` | each of the six documented exclusions in order, the half-open window test `[beg, end)`, chunk-wide admission under `joint_het_orientation`, and that a homopolymer indel at a textbook fraction still passes -- which is why this predicate re-admits chance-level repeat sites to the link list |

## The second case bug, found by writing the test

`var_is_repeat_region_pg`'s insertion branch compared byte-exact strings:

```cpp
std::string ref_b = ref_seq.substr(off, len);   // raw reference, soft-masked
...
alt_b[j] = var.alt[j];                          // uppercased at bam_digar.cpp:323
return ref_b == alt_b;
```

So an insertion inside a lowercase tandem repeat never matched. Its sibling
`var_is_homopolymer_pg` reads the reference through `nt4_from_ref_char`, which
uppercases, so the two complementary predicates disagreed on the same locus for
a reason unrelated to repeat structure. The deletion branch used `memcmp`
between two reference windows, which is case-consistent unless the windows
straddle a mask boundary. Both branches now compare through `nt4`, and an
ambiguous base is no longer a match, so a run of N cannot pass as a tandem
repeat.

## It is latent: chr20 output is byte-identical

| chr20 insertions | |
|---|---:|
| INS candidates | 11,258 |
| sampled | 400 |
| lowercase (soft-masked) context | 278 (69.5%) |
| `var_is_repeat_region_pg` verdict changed by the fix | 196 (49%) |
| of those, where `var_is_homopolymer_pg` does not already fire | **0** |

The classification at `collect_var.cpp:1490-1492` is
`var_is_homopolymer_pg(...) || var_is_repeat_region_pg(...)`, and the first is
itself a unit-1..6 three-copy STR test read case-insensitively. It fires at
every one of the 196 loci, so the OR never changed and the chr20 VCFs from the
two arms are byte-identical (116,292 records each). The sampling check
reimplements the sibling forward-only while the real one also scans backward, so
it can only fire more often -- the zero holds a fortiori.

The fix is therefore correctness hygiene, not a behaviour change, and the
`repeat_region` insertion branch is currently redundant on this data. A test case
pins that redundancy explicitly, so narrowing `var_is_homopolymer_pg` later
cannot silently expose the old bug.

## A process note

The first A/B was invalid and said so out loud: reverting the whole
`collect_var.cpp` to HEAD also reverted the `static` removal the header now
depends on, the base build failed with 2 errors, and `cp pgphase` copied the
still-fixed binary -- both arms ran identical code, caught only by an explicit
`cmp` guard printing "WARNING identical". The valid base arm reverts the two
logic blocks by string replacement and keeps the signatures.

Suites: unit 4/4, window 66/66, predicate 105/105.
