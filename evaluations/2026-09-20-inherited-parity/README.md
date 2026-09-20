# Should a site that only inherits a phase set inherit its parity? Measured: no

## The hypothesis

`iter_update_var_hap_cons_phase_set` (`collect_phase.cpp:670-692`) resolves
orientation for link-list members only:

```
if (hi >= 0) { phase_set = het_ps[hi]; if (parity[hi]) swap(cons[1], cons[2]); }
var.phase_set = phase_set;                  // every candidate, list member or not
```

So a candidate kept off the link list -- a homopolymer indel, for instance --
joins the block with whatever orientation its own consensus produced. longcallD
has exactly the same shape: `assign_hap.c:409` applies its running `flip` inside
the `is_het` guard while `:418` assigns the phase set to every candidate. Two
comments in our source asserted that this asymmetry is what put a
maternal-on-hap1 site inside a paternal-on-hap1 block at `chr20:48,225,786`.

The obvious fix is to carry the running parity onto inherited sites too.

## The measurement

Implemented exactly that: track the parity of the block the walk is inside and
apply it to any candidate whose `hap_to_cons_alle[1] != hap_to_cons_alle[2]`
that is not on the link list. Whole chr20, alignment arm:

| | before | after |
|---|---:|---:|
| phased VCF records | 116,896 | 116,895 |
| shared phased records | | 81,277 |
| genotypes changed | | **128** (109 SNPs, 19 indels) |
| phase sets changed | | 136 |

Then each changed site was judged against **its own block's orientation**: take
the site's parent-to-hap1 mapping from read truth (>= 5 alternate reads, >= 90%
pure), take the block's mapping from up to 12 unchanged single-base peers within
300 kb in the same phase set, and require at least 3 concordant peer votes.

| verdict | sites |
|---|---:|
| now AGREES with its block (fix) | 18 |
| now DISAGREES with its block (regression) | **31** |
| unscorable or verdict unchanged | 79 |

Net negative, so the change was reverted; the tree is byte-identical to the
pre-change output over whole chr20 and all three suites pass (unit 3/3,
predicate 151, window 125).

## Why it fails, and what it corrects

The solve is **iterative**. The caller re-derives read labels from the flipped
consensus (`iter_update_var_hap_to_cons_alle`) and calls this function again
until nothing changes. By convergence an inherited site's consensus is already
expressed in the block's gauge, so applying the block parity to it a second time
inverts a site that was already right -- which is what the 31 regressions are,
and the examples show precisely that shape (`hap1 PAT -> MAT` where the block is
`PAT`).

So the two source comments were wrong about the mechanism, not just incomplete,
and both are now corrected in place. The switch at `chr20:48,225,786` is real --
it is what motivated `allele_depths_call_het` -- but "parity is not applied to
inherited sites" is not its cause, and the cause is open again.

Caveat on strength: 79 of the 128 changed sites could not be judged (no usable
truth at the site, or fewer than 3 concordant peers), and the peer votes that
did decide were thin (3 votes in the examples shown). The direction is
consistent and there is no evidence of benefit, which is enough to reject a
change that moves 128 genotypes, but it is not a precise error rate.
