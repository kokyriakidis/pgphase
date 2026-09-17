# Can we just use the clean sites the alignment channel finds in the gap?

Counted every candidate the alignment channel places **strictly inside** each
panel gap (`gap_left < pos < gap_right`, so neither boundary counts), by
category, against the hybrid's retry arm over the same intervals.

| category | alignment channel | phase information |
|---|---:|---|
| `CLEAN_HOM` | **266** | none by construction |
| `NOISY_CAND_HOM` | 49 | none |
| `NOISY_CAND_HET` | **24** | the only interior het evidence |
| `CLEAN_HET_SNP` | **0** | -- |
| `CLEAN_HET_INDEL` | **1** | already phased |

Per window the het counts are 4, 7, 2, 3, 2 and 6 -- all `NOISY_CAND_HET` -- and
the hybrid's composition is the same to within one or two sites, which are
`CLEAN_HOM`/`NOISY_CAND_HOM` reclassifications.

**There is exactly one clean heterozygote inside any of the six gaps:**

```
48,183,977   CLEAN_HET_INDEL   DP=75   42/33   PS=48,147,227
```

It is already in a phase set, and its window still does not span.

## Why that settles it

The clean class needs no new machinery: `assign_hap_based_on_germline_het_vars_kmeans`
with `kCandGermlineClean` is stage 1, it runs in both channels by default, and it
already consumes every clean het there is. So "use the clean sites from the
alignment region" is not an available lever -- inside these gaps the clean class
holds 266 homozygotes and one heterozygote.

That is also the reason every approach in this investigation has had to engage
the noisy class: with 24 `NOISY_CAND_HET` records against 1 clean het across
6 gaps, the interior evidence is noisy by composition, not by choice. And it is
why the phantom problem is unavoidable rather than incidental -- the same 24
records contain both the informative sites (`24,121,714`, whose two alleles split
the haplotypes perfectly) and the bridges that invert flanks (`24,131,708` at
9-against-5, `55,883,020` at 9-against-8).

A gap in this pipeline is, by this census, precisely an interval where the clean
class has nothing to say.

## Are both sources' clean sites used? Yes, and all of them participate

Counted over the whole of each panel window, not just the gap interior:

| source of the clean het | candidates | phased (`PHASE_SET != 0`) |
|---|---:|---:|
| alignment-only | 65 | **65** |
| catalog-matched | 450 | **450** |
| **total** | **515** | **515** |

Every clean heterozygote from either source carries a phase set. Nothing is
excluded by provenance -- the stage-1 mask is a property of the category, and a
candidate's category does not record where it came from.

And the union is doing real work rather than duplicating: **37 sites are
`CLEAN_HET` only because the graph claimed them.** The alignment channel calls
each of them `NOISY_CAND_HET` at the *same* depth -- 64 against 64 at
`48,243,387`, 66 against 66 at `48,243,748`, 73 against 73 at `55,895,053` -- so
the promotion comes from the catalog claim, not from better counts, and all 37
end up phased. Without the claim they would sit in the noisy class, which the
hybrid keeps out of its solve.

The qualification that matters: **all 37 are in the flanks, none inside a gap.**
That is the same census result from the other direction. The union strengthens
the blocks that already exist on either side; it adds nothing to the intervals
between them, because inside those intervals there is one clean heterozygote in
total.

## What class does the competitor's in-gap evidence fall into?

Every heterozygote hiphase phases strictly inside a panel gap, matched to our
candidate at the same locus (indel anchors differ by a base or two, so matching
is by locus):

| our category at that locus | n |
|---|---:|
| `NOISY_CAND_HET` | **14** |
| `NOISY_CAND_HOM` | **3** |
| `CLEAN_HET_INDEL` | 1 |
| `CLEAN_HET_SNP` | **0** |
| no candidate at all | **0** |

**18 sites, and not one of them is a clean heterozygous substitution.** The
single clean record is `48,183,977`, the one clean in-gap heterozygote the census
found, and it sits on a gap boundary. Discovery is not the difference either --
we hold a candidate at all 18.

So the competitor's interior evidence *is* the class this pipeline excludes. That
is the same conclusion the census reaches from the supply side: there is nothing
clean in these intervals to use, so anything that crosses them crosses on noisy
candidates.

### Three of them we call homozygous

| locus | hiphase | our category | our DP |
|---|---|---|---:|
| `55,883,019` | `A,AAT` at `2\|1` | `NOISY_CAND_HOM` | 37 |
| `5,339,363` | `GTGTGTGT...` at `2\|1` | `NOISY_CAND_HOM` | 40 |
| `12,735,894` | `TA>T` at `0\|1` | `NOISY_CAND_HOM` | **8** |

The first two are multiallelic in hiphase and genotyped `2|1` -- both haplotypes
non-reference, the case a single biallelic record cannot describe and the one the
merge work addresses. The third is the depth signature: DP 8 where the locus
carries roughly 70-fold coverage.

And `24,121,713` is the only one of the 18 still at `PHASE_SET = 0` in the retry
arm -- hiphase emits it with three ALTs at `3|1`, we hold it as one merged record
with two, and it is excluded from the solve rather than mis-measured.
