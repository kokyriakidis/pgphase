# One arm, iterated

`arm.sh` is the only configuration. It is the alignment machinery run over the
union of the alignment's own candidates and the graph's catalog sites, one solve
per chunk, with nothing added on top: no `--recover-gaps` (the second pass zeroes
every non-graph candidate's category before phasing and skips the noisy-region
MSA), no retry, no margin filter (`min_read_hap_margin` defaults to 0 and the
hybrid subcommand does not override it), `--keep-noisy-kmeans` because
`skip_noisy_kmeans = true` is a hybrid-specific override the alignment pipeline
does not apply, and `-q 1`.

```sh
OUT=/tmp/arm REGION='CHM13#0#chr20:48126830-48279446' \
  ./evaluations/2026-09-16-single-arm/arm.sh
```

## Baseline on chr20:48,176,830-48,229,446

2 het blocks, no span, 6 in-gap phased hets, 393 reads tagged, **100.00%**
concordant against read truth. Blocks `48,149,548-48,229,226` (8 sites) and
`48,229,446-48,279,445` (49 sites).

Both site defects diagnosed under the recovery arm are **absent** here:
`48,177,780` is emitted `0|1` (the recovery arm emitted `NOISY_CAND_HOM 1|1`),
and the two sites carrying no haplotype information are not emitted at all.

## Iteration 1: the crossing site was genotyped homozygous

The block boundary sits where the chain has nothing to link across. Inside the
41.8 kb interval `48,183,976-48,225,787` the arm holds:

| source | in the interval |
|---|---|
| graph catalog | **500 sites** |
| arm candidates | 59 `CLEAN_HOM`, 1 `CLEAN_HET_INDEL`, 2 `NOISY_CAND_HET`, 6 `NOISY_CAND_HOM` |
| arm candidate at `48,204,384` | DEL, DP 71, **30 ref / 41 alt, AF 0.5774**, `NOISY_CAND_HET` |
| arm emitted record at `48,204,383` | `AT>A` **`GT=1|1`** |
| hiphase record at `48,204,383` | `AT>A` **`GT=0|1`** |

So the signal is not missing from discovery. The site hiphase crosses on is
found, classified a noisy het, and measured at a textbook heterozygous depth --
and then emitted homozygous, which is why the chain has no phased heterozygote
inside the interval and links across it with no spanning read.

`iter_update_var_hap_to_cons_alle` recomputes each haplotype's consensus
**independently by majority**, so both haplotypes select the deeper allele. The
multi-allele branch beside it already guards this case, with the comment
*"Independent haplotype majorities can select the same allele twice."* The
predicate that applies the joint orientation existed but was scoped to retry
windows; `--joint-het-orientation` applies it wherever a biallelic candidate's
own allele depths call it heterozygous, under the same conditions the predicate
already tested (category het, both haplotype profiles present, both depths above
`min_alt_depth`, allele fraction within `[min_af, max_af]`).

| | baseline | `--joint-het-orientation` |
|---|---|---|
| `48,204,383` | `1\|1` | **`0\|1`** (hiphase agrees) |
| `48,202,056` | absent | `0\|1` |
| in-gap phased hets | 6 | **9** |
| blocks / spans | 2 / no | 2 / no |
| reads tagged | 393 | 393 |
| concordant | **100.00%** | **100.00%** |

Per-site orientation shows the interval is now crossed **correctly**: the chain
runs `48,183,976` (1.000) to `48,204,383` (0.958), both PATERNAL-on-hap1 like the
rest of the body, over links of 18.1 kb / 7 reads, 2.3 kb / 58 reads and 21.4 kb
/ 7 reads. Off by default; measured on this window only.

## Iteration 2, the next target

One switch remains, and it is now isolated: `48,204,383` is PATERNAL at 0.958
and `48,225,786` (`CAAAA>C`, consistency **1.000**) is MATERNAL, so the
orientation flips across the 21.4 kb link carried by 7 reads. Everything from
`48,225,786` onward, including the right block, is MATERNAL and mutually
consistent.

Separately, the 220 bp boundary `48,229,226 -> 48,229,446` reports chain link
evidence `agree=0 conflict=0` despite **60** reads spanning it, because both of
the left block's terminal sites are homopolymer indels and those are excluded
from read scoring. That exclusion must not be lifted before the switch above is
fixed: doing so earlier produced a spanning block at 59.29% accuracy.

## Recheck: do we hold all of hiphase's signal, and encode it usably?

Every phased heterozygote hiphase emits in `48,145,000-48,240,000`, matched to
the arm's candidate and emitted record at the same locus (indel anchors differ by
a base or two, so matching is by locus, not by exact position), with each site
genotyped from the alignment and scored against read truth.

**Hiphase emits 19, all in a single phase set `48149548` spanning
`48,149,548-48,235,309`. The arm holds 15 of them and splits the window into
two blocks.**

| verdict | n |
|---|---:|
| used, alleles identical | 12 |
| used, re-represented | 3 |
| candidate only, screened out | 2 |
| absent from our candidates | 2 |

The four we do not use all carry real signal:

| site | hiphase allele | our status | n | segregation |
|---|---|---|---:|---:|
| 48,149,567 | `A>AGAG` | **present, placed 19 bp away** (see correction) | 64 | **0.984** |
| 48,173,317 | `TGGG>T` | **absent** | 66 | **1.000** |
| 48,177,725 | `T>TA` | candidate, `REP_HET_INDEL` | 74 | 0.865 |
| 48,234,100 | `C>CCT` | candidate, `REP_HET_INDEL` | 54 | **0.981** |

So two informative sites are never discovered, and two more are discovered and
then screened out by repeat demotion while hiphase phases both.

The three re-represented sites are encoded differently but remain usable:
`48,177,780` is multi-allelic for hiphase (`GAGA>G,GA`, `2|1`) and biallelic for
us (`0|1`, segregation 0.946), so one allele is dropped; `48,225,788`
(`AAAA>A`, `1|0`) is our nested-deletion split into a common `CA>C` emitted
`1|1`, which carries nothing, plus a residual `CAAAA>C` emitted `0|1` at
segregation **1.000**; `48,149,548` differs in anchor and allele only.

### This does not explain the switch

**None of the four unused sites lies inside the 21.4 kb link where our
orientation flips, and hiphase has no heterozygote there either.** It crosses
`48,204,383 -> 48,225,786` on the same 7 spanning reads we have, between two
sites we also have and which are both informative in our own data (0.958 and
**1.000**), and it gets the orientation right where we do not.

So the two findings are separate, and in this order:

1. **The switch is a decision, not a data deficit** -- same sites, same reads,
   different outcome. That is iteration 2.
2. **We do discard usable signal**, but recovering it lengthens and strengthens
   the chain rather than fixing the switch: two sites to discover
   (`48,149,567`, `48,173,317`) and two the repeat screen removes at 0.865 and
   0.981 (`48,177,725`, `48,234,100`).

## Iteration 2: a site excluded from the link list still inherited a phase set

The switch was not a decision made on thin evidence -- it was no decision at
all. `48,225,786` (`CAAAA>C`, in an A run) is `is_homopolymer_indel`, and the
het **link** list excludes such a site unless a recovery homopolymer window is
active:

```cpp
if (... && (var.msa_insertion_alts.empty() || var.gap_link_supported) &&
    (!var.is_homopolymer_indel || gap_hp_link) && !unsupported_gap_indel) {
    is_het[_vi] = true;
    het_var_idx.push_back(_vi);
}
```

The arm has no such window, so the site never entered the list and
`het_rank[_vi]` stayed `-1`. But the emit loop still hands it the running phase
set:

```cpp
const int hi = het_rank[_vi];
if (hi >= 0) {
    phase_set = het_ps[hi];
    if (parity[hi] == 1) std::swap(var.hap_to_cons_alle[1], var.hap_to_cons_alle[2]);
}
var.phase_set = phase_set;      // reached with hi < 0 too
```

So it joined the block carrying whatever orientation its own consensus produced,
with `parity` never applied -- never reconciled against the block it was labelled
into. That is why the switch was invisible in read space: the reads covering the
site belong to the next block, so per-block read concordance stayed at 100.00%
while the site itself sat on the opposite haplotype.

The fix admits a site to the link list when its own allele depths call it a clear
heterozygote (`allele_depths_call_het`, the predicate iteration 1 generalised),
so its orientation is decided by spanning reads like any other het instead of
inherited.

| | before | after |
|---|---|---|
| switches within a block | **1** | **0** |
| blocks over the window | 2 | 3 |
| reads tagged | 393 | 393 |
| concordant | **100.00%** | **100.00%** |
| `48,229,446` onward | own block `PS=48229446` | **merged into `PS=48225786`** |

Blocks after the fix: `48,147,225-48,202,056` (8 sites), `48,204,383` alone, and
`48,225,786-48,279,445` (52 sites). The 220 bp boundary that no merge could close
is now joined by the ordinary chain, because the site on its left finally has a
link.

**Correction.** An earlier version of this paragraph said the first block was
"8 sites, all paternal-on-hap1 at 0.946-1.000" and the second "all
maternal-on-hap1 at 0.933-1.000". Both ranges were taken from the reliably
scored sites only, and stated as if they covered every site in the block. The
full picture, from the same table:

| block | sites scoring >= 0.90 | sites below that floor |
|---|---|---|
| `48147225` | 5, **all paternal-on-hap1**: 48,149,548 (0.952), 48,162,480 (1.000), 48,176,830 (0.986), 48,177,780 (0.946), 48,183,976 (1.000) | 3: 48,147,225 (0.579), 48,173,989 (0.642), 48,202,056 (0.507) -- all nominally maternal |
| `48225786` | 48,225,786 `CAAAA>C` (1.000), 48,229,226 (0.933), 48,229,446 (1.000), 48,230,918 (0.983), 48,232,579 (1.000), 48,232,790 (1.000), 48,234,732 (1.000), 48,235,242 (1.000), **all maternal-on-hap1** | 2: 48,225,786 `CA>C` (0.733), 48,235,309 (0.610) -- also maternal |

So the true per-block ranges are 0.507-1.000 and 0.610-1.000, and the defensible
claim is narrower than the one made: **every site that carries usable haplotype
signal agrees on its block's orientation, and no site above the floor
contradicts it.** The three sub-floor sites in the first block reading the
opposite way is what an uninformative site does -- 0.507 is chance -- not
evidence of a switch.

## Iteration 3, the next target

The window now fragments where it used to switch. `48,204,383` is a singleton:
its links to `48,202,056` (2.3 kb, 58 spanning reads) and to `48,225,786`
(21.4 kb, 7 spanning reads) are both refused. `48,202,056` -- the site carrying
no haplotype information, segregation 0.507 -- sits between it and the left
block, so the chain's only route from the body to `48,204,383` runs through a
site whose alleles are noise. Removing that site, or refusing to link through it,
is the next step; hiphase holds no het there at all.

Measured on this window only; `--joint-het-orientation` remains off by default.

## Iteration 3: the arm was running the block-link vote in the wrong mode

`--link-by-alleles` was omitted from the arm on the stated grounds that it is
already default-on. That was wrong, and it is worth stating precisely because
the two options read almost identically:

| field | default | set by |
|---|---|---|
| `gap_link_by_alleles` (`phasing_types.hpp:310`) | **true** | cleared by `--no-gap-link-by-alleles` |
| `link_by_alleles` (`phasing_types.hpp:430`) | **false** | set by `--link-by-alleles` |

The block-link vote reads the **second**:

```cpp
const int agree = opts.link_by_alleles
        ? check_agree_alleles(chunk, read_i, vj, vi)
        : check_agree_haps(chunk, read_i, chunk.haps[read_i], vj, vi);
```

With it off, `check_agree_haps` requires the read to already carry a haplotype
**and** its allele at the left site to match that haplotype's own consensus, so
at a noisy site nearly every read returns `-1`. Probed at the two links in
question: `48,202,056 -> 48,204,383` has 82 overlapping reads, **57 with a usable
allele at both sites**, and scored **agree=1 conflict=0**. The allele histograms
hold only `-1/0/1`, so nothing exotic is being discarded -- the vote was simply
asking the wrong question.

| arm | blocks | spans | tagged | concordant |
|---|---:|---|---:|---|
| without `--link-by-alleles` | 3 | no | 393 | **100.00%** |
| with `--link-by-alleles` | **2** | no | 393 | **100.00%** |

The singleton rejoins the left block, which becomes `48,147,225-48,204,383`
(9 sites), at no cost in accuracy. `arm.sh` now passes the flag, with the
distinction between the two options written next to it.

## Iteration 4, the next target

One split remains, `48,204,383 -> 48,225,786`, and the probe says why it cannot
be linked: of 173 reads overlapping the pair, **exactly 1 has a usable allele at
both sites** -- 102 lack one at the left, 70 at the right. Positionally 7 reads
span both. So the link is starved of read-level observations, not of reads.

At `48,225,786` the allele is `-1` for 29 to 52 of the overlapping reads
depending on which of the two split records is used. Hiphase crosses the same
21.4 kb using its own `48,225,788 AAAA>A`, so the difference is that it obtains
an allele call for reads where we record none. That is the next thing to fix:
why the per-read allele is unset at that site.

## Correction: one of the two "absent" sites is a matching artifact

`48,149,567` was reported absent because nothing matched within 6 bp. That window
is too tight for a tandem repeat, where the aligner places the same event
anywhere within its tract. The arm does hold this variation, as **two records
19 bp upstream at `48,149,548`**:

```
48,149,548  A>AAGAAGAGAAG    GT=0|1      (+10 insertion)
48,149,548  AAGAAG>A         GT=1|0      (-5 deletion)
```

and that is exactly what the reads carry. Net length change over a 25 bp window
at `48,149,567`, split by truth haplotype:

| haplotype | net length |
|---|---|
| PATERNAL | **+10 (29 reads)**, +5..+13 scatter, 0 (2) |
| MATERNAL | **-5 (10 reads)**, 0 (11), -6 (3), -10 (1) |

So the locus is discovered, genotyped and phased, in the opposite orientation on
the two records, matching the reads. Hiphase writes the same event as
`48,149,548 AAGA>A` plus `48,149,567 A>AGAGA`; we write it as one position with
two records. **Matching competitor sites must be by event within the repeat
tract, not by position within a few bases** -- the same rule the project already
recorded for anchor offsets, applied at repeat scale rather than 1-2 bp.

The corrected accounting for this window is therefore **1 absent, not 2**.

## `48,173,317` is a genuine miss, and both channels fail on it

```
net length change over +/-25 bp, by truth haplotype:
  PATERNAL :  -7 (39 reads),  -8 (3),  -9 (1)
  MATERNAL :   0 (21 reads),  -1 (2)
reference context: ggttccttggggatgggggattgggggatggga     (GGGGATG tandem repeat)
```

A clean 7 bp deletion, 43 reads against 23 at roughly 66x, which hiphase calls
`TGGGGATG>T`. The arm has **no candidate and no emitted record within 120 bp**.

**Retraction.** An earlier version of this section said `collect-bam-variation`
alone does not call it either, and concluded that both channels fail. That was
never tested -- the check queried the *hybrid* run's outputs. Run standalone, the
alignment channel discovers **and phases** it:

```
POS=48173318  DEL  GGGGATG>.  DP=64  23/41  AF=0.6406  NOISY_CAND_HET
record: 48173317  TGGGGATG>T  GT=0|1
```

So only the hybrid path loses it, and the mechanism is the catalog claiming the
locus.

### Localized: the catalog claim routes it to a stricter classifier

Probing the injection path, nothing fails there -- the site is eligible, present
in the chunk's catalog view, and a graph-only candidate **is** created (index
586, key position 48,173,317). It is then dropped downstream. The counts after
`backfill_graph_candidate_counts` explain why:

```
pos=48173317  graph_site=1  ref_cov=23  alt_cov=42  total=65  af=0.000  ref_len=7
```

The coverage is right and matches the reads, and `allele_fraction` is **0.000** --
the backfill accumulates the three coverage fields and never publishes the
fraction they imply. That is fixed here (same expression as
`gap_evidence.cpp:256` and `graph_bam_adapter.cpp:693`), though it is not the
whole story, because `classify_graph_only_candidates` recomputes the fraction
itself.

The decisive difference is *which classifier judges the site*. Because the
catalog claims the locus, it goes through `classify_graph_only_candidates`,
whose het-indel test is `|AF - 0.5| <= graph_indel_af_margin`, i.e. AF in
**[0.39, 0.61]** at the default 0.11. The site's AF is **0.652**, so it falls to
the `else` branch -- and that branch assigns `LowCoverage`, which
`prune_not_candidate_variants` **deletes**. The identical site without a catalog
claim is judged by the BAM classifier's `[min_af, max_af]` = **[0.20, 0.80]**
and comes out `NOISY_CAND_HET`, phased `0|1`. So a catalog-claimed site is
judged more harshly than the same site would be without the catalog, and the
penalty is deletion rather than demotion -- while the branch's own comment says
its purpose is only to "keep out of k-means".

### Why simply admitting it is not the fix: the verification asymmetry

Changing that branch to `NoisyCandHet` does recover the site -- hiphase sites
used goes 15/19 to 16/19 and both blocks grow -- but read concordance on the arm
falls from **100.00% (393/393)** to **88.39% (449/508)**, with *both* blocks
turning internally inconsistent (87.06% and 90.45%). So it is not one flip.

The cause is that a catalog-claimed site cannot be verified. `msa_verified` is
set in exactly **one** place, inside the factory that *constructs* a candidate
from the MSA consensus (`collect_phase_noisy.cpp:213`). It is a property of
MSA-created candidates, not a stamp applied to existing ones -- and a
catalog-claimed candidate exists *before* that pass. Probed across
`48,145,000-48,240,000`:

| | sites | `msa_verified` |
|---|---:|---|
| BAM-discovered `NoisyCandHet` | 14 | **1 for all of them** |
| catalog-claimed `NoisyCandHet` | 5 | **0 for all of them** |

So admitting them puts *unverified* sites into noisy k-means -- three of the five
at AF 0.725-0.791 -- which is what corrupts both blocks. The BAM-discovered
sites in the same window are safe precisely because they arrive verified, having
been constructed by the MSA from the reads.

The fix is therefore to let the noisy MSA pass construct the verified version of
a catalog-claimed site, exactly as it already does for the identical site when
the catalog does not claim it. Admitting it raw is measurably wrong, and hiding
the admission behind a flag would only make the site unreachable by default. The
branch now carries this reasoning so the drop is a named decision rather than an
accident.

## The fix: an MSA call must not lose to a candidate that is about to be deleted

The noisy MSA pass was never the problem. Both the hybrid and the standalone
alignment runs build the **same 41 noisy regions**, and
`48,173,237-48,173,436` covers the site in both -- so the MSA constructs its
verified call at `48,173,317` in the hybrid run too. It is then thrown away.

`merge_var_profile` decides what happens when the MSA's new variant has the same
key as an existing candidate:

```cpp
const bool replace_repeat = (admit_all_in_region || whitelisted) &&
    old_vars[old_i].counts.category == VariantCategory::RepeatHetIndel &&
    admissible_type(new_vars[new_i]);
const bool replace_selected = replace_sites != nullptr && ...;
if (replace_repeat || replace_selected) { /* take the MSA variant */ }
else { merged_vars.push_back(old_vars[old_i]); }   // keep the old one
```

The MSA's call is preferred only when the old candidate is specifically
`RepeatHetIndel`, or explicitly listed. Our catalog-claimed candidate is
`LowCoverage` -- the category `prune_not_candidate_variants` **deletes** -- so
the pipeline kept a candidate it was about to throw away, in preference to the
verified call that would have survived. The site then vanished at prune time.

The fix adds one condition: an old candidate in a pruned category loses to the
MSA's call at the same key.

```cpp
const bool replace_pruned =
    (old_vars[old_i].counts.category == VariantCategory::LowCoverage ||
     old_vars[old_i].counts.category == VariantCategory::LowAlleleFraction) &&
    admissible_type(new_vars[new_i]);
```

This is safe by construction rather than by threshold: the alternative to
replacing a pruned candidate is **no site at all**, and what replaces it is the
MSA-constructed version, which carries `msa_verified` for the same reason the 14
BAM-discovered sites in this window do.

| arm | blocks | sites used (+/-6 bp) | sites used (event) | reads tagged | concordance |
|---|---:|---:|---:|---:|---|
| before, site pruned | 2 | 15/19 | 16/19 | 393 | **100.00%** |
| raw admission (`NoisyCandHet` in the classifier) | 2 | 16/19 | 17/19 | 508 | 88.39% |
| **MSA replaces the pruned candidate** | 2 | 16/19 | **17/19** | 393 | **100.00%** |

Two columns because the naive matcher is the one this document already
retracted: counting a hiphase site as used only when one of ours sits within
6 bp marks `48,149,567` unused, although the arm emits that event 19 bp away at
`48,149,548`. **The event column is the correct criterion** -- a tandem-repeat
indel may be placed anywhere in its tract -- and the narrow column is kept only
because the earlier tables in this file quote it.

The site comes back with the alignment channel's own numbers -- `48,173,318 DEL
GGGGATG>. DP 64, 23/41, AF 0.6406, NOISY_CAND_HET`, emitted `TGGGGATG>T` -- and
read tagging is untouched: 393 reads, both blocks still 100.00% internally
consistent. So the recovery is purely site-level, with no read-level cost, which
is what admitting an unverified site could not achieve.

`classify_graph_only_candidates` is left alone. Its `LowCoverage` verdict on an
off-centre-AF graph indel is now harmless, because the MSA's verified call
overrides it where one exists -- and where no MSA call exists, the strict verdict
is the conservative one.

## The same bug one category over: a repeat indel also outranked its MSA call

`48,177,725` and `48,234,100` were the last two hiphase sites the arm did not
use. Both sat in the candidate table as `REP_HET_INDEL` with `HAP_ALT = 0`,
`HAP_REF = 0` and `PS = 0` -- screened out of k-means by construction, so never
assigned a haplotype and never emitted. `48,177,716-48,177,915` is one of the
41 noisy regions, so the MSA runs over them.

The merge preferred the screened candidate for the same structural reason as
before, one condition earlier:

```cpp
const bool replace_repeat =
    (admit_all_in_region ||
     (site_whitelist != nullptr && site_whitelist->find(key) != site_whitelist->end())) &&
    old_vars[old_i].counts.category == VariantCategory::RepeatHetIndel && ...;
```

The swap required region-trust mode or a whitelist hit, so a plain run always
kept the screened, unusable version and discarded the MSA's verified call. The
`admit_all_in_region || whitelisted` requirement is dropped: a `RepeatHetIndel`
carries no phase information by construction, so the MSA's call at the same key
strictly dominates it.

| arm | blocks | sites used (+/-6 bp) | sites used (event) | reads tagged | concordance |
|---|---:|---:|---:|---:|---|
| site pruned | 2 | 15/19 | 16/19 | 393 | **100.00%** |
| MSA replaces pruned | 2 | 16/19 | 17/19 | 393 | **100.00%** |
| **+ MSA replaces repeat** | 2 | **18/19** | **19/19** | 393 | **100.00%** |

Both sites come back matching hiphase's alleles *and* genotypes exactly --
`48,177,725 T>TA` at `1|0` and `48,234,100 C>CCT` at `0|1`, both
`NOISY_CAND_HET`, joined to the blocks at `48147225` and `48225786`. Reads
tagged and read concordance are unchanged at 393 and 100.00%, both blocks still
internally consistent.

**The arm now uses every heterozygous site hiphase phases in this window.** The
only position the narrow matcher still reports missing is `48,149,567`, which
the arm emits as two records at `48,149,548`.

What remains is not a site deficit. The window is still two blocks where hiphase
emits one, and the split is the link `48,204,383 -> 48,225,786`: of 173 reads
overlapping the pair exactly 1 carries a usable allele at both ends, while 7
span the positions.

## Why the window is still two blocks: the overlap test, not the site set

With every hiphase site now used, the split `48,204,383 -> 48,225,786` is the
only thing left, and it is a read-observation problem. The seven reads that span
both positions separate **perfectly** by haplotype at the right-hand site:

| truth | net length at 48,204,383 | net length at 48,225,786 |
|---|---|---|
| PATERNAL | -1, 0, -1, 0 | **-7, -6, -6, -5** |
| MATERNAL | 0, -2, -2 | **-1, 0, 0** |

Our records there are `CA>C` (1 bp) and `CAAAA>C` (4 bp), so no paternal read
matches either length exactly, and the link scored one usable read out of seven.

The obvious repair -- resolve a non-exact overlapping indel by the closer
hypothesis, which is the rule this project's own genotyping tooling validated --
was implemented and measured, and it is **inert on this window**. Instrumented,
the branch is reached 1000 times, 861 of those inside a repeat tract, and it
decides 538 observations that were previously `-1`; the arm's VCF and BAM come
out byte-identical, and there is **no hit at all** at `48,225,78x`. It was
reverted rather than shipped unmeasured.

That absence is the finding. Those reads never reach the non-exact branch,
because `profile_ovlp_var_site` does not consider their deletion to overlap the
candidate in the first place -- the aligner has placed the event at a different
offset inside the A run, not merely at a different length. So the fix has to be
in the **overlap test**: inside a repeat tract, a same-kind indel anywhere in
the tract is the same event and should be compared by net length over the tract
window, which is exactly the rule already recorded for genotyping these sites by
hand ("never at the anchor position"). That is a change to how reads are matched
to candidates, wider in blast radius than anything in this file so far, and it
is the next thing to do.

## Not a mapping-quality gate: the observation is lost in the MSA allele call

Mapping quality is ruled out at the read level. All seven reads spanning
`48,204,383 -> 48,225,786` are **MAPQ 60**, primary, and none is skipped
(`is_skipped` is written in exactly one place, `gap_recovery.cpp:50`, which this
arm never reaches). `min_read_hap_margin` is 0 and the arm does not pass a
margin, so no read-tagging gate applies either.

The chain, probed end to end:

1. **The site does not exist when read profiles are built.** At
   `collect_var_build_profiles` time the candidates near there are `48,204,379`,
   `48,204,385` and `48,225,795` -- `48,204,383` and `48,225,787` are absent.
   They are **MSA-created**, so their per-read alleles never come from the digar
   comparison path at all, which is why the closer-hypothesis change to that
   path was inert.
2. **In the final profile the observation is simply missing.** At the two split
   records at `48,225,787` (`ref_len` 1 and 6, both `cate=0x100`), **six of the
   seven reads hold `allele = -1`**, and those same six hold `hap = 0`. The one
   read with an observation is the only one that is phased.
3. **The rescue for that case never runs.** `add_msa_site_observations` -- which
   re-calls a site against both consensuses and accepts it when two independently
   composed paths agree -- is invoked unconditionally, but probed at this site it
   reports `unassigned = 0`. The list is empty.
4. **Because this region takes the branch that cannot fill it.**
   `collect_noisy_reg_aln_strs` picks between two paths, and with
   `--verbose 1` the region reports `BranchSelect ps=-1 n_full_reads=72
   n_reads=76`. `ps = -1` means no phase set here has reads from both
   haplotypes, so the hap-aware path is skipped and
   `wfa_collect_noisy_aln_str_no_ps_hap` runs -- and that function **has no
   `unassigned` parameter at all**. Only the hap-aware path takes one
   (`align.cpp:1847`).
5. **So a clustered read whose allele call fails has no fallback.** In
   `update_cand_var_profile_from_cons_aln_str2`, a full-cover read gets
   `allele_i = get_var_allele_i_from_cons_aln_str(...)` for a variant from its
   own cluster's consensus, `0` for one from the other cluster, and whatever the
   first returns -- including `-1` -- is written straight into the profile. For a
   deletion the aligner placed elsewhere inside the A run, that call fails, and
   nothing re-tries it.

The consequence is circular, which is why it is stable: these reads are unphased
because they have no observation at the site, and they have no observation
because the rescue is reserved for reads that failed cluster assignment in a
branch this region does not take.

Two fixes follow, and they are not equivalent:

- **Reach the existing evidence standard from this branch.** Give
  `wfa_collect_noisy_aln_str_no_ps_hap` the same `unassigned` output, and route a
  clustered read whose allele call returns `-1` through it, so
  `call_msa_site_with_context` and the two-path agreement test decide the site.
  This adds no new standard -- it applies the one already used for ambiguous
  reads to a case that currently skips it.
- **Fix the allele call itself**, by reading the read's alignment over the
  variant's repeat tract by net length rather than at the key. Stronger, but it
  changes how every MSA site is called, so it wants the read-level gate on both
  windows before it is trusted.

The first is the smaller change and is the next one to make.
