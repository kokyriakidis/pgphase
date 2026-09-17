# Retry the window the solve could not phase, with the BAM's sites admitted

Instead of leaving an unphasable region to a second recovery pass, detect it in
the normal pass -- the solve assigned no phase set there -- admit the BAM's own
sites **in that window only**, and solve again. `--retry-unphased-with-bam`.

## How the pipeline signals the failure

`update_read_phase_set` (`collect_phase.cpp`) grants a read a phase set only from
an eligible **heterozygous** candidate: homopolymer indels, `kCandNoisyCandHom`
and unsupported MSA insertions are skipped, and a read with no such candidate
gets `ps = -1`. The VCF writer emits `PS` only for phased hets
(`collect_output.cpp:582`). So an unphasable window is not labelled, and inside
`chr20:48,176,830-48,229,446` **170 reads lie wholly within the gap and none
carries a phase set**.

Two ways of measuring that are wrong, both of which I used before getting here:

- **Counting `|` in the genotype.** A homozygous record is written `1|1`, so it
  reads as phased while carrying no linkage. This gap has 2 phased het records
  and 64 hom ones; counting separators reported 66.
- **Asking whether reads in the window carry a phase set.** A read reaching in
  from a flanking block carries one earned outside the window, which masks the
  failure. The detector therefore marks a bin as phased only where an eligible
  het candidate with a phase set exists, using the same test as
  `update_read_phase_set` (including the `hap_to_cons_alle[1]`/`[2]` indices --
  `[0]` is not a haplotype allele).

## What it does

`collect_unphased_windows` bins the chunk at 1 kb, marks bins holding no eligible
phasing candidate, and returns runs of them covered by at least
`--retry-min-unphased-reads` reads and at least `--retry-min-window-bp` wide.
For candidates inside those windows only, the saved category is restored -- the
categories are already saved in `discovery_flags` for the recovery pass, so
nothing new is stored -- and `collect_var_run_phasing` runs again on the chunk,
with `skip_noisy_kmeans` disabled for that call alone.

The confinement is the point. Admitting these sites chromosome-wide on chr20
doubles the read Hamming error, 0.878% -> 1.837% (1,831 -> 3,415 discordant
reads), fragments blocks 238 -> 408, and spans only 14 of the 196 gaps.

## Where it stands on the first gap

`chr20:48,176,830-48,229,446`, flag on, `--verbose 1`:

```
retry window 48107000-48162000 (55.0 kb): 360 candidates, 360 graph sites
retry window 48163000-48176000 (13.0 kb): 170 candidates, 170 graph sites
retry window 48177000-48229000 (52.0 kb): 662 candidates, 661 graph sites
retry window 48272000-48300000 (28.0 kb):  94 candidates,  94 graph sites
retry: 4 unphased window(s), 1 site(s) admitted
```

The detector lands exactly on the gap -- `48,177,000-48,229,000`, starting just
after the single phased het at 48,176,831. The retry admits `48,183,977`
(`CLEAN_HET_INDEL`) and the left block grows from 2 sites
(`48,162,480-48,176,830`) to 3 (`48,162,480-48,183,976`, +7.1 kb), with in-gap
phased hets 2 -> 3. Gate clean: 0 concordant->discordant, 0 tags lost, 0 newly
discordant.

## The bug: the noisy-region MSA is skipped whenever recovery is on

`collect_var_run_phasing` (`collect_var.cpp`):

```cpp
    if (!opts.recover_gaps)
        collect_noisy_vars_step4(chunk, opts, noisy_site_whitelist);
```

So with `--recover-gaps` the noisy-region MSA never runs in the normal pass; it
is deferred to the recovery pass. That is why the failed window holds **no
`NoisyCandHet` candidate at all** -- 1 clean het SNP, 2 clean het indels and 606
low-coverage catalog sites -- while the BAM channel run on the same interval
calls **8 noisy het sites** and phases the whole window into one block. The
candidates were never created, so nothing was there to admit, and the first
version of this retry could only re-admit the single non-catalog site.

Three explanations were checked and ruled out first: the 50 kb
`max_noisy_reg_len` cap (the BAM channel produces the same 8 sites whether run
over 52 kb or 152 kb), the `private_keys` branch that zeroes that cap (it is
entered only with `--private-sites-vcf`, which these runs do not pass), and
`skip_noisy_kmeans` (it is read *inside* the step that never ran).

The retry is the second try that deferral assumes, so it now runs that step for
its own call: `retry_opts.recover_gaps = false` plus
`skip_noisy_kmeans = false`.

| | retry off | retry on |
|---|---:|---:|
| phased het records in the gap | 2 | **8** |
| left block | `48,162,480-48,176,830` (2 sites, 14.3 kb) | `48,147,227-48,229,226` (**13 sites, 82.0 kb**) |
| gate: concordant -> discordant | -- | **0** |
| gate: tags lost / newly discordant | -- | 0 / 0 |

The evidence recovery works: 8 in-gap phased hets matches what the BAM channel
finds on its own, and the left block now reaches to within 220 bp of the right
one.

## But the block it forms is switched, and the read gate cannot see it

That 82 kb block spans `48,183,976 -> 48,225,786`, which is **41.8 kb with zero
reads covering both sites** (the neighbouring spacings carry 41 and 51; the
longest read there is 29.9 kb). Scoring its sites against read truth:

| site | side of the hole | hap1 carries | confidence |
|---|---|---|---:|
| 48,162,480 | left | PAT | 1.00 |
| 48,176,830 | left | PAT | 0.96 |
| 48,183,976 | left | PAT | 1.00 |
| 48,225,786 | right | **MAT** | 0.99 |
| 48,229,226 | right | **MAT** | 0.93 |

Left and right are on opposite haplotypes: the block asserts a phase across a
hole nothing supports. The read-level gate reports zero flips because no read
spans the hole, so it is structurally blind to this class of error -- only the
site-level truth check sees it.

## Next, and it is a rule already validated here

The retry must not join across a spacing with no spanning reads. The correct
output for this window is two blocks -- `48,147,227-48,183,976` and
`48,225,786-48,229,226` -- each anchored to its own flank, with the 41.8 kb
between them left open. That is the same conclusion the window-size comparison
reached (`evaluations/2026-09-16-window-stitch/`), where 20 kb chunks produced
exactly those two blocks and a 500 kb chunk produced this switched one.

Defaults unchanged (`retry_unphased_with_bam = false`); all five unit-test
binaries pass.

## The genotype collapse, and the fix

The retry recovered the evidence but the window still would not phase, because
the site that bridges it was emitted homozygous. `iter_update_var_hap_to_cons_alle`
recomputes each haplotype's consensus allele **independently by majority**,
except for verified multi-allele MSA insertions, which get a joint orientation
with this comment on it:

> Independent haplotype majorities can select the same allele twice.

That is exactly what happens to a plain biallelic site inside a window the first
solve could not phase: the reads there carry no haplotype labels, so both
majorities are the deeper allele, `hap_to_cons_alle[1] == hap_to_cons_alle[2]`,
and the record is written `1|1`. It is self-sustaining -- a hom site links
nothing, so the window stays unphasable and the collapse repeats. A probe inside
the iteration caught it in the act at `chr20:48,204,383`: `cons1=1 cons2=0` on
one round, `cons1=1 cons2=1` on the next.

Both of our channels did this, not just the hybrid one: `collect-bam-variation`
on the same interval also emitted `1|1` with `HAP_ALT=3`. Only hiphase called it
`0|1`, and it is the only heterozygote between 48,183,976 and 48,225,786 -- the
difference between a read-chained bridge and a 41.8 kb jump with no spanning read.

The fix applies the same joint orientation to a biallelic candidate whose allele
depths call it het (`ref_cov`/`alt_cov` over `min_alt_depth`, AF within
`min_af`..`max_af`), confined to the windows the retry is re-solving. Where the
labels carry no preference at all it seeds a het rather than a hom: an arbitrary
orientation is resolvable by the link votes, a collapsed one is not.

### Verified with `verify_retry.py` on two windows

`chr20:48,176,830-48,229,446`:

| | retry off | retry on, before this fix | retry on, now |
|---|---:|---:|---:|
| usable het sites in region | 2 | 8 | **10** |
| unsupported links (0 spanning reads) | 0 | **1 (41.8 kb)** | **0** |
| category-het / genotype-hom | 0 | **1** | **0** |
| competitor sites we call hom | 0 | **1** | **0** |
| competitor sites absent from our VCF | 34 | 21 | **21** |
| gate: conc -> disc / tags lost | -- | 0 / 0 | **0 / 0** |

The 41.8 kb unsupported link is gone because `48,204,383` now sits inside it as a
usable het and both halves carry spanning reads. The verdict is still FAIL, on a
switch between `48,147,227` and `48,149,548` -- 2.3 kb apart with 59 spanning
reads, so a supported link that is oriented wrong, and left of the gap rather
than inside it -- plus three sites below the confidence floor.

`chr20:36,217,274-36,268,291`, the same two arms:

| | retry off | retry on |
|---|---:|---:|
| usable het sites in region | 4 | **15** |
| switches | 1 | **0** |
| competitor sites absent from our VCF | 14 | **7** |
| gate: tagged | 358 | **240** |
| gate: concordant | 278 | 239 |
| gate: conc -> disc / lost / new | -- | 0 / **75** / 12c |

Accuracy on what remains improves sharply (77.65% -> 99.6%) and the pre-existing
switch disappears, but **75 correct read tags are lost against 12 gained**. So
this is not a default: the re-solve judges every read in the chunk against the
new site set, and on this window that costs coverage. The same shape as the
earlier filter-reorder regression, and the next thing to measure.

## Checking the fixed window for further bugs

Two things the verifier still fails on, and what each turned out to be.

### The duplicate-locus records are the MSA path's, not the retry's

At `48,147,227` two deletion records sit at one position in the same phase set
with **opposite** haplotype assignments and complementary depths: `GGT>G`
(`0|1`, 10 ref / 5 alt) and `GGTGTGTGTG>G` (`1|0`, 5 ref / 10 alt). They are the
2 bp and 9 bp forms of one GT tandem repeat (`GGT>G` deletes 2 bases,
`GGTGTGTGTG>G` deletes 9; an earlier version of this paragraph said 8 bp), decomposed into separate biallelic
candidates and then phased independently, which is contradictory as a genotype:
one locus cannot carry both.

Duplicated positions among phased hets, by arm:

| arm | duplicated / phased hets |
|---|---|
| retry off (noisy MSA not run) | 0 / 44 |
| `collect-bam-variation` (MSA on) | 2 / 8 |
| hybrid `--keep-noisy-kmeans` (MSA on) | **4 / 63**, same cluster, same opposite haps |
| retry on | 4 / 66 |

So this is a pre-existing defect of the noisy-MSA path, present wherever that
path runs; the retry arm shows it only because it runs that path. An earlier note
here attributed it to the het-seeding change -- that was wrong, and the
correction matters because the seeding was not what made these records: they are
phased by the ordinary per-haplotype majority once the MSA has created them.

A rule that seeded at most one candidate per overlapping locus was written and
then **reverted**: the output was byte-identical with and without it, because
suppressed candidates still reach the ordinary path and are phased there. Fixing
this belongs where the alternatives are created, not where they are oriented.

### The switch is anchored on that cluster

`verify_retry.py` reports one switch, `48,147,227 (MAT) -> 48,149,548 (PAT)`,
2.3 kb apart with 59 spanning reads. The left side of it is the repeat cluster
above, whose three records disagree with each other: `MAT` at 0.86, a tie at
0.50, and `MAT` at 0.57 (below the confidence floor). So the block's left anchor
is a contradictory cluster at one locus, and the switch flag follows from it
rather than from a long-range orientation error. It is also outside the gap
interval, left of it.

Both findings point at the same place: the noisy MSA should emit one multiallelic
record per locus, or mark its decomposed alternatives so they cannot be oriented
as independent heterozygotes.

## The retry was disabling five correctness guards, one of which mattered here

`recover_gaps = false` on the retry's `Options` copy was a proxy for "run the
noisy-region MSA". It is not one: five guards inside the phasing path are gated
on that field, and the retry cleared all of them --
`split_nested_msa_deletions` (`collect_phase_noisy.cpp:1145`), the paths at
`collect_phase_noisy.cpp:1151` and `:1624`, `collect_var.cpp:2025`, and the
`recovery_graph` composite.

`split_nested_msa_deletions` is the one that produced the duplicate-locus
records, and its own comment states the consequence:

> Both haplotypes delete the common prefix. Treating each full deletion as an
> independent het lets a longer-deletion read support both ALT rows,
> manufacturing a false phase bridge.

At `48,147,227` the two haplotype consensuses carry a 2 bp and a 24 bp deletion
of the same GT tandem repeat. With the splitter disabled both full deletions were
emitted as independent hets with opposite haplotype assignments. With it
restored they decompose correctly into the common deletion, homozygous because
both haplotypes carry it (`GGT>G`, `1|1`), and the residual het
(`TGTGTGTGTGTGTGTGTGTGTGT>T`, `1|0`).

The fix replaces the proxy with a named flag, `force_noisy_msa`, so the retry
asks for that one step and keeps every other recovery guard in force.

| `chr20:48,176,830-48,229,446` | before | after |
|---|---:|---:|
| duplicated positions among phased hets | 4 | **0** |
| switches within a phase set | 1 | **0** |
| unsupported links | 0 | 0 |
| bridging site `48,204,383` | `0\|1` | `0\|1` |
| gate: conc -> disc / tags lost | 0 / 0 | 0 / 0 |

The verdict is still FAIL, now only because two sites fall below the verifier's
confidence floor (`48,177,788` and `48,202,056`, the second inside the region),
which sets `gate_blind`. Also still open, and both now visible rather than
hidden behind the duplicates: `48,177,780` is called `NOISY_CAND_HOM` where the
competitor calls a het, and `48,173,990` is a category-het / genotype-hom site
outside the region.
