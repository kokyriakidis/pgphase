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

**Open, and the reason the gain is one site rather than several.** Of the 662
candidates in that window, 661 are graph sites, so the category zeroing withholds
only one of them. The BAM channel run alone on the same interval calls 11 het
sites there, so the rest are present at catalog positions and held back by
something other than the zeroing -- the noisy-candidate class and the
`skip_noisy_kmeans` default are the next thing to check, since the retry enables
that only for its own call and the earlier chromosome-wide arm with it enabled
reached 7 in-gap phased hets rather than 3.

Defaults unchanged (`retry_unphased_with_bam = false`); all five unit-test
binaries pass.
