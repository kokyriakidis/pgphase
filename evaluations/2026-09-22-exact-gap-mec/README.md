# Exact injected-site gap recovery (2026-09-22)

The retained recovery adds a bounded exact MEC pass after the ordinary
left-to-right seam stitch. It optimizes only exact BAM-injected sites in a
read-connected component joining both graph flanks. Existing imported phase
sets are atomic variables, and the same and flipped right-block orientations
are solved separately. A tied optimum abstains. Exact MEC is exponential, so a
seam with more than 20 variables abstains before search; the chr20 target panel
uses at most 13, while an unrelated 52.8 kb seam contained 131.

The final ordinary-pass fix also permits a new component to extend from a
BAM-injected upstream site when its strongest pair has the aggregate fallback's
eight-read net margin. This closes the high-support 58,834,248--58,836,037 gap.
The resulting phase set is 587/610 (96.23%) consistent with parental truth.
The exact MEC pass closes 18,194,808--18,218,259 and
31,886,648--31,901,501; the available crossing truth reads are 100% consistent
at both targets (one and two reads respectively).

A broader adjacent-pair MEC experiment closed 40/48 coordinate cases, but it
created a confirmed parental-orientation switch across
32,035,459--32,050,364 and raised discordant reads to 40,515. Allowing every
post-break strongest graph edge similarly flipped a 6,418-read block at
1.907 Mb. Both arms were rejected. The retained post-break edge must originate
at an injected BAM site; the equivalent graph-to-graph shortcut abstains.

The full chr20 comparison uses the same 500 kb chunks, 16 workers, input BAM,
graph catalog, GAF, and diplinator truth BAM. The retained arm spans 31/48
coordinate cases (30/47 distinct regions because 4.78 Mb occurs under two
coordinate conventions), up from 28/48 in the exact current-source control.
The remaining 17 coordinate cases have tied/disconnected injected observation
graphs or exceed the exact variable bound. Forcing them would choose a block
parity without production evidence. Read discordance rises by 860 while phased
reads rise by 405; this is accepted because the three new joins preserve the
observed parental orientation and the project preference allows noisier read
assignments when haplotype-to-haplotype orientation remains correct.

`results.tsv` records the chromosome totals. `no_mec_control` is the current
source with the MEC pass disabled. `adjacent_mec_rejected` is diagnostic only.

## Existing switch-window audit

The integration panel still reports wrong-orientation joins at 6.578, 22.981,
and 48.226 Mb. An exact A/B run with the pre-MEC binary produced the same
failures byte for byte, so the retained MEC and injected-endpoint rule did not
introduce them. Raising every ordinary recovery edge to margin eight corrected
the 6.578 Mb window, but left 22.981 and 48.226 Mb switched and lost valid 5.31
and 60.03 Mb closures; that global threshold was rejected. These older windows
need a separate orientation-quality criterion and are not counted among the 48
HiPhase-correct gaps targeted by this experiment.
