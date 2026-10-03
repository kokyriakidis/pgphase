# Close an NG50-priority gap with preserved physical boundary calls

## Result

Close **chr20:5,511,231–5,531,924**, joining the complete VCF phase blocks
5,393,615–5,511,231 and 5,531,924–6,148,782 into a **755,168 bp** block.
The full chromosome confirms the predicted NG50 gain. Every read's prior
truth-correctness state survives; the join introduces no newly incorrect read.

| Full chr20 metric | Before | After |
|---|---:|---:|
| VCF block NG50 | 633,077 bp | **643,699 bp** |
| VCF block N50 | 774,189 bp | 774,189 bp |
| Truth-scored phased reads | 237,199 | 237,199 |
| Truth-correct assignments | 230,072 | 230,072 |
| Discordant assignments | 7,127 | 7,127 |
| Conditional read accuracy | 96.995350% | 96.995350% |
| Read phase sets with truth-scored reads | 661 | 659 |
| VCF phase blocks | 333 | 332 |
| VCF variant keys | 63,630 | 63,630 |
| Coordinate panel spans | 86/100 | **87/101** |

NG50 increases by **10,622 bp (1.68%)**. It uses half the reference chr20
length, **66,210,255 bp**, as its cumulative-span threshold. N50 uses half the
sum of the emitted phase-block spans, so these metrics can change differently.
Both use inclusive extents of positive-PS phased heterozygotes. Overlapping
extents are counted separately; this is span NG50, not error-corrected NGC50.

Only the targeted pair of old VCF blocks merges. No phased SNP disappears or
old SNP block acquires a mixed allele gauge. No variant key or unordered
genotype changes. The 584 changed read HP/PS tuples and 124 changed VCF sample
rows reflect the connection; truth scoring retains all 230,072 previously
correct and all 7,127 previously discordant reads.

## Choosing the target

Audit the accepted chromosome against saved same-BAM HiPhase results. Exclude
25–30 Mb and require at least 95% competitor block purity on ten local overlaps,
plus at least five truth-scored reads in each disjoint 10 kb flank with matching
parental orientation. Deduplicate the resulting nominations and simulate joining
the **whole** adjacent phase-set extents. Rank by the resulting NG50 gain.

The audit supplies 38 distinct nominations, not an exhaustive census of all
unclosed gaps. Nine have a positive modeled single-join NG50 gain. Four tie for
the largest gain of 10,622 bp: 5.511, 57.085, 61.738 and 62.623 Mb. The selected
5.511 Mb gap has coherent boundary representations and flank orientations.
The other candidates have representation or source-gauge issues; the 62 Mb
counterexample includes an internal graph-block switch and remains protected.
A large projected block does not establish a safe haplotype connection.

`ranking.json` records the starting VCF hash, exact extents and competitor
evidence. `after-ranking.json` confirms the measured NG50 and ranks the 37
remaining baseline nominations. The next modeled gains at 57.085, 61.738 and
62.623 Mb are 29,299 bp each; none of those joins is certified by ranking alone.
The saved competitor outputs are from the October 1 same-BAM comparison under
`/tmp/pgphase-hiphase-comparison-2026-10-01/`; no new competitor runtime is claimed.

## Defect and retained repair

Two MAPQ60 physical reads call both boundary deletions:

| Read | Left `CAAT→C` | Right `TACAC→T` |
|---|---|---|
| `m84031_231217_034919_s2/213455387/ccs` | ALT | REF |
| `m84031_231217_062403_s3/126357809/ccs` | REF | ALT |

The graph read profiles end before the right boundary, but the BAM alignments
reach it. Targeted recovery retains both paired calls. The larger whole-flank
BAM validation solve omits the right MSA deletion observation on these partial
reads. Its original boundary rule also excludes this recalled deletion and
screens for exactly one spanning molecule. Thus increasing the solve context
does not by itself preserve the decisive physical evidence.

Keep the established physical certificate first, using its original boundary
eligibility and unchanged source solve. Only when it fails, consider a
downstream MSA-verified heterozygous deletion before the first shared graph
site. Screen physical coverage before the expensive whole-flank solve and
reuse that solve for both certificates. In the deletion fallback, apply the
existing targeted-recovery observation backfill over the seam, preserving
discovery counts, genotypes, source phase sets and source orientations.

The fallback requires multiple MAPQ30-or-higher callable molecules, evidence
on **both haplotypes**, and agreement of every callable allele pair. The exact
CIGAR alleles must agree with the BAM observations; deletion REF also requires
matching reference bases. Keep the full-block graph/BAM orientation and
matched-site coverage checks, source weak-cut rejection and normal stitch
conflict checks. The upstream eligibility and ordinary multi-read SNP stitch
path stay unchanged. No truth, competitor call, fixture coordinate, MSA row
merge or new allele realignment method enters production.

The final boundary genotypes are `1|0` and `0|1`, sharing PS 5,531,924. They
preserve the opposite-ALT relationship supported by the physical reads.

## Local accuracy and limits

The owning 5–6 Mb replay preserves **4,024 scored / 4,001 correct / 23
discordant** reads while its VCF blocks fall from four to three. The separate
upstream 5.31/5.39 Mb connection and source-cut regression remain intact.

There are 166 truth-scorable physical gap overlaps. Pgphase retains 139
phased overlaps, 136 correct and three discordant, in the same parental
assignments as before. The main joined block contains **131/131 correct**
local reads. The other eight tagged overlaps belong to an independent fallback
block (five correct, three discordant). Saved native-DV HiPhase has 136/136
correct in its main block. This repair certifies continuity without claiming
that pgphase has matched HiPhase's local block coverage or aggregate purity.

## Rejected experiments

- Filtering nongermline rows out of the trusted MEC matrix reopens the
  required 5.31 Mb gap and loses 15 correct owning-chunk assignments; the target
  remains split. Restore the original solver admission rule.
- Broadening both boundary sides to recalled indels, or accepting general
  multiple-read physical links, displaces the upstream source boundary. Trials
  lose correct reads or merge across the existing upstream cut. Restore them.
- Prioritizing the nearer downstream deletion **before** the original
  certificate closes this target locally but reopens the required 23.008 and
  47.714 Mb panel connections. Full chr20 loses two scored reads and three
  correct assignments, gains one discordant assignment, loses a variant key
  and has no NG50 gain. Preserve the original certificate first instead.

`rejected-boundary-first-parity.json` and
`rejected-boundary-first-transfer.json` record that failed full run. The retained
fallback leaves the 23–24 and 47–48 Mb owning read tags and VCF rows exactly
unchanged. No existing expectation is lowered to accommodate a trial.

## Regression and validation

Add the new coordinate case to the permanent panel with exact `spans=1`, an
in-gap heterozygote requirement and the measured local coverage floor. Extend
the existing 5 Mb owning-chunk test with opposite boundary genotypes, parental
orientation, whole-chunk correct/error limits and at least 131 correct local
reads in the main block. The starting executable fails three of these
assertions: shared PS, spanning and main-block read count. Existing protected
separations remain asserted. The panel total moves only for this measured
connection: 86 to 87 spans and 1,290 to 1,291 in-gap heterozygotes.

- Build and standalone units pass without new warnings.
- Predicate tests: **1,020 assertions / 42 cases**, all pass.
- Complete window suite: **7,092 assertions / 69 cases**, all pass.
- HiFi/ONT TSV and VCF goldens pass; HiFi one/four-thread output is identical.
- Run 105 fresh native panel requests with the final executable, then score
  verified outputs. `native-cache.json` and `replay_cached_panel.py` check the
  executable hash, every normalized CLI argument, input stats and output
  hashes. Requests absent from the cache execute natively.
- Whole chr20: 291.37 s at eight threads; native panel: 394.17 s with four
  concurrent workers. These concurrent validation timings are not a controlled
  runtime comparison.

Starting executable SHA256:
`b36dd6d8ea32b1f1707b64efc51878b2d8641344ef1b0f009f08eee94b1bac0b`.
Final executable SHA256:
`aa3b89bb6f6f5a49a4fc72e0c83577bfade869165e60f2c8556a48bf5e0026d4`.
Starting outputs: `test_data/tmp_gap_next34/full-contradictory/`.
Final outputs: `test_data/tmp_gap_next36/full-preserved/`.

Reproduce ranking from the repository root:

```bash
python3 evaluations/2026-10-02-ng50-gap-priority/rank_ng50_joins.py \
  --vcf test_data/tmp_gap_next34/full-contradictory/phased.vcf \
  --fai test_data/chm13v2.0.chr20.renamed.fa.fai \
  --chrom 'CHM13#0#chr20' \
  --targets evaluations/2026-10-02-ng50-gap-priority/competitor-targets.json \
  --output /tmp/chr20-ng50-priority.json
```

Replace the VCF with the final output to measure the retained NG50 and rank
the remaining nominations. `chromosome-parity.json` records read truth
transitions and variant parity; `block-transfer.json` records the sole new
block connection and every panel span. `validation-detail.json` confirms that
only VCF PS fields changed and records local read counts against the whole
phase-set parental orientation.
