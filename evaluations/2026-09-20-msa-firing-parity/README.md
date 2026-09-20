# Does the noisy MSA fire the same times, on the same regions, in both tools?

Answered by instrumenting both sides and diffing the per-region firing lists.
Upstream already logs every firing at `-V 1`: `Hap`, `NoHap` or `Skipped` with
the region, read count and fully-covering count (`align.c:1792`, `:1798`,
`:1800`). Ours was given the same three-way label per region.

## The answer: 53 of 53 regions, 52 of 53 branches

`CHM13#0#chr20:30,700,000-30,900,000`, both tools run directly:

| | regions | `Hap` | `NoHap` | `Skipped` |
|---|---:|---:|---:|---:|
| upstream longcallD | 53 | 43 | 1 | 9 |
| pgphase | 53 | 42 | 1 | 10 |

**Identical region bounds: 53. Upstream-only regions: 0. Ours-only: 0.** Region
formation is at parity -- same count, same boundaries.

**One branch disagreement**, and it is the cluster from the previous record:

| region | reads | fully covering | upstream | ours |
|---|---:|---:|---|---|
| 30,802,901-30,804,117 | 16 | 3 | **`Hap`**, n_cons 2 | **`Skipped`** |

## Why that one region differs -- and it is not the MSA logic

The branch is chosen the same way on both sides (`align.cpp:1942-1952` against
`align.c:1789-1801`): guided when a phase set carries both haplotypes,
otherwise unguided when `n_full_reads >= min_dp`, otherwise the region is
skipped. Dumping the phase-set tallies on the SAME span:

| | full (hap1/hap2) | all (hap1/hap2) | minor | outcome |
|---|---|---|---:|---|
| upstream | **1 / 2** | 3 / 13 | 1 | ps valid -> `Hap` |
| ours | **0 / 3** | 0 / 5 | 0 | ps = -1 -> `Skipped` |

Both require the minor haplotype to hold at least one fully-covering read, so
with `0 / 3` we correctly return -1 and with `1 / 2` upstream correctly does
not. The tallies differ because the READ PARTITION differs. By read length in
the region:

| | haplotype A | haplotype B |
|---|---|---|
| upstream | 1053, 5692, 8805 | 1048, 1052 + the eleven long reads |
| ours | the eleven long reads | 1048, 1052, 1053, 5692, 8805 |

Which label is A and which is B is a free gauge and does not matter. What
matters is that upstream places the short reads `1048` and `1052` with the long
reads, and we place them with `1053`, `5692` and `8805`. The three
fully-covering reads are exactly `1048`, `1052` and `1053`, so upstream's
partition splits them 1/2 across haplotypes and ours puts all three on one side.

That single grouping decision is the whole difference between the MSA firing and
not firing on this region.

## Verified identical, so ruled out as causes

Checked line for line or by measurement: the branch conditions themselves; the
phase-set derivation (`collect_phase_set_with_both_haps` against
`align.c:1225-1272`, same full/partial accounting, same minor-haplotype checks,
thresholds 1 and 2 -- and on a 15 kb span the two produce byte-identical
tallies, `full=0/3 all=0/5`, and both skip the region); the abPOA setup (`wb -1`,
`inc_path_score 1`, `cons_algrm ABPOA_MF`, `min_freq = min_af`); the scoring
defaults (match 2, mismatch 6, gaps 6/2/24/1); the noisy-region skip rule; all
seventeen pipeline thresholds; and the 500 kb chunk size.

So achieving firing parity does not require changing anything in the MSA path.
It requires the read partition to agree, which is a solver-level question.

## Correction to the previous record

Commit `6cea093` states that our branch gate diverges from upstream, "which
bails only on `n_full_reads == 0` (align.c:1175)". That citation is the INNER
function's own guard, one level below the branch. Upstream's branch gate is
`align.c:1794`, `else if (ps_with_both_haps <= 0 && n_full_reads >=
min_no_hap_full_read_count)`, with `min_no_hap_full_read_count = opt->min_dp`
(`align.c:1777`) -- identical to ours. There is no gate divergence; the
measurement in that record (the relaxation is neutral on GIAB and costs read
accuracy) stands and is now doubly justified, since relaxing our gate would fire
the MSA on regions upstream skips.

## Next experiment

Whether upstream's partition or ours is correct at this locus is testable
against parental truth, but it needs the read NAMES: upstream's `names[]` array
is unpopulated at that point in its code path, so the dump returned `?` for
every read and the truth join could not be made. Dumping names from the read ids
on both sides is the next step, and it decides whether this is a firing-parity
gap to close or a place where our partition is the better one -- as was already
shown at `26,591,362`, where our labelling was 8/8 truth-consistent against
upstream's 7/8.

## Tooling

The instrumented upstream lives at `/tmp/lcd-instr`, a copy -- the user's
`~/Downloads/longcallD` checkout is untouched. Probes there are gated on
`LCD_PSDUMP` and `LCD_READS`.

## Whose partition is right? Truth cannot say here -- claim withdrawn

An earlier reading of this record said our labelling was wrong, because the
three fully-covering reads are `062403_s3/42341906` (MATERNAL),
`034919_s2/144706284` (PATERNAL) and `034919_s2/218562686` (PATERNAL) and we put
all three on one haplotype while upstream splits them 1/2 -- matching that
composition. That conclusion does not hold, and the contradiction that broke it
is worth recording.

**Our labels are internally consistent.** Scoring all reads against our own 165
phased clean SNPs in the enclosing block, the two supposedly misassigned reads
agree with the haplotype we gave them at **152 of 153 sites**. A locally
misassigned read does not look like that.

**The parental truth map is self-inconsistent in this window.** Assigning each
site's alternate allele to a parent (>= 5 alternate reads, >= 90% pure) and
walking the block, the parent attached to hap1 alternates **42 times across 147
sites** -- often site to site. No real phasing produces that, and it cannot
coexist with reads matching one haplotype at 152 of 153 sites. This is the
low-mapping-quality, duplicated window already on record (21 of 28 reads at
MAPQ 1-19), and it is outside the GIAB confident regions, so neither the
benchmark nor read truth can adjudicate the partition here.

So the honest status of the single firing divergence: its mechanism is
established (the read partition decides the minor-haplotype count, which decides
the branch), but which partition is correct is **not** established, and this
locus cannot establish it.

## The harness

`compare_firing.sh REF BAM REGION` (with `LCD` pointing at the longcallD
binary) diffs the two firing lists directly. Ours now prints one comparable line
per region at `--verbose 1`:

    MsaFire <Hap|NoHap|Skipped> <beg>-<end> <len> <n> reads (<f> full) ps=<ps>

against upstream's `Hap` / `NoHap` / `Skipped` lines. This is a verbose-only
change; no behaviour is affected, and all three suites are unchanged (unit 3/3,
predicate 151, window 125).

Any future change that moves where the MSA fires is now visible in one command
rather than by inference.
