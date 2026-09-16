# Why a gap with complete evidence still reports `split`

Single-gap diagnosis of `chr20:36,332,599-36,381,019` (48.4 kb), driven off
`--gap-decision-audit`, which exports the frozen proposal: every site with its
category, `msa_verified` and homopolymer flags, every read, and every per-read
allele observation from both channels.

```sh
./pgphase collect-hybrid-variation ... -q 1 --recover-gaps \
  --gap-decision-audit DIR --gap-recovery-report DIR/tiers.tsv
```

## The evidence to join this gap is complete

Injection is not the problem. Inside the gap the proposal holds 581 sites, of
which 16 are phase-informative categories, and every one carries read
observations:

- 2 `CleanHetSnp` at the gap edges (68 and 48 observations), segregation 0.971
  and 1.000 against the diplinator read truth.
- **14 `NoisyMsaHet` sites, all `msa_verified=1`, none flagged
  `is_homopolymer_indel`**, with 27-84 BAM observations each.

10 of the 16 are informative (segregation >= 0.90) and they form an unbroken
read-backed chain across the gap: consecutive spacings 11.4, 10.9, 10.9, 8.0 and
7.2 kb, with **17-81 reads spanning every consecutive pair**. Nothing about the
data prevents a join.

## The proposal still splits, and the vote matrix shows why

`stitch_gap_proposal` returns `left_linked=1, right_linked=1, joined=0`, because
`joined = left_linked && right_linked && links[0].ps == links[1].ps`
(`gap_recovery.cpp:376`). The audit's vote matrix shows the two ids in every
tier:

| tier | proposal PS | reads | left votes | right votes |
|---|---|---:|---|---|
| 1 | 36317511 | 20 | L11=5, L22=5 | — |
| 1 | 36381019 | 224 | — | R11=110, R22=91 |
| 2 | 36317511 | 20 | L11=5, L22=5 | — |
| 2 | 36354935 | 438 | — | R12=110, R21=91 |
| 3 | 36317511 | 20 | L11=5, L22=5 | — |
| 3 | 36354890 | 388 | — | R12=110, R21=91 |

Only **10 of the 637 reads in the window carry the left block's phase set at
all**, and they never merge into the substantive block. Cross-tabulating each
read's original assignment against its proposal block: 201 right-block reads
join the big block in every tier, 82-107 previously unphased reads join it too,
and the 10 left-block reads stay alone. **Zero reads in the big block overlap the
left anchor**, so the big block carries no left vote and `joined` cannot become
true regardless of tier.

## The gate: verified SNPs are refused where verified indels are admitted

`select_gap_link_sites` (`collect_phase.cpp:1188-1193`, pre-fix line numbers) decides which in-gap
sites may carry link support:

```cpp
const bool verified_indel = var.msa_verified && !var.is_homopolymer_indel &&
                            var.key.type != VariantType::Snp &&
                            var.lcd_var_i_to_cate == kCandNoisyCandHet && in_gap;
const bool verified_snp = opts.gap_bridge_private_snps && var.msa_verified &&
                          var.key.type == VariantType::Snp &&
                          var.lcd_var_i_to_cate == kCandNoisyCandHet && in_gap;
```

An MSA-verified **indel** in the gap may link unconditionally. An MSA-verified
**SNP** may link only when `--gap-bridge-private-snps` is passed, which is off by
default. Both classes come from the same MSA verification and carry the same
`msa_verified` flag.

Measured against truth, the policy is inverted relative to informativeness:

| class | n | median segregation | range |
|---|---:|---:|---|
| verified SNPs (refused by default) | 7 | **0.988** | 0.852-1.000, three at 1.000 |
| verified indels (admitted by default) | 7 | 0.895 | 0.657-1.000 |

The default-admitted class contains the two worst sites in the window (0.657 at
36,380,739 and 0.792 at 36,373,836); the flag-gated class contains the three
best.

Independent corroboration: the whole-chr20 `--gap-bridge-private-snps`
experiment joined exactly one gap, and it was this one --
`['CHM13#0#chr20', '36332599', '36381019']` -- with zero
concordant-to-discordant reads. Two routes, same locus, same mechanism.

## Correction: the gate is two chained gates, not one

The first version of this note named only `select_gap_link_sites` as the gate.
Dropping `opts.gap_bridge_private_snps &&` there alone is **inert** -- it was
measured and produced byte-identical tier statuses on this gap. There are two
gates with the same asymmetry, chained:

| location | decides | verified indel | verified SNP |
|---|---|---|---|
| `collect_phase.cpp:1191` (`verified_snp` in `select_gap_link_sites`) | may the site earn `gap_link_supported` | unconditional | flag |
| `collect_phase.cpp:896` (`recovered_snp` in the block-vote loop) | may the site cast a bridge vote | unconditional | flag |

The second gate additionally requires `var.gap_link_supported`, so the first is a
precondition for the second. Opening only the first leaves the site eligible but
voteless; both must open for a verified SNP to bridge.

## Fix applied and verified on this gap

Both flag terms removed, so a verified SNP is treated exactly as a verified indel
already was. Build clean with no new warnings, all five unit-test binaries pass.
`test_private_snp_bridge_anchor_requires_flag` pinned the old behaviour and was
updated to `test_verified_msa_snp_bridges_without_flag`, asserting the bridge
verdict holds with the flag off.

On `chr20:36,332,599-36,381,019` with the flag **off**, tier 3 now reports
`joined` with `LEFT_LINK_PS = RIGHT_LINK_PS = 36317511`, identical to what the
flag produced before:

| | before | after |
|---|---:|---:|
| phase sets in the window | 2 | **1** |
| single set spanning the gap | no | **yes** (55 variants) |
| reads tagged | 462 | 437 |
| reads concordant / scored | 456 / 462 | 433 / 437 |
| accuracy of the spanning set | — | 99.08% |

Read-level regression gate, matched by read name: **0 concordant -> discordant**.
433 concordant reads stay concordant, 4 discordant stay discordant, and **23
concordant plus 2 discordant reads lose their tags** -- the join costs 23 correct
read tags to close a 48.4 kb gap. That cost is the same kind the whole-chr20 flag
arm showed (61 reads), so this is verified on one gap only and is **not yet a
default**: it needs the chr20 gate run before the flag is retired.

## It does not generalize: 1 of 8 gaps

Re-ran the whole eight-window panel with the patched build, flag off, against the
matched pre-change runs (`compare_panel.py`).

| gap | gap bp | before | after | spans after | gate c->d | concordant tags lost | in-gap MSA het SNPs |
|---|---:|---|---|---|---:|---:|---:|
| 7,073,919 | 44,697 | split | split | no | 0 | 0 | 1 |
| 7,163,303 | 43,347 | split | split | no | 0 | 0 | 6 |
| 10,296,487 | 47,108 | split | split | no | 0 | 0 | 1 |
| 13,429,829 | 201,996 | partial | partial | no | 0 | 0 | 0 |
| 16,673,759 | 41,482 | split | split | no | 0 | 0 | 0 |
| 35,919,404 | 236,915 | split | split | no | 0 | 0 | **12** |
| 36,332,599 | 48,420 | split | **joined** | **yes** | 0 | 23 | 7 |
| 60,084,884 | 47,398 | partial | partial | no | 0 | 0 | 0 |

**One gap joins; the other seven are unchanged down to the individual read tag**
(no verdict change, no tag gained or lost). Panel totals: gate 0 concordant ->
discordant, 23 concordant tags lost, tagged reads 3832 -> 3807, concordant
3808 -> 3785.

The no-ops are **not** for lack of material. `35,919,404` carries 12 in-gap
MSA-verified het SNPs -- more than the gap that joined -- plus 15 in-gap clean
het SNPs available as anchors, and `7,163,303` carries 6, and both still report
`split`. So the flag asymmetry was a real blocker but only one of at least two,
and it is not the explanation for the `split` population as a whole. The next
diagnosis is the same audit on `35,919,404`: with anchors and verified SNPs both
present and the gates now open, whatever stops it is a third mechanism --
candidates are the bridge-vote thresholds (`min_block_link_reads` and the
`strong`/`snp_strong` vote requirements) and the anchor-eligibility restriction
below.

## Remaining fix direction

1. Treat the two verified classes identically in `select_gap_link_sites`: drop
   `opts.gap_bridge_private_snps &&` from `verified_snp`. This is a strict subset
   of what the flag does today -- the flag also admits private sites into
   `select_graph_gap_bam_reads` and requires BAM confirmation -- so the 61
   concordant reads that flag cost are not necessarily attributable to the link
   change, and the narrow version needs its own whole-chr20 gate measurement.
2. Second-order limitation, worth measuring after: the anchor loop in the same
   function only accepts anchors with `(anchor.lcd_var_i_to_cate &
   kCandGermlineClean) != 0`, so an MSA-verified site can be the *linked* site
   but never the *anchor*. MSA sites therefore cannot chain to each other -- in a
   gap whose only informative sites are MSA-derived, every link must reach one of
   the flanking clean anchors, which caps the bridgeable distance regardless of
   how dense the MSA evidence is.
3. `kCandAnchorClean` (`collect_phase.hpp:47`) is defined, documented as "the
   categories allowed to anchor k-means read assignment", and **never
   referenced**. Either enforce it or delete it; as it stands the comment
   describes a policy the code does not implement.

## Why the patch cannot fix `35,919,404`: no read crosses the interior

Audited the largest non-joined panel gap (`chr20:35,919,404-36,156,319`, 236.9 kb)
with the patched build. It has more of the material the fix unlocks than the gap
that joined -- 10 MSA het SNPs and 27 MSA het indels reported per tier, and 48
in-gap sites in phase-informative categories in the audit export, 19 of them
truth-informative (11 clean, 8 MSA) --
and it still reports `split` at every tier, `rejected` at the homopolymer tier.

The vote matrix shows a different shape from the joined gap: **both flanks link
strongly** (left PS 35873360 with 124/126 votes, right PS 36156319 with 54/60)
into two blocks with nothing between them, where the joined gap had a 10-vote
left remnant. Nothing is being refused here; the two halves simply never meet.

The reason is read coverage, not policy. Chaining the truth-informative sites and
counting reads that span each consecutive pair: 18 steps, of which **3 have zero
spanning reads** -- 22.9 kb, 26.5 kb and 43.5 kb, together 92.9 kb, 39% of the
gap. Repeating the step test over the wider 60-site linking selection from the
candidate TSV (a different set, not truth-scored), one 21.2 kb step still has
zero spanning reads. No site-admission
policy can bridge a step no read crosses.

## Panel-wide: which gaps a site-admission fix could ever reach

`classify_blockers.py` tests, per gap, whether a positional chain of called
linking sites runs flank to flank with reads spanning every step
(`blocker_classes.tsv`).

| class | gaps | span | note |
|---|---:|---:|---|
| `read_coverage` (a step no read crosses) | 2 | 438.9 kb | `13,429,829` (16.8 kb break) and `35,919,404` (21.2 kb break) |
| `chain_present` | 6 | 272.5 kb | coverage does not exclude a join |

Both large gaps are coverage-limited, which is why neither the phase-transfer nor
the injection experiment moved them either. But of the six where a chain exists,
**the patch fixed one**. So the five remaining have a read-backed chain of called
sites and no flag refusing them, and still split -- a third blocker, neither site
admission nor read coverage. Candidates, in order: the sites in those chains may
not be informative (at `35,919,404` truth scoring left 19 of the 48 in-gap sites
in phase-informative categories informative; the 60 in `blocker_classes.tsv` is a
different, wider selection from the candidate TSV and was never truth-scored, so
the two counts must not be chained -- running the truth scoring on those five
would settle it), the bridge
vote thresholds (`min_block_link_reads`, `strong`/`snp_strong`), and the
anchor-eligibility restriction that stops MSA sites anchoring each other.

Note on the label: `chain_present` means only that coverage does not rule a join
out. It is a necessary condition, not evidence that a gate fix will produce one.
