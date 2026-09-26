# Chr20 recovery boundary evidence, 2026-09-24

The same pgphase graph, surjected alignment and truth inputs were used for the
before and after full-chromosome runs. These are read-truth scores on reads for
which the parental map provides a label; they do not establish correctness of
an individual phase-set join by themselves.

| Run | Phased truth reads | Correct | Discordant | Purity | VCF blocks | N50 |
|---|---:|---:|---:|---:|---:|---:|
| Prior one-haplotype source-path rule | 236,675 | 228,700 | 7,975 | 96.6304% | 410 | 486,139 bp |
| Verified clean indel/SNP boundary evidence | 236,675 | 228,700 | 7,975 | 96.6304% | 408 | 486,139 bp |
| Complete source path plus split-stable aggregate allele evidence | 236,675 | 228,701 | 7,974 | 96.6308% | 406 | 491,812 bp |
| One-neighbor weak-cut run attachment | 236,656 | 228,729 | 7,927 | 96.6504% | 398 | 491,812 bp |

The clean boundary test joins 14,235,594–14,239,930 and
36,016,138–36,018,966. The aggregate vote joins
32,490,058–32,490,150. These three boundaries previously lacked an exact
sequence-identical graph/BAM anchor. The 19,395,544–19,414,720 wrong-join
control remains split. The rule requires an intact original BAM source path,
exact MEC parity, a significant allele vote and agreement in both independent
read halves. An existing conflicting graph/BAM read gauge still vetoes it.

At the remaining four boundaries, 70, 58, 33 and 54 primary BAM reads
physically cross 0.528, 0.542, 17.62 and 37.46 Mb respectively; every one has
MAPQ >= 30. The recovery subsolve uses MAPQ 1, so read MAPQ filtering is not
the cause. At 17.62 Mb, only five of 33 spanning alignments have callable
alleles at both BAM boundary rows; 22 have an ambiguous left deletion call, and
the callable direct parity is 3:2. At 0.528 Mb, 49 of 70 spanning reads have an
ambiguous insertion call; only eight have both boundary alleles, split 4:4.
Nearby graph deletion rows at 17.62 Mb are allele-imbalanced and do not give a
consistent chain from the left phased SNPs. These are observation losses, not
molecule-length gaps.

The 0.542 and 37.46-Mb BAM source phase sets contain weak internal cuts. A
split-all prototype preserved the 19-Mb negative control but did not close
those gaps, so it was removed. Instead, the retained postpass identifies runs
between the original BAM weak cuts. It preselects the nearest oriented boundary
pair on each side and tests clean or MSA/alignment-verified candidates. Both
source alleles must occur, the full read vote must pass a one-sided exact
binomial `p <= 0.05`, and deterministic disjoint read halves must agree. When
only one neighboring phase set passes, only that run's rows and reads observing
no other source run move into its gauge. A root already used by a different run
of the same source PS cannot be reused. If both sides pass, the run remains
independent. This closes 0.542 and 37.46 Mb without changing the distant side
of their weak cut.

A simple two-sided local vote is unsafe: in the 19-Mb replay the left
graph/source pair has 8 same versus 49 cross reads, while the source/right
graph pair has 20 same versus 0 cross. Both are locally decisive, yet attaching
both ends recreates the known wrong whole-block join. The retained rule does
not make that second attachment. The 19-Mb wrong-join control remains split in
the full chr20 output.

The graph-arm window panel also newly spans 6,578,161–6,582,248,
12,269,536–12,277,077 and 48,225,787–48,229,445 bp. Each interval has
molecules physically crossing it; the panel observed no parental flank switch
and read concordances of 100%, 99.06% and 100% respectively. Only those three
span expectations and the graph-arm total were changed. In the full chr20
read-group audit, four new PS labels contain at least two prior PS groups with
at least ten truth-labeled reads each; none combines two old groups with
opposed parental orientation and >=90% purity.

The remaining 0.528 and 17.62-Mb gaps need an allele-observation audit before
another stitch policy change. Their candidate profiles contain many `-1`
indel calls even on MAPQ 60 spanning alignments. If those calls can be
recovered from the original CIGAR without realignment, the same statistical
link test can revisit them. Otherwise they should remain independent phase
sets rather than receive a forced join.

## Matched HiPhase realignment control

The frozen HiPhase 1.6.0-ac3f399 binary was run twice on each identical cropped
pgphase VCF and the same BAM/reference. Both runs used four threads and
`--ignore-read-groups`; the only changed option was
`--disable-global-realignment`. The default still has local realignment. The
exact boundary records were checked for a shared PS:

| Gap (Mb) | Default HiPhase | Global realignment disabled |
|---|---|---|
| 0.528827–0.528828 | joined | split (left SNP unphased) |
| 0.542052–0.545002 | joined | split (left deletion unphased) |
| 17.616778–17.625527 | joined | split (right insertion unphased) |
| 37.461999–37.466820 | joined | joined |

This controlled comparison shows that the first three **HiPhase** joins depend
on its global realignment setting for these cropped inputs. It does not prove
that pgphase must realign to solve them: a corrected CIGAR-to-allele projection
or another independent source of allele calls might recover the same signal.
Under the current no-realignment requirement, simply relaxing the stitch
threshold cannot reproduce those three joins from pgphase's observed matrix.
The 37.46-Mb case is a genuine stitch/transfer target that remains solvable
without HiPhase global realignment.

## Exact two-site HiPhase mechanism at the two remaining gaps

To isolate the reason for HiPhase's joins, its temporary diagnostic build
printed the per-read alleles handed to A* at the two exact boundary variants.
For each gap, the same two-site VCF and original BAM were run twice; only
`--disable-global-realignment` changed. The pgphase repository and production
HiPhase binary were not modified. The two-site result matches the larger
regional result, so no intermediate variant is required for either HiPhase
join.

| Gap | HiPhase global paired calls | HiPhase local paired calls | A* global | A* local |
|---|---|---|---|---|
| 528,827–528,828 | 53: 39 alt/ref, 14 ref/alt | 28: 27 ref/alt, 1 ref/ref | 2 phased hets, one PS | left SNP treated as homozygous; insertion alone phased |
| 17,616,778–17,625,527 | 32: 25 cross, 7 same | 32: 14 cross, 18 same | 2 phased hets, one PS | right insertion treated as homozygous; deletion alone phased |

At 0.528 Mb, the global SNP and insertion allele calls are each 100% concordant
with parental truth among the truth-labelled reads where HiPhase calls them.
Global graph alignment recovers both-site alleles on 53 reads, whereas local
alignment produces only 28 paired calls and calls the left SNP reference on
27 of them. Among those same 53 globally paired read names, pgphase has only
eight callable allele pairs at its boundary candidates; its insertion is
ambiguous on 39 and absent on one. Of the 13 reads with a callable insertion,
five have an ambiguous left SNP call.

At 17.62 Mb, local and global HiPhase each produce 32 paired calls, but their
**allele assignments** differ. Global alignment changes a nearly balanced
18:14 same/cross local relation into 7:25 same/cross. In the two-site
control, global deletion and insertion calls have 84.4% and 87.5% per-site
parental concordance, versus 71.9% and 65.6% locally. The boundary alleles
remain noisy, but the extra polarity lets A* keep both as phased
heterozygotes. Pgphase has only five callable pairs on those 32 globally paired
read names because the left deletion is ambiguous or absent on most of them. Its direct pair vote is 2 same / 3 cross.

HiPhase's source path explains the difference: `load_full_read_segments` calls
`global_realignment`, which constructs a reference/variant graph for each
aligned read span, aligns the read to that graph, and takes allele assignments
from traversed nodes before A* phasing. `--disable-global-realignment` calls
`local_realignment` independently at the sites. The phase-set builder connects
only sites that A* kept heterozygous and that have a valid co-observing read.
Thus the missing pgphase joins are **allele assignment/representation** issues,
not a failure to search a longer chain or a MAPQ filter. A stitch threshold
change cannot recreate the 45 missing paired observations at 0.528 Mb or the
changed 17.62-Mb allele polarization.


## Why the boundary calls are missing (original BAM CIGAR audit)

The paired names above were inspected directly in the original BAM, without
changing the pipeline. At 0.528 Mb, all 39 reads HiPhase globally calls
SNP-alt/insertion-ref have a **9–11 bp CIGAR deletion covering the SNP** at
528,827, so the BAM has no aligned base from which pgphase can call that SNP.
The other 14 reads have an aligned A at the SNP and an insertion at 528,829,
but its CIGAR length varies from 34 to 39 bp; the emitted insertion allele is
40 bp. None of those 14 alignments describes that emitted allele exactly.
The graph recovery matrix additionally represents the left boundary as an
8-base graph candidate at 528,820, while the emitted VCF has a SNP at 528,827;
a graph-walk allele at the broad site cannot simply be treated as an observed
base at the SNP. These are overlapping repeat-haplotype descriptions, not
independent exact observations of both VCF rows.

At 17.62 Mb, the boundaries lie in poly-A and poly-T tracts. Of the 32 reads
HiPhase globally calls at both sites, only 12 have the exact 1-bp deletion at
17,616,779 in the CIGAR; another eight globally called deletion-alt reads have
an aligned A there. At the right insertion, only two globally called
insertion-alt reads have the exact 1-bp T insertion at 17,625,528; most have an
aligned T without that exact insertion. Exact CIGAR projection would give 25
callable pairs, but its left deletion calls agree with parental labels on only
23/32 reads (71.9%), and its right insertion calls agree on only 14/25 (56.0%)
under the favorable haplotype orientation. HiPhase's global calls reach 27/32
(84.4%) and 28/32 (87.5%) respectively. Thus broadening the existing exact
CIGAR backfill to these clean MSA candidates would add weak or misleading
edges; the current backfill only admits noisy MSA candidates and would not
solve these two joins even if that guard changed.

This original-CIGAR audit does not establish that pgphase needs a new
realigner: the BAM recovery pipeline already runs a WFA/abPOA noisy-region
MSA. The decisive question is which reads that MSA admits and whether its
observations survive graph transfer; see the following controlled audit.


## Recovery MSA read-selection control

The standard BAM subsolve inherits `use_longcalld_bam_options`, which sets
`add_unplaced_msa_observations=false`. In a haplotype-aware noisy-region MSA,
`wfa_collect_noisy_aln_str_with_ps_hap` builds each consensus only from reads
already assigned to the **same selected phase set** and haplotype. With a null
`unassigned` pointer, it neither aligns other spanning reads to both
consensuses nor assigns their alleles at the MSA site. This is a source-MSA
observation loss, upstream of graph candidate transfer. The option is required
for the standalone BAM/longcalld parity mode; the following trial changed it
only in the graph recovery subsolve, then restored the production source.

The same narrow graph regions were run twice with identical inputs and only
that recovery option changed. Comparing the exact HiPhase-global paired read
names to the **BAM source** matrix gives:

| Gap | Standard source MSA paired calls | Admit unplaced source MSA paired calls | Trial source pair pattern |
|---|---:|---:|---|
| 0.528 Mb, deletion 528,826 / insertion 528,829 | 13/53 | 53/53 | 39 deletion-alt/insertion-ref, 14 deletion-ref/insertion-alt |
| 17.62 Mb, deletion 17,616,779 / insertion 17,625,528 | 5/32 | 32/32 | 26 cross, 6 same |

The trial's 17.62-Mb MSA pair pattern closely matches HiPhase global's 25
cross / 7 same. Across all truth-labeled source MSA reads in that narrow
region, the deletion agrees with parental read labels on 52/55 and the
insertion on 65/76. At 0.528 Mb, the trial's deletion and insertion each agree
on 70/72 source reads. Thus pgphase already has a realignment method capable
of extracting the missing allele information. The default selection rule
withholds many informative reads because they lack the particular prior PS
chosen for the MSA.

The 1-Mb controls confirm the downstream limits. At 17.62 Mb, the trial joins
the boundary VCF records (three phase sets in the region versus six baseline),
and the merged block's read labels retain the correct parental orientation.
At 0.528 Mb, the trial still splits the two VCF records even though its BAM
source MSA has all 53 paired calls. The graph left boundary is a broad 8-base
candidate at 528,820 rather than the BAM source's 9-bp deletion at 528,826;
its observations on the 53 reads do not polarize the BAM haplotypes cleanly.
The baseline source had a previously audited weak cut at 528,850. With
unplaced reads admitted, every local source edge has support; the remaining
blocker in that trial is transfer and stitching, not that cut. In a smaller 117-kb control,
the standard run joins the 0.528-Mb pair but the blanket unplaced-read trial
splits it. This is why enabling the option everywhere is not a validated
production change.

A robust implementation would keep standalone BAM parity unchanged and use
unplaced-read MSA evidence only for unresolved graph seams, without allowing
that evidence to silently rewrite already-correct graph block orientation.
At 0.528 Mb it must additionally reconcile the graph's broad left allele with
the BAM deletion/insertion haplotype and test the source weak cut using
independent spanning reads. The production join must be tested against both
flanks and independent read halves; parental labels remain evaluation-only.

## Follow-up: separate MSA admission, transfer and stitching

A same-command 1-Mb A/B used chr20:1–1,000,000, 17,000,000–18,000,000,
and 19,000,000–20,000,000 with one worker and the same BAM/GAF/catalog.
The parental read score counts each truth-labeled read once within its output
phase set, choosing the better whole-block orientation. These are local
controls, not whole-chromosome accuracy estimates.

| Trial | 0.528-Mb region tagged / discordant | 17.62-Mb region tagged / discordant | 19-Mb negative control tagged / discordant |
|---|---:|---:|---:|
| Existing recovery | 3,699 / 35 | 4,071 / 39 | 4,071 / 51 |
| Admit unplaced reads throughout targeted BAM MSA | 3,784 / 53 | 4,124 / 606 | 4,083 / 53 |
| Same admission plus a graph–BAM–graph transitive join when the direct outer vote is absent | 3,781 / 33 | 4,118 / 604 | 4,083 / 53 |
| Baseline BAM solve plus a second MSA pass used only to fill missing exact-site observations | 3,710 / 37 | 4,081 / 47 | 4,105 / 64 |

The blanket MSA trial supplies the 53/53 and 32/32 paired source calls, but
also rephases other source sites. The 17-Mb region then has about 600
discordant reads; a correct local MSA pair does not prove that the surrounding
phase block is correct. The observation-only second pass keeps original
source genotypes and phase sets, but neither target gap closes: the extra
calls do not by themselves turn two original BAM blocks into one block.
It also doubles the targeted subsolve. That experiment was removed.

At 0.528 Mb, a single imported BAM block can join both graph flanks in the
stitch transaction, but the transaction rolls it back because there is no
direct outer graph-flank allele vote. Allowing that transitive join on source
path and graph/BAM gauge votes closes this local seam, yet the same rule fails
the existing weak-cut regression near 8.64 Mb and does not cure the 17-Mb
read switch. An added end-to-end molecule check protects those focused
regressions, but with the observation-only input it cannot close the targets.
The join rule was restored.

A separate transfer check found that the 9-bp BAM deletion at 528,826 spans
the seam's 528,827 left boundary, but transfer tests candidate *start*
coordinates strictly inside the seam. It therefore omits this source row.
Naively admitting all boundary-spanning deletions also transfers long
repeat-region events near 0.864 Mb and provides no read-label gain in the
same-command A/B. That trial was restored. The catalog's 528,827 A>T row
comes from an eight-base multiallelic graph snarl and lies inside the BAM
deletion; graph calls on the 53 paired reads cover almost exclusively the
deletion-supporting read group. Thus a two-site HiPhase SNP/insertion phase
cannot alone establish the physical haplotype of this overlapping
SNP/deletion/insertion representation.

Focused baseline tests for the 3.85-Mb and 4.76-Mb windows and the supported
graph-flank attachment pass. Blanket MSA admission fails all three; the
transitive join also fails the 8.64-Mb weak-cut control. All experimental
production edits from this follow-up were removed. The next implementation
needs to build an independent BAM gap block from the added observations,
validate its internal path, then test its two graph attachments without
changing established source or graph flanks.

## Source-block label consistency bug and repair

A same-run matrix dump around 17-18 Mb localized the 606-discordant-read
trial regression to the in-chunk stitch, before BAM source attachment. The
original graph blocks at PS 17679082 and 17883198 were each truth-pure
(533/0 and 540/0 tagged/discordant). Recovery merged them through two
different BAM phase blocks. Their boundary allele profiles produced a
significant source-block link, but among spanning reads already tagged to
source PS 17774311, 25 inferred the opposite HP from its boundary alleles
and only 10 agreed. The former stitcher accepted that locally inconsistent
edge and merged the graph blocks into PS 17667022 with 587 discordant reads.

The imported/imported edge now collects one agreement or conflict per
spanning read against its own saved source HP. A one-sided exact binomial
conflict at `p <= 0.01` vetoes the edge and restores the atomic transaction.
In the admitted-read trial, final BAM output scores 4,133/43 for the 17-18
Mb control, versus 4,124/606 before. The normal recovery setting remains
4,071/39. The useful 4.76-Mb join remains present, and the 3.85-Mb floor
passes. Blanket MSA admission still fails those two window expectations and
remains disabled. A synthetic source-label regression failed on the original
stitcher and passed after this change.

## Scoped unplaced-MSA admission and the remaining representation gap

The ordinary BAM source solve at 17,616,778–17,625,527 has five callable
indel pairs among 33 MAPQ-30 reads physically crossing the two source blocks.
After exact-CIGAR backfill, a one-sided 50:50 missing-call test gives
`p=3.309e-05`. Re-solving that targeted group with unplaced reads admitted
to its existing MSA clusters makes both VCF boundary rows share PS 17559224
in the isolated replay; the dedicated integration test also finds no parental
flank switch and at least 80% local read concordance. In the same-command
17–18 Mb replay, this mode tags 4,133 truth-evaluable reads with 43 discordant,
versus 4,071/39 under ordinary recovery.

The trigger is evaluated after exact-CIGAR backfill. Evaluating before it
incorrectly selected 4.76 Mb: zero source pairs were initially callable, but
19 of 36 were callable after backfill; the needless MSA re-solve lost its
useful graph join. The 3.85-Mb seam has three source blocks and remains on the
ordinary solve. Requiring 20 MAPQ-30 crossing reads avoids a seven-read
14.66-Mb trigger that disrupted the established 14.58-Mb positive attachment.
A 15.096-Mb negative control had 66 crossing reads and only 12 callable
pairs, but both boundary coordinates had co-located alternatives. Re-solving
that multi-allelic pair made a wrong graph join. The retained trigger requires
a single candidate row at each boundary, so that control stays split.

The 0.528-Mb gap is a distinct representation problem. At its graph SNP
position, the original BAM has 27 maternal A bases, 42 paternal deletions,
one paternal A and no T base. Graph SNP observations show 24 T calls and 19
A calls, all on paternal reads; maternal graph calls are missing. The BAM
source has a deletion and insertion phase path spanning the position, but it
does not contain a source SNP row equivalent to the graph A>T representation.
Blanket MSA admission increases paired source calls without providing a
valid graph-SNP-to-source allele mapping. The graph SNP remains a separate
phase set until this overlapping representation is reconciled.

## Statistical source/graph approval and the 0.528-Mb allele trap

The 528,827–528,828 coordinate split exposed a zero-error approval rule:
29 exact shared-site candidate checks agree with the source gauge, and its
2x2 read vote has 173 agreeing reads versus one conflict, with both source
haplotypes represented. A one-sided exact binomial `p <= 0.01` vote admits
that graph flank. In a trial applying the vote to weak-cut source blocks,
pgphase joined all 35 tracked short gaps and improved full-chr20 read purity.
That apparent last closure was **wrong at the variant level**: pgphase emitted
`1|0` at both the 528,827 A>T SNP and the 528,828 insertion, while HiPhase
and its 53 phase-consistent globally called read pairs support opposite ALT
haplotypes (`0|1`, `1|0`). The graph SNP has 24 ALT and 19 REF GAF calls, but
24 ALT and 18 REF callers have a CIGAR deletion across the SNP base; only one
REF caller has a physical A. Neither allele group separates the BAM source
haplotypes. The source block also has a weak internal cut. Thus a strong
source-to-graph vote elsewhere in the block cannot validate this particular
SNP/insertion allele connection.

The retained statistical approval applies only when the original BAM source
path has no weak cut. Weak-cut local runs keep their previous zero-conflict
approval. A regression now rejects a shared PS if the two boundary ALT alleles
are on the same haplotype, and retains the 19,395,544–19,414,720 wrong-join
control. The 17,616,778–17,625,527 indel pair from scoped MSA admission and
the 542,052–545,002 join remain present.

| Matched full chr20 run | Phased truth reads | Correct | Discordant | Purity | VCF blocks | N50 | Correctly oriented tracked joins |
|---|---:|---:|---:|---:|---:|---:|---:|
| Scoped MSA and source-label veto | 236,827 | 228,847 | 7,980 | 96.6305% | 395 | 491,812 bp | 34/35 |
| Statistical approval restricted to complete source paths | 236,834 | 228,912 | 7,922 | 96.6550% | 392 | 491,812 bp | 34/35 |
| Rejected weak-cut statistical trial | 236,851 | 228,925 | 7,926 | 96.6536% | 393 | 506,103 bp | 34/35; 35/35 share a PS but the last has wrong allele phase |

The fixed inventory is the 35 noncentromeric, <10-kb intervals in
`evaluations/2026-09-24-hiphase-on-pgphase-sites/short_gap_replay.tsv` with
at least 20 truth reads and at least 98% HiPhase local purity. All 34
retained joins have the same relative biallelic genotype as HiPhase. The
full-run read-group audit found one new phase set combining two prior
truth-pure groups of at least ten reads; those groups have the same parental
orientation. `make unit-tests` passed, and the full `make window-tests`
panel passed 432 assertions in 16 cases with the added 528-kb
allele-orientation assertion.

A correct last join needs a sequence-aware translation between the broad
multiallelic graph snarl and the overlapping BAM deletion/insertion source
haplotypes. Current graph SNP observations do not identify both alleles on
both source haplotypes; a phase-set label or whole-block read gauge alone
cannot supply the missing mapping.

## Exact overlap projection closes the last tracked gap

The 528,827 graph A>T row is a terminal child-snarl site. Its REF and ALT
GAF observations mostly describe reads whose original BAM CIGAR deletes the
reference base. The BAM source contains a phased 9-bp deletion beginning at
528,826 and a 40-bp insertion at 528,829. Fifteen source reads call both rows:
eight deletion-ALT/insertion-REF and seven deletion-REF/insertion-ALT, with no
opposing pair. The source deletion call agrees with the physical CIGAR on all
33 callable source reads (25 reference, eight deletion).

Recovery now checks the original BAM CIGAR only at a terminal graph SNP inside
a phased source deletion. At this site, the 24 graph ALT calls are all on
physical deletion reads; the 18 graph REF calls on deletions and 28 physical
reference calls (27 lacking a graph child call) give Fisher two-sided
`p=2.04e-7` for deletion enrichment of graph ALT. No physical T occurs, and
no graph ALT occurs on a physical A read. Among source-phased reads, reference
versus deletion associates with the two source haplotypes at Fisher two-sided
`p=0.00976`. The source path from deletion to its next phased site has no weak
cut. Those independent checks map graph ALT to the BAM deletion haplotype and
move only this terminal SNP into that local BAM block. A one-site graph block's
read labels follow exact physical calls; a larger graph block keeps its other
read labels. The ordinary stitch and source transfer then set the final HP
orientation. The correction does not run a new sequence alignment.

The chr20:1-1,000,000 regression now requires the A>T SNP and adjacent
insertion to share a PS with opposite ALT haplotypes; the 19-Mb wrong-join
control remains split. The 100-kb and 1-Mb replays both meet that genotype
condition. The full chromosome run with the same inputs and eight workers
completed in 2:10.73 (peak RSS 20.5 GB). Its truth-evaluable read score and
VCF block statistics compare with the immediately preceding guarded run:

| Run | Phased truth reads | Correct | Discordant | Purity | VCF blocks | N50 | Tracked short gaps joined |
|---|---:|---:|---:|---:|---:|---:|---:|
| Guarded source transfer | 236,834 | 228,912 | 7,922 | 96.6550% | 392 | 491,812 bp | 34/35 |
| Terminal SNP overlap projection | 236,834 | 228,914 | 7,920 | 96.6559% | 392 | 491,812 bp | 35/35 |

The 35-gap inventory uses the committed short-gap replay with at least 20
truth reads, at least 98% HiPhase local purity, and the centromeric region
excluded. A coordinate counts as joined when any phased variant row at each
boundary shares a PS; this handles the two deletion alternatives at 37.462 Mb.
The newly joined 528,827/528,828 pair emits `0|1` and `1|0` in the full run.

After indexing each source block's phased coordinates once instead of scanning
all source candidates per deletion, the final eight-worker full chr20 rerun
produced a byte-identical phased VCF and phased BAM (`cmp` exit 0 for both)
to the scored run above. `make unit-tests` and `make window-tests` pass; the
window panel reports 433 assertions in 16 cases. `git diff --check` passes.

The before/after full VCF differs in exactly one record: 528,827 changes
`1|0:PS=130541` to `0|1:PS=487305`; the 528,828 insertion remains
`1|0:PS=487305`. Only two primary-read HP/PS tag pairs change in the full
BAM. This localizes the two-read truth improvement and shows the other 34
previously validated gap genotypes are unchanged.
