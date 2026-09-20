# Root cause of the largest remaining clean-region class: surplus reference reads turn a het into a tie, and a tie resolves to reference

## Where the remaining differences actually are

After `refresh_msa_observations` was turned off
(`evaluations/2026-09-20-msa-refresh-breaks-complementarity/`), 1,068 record
mismatches remained against longcallD on chr20. They are not spread evenly.
Binned at 100 kb and grouped by local mean candidate depth:

| local mean depth | bins | candidates | mismatches | rate |
|---|---:|---:|---:|---:|
| DP >= 40 | 613 | 106,828 | 283 | **0.26%** |
| DP 20-40 | 9 | 3,265 | 53 | 1.62% |
| DP < 20 | 26 | 10,036 | 732 | **7.29%** |

69% of all remaining mismatches sit in 26 low-depth bins -- chr20's
centromere/satellite band, where mean candidate depth is 13 against 59
chromosome-wide and both tools are working on degenerate input (in the worst
single bin, 26.6-26.7 Mb, we emit 802 records, upstream emits 689, and they
agree on 548). In the 613 normal-coverage bins, holding 89% of all
candidates, record identity is **99.74%**.

Outside the band, 296 mismatches remain, and they are one class: **218
records longcallD emits that we discover as candidates but leave at
`hap_alt = hap_ref = 0`**, so they never reach the VCF. All 218 are indels
(142 DEL, 76 INS), all in homopolymer or short-tandem-repeat context, at full
depth (mean DP 62) and clean heterozygous allele fraction (mean AF 0.32).

## Ruling out the obvious explanations

Diffed against upstream and found faithful, so none of these is the cause:

- `select_init_var` -- identical, `is_homopolymer_indel == 0` exclusion included.
- `init_assign_read_hap_based_on_cons_alle` -- identical, including the
  `is_homopolymer_indel == 1 || NOISY_CAND_HOM` skip (assign_hap.c:166).
- `update_var_hap_profile_based_on_read_hap` -- identical, including the
  `hap == 0` double-increment of both haplotype profiles.
- `update_var_hap_to_cons_alle` -- identical: `> max_cov` from `max_cov = 0`,
  so a tie resolves to the lowest allele index. Upstream's own comment on that
  line reads `// prefer ref allele`.
- Phase 1's outward sweep, its `hap = -1 -> hap = 1` seeding, and its HOM-var
  skip -- identical.
- `add_msa_site_observations` -- tested off at chr20:3,870,827 and produced
  byte-identical profiles, which ruled it out *at that locus* (wrongly, for
  the class -- see below).

`joint_het_orientation`, the existing mechanism that would force these sites
to het, was re-measured on the current baseline rather than trusted from
`evaluations/2026-09-19-joint-orientation-alignment-scope/` (whose numbers
predate this session's three fixes). It is still net-negative, now clearly:
it recovers 299 of the upstream-only records but invents 792 new ours-only
ones (identity 99.60% -> 98.94%) and drops read accuracy to 96.55%, below
longcallD's own 97.03%. Not shipped.

## The trace

chr20:3,997,065 (`CA>C`), where our allele counts and upstream's are
effectively identical (ours 67 ref / 22 alt, upstream 65 / 22), so everything
before genotyping agrees.

Both tools were run on the same 500 kb chunk with `-V 2` and their traces
parsed into the same table. longcallD prints `Hap<N>-Read:` with the **MSA
cluster** index (collect_var.c:2336) and `read: <q> hap: <N>` with the
**chunk-level haplotype** (assign_hap.c:276); pgphase prints the same two
under `-V 2`.

MSA cluster level -- our clustering is *cleaner* than upstream's:

| | cluster 1 | cluster 2 | cluster-2 argmax |
|---|---|---|---|
| ours | 45 reads, all ref (truth: 43 PAT) | 42 reads, 20 ref / **22 alt** (truth: **42 MAT, 0 PAT**) | **alt** |
| upstream | 45 reads, all ref (truth: 43 PAT) | 91 reads, 71 ref / 20 alt (truth: 77 MAT, 14 PAT) | ref |

Chunk level, which is what the genotype is read from:

| | hap1 | hap2 | argmax hap2 | emitted |
|---|---|---|---|---|
| longcallD | {ref: 45} | {ref: **20**, alt: 22} | alt, by 2 | `0\|1` |
| ours | {ref: 45} | {ref: **22**, alt: 22} | **exact tie -> ref** | dropped |

Everything else about the two runs is the same. Per-read alleles against the
parental-origin truth map are identical in both tools (MAT 22 ref / 22 alt,
PAT 43-45 ref / 0 alt). The read partition is the same and, in both tools,
**no read is left unassigned** (ours MAT hap1=2 hap2=42 hap0=0, PAT hap1=43
hap2=2 hap0=0; upstream the same but PAT hap2=0).

The whole difference is two reads. We admit two paternal reference reads that
upstream does not (DP 89 vs 87), the solve puts them in hap2 -- the maternal
haplotype -- and hap2's profile goes from 20 ref / 22 alt to 22 / 22. The
argmax tie-break then prefers reference, both haplotypes read reference,
`hap_alt = hap_ref = 0`, and the record is dropped.

## The cause generalises across the whole class

Comparing our depth with upstream's at all 218 clean-region dropped records:

| | value |
|---|---:|
| mean ours DP | 62.0 |
| mean upstream DP | 41.3 |
| ours deeper / equal / **shallower** | 141 / 77 / **0** |
| mean DP excess | +20.7 (median +16) |
| mean ALT excess | +6.7 |

We are deeper than upstream at every one of the 218 and shallower at none,
and the surplus is overwhelmingly reference (+20.7 depth against only +6.7
alt, so roughly 14 extra reference reads each). That is exactly the input
that turns an alt-majority haplotype profile into a tie or a ref majority.

The source is `add_msa_site_observations`, which adds observations for reads
the MSA could not place. It has no counterpart in longcallD: a noisy
candidate's depth upstream is exactly the reads its two cluster alignments
cover. The earlier single-locus test that appeared to exonerate it was run at
chr20:3,870,827, where both clusters already covered every read, so it had
nothing to add there -- a locus where the mechanism is inert, not a locus
that disproves it.

## Whole chr20, with it off

| | our records | identical | ours-only | upstream-only | total mismatch |
|---|---:|---:|---:|---:|---:|
| previous | 118,152 | 117,677 (99.60%) | 475 | 593 | 1,068 |
| **this fix** | 118,311 | 117,842 (99.60%) | **469** | **428** | **897** |

Read accuracy against the diplinator truth BAM improves as well, 98.34% ->
**98.45%** (3,658 -> 3,343 discordant), and phase-set count moves toward
upstream's 429 (324 -> 380). longcallD's own score on the same truth BAM is
97.03%.

chr20:3,997,065 now emits
`CA>C 0|1:87:65,22:0.253:60:3981464` -- byte-identical to longcallD's own
record, phase set included.

## Status: shipped

`opts.add_unplaced_msa_observations = false`, scoped to
`collect_bam_variation()` (`src/collect_pipeline.cpp`), alongside
`anchored_stage2`, `merge_colocated_msa_alleles` and
`refresh_msa_observations`. Struct default left `true`; the graph arm never
sets the field. `make window-tests` 125/125 (graph-arm windows included),
unit ALL PASS, predicate 151/151.

## What is left

897 mismatches. Roughly 700 of them are inside the low-depth centromeric
band, where the 7.29% mismatch rate reflects input both tools struggle on
rather than a divergence in the port. The clean-region remainder is now small
enough that the next step is to re-derive the depth comparison above on the
surviving records and see whether the residue is still a depth-surplus class
or something new.

Co-authored-by: Claude Opus 5 (1M context) <noreply@anthropic.com>
