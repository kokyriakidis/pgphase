# Depth attrition inside the guided MSA: traced past every ported function, into stage 1

## The symptom

Of the ~963 single-ALT record mismatches against longcallD that survive
after the `anchored_stage2` fix (`evaluations/2026-09-19-anchored-stage2-
alignment-arm/`), several show DP roughly half of upstream's at the same
locus, with neighboring candidates in the same reconstructed region
reporting wildly inconsistent depths -- not a smooth falloff, a jumble.

Two spot-checks: chr20:16,627,890 (`T>A`) reports DP 41; longcallD reports
DP 69 at the same position. chr20:22,823,437 (`T>A`): ours DP 30,
upstream DP 61. Both loci: 69 (61) reads physically overlap the position
in the BAM at MAPQ >= 47, matching upstream's DP exactly -- so upstream
uses essentially every available read, and something in our path does
not.

## What was ruled out, each checked against longcallD's actual source

- **Region/read discovery.** `MsaFire Hap 16627747-16627927 ... 69 full`
  (our own verbose log) -- the noisy region finds and full-covers exactly
  69 reads, matching upstream's DP. Not a discovery gap.
- **Phase-set gating into the MSA.** Instrumented directly: all 69
  full-covering reads carry a valid `hap` (1 or 2) and the SAME
  `phase_set` (16,502,593) entering `wfa_collect_noisy_aln_str_with_ps_hap`.
  Zero excluded here.
- **abPOA clustering.** Every `assign_hap_based_on_germline_het_vars_kmeans`
  cluster call in this region logs `n_ps == out_clu_n` -- reads handed to
  abPOA come back with none dropped.
- **The full-cover similarity check itself
  (`is_match_aln_str`/`collect_phase_noisy.cpp:283`).** Byte-for-byte
  identical to longcallD's `collect_var.c:1960-1993`, same 0.9 threshold,
  same `n_eq >= len * cons_sim_thres` rule, same `cover_start`/`cover_end`
  bookkeeping.
- **The boundary/trim logic feeding that check.** For a BOTH_COVER read
  (all 69 here are), `wfa_trim_aln_str` returns immediately
  (`collect_var.c`-equivalent guard, `align.cpp:951`) -- `query_beg=0,
  query_end=aln_len-1` span the whole alignment, so the per-read boundary
  window cannot be excluding anyone.
- **WFA is not even in this path for clustered reads.** Their
  cons-vs-read alignment string comes straight from abPOA's own MSA
  columns (`make_cons_read_aln_str`), not a fresh WFA realignment; WFA is
  only used for the ref-vs-consensus row and for reads NOT placed in
  either cluster.

## Stage 1 ruled out directly -- correcting the initial hypothesis

The first hypothesis here was that stage 1 (clean-sites-only phasing)
produces a different read partition than longcallD's equivalent step,
matching `d61e98a`'s unresolved conclusion at a different locus ("the
read partition... needs to agree"). Tested directly rather than assumed:
every CLEAN candidate in chr20:16,600,000-16,660,000 (79 records) was
diffed field-for-field against longcallD's own output. DP, ref count,
alt count and AF are **identical at every single position**. GT is a
perfect mirror image throughout (0|1 wherever we say 0|1, upstream says
1|0, and vice versa, at every one of the 79 sites) -- the expected,
harmless global gauge flip, not a partition difference. Stage 1's read
partition matches longcallD's exactly. That hypothesis is wrong;
retracted.

## What remains, narrowed one level further

With stage 1 cleared, the loss is entirely inside stage 2's guided MSA
for this specific region. Re-instrumented precisely: `clu0=36, clu1=35`
(71 reads split into the two abPOA-discovered clusters, summing to
`info.n_reads`, zero lost at clustering -- matches the earlier
`n_ps == out_clu_n` finding). The two clusters together hold 71 reads;
the final candidate at this position counts 41. Every mechanism between
"clustered" and "counted" -- the 0.9 full-cover similarity check, its
boundary/trim logic, the sort feeding read order into abPOA -- is now
verified byte-for-byte identical to longcallD's source
(`collect_var.c:1960-1993`, `align.c:954-984`). `sort_noisy_region_reads`
in particular: identical bubble sort, identical three-key comparison
(full-cover, error-rate, length), identical "stable on ties" behavior on
both sides.

What is NOT verified, and is the last remaining candidate: the read
*insertion order* feeding that stable sort's tie-break, and by extension
the exact read order abPOA receives. If that order differs even slightly
between the two tools -- both sides read the same coordinate-sorted BAM,
so this would have to be some other ordering detail, not gross read
selection -- a clustering/consensus library that is sensitive to input
order could produce a different (but individually valid) abPOA consensus
per haplotype, which is sufficient to explain the observed depth loss
without any of the surrounding C++ being wrong.

## The read-order/abPOA angle, opened and closed

Chased the "exact read-insertion-order" candidate to ground. Found and
compared the actual origin of `ordered_read_ids` on both sides:

- **Chunk-level read order.** longcallD's `sort_chunk_reads`
  (`bam_utils.c:1636-1656`, `qsort` over `pos asc, end desc, NM asc,
  qname asc`) against our `populate_chunk_read_indexes`
  (`collect_var.cpp:1163-1173`, `std::sort` over the identical four keys
  in the identical directions). Both read `NM` the same way -- directly
  from the BAM's own `NM` aux tag, default 0
  (`bam_get_NM`/`bam_digar.cpp:1273`). Because `qname` is the final,
  almost-always-unique tiebreak, `qsort`'s and `std::sort`'s differing
  stability guarantees have nothing to act on.
- **Noisy-region read order and extraction.** longcallD's read-coverage
  classification and per-read sequence extraction
  (`align.c:1391-1429`) against ours (`align.cpp:225-263`) --
  line-for-line identical: same DIGAR walk, same DEL-boundary special
  case, same `LONGCALLD_NOISY_REG_FLANK_LEN` / `kNoisyRegFlankLen`
  constant (10 on both sides), same four-way cover classification.
- **The noisy region boundary itself**, confirmed exact match at this
  locus: both tools log `16627747-16627927` for this cluster (upstream's
  own `-V 2` output additionally logs several narrower, overlapping
  sub-region lines nested inside that same span -- a per-read-window
  artifact of its logging, not a different reconstruction boundary; the
  outer, merged region used for reconstruction is the same span both
  sides use).

Every mechanism from BAM read ingestion through to the exact bytes
handed to abPOA is now verified identical, function by function. There
is no remaining C++/C glue code between the two tools left to compare
for this symptom.

## Status: open -- past the limit of source-level comparison

No fix attempted, none is possible at this level. `d61e98a` stopped at
"the read partition needs to agree" without identifying which stage.
This record rules out, with direct evidence rather than inference,
every stage between BAM ingestion and abPOA's own internals: stage 1
(byte-identical clean-site phasing), noisy-region discovery and
boundary computation, per-read coverage classification and sequence
extraction, chunk- and region-level read ordering, phase-set gating,
clustering, and the full-cover similarity check with its boundary
logic. What's left is abPOA's own behavior given input this project's
code constructs identically on both sides -- resolvable only by
instrumenting or diffing abPOA's internal state directly, not by
comparing this repository against longcallD's.

Co-authored-by: Claude Sonnet 5 <noreply@anthropic.com>
