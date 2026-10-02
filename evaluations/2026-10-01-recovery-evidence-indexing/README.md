# Recovery evidence indexing and corroboration (2026-10-01)

## Retained fixes

Two defects were reproduced independently of parental truth:

1. `update_read_var_profile_with_allele` grew primary alleles, graph/BAM calls,
   and query indices, but omitted `bam_base_qualities`. Prepending two sites
   left `{40, 5}` at offsets zero and one instead of `{0, 0, 40, 5}`. The Q40
   observation then belonged to the wrong candidate, and the original SNP
   could lose its quality-bearing vote. Appending sites left a short quality
   channel. The updater now pads each populated quality channel with zero at
   the same offsets as the alleles. A post-transfer invariant checks its length.
2. `corroborated_bam_block_flip` retained only the last visited clean SNP and
   MSA indel positions. With agreeing SNPs at 1,000 and 1,950 and an indel
   anchored at 1,999, it rejected corroboration because the last pair was only
   49 bp apart, despite the independent earlier SNP. It now keeps the extrema
   of all agreeing SNP/indel positions and checks the largest pair distance.
   This is linear in the observed sites with constant additional storage;
   the MAPQ, base-quality, 100 bp spacing, complete-source and conflict gates
   remain unchanged.

The existing updater and corroboration predicate have header declarations so
unit regressions exercise the actual implementations. Tests cover prepending,
appending, both-direction growth, absent qualities, preserved existing quality,
malformed quality channels, reversed visit order, corroboration on either
block, low-quality distant SNPs, conflicting SNPs, earlier indels, the exact
VCF-anchor spacing boundary, and opposite block parity. The old implementations
failed three assertions each before their fixes. No truth, competitor output,
or chr20 coordinate enters either implementation.

## Matched chromosome result

The comparison baseline includes the preceding accepted single-base insertion
repair (the uncommitted state after commit `74bc786`). These results must not be
compared against `74bc786` alone as if both fixes newly closed its four gaps.
Baseline output: `/tmp/pgphase-gap-next/fixed-full/`; final output:
`/tmp/pgphase-gap-next2/corroboration-full/`.

Both use the same annotated chr20 BAM, stripped graph catalog, coordinate GAF,
CHM13 reference, default recovery and eight threads. VCF files are byte-identical.
All 256,570 primary output read names and HP/PS pairs are identical, including
unphased reads. The retained defects do not trigger an output change on this
chr20 fixture.

| Measure | Before | After |
| --- | ---: | ---: |
| Truth-scored phased reads | 237,070 | 237,070 |
| Truth-correct assignments | 229,890 | 229,890 |
| Discordant assignments | 7,180 | 7,180 |
| Read concordance | 96.9714% | 96.9714% |
| Read phase sets | 692 | 692 |
| VCF variant keys | 62,368 | 62,368 |
| VCF blocks | 341 | 341 |
| Span N50 | 672,998 bp | 672,998 bp |

No new gap closes; the 13 tracked HiPhase-positive gaps remain open. Every
shared genotype and phase-set label is unchanged. Owning-chunk replays at
3, 10, 15, 19, 21, 23, 35, 36, 41, 50, 57 and 58 Mb likewise retain the
accepted coverage, read scores and VCF fields. Existing panel expectations
are not refreshed for these fixes.

## Rejected representation trials

A high-quality, shifted one-base deletion ALT repair before source phasing
closed no remaining gap. At 58–59 Mb it added 14 scored reads but only five
correct ones (discordant 57 to 66). It was discarded.

A post-solve veto against REF backfill when a nearby indel gives an identical
edited reference reduced whole-chromosome discordant assignments from 7,180
to 7,102 and added 122 phased reads, but reopened established connections at
1.180, 4.778, 24.581, 48.929 and 56.064 Mb; N50 fell to 643,699 bp. The
4.766 Mb window lost its required span. It was discarded. This confirms the
source-consistency pitfall already recorded in the September 30 shifted
insertion evaluation: correcting observations after solving can invalidate
the transfer/stitch relation based on the prior matrix.

An isolated pre-phase multi-base deletion ALT trial used only non-homopolymer
MSA-verified candidates and Q30 sequence-equivalent CIGAR events. At 10–11 Mb
it added three correct reads and removed one discordant read, but did not close
the target. At 36–37 Mb discordant assignments rose from 163 to 292, with
read phase sets dropping from 18 to 16. The production source and binary were
never changed by that isolated trial. It was discarded.

## Reproduction and validation

Use `../2026-10-01-single-base-msa-recovery/commands.sh` with binaries built
before and after these fixes for a matched full chr20 run. The scorer orients
each read PS independently using `test_data/derived/chr20_truth_hap.tsv`.
Truth is evaluation-only. Temporary replay, scoring scripts and logs are in
`/tmp/pgphase-gap-next2/`.

- Build, shared units, port parity and upstream parity pass.
- Predicate suite passes: 243 assertions in 24 cases, including 30 additional
  assertions beyond the accepted single-base insertion state.
- Full window panel passes: 3,868 assertions in 48 Catch2 cases over the
  unchanged 89-window panel, with all established closure and protected-split
  expectations retained.
- `git diff --check` passes; the final build introduces no warnings.
