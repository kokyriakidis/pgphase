# Automatic, cumulative gap recovery

`collect-hybrid-variation --recover-gaps` now detects and attempts recovery of
remaining internal phasing gaps automatically. No private-site whitelist,
competitor output, external graph merger, or PS-distance cap is needed.

## Use

Add these options to the existing hybrid command:

```bash
--recover-gaps --gap-recovery-report recovery.tsv
```

For example, the default-linker smoke test used:

```bash
./pgphase collect-hybrid-variation \
  --ref test_data/chm13v2.0.chr20.renamed.fa \
  --bam test_data/HG002_chr20_hifi_mapped_to_CHM13_chr20_annotated.bam \
  --graph-sites test_data/chr20.sites.striped.vcf.gz \
  --gaf test_data/HG002.chr20.annotated.coord.gaf.gz \
  -r 'CHM13#0#chr20:14950000-15180000' \
  --recover-gaps --gap-recovery-report recovery.tsv \
  --min-read-margin 2 -o candidates.tsv \
  --phased-vcf-out recovered.vcf -b recovered.bam
```

MSA admission retains its default margin of 24. Existing `--stitch-rule` and
`--stitch-min-margin` control endpoint acceptance; recovery uses the same vote
function as ordinary chunk stitching. The feature is opt-in and belongs to
the hybrid command, which has both BAM and GAF evidence. Its preserved scaffold
is that command's initial clean BAM+graph pass, not an externally supplied
graph-only BAM.

Manual private-site/tiering flags and `--gap-fill` cannot be combined with
automatic recovery, because they define different admission/output policies.
`--private-msa-margin` remains available, but neither regression needs it
lowered. Without `--recover-gaps`, the prior pipeline behavior is unchanged.

## Workflow and acceptance

1. Run clean hybrid phasing and normal inter-chunk stitching, then apply the
   configured read-margin and minimum-PS-size gates.
2. Detect internal coordinate gaps between distinct, read-supported phase
   blocks. Detection spans adjacent chunks, but never crosses an excluded
   interval between separate user-requested regions.
3. Reload each gap plus up to 50 kb of each flank, clipped to the contiguous
   requested interval. Try clean BAM+graph candidates first with the existing
   k-means core.
4. If the flanks are not connected, admit MSA SNPs and rerun the same core.
5. If still unresolved, retain the SNPs and their observations, admit MSA
   indels, and rerun it again.
6. After each tier, filter proposal reads and apply the ordinary overlap-read
   stitching vote rule independently at each endpoint. Stop once one proposal
   PS has accepted links to both endpoints. Otherwise preserve supported
   one-sided additions and continue on the remaining break.

MSA covers the remaining gap plus 5 kb near its boundaries. Existing noisy
intervals are retained intact. Stretches without detected noise are covered
by 512 bp MSA windows with 64 bp overlap, so lack of a noise trigger no longer
prevents examining a gap. Existing length, coverage, consensus and read-score
gates still apply; an unsuccessful MSA attempt is not evidence for a join.

The same local `PhasingChunk` survives the three tiers. Indels cannot enter
the SNP tier, and SNP observations are not discarded when escalating. A
complete two-sided proposal is preferred over independently stronger matches
to two disconnected proposal blocks. An isolated new block does not count as
a closed gap.

Initial read assignments are never independently overwritten. An accepted
join only flips/relabels a whole original PS uniformly. Unphased reads may be
added through an accepted endpoint, and new candidate orientations/profiles
are transferred together. Unsupported exact candidate matches can be
replaced, including previously filtered rows; trusted core candidates retain
their orientation. Raw observations from untagged proposal reads survive the
transfer without creating unsupported read tags.

The normal initial `stitch_chunk_haps` still runs. Recovery adds an adapter for
matching the reloaded window to specific existing phase sets and reuses the
factored normal overlap vote rule; it does not blindly apply chunk-wide
relabeling to all blocks in the recovery window. Reference/BAM/GAF handles are
owned by the recovery context. Recovery is sequential after initial parallel
chunk workers finish, so overlapping gap proposals do not race.

## Reproduce the regression

```bash
make -j"$(nproc)"
make unit-tests
bash evaluations/2026-09-13-auto-gap-recovery/run.sh
```

Environment overrides: `PGPHASE`, `OUT` (default
`/tmp/pgphase-auto-gap-regression`), and `DATA_ROOT` (default
`~/Downloads/pgphase-eval-data`). The runner writes exact commands, a binary
hash, BAM/VCF/candidate outputs, tier reports, read-truth results, and block
preservation assertions. It never uses truth to choose gaps or admit sites.

The main comparison retains `--link-by-alleles --block-link-window 8
--min-read-margin 2` from the earlier audits. Two worker threads are used.
The 15m split case uses 80 kb chunks, so recovery must handle a break crossing
an original chunk boundary. An additional 15m smoke test omitted the two
linker options and also produced 555 reads in one error-free block; its
results are in `default_linker.summary.json`.

## Measured results

These are assembly read-truth metrics over the selected windows, not NGC50 or
shared-VCF variant-switch metrics.

| window | configuration | evaluated reads | discordant | phase sets | read switch/flips |
|---|---|---:|---:|---:|---:|
| chr20:14.95–15.18 Mb | clean | 555 | 0 | 3 | 0 |
| chr20:14.95–15.18 Mb | automatic, margin 24 | 555 | 0 | **1** | 0 |
| chr20:14.95–15.18 Mb | automatic, margin 1 | 555 | 0 | **1** | 0 |
| chr20:47.60–47.82 Mb | clean | 510 | 0 | 3 | 0 |
| chr20:47.60–47.82 Mb | automatic, margin 24 | 510 | 0 | **1** | 0 |
| chr20:47.60–47.82 Mb | automatic, margin 1 | 510 | 0 | **1** | 0 |
| chr20:14.95–15.18 Mb, 80 kb chunks | clean | 534 | 0 | 3 | 0 |
| chr20:14.95–15.18 Mb, 80 kb chunks | automatic, margin 24 | 534 | 0 | **1** | 0 |
| chr20:14.95–15.18 Mb, 80 kb chunks | automatic, margin 1 | 534 | 0 | **1** | 0 |

Each case discovers two residual gaps after the clean pass; both close at tier
3. These runs improve connectivity without adding tagged reads. Margin 1
provides no extra benefit here. Unlike the earlier fixed-whitelist experiment,
automatic recovery examines all the local gap evidence rather than only
whitelist-triggered noisy intervals.

`invariants.json` records the original-block transformations. The verifier
asserts that every initially tagged read survives, each original block has
one consistent final PS/orientation, tier ordering is cumulative, no MSA het
indel appears before tier 3, escalation stops on success, and these regional
regressions retain zero read errors and end in one PS.

`LEFT_LINK=1` and `RIGHT_LINK=1` can coexist with `STATUS=partial`: the two
endpoints may link to different local proposal phase sets. Only `joined`
means the proposal establishes a connection between them. `NEW_SITES` counts
rows added within a tier, while `MSA_HET_SNPS/INDELS` count the accumulated
MSA heterozygous rows in the local proposal. Tier 1 has no MSA rows.

## Validation and remaining limits

Build and all five unit-test binaries pass. New tests cover cross-chunk gap
detection, partial extensions, complete-bridge selection, duplicate votes,
whole-block orientation preservation, unphased exact-candidate replacement,
untagged-read observations, cumulative SNP/indel admission, non-anchor masks,
requested-region boundaries, and gap-driven MSA coverage.

`make check` still fails before golden comparisons because its existing runner
uses obsolete `--phased-vcf-output` and positional input arguments. It was not
modified as part of recovery.

This implementation targets internal non-overlapping block gaps. It does not
invent a missing terminal anchor, count an isolated proposal as a bridge, or
attempt to merge overlapping blocks solely because their genomic spans
intersect. Failed evidence gates leave the connection open. Original clean
phasing errors are preserved rather than corrected by recovery. Whole-
chromosome accuracy, variant-switch metrics, and runtime scaling remain to be
validated before enabling the mode by default.
