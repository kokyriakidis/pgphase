# Why the largest block stops 12,277 bp before HiPhase

The current pgphase endpoint is 66,194,203; HiPhase ends at 66,206,480.
The production binary is unchanged, SHA256
`477aaecb8adadf8f6f20bc380676de81f15d51f8eacbdb787c2cf5a8bf18140e`.
The earlier explanation based only on GQ2 and MAPQ3–15 was incomplete.
This is a candidate/recovery limitation in a repeat, with a separate read
phasing problem. It is not a missing catalog site or a hard region cutoff.

## Exact filtering mechanism

The graph catalog contains the exact terminal T>A SNP, including both allele
walks. At default graph MAPQ5, the production walk matcher finds four reference
walks and zero alternate walks. The site is dropped as `ref_only`, before
phasing. Lowering MAPQ to 1 finds eight reference walks and still no alternate.
The trace explicitly includes forward and reverse traversals: this is not an
orientation-blind substring check. Three other original reads do not traverse
this graph site. See `catalog-terminal.vcf.txt`, `graph-observations.json`,
`walk-matcher-trace.txt`, and `filtered-sites.json`.

At that same physical position, the original BAM has eight callable Q40 bases:
six A and two T; three alignments have no base there. This reproduces HiPhase's
DeepVariant AD=2,6, DP8, GT1|0, GQ2. All six BAM A observations traverse
the graph's reference walk. Thus BAM placement and GAF walk placement disagree
in this repetitive sequence; the catalog already provides the alternative.
This audit establishes the disagreement, not which repeat-copy placement is
biologically correct. All eleven truth-scorable original reads are paternal,
so these observed T/A classes do not establish parental segregation.

The first graph SNP after pgphase's endpoint, 66,194,253, likewise has five
reference observations and zero alternate at MAPQ5. At MAPQ1 it has 28 reference
and 20 alternate observations. Low mapping qualities therefore really do
explain loss of usable graph candidates in the earlier tail. Merely reducing
depth leaves the reference-only terminal SNP excluded.

The BAM route has two additional limitations:

- The terminal noisy interval is 66,206,004–66,210,255. BAM classification masks
  candidates in noisy intervals for MSA recall (`collect_var.cpp`,
  `classify_cand_vars_pgphase`). In the MAPQ1 owning BAM replay, recall logs
  `MsaFire Skipped ... 11 reads (3 full) ps=-1`. There is no qualifying phase set
  with both full-cover haplotypes, and three full-cover reads are below the
  minimum depth of five (`align.cpp`, `collect_noisy_reg_aln_strs`). Reducing
  minimum depth to three changes the branch to `NoHap`, but does not emit the
  exact terminal T>A call. Shorter-window/partial-coverage recall is a possible
  future investigation; this experiment does not prove that it will recover
  a valid diploid terminal marker.
- `collect_phase_set_seams` only nominates intervals between differing phased
  anchors. A terminal tail has no downstream anchor. Consequently targeted
  recovery does not nominate it. The existing whole-chunk BAM fallback runs,
  but only attaches observations to existing graph candidates and rescues
  read assignments; it does not import every private terminal BAM variant.

## Same original reads, including abstentions

The additional interval is 66,194,204–66,206,480. It overlaps 54 original primary
alignments, all truth-scorable and paternal. Every HiPhase alignment was checked
against the corresponding original start/end, CIGAR and sequence. Pgphase's
synthetic output BAM is used for tags only. Each pgphase phase set is oriented
on its whole output cohort, and HiPhase uses the independently measured whole
PS63182011 orientation from the preceding largest-block audit. Rescue PS values
at or above 1e9 are excluded from the connected correct core.

| Output | Correct tail reads / 54 | Correct endpoint reads / 11 | Correct tail core |
|---|---:|---:|---:|
| Current frozen chromosome / default owner | 25 (46.3%) | 5 (45.5%) | 5 |
| HiPhase | 3 (5.6%) | 1 (9.1%) | 3 |
| Owning graph replay, MAPQ1 | 32 (59.3%) | 9 (81.8%) | 32 |
| Owning graph replay, MAPQ0/depth1/alt1 | 31 (57.4%) | 8 (72.7%) | 31 |

The default output has 29 discordant tail assignments and no abstentions;
HiPhase has six discordant and 45 unphased. HiPhase's whole phase set has an
internal parental switch, so its most favorable local flip would yield six
correct tail reads, or four correct endpoint reads, with the same abstentions.
Even that favorable flip is below 80%. A raw terminal VCF phase label therefore
does not demonstrate that HiPhase closes this tail correctly.

MAPQ1 extends the owning block to 66,207,824, beyond HiPhase's endpoint, without
calling its exact T>A marker. The permissive replay reaches 66,207,887 in the
original owner PS and creates additional short terminal blocks. These are real
span extensions, but neither reaches 80% correctness across the added interval.
Both lose 74 formerly correct assignments from the default owning replay;
these losses are in BAM rescue phase sets, with no old connected-core correct
read lost. Total correct owner assignments fall from 584 to 566/565, despite
56/55 newly correct assignments. Checking only the final coordinate would
miss both the poor intervening accuracy and these losses.

The next useful fix would address terminal recovery and reconcile repeat
placements with established flanking haplotypes. Blanket lowering of MAPQ or
MSA depth is not a validated fix. This investigation changes evidence and
explanation only; it does not change pipeline behavior or claim a new closure.

## Reproduction and checks

The three graph owning replays take about 2–4 seconds each on this dataset;
no full-chromosome rerun is needed for this investigation. Run from the repo:

```bash
bash evaluations/2026-10-05-terminal-endpoint/replay.sh
/home/kokyriakidis/micromamba/envs/bench-phasers/bin/python \
  evaluations/2026-10-05-terminal-endpoint/audit_terminal.py
/home/kokyriakidis/micromamba/envs/bench-phasers/bin/python \
  evaluations/2026-10-05-terminal-endpoint/audit_replays.py
```

The replay writes ignored state in `test_data/tmp_gap_fix70/`. The read audit
uses frozen production output `tmp_gap_fix69/frozen_final/0`, original input
alignments, parental truth, and the existing HiPhase comparator. The replay
audit uses the preceding `audit_reads.py` helper for consistent scoring.
All counterfactual-count assertions, exact HiPhase alignment comparisons, and
GAF/production matcher count checks pass. `bash -n` and Python compilation pass.
No production source or executable changed; existing native gap development
checks pass.
