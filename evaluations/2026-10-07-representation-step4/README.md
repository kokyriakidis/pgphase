# Representation repair: step 4, joint candidate loci

`GraphChunkBuildResult::joint_candidate_loci` now groups reference-equivalent
graph/BAM descriptions by complete normalized REF/ALT identity. Original rows
remain intact: topology IDs, physical edit coordinates, genotype/classification,
depths, observation channels, quality certificates and recovery source indices
keep their existing contracts. A canonical identity does not supply a new
genotype, quality score or phase gauge.

Both graph worker paths use `phase_joint_graph_candidates` for the initial clean
solve. Compatible clean heterozygous aliases temporarily share one representative
and one observation per read. Agreeing calls cast one vote; conflicting calls
abstain. The wrapper excludes duplicate category masks during k-means, shares
its result with compatible aliases, restores all original profiles/masks, and
rebuilds the interval tree. It preserves measured counts rather than summing
overlapping depths. Sparse alias-only coverage can temporarily reach the
representative without rewriting the original profile extent.

Projection requires matching admission/category contracts, no existing PS or
haplotype consensus, no MSA ALT list or homopolymer flag, and coverage of both
original physical edits/outside flanks for informative shifted observations.
One unsupported shifted observation keeps the entire cohort on the existing
solve. Different complete POS/REF/selected-ALT contexts also require an existing
independent graph observation at the representative for every informative alias
call. A primitive marker's REF cannot establish the full REF of a containing
snarl; exact complete contrasts may retain alias-only observations. In graph
modes where another ALT is the non-selected class, only
single-ALT sites qualify. Different DNA sequences and ALT/ALT contrasts remain
separate. Mixed graph/BAM classifications and recovered gauges are retained in
the locus view but are not coalesced into the clean solve. This stage does not
compact output rows or deduplicate late rescue votes; those consumers still need
complete allele contrasts and explicit source-confidence handling.

Graph read bounds are the first/last observed site positions, not the complete
GAF alignment span. The coverage test therefore uses a conservative extent;
it can reject otherwise valid projections and must not be widened by substituting
canonical coordinates. It does not certify a containing snarl's reference class.

The whole-chunk overlay rebuilds the locus table after recovery reindexing.
With `--phase-matrix-dump`, `.chunkN.joint-loci.tsv` records the canonical edit,
original member index, graph site ID, BAM provenance, raw edit, category and
phase/haplotype state. The initial solve reports duplicate/projected counts.

## Verification

The starting executable is the verified step-3 binary
`fe0ae0862c1f3596921218fee6b3f4b5f6cdcef9d6ee674df1cd687cc3d69b50`,
preserved with source snapshots under `test_data/tmp_representation_step4/`.
The final production executable is
`7bccceb86f200cb11ae7526d215d368c9ddc47255945747101fad008a024758f`.

The focused production-adapter fixtures run with
`./test_graph_bam_adapter --joint-loci` and average 2.5 ms including startup.
They compare two/four equivalent descriptions with one description, checking
read labels, vote count and margin, conflict abstention, sparse alias-only
coverage, all restored channels/qualities, original depths/provenance, consistent
alias orientation, stable source membership/certificates, different repeat
lengths, invalid reference, ALT/ALT, noisy exclusion and non-reference ALT classes.
They also verify that incomplete physical coverage prevents projection. Mutants
disabling projection, choosing a conflicting allele, losing original profiles,
ignoring ALT sequence, accepting non-reference ALT classes or bypassing physical
coverage fail the focused checks. A seventh mutant borrowing the parent
reference class also fails. Exact failure counts and timings are recorded in
the JSON results.
See `verify_fast.py`, `fast-checks.json` and `mutation-checks.json`.

The 4–5 Mb owning chunk recognizes **31 loci with 63 descriptions**, including
26 graph/BAM pairs with different classification contracts. The 5–6 Mb chunk
recognizes **seven loci with 14 descriptions**, including five mixed graph/BAM
classes. Their initial graph-only tables contain five and two duplicate loci,
respectively, with **zero** eligible clean projections. Consequently, both
owning matrices and emitted candidate/VCF rows, primary tags and parental statuses
are unchanged. The preserved mixed contracts are the next confidence/contrast
problem, not a reason to broaden clean admission. Parental results remain
3,862 correct/26 discordant/238 unphased in 4–5 Mb and
4,001/23/232 in 5–6 Mb. See `owner-checks.json` and `audit_owners.py`.

The 65–66 Mb owning chunk recognizes **11 loci with 22 descriptions** and safely
projects one insertion cohort. Its complete saved matrix, candidate/VCF rows,
primary tags and parental states remain identical: 2,962 correct, six discordant,
three unphased. The unchanged largest-block regression passes 742 assertions.

The accepted full chr20 replay preserves all 82,287 candidate rows, 64,483 VCF
records and 256,612 primary HP/PS assignments and parental states: 230,965 correct,
6,620 discordant and 19,027 unphased. N50 remains **955,496 bp**, largest block
**3,012,193 bp** across 250 variant blocks. See `full-checks.json`.

Build, unit tests, 47 phasing predicate cases (1,551 assertions), HiFi/ONT
TSV/VCF goldens and HiFi t1/t4 determinism pass. The identity helper still agrees
with `bcftools norm` on all 64,483 frozen chromosome records, including 504
realigned records. HiPhase measurement state validates the unchanged original
alignments, truth, panel, helper and competitor identities before reusing its
measurements. No expectation or read floor has been relaxed.

All **107/107 registered window/mechanism checks** pass in 714.1 seconds cold.
The cached standard `make window-tests` passes **18,138 assertions in five cases**,
including all four replay-cache checks.
All **141 native and full panel measurements** and HiPhase contract classifications
match step 3. Among 118 windows where HiPhase correctly phases at least 80% of
original primary truth-scorable reads, the full closure contract remains 59 passing
and 59 failing; total correctness and dominant connected-core comparisons are
unchanged. No new closure is claimed. See `window-checks.json`, `panel-audit.json`
and `panel-comparison.json`.

The provisional executable `a72c2fe268faa5417fd7c6795b474003de644e19748a5439af2f30bfec39df14`
passed the initial unit/owning checks. Review identified missing shifted-read
coverage handling; its in-progress broad runs were stopped and preserved under
`provisional-windows/` and `provisional-full/`. Their partial checks are not
accepted evaluation results. The final fixtures/coverage guard were added before
fresh final-binary replays.

The intermediate `52348f07d64b79615a10aaac781200a85dd0d56ba60a8f082e6219cdde03a6df`
passed 106/107 named checks. Its largest-block owner regressed from at least
2,962 correct/at most six discordant to 2,913 correct/55 discordant. The saved
65 Mb matrix explains the confidence defect: at the short deletion alias near
65,996,195, 24 source REF calls and one ALT call lack a callable containing-snarl
parent near 65,996,178. Normalized edit equivalence does not make their full
reference classes equivalent. The independent-parent guard restores the owning
check's 742 unchanged assertions; this rejected broad projection is not the
accepted result. Interrupted runners were also found still finishing children;
reports from shared work paths were discarded. Final acceptance uses the isolated
`accepted/` directories and the fixed production checksum above.

## Reproduction

Use the benchmark Python containing `pysam` for output audits:
`/home/kokyriakidis/micromamba/envs/bench-phasers/bin/python`.

```sh
make -j10
make unit-tests gap-dev-check
./test_graph_bam_adapter --joint-loci
python3 evaluations/2026-10-07-representation-step4/verify_fast.py
bash evaluations/2026-10-07-representation-step4/replay_matrix.sh \
  test_data/tmp_representation_step4/matrix4-before \
  test_data/tmp_representation_step4/pgphase.before
bash evaluations/2026-10-07-representation-step4/replay_matrix.sh \
  test_data/tmp_representation_step4/matrix4-final
bash evaluations/2026-10-07-representation-step4/replay_matrix.sh \
  test_data/tmp_representation_step4/matrix5-final ./pgphase CHM13#0#chr20:5000001-6000000
bash evaluations/2026-10-07-representation-step4/replay_matrix.sh \
  test_data/tmp_representation_step4/matrix65-before \
  test_data/tmp_representation_step4/pgphase.before CHM13#0#chr20:65000001-66000000
bash evaluations/2026-10-07-representation-step4/replay_matrix.sh \
  test_data/tmp_representation_step4/matrix65-final ./pgphase CHM13#0#chr20:65000001-66000000
make check
python3 evaluations/2026-10-06-representation-step1/check_windows.py \
  --binary ./pgphase --out test_data/tmp_representation_step4/accepted/windows
make window-tests
bash evaluations/2026-10-06-shared-insertion-source/replay.sh \
  test_data/tmp_representation_step4/accepted/full-current
```

After the corresponding complete replays, run `audit_owners.py`, `audit_full.py`,
`audit_panel.py` and `compare_panel.py`. Panel comparisons use original
overlapping primary alignments, parental truth and the fixed HiPhase result;
abstentions remain in the 80% denominator, and rescue phase sets are excluded
from the dominant connected-core count.
