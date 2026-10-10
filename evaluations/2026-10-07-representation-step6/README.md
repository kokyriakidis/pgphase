# Representation repair: step 6, original full-context sequence evidence

The joint evidence now has a source-independent physical sequence scorer.
It reconstructs both normalized selected alleles across the union of complete
parent REF spans and edits, with 16 outside reference bases per side. Every
other parent ALT and literal REF outside an ALT/ALT pair remain competing
hypotheses. A primitive REF call cannot establish full parent REF when another
branch changes elsewhere in that parent.

The original primary BAM records are already loaded by the whole-chunk solve;
the diagnostic extractor uses those records without another BAM query. Both
outer flanks must be aligned, internal I/D remain in the query slice, and
reference skips, unknown query bases, incomplete coverage and ambiguous primary
alignment identities abstain. SEQ already follows reference orientation for
reverse alignments. Original MAPQ, all base qualities (including unknown 255)
and physical query coordinates survive unchanged. Context coordinates are
one-based inclusive; query coordinates are zero-based half-open.

The score is global unit-cost edit distance, **not a calibrated posterior**.
A selected allele is uniquely nearest only if strictly closer than the other
selected allele and all competing parent sequences. A tie or closer unselected
allele leaves it unresolved. No threshold based on MSA provenance, source HP,
truth labels or source channel preference is introduced. Qualities are saved
for subsequent calibration; they do not become invented confidence here.
Contexts or sequences longer than 4,096 bases and unavailable reference flanks
are unsupported, bounding diagnostic work without new phasing admission rules.

This stage runs under `--phase-matrix-dump` for duplicate complete pairs.
`.chunkN.joint-contexts.tsv` saves reference and complete sequence hypotheses;
`.joint-parents.tsv` saves the original parent descriptions; `.joint-sequences.tsv`
saves one query slice, original quality bytes, source reduction and scores per
molecule/contrast. The existing original-channel matrix and molecule provenance
remain intact. With diagnostics disabled, the scorer does no work. Existing
phasing admission, graph/BAM profiles, parent guards, component joins, recovery
and output remain unchanged. Truth is evaluation-only.

## What the original reads show

| Owning replay | Duplicate pairs | Molecule records | Complete scored slices | Uncovered |
|---|---:|---:|---:|---:|
| Default 4–5 Mb | 31 | 1,749 | 1,714 | 35 |
| Default 5–6 Mb | 7 | 413 | 412 | 1 |
| Default 65–66 Mb | 11 | 391 | 386 | 5 |
| Existing whole-snarl 4–5 Mb mode | 23 | 1,370 | 1,337 | 33 |

Every context in these replays is supported. All source conflicts have complete
query coverage. The two default repeat-heavy owners contain 111 conflicts:
**55 uniquely favor a selected allele, 47 favor another parent allele, and nine
tie**. The split is 53/43/eight at 4–5 Mb and two/four/one at 5–6 Mb.
The 65–66 Mb owner has no source conflicts. Whole-snarl mode has 51:
32 favor a selected allele, 18 another parent allele, and one ties.

This explains why neither MSA verification nor choosing one channel resolves
the representation problem: the current two selected hypotheses can omit the
sequence that fits the molecule best. These are measured sequence rankings,
not assertions of genotype truth or new read accuracy. The next inference stage
must retain competing parent alleles and uncertainty instead of certifying a
primitive REF or promoting all nearest-selected calls blindly. It still needs
a validated quality/error model and a component-level diploid allele choice.

## Verification and fast replay

`make gap-dev-check` and `make unit-tests` include the new context tests.
There are **545 checks**, including 500 independent scalar dynamic-programming
score fixtures, complete REF/ALT and ALT/ALT contexts, source order reversal,
padded aliases, other parent alleles changing outside a primitive marker,
repeat-shifted I/D, reverse alignments, clips, incomplete coverage, skips,
unknown sequence and original/unknown qualities.

`audit_sequences.py` independently reconstructs every real context from original
parent descriptions and canonical keys against FASTA, validates each query
slice and quality byte against the original BAM via pysam, and reproduces all
**3,849 real scores** with scalar dynamic programming. It also checks exact
membership/provenance against the step-5 molecule table. This one-time independent
real-data audit takes 126.1 seconds; it is separate from the development loop.

Seven production mutations fail cleanly: omitting parent ALTs, omitting literal
REF competition, accepting selected ties, ignoring reference skips, erasing base
qualities, shifting query extraction and replacing global with local alignment.
They fail four, one, one, one, one, eight and six checks respectively.

Over 100 executions including startup, focused fixtures average **6.44 ms**.
Saved-state production scoring averages **26.78 ms**, **12.81 ms**, **48.27 ms**
and **16.92 ms** for the four owners above. Saved-state scoring accesses neither
BAM nor the pipeline. These measurements include concurrent broad verification
load. See `fast-checks.json`, `mutation-checks.json` and `sequence-checks.json`.

All four owning replays preserve original candidate/read-channel/quality matrices,
candidate rows, VCF calls, primary HP/PS and parental states. Default owner
correct/discordant/unphased counts remain 3,862/26/238, 4,001/23/232 and 2,962/6/3;
whole-snarl mode remains 3,866/27/235 under its own unchanged options.
`audit_contrasts.py` retains all step-5 original-index and channel checks, including
the three real whole-snarl ALT/ALT pairs and 99 collapsed other-ALT flags.

## Reproduction

Use the benchmark Python containing pysam for independent BAM/FASTA audits:
`/home/kokyriakidis/micromamba/envs/bench-phasers/bin/python`.
Run from the repository root:

```sh
make -j10
make unit-tests gap-dev-check
make check
bash evaluations/2026-10-07-representation-step5/replay_matrix.sh \
  test_data/tmp_representation_step6/accepted/matrix4-final
bash evaluations/2026-10-07-representation-step5/replay_matrix.sh \
  test_data/tmp_representation_step6/accepted/matrix5-final ./pgphase CHM13#0#chr20:5000001-6000000
bash evaluations/2026-10-07-representation-step5/replay_matrix.sh \
  test_data/tmp_representation_step6/accepted/matrix65-final ./pgphase CHM13#0#chr20:65000001-66000000
bash evaluations/2026-10-07-representation-step5/replay_whole_snarl.sh \
  test_data/tmp_representation_step6/accepted/whole-final
./test_allele_context --replay \
  test_data/tmp_representation_step6/accepted/matrix4-final/matrix.chunk0.joint-contexts.tsv \
  test_data/tmp_representation_step6/accepted/matrix4-final/matrix.chunk0.joint-sequences.tsv
python3 evaluations/2026-10-07-representation-step6/verify_fast.py
/home/kokyriakidis/micromamba/envs/bench-phasers/bin/python \
  evaluations/2026-10-07-representation-step6/audit_sequences.py
/home/kokyriakidis/micromamba/envs/bench-phasers/bin/python \
  evaluations/2026-10-07-representation-step6/audit_owners.py
python3 evaluations/2026-10-07-representation-step6/audit_contrasts.py
python3 evaluations/2026-10-06-representation-step1/check_windows.py \
  --binary ./pgphase --out test_data/tmp_representation_step6/accepted/windows
make window-tests
bash evaluations/2026-10-06-shared-insertion-source/replay.sh \
  test_data/tmp_representation_step6/accepted/full-current
```

After completed replays, run the local `audit_full.py`, `audit_panel.py` and
`compare_panel.py`. HiPhase measurements reuse the step-5 identity-validated
state only if input/index, panel, truth, competitor and helper identities still
match. `manifest.json` covers final production/extraction sources, original input
metadata, evidence-state file hashes and unchanged test floors.

## Full-chromosome verification

The final chr20 replay preserves all 82,287 candidate rows, 64,483 VCF records
and 256,612 primary HP/PS and parental states: 230,965 correct, 6,620 discordant
and 19,027 unphased. N50 remains **955,496 bp**, largest block **3,012,193 bp**,
across 250 variant blocks. See `full-checks.json`. No new gap closure or
biological accuracy gain is claimed for this sequence-evidence stage.

Build, units, all 47 phasing predicate cases (1,551 assertions), HiFi/ONT
TSV/VCF goldens and HiFi t1/t4 determinism pass without new compiler warnings.
Starting executable SHA256:
`cf849f5fa6fac264f2b5d7b01f6c768590da9cf8ac578e478264014ea53b1750`.
Final executable SHA256:
`e6a9aae89bdbe011be7d7759f67efcfe1f29b9c7c77d2c12298a867f37fdde68`.

All **107/107 registered window/mechanism checks** pass in 719.9 seconds cold.
All **141 native/full panel measurements and HiPhase classifications** are
identical to step 5. Among 118 windows where HiPhase correctly phases at least
80% of original primary truth-scorable reads, 59 still pass and 59 still fail
the full closure contract. Its original primary denominator includes abstentions;
total correct and dominant connected-core correct must meet HiPhase, excluding
rescue PS from the core. No expectation, closure certificate, panel or read floor
changed. See `window-checks.json`, `panel-audit.json`, `panel-comparison.json`
and the identity-validated `hiphase-verified.tsv(.json)`.

The cached standard `make window-tests` passes **18,138 assertions in five
cases**, plus all four replay-cache tests. The final build retains the evaluated
production checksum.
