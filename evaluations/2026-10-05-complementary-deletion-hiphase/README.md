# Recheck the complementary-deletion gap against HiPhase

Target: chr20:50,548,245–50,562,066, outside the centromere.

## Current status: connected, including both deletion rows

This is the interval in `2026-10-03-complementary-deletion-backfill/README.md`.
It was closed by that fix, not one of the ten remaining competitor nominations.
A new owning 50–51 Mb graph/BAM run and the latest full chr20 output both have:

| VCF row | GT | PS |
| --- | --- | --- |
| 50,548,245 CA→C | 1\|0 | 50,127,297 |
| 50,548,245 CAA→C | 0\|1 | 50,127,297 |
| 50,562,066 A→G | 0\|1 | 50,127,297 |

The two deletion ALTs are opposite. The two-base deletion ALT and SNP ALT
are on the same haplotype. The existing owning regression passes all 27
assertions, including these exact representations, connection, parental
orientation, 3,887 scored / 3,876 correct reads, and the 11-error ceiling.
No pipeline code, threshold, test expectation or genotype changes this turn.

## Why the earlier source split, and what HiPhase uses

Run HiPhase 1.6.0-ac3f399 on the current pgphase sites, with genotypes unphased,
and complete original annotated BAM alignments overlapping either the 114 kb
short interval or the owning Mb. The short default run has global alignment
enabled; repeat it with `--disable-global-realignment` for the local-only
comparison. Decode explicit named read segments and their site ranges from
serial `-vvv` output. Do not match anonymous vectors by BAM iteration order.
The VCF heterozygote count and vector bounds are checked.

There are 12 distinct MAPQ60 primary molecules physically spanning both
boundaries. Their right SNP calls are present from the start (six REF, six ALT).
Both left deletion calls are missing in the initial BAM solve. The signal loss
therefore precedes graph transfer. It is not missing alignment coverage.

| Observation stage | One-base deletion/SNP pairs | Two-base deletion/SNP pairs |
| --- | ---: | ---: |
| pgphase initial BAM MSA matrix | 0 | 0 |
| pgphase accepted full-block MSA retry | 10 | 8 |
| pgphase transferred matrix | 10 | 8 |
| pgphase final matrix | 10 | 8 |
| HiPhase default short solve | 5 | 10 |
| HiPhase local-only short solve | 11 | 10 |

These are callable row-specific binary observations, not an accuracy measure.
BAM complementary MSA rows encode ALT absence on the other row; this must not
be equated with literal reference sequence in a separate alignment caller.
`boundary-evidence.json` preserves every molecule and allele at every stage.
The accepted retry's calls reach the working matrix unchanged.

The original ordinary backfill excluded these homopolymer deletion rows:
literal CIGAR REF cannot identify their complementary ALT/other contrast.
Missing calls also hid the source conflict needed to select the focused retry.
The retained fix uses only exact or equivalent complementary ALT observations
as temporary retry-admission evidence, restores them before transfer, and runs
the existing full-block BAM/MSA retry. Validated MSA calls supply the actual
imported block. Full adjacent blocks provide shared SNP gauge context; a short
replay lacks it. The new owning solve has 455 and 459 agreeing shared candidate
orientations on its two graph flanks. A direct vote on just the 12 boundary
molecules is not the full source evidence or the whole-block certificate.

Two molecules still have a missing two-base deletion call in pgphase's retry,
although both HiPhase modes call that row ALT:

- `m84031_231217_062403_s3/125176998/ccs`: a three-base CIGAR deletion,
  shifted relative to the candidate, with Q40 surviving flanks.
- `m84031_231217_034919_s2/98238767/ccs`: a three-base CIGAR deletion,
  with surviving flank qualities 40 and 17.

Their right SNP ALT calls are present in all stages. HiPhase can assign these
noisy sequences to a nearby two-base allele; pgphase's fixed MSA projection
abstains. This is an allele-assignment difference, not lost verified transfer.
The same difference occurs with global alignment disabled, so it cannot be
attributed solely to HiPhase's global alignment. The CIGAR audit is diagnostic,
not replacement phasing evidence. No third length is silently converted to a
candidate ALT. `missing-calls-physical.json` records the actual events.

## Read truth and what “HiPhase closes it” means

Use original primary-read physical spans to select gap-overlapping molecules.
Choose each output PS's parental orientation from all of its truth-scored reads,
then apply that orientation to the gap subset. Truth is evaluation-only.

| Solve | Scored gap reads | Correct | Errors | Accuracy |
| --- | ---: | ---: | ---: | ---: |
| pgphase owning Mb | 112 | 105 | 7 | 93.75% |
| pgphase latest full chr20 | 112 | 105 | 7 | 93.75% |
| HiPhase owning Mb | 109 | 103 | 6 | 94.50% |
| HiPhase default short | 109 | 103 | 6 | 94.50% |
| HiPhase local-only short | 107 | 100 | 7 | 93.46% |

Pgphase currently tags three more gap reads, yielding two more correct and one
more discordant assignment. Its core boundary block is joined; a separate
output-only rescue group is not another VCF gap.

In both default HiPhase runs the one-base deletion and right SNP share a PS,
but the two-base deletion row has **no PS**. A connected read block is not proof
that both complementary VCF rows were phased. Pgphase phases all three rows.

The wider Mb comparison also cautions against treating connectivity as accuracy:
pgphase has 3,887 scored / 3,876 correct / 11 errors (99.7170%) in three scored
read groups; HiPhase has 3,898 / 3,438 / 460 (88.1991%) in one group. Of those
HiPhase errors, 442 occur on molecules that pgphase assigns to the separate
left block PS 50,002,195; their parental orientation differs from the larger
joined HiPhase block. This is outside the requested 13.8 kb gap and is not a
full-chromosome competitor benchmark. It explains why admitting a broader
whole-region join is not automatically an improvement.

## Reproduction and provenance

Work directory: `test_data/tmp_deletion_hiphase_audit`.
Fresh graph source dumps: `test_data/tmp_gap_fix50/deletion-hiphase/50`.
The graph run uses the committed reference, BAM, catalog and GAF, default
settings, `-t 4 -r CHM13#0#chr20:50000001-51000000`, and
`--phase-matrix-dump PREFIX`. HiPhase uses the same reference and original
BAM alignments, current pgphase VCF sites with PS cleared and GT unphased,
`--ignore-read-groups`, and one thread for named trace collection or four
threads for the owning solve. `prepare.py`, `audit.py` and `score.py` reproduce
input extraction and scoring after those solves.

The native existing regression is:

```bash
PGPHASE_TEST_WORKDIR="$PWD/test_data/tmp_deletion_hiphase_audit/window" \
LD_LIBRARY_PATH=/usr/lib/x86_64-linux-gnu \
./test_gap_windows 'complementary deletion recovery closes the owning 50 Mb gap'
```

Pgphase SHA256:
`9bda79a94b69241f8b9083ee902dceea6d0c4560f583ff09cfce2543afcb45df`.
HiPhase SHA256:
`4a949e4b3de9bea21d24c656f8b7ab1a8b4cec87e87be707099eec56eac63996`.
`manifest.json` hashes the source matrices, input BAM/VCFs and traces.
