# Test independent BAM phase-block transfer in an open graph gap

## Question and target

Can a BAM subregion be phased independently, its phase blocks injected intact,
and consecutive blocks stitched using all their site/read evidence? Test the
noncentromeric chr20:61,738,239–61,747,506 nomination. The detector's surrounding
graph seam is 61,732,321–61,757,551. Use the same BAM and reference for every
standalone arm, and existing saved same-BAM HiPhase for evaluation only.

The model is appropriate: source blocks have independent HP gauges, so retain
candidate alleles, per-read observations and source membership; infer relative
orientation before merging labels. A numeric BAM PS is not a certificate that
every internal allele edge is correct. Production has no truth or competitor
input. Experiments do not change source allele representations or realign reads.

## Native BAM subregions

All local scores use the same 96 physical gap-overlapping read names and a
majority parental orientation per output PS. These are local conditional scores,
not whole-chromosome accuracy or runtime comparisons.

| Solve | Phased / correct / discordant local reads | Local PS |
|---|---:|---:|
| Exact nominated gap, MAPQ30 | 69 / 69 / 0 | 1 |
| Gap plus 100 bp, MAPQ5 | 69 / 69 / 0 | 1 |
| Gap plus 1 kb, MAPQ5 | 69 / 69 / 0 | 1 |
| Gap plus 10 kb, MAPQ5 | 63 / 61 / 2 | 2 |
| Gap plus 50 kb, MAPQ30 | 53 / 41 / 12 | 1 |
| Gap plus 50 kb, MAPQ5 | 53 / 41 / 12 | 1 |
| Actual graph seam plus 1 kb, MAPQ5 | 55 / 41 / 14 | 1 |
| Current graph+recovery | 83 / 79 / 4 | 4 |
| Saved HiPhase, native DV | 92 / 91 / 1 | 1 |

The exact/100 bp/1 kb BAM arms emit only one phased heterozygote,
T→G at 61,738,239. They assign useful local reads but provide no diploid chain
to the right boundary. Adding context exposes the right deletion but also the
source's conflicting orientation. A shorter region can avoid that conflict
without solving the boundary connection. This is why local purity alone is an
insufficient criterion for selecting an injected source.

Each native BAM arm completes in approximately 0.3–0.6 s on this fixture.
The original command vectors are saved in `test_data/tmp_gap_next38/*/args.json`.
`local61-bam-comparison.json`, `local61-small-padding.json` and
`local61-seam1000.json` retain the local measurements.

## Trace existing transfer

The accepted targeted BAM source has 14 phased sites in PS 61,690,751 spanning
61,690,751–61,747,507 internally. All seven private phased rows in the graph
seam reach the graph chunk in one new PS, 61,690,752, with unchanged source GT.
The two shared clean SNPs retain graph representation and graph ownership.
The source provenance retains their BAM phase relation for stitching.

| Private source row (internal coordinate) | Source GT / PS | Imported GT / PS | After stitch GT / PS |
|---|---|---|---|
| SNP 61,738,239 | 1\|0 / 61,690,751 | 1\|0 / 61,690,752 | 0\|1 / 61,487,312 |
| Deletion 61,747,507 | 0\|1 / 61,690,751 | 0\|1 / 61,690,752 | 0\|1 / 61,755,063 |

The split happens after successful injection. Both original callable matrix
pairs between these rows are ALT/ALT. Both reads belong to the source PS;
under the source's opposite GT orientations, both pairs conflict with that PS.
There is no source-consistent callable pair at this edge. Earlier independent
physical clean-SNP/deletion checks also detect a source-gauge contradiction.
Blindly preserving the nominal PS would preserve this internal reversal.

`source-transfer.json` records every phased source row and its exact raw-key
match before and after stitching. Shared graph-walk rows need sequence metadata
and intentionally do not have an exact raw-key match. `audit_transfer.py`
recomputes the two boundary calls and their source-gauge consistency.

## Trials

### Trust the nominal source block

An evaluation-only control clears source weak/quality cuts before transfer,
then runs the normal stitcher. This preserves the nominated gap's private
rows in one BAM block, but the next graph block remains separate. The SNP and
deletion retain opposite GTs in that common PS, contradicting their supported
ALT/ALT relation. Local read counts improve to 83 / 80 / 3 and the owning chunk
improves one correct read, illustrating that aggregate read purity can miss
an incorrect allele connection. Reject this control; the temporary source change
is absent from production. `trust-block-parity.json` and `trust-block-local.json`
record its output. The permanent owning test is extended to reject the wrong
SNP/deletion relationship whenever those rows acquire a common PS; a future
correctly oriented join remains allowed.

### Retain every phased source flank in the live table

Remove both the focused-retry and weak-path limitations on transferring private
flank rows and adopting shared unphased repeats. Four initial owning-chunk tests
retain all read truth states, while adding 92/28/16/81 VCF keys at 61/57/48/56 Mb.
The selected gap stays open.

Full chromosome exposes the failure: 3,622 previously correct reads become
incorrect, four protected joins reopen, and 13 old phased SNPs disappear.
Conditional read accuracy falls 96.995350% to 95.509795%; correct reads fall
230,072 to 226,618 and discordant reads rise 7,127 to 10,654. The table gains
4,131 keys and loses 14, but these gains cannot justify the wrong connections.
Candidate insertion changes later live block ownership and stitching, even
though the original observations and genotype calls are unchanged. Reject.
`rejected-all-context-parity.json` and `rejected-all-context-transfer.json`
record this whole-chromosome counterexample. No existing gate is weakened.

### Retain complete internally supported private source blocks

Keep existing shared graph-site adoption rules and weak-cut checks. Broaden
only private flank transfer from selected BAM blocks with at least two phased
sites and no weak internal path cut. This retains all usable phased private
members of those blocks; it does not grant new orientation or override an
established graph genotype. Initial owning 61/23/53 Mb checks preserve every
read truth state and all old VCF keys, adding 72/7/0 keys respectively.
An experimental 61 Mb completeness guard confirms two complementary private
flank insertion rows were omitted by the starting transfer and survive this
trial. That extra expectation is not retained because the trial fails the full
chromosome gates. The permanent test retains the additional SNP/deletion
orientation guard, which remains relevant to future transfer designs.

Full chromosome retains aggregate conditional read accuracy (96.997065%), but
loses 31 scored and 26 correct assignments. Thirty-five previously correct
reads become incorrect and 40 become unphased. The table adds 2,652 keys and
loses 11, ten old phased SNPs disappear, and the protected 54.547514–54.569072 Mb
gap reopens. N50 falls 774,189 to 739,888 bp. Reject this version as well;
`cut-free-parity.json` and `cut-free-transfer.json` record the failure.

## Retained state and design implication

Restore the accepted production source and executable exactly. No new gap is
closed or runtime/accuracy improvement claimed. The executable SHA256 remains
`244d3e8f39b91c5ac17d9ad3c5c05a28ad42223a9cb12db681c77e0c571ef3a8`.
Keep the stronger permanent 61 Mb orientation regression. It allows a future
correct join but rejects the negative control's opposite SNP/deletion GTs
inside a single phase set. No existing span, coverage or accuracy gate changes.

The implementation needs to distinguish full source evidence from live graph
ownership: retain immutable complete BAM block/site/read records for pairwise
orientation, while graph-owned sites retain their established membership.
Independent source labels need source identity as well as PS. Use all usable
allele observations from the two complete block records when testing a join;
count each independent spanning molecule once per block pair and keep clean-SNP
and conflict checks. Apply an accepted whole-block flip and label union only
after those checks pass. A conflicting or internally unsupported source edge
cannot be repaired by trusting its numeric PS. This is the proposed design,
not a claim that these rejected live-table trials implement it safely.

Full unrestricted and cut-free trials complete in 254.98 and 243.35 seconds at
eight threads. They are concurrent validation timings, not controlled runtime
comparisons. The unrestricted native-panel replay is aborted after the full
chromosome gate fails; no complete passing panel result is claimed for a rejected
candidate. The restored binary can reuse its previously verified 105 native
panel outputs after executable/CLI/input/output hash checks; the newly expanded
owning test executes natively. Give each test runner a fresh
`PGPHASE_TEST_WORKDIR`: the runner's hard-linked replay cache can contaminate
other runs when they share output directories. A discarded shared-directory
validation is not counted as a product regression or a passing gate. The final
isolated suite passes **7,121 assertions / 70 cases**. All standalone units and
**1,031 predicate assertions / 42 cases** pass. The negative control fails exactly
the new SNP/deletion common-PS genotype relation; no read, coverage or span gate
is changed to accept it. The restored executable is byte-identical to the
accepted one, so the earlier verified chromosome accuracy and NG50 remain the
accepted result. No new chromosome run of the restored binary is needed.

Reproduce the source/transfer audit from the repository root:

```bash
python3 evaluations/2026-10-03-bam-block-transfer/audit_transfer.py \
  --source test_data/tmp_gap_next37/gauge-preserved/61/matrix.recovery.chunk0.window2.chunk-1.recovery-source.tsv \
  --input test_data/tmp_gap_next37/gauge-preserved/61/matrix.chunk0.recovery-input.tsv \
  --final test_data/tmp_gap_next37/gauge-preserved/61/matrix.chunk0.recovery-final.tsv \
  --phase-set 61690751 --left 61738239 --right 61747507 \
  --output /tmp/bam-block-source-transfer.json
PGPHASE_TEST_WORKDIR=/tmp/pgphase-independent-block-check make window-tests
```
