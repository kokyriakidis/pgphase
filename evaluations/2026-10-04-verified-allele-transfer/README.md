# Verified BAM allele transfer audit

## Finding

Two additional evidence-loss defects were confirmed and fixed:

1. **Late whole-chunk BAM fallback overwrote known targeted BAM alleles.**
   At VCF 61,802,959, six reads lost their paired calls at the separate
   19/20-base insertion rows (12 changed observations). One complex insertion
   observation at VCF 61,806,160 also changed. The primary working alleles were
   unchanged, hiding the error from output-only comparisons. The BAM channel
   later used for read rescue no longer carried the selected source calls.
   Fallback now fills only missing observations and preserves an existing
   allele's query position and quality certificate.
2. **Graph-anchor admission discarded fixed-consensus calls from a private
   BAM block.** At VCF 21,594,343, a BAM-private source block has no shared graph
   SNP. The shared-SNP requirement nevertheless rejected 102 verified calls
   (51 complementary pairs). Both separate insertion rows had depth 11 despite
   62 callable reads. Fixed-consensus MSA observations now survive in a block
   whose heterozygous anchors are entirely BAM-private. Its independent gauge,
   source genotype and normal stitch requirements are preserved. Here ownership
   means an exact source heterozygote matching the selected graph allele
   representation. This exception
   does not admit supplementary physical corrections or join blocks by PS alone.

The preceding cross-block fix is retained: a deferred call in another BAM phase
set cannot be compared with the read's numerical HP in its current block.
The within-block contradiction veto remains. See
[the preceding evaluation](../2026-10-04-cross-block-msa-observations/README.md).

No new gap closure is claimed. Restoring evidence and accepting a haplotype
connection are separate operations.

## What was checked

The audit follows these boundaries:

| Stage | Check and result |
|---|---|
| Fixed-consensus recall → source admission | Trace each pending observation and its explicit admission reason. Distinguish MSA recall from supplementary physical correction. |
| Selected source → candidate merge | Resolve exact sequence keys through the same canonical aliases as production. Raw graph-walk keys are not sequence alleles. |
| Candidate/read reindexing | Compare source calls with the destination BAM channel after candidate sorting and read qname sorting. |
| Targeted matrix → whole-BAM fallback | Compare every already known clean/MSA BAM allele before and after attachment. |
| Completed chromosome | Compare every VCF row and read tag, score parental truth per read phase set, and audit old SNP block gauges and all panel connections. |

Eighteen owning 1 Mb chunks were replayed, beginning at 3, 7, 10, 11, 13,
15, 17, 20, 21, 24, 34, 35, 37, 41, 50, 57, 58 and 61 Mb. They exclude
the 25–30 Mb centromeric interval. All comparisons use the existing reference,
graph catalog, GAF and surjected BAM fixture. Production uses no truth labels,
competitor calls or genomic-coordinate special cases.

The fill-only fix's audit finds:

- **506,149 unambiguous mapped source observations retained**, with zero missing
  or changed alleles after transfer/reindexing.
- **481,431 existing clean/MSA BAM observations retained** across late fallback;
  the old code overwrote 13 of them.
- **76 conflicting read/site records** from overlapping source solves. A working
  matrix has one allele slot; the existing first-available-call policy remains.
  Each accepted solve also retains its own independent source matrix for stitch
  evaluation. These conflicts are listed separately, not counted as loss of an
  unambiguous observation.

These totals count records across replay stages, not unique chromosome alleles.
The audit covers selected mapped calls, not unconditional injection of every
candidate from every padded solve. Its raw census also includes 114,277 call
records whose candidates are absent from the working table and 1,629 owned by
another seam. Existing admission restricts new rows to the seam or eligible
complete source flanks and applies the source phasing category mask. The
immutable stitch matrix separately retains eligible phased source heterozygotes,
including flank-only blocks. A working-table omission is therefore distinct
from disappearance of the full source block's stitch evidence. This audit does
not prove that every such exclusion is optimal.

An additional seven owning replays with private-block admission retain all
157,089 unambiguous mapped calls and all 149,316 preexisting overlay calls.
Both counts overlap the first panel and must not be added to it.
The final binary's additional 32 Mb owning replay retains all 55,703
unambiguous mapped calls and all 55,071 known overlay calls, with no source
conflict or overwrite. It covers the other two insertion pairs changed in the
full chromosome output.

### Admission still has deliberate limits

Detailed admission traces in chunks 11, 15, 21 and 34 exposed 124 withheld
fixed-consensus observations: the 102 private-block calls above, 12 at internal
11,862,430 and 10 at internal 34,835,167. The latter two have inconsistent graph
gauges and remain rejected. A broader seven-chunk sample also exposes such
conflicts in other solves. Counts refer to their respective replay scopes.

Removing graph-gauge checks broadly was previously unsafe: it joined the
19.374 Mb counterexample and changed 368 previously correct reads into errors.
See [that rejected experiment](../2026-10-03-tandem-insertion-recall/README.md).
The new exception applies only when no source heterozygous anchor is graph-owned.

Verified site discovery does not certify every read's original CIGAR allele.
MSA calls and CIGAR backfill can disagree, particularly at shifted repeat
insertions. The fix preserves selected verified calls instead of replacing them
with a different solve's interpretation or adding another CIGAR prerequisite.
Physical quality and stitching checks retain their distinct roles.

## Full chr20 result

Baseline is the preceding cross-block fix, not the earlier committed output.
The overwrite fix alone leaves all output read tags and VCF rows identical.
Adding private-block admission produces:

| Metric | Before | After |
|---|---:|---:|
| Truth-scored phased reads | 237,308 | 237,308 |
| Truth-correct reads | 230,211 | 230,212 |
| Discordant reads | 7,097 | 7,096 |
| Read concordance | 97.009372% | 97.009793% |
| Unique primary output reads | 256,609 | 256,610 |
| Scored read phase sets | 653 | 653 |
| VCF variant keys | 64,188 | 64,188 |
| VCF blocks | 331 | 331 |
| VCF span N50 | 806,449 bp | 806,449 bp |
| Connected panel coordinates | 91 | 91 |

All 230,211 old correct reads remain correct. One old erroneous read becomes
unphased and one formerly unphased read becomes correct. Two read tags and six
VCF rows change; no keys disappear, no old SNP block develops a mixed gauge,
and no panel connection changes. Thirteen panel coordinates remain open,
including three controls. No fresh competitor comparison was run here.
The other four changed rows are two complementary insertion pairs: at
32,364,883 their depth increases 7→8, and at 32,725,929 it increases 51→60.
All six rows preserve GT/PS.
The chromosome runs shared resources with tests and are not a controlled
runtime benchmark.

## Reproduction and evidence

Request `--phase-matrix-dump PREFIX` on an owning graph replay. It writes:

- `PENDING`/`REJECTED` provenance in matrix dumps;
- source `msa-transfer-pending`, `msa-transfer-selected` and `msa-admission.tsv`;
- canonical `chunkN.transfer.tsv` and optional `transfer-retry.tsv`;
- `chunkN.bam-overlay-input.tsv` and `bam-overlay-output.tsv`.

Admission status `eligible` means the graph-level filter passed. The normal
source commit still rejects contradictory duplicate recalls and same-block
allele contradictions. The selected-source matrix and canonical transfer trace
establish which calls actually survived those checks.

Run the independent, truth-free conservation check:

```sh
python3 evaluations/2026-10-04-verified-allele-transfer/audit_transfer.py \
  --root test_data/tmp_gap_fix48/transfer-fixed \
  --output /tmp/transfer-audit.json --strict
```

`before-overlay-audit.json` and `after-audit.json` retain the 18-chunk comparison;
`admission-audit.json` records the four detailed admission replays;
`overlap-conflicts.json` records competing source calls;
`private-block-transfer-audit.json` checks the seven private-admission replays.
`private-block-32-transfer-audit.json` checks the additional final-binary replay.
`full-parity.json`/`block-audit.json` check the overwrite fix;
`private-block-full-parity.json`/`private-block-block-audit.json` check the final
behavior. `private-block-owning-parity.json` isolates the changed 21 Mb chunk.

The two permanent native regression cases cover both transfer boundaries:
the old code changes 13 preexisting observations and fails the conservation
case; it also reports depth 11 instead of 62 and fails four private-block
assertions. The fixed cases pass 10 and 23 assertions respectively. Their logs
are included with this report. Existing span expectations are unchanged.

Final validation passes: build without new warnings, all unit tests,
1,517 predicate assertions across 47 cases, HiFi/ONT TSV and VCF goldens,
HiFi one/four-thread determinism, and the full window suite's **8,006 assertions
across 76 cases** using 110 fresh native requests. Build, unit, golden and
window logs are included. `input-manifest.json` records fixture hashes and
the final binary SHA256
`b4614c555ca06d411e68fc5f6f4e79bb9df8507c755d2172ff13557b6665fb47`.
