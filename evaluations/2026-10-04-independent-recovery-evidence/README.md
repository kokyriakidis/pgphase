# Preserve independent verified recovery evidence

## Defects

The starting binary preserved selected targeted calls against the late whole-BAM
overlay, but still rejected fixed-consensus calls when graph ownership had no
coherent clean-SNP gauge or a prior same-block HP disagreed. Overlapping selected
sources also used the first available working allele. Their immutable matrices
borrowed physical qualities from that merged table, so an allele from one source
could acquire another source's certificate or lose its own flank certificate.

Starting binary SHA256:
`b4614c555ca06d411e68fc5f6f4e79bb9df8507c755d2172ff13557b6665fb47`.
Its saved full chromosome output is
`test_data/tmp_transfer_audit/full-private-admission/`.

## Implementation

1. After selecting a targeted source, admit every unambiguous fixed-consensus
   recall at an existing verified phased heterozygote. Graph ownership and a
   read's prior HP no longer veto that allele. Rows, genotypes and source HP/PS
   are not relabeled by admission. Accepted focused source solves now consume
   their pending calls too.
2. When graph ownership has an unresolved gauge, cache the source's original
   internal-path certificate before adding recalled calls. Those calls cannot
   remove an original weak cut or certify a whole-block connection. Normal
   stitch evidence still establishes relative orientation.
3. A BAM-only working row takes its observations from the selected source
   owning its genotype. Every other source interpretation survives independently.
   Shared/unowned disagreements produce an explicit working BAM conflict
   (`-2`). The independent graph allele remains available. Conflict markers
   survive reordering, recovery retry and whole-chunk fallback; ordinary missing
   calls (`-1`) may still be filled.
4. Every immutable source matrix computes qualities from its own alignment and
   allele, including context-only sites. A working conflict has no physical
   certificate. Verified MSA calls take precedence over supplementary physical
   projections. Contradictions within that verified tier abstain, and distinct
   proposals survive as provenance rather than extra molecule votes.

Supplementary physical projections still require their original coherent graph
gauge and fixed-consensus context. These are not independently verified MSA
calls. No production branch reads truth, competitors, or fixture coordinates.
There is no new alignment algorithm or relaxed stitch threshold.

## Experiments

The first trial admitted verified calls but made every overlapping disagreement
unknown, even on a BAM row whose selected genotype already identified its
source. It retained the known wrong 19.374 Mb split but lost two correct read
assignments in owning61. This trial was rejected. The `trial-*-parity.json`
files retain its results.

The source-owner refinement preserved all old correct assignments in seven
noncentromeric owning chunks (3, 11, 19, 21, 34, 41 and 61 Mb). Owning19 gained
two correct assignments without increasing its 41 errors; the other six had
unchanged correct/error totals. All keys survived. See `owner-*-parity.json`.
These probes precede the final additional conflict persistence guard across
recovery retries; the final native regression audit covers the completed code.

Restored complementary MSA rows at 11,862,429 rise DP62→68 and at 34,835,166
rise DP51→56, retaining complementary GT and the independent source PS. The
new native admission regression fails on the starting binary in all four depth
assertions (two rows at each locus). It passes after the fix. Truth gates retain
owning11's 4,004 correct/≤90 errors and owning34's 3,539 correct/≤28 errors.

At 3,597,791, read `m84031_231217_034919_s2/16389926/ccs` has REF/Q40 in one
source and ALT/Q0 in another. Both original calls remain independently
available. The working BAM slot is `-2` with no certificate; the independent
graph REF survives. A permanent native regression checks these source calls
and the conflict before/after whole-chunk fallback. The transfer conservation
regression now checks every phased source call in the immutable snapshot,
including omitted context sites, as well as known and blocked working calls
across late fallback. It still requires fallback to fill missing calls.

`audit_evidence.py` checks canonical source preservation, selected working
observations, admission and overlay persistence without truth. The final
targeted audit retains all 66,029 checked independent source calls; 64 overlapping
disagreements consist of 63 selected-owner working calls and one shared
abstention. All 71,612 existing overlay calls and that explicit conflict survive.
The earlier seven-chunk audit likewise reports no loss, overwrite or rejected
verified recall. Counts include solve/stage records, not unique chromosome
alleles. Reports are `regression-evidence-audit.json` and
`owner-evidence-audit.json`.

## Validation

Targeted completed-code regressions: 90 assertions in four cases pass,
including the known wrong 19.374 Mb whole-block join guard. The native
admission test uses the existing panel boundaries for its owning replays.
Unit tests and 1,532 predicate assertions in 47 cases pass. HiFi and ONT
candidate/phased-VCF goldens match, and HiFi one/four-thread output is
deterministic. No regression floor, ceiling, required site or span expectation
was relaxed.

## Full chr20 result

The final eight-thread native graph run uses the same BAM, catalog, GAF,
reference and options as the starting full output. It produces:

| Metric | Before | After |
| --- | ---: | ---: |
| Scored phased reads | 237,308 | 237,328 |
| Truth-correct reads | 230,212 | 230,232 |
| Discordant reads | 7,096 | 7,096 |
| Truth concordance | 97.009793% | 97.010045% |
| Scored read phase sets | 653 | 654 |
| Variant keys | 64,188 | 64,188 |
| VCF phase blocks | 331 | 331 |
| Span N50 | 806,449 bp | 806,449 bp |
| Connected tracked coordinates | 91/104 | 91/104 |

All 230,212 previously correct read assignments remain correct, and all
7,096 previously incorrect ones have unchanged truth status. Exactly 20
previously unphased reads become correct.
The additional reads occupy two independent read blocks, with 10 paternal/8
maternal and 1 maternal/1 paternal respectively; both source haplotypes are
represented in each PS. Exact changed tags are in `changed-read-tags.json`.
All VCF GT/PS values are unchanged;
ten rows at five insertion loci change only depth/allele fractions. The
restored sites are 10.325039, 11.862429, 19.373922, 34.835166 and 60.090437 Mb.
See `full-parity.json`, `changed-sites.tsv` and `block-audit.json`. No old SNP
block has a mixed new gauge, no block join is introduced, and no key or phased
SNP is lost. There is no new gap closure or fresh competitor comparison.
Thirteen tracked coordinates stay open: ten competitor nominations and three
controls. Evidence preservation is fixed without forcing those connections.

The final full run takes 376.23 s during concurrent native regression work;
this is not a controlled runtime benchmark. Exact input/binary and output
hashes are in `input-manifest.json` and `output-manifest.json`.

Final production binary SHA256:
`54e67e0911beb6402e2558b44ac5dbf31d04da8da8f26b93a66e9acd8a91348f`.

The complete native window suite passes all 8,060 assertions in 77 cases,
covering 104 panel coordinates through 110 freshly generated pipeline requests.
The rebuilt targeted suite additionally passes 90 assertions in four cases
after aligning the admission replay's 34 Mb boundary with the existing panel;
that correction leaves its native owning-chunk command and truth floors
unchanged. The final diagnostic panel audit also verifies independent source
preservation and late conflict persistence. See `window-tests.txt`,
`regression-after.txt` and `panel-evidence-audit.json`.
