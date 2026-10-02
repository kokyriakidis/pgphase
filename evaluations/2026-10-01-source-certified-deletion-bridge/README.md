# Original-source certification closes the 11.573 Mb gap

## Retained changes

The physical SNP-to-deletion helper looked up a BAM source-path certificate
using the deletion's **current** PS. A row transferred to another block could
therefore borrow that block's certificate even though its original source had
an unsupported cut or its allele orientation did not match the destination.

`bam_source_site_path_supported` now requires unique adoptable provenance, a
complete original BAM source without weak or quality cuts, and another oriented
site at a distinct coordinate from that source in the current block. Every
such current source row must have the same source-to-current allele orientation.
Homozygotes do not orient the block. The existing graph path checks remain.

The clean SNP pair remains first. When it has zero callable pairs, its nearer
right MSA deletion can now use the existing physical bridge. MAPQ30, Q10 SNP,
Q30 deletion, independent source/graph paths, and the existing quality-weighted
wrong-parity bound of 0.05 are unchanged. Candidate rows and allele meanings
remain separate. Production uses neither truth nor competitor data and has no
fixture coordinate exception or new realignment.

Matrix dumps append the independent graph and BAM alleles after the existing
working-allele column. This lets a local audit distinguish a missing working
call from conflicting channels without rerunning the chromosome.

## New connection and its evidence

**chr20:11,573,074–11,586,531, 13,457 bp.** The downstream independent BAM
source contains complementary insertions at 11,586,531, the verified `CT>C`
deletion at 11,591,585, and clean SNPs beginning at 11,599,138. Its complete
source path and consistent current allele gauge certify the nearer deletion.

No high-quality BAM molecule calls both preselected clean SNPs at 11,573,074
and 11,599,138. One primary MAPQ60 molecule instead calls the left SNP REF at
Q40 and the deletion REF across Q30 verified flanks; it ends before the right
SNP. Both REF calls imply the same haplotype in the retained source gauges.
`bridge-read-calls.tsv` records the physical calls; `audit_bridge.cpp` uses the
production callers against the same reference and BAM. It is a fixture audit,
not production phasing code. Its arguments are BAM, reference, contig, left SNP
position, right SNP position, and internal one-base deletion position; the
fixture SNP alleles are `G>C` and `G>A`.

| Window read measurement | Accepted baseline | Retained |
|---|---:|---:|
| Truth-scorable input reads | 134 | 134 |
| Phased / correct / discordant | 111 / 108 / 3 | 111 / 108 / 3 |
| Local read concordance | 97.2973% | 97.2973% |
| Correct fraction placed in one block | 47.7612% | 64.9254% |

The archived HiPhase run phases 125 local reads, 118 correct and seven
incorrect (94.4% local concordance; 88.0597% correctly separated). This is an
existing competitor output, not a new run on the current pgphase call set.
The new case is below the strict 98% HiPhase target selection threshold and
does not reduce the fourteen previously tracked unresolved targets.

The 100 kb padded replay stays split: it lacks the earlier BAM source and full
right-block context. The panel therefore replays the owning 11–12 Mb chunk.
That chunk retains 4,092 scored reads, 4,002 correct and 90 discordant. Its
absolute error ceiling covers the whole chunk; the local three-error tally
above is measured separately in `local-truth.json`.

## Counterexample and regressions

- Add **20 adapter checks** for original versus destination certificates,
  missing/duplicate provenance, weak and quality cuts, consistent/reversed
  gauges, homozygotes, co-located rows, and independent source anchors.
- The actual old certificate lookup, compiled in a scratch binary with the
  zero-pair fallback, joins **58,366,458–58,385,702 incorrectly**. Its owning
  correct/discordant counts change 3,492/57 to 3,428/122. The existing owning
  regression fails the forbidden join, parental switch, and concordance checks.
  The retained original-source check keeps it split.
- Add the newly closed gap to the committed panel with **spans=1** and measured
  0.97 concordance / 0.64 separation floors. Update the existing 11.599 Mb
  owning regression for this reviewed physical join: check six exact allele
  rows, their original genotypes, shared PS, no parental reversal on both
  intervals, and the unchanged 90-error owning ceiling. The original baseline
  fails five connection/separation assertions. This does not admit or replace
  the previously rejected incomplete focused MSA solve.
- Every old panel span and accuracy floor remains unchanged. The panel grows
  **92 to 93** windows; its total connection floor grows **74 to 75**. The old
  required-site entries remain unchanged; the new gap has no strictly interior
  heterozygote, and its external bridge/anchors are checked by exact VCF key.
- All **105 distinct native window commands** run on the final binary with four
  workers in **231.96 s**. The new owning panel row reuses an existing identical
  owning-chunk command. Exact command and binary SHA256 verified scoring passes
  **4,082 assertions / 53 cases**. Unit, predicate, port-parity, and original
  upstream-parity gates pass: 491, 27, and 168,696 framework assertions,
  respectively, plus the standalone units. No new compiler warning.

An audit of working-channel graph-indel paths across all fourteen remaining
cases, and separate graph-channel calls in seven cases in the 21, 35, 36 and
58 Mb chunks, found no path passing the inspected two-haplotype statistical
edge criteria. No relaxed conflict or MAPQ threshold was retained.

## Default chromosome comparison

Accepted baseline: `/tmp/pgphase-gap-next9/final-full`.
Retained native run: `/tmp/pgphase-gap-next10/source-full`, eight threads,
178.09 s, identical annotated BAM, catalog, GAF, reference and defaults.

| Metric | Baseline | Retained |
|---|---:|---:|
| Primary read names | 256,586 | 256,586 |
| Phased/truth-scored reads | 237,101 | 237,101 |
| Correct / discordant | 229,921 / 7,180 | 229,921 / 7,180 |
| Read concordance | 96.971755% | 96.971755% |
| Read PS | 686 | 684 |
| VCF records / blocks | 62,798 / 340 | 62,798 / 339 |
| VCF span N50 | 739,888 bp | 739,888 bp |
| Spanned expanded panel windows | 74/93 | 75/93 |

Every primary read name, HP, VCF key and genotype is unchanged. Exactly 609
read PS labels and 60 VCF PS labels change at this connection. No variant key
is added or removed. Fourteen tracked competitor targets remain open.

Reproduce the focused checks with:

```bash
make unit-tests
./test_gap_windows 'certified deletion bridge preserves the 11.599 Mb source alleles'
./test_gap_windows '[moved-deletion]'
```

Run all integration gates with `make window-tests`. Validation, exact output
comparisons, physical-call evidence, both failing counterexamples, and the
native command inventory accompany this report.
