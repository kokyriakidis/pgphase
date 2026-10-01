# A complete BAM source fills tied and absent graph SNP edges

At chr20:12,256,072 `A>G` to 12,269,535 `TA>T`, HiPhase joins the
flanks correctly but the graph recovery output split them. A primary
MAPQ-60 read calls the clean left SNP and the right clean SNP at
12,277,078 at Q40 on both bases. The imported right BAM source
(PS 12,361,654) has a complete path with no weak or quality cuts.
The graph path checker nevertheless rejected the whole right block:
it encounters balanced two-haplotype ties and later a graph SNP pair
at 12,411,075 to 12,458,901 with no callable GAF pair. A complete
BAM source can supply those absent edges, but the clean SNP pair still
has to establish the interblock orientation independently.

An unrestricted complete-source shortcut was unsafe. It joined three
other chr20 blocks at 7.9, 54.5 and 62.6 Mb across graph edges with
one-haplotype support or reversal. The trial lost 713 truth-correct
reads (229,131 to 228,418), and was discarded. The accepted rule lets
an independently complete right BAM source fill balanced ties with at
least two votes from each haplotype or zero-vote graph edges only after
a supported graph SNP edge with at least two votes from each haplotype.
One-haplotype and reversal edges remain vetoes. The normal physical
clean-SNP likelihood and left-block path checks still apply. No truth
label participates in the runtime decision.

The 12–13 Mb owning-chunk regression now has the left SNP, recovered
deletion, and right SNP in one phase set with the same ALT haplotype.
Its earlier clean SNP at 12,243,047 is in the same phase set with the
opposite ALT haplotype. The window spans, with 4,147/4,168 local
truth-scored reads correct (99.50%); 121/151 separated reads are
correct (80.13%). The exact VCF rows are checked by
`src/test_gap_windows.cpp`, and the panel expectation requires a span.

The full chr20 comparison uses the preceding accepted output in
`/tmp/pgphase-left-deletion-final-full/` as baseline and the new output
in `/tmp/pgphase-source-anchored-final-full/`:

| Metric | Baseline | Source-backed path |
|---|---:|---:|
| VCF variant keys | 62,361 | 62,361 |
| Truth-scored tagged reads | 236,866 | 236,866 |
| Truth-correct / discordant reads | 229,131 / 7,735 | 229,131 / 7,735 |
| Read phase sets | 697 | 695 |

No variant keys are added or removed. Exactly 760 VCF sample fields from
the right block (12,269,535–12,717,796) change PS 12,361,654 to
11,813,446, with their GT orientation flipped together. The clean left
SNP, recovered deletion, and right clean SNP now have the same ALT-haplotype
orientation. In the phased BAM, 2,004 source-block reads and 23
auxiliary read-only assignments move to the corresponding left PS with
HP flipped; no read gains or loses a tag. The two former read blocks have
1,906/1,914 and 1,997/2,004 truth-correct reads independently; the joined
block has
3,903/3,918 correct, exactly their sum. This directly excludes a
whole-block parental switch at the new join. The final full VCF and BAM
are byte-identical to the initial narrow-rule trial, including after
tightening the minimum graph anchor to two votes per haplotype.

Reproduce the focused regression with:

```bash
make -j4 test_gap_windows
./test_gap_windows 'chr20 gap windows' \
  -c 'graph / window 12256072-12269535' -r compact
```

The full comparison was run with the same chr20 input fixture:

```bash
./pgphase collect-graph-variation \
  --ref test_data/chm13v2.0.chr20.renamed.fa \
  --bam test_data/HG002_chr20_hifi_mapped_to_CHM13_chr20_annotated.bam \
  --sites test_data/chr20.sites.striped.vcf.gz \
  --gaf test_data/HG002.chr20.annotated.coord.gaf.gz \
  -r 'CHM13#0#chr20' -t 4 \
  -o /tmp/pgphase-source-anchored-final-full/candidates.tsv \
  --phased-vcf-out /tmp/pgphase-source-anchored-final-full/phased.vcf \
  --phased-bam-out /tmp/pgphase-source-anchored-final-full/phased.bam
```

Validation: `make -j4 pgphase test_gap_windows`, `make unit-tests`,
`make check`, and the full `make window-tests` panel passed. The panel
reported 3,557 assertions in 45 cases; one additional exact-allele
assertion was then added and passed in the focused 31-assertion replay.
`git diff --check` is clean.
