# Reuse a certified BAM graph-block join (2026-09-28)

The 3,529,324–3,542,977 chr20 SNP gap has 11 MAPQ-60 original BAM reads
calling both boundary bases at quality >=30. All 11 support one phase
relation. pgphase nevertheless left it split in the 3–4 Mb owning chunk:
a subsequent BAM recovery stitch had already joined the right block across
a 3,573,979–3,597,791 graph-only path gap, and the final physical SNP
check demanded GAF calls from both haplotypes across that gap.

The recovery matrix shows two original graph blocks on either side of the
3.573–3.597 Mb gap. The main BAM stitch joins them through its existing
validated SNP-to-deletion path. A MAPQ-60 molecule reaches the first
deletion; its reference base there has quality 10, so the final SNP stitch
must not reinterpret that single base as a fresh high-quality bridge. It
can instead reuse the join the main stitch already accepted. A later
3,786,668–3,806,891 graph edge has only one haplotype's GAF call; it too
connects two original graph blocks joined by the main BAM stitch.

The final SNP path check now retains original graph phase-set identities by
stable site ID across the two recovery passes. It also records each graph
site's phase-set identity immediately after the *first main BAM stitch*,
before source attachments. An edge with insufficient two-haplotype GAF
support can pass only when its two sites came from different original graph
blocks, the main BAM stitch had already put those blocks in one phase set,
and no GAF read votes for the reversing relation. Edges inside an original
block still require their usual two-haplotype read path. The unsafe 62.623–
62.645 Mb whole-block control remains split.

The 3–4 Mb owning chunk now puts 3,529,324 and 3,542,977 in the same phase
set with opposite ALT haplotypes. The established 3,573,979-to-3,596,663
deletion bridge stays in that same set. Local truth scoring remains
3,756/3,832 correct (76 discordant) before and after. Its regression tests
both joins and the 3.529 Mb gap's read concordance. Full chr20 with the
same annotated BAM, striped graph catalog, GAF, reference, and parental
truth map:

| State | Read phase sets | Scored reads | Correct | Discordant | Accuracy | Tracked gaps closed |
|---|---:|---:|---:|---:|---:|---:|
| Guarded SNP/prefix stitch | 722 | 236,847 | 229,090 | 7,757 | 96.7249% | 9/23 |
| Certified graph path | 717 | 236,847 | 229,090 | 7,757 | 96.7249% | 10/23 |

Both runs emit the same 62,154 VCF variant keys. The same 236,847 truth-
scored read names are present in both BAM outputs; no read changes from
truth-correct to discordant or vice versa. The only newly closed tracked
target is 3,529,324–3,542,977. Thirteen tracked targets remain
open. Outputs: `/tmp/pgphase-guarded-final-full/` and
`/tmp/pgphase-certified-join-full/`. The target manifest is
`evaluations/2026-09-27-remaining-hiphase-correct-gaps/targets.tsv`.

Reproduce from repository root:

```bash
mkdir -p /tmp/pgphase-certified-join-full
./pgphase collect-graph-variation \
  --ref test_data/chm13v2.0.chr20.renamed.fa \
  --bam test_data/HG002_chr20_hifi_mapped_to_CHM13_chr20_annotated.bam \
  --sites test_data/chr20.sites.striped.vcf.gz \
  --gaf test_data/HG002.chr20.annotated.coord.gaf.gz \
  -r 'CHM13#0#chr20:1-66210255' -t 8 \
  -o /tmp/pgphase-certified-join-full/candidates.tsv \
  --phased-vcf-out /tmp/pgphase-certified-join-full/phased.vcf \
  --phased-bam-out /tmp/pgphase-certified-join-full/phased.bam
make unit-tests
make check
make window-tests
```

The full window panel passes 1,728 assertions in 34 test cases; unit tests
and BAM/ONT validation gates pass. Local truth scoring helper:
`/tmp/score_pgphase_weakcut.py`.
