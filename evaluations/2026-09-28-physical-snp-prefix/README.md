# Physical SNP seams and disconnected graph prefixes (2026-09-28)

The 3,529,324–3,542,977 chr20 gap remained split although 11 distinct
MAPQ-60 original BAM alignments call both clean SNPs at base quality >=30.
All 11 choose the same parental phase relation (seven G/T, four A/A pairs).
The targeted BAM subsolve demoted the right SNP as `NoisyCandHet`/`LOW_COV`,
so the transferred read profile lacked most of these physical calls. The
right SNP also came from a multiallelic snarl whose two biallelic SNP rows
shared a graph locus. Treating the rows as two consecutive independent graph
edges falsely rejected the right block's SNP path.

The final within-chunk stitch reads the physical boundary bases from the
original indexed BAM, normalizes padded graph SNP alleles, and chooses the
phase relation by base/MAPQ-weighted likelihood. For multiple paired reads,
both boundary haplotypes must occur. The wrong-parity posterior must be
<=0.001. Original graph SNPs are checked as a path, with rows from the same
multiallelic snarl counted as one locus. The right block at 3.529 Mb has a
zero-vote graph-SNP path cut at 3,573,979–3,597,791, but two phased deletion
candidates lie inside that cut. The earlier stitch had already joined the
3,573,979 SNP to the 3,596,663 deletion using a single validated BAM
molecule; the subsequent deletion links are also supported. An initial
prefix-only trial closed 3.529 Mb but broke that established 3.573 Mb join.
The retained rule refuses a partial join if the cut contains a candidate of
the right block, any primary BAM read spans it, or any profile calls right-
block candidates on both sides. Thus the 3.529 Mb owning chunk stays split,
while its narrower panel replay spans because the later cut falls outside
that replay. Physical graph alleles are normalized before testing the cut.

The retained implementation joins the 8,638,940–8,662,670 and
54,547,514–54,569,072 owning-chunk gaps. The former merged block has
2,145/2,145 truth-correct reads on the left, 71/71 near the gap, and 297/297
on the right in the full chromosome replay. The 54–55 Mb owning chunk
has nine discordant reads before and after the change. The unsafe
62,623,253–62,645,168 whole-block control remains split; its right graph
path is not continuous.

The same SNP stitch also joins 5,345,085–5,350,509 in the 5–6 Mb owning
chunk. The upstream 5,309,406–5,345,085 bridge remains joined; the distant
5,393,615 graph block beyond the known weak source cut stays independent.
The newly joined local edge has at least 99% parental read concordance in its
owning-chunk regression, with no block switch. The full chromosome's read
truth counts above include this join.

Full chr20 with the same annotated BAM, striped graph catalog, GAF, reference,
and parental truth map:

| State | Read phase sets | Scored reads | Correct | Discordant | Accuracy | Tracked gaps closed |
|---|---:|---:|---:|---:|---:|---:|
| Prior certified seam | 727 | 236,849 | 229,092 | 7,757 | 96.7249% | 7/23 |
| Guarded SNP/prefix stitch | 722 | 236,847 | 229,090 | 7,757 | 96.7249% | 9/23 |

The two fewer truth-scored reads are both truth-correct; no new discordant
reads were introduced. VCF variant keys remain 62,150 in the prior run;
the guarded run has four additional recovered heterozygous keys in the 54 Mb
region after the earlier join opens a second-pass recovery seam.

Results: `/tmp/pgphase-singleton-full/` and
`/tmp/pgphase-guarded-final-full/`. The target manifest is
`evaluations/2026-09-27-remaining-hiphase-correct-gaps/targets.tsv`.
The owning-chunk assertions preserve the 3.573 Mb join, reject the unsafe
3.529 Mb split, and require the 5.31, 8.638, and 54.547 Mb joins. The 3.529 Mb panel expectation
records only its narrower replay. These tests are in `src/test_gap_windows.cpp`
and `src/test_gap_windows_expect.tsv`.

Reproduce from repository root:

```bash
mkdir -p /tmp/pgphase-guarded-final-full
./pgphase collect-graph-variation \
  --ref test_data/chm13v2.0.chr20.renamed.fa \
  --bam test_data/HG002_chr20_hifi_mapped_to_CHM13_chr20_annotated.bam \
  --sites test_data/chr20.sites.striped.vcf.gz \
  --gaf test_data/HG002.chr20.annotated.coord.gaf.gz \
  -r 'CHM13#0#chr20:1-66210255' -t 8 \
  -o /tmp/pgphase-guarded-final-full/candidates.tsv \
  --phased-vcf-out /tmp/pgphase-guarded-final-full/phased.vcf \
  --phased-bam-out /tmp/pgphase-guarded-final-full/phased.bam
make unit-tests
make check
make window-tests
```

The final full window panel passes with 1,723 assertions in 34 test cases;
`make unit-tests` and `make check` pass as well. The local truth scoring
helper is `/tmp/score_pgphase_weakcut.py`.
