# Certified single-molecule graph seam (2026-09-28)

The 54.894–54.912 Mb gap stayed split in its 54–55 Mb owning chunk although
a short replay closed it. One annotated BAM molecule spans both boundary SNPs
(MAPQ 60; base quality 40 at each) and calls the same phase relation as
HiPhase. Both chunk-local graph blocks retain consistent two-haplotype clean-SNP
paths after recovery transfer. The right graph block of the unsafe 62.623–
62.645 Mb one-read seam fails that whole-block path test.

A global `--min-block-link-reads 1` trial joined two 54 Mb gaps but changed
chr20 truth scoring to 227,821 correct / 9,026 discordant of 236,847 scored
reads and reopened two previously solved tracked gaps. Applying the one-read
setting only to targeted BAM sub-solves still produced 228,288 correct /
8,781 discordant of 237,069 scored reads. Both trials were discarded.

The retained post-transfer stitch considers only the closest phased clean
SNP pair in a recovery seam, exactly one distinct BAM read with callable
alleles on both sides, and a conservative sum of allele and mapping error
probabilities <=0.001. It joins only when both chunk-local original graph blocks
have consistent two-haplotype clean-SNP read paths. The 54–55 Mb owner now
puts `54894127 G>C` and `54912022 G>A` in the same phase set. Its truth score
is unchanged: 4,305 correct / 9 discordant of 4,314 scored reads. The 62 Mb
control remains split.

Full chr20 replay (`CHM13#0#chr20:1-66210255`) against the annotated BAM,
striped graph catalog, GAF, and parental truth map:

| State | Read phase sets | Scored reads | Correct | Discordant | Accuracy | Tracked gaps closed |
|---|---:|---:|---:|---:|---:|---:|
| Prior | 728 | 236,849 | 229,092 | 7,757 | 96.7249% | 6/23 |
| Certified seam | 727 | 236,849 | 229,092 | 7,757 | 96.7249% | 7/23 |

Results: `/tmp/pgphase-snp-insertion-final-full/` (prior),
`/tmp/pgphase-singleton-full/` (certified seam). Both outputs contain the same 62,150 VCF rows. Only 458 right-block
VCF rows change their phase-set label; none changes genotype. The new owning-chunk
regression is in `src/test_gap_windows.cpp`; its unsafe 62 Mb control is
retained.

Reproduction (from repository root):

```bash
mkdir -p /tmp/pgphase-singleton-full
./pgphase collect-graph-variation \
  --ref test_data/chm13v2.0.chr20.renamed.fa \
  --bam test_data/HG002_chr20_hifi_mapped_to_CHM13_chr20_annotated.bam \
  --sites test_data/chr20.sites.striped.vcf.gz \
  --gaf test_data/HG002.chr20.annotated.coord.gaf.gz \
  -r 'CHM13#0#chr20:1-66210255' -t 8 \
  -o /tmp/pgphase-singleton-full/candidates.tsv \
  --phased-vcf-out /tmp/pgphase-singleton-full/phased.vcf \
  --phased-bam-out /tmp/pgphase-singleton-full/phased.bam
make unit-tests
make check
make window-tests
```

The scoring helper used for this fixture is
`/tmp/score_pgphase_weakcut.py`; the committed gap target manifest is
`evaluations/2026-09-27-remaining-hiphase-correct-gaps/targets.tsv`.
