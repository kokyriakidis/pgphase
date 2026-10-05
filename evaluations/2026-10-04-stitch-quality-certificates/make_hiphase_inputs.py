from pathlib import Path
import pysam
for mb,beg,end in [(35,35445000,35560000),(57,57810000,57900000)]:
 root=Path(f'test_data/tmp_gap_fix50/hiphase/{mb}');root.mkdir(parents=True,exist_ok=True)
 with pysam.AlignmentFile('test_data/HG002_chr20_hifi_mapped_to_CHM13_chr20_annotated.bam') as src, pysam.AlignmentFile(str(root/'input.bam'),'wb',header=src.header) as dst:
  for r in src.fetch('CHM13#0#chr20',beg,end):dst.write(r)
 pysam.index(str(root/'input.bam'))
 with pysam.VariantFile(f'test_data/tmp_gap_fix50/baseline/{mb}/phased.vcf') as src,pysam.VariantFile(str(root/'input.vcf.gz'),'wz',header=src.header) as dst:
  for r in src:
   if not beg<=r.start<end:continue
   for s in r.samples.values():s.phased=False;s['PS']=None
   dst.write(r)
 pysam.tabix_index(str(root/'input.vcf.gz'),preset='vcf')
