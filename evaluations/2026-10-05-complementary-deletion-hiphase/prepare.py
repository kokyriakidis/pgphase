from pathlib import Path
import pysam
root=Path('test_data/tmp_deletion_hiphase_audit')
for label,beg,end in [('short',50498244,50612066),('owner',50000000,51000000)]:
 p=root/label;p.mkdir(exist_ok=True)
 with pysam.AlignmentFile('test_data/HG002_chr20_hifi_mapped_to_CHM13_chr20_annotated.bam') as src,pysam.AlignmentFile(str(p/'input.bam'),'wb',header=src.header) as dst:
  for r in src.fetch('CHM13#0#chr20',beg,end):dst.write(r)
 pysam.index(str(p/'input.bam'))
 with pysam.VariantFile('test_data/tmp_gap_fix50/full-final/phased.vcf') as src,pysam.VariantFile(str(p/'input.vcf.gz'),'wz',header=src.header) as dst:
  for r in src:
   if not beg<=r.start<end:continue
   for sample in r.samples.values():sample.phased=False;sample['PS']=None
   dst.write(r)
 pysam.tabix_index(str(p/'input.vcf.gz'),preset='vcf')
