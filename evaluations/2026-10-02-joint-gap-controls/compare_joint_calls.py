from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
import pysam,subprocess,json,time
root=Path('test_data/tmp_gap_next28').resolve()
source=Path('test_data/tmp_gap_next27/accepted/phased.vcf')
bam_path='test_data/HG002_chr20_hifi_mapped_to_CHM13_chr20_annotated.bam'
reference='test_data/chm13v2.0.chr20.renamed.fa'
binary='/tmp/pgphase-gap-next12/HiPhase/target/release/hiphase'
contig='CHM13#0#chr20'
windows=[(7264321,7280346),(24121713,24131707)]
def run(pair):
 lo,hi=pair;d=root/f'{lo}-{hi}';d.mkdir(exist_ok=True)
 with pysam.AlignmentFile(bam_path) as src, pysam.AlignmentFile(str(d/'input.bam'),'wb',template=src) as out:
  for read in src.fetch(contig,lo-100001,hi+100000):out.write(read)
 pysam.index(str(d/'input.bam'))
 with pysam.VariantFile(str(source)) as src,pysam.VariantFile(str(d/'input.vcf.gz'),'wz',header=src.header) as out:
  for record in src:
   if lo-100000<=record.pos<=hi+100000:out.write(record)
 pysam.tabix_index(str(d/'input.vcf.gz'),preset='vcf',force=True)
 for mode in ['default','local']:
  out=d/mode;out.mkdir(exist_ok=True);t=time.monotonic()
  cmd=[binary,'--bam',str(d/'input.bam'),'--vcf',str(d/'input.vcf.gz'),'--reference',reference,'--ignore-read-groups','--min-mapq','5','--threads','1','--output-vcf',str(out/'phased.vcf.gz'),'--output-bam',str(out/'phased.bam')]
  if mode=='local':cmd+=['--disable-global-realignment']
  with (out/'run.log').open('w') as log:r=subprocess.run(cmd,stdout=log,stderr=log)
  print(lo,hi,mode,r.returncode,round(time.monotonic()-t,3),flush=True)
  if r.returncode:continue
  blocks={}
  with pysam.VariantFile(str(out/'phased.vcf.gz')) as v:
   for x in v:
    s=list(x.samples.values())[0];gt=s.get('GT');ps=s.get('PS')
    if s.phased and ps and gt and len(set(gt))>1:blocks.setdefault(ps,[]).append(x.pos)
  spans=[ps for ps,p in blocks.items() if min(p)<=lo and max(p)>=hi]
  print('SPANS',mode,lo,spans,flush=True)
with ThreadPoolExecutor(max_workers=2) as pool:list(pool.map(run,windows))
