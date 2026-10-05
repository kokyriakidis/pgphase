from pathlib import Path
from collections import Counter,defaultdict
import json,re
import pysam
root=Path('test_data/tmp_deletion_hiphase_audit'); pg=Path('test_data/tmp_gap_fix50/deletion-hiphase/50')
left,right=50548245,50562066
with pysam.AlignmentFile(str(root/'owner/input.bam')) as bam:
 spanning={r.query_name for r in bam.fetch('CHM13#0#chr20',left-1,right) if not(r.flag&(4|256|512|1024|2048)) and r.mapping_quality>=30 and r.reference_start<=left-1 and r.reference_end>=right}
print('physical spanning',len(spanning))
def load_pg(path):
 variants={};calls=defaultdict(dict)
 for line in path.read_text().splitlines():
  f=line.split('\t')
  if f[0]=='VAR':variants[int(f[1])]=f
  elif f[0]=='OBS':calls[f[1]][int(f[2])]=int(f[3])
 indices={str(left)+(':DEL1' if int(f[6])==1 else ':DEL2'):i for i,f in variants.items() if f[2]=='50548246' and f[3]=='D' and f[6] in ['1','2']}
 indices.update({str(right)+':SNP':i for i,f in variants.items() if f[2]==str(right) and f[3]=='X'})
 assert len(indices)==3,(path,indices)
 return {q:{key:calls[q].get(i,-1) for key,i in indices.items()} for q in spanning}
pg_sources={
 'pg_initial':load_pg(pg/'matrix.recovery.chunk0.window2.initial.chunk-1.flags5004.tsv'),
 'pg_retry':load_pg(pg/'matrix.recovery.chunk0.window2.focused.seam2.chunk-1.recovery-source.tsv'),
 'pg_transfer':load_pg(pg/'matrix.chunk0.recovery-input.tsv'),
 'pg_final':load_pg(pg/'matrix.chunk0.recovery-final.tsv')}
convert={'Reference':0,'Alternate':1,'Ambiguous':-1,'NoOverlap':-1}
def parse_hi(path):
 blocks=[]
 for line in path.read_text().splitlines():
  if 'Solving problem: PhaseBlock' in line:
   m=re.search(r'coordinates: "(.+):(\d+)-(\d+)".*num_variants: (\d+)',line);assert m
   blocks.append({'beg':int(m[2]),'end':int(m[3]),'nvars':int(m[4]),'rows':[]})
  elif 'read segment #' in line and 'read_name:' in line:
   m=re.search(r'read_name: "([^"]+)", alleles: \[(.*?)\], quals: \[(.*?)\], region: (\d+)\.\.(\d+)',line);assert m,line
   aa=m[2].split(', '); beg,end=int(m[4]),int(m[5]);assert len(aa)==end-beg
   blocks[-1]['rows'].append((m[1],beg,end,aa))
 with pysam.VariantFile(str(root/'short/input.vcf.gz')) as vcf:
  variants=[r for r in vcf if len(set(next(iter(r.samples.values()))['GT']))==2]
 out={q:{str(left)+':DEL1':-1,str(left)+':DEL2':-1,str(right)+':SNP':-1} for q in spanning}
 for b in blocks:
  vs=[r for r in variants if b['beg']<=r.start<=b['end']];assert len(vs)==b['nvars'],(len(vs),b['nvars'])
  idx={i:str(left)+(':DEL1' if len(r.ref)==2 else ':DEL2') for i,r in enumerate(vs) if r.pos==left and len(r.ref)>len(r.alts[0])}
  idx.update({i:str(right)+':SNP' for i,r in enumerate(vs) if r.pos==right})
  for name,beg,end,aa in b['rows']:
   assert end<=len(vs)
   if name not in out:continue
   for i,key in idx.items():
    if beg<=i<end:out[name][key]=convert[aa[i-beg]]
 return out
hi=parse_hi(root/'short/trace.log'); hlocal=parse_hi(root/'short/local-only.trace.log')
all_stages={**pg_sources,'hi_default':hi,'hi_local':hlocal};keys=[str(left)+':DEL1',str(left)+':DEL2',str(right)+':SNP'];summaries={}
for stage,reads in all_stages.items():
 summary={}
 for key in keys[:2]:
  table=Counter((row[key],row[keys[2]]) for row in reads.values());summary[key]={f'{a},{b}':n for (a,b),n in sorted(table.items())}; summary[key+'_paired']=sum(n for (a,b),n in table.items() if a in [0,1] and b in [0,1])
 summaries[stage]=summary;print(stage,summary)
result={'coordinates':[left,right],'physical_mapq30_spanners':len(spanning),'summary':summaries,'reads':[{ 'name':q,**{stage:reads[q] for stage,reads in all_stages.items()}} for q in sorted(spanning)]}
(root/'evidence.json').write_text(json.dumps(result,indent=2)+'\n')
