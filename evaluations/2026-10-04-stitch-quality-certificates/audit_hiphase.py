import json,re
from pathlib import Path
from collections import Counter,defaultdict
import pysam
root=Path('test_data/tmp_gap_fix50')
result=[]
for mb,lo,hi in [(35,35498368,35516845),(57,57854341,57866713)]:
 p=root/'hiphase'/str(mb);rows=[];v=[];blocks=[]
 for l in (p/'trace.log').read_text().splitlines():
  if 'Solving problem: PhaseBlock' in l:
   m=re.search(r'coordinates: "(.+):(\d+)-(\d+)".*num_variants: (\d+)',l);assert m
   blocks.append({'chrom':m[1],'beg':int(m[2]),'end':int(m[3]),'nvars':int(m[4]),'rows':[]})
  elif 'read segment #' in l and 'read_name:' in l:
   m=re.search(r'read_name: "([^"]+)", alleles: \[(.*?)\], quals: \[(.*?)\], region: (\d+)\.\.(\d+)',l);assert m,l
   aa=m[2].split(', ');beg,end=int(m[4]),int(m[5]);assert len(aa)==end-beg
   blocks[-1]['rows'].append((m[1],beg,end,aa))
 with pysam.VariantFile(str(p/'input.vcf.gz')) as f:
  candidates=[r for r in f if len(set(next(iter(r.samples.values()))['GT']))==2]
 with pysam.AlignmentFile(str(p/'input.bam')) as f:
  for b in blocks:
   vs=[r for r in candidates if b['beg']<=r.start<=b['end']]
   rs={r.query_name:r for r in f.fetch(b['chrom'],b['beg'],b['end']+1) if not r.flag & (4|256|512|1024) and r.mapping_quality>=5}
   assert len(vs)==b['nvars'] and all(row[0] in rs for row in b['rows']),(mb,len(vs),b['nvars'],len(rs),len(b['rows']))
   assert all(row[2]<=len(vs) for row in b['rows'])
   indices=[i for i,r in enumerate(vs) if r.pos in (lo,hi)]
   calls={str(vs[i].pos)+':'+'>'.join(vs[i].alleles):Counter() for i in indices};pairs=[]
   for name,start,end,known in b['rows']:
    r=rs[name];aa=['NoOverlap']*len(vs);aa[start:end]=known
    for i in indices:calls[str(vs[i].pos)+':'+'>'.join(vs[i].alleles)][aa[i]]+=1
    left=[(vs[i].pos,vs[i].alleles,aa[i]) for i in indices if vs[i].pos==lo and aa[i] in ('Reference','Alternate')];right=[(vs[i].pos,vs[i].alleles,aa[i]) for i in indices if vs[i].pos==hi and aa[i] in ('Reference','Alternate')]
    if left and right:pairs.append({'read':r.query_name,'mapq':r.mapping_quality,'left':left,'right':right})
   print(mb,dict(calls),'paired',len(pairs))
   result.append({'chunk':mb,'beg':b['beg'],'end':b['end'],'variant_count':len(vs),'read_count':len(rs),'boundary_calls':calls,'bridges':pairs})
(root/'hiphase'/'boundary-calls.json').write_text(json.dumps(result,indent=2)+'\n')
