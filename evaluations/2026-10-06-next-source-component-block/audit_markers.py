from pathlib import Path
from collections import Counter,defaultdict
import pysam,json,importlib.util
P=Path(__file__).resolve().parent
spec=importlib.util.spec_from_file_location('a',P/'read_audit.py');a=importlib.util.module_from_spec(spec);spec.loader.exec_module(a)
truth={f[0]:f[1] for l in open('test_data/derived/chr20_truth_hap.tsv') if len(f:=l.rstrip().split('\t'))==2}
_,status=a.assignments(Path('test_data/tmp_gap_fix89/baseline/46/phased.bam'),{q:v=='PATERNAL' for q,v in truth.items()})
fa=pysam.FastaFile('test_data/chm13v2.0.chr20.renamed.fa');chrom='CHM13#0#chr20'
def distance(a,b):
 row=list(range(len(b)+1))
 for i,x in enumerate(a,1):
  nxt=[i]
  for j,y in enumerate(b,1):nxt.append(min(nxt[-1]+1,row[j]+1,row[j-1]+(x!=y)))
  row=nxt
 return row[-1]
def base(r,p):
 rp,qp=r.reference_start+1,0
 for op,n in r.cigartuples:
  if op in(0,7,8) and rp<=p<rp+n:
   q=qp+p-rp;return r.query_sequence[q].upper(),r.query_qualities[q],q
  if op in(0,1,4,7,8):qp+=n
  if op in(0,2,3,7,8):rp+=n
 return '',255,-1
markers={46389865:(0,'A'),46402980:(1,''),46405049:(0,'A')}
report={}
with pysam.AlignmentFile('test_data/HG002_chr20_hifi_mapped_to_CHM13_chr20_annotated.bam') as bam:
 reads=[r for r in bam.fetch(chrom,46378000,46413000) if not r.is_secondary and not r.is_supplementary and not r.is_duplicate and not r.is_qcfail and r.mapping_quality!=255]
 for pos,(dl,alt) in markers.items():
  motif=alt or fa.fetch(chrom,pos-1,pos-1+dl).upper();beg=pos;end=pos+dl
  while pos-beg<64 and fa.fetch(chrom,beg-2,beg-1).upper()==motif[(beg-pos-1)%len(motif)]:beg-=1
  while end-pos<64 and fa.fetch(chrom,end-1,end).upper()==motif[(end-pos)%len(motif)]:end+=1
  beg-=16;end+=16
  ref=fa.fetch(chrom,beg-1,end).upper();alts=[ref,ref[:pos-beg]+alt+ref[pos-beg+dl:]]
  calls={};counts=Counter();gauge=defaultdict(Counter)
  for r in reads:
   _,bq,qb=base(r,beg);_,eq,qe=base(r,end)
   if qb<0 or qe<=qb or min(bq,eq)<10 or max(bq,eq)==255:continue
   query=r.query_sequence[qb:qe+1].upper();ds=[distance(query,s) for s in alts];net=len(query)-len(ref)
   lengths=[0,len(alt)-dl];near=[abs(net-n) for n in lengths];winner=ds.index(min(ds))
   if ds.count(min(ds))>1 or near.count(min(near))>1 or near.index(min(near))!=winner or min(near)>2 or min(ds)>min(near)+1:continue
   err=2*10**(-r.mapping_quality/10)+10**(-bq/10)+10**(-eq/10)
   q=r.query_name;parent=truth.get(q)
   calls[q]={'allele':winner,'distances':ds,'net':net,'error':err,'MQ':r.mapping_quality,'parent':parent,'status':status.get(q),'span':[r.reference_start+1,r.reference_end]}
   counts[(winner,parent)]+=1
   if r.mapping_quality<30 or err>.01:continue
   for sp in [46384021,46407462]:
    bb,sq,_=base(r,sp)
    if sq>=20 and sq!=255 and err+10**(-sq/10)<=.01 and bb:
     if (sp<pos and r.reference_end<46407462) or (sp>pos and r.reference_start+1>46384021):gauge[sp][(winner,bb)]+=1
  print(pos,'parents',dict(counts),'gauges',{sp:dict(v) for sp,v in gauge.items()})
  report[str(pos)]={'region':[beg,end],'calls':calls,'gauges':{str(sp):{str(k):v for k,v in c.items()} for sp,c in gauge.items()}}
(P/'physical-calls.json').write_text(json.dumps(report,indent=2)+'\n')
