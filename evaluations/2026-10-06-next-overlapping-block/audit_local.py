#!/usr/bin/env python3
"""Audit complete local sequence contrasts on original primary alignments."""
from collections import Counter,defaultdict
from pathlib import Path
import json,pysam,importlib.util
P=Path(__file__).resolve().parent
spec=importlib.util.spec_from_file_location('a',P/'read_audit.py');a=importlib.util.module_from_spec(spec);spec.loader.exec_module(a)
truth={f[0]:f[1] for l in open('test_data/derived/chr20_truth_hap.tsv') if len(f:=l.rstrip().split('\t'))==2}
_,status=a.assignments(Path('test_data/tmp_gap_fix87/final/0/phased.bam'),{q:v=='PATERNAL' for q,v in truth.items()})
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
markers={19365841:['','TTTAT'],19373923:['','TTCC','TTCCTTCC'],19377346:['','ATATATAGAG','ATATATATAGAGAG']}
calls=defaultdict(dict);counts=defaultdict(Counter)
with pysam.AlignmentFile('test_data/HG002_chr20_hifi_mapped_to_CHM13_chr20_annotated.bam') as bam:
 for r in bam.fetch(chrom,19358000,19405000):
  if r.is_secondary or r.is_supplementary or r.is_duplicate or r.is_qcfail or not 30<=r.mapping_quality<255:continue
  for pos,alleles in markers.items():
   beg,end=pos-40,pos+160
   _,bq,qb=base(r,beg);_,eq,qe=base(r,end)
   if qb<0 or qe<=qb or min(bq,eq)<10 or max(bq,eq)==255:continue
   ref=fa.fetch(chrom,beg-1,end).upper();query=r.query_sequence[qb:qe+1].upper()
   ds=[distance(query,ref[:pos-beg]+alt+ref[pos-beg:]) for alt in alleles]
   net=len(query)-len(ref);winner=ds.index(min(ds))
   calls[pos][r.query_name]={'ds':ds,'net':net,'parent':truth.get(r.query_name),'status':status.get(r.query_name),'span':[r.reference_start+1,r.reference_end]}
   counts[(pos,'lengths')][net]+=1
   if ds.count(min(ds))>1 or min(ds)>3:continue
   counts[(pos,'parents')][(winner,truth.get(r.query_name))]+=1
   for sp in [19341833,19358995,19395544,19403172]:
    b,q,_=base(r,sp)
    if q>=20 and q!=255 and b:counts[(pos,sp)][(winner,b)]+=1
for key,count in counts.items():print(key,dict(count))
(P/'local-calls.json').write_text(json.dumps(calls,indent=2)+'\n')
