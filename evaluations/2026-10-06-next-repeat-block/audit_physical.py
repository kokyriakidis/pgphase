#!/usr/bin/env python3
"""Reconstruct the production cohorts from original CIGARs, without parental truth."""
import json,math
from pathlib import Path
import pysam
ROOT=Path(__file__).resolve().parents[2];OUT=Path(__file__).resolve().parent
chrom='CHM13#0#chr20';p=45866905;q=45883707
ref=pysam.FastaFile(str(ROOT/'test_data/chm13v2.0.chr20.renamed.fa'))
bam=pysam.AlignmentFile(str(ROOT/'test_data/HG002_chr20_hifi_mapped_to_CHM13_chr20_annotated.bam'))
alt=['ACAGACAGACACACACAC','ACAGACAGACACACACACAC']
profiles={}
for line in (ROOT/'test_data/tmp_gap_fix84/baseline/45/matrix.chunk0.recovery-final.tsv').open():
 f=line.rstrip().split('\t')
 if f[0]=='READ':profiles[f[1]]=(int(f[6]),int(f[7]))
def e(quality):return 10**(-quality/10)
def base(pos):return ref.fetch(chrom,pos-1,pos).upper()
def distance(a,b):
 d=list(range(len(b)+1))
 for i,x in enumerate(a,1):
  n=[i]
  for j,y in enumerate(b,1):n.append(min(n[-1]+1,d[j]+1,d[j-1]+(x!=y)))
  d=n
 return d[-1]
def left_call(r):
 beg,end=p-16,p+16;rp,qp=r.reference_start+1,0;qb=qe=None;net=0;error=0
 for op,n in r.cigartuples:
  if op in (0,7,8):
   if rp<=beg<rp+n:qb=qp+beg-rp
   if rp<=end<rp+n:qe=qp+end-rp+1
  if op==1 and beg<=rp<=end:
   net+=n
   if any(v==255 or (i>=len(alt[0]) and v<10) for i,v in enumerate(r.query_qualities[qp:qp+n])):return None
   error+=sum(e(v) for v in r.query_qualities[qp+len(alt[0]):qp+n])
  if op==2 and rp<end and rp+n>beg:net-=min(rp+n,end)-max(rp,beg)
  if op==3 and rp<end and rp+n>beg:return None
  if op in (0,1,4,7,8):qp+=n
  if op in (0,2,3,7,8):rp+=n
 if min(abs(net-len(a)) for a in alt)>2:return None
 if abs(net)>64 or qb is None or qe is None or qe<=qb:return None
 aq,bq=r.query_qualities[qb],r.query_qualities[qe-1]
 if min(aq,bq)<20 or 255 in (aq,bq):return None
 reference=ref.fetch(chrom,beg-1,end).upper();query=r.query_sequence[qb:qe].upper()
 ds=[distance(query,reference[:16]+a+reference[16:]) for a in alt]
 if ds[0]==ds[1]:return None
 allele=int(ds[1]<ds[0]);lengths=[abs(net-len(a)) for a in alt]
 if lengths[0]!=lengths[1] and int(lengths[1]<lengths[0])!=allele:return None
 if ds[allele]>lengths[allele]+1:return None
 return 1-allele,error+e(aq)+e(bq),net
beg=q
while base(beg-1)=='A' and beg>q-64:beg-=1
end=q
while base(end)=='A' and end<q+64:end+=1
def right_call(r):
 rp,qp=r.reference_start+1,0;net=events=0;eq=255;flanks=[0,0];pure=True
 for op,n in r.cigartuples:
  if op==1 and beg<=rp<=end:
   net+=n;events+=1;pure&=set(r.query_sequence[qp:qp+n])=={'A'};eq=min(eq,min(r.query_qualities[qp:qp+n]))
  if op==2 and rp<end and rp+n>beg:
   if rp<beg or rp+n>end or qp==0 or qp>=r.query_length:return None
   net-=n;events+=1;eq=min(eq,r.query_qualities[qp-1],r.query_qualities[qp])
  if op==3 and rp<end+16 and rp+n>beg-16:return None
  if op in (0,7,8):
   for pos in range(max(rp,beg),min(rp+n,end)):pure&=r.query_sequence[qp+pos-rp]=='A'
   for side,(a,b) in enumerate(((beg-16,beg),(end,end+16))):
    for pos in range(max(rp,a),min(rp+n,b)):
     qi=qp+pos-rp;quality=r.query_qualities[qi]
     if quality!=255 and r.query_sequence[qi]==base(pos):flanks[side]=max(flanks[side],quality)
  if op in (0,1,4,7,8):qp+=n
  if op in (0,2,3,7,8):rp+=n
 allele=int(net>0)
 if not pure or events>1 or abs(net-allele)>1 or min(flanks)<20:return None
 quality=min(flanks) if events==0 else eq
 return allele,sum(e(v) for v in flanks),quality
cohorts=[[[0,0],[0,0]],[[0,0],[0,0]]];classes=[[0,0],[0,0]];errors=[1.,1.];log_odds=0;records=[];seen=set()
for r in bam.fetch(chrom,p-1,q+1):
 name=r.query_name
 if r.flag&(4|256|2048|1024|512) or r.mapping_quality<30 or r.mapping_quality==255 or name in seen:continue
 seen.add(name)
 if name not in profiles:continue
 a,b=left_call(r),right_call(r);mapping=2*e(r.mapping_quality)
 ah=a[0] if a and a[1]+mapping<=.01 else None
 bh=b[0] if b and 20<=b[2]<255 and b[1]+mapping+e(b[2])<=.01 else None
 h,ps=profiles[name];kind=None
 if ah is not None and bh is not None:
  parity=int(ah!=bh);classes[bh][parity]+=1;error=a[1]+b[1]+mapping+e(b[2]);errors[bh]*=error;log_odds+=(1 if parity else -1)*math.log((1-error)/error);kind='bridge'
 elif ah is not None and r.reference_end<q and ps==45854428 and h in(1,2):cohorts[0][ah][h-1]+=1;kind='left_calibration'
 elif bh is not None and b[2]>=30 and r.reference_start+1>p and ps==45883706 and h in(1,2):cohorts[1][bh][h-1]+=1;kind='right_calibration'
 if kind:records.append({'qname':name,'cohort':kind,'left':a,'right':b,'mapq':r.mapping_quality})
def wilson(d,n):
 z=1.6448536269514722;rate=d/n
 return (rate+z*z/(2*n)+z*math.sqrt(rate*(1-rate)/n+z*z/(4*n*n)))/(1+z*z/n)
result={'left_joint_error_bound':wilson(3,44)+.01+1/(1+math.exp(log_odds)),
        'right_diploid_joint_error_bound':(wilson(1,16)+.01+errors[0])*(wilson(1,5)+.01+errors[1]),
        'gauges':cohorts,'class_parity':classes,'class_error_products':errors,'log_odds':log_odds,'molecules':records,'zero_length_compound_calibrators':[x for x in records if x['left'] and x['left'][2]==0]}
assert cohorts==[[[27,0],[3,14]],[[15,1],[1,4]]],cohorts
assert classes==[[0,1],[0,3]],classes
(OUT/'physical-evidence.json').write_text(json.dumps(result,indent=2)+'\n');print({k:v for k,v in result.items() if k!='molecules'})
