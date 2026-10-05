import pysam,collections,json
v=pysam.VariantFile('test_data/tmp_gap_fix55/full-final/phased.vcf');blocks=collections.defaultdict(list)
for r in v:
 s=next(iter(r.samples.values()));p=s.get('PS')
 if p and s.phased and len(set(s['GT']))>1:blocks[p].append((r.pos,r.ref,','.join(r.alts),s['GT']))
a=sorted((min(r[0] for r in rr),max(r[0] for r in rr),p) for p,rr in blocks.items());bam=pysam.AlignmentFile('test_data/HG002_chr20_hifi_mapped_to_CHM13_chr20_annotated.bam');out=[]
frontier=a[0];pairs=[]
for current in a[1:]:
 if current[0]>frontier[1]:pairs.append((frontier,current))
 if current[1]>frontier[1]:frontier=current
for (_,x,p),(y,_,q) in pairs:
 if not 0<y-x<=50000:continue
 ls=[z for z in blocks[p] if z[0]==x];rs=[z for z in blocks[q] if z[0]==y]
 if len(ls)!=1 or len(rs)!=1:continue
 l,r=ls[0],rs[0]
 if all(len(z[1])==len(z[2])==1 for z in (l,r)):continue
 if any(',' in z[2] or (len(z[1])>1 and len(z[2])>1) for z in (l,r)):continue
 counts=[0]*4
 for read in bam.fetch('CHM13#0#chr20',x-1,y):
  if read.is_secondary or read.is_supplementary or read.is_duplicate or read.is_qcfail or read.mapping_quality<30 or read.mapping_quality==255 or read.reference_start>x-1 or read.reference_end<y:continue
  rp=read.reference_start;qp=0;bases={};events={}
  for op,n in read.cigartuples:
   if op in (0,7,8):
    for z in (l,r):
     for pos in range(z[0]-1,z[0]+len(z[1])):
      if rp<=pos<rp+n:
       qi=qp+pos-rp;bases[pos]=(read.query_sequence[qi],read.query_qualities[qi])
    rp+=n;qp+=n
   elif op==1:
    if rp in (x,y):events[rp]=('I',read.query_sequence[qp:qp+n],min(read.query_qualities[max(0,qp-1):qp+n+1]))
    qp+=n
   elif op in (2,3):
    if rp in (x,y):events[rp]=('D',n, min(read.query_qualities[max(0,qp-1):qp+1]))
    rp+=n
   elif op==4:qp+=n
  calls=[]
  for z in (l,r):
   pos,ref,alt,gt=z;call=-1
   if len(ref)==len(alt)==1:
    b,bq=bases.get(pos-1,('',0));call=int(b==alt) if bq>=30 and b in (ref,alt) else -1
   elif pos in events:
    op,seq,bq=events[pos]
    if bq>=30 and ((op=='I' and len(ref)==1 and seq==alt[1:]) or (op=='D' and len(alt)==1 and seq==len(ref)-1)):call=1
   elif all(bases.get(pos-1+i,('',0))[0]==b and bases.get(pos-1+i,('',0))[1]>=30 for i,b in enumerate(ref)) and bases.get(pos,('',0))[1]>=30:call=0
   calls.append(call)
  if min(calls)<0:continue
  counts[calls[0]*2+calls[1]]+=1
 same=counts[0]+counts[3];cross=counts[1]+counts[2]
 d=dict(left=l,right=r,left_ps=p,right_ps=q,counts=counts);out.append(d)
 if min(same,cross)==0 and max(same,cross)>=4 and ((min(counts[0],counts[3])>=2) or min(counts[1],counts[2])>=2):print(d,flush=True)
open('test_data/tmp_gap_fix56/indel-screen.json','w').write(json.dumps(out,indent=2)+'\n')
