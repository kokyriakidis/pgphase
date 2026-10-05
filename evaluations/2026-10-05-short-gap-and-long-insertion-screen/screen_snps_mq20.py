import collections,json,pysam
blocks=collections.defaultdict(list)
with pysam.VariantFile('test_data/tmp_gap_fix55/full-final/phased.vcf') as v:
 for r in v:
  s=next(iter(r.samples.values()));g=s['GT'];p=s.get('PS')
  if s.phased and p and g and len(set(g))>1:
   blocks[p].append((r.pos,r.ref,','.join(r.alts),g,r.info.get('CAT')))
ordered=sorted((min(x[0] for x in rr),max(x[0] for x in rr),p) for p,rr in blocks.items())
bam=pysam.AlignmentFile('test_data/HG002_chr20_hifi_mapped_to_CHM13_chr20_annotated.bam');results=[]
frontier=ordered[0];pairs=[]
for current in ordered[1:]:
 if current[0]>frontier[1]:pairs.append((frontier,current))
 if current[1]>frontier[1]:frontier=current
for (_,a,p),(z,_,q) in pairs:
 if not 0<z-a<=50000:continue
 l=[x for x in blocks[p] if x[0]==a];r=[x for x in blocks[q] if x[0]==z]
 if len(l)!=1 or len(r)!=1:continue
 l,r=l[0],r[0]
 if any(len(x[1])!=1 or len(x[2])!=1 for x in(l,r)):continue
 counts=[0,0,0,0];names=set()
 for read in bam.fetch('CHM13#0#chr20',a-1,z):
  if (read.is_secondary or read.is_supplementary or read.is_duplicate or read.is_qcfail or read.mapping_quality<20 or read.mapping_quality==255 or read.reference_start>a-1 or read.reference_end<z or read.query_name in names):continue
  rp=read.reference_start;qp=0;calls={}
  for op,n in read.cigartuples:
   if op in(0,7,8):
    for x in(l,r):
     if rp<=x[0]-1<rp+n:
      qi=qp+x[0]-1-rp;b=read.query_sequence[qi];qual=read.query_qualities[qi]
      if 20<=qual<255 and b in(x[1],x[2]):calls[x[0]]=int(b==x[2])
    rp+=n;qp+=n
   elif op in(2,3):rp+=n
   elif op in(1,4):qp+=n
   if rp>=z:break
  if len(calls)!=2:continue
  names.add(read.query_name);counts[calls[a]*2+calls[z]]+=1
 same=counts[0]+counts[3];cross=counts[1]+counts[2]
 result=dict(left=a,right=z,left_row=l,right_row=r,left_ps=p,right_ps=q,counts=counts,left_sites=len(blocks[p]),right_sites=len(blocks[q]))
 results.append(result)
 if min(same,cross)==0 and max(same,cross)>=2 and ((counts[0]>0 and counts[3]>0) or (counts[1]>0 and counts[2]>0)):print(result,flush=True)
open('test_data/tmp_gap_fix57/snp-screen-mq20.json','w').write(json.dumps(results,indent=2)+'\n')
