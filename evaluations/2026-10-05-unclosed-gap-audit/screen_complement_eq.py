import pysam,collections,json
v=pysam.VariantFile('test_data/tmp_gap_fix55/full-final/phased.vcf');blocks=collections.defaultdict(list)
for r in v:
 s=next(iter(r.samples.values()));p=s.get('PS')
 if p and s.phased and len(set(s['GT']))>1:blocks[p].append((r.pos,r.ref,','.join(r.alts),s['GT']))
a=sorted((min(r[0] for r in rr),max(r[0] for r in rr),p) for p,rr in blocks.items());bam=pysam.AlignmentFile('test_data/HG002_chr20_hifi_mapped_to_CHM13_chr20_annotated.bam');fa=pysam.FastaFile('test_data/chm13v2.0.chr20.renamed.fa');out=[]
frontier=a[0];pairs=[]
for current in a[1:]:
 if current[0]>frontier[1]:pairs.append((frontier,current))
 if current[1]>frontier[1]:frontier=current
for (_,x,p),(y,_,q) in pairs:
 if not 0<y-x<=50000:continue
 ls=[z for z in blocks[p] if z[0]==x];rs=[z for z in blocks[q] if z[0]==y]
 if len(ls)==2 and len(rs)==1:pair,snp,side=ls,rs[0],0
 elif len(rs)==2 and len(ls)==1:pair,snp,side=rs,ls[0],1
 else:continue
 if len(snp[1])!=1 or len(snp[2])!=1 or any(len(z[1])!=1 or len(z[2])<=1 or ',' in z[2] for z in pair):continue
 pos=pair[0][0];lo=pos-17;hi=pos+17;ref=fa.fetch('CHM13#0#chr20',lo,hi);expected=[ref[:pos-lo]+z[2][1:]+ref[pos-lo:] for z in pair];counts=[0]*4;names=set()
 for read in bam.fetch('CHM13#0#chr20',x-1,y):
  if read.is_secondary or read.is_supplementary or read.is_duplicate or read.is_qcfail or read.mapping_quality<30 or read.mapping_quality==255 or read.reference_start>x-1 or read.reference_end<y or read.query_name in names:continue
  rp=read.reference_start;qp=0;call=-1;events=[]
  for op,n in read.cigartuples:
   if op in (0,7,8):
    if rp<=snp[0]-1<rp+n:
     qi=qp+snp[0]-1-rp;b=read.query_sequence[qi]
     if read.query_qualities[qi]>=30 and b in (snp[1],snp[2]):call=int(b==snp[2])
    rp+=n;qp+=n
   elif op==1:
    if lo<=rp<hi and min(read.query_qualities[max(0,qp-1):qp+n+1])>=30:events.append((rp,read.query_sequence[qp:qp+n]))
    qp+=n
   elif op in (2,3):rp+=n
   elif op==4:qp+=n
  if call<0 or len(events)!=1:continue
  ep,seq=events[0];actual=ref[:ep-lo]+seq+ref[ep-lo:];matches=[i for i,z in enumerate(expected) if z==actual]
  if len(matches)!=1:continue
  names.add(read.query_name);counts[matches[0]*2+call]+=1
 d=dict(left=x,right=y,left_ps=p,right_ps=q,pair=pair,snp=snp,pair_side=side,counts=counts);out.append(d)
 if min(counts[0]+counts[3],counts[1]+counts[2])==0 and ((min(counts[0],counts[3])>=1) or min(counts[1],counts[2])>=1):print(d,flush=True)
open('test_data/tmp_gap_fix56/complement-eq-screen.json','w').write(json.dumps(out,indent=2)+'\n')
