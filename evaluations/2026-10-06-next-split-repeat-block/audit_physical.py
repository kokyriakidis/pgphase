#!/usr/bin/env python3
"""Reproduce complete-repeat calls and independent SNP calibration from original alignments."""
import pysam,collections,json
bam=pysam.AlignmentFile('test_data/HG002_chr20_hifi_mapped_to_CHM13_chr20_annotated.bam');fa=pysam.FastaFile('test_data/chm13v2.0.chr20.renamed.fa');chrom='CHM13#0#chr20'
truth={f[0]:f[1] for l in open('test_data/derived/chr20_truth_hap.tsv') if len(f:=l.rstrip().split('\t'))==2}
def ref(p,n=1):return fa.fetch(chrom,p-1,p-1+n).upper()
def base(r,p):
 rp,qp=r.reference_start+1,0
 for op,n in r.cigartuples:
  if op in(0,7,8) and rp<=p<rp+n:
   q=qp+p-rp;return r.query_sequence[q],r.query_qualities[q],q
  if op in(0,1,4,7,8):qp+=n
  if op in(0,2,3,7,8):rp+=n
 return '',255,-1
# Sequence and nearest signed-length agreement; complete bounded repeat.
def distance(a,b):
 row=list(range(len(b)+1))
 for i,x in enumerate(a,1):
  nxt=[i]
  for j,y in enumerate(b,1):nxt.append(min(nxt[-1]+1,row[j]+1,row[j-1]+(x!=y)))
  row=nxt
 return row[-1]
def deletion(r,p,lengths):
 s=ref(p,max(lengths));period=next((i for i in(1,2) if len(s)%i==0 and s==s[:i]*(len(s)//i)),None)
 if period is None:return
 left,right=p,p+max(lengths)
 while p-left<64 and ref(left-1)==ref(left-1+period):left-=1
 while right-p<64 and ref(right)==ref(right-period):right+=1
 if p-left==64 or right-p==64:return
 beg,end=left-16,right+16;_,bq,qb=base(r,beg);_,eq,qe=base(r,end)
 if qb<0 or qe<=qb or min(bq,eq)<20 or max(bq,eq)==255:return
 rp,qp,net=r.reference_start+1,0,0
 for op,n in r.cigartuples:
  if op==1 and beg<=rp<=end:net-=n
  if op==2 and rp<end and rp+n>beg:net+=min(rp+n,end)-max(rp,beg)
  if op==3 and rp<end and rp+n>beg:return
  if op in(0,1,4,7,8):qp+=n
  if op in(0,2,3,7,8):rp+=n
 near=[abs(net-n) for n in lengths]
 if near[0]==near[1]:return
 a=int(near[1]<near[0]);reference=ref(beg,end-beg+1);query=r.query_sequence[qb:qe+1]
 ds=[distance(query,reference[:p-beg]+reference[p-beg+n:]) for n in lengths]
 if ds[a]>=ds[1-a] or ds[a]>near[a]+1 or near[a]>2:return
 return a,10**(-bq/10)+10**(-eq/10),net,ds
if __name__=='__main__':
 pairs=[(56007501,[4,6]),(56026465,[1,2]),(56027379,[1,2])];snps=[55999194,55999561,56040612]
 calls={p:{} for p,ls in pairs};counts=collections.defaultdict(collections.Counter)
 for r in bam.fetch(chrom,55998000,56042000):
  if r.is_secondary or r.is_supplementary or r.is_duplicate or r.is_qcfail or not 30<=r.mapping_quality<255:continue
  for p,ls in pairs:
   c=deletion(r,p,ls)
   if c is None:continue
   a,e,net,ds=c;e+=2*10**(-r.mapping_quality/10)
   if e>.01:continue
   calls[p][r.query_name]=(a,e,net,ds)
   if r.query_name in truth:counts[(p,'parents')][(a,truth[r.query_name])]+=1
   for snp in snps:
    allele,q,_=base(r,snp)
    if q>=20 and q!=255 and allele and e+10**(-q/10)<=.01:counts[(p,snp)][(a,allele)]+=1
 for (p,l),(q,m) in zip(pairs,pairs[1:]):
  for name in calls[p].keys()&calls[q].keys():counts[(p,q)][(calls[p][name][0],calls[q][name][0])]+=1
 for k,v in counts.items():print(k,dict(v))
