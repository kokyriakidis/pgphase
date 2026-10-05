#!/usr/bin/env python3
"""Reproduce the primary-alignment nomination screen recorded here."""
import collections,json,pysam,math
blocks=collections.defaultdict(list);snps=collections.defaultdict(list)
with pysam.VariantFile('test_data/tmp_gap_fix58/full-projected-final/phased.vcf') as v:
 for r in v:
  s=next(iter(r.samples.values()));g=s['GT'];p=s.get('PS')
  if not(s.phased and p and g and len(set(g))>1):continue
  blocks[p].append(r.pos)
  if r.info.get('CAT')=='CLEAN_HET_SNP' and len(r.ref)==1 and len(r.alts)==1 and len(r.alts[0])==1:snps[p].append((r.pos,r.ref,r.alts[0],g[0]))
ordered=sorted((min(rr),max(rr),p) for p,rr in blocks.items());frontier=ordered[0];pairs=[]
for cur in ordered[1:]:
 if cur[0]>frontier[1]:pairs.append((frontier,cur))
 if cur[1]>frontier[1]:frontier=cur
out=[]
with pysam.AlignmentFile('test_data/HG002_chr20_hifi_mapped_to_CHM13_chr20_annotated.bam') as bam:
 for (_,a,p),(z,_,q) in pairs:
  if not 0<z-a<=50000:continue
  ll=[x for x in snps[p] if a-10000<=x[0]<=a];rr=[x for x in snps[q] if z<=x[0]<=z+10000]
  if not ll or not rr:continue
  l=max(ll);r=min(rr)
  if l[0]==a and r[0]==z:continue
  counts=[0,0,0,0];lod=0.
  for read in bam.fetch('CHM13#0#chr20',l[0]-1,r[0]):
   if read.flag&(4|256|2048|512|1024) or not 30<=read.mapping_quality<255:continue
   aligned={pos:qi for qi,pos in read.get_aligned_pairs() if pos is not None and qi is not None};calls=[];quals=[]
   for row in(l,r):
    qi=aligned.get(row[0]-1)
    if qi is None or not 30<=read.query_qualities[qi]<255 or read.query_sequence[qi] not in row[1:3]:break
    calls.append(int((read.query_sequence[qi]==row[2])==(row[3]==1)));quals.append(read.query_qualities[qi])
   if len(calls)!=2:continue
   counts[calls[0]*2+calls[1]]+=1
   error=sum(10**(-v/10) for v in quals)+2*10**(-read.mapping_quality/10)
   lod+=(1 if calls[0]!=calls[1] else -1)*math.log((1-error)/error)
  x=dict(gap=(a,z),left_ps=p,right_ps=q,left=l,right=r,counts=counts,lod=lod);out.append(x)
  if sum(counts) and (counts[0]+counts[3]==0 or counts[1]+counts[2]==0) and abs(lod)>=math.log(999):print(x,flush=True)
open('test_data/tmp_gap_fix59/outer-snp-evidence.json','w').write(json.dumps(out,indent=2)+'\n')
