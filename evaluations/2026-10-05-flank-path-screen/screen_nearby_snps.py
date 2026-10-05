#!/usr/bin/env python3
"""Reproduce the primary-alignment nomination screen recorded here."""
import collections,json,pysam,math
blocks=collections.defaultdict(list);snps=collections.defaultdict(list)
with pysam.VariantFile('test_data/tmp_gap_fix58/full-projected-final/phased.vcf') as v:
 for r in v:
  s=next(iter(r.samples.values()));g=s['GT'];p=s.get('PS')
  if not(s.phased and p and g and len(set(g))>1):continue
  blocks[p].append(r.pos)
  if r.info.get('CAT') in ('CLEAN_HET_SNP','NOISY_CAND_HET') and len(r.ref)==1 and len(r.alts)==1 and len(r.alts[0])==1:snps[p].append((r.pos,r.ref,r.alts[0],g[0]))
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
  ll=sorted(ll,reverse=True)[:3];rr=sorted(rr)[:3]
  pair_counts={(l,r):[0,0,0,0] for l in ll for r in rr};pair_lod=collections.Counter()
  for read in bam.fetch('CHM13#0#chr20',min(l[0] for l in ll)-1,max(r[0] for r in rr)):
   if read.flag&(4|256|2048|512|1024) or not 30<=read.mapping_quality<255:continue
   aligned={pos:qi for qi,pos in read.get_aligned_pairs() if pos is not None and qi is not None};calls={};quals={}
   for row in ll+rr:
    qi=aligned.get(row[0]-1)
    if qi is None or not 30<=read.query_qualities[qi]<255 or read.query_sequence[qi] not in row[1:3]:continue
    calls[row]=int((read.query_sequence[qi]==row[2])==(row[3]==1));quals[row]=read.query_qualities[qi]
   for (l,r),counts in pair_counts.items():
    if l not in calls or r not in calls:continue
    counts[calls[l]*2+calls[r]]+=1
    error=10**(-quals[l]/10)+10**(-quals[r]/10)+2*10**(-read.mapping_quality/10)
    pair_lod[(l,r)]+=(1 if calls[l]!=calls[r] else -1)*math.log((1-error)/error)
  for (l,r),counts in pair_counts.items():
   if l[0]==a and r[0]==z:continue
   x=dict(gap=(a,z),left_ps=p,right_ps=q,left=l,right=r,counts=counts,lod=pair_lod[(l,r)])
   if (not(counts[1]+counts[2]) and min(counts[0],counts[3])>=2 or not(counts[0]+counts[3]) and min(counts[1],counts[2])>=2) and math.ldexp(2.,-sum(counts))<=.01:
    out.append(x);print(x,flush=True)
open('test_data/tmp_gap_fix59/nearby-snp-evidence.json','w').write(json.dumps(out,indent=2)+'\n')
