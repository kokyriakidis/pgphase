#!/usr/bin/env python3
"""Reproduce the original-base deletion certificate without parental truth."""
import pysam,collections,math,json
from pathlib import Path
ref=pysam.FastaFile('test_data/chm13v2.0.chr20.renamed.fa');chrom='CHM13#0#chr20';pos=52715882;length=4
OUT=Path(__file__).resolve().parent
calibration=[]
sites=[]
for r in pysam.VariantFile('test_data/tmp_gap_fix74/baseline/52/phased.vcf'):
 s=next(iter(r.samples.values()))
 if s.get('PS')==52711825 and r.info.get('CAT')=='CLEAN_HET_SNP' and len(r.ref)==1 and len(r.alts or())==1 and len(r.alts[0])==1 and s.phased and len(set(s['GT']))==2:sites.append((r.pos,r.ref,r.alts[0],s['GT'][0]))
def rb(p):return ref.fetch(chrom,p-1,p).upper()
rows=[];gauge=collections.Counter();parity=collections.Counter()
bam=pysam.AlignmentFile('test_data/HG002_chr20_hifi_mapped_to_CHM13_chr20_annotated.bam')
for r in bam.fetch(chrom,52696939,pos+length):
 if r.flag&(4|256|2048|1024|512) or r.mapping_quality<30 or r.mapping_quality==255:continue
 bases={};events=[];rp=r.reference_start+1;qp=0
 for op,n in r.cigartuples:
  if op in(0,7,8):
   for p in range(max(rp,pos-33),min(rp+n,pos+length+34)):
    q=qp+p-rp;bases[p]=(r.query_sequence[q],r.query_qualities[q])
   for p in [52696940]+[x[0] for x in sites if rp<=x[0]<rp+n]:
    if rp<=p<rp+n:q=qp+p-rp;bases[p]=(r.query_sequence[q],r.query_qualities[q])
  if op in(1,2):events.append((rp,op,n))
  if op in(0,7,8,2,3):rp+=n
  if op in(0,7,8,1,4):qp+=n
 selected=[]
 for p,op,n in events:
  if op!=2 or n!=length or abs(p-pos)>32:continue
  beg,end=min(pos,p),max(pos,p)+length;s=ref.fetch(chrom,beg-1,end-1)
  if s[:pos-beg]+s[pos-beg+length:]==s[:p-beg]+s[p-beg+length:]:selected.append(p)
 if len(selected)>1:continue
 observed=selected[0] if selected else pos;beg,end=min(pos,observed),max(pos,observed)+length
 if any((op==1 and beg-1<=p<=end) or (op==2 and p<=end and p+n>beg-1) for p,op,n in events if not(selected and p==observed and op==2 and n==length)):continue
 checks=[p for p in range(beg-1,end+1) if not(selected and observed<=p<observed+length)]
 if any(p not in bases or bases[p][0]!=rb(p) or not 10<=bases[p][1]<255 for p in checks):continue
 error=sum(10**(-bases[p][1]/10) for p in checks)+2*10**(-r.mapping_quality/10);allele=int(bool(selected));hap=1 if allele else 2
 votes=[];snperror=0
 for sp,sr,sa,sh in sites:
  if sp not in bases:continue
  base,q=bases[sp]
  if 20<=q<255 and base in(sr,sa):votes.append(1 if int(base==sa)==sh else 2);snperror+=10**(-q/10)
 if r.reference_start>=52696940 and votes and error+snperror<=.01 and all(bases[p][1]>=20 for p in checks):
  observed=votes[0] if len(set(votes))==1 else 3-hap
  gauge[hap,observed]+=1
  calibration.append(dict(qname=r.query_name,deletion_hap=hap,snp_hap=observed,error=error+snperror))
 left=bases.get(52696940)
 if left and left[0] in('G','A') and 10<=left[1]<255:
  le=error+10**(-left[1]/10);flip=(left[0]=='A')!=(hap==1)
  row=dict(qname=r.query_name,left=left,allele=allele,anchors=[bases[p] for p in checks],error=le,flip=flip)
  rows.append(row)
names={r['qname'] for r in rows}
assert names.isdisjoint(r['qname'] for r in calibration)
assert gauge[1,1]>=2 and gauge[2,2]>=2 and gauge[1,2]==gauge[2,1]==0
association=2*0.5**len(calibration)
quality_error=sum(r['error'] for r in rows+calibration)
assert association<=.01 and rows and len({r['flip'] for r in rows})==1
assert quality_error<=.20
report={'bridge_pairs':rows,'calibration':calibration,'gauge':[[gauge[1,1],gauge[1,2]],[gauge[2,1],gauge[2,2]]],
        'association_p':association,'joint_base_mapping_error_bound':quality_error,
        'disjoint_original_molecules':True,'parental_truth_loaded':False}
(OUT/'physical-certificate.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps(report,indent=2))
