#!/usr/bin/env python3
"""Measure original physical repeat classes and record independent HiPhase calls."""
from collections import Counter
import json
from pathlib import Path
import pysam

ROOT = Path(__file__).resolve().parents[2]
OUT = Path(__file__).resolve().parent
b = pysam.AlignmentFile(str(ROOT/'test_data/HG002_chr20_hifi_mapped_to_CHM13_chr20_annotated.bam'))
f = pysam.FastaFile(str(ROOT/'test_data/chm13v2.0.chr20.renamed.fa'))
c = 'CHM13#0#chr20'
truth = {z[0]: z[1] for line in (ROOT/'test_data/derived/chr20_truth_hap.tsv').read_text().splitlines()
         if len(z := line.split('\t')) == 2}

def call(r,p,delta):
 base=f.fetch(c,p-1,p).upper();beg=p;end=p
 while beg>p-64 and f.fetch(c,beg-2,beg-1).upper()==base:beg-=1
 while end<p+64 and f.fetch(c,end-1,end).upper()==base:end+=1
 if end==p+64 or beg==p-64:return
 rp=r.reference_start+1;qi=0;net=0;ev=0;seq=r.query_sequence;q=r.query_qualities;fl=[0,0];pure=True
 for op,n in r.cigartuples:
  if op==1 and beg<=rp<=end:
   net+=n;ev+=1;pure&=all(z==base for z in seq[qi:qi+n])
  if op==2 and rp<end and rp+n>beg:
   if rp<beg or rp+n>end:return
   net-=n;ev+=1
  if op==3 and rp<end+16 and rp+n>beg-16:return
  if op in [0,7,8]:
   for x in range(max(rp,beg),min(rp+n,end)):pure&=seq[qi+x-rp]==base
   for side,(lo,hi) in enumerate([(beg-16,beg),(end,end+16)]):
    for x in range(max(rp,lo),min(rp+n,hi)):
     if seq[qi+x-rp]==f.fetch(c,x-1,x).upper() and q[qi+x-rp]!=255:fl[side]=max(fl[side],q[qi+x-rp])
  if op in [0,1,4,7,8]:qi+=n
  if op in [0,2,3,7,8]:rp+=n
 allele=int(net*delta>0)
 if ev>1 or not pure or min(fl)<20 or abs(net-allele*delta)>1:return
 return allele

counts = [Counter(), Counter()]
bridge = Counter()
seen = set()
for r in b.fetch(c, 12717000, 12754000):
    if r.is_secondary or r.is_supplementary or r.is_unmapped or r.is_duplicate or r.is_qcfail or r.mapping_quality < 20 or r.mapping_quality == 255:
        continue
    if r.query_name in seen:
        continue
    seen.add(r.query_name)
    a, d = call(r, 12721113, 1), call(r, 12735895, -1)
    for i, allele in enumerate((a, d)):
        if allele is not None:
            counts[i][(truth.get(r.query_name, 'UNSCORABLE'), allele)] += 1
    if a is not None and d is not None:
        bridge[a, d] += 1
hi = []
v = pysam.VariantFile(str(ROOT/'test_data/tmp_gap_fix48/competitor/hiphase_dv/phased.vcf.gz'))
for r in v.fetch('chr20', 12721000, 12752500):
    if r.pos in (12721112, 12735894, 12740002, 12752291):
        hi.append({'pos': r.pos, 'ref': r.ref, 'alt': r.alts, 'sample': dict(next(iter(r.samples.values())))})
result = {'original_repeat_calls': [{f'{parent}:{allele}': n for (parent, allele), n in x.items()} for x in counts],
          'interior_pairs': {f'{a}:{d}': n for (a, d), n in bridge.items()},
          'hiphase_calls': hi, 'truth_used_by_production': False}
for name, path in [('native', 'final_owner/12'), ('two_chunk', 'continuation4/11')]:
    result[name + '_trace'] = [line for line in (ROOT/'test_data/tmp_gap_fix79'/path/'stderr.log').read_text().splitlines()
                             if '[repeat-indel-chain]' in line]
(OUT/'physical-evidence.json').write_text(json.dumps(result, indent=2) + '\n')
print(json.dumps(result, indent=2))
