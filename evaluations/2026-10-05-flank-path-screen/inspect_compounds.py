#!/usr/bin/env python3
"""Reproduce the primary-alignment nomination screen recorded here."""
import pysam,collections,json
positions=((40633644,40636354),(45416750,45439920),(49800887,49845376),(57584041,57602614));out=[]
bam=pysam.AlignmentFile('test_data/HG002_chr20_hifi_mapped_to_CHM13_chr20_annotated.bam')
for a,z in positions:
 counts=collections.Counter();details=[]
 for read in bam.fetch('CHM13#0#chr20',a-1,z):
  if read.flag&(4|256|2048|512|1024) or not 30<=read.mapping_quality<255 or read.reference_start>=a or read.reference_end<z:continue
  aligned={p:q for q,p in read.get_aligned_pairs() if p is not None and q is not None}
  call=[];events=[];rp=read.reference_start+1;qp=0
  for op,n in read.cigartuples:
   if op==1 and any(abs(rp-(pos+1))<=16 for pos in(a,z)):
    events.append(('I',rp,read.query_sequence[qp:qp+n],min(read.query_qualities[max(0,qp-1):min(read.query_length,qp+n+1)])))
   if op==2 and any(abs(rp-(pos+1))<=16 for pos in(a,z)):
    events.append(('D',rp,n,min(read.query_qualities[max(0,qp-1):min(read.query_length,qp+1)])))
   if op in(0,7,8):rp+=n;qp+=n
   elif op in(2,3):rp+=n
   elif op in(1,4):qp+=n
  for p in(a,z):
   q=aligned.get(p-1)
   call.append((read.query_sequence[q],read.query_qualities[q]) if q is not None else None)
  details.append(dict(name=read.query_name,bases=call,events=events))
 out.append(dict(left=a,right=z,reads=details));print(a,z,len(details));
 for r in details:print(r)
open('test_data/tmp_gap_fix59/compound-evidence.json','w').write(json.dumps(out,indent=2)+'\n')
