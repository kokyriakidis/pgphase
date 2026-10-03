from collections import Counter
import json,pysam
truth={f[0]:f[1] for line in open('test_data/derived/chr20_truth_hap.tsv') if len(f:=line.rstrip().split('\t'))==2}
out={}
with pysam.AlignmentFile('test_data/HG002_chr20_hifi_mapped_to_CHM13_chr20_annotated.bam') as bam:
 for pos in (50562679,50719983):
  counts=Counter();examples=[]
  for read in bam.fetch('CHM13#0#chr20',pos-1,pos):
   if read.is_secondary or read.is_supplementary or read.mapping_quality<30:continue
   ref=read.reference_start;query=0;base=None;insertion=''
   for op,n in read.cigartuples:
    if op in (0,7,8):
     if ref<=pos-1<ref+n:
      qi=query+pos-1-ref
      if read.query_qualities[qi]>=30:base=read.query_sequence[qi]
     ref+=n;query+=n
    elif op==1:
     if ref==pos:insertion=read.query_sequence[query:query+n]
     query+=n
    elif op in (2,3):ref+=n
    elif op==4:query+=n
   if base and read.query_name in truth:counts[(truth[read.query_name],base,insertion)]+=1
  out[pos]=[{'parent':p,'base':b,'following_insertion':i,'reads':n} for (p,b,i),n in counts.items()]
open('evaluations/2026-10-03-complementary-deletion-backfill/changed-site-cigar-audit.json','w').write(json.dumps(out,indent=2)+'\n')
print(json.dumps(out,indent=2))
