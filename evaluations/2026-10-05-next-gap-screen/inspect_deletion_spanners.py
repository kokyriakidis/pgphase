import pysam
b=pysam.AlignmentFile('test_data/HG002_chr20_hifi_mapped_to_CHM13_chr20_annotated.bam')
for r in b.fetch('CHM13#0#chr20',11235200,11255400):
 if r.is_secondary or r.is_supplementary or r.reference_start>11235277 or r.reference_end<11255370:continue
 print(r.query_name,r.mapping_quality,r.reference_start,r.reference_end)
 pos=r.reference_start;q=0
 for op,n in r.cigartuples:
  if op in (1,2) and any(abs(pos-t)<12 for t in (11235277,11255369)):print('event',pos,op,n,r.query_sequence[max(0,q-5):q+n+5])
  if op in (0,2,3,7,8):pos+=n
  if op in (0,1,4,7,8):q+=n
