import pysam,collections
b=pysam.AlignmentFile('test_data/HG002_chr20_hifi_mapped_to_CHM13_chr20_annotated.bam'); c=collections.Counter(); cq=collections.Counter()
for r in b.fetch('CHM13#0#chr20',7280345,7280357):
 if r.is_secondary or r.is_supplementary or r.mapping_quality<30: continue
 a={rp:q for q,rp in r.get_aligned_pairs() if rp in (7280345,7280355) and q is not None}
 if len(a)<2:continue
 bases=tuple(r.query_sequence[a[p]] for p in (7280345,7280355)); c[bases]+=1
 if min(r.query_qualities[a[p]] for p in (7280345,7280355))>=30:cq[bases]+=1
print(c,cq)
