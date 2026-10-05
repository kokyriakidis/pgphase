import pysam,json,collections
b=pysam.AlignmentFile('test_data/HG002_chr20_hifi_mapped_to_CHM13_chr20_annotated.bam')
for x in json.load(open('test_data/tmp_gap_fix51/competitor-gaps.json')):
 a,z=x['left'],x['right'];n=sum(not r.is_secondary and not r.is_supplementary and r.mapping_quality>=30 and r.reference_start<=a-1 and r.reference_end>=z for r in b.fetch('CHM13#0#chr20',a-1,z));print(a,z,n)
