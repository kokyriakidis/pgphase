from pathlib import Path
import pysam,json
root=Path('test_data/tmp_gap_fix61');c=json.loads((Path('evaluations/2026-10-05-short-insertion-source-path')/'hiphase-comparison.json').read_text());names=set(c['hiphase_correct_pgphase_not']+c['pgphase_correct_hiphase_not']);result={q:{} for q in names}
for label,path in [('pg','test_data/tmp_gap_fix61/bam-only/33/phased.bam'),('hp','test_data/tmp_gap_fix48/competitor/hiphase_dv/phased.bam')]:
 with pysam.AlignmentFile(path,check_sq=False) as b:
  it=b.fetch('chr20',33747590,33749688) if label=='hp' else b
  for r in it:
   if r.query_name not in names or r.is_secondary or r.is_supplementary:continue
   result[r.query_name][label]={'hp':r.get_tag('HP') if r.has_tag('HP') else 0,'ps':r.get_tag('PS') if r.has_tag('PS') else 0}
with pysam.AlignmentFile('test_data/HG002_chr20_hifi_mapped_to_CHM13_chr20_annotated.bam') as b:
 for r in b.fetch('CHM13#0#chr20',33747590,33749688):
  if r.query_name not in names or r.is_secondary or r.is_supplementary:continue
  result[r.query_name]['alignment']={'beg':r.reference_start+1,'end':r.reference_end,'mapq':r.mapping_quality}
for q in names:
 result[q]['final_matrix']=[l for l in (root/'bam-only/33/matrix.chunk0.recovery-final.tsv').read_text().splitlines() if '\t'+q+'\t' in l]
 result[q]['source_matrix']=[l for l in (root/'bam-only/33/matrix.recovery.chunk0.window2.chunk-1.recovery-source.tsv').read_text().splitlines() if '\t'+q+'\t' in l]
 result[q]['source_calls']=[l for l in (root/'bam-only/33/matrix.chunk0.bam-source-evidence.tsv').read_text().splitlines() if '\t'+q+'\t' in l]
(Path('evaluations/2026-10-05-short-insertion-source-path')/'read-disagreements.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result,indent=2))
