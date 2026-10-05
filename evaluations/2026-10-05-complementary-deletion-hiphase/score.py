"""Compare read truth and actual boundary rows in matched owning-region solves."""
from pathlib import Path
from collections import Counter, defaultdict
import json
import pysam
root=Path('test_data/tmp_deletion_hiphase_audit')
truth={f[0]:f[1] for line in open('test_data/derived/chr20_truth_hap.tsv') if len(f:=line.rstrip().split('\t'))==2}
left,right=50548245,50562066
with pysam.AlignmentFile(str(root/'owner/input.bam')) as bam:
 eligible={r.query_name for r in bam if not(r.is_secondary or r.is_supplementary) and r.reference_start<=right-1 and r.reference_end>left-1}
outputs=[('pgphase_owner',Path('test_data/tmp_gap_fix50/deletion-hiphase/50/phased.bam'),Path('test_data/tmp_gap_fix50/deletion-hiphase/50/phased.vcf')),
         ('pgphase_chr20',Path('test_data/tmp_gap_fix50/full-final/phased.bam'),Path('test_data/tmp_gap_fix50/full-final/phased.vcf')),
         ('hiphase_owner',root/'owner/phased.bam',root/'owner/phased.vcf.gz'),
         ('hiphase_short',root/'short/phased.bam',root/'short/phased.vcf.gz'),
         ('hiphase_local',root/'short/local-only.bam',root/'short/local-only.vcf.gz')]
result=[]
for stage,bam_path,vcf_path in outputs:
 all_votes=defaultdict(Counter);gap_votes=defaultdict(Counter);tags={}
 with pysam.AlignmentFile(str(bam_path),check_sq=False) as bam:
  for r in bam:
   q=r.query_name
   if r.is_secondary or r.is_supplementary or q in tags:continue
   hp=r.get_tag('HP') if r.has_tag('HP') else 0;ps=r.get_tag('PS') if r.has_tag('PS') else 0
   tags[q]=(hp,ps)
   if hp not in [1,2] or ps<=0 or q not in truth:continue
   relation=(hp==1)==(truth[q]=='PATERNAL');all_votes[ps][relation]+=1
   if q in eligible:gap_votes[ps][relation]+=1
 majority={ps:v[True]>v[False] for ps,v in all_votes.items()}
 scored=sum(sum(v.values()) for v in gap_votes.values());correct=sum(v[majority[ps]] for ps,v in gap_votes.items())
 rows=[]
 with pysam.VariantFile(str(vcf_path)) as vcf:
  for r in vcf:
   if r.pos not in [left,right]:continue
   sample=next(iter(r.samples.values()));rows.append({'pos':r.pos,'alleles':r.alleles,'gt':sample['GT'],'phased':sample.phased,'ps':sample.get('PS')})
 result.append(dict(stage=stage,scored=scored,correct=correct,discordant=scored-correct,accuracy=correct/scored if scored else None,phase_sets=len(gap_votes),boundary_rows=rows))
p=Path('evaluations/2026-10-05-complementary-deletion-hiphase/read-truth.json');p.write_text(json.dumps(result,indent=2)+'\n')
for r in result: print(r['stage'],r['scored'],r['correct'],r['discordant'],r['accuracy'],r['phase_sets'])
