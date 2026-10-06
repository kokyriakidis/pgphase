from pathlib import Path
import json,pysam,importlib.util,collections
spec=importlib.util.spec_from_file_location('audit','evaluations/2026-10-05-largest-hiphase-block/audit_reads.py');a=importlib.util.module_from_spec(spec);spec.loader.exec_module(a)
truth={f[0]:f[1]=='PATERNAL' for l in open('test_data/derived/chr20_truth_hap.tsv') if len(f:=l.rstrip().split('\t'))==2}
p=Path('test_data/tmp_gap_fix72');before=p/'baseline/24';after=p/'final_checked/24'
r=a.audit(before,after,truth)
with pysam.AlignmentFile('test_data/HG002_chr20_hifi_mapped_to_CHM13_chr20_annotated.bam') as bam:
 names={x.query_name for x in bam.fetch('CHM13#0#chr20',24121712,24131707) if not x.is_secondary and not x.is_supplementary and x.query_name in truth}
for tool,path in [('before',before/'phased.bam'),('after',after/'phased.bam')]:
 tags,states=a.assignments(path,truth);counts=collections.Counter(states.get(q,'unphased') for q in names);cores=collections.Counter(tags[q][1] for q in names if states.get(q)=='correct' and tags[q][1]<1000000000)
 r[tool]={'scorable':len(names),'counts':dict(counts),'core':cores.most_common(),'correct_all':counts['correct']/len(names)}
benchmark=json.loads(Path('evaluations/2026-10-05-third-largest-block-target/next-block.json').read_text())['tools']['hiphase']['regions']['gap_1']
r['hiphase']=benchmark
r['benchmark_source']='evaluations/2026-10-05-third-largest-block-target/next-block.json (identical original alignments verified)'
assert r['after']['counts']['correct'] >= benchmark['counts']['correct']
assert r['after']['core'][0][1] >= benchmark['dominant_connected_correct_core']
assert r['after']['correct_all'] >= 0.8
Path('evaluations/2026-10-05-complementary-repeat-insertion/owner-results.json').write_text(json.dumps(r,indent=2)+'\n')
print(json.dumps(r,indent=2))
