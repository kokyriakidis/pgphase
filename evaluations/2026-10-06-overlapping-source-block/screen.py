from pathlib import Path
from collections import Counter,defaultdict
import json,importlib.util,pysam
spec=importlib.util.spec_from_file_location('a','evaluations/2026-10-06-overlapping-source-block/read_audit.py');a=importlib.util.module_from_spec(spec);spec.loader.exec_module(a)
truth={f[0]:f[1]=='PATERNAL' for l in open('test_data/derived/chr20_truth_hap.tsv') if len(f:=l.rstrip().split('\t'))==2}
tools={'pgphase':a.assignments(Path('test_data/tmp_gap_fix86/final/0/phased.bam'),truth), 'hiphase': json.load(open('test_data/tmp_gap_fix78/hi-tags.json'))}
transitions=[]
last_pos,last_sets=None,set()
rows={}
with pysam.VariantFile('test_data/tmp_gap_fix86/final/0/phased.vcf') as vcf:
 for row in vcf:
  sample=next(iter(row.samples.values()));gt,ps=sample.get('GT'),sample.get('PS')
  if sample.phased and gt and None not in gt and len(set(gt))>1 and ps and ps>0:rows.setdefault(row.pos,set()).add(ps)
for pos,sets in sorted(rows.items()):
 if last_pos is not None and last_sets.isdisjoint(sets):transitions.append((last_pos,pos))
 last_pos,last_sets=pos,sets
original=pysam.AlignmentFile('test_data/HG002_chr20_hifi_mapped_to_CHM13_chr20_annotated.bam');results=[]
for target in json.load(open('evaluations/2026-10-06-overlapping-source-block/ranked-blocks.json')):
 if target['exact_span_covered']:continue
 bs=target['pgphase_blocks'];intervals=[];frontier=target['left']
 for b in bs:
  if b['left']>frontier:intervals.append((frontier,b['left']))
  frontier=max(frontier,b['right'])
 if frontier<target['right']:intervals.append((frontier,target['right']))
 intervals=sorted(set(intervals)|{(l,r) for l,r in transitions if target['left']<=l<r<=target['right']})
 for left,right in intervals:
  names={r.query_name for r in original.fetch('CHM13#0#chr20',left-1,right) if not r.is_secondary and not r.is_supplementary and r.query_name in truth}
  row={'rank':target['hiphase_rank'],'target':[target['left'],target['right']],'span':target['length'],'gap':[left,right],'scorable':len(names),'tools':{}}
  for tool,(tags,status) in tools.items():
   votes=defaultdict(Counter);cores=Counter();counts=Counter(status.get(q,'unphased') for q in names)
   for q in names:
    h,ps=tags.get(q,(0,0))
    if h in(1,2) and ps>0:votes[ps][(h==1)!=truth[q]]+=1
    if status.get(q)=='correct' and 0<ps<1000000000:cores[ps]+=1
   row['tools'][tool]={'counts':dict(counts),'core':max(cores.values(),default=0),'best_local_orientation_correct':sum(max(v.values()) for v in votes.values()),'orientation_votes':{str(ps):dict(v) for ps,v in votes.items()}}
  row['hiphase_global_80']=counts['correct']>=.8*len(names)
  row['hiphase_local_80']=row['tools']['hiphase']['best_local_orientation_correct']>=.8*len(names)
  results.append(row)
  print(json.dumps({k:v for k,v in row.items() if k!='tools'}),'HI',row['tools']['hiphase']['counts'],'PG',row['tools']['pgphase']['counts'],flush=True)
Path('evaluations/2026-10-06-overlapping-source-block/screened-seams.json').write_text(json.dumps(results,indent=2)+'\n')
