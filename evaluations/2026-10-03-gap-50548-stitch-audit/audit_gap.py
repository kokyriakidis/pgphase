#!/usr/bin/env python3
"""Evaluate one gap's actual matrices and parental read assignments."""
from collections import Counter, defaultdict
from pathlib import Path
import hashlib, json
import pysam
root=Path('test_data/tmp_gap_inspect40/probe/50')
report=Path('evaluations/2026-10-03-gap-50548-stitch-audit')
left,right=50548245,50562066
bam_path='test_data/HG002_chr20_hifi_mapped_to_CHM13_chr20_annotated.bam'
truth={f[0]:f[1]=='PATERNAL' for l in open('test_data/derived/chr20_truth_hap.tsv') if len(f:=l.rstrip().split('\t'))==2}
def matrix(path):
 sites,reads,obs,qualities={}, {}, defaultdict(dict),defaultdict(dict)
 for l in path.open():
  f=l.rstrip('\n').split('\t')
  if f[0]=='VAR':sites[int(f[1])]={'pos':int(f[2]),'type':f[3],'mask':int(f[4]),'ref_len':int(f[6]),'alt':f[7],'ps':int(f[8]),'gt':list(map(int,f[9:11])),'injected':bool(int(f[11]))}
  elif f[0]=='READ':reads[f[1]]={'beg':int(f[2]),'end':int(f[3]),'mapq':int(f[4]),'skip':bool(int(f[5])),'hp':int(f[6]),'ps':int(f[7])}
  elif f[0]=='OBS':obs[f[1]][int(f[2])]=list(map(int,f[3:6]))
  elif f[0]=='BAMQ':qualities[f[1]][int(f[2])]=int(f[3])
 return sites,reads,obs,qualities
live,live_reads,live_obs,live_q=matrix(root/'matrix.chunk0.recovery-final.tsv')
source,source_reads,source_obs,source_q=matrix(root/'matrix.recovery.chunk0.window2.chunk-1.recovery-source.tsv')
saved_sites,saved_obs={},defaultdict(dict)
for l in (root/'stderr.log').open():
 f=l.rstrip().split('\t')
 if f[0]=='TRACE_SOURCE_SITE' and f[1]=='1':saved_sites[int(f[2])]={'pos':int(f[3]),'type':int(f[4]),'ps':int(f[7]),'gt':list(map(int,f[8:10])),'clean_snp':bool(int(f[10]))}
 elif f[0]=='TRACE_SOURCE_OBS' and f[1]=='1':saved_obs[f[2]][int(f[4])]={'allele':int(f[5]),'quality':int(f[6])}
local,spans={},{}
with pysam.AlignmentFile(bam_path) as bam:
 for r in bam.fetch('CHM13#0#chr20',left-1,right):
  if r.is_secondary or r.is_supplementary or r.is_unmapped:continue
  local[r.query_name]=r
  if r.reference_start<=left-1 and r.reference_end>=right:spans[r.query_name]=r
high={n:r for n,r in spans.items() if 30<=r.mapping_quality<255}
source_blocks=defaultdict(list)
for s in source.values():
 if s['ps']>0 and min(s['gt'])>=0 and len(set(s['gt']))==2:source_blocks[s['ps']].append(s['pos'])
saved_blocks=defaultdict(list)
for s in saved_sites.values():saved_blocks[s['ps']].append(s['pos'])
boundary_live=[i for i,s in live.items() if s['pos']==left+1 and s['type']=='D' and s['ref_len'] in [1,2] and s['alt']=='']
boundary_source=[i for i,s in source.items() if s['pos']==left+1 and s['type']=='D' and s['ref_len'] in [1,2] and s['alt']=='']
read_details=[]
for n,r in high.items():
 calls=[{'pos':saved_sites[i]['pos'],'gt':saved_sites[i]['gt'],**o} for i,o in saved_obs.get(n,{}).items() if o['allele'] in saved_sites[i]['gt']]
 read_details.append({'qname':n,'beg':r.reference_start+1,'end':r.reference_end,'mapq':r.mapping_quality,'source_tag':source_reads.get(n),'saved_callable_sites':calls,'live_deletion_calls':{live[i]['ref_len']:live_obs.get(n,{}).get(i,[-99]*3)[0] for i in boundary_live},'source_deletion_calls':{source[i]['ref_len']:source_obs.get(n,{}).get(i,[-99]*3)[0] for i in boundary_source}})
def tags_and_score(path):
 tags={};votes=defaultdict(Counter)
 with pysam.AlignmentFile(str(path),check_sq=False) as bam:
  for r in bam:
   if r.is_secondary or r.is_supplementary or not r.has_tag('HP') or not r.has_tag('PS'):continue
   hp,ps=r.get_tag('HP'),r.get_tag('PS')
   if hp not in [1,2] or ps<=0:continue
   tags[r.query_name]=(hp,ps)
   if r.query_name in truth:votes[ps][(hp==1)!=truth[r.query_name]]+=1
 majority={ps:v[True]>v[False] for ps,v in votes.items()}
 names=[n for n in local if n in tags and n in truth]
 correct=[n for n in names if ((tags[n][0]==1)!=truth[n])==majority[tags[n][1]]]
 return {'phased_truth_reads_overlapping_gap':len(names),'correct_using_whole_block_orientation':len(correct),'discordant':len(names)-len(correct),'spanning_reads_phased':sum(n in tags for n in high),'spanning_read_tags':{n:tags[n] for n in high if n in tags}}
competitor_root=Path('/tmp/pgphase-hiphase-comparison-2026-10-01')
result={'gap':[left,right],'production_sha256':hashlib.sha256(Path('pgphase').read_bytes()).hexdigest(),'input_primary_overlaps':len(local),'physical_primary_spanners':len(spans),'spanner_mapq_distribution':dict(Counter(r.mapping_quality for r in spans.values())),'mapq30_spanners':len(high),'mapq30_spanners_present_in_live_matrix':sum(n in live_reads for n in high),'mapq30_spanners_present_in_source_matrix':sum(n in source_reads for n in high),'mapq30_spanners_present_in_saved_matrix':sum(n in saved_obs for n in high),'all_raw_source_blocks':{ps:{'n_sites':len(p),'first':min(p),'last':max(p)} for ps,p in source_blocks.items()},'retained_source_blocks':{ps:{'n_sites':len(p),'first':min(p),'last':max(p)} for ps,p in saved_blocks.items()},'mapq30_spanners_with_callable_left_source':sum(bool(d['saved_callable_sites']) for d in read_details),'read_details':read_details,'outcomes':{'pgphase':tags_and_score(Path('test_data/tmp_gap_next39/full-retained/phased.bam')),'hiphase_dv':tags_and_score(competitor_root/'hiphase_dv/phased.bam'),'hiphase_pg_calls':tags_and_score(competitor_root/'hiphase_pg_calls/phased.bam')}}
(report/'audit.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps({k:v for k,v in result.items() if k not in ['read_details','outcomes']},indent=2))
print(json.dumps({a:{k:v for k,v in o.items() if k!='spanning_read_tags'} for a,o in result['outcomes'].items()},indent=2))
