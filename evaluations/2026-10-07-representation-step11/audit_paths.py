#!/usr/bin/env python3
"""Independent anchored-edit pruning, complete DFS census and path rendering."""
from collections import Counter,defaultdict
import csv
import hashlib
import itertools
import json
from pathlib import Path
import time

ROOT=Path(__file__).resolve().parents[2];OUT=Path(__file__).resolve().parent
WORK=ROOT/'test_data/tmp_representation_step11';BASE=ROOT/'test_data/tmp_representation_step10/accepted'
previous=json.loads((ROOT/'evaluations/2026-10-07-representation-step10/manifest.json').read_text())
known={r['path']:r['sha256'] for r in previous['files']}

def rows(path):
 with path.open() as handle:return list(csv.DictReader(handle,delimiter='\t'))

def placements(pos,ref,alt):
 if ref==alt:return [None]
 left=0
 while left<min(len(ref),len(alt)) and ref[left]==alt[left]:left+=1
 right=0
 while right<min(len(ref),len(alt)) and ref[-1-right]==alt[-1-right]:right+=1
 trimmed=min(len(ref),len(alt),left+right)
 return [(pos+p,len(ref)-trimmed,alt[p:len(alt)-(trimmed-p) if trimmed-p else len(alt)])
         for p in range(max(0,trimmed-right),min(left,trimmed)+1)]

def render(reference,beg,events):
 cursor=beg;parts=[];insertions=set()
 for pos,length,alt in sorted(set(e for e in events if e is not None)):
  if pos<cursor or not length and pos in insertions:return None
  parts.extend([reference[cursor-beg:pos-beg],alt]);cursor=pos+length
  if not length:insertions.add(pos)
 return ''.join(parts)+reference[cursor-beg:]

class Oracle:
 def __init__(self,beg,reference):self.beg=beg;self.reference=reference;self.cache={}
 def prepare(self,raw):
  if raw in self.cache:return self.cache[raw]
  pos,ref,alt=raw;ref=ref.upper();alt=alt.upper()
  if pos<self.beg or pos+len(ref)>self.beg+len(self.reference) or len(ref)>4096 or len(alt)>4096 or set(ref+alt)-set('ACGT') or self.reference[pos-self.beg:pos-self.beg+len(ref)]!=ref:
   result=('unsupported',None,None,None,None)
  else:
   choices=placements(pos,ref,alt);key=choices[0]
   if key is None:result=('noop',None,None,None,choices)
   else:
    effect=render(self.reference,self.beg,[key]);region=None
    if choices[0][0]!=choices[-1][0]:region=(choices[0][0],choices[-1][0]+key[1])
    result=('edit',key,region,effect,choices)
  self.cache[raw]=result;return result
 def compose(self,raw_edits):
  canonical=[(pos,ref.upper(),alt.upper()) for pos,ref,alt in raw_edits]
  data=[self.prepare(raw) for raw in dict.fromkeys(canonical)]
  if any(r[0]=='unsupported' for r in data):return 'unsupported',None
  data=[r for r in data if r[0]=='edit'];counts=Counter(r[1] for r in data)
  for _,key,region,_,_ in data:
   if region:
    if counts[key]>1:return 'overlap',None
    if any(k!=key and k[0]<=region[1] and k[0]+k[1]>=region[0] for k in counts):return 'overlap',None
  effects={}
  for _,key,_,effect,_ in data:
   if effect in effects and effects[effect]!=key:return 'overlap',None
   effects[effect]=key
  sequence=render(self.reference,self.beg,list(counts))
  if sequence is None:return 'overlap',None
  if len(sequence)>4096:return 'unsupported',None
  return 'valid',sequence

def physical_groups(sites):
 grouped={}
 for candidate,pos,ref,alts in sites:
  grouped.setdefault((pos,ref,tuple(sorted(set(alts)))),[]).append((candidate,alts))
 return list(grouped.items())

def group_selection(group,choice):
 (pos,ref,alts),members=group
 if not choice:return [],[]
 alt=alts[choice-1]
 return [(pos,ref,alt)],[(candidate,original.index(alt)+1) for candidate,original in members]

def enumerate_paths(reference,beg,parents,sites,budget=65536):
 oracle=Oracle(beg,reference);paths=[];visits=overlaps=unsupported=0;limited=False
 groups=physical_groups(sites)
 def search(parent,next_,edits,selected):
  nonlocal visits,overlaps,unsupported,limited
  if limited:return
  if visits==budget:limited=True;paths.clear();return
  visits+=1;status,sequence=oracle.compose(edits)
  if status=='overlap':overlaps+=1;return
  if next_==len(groups):
   if status=='valid':paths.append((parent,tuple(sorted(selected)),sequence))
   else:unsupported+=1
   return
  group=groups[next_]
  search(parent,next_+1,edits,selected)
  for choice in range(1,len(group[0][2])+1):
   edit,aliases=group_selection(group,choice)
   search(parent,next_+1,edits+edit,selected+aliases)
   if limited:break
 for parent,sequence in enumerate(parents):
  search(parent,0,[(beg,reference,sequence)],[])
  if limited:break
 return dict(status='limited' if limited else 'complete',visited_prefixes=visits,overlap_prefixes=overlaps,unsupported_paths=unsupported,paths=len(paths)),paths,oracle

start=time.monotonic();reports=[]
for owner in ['matrix4-final','matrix5-final','matrix65-final','whole-final']:
 after=WORK/'accepted-final'/owner;before=BASE/owner
 for p in before.glob('matrix.chunk0.*.tsv'):
  assert hashlib.sha256(p.read_bytes()).hexdigest()==known[str(p.relative_to(ROOT))],p
  assert p.read_bytes()==(after/p.name).read_bytes(),(owner,p.name,'existing state changed')
 contexts={int(r['physical']):r for r in rows(after/'matrix.chunk0.composition-contexts.tsv')}
 parents=defaultdict(list)
 for r in rows(after/'matrix.chunk0.physical-alleles.tsv'):parents[int(r['physical'])].append(r['sequence'])
 old=defaultdict(set)
 for r in rows(after/'matrix.chunk0.composed-alleles.tsv'):old[int(r['physical'])].add(r['sequence'])
 scoped=defaultdict(dict)
 for r in rows(after/'matrix.chunk0.composition-sites.tsv'):
  id_,candidate=int(r['physical']),int(r['candidate']);sites=scoped[id_]
  if candidate not in sites:sites[candidate]=[candidate,int(r['pos']),r['ref'],[]]
  assert int(r['alt_index'])==len(sites[candidate][3])+1
  sites[candidate][3].append(r['alt'])
 emitted=defaultdict(list)
 for r in rows(after/'matrix.chunk0.sequence-paths.tsv'):
  id_=int(r['physical']);assert int(r['path'])==len(emitted[id_])
  selected=tuple(tuple(map(int,choice.split(':'))) for choice in r['neighbors'].split(',')) if r['neighbors']!='.' else ()
  assert len({c for c,a in selected})==len(selected)
  emitted[id_].append((int(r['parent_allele']),selected,r['sequence']))
 status_rows=rows(after/'matrix.chunk0.path-status.tsv');assert len(status_rows)==len(contexts)
 statuses={int(r['physical']):r for r in status_rows}
 hypotheses=defaultdict(list)
 for r in rows(after/'matrix.chunk0.path-alleles.tsv'):
  id_=int(r['physical']);assert int(r['allele'])==len(hypotheses[id_]);hypotheses[id_].append(r['sequence'])
 counts=Counter();new=multi=source_multi=checked=exhaustive_leaves=0;details=[]
 for id_,ctx in contexts.items():
  beg=int(ctx['beg']);ref=ctx['reference'];sites=[v for k,v in sorted(scoped[id_].items())]
  expected,paths,oracle=enumerate_paths(ref,beg,parents[id_],sites)
  actual={k:v if k=='status' else int(v) for k,v in statuses[id_].items() if k!='physical'}
  assert expected==actual,(owner,id_,expected,actual)
  assert paths==emitted[id_],(owner,id_,'path census/provenance mismatch')
  assert sorted(old[id_]|{s for p,n,s in paths})==hypotheses[id_]
  if actual['status']=='limited':assert not paths and sorted(old[id_])==hypotheses[id_]
  by_candidate={s[0]:s for s in sites}
  if actual['status']=='complete':
   # No pruning: enumerate every full raw assignment to prove no valid path
   # disappeared through an intermediate overlap or unsupported prefix.
   exhaustive=[]
   for parent,sequence in enumerate(parents[id_]):
    for assignment in itertools.product(*(range(len(group[0][2])+1) for group in physical_groups(sites))):
     exhaustive_leaves+=1
     selections=[group_selection(group,choice) for group,choice in zip(physical_groups(sites),assignment)]
     selected=tuple(sorted(alias for edit,aliases in selections for alias in aliases))
     raw=[(beg,ref,sequence)]+[edit for edits,aliases in selections for edit in edits]
     state,rendered=oracle.compose(raw)
     if state=='valid':exhaustive.append((parent,selected,rendered))
   assert exhaustive==paths,(owner,id_,'pruning lost a valid full assignment')
  for parent,selected,sequence in paths:
   raw=[(beg,ref,parents[id_][parent])]+[(by_candidate[c][1],by_candidate[c][2],by_candidate[c][3][a-1]) for c,a in selected]
   choices=[]
   for edit in dict.fromkeys(raw):
    positions=oracle.prepare(edit)[4];choices.append(list(dict.fromkeys([positions[0],positions[-1]])))
   assert {render(ref,beg,events) for events in itertools.product(*choices)}=={sequence},(owner,id_,selected,'padding placements differ')
   checked+=1;source_multi+=len(selected)>=2
   selected_groups={(by_candidate[c][1],by_candidate[c][2],tuple(sorted(set(by_candidate[c][3])))) for c,a in selected}
   multi+=len(selected_groups)>=2
  counts[actual['status']]+=1;added=len(set(hypotheses[id_])-old[id_]);new+=added
  details.append(dict(physical=id_,neighbor_sites=len(sites),physical_neighbor_sites=len(physical_groups(sites)),**actual,added_sequences=added))
 for kind in ['path-status','sequence-paths','path-alleles','composition-sites','compositions','composed-alleles']:
  assert (after/f'matrix.chunk0.{kind}.tsv').read_bytes()==(WORK/'cached'/f'{owner}.{kind}.tsv').read_bytes(),(owner,kind,'cached replay mismatch')
 report=dict(owner=owner,contexts=len(contexts),statuses=dict(counts),paths=checked,unpruned_full_assignments_verified=exhaustive_leaves,multiple_physical_neighbor_paths=multi,multiple_source_description_paths=source_multi,new_unique_sequences=new,old_hypotheses=sum(map(len,old.values())),path_hypotheses=sum(map(len,hypotheses.values())),previous_tables_byte_identical=True,physical_contexts=details)
 reports.append(report);print(json.dumps({k:v for k,v in report.items() if k!='physical_contexts'}),flush=True)
for r in previous['input_metadata']:
 stat=(ROOT/r['path']).stat();assert stat.st_size==r['size'] and stat.st_mtime_ns==r['mtime_ns']
(OUT/'path-checks.json').write_text(json.dumps(dict(seconds=time.monotonic()-start,owners=reports),indent=2)+'\n')
