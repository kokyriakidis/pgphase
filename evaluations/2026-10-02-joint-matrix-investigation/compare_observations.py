"""Compare actual HiPhase local calls to unchanged graph/BAM working calls."""
from pathlib import Path
from collections import Counter,defaultdict
import json
root=Path('test_data/tmp_gap_next29');out=[]
for c,lo,hi,site_map,left_idx,right_idx in [(7,7264321,7280346,{114:7250571,116:7264321,123:7280346,124:7280356,125:7290365},[114,116],[123,124,125]),(24,24121713,24131707,{61:24103779,62:24105188,70:24121713,71:24121713,78:24131707,83:24142287,84:24142446,85:24143051},[61,62,70,71],[78,83,84,85])]:
 pg={};v={};prof=defaultdict(dict);channels=defaultdict(dict)
 for line in Path(f'test_data/tmp_gap_next27/joint-control-{c}/matrix.chunk0.recovery-final.tsv').read_text().splitlines():
  f=line.split('\t')
  if f[0]=='VAR':v[int(f[1])]=f
  elif f[0]=='READ':pg[f[1]]=f
  elif f[0]=='OBS':prof[f[1]][int(f[2])]=int(f[3]);channels[f[1]][int(f[2])]=list(map(int,f[3:6]))
 b=next(b for b in json.loads((root/f'hiphase{c}.matrix.json').read_text()) if b['beg']<=lo and b['end']>=hi)
 hipidx={}
 for i,pos in site_map.items():
  match=[j for j,h in enumerate(b['variants']) if h['pos']==pos-1]
  if len(match)>1:match=[j for j in match if b['variants'][j]['alleles'][1][1:]==v[i][7]]
  assert len(match)==1,(c,i,match)
  hipidx[i]=match[0]
 hip={};convert={'Reference':0,'Alternate':1,'Ambiguous':-1,'NoOverlap':-2}
 for r in b['reads']:
  calls={i:convert[r['alleles'][j]] for i,j in hipidx.items()}
  assert r['name'] not in hip,'No duplicated segments in these controls'
  hip[r['name']]=calls
 stat=[]
 for i,pos in site_map.items():
  counts=Counter();pairs=[]
  for n,p in hip.items():
   a=p[i];g=prof[n].get(i,-2)
   counts[(a,g)]+=1
   if a in [0,1] and g!=a:pairs.append({'read':n,'hi':a,'pg':g,'channels':channels[n].get(i)})
  stat.append({'index':i,'pos':pos,'comparison':{f'{a},{g}':n for(a,g),n in counts.items()},'different_calls':pairs})
 bridges=[]
 for n,h in hip.items():
  ha={i:h[i] for i in left_idx if h[i] in [0,1]};hb={i:h[i] for i in right_idx if h[i] in [0,1]}
  if ha and hb:
   p=prof[n];pa={i:p[i] for i in left_idx if p.get(i) in [0,1]};pb={i:p[i] for i in right_idx if p.get(i) in [0,1]}
   bridges.append({'read':n,'hip_left':ha,'hip_right':hb,'pg_left':pa,'pg_right':pb,'pg_channels':{i:channels[n].get(i) for i in site_map}})
 result={'chunk':c,'sites':stat,'hi_crossing_reads':len(bridges),'pg_crossing_among_hi':sum(bool(x['pg_left'] and x['pg_right']) for x in bridges),'bridges':bridges}
 out.append(result)
 print(c,'hipaired',len(bridges),'pgpaired',result['pg_crossing_among_hi'])
 for s in stat:
  print(s['index'],s['pos'],s['comparison'])
(root/'observation-comparison.json').write_text(json.dumps(out,indent=2)+'\n')
