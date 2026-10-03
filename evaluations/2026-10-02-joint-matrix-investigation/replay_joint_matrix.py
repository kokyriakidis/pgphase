from pathlib import Path
from collections import Counter,defaultdict
import argparse,json,subprocess,time
parser=argparse.ArgumentParser(description='Evaluation-only exact MEC on unchanged pgphase observations')
parser.add_argument('--chunk',type=int,required=True)
parser.add_argument('--free-flank-indels', action='store_true', help='Hold complete SNP flank gauges; rephase boundary indels as gap variables')
parser.add_argument('--hi-local-calls', action='store_true', help='Diagnostic override of shared calls; not a same-matrix comparison')
parser.add_argument('--left-ps',type=int,required=True)
parser.add_argument('--right-ps',type=int,required=True)
parser.add_argument('--gap-left',type=int,required=True)
parser.add_argument('--gap-right',type=int,required=True)
args=parser.parse_args()
root=Path('test_data/tmp_gap_next29')
path=Path('test_data/tmp_gap_next27')/f'joint-control-{args.chunk}'/'matrix.chunk0.recovery-final.tsv'
variants={};reads={};profiles=defaultdict(dict);seen=defaultdict(Counter)
for line in path.read_text().splitlines():
 f=line.split('\t')
 if f[0]=='VAR':
  variants[int(f[1])]={'pos':int(f[2]),'type':f[3],'flags':int(f[4]),'weight':int(f[5]),'ref_len':int(f[6]),'alt':f[7],'ps':int(f[8]),'h1':int(f[9]),'h2':int(f[10]),'bam':int(f[11]),'verified':int(f[12]),'af':float(f[16])}
 elif f[0]=='READ':reads[f[1]]={'mapq':int(f[4]),'skipped':int(f[5])!=0,'hap':int(f[6]),'ps':int(f[7])}
 elif f[0]=='OBS':
  a=int(f[3]);i=int(f[2]);profiles[f[1]][i]=a
  if not reads[f[1]]['skipped'] and a>=0:seen[i][a]+=1

if args.hi_local_calls:
 comparison=next(x for x in json.loads((root/'observation-comparison.json').read_text()) if x['chunk']==args.chunk)
 block=next(b for b in json.loads((root/f'hiphase{args.chunk}.matrix.json').read_text()) if b['beg']<=args.gap_left and b['end']>=args.gap_right)
 idx={}
 for site in comparison['sites']:
  i=site['index'];matches=[j for j,v in enumerate(block['variants']) if v['pos']==site['pos']-1]
  if len(matches)>1:matches=[j for j in matches if block['variants'][j]['alleles'][1][1:]==variants[i]['alt']]
  assert len(matches)==1
  idx[i]=matches[0]
 convert={'Reference':0,'Alternate':1,'Ambiguous':-1,'NoOverlap':-1}
 for read in block['reads']:
  if read['name'] not in reads:continue
  for i,j in idx.items():profiles[read['name']][i]=convert[read['alleles'][j]]
 seen=defaultdict(Counter)
 for n,p in profiles.items():
  if not reads[n]['skipped']:
   for i,a in p.items():
    if a>=0:seen[i][a]+=1

def oriented(v):return {v['h1'],v['h2']}=={0,1}
fixed={i for i,v in variants.items() if v['ps'] in [args.left_ps,args.right_ps] and oriented(v) and(not args.free_flank_indels or v['type']=='X')}
left={i for i in fixed if variants[i]['ps']==args.left_ps};right=fixed-left
beg=min(variants[i]['pos'] for i in left)
end=min(variants[i]['pos'] for i in right if variants[i]['type']=='X' )
# Genotyped homozygotes supply alignment context, but cannot orient haplotypes.
# This experiment treats each other binary site with both alleles as a variable.
extra={i for i,v in variants.items() if i not in fixed and beg<=v['pos']<=end and seen[i][0] and seen[i][1] and not any(a>1 for a in seen[i]) and not(v['ps']>0 and v['h1']==v['h2'] and v['h1']>=0)}
truth={f[0]:f[1]=='PATERNAL' for line in open('test_data/derived/chr20_truth_hap.tsv') if len(f:=line.rstrip().split('\t'))==2}
def fnv(s):
 h=14695981039346656037
 for b in s.encode():h=((h^b)*1099511628211)&((1<<64)-1)
 return h&1
out=[]
for mode in ['centered','all']:
 selected=fixed|{i for i in extra if mode=='all' or abs(variants[i]['af']-.5)<=.12}
 # Preserve the complete flanks; only candidate paths actually attached to
 # them need free variables. Disconnected components cannot choose parity.
 parent={i:i for i in selected}
 def find(i):
  while parent[i]!=i:parent[i]=parent[parent[i]];i=parent[i]
  return i
 for n,p in profiles.items():
  if reads[n]['skipped']:continue
  known=[i for i,a in p.items() if i in selected and a in [0,1]]
  for i in known[1:]:parent[find(i)]=find(known[0])
 lroots={find(i) for i in left};rroots={find(i) for i in right}
 bridge=lroots&rroots
 retained={i for i in selected if find(i) in bridge}
 var_ids=sorted(retained-fixed);var_of={i:j for j,i in enumerate(var_ids)}
 record={'mode':mode,'left_ps':args.left_ps,'right_ps':args.right_ps,'variable_scope':[beg,end],'extra_candidates':[dict(index=i,**variants[i],calls=dict(seen[i])) for i in sorted(extra)],'selected_extra':len(selected-fixed),'connected':bool(bridge),'variables':len(var_ids),'retained_fixed_sites':len(retained&fixed),'variable_sites':[dict(index=i,**variants[i]) for i in var_ids]}
 if len(var_ids)>20:
  record['resource_limit']=True;out.append(record);continue
 for weights in ['uniform','snp_priority']:
  calls={n:{i:a for i,a in p.items() if i in retained and a in [0,1]} for n,p in profiles.items() if not reads[n]['skipped']}
  indel_bound=sum(1 for p in calls.values() for i in p if variants[i]['type']!='X')
  weight=lambda i:indel_bound+1 if weights=='snp_priority' and variants[i]['type']=='X' else 1
  folds=[]
  for fold in [-1,0,1]:
   rows={n:p for n,p in calls.items() if len(p)>=2 and(fold<0 or fnv(n)==fold)}
   optimum=[];started=time.monotonic()
   for flip in [0,1]:
    encoded=[f'{len(rows)} {len(var_ids)}']
    for n,p in rows.items():
     mismatches=observations=0;costs={}
     for i,a in p.items():
      w=weight(i)
      if i in fixed:
       h=variants[i]['h1']^(flip if i in right else 0);mismatches+=w*(a!=h);observations+=w
      else:costs[var_of[i]]=[w*a,w*(1-a)]
     encoded.append(f'{mismatches} {observations} {len(costs)} '+' '.join(f'{vi} {c[0]} {c[1]}' for vi,c in costs.items()))
    t=subprocess.run([str(root/'mec_replay')],input='\n'.join(encoded)+'\n',text=True,capture_output=True,check=True)
    values=list(map(int,t.stdout.split()));optimum.append({'score':values[0],'bits':values[1:]})
   folds.append({'fold':fold,'rows':len(rows),'same':optimum[0],'flip':optimum[1],'elapsed_seconds':round(time.monotonic()-started,6),'winner':None if optimum[0]['score']==optimum[1]['score'] else int(optimum[1]['score']<optimum[0]['score'])})
  result=dict(record,weights=weights,folds=folds)
  # Truth is used below only to score the externally computed solution and
  # identify the correct connection between the pre-existing flank gauges.
  side_votes=[]
  for side in [left,right]:
   votes=Counter()
   for n,p in profiles.items():
    if n not in truth or reads[n]['skipped']:continue
    hp_votes=Counter()
    for i,a in p.items():
     if i in side and a in [0,1]:hp_votes[a==variants[i]['h1']]+=weight(i)
    if hp_votes[True]==hp_votes[False]:continue
    hp1=hp_votes[True]>hp_votes[False];votes[hp1!=truth[n]]+=1
   side_votes.append(dict(votes))
  result['flank_truth_votes']=side_votes
  result['expected_flip']=int(max(side_votes[0],key=side_votes[0].get)!=max(side_votes[1],key=side_votes[1].get))
  out.append(result)
p=root/f"matrix-replay-{args.chunk}{'-hi-local' if args.hi_local_calls else ''}{'-free-indels' if args.free_flank_indels else ''}.json";p.write_text(json.dumps(out,indent=2)+'\n')
for r in out:print(r['mode'],r.get('weights'),'connected',r['connected'],'variables',r['variables'],'scores',[(x['fold'],x['same']['score'],x['flip']['score'],x['winner']) for x in r.get('folds',[])],'expected',r.get('expected_flip'))
