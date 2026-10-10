#!/usr/bin/env python3
"""Reuse frozen old costs, independently check novel costs and excluded discovery."""
import csv
import hashlib
import json
from collections import defaultdict
from pathlib import Path
import subprocess
ROOT=Path(__file__).resolve().parents[2]; OUT=Path(__file__).resolve().parent
WORK=ROOT/'test_data/tmp_representation_step14'
def rows(path):
    with path.open() as f:return list(csv.DictReader(f,delimiter='\t'))
def digest(path):return hashlib.sha256(path.read_bytes()).hexdigest()
def distance(a,b):
    previous=list(range(len(b)+1))
    for i,x in enumerate(a,1):
        current=[i]
        for j,y in enumerate(b,1):current.append(min(previous[j]+1,current[-1]+1,previous[j-1]+(x!=y)))
        previous=current
    return previous[-1]
known={r['path']:r['sha256'] for r in json.loads((ROOT/'evaluations/2026-10-07-representation-step11/manifest.json').read_text())['files']}
reports=[];target=None
for owner in ['matrix4-final','matrix5-final','matrix65-final','whole-final']:
    folder=WORK/'accepted-final'/owner
    oldfolder=ROOT/'test_data/tmp_representation_step11/accepted-final'/owner
    for kind in ['physical-cohort','path-alleles']:
        oldpath=oldfolder/f'matrix.chunk0.{kind}.tsv'
        assert digest(oldpath)==known[str(oldpath.relative_to(ROOT))]
        assert oldpath.read_bytes()==(folder/oldpath.name).read_bytes()
    assert (folder/'matrix.chunk0.path-alleles.tsv').read_bytes()==(folder/'matrix.chunk0.nested-alleles.tsv').read_bytes()
    baseline=ROOT/f'test_data/tmp_representation_step11/cached/{owner}.residuals.tsv'
    assert digest(baseline)==known[str(baseline.relative_to(ROOT))]
    old={(r['physical'],r['read']):int(r['path_min']) for r in rows(baseline)}
    cohort={(r['physical'],r['read']):r['sequence'] for r in rows(folder/'matrix.chunk0.physical-cohort.tsv')}
    assert old.keys()==cohort.keys()
    groups=defaultdict(lambda:defaultdict(set));novel=defaultdict(set);seeds=defaultdict(set)
    for r in rows(folder/'matrix.chunk0.nested-alleles.tsv'):seeds[r['physical']].add(r['sequence'])
    for (physical,name),sequence in cohort.items():groups[physical][sequence].add(name)
    for physical,seqs in groups.items():novel[physical]={s for s,names in seqs.items() if len(names)>=2 and s not in seeds[physical]}
    costs=WORK/'cached'/f'{owner}.read-costs.tsv'
    with costs.open('w') as f:subprocess.run([str(WORK/'score_read_catalog'),str(folder)],check=True,stdout=f)
    full=dict(old);heldout=dict(old);seen=set()
    for r in rows(costs):
        key=r['physical'],r['read'];sequence=r['sequence'];seen.add((*key,sequence))
        d=distance(cohort[key],sequence);assert d==int(r['distance'])
        support=groups[key[0]][sequence]-{key[1]}
        admitted=len(support)>=2;assert admitted==bool(int(r['heldout_admitted']))
        full[key]=min(full[key],d)
        if admitted:heldout[key]=min(heldout[key],d)
    assert seen=={(*key,s) for key in cohort for s in novel[key[0]]}
    residual_path=WORK/'cached'/f'{owner}.read-residuals.tsv'
    with residual_path.open('w') as f:
        f.write('physical\tread\told_min\tread_min\theldout_read_min\n')
        for key in old:f.write('\t'.join([*key,str(old[key]),str(full[key]),str(heldout[key])])+'\n')
    report=dict(owner=owner,rows=len(old),novel_distance_checks=len(seen),old_residual=sum(old.values()),read_residual=sum(full.values()),discovery_heldout_residual=sum(heldout.values()),improved_reads=sum(heldout[k]<v for k,v in old.items()),old_costs_reused_after_hash_and_input_identity=True,independent_scalar_dp_matches=True)
    reports.append(report)
    if owner=='matrix65-final':
        keys=[k for k in old if k[0]=='3']
        target=dict(physical=3,rows=len(keys),old_residual=sum(old[k] for k in keys),read_residual=sum(full[k] for k in keys),discovery_heldout_residual=sum(heldout[k] for k in keys),improved_reads=sum(heldout[k]<old[k] for k in keys),old_max_allele_length=max(map(len,seeds['3'])),new_max_allele_length=max(map(len,novel['3'])),exact_long_molecules=len(groups['3'][max(novel['3'],key=len)]))
(OUT/'residual-checks.json').write_text(json.dumps(dict(owners=reports,target=target),indent=2)+'\n')
print(json.dumps(dict(owners=reports,target=target),indent=2))
