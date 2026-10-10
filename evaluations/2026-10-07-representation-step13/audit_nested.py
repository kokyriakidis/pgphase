#!/usr/bin/env python3
"""Independent all-optimal mapping, anchored-edit composition and raw-state audit."""
from collections import Counter, defaultdict
import csv
import hashlib
import json
from pathlib import Path
import subprocess
import time
from composition_oracle import Oracle
from mapping_oracle import oracle

ROOT=Path(__file__).resolve().parents[2];OUT=Path(__file__).resolve().parent
WORK=ROOT/'test_data/tmp_representation_step13';BASE=ROOT/'test_data/tmp_representation_step12/accepted-final'
previous=json.loads((ROOT/'evaluations/2026-10-07-representation-step12/manifest.json').read_text())
known={r['path']:r['sha256'] for r in previous['files']}
def rows(path):
    with path.open() as stream:return list(csv.DictReader(stream,delimiter='\t'))
reports=[];start=time.monotonic()
for owner in ['matrix4-final','matrix5-final','matrix65-final','whole-final']:
    folder=WORK/'accepted-final'/owner;prefix=WORK/'cached'/owner
    for path in (BASE/owner).glob('matrix.chunk0.*.tsv'):
        assert hashlib.sha256(path.read_bytes()).hexdigest()==known[str(path.relative_to(ROOT))],path
        assert path.read_bytes()==(folder/path.name).read_bytes(),(owner,path.name)
    subprocess.run([str(ROOT/'test_allele_context'),'--nested',str(folder),str(prefix)],check=True)
    for kind in ['composition-sites','compositions','composed-alleles','path-status','sequence-paths','path-alleles','parent-maps','matched-subpaths','nested-compositions','nested-alleles']:
        assert (folder/f'matrix.chunk0.{kind}.tsv').read_bytes()==Path(str(prefix)+f'.{kind}.tsv').read_bytes(),(owner,kind)
    contexts={r['physical']:r for r in rows(folder/'matrix.chunk0.composition-contexts.tsv')}
    parents={(r['physical'],r['allele']):r['sequence'] for r in rows(folder/'matrix.chunk0.physical-alleles.tsv')}
    maps={};reference_oracles={}
    exported=defaultdict(list)
    for r in rows(folder/'matrix.chunk0.matched-subpaths.tsv'):
        exported[(r['physical'],r['parent_allele'])].append([int(r[k]) for k in ['ref_beg','alt_beg','length']])
    for key,allele in parents.items():
        reference=contexts[key[0]]['reference']
        cost,spans,count=oracle(reference,allele)
        assert spans==exported[key],(owner,key,'map')
        maps[key]=spans
    sites={(r['physical'],r['candidate'],r['alt_index']):r for r in rows(folder/'matrix.chunk0.composition-sites.tsv')}
    old=rows(folder/'matrix.chunk0.compositions.tsv');actual=rows(folder/'matrix.chunk0.nested-compositions.tsv')
    assert len(old)==len(actual)
    changes=Counter();resolved=[];complete=defaultdict(set);before=defaultdict(set)
    for r in rows(folder/'matrix.chunk0.path-alleles.tsv'):
        before[r['physical']].add(r['sequence']);complete[r['physical']].add(r['sequence'])
    for legacy,r in zip(old,actual):
        for k in ['physical','parent_allele','candidate','alt_index']:assert legacy[k]==r[k]
        key=(r['physical'],r['parent_allele']);ctx=contexts[key[0]];beg=int(ctx['beg']);reference=ctx['reference'];parent=parents[key]
        site=sites[(key[0],r['candidate'],r['alt_index'])];pos=int(site['pos']);ref=site['ref'];alt=site['alt']
        original=reference_oracles.setdefault(key[0],Oracle(beg,reference))
        status,sequence=original.compose([(beg,reference,parent),(pos,ref,alt)])
        assert (status,sequence or '.')==(legacy['status'],legacy['sequence']),(owner,key,'legacy oracle')
        if status=='overlap':
            independent,_=original.compose([(pos,ref,alt)])
            if independent!='valid':status,sequence=independent,None
            else:
                offset=pos-beg
                spans=[s for s in maps[key] if s[0]<=offset and offset+len(ref)<=s[0]+s[2] and
                       (len(ref)==len(alt) or s[0]<offset and offset+len(ref)<s[0]+s[2])]
                if not spans:status,sequence='overlap',None
                else:
                    i,j,length=spans[0]
                    status,sequence=Oracle(1,parent).compose([(j+offset-i+1,ref,alt)])
        assert (status,sequence or '.')==(r['status'],r['sequence']),(owner,key,r['candidate'],'nested oracle')
        if legacy['status']=='valid':assert r['status']=='valid' and r['sequence']==legacy['sequence']
        if (r['status'],r['sequence'])!=(legacy['status'],legacy['sequence']):
            changes[legacy['status']+'->'+r['status']]+=1;resolved.append(r)
        if status=='valid':complete[key[0]].add(sequence)
    hypotheses=defaultdict(list)
    for r in rows(folder/'matrix.chunk0.nested-alleles.tsv'):
        assert int(r['allele'])==len(hypotheses[r['physical']]);hypotheses[r['physical']].append(r['sequence'])
    assert all(hypotheses[k]==sorted(v) and before[k]<=v for k,v in complete.items())
    reports.append(dict(owner=owner,contexts=len(contexts),parents=len(parents),attempts=len(actual),changes=dict(changes),
        new_unique_sequences=sum(len(complete[k]-before[k]) for k in complete),final_hypotheses=sum(map(len,complete.values())),
        previous_tables_byte_identical=True,all_original_valid_compositions_retained=True,raw_replay_byte_identical=True,resolved=resolved))
result=dict(oracle='arbitrary-precision counts of all optimal alignments plus independent anchored-edit oracle',wall_seconds=time.monotonic()-start,owners=reports)
(OUT/'nested-checks.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result,indent=2))
