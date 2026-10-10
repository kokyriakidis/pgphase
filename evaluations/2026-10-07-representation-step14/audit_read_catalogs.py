#!/usr/bin/env python3
"""Independently group physical read sequences and verify raw catalog replay."""
from collections import defaultdict
import csv
import hashlib
import json
from pathlib import Path
import subprocess

ROOT=Path(__file__).resolve().parents[2];OUT=Path(__file__).resolve().parent
WORK=ROOT/'test_data/tmp_representation_step14';BASE=ROOT/'test_data/tmp_representation_step13/accepted-final'
previous=json.loads((ROOT/'evaluations/2026-10-07-representation-step13/manifest.json').read_text())
known={r['path']:r['sha256'] for r in previous['files']}
def rows(path):
    with path.open() as f:return list(csv.DictReader(f,delimiter='\t'))
def oracle(observations,excluded=None):
    by_name=defaultdict(set);unsupported=set()
    for name,sequence in observations:
        if name==excluded:continue
        sequence=sequence.upper()
        if not name or not sequence or len(sequence)>4096 or set(sequence)-set('ACGT'):unsupported.add(name)
        else:by_name[name].add(sequence)
    groups=defaultdict(set);conflicts=set()
    for name,sequences in by_name.items():
        if name in unsupported:continue
        if len(sequences)!=1:conflicts.add(name)
        else:groups[next(iter(sequences))].add(name)
    return {seq:sorted(names) for seq,names in sorted(groups.items()) if len(names)>=2},conflicts,unsupported
reports=[]
for owner in ['matrix4-final','matrix5-final','matrix65-final','whole-final']:
    folder=WORK/'accepted-final'/owner
    for path in (BASE/owner).glob('matrix.chunk0.*.tsv'):
        assert hashlib.sha256(path.read_bytes()).hexdigest()==known[str(path.relative_to(ROOT))],path
        assert path.read_bytes()==(folder/path.name).read_bytes(),(owner,path.name)
    subprocess.run([str(ROOT/'test_allele_context'),'--read-catalogs',str(folder),str(WORK/'cached'/owner)],check=True)
    for kind in ['composition-sites','compositions','composed-alleles','path-status','sequence-paths','path-alleles','parent-maps','matched-subpaths','nested-compositions','nested-alleles','read-catalog-status','read-hypotheses','read-support','read-exclusions','read-alleles']:
        assert (folder/f'matrix.chunk0.{kind}.tsv').read_bytes()==Path(str(WORK/'cached'/owner)+f'.{kind}.tsv').read_bytes(),(owner,kind)
    cohort=defaultdict(list)
    for r in rows(folder/'matrix.chunk0.physical-cohort.tsv'):cohort[r['physical']].append((r['read'],r['sequence']))
    actual=defaultdict(dict);support=defaultdict(list)
    for r in rows(folder/'matrix.chunk0.read-support.tsv'):support[r['physical'],r['hypothesis']].append(r['molecule'])
    for r in rows(folder/'matrix.chunk0.read-hypotheses.tsv'):
        actual[r['physical']][r['sequence']]=support[r['physical'],r['hypothesis']]
        assert len(support[r['physical'],r['hypothesis']])==int(r['molecules'])
    old=defaultdict(set);expanded=defaultdict(list)
    for r in rows(folder/'matrix.chunk0.nested-alleles.tsv'):old[r['physical']].add(r['sequence'])
    for r in rows(folder/'matrix.chunk0.read-alleles.tsv'):
        assert int(r['allele'])==len(expanded[r['physical']]);expanded[r['physical']].append(r['sequence'])
    statuses=rows(folder/'matrix.chunk0.read-catalog-status.tsv')
    assert len(statuses)==len(old)
    heldout=0
    for r in statuses:
        key=r['physical'];expected,conflicts,unsupported=oracle(cohort[key])
        assert expected==actual[key] and not conflicts and not unsupported
        assert r['status']=='complete' and int(r['molecules'])==len(cohort[key]) and int(r['hypotheses'])==len(expected)
        assert r['conflicting']==r['unsupported']=='0'
        assert expanded[key]==sorted(old[key]|set(expected))
        for name,_ in cohort[key]:
            withheld,_,_=oracle(cohort[key],name)
            assert all(name not in members for members in withheld.values())
            assert withheld=={s:[n for n in names if n!=name] for s,names in expected.items() if len(names)-(name in names)>=2}
            heldout+=1
    assert not rows(folder/'matrix.chunk0.read-exclusions.tsv')
    reports.append(dict(owner=owner,contexts=len(old),read_context_rows=sum(map(len,cohort.values())),
        hypotheses=sum(map(len,actual.values())),new_unique_sequences=sum(len(set(actual[k])-old[k]) for k in old),
        full_hypotheses=sum(map(len,expanded.values())),statuses_complete=True,previous_tables_byte_identical=True,
        independent_support_and_catalogs_verified=True,heldout_discovery_censuses=heldout,raw_replay_byte_identical=True))
(OUT/'read-catalog-checks.json').write_text(json.dumps(reports,indent=2)+'\n')
print(json.dumps(reports,indent=2))
