#!/usr/bin/env python3
"""Independently verify exact physical grouping, unique cohorts and all fits."""
from collections import defaultdict
import csv
import hashlib
import json
from pathlib import Path
import time

import numpy as np

ROOT=Path(__file__).resolve().parents[2]
OUT=Path(__file__).resolve().parent
WORK=ROOT/'test_data/tmp_representation_step9'
BASE=ROOT/'test_data/tmp_representation_step8/accepted'
previous=json.loads((ROOT/'evaluations/2026-10-07-representation-step8/manifest.json').read_text())
known={r['path']:r['sha256'] for r in previous['files']}

def rows(path):
    with path.open() as h:return list(csv.DictReader(h,delimiter='\t'))

def fit(matrix):
    if len(matrix)==0:return None,None,None,0
    pairs=[(a,b) for a in range(matrix.shape[1]) for b in range(a,matrix.shape[1])]
    totals=[int(np.minimum(matrix[:,a],matrix[:,b]).sum(dtype=np.int64)) for a,b in pairs]
    ranked=sorted(zip(totals,pairs))
    best=ranked[0][0];ties=sum(c==best for c in totals)
    return ranked[0][1] if ties==1 else None,best,ranked[1][0] if len(ranked)>1 else None,ties

def saved(r):
    return (None if r['allele0']=='.' else (int(r['allele0']),int(r['allele1'])),
            None if r['cost']=='.' else int(r['cost']),None if r['runner_up']=='.' else int(r['runner_up']),int(r['tied_pairs']))

reports=[];started=time.monotonic()
for owner in ['matrix4-final','matrix5-final','matrix65-final','whole-final']:
    base=BASE/owner;prod=WORK/'accepted'/owner
    for p in base.glob('matrix.chunk0.joint-*.tsv'):
        assert hashlib.sha256(p.read_bytes()).hexdigest()==known[str(p.relative_to(ROOT))],p
        assert p.read_bytes()==(prod/p.name).read_bytes(),(owner,p.name,'legacy projection changed')
    contexts={int(r['locus']):r for r in rows(base/'matrix.chunk0.joint-contexts.tsv')}
    expected_groups=defaultdict(list)
    for locus,c in contexts.items():
        if c['status']=='valid':
            alleles=tuple(sorted({c['allele0'],c['allele1']}|set(c['other_alleles'].split(','))-{''}))
            expected_groups[int(c['beg']),int(c['end']),alleles].append(locus)
    expected_groups=sorted(expected_groups.items())
    phys=rows(prod/'matrix.chunk0.physical-contexts.tsv')
    assert len(phys)==len(expected_groups)
    members=rows(prod/'matrix.chunk0.physical-members.tsv')
    hypotheses=defaultdict(list)
    for r in rows(prod/'matrix.chunk0.physical-alleles.tsv'):
        id_=int(r['physical']);assert int(r['allele'])==len(hypotheses[id_]);hypotheses[id_].append(r['sequence'])
    cohort=defaultdict(dict)
    for r in rows(prod/'matrix.chunk0.physical-cohort.tsv'):
        id_=int(r['physical']);name=r['read'];assert name not in cohort[id_]
        cohort[id_][name]={k:r[k] for k in ['query_beg','query_end','mapq','sequence','qualities_hex']}
    old_cohort=defaultdict(dict)
    for r in rows(base/'matrix.chunk0.joint-cohort.tsv'):
        old_cohort[int(r['locus'])][r['read']]={k:r[k] for k in ['query_beg','query_end','mapq','sequence','qualities_hex']}
    costs=defaultdict(dict)
    for r in rows(prod/'matrix.chunk0.physical-costs.tsv'):
        id_=int(r['physical']);assert r['read'] not in costs[id_]
        costs[id_][r['read']]=list(map(int,r['costs'].split(',')))
    old_costs=defaultdict(dict)
    for r in rows(base/'matrix.chunk0.joint-costs.tsv'):old_costs[int(r['locus'])][r['read']]=list(map(int,r['costs'].split(',')))
    full={int(r['physical']):r for r in rows(prod/'matrix.chunk0.physical-genotypes.tsv')}
    heldout={(int(r['physical']),r['read']):r for r in rows(prod/'matrix.chunk0.physical-heldout.tsv')}
    assert len(full)==len(expected_groups) and len(heldout)==sum(map(len,cohort.values()))
    membership=[];pair_checks=individual_cost_checks=0
    for id_,((beg,end,alleles),aliases) in enumerate(expected_groups):
        assert phys[id_]==dict(physical=str(id_),beg=str(beg),end=str(end))
        assert hypotheses[id_]==list(alleles)
        for alias in aliases:
            c=contexts[alias]
            membership.append(dict(physical=str(id_),locus=str(alias),allele0=str(alleles.index(c['allele0'])),allele1=str(alleles.index(c['allele1']))))
            assert cohort[id_]==old_cohort[alias] and costs[id_]==old_costs[alias]
        assert cohort[id_].keys()==costs[id_].keys()
        names=list(costs[id_]);matrix=np.asarray([costs[id_][n] for n in names],dtype=np.int64)
        expected=fit(matrix);assert saved(full[id_])==expected
        assert int(full[id_]['molecules'])==len(names) and full[id_]['conflicting_molecules']=='0'
        individual_cost_checks+=sum(map(len,costs[id_].values()));pair_checks+=1
        for i,name in enumerate(names):
            training=fit(np.delete(matrix,i,axis=0));actual=heldout[id_,name]
            assert saved(actual)==training
            winner=int(np.argmin(matrix[i]));unique=int(sum(matrix[i]==matrix[i,winner]))==1
            allele=winner if training[0] is not None and unique and winner in training[0] else -1
            same=training[0] is not None and training[0]==expected[0]
            assert int(actual['allele'])==allele and int(actual['same_pair'])==same
            pair_checks+=1
    assert members==membership
    for kind in ['contexts','members','alleles','costs','genotypes','heldout']:
        assert (prod/f'matrix.chunk0.physical-{kind}.tsv').read_bytes()==(WORK/'cached'/f'{owner}.physical-{kind}.tsv').read_bytes(),(owner,kind,'physical replay mismatch')
    for kind in ['alleles','costs','genotypes','heldout']:
        assert (prod/f'matrix.chunk0.joint-{kind}.tsv').read_bytes()==(WORK/'cached'/f'{owner}.{kind}.tsv').read_bytes(),(owner,kind,'projection replay mismatch')
    old_records=sum(map(len,old_cohort.values()));unique_records=sum(map(len,cohort.values()))
    report=dict(owner=owner,descriptions=sum(len(v) for _,v in expected_groups),physical_contexts=len(phys),
                alias_groups=[v for _,v in expected_groups if len(v)>1],legacy_molecule_records=old_records,
                unique_physical_molecule_records=unique_records,redundant_records_removed=old_records-unique_records,
                independently_verified_full_and_heldout_fits=pair_checks,unchanged_verified_individual_costs=individual_cost_checks)
    reports.append(report);print(json.dumps(report),flush=True)
for record in previous['input_metadata']:
    stat=(ROOT/record['path']).stat();assert stat.st_size==record['size'] and stat.st_mtime_ns==record['mtime_ns']
assert sum(r['redundant_records_removed'] for r in reports)==129
(OUT/'physical-checks.json').write_text(json.dumps(dict(seconds=time.monotonic()-started,owners=reports,
    physical_contexts=sum(r['physical_contexts'] for r in reports),physical_molecule_records=sum(r['unique_physical_molecule_records'] for r in reports),
    legacy_projections_identical=True,raw_cohorts_and_costs_match_independently_verified_bam_state=True),indent=2)+'\n')
