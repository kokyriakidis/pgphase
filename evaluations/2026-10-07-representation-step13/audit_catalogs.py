#!/usr/bin/env python3
"""Prove identical catalogs/cohorts imply unchanged per-read costs without rescoring."""
import csv
import hashlib
import json
from pathlib import Path

ROOT=Path(__file__).resolve().parents[2];OUT=Path(__file__).resolve().parent
WORK=ROOT/'test_data/tmp_representation_step13/accepted-final'
previous=json.loads((ROOT/'evaluations/2026-10-07-representation-step11/manifest.json').read_text())
known={r['path']:r['sha256'] for r in previous['files']}
def digest(path):return hashlib.sha256(path.read_bytes()).hexdigest()
def rows(path):
    with path.open() as f:return list(csv.DictReader(f,delimiter='\t'))
reports=[]
for owner in ['matrix4-final','matrix5-final','matrix65-final','whole-final']:
    folder=WORK/owner
    old=folder/'matrix.chunk0.path-alleles.tsv';new=folder/'matrix.chunk0.nested-alleles.tsv'
    assert old.read_bytes()==new.read_bytes()
    baseline=ROOT/'test_data/tmp_representation_step11/accepted-final'/owner
    for kind in ['physical-alleles','physical-cohort','physical-costs','physical-genotypes','physical-heldout','path-alleles']:
        path=baseline/f'matrix.chunk0.{kind}.tsv'
        assert digest(path)==known[str(path.relative_to(ROOT))]
        assert path.read_bytes()==(folder/path.name).read_bytes()
    saved=ROOT/'test_data/tmp_representation_step11/cached'/f'{owner}.residuals.tsv'
    assert digest(saved)==known[str(saved.relative_to(ROOT))]
    cohort={(r['physical'],r['read']):r for r in rows(folder/'matrix.chunk0.physical-cohort.tsv')}
    costs={(r['physical'],r['read']):r['costs'] for r in rows(folder/'matrix.chunk0.physical-costs.tsv')}
    residuals={(r['physical'],r['read']):r for r in rows(saved)}
    assert cohort.keys()==costs.keys()==residuals.keys()
    assert all(costs[k]==r['old_costs'] for k,r in residuals.items())
    reports.append(dict(owner=owner,catalog_sha256=digest(new),catalogs_byte_identical=True,original_cohorts_costs_and_fits_byte_identical=True,
        baseline_residual_sha256=digest(saved),physical_read_context_rows=len(cohort),
        old_path_residual=sum(int(r['path_min']) for r in residuals.values()),nested_residual=sum(int(r['path_min']) for r in residuals.values()),
        costs_equal_by_identical_catalogs_and_queries=True,new_scoring_required=False))
(OUT/'catalog-checks.json').write_text(json.dumps(reports,indent=2)+'\n')
print(json.dumps(reports,indent=2))
