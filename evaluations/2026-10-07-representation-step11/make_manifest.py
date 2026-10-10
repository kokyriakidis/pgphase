#!/usr/bin/env python3
"""Freeze accepted step-11 evidence and require unchanged regression contracts."""
import hashlib
import json
from pathlib import Path

ROOT=Path(__file__).resolve().parents[2]
OUT=Path(__file__).resolve().parent
WORK=ROOT/'test_data/tmp_representation_step11'

def digest(path):
    h=hashlib.sha256()
    with path.open('rb') as handle:
        while data:=handle.read(1024*1024):h.update(data)
    return h.hexdigest()

previous=json.loads((ROOT/'evaluations/2026-10-07-representation-step10/manifest.json').read_text())
for record in previous['input_metadata']:
    stat=(ROOT/record['path']).stat()
    assert stat.st_size==record['size'] and stat.st_mtime_ns==record['mtime_ns'],record['path']
unchanged=[]
for line in (WORK/'baseline/floors.sha256').read_text().splitlines():
    sha,name=line.split(maxsplit=1)
    if name=='pgphase':continue
    assert digest(ROOT/name)==sha,name
    unchanged.append(name)
sha=digest(ROOT/'pgphase')
assert sha==(WORK/'final-binary.sha256').read_text().split()[0]
checks=json.loads((WORK/'accepted-final/windows/checks.json').read_text())
assert checks['registered_checks']==checks['passed']==107 and checks['binary_sha256']==sha
full=json.loads((OUT/'full-checks.json').read_text())
assert full['candidates_byte_identical'] and full['vcf_rows_identical']
assert full['changed_primary_hp_ps']==full['changed_parental_status']==0
assert full['before_blocks']==full['after_blocks']
assert full['connected_core_tags_unchanged'] and full['previously_phased_tags_unchanged']
assert full['before_parental_counts']==full['after_parental_counts']
panel=json.loads((OUT/'panel-comparison.json').read_text())
assert panel['panel_windows']==141 and panel['changed_windows']==0 and panel['summary_identical']
for log in ['build-final-fast.log']:
    assert 'warning:' not in (OUT/log).read_text() and 'error:' not in (OUT/log).read_text(),log
for path in ['path-checks.json','residual-checks.json','owner-checks.json','mutation-checks.json','fast-checks.json','raw-state-checks.json']:
    assert (OUT/path).is_file()
paths=json.loads((OUT/'path-checks.json').read_text())['owners']
assert sum(r['contexts'] for r in paths)==70
assert all(r['statuses']=={'complete':r['contexts']} for r in paths)
assert sum(r['paths'] for r in paths)==24260
assert sum(r['new_unique_sequences'] for r in paths)==11018
assert sum(r['unpruned_full_assignments_verified'] for r in paths)==62577
assert all(r['previous_tables_byte_identical'] for r in paths)
assert len(json.loads((OUT/'mutation-checks.json').read_text()))==12
assert all(r['exit']!=0 for r in json.loads((OUT/'mutation-checks.json').read_text()).values())
assert 'All tests passed' in (OUT/'standard-window-tests-final.log').read_text()
files=[]
selected=[ROOT/'pgphase',ROOT/'test_allele_context',ROOT/'test_allele_genotype',ROOT/'docs/IMPLEMENTATION.md',ROOT/'CHECKPOINT.md']
selected+=list((ROOT/'src').glob('allele_*.?pp'))+list((ROOT/'src').glob('test_allele_*.cpp'))+[ROOT/'src/collect_pipeline.cpp']
selected += [ROOT/p for p in unchanged]
for directory in [OUT,WORK/'accepted-final',WORK/'cached',WORK/'mutants',WORK/'timing']:
    selected += [p for p in directory.rglob('*') if p.is_file() and '__pycache__' not in p.parts and p.name not in {'manifest.json','standard-window-tests.log'}]
selected += [WORK/'final-binary.sha256',WORK/'baseline/floors.sha256',WORK/'score_catalog',WORK/'candidate-pgphase',WORK/'candidate-binary.sha256']
selected += [p for p in (WORK/'raw-state-only').rglob('*') if p.is_file()]
for p in sorted(set(selected)):
    files.append(dict(path=str(p.relative_to(ROOT)),size=p.stat().st_size,sha256=digest(p)))
manifest=dict(production_sha256=sha,starting_sha256=digest(WORK/'baseline/pgphase'),
              accepted_output='test_data/tmp_representation_step11/accepted-final/full-frozen',
              unchanged_floor_and_panel_files=unchanged,input_metadata=previous['input_metadata'],
              competitor_panel_signature=json.loads((OUT/'hiphase-verified.tsv.json').read_text())['signature'],
              baseline_independent_cohort_and_fit_manifest_sha256=digest(ROOT/'evaluations/2026-10-07-representation-step10/manifest.json'),
              files=files)
(OUT/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
print(f'Frozen {len(files)} files; all 107 checks, 141 panel contracts and phase outputs unchanged')
