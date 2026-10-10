#!/usr/bin/env python3
"""Freeze accepted step-9 evidence and require unchanged regression contracts."""
import hashlib
import json
from pathlib import Path

ROOT=Path(__file__).resolve().parents[2]
OUT=Path(__file__).resolve().parent
WORK=ROOT/'test_data/tmp_representation_step9'

def digest(path):
    h=hashlib.sha256()
    with path.open('rb') as handle:
        while data:=handle.read(1024*1024):h.update(data)
    return h.hexdigest()

previous=json.loads((ROOT/'evaluations/2026-10-07-representation-step8/manifest.json').read_text())
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
checks=json.loads((WORK/'accepted/windows/checks.json').read_text())
assert checks['registered_checks']==checks['passed']==107 and checks['binary_sha256']==sha
full=json.loads((OUT/'full-checks.json').read_text())
assert full['candidates_byte_identical'] and full['vcf_rows_identical']
assert full['changed_primary_hp_ps']==full['changed_parental_status']==0
assert full['before_blocks']==full['after_blocks']
assert full['connected_core_tags_unchanged'] and full['previously_phased_tags_unchanged']
assert full['before_parental_counts']==full['after_parental_counts']
panel=json.loads((OUT/'panel-comparison.json').read_text())
assert panel['panel_windows']==141 and panel['changed_windows']==0 and panel['summary_identical']
for log in ['build.log','fast.log']:
    assert 'warning:' not in (OUT/log).read_text() and 'error:' not in (OUT/log).read_text(),log
for path in ['physical-checks.json','parental-checks.json','owner-checks.json','group-mutation-checks.json']:
    assert (OUT/path).is_file()
physical=json.loads((OUT/'physical-checks.json').read_text())
assert physical['physical_contexts']==70 and physical['physical_molecule_records']==3971
assert physical['legacy_projections_identical'] and physical['raw_cohorts_and_costs_match_independently_verified_bam_state']
assert 'All tests passed' in (OUT/'standard-window-tests.log').read_text()
files=[]
selected=[ROOT/'pgphase',ROOT/'test_allele_context',ROOT/'test_allele_genotype',ROOT/'docs/IMPLEMENTATION.md',ROOT/'CHECKPOINT.md']
selected+=list((ROOT/'src').glob('allele_*.?pp'))+list((ROOT/'src').glob('test_allele_*.cpp'))+[ROOT/'src/collect_pipeline.cpp']
selected += [ROOT/p for p in unchanged]
for directory in [OUT,WORK/'accepted',WORK/'cached',WORK/'mutants',WORK/'group-mutants',WORK/'timing']:
    selected += [p for p in directory.rglob('*') if p.is_file() and '__pycache__' not in p.parts and p.name!='manifest.json']
selected += [WORK/'final-binary.sha256',WORK/'baseline/floors.sha256']
for p in sorted(set(selected)):
    files.append(dict(path=str(p.relative_to(ROOT)),size=p.stat().st_size,sha256=digest(p)))
manifest=dict(production_sha256=sha,starting_sha256=digest(WORK/'baseline/pgphase'),
              accepted_output='test_data/tmp_representation_step9/accepted/full-current',
              unchanged_floor_and_panel_files=unchanged,input_metadata=previous['input_metadata'],
              competitor_panel_signature=json.loads((OUT/'hiphase-verified.tsv.json').read_text())['signature'],
              baseline_independent_cohort_and_fit_manifest_sha256=digest(ROOT/'evaluations/2026-10-07-representation-step8/manifest.json'),
              files=files)
(OUT/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
print(f'Frozen {len(files)} files; all 107 checks, 141 panel contracts and phase outputs unchanged')
