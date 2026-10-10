#!/usr/bin/env python3
"""Freeze step-12 source, accepted evidence and unchanged regression contracts."""
import hashlib
import json
from pathlib import Path

ROOT=Path(__file__).resolve().parents[2]
OUT=Path(__file__).resolve().parent
WORK=ROOT/'test_data/tmp_representation_step12'

def digest(path):
    h=hashlib.sha256()
    with path.open('rb') as stream:
        while block:=stream.read(1024*1024):h.update(block)
    return h.hexdigest()

previous=json.loads((ROOT/'evaluations/2026-10-07-representation-step11/manifest.json').read_text())
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
assert sha==digest(WORK/'candidate-pgphase')==(WORK/'final-binary.sha256').read_text().split()[0]
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
for log in ['build-fast.log','unit-final.log','golden-checks.log']:
    text=(OUT/log).read_text()
    assert 'warning:' not in text and 'error:' not in text,log
maps=json.loads((OUT/'map-checks.json').read_text())['owners']
assert sum(r['contexts'] for r in maps)==70 and sum(r['parents'] for r in maps)==734
assert all(r['parents']==r['complete'] and r['previous_tables_byte_identical'] and r['raw_replay_byte_identical'] for r in maps)
assert sum(r['multiple_optimal_alignments'] for r in maps)==629
assert sum(r['certified_matched_bases'] for r in maps)==26017
mutants=json.loads((OUT/'mutation-checks.json').read_text())
assert len(mutants)==8 and all(r['exit']!=0 and r['failed_checks']>0 for r in mutants.values())
raw=json.loads((OUT/'raw-state-checks.json').read_text())
assert raw['missing_raw_catalog_rejected'] and all(r['byte_identical'] and r['precomputed_outputs_absent'] for r in raw['owners'])
deficit=json.loads((OUT/'length-deficit-checks.json').read_text())
assert deficit['molecules']==39 and deficit['same_hiphase_input_alignments'] and deficit['raw_query_slices_and_cigar_length_balance_verified']
assert deficit['unavoidable_length_residual']==4327 and deficit['path_residual']==4383
assert 'All tests passed' in (OUT/'standard-window-tests-final.log').read_text()
hi=json.loads((OUT/'hiphase-verified.tsv.json').read_text())
assert hi['signature']==previous['competitor_panel_signature']
assert digest(OUT/'hiphase-verified.tsv')==hi['output_sha256']
selected=[ROOT/'pgphase',ROOT/'test_allele_context',ROOT/'test_allele_genotype',ROOT/'docs/IMPLEMENTATION.md',ROOT/'CHECKPOINT.md',ROOT/'src/collect_pipeline.cpp']
selected+=list((ROOT/'src').glob('allele_*.?pp'))+list((ROOT/'src').glob('test_allele_*.cpp'))+[ROOT/p for p in unchanged]
for directory in [OUT,WORK/'accepted-final',WORK/'cached',WORK/'raw-state-only',WORK/'mutants',WORK/'timing',WORK/'baseline',WORK/'sources']:
    selected += [p for p in directory.rglob('*') if p.is_file() and '__pycache__' not in p.parts and p.name!='manifest.json']
selected += [WORK/'candidate-pgphase',WORK/'final-binary.sha256']
files=[dict(path=str(p.relative_to(ROOT)),size=p.stat().st_size,sha256=digest(p)) for p in sorted(set(selected))]
manifest=dict(production_sha256=sha,starting_sha256=digest(WORK/'baseline/pgphase'),
    accepted_output=str((WORK/'accepted-final/full-frozen').relative_to(ROOT)),
    unchanged_floor_and_panel_files=unchanged,input_metadata=previous['input_metadata'],
    competitor_panel_signature=hi['signature'],baseline_manifest_sha256=digest(ROOT/'evaluations/2026-10-07-representation-step11/manifest.json'),files=files)
(OUT/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
print(f'Frozen {len(files)} files; 107 checks, 141 panel contracts and phase outputs unchanged')
