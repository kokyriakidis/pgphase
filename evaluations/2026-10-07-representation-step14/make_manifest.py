#!/usr/bin/env python3
"""Freeze step-14 source, accepted evidence and unchanged regression contracts."""
import hashlib
import json
from pathlib import Path

ROOT=Path(__file__).resolve().parents[2]
OUT=Path(__file__).resolve().parent
WORK=ROOT/'test_data/tmp_representation_step14'

def digest(path):
    h=hashlib.sha256()
    with path.open('rb') as stream:
        while block:=stream.read(1024*1024):h.update(block)
    return h.hexdigest()

previous=json.loads((ROOT/'evaluations/2026-10-07-representation-step13/manifest.json').read_text())
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
catalogs=json.loads((OUT/'read-catalog-checks.json').read_text())
assert sum(r['contexts'] for r in catalogs)==70
assert sum(r['read_context_rows'] for r in catalogs)==3971
assert sum(r['hypotheses'] for r in catalogs)==255
assert sum(r['new_unique_sequences'] for r in catalogs)==63
assert all(r['previous_tables_byte_identical'] and r['raw_replay_byte_identical'] and r['independent_support_and_catalogs_verified'] and r['statuses_complete'] for r in catalogs)
residual=json.loads((OUT/'residual-checks.json').read_text())
assert all(r['discovery_heldout_residual']<=r['old_residual'] and r['independent_scalar_dp_matches'] for r in residual['owners'])
assert residual['target']['old_residual']==4383 and residual['target']['discovery_heldout_residual']==715
msa=json.loads((OUT/'msa-checks.json').read_text())
assert msa['reversed_order_consensuses_and_membership_identical'] and msa['long_supporter_exclusions_verified']==6
assert msa['all_msa_input_rows_preserved'] and msa['long_consensus_length']==312
assert not msa['production_uses_msa_results']
mutants=json.loads((OUT/'mutation-checks.json').read_text())
assert len(mutants)==8 and all(r['exit']!=0 and r['failed_checks']>0 for r in mutants.values())
raw=json.loads((OUT/'raw-state-checks.json').read_text())
assert raw['missing_raw_catalog_rejected'] and all(r['byte_identical'] and r['precomputed_outputs_absent'] for r in raw['owners'])
assert 'All tests passed' in (OUT/'standard-window-tests-final.log').read_text()
hi=json.loads((OUT/'hiphase-verified.tsv.json').read_text())
assert hi['signature']==previous['competitor_panel_signature']
assert digest(OUT/'hiphase-verified.tsv')==hi['output_sha256']
selected=[ROOT/'pgphase',ROOT/'test_allele_context',ROOT/'test_allele_genotype',ROOT/'docs/IMPLEMENTATION.md',ROOT/'CHECKPOINT.md',ROOT/'src/collect_pipeline.cpp']
selected+=list((ROOT/'src').glob('allele_*.?pp'))+list((ROOT/'src').glob('test_allele_*.cpp'))+[ROOT/p for p in unchanged]
for directory in [OUT,WORK/'accepted-final',WORK/'cached',WORK/'raw-state-only',WORK/'mutants',WORK/'timing',WORK/'baseline',WORK/'sources']:
    selected += [p for p in directory.rglob('*') if p.is_file() and '__pycache__' not in p.parts and p.name!='manifest.json']
selected += [WORK/'candidate-pgphase',WORK/'final-binary.sha256',WORK/'msa_verify',WORK/'score_read_catalog']
files=[dict(path=str(p.relative_to(ROOT)),size=p.stat().st_size,sha256=digest(p)) for p in sorted(set(selected))]
manifest=dict(production_sha256=sha,starting_sha256=digest(WORK/'baseline/pgphase'),
    accepted_output=str((WORK/'accepted-final/full-frozen').relative_to(ROOT)),
    unchanged_floor_and_panel_files=unchanged,input_metadata=previous['input_metadata'],
    competitor_panel_signature=hi['signature'],baseline_manifest_sha256=digest(ROOT/'evaluations/2026-10-07-representation-step13/manifest.json'),files=files)
(OUT/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
print(f'Frozen {len(files)} files; 107 checks, 141 panel contracts and phase outputs unchanged')
