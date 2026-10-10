#!/usr/bin/env python3
"""Prove composition replay needs raw state, not precomputed composition outputs."""
import json
from pathlib import Path
import resource
import shutil
import subprocess

ROOT=Path(__file__).resolve().parents[2];OUT=Path(__file__).resolve().parent
WORK=ROOT/'test_data/tmp_representation_step12'
inputs=['composition-contexts','site-catalog','physical-alleles','physical-members','joint-parents']
outputs=['composition-sites','compositions','composed-alleles','path-status','sequence-paths','path-alleles','parent-maps','matched-subpaths']
reports=[]
for owner in ['matrix4-final','matrix5-final','matrix65-final','whole-final']:
 source=WORK/'accepted-final'/owner;state=WORK/'raw-state-only'/owner;state.mkdir(parents=True,exist_ok=True)
 for kind in inputs:shutil.copyfile(source/f'matrix.chunk0.{kind}.tsv',state/f'matrix.chunk0.{kind}.tsv')
 assert {p.name for p in state.glob('matrix.chunk0.*')}=={f'matrix.chunk0.{kind}.tsv' for kind in inputs}
 run=subprocess.run(['./test_allele_context','--maps',str(state),str(state/'replay')],cwd=ROOT,capture_output=True,text=True,check=True)
 for kind in outputs:assert (state/f'replay.{kind}.tsv').read_bytes()==(source/f'matrix.chunk0.{kind}.tsv').read_bytes()
 reports.append(dict(owner=owner,raw_tables=len(inputs),result=run.stdout.strip(),precomputed_outputs_absent=True,byte_identical=True))
negative=WORK/'raw-state-only/missing-catalog';negative.mkdir(exist_ok=True)
for kind in inputs:
 if kind!='site-catalog':shutil.copyfile(WORK/'accepted-final/matrix4-final'/f'matrix.chunk0.{kind}.tsv',negative/f'matrix.chunk0.{kind}.tsv')
def disable_core():resource.setrlimit(resource.RLIMIT_CORE,(0,0))
run=subprocess.run(['./test_allele_context','--maps',str(negative),str(negative/'replay')],cwd=ROOT,capture_output=True,text=True,preexec_fn=disable_core)
assert run.returncode!=0 and 'cannot read composition state:' in run.stderr and 'site-catalog' in run.stderr
(negative/'rejection.log').write_text(run.stderr)
(OUT/'raw-state-checks.json').write_text(json.dumps(dict(owners=reports,missing_raw_catalog_rejected=True,missing_catalog_exit=run.returncode),indent=2)+'\n')
print(json.dumps(reports,indent=2));print('Missing raw catalog rejected')
