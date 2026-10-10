#!/usr/bin/env python3
"""Measure best-sequence residual lower bounds; these are not diploid fits."""
import csv
import json
from pathlib import Path

ROOT=Path(__file__).resolve().parents[2];OUT=Path(__file__).resolve().parent
WORK=ROOT/'test_data/tmp_representation_step10'
reports=[]
for owner in ['matrix4-final','matrix5-final','matrix65-final','whole-final']:
 def rows(path):return list(csv.DictReader(path.open(),delimiter='\t'))
 base=WORK/'accepted'/owner
 original={(r['physical'],r['read']):r['costs'] for r in rows(base/'matrix.chunk0.physical-costs.tsv')}
 measured=rows(WORK/'cached'/f'{owner}.residuals.tsv')
 assert len(measured)==len(original) and {(r['physical'],r['read']) for r in measured}==set(original)
 totals={}
 for r in measured:
  assert r['old_costs']==original[r['physical'],r['read']]
  old,new=int(r['old_min']),int(r['composed_min'])
  assert old==min(map(int,r['old_costs'].split(','))) and 0<=new<=old
  ctx=totals.setdefault(r['physical'],dict(molecules=0,old_residual=0,composed_residual=0,improved_molecules=0))
  ctx['molecules']+=1;ctx['old_residual']+=old;ctx['composed_residual']+=new;ctx['improved_molecules']+=new<old
 report=dict(owner=owner,**{k:sum(v[k] for v in totals.values()) for k in ['molecules','old_residual','composed_residual','improved_molecules']},physical_contexts=totals)
 reports.append(report);print(json.dumps({k:v for k,v in report.items() if k!='physical_contexts'}),flush=True)
(OUT/'residual-checks.json').write_text(json.dumps(dict(interpretation='Unconstrained per-read best-sequence edit-cost lower bound; no genotype, error confidence, phase decision or closure claim.',owners=reports),indent=2)+'\n')
