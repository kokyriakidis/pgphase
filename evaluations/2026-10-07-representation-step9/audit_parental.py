#!/usr/bin/env python3
"""Project verified cached parental separability onto unique physical units."""
from collections import Counter, defaultdict
import csv
import hashlib
import json
from pathlib import Path

ROOT=Path(__file__).resolve().parents[2]
OUT=Path(__file__).resolve().parent
BASE=ROOT/'evaluations/2026-10-07-representation-step8'
manifest=json.loads((BASE/'manifest.json').read_text());known={r['path']:r['sha256'] for r in manifest['files']}
paths=[BASE/'heldout-parental-audit.tsv',ROOT/'test_data/tmp_representation_step8/hiphase-molecule-tags.json']
for path in paths:assert hashlib.sha256(path.read_bytes()).hexdigest()==known[str(path.relative_to(ROOT))],path
for record in manifest['input_metadata']:
 stat=(ROOT/record['path']).stat();assert stat.st_size==record['size'] and stat.st_mtime_ns==record['mtime_ns']
_,hi=json.loads(paths[1].read_text())
with paths[0].open() as handle:cached={(r['owner'],int(r['locus']),r['read']):r for r in csv.DictReader(handle,delimiter='\t')}
reports=[];names=set()
for owner in ['matrix4-final','matrix5-final','matrix65-final','whole-final']:
 folder=ROOT/'test_data/tmp_representation_step9/accepted'/owner
 groups=defaultdict(list)
 with (folder/'matrix.chunk0.physical-members.tsv').open() as h:
  for r in csv.DictReader(h,delimiter='\t'):groups[int(r['physical'])].append(int(r['locus']))
 outcomes=Counter();partitions=Counter();hiphase_correct=0;count=0
 with (folder/'matrix.chunk0.physical-cohort.tsv').open() as h:
  for r in csv.DictReader(h,delimiter='\t'):
   aliases=groups[int(r['physical'])];name=r['read'];names.add(name)
   reference=cached[owner,aliases[0],name]
   for alias in aliases[1:]:
    row=cached[owner,alias,name];assert (row['outcome'],row['partition'])==(reference['outcome'],reference['partition'])
   outcomes[reference['outcome']]+=1;partitions[reference['partition']]+=1
   hiphase_correct+=hi[name]=='correct';count+=1
 report=dict(owner=owner,physical_molecule_records=count,heldout_evaluation_only_outcomes=dict(outcomes),partitions=dict(partitions),identical_physical_records_hiphase_correct=hiphase_correct)
 reports.append(report);print(json.dumps(report),flush=True)
assert names==set(hi)
(OUT/'parental-checks.json').write_text(json.dumps(dict(interpretation='Cached independent step-8 parental separability projected onto exactly equal deduplicated physical input and full/held-out fits; not production phasing or independent folds',verified_distinct_molecules=len(names),owners=reports),indent=2)+'\n')
