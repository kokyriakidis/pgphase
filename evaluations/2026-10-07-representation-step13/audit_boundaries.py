#!/usr/bin/env python3
"""Classify the 30 whole-REF-contained overlap attempts identified in step 12."""
import csv
import json
from collections import Counter, defaultdict
from pathlib import Path

ROOT=Path(__file__).resolve().parents[2];OUT=Path(__file__).resolve().parent
reports=[]
for owner in ['matrix4-final','matrix5-final','matrix65-final','whole-final']:
    folder=ROOT/'test_data/tmp_representation_step13/accepted-final'/owner
    def rows(kind):
        with (folder/f'matrix.chunk0.{kind}.tsv').open() as f:return list(csv.DictReader(f,delimiter='\t'))
    contexts={r['physical']:r for r in rows('composition-contexts')}
    spans=defaultdict(list)
    for r in rows('matched-subpaths'):
        spans[r['physical'],r['parent_allele']].append(tuple(int(r[k]) for k in ['ref_beg','alt_beg','length']))
    sites={(r['physical'],r['candidate'],r['alt_index']):r for r in rows('composition-sites')}
    outcomes=Counter();details=[]
    for old,new in zip(rows('compositions'),rows('nested-compositions')):
        if old['status']!='overlap':continue
        site=sites[old['physical'],old['candidate'],old['alt_index']]
        offset=int(site['pos'])-int(contexts[old['physical']]['beg']);length=len(site['ref'])
        contained=[s for s in spans[old['physical'],old['parent_allele']] if s[0]<=offset and offset+length<=s[0]+s[2]]
        if not contained:continue
        reason='resolved' if new['status']=='valid' else 'unguarded_length_change'
        assert reason=='resolved' or len(site['ref'])!=len(site['alt']) and not any(s[0]<offset and offset+length<s[0]+s[2] for s in contained)
        outcomes[reason]+=1
        details.append(dict(physical=old['physical'],parent=old['parent_allele'],candidate=old['candidate'],
            pos=site['pos'],ref=site['ref'],alt=site['alt'],reason=reason))
    reports.append(dict(owner=owner,outcomes=dict(outcomes),details=details))
assert sum(sum(r['outcomes'].values()) for r in reports)==30
assert sum(r['outcomes'].get('resolved',0) for r in reports)==7
(OUT/'boundary-checks.json').write_text(json.dumps(reports,indent=2)+'\n')
print([(r['owner'],r['outcomes']) for r in reports])
