#!/usr/bin/env python3
"""Summarize independent regions without assuming successful recovery."""
import csv, json, shutil
from collections import defaultdict
from pathlib import Path
import pysam
REPORT=Path(__file__).resolve().parent
OUT=Path('/tmp/pgphase-site-observation-recovery')
def phase(path):
    with pysam.AlignmentFile(path) as bam:
        return {r.query_name:(r.get_tag('HP'),r.get_tag('PS')) for r in bam if r.has_tag('HP') and r.has_tag('PS') and r.get_tag('HP') in (1,2)}
results=[]; invariants={}
for row in json.loads((REPORT/'manifest.json').read_text()):
    name=row['name']; case=OUT/name; dest=REPORT/name; dest.mkdir(exist_ok=True)
    initial=phase(case/'clean.bam'); after=phase(case/'auto.bam'); transforms=defaultdict(set)
    for read,(hp,ps) in initial.items():
        assert read in after,(name,'lost read',read)
        nh,np=after[read]; transforms[ps].add((np,int(nh!=hp)))
    assert all(len(x)==1 for x in transforms.values()),(name,'mixed block')
    tiers=list(csv.DictReader((case/'auto.tiers.tsv').open(),delimiter='\t')); gaps=defaultdict(list)
    for tier in tiers: gaps[tier['GAP_LEFT'],tier['GAP_RIGHT']].append(tier)
    for gap,records in gaps.items():
        assert [int(r['TIER']) for r in records]==list(range(1,len(records)+1))
        assert all(r['STATUS']!='joined' for r in records[:-1])
        assert all(int(r['TIER'])==3 or int(r['MSA_HET_INDELS'])==0 for r in records)
        assert [int(r['MSA_HET_SNPS']) for r in records]==sorted(int(r['MSA_HET_SNPS']) for r in records)
    invariants[name]={'preserved_reads':len(initial),'original_blocks':len(transforms),'uniform_transforms':True,'gaps_attempted':len(gaps),'gaps_joined':sum(r[-1]['STATUS']=='joined' for r in gaps.values()),'new_tagged_reads':len(after.keys()-initial.keys())}
    shutil.copy(case/'auto.tiers.tsv',dest/'auto.tiers.tsv')
    for arm in ('clean','auto'):
        s=json.loads((case/f'{arm}.eval/summary.json').read_text())
        shutil.copy(case/f'{arm}.eval/summary.json',dest/f'{arm}.summary.json')
        shutil.copy(case/f'{arm}.command.sh',dest/f'{arm}.command.sh')
        results.append([name,arm]+[s[k] for k in ('total_reads_evaluated','discordant_reads','total_phase_sets','switchflip_errors')])
(REPORT/'invariants.json').write_text(json.dumps(invariants,indent=2)+'\n')
with (REPORT/'results.tsv').open('w') as f:
    w=csv.writer(f,delimiter='\t',lineterminator='\n'); w.writerow(['region','arm','evaluated_reads','discordant_reads','phase_sets','read_switchflips']); w.writerows(results)
print((REPORT/'results.tsv').read_text())
print(json.dumps(invariants,indent=2))

baseline = REPORT.parent / '2026-09-13-auto-gap-validation'
comparison = []
for row in json.loads((REPORT/'manifest.json').read_text()):
    name = row['name']
    before = json.loads((baseline/name/'auto.summary.json').read_text())
    after = json.loads((REPORT/name/'auto.summary.json').read_text())
    for key in ('discordant_reads', 'switchflip_errors'):
        assert after[key] <= before[key], (name, key, 'truth regression')
    comparison.append([name, before['total_phase_sets'], after['total_phase_sets'],
                       before['total_reads_evaluated'], after['total_reads_evaluated'],
                       before['discordant_reads'], after['discordant_reads']])
with (REPORT/'comparison.tsv').open('w') as f:
    w=csv.writer(f,delimiter='\t',lineterminator='\n')
    w.writerow(['region','previous_blocks','new_blocks','previous_reads','new_reads',
                'previous_discordant','new_discordant'])
    w.writerows(comparison)
print((REPORT/'comparison.tsv').read_text())
