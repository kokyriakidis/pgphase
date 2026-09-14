#!/usr/bin/env python3
"""Compare one-read linking against the fixed two-read recovery baseline."""
import csv,json,shutil
from pathlib import Path
from collections import defaultdict
import pysam
REPORT=Path(__file__).resolve().parent
OUT=Path('/tmp/pgphase-singleton-link-validation')
BASE=REPORT.parent/'2026-09-13-site-observation-recovery'
def tags(path):
    with pysam.AlignmentFile(path) as f:
        return {r.query_name:(r.get_tag('HP'),r.get_tag('PS')) for r in f
                if r.has_tag('HP') and r.has_tag('PS') and r.get_tag('HP') in (1,2)}
results=[]; invariants={}
for m in json.loads((REPORT/'manifest.json').read_text()):
    name=m['name'];d=REPORT/name;d.mkdir(exist_ok=True);case=OUT/name
    old=json.loads((BASE/name/'auto.summary.json').read_text())
    new=json.loads((case/'auto.eval/summary.json').read_text())
    for src,dst in [(case/'auto.eval/summary.json',d/'summary.json'),(case/'auto.tiers.tsv',d/'tiers.tsv'),(case/'auto.command.sh',d/'command.txt')]:shutil.copy(src,dst)
    a=tags(Path('/tmp/pgphase-site-observation-recovery')/name/'auto.bam');b=tags(case/'auto.bam');transforms=defaultdict(set)
    lost=[]
    for n,(hp,ps) in a.items():
        if n not in b:lost.append(n);continue
        nh,np=b[n];transforms[ps].add((np,int(nh!=hp)))
    invariants[name]={'original_reads':len(a),'lost_reads':lost,'block_transforms':{k:sorted(v) for k,v in transforms.items()}}
    results.append([name]+[s[k] for k in ('total_phase_sets','total_reads_evaluated','discordant_reads','switchflip_errors') for s in (old,new)])
with (REPORT/'results.tsv').open('w') as f:
    w=csv.writer(f,delimiter='\t',lineterminator='\n');w.writerow(['region','blocks_link2','blocks_link1','reads_link2','reads_link1','discordant_link2','discordant_link1','switchflips_link2','switchflips_link1']);w.writerows(results)
(REPORT/'invariants.json').write_text(json.dumps(invariants,indent=2)+'\n')
print((REPORT/'results.tsv').read_text())
print(json.dumps(invariants,indent=2))
