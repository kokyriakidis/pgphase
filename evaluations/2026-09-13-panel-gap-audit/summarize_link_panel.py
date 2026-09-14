#!/usr/bin/env python3
"""Compare exact target endpoints, read truth, and preserved block orientations."""
import argparse,csv,gzip,json
from collections import Counter,defaultdict
from pathlib import Path
import pysam
p=argparse.ArgumentParser(description=__doc__)
p.add_argument('--label',required=True)
p.add_argument('--baseline-label')
a=p.parse_args()
REPORT=Path(__file__).resolve().parent/a.label
OUT=Path('/tmp')/('pgphase-'+a.label)

def tags(path):
    with pysam.AlignmentFile(str(path)) as bam:
        return {r.query_name:(r.get_tag('PS'),r.get_tag('HP')) for r in bam if r.has_tag('PS') and r.has_tag('HP') and r.get_tag('HP') in (1,2)}

def target(path,left,right):
    ps={}
    with pysam.VariantFile(str(path)) as vcf:
        for r in vcf:
            s=r.samples[0]
            if r.pos in (left,right) and s.phased:ps[r.pos]=s.get('PS')
    if ps.get(left) is not None and ps.get(left)==ps.get(right):return 'joined'
    return 'split' if len(ps)==2 else 'endpoint_unphased'
rows=[];tiers=[];baseline_changes=[];block_orientations=[]
for r in json.loads((REPORT/'manifest.json').read_text()):
    case=OUT/r['name'];row={'region':r['name'],'left':r['left'],'right':r['right']}
    for arm in ('clean','recovery2'):
        d=case/arm;s=json.loads((d/'read_eval/summary.json').read_text())
        row[arm+'_target']=target(d/'shared.vcf',r['left'],r['right'])
        for k in ('total_reads_evaluated','discordant_reads','total_phase_sets','switchflip_errors'):row[arm+'_'+k]=s[k]
    initial=tags(case/'clean/phased.bam');after=tags(case/'recovery2/phased.bam');transforms=defaultdict(set)
    for name,(ps,hp) in initial.items():
        assert name in after,(r['name'],'lost read',name)
        nps,nhp=after[name];transforms[ps].add((nps,int(nhp!=hp)))
    assert all(len(t)==1 for t in transforms.values()),(r['name'],'nonuniform original block')
    # Separate a wrong relative block orientation from uncertain individual
    # reads. Truth labels evaluate the stitch; they never select its evidence.
    def read_truth(arm):
        with gzip.open(case/arm/'read_eval/per_read.tsv.gz','rt') as f:
            return {v['read_name']:v for v in csv.DictReader(f,delimiter='\t')}
    before_truth=read_truth('clean');after_truth=read_truth('recovery2')
    original_blocks=defaultdict(list)
    for name,(ps,hp) in initial.items():
        if name not in before_truth or name not in after_truth:continue
        before=before_truth[name]['status'];after_status=after_truth[name]['status']
        if before=='skipped' or after_status=='skipped':continue
        original_blocks[ps].append((before=='DISCORDANT',after_status=='DISCORDANT'))
    for ps,observations in original_blocks.items():
        n=len(observations);before=sum(v[0] for v in observations);after_errors=sum(v[1] for v in observations)
        final_ps,flipped=next(iter(transforms[ps]))
        block_orientations.append({'region':r['name'],'original_ps':ps,'final_ps':final_ps,
            'truth_reads':n,'before_discordant':before,'after_discordant':after_errors,
            'uniform_label_flip':flipped,'majority_reversed':before*2<n and after_errors*2>n})
    row['original_reads_preserved']=len(initial);row['added_tagged_reads']=len(after.keys()-initial.keys())
    for t in csv.DictReader((case/'recovery2/tiers.tsv').open(),delimiter='\t'):tiers.append(dict(region=r['name'],**t))
    rows.append(row)
    if a.baseline_label:
        previous=Path('/tmp')/('pgphase-'+a.baseline_label)/r['name']/'recovery2'
        s=json.loads((previous/'read_eval/summary.json').read_text())
        baseline_changes.append({'region':r['name'],'before_target':target(previous/'shared.vcf',r['left'],r['right']),'after_target':row['recovery2_target'],'before_reads':s['total_reads_evaluated'],'after_reads':row['recovery2_total_reads_evaluated'],'before_discordant':s['discordant_reads'],'after_discordant':row['recovery2_discordant_reads'],'before_switchflips':s['switchflip_errors'],'after_switchflips':row['recovery2_switchflip_errors']})
for filename,data in [('results.tsv',rows),('tiers.tsv',tiers),('baseline_comparison.tsv',baseline_changes),('block_orientations.tsv',block_orientations)]:
    if not data:continue
    with (REPORT/filename).open('w') as f:
        w=csv.DictWriter(f,fieldnames=list(data[0]),delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(data)
summary={'cases':len(rows),'target_transitions':dict(Counter(r['clean_target']+' -> '+r['recovery2_target'] for r in rows)),'original_reads_preserved':sum(r['original_reads_preserved'] for r in rows),'added_tagged_reads':sum(r['added_tagged_reads'] for r in rows),'discordance_increases':[r['region'] for r in rows if r['recovery2_discordant_reads']>r['clean_discordant_reads']]}
summary['majority_reversed_original_blocks']=[r for r in block_orientations if r['majority_reversed']]
if baseline_changes:
    summary['baseline_target_transitions']=dict(Counter(r['before_target']+' -> '+r['after_target'] for r in baseline_changes))
    summary['baseline_discordance_increases']=[r for r in baseline_changes if r['after_discordant']>r['before_discordant']]
    summary['baseline_switchflip_increases']=[r for r in baseline_changes if r['after_switchflips']>r['before_switchflips']]
    summary['baseline_discordance_decreases']=[r for r in baseline_changes if r['after_discordant']<r['before_discordant']]
(REPORT/'summary.json').write_text(json.dumps(summary,indent=2)+'\n');print(json.dumps(summary,indent=2))
