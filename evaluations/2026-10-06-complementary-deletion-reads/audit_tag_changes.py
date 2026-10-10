#!/usr/bin/env python3
"""Locate every new full-chromosome assignment in the original primary BAM."""
from pathlib import Path
from collections import Counter
import importlib.util
import json
import pysam

ROOT = Path(__file__).resolve().parents[2]
OUT = Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location('audit', OUT/'read_audit.py')
a = importlib.util.module_from_spec(spec)
spec.loader.exec_module(a)
truth = {f[0]: f[1] == 'PATERNAL' for line in (ROOT/'test_data/derived/chr20_truth_hap.tsv').open()
         if len(f := line.rstrip().split('\t')) == 2}
x, xs = a.assignments(ROOT/'test_data/tmp_gap_fix89/final/0/phased.bam', truth)
y, ys = a.assignments(ROOT/'test_data/tmp_gap_fix90/final/0/phased.bam', truth)
changed = {q for q in x if x[q] != y[q]}
rows = []
regions = [(3000000,4000000),(4000000,5000000),(11000000,12000000),(19000000,20000000),
           (22000000,24000000),(41000000,42000000),(48000000,49000000)]
with pysam.AlignmentFile(str(ROOT/'test_data/HG002_chr20_hifi_mapped_to_CHM13_chr20_annotated.bam')) as bam:
    for left, right in regions:
        for r in bam.fetch('CHM13#0#chr20', left, right):
            q = r.query_name
            if q not in changed or r.is_secondary or r.is_supplementary: continue
            assert xs[q] == 'unphased' and ys[q] == 'correct'
            rows.append({'qname':q,'span':[r.reference_start+1,r.reference_end],
                         'before':x[q],'after':y[q],'parent':'PATERNAL' if truth[q] else 'MATERNAL',
                         'before_status':xs[q],'after_status':ys[q],
                         'target_neighborhood':r.reference_start < 4645000 and r.reference_end > 4618000})
assert {r['qname'] for r in rows} == changed
assert len(rows) == len(changed) == 40
(OUT/'tag-changes.json').write_text(json.dumps(rows,indent=2)+'\n')
print(Counter(('target' if r['target_neighborhood'] else (r['span'][0]//1000000)) for r in rows))
