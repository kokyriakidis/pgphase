#!/usr/bin/env python3
"""Preserve every previous correct/phased read and unrelated variant evidence."""
from collections import Counter
import importlib.util
import json
from pathlib import Path
ROOT = Path(__file__).resolve().parents[2]
OUT = Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location('audit', ROOT/'evaluations/2026-10-05-largest-hiphase-block/audit_reads.py')
a = importlib.util.module_from_spec(spec)
spec.loader.exec_module(a)
truth = {f[0]: f[1] == 'PATERNAL' for line in (ROOT/'test_data/derived/chr20_truth_hap.tsv').open()
         if len(f := line.rstrip().split('\t')) == 2}
before = ROOT/'test_data/tmp_gap_fix76/final/0'
after = ROOT/'test_data/tmp_gap_fix77/final/0'
x, xs = a.assignments(before/'phased.bam', truth)
y, ys = a.assignments(after/'phased.bam', truth)
assert x.keys() == y.keys()
assert all(ys[q] == 'correct' for q in xs if xs[q] == 'correct')
assert all(y[q][0] in (1, 2) and y[q][1] > 0 for q, (hp, ps) in x.items() if hp in (1, 2) and ps > 0)
vx, vy = a.variants(before/'phased.vcf'), a.variants(after/'phased.vcf')
changed_evidence = []
for key in vx.keys() & vy.keys():
    info, sample = vx[key]
    ninfo, nsample = vy[key]
    if info != ninfo or {k: v for k, v in sample.items() if k not in ('GT', 'PS')} != {
            k: v for k, v in nsample.items() if k not in ('GT', 'PS')}:
        changed_evidence.append(key)
assert all(k[1] == '40636354' and k[3:5] == ('G', 'T') for k in changed_evidence), changed_evidence
removed = sorted(vx.keys()-vy.keys())
added = sorted(vy.keys()-vx.keys())
assert all(k[1] == '60171648' and k[3:5] == ('T', 'C') for k in removed), removed
assert all(k[1] == '40636353' and k[3:5] == ('TG', 'T') for k in added), added
r = {'primary_reads': len(y), 'before_counts': dict(Counter(xs.values())), 'after_counts': dict(Counter(ys.values())),
     'lost_correct_assignments': 0, 'lost_phased_assignments': 0,
     'changed_read_tags': sum(x[q] != y[q] for q in x), 'changed_evidence': changed_evidence,
     'removed_variant_keys': removed, 'added_variant_keys': added,
     'newly_correct': [q for q in xs if xs[q] != 'correct' and ys[q] == 'correct']}
(OUT/'full-preservation.json').write_text(json.dumps(r, indent=2)+'\n')
print(json.dumps(r, indent=2))
