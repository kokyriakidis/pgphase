#!/usr/bin/env python3
"""Audit full-chromosome read transitions and unrelated variant evidence."""
from collections import Counter
import importlib.util
import json
from pathlib import Path
import pysam

ROOT = Path(__file__).resolve().parents[2]
OUT = Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location('audit', OUT/'read_audit.py')
audit = importlib.util.module_from_spec(spec)
spec.loader.exec_module(audit)
truth = {f[0]: f[1] == 'PATERNAL' for line in (ROOT/'test_data/derived/chr20_truth_hap.tsv').open()
         if len(f := line.rstrip().split('\t')) == 2}
before = ROOT/'test_data/tmp_gap_fix89/final/0'
after = ROOT/'test_data/tmp_gap_fix90/final/0'
x, xs = audit.assignments(before/'phased.bam', truth)
y, ys = audit.assignments(after/'phased.bam', truth)
assert x.keys() <= y.keys()
with pysam.AlignmentFile(str(ROOT/'test_data/HG002_chr20_hifi_mapped_to_CHM13_chr20_annotated.bam')) as bam:
    gap = {read.query_name for read in bam.fetch('CHM13#0#chr20', 4618000, 4645000)
           if not read.is_secondary and not read.is_supplementary}
lost_correct = [q for q in xs if xs[q] == 'correct' and ys[q] != 'correct']
assert not lost_correct, lost_correct
changed_status = {q for q in xs if xs[q] != ys[q]}
outside_gap = sorted(changed_status-gap)
assert all(x[q][0] not in (1, 2) or x[q][1] <= 0 for q in outside_gap), outside_gap
assert all(y[q] == x[q] for q, (hp, ps) in x.items() if hp in (1, 2) and ps > 0)
assert all(y[q][0] in (1, 2) and y[q][1] > 0 for q, (hp, ps) in x.items() if hp in (1, 2) and ps > 0)
vx, vy = audit.variants(before/'phased.vcf'), audit.variants(after/'phased.vcf')
changed_evidence = []
for key in vx.keys() & vy.keys():
    info, sample = vx[key]
    ninfo, nsample = vy[key]
    if info != ninfo or {k: v for k, v in sample.items() if k not in ('GT', 'PS')} != {
            k: v for k, v in nsample.items() if k not in ('GT', 'PS')}:
        changed_evidence.append(key)
assert not changed_evidence, changed_evidence
removed, added = sorted(vx.keys()-vy.keys()), sorted(vy.keys()-vx.keys())
assert not removed, removed
assert not added, added
changed_records = [key for key in vx if vx[key] != vy[key]]
assert not changed_records, changed_records
result = {'primary_reads': len(y), 'before_counts': dict(Counter(xs.values())), 'after_counts': dict(Counter(ys.values())),
          'lost_correct_assignments': lost_correct, 'lost_correct_changed_tags': [], 'lost_phased_assignments': 0,
          'changed_read_tags': sum(x[q] != y[q] for q in x), 'changed_evidence': changed_evidence,
          'changed_variant_records': changed_records, 'variant_records': len(vx),
          'changed_tags_with_status_change_outside_gap': outside_gap, 'outside_gap_transitions': dict(Counter(f'{xs[q]}->{ys[q]}' for q in outside_gap)), 'previous_phased_tags_identical': True, 'removed_variant_keys': removed, 'added_variant_keys': added,
          'transitions': dict(Counter(f'{xs[q]}->{ys[q]}' for q in xs)),
          'newly_correct': sorted(q for q in xs if xs[q] != 'correct' and ys[q] == 'correct')}
(OUT/'full-preservation.json').write_text(json.dumps(result, indent=2)+'\n')
print(json.dumps({k: v for k, v in result.items() if k != 'newly_correct'}, indent=2))
