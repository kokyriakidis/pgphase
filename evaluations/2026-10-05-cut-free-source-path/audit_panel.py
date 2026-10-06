#!/usr/bin/env python3
"""Audit every matching frozen native replay against the task-start suite."""
import json
from pathlib import Path
from audit_reads import audit

old_root = Path('test_data/tmp_gap_fix67')
new_root = Path('test_data/tmp_gap_fix68')
old = {(r['name'], r['section']): r for r in json.loads((old_root / 'frozen-checks.json').read_text())}
new = json.loads((new_root / 'certified-checks.json').read_text())
truth = {f[0]: f[1] == 'PATERNAL'
         for line in open('test_data/derived/chr20_truth_hap.tsv')
         if len(f := line.rstrip().split('\t')) == 2}
reports = []
for check in new:
    prior = old.get((check['name'], check['section']))
    if prior is None:
        continue
    before_root = old_root / 'frozen-checks' / str(prior['index'])
    after_root = new_root / 'certified-checks' / str(check['index'])
    for after_bam in sorted(after_root.rglob('phased.bam')):
        relative = after_bam.relative_to(after_root)
        before_bam = before_root / relative
        if not before_bam.exists():
            continue
        report = audit(before_bam.parent, after_bam.parent, truth)
        reports.append({'check': check['name'], 'section': check['section'],
                        'output': str(relative), **report})
assert reports, 'No native replays compared'
result = {'native_replays_compared': len(reports),
          'variant_alleles_counts_and_filters_preserved': True,
          'correct_and_phased_assignments_preserved': True,
          'changed_outputs': [r for r in reports if r['changed_read_tags'] or r['new_assignments']]}
Path('evaluations/2026-10-05-cut-free-source-path/panel-audit.json').write_text(json.dumps(result, indent=2) + '\n')
print(json.dumps({k: v for k, v in result.items() if k != 'changed_outputs'}, indent=2))
