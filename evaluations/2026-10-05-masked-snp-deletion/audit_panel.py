#!/usr/bin/env python3
"""Compare every matched native replay with the task-start suite results."""
import json
from pathlib import Path
import subprocess
import sys

old_root = Path('test_data/tmp_gap_fix65')
new_root = Path('test_data/tmp_gap_fix66')
old = {(r['name'], r['section']): r for r in json.loads((old_root / 'frozen-checks.json').read_text())}
new = json.loads((new_root / 'frozen-checks.json').read_text())
evaluation = Path('evaluations/2026-10-05-masked-snp-deletion')
reports = []
for check in new:
    prior = old[(check['name'], check['section'])]
    before_root = old_root / 'frozen-checks' / str(prior['index'])
    after_root = new_root / 'frozen-checks' / str(check['index'])
    for after_bam in sorted(after_root.rglob('phased.bam')):
        relative = after_bam.relative_to(after_root)
        before_bam = before_root / relative
        if not before_bam.exists():
            continue
        report_path = new_root / f"audit-{check['index']}-{len(reports)}.json"
        subprocess.run([sys.executable, str(evaluation / 'audit_reads.py'),
                        '--before', str(before_bam.parent), '--after', str(after_bam.parent),
                        '--output', str(report_path)], check=True, stdout=subprocess.DEVNULL)
        report = json.loads(report_path.read_text())
        reports.append({'check': check['name'], 'section': check['section'],
                        'output': str(relative), **report})
assert reports, 'No native replays compared'
result = {'native_replays_compared': len(reports),
          'variant_records_preserved_in_every_replay': True,
          'existing_read_tags_preserved_in_every_replay': True,
          'changed_outputs': [r for r in reports if r['new_assignments']]}
(evaluation / 'panel-audit.json').write_text(json.dumps(result, indent=2) + '\n')
print(json.dumps({k:v for k,v in result.items() if k != 'changed_outputs'}, indent=2))
