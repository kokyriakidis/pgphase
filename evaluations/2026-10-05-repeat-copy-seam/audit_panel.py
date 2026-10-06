#!/usr/bin/env python3
"""Preserve prior native replay assignments and variants across the gap fix."""
import argparse
import importlib.util
import hashlib
import re
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
OUT = Path(__file__).resolve().parent
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--after', type=Path, default=Path('/tmp/pgphase-window-tests'))
args = parser.parse_args()
spec = importlib.util.spec_from_file_location('audit', OUT/'audit_preservation.py')
audit = importlib.util.module_from_spec(spec)
spec.loader.exec_module(audit)
truth = {f[0]: f[1] == 'PATERNAL'
         for line in (ROOT/'test_data/derived/chr20_truth_hap.tsv').open()
         if len(f := line.rstrip().split('\t')) == 2}


reports = []
missing = []
for before_bam in sorted(Path('/tmp/pgphase72-final').rglob('phased.bam')):
    if '.replay-cache' in before_bam.parts: continue
    relative = before_bam.relative_to(Path('/tmp/pgphase72-final'))
    after_bam = args.after/relative
    if not after_bam.exists():
        missing.append(str(relative))
        continue
    match = re.search(r'(/\S*\.gap-replay-cache/[a-f0-9]{64})/candidates.tsv', (after_bam.parent/'stderr.log').read_text())
    assert match, after_bam
    state = json.loads((Path(match.group(1))/'state.json').read_text())
    assert state['signature']['binary'] == hashlib.sha256((ROOT/'pgphase').read_bytes()).hexdigest(), after_bam
    result = audit.compare(before_bam.parent, after_bam.parent, truth)
    reports.append({'baseline': str(before_bam), 'output': str(relative), **result})
assert len(reports) == 220, f'Expected 220 task-start native outputs, got {len(reports)}'
assert not missing, missing
result = {'native_replays_compared': len(reports), 'unmatched_prior_replays': missing,
          'retained_variant_alleles_counts_and_filters_preserved': True,
          'correct_assignments_preserved': True,
          'phased_assignments_preserved': not any(r['lost_phased_assignments'] for r in reports),
          'changed_outputs': [r for r in reports if r['changed_read_tags'] or r['added_recovery_variants']]}
(OUT/'panel-audit.json').write_text(json.dumps(result, indent=2)+'\n')
print(json.dumps({k: v for k, v in result.items() if k != 'changed_outputs'}, indent=2))
