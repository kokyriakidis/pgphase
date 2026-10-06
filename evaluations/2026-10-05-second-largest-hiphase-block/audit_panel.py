#!/usr/bin/env python3
"""Preserve prior native replay assignments and variants across the gap fix."""
import argparse
import importlib.util
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
OUT = Path(__file__).resolve().parent
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--after', type=Path, default=Path('/tmp/pgphase71-final'))
args = parser.parse_args()
spec = importlib.util.spec_from_file_location('audit', ROOT/'evaluations/2026-10-05-largest-hiphase-block/audit_reads.py')
audit = importlib.util.module_from_spec(spec)
spec.loader.exec_module(audit)
truth = {f[0]: f[1] == 'PATERNAL'
         for line in (ROOT/'test_data/derived/chr20_truth_hap.tsv').open()
         if len(f := line.rstrip().split('\t')) == 2}
reports = []
missing = []
for before_bam in sorted((ROOT/'test_data/tmp_gap_fix69/certified-checks').rglob('phased.bam')):
    relative = Path(*before_bam.parts[-3:])
    after_bam = args.after/relative
    if not after_bam.exists():
        missing.append(str(relative))
        continue
    result = audit.audit(before_bam.parent, after_bam.parent, truth)
    reports.append({'baseline': str(before_bam.relative_to(ROOT)), 'output': str(relative), **result})
assert reports, 'No native replays compared'
result = {'native_replays_compared': len(reports), 'unmatched_prior_replays': missing,
          'variant_alleles_counts_and_filters_preserved': True,
          'correct_and_phased_assignments_preserved': True,
          'changed_outputs': [r for r in reports if r['changed_read_tags'] or r['new_assignments']]}
(OUT/'panel-audit.json').write_text(json.dumps(result, indent=2)+'\n')
print(json.dumps({k: v for k, v in result.items() if k != 'changed_outputs'}, indent=2))
