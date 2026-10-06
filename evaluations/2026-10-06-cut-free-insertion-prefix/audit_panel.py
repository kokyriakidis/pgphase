#!/usr/bin/env python3
"""Report changes to every task-start native replay, without hiding label losses."""
from collections import Counter
import hashlib
import importlib.util
import json
from pathlib import Path
import re

ROOT = Path(__file__).resolve().parents[2]
OUT = Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location('audit', ROOT/'evaluations/2026-10-05-largest-hiphase-block/audit_reads.py')
audit = importlib.util.module_from_spec(spec)
spec.loader.exec_module(audit)
truth = {f[0]: f[1] == 'PATERNAL' for line in (ROOT/'test_data/derived/chr20_truth_hap.tsv').open()
         if len(f := line.rstrip().split('\t')) == 2}
binary = hashlib.sha256((ROOT/'pgphase').read_bytes()).hexdigest()
reports = []
compared = {}
for before in sorted(Path('/tmp/pgphase75-final').rglob('phased.bam')):
    if '.replay-cache' in before.parts:
        continue
    relative = before.relative_to('/tmp/pgphase75-final')
    after = Path('/tmp/pgphase-window-tests')/relative
    assert after.exists(), after
    paths = []
    for path in (before, after):
        match = re.search(r'(/\S*\.gap-replay-cache/[a-f0-9]{64})/candidates.tsv', (path.parent/'stderr.log').read_text())
        assert match, path
        paths.append(Path(match.group(1)))
    state = json.loads((paths[1]/'state.json').read_text())
    assert state['signature']['binary'] == binary, after
    key = tuple(map(str, paths))
    if key not in compared:
        old_tags, old = audit.assignments(paths[0]/'phased.bam', truth)
        new_tags, new = audit.assignments(paths[1]/'phased.bam', truth)
        assert old_tags.keys() == new_tags.keys(), after
        old_vcf = audit.variants(paths[0]/'phased.vcf')
        new_vcf = audit.variants(paths[1]/'phased.vcf')
        changed = sum(old_tags[q] != new_tags[q] for q in old_tags)
        compared[key] = {
            'before_counts': dict(Counter(old.values())), 'after_counts': dict(Counter(new.values())),
            'changed_read_tags': changed,
            'lost_correct_assignments': [q for q in old if old[q] == 'correct' and new[q] != 'correct'],
            'changed_variant_records': sum(old_vcf[k] != new_vcf.get(k) for k in old_vcf),
            'missing_variant_keys': sorted(old_vcf.keys()-new_vcf.keys()),
            'added_variant_keys': sorted(new_vcf.keys()-old_vcf.keys())}
    reports.append({'output': str(relative), **compared[key]})
assert len(reports) == 223, len(reports)
result = {'native_output_labels': len(reports), 'independent_requests': len(compared),
          'binary_sha256': binary, 'verified_current_binary_states': True,
          'changed_outputs': [r for r in reports if r['changed_read_tags'] or r['changed_variant_records'] or r['added_variant_keys']]}
(OUT/'panel-audit.json').write_text(json.dumps(result, indent=2)+'\n')
print(json.dumps({'labels': len(reports), 'requests': len(compared), 'changed_labels': len(result['changed_outputs'])}))
