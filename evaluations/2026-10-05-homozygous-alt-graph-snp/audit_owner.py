#!/usr/bin/env python3
"""Audit the native owner, label transitions, and one/four-thread identity."""
from collections import Counter
import hashlib
import importlib.util
import json
from pathlib import Path
import pysam

ROOT = Path(__file__).resolve().parents[2]
OUT = Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location('audit', ROOT/'evaluations/2026-10-05-largest-hiphase-block/audit_reads.py')
audit = importlib.util.module_from_spec(spec)
spec.loader.exec_module(audit)
truth = {f[0]: f[1] == 'PATERNAL' for line in (ROOT/'test_data/derived/chr20_truth_hap.tsv').open()
         if len(f := line.rstrip().split('\t')) == 2}
base = ROOT/'test_data/tmp_gap_fix75/baseline/36'
fixed = Path('/tmp/pgphase-window-tests/graph/w36614185-36623545')
single = ROOT/'test_data/tmp_gap_fix75/retired_single/36'
old_tags, old = audit.assignments(base/'phased.bam', truth)
new_tags, new = audit.assignments(fixed/'phased.bam', truth)
assert old_tags.keys() == new_tags.keys()
assert (fixed/'candidates.tsv').read_bytes() == (single/'candidates.tsv').read_bytes()
assert (fixed/'native.vcf').read_bytes() == (single/'phased.vcf').read_bytes()
assert (new_tags, new) == audit.assignments(single/'phased.bam', truth)
with pysam.AlignmentFile(str(ROOT/'test_data/HG002_chr20_hifi_mapped_to_CHM13_chr20_annotated.bam')) as bam:
    names = {r.query_name for r in bam.fetch('CHM13#0#chr20', 36614184, 36623545)
             if not r.is_secondary and not r.is_supplementary and r.query_name in truth}
result = {'binary_sha256': hashlib.sha256((ROOT/'pgphase').read_bytes()).hexdigest(),
          'one_four_thread_identity': True, 'primary_reads': len(new_tags),
          'before_counts': dict(Counter(old.values())), 'after_counts': dict(Counter(new.values())),
          'transitions': dict(Counter(f'{old[q]}->{new[q]}' for q in old)),
          'lost_correct_assignments': [q for q in old if old[q] == 'correct' and new[q] != 'correct'],
          'gap': {}}
for label, tags, status in [('before', old_tags, old), ('pgphase', new_tags, new)]:
    counts = Counter(status.get(q, 'unphased') for q in names)
    core = Counter(tags[q][1] for q in names if status.get(q) == 'correct' and tags[q][1] < 1_000_000_000)
    result['gap'][label] = {'scorable': len(names), 'counts': dict(counts),
                          'core_correct': max(core.values(), default=0)}
v, w = audit.variants(base/'phased.vcf'), audit.variants(fixed/'native.vcf')
result['variants'] = {'before': len(v), 'after': len(w), 'missing_keys': sorted(v.keys()-w.keys()),
                      'added_keys': sorted(w.keys()-v.keys())}
hom = w[('CHM13#0#chr20', '36620864', '.', 'G', 'A')]
assert hom[1]['GT'] == '1/1' and 'PS' not in hom[1] and hom[1]['AD'] == '11,44'
assert 'CAT=CLEAN_HOM' in hom[0][2]
assert result['gap']['pgphase']['counts']['correct'] >= 91
assert result['gap']['pgphase']['core_correct'] >= 91
assert result['gap']['pgphase']['counts']['correct']/len(names) >= .8
with pysam.AlignmentFile(str(ROOT/'test_data/HG002_chr20_hifi_mapped_to_CHM13_chr20_annotated.bam')) as bam:
    control = {r.query_name for r in bam.fetch('CHM13#0#chr20', 36611592, 36620864)
               if not r.is_secondary and not r.is_supplementary and r.query_name in truth}
control_counts = Counter(new.get(q, 'unphased') for q in control)
control_cores = Counter(new_tags[q][1] for q in control if new.get(q) == 'correct' and new_tags[q][1] < 1_000_000_000)
result['control_gap'] = {'scorable': len(control), 'counts': dict(control_counts),
                         'core_correct': max(control_cores.values(), default=0)}
assert control_counts['correct'] >= 82 and max(control_cores.values()) >= 82

(OUT/'owner-results.json').write_text(json.dumps(result, indent=2)+'\n')
print(json.dumps(result, indent=2))
