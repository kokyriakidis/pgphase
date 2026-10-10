#!/usr/bin/env python3
"""Preserve all native owning read labels and candidate evidence across repair."""
from collections import Counter
import hashlib
import importlib.util
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
OUT = Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location('audit', OUT/'read_audit.py')
a = importlib.util.module_from_spec(spec)
spec.loader.exec_module(a)
truth = {f[0]: f[1] == 'PATERNAL' for line in (ROOT/'test_data/derived/chr20_truth_hap.tsv').open()
         if len(f := line.rstrip().split('\t')) == 2}
before = ROOT/'test_data/tmp_gap_fix90/baseline/4'
after = ROOT/'test_data/tmp_gap_fix90/final_owner/4'
single = ROOT/'test_data/tmp_gap_fix90/final_single/4'
x, xs = a.assignments(before/'phased.bam', truth)
y, ys = a.assignments(after/'phased.bam', truth)
z, zs = a.assignments(single/'phased.bam', truth)
assert x.keys() == y.keys() == z.keys()
assert all(y[q] == x[q] for q in x if x[q][0] in (1, 2) and x[q][1] > 0)
assert all(ys[q] == 'correct' for q in xs if xs[q] == 'correct')
assert y == z and ys == zs
assert before.joinpath('candidates.tsv').read_bytes() == after.joinpath('candidates.tsv').read_bytes() == single.joinpath('candidates.tsv').read_bytes()
assert a.variants(before/'phased.vcf') == a.variants(after/'phased.vcf') == a.variants(single/'phased.vcf')
changed = sorted(q for q in x if x[q] != y[q])
assert len(changed) == 18
assert all(xs[q] == 'unphased' and ys[q] == 'correct' for q in changed)
result = {'binary_sha256': hashlib.sha256((ROOT/'pgphase').read_bytes()).hexdigest(),
          'before_counts': dict(Counter(xs.values())), 'after_counts': dict(Counter(ys.values())),
          'changed_tags': len(changed), 'changed_reads': changed,
          'previous_phased_tags_identical': True, 'lost_correct_assignments': [],
          'candidate_bytes_identical': True, 'complete_vcf_records_identical': True,
          'one_four_thread_tags_candidates_vcf_identical': True}
(OUT/'owner-results.json').write_text(json.dumps(result, indent=2)+'\n')
print({k:v for k,v in result.items() if k != 'changed_reads'})
