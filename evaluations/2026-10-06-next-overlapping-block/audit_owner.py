#!/usr/bin/env python3
"""Audit the native owner and thread determinism, preserving unrelated evidence."""
from collections import Counter
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
before = ROOT/'test_data/tmp_gap_fix88/baseline/19'
after = ROOT/'test_data/tmp_gap_fix88/final_owner/19'
single = ROOT/'test_data/tmp_gap_fix88/final_single/19'
old_tags, old = a.assignments(before/'phased.bam', truth)
new_tags, new = a.assignments(after/'phased.bam', truth)
one_tags, one = a.assignments(single/'phased.bam', truth)
assert old_tags.keys() <= new_tags.keys()
assert new_tags == one_tags
lost = [q for q in old if old[q] == 'correct' and new[q] != 'correct']
assert not lost, lost
assert all(new_tags[q][0] in (1, 2) and new_tags[q][1] > 0
           for q, (hp, ps) in old_tags.items() if hp in (1, 2) and ps > 0)
vold, vnew, vone = a.variants(before/'phased.vcf'), a.variants(after/'phased.vcf'), a.variants(single/'phased.vcf')
assert vnew == vone
assert vold == vnew
assert (before/'candidates.tsv').read_bytes() == (after/'candidates.tsv').read_bytes()
assert (after/'candidates.tsv').read_bytes() == (single/'candidates.tsv').read_bytes()
assert vold.keys() == vnew.keys()
changed_evidence = []
for key, (info, sample) in vold.items():
    if key not in vnew: continue
    ninfo, nsample = vnew[key]
    if info != ninfo or {k: v for k, v in sample.items() if k not in ('GT', 'PS')} != {
            k: v for k, v in nsample.items() if k not in ('GT', 'PS')}:
        changed_evidence.append(key)

assert not changed_evidence
result = {'primary_owner_reads': len(new_tags), 'before_counts': dict(Counter(old.values())),
          'after_counts': dict(Counter(new.values())),
          'transitions': dict(Counter(f'{old[q]}->{new[q]}' for q in old)),
          'lost_correct_assignments': lost, 'lost_correct_changed_tags': [], 'lost_phased_assignments': 0,
          'changed_evidence': changed_evidence, 'added_variant_keys': sorted(vnew.keys()-vold.keys()), 'removed_variant_keys': sorted(vold.keys()-vnew.keys()),
          'changed_read_tags': sum(old_tags[q] != new_tags[q] for q in old_tags),
          'single_vs_four_threads': 'identical candidates, VCF records and BAM assignments'}
(OUT/'owner-results.json').write_text(json.dumps(result, indent=2)+'\n')
print(json.dumps(result, indent=2))
