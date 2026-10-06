#!/usr/bin/env python3
"""Verify original-label preservation and identical native one/four-thread results."""
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
base = ROOT/'test_data/tmp_gap_fix76/baseline/36'
fixed = ROOT/'test_data/tmp_gap_fix76/stable/36'
single = ROOT/'test_data/tmp_gap_fix76/stable_single/36'
old_tags, old = audit.assignments(base/'phased.bam', truth)
new_tags, new = audit.assignments(fixed/'phased.bam', truth)
assert (fixed/'candidates.tsv').read_bytes() == (single/'candidates.tsv').read_bytes()
assert (fixed/'phased.vcf').read_bytes() == (single/'phased.vcf').read_bytes()
assert (new_tags, new) == audit.assignments(single/'phased.bam', truth)
result = audit.audit(base, fixed, truth)
result.update({'binary_sha256': hashlib.sha256((ROOT/'pgphase').read_bytes()).hexdigest(),
               'one_four_thread_identity': True,
               'transitions': dict(Counter(f'{old[q]}->{new[q]}' for q in old))})
with pysam.AlignmentFile(str(ROOT/'test_data/HG002_chr20_hifi_mapped_to_CHM13_chr20_annotated.bam')) as bam:
    names = {r.query_name for r in bam.fetch('CHM13#0#chr20', 36268557, 36286778)
             if not r.is_secondary and not r.is_supplementary and r.query_name in truth}
counts = Counter(new.get(q, 'unphased') for q in names)
cores = Counter(new_tags[q][1] for q in names if new.get(q) == 'correct' and new_tags[q][1] < 1_000_000_000)
result['gap'] = {'scorable': len(names), 'counts': dict(counts), 'core_correct': max(cores.values())}
assert len(names) == 109 and counts['correct'] >= 100 and max(cores.values()) >= 100
assert counts['correct']/len(names) >= .8
(OUT/'owner-results.json').write_text(json.dumps(result, indent=2)+'\n')
print(json.dumps(result, indent=2))
