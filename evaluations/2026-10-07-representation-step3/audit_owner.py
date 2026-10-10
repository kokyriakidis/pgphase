#!/usr/bin/env python3
"""Check that the unchanged 5–6 Mb matrix also preserves output and parental calls."""
from collections import Counter
import importlib.util
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
spec = importlib.util.spec_from_file_location('read_audit',
    ROOT / 'evaluations/2026-10-06-shared-insertion-source/read_audit.py')
audit = importlib.util.module_from_spec(spec)
spec.loader.exec_module(audit)
truth = {f[0]: f[1] == 'PATERNAL'
         for line in (ROOT / 'test_data/derived/chr20_truth_hap.tsv').read_text().splitlines()
         if len(f := line.split('\t')) == 2 and f[1] in ('PATERNAL', 'MATERNAL')}
before = ROOT / 'test_data/tmp_representation_step2/matrix5'
after = ROOT / 'test_data/tmp_representation_step3/matrix5'
bt, bs = audit.assignments(before / 'phased.bam', truth)
at, ast = audit.assignments(after / 'phased.bam', truth)
assert bt == at and bs == ast
assert audit.variants(before / 'phased.vcf') == audit.variants(after / 'phased.vcf')
assert (before / 'candidates.tsv').read_bytes() == (after / 'candidates.tsv').read_bytes()
report = dict(owning_region='CHM13#0#chr20:5000001-6000000', candidate_rows_identical=True,
              vcf_rows_identical=True, primary_tags_identical=True, parental_status_identical=True,
              primary_reads=len(at), parental_counts=dict(Counter(ast.values())))
(Path(__file__).resolve().parent / 'owner-output-checks.json').write_text(
    json.dumps(report, indent=2) + '\n')
print(json.dumps(report, indent=2))
