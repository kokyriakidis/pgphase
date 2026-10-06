#!/usr/bin/env python3
"""Check the union in its preceding-chunk context without weakening evidence."""
from collections import Counter
import importlib.util
import json
from pathlib import Path
import pysam

ROOT = Path(__file__).resolve().parents[2]
OUT = Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location('audit', OUT/'read_audit.py')
audit = importlib.util.module_from_spec(spec)
spec.loader.exec_module(audit)
truth = {f[0]: f[1] == 'PATERNAL' for line in (ROOT/'test_data/derived/chr20_truth_hap.tsv').open()
         if len(f := line.rstrip().split('\t')) == 2}
before = ROOT/'test_data/tmp_gap_fix83/baseline_continuation/54'
after = ROOT/'test_data/tmp_gap_fix83/final_continuation/54'
x, xs = audit.assignments(before/'phased.bam', truth)
y, ys = audit.assignments(after/'phased.bam', truth)
assert x.keys() == y.keys()
lost = [q for q in xs if xs[q] == 'correct' and ys[q] != 'correct']
assert all(x[q] == y[q] and x[q][1] >= 1_000_000_000 for q in lost), lost
assert all(y[q][0] in (1, 2) and y[q][1] > 0 for q, (hp, ps) in x.items() if hp in (1, 2) and ps > 0)
vx, vy = audit.variants(before/'phased.vcf'), audit.variants(after/'phased.vcf')
assert not (vx.keys()-vy.keys())
assert vy.keys()-vx.keys() == {('CHM13#0#chr20', '55309794', '.', 'T', 'TT'), ('CHM13#0#chr20', '55309794', '.', 'TTT', 'T')}
for key in vx:
    info, sample = vx[key]
    ninfo, nsample = vy[key]
    assert info == ninfo
    assert {k:v for k,v in sample.items() if k not in ('GT','PS')} == {k:v for k,v in nsample.items() if k not in ('GT','PS')}
markers = [(key, row[1]) for key, row in vy.items() if int(key[1]) in (54453020, 55309475, 55309789, 55309794)]
assert len({sample['PS'] for key, sample in markers}) == 1
with pysam.AlignmentFile(str(ROOT/'test_data/HG002_chr20_hifi_mapped_to_CHM13_chr20_annotated.bam')) as bam:
    names = {r.query_name for r in bam.fetch('CHM13#0#chr20', 55309788, 55309794)
             if not r.is_secondary and not r.is_supplementary and r.query_name in truth}
counts = Counter(ys[q] for q in names)
core = Counter(y[q][1] for q in names if ys[q] == 'correct' and y[q][1] < 1_000_000_000)
assert len(names) == 67 and counts['correct'] >= 66 and max(core.values()) >= 66
result = {'primary_reads': len(y), 'before_counts': dict(Counter(xs.values())), 'after_counts': dict(Counter(ys.values())),
          'lost_correct_assignments': lost, 'lost_correct_changed_tags': [], 'changed_evidence': [],
          'gap': {'scorable': len(names), 'counts': dict(counts), 'correct_cores': dict(core)},
          'markers': markers}
(OUT/'continuation-results.json').write_text(json.dumps(result, indent=2)+'\n')
print(json.dumps(result, indent=2))
