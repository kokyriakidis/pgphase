#!/usr/bin/env python3
"""Audit read transitions and unrelated variant evidence after the repeat-SNP repair."""
from collections import Counter
import importlib.util
import json
from pathlib import Path
import pysam
ROOT = Path(__file__).resolve().parents[2]
OUT = Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location('audit', ROOT/'evaluations/2026-10-05-largest-hiphase-block/audit_reads.py')
a = importlib.util.module_from_spec(spec)
spec.loader.exec_module(a)
truth = {f[0]: f[1] == 'PATERNAL' for line in (ROOT/'test_data/derived/chr20_truth_hap.tsv').open()
         if len(f := line.rstrip().split('\t')) == 2}
before = ROOT/'test_data/tmp_gap_fix77/final/0'
after = ROOT/'test_data/tmp_gap_fix78/final/0'
x, xs = a.assignments(before/'phased.bam', truth)
y, ys = a.assignments(after/'phased.bam', truth)
assert x.keys() == y.keys()
lost_correct = [q for q in xs if xs[q] == 'correct' and ys[q] != 'correct']
assert lost_correct == ['m84031_231217_062403_s3/74975318/ccs'], lost_correct
with pysam.AlignmentFile(str(ROOT/'test_data/HG002_chr20_hifi_mapped_to_CHM13_chr20_annotated.bam')) as bam:
    gap_reads = {read.query_name for read in bam.fetch('CHM13#0#chr20', 11796978, 11813446)
                 if not read.is_secondary and not read.is_supplementary}
changed_status = {q for q in xs if xs[q] != ys[q]}
assert changed_status <= gap_reads, changed_status-gap_reads
assert all(y[q][0] in (1, 2) and y[q][1] > 0 for q, (hp, ps) in x.items() if hp in (1, 2) and ps > 0)
vx, vy = a.variants(before/'phased.vcf'), a.variants(after/'phased.vcf')
changed_evidence = []
for key in vx.keys() & vy.keys():
    info, sample = vx[key]
    ninfo, nsample = vy[key]
    if info != ninfo or {k: v for k, v in sample.items() if k not in ('GT', 'PS')} != {
            k: v for k, v in nsample.items() if k not in ('GT', 'PS')}:
        changed_evidence.append(key)
assert not changed_evidence, changed_evidence
removed = sorted(vx.keys()-vy.keys())
added = sorted(vy.keys()-vx.keys())
expected_removed = {11813622: ('T', 'C'), 2282198: ('T', 'A'), 42290009: ('A', 'G'),
                    47003854: ('A', 'G'), 50431256: ('A', 'T'), 62940433: ('A', 'G'),
                    63734075: ('C', 'T'), 63793665: ('C', 'T')}
assert removed == sorted(('CHM13#0#chr20', str(pos), '.', ref, alt)
                         for pos, (ref, alt) in expected_removed.items()), removed
assert added == [('CHM13#0#chr20', '11813622', '.', 'T', 'TAC'),
                 ('CHM13#0#chr20', '11813668', '.', 'T', 'C'),
                 ('CHM13#0#chr20', '47003826', '.', 'G', 'GAGAA'),
                 ('CHM13#0#chr20', '50880007', '.', 'T', 'TA')], added
r = {'primary_reads': len(y), 'before_counts': dict(Counter(xs.values())), 'after_counts': dict(Counter(ys.values())),
     'lost_correct_assignments': lost_correct, 'lost_phased_assignments': 0,
     'changed_read_tags': sum(x[q] != y[q] for q in x), 'changed_evidence': changed_evidence,
     'read_status_changes_confined_to_gap': True,
     'removed_variant_keys': removed, 'added_variant_keys': added,
     'transitions': dict(Counter(f'{xs[q]}->{ys[q]}' for q in xs)),
     'newly_correct': [q for q in xs if xs[q] != 'correct' and ys[q] == 'correct']}
(OUT/'full-preservation.json').write_text(json.dumps(r, indent=2)+'\n')
print(json.dumps(r, indent=2))
