#!/usr/bin/env python3
"""Locate chromosome label losses and audit the other withdrawn singleton."""
from collections import Counter
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
before, old = audit.assignments(ROOT/'test_data/tmp_gap_fix74/frozen_final/0/phased.bam', truth)
after, new = audit.assignments(ROOT/'test_data/tmp_gap_fix75/retired_final/0/phased.bam', truth)
lost = {q for q in old if old[q] == 'correct' and new[q] != 'correct'}
rows, names, physical = [], set(), Counter()
pos = 45859664
with pysam.AlignmentFile(str(ROOT/'test_data/HG002_chr20_hifi_mapped_to_CHM13_chr20_annotated.bam')) as bam:
    for left, right in [(36600000,36630000),(45835000,45875000)]:
        for read in bam.fetch('CHM13#0#chr20',left-1,right):
            if read.is_secondary or read.is_supplementary:
                continue
            q = read.query_name
            if q in lost:
                rows.append({'read':q,'start':read.reference_start+1,'end':read.reference_end,
                             'before':before[q],'after':after[q]})
            if left == 45835000 and q in truth:
                names.add(q)
    for read in bam.fetch('CHM13#0#chr20',pos-1,pos):
        if read.flag & (4 | 256 | 2048 | 512 | 1024) or read.mapping_quality < 30 or read.mapping_quality == 255:
            continue
        for qi, rp in read.get_aligned_pairs():
            if rp != pos-1:
                continue
            if qi is None:
                physical['DEL'] += 1
            elif read.query_qualities[qi] >= 20:
                physical[read.query_sequence[qi]] += 1
            break
assert {r['read'] for r in rows} == lost
assert all(r['before'][1] in (36620864,45859664) for r in rows)
assert set(physical) == {'G'}
result = {'lost_correct_rows':rows,'lost_by_old_phase_set':dict(Counter(r['before'][1] for r in rows)),
          'other_withdrawn_singleton':{'position':pos,'ref':'A','alt':'G','physical_calls':dict(physical),
               'original_region':[45835000,45875000],'scorable':len(names),
               'before_counts':dict(Counter(old.get(q,'unphased') for q in names)),
               'after_counts':dict(Counter(new.get(q,'unphased') for q in names))}}
(OUT/'side-effects.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps({k:v for k,v in result.items() if k!='lost_correct_rows'},indent=2))
