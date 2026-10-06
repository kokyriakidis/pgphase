#!/usr/bin/env python3
"""Audit both seams and the target block on the same original primary alignments."""
from collections import Counter, defaultdict
import hashlib
import importlib.util
import json
from pathlib import Path
import pysam

ROOT = Path(__file__).resolve().parents[2]
OUT = Path(__file__).resolve().parent
BEFORE = ROOT/'test_data/tmp_gap_fix84/final3/0'
AFTER = ROOT/'test_data/tmp_gap_fix85/final15/0'
HI = ROOT/'test_data/tmp_gap_fix48/competitor/hiphase_dv'
spec = importlib.util.spec_from_file_location('audit', OUT/'read_audit.py')
a = importlib.util.module_from_spec(spec)
spec.loader.exec_module(a)
truth = {f[0]: f[1] == 'PATERNAL' for line in (ROOT/'test_data/derived/chr20_truth_hap.tsv').open()
         if len(f := line.rstrip().split('\t')) == 2}
regions = {'gap1': (15023123, 15039543), 'gap2': (15100456, 15101262),
           'control': (15095642, 15101261), 'block': (14719378, 15335938)}
groups, geometry = {}, {}
def original_geometry(read):
    return (read.reference_start, read.reference_end, read.cigarstring, read.flag, read.mapping_quality,
            hashlib.sha256(read.query_sequence.encode()).hexdigest(), hashlib.sha256(bytes(read.query_qualities)).hexdigest())
with pysam.AlignmentFile(str(ROOT/'test_data/HG002_chr20_hifi_mapped_to_CHM13_chr20_annotated.bam')) as bam:
    for name, (left, right) in regions.items():
        groups[name] = set()
        for read in bam.fetch('CHM13#0#chr20', left-1, right):
            if read.is_secondary or read.is_supplementary or read.query_name not in truth: continue
            groups[name].add(read.query_name)
            geometry[read.query_name] = original_geometry(read)
with pysam.AlignmentFile(str(HI/'phased.bam')) as bam:
    matched = set()
    for read in bam:
        if read.is_secondary or read.is_supplementary or read.query_name not in geometry: continue
        assert original_geometry(read) == geometry[read.query_name], read.query_name
        matched.add(read.query_name)
assert matched == geometry.keys()
report = {'target': {'left': 14719378, 'right': 15335938, 'length': 616561, 'hiphase_rank': 36},
          'production_sha256': hashlib.sha256((ROOT/'pgphase').read_bytes()).hexdigest(), 'tools': {}}
for tool, path in [('before', BEFORE), ('pgphase', AFTER), ('hiphase', HI)]:
    tags, status = a.assignments(path/'phased.bam', truth)
    measurements = {}
    for name, names in groups.items():
        counts = Counter(status.get(q, 'unphased') for q in names)
        cores = Counter(tags[q][1] for q in names if status.get(q) == 'correct' and 0 < tags[q][1] < 1000000000)
        measurements[name] = {'scorable': len(names), 'counts': dict(counts), 'correct_fraction': counts['correct']/len(names),
                              'dominant_core_correct': max(cores.values(), default=0), 'core_counts': dict(cores)}
    report['tools'][tool] = measurements
old_rows, new_rows = a.variants(BEFORE/'phased.vcf'), a.variants(AFTER/'phased.vcf')
old_extents, new_extents = a.block_extents(old_rows), a.block_extents(new_rows)
assert all(any(lo <= left and hi >= right for lo, hi in new_extents) for left, right in old_extents)
covering = [b for b in new_extents if b[0] <= 14719378 and b[1] >= 15335938]
assert covering
lengths = sorted((hi-lo+1 for lo, hi in new_extents), reverse=True)
cumulative = 0
for length in lengths:
    cumulative += length
    if cumulative >= sum(lengths)/2:
        n50 = length; break
report['blocks'] = {'count': len(lengths), 'n50_bp': n50, 'largest_block_bp': max(lengths),
                    'covering_target': covering, 'previous_extents_preserved': True}
for name in ('gap1', 'gap2', 'control'):
    pg, hi = report['tools']['pgphase'][name], report['tools']['hiphase'][name]
    assert pg['correct_fraction'] >= .8
    assert pg['counts']['correct'] >= hi['counts']['correct']
    assert pg['dominant_core_correct'] >= hi['dominant_core_correct']
assert report['tools']['pgphase']['gap1']['scorable'] == 126
assert report['tools']['pgphase']['gap2']['scorable'] == 73
report['whole_block_deficits'] = {field: max(0, report['tools']['hiphase']['block'][field]-report['tools']['pgphase']['block'][field])
    for field in ('dominant_core_correct',)}
report['whole_block_deficits']['correct'] = max(0, report['tools']['hiphase']['block']['counts'].get('correct', 0)-report['tools']['pgphase']['block']['counts'].get('correct', 0))
(OUT/'results.json').write_text(json.dumps(report, indent=2)+'\n')
print(json.dumps(report, indent=2))
