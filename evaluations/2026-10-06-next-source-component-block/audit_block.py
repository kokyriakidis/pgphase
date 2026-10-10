#!/usr/bin/env python3
"""Audit both seams and the target block on the same original primary alignments."""
from collections import Counter
import hashlib
import importlib.util
import json
from pathlib import Path
import pysam

ROOT = Path(__file__).resolve().parents[2]
OUT = Path(__file__).resolve().parent
BEFORE = ROOT/'test_data/tmp_gap_fix88/final/0'
AFTER = ROOT/'test_data/tmp_gap_fix89/final/0'
HI = ROOT/'test_data/tmp_gap_fix48/competitor/hiphase_dv'
spec = importlib.util.spec_from_file_location('audit', OUT/'read_audit.py')
a = importlib.util.module_from_spec(spec)
spec.loader.exec_module(a)
truth = {f[0]: f[1] == 'PATERNAL' for line in (ROOT/'test_data/derived/chr20_truth_hap.tsv').open()
         if len(f := line.rstrip().split('\t')) == 2}
regions = {'gap1': (46389864, 46402979), 'gap2': (46402979, 46405048),
           'block': (45866904, 46636707)}
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
    chrom = "CHM13#0#chr20" if "CHM13#0#chr20" in bam.references else "chr20"
    for read in bam.fetch(chrom, 45866903, 46636707):
        if read.is_secondary or read.is_supplementary or read.query_name not in geometry: continue
        assert original_geometry(read) == geometry[read.query_name], read.query_name
        matched.add(read.query_name)
assert matched == geometry.keys()
report = {'target': {'left': 45866904, 'right': 46636707, 'length': 769804, 'hiphase_rank': 30},
          'production_sha256': hashlib.sha256((ROOT/'pgphase').read_bytes()).hexdigest(), 'tools': {}}
for tool, path in [('before', BEFORE), ('pgphase', AFTER), ('hiphase', HI)]:
    tags, status = (json.loads((ROOT/'test_data/tmp_gap_fix78/hi-tags.json').read_text())
                    if tool == 'hiphase' else a.assignments(path/'phased.bam', truth))
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
covering = [b for b in new_extents if b[0] <= 45866904 and b[1] >= 46636707]
assert covering
covering_gap = [b for b in new_extents if b[0] <= 46389864 and b[1] >= 46402979]
assert covering_gap
lengths = sorted((hi-lo+1 for lo, hi in new_extents), reverse=True)
cumulative = 0
for length in lengths:
    cumulative += length
    if cumulative >= sum(lengths)/2:
        n50 = length; break
report['blocks'] = {'count': len(lengths), 'n50_bp': n50, 'largest_block_bp': max(lengths),
                    'covering_target': covering, 'covering_gap': covering_gap, 'previous_extents_preserved': True}
for name in ('gap1', 'gap2'):
    pg, hi = report['tools']['pgphase'][name], report['tools']['hiphase'][name]
    assert pg['correct_fraction'] >= .8
    assert pg['counts']['correct'] >= hi['counts']['correct']
    assert pg['dominant_core_correct'] >= hi['dominant_core_correct']
assert report['tools']['pgphase']['gap1']['scorable'] == 118
assert report['tools']['pgphase']['gap1']['counts']['correct'] >= 107
assert report['tools']['pgphase']['gap1']['dominant_core_correct'] >= 107
assert report['tools']['pgphase']['gap2']['scorable'] == 75
assert report['tools']['pgphase']['gap2']['counts']['correct'] >= 72
assert report['tools']['pgphase']['gap2']['dominant_core_correct'] >= 72
report['whole_block_deficits'] = {field: max(0, report['tools']['hiphase']['block'][field]-report['tools']['pgphase']['block'][field])
    for field in ('dominant_core_correct',)}
report['whole_block_deficits']['correct'] = max(0, report['tools']['hiphase']['block']['counts'].get('correct', 0)-report['tools']['pgphase']['block']['counts'].get('correct', 0))
(OUT/'results.json').write_text(json.dumps(report, indent=2)+'\n')
print(json.dumps(report, indent=2))
