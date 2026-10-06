#!/usr/bin/env python3
"""Measure block spans and correctness on the largest HiPhase interval."""
from collections import Counter, defaultdict
import json
from pathlib import Path
import pysam
from audit_reads import assignments

root = Path('evaluations/2026-10-05-largest-hiphase-block')
truth = {f[0]: f[1] == 'PATERNAL' for line in open('test_data/derived/chr20_truth_hap.tsv')
         if len(f := line.rstrip().split('\t')) == 2}
left, right = 63182011, 66206480
eligible = {}
with pysam.AlignmentFile('test_data/HG002_chr20_hifi_mapped_to_CHM13_chr20_annotated.bam') as bam:
    for read in bam.fetch('CHM13#0#chr20', left - 1, right):
        if not read.is_secondary and not read.is_supplementary and read.query_name in truth:
            eligible[read.query_name] = read.reference_start, read.reference_end, read.cigarstring, read.query_sequence


def spans(path):
    blocks = defaultdict(list)
    with pysam.VariantFile(str(path)) as vcf:
        for row in vcf:
            sample = next(iter(row.samples.values()))
            gt = sample.get('GT')
            ps = sample.get('PS')
            if sample.phased and gt and None not in gt and len(set(gt)) == 2 and ps and ps > 0:
                blocks[ps].append(row.pos)
    lengths = sorted((max(pos) - min(pos) + 1 for pos in blocks.values() if len(pos) > 1), reverse=True)
    accumulated = 0
    n50 = 0
    for length in lengths:
        accumulated += length
        if accumulated * 2 >= sum(lengths):
            n50 = length
            break
    return {'n50_bp': n50, 'largest_block_bp': max(lengths), 'blocks_with_two_hets': len(lengths),
            'overlapping_blocks': [{'ps': ps, 'left': min(pos), 'right': max(pos), 'rows': len(pos)}
                                   for ps, pos in blocks.items() if max(pos) >= left and min(pos) <= right]}

report = {}
for name, directory in [('before', 'test_data/tmp_gap_fix68/frozen_final/0'),
                        ('after', 'test_data/tmp_gap_fix69/frozen_final/0')]:
    directory = Path(directory)
    tags, status = assignments(directory / 'phased.bam', truth)
    counts = Counter(status.get(q, 'unphased') for q in eligible)
    core = Counter(tags[q][1] for q in eligible if status.get(q) == 'correct' and tags[q][1] < 1_000_000_000)
    report[name] = {**spans(directory / 'phased.vcf'), 'scorable': len(eligible), 'counts': dict(counts),
                    'correct_all_reads': counts['correct'] / len(eligible), 'dominant_core_correct': max(core.values())}

path = 'test_data/tmp_gap_fix48/competitor/hiphase_dv/phased.bam'
tags, status = assignments(Path(path), truth)
matched = 0
with pysam.AlignmentFile(path) as bam:
    contig = next(c for c in bam.references if c.endswith('chr20'))
    for read in bam.fetch(contig, left - 1, right):
        if read.is_secondary or read.is_supplementary or read.query_name not in eligible:
            continue
        assert eligible[read.query_name] == (read.reference_start, read.reference_end,
                                             read.cigarstring, read.query_sequence)
        matched += 1
assert matched == len(eligible)
counts = Counter(status.get(q, 'unphased') for q in eligible)
core = Counter(tags[q][1] for q in eligible if status.get(q) == 'correct')
report['hiphase'] = {**spans(Path(path).with_name('phased.vcf.gz')), 'scorable': len(eligible),
                     'counts': dict(counts), 'correct_all_reads': counts['correct'] / len(eligible),
                     'dominant_core_correct': max(core.values()), 'identical_input_alignments': matched}
assert report['after']['correct_all_reads'] >= .8
assert report['after']['counts']['correct'] >= report['hiphase']['counts']['correct']
assert report['after']['dominant_core_correct'] >= report['hiphase']['dominant_core_correct']
(root / 'block-metrics.json').write_text(json.dumps(report, indent=2) + '\n')
print(json.dumps(report, indent=2))
