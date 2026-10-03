#!/usr/bin/env python3
"""Evaluate old SNP gauges and coordinate coverage; never used by phasing."""
import argparse
from collections import Counter, defaultdict
import csv
import json
from pathlib import Path

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--before', type=Path, required=True)
parser.add_argument('--after', type=Path, required=True)
parser.add_argument('--panel', type=Path, required=True)
parser.add_argument('--output', type=Path, required=True)
args = parser.parse_args()


def load(path):
    rows = {}
    blocks = defaultdict(list)
    for line in path.open():
        if line.startswith('#'):
            continue
        fields = line.rstrip().split('\t')
        values = dict(zip(fields[8].split(':'), fields[9].split(':')))
        gt, ps = values['GT'], values.get('PS', '.')
        if '|' not in gt or len(set(gt.split('|'))) != 2 or ps in ('.', '0', '-1'):
            continue
        key = (int(fields[1]), fields[3], fields[4])
        rows[key] = (ps, gt)
        blocks[ps].append(key[0])
    return rows, {ps: (min(p), max(p)) for ps, p in blocks.items()}


before, before_blocks = load(args.before)
after, after_blocks = load(args.after)
gauges = defaultdict(Counter)
lost = []
for key, (ps, gt) in before.items():
    if len(key[1]) != 1 or len(key[2]) != 1:
        continue
    if key not in after:
        lost.append(key)
        continue
    new_ps, new_gt = after[key]
    gauges[ps][(new_ps, gt != new_gt)] += 1
mixed = {ps: [{'new_ps': new_ps, 'flipped': flip, 'snps': count}
              for (new_ps, flip), count in counter.items()]
         for ps, counter in gauges.items() if len(counter) > 1}
all_block_gauges = defaultdict(Counter)
for key, (old_ps, gt) in before.items():
    if key in after:
        new_ps, new_gt = after[key]
        all_block_gauges[old_ps][(new_ps, gt != new_gt)] += 1
joined = defaultdict(dict)
for old_ps, counter in all_block_gauges.items():
    for (new_ps, flip), count in counter.items():
        joined[new_ps][old_ps] = {'flipped': flip, 'variants': count,
                                 'extent': before_blocks[old_ps]}
joined = {ps: members for ps, members in joined.items() if len(members) > 1}


def spans(blocks, left, right):
    return any(beg <= left and end >= right for beg, end in blocks.values())


panel = list(csv.DictReader((line for line in args.panel.open()
                            if not line.startswith('#')), delimiter='\t'))
changes = []
open_windows = []
before_count = after_count = 0
for row in panel:
    left, right = int(row['gap_left']), int(row['gap_right'])
    old, new = spans(before_blocks, left, right), spans(after_blocks, left, right)
    before_count += old
    after_count += new
    if old != new:
        changes.append({'left': left, 'right': right, 'before': old, 'after': new})
    if not new:
        open_windows.append({'left': left, 'right': right,
                             'competitor': row['competitor']})
result = dict(lost_phased_snps=lost, old_snp_blocks_with_mixed_gauge=mixed,
              new_blocks_containing_multiple_old_vcf_blocks=joined,
              panel_spanned_before=before_count, panel_spanned_after=after_count,
              changed_panel_spans=changes, open_panel_windows=open_windows)
args.output.write_text(json.dumps(result, indent=2) + '\n')
print(json.dumps(result, indent=2))
