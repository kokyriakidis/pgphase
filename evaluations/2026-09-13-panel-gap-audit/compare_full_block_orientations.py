#!/usr/bin/env python3
"""Evaluate original-block orientation separately from added-read discordance."""
import argparse
import csv
import gzip
from collections import defaultdict
from pathlib import Path

p = argparse.ArgumentParser(description=__doc__)
p.add_argument('--chromosome', default='chr20')
p.add_argument('--arms', nargs='+', default=['recovery2', 'graph_support2'])
a = p.parse_args()
root = Path('/tmp/pgphase-panel-gap-audit') / a.chromosome
report = Path(__file__).resolve().parent

def reads(arm):
    with gzip.open(root / arm / 'read_eval/per_read.tsv.gz', 'rt') as f:
        return {r['read_name']: r for r in csv.DictReader(f, delimiter='\t')}

before = reads('clean')
for arm in a.arms:
    after = reads(arm)
    groups = defaultdict(list)
    for name, old in before.items():
        if name not in after or old['status'] == 'skipped' or after[name]['status'] == 'skipped':
            continue
        new = after[name]
        groups[old['PS'], new['PS']].append((old, new))
    destinations = defaultdict(set)
    for old_ps, new_ps in groups:
        destinations[old_ps].add(new_ps)
    assert all(len(v) == 1 for v in destinations.values()), (arm, 'split original block')
    rows = []
    for (old_ps, new_ps), observations in groups.items():
        flips = {old['HP'] != new['HP'] for old, new in observations}
        assert len(flips) == 1, (arm, old_ps, 'nonuniform original block')
        n = len(observations)
        old_errors = sum(old['status'] == 'DISCORDANT' for old, _ in observations)
        new_errors = sum(new['status'] == 'DISCORDANT' for _, new in observations)
        rows.append(dict(original_ps=old_ps, final_ps=new_ps, truth_reads=n,
                         before_discordant=old_errors, after_discordant=new_errors,
                         uniform_label_flip=next(iter(flips)),
                         majority_reversed=old_errors * 2 < n and new_errors * 2 > n))
    rows.sort(key=lambda r: int(r['original_ps']))
    with (report / f'{a.chromosome}.{arm}.block_orientations.tsv').open('w') as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]), delimiter='\t', lineterminator='\n')
        w.writeheader(); w.writerows(rows)
    print(arm, 'original blocks with reversed truth majority:',
          sum(r['majority_reversed'] for r in rows))
