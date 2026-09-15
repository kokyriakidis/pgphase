#!/usr/bin/env python3
"""Attribute recovery joins using frozen original blocks and parental read truth."""
import argparse
import collections
import csv
import gzip
import json
from pathlib import Path

p = argparse.ArgumentParser(description=__doc__)
p.add_argument('--baseline', type=Path, required=True)
p.add_argument('--run', type=Path, required=True)
p.add_argument('--audit', type=Path, required=True)
p.add_argument('--output', type=Path, required=True)
a = p.parse_args()

def reads(root):
    with gzip.open(root / 'read_eval/per_read.tsv.gz', 'rt') as f:
        return {r['read_name']: r for r in csv.DictReader(f, delimiter='\t')}

before, after = reads(a.baseline), reads(a.run)
groups = collections.defaultdict(list)
transitions = collections.Counter()
for name, old in before.items():
    new = after.get(name)
    if new is None:
        transitions['not_evaluated_after'] += 1
        continue
    transitions[old['status'] + '->' + new['status']] += 1
    if old['status'] not in ('concordant', 'DISCORDANT') or new['status'] not in ('concordant', 'DISCORDANT'):
        continue
    groups[int(old['PS'])].append((old, new))
block_map = {}
block_rows = []
for ps, pairs in sorted(groups.items()):
    mapping = {(int(new['PS']), old['HP'] != new['HP']) for old, new in pairs}
    if len(mapping) != 1:
        raise ValueError(f'Original block {ps} has nonuniform transformation: {mapping}')
    destination, flip = next(iter(mapping))
    old_errors = sum(old['status'] == 'DISCORDANT' for old, _ in pairs)
    new_errors = sum(new['status'] == 'DISCORDANT' for _, new in pairs)
    block_map[ps] = destination, flip, pairs[0][0]['orientation'], old_errors / len(pairs)
    block_rows.append([ps, destination, len(pairs), old_errors, new_errors, flip])
edge_rows = []
for path in sorted(a.audit.glob('*.evidence.tsv')):
    with path.open() as f:
        assert f.readline().rstrip() == 'SCHEMA\t1'
        _, tid, left, right, beg, end, *_ = f.readline().rstrip().split('\t')
    left, right = int(left), int(right)
    status = 'UNASSESSED'
    expected = observed = ''
    if left in block_map and right in block_map:
        l, r = block_map[left], block_map[right]
        if l[0] != r[0]:
            status = 'SPLIT'
        elif max(l[3], r[3]) > 0.1:
            status = 'UNCERTAIN_TRUTH'
        else:
            expected = l[2] != r[2]
            observed = l[1] != r[1]
            status = 'CORRECT' if expected == observed else 'WRONG'
    edge_rows.append([tid, beg, end, left, right, status, expected, observed])
a.output.mkdir(parents=True, exist_ok=True)
for name, header, rows in [
    ('blocks.tsv', ['original_ps', 'final_ps', 'common_reads', 'before_errors', 'after_errors', 'label_flip'], block_rows),
    ('edges.tsv', ['tid', 'left_end', 'right_beg', 'left_ps', 'right_ps', 'status', 'expected_flip', 'observed_flip'], edge_rows),
]:
    with (a.output / name).open('w') as f:
        w = csv.writer(f, delimiter='\t', lineterminator='\n')
        w.writerow(header)
        w.writerows(rows)
summary = dict(edges=dict(collections.Counter(row[5] for row in edge_rows)),
               common_read_transitions=dict(transitions),
               original_blocks=len(block_rows),
               uniform_original_block_transformations=True)
(a.output / 'summary.json').write_text(json.dumps(summary, indent=2) + '\n')
print(json.dumps(summary, indent=2))
