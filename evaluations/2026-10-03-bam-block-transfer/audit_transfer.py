#!/usr/bin/env python3
"""Audit BAM source observations and imported rows; truth is evaluation only."""
import argparse
import json
from collections import Counter
from pathlib import Path

p = argparse.ArgumentParser(description=__doc__)
p.add_argument('--source', type=Path, required=True)
p.add_argument('--input', type=Path, required=True)
p.add_argument('--final', type=Path, required=True)
p.add_argument('--phase-set', type=int, required=True)
p.add_argument('--left', type=int, required=True)
p.add_argument('--right', type=int, required=True)
p.add_argument('--output', type=Path, required=True)
a = p.parse_args()


def load(path):
    variants, reads, observations = {}, {}, {}
    for line in path.read_text().splitlines():
        f = line.split('\t')
        if f[0] == 'VAR':
            variants[int(f[1])] = dict(pos=int(f[2]), type=f[3], length=int(f[6]),
                                      alt=f[7], ps=int(f[8]),
                                      hap1=int(f[9]), hap2=int(f[10]))
        elif f[0] == 'READ':
            reads[f[1]] = dict(mapq=int(f[4]), skipped=int(f[5]),
                              hap=int(f[6]), ps=int(f[7]))
        elif f[0] == 'OBS':
            observations.setdefault(f[1], {})[int(f[2])] = int(f[3])
    return variants, reads, observations


source, reads, observations = load(a.source)
transferred, _, _ = load(a.input)
final, _, _ = load(a.final)
selected = {i: v for i, v in source.items()
            if v['ps'] == a.phase_set and v['hap1'] >= 0 and v['hap2'] >= 0
            and v['hap1'] != v['hap2']}
key = lambda v: (v['pos'], v['type'], v['length'], v['alt'])
index = {key(v): v for v in transferred.values()}
final_index = {key(v): v for v in final.values()}
rows = []
for v in selected.values():
    t, f = index.get(key(v)), final_index.get(key(v))
    rows.append(dict(source=v, transferred=t, final=f))
left = next(i for i, v in selected.items() if v['pos'] == a.left)
right = next(i for i, v in selected.items() if v['pos'] == a.right)
pairs = []
for name, obs in observations.items():
    if obs.get(left, -1) < 0 or obs.get(right, -1) < 0:
        continue
    x, y = selected[left], selected[right]
    lh = 1 if obs[left] == x['hap1'] else 2 if obs[left] == x['hap2'] else 0
    rh = 1 if obs[right] == y['hap1'] else 2 if obs[right] == y['hap2'] else 0
    pairs.append(dict(read=name, left_allele=obs[left], right_allele=obs[right],
                      source_consistent=lh == rh and lh != 0, **reads[name]))
result = dict(source_phase_set=a.phase_set, source_phased_sites=len(selected),
              source_extent=[min(v['pos'] for v in selected.values()),
                             max(v['pos'] for v in selected.values())],
              rows=rows, boundary_pairs=pairs,
              assigned_pair_votes=dict(Counter(
                  'consistent' if pair['source_consistent'] else 'conflicting'
                  for pair in pairs if pair['ps'] == a.phase_set)))
a.output.write_text(json.dumps(result, indent=2) + '\n')
