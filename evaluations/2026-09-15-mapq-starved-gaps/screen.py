#!/usr/bin/env python3
"""Screen unresolved chr20 gaps for MAPQ starvation.

For every gap the pipeline does not join, compare the coverage actually present
in the surjected BAM against the coverage that survives the min_mapq floor, and
report how many het sites a competitor phases inside the same window. Reads the
frozen inputs only -- no pgphase rerun.
"""
import argparse
import csv
import collections
import subprocess
from pathlib import Path

p = argparse.ArgumentParser(description=__doc__)
p.add_argument('--tiers', type=Path, required=True)
p.add_argument('--bam', type=Path, required=True)
p.add_argument('--candidates', type=Path, required=True)
p.add_argument('--competitor-vcf', type=Path, required=True)
p.add_argument('--contig', default='CHM13#0#chr20')
p.add_argument('--competitor-contig', default='chr20')
p.add_argument('--min-mapq', type=int, default=30)
p.add_argument('--output', type=Path, required=True)
a = p.parse_args()

RANK = {'joined': 0, 'split': 1, 'partial': 2, 'vetoed': 3, 'rejected': 4, 'open': 5}
best, reads_added = {}, {}
with a.tiers.open() as stream:
    for r in csv.DictReader(stream, delimiter='\t'):
        key = (int(r['GAP_LEFT']), int(r['GAP_RIGHT']))
        if key not in best or RANK.get(r['STATUS'], 9) < RANK.get(best[key], 9):
            best[key] = r['STATUS']
        reads_added[key] = max(reads_added.get(key, 0), int(r['READS_ADDED']))

cands = collections.defaultdict(int)
with a.candidates.open() as stream:
    for r in csv.DictReader(stream, delimiter='\t'):
        cands[int(r['POS'])] += 1
positions = sorted(cands)


def mean_depth(region, min_mapq):
    out = subprocess.run(['samtools', 'depth', '-a', '-Q', str(min_mapq), '-r', region, str(a.bam)],
                         capture_output=True, text=True).stdout
    total = n = 0
    for line in out.splitlines():
        total += int(line.rsplit('\t', 1)[1]); n += 1
    return total / n if n else 0.0


def competitor_hets(region):
    out = subprocess.run(['tabix', str(a.competitor_vcf), region], capture_output=True, text=True).stdout
    n = 0
    for line in out.splitlines():
        f = line.split('\t')
        gt = f[9].split(':')[0]
        if gt in ('0|1', '1|0'): n += 1
    return n


import bisect
rows = []
for (left, right), status in sorted(best.items()):
    if status == 'joined': continue
    span = right - left
    region = f'{a.contig}:{left}-{right}'
    d_all, d_pass = mean_depth(region, 0), mean_depth(region, a.min_mapq)
    lo, hi = bisect.bisect_left(positions, left), bisect.bisect_right(positions, right)
    rows.append({
        'gap_left': left, 'gap_right': right, 'gap_bp': span, 'status': status,
        'reads_added': reads_added[(left, right)],
        'depth_all': round(d_all, 1), 'depth_pass_mapq': round(d_pass, 1),
        'depth_lost_pct': round(100.0 * (1 - d_pass / d_all), 1) if d_all else 0.0,
        'our_candidates': hi - lo,
        'our_cand_per_kb': round(1000.0 * (hi - lo) / span, 2) if span else 0.0,
        'competitor_phased_hets': competitor_hets(f'{a.competitor_contig}:{left}-{right}'),
    })

with a.output.open('w', newline='') as stream:
    w = csv.DictWriter(stream, fieldnames=list(rows[0]), delimiter='\t')
    w.writeheader()
    w.writerows(rows)
print(f'wrote {len(rows)} unresolved gaps to {a.output}')
