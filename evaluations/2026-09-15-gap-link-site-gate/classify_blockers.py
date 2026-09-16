#!/usr/bin/env python3
"""Separate gate-limited gaps from read-coverage-limited ones.

A gap can only be joined if some chain of called sites runs from one flank to the
other with reads spanning every step. Where that chain is complete, admission
policy decides the outcome and a gate fix can help; where a step has no spanning
read, no site-admission policy can help, because the linking reads do not exist.
Tests the chain over the sites each run actually called, straight from the
alignment, so it needs neither truth nor the audit export.
"""
import argparse
import csv
from pathlib import Path

import pysam

p = argparse.ArgumentParser(description=__doc__)
p.add_argument('--panel-root', type=Path, required=True)
p.add_argument('--gaps', type=Path, required=True)
p.add_argument('--bam', type=Path, required=True)
p.add_argument('--contig', default='CHM13#0#chr20')
p.add_argument('--arm', default='union_recover')
p.add_argument('--flank', type=int, default=25000)
p.add_argument('--min-spanning', type=int, default=2)
p.add_argument('--output', type=Path, required=True)
a = p.parse_args()

bam = pysam.AlignmentFile(str(a.bam), 'rb')
gaps = {int(r['gap_left']): r for r in csv.DictReader(a.gaps.open(), delimiter='\t')}
LINKING = {'CLEAN_HET_SNP', 'CLEAN_HET_INDEL', 'NOISY_CAND_HET'}

rows = []
for gap_left, gap in sorted(gaps.items()):
    cand = a.panel_root / f'gap_{gap_left}' / a.arm / 'candidates.tsv'
    if not cand.exists():
        continue
    gap_right = int(gap['gap_right'])
    sites = sorted({int(r['POS']) for r in csv.DictReader(cand.open(), delimiter='\t')
                    if r['CATEGORY'] in LINKING and gap_left <= int(r['POS']) <= gap_right})
    chain = sorted(set(sites) | {gap_left, gap_right})
    spans = [(r.reference_start, r.reference_end) for r in
             bam.fetch(a.contig, max(0, gap_left - a.flank), gap_right + a.flank)
             if not r.is_unmapped and r.reference_end is not None]
    steps = []
    for x, y in zip(chain, chain[1:]):
        n = sum(1 for s, e in spans if s <= x and e >= y)
        steps.append((x, y, y - x, n))
    breaks = [s for s in steps if s[3] < a.min_spanning]
    rows.append({
        'gap_left': gap_left, 'gap_right': gap_right, 'gap_bp': int(gap['gap_bp']),
        'linking_sites': len(sites), 'steps': len(steps), 'unbridged_steps': len(breaks),
        'unbridged_bp': sum(s[2] for s in breaks),
        'widest_unbridged_bp': max((s[2] for s in breaks), default=0),
        # 'chain_present' says only that a positional chain of called sites exists
        # with reads spanning every step -- a necessary condition for any join, not
        # evidence that a gate fix will deliver one, since those sites still have to
        # be informative. 'read_coverage' is the harder verdict: a step no read
        # crosses cannot be bridged by any site-admission policy.
        'blocker': 'read_coverage' if breaks else 'chain_present'})

with a.output.open('w', newline='') as stream:
    w = csv.DictWriter(stream, fieldnames=list(rows[0]), delimiter='\t', lineterminator='\n')
    w.writeheader()
    w.writerows(rows)

print('%-11s %9s %7s %6s %10s %12s %14s %s' % (
    'gap_left', 'gap_bp', 'sites', 'steps', 'unbridged', 'unbridged bp', 'widest break', 'blocker'))
for r in rows:
    print('%-11d %9d %7d %6d %10d %12d %14d %s' % (
        r['gap_left'], r['gap_bp'], r['linking_sites'], r['steps'], r['unbridged_steps'],
        r['unbridged_bp'], r['widest_unbridged_bp'], r['blocker']))
pol = [r for r in rows if r['blocker'] == 'chain_present']
cov = [r for r in rows if r['blocker'] == 'read_coverage']
print()
print('chain present (join not excluded by coverage): %d gaps, %.1f kb' % (
    len(pol), sum(r['gap_bp'] for r in pol) / 1e3))
print('read-coverage-limited (no site policy can help): %d gaps, %.1f kb, %.1f kb of it unbridgeable' % (
    len(cov), sum(r['gap_bp'] for r in cov) / 1e3, sum(r['unbridged_bp'] for r in cov) / 1e3))
