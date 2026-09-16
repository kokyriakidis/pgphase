#!/usr/bin/env python3
"""Compare graph-only arms: coverage, accuracy, gaps, and the read-level gate.

Aggregate accuracy can improve while individual reads are corrupted, so the
deciding number here is the project's gate -- reads that were concordant in the
baseline arm and become discordant in the arm under test, counted over the reads
both arms phased.
"""
import argparse
import collections
import csv
import gzip
import json
from pathlib import Path

p = argparse.ArgumentParser(description=__doc__)
p.add_argument('--baseline', type=Path, required=True, help='arm directory (holds eval/ and gaps.tsv)')
p.add_argument('--baseline-label', default='baseline')
p.add_argument('--arm', type=Path, action='append', required=True)
p.add_argument('--label', action='append', required=True)
p.add_argument('--output', type=Path, required=True)
a = p.parse_args()
if len(a.arm) != len(a.label):
    raise SystemExit('--arm and --label must be paired')


def per_read(root):
    with gzip.open(root / 'eval/per_read.tsv.gz', 'rt') as stream:
        return {r['read_name']: r['status'].lower() for r in csv.DictReader(stream, delimiter='\t')}


def summary(root):
    return json.load((root / 'eval/summary.json').open())


def gaps(root):
    path = root / 'gaps.tsv'
    if not path.exists():
        return None
    rows = list(csv.DictReader(path.open(), delimiter='\t'))
    return {'gaps': len(rows), 'gap_span_bp': sum(int(r['gap_bp']) for r in rows),
            'untagged_in_gaps': sum(int(r['reads_untagged']) for r in rows),
            'gaps_over_200kb': sum(1 for r in rows if int(r['gap_bp']) >= 200000)}


base_reads, base_summary, base_gaps = per_read(a.baseline), summary(a.baseline), gaps(a.baseline)
METRICS = ('total_phased_reads', 'unphased_reads', 'total_phase_sets', 'concordant_reads',
           'discordant_reads', 'hamming_error_rate', 'phase_block_n50_bp', 'phase_block_aun_bp')
rows = []
for root, label in [(a.baseline, a.baseline_label)] + list(zip(a.arm, a.label)):
    s, g, reads = summary(root), gaps(root), per_read(root)
    shared = [n for n in base_reads if n in reads]
    corrupted = sum(1 for n in shared if base_reads[n] == 'concordant' and reads[n] == 'discordant')
    repaired = sum(1 for n in shared if base_reads[n] == 'discordant' and reads[n] == 'concordant')
    row = {'arm': label}
    row.update({m: s[m] for m in METRICS})
    row.update(g or {})
    row.update({'reads_shared_with_baseline': len(shared),
                'gate_concordant_to_discordant': corrupted,
                'discordant_to_concordant': repaired,
                'reads_gained_vs_baseline': len([n for n in reads if n not in base_reads]),
                'reads_lost_vs_baseline': len([n for n in base_reads if n not in reads])})
    rows.append(row)

with a.output.open('w', newline='') as stream:
    w = csv.DictWriter(stream, fieldnames=list(rows[0]), delimiter='\t')
    w.writeheader()
    w.writerows(rows)

print('%-9s %9s %8s %5s %8s %8s %7s %9s %6s %8s %7s %7s' % (
    'arm', 'phased', 'unphased', 'PS', 'concord', 'discord', 'hamm%', 'N50',
    'gaps', 'gapSpan', 'GATE', 'fixed'))
for r in rows:
    print('%-9s %9d %8d %5d %8d %8d %7.3f %9d %6s %8s %7d %7d' % (
        r['arm'], r['total_phased_reads'], r['unphased_reads'], r['total_phase_sets'],
        r['concordant_reads'], r['discordant_reads'], 100 * r['hamming_error_rate'],
        r['phase_block_n50_bp'], r.get('gaps', '-'),
        ('%.2f Mb' % (r['gap_span_bp'] / 1e6)) if 'gap_span_bp' in r else '-',
        r['gate_concordant_to_discordant'], r['discordant_to_concordant']))
print()
print('GATE = reads concordant in %s that become discordant in the arm' % a.baseline_label)
