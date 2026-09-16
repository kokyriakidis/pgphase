#!/usr/bin/env python3
"""Are the competitor spans across our deficit gaps actually correct?

A gap a competitor spans is only an opportunity if its span is right. For each
deficit gap this takes the reads overlapping it in the surjected BAM, looks each
one up in the frozen per-tool read evaluations (scored against the diplinator
read truth), and reports the spanning tools' read-level accuracy inside that gap
alongside our own. A tool that spans a gap at chance is not showing us signal we
are missing; it is guessing where we abstain.
"""
import argparse
import collections
import csv
import gzip
import subprocess
from pathlib import Path

p = argparse.ArgumentParser(description=__doc__,
                            formatter_class=argparse.RawDescriptionHelpFormatter)
p.add_argument('--deficit', type=Path, required=True)
p.add_argument('--bam', type=Path, required=True)
p.add_argument('--contig', default='CHM13#0#chr20')
p.add_argument('--eval-root', type=Path, required=True)
p.add_argument('--tools', default='hiphase,longphase,whatshap_opt')
p.add_argument('--output', type=Path, required=True)
a = p.parse_args()

tools = a.tools.split(',')
status = {}
for tool in tools:
    path = a.eval_root / f'{tool}_reads' / 'per_read.tsv.gz'
    with gzip.open(path, 'rt') as f:
        status[tool] = {r['read_name']: r['status'].lower()
                        for r in csv.DictReader(f, delimiter='\t')}

rows = []
with a.deficit.open() as f:
    deficit = [r for r in csv.DictReader(f, delimiter='\t')
               if int(r['n_tools_spanning']) >= 1]

for r in deficit:
    lo, hi = int(r['gap_left']), int(r['gap_right'])
    names = subprocess.run(['samtools', 'view', str(a.bam), f'{a.contig}:{lo}-{hi}'],
                           capture_output=True, text=True, check=True).stdout
    reads = {line.split('\t', 1)[0] for line in names.splitlines()}
    out = {'gap_left': lo, 'gap_right': hi, 'gap_bp': hi - lo,
           'reads_in_gap': len(reads),
           'n_tools_spanning': int(r['n_tools_spanning'])}
    for tool in tools:
        scored = [status[tool][q] for q in reads if q in status[tool]]
        conc = sum(1 for s in scored if s == 'concordant')
        out[f'{tool}_spans'] = int(r[f'{tool}_spans'])
        out[f'{tool}_scored'] = len(scored)
        out[f'{tool}_concordant'] = conc
        out[f'{tool}_acc'] = round(conc / len(scored), 4) if scored else 0.0
    rows.append(out)

with a.output.open('w', newline='') as stream:
    w = csv.DictWriter(stream, fieldnames=list(rows[0]), delimiter='\t',
                       lineterminator='\n')
    w.writeheader()
    w.writerows(rows)

print(f'deficit gaps scored: {len(rows)}')
for tool in tools:
    sel = [r for r in rows if r[f'{tool}_spans'] and r[f'{tool}_scored'] >= 20]
    if not sel:
        print(f'  {tool}: no spanned gap with >=20 scored reads')
        continue
    sc = sum(r[f'{tool}_scored'] for r in sel)
    co = sum(r[f'{tool}_concordant'] for r in sel)
    good = [r for r in sel if r[f'{tool}_acc'] >= 0.98]
    poor = [r for r in sel if r[f'{tool}_acc'] < 0.90]
    print(f'  {tool}: spans {len(sel)} scorable deficit gaps, '
          f'{co}/{sc} reads concordant ({100 * co / sc:.2f}%)')
    print(f'      >=98% accurate in {len(good)} gaps ({sum(r["gap_bp"] for r in good) / 1e6:.2f} Mb), '
          f'<90% in {len(poor)} ({sum(r["gap_bp"] for r in poor) / 1e6:.2f} Mb)')

print()
print('per-gap detail, worst competitor accuracy first:')
hdr = '%-11s %-11s %8s %6s ' % ('gap_left', 'gap_right', 'gap_bp', 'reads')
hdr += ' '.join('%18s' % t for t in tools)
print('  ' + hdr)
def worst(r):
    vals = [r[f'{t}_acc'] for t in tools if r[f'{t}_spans'] and r[f'{t}_scored'] >= 20]
    return min(vals) if vals else 1.0
for r in sorted(rows, key=worst)[:20]:
    line = '%-11d %-11d %8d %6d ' % (r['gap_left'], r['gap_right'], r['gap_bp'],
                                     r['reads_in_gap'])
    for t in tools:
        line += '%18s' % (('SPAN %5.1f%% n=%-4d' % (100 * r[f'{t}_acc'], r[f'{t}_scored']))
                          if r[f'{t}_spans'] else '-')
    print('  ' + line)
print(f'wrote {a.output}')
