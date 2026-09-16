#!/usr/bin/env python3
"""Find the gaps a competitor closes and we do not.

Our unresolved gaps are only a competitive deficit where some competitor actually
phases across them. For each gap in the pass-1 inventory, test whether a
competitor phase set spans it end to end and carries sites in its interior --
bracketing a gap without phasing inside it is not a join. Competitors and our
pipeline share CHM13 coordinates, so the gap bounds are used directly.
"""
import argparse
import csv
import gzip
from pathlib import Path

p = argparse.ArgumentParser(description=__doc__)
p.add_argument('--gaps', type=Path, required=True)
p.add_argument('--results-root', type=Path, required=True)
p.add_argument('--tools', nargs='+', default=['hiphase', 'longphase', 'whatshap_opt'])
p.add_argument('--contig', default='chr20')
p.add_argument('--output', type=Path, required=True)
a = p.parse_args()


def phase_sets(vcf, contig):
    """phase set id -> (first site, last site, n sites), merged across the file."""
    spans = {}
    with gzip.open(vcf, 'rt') as stream:
        for line in stream:
            if line.startswith('#'):
                continue
            c = line.rstrip('\n').split('\t')
            if len(c) < 10 or c[0] != contig or '|' not in c[9].split(':')[0]:
                continue
            keys, values = c[8].split(':'), c[9].split(':')
            if 'PS' not in keys:
                continue
            # The PS value can be the final FORMAT field; strip stray whitespace so
            # one phase set does not split into two keys.
            ps = values[keys.index('PS')].strip()
            pos = int(c[1])
            lo, hi, n = spans.get(ps, (pos, pos, 0))
            spans[ps] = (min(lo, pos), max(hi, pos), n + 1)
    return spans


tools = {}
for tool in a.tools:
    vcf = a.results_root / tool / 'phased.vcf.gz'
    if not vcf.exists():
        print('missing: %s' % vcf)
        continue
    tools[tool] = phase_sets(vcf, a.contig)
    print('%-13s %d phase sets' % (tool, len(tools[tool])))

sites_by_tool = {}
for tool in tools:
    positions = []
    with gzip.open(a.results_root / tool / 'phased.vcf.gz', 'rt') as stream:
        for line in stream:
            if line.startswith('#'):
                continue
            c = line.rstrip('\n').split('\t')
            if len(c) < 10 or c[0] != a.contig or '|' not in c[9].split(':')[0]:
                continue
            positions.append(int(c[1]))
    sites_by_tool[tool] = sorted(positions)

import bisect
rows = []
for r in csv.DictReader(a.gaps.open(), delimiter='\t'):
    gl, gr = int(r['gap_left']), int(r['gap_right'])
    row = {'gap_left': gl, 'gap_right': gr, 'gap_bp': int(r['gap_bp'])}
    n_span = 0
    for tool, spans in tools.items():
        best = None
        for ps, s in spans.items():
            if s[0] <= gl and s[1] >= gr:
                interior = bisect.bisect_right(sites_by_tool[tool], gr) - \
                           bisect.bisect_left(sites_by_tool[tool], gl)
                if best is None or s[2] > best[0]:
                    best = (s[2], interior, ps)
        row[f'{tool}_spans'] = int(best is not None)
        row[f'{tool}_interior_sites'] = best[1] if best else 0
        n_span += int(best is not None)
    row['n_tools_spanning'] = n_span
    rows.append(row)

with a.output.open('w', newline='') as stream:
    w = csv.DictWriter(stream, fieldnames=list(rows[0]), delimiter='\t', lineterminator='\n')
    w.writeheader()
    w.writerows(rows)

deficit = [r for r in rows if r['n_tools_spanning'] > 0]
print()
print('our pass-1 gaps: %d spanning %.2f Mb' % (len(rows), sum(r['gap_bp'] for r in rows) / 1e6))
print('DEFICIT (at least one competitor spans it): %d gaps, %.2f Mb' % (
    len(deficit), sum(r['gap_bp'] for r in deficit) / 1e6))
for k in (3, 2, 1):
    sel = [r for r in rows if r['n_tools_spanning'] == k]
    print('  spanned by %d of %d tools: %d gaps, %.2f Mb' % (
        k, len(tools), len(sel), sum(r['gap_bp'] for r in sel) / 1e6))
print('  spanned by none (nobody phases these): %d gaps, %.2f Mb' % (
    len(rows) - len(deficit), sum(r['gap_bp'] for r in rows if r['n_tools_spanning'] == 0) / 1e6))
print()
print('largest deficit gaps (all tools spanning first):')
print('%-11s %-11s %9s %6s %8s %9s %10s' % (
    'gap_left', 'gap_right', 'gap_bp', 'tools', 'hiphase', 'longphase', 'whatshap'))
for r in sorted(deficit, key=lambda x: (-x['n_tools_spanning'], -x['gap_bp']))[:20]:
    print('%-11d %-11d %9d %6d %8d %9d %10d' % (
        r['gap_left'], r['gap_right'], r['gap_bp'], r['n_tools_spanning'],
        r['hiphase_interior_sites'], r['longphase_interior_sites'], r['whatshap_opt_interior_sites']))
