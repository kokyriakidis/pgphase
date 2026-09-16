#!/usr/bin/env python3
"""Which gaps does the CURRENT pipeline leave open that a competitor closes?

The earlier deficit list was derived from the graph-only pass-1 gap inventory,
which overstates the remaining work: the hybrid path already spans many of those
gaps outright (chr20:30,794,962-30,814,005, which all three competitors close, is
covered by a single 76.1 kb phase set at 99.6% read accuracy, and gap recovery is
never even asked about it). This measures the deficit against what the pipeline
actually emits.

Blocks come from the phased VCF: a phase set spans its first to its last phased
variant. Blocks may overlap, so a gap is only counted where the next block starts
beyond the running maximum end. For each gap, a competitor "spans" it when one of
its phase sets covers the whole interval; interior sites are its phased
heterozygous records inside the gap.
"""
import argparse
import collections
import csv
import gzip
from pathlib import Path


def read_blocks(path, contig=None):
    """phase set -> (first phased variant, last phased variant, count)."""
    opener = gzip.open if str(path).endswith('.gz') else open
    blocks = {}
    with opener(path, 'rt') as stream:
        for line in stream:
            if line.startswith('#'):
                continue
            c = line.rstrip('\n').split('\t')
            if len(c) < 10:
                continue
            if contig is not None and c[0] != contig:
                continue
            gt = c[9].split(':')[0]
            if '|' not in gt:
                continue
            keys = c[8].split(':')
            if 'PS' not in keys:
                continue
            ps = c[9].split(':')[keys.index('PS')].strip()
            if ps in ('.', '', '0'):
                continue
            pos = int(c[1])
            lo, hi, n = blocks.get(ps, (pos, pos, 0))
            blocks[ps] = (min(lo, pos), max(hi, pos), n + 1)
    return blocks


def het_records(path, contig=None):
    """Sorted (pos, is_snp, phase_set) for every phased het record."""
    opener = gzip.open if str(path).endswith('.gz') else open
    out = []
    with opener(path, 'rt') as stream:
        for line in stream:
            if line.startswith('#'):
                continue
            c = line.rstrip('\n').split('\t')
            if len(c) < 10:
                continue
            if contig is not None and c[0] != contig:
                continue
            gt = c[9].split(':')[0]
            if '|' not in gt:
                continue
            a, b = gt.split('|')[:2]
            if a == b:
                continue
            keys = c[8].split(':')
            ps = (c[9].split(':')[keys.index('PS')].strip()
                  if 'PS' in keys else '')
            out.append((int(c[1]), len(c[3]) == 1 and len(c[4]) == 1, ps))
    out.sort()
    return out


p = argparse.ArgumentParser(description=__doc__,
                            formatter_class=argparse.RawDescriptionHelpFormatter)
p.add_argument('--ours', type=Path, required=True, help='our phased VCF')
p.add_argument('--our-contig', default=None)
p.add_argument('--competitor', action='append', default=[],
               help='name=path/to/phased.vcf.gz, repeatable')
p.add_argument('--competitor-contig', default='chr20')
p.add_argument('--min-gap-bp', type=int, default=1,
               help='ignore gaps narrower than this')
p.add_argument('--output', type=Path, required=True)
a = p.parse_args()

ours = read_blocks(a.ours, a.our_contig)
ordered = sorted(ours.values())
gaps = []
run_hi = None
for lo, hi, _ in ordered:
    if run_hi is not None and lo > run_hi and lo - run_hi >= a.min_gap_bp:
        gaps.append((run_hi, lo))
    run_hi = hi if run_hi is None else max(run_hi, hi)
print(f'our phase blocks: {len(ours)}   phased span '
      f'{sum(hi - lo for lo, hi, _ in ordered) / 1e6:.2f} Mb')
print(f'gaps between consecutive blocks: {len(gaps)}   '
      f'total gap span {sum(hi - lo for lo, hi in gaps) / 1e6:.2f} Mb')

tools = {}
for spec in a.competitor:
    name, _, path = spec.partition('=')
    tools[name] = (read_blocks(path, a.competitor_contig),
                   het_records(path, a.competitor_contig))
    print(f'  {name}: {len(tools[name][0])} phase sets, '
          f'{len(tools[name][1])} phased hets')

rows = []
for lo, hi in gaps:
    row = {'gap_left': lo, 'gap_right': hi, 'gap_bp': hi - lo}
    spanning = 0
    for name, (blocks, hets) in tools.items():
        spans = any(blo <= lo and bhi >= hi for blo, bhi, _ in blocks.values())
        inside = [r for r in hets if lo <= r[0] <= hi]
        row[f'{name}_spans'] = int(spans)
        row[f'{name}_interior_sites'] = len(inside)
        row[f'{name}_interior_snps'] = sum(1 for r in inside if r[1])
        spanning += int(spans)
    row['n_tools_spanning'] = spanning
    rows.append(row)

with a.output.open('w', newline='') as stream:
    w = csv.DictWriter(stream, fieldnames=list(rows[0]), delimiter='\t',
                       lineterminator='\n')
    w.writeheader()
    w.writerows(rows)

deficit = [r for r in rows if r['n_tools_spanning'] >= 1]
print()
print(f'DEFICIT: {len(deficit)} of {len(rows)} gaps are spanned by at least one '
      f'competitor, {sum(r["gap_bp"] for r in deficit) / 1e6:.2f} Mb')
by = collections.Counter(r['n_tools_spanning'] for r in deficit)
print(f'  by number of tools spanning: {dict(sorted(by.items()))}')
none = [r for r in rows if r['n_tools_spanning'] == 0]
print(f'  spanned by nobody: {len(none)} gaps, '
      f'{sum(r["gap_bp"] for r in none) / 1e6:.2f} Mb')
print()
print('largest deficit gaps:')
hdr = ['gap_left', 'gap_right', 'gap_bp', 'n_tools_spanning'] + \
      [f'{n}_interior_sites' for n in tools]
print('  ' + ' '.join(f'{h:>12s}' for h in hdr))
for r in sorted(deficit, key=lambda r: -r['gap_bp'])[:15]:
    print('  ' + ' '.join(f'{r[h]:>12d}' for h in hdr))
print()
print('deficit gaps with the most competitor interior sites:')
key = f'{list(tools)[0]}_interior_sites' if tools else 'gap_bp'
for r in sorted(deficit, key=lambda r: -r[key])[:15]:
    print('  ' + ' '.join(f'{r[h]:>12d}' for h in hdr))
print(f'wrote {a.output}')
