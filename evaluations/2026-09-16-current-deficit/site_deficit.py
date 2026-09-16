#!/usr/bin/env python3
"""In the gaps a competitor closes correctly, what do we call at its sites?

Restricted to the deficit gaps where a competitor spans at high read-level
accuracy, so the comparison is against spans that are demonstrably right. For
each phased heterozygous record the competitor places inside such a gap, this
finds our own candidate at that position and reports its category, or records the
site as absent from our table entirely. The aggregate says which admission
decision costs us these gaps: a classifier that demotes the site, a filter that
drops it, or discovery that never proposes it.
"""
import argparse
import collections
import csv
import gzip
from pathlib import Path

p = argparse.ArgumentParser(description=__doc__,
                            formatter_class=argparse.RawDescriptionHelpFormatter)
p.add_argument('--scored', type=Path, required=True, help='deficit_scored.tsv')
p.add_argument('--candidates', type=Path, required=True, help='our candidates.tsv')
p.add_argument('--competitor-vcf', type=Path, required=True)
p.add_argument('--competitor-contig', default='chr20')
p.add_argument('--tool', default='hiphase')
p.add_argument('--min-acc', type=float, default=0.98)
p.add_argument('--min-scored', type=int, default=50)
p.add_argument('--tolerance', type=int, default=10,
               help='bp window for matching our candidate to their position')
p.add_argument('--output', type=Path, required=True)
a = p.parse_args()

gaps = []
with a.scored.open() as f:
    for r in csv.DictReader(f, delimiter='\t'):
        if (int(r[f'{a.tool}_spans']) and int(r[f'{a.tool}_scored']) >= a.min_scored
                and float(r[f'{a.tool}_acc']) >= a.min_acc):
            gaps.append((int(r['gap_left']), int(r['gap_right'])))
gaps.sort()
print(f'{a.tool} spans {len(gaps)} deficit gaps at >= {100 * a.min_acc:.0f}% '
      f'({sum(hi - lo for lo, hi in gaps) / 1e6:.2f} Mb)')

ours = collections.defaultdict(list)
with a.candidates.open() as f:
    reader = csv.DictReader(f, delimiter='\t')
    for r in reader:
        pos = int(r['POS'])
        if any(lo - a.tolerance <= pos <= hi + a.tolerance for lo, hi in gaps):
            ours[pos].append(r)

def ours_at(pos):
    for d in range(a.tolerance + 1):
        for q in (pos - d, pos + d):
            if q in ours:
                return ours[q][0]
    return None

rows = []
kinds = collections.Counter()
by_gap = collections.Counter()
with gzip.open(a.competitor_vcf, 'rt') as f:
    for line in f:
        if line.startswith('#'):
            continue
        c = line.rstrip('\n').split('\t')
        if c[0] != a.competitor_contig:
            continue
        pos = int(c[1])
        hit = next(((lo, hi) for lo, hi in gaps if lo <= pos <= hi), None)
        if hit is None:
            continue
        gt = c[9].split(':')[0]
        if '|' not in gt:
            continue
        alt_a, alt_b = gt.split('|')[:2]
        if alt_a == alt_b:
            continue
        ref, alt = c[3], c[4].split(',')[0]
        their_kind = 'SNP' if len(ref) == 1 and len(alt) == 1 else (
            'INS' if len(alt) > len(ref) else 'DEL')
        mine = ours_at(pos)
        cat = mine['CATEGORY'] if mine else 'ABSENT'
        rows.append({'gap_left': hit[0], 'gap_right': hit[1], 'pos': pos,
                     'their_type': their_kind, 'ref': ref[:20], 'alt': alt[:20],
                     'our_category': cat,
                     'our_type': mine['TYPE'] if mine else '',
                     'our_dp': mine.get('DP', '') if mine else '',
                     'our_af': mine.get('AF', '') if mine else ''})
        kinds[(their_kind, cat)] += 1
        by_gap[hit] += 1

with a.output.open('w', newline='') as stream:
    w = csv.DictWriter(stream, fieldnames=list(rows[0]), delimiter='\t',
                       lineterminator='\n')
    w.writeheader()
    w.writerows(rows)

print(f'their phased het sites inside those gaps: {len(rows)}')
print(f'  by their variant type: '
      f'{dict(collections.Counter(r["their_type"] for r in rows))}')
print()
print('what WE call at each of those positions:')
print('  %-8s %-22s %6s' % ('theirs', 'our category', 'n'))
for (kind, cat), n in kinds.most_common():
    print('  %-8s %-22s %6d' % (kind, cat, n))
print()
usable = {'CLEAN_HET_SNP', 'CLEAN_HET_INDEL'}
tot = len(rows)
absent = sum(1 for r in rows if r['our_category'] == 'ABSENT')
demoted = sum(1 for r in rows if r['our_category'] not in usable
              and r['our_category'] != 'ABSENT')
print(f'  usable to us today (clean het): '
      f'{sum(1 for r in rows if r["our_category"] in usable)} of {tot}')
print(f'  present but not clean het     : {demoted}')
print(f'  absent from our table entirely : {absent}')
print()
print('per-gap: how many of their sites we hold as clean het')
print('  %-11s %-11s %7s %7s %7s %7s' % ('gap_left', 'gap_right', 'theirs',
                                          'cleanHet', 'notClean', 'absent'))
for (lo, hi), n in sorted(by_gap.items()):
    sel = [r for r in rows if r['gap_left'] == lo]
    print('  %-11d %-11d %7d %7d %7d %7d' % (
        lo, hi, n,
        sum(1 for r in sel if r['our_category'] in usable),
        sum(1 for r in sel if r['our_category'] not in usable
            and r['our_category'] != 'ABSENT'),
        sum(1 for r in sel if r['our_category'] == 'ABSENT')))
print(f'wrote {a.output}')
