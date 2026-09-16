#!/usr/bin/env python3
"""Do the competitor's in-gap sites carry real haplotype signal in OUR alignment?

Genotypes each site directly from the surjected BAM and scores how well its
allele partition follows the read-level truth. Substitutions are read as the base
at the position; indels are read as the net insert-minus-delete length over a
window, never at the anchor position, because the aligner places a repeat-context
indel arbitrarily within its run.

The classes we already trust (clean het SNP, clean het indel) act as the control:
if they do not come out informative, the genotyping is wrong rather than the
sites. Every measurement here should be read against that control.
"""
import argparse
import collections
import csv
import subprocess
from pathlib import Path

CONSUME_REF = frozenset('MDN=X')

p = argparse.ArgumentParser(description=__doc__,
                            formatter_class=argparse.RawDescriptionHelpFormatter)
p.add_argument('--sites', type=Path, required=True, help='site_deficit_*.tsv')
p.add_argument('--bam', type=Path, required=True)
p.add_argument('--contig', default='CHM13#0#chr20')
p.add_argument('--truth-map', type=Path, required=True)
p.add_argument('--flank', type=int, default=25, help='window for the net-length test')
p.add_argument('--min-reads', type=int, default=10)
p.add_argument('--output', type=Path, required=True)
a = p.parse_args()

truth = dict(line.rstrip('\n').split('\t') for line in a.truth_map.open())


def cigar_ops(cig):
    n = ''
    for ch in cig:
        if ch.isdigit():
            n += ch
        else:
            yield int(n), ch
            n = ''


def fetch(pos, span):
    lo = max(1, pos - a.flank - 50)
    hi = pos + span + a.flank + 50
    out = subprocess.run(['samtools', 'view', '-q', '1', str(a.bam),
                          f'{a.contig}:{lo}-{hi}'],
                         capture_output=True, text=True, check=True).stdout
    for line in out.splitlines():
        f = line.split('\t')
        if len(f) < 10:
            continue
        yield f[0], int(f[3]), f[5], f[9]


def net_indel(start, cig, lo, hi):
    ref = start
    delta = 0
    for length, op in cigar_ops(cig):
        if op == 'I':
            if lo <= ref <= hi:
                delta += length
        elif op in 'DN':
            if ref + length >= lo and ref <= hi:
                delta -= length
            ref += length
        elif op in CONSUME_REF:
            ref += length
        elif op == 'S' or op == 'H':
            continue
    return delta, start, ref


def base_at(start, cig, seq, pos):
    ref = start
    qi = 0
    for length, op in cigar_ops(cig):
        if op in 'M=X':
            if ref <= pos < ref + length:
                return seq[qi + (pos - ref)]
            ref += length
            qi += length
        elif op == 'I':
            qi += length
        elif op in 'DN':
            if ref <= pos < ref + length:
                return None
            ref += length
        elif op == 'S':
            qi += length
    return None


rows = []
with a.sites.open() as f:
    sites = list(csv.DictReader(f, delimiter='\t'))

for s in sites:
    pos = int(s['pos'])
    ref, alt = s['ref'], s['alt']
    kind = s['their_type']
    span = max(len(ref), 1)
    calls = {}
    if kind == 'SNP':
        for name, start, cig, seq in fetch(pos, span):
            # SAM POS and VCF POS are both 1-based, and the reference walk
            # starts at SAM POS, so the site's base sits at ref == pos. Passing
            # pos - 1 reads the neighbouring base and scored the clean-het-SNP
            # control at 0.542 with 0% informative -- which is what that control
            # is for.
            b = base_at(start, cig, seq, pos)
            if b is None:
                continue
            if b.upper() == alt.upper():
                calls[name] = 1
            elif b.upper() == ref.upper():
                calls[name] = 0
    else:
        expect = len(alt) - len(ref)
        lo, hi = pos - a.flank, pos + span + a.flank
        for name, start, cig, seq in fetch(pos, span):
            delta, rs, re = net_indel(start, cig, lo, hi)
            if rs > lo or re < hi:
                continue
            if abs(delta - expect) < abs(delta):
                calls[name] = 1
            elif abs(delta) < abs(delta - expect):
                calls[name] = 0
    scored = [(v, truth[q]) for q, v in calls.items() if q in truth]
    row = dict(s)
    row['reads_called'] = len(calls)
    row['reads_scored'] = len(scored)
    if len(scored) >= a.min_reads:
        c = collections.Counter(scored)
        straight = c[(1, 'MATERNAL')] + c[(0, 'PATERNAL')]
        flipped = c[(1, 'PATERNAL')] + c[(0, 'MATERNAL')]
        row['segregation'] = round(max(straight, flipped) / len(scored), 4)
        alt_n = sum(1 for v, _ in scored if v == 1)
        row['alt_fraction'] = round(alt_n / len(scored), 4)
    else:
        row['segregation'] = ''
        row['alt_fraction'] = ''
    rows.append(row)

with a.output.open('w', newline='') as stream:
    w = csv.DictWriter(stream, fieldnames=list(rows[0]), delimiter='\t',
                       lineterminator='\n')
    w.writeheader()
    w.writerows(rows)

groups = collections.defaultdict(list)
for r in rows:
    if r['segregation'] == '':
        continue
    groups[r['our_category']].append(float(r['segregation']))
print(f'sites: {len(rows)}   scorable (>= {a.min_reads} reads with truth): '
      f'{sum(1 for r in rows if r["segregation"] != "")}')
print()
print('%-22s %6s %10s %12s %12s' % ('our category', 'n', 'median seg',
                                     'informative', 'phantom'))
for cat, vals in sorted(groups.items(), key=lambda kv: -len(kv[1])):
    vals.sort()
    med = vals[len(vals) // 2]
    inf = sum(1 for v in vals if v >= 0.90)
    ph = sum(1 for v in vals if v < 0.70)
    print('%-22s %6d %10.3f %11.0f%% %11.0f%%' % (
        cat, len(vals), med, 100 * inf / len(vals), 100 * ph / len(vals)))
unscorable = [r for r in rows if r['segregation'] == '']
print()
print(f'unscorable: {len(unscorable)}  '
      f'{dict(collections.Counter(r["our_category"] for r in unscorable))}')
print(f'wrote {a.output}')
