#!/usr/bin/env python3
"""Diagnose one deficit gap: what we do, what a competitor does, and why we stop.

Runs against the outputs of one windowed `collect-hybrid-variation --recover-gaps
--gap-decision-audit` invocation. Reports our residual break inside the gap, the
tier verdicts, the competitor's site inventory matched against our candidates
(allowing for indel anchor offsets), and for every site inside the break its
category, verification flags, truth segregation and which recovery tier may admit
it -- plus whether reads chain across the break and whether the competitor's own
phasing there is correct. Site-admission questions and read-coverage questions
look identical in a tier report and need separating before any fix is proposed.
"""
import argparse
import bisect
import collections
import csv
import glob
import gzip
from pathlib import Path

p = argparse.ArgumentParser(description=__doc__)
p.add_argument('--run-dir', type=Path, required=True)
p.add_argument('--gap-left', type=int, required=True)
p.add_argument('--gap-right', type=int, required=True)
p.add_argument('--truth-map', type=Path, required=True)
p.add_argument('--competitor-vcf', type=Path, required=True)
p.add_argument('--competitor-per-read', type=Path, required=True)
p.add_argument('--competitor-name', default='hiphase')
p.add_argument('--contig', default='chr20')
p.add_argument('--offset-tolerance', type=int, default=15)
p.add_argument('--output', type=Path)
a = p.parse_args()

GL, GR = a.gap_left, a.gap_right
truth = {}
with a.truth_map.open() as stream:
    for line in stream:
        name, hap = line.rstrip('\n').split('\t')
        truth[name] = hap

CAT = {0x004: 'CleanHetSnp', 0x008: 'CleanHetIndel', 0x010: 'RepHetIndel',
       0x080: 'CleanHom', 0x100: 'NoisyMsaHet', 0x200: 'NoisyMsaHom',
       0x1000: 'NonAnchorHet'}
TYPE = {8: 'SNP', 1: 'INS', 2: 'DEL'}


def ps_spans(vcf):
    spans = {}
    with vcf.open() as stream:
        for line in stream:
            if line.startswith('#'):
                continue
            c = line.rstrip('\n').split('\t')
            if len(c) < 10 or '|' not in c[9].split(':')[0]:
                continue
            keys, values = c[8].split(':'), c[9].split(':')
            if 'PS' not in keys:
                continue
            ps, pos = values[keys.index('PS')].strip(), int(c[1])
            lo, hi, n = spans.get(ps, (pos, pos, 0))
            spans[ps] = (min(lo, pos), max(hi, pos), n + 1)
    return spans


def residual_break(spans):
    """Largest stretch inside the gap covered by no single phase set."""
    iv = sorted((max(s[0], GL), min(s[1], GR)) for s in spans.values()
                if s[1] >= GL and s[0] <= GR)
    if not iv:
        return GL, GR
    merged = [list(iv[0])]
    for b, e in iv[1:]:
        if b <= merged[-1][1]:
            merged[-1][1] = max(merged[-1][1], e)
        else:
            merged.append([b, e])
    best = (GL, merged[0][0]) if merged[0][0] > GL else None
    for x, y in zip(merged, merged[1:]):
        if best is None or (y[0] - x[1]) > (best[1] - best[0]):
            best = (x[1], y[0])
    if merged[-1][1] < GR and (best is None or (GR - merged[-1][1]) > (best[1] - best[0])):
        best = (merged[-1][1], GR)
    return best if best else (GL, GR)


spans = ps_spans(a.run_dir / 'native.vcf')
inside = sorted(((s[0], s[1], s[2]) for s in spans.values() if s[1] >= GL and s[0] <= GR))
spanning = [s for s in spans.values() if s[0] <= GL and s[1] >= GR]
BL, BR = residual_break(spans)
print('gap %d-%d (%.1f kb)' % (GL, GR, (GR - GL) / 1e3))
print('  our blocks overlapping it: %s' % ' | '.join(
    '%d-%d (%d sites)' % (lo, hi, n) for lo, hi, n in inside))
print('  spanning phase set: %s' % ('yes' if spanning else 'NONE'))
print('  residual break: %d-%d (%.1f kb)' % (BL, BR, (BR - BL) / 1e3))

tiers = a.run_dir / 'tiers.tsv'
if tiers.exists():
    rows = list(csv.reader(tiers.open(), delimiter='\t'))[1:]
    rel = [r for r in rows if len(r) > 19 and int(r[2]) >= BL - 1 and int(r[1]) <= BR + 1]
    print('  tier verdicts for the break:')
    for r in rel:
        print('    gap %s-%s tier %s pass %s leftPS=%-11s rightPS=%-11s msaSNP=%s msaIND=%s reads=%-4s %s'
              % (r[1], r[2], r[3], r[17], r[11], r[12], r[7], r[8], r[13], r[19]))

comp = []
op = gzip.open if str(a.competitor_vcf).endswith('.gz') else open
with op(a.competitor_vcf, 'rt') as stream:
    for line in stream:
        if line.startswith('#'):
            continue
        c = line.rstrip('\n').split('\t')
        if len(c) < 10 or c[0] != a.contig or '|' not in c[9].split(':')[0]:
            continue
        pos = int(c[1])
        if GL <= pos <= GR:
            comp.append((pos, c[3], c[4]))
cand = {int(r['POS']): r for r in csv.DictReader((a.run_dir / 'candidates.tsv').open(), delimiter='\t')}
cand_pos = sorted(cand)
exact = off = absent = 0
absent_rows = []
for pos, ref, alt in comp:
    if pos in cand:
        exact += 1
    else:
        i = bisect.bisect_left(cand_pos, pos)
        near = [q for q in cand_pos[max(0, i - 2):i + 2] if abs(q - pos) <= a.offset_tolerance]
        if near:
            off += 1
        else:
            absent += 1
            absent_rows.append((pos, ref, alt))
print('  %s sites in the gap: %d   we hold %d exactly, %d within %d bp, %d absent'
      % (a.competitor_name, len(comp), exact, off, a.offset_tolerance, absent))
for pos, ref, alt in absent_rows[:8]:
    print('    absent: %d  %s > %s' % (pos, ref[:14], alt[:24]))

files = glob.glob(str(a.run_dir / 'tid*.evidence.tsv'))
sites, reads, obs = {}, {}, collections.defaultdict(dict)
for path in files:
    for line in open(path):
        f = line.rstrip('\n').split('\t')
        if f[0] == 'SITE':
            sites[(path, int(f[1]))] = dict(pos=int(f[2]), vtype=int(f[3]), cate=int(f[7]),
                                            msa=int(f[8]), hp=int(f[9]))
        elif f[0] == 'READ':
            reads[(path, int(f[1]))] = dict(q=f[3], beg=int(f[4]), end=int(f[5]), skip=int(f[7]))
        elif f[0] == 'OBS':
            if int(f[3]) >= 0:
                obs[(path, int(f[2]))][(path, int(f[1]))] = int(f[3])


def segregation(key):
    pairs = [(al, truth[reads[r]['q']]) for r, al in obs[key].items()
             if reads[r]['q'] in truth and not reads[r]['skip']]
    if len(pairs) < 10:
        return None, len(pairs)
    c = collections.Counter(pairs)
    alleles = sorted({x for x, _ in pairs})
    if len(alleles) < 2:
        return 0.0, len(pairs)
    s = c[(alleles[0], 'MATERNAL')] + c[(alleles[1], 'PATERNAL')]
    f_ = c[(alleles[0], 'PATERNAL')] + c[(alleles[1], 'MATERNAL')]
    return (max(s, f_) / (s + f_) if s + f_ else 0.0), s + f_


def admitted_by(s):
    """Which recovery tier may admit this site as link evidence."""
    if s['cate'] & 0x004 or s['cate'] & 0x008 or s['cate'] & 0x080:
        return 'all (clean germline)'
    if s['cate'] == 0x100 and s['msa']:
        if s['hp']:
            return '4 only (homopolymer)'
        return '2-4' if s['vtype'] == 8 else '3-4'
    if s['cate'] & 0x010:
        return 'never (repeat-demoted)'
    return 'never'


break_sites = [(k, s) for k, s in sites.items()
               if BL <= s['pos'] <= BR and s['cate'] & (0x004 | 0x008 | 0x010 | 0x100)]
out_rows = []
print('  sites inside the break:')
print('    %-11s %-5s %-4s %-4s %-14s %7s %12s %s'
      % ('pos', 'type', 'msa', 'hp', 'category', 'reads', 'segregation', 'admitted by'))
seen = set()
for k, s in sorted(break_sites, key=lambda kv: kv[1]['pos']):
    if s['pos'] in seen:
        continue
    seen.add(s['pos'])
    v, n = segregation(k)
    cat = '+'.join(nm for b, nm in CAT.items() if s['cate'] & b) or hex(s['cate'])
    adm = admitted_by(s)
    print('    %-11d %-5s %-4d %-4d %-14s %7d %12s %s'
          % (s['pos'], TYPE.get(s['vtype'], s['vtype']), s['msa'], s['hp'], cat, n,
             '%.3f' % v if v is not None else 'too few', adm))
    out_rows.append({'gap_left': GL, 'gap_right': GR, 'break_beg': BL, 'break_end': BR,
                     'pos': s['pos'], 'type': TYPE.get(s['vtype'], s['vtype']),
                     'msa_verified': s['msa'], 'is_homopolymer_indel': s['hp'],
                     'category': cat, 'reads': n,
                     'segregation': round(v, 4) if v is not None else '',
                     'admitted_by_tier': adm})

chain = sorted(seen | {BL, BR})
print('  read spanning across the break chain:')
for x, y in zip(chain, chain[1:]):
    n = sum(1 for r in reads.values() if not r['skip'] and r['beg'] <= x and r['end'] >= y)
    print('    %d -> %d  (%6.1f kb)  spanning reads=%d' % (x, y, (y - x) / 1e3, n))

qn = {reads[k]['q'] for k in reads if not reads[k]['skip']
      and reads[k]['end'] >= BL and reads[k]['beg'] <= BR}
with gzip.open(a.competitor_per_read, 'rt') as stream:
    status = collections.Counter(r['status'].lower() for r in csv.DictReader(stream, delimiter='\t')
                                 if r['read_name'] in qn)
tot = status['concordant'] + status['discordant']
print('  %s on the %d reads overlapping our break: %s%s' % (
    a.competitor_name, len(qn), dict(status),
    '  -> %.2f%% accurate' % (100 * status['concordant'] / tot) if tot else ''))

if a.output and out_rows:
    with a.output.open('w', newline='') as stream:
        w = csv.DictWriter(stream, fieldnames=list(out_rows[0]), delimiter='\t', lineterminator='\n')
        w.writeheader()
        w.writerows(out_rows)
    print('  wrote %s' % a.output)
