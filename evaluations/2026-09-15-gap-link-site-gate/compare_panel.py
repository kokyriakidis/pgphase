#!/usr/bin/env python3
"""Compare gap recovery before and after removing the private-SNP flag gates.

For each gap window, report the tier verdict, whether one phase set spans the
gap, and the read-level regression gate matched by read name. A join that turns
concordant reads discordant is a regression whatever it does to block sizes, so
the gate is reported per gap rather than only in aggregate.
"""
import argparse
import collections
import csv
import subprocess
from pathlib import Path

p = argparse.ArgumentParser(description=__doc__)
p.add_argument('--before-root', type=Path, required=True)
p.add_argument('--after-root', type=Path, required=True)
p.add_argument('--gaps', type=Path, required=True)
p.add_argument('--truth-map', type=Path, required=True)
p.add_argument('--arm', default='union_recover')
p.add_argument('--output', type=Path, required=True)
a = p.parse_args()

truth = {}
with a.truth_map.open() as stream:
    for line in stream:
        name, hap = line.rstrip('\n').split('\t')
        truth[name] = hap
gaps = {int(r['gap_left']): r for r in csv.DictReader(a.gaps.open(), delimiter='\t')}


def tags(bam):
    done = subprocess.run(['samtools', 'view', str(bam)], capture_output=True, text=True)
    if done.returncode != 0:
        raise SystemExit('samtools view failed on %s: %s' % (bam, done.stderr.strip()[:200]))
    out = {}
    for line in done.stdout.splitlines():
        f = line.split('\t')
        hp = ps = None
        for x in f[11:]:
            if x.startswith('HP:i:'):
                hp = x[5:]
            elif x.startswith('PS:i:'):
                ps = x[5:]
        if hp is not None and ps is not None:
            out[f[0]] = (hp, ps)
    return out


def correctness(t):
    """Per-read correct/incorrect under each phase set's majority orientation."""
    by_ps = collections.defaultdict(list)
    for q, (hp, ps) in t.items():
        if q in truth:
            by_ps[ps].append((q, hp))
    out = {}
    for ps, members in by_ps.items():
        c = collections.Counter((hp, truth[q]) for q, hp in members)
        straight = c[('1', 'MATERNAL')] + c[('2', 'PATERNAL')]
        flip = c[('1', 'PATERNAL')] + c[('2', 'MATERNAL')]
        exp = ({'1': 'MATERNAL', '2': 'PATERNAL'} if straight >= flip
               else {'1': 'PATERNAL', '2': 'MATERNAL'})
        for q, hp in members:
            out[q] = exp[hp] == truth[q]
    return out


def ps_spans(vcf):
    spans = {}
    with vcf.open() as stream:
        for line in stream:
            if line.startswith('#'):
                continue
            f = line.rstrip('\n').split('\t')
            if len(f) < 10 or '|' not in f[9].split(':')[0]:
                continue
            keys, values = f[8].split(':'), f[9].split(':')
            if 'PS' not in keys:
                continue
            ps, pos = values[keys.index('PS')], int(f[1])
            lo, hi, n = spans.get(ps, (pos, pos, 0))
            spans[ps] = (min(lo, pos), max(hi, pos), n + 1)
    return spans


def verdict(tiers):
    if not tiers.exists():
        return 'missing'
    rows = list(csv.reader(tiers.open(), delimiter='\t'))[1:]
    states = [r[19] for r in rows if len(r) > 19]
    return 'joined' if 'joined' in states else (states[-1] if states else 'none')


rows = []
for gap_left, gap in sorted(gaps.items()):
    bd, ad = a.before_root / f'gap_{gap_left}' / a.arm, a.after_root / f'gap_{gap_left}' / a.arm
    if not (bd / 'phased.bam').exists() or not (ad / 'phased.bam').exists():
        continue
    gap_right = int(gap['gap_right'])
    cb, ca = correctness(tags(bd / 'phased.bam')), correctness(tags(ad / 'phased.bam'))
    sb, sa = ps_spans(bd / 'native.vcf'), ps_spans(ad / 'native.vcf')
    spans_b = any(s[0] <= gap_left and s[1] >= gap_right for s in sb.values())
    spans_a = any(s[0] <= gap_left and s[1] >= gap_right for s in sa.values())
    x = collections.Counter()
    for q in set(cb) | set(ca):
        x[(cb.get(q), ca.get(q))] += 1
    rows.append({
        'gap_left': gap_left, 'gap_bp': int(gap['gap_bp']),
        'verdict_before': verdict(bd / 'tiers.tsv'), 'verdict_after': verdict(ad / 'tiers.tsv'),
        'phase_sets_before': len(sb), 'phase_sets_after': len(sa),
        'spans_before': int(spans_b), 'spans_after': int(spans_a),
        'tagged_before': len(cb), 'tagged_after': len(ca),
        'correct_before': sum(cb.values()), 'correct_after': sum(ca.values()),
        'gate_conc_to_disc': x[(True, False)], 'lost_tag_was_concordant': x[(True, None)],
        'gained_tag_concordant': x[(None, True)], 'gained_tag_discordant': x[(None, False)]})

with a.output.open('w', newline='') as stream:
    w = csv.DictWriter(stream, fieldnames=list(rows[0]), delimiter='\t', lineterminator='\n')
    w.writeheader()
    w.writerows(rows)

print('%-11s %8s %-9s %-9s %6s %6s %9s %8s %7s %7s' % (
    'gap_left', 'gap_bp', 'before', 'after', 'spanB', 'spanA', 'GATE c->d', 'lostConc', 'gainC', 'gainD'))
for r in rows:
    print('%-11d %8d %-9s %-9s %6s %6s %9d %8d %7d %7d' % (
        r['gap_left'], r['gap_bp'], r['verdict_before'], r['verdict_after'],
        'yes' if r['spans_before'] else 'no', 'yes' if r['spans_after'] else 'no',
        r['gate_conc_to_disc'], r['lost_tag_was_concordant'],
        r['gained_tag_concordant'], r['gained_tag_discordant']))
print()
print('gaps: %d   newly joined: %d   newly spanning: %d' % (
    len(rows),
    sum(1 for r in rows if r['verdict_after'] == 'joined' and r['verdict_before'] != 'joined'),
    sum(1 for r in rows if r['spans_after'] and not r['spans_before'])))
print('GATE concordant -> discordant, total: %d' % sum(r['gate_conc_to_disc'] for r in rows))
print('concordant reads that lost their tag: %d   newly tagged: %d concordant / %d discordant' % (
    sum(r['lost_tag_was_concordant'] for r in rows),
    sum(r['gained_tag_concordant'] for r in rows), sum(r['gained_tag_discordant'] for r in rows)))
print('tagged reads: %d -> %d   concordant: %d -> %d' % (
    sum(r['tagged_before'] for r in rows), sum(r['tagged_after'] for r in rows),
    sum(r['correct_before'] for r in rows), sum(r['correct_after'] for r in rows)))
