#!/usr/bin/env python3
"""Did injecting BAM-derived sites let one phase set span the gap?

For each gap window and arm, take the phase sets from the phased VCF (first to
last phased variant defines a set's span) and ask whether any single set covers
the whole gap interval the graph-only pass left open. Reads are then scored
against the read-level truth so a spanned gap can be distinguished from a
mis-joined one -- spanning the gap wrongly is worse than leaving it open.
"""
import argparse
import collections
import csv
import subprocess
from pathlib import Path

p = argparse.ArgumentParser(description=__doc__)
p.add_argument('--root', type=Path, required=True)
p.add_argument('--gaps', type=Path, required=True)
p.add_argument('--truth-map', type=Path, required=True)
p.add_argument('--arms', nargs='+', default=['union', 'union_recover'])
p.add_argument('--output', type=Path, required=True)
a = p.parse_args()

truth = {}
with a.truth_map.open() as stream:
    for line in stream:
        name, hap = line.rstrip('\n').split('\t')
        truth[name] = hap
gaps = {int(r['gap_left']): r for r in csv.DictReader(a.gaps.open(), delimiter='\t')}


def ps_spans(vcf):
    """phase set -> (first phased variant, last phased variant, count)."""
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


def read_tags(bam):
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


rows = []
for window_dir in sorted(a.root.iterdir()):
    if not window_dir.is_dir() or not window_dir.name.startswith('gap_'):
        continue
    gap_left = int(window_dir.name.split('_')[1])
    gap = gaps.get(gap_left)
    if gap is None:
        continue
    gap_right = int(gap['gap_right'])
    for arm in a.arms:
        d = window_dir / arm
        if not (d / 'native.vcf').exists() or not (d / 'phased.bam').exists():
            continue
        spans = ps_spans(d / 'native.vcf')
        spanning = [(ps, s) for ps, s in spans.items() if s[0] <= gap_left and s[1] >= gap_right]
        tags = read_tags(d / 'phased.bam')

        # Accuracy of the spanning set (or of the largest set when none spans).
        target = max(spanning, key=lambda kv: kv[1][2])[0] if spanning else (
            max(((ps, sum(1 for _, q in tags.values() if q == ps)) for ps in spans),
                key=lambda kv: kv[1])[0] if spans else None)
        members = [r for r, (hp, ps) in tags.items() if ps == target and r in truth]
        counts = collections.Counter((tags[r][0], truth[r]) for r in members)
        straight = counts[('1', 'MATERNAL')] + counts[('2', 'PATERNAL')]
        flipped = counts[('1', 'PATERNAL')] + counts[('2', 'MATERNAL')]
        correct = max(straight, flipped)
        rows.append({
            'gap_left': gap_left, 'gap_bp': int(gap['gap_bp']), 'arm': arm,
            'phase_sets': len(spans), 'spans_gap': int(bool(spanning)),
            'target_ps': target or '', 'target_variants': spans[target][2] if target in spans else 0,
            'reads_tagged': sum(1 for _, q in tags.values() if q == target),
            'reads_scored': len(members), 'reads_correct': correct,
            'accuracy': round(correct / len(members), 4) if members else 0.0})

with a.output.open('w', newline='') as stream:
    w = csv.DictWriter(stream, fieldnames=list(rows[0]), delimiter='\t', lineterminator='\n')
    w.writeheader()
    w.writerows(rows)

print('%-11s %9s %-14s %7s %6s %9s %8s %8s %9s' % (
    'gap_left', 'gap_bp', 'arm', 'phasePS', 'spans', 'targetVar', 'tagged', 'scored', 'accuracy'))
for r in rows:
    print('%-11d %9d %-14s %7d %6s %9d %8d %8d %8.1f%%' % (
        r['gap_left'], r['gap_bp'], r['arm'], r['phase_sets'],
        'YES' if r['spans_gap'] else 'no', r['target_variants'], r['reads_tagged'],
        r['reads_scored'], 100 * r['accuracy']))
print()
for arm in a.arms:
    sel = [r for r in rows if r['arm'] == arm]
    if not sel:
        continue
    spanned = [r for r in sel if r['spans_gap']]
    print('  %-14s spans %d of %d gaps' % (arm, len(spanned), len(sel)), end='')
    if spanned:
        sc = sum(r['reads_scored'] for r in spanned)
        co = sum(r['reads_correct'] for r in spanned)
        print('   reads in spanning sets %d, correct %d (%.2f%%)' % (sc, co, 100 * co / sc))
    else:
        print()
