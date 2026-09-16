#!/usr/bin/env python3
"""Score one arm over the test panel, with an optional read-level gate.

Per window this reports what has to move for the window to count as fixed --
whether a single block spans the gap, how many heterozygous sites are phased
inside it, how many reads are tagged and what fraction agree with the read-level
truth -- and, when a baseline arm is given, the gate that decides whether a
change is admissible: a read that was concordant and is now discordant is a
regression no coverage gain offsets.

Read concordance is computed per phase set, orienting each block by its own
majority, because a block's hap labels are arbitrary; a switch inside one block
therefore shows up as reduced accuracy rather than being hidden.
"""
import argparse
import collections
import csv
import subprocess
from pathlib import Path

CONTIG = 'CHM13#0#chr20'

p = argparse.ArgumentParser(description=__doc__)
p.add_argument('--panel', type=Path, required=True)
p.add_argument('--arm', type=Path, required=True, help='OUT dir of run_panel.sh')
p.add_argument('--baseline', type=Path, help='OUT dir of the arm to gate against')
p.add_argument('--truth-map', type=Path, required=True)
p.add_argument('--flank', type=int, default=60000)
p.add_argument('--output', type=Path, required=True)
a = p.parse_args()

truth = dict(line.split('\t') for line in a.truth_map.read_text().splitlines())


def tags(bam, lo, hi):
    out = subprocess.run(['samtools', 'view', str(bam), f'{CONTIG}:{lo}-{hi}'],
                         capture_output=True, text=True)
    if out.returncode != 0:
        raise RuntimeError(f'samtools view failed on {bam}: {out.stderr[:200]}')
    d = {}
    for line in out.stdout.splitlines():
        f = line.split('\t')
        hap = ps = None
        for x in f[11:]:
            if x.startswith('HP:i:'):
                hap = int(x[5:])
            elif x.startswith('PS:i:'):
                ps = x[5:]
        if hap is not None and ps is not None:
            d[f[0]] = (hap, ps)
    return d


def concordance(d):
    """read -> True/False, orienting each phase set by its own majority."""
    per = collections.defaultdict(list)
    for q, (hap, ps) in d.items():
        if q in truth:
            per[ps].append((q, hap))
    out = {}
    for ps, members in per.items():
        c = collections.Counter((hap, truth[q]) for q, hap in members)
        straight = c[(1, 'MATERNAL')] + c[(2, 'PATERNAL')]
        flipped = c[(1, 'PATERNAL')] + c[(2, 'MATERNAL')]
        flip = flipped > straight
        for q, hap in members:
            expect = 'MATERNAL' if (hap == 1) != flip else 'PATERNAL'
            out[q] = truth[q] == expect
    return out


def blocks(vcf):
    spans = {}
    with vcf.open() as stream:
        for line in stream:
            if line.startswith('#'):
                continue
            c = line.rstrip('\n').split('\t')
            if len(c) < 10 or c[9].split(':')[0] not in ('0|1', '1|0'):
                continue
            keys = c[8].split(':')
            if 'PS' not in keys:
                continue
            ps = c[9].split(':')[keys.index('PS')].strip()
            pos = int(c[1])
            e = spans.setdefault(ps, [pos, pos, 0])
            e[0] = min(e[0], pos)
            e[1] = max(e[1], pos)
            e[2] += 1
    return spans


rows = []
with a.panel.open() as stream:
    panel = list(csv.DictReader(stream, delimiter='\t'))

for w in panel:
    gl, gr = int(w['gap_left']), int(w['gap_right'])
    lo, hi = gl - a.flank, gr + a.flank
    d = a.arm / f'w{gl}'
    if not (d / 'native.vcf').exists():
        print(f'  w{gl}: MISSING')
        continue
    spans = blocks(d / 'native.vcf')
    spanning = [e for e in spans.values() if e[0] <= gl and e[1] >= gr]
    in_gap = sum(1 for line in (d / 'native.vcf').read_text().splitlines()
                 if not line.startswith('#')
                 and gl <= int(line.split('\t')[1]) <= gr
                 and line.split('\t')[9].split(':')[0] in ('0|1', '1|0'))
    after = tags(d / 'phased.bam', lo, hi)
    ca = concordance(after)
    rec = dict(gap_left=gl, gap_right=gr, gap_bp=int(w['gap_bp']),
               competitor=w['competitor'], competitor_acc=w['competitor_acc'],
               blocks=len(spans), spans_gap=int(bool(spanning)),
               in_gap_hets=in_gap, tagged=len(after), scored=len(ca),
               concordant=sum(1 for v in ca.values() if v))
    rec['accuracy'] = round(rec['concordant'] / max(rec['scored'], 1), 5)
    if a.baseline is not None:
        bd = a.baseline / f'w{gl}'
        if (bd / 'phased.bam').exists():
            cb = concordance(tags(bd / 'phased.bam', lo, hi))
            x = collections.Counter()
            for q in set(cb) | set(ca):
                x[(cb.get(q), ca.get(q))] += 1
            rec.update(gate_conc_to_disc=x[(True, False)],
                       lost_concordant=x[(True, None)],
                       gained_concordant=x[(None, True)],
                       gained_discordant=x[(None, False)])
    rows.append(rec)

with a.output.open('w', newline='') as stream:
    w_ = csv.DictWriter(stream, fieldnames=list(rows[0]), delimiter='\t',
                        lineterminator='\n')
    w_.writeheader()
    w_.writerows(rows)

gated = 'gate_conc_to_disc' in rows[0]
hdr = '%-11s %8s %7s %6s %9s %8s %9s' % (
    'gap_left', 'gap_bp', 'blocks', 'spans', 'inGapHet', 'tagged', 'accuracy')
if gated:
    hdr += ' %9s %8s %14s' % ('c->d', 'lostC', 'new c/d')
print(hdr)
for r in rows:
    line = '%-11d %8d %7d %6s %9d %8d %8.2f%%' % (
        r['gap_left'], r['gap_bp'], r['blocks'], 'YES' if r['spans_gap'] else 'no',
        r['in_gap_hets'], r['tagged'], 100 * r['accuracy'])
    if gated:
        line += ' %9d %8d %7d/%-6d' % (r['gate_conc_to_disc'], r['lost_concordant'],
                                       r['gained_concordant'], r['gained_discordant'])
    print(line)
print()
print('panel totals: %d windows, %d spanned, %d in-gap hets, %d tagged, %.2f%% concordant' % (
    len(rows), sum(r['spans_gap'] for r in rows), sum(r['in_gap_hets'] for r in rows),
    sum(r['tagged'] for r in rows),
    100 * sum(r['concordant'] for r in rows) / max(sum(r['scored'] for r in rows), 1)))
if gated:
    print('panel gate: %d concordant->discordant, %d concordant tags lost, %d new concordant / %d new discordant' % (
        sum(r['gate_conc_to_disc'] for r in rows), sum(r['lost_concordant'] for r in rows),
        sum(r['gained_concordant'] for r in rows), sum(r['gained_discordant'] for r in rows)))
print(f'wrote {a.output}')
