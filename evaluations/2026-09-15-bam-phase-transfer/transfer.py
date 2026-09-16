#!/usr/bin/env python3
"""Transfer BAM-channel phasing into the graph's gauge and try to stitch the gap.

The two channels phase the same reads, so read identity is an exact bridge: a
read that carries both a BAM haplotype and a graph phase set says how the two
gauges line up. For each gap this asks whether one BAM phase set reaches both
adjacent graph blocks, whether the orientation it implies is consistent on each
side, and how many reads the graph left unphased would inherit a haplotype. The
transferred assignments are then scored against the read-level truth.

Refusing is a valid outcome and is reported as such: without a consistent anchor
on both flanks the gap must stay open rather than be merged on one side's word.
"""
import argparse
import collections
import csv
import subprocess
from pathlib import Path

p = argparse.ArgumentParser(description=__doc__)
p.add_argument('--transfer-root', type=Path, required=True, help='per-window BAM phasing output')
p.add_argument('--graph-reads', type=Path, required=True, help='graph-only --phase-reads-out TSV')
p.add_argument('--gaps', type=Path, required=True, help='graph phase-block gap inventory')
p.add_argument('--truth-map', type=Path, required=True)
p.add_argument('--min-anchor-reads', type=int, default=5)
p.add_argument('--min-anchor-agreement', type=float, default=0.90)
p.add_argument('--output', type=Path, required=True)
a = p.parse_args()

truth = {}
with a.truth_map.open() as stream:
    for line in stream:
        name, hap = line.rstrip('\n').split('\t')
        truth[name] = hap

graph = {}
with a.graph_reads.open() as stream:
    for r in csv.DictReader(stream, delimiter='\t'):
        if r['HAP'] in ('1', '2') and r['PHASE_SET'] not in ('', '0', '-1', '.'):
            graph[r['READ']] = (r['HAP'], r['PHASE_SET'])

gaps = {}
with a.gaps.open() as stream:
    for r in csv.DictReader(stream, delimiter='\t'):
        gaps[int(r['gap_left'])] = r


positions = {}


def bam_tags(path):
    """read name -> (HP, PS) from a phased BAM, streamed (no index needed)."""
    done = subprocess.run(['samtools', 'view', str(path)], capture_output=True, text=True)
    if done.returncode != 0:
        raise SystemExit('samtools view failed on %s: %s' % (path, done.stderr.strip()[:200]))
    out = {}
    for line in done.stdout.splitlines():
        fields = line.split('\t')
        hp = ps = None
        for x in fields[11:]:
            if x.startswith('HP:i:'):
                hp = x[5:]
            elif x.startswith('PS:i:'):
                ps = x[5:]
        if hp is not None and ps is not None:
            out[fields[0]] = (hp, ps)
            positions[fields[0]] = int(fields[3])
    return out


def orientation(shared):
    """Majority mapping from BAM HP to graph HAP, with its agreement fraction."""
    counts = collections.Counter((hp, hap) for hp, hap in shared)
    straight = counts[('1', '1')] + counts[('2', '2')]
    flipped = counts[('1', '2')] + counts[('2', '1')]
    total = straight + flipped
    if total == 0:
        return None, 0.0, 0
    if straight >= flipped:
        return 'straight', straight / total, total
    return 'flipped', flipped / total, total


rows = []
for window_dir in sorted(a.transfer_root.iterdir()):
    if not (window_dir / 'phased.bam').exists():
        continue
    name = window_dir.name
    gap_left = int(name.split('_')[1]) if name.startswith('gap_') else None
    gap = gaps.get(gap_left)
    if gap is None:
        continue
    left_ps, right_ps = gap['left_ps'], gap['right_ps']
    tags = bam_tags(window_dir / 'phased.bam')

    # Which BAM phase set reaches both graph blocks?
    reach = collections.defaultdict(lambda: {'left': [], 'right': []})
    for read, (hp, ps) in tags.items():
        g = graph.get(read)
        if g is None:
            continue
        hap, gps = g
        if gps == left_ps:
            reach[ps]['left'].append((hp, hap))
        elif gps == right_ps:
            reach[ps]['right'].append((hp, hap))

    # Anchor each side independently. A single BAM phase set spanning both graph
    # blocks would stitch the gap outright; when none does, each side's anchored
    # phase set can still be transferred inward, and what remains between them is
    # the residual break -- the only interval that still needs new evidence.
    def pick(side, other):
        best, best_n, best_agree = None, 0, 0.0
        for ps, sides in reach.items():
            direction, agree, n = orientation(sides[side])
            if n >= a.min_anchor_reads and agree >= a.min_anchor_agreement and n > best_n:
                best, best_n, best_agree = (ps, direction), n, agree
        return best, best_n, best_agree

    left_pick, left_n, left_agree = pick('left', 'right')
    right_pick, right_n, right_agree = pick('right', 'left')
    spanning = [ps for ps, sides in reach.items()
                if orientation(sides['left'])[2] >= a.min_anchor_reads
                and orientation(sides['right'])[2] >= a.min_anchor_reads]

    # A phase set's span is first to last phased variant, taken from the phased
    # VCF. Read start positions understate it by up to a read length at each end.
    extent = {}
    vcf = window_dir / 'native.vcf'
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
            ps = values[keys.index('PS')]
            pos = int(f[1])
            lo_hi = extent.setdefault(ps, [10 ** 12, 0])
            lo_hi[0] = min(lo_hi[0], pos)
            lo_hi[1] = max(lo_hi[1], pos)

    transferred, concordant, scored = 0, 0, 0
    for pick_, side_ps in ((left_pick, left_ps), (right_pick, right_ps)):
        if pick_ is None:
            continue
        ps, direction = pick_
        expect = {'1': '1', '2': '2'} if direction == 'straight' else {'1': '2', '2': '1'}
        gained = [(read, hp) for read, (hp, q) in tags.items()
                  if q == ps and read not in graph]
        transferred += len(gained)
        # Truth check: the graph block's own orientation to truth is read off the
        # reads it already phased, so the transfer is scored in the graph's gauge.
        block_map = collections.Counter()
        for read, (hap, gps) in graph.items():
            if gps == side_ps and read in truth:
                block_map[(hap, truth[read])] += 1
        straight = block_map[('1', 'MATERNAL')] + block_map[('2', 'PATERNAL')]
        flipped = block_map[('1', 'PATERNAL')] + block_map[('2', 'MATERNAL')]
        hap_to_truth = {'1': 'MATERNAL', '2': 'PATERNAL'} if straight >= flipped else \
                       {'1': 'PATERNAL', '2': 'MATERNAL'}
        for read, hp in gained:
            if read not in truth:
                continue
            scored += 1
            if truth[read] == hap_to_truth[expect[hp]]:
                concordant += 1

    # The residual break is an interval, not just a length: it runs from where the
    # left-anchored phase set stops to where the right-anchored one starts. Its
    # endpoints are reported so downstream work targets the actual interval
    # rather than assuming it sits at the gap midpoint.
    residual = 0
    residual_beg = residual_end = 0
    if left_pick and right_pick:
        left_end = extent.get(left_pick[0], [0, 0])[1]
        right_beg = extent.get(right_pick[0], [0, 0])[0]
        residual = max(0, right_beg - left_end)
        if residual:
            residual_beg, residual_end = left_end, right_beg
    if spanning:
        verdict = 'one_bam_ps_spans_gap'
    elif left_pick and right_pick:
        verdict = 'both_sides_transfer'
    elif left_pick or right_pick:
        verdict = 'one_side_transfer'
    else:
        verdict = 'no_anchor'

    rows.append({
        'gap_left': gap_left, 'gap_bp': int(gap['gap_bp']), 'verdict': verdict,
        'bam_phase_sets': len(set(ps for _, ps in tags.values())),
        'anchor_left': left_n, 'anchor_right': right_n,
        'agree_left': round(left_agree, 4), 'agree_right': round(right_agree, 4),
        'residual_break_bp': residual,
        'residual_beg': residual_beg, 'residual_end': residual_end,
        'reads_transferred': transferred, 'transferred_concordant': concordant,
        'transfer_accuracy': round(concordant / scored, 4) if scored else 0.0})

with a.output.open('w', newline='') as stream:
    w = csv.DictWriter(stream, fieldnames=list(rows[0]), delimiter='\t')
    w.writeheader()
    w.writerows(rows)

print('%-11s %9s %-21s %6s %7s %7s %7s %10s %8s %9s' % (
    'gap_left', 'gap_bp', 'verdict', 'bamPS', 'anchorL', 'anchorR', 'agree',
    'residual', 'gained', 'accuracy'))
for r in rows:
    print('%-11s %9d %-21s %6d %7d %7d %7.3f %10d %8d %8.1f%%' % (
        r['gap_left'], r['gap_bp'], r['verdict'], r['bam_phase_sets'],
        r['anchor_left'], r['anchor_right'], min(r['agree_left'], r['agree_right']),
        r['residual_break_bp'], r['reads_transferred'], 100 * r['transfer_accuracy']))
print()
for verdict, n in collections.Counter(r['verdict'] for r in rows).most_common():
    print('  %-22s %d gaps' % (verdict, n))
tot = sum(r['reads_transferred'] for r in rows)
con = sum(r['transferred_concordant'] for r in rows)
print('  reads transferred: %d, concordant against truth %d (%.2f%%)' % (
    tot, con, 100 * con / tot if tot else 0.0))
print('  gap span before: %.1f kb   residual break after transfer: %.1f kb (%.1f%%)' % (
    sum(r['gap_bp'] for r in rows) / 1e3, sum(r['residual_break_bp'] for r in rows) / 1e3,
    100 * sum(r['residual_break_bp'] for r in rows) / sum(r['gap_bp'] for r in rows)))
