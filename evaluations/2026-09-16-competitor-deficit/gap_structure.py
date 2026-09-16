#!/usr/bin/env python3
"""Why does a gap not join? Print the structure that decides it.

Reads a `--gap-decision-audit` export plus the run's outputs and reports, in the
order the join actually depends on:

  1. the proposal's blocks per pass/tier -- span, reads, truth accuracy, and
     whether each reaches the gap's left or right edge (a join needs ONE block
     holding both, so a fragmented proposal cannot join however good it is);
  2. the interior site inventory by class, and the site spacings that break the
     read chain (zero reads covering both flanking sites is a hard linkage
     break: no read-based phaser crosses it and abstaining is correct);
  3. what survives into the output, so interior blocks the pipeline computed and
     then discarded are visible.

Truth accuracy is per block under its own best orientation, from a read-name to
haplotype map.
"""
import argparse
import collections
import csv
import glob
import subprocess
from pathlib import Path

CLEAN_SNP, CLEAN_IND, REP, MSA = 0x004, 0x008, 0x010, 0x100

p = argparse.ArgumentParser(description=__doc__,
                            formatter_class=argparse.RawDescriptionHelpFormatter)
p.add_argument('--audit', type=Path, required=True, help='--gap-decision-audit directory')
p.add_argument('--gap-left', type=int, required=True)
p.add_argument('--gap-right', type=int, required=True)
p.add_argument('--truth-map', type=Path, required=True, help='read name -> MATERNAL/PATERNAL')
p.add_argument('--bam', type=Path, help='surjected input BAM, for raw coverage at a break')
p.add_argument('--contig', default='CHM13#0#chr20')
p.add_argument('--phased-bam', type=Path, help='the run\'s output BAM')
p.add_argument('--tiers', help='only these tiers (default 1,4)', default='1,4')
a = p.parse_args()

truth = dict(line.rstrip('\n').split('\t') for line in a.truth_map.open())
GL, GR = a.gap_left, a.gap_right
hits = [f for f in glob.glob(str(a.audit / 'tid*.evidence.tsv')) if f'{GL}.{GR}' in f]
if not hits:
    raise SystemExit(f'no audit evidence for {GL}-{GR} in {a.audit}')
ev = hits[0]

sites, reads, obs = {}, {}, collections.defaultdict(dict)
for line in open(ev):
    f = line.rstrip('\n').split('\t')
    if f[0] == 'SITE':
        sites[int(f[1])] = dict(pos=int(f[2]), vtype=int(f[3]), cate=int(f[7]),
                                msa=int(f[8]), hp=int(f[9]))
    elif f[0] == 'READ':
        reads[int(f[1])] = dict(q=f[3], beg=int(f[4]), end=int(f[5]), skip=int(f[7]))
    elif f[0] == 'OBS' and int(f[3]) >= 0:
        obs[int(f[2])][int(f[1])] = int(f[3])

blocks = collections.defaultdict(lambda: collections.defaultdict(dict))
for line in open(ev.replace('.evidence.tsv', '.reads.tsv')):
    f = line.rstrip('\n').split('\t')
    if f[0] == 'VIEW':
        continue
    view, tier, ri, hp, ps, sk = (int(x) for x in f[:6])
    if sk or hp == 0 or ps < 0:
        continue
    blocks[(view, tier)][ps][ri] = hp


def accuracy(part):
    sel = [(r, h) for r, h in part.items() if reads[r]['q'] in truth]
    if not sel:
        return None
    c = collections.Counter((h, truth[reads[r]['q']]) for r, h in sel)
    good = c[(1, 'MATERNAL')] + c[(2, 'PATERNAL')]
    bad = c[(1, 'PATERNAL')] + c[(2, 'MATERNAL')]
    return max(good, bad) / (good + bad), good + bad


want = {int(t) for t in a.tiers.split(',')}
print(f'gap {GL}-{GR}  ({(GR - GL) / 1e3:.1f} kb)   audit reads {len(reads)}')
for key in sorted(blocks):
    if key[1] not in want:
        continue
    bl = blocks[key]
    print(f'\npass {key[0]} tier {key[1]}: {len(bl)} proposal blocks')
    rows = []
    for ps, part in bl.items():
        rr = [reads[r] for r in part]
        acc = accuracy(part)
        rows.append((min(r['beg'] for r in rr), max(r['end'] for r in rr), ps, len(part), acc))
    for lo, hi, ps, n, acc in sorted(rows):
        touch_l = sum(1 for r in bl[ps] if reads[r]['beg'] <= GL <= reads[r]['end'])
        touch_r = sum(1 for r in bl[ps] if reads[r]['beg'] <= GR <= reads[r]['end'])
        shown = f'{100 * acc[0]:.1f}% (n={acc[1]})' if acc else '--'
        print(f'   ps={ps:<11d} reads={n:<4d} span {lo}-{hi} ({(hi - lo) / 1e3:6.1f} kb) '
              f'acc={shown:<16s} leftEdge={touch_l:<3d} rightEdge={touch_r}')
    spanning = [ps for ps in bl
                if any(reads[r]['beg'] <= GL <= reads[r]['end'] for r in bl[ps])
                and any(reads[r]['beg'] <= GR <= reads[r]['end'] for r in bl[ps])]
    print(f'   blocks reaching BOTH edges (a join needs one): {len(spanning)}')

live = [r for r in reads.values() if not r['skip']]
interior = {}
for vi in sorted((v for v, s in sites.items()
                  if GL <= s['pos'] <= GR
                  and s['cate'] & (CLEAN_SNP | CLEAN_IND | REP | MSA)
                  and len(obs[v]) >= 10),
                 key=lambda v: (sites[v]['pos'], -len(obs[v]))):
    interior.setdefault(sites[vi]['pos'], vi)
order = sorted(interior.values(), key=lambda v: sites[v]['pos'])


def klass(s):
    if s['cate'] & CLEAN_SNP:
        return 'CLEAN_SNP'
    if s['cate'] & CLEAN_IND:
        return 'CLEAN_IND'
    if s['cate'] & MSA:
        return 'MSA'
    return 'REP'


print(f'\ninterior sites with >=10 observations: {len(order)}  '
      f'{dict(collections.Counter(klass(sites[v]) for v in order))}')
pos = [sites[v]['pos'] for v in order]
breaks = []
if len(pos) > 1:
    spacing = sorted(((pos[i + 1] - pos[i], pos[i], pos[i + 1]) for i in range(len(pos) - 1)),
                     reverse=True)
    print('   largest spacings, with reads covering both flanking sites:')
    for d, lo, hi in spacing[:6]:
        n = sum(1 for r in live if r['beg'] <= lo and r['end'] >= hi)
        print(f'      {d / 1e3:8.1f} kb  {lo} -> {hi}   spanning reads: {n}'
              + ('   <-- LINKAGE BREAK' if n == 0 else ''))
        if n == 0:
            breaks.append((lo, hi))
    print(f'   site range {pos[0]}-{pos[-1]} covers {(pos[-1] - pos[0]) / 1e3:.1f} kb '
          f'of the {(GR - GL) / 1e3:.1f} kb gap')

if breaks and a.bam:
    print('\n   raw coverage at each linkage break (breaks are not coverage holes):')
    for lo, hi in breaks:
        out = subprocess.run(['samtools', 'depth', '-a', '-Q', '1', '-r',
                              f'{a.contig}:{lo}-{hi}', str(a.bam)],
                             capture_output=True, text=True, check=True).stdout
        d = [int(l.split('\t')[2]) for l in out.splitlines()]
        reads_out = subprocess.run(['samtools', 'view', '-q', '1', str(a.bam),
                                    f'{a.contig}:{lo}-{hi}'],
                                   capture_output=True, text=True, check=True).stdout.splitlines()
        longest = max((len(l.split('\t')[9]) for l in reads_out), default=0)
        print(f'      {lo}-{hi} ({(hi - lo) / 1e3:.1f} kb): mean depth '
              f'{sum(d) / max(len(d), 1):.1f}, {len(reads_out)} reads, longest {longest / 1e3:.1f} kb')

if a.phased_bam:
    out = subprocess.run(['samtools', 'view', str(a.phased_bam)],
                         capture_output=True, text=True, check=True).stdout
    got = collections.defaultdict(list)
    for line in out.splitlines():
        f = line.split('\t')
        hp = ps = None
        for x in f[11:]:
            if x.startswith('HP:i:'):
                hp = x[5:]
            elif x.startswith('PS:i:'):
                ps = x[5:]
        if hp and ps:
            got[ps].append((int(hp), f[0], int(f[3])))
    print('\noutput blocks (>=5 reads):')
    for ps, v in sorted(got.items(), key=lambda kv: min(z[2] for z in kv[1])):
        if len(v) < 5:
            continue
        sel = [(h, q) for h, q, _ in v if q in truth]
        c = collections.Counter((h, truth[q]) for h, q in sel)
        good = c[(1, 'MATERNAL')] + c[(2, 'PATERNAL')]
        bad = c[(1, 'PATERNAL')] + c[(2, 'MATERNAL')]
        inside = 'inside gap' if GL <= min(z[2] for z in v) <= GR else ''
        print(f'   ps={ps:<11s} reads={len(v):<4d} span {min(z[2] for z in v)}-'
              f'{max(z[2] for z in v)} acc={100 * max(good, bad) / max(good + bad, 1):.1f}% {inside}')
