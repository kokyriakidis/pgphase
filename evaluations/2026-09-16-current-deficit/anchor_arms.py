#!/usr/bin/env python3
"""Would admitting MSA-verified sites to the INITIAL partition connect a gap?

The hybrid solve already carries BAM variation -- the BAM chunk is the base and
graph sites are injected into it before one unified k-means. What differs between
the channels is which candidates may act as anchors: the initial k-means votes
only on the clean germline classes, and an MSA-verified noisy het indel is
admitted later, by a gap-recovery tier, if at all.

This tests the alternative offline, on one window, from the audit export. Two
arms: clean classes only (what the initial solve uses) and clean plus
MSA-verified. For each, sites are oriented and reads assigned by alternating
majority vote, blocks are defined by the site chain (two consecutive sites are
linked when enough reads observe both), and every block is scored against read
truth.

The clean-only arm is the control: it must reproduce the pipeline's own
fragmentation for the second arm to mean anything. An offline solver that
disagrees with the pipeline on the control is measuring itself.
"""
import argparse
import collections
import glob
from pathlib import Path

CLEAN_SNP, CLEAN_IND, REP, MSA = 0x004, 0x008, 0x010, 0x100

p = argparse.ArgumentParser(description=__doc__,
                            formatter_class=argparse.RawDescriptionHelpFormatter)
p.add_argument('--audit', type=Path, required=True)
p.add_argument('--gap-left', type=int, required=True)
p.add_argument('--gap-right', type=int, required=True)
p.add_argument('--target-left', type=int, required=True,
               help='interval the competitor spans, which is what must be covered')
p.add_argument('--target-right', type=int, required=True)
p.add_argument('--truth-map', type=Path, required=True)
p.add_argument('--min-obs', type=int, default=10)
p.add_argument('--min-link', type=int, default=2,
               help='reads observing both of two consecutive sites before they chain')
p.add_argument('--margin', type=int, default=2, help='read assignment margin')
a = p.parse_args()

truth = dict(line.rstrip('\n').split('\t') for line in a.truth_map.open())
hits = [f for f in glob.glob(str(a.audit / 'tid*.evidence.tsv'))
        if f'{a.gap_left}.{a.gap_right}' in f]
if not hits:
    raise SystemExit(f'no audit evidence for {a.gap_left}-{a.gap_right}')

sites, reads, obs = {}, {}, collections.defaultdict(dict)
for line in open(hits[0]):
    f = line.rstrip('\n').split('\t')
    if f[0] == 'SITE':
        sites[int(f[1])] = dict(pos=int(f[2]), cate=int(f[7]), msa=int(f[8]), hp=int(f[9]))
    elif f[0] == 'READ':
        reads[int(f[1])] = dict(q=f[3], beg=int(f[4]), end=int(f[5]), skip=int(f[7]))
    elif f[0] == 'OBS' and int(f[3]) >= 0:
        obs[int(f[2])][int(f[1])] = int(f[3])


def solve(site_ids):
    """Alternating orientation/assignment, then blocks from the site chain."""
    site_ids = sorted(site_ids, key=lambda v: sites[v]['pos'])
    orient = {}
    seed = max(site_ids, key=lambda v: len(obs[v]))
    orient[seed] = {al: (1 if i == 0 else 2)
                    for i, al in enumerate(sorted({v for v in obs[seed].values()}))}
    assign = {}
    for it in range(24):
        # The first passes must run at margin 1: with only the seed site
        # oriented, no read can reach a margin of 2, so the alternating loop
        # never leaves the seed and reports zero blocks -- a property of the
        # solver, not of the evidence.
        margin = 1 if it < 3 else a.margin
        new_assign = {}
        for ri in reads:
            if reads[ri]['skip']:
                continue
            tab = collections.Counter()
            for vi, mapping in orient.items():
                al = obs[vi].get(ri)
                if al is not None and al in mapping:
                    tab[mapping[al]] += 1
            if tab and abs(tab[1] - tab[2]) >= margin:
                new_assign[ri] = 1 if tab[1] > tab[2] else 2
        new_orient = dict(orient)
        for vi in site_ids:
            tab = collections.Counter()
            for ri, al in obs[vi].items():
                if ri in new_assign:
                    tab[(al, new_assign[ri])] += 1
            alleles = sorted({al for al, _ in tab})
            if len(alleles) != 2:
                continue
            mapping = {}
            for al in alleles:
                if tab[(al, 1)] == tab[(al, 2)]:
                    mapping = {}
                    break
                mapping[al] = 1 if tab[(al, 1)] > tab[(al, 2)] else 2
            if len(mapping) == 2 and mapping[alleles[0]] != mapping[alleles[1]]:
                new_orient[vi] = mapping
        if new_assign == assign and new_orient == orient:
            break
        assign, orient = new_assign, new_orient

    used = [v for v in site_ids if v in orient]
    chain = []
    cur = [used[0]] if used else []
    for left, right in zip(used, used[1:]):
        shared = sum(1 for ri in obs[left] if ri in obs[right])
        if shared >= a.min_link:
            cur.append(right)
        else:
            chain.append(cur)
            cur = [right]
    if cur:
        chain.append(cur)

    blocks = []
    for comp in chain:
        members = {}
        for vi in comp:
            for ri in obs[vi]:
                if ri in assign:
                    members[ri] = assign[ri]
        if not members:
            continue
        lo = min(sites[v]['pos'] for v in comp)
        hi = max(sites[v]['pos'] for v in comp)
        scored = [(h, truth[reads[r]['q']]) for r, h in members.items()
                  if reads[r]['q'] in truth]
        acc = None
        if scored:
            c = collections.Counter(scored)
            good = c[(1, 'MATERNAL')] + c[(2, 'PATERNAL')]
            bad = c[(1, 'PATERNAL')] + c[(2, 'MATERNAL')]
            acc = (max(good, bad) / len(scored), len(scored))
        blocks.append(dict(lo=lo, hi=hi, sites=len(comp), reads=len(members), acc=acc))
    return blocks, len(used)


interior = {}
for vi in sorted(sites, key=lambda v: (sites[v]['pos'], -len(obs[v]))):
    s = sites[vi]
    if a.gap_left <= s['pos'] <= a.gap_right and len(obs[vi]) >= a.min_obs:
        interior.setdefault(s['pos'], vi)
pool = list(interior.values())
clean = [v for v in pool if sites[v]['cate'] & (CLEAN_SNP | CLEAN_IND)]
verified = [v for v in pool if sites[v]['msa'] == 1
            and not sites[v]['cate'] & (CLEAN_SNP | CLEAN_IND)]

print(f'window {a.gap_left}-{a.gap_right}   competitor-spanned target '
      f'{a.target_left}-{a.target_right} ({(a.target_right - a.target_left) / 1e3:.1f} kb)')
print(f'interior sites with >= {a.min_obs} observations: {len(pool)}  '
      f'clean {len(clean)}  MSA-verified non-clean {len(verified)}')

for name, ids in (('clean only (control)', clean),
                  ('clean + MSA-verified', clean + verified)):
    blocks, used = solve(ids)
    print(f'\n{name}: {len(ids)} sites offered, {used} oriented, {len(blocks)} blocks')
    covering = 0
    for b in sorted(blocks, key=lambda b: b['lo']):
        shown = (f"{100 * b['acc'][0]:.1f}% (n={b['acc'][1]})" if b['acc'] else '--')
        covers = b['lo'] <= a.target_left and b['hi'] >= a.target_right
        covering += covers
        print(f"   {b['lo']}-{b['hi']} ({(b['hi'] - b['lo']) / 1e3:6.1f} kb) "
              f"sites={b['sites']:<3d} reads={b['reads']:<4d} acc={shown:<16s}"
              + ('   <-- covers the target' if covers else ''))
    print(f'   blocks covering the competitor-spanned target: {covering}')
