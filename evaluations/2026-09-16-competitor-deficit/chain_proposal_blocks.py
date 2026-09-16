#!/usr/bin/env python3
"""Can a fragmented gap proposal be closed by chaining its own blocks?

Gap recovery asks one question: does any single proposal block link BOTH flanks?
When the interior fragments -- which is what happens in wide gaps -- the answer is
no however accurate the fragments are, and the gap is abandoned. But the fragments
chain: reads cross the holes between consecutive blocks.

This scores that chain offline, from the audit export alone. For each block it
derives a per-site allele->haplotype consensus from that block's own reads, then
for each adjacent pair tabulates the reads that observe sites in both, exactly as
the pipeline's own shared-read orientation vote does, and composes the parity from
the left-flank block through to the right-flank block. Every accepted composition
is scored against read truth, so a chain that closes the gap incorrectly is
distinguishable from one that closes it correctly.
"""
import argparse
import collections
import glob
from pathlib import Path

p = argparse.ArgumentParser(description=__doc__,
                            formatter_class=argparse.RawDescriptionHelpFormatter)
p.add_argument('--audit', type=Path, required=True)
p.add_argument('--gap-left', type=int, required=True)
p.add_argument('--gap-right', type=int, required=True)
p.add_argument('--truth-map', type=Path, required=True)
p.add_argument('--tier', type=int, default=4)
p.add_argument('--view', type=int, default=0)
p.add_argument('--min-link-reads', type=int, default=2,
               help='minimum reads observing both blocks before a pair may link')
p.add_argument('--min-margin', type=int, default=2,
               help='net-margin threshold on the pair vote, as select_stitch_orientation uses')
a = p.parse_args()

truth = dict(line.rstrip('\n').split('\t') for line in a.truth_map.open())
GL, GR = a.gap_left, a.gap_right
hits = [f for f in glob.glob(str(a.audit / 'tid*.evidence.tsv')) if f'{GL}.{GR}' in f]
if not hits:
    raise SystemExit(f'no audit evidence for {GL}-{GR}')
ev = hits[0]

sites, reads, obs = {}, {}, collections.defaultdict(dict)
for line in open(ev):
    f = line.rstrip('\n').split('\t')
    if f[0] == 'SITE':
        sites[int(f[1])] = dict(pos=int(f[2]))
    elif f[0] == 'READ':
        reads[int(f[1])] = dict(q=f[3], beg=int(f[4]), end=int(f[5]), skip=int(f[7]))
    elif f[0] == 'OBS' and int(f[3]) >= 0:
        obs[int(f[2])][int(f[1])] = int(f[3])

part = collections.defaultdict(dict)
for line in open(ev.replace('.evidence.tsv', '.reads.tsv')):
    f = line.rstrip('\n').split('\t')
    if f[0] == 'VIEW':
        continue
    view, tier, ri, hp, ps, sk = (int(x) for x in f[:6])
    if view != a.view or tier != a.tier or sk or hp == 0 or ps < 0:
        continue
    part[ps][ri] = hp
if not part:
    raise SystemExit(f'no blocks for view {a.view} tier {a.tier}')

# Per-block, per-site allele -> haplotype consensus from that block's own reads.
consensus = {}
for ps, members in part.items():
    per_site = {}
    for vi, per_read in obs.items():
        tab = collections.Counter()
        for ri, allele in per_read.items():
            if ri in members:
                tab[(allele, members[ri])] += 1
        alleles = sorted({al for al, _ in tab})
        if len(alleles) != 2:
            continue
        mapping = {}
        for al in alleles:
            h1, h2 = tab[(al, 1)], tab[(al, 2)]
            if h1 == h2:
                mapping = {}
                break
            mapping[al] = 1 if h1 > h2 else 2
        if len(mapping) == 2 and mapping[alleles[0]] != mapping[alleles[1]]:
            per_site[vi] = mapping
    consensus[ps] = per_site

def implied(ri, ps):
    """Haplotype this read implies under block `ps`, and how many sites voted."""
    tab = collections.Counter()
    for vi, mapping in consensus[ps].items():
        allele = obs[vi].get(ri)
        if allele is None or allele not in mapping:
            continue
        tab[mapping[allele]] += 1
    if not tab:
        return 0, 0
    total = sum(tab.values())
    if tab[1] == tab[2]:
        return 0, total
    return (1 if tab[1] > tab[2] else 2), total

extent = {ps: (min(reads[r]['beg'] for r in m), max(reads[r]['end'] for r in m))
          for ps, m in part.items()}
order = sorted(part, key=lambda ps: extent[ps][0])
print(f'gap {GL}-{GR} ({(GR - GL) / 1e3:.1f} kb)  view {a.view} tier {a.tier}: '
      f'{len(order)} blocks')
for ps in order:
    print(f'   ps={ps:<11d} reads={len(part[ps]):<4d} span {extent[ps][0]}-{extent[ps][1]} '
          f'consensus sites={len(consensus[ps])}')

print('\nadjacent-pair link votes (reads observing sites in both blocks):')
edges = {}
for left, right in zip(order, order[1:]):
    tab = collections.Counter()
    for ri in reads:
        if reads[ri]['skip']:
            continue
        hl, nl = implied(ri, left)
        hr, nr = implied(ri, right)
        if hl == 0 or hr == 0 or nl < 1 or nr < 1:
            continue
        tab[(hl, hr)] += 1
    straight = tab[(1, 1)] + tab[(2, 2)]
    flipped = tab[(1, 2)] + tab[(2, 1)]
    n = straight + flipped
    margin = abs(straight - flipped)
    ok = n >= a.min_link_reads and margin > a.min_margin
    flip = flipped > straight
    verdict = ('LINK flip' if flip else 'LINK same') if ok else 'no link'
    print(f'   {left} -> {right}: n={n:<4d} straight={straight:<4d} flipped={flipped:<4d} '
          f'margin={margin:<4d} {verdict}')
    if ok:
        edges[(left, right)] = flip

print('\nchain composition:')
parity = {order[0]: False}
for left, right in zip(order, order[1:]):
    if (left, right) not in edges or left not in parity:
        break
    parity[right] = parity[left] ^ edges[(left, right)]
reached = [ps for ps in order if ps in parity]
left_edge = [ps for ps in order if any(reads[r]['beg'] <= GL <= reads[r]['end'] for r in part[ps])]
right_edge = [ps for ps in order if any(reads[r]['beg'] <= GR <= reads[r]['end'] for r in part[ps])]
print(f'   blocks in one composed frame: {len(reached)} of {len(order)}')
print(f'   left-edge block(s) {left_edge}  right-edge block(s) {right_edge}')
closes = bool(left_edge) and bool(right_edge) and left_edge[0] in parity and right_edge[-1] in parity
print(f'   composed frame spans both gap edges: {"YES" if closes else "no"}')

scored = []
for ps in reached:
    for ri, hp in part[ps].items():
        q = reads[ri]['q']
        if q not in truth:
            continue
        scored.append((3 - hp if parity[ps] else hp, truth[q]))
if scored:
    c = collections.Counter(scored)
    good = c[(1, 'MATERNAL')] + c[(2, 'PATERNAL')]
    bad = c[(1, 'PATERNAL')] + c[(2, 'MATERNAL')]
    print(f'   composed frame scored against truth: {max(good, bad)}/{good + bad} '
          f'({100 * max(good, bad) / (good + bad):.2f}%)')
