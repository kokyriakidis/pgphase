#!/usr/bin/env python3
"""Is there a chain of sites that phases a panel gap, and is it findable?

Takes the union of every site available in a window -- the pipeline's own
candidates (alignment-derived plus injected catalog sites) and every graph
catalog site, most of which never become candidates -- genotypes each one
directly from the alignment, and searches for a path of sites linking the left
flank to the right flank.

Two rules make the search answer a useful question rather than an oracle one.
Sites are genotyped and partitioned WITHOUT truth: a substitution splits on the
base, an indel splits on the net length change across its tract, which is what
separates the haplotypes in a repeat where matching the emitted allele exactly
does not. And the path is chosen WITHOUT truth, maximising the weakest link's
shared-read count. Truth is used only afterwards, to say whether the chain the
search found is the correct one -- so a chain reported here is one the pipeline
could in principle find, not one only hindsight can see.
"""
import argparse
import bisect
import collections
import csv
import gzip
import re
import subprocess
from pathlib import Path

CONTIG = 'CHM13#0#chr20'
CIGAR = re.compile(r'(\d+)([MIDNSHP=X])')
TRACT_FLANK = 25
MIN_GROUP = 5          # reads needed on each side for a site to define a partition
MIN_SHARED = 5         # reads needed at both ends for an edge (overridden by --min-shared)
MIN_CONSISTENCY = 0.90  # of the shared reads, the fraction the orientation must explain
MAX_EDGE_BP = 40000    # longest link considered


class Read:
    __slots__ = ('name', 'beg', 'end', 'blocks', 'indels', 'seq', '_starts', '_ipos')

    def __init__(self, name, beg, cigar, seq):
        self.name = name
        self.beg = beg
        self.seq = seq
        self.blocks = []   # (ref_beg, ref_end, query_beg) for aligned stretches
        self.indels = []   # (ref_pos, +len for insertion / -len for deletion)
        ref = beg
        qry = 0
        for n, op in CIGAR.findall(cigar):
            n = int(n)
            if op in 'M=X':
                self.blocks.append((ref, ref + n, qry))
                ref += n
                qry += n
            elif op == 'I':
                self.indels.append((ref, n))
                qry += n
            elif op in 'DN':
                self.indels.append((ref, -n))
                ref += n
            elif op == 'S':
                qry += n
        self.end = ref
        self._starts = [b[0] for b in self.blocks]
        self._ipos = [i[0] for i in self.indels]

    def base_at(self, pos):
        i = bisect.bisect_right(self._starts, pos) - 1
        if i < 0:
            return None
        rb, re_, qb = self.blocks[i]
        if not (rb <= pos < re_):
            return None
        return self.seq[qb + (pos - rb)]

    def net_length(self, lo, hi):
        total = 0
        i = bisect.bisect_left(self._ipos, lo - 1)
        while i < len(self.indels) and self.indels[i][0] <= hi:
            pos, delta = self.indels[i]
            if delta > 0:
                if lo <= pos <= hi:
                    total += delta
            else:
                if pos + (-delta) >= lo and pos <= hi:
                    total += delta
            i += 1
        return total


def load_reads(bam, lo, hi, min_mapq):
    out = subprocess.run(['samtools', 'view', '-q', str(min_mapq), bam,
                          f'{CONTIG}:{lo}-{hi}'], capture_output=True, text=True)
    if out.returncode != 0:
        raise RuntimeError(out.stderr[:300])
    reads = []
    for line in out.stdout.splitlines():
        f = line.split('\t')
        if f[5] == '*':
            continue
        reads.append(Read(f[0], int(f[3]), f[5], f[9]))
    return reads


def load_sites(candidates, catalog, lo, hi):
    """Union of pipeline candidates and catalog sites, with provenance."""
    cand = {}
    with open(candidates) as stream:
        for r in csv.DictReader(stream, delimiter='\t'):
            pos = int(r['POS'])
            if not (lo <= pos <= hi):
                continue
            cand[(pos, r['REF'], r['ALT'])] = r['CATEGORY']
    cat = set()
    with gzip.open(catalog, 'rt') as stream:
        for line in stream:
            if line.startswith('#'):
                continue
            c = line.split('\t', 6)
            pos = int(c[1])
            if lo <= pos <= hi:
                for alt in c[4].split(','):
                    cat.add((pos, c[3], alt))
    sites = []
    for key in set(cand) | cat:
        pos, ref, alt = key
        in_c, in_g = key in cand, key in cat
        sites.append(dict(pos=pos, ref=ref, alt=alt,
                          provenance='both' if (in_c and in_g) else
                                     ('candidate' if in_c else 'catalog'),
                          category=cand.get(key, '-')))
    sites.sort(key=lambda s: (s['pos'], s['ref'], s['alt']))
    return sites


def partition_site(site, reads):
    """Label reads 0/1 at this site without using truth. Returns {read: label}."""
    pos, ref, alt = site['pos'], site['ref'], site['alt']
    labels = {}
    if len(ref) == 1 and len(alt) == 1 and alt != '.':
        for rd in reads:
            if rd.beg > pos or rd.end <= pos:
                continue
            b = rd.base_at(pos)
            if b is None:
                continue
            if b.upper() == ref.upper():
                labels[rd.name] = 0
            elif b.upper() == alt.upper():
                labels[rd.name] = 1
        return labels, 'base'

    # Indel: split on net length across the tract. The two haplotypes of a
    # repeat differ by net length even when neither matches the emitted allele.
    lo = pos - TRACT_FLANK
    hi = pos + max(len(ref), 1) + TRACT_FLANK
    nets = {}
    for rd in reads:
        if rd.beg > lo or rd.end < hi:
            continue
        nets[rd.name] = rd.net_length(lo, hi)
    if len(nets) < 2 * MIN_GROUP:
        return {}, 'net'
    values = sorted(nets.values())
    # Widest gap between consecutive observed lengths, with both sides populated.
    best = None
    for i in range(MIN_GROUP - 1, len(values) - MIN_GROUP):
        if values[i + 1] == values[i]:
            continue
        width = values[i + 1] - values[i]
        if best is None or width > best[0]:
            best = (width, (values[i] + values[i + 1]) / 2.0)
    if best is None or best[0] < 1:
        return {}, 'net'
    cut = best[1]
    for name, v in nets.items():
        labels[name] = 0 if v < cut else 1
    return labels, 'net'


def orient(a_labels, b_labels):
    """Shared reads and the fraction the better orientation explains."""
    shared = set(a_labels) & set(b_labels)
    if len(shared) < MIN_SHARED:
        return len(shared), 0.0, None
    same = sum(1 for q in shared if a_labels[q] == b_labels[q])
    diff = len(shared) - same
    keep = max(same, diff)
    return len(shared), keep / len(shared), (0 if same >= diff else 1)


def best_chain(nodes, edges, left_anchor, right_anchor):
    """Max-bottleneck path: maximise the weakest edge's shared-read count."""
    parent = {n: n for n in nodes}

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    order = sorted(edges, key=lambda e: -e['shared'])
    kept = []
    for e in order:
        kept.append(e)
        ra, rb = find(e['a']), find(e['b'])
        if ra != rb:
            parent[ra] = rb
        if find(left_anchor) == find(right_anchor):
            bottleneck = e['shared']
            adj = collections.defaultdict(list)
            for k in kept:
                adj[k['a']].append(k)
                adj[k['b']].append(k)
            # Shortest path in edge count among edges at or above the bottleneck.
            prev = {left_anchor: None}
            queue = collections.deque([left_anchor])
            while queue:
                cur = queue.popleft()
                if cur == right_anchor:
                    break
                for k in adj[cur]:
                    nxt = k['b'] if k['a'] == cur else k['a']
                    if nxt not in prev:
                        prev[nxt] = (cur, k)
                        queue.append(nxt)
            if right_anchor not in prev:
                continue
            path = []
            cur = right_anchor
            while prev[cur] is not None:
                p, k = prev[cur]
                path.append(k)
                cur = p
            return bottleneck, list(reversed(path))
    return 0, []


def main():
    global MIN_SHARED
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--panel', type=Path, required=True)
    p.add_argument('--arm', type=Path, required=True,
                   help='run_panel.sh OUT dir, for each window candidates.tsv')
    p.add_argument('--catalog', type=Path, required=True)
    p.add_argument('--bam', type=Path, required=True)
    p.add_argument('--truth-map', type=Path, required=True)
    p.add_argument('--min-mapq', type=int, default=1)
    p.add_argument('--flank', type=int, default=50000)
    p.add_argument('--outdir', type=Path, required=True)
    p.add_argument('--min-shared', type=int, default=MIN_SHARED,
                   help='reads needed at both ends of a link; hiphase links some of '
                        'these breaks on 2, so the default of 5 can refuse a link it '
                        'would make')
    p.add_argument('--corridor', type=int, default=5000,
                   help='sites are restricted to [gap_left - corridor, gap_right + '
                        'corridor] so the chain must cross the gap rather than '
                        'detour through the flanking blocks')
    p.add_argument('--tag', default='')
    a = p.parse_args()
    MIN_SHARED = a.min_shared
    a.outdir.mkdir(parents=True, exist_ok=True)
    truth = dict(line.split('\t') for line in a.truth_map.read_text().splitlines())

    summary = []
    with a.panel.open() as stream:
        panel = list(csv.DictReader(stream, delimiter='\t'))

    for w in panel:
        gl, gr = int(w['gap_left']), int(w['gap_right'])
        lo, hi = gl - a.flank, gr + a.flank
        reads = load_reads(str(a.bam), lo, hi, a.min_mapq)
        sites = load_sites(a.arm / f'w{gl}' / 'candidates.tsv', a.catalog, lo, hi)

        nodes = {}
        corridor = (gl - a.corridor, gr + a.corridor)
        for s in sites:
            if not (corridor[0] <= s['pos'] <= corridor[1]):
                continue
            labels, how = partition_site(s, reads)
            groups = collections.Counter(labels.values())
            if len(groups) < 2 or min(groups.values()) < MIN_GROUP:
                continue
            key = (s['pos'], s['ref'], s['alt'])
            nodes[key] = dict(site=s, labels=labels, how=how)

        keys = sorted(nodes)
        edges = []
        for i, ka in enumerate(keys):
            for kb in keys[i + 1:]:
                if kb[0] - ka[0] > MAX_EDGE_BP:
                    break
                shared, cons, flip = orient(nodes[ka]['labels'], nodes[kb]['labels'])
                if flip is None or cons < MIN_CONSISTENCY:
                    continue
                edges.append(dict(a=ka, b=kb, shared=shared, cons=cons, flip=flip))

        left = [k for k in keys if k[0] <= gl]
        right = [k for k in keys if k[0] >= gr]
        if not left or not right or not edges:
            summary.append(dict(gap_left=gl, gap_bp=int(w['gap_bp']), usable_sites=len(nodes),
                                edges=len(edges), chain_found=0, bottleneck=0, chain_sites=0,
                                in_gap_chain_sites=0, from_catalog_only=0, min_segregation=0.0,
                                orientation='no chain'))
            continue
        la, ra_ = left[-1], right[0]
        bottleneck, path = best_chain(set(keys), edges, la, ra_)

        chain_keys = []
        if path:
            chain_keys = [path[0]['a']]
            for e in path:
                chain_keys.append(e['b'] if e['a'] == chain_keys[-1] else e['a'])

        # Truth check, after the fact: propagate one orientation and score.
        seg = []
        rows = []
        flip_state = 0
        for i, k in enumerate(chain_keys):
            if i > 0:
                flip_state ^= path[i - 1]['flip']
            labels = nodes[k]['labels']
            c = collections.Counter()
            for q, lab in labels.items():
                if q in truth:
                    c[((lab ^ flip_state), truth[q])] += 1
            straight = c[(0, 'MATERNAL')] + c[(1, 'PATERNAL')]
            crossed = c[(0, 'PATERNAL')] + c[(1, 'MATERNAL')]
            n = straight + crossed
            s_ = max(straight, crossed) / n if n else 0.0
            seg.append(s_)
            rows.append(dict(pos=k[0], ref=k[1][:20], alt=k[2][:20],
                             provenance=nodes[k]['site']['provenance'],
                             category=nodes[k]['site']['category'],
                             how=nodes[k]['how'], reads=len(labels),
                             scored=n, segregation=round(s_, 4),
                             hap1_is='MATERNAL' if straight >= crossed else 'PATERNAL',
                             in_gap=int(gl <= k[0] <= gr)))
        # Does the chain carry the two flanks at a consistent orientation?
        verdict = 'no chain'
        if rows:
            first, last = rows[0], rows[-1]
            reliable = [r for r in rows if r['segregation'] >= 0.90]
            verdict = 'CORRECT' if (first['hap1_is'] == last['hap1_is'] and
                                    len(reliable) >= 2 and
                                    all(r['hap1_is'] == first['hap1_is'] for r in reliable)) \
                else 'SWITCHED'
        if rows:
            out = a.outdir / f'chain_{gl}{a.tag}.tsv'
            with out.open('w', newline='') as stream:
                wr = csv.DictWriter(stream, fieldnames=list(rows[0]), delimiter='\t',
                                    lineterminator='\n')
                wr.writeheader()
                wr.writerows(rows)
        # Truth-free chain statistics: the weakest link by each measure, and how
        # far it reaches. These are the quantities a selection rule could use.
        min_cons = min((e['cons'] for e in path), default=0.0)
        weak = min(path, key=lambda e: (e['cons'], e['shared'])) if path else None
        longest = max((abs(e['b'][0] - e['a'][0]) for e in path), default=0)
        summary.append(dict(
            gap_left=gl, gap_bp=int(w['gap_bp']), usable_sites=len(nodes), edges=len(edges),
            chain_found=int(bool(chain_keys)), bottleneck=bottleneck,
            chain_sites=len(chain_keys),
            in_gap_chain_sites=sum(r['in_gap'] for r in rows),
            from_catalog_only=sum(1 for r in rows if r['provenance'] == 'catalog'),
            min_consistency=round(min_cons, 4),
            weakest_link=('%d-%d' % (weak['a'][0], weak['b'][0])) if weak else '-',
            weakest_shared=weak['shared'] if weak else 0,
            weakest_cons=round(weak['cons'], 4) if weak else 0.0,
            longest_link_bp=longest,
            min_segregation=round(min(seg), 4) if seg else 0.0,
            orientation=verdict))
        print('  w%-10d sites=%-3d inGap=%-3d catOnly=%-3d bottle=%-3d minCons=%.3f '
              'weakest=%-19s %2dr/%.2f  maxLink=%5.1fkb minSeg=%.3f  %s' % (
                  gl, len(chain_keys), sum(r['in_gap'] for r in rows),
                  sum(1 for r in rows if r['provenance'] == 'catalog'), bottleneck,
                  min_cons, weak and ('%d-%d' % (weak['a'][0], weak['b'][0])) or '-',
                  weak['shared'] if weak else 0, weak['cons'] if weak else 0.0,
                  longest / 1e3, min(seg) if seg else 0.0, verdict))

    out = a.outdir / f'chain_results{a.tag}.tsv'
    with out.open('w', newline='') as stream:
        wr = csv.DictWriter(stream, fieldnames=list(summary[0]), delimiter='\t',
                            lineterminator='\n')
        wr.writeheader()
        wr.writerows(summary)
    print(f'wrote {out}')


if __name__ == '__main__':
    main()
