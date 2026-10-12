#!/usr/bin/env python3
"""Prototype: weighted-MEC beam search (HiPhase-style) on the EM dump's read x
site matrices, compared with the EM's own in-chunk labels (evaluation only).

Usage: proto_mec.py EM_DUMP_DIR [--beam B] [--cut T] [--margin M] [--jobs J]

Per chunk: sites in position order; each takes one of four genotypes
(x|y, y|x, x/x, y/y). A read costs min(mismatches vs hap1, vs hap2); the
search keeps the B cheapest partial solutions site by site (no per-site error
is learned). A site whose reads are all new starts a component and takes only
one het orientation. Blocks are cut where flipping the rest would cost fewer
than T mismatches; a read is labelled in the block where |cost difference| is
largest, if that is at least M. Truth only scores the result.
"""
import argparse
import collections
import glob
import multiprocessing as mp
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import score  # noqa: E402

A = None


def parse(path):
    names, mapq, em, obs = [], [], [], []
    for line in open(path):
        p = line.rstrip("\n").split("\t")
        names.append(p[0])
        mapq.append(int(p[1]))
        em.append((int(p[2]), int(p[3])))
        o = []
        first = 6 if len(p) > 5 and ":" not in p[5] else 5  # v2 dumps carry a start label
        for tok in p[first:]:
            pos, kind, side, err = tok.split(":")[:4]
            o.append((int(pos), kind, int(side), float(err)))
        obs.append(o)
    return names, mapq, em, obs


def solve(path):
    names, mapq, em, obs = parse(path)
    chunk = os.path.basename(path)
    keys = sorted({(o[0], o[1]) for r in obs for o in r})
    site_err = {}
    for r in obs:
        for o in r:
            site_err[(o[0], o[1])] = o[3]
    sid = {k: i for i, k in enumerate(keys)}
    n_sites, n_reads = len(keys), len(names)
    site_reads = [[] for _ in range(n_sites)]
    site_sides = [[] for _ in range(n_sites)]
    read_sites = []
    for ri, r in enumerate(obs):
        seen = {}
        for pos, kind, side, _ in r:
            seen[sid[(pos, kind)]] = side
        items = sorted(seen.items())
        read_sites.append(items)
        for s, side in items:
            site_reads[s].append(ri)
            site_sides[s].append(side)
    weight = np.ones(n_sites)
    for k, i in sid.items():
        if A.weights == "emerr":
            e = min(max(site_err[k], 0.01), 0.49)
            weight[i] = 0.0 if e >= A.unreliable else np.log((1 - e) / e)
        elif A.weights != "uniform":
            table = dict(kv.split("=") for kv in A.weights.split(","))
            weight[i] = float(table.get(k[1], 1.0))
    first_site = np.array([rs[0][0] if rs else -1 for rs in read_sites])
    out = {}
    if n_sites == 0:
        return chunk, out, em, names, [len(x) for x in read_sites]
    B = A.beam
    # beam: arrays a, b (per read mismatches vs hap1 / hap2), cost, history
    a = np.zeros((1, n_reads))
    b = np.zeros((1, n_reads))
    cost = np.zeros(1)
    back = []  # per site: (parent index array, option array)
    opts = np.array([[0, 1], [1, 0], [0, 0], [1, 1]])  # (hap1 allele, hap2 allele)
    for s in range(n_sites):
        rd = np.array(site_reads[s], dtype=np.int64)
        sd = np.array(site_sides[s], dtype=np.int32)
        new_component = bool(np.all(first_site[rd] == s))
        allowed = [0, 2, 3] if new_component else [0, 1, 2, 3]
        old = np.minimum(a[:, rd], b[:, rd]).sum(axis=1)
        cand_cost, cand_parent, cand_opt = [], [], []
        for o in allowed:
            h1, h2 = opts[o]
            da = (sd != h1) * weight[s]
            db = (sd != h2) * weight[s]
            new = np.minimum(a[:, rd] + da, b[:, rd] + db).sum(axis=1)
            c = cost + (new - old)
            # prefer heterozygous genotypes on ties (as HiPhase does)
            cand_cost.append(c + (0.0 if o < 2 else 1e-6))
            cand_parent.append(np.arange(len(cost)))
            cand_opt.append(np.full(len(cost), o))
        cc = np.concatenate(cand_cost)
        cp = np.concatenate(cand_parent)
        co = np.concatenate(cand_opt)
        keep = np.argsort(cc, kind="stable")[:B]
        pa, po = cp[keep], co[keep]
        na, nb = a[pa].copy(), b[pa].copy()
        h1 = opts[po][:, 0][:, None]
        h2 = opts[po][:, 1][:, None]
        na[:, rd] += (sd[None, :] != h1) * weight[s]
        nb[:, rd] += (sd[None, :] != h2) * weight[s]
        a, b, cost = na, nb, np.round(cc[keep], 6)
        back.append((pa, po))
    # backtrack the best solution
    geno = np.zeros(n_sites, dtype=np.int64)
    idx = int(np.argmin(cost))
    for s in range(n_sites - 1, -1, -1):
        pa, po = back[s]
        geno[s] = po[idx]
        idx = int(pa[idx])
    het = geno < 2
    hap1 = opts[geno][:, 0]
    # boundary flip costs over every read (like the EM's every-read cut)
    diff = np.zeros(n_sites + 1)
    for ri, items in enumerate(read_sites):
        hs = [(s, side) for s, side in items if het[s]]
        if len(hs) < 2:
            continue
        w = np.array([weight[s] for s, _ in hs])
        m1 = np.array([side != hap1[s] for s, side in hs]) * w  # weighted mismatch vs hap1
        m2 = w - m1
        ta, tb = m1.sum(), m2.sum()
        ca, cb = np.cumsum(m1)[:-1], np.cumsum(m2)[:-1]
        now = min(ta, tb)
        flipped = np.minimum(ca + (tb - cb), cb + (ta - ca))
        gain = flipped - now
        for j in range(len(hs) - 1):
            diff[hs[j][0]] += gain[j]
            diff[hs[j + 1][0]] -= gain[j]
    delta = np.cumsum(diff)[:n_sites]
    block = np.zeros(n_sites, dtype=np.int64)
    cur = 0
    het_idx = np.nonzero(het)[0]
    for i, s in enumerate(het_idx):
        block[s] = cur
        if i + 1 < len(het_idx) and delta[s] < A.cut:
            cur = het_idx[i + 1]
    for ri, items in enumerate(read_sites):
        per = collections.defaultdict(lambda: [0, 0])
        for s, side in items:
            if het[s]:
                per[block[s]][0 if side == hap1[s] else 1] += weight[s]
        best = None
        for blk, (agree1, agree2) in per.items():
            m = agree1 - agree2
            if abs(m) >= A.margin and (best is None or abs(m) > abs(best[1])):
                best = (blk, m)
        if best:
            out[names[ri]] = (1 if best[1] > 0 else 2, (chunk, int(best[0])))
    return chunk, out, em, names, [len(x) for x in read_sites]


def classes(tags, scope):
    votes = collections.defaultdict(lambda: [0, 0])
    for q, (hp, ps) in tags.items():
        if q in scope:
            votes[ps][(hp == 1) == scope[q][2]] += 1
    orient = {ps: v[1] >= v[0] for ps, v in votes.items()}
    res = {}
    for q, info in scope.items():
        t = tags.get(q)
        res[q] = "unphased" if t is None else ("correct" if ((t[0] == 1) == info[2]) == orient[t[1]] else "discordant")
    return res


def init(args):
    global A
    A = args


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("dump")
    ap.add_argument("--beam", type=int, default=32)
    ap.add_argument("--cut", type=float, default=3)
    ap.add_argument("--margin", type=float, default=1)
    ap.add_argument("--weights", default="uniform", help="uniform | emerr | C=1,O=0.5,I=0.5,L=0.5")
    ap.add_argument("--unreliable", type=float, default=0.3)
    ap.add_argument("--jobs", type=int, default=10)
    ap.add_argument("--chunks", default="")
    args = ap.parse_args()
    files = sorted(glob.glob(os.path.join(args.dump, "em.*.tsv")))
    if args.chunks:
        files = [f for f in files if any(c in f for c in args.chunks.split(","))]
    with mp.Pool(args.jobs, initializer=init, initargs=(args,)) as pool:
        results = pool.map(solve, files, chunksize=1)
    # one visit per read: the chunk where it has most observations
    best_visit = {}
    for chunk, out, em, names, nobs in results:
        for i, q in enumerate(names):
            if q not in best_visit or nobs[i] > best_visit[q][0]:
                best_visit[q] = (nobs[i], chunk)
    mec, emt = {}, {}
    for chunk, out, em, names, nobs in results:
        for i, q in enumerate(names):
            if best_visit[q][1] != chunk:
                continue
            if q in out:
                mec[q] = out[q]
            if em[i][0] in (1, 2):
                emt[q] = (em[i][0], (chunk, em[i][1]))
    full = score.read_scope()
    seen = set(best_visit)
    scope = {q: v for q, v in full.items() if q in seen}
    groups = ("arms_phaseable", "arms_nohet", "cen_phaseable")
    print(f"weights {args.weights} beam {args.beam} cut {args.cut} margin {args.margin}; reads in solve with truth: {len(scope):,}")
    for name, tags in (("EM (dump)", emt), ("MEC beam", mec)):
        cl = classes(tags, scope)
        cells = []
        for g in groups:
            c = collections.Counter(cl[q] for q, v in scope.items() if v[1] == g)
            cells.append(f"{g} {c['correct']:,}/{c['discordant']:,}/{c['unphased']:,}")
        n_ps = len({t[1] for t in tags.values()})
        print(f"  {name:10} " + "   ".join(cells) + f"   phase sets {n_ps:,}")


if __name__ == "__main__":
    main()
