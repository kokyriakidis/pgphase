#!/usr/bin/env python3
"""Prototype: EM whose site errors are estimated leave-one-out, on the EM dump's
read x site matrices, against the EM's own in-chunk labels (evaluation only).

Usage: proto_loo.py EM_DUMP_DIR [--iters N] [--prior E] [--prior-weight W] [--cut T]

Model as the production EM: two haplotypes, per-site error e_s, read posterior
from all its observations. Difference: a site's phase and error are learned
from each read's posterior WITHOUT that site's own observation, so a site
cannot confirm itself; a read with no other evidence carries no information
about the site (its leave-one-out posterior is 0.5). A weak prior (W pseudo
reads at error E) keeps sites with little independent evidence near E.
Blocks are cut where flipping the rest costs less than T log-likelihood over
every read; reads are labelled in their strongest block at posterior >= 0.8
from sites with error < 0.3.
"""
import argparse
import collections
import glob
import multiprocessing as mp
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import proto_mec as P  # noqa: E402  (parse, classes)
import score  # noqa: E402

A = None
EMIN, EMAX, UNREL, POST = 0.01, 0.45, 0.3, 0.8


def solve(path):
    names, mapq, em, obs = P.parse(path)
    chunk = os.path.basename(path)
    keys = sorted({(o[0], o[1]) for r in obs for o in r})
    sid = {k: i for i, k in enumerate(keys)}
    n = len(keys)
    err0 = np.full(n, 0.2)
    rr, ss, bb = [], [], []
    for ri, r in enumerate(obs):
        seen = {}
        for pos, kind, side, e in r:
            seen[sid[(pos, kind)]] = side
            err0[sid[(pos, kind)]] = e
        for s, side in seen.items():
            rr.append(ri); ss.append(s); bb.append(side)
    out = {}
    if n == 0 or not rr:
        return chunk, out, em, names, [len(r) for r in obs]
    rr, ss, bb = np.array(rr), np.array(ss), np.array(bb)
    learns = np.array(mapq) >= 20
    err = np.clip(err0, EMIN, EMAX)
    # starting read posteriors: the EM's labels in its own gauge are per block,
    # so start from the EM's site orientation implied by its labelled reads
    hap = np.array([e[0] for e in em])
    p = np.where(hap == 1, 0.95, np.where(hap == 2, 0.05, 0.5))
    phase = np.zeros(n, dtype=int)  # 0: hap1 carries side 0

    def site_ll(phase, err):
        e = err[ss]
        match1 = bb == phase[ss]
        return np.where(match1, np.log(1 - e), np.log(e)), np.where(match1, np.log(e), np.log(1 - e))

    for it in range(A.iters):
        # leave-one-out posterior of each observation's read
        l1, l2 = site_ll(phase, err) if it > 0 else (None, None)
        if it == 0:
            q = p[rr]
        else:
            t1 = np.bincount(rr, weights=l1, minlength=len(names))
            t2 = np.bincount(rr, weights=l2, minlength=len(names))
            d = (t2[rr] - l2) - (t1[rr] - l1)
            q = 1.0 / (1.0 + np.exp(np.clip(d, -700, 700)))
        w = learns[rr].astype(float)
        # Each read votes with its leave-one-out label, discounted by how sure
        # that label is (c = |2q - 1|). With label accuracy (1 + c) / 2 the
        # expected disagreement is (1 - c) / 2 + e * c, so
        # e = sum(d - (1 - c) / 2) / sum(c).
        c = np.abs(2 * q - 1) * w
        lab1 = q >= 0.5                       # read labelled hap1 without this site
        side_on_hap1 = np.where(lab1, bb, 1 - bb)  # the side this read puts on hap1
        v0 = np.bincount(ss, weights=c * (side_on_hap1 == 0), minlength=n)
        v1 = np.bincount(ss, weights=c * (side_on_hap1 == 1), minlength=n)
        phase = np.where(v0 >= v1, 0, 1)
        d = (side_on_hap1 != phase[ss]).astype(float)
        num = np.bincount(ss, weights=w * (d - (1 - np.abs(2 * q - 1)) / 2), minlength=n)
        den = np.bincount(ss, weights=c, minlength=n)
        err = np.clip((num + A.prior_weight * A.prior) / (den + A.prior_weight), EMIN, EMAX)
    l1, l2 = site_ll(phase, err)
    # block cuts over every read
    diff = np.zeros(n + 1)
    order = np.lexsort((ss, rr))
    rr_o, ss_o, l1_o, l2_o = rr[order], ss[order], l1[order], l2[order]
    bounds = np.flatnonzero(np.diff(rr_o)) + 1
    for seg in np.split(np.arange(len(rr_o)), bounds):
        if len(seg) < 2:
            continue
        a1, a2 = l1_o[seg], l2_o[seg]
        T1, T2 = a1.sum(), a2.sum()
        c1, c2 = np.cumsum(a1)[:-1], np.cumsum(a2)[:-1]
        now = np.logaddexp(T1, T2)
        flip = np.logaddexp(c1 + (T2 - c2), c2 + (T1 - c1))
        gain = now - flip
        s_here = ss_o[seg]
        np.add.at(diff, s_here[:-1], gain)
        np.add.at(diff, s_here[1:], -gain)
    delta = np.cumsum(diff)[:n]
    block = np.zeros(n, dtype=int)
    cur = 0
    for s in range(n):
        block[s] = cur
        if s + 1 < n and delta[s] < A.cut:
            cur = s + 1
    ok = err[ss] < UNREL
    per = collections.defaultdict(lambda: [0.0, 0.0])
    for ri, s, x1, x2 in zip(rr[ok], ss[ok], l1[ok], l2[ok]):
        v = per[(ri, block[s])]
        v[0] += x1; v[1] += x2
    best = {}
    for (ri, blk), (x1, x2) in per.items():
        m = x1 - x2
        if abs(m) >= np.log(POST / (1 - POST)) and (ri not in best or abs(m) > abs(best[ri][1])):
            best[ri] = (blk, m)
    for ri, (blk, m) in best.items():
        out[names[ri]] = (1 if m > 0 else 2, (chunk, int(blk)))
    return chunk, out, em, names, [len(r) for r in obs]


def init(args):
    global A
    A = args


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("dump")
    ap.add_argument("--iters", type=int, default=15)
    ap.add_argument("--prior", type=float, default=0.2)
    ap.add_argument("--prior-weight", type=float, default=2.0)
    ap.add_argument("--cut", type=float, default=4.0)
    ap.add_argument("--jobs", type=int, default=10)
    ap.add_argument("--chunks", default="")
    args = ap.parse_args()
    files = sorted(glob.glob(os.path.join(args.dump, "em.*.tsv")))
    if args.chunks:
        files = [f for f in files if any(c in f for c in args.chunks.split(","))]
    with mp.Pool(args.jobs, initializer=init, initargs=(args,)) as pool:
        results = pool.map(solve, files, chunksize=1)
    best_visit = {}
    for chunk, out, em, names, nobs in results:
        for i, q in enumerate(names):
            if q not in best_visit or nobs[i] > best_visit[q][0]:
                best_visit[q] = (nobs[i], chunk)
    loo, emt = {}, {}
    for chunk, out, em, names, nobs in results:
        for i, q in enumerate(names):
            if best_visit[q][1] != chunk:
                continue
            if q in out:
                loo[q] = out[q]
            if em[i][0] in (1, 2):
                emt[q] = (em[i][0], (chunk, em[i][1]))
    full = score.read_scope()
    scope = {q: v for q, v in full.items() if q in best_visit}
    groups = ("arms_phaseable", "arms_nohet", "cen_phaseable")
    print(f"iters {args.iters} prior {args.prior}x{args.prior_weight} cut {args.cut}")
    for name, tags in (("EM (dump)", emt), ("LOO EM", loo)):
        cl = P.classes(tags, scope)
        cells = []
        for g in groups:
            c = collections.Counter(cl[q] for q, v in scope.items() if v[1] == g)
            cells.append(f"{g} {c['correct']:,}/{c['discordant']:,}/{c['unphased']:,}")
        print(f"  {name:10} " + "   ".join(cells) + f"   phase sets {len({t[1] for t in tags.values()}):,}")


if __name__ == "__main__":
    main()
