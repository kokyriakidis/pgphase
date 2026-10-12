#!/usr/bin/env python3
"""A numpy replica of the production union EM (src/union_em.cpp), run on a v2
EM dump, so that algorithm variants can be tested one change at a time
(evaluation only).

Usage: proto_em.py EM_DUMP_DIR [--variant NAME] [--chunks a,b] [--jobs J]

Variants:
  em        the production algorithm (validation: should match the dump)
  loo       site phase and error learned leave-one-out (a site cannot
            confirm itself)
  daem      deterministic annealing: posteriors tempered (beta 0.2 -> 1)
            during the first round, so starting labels cannot lock a site in
Reports, for the dump's own labels and the variant's: arms_phaseable,
arms_nohet and cen_phaseable correct / discordant / unphased.
"""
import argparse
import collections
import glob
import multiprocessing as mp
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import proto_mec as P  # noqa: E402  (classes)
import score  # noqa: E402

dedup = []
EMIN, EMAX, CLEAN0, OTHER0, UNREL = 0.01, 0.45, 0.05, 0.2, 0.3
ITERS, ROUNDS, POST, CUT, WEAK_CUT, LEARN_MAPQ = 15, 30, 0.8, 4.0, 10.0, 20
A = None


def load(path):
    names, mapq, final, start, obs = [], [], [], [], []
    for line in open(path):
        p = line.rstrip("\n").split("\t")
        names.append(p[0]); mapq.append(int(p[1])); final.append((int(p[2]), int(p[3])))
        start.append(int(p[5]))
        o = []
        for tok in p[6:]:
            f = tok.split(":")
            o.append((int(f[0]), f[1], int(f[2]), int(f[5]), float(f[6])))
        obs.append(o)
    return names, np.array(mapq), final, np.array(start), obs


def solve(path):
    names, mapq, final, start, obs = load(path)
    chunk = os.path.basename(path)
    # global site order: dump observations are in site-index order per read
    keys = sorted({(o[0], o[1]) for r in obs for o in r}, key=lambda k: (k[0], k[1] == "L"))
    sid = {k: i for i, k in enumerate(keys)}
    n = len(keys)
    out = {}
    if n < 2:
        return chunk, out, final, names, [len(r) for r in obs]
    kind = np.array([k[1] for k in keys])
    weak = np.zeros(n, dtype=bool)
    rr, ss, bb, ee = [], [], [], []
    for ri, r in enumerate(obs):
        for pos, k, side, w, e in r:
            s = sid[(pos, k)]
            weak[s] |= bool(w)
            rr.append(ri); ss.append(s); bb.append(side); ee.append(e)
    rr, ss, bb, ee = map(np.array, (rr, ss, bb, ee))
    if A.dedup_window > 0:
        # a window and the sites inside it are one locus: keep the window's call
        pos = np.array([k[0] for k in keys])
        is_l = kind == "L"
        win_reads = collections.defaultdict(list)
        for i in np.flatnonzero(is_l[ss]):
            win_reads[rr[i]].append(pos[ss[i]])
        keep = np.ones(len(rr), dtype=bool)
        for i in range(len(rr)):
            if is_l[ss[i]]:
                continue
            for wp in win_reads.get(rr[i], ()):
                if wp <= pos[ss[i]] < wp + A.dedup_window:
                    keep[i] = False
                    break
        dedup.append(int((~keep).sum()))
        rr, ss, bb, ee = rr[keep], ss[keep], bb[keep], ee[keep]
    order = np.lexsort((ss, rr))
    rr, ss, bb, ee = rr[order], ss[order], bb[order], ee[order]
    nr = len(names)
    learns = mapq >= LEARN_MAPQ
    phase = np.zeros(n, dtype=int)
    err = np.where(kind == "C", CLEAN0, OTHER0).astype(float)
    post = np.where(start == 1, 0.95, np.where(start == 2, 0.05, 0.5))
    has_obs = np.bincount(rr, minlength=nr) > 0
    first_s = np.full(nr, n)
    np.minimum.at(first_s, rr, ss)
    beta = 1.0

    def obs_err():
        return np.minimum(EMAX, 1 - (1 - err[ss]) * (1 - ee))

    def lls():
        e = obs_err()
        m1 = bb == phase[ss]
        return np.where(m1, np.log1p(-e), np.log(e)), np.where(m1, np.log(e), np.log1p(-e))

    def m_step():
        nonlocal phase, err
        w = learns[rr].astype(float)
        if A.variant == "loo" and m_step.started:
            l1, l2 = lls()
            t1 = np.bincount(rr, weights=l1, minlength=nr)
            t2 = np.bincount(rr, weights=l2, minlength=nr)
            d = np.clip((t2[rr] - l2) - (t1[rr] - l1), -700, 700)
            q = 1 / (1 + np.exp(d))
        else:
            q = post[rr]
        if A.variant == "loo2" and m_step.started:
            l1, l2 = lls()
            t1 = np.bincount(rr, weights=l1, minlength=nr)
            t2 = np.bincount(rr, weights=l2, minlength=nr)
            d = np.clip((t2[rr] - l2) - (t1[rr] - l1), -700, 700)
            q = 1 / (1 + np.exp(d))
            c = np.abs(2 * q - 1) * w
            side_h1 = np.where(q >= 0.5, bb, 1 - bb)
            v0 = np.bincount(ss, weights=c * (side_h1 == 0), minlength=n)
            v1 = np.bincount(ss, weights=c * (side_h1 == 1), minlength=n)
            den = v0 + v1
            upd = den >= 1e-6
            phase = np.where(upd, np.where(v0 >= v1, 0, 1), phase)
            dis = (side_h1 != phase[ss]).astype(float)
            num = np.bincount(ss, weights=w * (dis - (1 - np.abs(2 * q - 1)) / 2), minlength=n)
            err = np.clip((num + A.prior_weight * A.prior) / (den + A.prior_weight), EMIN, EMAX)
            return
        w0 = np.bincount(ss, weights=w * np.where(bb == 0, q, 1 - q), minlength=n)
        w1 = np.bincount(ss, weights=w * np.where(bb == 1, q, 1 - q), minlength=n)
        tot = w0 + w1
        upd = tot >= 1e-6
        phase = np.where(upd, np.where(w0 >= w1, 0, 1), phase)
        err = np.where(upd, np.clip(np.minimum(w0, w1) / np.where(upd, tot, 1), EMIN, EMAX), err)
        m_step.started = True
    m_step.started = False

    def e_step():
        nonlocal post
        l1, l2 = lls()
        t1 = np.bincount(rr, weights=l1, minlength=nr)
        t2 = np.bincount(rr, weights=l2, minlength=nr)
        d = np.clip(beta * (t2 - t1), -700, 700)
        post = np.where(has_obs, 1 / (1 + np.exp(d)), post)

    def boundary_delta(every_read, strong_only=False):
        l1, l2 = lls()
        keep = np.ones(len(rr), dtype=bool)
        if strong_only:
            keep = ~weak[ss]
        if not every_read:
            keep &= learns[rr]
        r_, s_, a1, a2 = rr[keep], ss[keep], l1[keep], l2[keep]
        if len(r_) < 2:
            return np.zeros(n)
        starts = np.r_[0, np.flatnonzero(np.diff(r_)) + 1]
        ends = np.r_[starts[1:], len(r_)]
        c1, c2 = np.cumsum(a1), np.cumsum(a2)
        base1 = np.repeat(np.r_[0, c1][starts], ends - starts)
        base2 = np.repeat(np.r_[0, c2][starts], ends - starts)
        p1, p2 = c1 - base1, c2 - base2              # prefix sums within read
        T1 = np.repeat(p1[ends - 1], ends - starts)
        T2 = np.repeat(p2[ends - 1], ends - starts)
        nxt = np.r_[np.arange(1, len(r_)), 0]
        same_read = np.r_[r_[1:] == r_[:-1], False]
        valid = same_read & (s_[nxt] != s_)
        # reads whose observations sit on one site only do not count
        span_ok = np.repeat(s_[starts] != s_[ends - 1], ends - starts)
        valid &= span_ok
        gain = np.logaddexp(T1, T2) - np.logaddexp(p1 + (T2 - p2), p2 + (T1 - p1))
        diff = np.zeros(n + 1)
        np.add.at(diff, s_[valid], gain[valid])
        np.add.at(diff, s_[nxt][valid], -gain[valid])
        return np.cumsum(diff)[:n]

    for rnd in range(ROUNDS):
        for it in range(ITERS):
            if A.variant == "daem" and rnd == 0:
                beta = 0.2 + 0.8 * it / (ITERS - 1)
            else:
                beta = 1.0
            m_step()
            e_step()
        delta = boundary_delta(False)
        neg = np.flatnonzero(delta[:-1] < -1e-6)
        if len(neg) == 0:
            break
        worst = neg[np.argmin(delta[neg])]
        phase[worst + 1:] ^= 1
        flip = has_obs & (first_s > worst)
        post = np.where(flip, 1 - post, post)
    delta = boundary_delta(True)
    strong = boundary_delta(True, True)
    # reads with observations on both sides of each boundary
    span = np.zeros(n + 1)
    starts = np.r_[0, np.flatnonzero(np.diff(rr)) + 1]
    ends = np.r_[starts[1:], len(rr)]
    np.add.at(span, ss[starts], 1)
    np.add.at(span, ss[ends - 1], -1)
    span = np.cumsum(span)[:n]
    block = np.zeros(n, dtype=int)
    cur = 0
    for k in range(n):
        block[k] = cur
        if k + 1 < n and (delta[k] < CUT or (strong[k] < CUT and delta[k] < WEAK_CUT) or span[k] < A.min_join_reads):
            cur = k + 1
    l1, l2 = lls()
    ok = err[ss] < A.unrel
    m = collections.defaultdict(float)
    for ri, s, x1, x2 in zip(rr[ok], ss[ok], l1[ok], l2[ok]):
        m[(ri, block[s])] += x1 - x2
    best = {}
    lim = np.log(A.post / (1 - A.post))
    for (ri, blk), v in m.items():
        if abs(v) >= lim and (ri not in best or abs(v) > abs(best[ri][1])):
            best[ri] = (blk, v)
    for ri, (blk, v) in best.items():
        out[names[ri]] = (1 if v > 0 else 2, (chunk, int(blk)))
    return chunk, out, final, names, [len(r) for r in obs]


def init(args):
    global A
    A = args


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("dump")
    ap.add_argument("--variant", default="em")
    ap.add_argument("--chunks", default="")
    ap.add_argument("--jobs", type=int, default=10)
    ap.add_argument("--post", type=float, default=POST, help="read label posterior")
    ap.add_argument("--unrel", type=float, default=UNREL, help="site error at which a site stops labelling")
    ap.add_argument("--min-join-reads", type=int, default=0)
    ap.add_argument("--dedup-window", type=int, default=0)
    ap.add_argument("--prior", type=float, default=0.2)
    ap.add_argument("--prior-weight", type=float, default=2.0)
    args = ap.parse_args()
    files = sorted(glob.glob(os.path.join(args.dump, "em.*.tsv")))
    if args.chunks:
        files = [f for f in files if any(c in f for c in args.chunks.split(","))]
    with mp.Pool(args.jobs, initializer=init, initargs=(args,)) as pool:
        results = pool.map(solve, files, chunksize=1)
    best_visit = {}
    for chunk, out, final, names, nobs in results:
        for i, q in enumerate(names):
            if q not in best_visit or nobs[i] > best_visit[q][0]:
                best_visit[q] = (nobs[i], chunk)
    mine, emt = {}, {}
    for chunk, out, final, names, nobs in results:
        for i, q in enumerate(names):
            if best_visit[q][1] != chunk:
                continue
            if q in out:
                mine[q] = out[q]
            if final[i][0] in (1, 2):
                emt[q] = (final[i][0], (chunk, final[i][1]))
    scope = {q: v for q, v in score.read_scope().items() if q in best_visit}
    agree = sum(1 for q in scope if (q in mine) == (q in emt))
    print(f"variant {args.variant} post {args.post} unrel {args.unrel} min_join {args.min_join_reads}: reads labelled the same way as the dump (labelled/unlabelled): "
          f"{agree:,} / {len(scope):,}")
    for name, tags in (("EM (dump)", emt), (f"replica:{args.variant}", mine)):
        cl = P.classes(tags, scope)
        cells = []
        for g in ("arms_phaseable", "arms_nohet", "cen_phaseable"):
            c = collections.Counter(cl[q] for q, v in scope.items() if v[1] == g)
            cells.append(f"{g} {c['correct']:,}/{c['discordant']:,}/{c['unphased']:,}")
        print(f"  {name:14} " + "   ".join(cells) + f"   phase sets {len({t[1] for t in tags.values()}):,}")


if __name__ == "__main__":
    main()
