#!/usr/bin/env python3
"""Simulate combining a primary phasing (EM) with a secondary one (k-means) per
read, without truth, and score the result with truth (evaluation only).

Usage: combine.py PRIMARY_DIR SECONDARY_DIR [--min-shared N] [--min-agree F]

Each secondary phase set is mapped onto the primary phase set it shares the
most commonly labelled reads with, when at least N shared reads agree on one
relative orientation at fraction >= F. Then:
  fill   a read the primary leaves unphased takes the mapped secondary label;
  veto   a read both label, but in opposite mapped orientation, is unlabelled.
Scored like score.py: per phase set majority parental orientation.
"""
import argparse
import collections
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import score  # noqa: E402


def classes(tags, scope):
    votes = collections.defaultdict(lambda: [0, 0])
    for q, (hp, ps) in tags.items():
        votes[ps][(hp == 1) == scope[q][2]] += 1
    orient = {ps: v[1] >= v[0] for ps, v in votes.items()}
    out = {}
    for q, info in scope.items():
        t = tags.get(q)
        out[q] = "unphased" if t is None else ("correct" if ((t[0] == 1) == info[2]) == orient[t[1]] else "discordant")
    return out


def summary(cls, scope, groups=("arms_phaseable", "cen_phaseable", "arms_nohet")):
    out = {}
    for g in groups:
        c = collections.Counter(cls[q] for q, v in scope.items() if v[1] == g)
        out[g] = (c["correct"], c["discordant"], c["unphased"])
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("primary")
    ap.add_argument("secondary")
    ap.add_argument("--min-shared", type=int, default=10)
    ap.add_argument("--min-agree", type=float, default=0.95)
    a = ap.parse_args()
    scope = score.read_scope()
    p = score.read_tags(a.primary, scope)
    s = score.read_tags(a.secondary, scope)
    # Map secondary PS -> (primary PS, flip) by commonly labelled reads.
    pair = collections.defaultdict(lambda: [0, 0])  # (sps, pps) -> [same, opposite]
    for q, (shp, sps) in s.items():
        if q in p:
            php, pps = p[q]
            pair[(sps, pps)][shp != php] += 1
    best = {}
    for (sps, pps), (same, opp) in pair.items():
        n = same + opp
        if n < a.min_shared or max(same, opp) < a.min_agree * n:
            continue
        if sps not in best or n > best[sps][2]:
            best[sps] = (pps, opp > same, n)
    fill, veto = dict(p), dict(p)
    both = dict(p)
    n_fill = n_veto = 0
    for q, (shp, sps) in s.items():
        m = best.get(sps)
        if m is None:
            continue
        pps, flip, _ = m
        hp = (3 - shp) if flip else shp
        if q not in p:
            fill[q] = (hp, pps)
            both[q] = (hp, pps)
            n_fill += 1
        elif p[q] != (hp, pps) and p[q][1] == pps:
            del veto[q]
            del both[q]
            n_veto += 1
    print(f"secondary phase sets mapped: {len(best):,} of {len({v[1] for v in s.values()}):,}; "
          f"reads filled {n_fill:,}, vetoed {n_veto:,}")
    rows = [("primary", p), ("secondary", s), ("fill", fill), ("veto", veto), ("fill+veto", both)]
    base = summary(classes(p, scope), scope)
    print(f"{'':12}" + "".join(f"{g:>34}" for g in base))
    for name, tags in rows:
        sm = summary(classes(tags, scope), scope)
        cells = []
        for g, (c, d, u) in sm.items():
            bc, bd, bu = base[g]
            cells.append(f"{c:>8,} {d:>6,} {u:>6,} ({c - bc:+,}/{d - bd:+,})")
        print(f"{name:12}" + "".join(f"{x:>34}" for x in cells))


if __name__ == "__main__":
    main()
