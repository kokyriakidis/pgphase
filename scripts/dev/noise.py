#!/usr/bin/env python3
"""EM-reliable sites the parents do not split, and the reads they mislabel (evaluation only).

Usage: noise.py DIAG_DIR [--min-parent N] [--examples N] [--rival NAME]

Every EM site the solve keeps (learned error < 0.3, not U) with at least
--min-parent scored reads from each parent is classed by the parents' sides:
  real   the parents' majorities differ, purity >= 0.9
  mid    the parents' majorities differ, purity in (0.65, 0.9)
  noise  the majorities agree, or purity <= 0.65
Each discordant arms_phaseable read is classed by its own reliable calls (as
readaudit.py). Sites are then tabled by production features (origin, kind,
learned error, minor allele fraction, reads, mean MAPQ of the callers, the
share of callers the EM labelled, site density around the site, repeat run),
real against noise, with the discordant reads each noise class labels.
--rival NAME counts only the discordant reads NAME phases correctly.
"""
import argparse
import bisect
import collections
import os
import sys

import pysam

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import diaglib as D  # noqa: E402
import score  # noqa: E402
from islands import repeat_run  # noqa: E402


def bin_of(edges, v, fmt="{:g}"):
    i = max(0, bisect.bisect_right(edges, v) - 1)
    hi = edges[i + 1] if i + 1 < len(edges) else None
    return f"[{fmt.format(edges[i])}, {fmt.format(hi) if hi is not None else 'inf'})"


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("diag")
    ap.add_argument("--min-parent", type=int, default=3)
    ap.add_argument("--examples", type=int, default=0)
    ap.add_argument("--rival", default="")
    a = ap.parse_args()
    contig = os.environ["PGDEV_CONTIG"]
    fa = pysam.FastaFile(os.environ["PGDEV_REF"])
    ref_name = contig if contig in fa.references else contig.split("#")[-1]
    scope = score.read_scope()
    classes = score.read_classes(a.diag)
    rival = score.read_classes(D.competitors()[a.rival]) if a.rival else None
    # per read: the longest dump line (a read can appear in two chunks)
    reads = {}
    for f in D.chunk_files(os.path.join(a.diag, "em"), "em"):
        with open(f) as fh:
            for line in fh:
                p = line.rstrip("\n").split("\t")
                first = 6 if len(p) > 5 and ":" not in p[5] else 5
                obs = []
                for t in p[first:]:
                    x = t.split(":")
                    if x[1] in ("U", "L") or x[2] not in ("0", "1"):
                        continue
                    obs.append((int(x[0]), x[1], int(x[2]), float(x[3]), int(x[4])))
                old = reads.get(p[0])
                if old is None or len(obs) > len(old[2]):
                    reads[p[0]] = (int(p[1]), p[2] in ("1", "2"), obs)
    sites = {}  # (pos, kind) -> dict
    for q, (mq, lab, obs) in reads.items():
        sc = scope.get(q)
        for pos, kind, side, err, blk in obs:
            s = sites.get((pos, kind))
            if s is None:
                s = sites[(pos, kind)] = dict(err=err, par=[0, 0, 0, 0], n=0, ones=0, mq=0, lab=0, blk=blk)
            s["n"] += 1
            s["ones"] += side
            s["mq"] += mq
            s["lab"] += lab
            if sc is not None:
                s["par"][(0 if sc[2] else 2) + side] += 1
    for key, s in sites.items():
        c = s["par"]
        if c[0] + c[1] < a.min_parent or c[2] + c[3] < a.min_parent:
            s["cls"], s["mat"] = None, None
            continue
        ms, ps = (0 if c[0] >= c[1] else 1), (0 if c[2] >= c[3] else 1)
        if ms == ps:
            s["cls"], s["mat"] = "noise", None
            continue
        pur = (c[ms] + c[2 + ps]) / sum(c)
        s["cls"] = "real" if pur >= 0.9 else "noise" if pur <= 0.65 else "mid"
        s["mat"] = ms
    # discordant arms_phaseable reads by their own reliable calls
    disc = collections.Counter()
    noise_reads = collections.defaultdict(set)  # site key -> discordant reads it labels (noise-labelled class)
    for q, (mq, lab, obs) in reads.items():
        sc = scope.get(q)
        if sc is None or sc[1] != "arms_phaseable" or classes.get(q) != "discordant":
            continue
        if rival is not None and rival.get(q) != "correct":
            continue
        right = wrong = noise = 0
        nkeys = []
        for pos, kind, side, err, blk in obs:
            if err >= 0.3:
                continue
            s = sites[(pos, kind)]
            if s["cls"] == "real":
                if (side == s["mat"]) == sc[2]:
                    right += 1
                else:
                    wrong += 1
            elif s["cls"] == "noise":
                noise += 1
                nkeys.append((pos, kind))
        sub = ("no clear site" if right + wrong + noise == 0 else
               "wrong calls at real sites" if wrong > right else
               "labelled by noise sites" if noise > right + wrong else
               "right calls (block)" if wrong == 0 else "mixed")
        disc[sub] += 1
        if sub == "labelled by noise sites":
            for k in nkeys:
                noise_reads[k].add(q)
    n_disc = sum(disc.values())
    print(f"{contig}: discordant arms_phaseable reads{' ' + a.rival + ' phases correctly' if a.rival else ''} {n_disc:,}")
    for k, n in disc.most_common():
        print(f"  {k:28} {n:6,}")
    rel = {k: s for k, s in sites.items() if s["err"] < 0.3 and s["cls"] is not None}
    cnt = collections.Counter(s["cls"] for s in rel.values())
    print(f"reliable EM sites with >= {a.min_parent} reads per parent: " + ", ".join(f"{c} {n:,}" for c, n in cnt.most_common()))
    print(f"noise sites labelling a noise-labelled discordant read: {len(noise_reads):,}")
    # site density: reliable sites within +-1 kb
    allpos = sorted(p for p, k in sites if sites[(p, k)]["err"] < 0.3)
    bad = set(noise_reads)
    feats = {}
    for key, s in rel.items():
        pos, kind = key
        maf = min(s["ones"], s["n"] - s["ones"]) / s["n"]
        dens = bisect.bisect_right(allpos, pos + 1000) - bisect.bisect_left(allpos, pos - 1000)
        run = repeat_run(fa, ref_name, pos) if kind != "C" or True else 0
        feats[key] = {
            "origin/kind": kind,
            "learned error": bin_of([0, 0.02, 0.05, 0.1, 0.15, 0.2, 0.3], s["err"]),
            "minor allele fraction": bin_of([0, 0.15, 0.25, 0.3, 0.35, 0.4, 0.45], maf),
            "reads": bin_of([0, 10, 20, 30, 40, 60, 100], s["n"], "{:d}"),
            "mean MAPQ": bin_of([0, 30, 50, 59, 60.01], s["mq"] / s["n"]),
            "labelled share": bin_of([0, 0.5, 0.8, 0.9, 0.95, 1.0001], s["lab"] / s["n"]),
            "sites within 1 kb": bin_of([0, 2, 4, 6, 10, 20], dens, "{:d}"),
            "repeat run": bin_of([0, 2, 4, 6, 10, 20], run, "{:d}"),
        }
    for g in ("origin/kind", "learned error", "minor allele fraction", "reads", "mean MAPQ", "labelled share",
              "sites within 1 kb", "repeat run"):
        tab = collections.defaultdict(collections.Counter)
        for key, s in rel.items():
            v = feats[key][g]
            tab[v][s["cls"]] += 1
            if key in bad:
                tab[v]["bad"] += 1
                tab[v]["badreads"] += len(noise_reads[key])
        print(f"by {g}:  {'real':>8} {'mid':>7} {'noise':>7} {'noise %':>8} {'bad sites':>10} {'their disc reads':>17}")
        for v in sorted(tab, key=lambda v: (float(v.split(",")[0].strip("[")) if v[:1] == "[" else 0, v)):
            t = tab[v]
            tot = t["real"] + t["mid"] + t["noise"]
            print(f"  {v:22} {t['real']:8,} {t['mid']:7,} {t['noise']:7,} {100 * t['noise'] / max(1, tot):7.1f}% "
                  f"{t['bad']:10,} {t['badreads']:17,}")
    if a.examples:
        print("worst noise sites (most noise-labelled discordant reads): pos kind err n maf mean-MAPQ labelled parents(m0 m1 p0 p1)")
        for key in sorted(bad, key=lambda k: -len(noise_reads[k]))[: a.examples]:
            s = sites[key]
            print(f"  {key[0]:>11,} {key[1]} {s['err']:.3f} {s['n']:4} {min(s['ones'], s['n'] - s['ones']) / s['n']:.2f} "
                  f"{s['mq'] / s['n']:5.1f} {s['lab'] / s['n']:.2f} {s['par']} disc {len(noise_reads[key])}")


if __name__ == "__main__":
    main()
