#!/usr/bin/env python3
"""Accuracy of the EM's per-read calls by read and site features (evaluation only).

Usage: calls.py DIAG_DIR [--sites rejected|all] [--min-err E]

Sites: EM candidate sites (not windows, not unverified), by default those the
EM rejects (learned error >= 0.3). Each site's repeat run on the reference
and whether a truth het lies within max(3, run) bp.
At truth-het sites a call is right when the read carries its parent's
majority side there (only sites where the two parents' majorities differ);
at sites without a truth het every read should carry the same allele, so a
call is wrong when it is on the site's minority side. Calls are binned by
  end    distance from the repeat tract to the nearer end of the read's alignment
  label  whether the EM labelled the read
  run    the site's repeat run
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

END_BINS = [0, 25, 50, 100, 200, 500, 1000, 2000, 10 ** 9]
RUN_BINS = [1, 4, 8, 13, 20, 30, 10 ** 6]


def bin_of(edges, v):
    i = bisect.bisect_right(edges, v) - 1
    return f"[{edges[i]}, {edges[i + 1] if edges[i + 1] < 10 ** 6 else 'inf'})"


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("diag")
    ap.add_argument("--sites", default="rejected", choices=["rejected", "all"])
    ap.add_argument("--min-err", type=float, default=0.3)
    a = ap.parse_args()
    contig = os.environ["PGDEV_CONTIG"]
    fa = pysam.FastaFile(os.environ["PGDEV_REF"])
    ref_name = contig if contig in fa.references else contig.split("#")[-1]
    scope = score.read_scope()
    min_err = a.min_err if a.sites == "rejected" else -1.0
    calls = collections.defaultdict(dict)  # pos -> read -> (side, labelled)
    errs = {}
    extra = {}  # (pos, read) -> (read error, MAPQ)
    for f in D.chunk_files(os.path.join(a.diag, "em"), "em"):
        with open(f) as fh:
            for line in fh:
                p = line.rstrip("\n").split("\t")
                if p[0] not in scope:
                    continue
                lab = p[2] in ("1", "2")
                first = 6 if len(p) > 5 and ":" not in p[5] else 5
                for t in p[first:]:
                    x = t.split(":")
                    if x[1] in ("L", "U") or x[2] not in ("0", "1"):
                        continue
                    e = float(x[3])
                    if e < min_err:
                        continue
                    pos = int(x[0])
                    errs[pos] = e
                    if p[0] not in calls[pos]:
                        calls[pos][p[0]] = (int(x[2]), lab)
                        extra[(pos, p[0])] = (float(x[6]) if len(x) > 6 else -1.0, int(p[1]))
    truth_gt = {}
    for p, recs in D.vcf_records(os.environ["PGDEV_TRUTH_VCF"], contig, 1, 2 ** 31 - 1).items():
        for ref, alts, gt, _, _ in recs:
            if len(gt) == 2 and None not in gt and gt[0] != gt[1]:
                truth_gt[p] = "ref/alt" if 0 in gt else "alt/alt"
    truth = sorted(truth_gt)
    het_tab = collections.defaultdict(lambda: [0, 0])    # key -> [right, wrong]
    nohet_tab = collections.defaultdict(lambda: [0, 0])  # key -> [majority, minority]
    n_het = n_nohet = 0
    het_sites = []
    for pos, rc in calls.items():
        run = repeat_run(fa, ref_name, pos)
        tol = max(3, run)
        i = bisect.bisect_left(truth, pos - tol)
        het = i < len(truth) and truth[i] <= pos + tol
        near = [truth[j] for j in range(i, min(i + 6, len(truth))) if truth[j] <= pos + tol]
        tg = ("2+ truth hets" if len(near) > 1 else truth_gt[near[0]] + (" same pos" if near[0] == pos else " offset")) if near else "none"
        rb = bin_of(RUN_BINS, run)
        if het:
            par = [0, 0, 0, 0]
            for q, (side, lab) in rc.items():
                par[(0 if scope[q][2] else 2) + side] += 1
            if not (par[0] + par[1] and par[2] + par[3]):
                continue
            ms, ps = (0 if par[0] >= par[1] else 1), (0 if par[2] >= par[3] else 1)
            if ms == ps:
                continue
            n_het += 1
            het_sites.append((pos, tol, rc))
        else:
            ones = sum(s for s, _ in rc.values())
            major = 1 if 2 * ones >= len(rc) else 0
            n_nohet += 1
        for q, (side, lab) in rc.items():
            s, e = scope[q][3], scope[q][4]
            dist = max(0, min(pos - s, e - (pos + run)))
            rerr, mq = extra[(pos, q)]
            keys = [("end", bin_of(END_BINS, dist)), ("label", "labelled" if lab else "unlabelled"), ("run", rb),
                    ("read error", f"{rerr:.4f}" if rerr in (0.0, 0.0474, 0.0025, 0.0001, 0.25, 0.35) else
                     ("< 0.01" if rerr < 0.01 else "other")),
                    ("mapq", "60" if mq >= 60 else "30-59" if mq >= 30 else "<30"), ("truth", tg),
                    ("end x label", f"{bin_of(END_BINS, dist)} {'labelled' if lab else 'unlabelled'}")]
            if het:
                right = side == (ms if scope[q][2] else ps)
                for k in keys:
                    het_tab[k][0 if right else 1] += 1
            else:
                for k in keys:
                    nohet_tab[k][0 if side == major else 1] += 1
    print(f"{contig}: {a.sites} sites with calls {len(calls):,}: truth het (parents split) {n_het:,}, no truth het {n_nohet:,}")
    for group in ("truth", "end", "label", "run", "read error", "mapq", "end x label"):
        keys = sorted({k for k in list(het_tab) + list(nohet_tab) if k[0] == group},
                      key=lambda k: (int(k[1].split(",")[0].strip("[")) if k[1][:1] == "[" else 0, k[1]))
        print(f"by {group}:   truth-het sites: calls, wrong %   |   no-het sites: calls, minority %")
        for k in keys:
            h, n = het_tab.get(k, [0, 0]), nohet_tab.get(k, [0, 0])
            hp = f"{100 * h[1] / sum(h):5.1f}" if sum(h) else "  n/a"
            npct = f"{100 * n[1] / sum(n):5.1f}" if sum(n) else "  n/a"
            print(f"  {k[1]:32} {sum(h):8,} {hp}%   |   {sum(n):8,} {npct}%")
    # within-site comparison: Simpson check for the labelled / unlabelled gap
    within = collections.Counter()
    nsite = 0
    per_site_err = []
    for pos, tol, rc in het_sites:
        par = [0, 0, 0, 0]
        for q, (side, lab) in rc.items():
            par[(0 if scope[q][2] else 2) + side] += 1
        ms, ps = (0 if par[0] >= par[1] else 1), (0 if par[2] >= par[3] else 1)
        c = collections.Counter()
        for q, (side, lab) in rc.items():
            c[(lab, side == (ms if scope[q][2] else ps))] += 1
        nl, nu = c[(True, True)] + c[(True, False)], c[(False, True)] + c[(False, False)]
        per_site_err.append(((c[(True, False)] + c[(False, False)]) / (nl + nu), nl, nu))
        if nl >= 3 and nu >= 3:
            nsite += 1
            for k, v in c.items():
                within[k] += v
    if nsite:
        wl = within[(True, False)] / max(1, within[(True, False)] + within[(True, True)])
        wu = within[(False, False)] / max(1, within[(False, False)] + within[(False, True)])
        print(f"truth-het sites with >= 3 labelled and >= 3 unlabelled calls: {nsite}; wrong % labelled "
              f"{100 * wl:.1f}, unlabelled {100 * wu:.1f}")
    print("truth-het sites by their own call error (all reads): sites, labelled calls, unlabelled calls")
    for lo, hi in ((0, 0.1), (0.1, 0.2), (0.2, 0.3), (0.3, 0.4), (0.4, 0.51)):
        sub = [x for x in per_site_err if lo <= x[0] < hi]
        print(f"  error [{lo}, {hi})  sites {len(sub):4}  labelled {sum(x[1] for x in sub):7,}  unlabelled {sum(x[2] for x in sub):6,}")
    if D.chunk_files(os.path.join(a.diag, "aln"), "alnreads"):
        aln_compare(a.diag, scope, het_sites)


def aln_compare(diag, scope, het_sites):
    """The same truth-het sites in the alignment solve: each read's allele there
    (nearest solve position within the site's tolerance with the most calls),
    judged by the parents' majorities, against the read's EM call."""
    want = {}
    for pos, tol, rc in het_sites:
        for p in range(pos - tol, pos + tol + 1):
            want.setdefault(p, []).append(pos)
    aln = collections.defaultdict(dict)  # solve pos -> read -> allele
    for f in D.chunk_files(os.path.join(diag, "aln"), "alnreads"):
        with open(f) as fh:
            for line in fh:
                p = line.rstrip("\n").split("\t")
                if p[0] not in scope:
                    continue
                for t in p[5:]:
                    x = t.split(":")
                    if len(x) < 3 or x[2] not in ("0", "1"):
                        continue
                    ap = int(x[0])
                    if ap in want and p[0] not in aln[ap]:
                        aln[ap][p[0]] = int(x[2])
    tab = collections.defaultdict(lambda: [0, 0])
    both = collections.Counter()
    for pos, tol, rc in het_sites:
        cands = [p for p in range(pos - tol, pos + tol + 1) if aln.get(p)]
        if not cands:
            tab[("no solve position", "")][0] += 0
            continue
        ap = max(cands, key=lambda p: len(aln[p]))
        par = [0, 0, 0, 0]
        for q, al in aln[ap].items():
            par[(0 if scope[q][2] else 2) + al] += 1
        if not (par[0] + par[1] and par[2] + par[3]):
            continue
        ms, ps = (0 if par[0] >= par[1] else 1), (0 if par[2] >= par[3] else 1)
        if ms == ps:
            continue
        epar = [0, 0, 0, 0]
        for q, (side, lab) in rc.items():
            epar[(0 if scope[q][2] else 2) + side] += 1
        ems, eps = (0 if epar[0] >= epar[1] else 1), (0 if epar[2] >= epar[3] else 1)
        for q, al in aln[ap].items():
            lab = rc[q][1] if q in rc else None
            right = al == (ms if scope[q][2] else ps)
            g = "EM labelled" if lab else "EM unlabelled" if lab is False else "no EM call"
            tab[(g, "")][0 if right else 1] += 1
            if q in rc:
                em_right = rc[q][0] == (ems if scope[q][2] else eps)
                both[(g, "solve " + ("right" if right else "wrong"), "EM " + ("right" if em_right else "wrong"))] += 1
    print("the same truth-het sites in the alignment solve: calls, wrong %")
    for k, v in sorted(tab.items()):
        if sum(v):
            print(f"  {k[0]:32} {sum(v):8,} {100 * v[1] / sum(v):5.1f}%")
    print("  solve call vs EM call, same read:")
    for k, v in sorted(both.items()):
        print(f"    {k[0]:16} {k[1]:12} {k[2]:10} {v:8,}")


if __name__ == "__main__":
    main()
