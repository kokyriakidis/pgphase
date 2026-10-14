#!/usr/bin/env python3
"""How many haplotypes do the reads at a region carry? (evaluation only)

Usage: copies.py RUN_DIR REGION [--min-af F] [--k K] [--min-shared N]

Every primary read overlapping REGION (any MAPQ) is piled up; a site is every
SNP position where the minor base reaches --min-af of the reads. Reads are
clustered on their bases at those sites: distance is the mismatch rate over
shared sites (at least --min-shared), clusters by k-medoids for K = 2..--k.
Per cluster (each K): reads, mean depth contribution, MAPQ, parents (truth,
evaluation only), our class for the read, and how many reads reach past the
region into flank sequence. Sites are also listed by how the parents split
them (real het / not) and by cluster.
"""
import argparse
import collections
import os
import random
import sys

import pysam

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import diaglib as D  # noqa: E402
import score  # noqa: E402


def kmedoids(n, dist, k, rng, rounds=30):
    """k-medoids on a precomputed n x n distance matrix; returns labels, cost."""
    import numpy as np
    med = rng.sample(range(n), k)
    for _ in range(rounds):
        lab = np.argmin(dist[:, med], axis=1)
        new = []
        for i in range(k):
            idx = np.where(lab == i)[0]
            if len(idx) == 0:
                new.append(med[i])
                continue
            sub = dist[np.ix_(idx, idx)].sum(axis=1)
            new.append(int(idx[np.argmin(sub)]))
        if new == med:
            break
        med = new
    lab = np.argmin(dist[:, med], axis=1)
    return lab, float(dist[np.arange(n), np.array(med)[lab]].sum())


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("run")
    ap.add_argument("region")
    ap.add_argument("--min-af", type=float, default=0.1)
    ap.add_argument("--k", type=int, default=5)
    ap.add_argument("--min-shared", type=int, default=5)
    a = ap.parse_args()
    contig, beg, end = D.parse_region(a.region)
    bam = pysam.AlignmentFile(os.environ["PGDEV_BAM"])
    bc = contig if contig in bam.references else contig.split("#")[-1]
    scope = score.read_scope()
    classes = score.read_classes(a.run)
    calls = collections.defaultdict(dict)  # read -> pos -> base
    mapq, span = {}, {}
    sites = []
    seqs = {}
    for col in bam.pileup(bc, beg - 1, end, truncate=True, min_base_quality=0, stepper="nofilter",
                          max_depth=100000):
        bases = {}
        for r in col.pileups:
            al = r.alignment
            if al.is_secondary or al.is_supplementary or r.is_del or r.is_refskip:
                continue
            q = al.query_name
            if q not in seqs:
                seqs[q] = al.query_sequence
            bases[q] = seqs[q][r.query_position]
            mapq[al.query_name] = al.mapping_quality
            span[al.query_name] = (al.reference_start, al.reference_end)
        if len(bases) < 6:
            continue
        cnt = collections.Counter(bases.values())
        if len(cnt) < 2:
            continue
        top = cnt.most_common(2)
        if top[1][1] < a.min_af * len(bases):
            continue
        pos = col.reference_pos + 1
        sites.append((pos, top[0][0], top[1][0]))
        for q, b in bases.items():
            if b in (top[0][0], top[1][0]):
                calls[q][pos] = 0 if b == top[0][0] else 1
    import numpy as np
    reads = {q: c for q, c in calls.items() if len(c) >= a.min_shared}
    names = sorted(reads)
    col = {pos: j for j, (pos, b0, b1) in enumerate(sites)}
    M = np.full((len(names), len(sites)), -1, dtype=np.int8)
    for i, q in enumerate(names):
        for pos, v in reads[q].items():
            M[i, col[pos]] = v
    has = (M >= 0).astype(np.float32)
    one = (M == 1).astype(np.float32)
    zero = (M == 0).astype(np.float32)
    shared = has @ has.T
    mism = one @ zero.T + zero @ one.T
    dist = np.where(shared >= a.min_shared, mism / np.maximum(shared, 1), 0.5)
    print(f"{contig}:{beg:,}-{end:,}: sites (minor >= {a.min_af}) {len(sites):,}; reads with >= {a.min_shared} "
          f"calls {len(reads):,} (MAPQ 0: {sum(mapq[q] == 0 for q in reads)}, <30: {sum(mapq[q] < 30 for q in reads)})",
          flush=True)
    # sites: parental split
    real = 0
    for pos, b0, b1 in sites:
        par = [0, 0, 0, 0]
        for q, c in calls.items():
            if pos in c and q in scope:
                par[(0 if scope[q][2] else 2) + c[pos]] += 1
        if par[0] + par[1] >= 3 and par[2] + par[3] >= 3:
            ms, ps = (0 if par[0] >= par[1] else 1), (0 if par[2] >= par[3] else 1)
            if ms != ps and (par[ms] + par[2 + ps]) / sum(par) >= 0.9:
                real += 1
    print(f"  sites splitting the parents at purity >= 0.9: {real} of {len(sites)}", flush=True)
    rng = random.Random(1)
    for k in range(2, a.k + 1):
        best = None
        for _ in range(5):
            lab, cost = kmedoids(len(names), dist, k, rng)
            if best is None or cost < best[1]:
                best = (lab, cost)
        lab, cost = best
        print(f"  K={k}: mean distance to medoid {cost / max(1, len(names)):.3f}", flush=True)
        for i in sorted(set(lab.tolist()), key=lambda i: -int((lab == i).sum())):
            g = [names[j] for j in np.where(lab == i)[0]]
            mat = sum(1 for q in g if q in scope and scope[q][2])
            pat = sum(1 for q in g if q in scope and not scope[q][2])
            cl = collections.Counter(classes.get(q, "-") for q in g)
            flank = sum(1 for q in g if span[q][0] < beg - 2000 or span[q][1] > end + 2000)
            mq = sorted(mapq[q] for q in g)
            print(f"    cluster {i}: reads {len(g):4}  maternal {mat:4}  paternal {pat:4}  "
                  f"MAPQ median {mq[len(mq) // 2]:2}  <30 {sum(m < 30 for m in mq):4}  "
                  f"reach flank {flank:4}  ours {dict(cl)}", flush=True)


if __name__ == "__main__":
    main()
