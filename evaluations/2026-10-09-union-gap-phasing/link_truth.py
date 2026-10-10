#!/usr/bin/env python3
"""Graph-only links a run breaks or flips, judged against the phased truth (eval only).
Usage: link_truth.py INV_DIR TRUTH_VCF RUN_DIR [RUN_DIR ...]"""
import collections, gzip, sys
sys.path.insert(0, sys.argv[1])
import audit_hybrid as A

def phased(path):
    out = {}
    op = gzip.open if path.endswith(".gz") else open
    for line in op(path, "rt"):
        if line[0] == "#": continue
        f = line.rstrip("\n").split("\t")
        fmt = dict(zip(f[8].split(":"), f[9].split(":")))
        gt = fmt.get("GT", "")
        if "|" not in gt: continue
        a, b = gt.split("|")
        if a == b or "." in (a, b): continue
        alleles = [f[3]] + f[4].split(",")
        ps = fmt.get("PS", "0")
        out[int(f[1])] = (alleles[int(a)], alleles[int(b)], ps)
    return out

inv, truth_vcf = sys.argv[1:3]
truth = phased(truth_vcf)
graph = phased(f"{inv}/graph/phased.vcf")
gpos = sorted(graph)
links = [(a, b) for a, b in zip(gpos, gpos[1:]) if graph[a][2] == graph[b][2]]

def rel(call, p):  # True if hap1 carries the same allele as truth hap1
    t = truth.get(p)
    if t is None or set(t[:2]) != set(call[:2]): return None
    return call[0] == t[0]

def graph_ok(a, b):
    x, y = rel(graph[a], a), rel(graph[b], b)
    return None if x is None or y is None else x == y

for d in sys.argv[3:]:
    run = phased(f"{d}/phased.vcf")
    ph, _ = A.load_vcf(f"{d}/phased.vcf")
    cover = A.covering(A.blocks_of(ph))
    tab = collections.Counter()
    for a, b in links:
        cen = "cen" if 26_000_000 <= a < 32_000_000 else "arms"
        g = graph_ok(a, b)
        gs = "graph_right" if g else "graph_wrong" if g is not None else "no_truth"
        if cover(a, b) is None:
            tab[(cen, "broken", gs)] += 1
        elif a in run and b in run and run[a][2] == run[b][2]:
            same = (run[a][0] == graph[a][0]) == (run[b][0] == graph[b][0])
            if not same: tab[(cen, "flipped", gs)] += 1
    print(d)
    for k in sorted(tab): print("  ", *k, tab[k])
