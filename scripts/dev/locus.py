#!/usr/bin/env python3
"""EM sites in a window: truth purity, learned error, block, and how cleanly
each site's calls agree with nearby sites' calls on shared reads (eval only).

Usage: locus.py EM_DUMP_DIR BEG END [--pairs K]
Pair concordance is orientation-free: max(same, opposite) / shared reads.
"""
import collections
import itertools
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import score  # noqa: E402
import why  # noqa: E402


def main():
    dump_dir, beg, end = sys.argv[1], int(sys.argv[2]), int(sys.argv[3])
    k_near = int(sys.argv[sys.argv.index("--pairs") + 1]) if "--pairs" in sys.argv else 3
    dump = why.load_dump(dump_dir)
    scope = score.read_scope()
    calls = collections.defaultdict(dict)  # site -> read -> side
    meta = {}
    for q, d in dump.items():
        for pos, kind, side, err, blk in d["obs"]:
            if beg <= pos <= end:
                calls[pos][q] = side
                meta[pos] = (kind, err, blk)
    sites = sorted(calls)
    print(f"{'pos':>10} kind  err   block       reads  purity  labelled(ours)  | concordance with next sites (n)")
    for i, p in enumerate(sites):
        pur = [0, 0]
        for q, side in calls[p].items():
            if q in scope:
                pur[(side == 0) == scope[q][2]] += 1
        purity = max(pur) / sum(pur) if sum(pur) else float("nan")
        labelled = sum(1 for q in calls[p] if dump[q]["hap"])
        pairs = []
        for p2 in sites[i + 1:i + 1 + k_near]:
            shared = calls[p].keys() & calls[p2].keys()
            if len(shared) < 3:
                continue
            same = sum(1 for q in shared if calls[p][q] == calls[p2][q])
            pairs.append(f"{p2}:{max(same, len(shared) - same) / len(shared):.2f}({len(shared)})")
        kind, err, blk = meta[p]
        print(f"{p:>10} {kind:4} {err:5.3f} {blk:>10} {len(calls[p]):6} {purity:7.2f} {labelled:8}        | {' '.join(pairs)}")


if __name__ == "__main__":
    main()
