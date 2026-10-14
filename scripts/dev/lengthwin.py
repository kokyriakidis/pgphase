#!/usr/bin/env python3
"""Length windows in a diag run against truth (evaluation only).

Usage: lengthwin.py DIAG_DIR

Per window (lengthwin.* dump): a truth het within it, and the parental purity
of its read sides (share of reads on their parent's majority side). Tallied by
het score and weak flag.
"""
import bisect, collections, os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import diaglib as D  # noqa: E402
import score  # noqa: E402

diag = sys.argv[1]
contig = os.environ["PGDEV_CONTIG"]
scope = score.read_scope()
truth = sorted(p for p, recs in D.vcf_records(os.environ["PGDEV_TRUTH_VCF"], contig, 1, 2 ** 31 - 1).items()
               if any(len(gt) == 2 and None not in gt and gt[0] != gt[1] for _, _, gt, _, _ in recs))
rows = []
seen = set()
for f in D.chunk_files(os.path.join(diag, "aln"), "lengthwin"):
    for line in open(f):
        p = line.rstrip("\n").split("\t")
        key = (p[0], p[1])
        if key in seen:
            continue
        seen.add(key)
        pos, end, llr, weak = int(p[0]), int(p[1]), float(p[6]), p[9] == "1"
        par = [0, 0, 0, 0]
        for t in p[10:]:
            q, side = t.rsplit(":", 1)
            if q in scope:
                par[(0 if scope[q][2] else 2) + int(side)] += 1
        pur = None
        if par[0] + par[1] and par[2] + par[3]:
            ms, ps = (0 if par[0] >= par[1] else 1), (0 if par[2] >= par[3] else 1)
            pur = 0.5 if ms == ps else (par[ms] + par[2 + ps]) / sum(par)
        i = bisect.bisect_left(truth, pos - 1)
        het = i < len(truth) and truth[i] <= end + 1
        rows.append((llr, weak, het, pur, sum(par)))
print(f"{contig}: length windows {len(rows):,}; truth het inside {sum(r[2] for r in rows):,}")
for lo, hi in ((6, 10), (10, 20), (20, 50), (50, 1e9)):
    for weak in (False, True):
        sub = [r for r in rows if lo <= r[0] < hi and r[1] == weak]
        if not sub:
            continue
        pur = sorted(r[3] for r in sub if r[3] is not None)
        good = sum(1 for r in sub if r[3] is not None and r[3] >= 0.9)
        print(f"  score [{lo:g}, {hi:g}) {'weak' if weak else 'strong'}: windows {len(sub):5,}  truth het {sum(r[2] for r in sub):5,}"
              f"  purity median {pur[len(pur) // 2] if pur else 0:.2f}  purity>=0.9 {good:5,}")
