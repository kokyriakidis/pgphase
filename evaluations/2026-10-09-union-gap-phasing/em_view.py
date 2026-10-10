"""Remaining HiPhase-correct unphased reads: what the final EM holds for them."""
import collections
import math
import sys

sys.path.insert(0, ".")
import audit_hybrid as A

run, dump_path, qlist = sys.argv[1:4]
tags = A.load_tags(f"{run}/phased.bam")
reads = {q for q in open(qlist).read().split() if q not in tags}
best = {}
for line in open(dump_path):
    f = line.rstrip("\n").split("\t")
    if f[0] not in reads:
        continue
    obs = [x.split(":") for x in f[6:]]
    rel = [o for o in obs if float(o[3]) < 0.3]
    per = collections.defaultdict(float)
    for o in rel:
        e = min(0.45, 1 - (1 - float(o[3])) * (1 - float(o[4])))
        per[o[6]] += (1 if o[5] == "1" else -1) * math.log((1 - e) / e)
    m = max((abs(v) for v in per.values()), default=0.0)
    kinds = collections.Counter("window" if o[1] == "-1" else ("inj" if o[2] == "1" else "graph") for o in obs)
    rec = (len(obs), len(rel), m, dict(kinds))
    if f[0] not in best or rec[0] > best[f[0]][0]:
        best[f[0]] = rec
tab = collections.Counter()
for q in reads:
    if q not in best:
        tab["not in any chunk"] += 1
        continue
    n, r, m, kinds = best[q]
    if n == 0:
        cls = "in chunk, no observations"
    elif r == 0:
        cls = "only unreliable observations (" + "+".join(sorted(kinds)) + ")"
    else:
        cls = "reliable but posterior < 0.8 (margin " + ("<0.7" if m < 0.7 else "0.7-1.39") + ")"
    tab[cls] += 1
for k, v in tab.most_common():
    print(v, k)
