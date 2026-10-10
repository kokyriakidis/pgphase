import collections
import sys

sys.path.insert(0, ".")
import audit_hybrid as A

tags = A.load_tags("full/af22/phased.bam")
rem = {q for q in open("hiphase_wins_af9.txt").read().split() if q not in tags}
seen = {}
for l in open("/tmp/cons.tsv"):
    q, n, a, b, s = l.rstrip("\n").split("\t")
    if q in rem:
        seen[q] = (int(n), int(a), int(b), int(s))
tab = collections.Counter()
for q in rem:
    if q not in seen:
        tab["not considered (no alignment in chunk, or not unlabelled there)"] += 1
        continue
    n, a, b, s = seen[q]
    if n == 0:
        tab["no markers covered"] += 1
        continue
    frac = max(a, b) / n
    nb = str(n) if n < 5 else "5+"
    fb = ">=0.9" if frac >= 0.9 else "0.7-0.9" if frac >= 0.7 else "<0.7"
    tab[f"votes {nb} agreement {fb}"] += 1
for k, v in tab.most_common():
    print(v, k)
