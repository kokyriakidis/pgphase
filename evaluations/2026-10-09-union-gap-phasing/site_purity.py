"""Unreliable injected sites behind the remaining reads: do the EM's calls
there separate maternal from paternal reads? (evaluation only)"""
import collections
import sys

sys.path.insert(0, ".")
import audit_hybrid as A

run, dump_path, qlist = sys.argv[1:4]
truth = A.load_truth("../derived/chr20_truth_hap.tsv")
tags = A.load_tags(f"{run}/phased.bam")
targets = {q for q in open(qlist).read().split() if q not in tags}
target_sites = set()
site_calls = collections.defaultdict(collections.Counter)  # (chunk, pos) -> (side-vs-phase, parent)
site_err = {}
for line in open(dump_path):
    f = line.rstrip("\n").split("\t")
    q = f[0]
    for o in (x.split(":") for x in f[6:]):
        if o[1] == "-1" or o[2] != "1":
            continue
        key = (f[4], int(o[0]))
        site_err[key] = float(o[3])
        if q in truth:
            site_calls[key][(o[5], truth[q])] += 1
        if q in targets and float(o[3]) >= 0.3:
            target_sites.add(key)
buckets = collections.Counter()
reads_behind = collections.Counter()
for key in target_sites:
    c = site_calls[key]
    cis = c[("1", "MATERNAL")] + c[("2", "PATERNAL")]
    trans = c[("2", "MATERNAL")] + c[("1", "PATERNAL")]
    n = cis + trans
    if n == 0:
        buckets["no truth-labelled calls"] += 1
        continue
    pur = max(cis, trans) / n
    buckets["purity >=0.9" if pur >= 0.9 else "0.75-0.9" if pur >= 0.75 else "0.6-0.75" if pur >= 0.6 else "<0.6"] += 1
print("unreliable injected sites behind remaining reads:", len(target_sites))
for k, v in buckets.most_common():
    print(v, k)
