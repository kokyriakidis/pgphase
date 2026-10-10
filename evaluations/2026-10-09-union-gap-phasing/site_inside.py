"""Inside the EM at one site: reads observing it, their labels, sides, truth."""
import collections
import sys

sys.path.insert(0, ".")
import audit_hybrid as A

dump_path, site = sys.argv[1], int(sys.argv[2])
truth = A.load_truth("../derived/chr20_truth_hap.tsv")
tab = collections.Counter()
blocks = collections.Counter()
other_obs = collections.Counter()
err = None
for line in open(dump_path):
    f = line.rstrip("\n").split("\t")
    obs = [x.split(":") for x in f[6:]]
    at = [o for o in obs if int(o[0]) == site]
    if not at:
        continue
    o = at[0]
    err = o[3]
    lab = f[2]
    t = truth.get(f[0], "?")[0]
    tab[(f"label {lab}", f"side {o[5]}", f"truth {t}")] += 1
    blocks[o[6]] += 1
    rel_other = sum(1 for x in obs if int(x[0]) != site and float(x[3]) < 0.3)
    other_obs["other reliable obs " + ("0" if rel_other == 0 else "1-2" if rel_other < 3 else "3+")] += 1
print("site", site, "learned err", err, "blocks", dict(blocks))
print(dict(other_obs))
for k, v in sorted(tab.items()):
    print(v, k)
