import collections, sys
sys.path.insert(0, ".")
import audit_hybrid as A
subset = set(open("full/clean5/phaseable_unphased.txt").read().split())
cls = {}
for line in open("full/diag/reads.tsv"):
    f = line.rstrip("\n").split("\t")
    q = f[0]
    if q not in subset: continue
    obs = [x.split(":") for x in f[6:]]
    c = "no_obs" if not obs else "only_unreliable" if all(float(o[3]) >= 0.3 for o in obs) else "low_posterior"
    rank = {"absent": 0, "no_obs": 1, "only_unreliable": 2, "low_posterior": 3}
    if q not in cls or rank[c] > rank[cls[q]]: cls[q] = c
for q in subset: cls.setdefault(q, "absent")
truth = A.load_truth("../derived/chr20_truth_hap.tsv")
for name, d in (a.split("=", 1) for a in sys.argv[1:]):
    tags = A.load_tags(f"{d}/phased.bam")
    _, orient = A.score_reads(tags, truth, set(truth))
    c = collections.Counter()
    for q, k in cls.items():
        t = tags.get(q)
        if t is None: c[(k, "unphased")] += 1; continue
        ok = ((t[0] == 1) == (truth[q] == "MATERNAL")) == orient[t[1]]
        c[(k, "correct" if ok else "wrong")] += 1
    print(name)
    for k in ("low_posterior", "only_unreliable", "no_obs", "absent"):
        print(f"   {k:16s} correct {c[(k,'correct')]:5d} wrong {c[(k,'wrong')]:4d} unphased {c[(k,'unphased')]:5d}")
