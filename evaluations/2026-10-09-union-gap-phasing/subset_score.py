#!/usr/bin/env python3
"""Score runs on a fixed read subset (per-PS majority orientation over all reads). Eval only.
Usage: subset_score.py QNAME_LIST NAME=DIR [...]"""
import collections, sys
sys.path.insert(0, ".")
import audit_hybrid as A
subset = set(open(sys.argv[1]).read().split())
truth = A.load_truth("../derived/chr20_truth_hap.tsv")
for name, d in (a.split("=", 1) for a in sys.argv[2:]):
    tags = A.load_tags(f"{d}/phased.bam")
    _, orient = A.score_reads(tags, truth, set(truth))
    c = collections.Counter()
    for q in subset:
        t = tags.get(q)
        if t is None: c["unphased"] += 1; continue
        ok = ((t[0] == 1) == (truth[q] == "MATERNAL")) == orient[t[1]]
        c["correct" if ok else "discordant"] += 1
    print(f"{name}\tcorrect {c['correct']}\tdiscordant {c['discordant']}\tunphased {c['unphased']}")
