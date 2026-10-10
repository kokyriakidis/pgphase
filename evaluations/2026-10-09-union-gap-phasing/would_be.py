#!/usr/bin/env python3
"""How accurate would labels be for reads the EM leaves unlabelled? (eval only)

Usage: would_be.py DIAG_READS TRUTH_TSV
Each dump line is one read in one chunk after the final EM pass: label, and
per observation (site pos, category bits, injected, site error, own error,
side relative to hap1 (1 match / 2 mismatch), block). Block orientation per
(chunk, block) is the truth majority of the reads the EM labelled there. For
unlabelled reads, the would-be label is the sign of the summed evidence in
their strongest block, under several evidence rules.
"""
import collections
import math
import sys

diag_path, truth_path = sys.argv[1:3]
only = set(open(sys.argv[3]).read().split()) if len(sys.argv) > 3 else None
truth = {}
for line in open(truth_path):
    q, h = line.rstrip("\n").split("\t")[:2]
    if h in ("MATERNAL", "PATERNAL"):
        truth[q] = h

rows = []
orient_votes = collections.Counter()
for line in open(diag_path):
    f = line.rstrip("\n").split("\t")
    q, mapq, hap, chunk = f[0], int(f[1]), int(f[2]), f[4]
    obs = []
    for item in f[6:]:
        pos, cat, inj, err, oe, match, block = item.split(":")
        obs.append((int(pos), int(cat), int(inj), float(err), float(oe), int(match), (chunk, int(block))))
    if hap in (1, 2) and q in truth:
        blocks = collections.Counter(o[6] for o in obs if o[3] < 0.3)
        if blocks:
            b = blocks.most_common(1)[0][0]
            orient_votes[b] += 1 if (hap == 1) == (truth[q] == "MATERNAL") else -1
    rows.append((q, mapq, hap, obs))

labelled = {q for q, _, hap, _ in rows if hap in (1, 2)}


def evidence(obs, max_err):
    per = collections.defaultdict(float)
    for pos, cat, inj, err, oe, match, block in obs:
        if err >= max_err:
            continue
        e = min(0.45, 1 - (1 - err) * (1 - oe))
        per[block] += (1 if match == 1 else -1) * math.log((1 - e) / e)
    if not per:
        return None
    b, m = max(per.items(), key=lambda kv: abs(kv[1]))
    return b, m


def cls_of(obs):
    if not obs:
        return "no_obs"
    if all(o[3] >= 0.3 for o in obs):
        return "only_unreliable"
    return "low_posterior"


tab = collections.Counter()
seen = set()
for q, mapq, hap, obs in rows:
    if q in labelled or q not in truth or q in seen or (only is not None and q not in only):
        continue
    cls = cls_of(obs)
    for rule, max_err in (("reliable", 0.3), ("all_sites", 1.0)):
        ev = evidence(obs, max_err)
        if ev is None:
            continue
        b, m = ev
        o = orient_votes.get(b, 0)
        if o == 0:
            tab[(cls, rule, "block_unoriented")] += 1
            continue
        p = 1 / (1 + math.exp(-abs(m)))
        bucket = "p<0.6" if p < 0.6 else "p0.6-0.8" if p < 0.8 else "p0.8-0.9" if p < 0.9 else "p>=0.9"
        ok = ((m > 0) == (truth[q] == "MATERNAL")) == (o > 0)
        tab[(cls, rule, bucket, "correct" if ok else "wrong")] += 1
    seen.add(q)

for cls in ("low_posterior", "only_unreliable"):
    for rule in ("reliable", "all_sites"):
        print(f"== {cls} / {rule}  (unoriented blocks: {tab[(cls, rule, 'block_unoriented')]})")
        for bucket in ("p<0.6", "p0.6-0.8", "p0.8-0.9", "p>=0.9"):
            c, w = tab[(cls, rule, bucket, "correct")], tab[(cls, rule, bucket, "wrong")]
            if c + w:
                print(f"   {bucket:9s} correct {c:5d} wrong {w:5d}  ({100 * c / (c + w):.0f}%)")
