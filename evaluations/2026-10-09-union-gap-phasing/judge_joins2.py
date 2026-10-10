#!/usr/bin/env python3
"""Joins run A makes that run B cuts: was A's join orientation right? (eval only)
Usage: judge_joins.py RUN_A RUN_B"""
import bisect
import collections
import sys

import pysam

sys.path.insert(0, ".")
import audit_hybrid as A

run_a, run_b = sys.argv[1:3]
truth = A.load_truth("../derived/chr20_truth_hap.tsv")


recs = {}


def blocks(path):
    ps_of = []
    for line in open(path):
        if line[0] == "#":
            continue
        f = line.rstrip("\n").split("\t")
        fmt = dict(zip(f[8].split(":"), f[9].split(":")))
        if "|" in fmt.get("GT", "") and fmt.get("PS", ".") not in (".", ""):
            ps_of.append((int(f[1]), int(fmt["PS"])))
            if path.startswith(run_a):
                recs[int(f[1])] = (f[3][:8], f[4][:12], fmt.get("AD", "?"), fmt.get("DP", "?"))
    ps_of.sort()
    return ps_of


ba = blocks(f"{run_a}/phased.vcf")
bb = blocks(f"{run_b}/phased.vcf")
b_ps = dict(bb)
# Boundaries: consecutive records in the same A block but different B blocks.
cuts = []
for (p1, a1), (p2, a2) in zip(ba, ba[1:]):
    if a1 == a2 and p1 in b_ps and p2 in b_ps and b_ps[p1] != b_ps[p2]:
        cuts.append((p1, p2, a1))
tags = A.load_tags(f"{run_a}/phased.bam")
bam = pysam.AlignmentFile("../HG002_chr20_hifi_mapped_to_CHM13_chr20_annotated.bam")
verdict = collections.Counter()
reads_at_stake = collections.Counter()
for p1, p2, ps in cuts:
    side = {"L": collections.Counter(), "R": collections.Counter()}
    for lo, hi, key in ((p1 - 30000, p1, "L"), (p2, p2 + 30000, "R")):
        for r in bam.fetch("CHM13#0#chr20", max(0, lo), hi):
            q = r.query_name
            if r.is_secondary or r.is_supplementary or q not in truth or q not in tags:
                continue
            t = tags[q]
            if t[1] != ps:
                continue
            side[key][(t[0] == 1) == (truth[q] == "MATERNAL")] += 1
    if sum(side["L"].values()) < 5 or sum(side["R"].values()) < 5:
        verdict["undetermined"] += 1
        continue
    left = side["L"][True] >= side["L"][False]
    right = side["R"][True] >= side["R"][False]
    v = "join right" if left == right else "join WRONG"
    verdict[v] += 1
    print(f"{v:10s} gap {p2 - p1:7d} p1 {p1} p2 {p2}  L {recs.get(p1)}  R {recs.get(p2)}")
    reads_at_stake[v] += sum(side["R"].values())
print("joins in A cut in B:", len(cuts))
for k, v in verdict.most_common():
    print(f"  {v:5d} {k}   (reads right of the cut, 30 kb: {reads_at_stake[k]})")
