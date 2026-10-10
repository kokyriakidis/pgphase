#!/usr/bin/env python3
"""Offline what-if: relabel reads from the final-EM dump with candidate site
filters, score arm reads split into phaseable / no-het. (evaluation only)

Usage: sim_filters.py RUN_DIR DIAG_READS TRUTH_TSV TRUTH_VCF
Approximation: site errors, phases and blocks are the run's; a filter only
removes a site's evidence from read labelling (the real fix also removes it
from learning). Block orientation per (chunk, block) comes from the truth
majority of the reads the run labelled there.
"""
import bisect
import collections
import gzip
import math
import sys

import pysam

run, diag_path, truth_path, truth_vcf = sys.argv[1:5]
BAM = "/home/kokyriakidis/Downloads/pgphase/test_data/HG002_chr20_hifi_mapped_to_CHM13_chr20_annotated.bam"
CEN = (26_000_000, 32_000_000)

truth = {}
for line in open(truth_path):
    q, h = line.rstrip("\n").split("\t")[:2]
    if h in ("MATERNAL", "PATERNAL"):
        truth[q] = h
thet = []
for line in gzip.open(truth_vcf, "rt"):
    if line[0] == "#":
        continue
    f = line.split("\t")
    gt = f[9].split(":")[0].replace("/", "|").split("|")
    if len(gt) == 2 and gt[0] != gt[1] and "." not in gt:
        thet.append(int(f[1]) - 1)
thet.sort()

read_class = {}
with pysam.AlignmentFile(BAM) as bam:
    for r in bam.fetch(until_eof=True):
        if r.is_unmapped or r.is_secondary or r.is_supplementary or r.is_duplicate or r.query_name not in truth:
            continue
        lo, hi = r.reference_start, r.reference_end
        if CEN[0] <= (lo + hi) // 2 < CEN[1]:
            continue
        read_class[r.query_name] = "phaseable" if bisect.bisect_left(thet, hi) > bisect.bisect_left(thet, lo) else "nohet"

# Site features from the run VCF, keyed by 1-based POS (= the dump's site position).
site = {}
for line in open(f"{run}/phased.vcf"):
    if line[0] == "#":
        continue
    f = line.rstrip("\n").split("\t")
    info = dict(kv.split("=", 1) for kv in f[7].split(";") if "=" in kv)
    fmt = dict(zip(f[8].split(":"), f[9].split(":")))
    snp = len(f[3]) == 1 and all(len(a) == 1 for a in f[4].split(","))
    site[int(f[1])] = dict(af=float(info.get("AF", 0.5)), cat=info.get("CAT", "?"), snp=snp,
                           ps=fmt.get("PS", "."))
by_ps = collections.defaultdict(list)
for p, s in site.items():
    if s["ps"] not in (".", ""):
        by_ps[s["ps"]].append(p)
for v in by_ps.values():
    v.sort()
for p, s in site.items():
    ps = by_ps.get(s["ps"])
    s["nbr"] = (bisect.bisect_left(ps, p + 15000) - bisect.bisect_left(ps, p - 15000) - 1) if ps else 0

rows = []
orient = collections.Counter()
for line in open(diag_path):
    f = line.rstrip("\n").split("\t")
    q, hap, chunk = f[0], int(f[2]), f[4]
    obs = []
    for item in f[6:]:
        pos, cat, inj, err, oe, match, block = item.split(":")
        obs.append((int(pos), int(cat), float(err), float(oe), int(match), (chunk, int(block))))
    rows.append((q, hap, obs))
    if hap in (1, 2) and q in truth:
        blocks = collections.Counter(o[5] for o in obs if o[2] < 0.3)
        if blocks:
            orient[blocks.most_common(1)[0][0]] += 1 if (hap == 1) == (truth[q] == "MATERNAL") else -1


def simulate(keep, unreliable=0.3, posterior=0.9):
    best = {}
    for q, hap, obs in rows:
        if q not in read_class:
            continue
        per = collections.defaultdict(float)
        for pos, cat, err, oe, match, block in obs:
            if err >= unreliable or not keep(pos, cat, err):
                continue
            e = min(0.45, 1 - (1 - err) * (1 - oe))
            per[block] += (1 if match == 1 else -1) * math.log((1 - e) / e)
        if not per:
            continue
        b, m = max(per.items(), key=lambda kv: abs(kv[1]))
        if 1 / (1 + math.exp(-abs(m))) < posterior:
            continue
        if q not in best or abs(m) > abs(best[q][1]):
            best[q] = (b, m)
    tab = collections.Counter()
    for q, (b, m) in best.items():
        o = orient.get(b, 0)
        if o == 0:
            tab[(read_class[q], "unoriented")] += 1
            continue
        ok = ((m > 0) == (truth[q] == "MATERNAL")) == (o > 0)
        tab[(read_class[q], "correct" if ok else "wrong")] += 1
    return tab


def kept_all(pos, cat, err):
    return True


def feat(pos):
    return site.get(pos)


FILTERS = {
    "baseline": kept_all,
    "af_0.25-0.75": lambda p, c, e: (s := feat(p)) is None or 0.25 <= s["af"] <= 0.75,
    "af_0.3-0.7": lambda p, c, e: (s := feat(p)) is None or 0.3 <= s["af"] <= 0.7,
    "nbr>=1": lambda p, c, e: (s := feat(p)) is None or s["nbr"] >= 1,
    "nbr>=2": lambda p, c, e: (s := feat(p)) is None or s["nbr"] >= 2,
    "no_noisy_snp": lambda p, c, e: (s := feat(p)) is None or not (s["snp"] and s["cat"] == "NOISY_CAND_HET"),
    "af_0.25-0.75+nbr>=1": lambda p, c, e: (s := feat(p)) is None or (0.25 <= s["af"] <= 0.75 and s["nbr"] >= 1),
    "af_0.25-0.75+nbr>=1+no_noisy_snp": lambda p, c, e: (s := feat(p)) is None or (
        0.25 <= s["af"] <= 0.75 and s["nbr"] >= 1 and not (s["snp"] and s["cat"] == "NOISY_CAND_HET")),
}
print("filter\tphaseable_correct\tphaseable_wrong\tnohet_correct\tnohet_wrong\tunoriented")
for name, keep in FILTERS.items():
    t = simulate(keep)
    print(f"{name}\t{t[('phaseable', 'correct')]}\t{t[('phaseable', 'wrong')]}\t{t[('nohet', 'correct')]}\t"
          f"{t[('nohet', 'wrong')]}\t{t[('phaseable', 'unoriented')] + t[('nohet', 'unoriented')]}")
for thr in (0.25, 0.2):
    t = simulate(kept_all, unreliable=thr)
    print(f"unreliable<{thr}\t{t[('phaseable', 'correct')]}\t{t[('phaseable', 'wrong')]}\t{t[('nohet', 'correct')]}\t{t[('nohet', 'wrong')]}\t-")
print("posterior sweeps")
af = FILTERS["af_0.25-0.75"]
for name, keep in (("baseline", kept_all), ("af_0.25-0.75", af)):
    for post in (0.9, 0.85, 0.8, 0.75, 0.7):
        t = simulate(keep, posterior=post)
        pc, pw = t[("phaseable", "correct")], t[("phaseable", "wrong")]
        nc, nw = t[("nohet", "correct")], t[("nohet", "wrong")]
        print(f"{name}\tpost>={post}\t{pc}\t{pw}\t{nc}\t{nw}")
