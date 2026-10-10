#!/usr/bin/env python3
"""Phased arm het records: which production-visible features separate truth
hets from non-truth hets? (evaluation only)
Usage: site_features.py RUN_DIR TRUTH_VCF DIAG_READS"""
import bisect
import collections
import gzip
import sys

run, truth_vcf, diag_path = sys.argv[1:4]
CEN = (26_000_000, 32_000_000)

thet = []
for line in gzip.open(truth_vcf, "rt"):
    if line[0] == "#":
        continue
    f = line.split("\t")
    gt = f[9].split(":")[0].replace("/", "|").split("|")
    if len(gt) == 2 and gt[0] != gt[1] and "." not in gt:
        thet.append(int(f[1]) - 1)
thet.sort()
tset = set(thet)


def is_truth(p, snp):
    if snp:
        return p in tset
    i = bisect.bisect_left(thet, p - 10)
    return i < len(thet) and thet[i] <= p + 10


# Learned site error and read count per site position (0-based sort pos ~ VCF pos) from the dump.
site_err = {}
site_reads = collections.Counter()
for line in open(diag_path):
    f = line.rstrip("\n").split("\t")
    for item in f[6:]:
        pos, cat, inj, err, oe, match, block = item.split(":")
        if cat == "-1":
            continue
        site_err[int(pos)] = float(err)
        site_reads[int(pos)] += 1

recs = []
for line in open(f"{run}/phased.vcf"):
    if line[0] == "#":
        continue
    f = line.rstrip("\n").split("\t")
    fmt = dict(zip(f[8].split(":"), f[9].split(":")))
    gt = fmt.get("GT", "")
    if "|" not in gt or fmt.get("PS", ".") in (".", "") or len(set(gt.split("|"))) < 2:
        continue
    p = int(f[1]) - 1
    if CEN[0] <= p < CEN[1]:
        continue
    info = dict(kv.split("=", 1) for kv in f[7].split(";") if "=" in kv)
    snp = len(f[3]) == 1 and all(len(a) == 1 for a in f[4].split(","))
    err = site_err.get(p + 1, site_err.get(p, site_err.get(p + 2)))
    recs.append(dict(p=p, ps=int(fmt["PS"]), snp=snp, cat=info.get("CAT", "?"), af=float(info.get("AF", 0)),
                     dp=int(info.get("DP", 0)), graph=f[2] != ".", err=err, truth=is_truth(p, snp)))
recs.sort(key=lambda r: r["p"])
# Production-visible linkage: phased neighbours in the same PS within 15 kb.
by_ps = collections.defaultdict(list)
for r in recs:
    by_ps[r["ps"]].append(r["p"])
for r in recs:
    ps = by_ps[r["ps"]]
    i = bisect.bisect_left(ps, r["p"])
    r["nbr15k"] = bisect.bisect_left(ps, r["p"] + 15000) - bisect.bisect_left(ps, r["p"] - 15000) - 1
    allp = [x["p"] for x in recs]
    r["block_size"] = len(ps)


def show(name, key):
    tab = collections.defaultdict(lambda: [0, 0])
    for r in recs:
        tab[key(r)][0 if r["truth"] else 1] += 1
    print(f"== {name}")
    for k in sorted(tab, key=str):
        t, n = tab[k]
        print(f"   {str(k):28s} truth {t:6d}  not {n:5d}  ({100 * n / (t + n):.1f}% not)")


show("source", lambda r: ("graph" if r["graph"] else "injected", "snp" if r["snp"] else "indel", r["cat"]))
show("phased neighbours within 15 kb", lambda r: min(r["nbr15k"], 5))
show("block size (records)", lambda r: "1" if r["block_size"] == 1 else "2-3" if r["block_size"] <= 3 else "4-10" if r["block_size"] <= 10 else ">10")
show("AF", lambda r: f"{min(int(r['af'] * 10), 9) / 10:.1f}")
show("learned err", lambda r: "?" if r["err"] is None else f"{min(int(r['err'] * 20), 9) / 20:.2f}")
show("injected & isolated (nbr15k==0)", lambda r: (not r["graph"], r["nbr15k"] == 0))
