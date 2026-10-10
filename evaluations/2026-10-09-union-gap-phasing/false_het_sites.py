#!/usr/bin/env python3
"""Run-phased heterozygous records by category: how many the truth confirms,
and how many sit inside reads that span no truth heterozygote. (eval only)
Usage: false_het_sites.py RUN_DIR TRUTH_VCF"""
import bisect
import collections
import gzip
import sys

import pysam

run, truth_vcf = sys.argv[1:3]
BAM = "/home/kokyriakidis/Downloads/pgphase/test_data/HG002_chr20_hifi_mapped_to_CHM13_chr20_annotated.bam"
CEN = (26_000_000, 32_000_000)

truth_het, truth_any = [], set()
for line in gzip.open(truth_vcf, "rt"):
    if line[0] == "#":
        continue
    f = line.split("\t")
    gt = f[9].split(":")[0].replace("/", "|").split("|")
    p = int(f[1]) - 1
    truth_any.add(p)
    if len(gt) == 2 and gt[0] != gt[1]:
        truth_het.append(p)
truth_het.sort()
het_set = set(truth_het)


def near_truth_het(p, tol):
    i = bisect.bisect_left(truth_het, p - tol)
    return i < len(truth_het) and truth_het[i] <= p + tol


recs = []
for line in open(f"{run}/phased.vcf"):
    if line[0] == "#":
        continue
    f = line.rstrip("\n").split("\t")
    fmt = dict(zip(f[8].split(":"), f[9].split(":")))
    gt = fmt.get("GT", "")
    if "|" not in gt or fmt.get("PS", ".") in (".", "") or len(set(gt.split("|"))) < 2:
        continue
    info = dict(kv.split("=", 1) for kv in f[7].split(";") if "=" in kv)
    p = int(f[1]) - 1
    if CEN[0] <= p < CEN[1]:
        continue
    snp = len(f[3]) == 1 and len(f[4]) == 1
    recs.append((p, info.get("CAT", "?"), "snp" if snp else "indel", "ID" if f[2] != "." else "noID"))
recs.sort()
rpos = [r[0] for r in recs]

# Records inside arm reads that span no truth het.
in_nohet = set()
with pysam.AlignmentFile(BAM) as bam:
    for r in bam.fetch(until_eof=True):
        if r.is_unmapped or r.is_secondary or r.is_supplementary or r.is_duplicate:
            continue
        lo, hi = r.reference_start, r.reference_end
        if CEN[0] <= (lo + hi) // 2 < CEN[1]:
            continue
        if bisect.bisect_left(truth_het, hi) != bisect.bisect_left(truth_het, lo):
            continue
        for i in range(bisect.bisect_left(rpos, lo), bisect.bisect_left(rpos, hi)):
            in_nohet.add(i)

tab = collections.Counter()
for i, (p, cat, kind, rid) in enumerate(recs):
    key = (cat, kind)
    tab[key + ("total",)] += 1
    ok = p in het_set if kind == "snp" else near_truth_het(p, 10)
    tab[key + ("truth_het" if ok else "not_truth_het",)] += 1
    if i in in_nohet:
        tab[key + ("in_nohet_reads",)] += 1
print("category\tkind\ttotal\ttruth_het\tnot_truth_het\tin_nohet_reads")
for cat, kind in sorted({(k[0], k[1]) for k in tab}):
    g = lambda x: tab[(cat, kind, x)]
    print(f"{cat}\t{kind}\t{g('total')}\t{g('truth_het')}\t{g('not_truth_het')}\t{g('in_nohet_reads')}")
