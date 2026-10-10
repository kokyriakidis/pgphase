#!/usr/bin/env python3
"""Why does union never call these truth indels? (evaluation only)
Usage: missing_indels.py RUN_DIR MISSING_TSV SITES_VCF"""
import bisect
import collections
import gzip
import sys

import pysam

run, missing_path, sites_vcf = sys.argv[1:4]
ROOT = "/home/kokyriakidis/Downloads/pgphase/test_data"
BAM = f"{ROOT}/HG002_chr20_hifi_mapped_to_CHM13_chr20_annotated.bam"
CONTIG = "CHM13#0#chr20"
TOL = 25


def load_tsv(path):
    rows = []
    with open(path) as fh:
        next(fh)
        for line in fh:
            f = line.rstrip("\n").split("\t")
            rows.append((int(f[1]), f[3], f[4], f[14], f[13]))
    rows.sort()
    return rows


bam_cand = load_tsv("bam/candidates.tsv")
run_cand = load_tsv(f"{run}/candidates.tsv")
bpos = [r[0] for r in bam_cand]
rpos = [r[0] for r in run_cand]
sites = []
for line in gzip.open(sites_vcf, "rt"):
    if line[0] == "#":
        continue
    f = line.split("\t", 6)
    sites.append((int(f[1]), int(f[1]) + len(f[3]), f[2]))
sites.sort()
spos = [s[0] for s in sites]


def near(rows, pos, p, lo, hi):
    return rows[bisect.bisect_left(pos, p - lo):bisect.bisect_right(pos, p + hi)]


bam = pysam.AlignmentFile(BAM)
summary = collections.Counter()
for line in open(missing_path):
    p, ref, alt = line.rstrip("\n").split("\t")
    p = int(p)
    dlen = len(alt) - len(ref)
    # Read support: primary reads with an indel of the same length change within +-10 bp.
    depth = support = support_mq30 = 0
    for r in bam.fetch(CONTIG, p - 1, p):
        if r.is_secondary or r.is_supplementary or r.is_duplicate:
            continue
        depth += 1
        rp = r.reference_start
        found = False
        for op, n in r.cigartuples:
            if op in (0, 7, 8):
                rp += n
            elif op == 2:
                if -n == dlen and abs(rp - p) <= 10:
                    found = True
                rp += n
            elif op == 1:
                if n == dlen and abs(rp - p) <= 10:
                    found = True
            if rp > p + 20:
                break
        if found:
            support += 1
            support_mq30 += r.mapping_quality >= 30
    bc = near(bam_cand, bpos, p, TOL, TOL)
    rc = near(run_cand, rpos, p, TOL, TOL)
    sc = [s for s in sites[max(0, bisect.bisect_left(spos, p - 2000)):bisect.bisect_right(spos, p + TOL)]
          if s[1] >= p - TOL]
    bdesc = ";".join(f"{x[0]}:{x[1][:6]}>{x[2][:6]}:{x[3]}" for x in bc if len(x[1]) != len(x[2])) or "-"
    rdesc = ";".join(f"{x[0]}:{x[3]}" for x in rc) or "-"
    print(f"{p}\t{ref[:14]}>{alt[:14]}\tdlen {dlen}\tdepth {depth}\tsupport {support} (mq30 {support_mq30})\t"
          f"catalog_sites {len(sc)}\tbam_solve: {bdesc}\tunion_rows: {rdesc}")
    summary["catalog_site" if sc else "no_catalog_site"] += 1
    summary["bam_called_indel" if any(len(x[1]) != len(x[2]) for x in bc) else "bam_no_indel"] += 1
print(dict(summary))
