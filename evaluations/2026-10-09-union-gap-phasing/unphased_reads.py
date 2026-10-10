#!/usr/bin/env python3
"""Why are arm reads unphased? (evaluation only)

Usage: unphased_reads.py INV_DIR RUN_DIR TRUTH_VCF SD_BED [OTHER=DIR ...]
For every truth-labelled primary read on the arms (midpoint outside 26-32 Mb)
that RUN_DIR leaves without HP/PS, report: segmental-duplication overlap,
truth heterozygotes spanned, MAPQ, and the run's own phased / emitted hets
spanned. OTHER runs (competitors) say whether they phase the same reads.
"""
import bisect
import collections
import gzip
import sys

import pysam

inv, run, truth_vcf, sd_bed = sys.argv[1:5]
others = [a.split("=", 1) for a in sys.argv[5:]]
sys.path.insert(0, inv)
import audit_hybrid as A  # noqa: E402

ROOT = "/home/kokyriakidis/Downloads/pgphase/test_data"
BAM = f"{ROOT}/HG002_chr20_hifi_mapped_to_CHM13_chr20_annotated.bam"
TRUTH_TSV = f"{ROOT}/derived/chr20_truth_hap.tsv"
CEN = (26_000_000, 32_000_000)


def het_positions(path, phased_only):
    out = []
    op = gzip.open if path.endswith(".gz") else open
    for line in op(path, "rt"):
        if line[0] == "#":
            continue
        f = line.rstrip("\n").split("\t")
        fmt = dict(zip(f[8].split(":"), f[9].split(":")))
        gt = fmt.get("GT", "").replace("/", "|")
        a = gt.split("|")
        if len(a) != 2 or a[0] == a[1] or "." in a:
            continue
        if phased_only and ("|" not in fmt.get("GT", "") or fmt.get("PS", ".") in (".", "")):
            continue
        out.append(int(f[1]) - 1)
    out.sort()
    return out


def count_in(sorted_pos, lo, hi):
    return bisect.bisect_left(sorted_pos, hi) - bisect.bisect_left(sorted_pos, lo)


sd = []
for line in open(sd_bed):
    f = line.split("\t")
    if f[0] == "chr20":
        sd.append((int(f[1]), int(f[2])))
sd.sort()
merged = []
for a, b in sd:
    if merged and a <= merged[-1][1]:
        merged[-1][1] = max(merged[-1][1], b)
    else:
        merged.append([a, b])
sd_starts = [m[0] for m in merged]


def sd_fraction(lo, hi):
    i = max(0, bisect.bisect_right(sd_starts, lo) - 1)
    cov = 0
    while i < len(merged) and merged[i][0] < hi:
        cov += max(0, min(hi, merged[i][1]) - max(lo, merged[i][0]))
        i += 1
    return cov / max(1, hi - lo)


truth = A.load_truth(TRUTH_TSV)
truth_hets = het_positions(truth_vcf, False)
run_phased = het_positions(f"{run}/phased.vcf", True)
run_any = het_positions(f"{run}/phased.vcf", False)
run_tags = A.load_tags(f"{run}/phased.bam")
other_tags = {name: A.load_tags(f"{d}/phased.bam") for name, d in others}

rows = collections.Counter()
by_other = collections.Counter()
examples = []
with pysam.AlignmentFile(BAM) as bam:
    for r in bam.fetch(until_eof=True):
        if r.is_unmapped or r.is_secondary or r.is_supplementary or r.is_duplicate:
            continue
        q = r.query_name
        if q not in truth:
            continue
        lo, hi = r.reference_start, r.reference_end
        if CEN[0] <= (lo + hi) // 2 < CEN[1]:
            continue
        rows["arm_reads"] += 1
        if q in run_tags:
            continue
        rows["unphased"] += 1
        sdf = sd_fraction(lo, hi)
        nt = count_in(truth_hets, lo, hi)
        if sdf >= 0.5:
            cls = "A_segdup>=50%"
        elif nt == 0:
            cls = "B_no_truth_het_spanned"
        elif r.mapping_quality < 20:
            cls = "C_mapq<20"
        elif count_in(run_phased, lo, hi) > 0:
            cls = "D1_run_has_phased_hets_in_span"
        elif count_in(run_any, lo, hi) > 0:
            cls = "D2_run_has_only_unphased_hets"
        else:
            cls = "D3_run_has_no_het_in_span"
        rows[cls] += 1
        if cls.startswith("D"):
            rows["D_truth_hets_" + ("1" if nt == 1 else "2-4" if nt <= 4 else ">=5")] += 1
            if len(examples) < 2000:
                examples.append((q, lo, hi, r.mapping_quality, nt, count_in(run_phased, lo, hi), cls))
        for name, tags in other_tags.items():
            if q in tags:
                by_other[(cls[0], name)] += 1

for k in sorted(rows):
    print(f"{k}\t{rows[k]}")
for k in sorted(by_other):
    print(f"phased_by\t{k[0]}\t{k[1]}\t{by_other[k]}")
with open(f"{run}/unphased_arm_examples.tsv", "w") as fh:
    fh.write("qname\tbeg\tend\tmapq\ttruth_hets\trun_phased_hets\tclass\n")
    for e in examples:
        fh.write("\t".join(map(str, e)) + "\n")
