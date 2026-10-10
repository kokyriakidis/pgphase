#!/usr/bin/env python3
"""Arm reads spanning no truth heterozygote: where are they, and are the ones
other runs phase phased correctly? (evaluation only)
Usage: nohet_reads.py INV_DIR TRUTH_VCF NAME=RUN_DIR [NAME=RUN_DIR ...]"""
import bisect
import collections
import gzip
import sys

import pysam

inv, truth_vcf = sys.argv[1:3]
runs = [a.split("=", 1) for a in sys.argv[3:]]
sys.path.insert(0, inv)
import audit_hybrid as A  # noqa: E402

ROOT = "/home/kokyriakidis/Downloads/pgphase/test_data"
BAM = f"{ROOT}/HG002_chr20_hifi_mapped_to_CHM13_chr20_annotated.bam"
CEN = (26_000_000, 32_000_000)

hets = []
for line in gzip.open(truth_vcf, "rt"):
    if line[0] == "#":
        continue
    f = line.split("\t")
    gt = f[9].split(":")[0].replace("/", "|").split("|")
    if len(gt) == 2 and gt[0] != gt[1]:
        hets.append(int(f[1]) - 1)
hets.sort()
truth = A.load_truth(f"{ROOT}/derived/chr20_truth_hap.tsv")

nohet, spans = set(), []
with pysam.AlignmentFile(BAM) as bam:
    for r in bam.fetch(until_eof=True):
        if r.is_unmapped or r.is_secondary or r.is_supplementary or r.is_duplicate:
            continue
        if r.query_name not in truth:
            continue
        lo, hi = r.reference_start, r.reference_end
        if CEN[0] <= (lo + hi) // 2 < CEN[1]:
            continue
        if bisect.bisect_left(hets, hi) == bisect.bisect_left(hets, lo):
            nohet.add(r.query_name)
            spans.append((lo, hi))
spans.sort()
merged = []
for lo, hi in spans:
    if merged and lo <= merged[-1][1]:
        merged[-1][1] = max(merged[-1][1], hi)
    else:
        merged.append([lo, hi])
lengths = sorted((b - a for a, b in merged), reverse=True)
print(f"arm reads with no truth het: {len(nohet)}")
print(f"merged intervals: {len(merged)}; total {sum(lengths):,} bp; "
      f"intervals >=100 kb: {sum(1 for x in lengths if x >= 100_000)} covering "
      f"{sum(x for x in lengths if x >= 100_000):,} bp")
big = [[a, b] for a, b in merged if b - a >= 100_000]
for a, b in sorted(big)[:15]:
    n = sum(1 for lo, hi in spans if a <= lo and hi <= b)
    print(f"  {a:,}-{b:,} ({(b - a) / 1e6:.2f} Mb) reads {n}")

for name, d in runs:
    tags = A.load_tags(f"{d}/phased.bam")
    # Orientation per phase set from all of that phase set's reads.
    _, orient = A.score_reads(tags, truth, set(truth))
    c = collections.Counter()
    for q in nohet:
        t = tags.get(q)
        if t is None:
            continue
        ok = ((t[0] == 1) == (truth[q] == "MATERNAL")) == orient[t[1]]
        c["correct" if ok else "discordant"] += 1
    print(f"{name}: phases {sum(c.values())} no-het reads: correct {c['correct']} discordant {c['discordant']}")
