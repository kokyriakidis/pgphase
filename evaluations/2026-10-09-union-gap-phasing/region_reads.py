#!/usr/bin/env python3
"""Read-level scores split into chromosome arms and the centromere (26-32 Mb),
and by whether the read crosses a truth heterozygote.

Usage: region_reads.py INPUT_BAM TRUTH_TSV TRUTH_VCF NAME=DIR [NAME=DIR ...]
Same rules as audit_hybrid.py: truth-labelled primary reads are the
denominator; each phase set takes its majority orientation over all of its
reads. A read belongs to the centromere when its alignment midpoint is in
[26, 32) Mb. A read is "phaseable" when its aligned span contains at least one
heterozygous record of TRUTH_VCF; otherwise "nohet" (the two haplotypes are
identical over it, so any label it gets is noise). Rows per (name, region)
with region in arms, cen, arms_phaseable, arms_nohet, cen_phaseable,
cen_nohet. Truth is evaluation-only.
"""
import bisect
import gzip
import multiprocessing as mp
import os
import sys

import pysam

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import audit_hybrid as A  # noqa: E402

CEN = (26_000_000, 32_000_000)
REGIONS = ("arms", "cen", "arms_phaseable", "arms_nohet", "cen_phaseable", "cen_nohet")
G = {}


def truth_hets(path):
    out = []
    for line in gzip.open(path, "rt"):
        if line[0] == "#":
            continue
        f = line.split("\t")
        gt = f[9].split(":")[0].replace("/", "|").split("|")
        if len(gt) == 2 and gt[0] != gt[1] and "." not in gt:
            out.append(int(f[1]) - 1)
    out.sort()
    return out


def classify(bam_path, truth, hets):
    region = {}
    with pysam.AlignmentFile(bam_path) as bam:
        for r in bam.fetch(until_eof=True):
            if r.is_unmapped or r.is_secondary or r.is_supplementary or r.is_duplicate:
                continue
            if r.query_name in truth:
                lo, hi = r.reference_start, r.reference_end
                arm = "cen" if CEN[0] <= (lo + hi) // 2 < CEN[1] else "arms"
                crosses = bisect.bisect_left(hets, hi) > bisect.bisect_left(hets, lo)
                region[r.query_name] = (arm, f"{arm}_{'phaseable' if crosses else 'nohet'}")
    return region


def score(item):
    name, d = item
    truth, region = G["truth"], G["region"]
    tags = A.load_tags(f"{d}/phased.bam")
    _, orient = A.score_reads(tags, truth, set(region))
    out = {r: {"correct": 0, "discordant": 0, "unphased": 0} for r in REGIONS}
    for q, regs in region.items():
        t = tags.get(q)
        if t is None:
            cls = "unphased"
        else:
            ok = ((t[0] == 1) == (truth[q] == "MATERNAL")) == orient[t[1]]
            cls = "correct" if ok else "discordant"
        for reg in regs:
            out[reg][cls] += 1
    return name, out


def main():
    in_bam, truth_path, truth_vcf = sys.argv[1:4]
    arms = [a.split("=", 1) for a in sys.argv[4:]]
    G["truth"] = A.load_truth(truth_path)
    G["region"] = classify(in_bam, G["truth"], truth_hets(truth_vcf))
    with mp.Pool(min(len(arms), 8)) as pool:
        results = pool.map(score, arms, chunksize=1)
    print("name\tregion\tcorrect\tdiscordant\tunphased\tdenominator")
    for name, out in results:
        for reg in REGIONS:
            o = out[reg]
            print(f"{name}\t{reg}\t{o['correct']}\t{o['discordant']}\t{o['unphased']}\t{sum(o.values())}")


if __name__ == "__main__":
    main()
