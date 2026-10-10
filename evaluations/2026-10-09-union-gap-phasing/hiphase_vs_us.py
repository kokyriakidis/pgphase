#!/usr/bin/env python3
"""For phaseable arm reads HiPhase phases correctly but RUN leaves unphased:
which HiPhase-phased variants lie in the read, and what did RUN do with each?
(evaluation only)
Usage: hiphase_vs_us.py RUN_DIR HIPHASE_DIR QNAME_LIST"""
import bisect
import collections
import gzip
import sys

import pysam

sys.path.insert(0, ".")
import audit_hybrid as A

run, hip, qlist = sys.argv[1:4]
ROOT = "/home/kokyriakidis/Downloads/pgphase/test_data"
BAM = f"{ROOT}/HG002_chr20_hifi_mapped_to_CHM13_chr20_annotated.bam"
TRUTH_VCF = "/home/kokyriakidis/Downloads/pgphase-eval-data/results/chr12-18-20-comparison/chr20/truth.vcf.gz"
subset = set(open(qlist).read().split())
truth = A.load_truth(f"{ROOT}/derived/chr20_truth_hap.tsv")


def opener(p):
    return gzip.open(p, "rt") if p.endswith(".gz") else open(p)


def records(path, phased_only):
    out = {}
    for line in opener(path):
        if line[0] == "#":
            continue
        f = line.rstrip("\n").split("\t")
        fmt = dict(zip(f[8].split(":"), f[9].split(":")))
        gt = fmt.get("GT", "")
        a = gt.replace("/", "|").split("|")
        if len(a) != 2 or a[0] == a[1] or "." in a:
            continue
        phased = "|" in gt and fmt.get("PS", ".") not in (".", "")
        if phased_only and not phased:
            continue
        out[int(f[1])] = (f[3], f[4], phased)
    return out


hip_rec = records(f"{hip}/phased.vcf.gz", True)
hip_pos = sorted(hip_rec)
run_rec = records(f"{run}/phased.vcf", False)
truth_rec = records(TRUTH_VCF, False)
run_cand = {}
with open(f"{run}/candidates.tsv") as fh:
    next(fh)
    for line in fh:
        f = line.split("\t")
        run_cand.setdefault(int(f[1]), f[14])

hip_tags = A.load_tags(f"{hip}/phased.bam")
_, orient = A.score_reads(hip_tags, truth, set(truth))
run_tags = A.load_tags(f"{run}/phased.bam")

fate = collections.Counter()
per_read = collections.Counter()
kinds = collections.Counter()
n = 0
with pysam.AlignmentFile(BAM) as bam:
    for r in bam.fetch(until_eof=True):
        q = r.query_name
        if q not in subset or r.is_secondary or r.is_supplementary or q in run_tags:
            continue
        t = hip_tags.get(q)
        if t is None:
            continue
        if (((t[0] == 1) == (truth[q] == "MATERNAL")) == orient[t[1]]) is False:
            continue
        n += 1
        lo, hi = r.reference_start, r.reference_end
        inside = hip_pos[bisect.bisect_left(hip_pos, lo + 1):bisect.bisect_left(hip_pos, hi + 1)]
        per_read[min(len(inside), 5)] += 1
        for p in inside:
            ref, alt, _ = hip_rec[p]
            kind = "snp" if len(ref) == 1 and len(alt) == 1 else "indel"
            tol = 0 if kind == "snp" else 25
            near = lambda d: [x for x in range(p - tol, p + tol + 1) if x in d]
            in_truth = bool(near(truth_rec))
            rr = near(run_rec)
            rc = near(run_cand)
            if rr:
                s = "run_phased" if any(run_rec[x][2] for x in rr) else "run_unphased_het"
            elif rc:
                s = "run_cand:" + run_cand[rc[0]]
            else:
                s = "run_absent"
            fate[(kind, "truth" if in_truth else "not_truth", s)] += 1
            kinds[kind] += 1
print(f"reads HiPhase phases correctly, run leaves unphased: {n}")
print("HiPhase phased variants per read:", sorted(per_read.items()))
for k, v in fate.most_common(20):
    print(f"{v}\t{k}")
