"""For HiPhase-correct reads RUN still leaves unphased: the HiPhase variants in
each read and what the alignment-only solve called near them. (eval only)
Usage: remaining_variants.py RUN_DIR QLIST"""
import bisect
import collections
import gzip
import sys

import pysam

sys.path.insert(0, ".")
import audit_hybrid as A

run, qlist = sys.argv[1:3]
C = "/home/kokyriakidis/Downloads/pgphase-eval-data/results/chr12-18-20-comparison/chr20"
tags = A.load_tags(f"{run}/phased.bam")
reads = {q for q in open(qlist).read().split() if q not in tags}

hip = {}
for line in gzip.open(f"{C}/hiphase/phased.vcf.gz", "rt"):
    if line[0] == "#":
        continue
    f = line.rstrip("\n").split("\t")
    fmt = dict(zip(f[8].split(":"), f[9].split(":")))
    if "|" in fmt.get("GT", "") and fmt.get("PS", ".") != ".":
        hip[int(f[1])] = (f[3], f[4], fmt["GT"])
hpos = sorted(hip)
bc = []
with open("bam/candidates.tsv") as fh:
    next(fh)
    for line in fh:
        f = line.rstrip("\n").split("\t")
        bc.append((int(f[1]), f[3], f[4], f[14], int(f[6]), int(f[7])))
bc.sort()
bpos = [b[0] for b in bc]

tab = collections.Counter()
per_read = collections.Counter()
ex = collections.defaultdict(list)
bam = pysam.AlignmentFile("../HG002_chr20_hifi_mapped_to_CHM13_chr20_annotated.bam")
for r in bam.fetch(until_eof=True):
    q = r.query_name
    if q not in reads or r.is_secondary or r.is_supplementary:
        continue
    vs = hpos[bisect.bisect_left(hpos, r.reference_start + 1):bisect.bisect_left(hpos, r.reference_end + 1)]
    per_read[min(len(vs), 4)] += 1
    for p in vs:
        ref, alt, gt = hip[p]
        alts = alt.split(",")
        kind = "snp" if len(ref) == 1 and all(len(a) == 1 for a in alts) else (
            "indel_1bp" if all(abs(len(a) - len(ref)) == 1 for a in alts) else "indel_big")
        multi = "1/2" if len(alts) > 1 else "0/1"
        tol = 0 if kind == "snp" else 15
        near = bc[bisect.bisect_left(bpos, p - tol):bisect.bisect_right(bpos, p + tol + 1)]
        if kind != "snp":
            near = [b for b in near if len(b[1]) != len(b[2]) or b[2] == "."]
        cat = ",".join(sorted(set(b[3] for b in near))) or "not_called"
        tab[(kind, multi, cat)] += 1
        if len(ex[(kind, multi, cat)]) < 2:
            ex[(kind, multi, cat)].append((p, ref[:10], alt[:16], [(b[0], b[1][:6], b[2][:6], b[4], b[5]) for b in near][:2]))
print("reads:", len(reads), "variants per read:", sorted(per_read.items()))
for k, v in tab.most_common(16):
    print(v, k, ex[k][:1])
