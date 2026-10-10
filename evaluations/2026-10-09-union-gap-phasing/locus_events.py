"""At remaining reads' HiPhase variants: events of labelled reads by haplotype."""
import bisect
import collections
import gzip
import sys

import pysam

sys.path.insert(0, ".")
import audit_hybrid as A

run, qlist, nshow = sys.argv[1], sys.argv[2], int(sys.argv[3])
C = "/home/kokyriakidis/Downloads/pgphase-eval-data/results/chr12-18-20-comparison/chr20"
tags = A.load_tags(f"{run}/phased.bam")
rem = [q for q in open(qlist).read().split() if q not in tags]
hip = {}
for line in gzip.open(f"{C}/hiphase/phased.vcf.gz", "rt"):
    if line[0] == "#":
        continue
    f = line.split("\t")
    if "|" in f[9].split(":")[0]:
        hip[int(f[1]) - 1] = (f[3], f[4], f[9].split(":")[0])
hpos = sorted(hip)
bam = pysam.AlignmentFile("../HG002_chr20_hifi_mapped_to_CHM13_chr20_annotated.bam")
first = {}
remset = set(rem)
for r in bam.fetch(until_eof=True):
    if r.query_name in remset and not r.is_secondary and not r.is_supplementary and r.query_name not in first:
        first[r.query_name] = (r.reference_start, r.reference_end)
shown = 0
for q in rem:
    if q not in first:
        continue
    lo, hi = first[q]
    vs = hpos[bisect.bisect_left(hpos, lo):bisect.bisect_left(hpos, hi)]
    if not vs:
        continue
    p = vs[0]
    ref, alt, gt = hip[p]
    print(f"== read {q} span {lo}-{hi}; HiPhase variant {p+1} {ref}>{alt} {gt}")
    by = collections.defaultdict(collections.Counter)
    for o in bam.fetch("CHM13#0#chr20", p - 5, p + 5):
        if o.is_secondary or o.is_supplementary or o.reference_start > p - 20 or o.reference_end < p + len(ref) + 20:
            continue
        t = tags.get(o.query_name)
        key = f"PS{t[1]}:HP{t[0]}" if t else ("TARGET" if o.query_name == q else "unlabelled")
        rp = o.reference_start
        ev = []
        for op, n in o.cigartuples:
            if op in (0, 7, 8):
                rp += n
            elif op == 2:
                if p - 15 <= rp <= p + len(ref) + 15:
                    ev.append(f"D{n}@{rp - p:+d}")
                rp += n
            elif op == 1:
                if p - 15 <= rp <= p + len(ref) + 15:
                    ev.append(f"I{n}@{rp - p:+d}")
            if rp > p + 40:
                break
        by[key][" ".join(ev) or "-"] += 1
    for k in sorted(by):
        print("   ", k, dict(by[k].most_common(4)))
    shown += 1
    if shown >= nshow:
        break
