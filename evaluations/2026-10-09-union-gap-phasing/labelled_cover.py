"""For remaining HiPhase-correct unphased reads: at their HiPhase variants, how
many reads does RUN label per haplotype (dominant phase set)? (eval only)"""
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
hip = []
for line in gzip.open(f"{C}/hiphase/phased.vcf.gz", "rt"):
    if line[0] == "#":
        continue
    f = line.split("\t")
    if "|" in f[9].split(":")[0]:
        hip.append(int(f[1]) - 1)
hip.sort()
bam = pysam.AlignmentFile("../HG002_chr20_hifi_mapped_to_CHM13_chr20_annotated.bam")
tab = collections.Counter()
done = 0
for r in bam.fetch(until_eof=True):
    q = r.query_name
    if q not in reads or r.is_secondary or r.is_supplementary:
        continue
    vs = hip[bisect.bisect_left(hip, r.reference_start):bisect.bisect_left(hip, r.reference_end)]
    if not vs:
        continue
    p = vs[len(vs) // 2]
    cov = collections.Counter()
    unl = 0
    for o in bam.fetch("CHM13#0#chr20", p - 10, p + 10):
        if o.is_secondary or o.is_supplementary or o.reference_start > p - 10 or o.reference_end < p + 10:
            continue
        t = tags.get(o.query_name)
        if t is None:
            unl += 1
        else:
            cov[(t[1], t[0])] += 1
    ps_tot = collections.Counter()
    for (ps, hp), n in cov.items():
        ps_tot[ps] += n
    if ps_tot:
        ps = ps_tot.most_common(1)[0][0]
        h1, h2 = cov[(ps, 1)], cov[(ps, 2)]
    else:
        h1 = h2 = 0
    lo = min(h1, h2)
    tab["min_side " + ("0" if lo == 0 else "1-2" if lo < 3 else "3-5" if lo < 6 else ">=6")] += 1
    tab["unlabelled_cover " + ("<10" if unl < 10 else "10-20" if unl < 20 else ">=20")] += 1
    done += 1
print("reads:", done)
for k, v in sorted(tab.items()):
    print(v, k)
