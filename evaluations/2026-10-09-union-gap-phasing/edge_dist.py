"""For reads spanning no site we phased: distance from the HiPhase variant to the read end."""
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
ours = sorted(int(l.split("\t")[1]) for l in open(f"{run}/phased.vcf") if l[0] != "#" and "|" in l.split("\t")[9].split(":")[0])
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
for r in bam.fetch(until_eof=True):
    q = r.query_name
    if q not in reads or r.is_secondary or r.is_supplementary:
        continue
    lo, hi = r.reference_start, r.reference_end
    if bisect.bisect_left(ours, hi + 1) > bisect.bisect_left(ours, lo + 1):
        continue  # spans a site we phased
    vs = hip[bisect.bisect_left(hip, lo):bisect.bisect_left(hip, hi)]
    if not vs:
        tab["no HiPhase variant in alignment span"] += 1
        continue
    d = max(min(p - lo, hi - p) for p in vs)  # best-placed variant
    tab["edge <=100bp" if d <= 100 else "edge 100-500bp" if d <= 500 else "interior >500bp"] += 1
for k, v in tab.most_common():
    print(v, k)
