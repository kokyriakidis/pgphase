import bisect
import collections
import sys

import pysam

sys.path.insert(0, ".")
import audit_hybrid as A

targets = set(open("full/clean5/phaseable_unphased.txt").read().split())
miss = sorted(int(l.split("\t")[0]) for l in open("full/af9/hiphase_missing_indels.tsv"))
loci = []
for l in open("full/diag2/loci.tsv"):
    a, b, why = l.split("\t")[:3]
    loci.append((int(a), int(b), why))
loci.sort()
ls = [x[0] for x in loci]
dump = {}
for line in open("full/diag2/reads.tsv"):
    f = line.split("\t", 6)
    dump[f[0]] = (f[2], f[3])
tags9 = A.load_tags("full/af9/phased.bam")
bam = pysam.AlignmentFile("../HG002_chr20_hifi_mapped_to_CHM13_chr20_annotated.bam")
tab = collections.Counter()
for p in miss:
    i = bisect.bisect_right(ls, p) - 1
    win = None
    for j in range(max(0, i - 3), min(len(loci), i + 3)):
        if loci[j][0] - 50 <= p <= loci[j][1] + 50:
            win = loci[j]
    for r in bam.fetch("CHM13#0#chr20", p - 1, p):
        q = r.query_name
        if q not in targets or r.is_secondary or r.is_supplementary or q in tags9:
            continue
        spans = win is not None and r.reference_start <= win[0] and r.reference_end >= win[1]
        d = dump.get(q)
        mq = "mapq>=20" if r.mapping_quality >= 20 else "mapq<20"
        nobs = int(d[1]) if d else -1
        tab[("window " + (win[2] if win else "none"), "spans" if spans else "partial",
             "in_chunk" if d else "absent", mq, "obs>0" if nobs > 0 else "no_obs")] += 1
for k, v in tab.most_common():
    print(v, k)
