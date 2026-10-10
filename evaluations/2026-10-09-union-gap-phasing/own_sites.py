"""HiPhase-correct reads RUN leaves unphased: our own phased sites in each read,
and the EM's observations of the read. (evaluation only)"""
import bisect
import collections
import sys

import pysam

sys.path.insert(0, ".")
import audit_hybrid as A

run, dump_path, qlist = sys.argv[1:4]
tags = A.load_tags(f"{run}/phased.bam")
reads = {q for q in open(qlist).read().split() if q not in tags}
ours = []
for line in open(f"{run}/phased.vcf"):
    if line[0] == "#":
        continue
    f = line.rstrip("\n").split("\t")
    fmt = dict(zip(f[8].split(":"), f[9].split(":")))
    if "|" in fmt.get("GT", "") and fmt.get("PS", ".") not in (".", "") and len(set(fmt["GT"].split("|"))) > 1:
        ours.append((int(f[1]), len(f[3]) == 1 and len(f[4]) == 1, int(fmt["PS"])))
ours.sort()
opos = [o[0] for o in ours]
dump = collections.defaultdict(list)
for line in open(dump_path):
    f = line.rstrip("\n").split("\t")
    if f[0] in reads:
        dump[f[0]].append((int(f[1]), int(f[3]), [x.split(":") for x in f[6:]], int(f[4]), int(f[5])))
bam = pysam.AlignmentFile("../HG002_chr20_hifi_mapped_to_CHM13_chr20_annotated.bam")
tab = collections.Counter()
examples = collections.defaultdict(list)
for r in bam.fetch(until_eof=True):
    q = r.query_name
    if q not in reads or r.is_secondary or r.is_supplementary:
        continue
    lo, hi = r.reference_start, r.reference_end
    inside = ours[bisect.bisect_left(opos, lo + 1):bisect.bisect_left(opos, hi + 1)]
    n_snp = sum(1 for o in inside if o[1])
    d = dump.get(q)
    if not inside:
        cls = "spans no site we phased"
    elif d is None:
        cls = "spans our phased sites; read in no chunk"
    else:
        obs_pos = set()
        for mapq, nobs, obs, cb, ce in d:
            for o in obs:
                obs_pos.add(int(o[0]))
        seen = sum(1 for o in inside if any(abs(o[0] - p) <= 1 for p in obs_pos))
        if seen == 0:
            cls = "spans our phased sites; EM has no observation at any"
        elif seen < len(inside):
            cls = "spans our phased sites; EM observes some"
        else:
            cls = "spans our phased sites; EM observes all"
    key = (cls, "has_snp_site" if n_snp else "indel_sites_only")
    tab[key] += 1
    if len(examples[key]) < 3:
        examples[key].append((q, lo, hi, r.mapping_quality, len(inside), n_snp,
                              [(m, n, cb, ce) for m, n, _, cb, ce in (d or [])][:2]))
print("reads:", sum(tab.values()))
for k, v in tab.most_common():
    print(v, k)
    for e in examples[k][:2]:
        print("     ", e)
