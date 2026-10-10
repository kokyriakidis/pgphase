"""HiPhase-correct reads we leave unphased whose HiPhase variant matches a site
we phase: does the EM observe the read there? (evaluation only)"""
import bisect
import collections
import gzip
import sys

import pysam

sys.path.insert(0, ".")
import audit_hybrid as A

run, diag, qlist = sys.argv[1:4]
C = "/home/kokyriakidis/Downloads/pgphase-eval-data/results/chr12-18-20-comparison/chr20"
reads = set(open(qlist).read().split())
tags = A.load_tags(f"{run}/phased.bam")
reads = {q for q in reads if q not in tags}

hip = {}
for line in gzip.open(f"{C}/hiphase/phased.vcf.gz", "rt"):
    if line[0] == "#":
        continue
    f = line.rstrip("\n").split("\t")
    fmt = dict(zip(f[8].split(":"), f[9].split(":")))
    if "|" in fmt.get("GT", "") and fmt.get("PS", ".") != ".":
        hip[int(f[1])] = (f[3], f[4])
hpos = sorted(hip)
ours = {}
for line in open(f"{run}/phased.vcf"):
    if line[0] == "#":
        continue
    f = line.rstrip("\n").split("\t")
    fmt = dict(zip(f[8].split(":"), f[9].split(":")))
    if "|" in fmt.get("GT", "") and fmt.get("PS", ".") not in (".", ""):
        ours[int(f[1])] = f[7].split("CAT=")[-1].split(";")[0]
opos = sorted(ours)

dump = collections.defaultdict(list)
for line in open(diag):
    f = line.rstrip("\n").split("\t")
    if f[0] in reads:
        obs = {}
        for item in f[6:]:
            pos, cat, inj, err, oe, match, block = item.split(":")
            obs[int(pos)] = (float(err), cat, inj)
        dump[f[0]].append(obs)

tab = collections.Counter()
bam = pysam.AlignmentFile("../HG002_chr20_hifi_mapped_to_CHM13_chr20_annotated.bam")
for r in bam.fetch(until_eof=True):
    q = r.query_name
    if q not in reads or r.is_secondary or r.is_supplementary:
        continue
    lo, hi = r.reference_start, r.reference_end
    for p in hpos[bisect.bisect_left(hpos, lo + 1):bisect.bisect_left(hpos, hi + 1)]:
        ref, alt = hip[p]
        if len(ref) == len(alt):
            continue
        j = bisect.bisect_left(opos, p - 25)
        match = [x for x in opos[j:bisect.bisect_right(opos, p + 25)]]
        if not match:
            continue
        cat = ours[match[0]]
        if q not in dump:
            tab[(cat, "read_not_in_chunk")] += 1
            continue
        seen = [o[x] for o in dump[q] for x in match if x in o]
        if not seen:
            tab[(cat, "no_observation_at_site")] += 1
        elif all(s[0] >= 0.3 for s in seen):
            tab[(cat, "observed_site_unreliable")] += 1
        else:
            tab[(cat, "observed_reliable_but_low_posterior")] += 1
for k, v in tab.most_common():
    print(v, k)
