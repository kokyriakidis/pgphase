#!/usr/bin/env python3
"""Full accounting: every read HiPhase phases correctly that RUN leaves
unphased, classified by why. (evaluation only)

Usage: why_missing.py RUN_DIR DIAG_READS QLIST
Per read, each HiPhase-phased variant in its span is located in our pipeline:
  A  the read observes a reliable EM site at that variant (label withheld:
     posterior below 0.8 or conflicting sites)
  B  the read observes an unreliable EM site there, and that site's calls do
     separate the parents (>= 90%): EM error drift
  C  the read observes an unreliable EM site there whose calls do not
     separate the parents: our representation / calls are wrong
  D  an EM site exists there, but this read has no call at it
  E  no EM site there: why (graph filter reason / alignment category / none)
A read takes its best variant (A..E order). Also reported: whether HiPhase's
variant is a truth het at all, and the read's mapping quality.
"""
import bisect
import collections
import gzip
import sys

import pysam

sys.path.insert(0, ".")
import audit_hybrid as A

run, dump_path, qlist = sys.argv[1:4]
ROOT = "/home/kokyriakidis/Downloads/pgphase/test_data"
C = "/home/kokyriakidis/Downloads/pgphase-eval-data/results/chr12-18-20-comparison/chr20"
TOL = 25

truth_lab = A.load_truth(f"{ROOT}/derived/chr20_truth_hap.tsv")
tags = A.load_tags(f"{run}/phased.bam")
reads = {q for q in open(qlist).read().split() if q not in tags}


def vcf_records(path, phased_only=False):
    out = []
    op = gzip.open if path.endswith(".gz") else open
    for line in op(path, "rt"):
        if line[0] == "#":
            continue
        f = line.rstrip("\n").split("\t")
        gt = f[9].split(":")[0]
        al = gt.replace("/", "|").split("|")
        if len(al) != 2 or al[0] == al[1] or "." in al:
            continue
        if phased_only and "|" not in gt:
            continue
        out.append((int(f[1]), f[3], f[4]))
    out.sort()
    return out


hip = vcf_records(f"{C}/hiphase/phased.vcf.gz", phased_only=True)
hpos = [h[0] for h in hip]
truth_het = vcf_records(f"{C}/truth.vcf.gz")
tpos = [t[0] for t in truth_het]


def near(sorted_pos, p, tol):
    return bisect.bisect_left(sorted_pos, p - tol) < bisect.bisect_right(sorted_pos, p + tol)


# EM view: per read, observations (pos, err, injected/locus); per site, calls by parent.
read_obs = collections.defaultdict(list)
site_calls = collections.defaultdict(collections.Counter)
site_positions = set()
for line in open(dump_path):
    f = line.rstrip("\n").split("\t")
    q = f[0]
    for o in (x.split(":") for x in f[6:]):
        p, err = int(o[0]), float(o[3])
        site_positions.add(p)
        kind = "window" if o[1] == "-1" else ("injected" if o[2] == "1" else "graph")
        if q in reads:
            read_obs[q].append((p, err, kind))
        if q in truth_lab:
            site_calls[p][(o[5], truth_lab[q])] += 1
spos = sorted(site_positions)


def purity(p):
    c = site_calls[p]
    cis = c[("1", "MATERNAL")] + c[("2", "PATERNAL")]
    tr = c[("2", "MATERNAL")] + c[("1", "PATERNAL")]
    return max(cis, tr) / (cis + tr) if cis + tr else 0.0


filtered = collections.defaultdict(set)
with open(f"{run}/filtered.tsv") as fh:
    next(fh)
    for line in fh:
        f = line.split("\t")
        filtered[int(f[1])].add(f[7].strip())
fpos = sorted(filtered)
bam_cat = collections.defaultdict(set)
with open("bam/candidates.tsv") as fh:
    next(fh)
    for line in fh:
        f = line.split("\t")
        bam_cat[int(f[1])].add(f[14])
bpos = sorted(bam_cat)


def nearby(d, sorted_keys, p, tol):
    out = set()
    for k in sorted_keys[bisect.bisect_left(sorted_keys, p - tol):bisect.bisect_right(sorted_keys, p + tol)]:
        out |= d[k]
    return out


ORDER = "ABCDE"
tab = collections.Counter()
sub = collections.Counter()
mapq_tab = collections.Counter()
examples = collections.defaultdict(list)
bam = pysam.AlignmentFile(f"{ROOT}/HG002_chr20_hifi_mapped_to_CHM13_chr20_annotated.bam")
for r in bam.fetch(until_eof=True):
    q = r.query_name
    if q not in reads or r.is_secondary or r.is_supplementary:
        continue
    reads.discard(q)
    lo, hi = r.reference_start + 1, r.reference_end
    vs = hip[bisect.bisect_left(hpos, lo):bisect.bisect_right(hpos, hi)]
    obs = read_obs.get(q, [])
    best = None
    for p, ref, alt in vs:
        tol = 0 if len(ref) == 1 and all(len(a) == 1 for a in alt.split(",")) else TOL
        is_truth = near(tpos, p, max(tol, 1) if tol else 0) if tol else (p in set(tpos[bisect.bisect_left(tpos, p):bisect.bisect_right(tpos, p)]))
        mine = [o for o in obs if abs(o[0] - p) <= max(tol, 1)]
        if any(o[1] < 0.3 for o in mine):
            cls, detail = "A", "reliable site observed; label withheld"
        elif mine:
            pur = max(purity(o[0]) for o in mine)
            kinds = "+".join(sorted(set(o[2] for o in mine)))
            cls, detail = ("B", f"unreliable {kinds} site, calls separate parents") if pur >= 0.9 else \
                          ("C", f"unreliable {kinds} site, calls do not separate parents")
        elif near(spos, p, max(tol, 1)):
            cls, detail = "D", "EM site exists, read has no call"
        else:
            g = nearby(filtered, fpos, p, max(tol, 1)) - {"ref_only"}
            b = nearby(bam_cat, bpos, p, max(tol, 1))
            why = ("graph filtered:" + ",".join(sorted(g))) if g else ("alignment:" + ",".join(sorted(b))) if b else "called by neither"
            cls, detail = "E", "no EM site; " + why
        detail += " | HiPhase variant " + ("is truth het" if is_truth else "NOT truth het")
        if best is None or ORDER.index(cls) < ORDER.index(best[0]):
            best = (cls, detail, p, ref[:10], alt[:16])
    if best is None:
        best = ("F", "no HiPhase phased variant in the read span", 0, "", "")
    tab[best[0]] += 1
    sub[(best[0], best[1])] += 1
    mapq_tab[(best[0], "mapq>=20" if r.mapping_quality >= 20 else "mapq<20")] += 1
    if len(examples[(best[0], best[1])]) < 2:
        examples[(best[0], best[1])].append((q.split("/")[1], best[2], best[3], best[4]))
names = {"A": "observed reliable site, label withheld", "B": "EM error drift (good calls marked unreliable)",
         "C": "our calls do not separate the parents", "D": "EM site exists, read has no call",
         "E": "no EM site at HiPhase's variant", "F": "no HiPhase variant in span"}
print("reads:", sum(tab.values()))
for k in ORDER + "F":
    if tab[k]:
        print(f"\n{k} {names[k]}: {tab[k]}  ({mapq_tab[(k, 'mapq>=20')]} MAPQ>=20)")
        for (c, d), v in sorted(sub.items(), key=lambda kv: -kv[1]):
            if c == k:
                print(f"    {v:4d}  {d}   e.g. {examples[(c, d)][:1]}")
