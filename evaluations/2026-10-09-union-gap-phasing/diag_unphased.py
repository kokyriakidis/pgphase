#!/usr/bin/env python3
"""Why are phaseable arm reads unphased? (evaluation only)

Usage: diag_unphased.py INV_DIR RUN_DIR DIAG_READS TRUTH_VCF SD_BED SITES_VCF
A phaseable read: truth-labelled primary arm read, MAPQ>=20, <50% in segdups,
spanning >=1 truth heterozygote. For those RUN_DIR leaves unphased, report the
EM's view of the read (DIAG_READS, written by a PGPHASE_DIAG_READS build) and,
for every truth het in its span, what the catalog, the alignment solve and the
run did with that het.
"""
import bisect
import collections
import gzip
import math
import sys

import pysam

inv, run, diag_path, truth_vcf, sd_bed, sites_vcf = sys.argv[1:7]
sys.path.insert(0, inv)
import audit_hybrid as A  # noqa: E402

ROOT = "/home/kokyriakidis/Downloads/pgphase/test_data"
BAM = f"{ROOT}/HG002_chr20_hifi_mapped_to_CHM13_chr20_annotated.bam"
CEN = (26_000_000, 32_000_000)
CAT_BITS = {0x004: "cleanSNP", 0x008: "cleanIndel", 0x100: "noisyHet"}


def opener(p):
    return gzip.open(p, "rt") if p.endswith(".gz") else open(p)


# Truth hets: 0-based pos -> is_snp
truth = {}
for line in opener(truth_vcf):
    if line[0] == "#":
        continue
    f = line.split("\t")
    gt = f[9].split(":")[0].replace("/", "|").split("|")
    if len(gt) == 2 and gt[0] != gt[1] and "." not in gt:
        truth[int(f[1]) - 1] = len(f[3]) == 1 and all(len(a) == 1 for a in f[4].split(","))
tpos = sorted(truth)


def table(path):
    """candidates.tsv: 0-based pos -> list of (category, phase_set)."""
    out = collections.defaultdict(list)
    with open(path) as fh:
        next(fh)
        for line in fh:
            f = line.split("\t")
            out[int(f[1]) - 1].append((f[14], f[16]))
    return out


def vcf_records(path):
    out = collections.defaultdict(list)
    for line in opener(path):
        if line[0] == "#":
            continue
        f = line.rstrip("\n").split("\t")
        fmt = dict(zip(f[8].split(":"), f[9].split(":")))
        gt = fmt.get("GT", "")
        alle = gt.replace("/", "|").split("|")
        state = "hom" if len(set(alle)) == 1 else ("phased" if "|" in gt and fmt.get("PS", ".") not in (".", "") else "unphased_het")
        out[int(f[1]) - 1].append(state)
    return out


def lookup(d, p, snp):
    if snp:
        return d.get(p, [])
    hits = []
    for q in range(p - 10, p + 11):
        hits += d.get(q, [])
    return hits


run_cand = table(f"{run}/candidates.tsv")
bam_cand = table(f"{inv}/bam/candidates.tsv")
run_vcf = vcf_records(f"{run}/phased.vcf")
sites = []
for line in opener(sites_vcf):
    if line[0] == "#":
        continue
    f = line.split("\t", 5)
    p = int(f[1]) - 1
    sites.append((p, p + len(f[3])))
sites.sort()
site_starts = [s[0] for s in sites]


def in_catalog(p):
    i = bisect.bisect_right(site_starts, p + 10)
    for j in range(max(0, i - 50), i):
        if sites[j][0] - 10 <= p < sites[j][1] + 10:
            return True
    return False


sd = sorted((int(f[1]), int(f[2])) for f in (l.split("\t") for l in open(sd_bed)) if f[0] == "chr20")
merged = []
for a, b in sd:
    if merged and a <= merged[-1][1]:
        merged[-1][1] = max(merged[-1][1], b)
    else:
        merged.append([a, b])
sd_starts = [m[0] for m in merged]


def sd_fraction(lo, hi):
    i = max(0, bisect.bisect_right(sd_starts, lo) - 1)
    cov = 0
    while i < len(merged) and merged[i][0] < hi:
        cov += max(0, min(hi, merged[i][1]) - max(lo, merged[i][0]))
        i += 1
    return cov / max(1, hi - lo)


labels = A.load_truth(f"{ROOT}/derived/chr20_truth_hap.tsv")
tags = A.load_tags(f"{run}/phased.bam")
targets = {}
with pysam.AlignmentFile(BAM) as bam:
    for r in bam.fetch(until_eof=True):
        if r.is_unmapped or r.is_secondary or r.is_supplementary or r.is_duplicate:
            continue
        q = r.query_name
        lo, hi = r.reference_start, r.reference_end
        if q not in labels or q in tags or CEN[0] <= (lo + hi) // 2 < CEN[1]:
            continue
        if r.mapping_quality < 20 or sd_fraction(lo, hi) >= 0.5:
            continue
        hets = tpos[bisect.bisect_left(tpos, lo):bisect.bisect_left(tpos, hi)]
        if hets:
            targets[q] = (lo, hi, hets)
print(f"phaseable unphased arm reads: {len(targets)}")
with open(f"{run}/phaseable_unphased.txt", "w") as fh:
    fh.write("\n".join(targets) + "\n")

dumps = collections.defaultdict(list)
for line in open(diag_path):
    f = line.rstrip("\n").split("\t")
    if f[0] in targets:
        obs = []
        for item in f[6:]:
            pos, cat, inj, err, oe, match, block = item.split(":")
            obs.append((int(pos), int(cat), int(inj), float(err), float(oe), int(match), int(block)))
        dumps[f[0]].append((int(f[2]), obs))


def best_margin(obs):
    per = collections.defaultdict(float)
    for pos, cat, inj, err, oe, match, block in obs:
        if err >= 0.3:
            continue
        e = min(0.45, 1 - (1 - err) * (1 - oe))
        per[block] += (1 if match == 1 else -1) * math.log((1 - e) / e)
    return max((abs(v) for v in per.values()), default=0.0)


read_cls = collections.Counter()
het_cls = collections.Counter()
margin_hist = collections.Counter()
single_err = collections.Counter()
for q, (lo, hi, hets) in targets.items():
    entries = dumps.get(q)
    if not entries:
        cls = "absent_from_chunks"
    elif any(h != 0 for h, _ in entries):
        cls = "labelled_in_chunk_then_lost"
    else:
        obs = max((o for _, o in entries), key=len)
        reliable = [o for o in obs if o[3] < 0.3]
        if not obs:
            cls = "in_chunk_no_observations"
        elif not reliable:
            cls = "only_unreliable_sites(err>=0.3)"
        else:
            m = best_margin(obs)
            cls = "posterior_below_0.9"
            margin_hist["<1" if m < 1 else "1-2" if m < 2 else "2-2.2" if m < 2.197 else ">=2.2"] += 1
            if len(reliable) == 1:
                single_err[f"{reliable[0][3]:.1f}"] += 1
    read_cls[cls] += 1
    # What happened to the truth hets in this read's span.
    for p in hets:
        snp = truth[p]
        rv = lookup(run_vcf, p, snp)
        if "phased" in rv:
            s = "run_vcf_phased"
        elif "unphased_het" in rv:
            s = "run_vcf_unphased_het"
        elif "hom" in rv:
            s = "run_vcf_hom"
        else:
            rc = lookup(run_cand, p, snp)
            bc = lookup(bam_cand, p, snp)
            cat = in_catalog(p)
            s = ("run_cand:" + rc[0][0]) if rc else ("bam_cand:" + bc[0][0]) if bc else "not_called"
            s += "|catalog" if cat else "|no_catalog"
        het_cls[(cls, "snp" if snp else "indel", s)] += 1

print("\n== reads")
for k, v in read_cls.most_common():
    print(f"{v}\t{k}")
print("\n== best block margin for posterior_below_0.9 (2.197 = posterior 0.9)")
for k in ("<1", "1-2", "2-2.2", ">=2.2"):
    print(f"{margin_hist[k]}\t{k}")
print("\n== learned site error when the read has one reliable observation")
for k, v in sorted(single_err.items()):
    print(f"{v}\terr {k}")
print("\n== truth hets inside these reads (read class, kind, fate)")
for k, v in het_cls.most_common(40):
    print(f"{v}\t{k[0]}\t{k[1]}\t{k[2]}")
