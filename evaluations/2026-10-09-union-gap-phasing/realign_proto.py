#!/usr/bin/env python3
"""Prototype: realignment allele calls at sites the EM learned as unreliable.

Usage: realign_proto.py RUN_DIR DIAG_READS [FLANK]
Sites: positions the final EM gave error >= 0.3 (dump), arms only, with their
REF/ALT from RUN_DIR/phased.vcf. For each primary MAPQ>=5 read spanning the
site's window (variant +- tandem-repeat extension + FLANK), the read's segment
is aligned (edlib NW) to the reference haplotype and to the ALT haplotype; the
closer one is its allele (ties abstain). Per site, the calls are scored against
parental truth: purity = fraction of called reads consistent with the best
allele<->parent pairing (evaluation only). Compared with the EM's learned error.
"""
import bisect
import collections
import sys

import edlib
import pysam

run, diag_path = sys.argv[1:3]
FLANK = int(sys.argv[3]) if len(sys.argv) > 3 else 20
MODE = sys.argv[4] if len(sys.argv) > 4 else "unreliable"
TRUTH_VCF = "/home/kokyriakidis/Downloads/pgphase-eval-data/results/chr12-18-20-comparison/chr20/truth.vcf.gz"
ROOT = "/home/kokyriakidis/Downloads/pgphase/test_data"
BAM = f"{ROOT}/HG002_chr20_hifi_mapped_to_CHM13_chr20_annotated.bam"
REF = f"{ROOT}/chm13v2.0.chr20.renamed.fa"
CONTIG = "CHM13#0#chr20"
CEN = (26_000_000, 32_000_000)

truth = {}
for line in open(f"{ROOT}/derived/chr20_truth_hap.tsv"):
    q, h = line.rstrip("\n").split("\t")[:2]
    if h in ("MATERNAL", "PATERNAL"):
        truth[q] = h
ref = pysam.FastaFile(REF).fetch(CONTIG).upper()

site_err = {}
for line in open(diag_path):
    f = line.rstrip("\n").split("\t")
    for item in f[6:]:
        pos, cat, inj, err, oe, match, block = item.split(":")
        if cat != "-1":
            site_err[int(pos)] = float(err)

vcf = {}
with open(f"{run}/candidates.tsv") as fh:
    next(fh)
    for line in fh:
        f = line.split("\t")
        if f[14] not in ("NOISY_CAND_HET", "CLEAN_HET_INDEL", "CLEAN_HET_SNP", "REP_HET_INDEL"):
            continue
        vcf.setdefault(int(f[1]), (f[3].upper(), f[4].split(",")[0].upper(), f[14]))


def repeat_extent(start):
    best = 0
    for period in range(1, 51):
        k = start + period
        while k - start < 500 and k < len(ref) and ref[k] == ref[k - period]:
            k += 1
        if k - start >= 2 * period:
            best = max(best, k - start)
    return best


def query_pos(read, ref_pos):
    if ref_pos < read.reference_start or ref_pos >= read.reference_end:
        return None
    r, q = read.reference_start, 0
    for op, length in read.cigartuples:
        if op in (0, 7, 8):
            if r + length > ref_pos:
                return q + (ref_pos - r)
            r += length
            q += length
        elif op in (2, 3):
            if r + length > ref_pos:
                return q
            r += length
        elif op in (1, 4):
            q += length
    return None


def dist(a, b):
    return edlib.align(a, b, mode="NW", task="distance")["editDistance"]


import gzip, random
thet = []
for line in gzip.open(TRUTH_VCF, "rt"):
    if line[0] == "#":
        continue
    f = line.split("\t")
    gt = f[9].split(":")[0].replace("/", "|").split("|")
    if len(gt) == 2 and gt[0] != gt[1] and "." not in gt:
        thet.append(int(f[1]))
thet.sort()
def truth_het_near(p):
    i = bisect.bisect_left(thet, p - 10)
    return i < len(thet) and thet[i] <= p + 10
if MODE == "rep":
    pool = sorted(p for p, v in vcf.items() if v[2] == "REP_HET_INDEL" and not CEN[0] <= p < CEN[1])
    random.seed(7)
    targets = sorted(random.sample(pool, min(600, len(pool))))
    site_err = {p: -1.0 for p in targets}
else:
    targets = sorted(p for p, e in site_err.items() if e >= 0.3 and p in vcf and not CEN[0] <= p < CEN[1])
print(f"unreliable sites on the arms with a VCF record: {len(targets)}")
bam = pysam.AlignmentFile(BAM)
hist = collections.Counter()
by_err = collections.defaultdict(list)
examples = []
for p in targets:
    r_allele, a_allele, info = vcf[p]
    p0 = p - 1
    end0 = p0 + len(r_allele)
    w0 = max(0, p0 - FLANK)
    w1 = end0 + repeat_extent(end0) + FLANK
    hap_ref = ref[w0:w1]
    hap_alt = ref[w0:p0] + a_allele + ref[end0:w1]
    if hap_ref == hap_alt:
        continue
    calls = collections.Counter()
    margins = collections.Counter()
    for r in bam.fetch(CONTIG, w0, w1):
        if r.is_secondary or r.is_supplementary or r.is_duplicate or r.mapping_quality < 5:
            continue
        if r.reference_start > w0 or r.reference_end < w1 or r.query_name not in truth:
            continue
        q0, q1 = query_pos(r, w0), query_pos(r, w1 - 1)
        if q0 is None or q1 is None or q1 < q0:
            continue
        seg = r.query_sequence[q0:q1 + 1]
        d0, d1 = dist(seg, hap_ref), dist(seg, hap_alt)
        if d0 == d1:
            calls[("tie", truth[r.query_name])] += 1
            continue
        calls[(0 if d0 < d1 else 1, truth[r.query_name])] += 1
        margins[min(abs(d0 - d1), 3)] += 1
    cis = calls[(0, "MATERNAL")] + calls[(1, "PATERNAL")]
    trans = calls[(1, "MATERNAL")] + calls[(0, "PATERNAL")]
    n = cis + trans
    ties = calls[("tie", "MATERNAL")] + calls[("tie", "PATERNAL")]
    if n == 0:
        hist["no_calls"] += 1
        continue
    purity = max(cis, trans) / n
    bucket = ">=0.95" if purity >= 0.95 else "0.85-0.95" if purity >= 0.85 else "0.7-0.85" if purity >= 0.7 else "<0.7"
    hist[bucket] += 1
    hist[("truth_het" if truth_het_near(p) else "not_truth", bucket)] += 1
    by_err[bucket].append(site_err[p])
    if len(examples) < 12:
        examples.append((p, r_allele[:12], a_allele[:12], round(site_err[p], 2), n, ties, round(purity, 2), dict(margins)))
print(f"realigned per-site purity (FLANK={FLANK}):")
for k in (">=0.95", "0.85-0.95", "0.7-0.85", "<0.7", "no_calls"):
    print(f"   {k:10s} {hist[k]}")
for t in ("truth_het", "not_truth"):
    print(t, {k: hist[(t, k)] for k in (">=0.95", "0.85-0.95", "0.7-0.85", "<0.7")})
for e in examples[:6]:
    print("  ", e)
