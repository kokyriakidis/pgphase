#!/usr/bin/env python3
"""Catalog sites whose decomposed rows were all filtered although two
non-reference alleles split the reads (a 1/2 heterozygote). (evaluation only)
Usage: lost_multiallelic.py FILTERED_TSV CANDIDATES_TSV TRUTH_VCF"""
import bisect
import collections
import gzip
import sys

filt_path, cand_path, truth_vcf = sys.argv[1:4]
CEN = (26_000_000, 32_000_000)

sites = collections.defaultdict(list)  # (pos, base id) -> rows
for line in open(filt_path):
    if line.startswith("CHROM"):
        continue
    f = line.rstrip("\n").split("\t")
    sid = f[2]
    base = sid.rsplit(":", 1)[0] if ":" in sid else sid
    sites[(int(f[1]), base)].append((sid, int(f[3]), int(f[4]), f[7]))

kept_pos = set()
with open(cand_path) as fh:
    next(fh)
    for line in fh:
        kept_pos.add(int(line.split("\t", 2)[1]))

thet = []
for line in gzip.open(truth_vcf, "rt"):
    if line[0] == "#":
        continue
    f = line.split("\t")
    gt = f[9].split(":")[0].replace("/", "|").split("|")
    if len(gt) == 2 and gt[0] != gt[1] and "." not in gt:
        thet.append(int(f[1]))
thet.sort()


def truth_het_in(lo, hi):
    return bisect.bisect_left(thet, hi + 1) > bisect.bisect_left(thet, lo)


tab = collections.Counter()
examples = []
for (pos, base), rows in sites.items():
    if len(rows) < 2:
        continue
    ref = max(r[1] for r in rows)
    alts = sorted((r[2] for r in rows), reverse=True)
    depth = ref + sum(alts)
    if depth < 10:
        continue
    # Two leading alleles (REF counts as one) each carrying >= 25% of the site's reads.
    top = sorted([ref] + alts, reverse=True)
    if top[1] < 0.25 * depth or ref >= top[1]:
        continue  # not a two-non-REF split
    region = "cen" if CEN[0] <= pos < CEN[1] else "arms"
    covered = any(p in kept_pos for p in range(pos - 25, pos + 26))
    truth = truth_het_in(pos - 5, pos + 60)
    tab[(region, "truth_het" if truth else "no_truth_het", "something_kept_nearby" if covered else "nothing_kept")] += 1
    if region == "arms" and truth and not covered and len(examples) < 8:
        examples.append((pos, base, ref, alts[:4], sorted(set(r[3] for r in rows))))
for k in sorted(tab):
    print(k, tab[k])
for e in examples:
    print("  ", e)
