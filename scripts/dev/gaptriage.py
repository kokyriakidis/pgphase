#!/usr/bin/env python3
"""Every gap a rival closes, put in one class by measured evidence (evaluation only).

Usage: gaptriage.py DIAG_DIR GAPSITES_TSV RIVAL [--list CLASS] [--top N]

GAPSITES_TSV is `gapsites.py --kind all --out`. For each locus: the gap reads
(arms, unphased by us, correct by RIVAL, inside the locus), their MAPQ, the
depth at the locus against the contig's median, and the rival's sites the gap
reads cover. Classes, first match wins:
  low_mapq        at least half the gap reads have MAPQ < 30
  high_depth      depth >= 1.5x the contig median (a collapsed duplication)
  no_truth_site   no truth heterozygote among the covered rival sites
  repeat_only     every covered truth-het rival site is a homopolymer / STR
  clean:<stage>   a covered truth-het SNP or non-repeat indel; by our furthest
                  stage of the best such site
Prints reads and loci per class, and the largest loci of --list CLASS.
"""
import argparse
import collections
import os
import random
import statistics
import sys

import pysam

sys.path.insert(0, "scripts/dev")
import diaglib as D  # noqa: E402
import score  # noqa: E402

RANK = ["em_reliable", "em_unreliable", "em_unverified", "graph", "inj", "aln", "absent"]


def stage_rank(st):
    for i, p in enumerate(RANK):
        if st.startswith(p):
            return i
    return 9


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("diag")
    ap.add_argument("tsv")
    ap.add_argument("rival")
    ap.add_argument("--list", default="")
    ap.add_argument("--top", type=int, default=15)
    a = ap.parse_args()
    contig = os.environ["PGDEV_CONTIG"]
    loci = collections.defaultdict(list)
    for line in open(a.tsv):
        p = line.rstrip("\n").split("\t")
        if p[0] == "locus_beg":
            continue
        loci[(int(p[0]), int(p[1]))].append(dict(pos=int(p[3]), ref=p[4], alt=p[5], kind=p[6], tract=int(p[7]),
                                                 truth=p[8], stage=p[9], detail=p[10], cover=int(p[12])))
    scope = {q: v for q, v in score.read_scope().items() if v[1] == "arms_phaseable"}
    ours = score.read_classes(a.diag)
    theirs = score.read_classes(D.competitors()[a.rival])
    target = [q for q in scope if ours[q] == "unphased" and theirs[q] == "correct"]
    starts = sorted((scope[q][3], q) for q in target)
    bam = pysam.AlignmentFile(os.environ["PGDEV_BAM"])
    # contig median depth from random positions on the arms
    random.seed(3)
    clen = bam.get_reference_length(contig)
    depths = [bam.count(contig, p, p + 1) for p in (random.randrange(1_000_000, clen - 1_000_000) for _ in range(200))]
    median_depth = statistics.median(d for d in depths if d > 0)
    mapq = {}
    rows = []
    for (beg, end), sites in loci.items():
        reads = [q for s, q in starts if beg <= s <= end]
        if not reads:
            continue
        want = set(reads)
        for r in bam.fetch(contig, beg, end):
            if r.query_name in want and not r.is_secondary and not r.is_supplementary:
                mapq[r.query_name] = r.mapping_quality
        low = sum(1 for q in reads if mapq.get(q, 60) < 30) / len(reads)
        mid = (beg + end) // 2
        depth = bam.count(contig, mid, mid + 1)
        covered = [s for s in sites if s["cover"] > 0]
        truth_sites = [s for s in covered if s["truth"] == "het"]
        clean = [s for s in truth_sites if s["kind"] in ("snp", "indel")]
        if low >= 0.5:
            cls = "low_mapq"
        elif depth >= 1.5 * median_depth:
            cls = "high_depth"
        elif not truth_sites:
            cls = "no_truth_site"
        elif not clean:
            cls = "repeat_only"
        else:
            best = min(clean, key=lambda s: stage_rank(s["stage"]))
            st = best["stage"]
            cls = "clean:" + (st.split(":")[0] if st.startswith("inj_category") is False else st[:28])
        rows.append((len(reads), beg, end, cls, low, depth, covered))
    tab = collections.Counter(); nloci = collections.Counter()
    for n, beg, end, cls, low, depth, covered in rows:
        tab[cls] += n; nloci[cls] += 1
    total = sum(tab.values())
    print(f"{a.rival}: {len(rows)} loci, {total:,} gap reads (contig median depth {median_depth:.0f})")
    for cls, n in sorted(tab.items(), key=lambda kv: -kv[1]):
        print(f"  {cls:40} reads {n:5} ({100 * n / total:4.1f}%)  loci {nloci[cls]:4}")
    if a.list:
        print(f"\nlargest loci of {a.list}:")
        for n, beg, end, cls, low, depth, covered in sorted(rows, reverse=True):
            if not cls.startswith(a.list):
                continue
            print(f"  {beg:,}-{end:,}  reads {n}  lowMQ {low:.2f}  depth {depth}  sites: " +
                  "; ".join(f"{s['pos']}:{s['kind']}{s['tract'] or ''}:{s['truth']}:{s['stage']}" for s in covered[:6]))
            a.top -= 1
            if a.top <= 0:
                break


if __name__ == "__main__":
    main()
