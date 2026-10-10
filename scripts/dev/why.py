#!/usr/bin/env python3
"""Why are these arm reads not phased correctly? (evaluation only)

Usage: why.py RUN_DIR EM_DUMP_DIR COMPETITOR_DIR [--regions BED] [--out TSV]

RUN_DIR is a --bam run made with PGPHASE_EM_DUMP=EM_DUMP_DIR. Two read sets,
both arms_phaseable (see score.py):
  A  the competitor phases the read correctly and we do not;
  B  we leave the read unphased although it carries a good heterozygote: a
     biallelic truth SNP with no other truth variant within 10 bp, outside a
     homopolymer or short tandem repeat, read at base quality >= 20 with one
     of its two alleles.
Each read takes the first class that applies:
  not_in_solve    absent from the final EM solve (no graph alignment, and not
                  added as an alignment-only read)
  no_observation  in the solve without one usable allele call
  unreliable_only every observation sits on a site the EM rejects (error >= 0.3)
  weak_posterior  reliable evidence, but no block reaches posterior 0.8
                  (split: one site / several agreeing / conflicting)
  em_unphased_out EM labelled it; the output has no label
  wrong_side      labelled, and the read's own evidence disagrees with its
                  block's majority orientation (bad calls or a switch)
For B reads each good het is also located: not an EM site; an EM site the read
has no call at; unreliable; or a reliable call.
"""
import argparse
import bisect
import collections
import glob
import os
import sys

import pysam

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import score  # noqa: E402

ENV = os.environ
LOG4 = 1.3862943611198906  # posterior 0.8


def load_dump(path):
    reads = {}
    for f in glob.glob(os.path.join(path, "em.*.tsv")):
        with open(f) as fh:
            for line in fh:
                p = line.rstrip("\n").split("\t")
                obs = []
                for tok in p[5:]:
                    pos, kind, side, err, block = tok.split(":")
                    obs.append((int(pos), kind, int(side), float(err), int(block)))
                margins = {}
                if p[4]:
                    for kv in p[4].split(","):
                        k, v = kv.split("=")
                        margins[int(k)] = float(v)
                # A read in two chunks: keep the visit with more evidence.
                prev = reads.get(p[0])
                if prev is None or len(obs) > len(prev["obs"]):
                    reads[p[0]] = {"mapq": int(p[1]), "hap": int(p[2]), "ps": int(p[3]),
                                   "margins": margins, "obs": obs}
    return reads


def good_hets(truth_phased, ref):
    """Biallelic truth SNPs, isolated and outside short repeats: pos -> (ref, alt)."""
    contig = ENV["PGDEV_CONTIG"]
    keys = sorted(truth_phased)
    pos_all = [k[0] for k in keys]
    out = {}
    for (pos, r, a) in keys:
        if len(r) != 1 or len(a) != 1:
            continue
        i = bisect.bisect_left(pos_all, pos - 10)
        j = bisect.bisect_right(pos_all, pos + 10)
        if j - i > 1:
            continue
        ctx = ref.fetch(contig, pos - 9, pos + 8).upper()  # 1-based pos at index 8
        left, right = ctx[:8], ctx[9:]
        if any(s * 4 in left[-4:] + right[:4] for s in "ACGT"):
            continue
        rep = False
        for period in (1, 2, 3):
            for seq in (left[-3 * period:], right[:3 * period]):
                if len(seq) == 3 * period and seq == seq[:period] * 3:
                    rep = True
        if not rep:
            out[pos] = (r.upper(), a.upper())
    return out


def read_alleles(bam_path, names, hets):
    """name -> [(het pos, 'ref'|'alt')] for good hets read at BQ >= 20 with one allele."""
    hpos = sorted(hets)
    out = collections.defaultdict(list)
    with pysam.AlignmentFile(bam_path) as fh:
        for r in fh.fetch(until_eof=True):
            if r.query_name not in names or r.is_secondary or r.is_supplementary:
                continue
            lo = bisect.bisect_left(hpos, r.reference_start + 1)
            hi = bisect.bisect_right(hpos, r.reference_end)
            if lo == hi:
                continue
            want = set(hpos[lo:hi])
            seq, qual = r.query_sequence, r.query_qualities
            for qp, rp in r.get_aligned_pairs(matches_only=True):
                p1 = rp + 1
                if p1 in want and qual[qp] >= 20:
                    b = seq[qp].upper()
                    if b == hets[p1][0]:
                        out[r.query_name].append((p1, "ref"))
                    elif b == hets[p1][1]:
                        out[r.query_name].append((p1, "alt"))
    return out


def obs_bucket(n):
    return "1 obs" if n == 1 else "2-4 obs" if n <= 4 else "5-19 obs" if n < 20 else "20+ obs"


def site_purity(dump, scope):
    """pos -> (reads, purity, err, kind): how well a site's calls separate the parents."""
    calls = collections.defaultdict(lambda: [0, 0])
    meta = {}
    for q, d in dump.items():
        if q not in scope:
            continue
        mat = scope[q][2]
        for pos, kind, side, err, _ in d["obs"]:
            calls[pos][(side == 0) == mat] += 1
            meta[pos] = (err, kind)
    out = {}
    for pos, (a, b) in calls.items():
        out[pos] = (a + b, max(a, b) / (a + b), meta[pos][0], meta[pos][1])
    return out


def classify(d, ours, maternal):
    if d is None:
        return "not_in_solve", ""
    if not d["obs"]:
        return "no_observation", ""
    reliable = [o for o in d["obs"] if o[3] < 0.3]
    if not reliable:
        return "unreliable_only", obs_bucket(len(d["obs"]))
    if d["hap"] == 0:
        best = max(abs(m) for m in d["margins"].values()) if d["margins"] else 0.0
        sides = collections.Counter()
        for o in reliable:
            sides[(o[4], o[2])] += 1
        blocks = collections.defaultdict(set)
        for (blk, side) in sides:
            blocks[blk].add(side)
        if len(reliable) == 1:
            sub = "one site"
        elif best < LOG4 and all(len(s) == 1 for s in blocks.values()):
            sub = "agreeing, low margin"
        else:
            sub = "conflicting"
        return "weak_posterior", f"{sub}; best margin {best:.2f}; {len(reliable)} reliable"
    if ours == "unphased":
        return "em_unphased_out", ""
    return "wrong_side", obs_bucket(len(reliable)) + " reliable"


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("run")
    ap.add_argument("dump")
    ap.add_argument("competitor")
    ap.add_argument("--regions")
    ap.add_argument("--out")
    a = ap.parse_args()
    regions = score.load_regions(a.regions)
    scope = {q: v for q, v in score.read_scope(regions).items() if v[1] == "arms_phaseable"}
    ours = score.classify_reads(a.run, scope)
    theirs = score.classify_reads(a.competitor, scope)
    t = score.truth_tables()
    dump = load_dump(a.dump)
    site_pos = {o[0] for d in dump.values() for o in d["obs"]}
    set_a = {q for q in scope if theirs[q] == "correct" and ours[q] != "correct"}
    unphased = {q for q in scope if ours[q] == "unphased"}
    ref = pysam.FastaFile(ENV["PGDEV_REF"])
    hets = good_hets(t["truth_phased"], ref)
    carried = read_alleles(ENV["PGDEV_BAM"], unphased, hets)
    set_b = {q for q in unphased if carried.get(q)}

    rows = []
    for name, reads in (("A", set_a), ("B", set_b)):
        tab = collections.Counter()
        sub = collections.Counter()
        mapq = collections.Counter()
        het_loc = collections.Counter()
        ex = collections.defaultdict(list)
        for q in sorted(reads, key=lambda q: scope[q][3]):
            d = dump.get(q)
            cls, detail = classify(d, ours[q], scope[q][2])
            tab[cls] += 1
            if detail:
                sub[(cls, detail.split(";")[0])] += 1
            m = d["mapq"] if d else -1
            mapq[(cls, "MAPQ>=20" if m >= 20 else ("MAPQ<20" if m >= 0 else "not in solve"))] += 1
            het_detail = ""
            if name == "B":
                obs_at = {o[0]: o for o in (d["obs"] if d else [])}
                locs = []
                for p, allele in carried[q]:
                    # EM keys a catalog SNP at its VCF position or one base left.
                    o = obs_at.get(p, obs_at.get(p - 1))
                    if o is None:
                        locs.append("no_call" if p in site_pos or p - 1 in site_pos else "not_em_site")
                    elif o[3] >= 0.3:
                        locs.append("unreliable")
                    else:
                        locs.append("reliable")
                het_detail = ",".join(f"{p}:{l}" for (p, _), l in zip(carried[q], locs))
                het_loc[min(locs, key=["reliable", "unreliable", "no_call", "not_em_site"].index)] += 1
            if len(ex[cls]) < 3:
                ex[cls].append(f"{q.split('/')[1]}@{scope[q][3]:,}")
            rows.append((name, q, scope[q][3], scope[q][4], ours[q], theirs[q], cls, detail, m, het_detail))
        total = sum(tab.values())
        print(f"\nSet {name}: {total:,} reads "
              + ("(competitor correct, ours not)" if name == "A" else "(ours unphased, carry a good het)"))
        if name == "A":
            print(f"  ours: {collections.Counter(ours[q] for q in reads)}")
        for cls, n in tab.most_common():
            mq = ", ".join(f"{k[1]} {v}" for k, v in mapq.items() if k[0] == cls)
            print(f"  {cls:16} {n:6,}  ({mq})  e.g. {' '.join(ex[cls])}")
            for (c, s), v in sub.most_common():
                if c == cls:
                    print(f"      {v:6,}  {s}")
        if name == "B":
            print("  best location of a carried good het: " + ", ".join(f"{k} {v:,}" for k, v in het_loc.most_common()))
    # Sites behind set-A reads that EM rejects or weighs too lightly.
    full_scope = score.read_scope(regions)
    purity = site_purity(dump, full_scope)
    hinge = collections.Counter()
    hinge_sites = collections.defaultdict(set)
    for q in set_a:
        d = dump.get(q)
        if d is None or d["hap"] != 0:
            continue
        for pos, kind, side, err, _ in d["obs"]:
            if err >= 0.2:
                n, pur, e, k = purity[pos]
                good = "separates parents (>=0.9)" if pur >= 0.9 else "does not separate parents"
                hinge[(k, "unreliable" if e >= 0.3 else "weak", good)] += 1
                hinge_sites[(k, "unreliable" if e >= 0.3 else "weak", good)].add(pos)
    print("\nSet A unphased reads: observations at sites with error >= 0.2, by site kind and truth purity")
    print("  (kind C clean het, O other catalog, I injected from the alignment, L window)")
    for key, n in hinge.most_common():
        print(f"  {key[0]} {key[1]:10} {key[2]:28} {n:5,} obs at {len(hinge_sites[key]):4,} sites")
    tight = [p for p, v in purity.items() if v[1] >= 0.95 and v[0] >= 10 and v[2] >= 0.3]
    print(f"  sites with >= 10 truth reads, purity >= 0.95, yet error >= 0.3: {len(tight):,}"
          + (f"  e.g. {sorted(tight)[:6]}" if tight else ""))

    if a.out:
        with open(a.out, "w") as fh:
            fh.write("set\tread\tstart\tend\tours\tcompetitor\tclass\tdetail\tmapq\tgood_hets\n")
            for r in rows:
                fh.write("\t".join(map(str, r)) + "\n")
        print(f"\nreads: {a.out}")


if __name__ == "__main__":
    sys.exit(main())
