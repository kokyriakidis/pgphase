#!/usr/bin/env python3
"""Gaps we leave that a rival closes (evaluation only).

Usage: gapclose.py DIAG_DIR [--top N] [--out TSV]

Arm reads crossing a truth het that we leave unphased and HiPhase or
LongPhase phases correctly are merged into loci (reads overlapping, 5 kb
slack). Each locus, ranked by reads:
  reads     reads recovered by a rival (per rival: correct)
  ours      our phase sets of the phased reads overlapping the locus's left
            and right ends: "break" when they differ (we cut a block there),
            "hole" when one block surrounds the locus, "edge" when a side has
            no phased read
  rival     per rival, the phase set holding most of the locus's reads, its
            span (first-last phased variant), and how many of the reads it
            holds
  why       why.py's reason per read on our side (unreliable_only,
            no_observation, weak_posterior, not_in_solve)
"""
import argparse
import bisect
import collections
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import diaglib as D  # noqa: E402
import score  # noqa: E402
import why  # noqa: E402


def ps_spans(vcf_path, contig):
    spans = {}
    for pos, recs in D.vcf_records(vcf_path, contig, 1, 2 ** 31 - 1).items():
        for ref, alts, gt, ps, c in recs:
            if ps:
                lo, hi = spans.get(ps, (pos, pos))
                spans[ps] = (min(lo, pos), max(hi, pos))
    return spans


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("diag")
    ap.add_argument("--top", type=int, default=25)
    ap.add_argument("--out")
    a = ap.parse_args()
    contig = os.environ["PGDEV_CONTIG"]
    comps = D.competitors()
    full = score.read_scope()
    scope = {q: v for q, v in full.items() if v[1] == "arms_phaseable"}
    ours = score.read_classes(a.diag)
    theirs = {n: score.read_classes(d) for n, d in comps.items()}
    target = sorted((q for q in scope if ours[q] == "unphased" and any(t[q] == "correct" for t in theirs.values())),
                    key=lambda q: scope[q][3])
    loci, cur = [], None
    for q in target:
        s, e = scope[q][3], scope[q][4]
        if cur and s <= cur["end"] + 5000:
            cur["end"] = max(cur["end"], e)
            cur["reads"].append(q)
        else:
            cur = {"beg": s, "end": e, "reads": [q]}
            loci.append(cur)
    # our tags for all scored reads (to find flanking blocks), rivals' for the targets
    our_tags = score.read_tags(a.diag, set(full))
    by_start = sorted((v[3], v[4], q) for q, v in full.items() if q in our_tags)
    starts = [x[0] for x in by_start]
    riv_tags = {n: score.read_tags(d, set(target)) for n, d in comps.items()}
    riv_spans = {n: ps_spans(os.path.join(d, "phased.vcf"), contig) for n, d in comps.items()}
    reads_em, _ = D.load_em(a.diag)

    def flank_ps(pos):
        i = bisect.bisect_right(starts, pos)
        c = collections.Counter()
        for s, e, q in by_start[max(0, i - 400):i]:
            if e >= pos:
                c[our_tags[q][1]] += 1
        return c.most_common(1)[0][0] if c else None

    rows = []
    for l in loci:
        rs = l["reads"]
        left, right = flank_ps(l["beg"]), flank_ps(l["end"])
        kind = "edge" if left is None or right is None else ("hole" if left == right else "break")
        rcell = []
        for n in comps:
            ps = collections.Counter(riv_tags[n][q][1] for q in rs if q in riv_tags[n] and theirs[n][q] == "correct")
            if ps:
                p, k = ps.most_common(1)[0]
                lo, hi = riv_spans[n].get(str(p), (0, 0))
                rcell.append(f"{n} {sum(theirs[n][q] == 'correct' for q in rs)} in PS {p} ({lo:,}-{hi:,}, holds {k})")
            else:
                rcell.append(f"{n} 0")
        whys = collections.Counter(why.classify(reads_em.get(q), "unphased", scope[q][2])[0] for q in rs)
        rows.append((len(rs), l["beg"], l["end"], kind, left, right, rcell, dict(whys)))
    rows.sort(reverse=True)
    tot = collections.Counter(r[3] for r in rows for _ in range(r[0]))
    print(f"unphased reads a rival phases correctly: {len(target):,} in {len(loci):,} loci "
          f"(break {tot['break']:,}, hole {tot['hole']:,}, edge {tot['edge']:,} reads)")
    for n, beg, end, kind, left, right, rcell, whys in rows[:a.top]:
        print(f"\n{beg:,}-{end:,} ({end - beg:,} bp)  reads {n}  {kind}  ours left PS {left} right PS {right}")
        for c in rcell:
            print(f"    {c}")
        print(f"    why: " + ", ".join(f"{k} {v}" for k, v in sorted(whys.items(), key=lambda kv: -kv[1])))
    if a.out:
        with open(a.out, "w") as fh:
            for r in rows:
                fh.write("\t".join(map(str, r)) + "\n")


if __name__ == "__main__":
    main()
