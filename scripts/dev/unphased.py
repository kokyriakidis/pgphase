#!/usr/bin/env python3
"""Every arms_phaseable read we leave unphased, in one primary cause (evaluation only).

Usage: unphased.py DIAG_DIR [--loci N] [--list CAUSE]

Evidence per read: its calls at every candidate as the final EM saw them
(emreads.*, candidates the EM left out included), each judged without
comparing allele spellings: at a position, the allele most maternal reads
carry and the one most paternal reads carry (all scored reads); a call is
right when it is its parent's allele, wrong when it is the other parent's,
and the position counts only where the two parents' alleles differ. The
position's state: reliable EM site (learned error < 0.3), rejected EM site,
or not an EM site. Primary cause, first match wins:
  skipped                    the read is skipped (filters)
  no call                    no allele at any candidate; split by the truth
                             hets in its span (snp / indel / repeat) and our
                             furthest stage at them
  right call, reliable site  a right call at a reliable EM site, yet no label:
                             split by whether it also has wrong reliable calls
  right call, rejected site  right calls only at rejected EM sites: island
                             claimed / too few unlabelled callers / site error
                             above the island cut
  right call, non-EM site    right calls only at candidates outside the EM
                             (repeat het indels, island-only, other categories)
  no clean parental split    calls only where the parents' alleles do not
                             differ, or split at under 80% purity (noise)
  wrong calls only
Each cause is crossed with HiPhase's and LongPhase's class for the read, and
its reads are grouped into loci (reads within 5 kb) to show concentration.
"""
import argparse
import bisect
import collections
import os
import sys

import pysam

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import diaglib as D  # noqa: E402
import score  # noqa: E402
from readaudit import site_kind  # noqa: E402


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("diag")
    ap.add_argument("--loci", type=int, default=4)
    ap.add_argument("--list", default="")
    a = ap.parse_args()
    contig = os.environ["PGDEV_CONTIG"]
    fa = pysam.FastaFile(os.environ["PGDEV_REF"])
    ref_name = contig if contig in fa.references else contig.split("#")[-1]
    scope_all = score.read_scope()
    scope = {q: v for q, v in scope_all.items() if v[1] == "arms_phaseable"}
    ours = score.read_classes(a.diag)
    U = {q for q in scope if ours[q] == "unphased"}
    rivals = {name: score.read_classes(d) for name, d in D.competitors().items()}

    # every scored read's calls (emreads), keeping the copy with most calls
    calls, mapq, skipped = {}, {}, {}
    for f in D.chunk_files(os.path.join(a.diag, "aln"), "emreads"):
        with open(f) as fh:
            for line in fh:
                p = line.rstrip("\n").split("\t")
                if p[0] not in scope_all:
                    continue
                mine = {}
                for t in p[3:]:
                    if ":" not in t:
                        continue
                    pos, al = t.rsplit(":", 1)
                    if al != "-1":
                        mine[int(pos)] = int(al)
                if p[0] not in calls or len(mine) > len(calls[p[0]]):
                    calls[p[0]] = mine
                    mapq[p[0]] = int(p[1])
                    skipped[p[0]] = p[2] == "1"
    want = set()
    for q in U:
        want.update(calls.get(q, {}))
    par = collections.defaultdict(lambda: collections.defaultdict(lambda: [0, 0]))  # pos -> allele -> [mat, pat]
    for q, mine in calls.items():
        mat = scope_all[q][2]
        for pos, al in mine.items():
            if pos in want:
                par[pos][al][0 if mat else 1] += 1

    def parent_alleles(pos):
        """(maternal allele, paternal allele) where the parents split cleanly:
        each parent's majority allele differs and together they carry at
        least 80% of the calls (a noisy split, e.g. homopolymer length noise
        at a real het, is not evidence a read could be phased by)."""
        d = par.get(pos)
        if not d:
            return None
        m = max(d, key=lambda al: d[al][0])
        pa = max(d, key=lambda al: d[al][1])
        if d[m][0] < 2 or d[pa][1] < 2 or m == pa:
            return None
        total = sum(v[0] + v[1] for v in d.values())
        if (d[m][0] + d[pa][1]) / total < 0.8:
            return None
        return m, pa

    # EM state per position and the EM's labels (unlabelled callers per site)
    em_err = {}
    em_label = {}
    em_margin = {}   # read -> (best |margin|, number of blocks with evidence)
    em_obs = {}      # read -> [(pos, kind, side, site err, read err)]
    for f in D.chunk_files(os.path.join(a.diag, "em"), "em"):
        with open(f) as fh:
            for line in fh:
                p = line.rstrip("\n").split("\t")
                if p[0] in U:
                    ms = [abs(float(x.split("=")[1])) for x in p[4].split(",") if "=" in x]
                    first0 = 6 if len(p) > 5 and ":" not in p[5] else 5
                    obs = []
                    for t in p[first0:]:
                        x = t.split(":")
                        obs.append((int(x[0]), x[1], x[2], float(x[3]), float(x[6]) if len(x) > 6 else 0.0))
                    if p[0] not in em_obs or len(obs) > len(em_obs[p[0]]):
                        em_obs[p[0]] = obs
                        em_margin[p[0]] = (max(ms) if ms else 0.0, len(ms))
                em_label.setdefault(p[0], int(p[2]))
                if int(p[2]) in (1, 2):
                    em_label[p[0]] = int(p[2])
                first = 6 if len(p) > 5 and ":" not in p[5] else 5
                for t in p[first:]:
                    x = t.split(":")
                    if x[1] not in ("L",):
                        em_err[int(x[0])] = (x[1], float(x[3]))
    unl_callers = collections.Counter()
    for q, mine in calls.items():
        if em_label.get(q, 0) == 0:
            for pos in mine:
                if pos in em_err:
                    unl_callers[pos] += 1
    claimed = {}
    for f in D.chunk_files(os.path.join(a.diag, "aln"), "islands"):
        with open(f) as fh:
            for line in fh:
                p = line.rstrip("\n").split("\t")
                claimed[p[0]] = p[10] if len(p) > 10 else "E"
    truth = []
    for pos, recs in D.vcf_records(os.environ["PGDEV_TRUTH_VCF"], contig, 1, 2 ** 31 - 1).items():
        for ref, alts, gt, ps, c in recs:
            if len(gt) == 2 and None not in gt and gt[0] != gt[1]:
                truth.append((pos, ref, alts[max(gt) - 1]))
    truth.sort()
    tpos = [t[0] for t in truth]
    graph, inj, aln = D.load_graph(a.diag), D.load_inj(a.diag), D.load_aln(a.diag)

    def stage(pos, tol):
        for d in range(tol + 1):
            for p in ((pos,) if d == 0 else (pos - d, pos + d)):
                if p in em_err:
                    k, e = em_err[p]
                    return "em_unverified" if k == "U" else "em_unreliable" if e >= 0.3 else "em_reliable"
        g = D.near(graph, pos, tol)
        if g:
            rows = g[0][1]
            return "graph_voting" if any(r["voting"] for r in rows) else "graph_" + rows[0]["cat"]
        i = D.near(inj, pos, tol)
        if i:
            r = min(i[0][1], key=lambda r: D.decision_rank(r["decision"]))
            return "inj_" + (f"category:{r['cat']}" if r["decision"] == "category" else r["decision"])
        al = D.near(aln, pos, tol)
        if al:
            r = max(al[0][1], key=lambda r: r["alt_cov"])
            return f"aln_{r['init']}>{r['final']}"
        return "absent"

    tab = collections.Counter()
    loci = collections.defaultdict(list)
    detail = collections.Counter()
    for q in U:
        s, e = scope[q][3], scope[q][4]
        mine = calls.get(q, {})
        right_rel = wrong_rel = right_rej = right_out = wrong_any = nosplit = 0
        rej_detail = set()
        for pos, al in mine.items():
            pa = parent_alleles(pos)
            if pa is None:
                nosplit += 1
                continue
            right = al == (pa[0] if scope[q][2] else pa[1])
            wrong = al == (pa[1] if scope[q][2] else pa[0])
            st = em_err.get(pos)
            if right:
                if st is None:
                    right_out += 1
                elif st[0] != "U" and st[1] < 0.3:
                    right_rel += 1
                else:
                    right_rej += 1
                    if st[1] >= 0.40:
                        rej_detail.add("site error >= 0.40")
                    elif unl_callers[pos] < 2:
                        rej_detail.add("lone caller")
                    else:
                        rej_detail.add("island candidate")
            elif wrong:
                wrong_any += 1
                if st is not None and st[0] != "U" and st[1] < 0.3:
                    wrong_rel += 1
        if skipped.get(q):
            cause = "skipped"
        elif not mine:
            i, j = bisect.bisect_left(tpos, s), bisect.bisect_right(tpos, e)
            kinds = collections.Counter()
            stages = collections.Counter()
            for pos, ref, alt in truth[i:j]:
                kind, tract = site_kind(fa, ref_name, pos, ref, alt)
                kinds["repeat" if kind in ("hp", "str") else kind] += 1
                st = stage(pos, 0 if kind == "snp" else max(3, tract))
                stages[st.split(":")[0] if st.startswith("inj_category") else st] += 1
            kk = "/".join(sorted(kinds))
            cause = f"no call ({kk})"
            detail[(cause, stages.most_common(1)[0][0] if stages else "-")] += 1
        elif right_rel:
            cause = "right call, reliable site" + (" + wrong reliable" if wrong_rel else "")
        elif right_rej:
            cause = "right call, rejected site: " + ("claimed by an island" if q in claimed else
                                                    " / ".join(sorted(rej_detail)))
        elif right_out:
            cause = "right call, non-EM site" + (" (claimed)" if q in claimed else "")
        elif wrong_any:
            cause = "wrong calls only"
        else:
            cause = "no clean parental split"
        hc = rivals.get("hiphase", {}).get(q, "?")
        lc = rivals.get("longphase", {}).get(q, "?")
        tab[(cause, hc, lc)] += 1
        loci[cause].append((s, e, q))
    # why reads with a right call at a reliable site stay unlabelled
    why = collections.Counter()
    margins = []
    for c, lst in loci.items():
        if not c.startswith("right call, reliable site"):
            continue
        for s, e, q in lst:
            m, nb = em_margin.get(q, (0.0, 0))
            margins.append(m)
            rel = [o for o in em_obs.get(q, []) if o[1] not in ("U",) and o[3] < 0.3]
            nrel = len({o[0] for o in rel})
            blocks = nb
            mb = "margin < 0.5" if m < 0.5 else "margin 0.5-1.39" if m < 1.386 else "margin >= 1.39 (posterior >= 0.8)"
            why[(mb, "1 reliable call" if nrel == 1 else "2-3 reliable calls" if nrel <= 3 else "4+ reliable calls",
                 "1 block" if blocks <= 1 else "2+ blocks")] += 1
    # every unphased read whose only reliable evidence is one call: right / wrong
    single = collections.Counter()
    for q in U:
        rel = {}
        for o in em_obs.get(q, []):
            if o[1] not in ("U", "L") and o[3] < 0.3:
                rel[o[0]] = o
        if len(rel) != 1:
            continue
        pos, o = next(iter(rel.items()))
        al = calls.get(q, {}).get(pos)
        pa = parent_alleles(pos)
        if al is None or pa is None:
            verdict = "no split"
        else:
            verdict = "right" if al == (pa[0] if scope[q][2] else pa[1]) else \
                      "wrong" if al == (pa[1] if scope[q][2] else pa[0]) else "other allele"
        m = em_margin.get(q, (0.0, 0))[0]
        mb = "<0.5" if m < 0.5 else "0.5-0.85" if m < 0.85 else "0.85-1.1" if m < 1.1 else "1.1-1.39" if m < 1.386 else ">=1.39"
        eb = "<0.1" if o[3] < 0.1 else "0.1-0.2" if o[3] < 0.2 else "0.2-0.3"
        single[(mb, eb, o[1], verdict)] += 1
    print(f"{contig}: arms_phaseable unphased reads {len(U):,}")
    causes = sorted({k[0] for k in tab}, key=lambda c: -sum(v for k, v in tab.items() if k[0] == c))
    print(f"  {'primary cause':58} {'reads':>6} {'HiPhase correct':>16} {'LongPhase correct':>18} {'either':>7}")
    for c in causes:
        n = sum(v for k, v in tab.items() if k[0] == c)
        hcor = sum(v for k, v in tab.items() if k[0] == c and k[1] == "correct")
        lcor = sum(v for k, v in tab.items() if k[0] == c and k[2] == "correct")
        either = sum(v for k, v in tab.items() if k[0] == c and "correct" in (k[1], k[2]))
        print(f"  {c:58} {n:6,} {hcor:16,} {lcor:18,} {either:7,}")
    if single:
        print("  unphased reads with exactly one reliable call: by margin, site error, site kind -> right / wrong / no split")
        keys = sorted({k[:3] for k in single})
        for k in keys:
            r, w, n = single[k + ("right",)], single[k + ("wrong",)], single[k + ("no split",)] + single[k + ("other allele",)]
            print(f"    margin {k[0]:9} site err {k[1]:8} kind {k[2]:2}  right {r:4}  wrong {w:4}  other {n:4}")
        tot = collections.Counter()
        for k, v in single.items():
            tot[(k[0], k[3])] += v
        print("    by margin: " + "; ".join(f"{mb}: {tot[(mb, 'right')]} right / {tot[(mb, 'wrong')]} wrong"
                                        for mb in ("<0.5", "0.5-0.85", "0.85-1.1", "1.1-1.39", ">=1.39")))
    if why:
        print("  'right call, reliable site' reads: best block |log-likelihood margin| (labelled at >= 1.39), reliable calls, blocks")
        for k, v in why.most_common():
            print(f"    {k[0]:36} {k[1]:20} {k[2]:10} {v:5,}")
    if detail:
        print("  no-call reads, our furthest stage at their truth hets (most common per read):")
        for (c, st), v in detail.most_common(12):
            print(f"    {c:32} {st:40} {v:5,}")
    print(f"  loci (reads within 5 kb), largest {a.loci} per cause:")
    for c in causes:
        rows = sorted(loci[c])
        groups, cur = [], None
        for s, e, q in rows:
            if cur and s <= cur[1] + 5000:
                cur[1] = max(cur[1], e)
                cur[2] += 1
            else:
                cur = [s, e, 1]
                groups.append(cur)
        groups.sort(key=lambda g: -g[2])
        print(f"    {c:58} " + ", ".join(f"{g[0]:,}-{g[1]:,} ({g[2]})" for g in groups[:a.loci]) +
              f"  [{len(groups)} loci]")
    if a.list:
        for c in causes:
            if c.startswith(a.list):
                for s, e, q in sorted(loci[c]):
                    print(f"  {q}\t{s}\t{e}\t{c}")


if __name__ == "__main__":
    main()
