#!/usr/bin/env python3
"""HiPhase's allele assignment against ours, read by read (evaluation only).

Usage: realign.py DIAG_DIR [--vcf DV_VCF] [--max-sites N] [--min-run BP] [--seed S]

HiPhase (default settings, global realignment; src/read_parsing.rs,
src/wfa_graph.rs) aligns each read's aligned span by unit-cost edit distance
to a graph of the reference with every heterozygous call (any non-missing,
non-hom-ref genotype, no FILTER check) and every homozygous-ALT call as a
bubble; a variant's allele is the branch the best path traverses, and when
equally good paths traverse both branches it is ambiguous (no call). This
reproduces that in a window: the site's repeat +- 50 bp, every call fully
inside as a bubble, every combination of non-overlapping bubbles as a path,
minimum edit distance (Biopython global alignment, unit costs) of the read's
segment between the window ends, and the target's allele from the best paths.

Sites: heterozygous indels in the rival's calls (default the DeepVariant VCF
HiPhase phased) with a reference repeat run >= --min-run bp, grouped by
truth (a truth het within max(3, run) bp) and our EM state there (rejected:
learned error >= 0.3, reliable, none). Reads: scored primary alignments
spanning the window. A call is right when it is the side the read's parent's
reads mostly carry (each method judged by its own calls, as calls.py does);
our EM calls on the same reads are judged the same way.
"""
import argparse
import bisect
import collections
import itertools
import os
import random
import sys

import pysam
from Bio import Align

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import diaglib as D  # noqa: E402
import score  # noqa: E402
from islands import repeat_run  # noqa: E402

FLANK = 50
MAX_BUBBLES = 10

EDIT = Align.PairwiseAligner()
EDIT.mode = "global"
EDIT.match_score, EDIT.mismatch_score = 0, -1
EDIT.open_gap_score = EDIT.extend_gap_score = -1


def default_vcf():
    root = os.environ.get("PGDEV_EVAL_ROOT", "")
    return os.path.join(os.path.dirname(os.path.dirname(root)), "shared_calls", os.environ["PGDEV_CHROM"],
                        "deepvariant.vcf.gz")


def qpos(read, ref_pos):
    """Query position aligned at ref_pos (0-based), or the next aligned base; None outside."""
    rp, qp = read.reference_start, 0
    for op, ln in read.cigartuples:
        if op in (0, 7, 8):
            if rp <= ref_pos < rp + ln:
                return qp + (ref_pos - rp)
            rp += ln
            qp += ln
        elif op in (1, 4):
            qp += ln
        elif op in (2, 3):
            if rp <= ref_pos < rp + ln:
                return qp
            rp += ln
    return None


def load_calls(path, contig):
    """Sorted (p0, end0, allele0, allele1, het) of HiPhase-eligible calls, 0-based."""
    out = []
    for pos, recs in D.vcf_records(path, contig, 1, 2 ** 31 - 1).items():
        for ref, alts, gt, ps, c in recs:
            if len(gt) != 2 or None in gt or gt == (0, 0):
                continue
            al = [ref] + list(alts)
            i0, i1 = sorted(gt)
            het = i0 != i1
            if not het:
                i0 = 0
            a0, a1 = al[i0].upper(), al[i1].upper()
            if set(a0 + a1) - set("ACGT"):
                continue
            out.append((pos - 1, pos - 1 + len(ref), a0, a1, het))
    out.sort()
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("diag")
    ap.add_argument("--vcf", default="")
    ap.add_argument("--max-sites", type=int, default=200)
    ap.add_argument("--min-run", type=int, default=8)
    ap.add_argument("--seed", type=int, default=1)
    a = ap.parse_args()
    contig = os.environ["PGDEV_CONTIG"]
    fa = pysam.FastaFile(os.environ["PGDEV_REF"])
    ref_name = contig if contig in fa.references else contig.split("#")[-1]
    bam = pysam.AlignmentFile(os.environ["PGDEV_BAM"])
    scope = score.read_scope()
    cen = tuple(int(x) for x in os.environ.get("PGDEV_CEN", "0-0").split("-"))
    vcf = a.vcf or default_vcf()
    calls = load_calls(vcf, contig)
    starts = [c[0] for c in calls]
    # our EM: learned error per site and every scored read's side
    em_err, em_calls = {}, collections.defaultdict(dict)
    for f in D.chunk_files(os.path.join(a.diag, "em"), "em"):
        with open(f) as fh:
            for line in fh:
                p = line.rstrip("\n").split("\t")
                first = 6 if len(p) > 5 and ":" not in p[5] else 5
                for t in p[first:]:
                    x = t.split(":")
                    if x[1] in ("L", "U"):
                        continue
                    pos = int(x[0])
                    em_err[pos] = float(x[3])
                    if x[2] in ("0", "1") and p[0] in scope and p[0] not in em_calls[pos]:
                        em_calls[pos][p[0]] = int(x[2])
    em_pos = sorted(em_err)
    truth = sorted(p for p, recs in D.vcf_records(os.environ["PGDEV_TRUTH_VCF"], contig, 1, 2 ** 31 - 1).items()
                   if any(len(gt) == 2 and None not in gt and gt[0] != gt[1] for _, _, gt, _, _ in recs))
    groups = collections.defaultdict(list)
    for idx, (p0, end0, a0, a1, het) in enumerate(calls):
        if not het or (len(a0) == 1 and len(a1) == 1) or cen[0] <= p0 < cen[1]:
            continue
        run = repeat_run(fa, ref_name, p0 + 1)
        if run < a.min_run:
            continue
        tol = max(3, run)
        i = bisect.bisect_left(truth, p0 + 1 - tol)
        is_truth = i < len(truth) and truth[i] <= p0 + 1 + tol
        j = bisect.bisect_left(em_pos, p0 + 1 - tol)
        near = [em_pos[k] for k in range(j, len(em_pos)) if em_pos[k] <= p0 + 1 + tol]
        ep = max(near, key=lambda q: len(em_calls.get(q, {}))) if near else None
        state = "no EM site" if ep is None else "EM rejected" if em_err[ep] >= 0.3 else "EM reliable"
        groups[("truth het" if is_truth else "not truth het", state)].append((idx, run, ep))
    graph, inj, aln = D.load_graph(a.diag), D.load_inj(a.diag), D.load_aln(a.diag)

    def stage(pos, tol):
        """Our furthest stage at a site without an EM site (as readaudit.py)."""
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
    stages = collections.Counter()
    msa_fate = collections.Counter()
    verified_lost = collections.Counter()
    verified_examples = []
    noisy = D.load_noisy(a.diag, 1, 2 ** 31 - 1)
    noisy_beg = [r[0] for r in noisy]
    msa = D.load_msa(a.diag)
    msa_pos = sorted(msa)
    for (truth_g, state), lst in groups.items():
        if truth_g == "truth het" and state == "no EM site":
            for idx, run, ep in lst:
                pos, tol = calls[idx][0] + 1, max(3, run)
                st = stage(pos, tol)
                stages[st] += 1
                # the alignment solve's MSA at the site
                i = bisect.bisect_right(noisy_beg, pos + tol) - 1
                reg = None
                while i >= 0 and noisy[i][1] >= pos - tol - 5000:
                    if noisy[i][0] <= pos + tol and noisy[i][1] >= pos - tol:
                        reg = noisy[i]
                        break
                    i -= 1
                if reg is None:
                    fate = "no noisy region"
                else:
                    j = bisect.bisect_left(msa_pos, reg[0])
                    in_reg = [msa_pos[k] for k in range(j, len(msa_pos)) if msa_pos[k] <= reg[1]]
                    at = [p for p in in_reg if abs(p - pos) <= tol]
                    rlen = reg[1] - reg[0] + 1
                    rb = "<200" if rlen < 200 else "<500" if rlen < 500 else "<1000" if rlen < 1000 else ">=1000"
                    if at:
                        r = msa[at[0]][0]
                        cl = "1 cluster" if min(r["clusters"]) == 0 else "2 clusters"
                        fate = f"MSA call at site ({cl}, {r['cat']})"
                    elif in_reg:
                        r = msa[in_reg[0]][0]
                        cl = "1 cluster" if min(r["clusters"]) == 0 else "2 clusters"
                        fate = f"MSA ran, no call at site ({cl})"
                    else:
                        fate = f"noisy region {rb} bp, no MSA call in it"
                msa_fate[(st.split(":")[0] if st.startswith("inj_category") else st.split("_")[0], fate)] += 1
                if fate.startswith("MSA call at site (2 clusters, NOISY_CAND_HET"):
                    g = D.near(graph, pos, tol)
                    gi = ("graph " + ",".join(sorted({r["cat"] + ("/vote" if r["voting"] else "") for r in g[0][1]}))) if g else "no graph row"
                    ij = D.near(inj, pos, tol)
                    ii = ("inj " + ",".join(sorted({r["decision"] + (":" + r["cat"] if r["decision"] == "category" else "") for r in ij[0][1]}))) if ij else "no inj row"
                    al = D.near(aln, pos, tol)
                    ai = ("aln " + ",".join(sorted({f"{r['init']}>{r['final']}" for r in al[0][1]}))) if al else "no aln row"
                    verified_lost[(gi, ii, ai)] += 1
                    if len(verified_examples) < 6:
                        verified_examples.append(pos)
    rng = random.Random(a.seed)
    tab = collections.defaultdict(collections.Counter)
    skipped = collections.Counter()
    for g, lst in sorted(groups.items()):
        for idx, run, ep in rng.sample(lst, min(a.max_sites, len(lst))):
            p0, end0, a0, a1, het = calls[idx]
            w0 = max(0, p0 - run - FLANK)
            w1 = end0 + run + FLANK
            lo, hi = bisect.bisect_left(starts, w0), bisect.bisect_left(starts, w1)
            bubbles = [k for k in range(lo, hi) if calls[k][0] >= w0 and calls[k][1] <= w1]
            if len(bubbles) > MAX_BUBBLES:
                skipped["too many calls in window"] += 1
                continue
            refseq = fa.fetch(ref_name, w0, w1).upper()
            paths = []  # (sequence, target allele or None)
            for choice in itertools.product((0, 1), repeat=len(bubbles)):
                used = [(calls[k], ch, k) for k, ch in zip(bubbles, choice) if ch == 1 or calls[k][2] != refseq[calls[k][0] - w0:calls[k][1] - w0]]
                spans = sorted((c[0], c[1]) for c, ch, k in used)
                if any(spans[i][1] > spans[i + 1][0] for i in range(len(spans) - 1)):
                    continue
                seq, cur, target = [], w0, None
                for c, ch, k in sorted(used, key=lambda u: u[0][0]):
                    seq.append(refseq[cur - w0:c[0] - w0])
                    seq.append(c[3] if ch == 1 else c[2])
                    cur = c[1]
                    if k == idx:
                        target = ch
                seq.append(refseq[cur - w0:])
                if target is None:
                    choice_of = dict(zip(bubbles, choice))
                    target = choice_of[idx]
                paths.append(("".join(seq), target))
            uniq = {}
            for s, t in paths:
                uniq.setdefault(s, set()).add(t)
            reads = []
            for r in bam.fetch(contig, w0, w1):
                if r.is_secondary or r.is_supplementary or r.query_name not in scope:
                    continue
                if r.reference_start > w0 or r.reference_end < w1:
                    continue
                q0, q1 = qpos(r, w0), qpos(r, w1 - 1)
                if q0 is None or q1 is None or q1 < q0:
                    continue
                seg = r.query_sequence[q0:q1 + 1].upper()
                best, alle = None, set()
                for s, ts in uniq.items():
                    d = -EDIT.score(seg, s)
                    if best is None or d < best:
                        best, alle = d, set(ts)
                    elif d == best:
                        alle |= ts
                call = next(iter(alle)) if len(alle) == 1 else None
                reads.append((r.query_name, call))
            # judge HiPhase calls by the parents' majority of HiPhase calls
            par = [0, 0, 0, 0]
            for q, call in reads:
                if call is not None:
                    par[(0 if scope[q][2] else 2) + call] += 1
            if not (par[0] + par[1] and par[2] + par[3]):
                skipped["no split among HiPhase calls"] += 1
                continue
            ms, ps = (0 if par[0] >= par[1] else 1), (0 if par[2] >= par[3] else 1)
            t = tab[g]
            t["sites"] += 1
            if ms == ps:
                t["sites without a parental split"] += 1
            for q, call in reads:
                t["reads"] += 1
                if call is None:
                    t["hiphase ambiguous"] += 1
                elif ms != ps:
                    t["hiphase right" if call == (ms if scope[q][2] else ps) else "hiphase wrong"] += 1
            if ep is not None:
                epar = [0, 0, 0, 0]
                for q, side in em_calls[ep].items():
                    epar[(0 if scope[q][2] else 2) + side] += 1
                ems, eps = (0 if epar[0] >= epar[1] else 1), (0 if epar[2] >= epar[3] else 1)
                for q, call in reads:
                    if q not in em_calls[ep]:
                        t["EM no call"] += 1
                    elif ems != eps:
                        t["EM right" if em_calls[ep][q] == (ems if scope[q][2] else eps) else "EM wrong"] += 1
    print(f"{contig}: HiPhase-exact allele assignment at {os.path.basename(vcf)} het indels in repeats >= {a.min_run} bp "
          f"(up to {a.max_sites} sites per group)" + (f"; skipped {dict(skipped)}" if skipped else ""))
    print("  sites per group (all, before sampling): " + ", ".join(f"{g[0]}/{g[1]} {len(v):,}" for g, v in sorted(groups.items())))
    print(f"  {'group':36} {'sites':>6} {'reads':>7}   {'HiPhase: called':>16} {'wrong %':>8}   {'EM: called':>11} {'wrong %':>8}")
    for g in sorted(tab):
        t = tab[g]
        hc = t["hiphase right"] + t["hiphase wrong"]
        ec = t["EM right"] + t["EM wrong"]
        n = max(1, t["reads"])
        print(f"  {g[0] + ', ' + g[1]:36} {t['sites']:6,} {t['reads']:7,}   "
              f"{hc:7,} ({100 * (hc + t['hiphase ambiguous']) / n if False else 100 * hc / n:4.0f}%) "
              f"{100 * t['hiphase wrong'] / max(1, hc):7.1f}%   "
              f"{ec:6,} ({100 * ec / n:3.0f}%) {100 * t['EM wrong'] / max(1, ec):7.1f}%")
    print("  truth-het calls without an EM site, our furthest stage: " +
          ", ".join(f"{k} {v}" for k, v in stages.most_common()))
    print("  ... and the alignment solve's MSA there (stage group, MSA fate):")
    for (sg, fate), v in sorted(msa_fate.items(), key=lambda kv: -kv[1]):
        print(f"    {sg:14} {fate:52} {v:5,}")
    print("  MSA-verified hets (2 clusters, NOISY_CAND_HET) that never reach the EM: graph row, injection, solve")
    for k, v in verified_lost.most_common(12):
        print(f"    {v:5,}  {k[0]}  |  {k[1]}  |  {k[2]}")
    print("    examples: " + ", ".join(f"{p:,}" for p in verified_examples))
    print(f"  (HiPhase ambiguous reads per group: " +
          ", ".join(f"{g[0][:3]}/{g[1]}: {tab[g]['hiphase ambiguous']}" for g in sorted(tab)) + ")")


if __name__ == "__main__":
    main()
