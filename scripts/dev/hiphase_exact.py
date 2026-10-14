#!/usr/bin/env python3
"""HiPhase's own per-read allele assignment against ours, exactly (evaluation only).

Usage: hiphase_exact.py DIAG_DIR [--sites N] [--min-run BP] [--workers W] [--seed S]

The allele assignment is hiphase_port.py, a line-for-line port of HiPhase's
default (global realignment) path, checked against brute force. It runs on
HiPhase's own inputs, read from the command line recorded in the rival's
resources.txt (--reference, --bam, --vcf): every het and hom-alt call in the
VCF, every read passing HiPhase's filter (MAPQ >= 5; not unmapped, secondary,
QC-fail or duplicate), the read's whole aligned span aligned to the graph of
all calls in it, falling back to local realignment when the edit distance
exceeds min(500, ceil(10% of the span)).

Sites: het indels in the VCF with a reference repeat run >= --min-run bp,
grouped by truth (a truth het within max(3, run) bp) and our EM state there,
--sites per group. A call is right when it is the allele the read's parent's
reads mostly carry (HiPhase's calls and our EM's each judged by their own
parental majority); reads are matched to our EM by name.
"""
import argparse
import bisect
import collections
import multiprocessing
import os
import random
import re
import shlex
import sys

import pysam

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import diaglib as D  # noqa: E402
import hiphase_port as H  # noqa: E402
import score  # noqa: E402
from islands import repeat_run  # noqa: E402

G = {}  # shared with forked workers


def hiphase_inputs():
    rdir = D.competitors()["hiphase"]
    with open(os.path.join(rdir, "resources.txt")) as fh:
        line = next(l for l in fh if "Command being timed" in l)
    cmd = shlex.split(re.search(r'"(.*)"', line).group(1))
    get = lambda flag: cmd[cmd.index(flag) + 1]  # noqa: E731
    return get("--reference"), get("--bam"), get("--vcf")


def work(site_i):
    hets, homs, seq, contig, scope = G["hets"], G["homs"], G["seq"], G["contig"], G["scope"]
    v = hets[site_i]
    bam = pysam.AlignmentFile(G["bam"])
    out, fallbacks = [], 0
    starts = G["starts"]
    for r in bam.fetch(contig, v.position, v.position + v.ref_len):
        if H.filter_out(r) or r.is_supplementary or r.query_name not in scope:
            continue
        try:
            alleles, _ = H.global_realignment(r, hets, homs, seq, 0)
            a = alleles[site_i]
        except H.MaxEditDistance:
            fallbacks += 1
            lo = bisect.bisect_left(starts, r.reference_start - 20)
            hi = bisect.bisect_right(starts, r.reference_end + 20)
            sub = hets[lo:hi]
            loc = H.local_realignment(r, sub)
            a = loc[site_i - lo] if lo <= site_i < hi else H.NOOVERLAP
        out.append((r.query_name, a))
    return site_i, out, fallbacks


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("diag")
    ap.add_argument("--sites", type=int, default=30)
    ap.add_argument("--min-run", type=int, default=8)
    ap.add_argument("--workers", type=int, default=3)
    ap.add_argument("--seed", type=int, default=1)
    a = ap.parse_args()
    ref_path, bam_path, vcf_path = hiphase_inputs()
    fa = pysam.FastaFile(ref_path)
    contig = os.environ["PGDEV_CONTIG"]
    hp_contig = contig if contig in fa.references else contig.split("#")[-1]
    seq = fa.fetch(hp_contig).upper()
    our_fa = pysam.FastaFile(os.environ["PGDEV_REF"])
    our_contig = contig if contig in our_fa.references else contig.split("#")[-1]
    scope = score.read_scope()
    cen = tuple(int(x) for x in os.environ.get("PGDEV_CEN", "0-0").split("-"))
    records = []
    for pos, recs in sorted(D.vcf_records(vcf_path, contig, 1, 2 ** 31 - 1).items()):
        for ref, alts, gt, ps, c in recs:
            records.append((pos - 1, ref.upper(), [x.upper() for x in alts], gt))
    hets, homs = H.load_variants(records, lambda x, y: seq[x:y], reference_buffer=15)
    starts = [v.position for v in hets]
    # our EM
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
    for i, v in enumerate(hets):
        if v.vtype == "snv" or cen[0] <= v.position < cen[1]:
            continue
        run = repeat_run(our_fa, our_contig, v.position + 1)
        if run < a.min_run:
            continue
        p1, tol = v.position + 1, max(3, run)
        j = bisect.bisect_left(truth, p1 - tol)
        is_truth = j < len(truth) and truth[j] <= p1 + tol
        k = bisect.bisect_left(em_pos, p1 - tol)
        near = [em_pos[m] for m in range(k, len(em_pos)) if em_pos[m] <= p1 + tol]
        ep = max(near, key=lambda q: len(em_calls.get(q, {}))) if near else None
        state = "no EM site" if ep is None else "EM rejected" if em_err[ep] >= 0.3 else "EM reliable"
        groups[("truth het" if is_truth else "not truth het", state)].append((i, ep))
    rng = random.Random(a.seed)
    chosen = {g: rng.sample(lst, min(a.sites, len(lst))) for g, lst in groups.items()}
    G.update(hets=hets, homs=homs, seq=seq, contig=hp_contig, scope=scope, bam=bam_path, starts=starts)
    todo = sorted({i for lst in chosen.values() for i, ep in lst})
    with multiprocessing.get_context("fork").Pool(a.workers) as pool:
        results = {i: (out, fb) for i, out, fb in pool.imap_unordered(work, todo)}
    print(f"{contig}: HiPhase inputs {os.path.basename(bam_path)}, {os.path.basename(vcf_path)}; "
          f"sites per group (all): " + ", ".join(f"{g[0]}/{g[1]} {len(v):,}" for g, v in sorted(groups.items())))
    print(f"  {'group':36} {'sites':>6} {'reads':>7}   {'HiPhase: called':>16} {'wrong %':>8} {'ambig':>6} {'fallback':>8}"
          f"   {'EM: called':>11} {'wrong %':>8}")
    for g in sorted(chosen):
        t = collections.Counter()
        for i, ep in chosen[g]:
            out, fb = results[i]
            t["fallback"] += fb
            calls = [(q, al) for q, al in out if al in (H.REF, H.ALT)]
            par = [0, 0, 0, 0]
            for q, al in calls:
                par[(0 if scope[q][2] else 2) + al] += 1
            t["sites"] += 1
            t["reads"] += len(out)
            t["ambig"] += sum(1 for q, al in out if al == H.AMBIG)
            if par[0] + par[1] and par[2] + par[3]:
                ms, ps = (0 if par[0] >= par[1] else 1), (0 if par[2] >= par[3] else 1)
                if ms != ps:
                    for q, al in calls:
                        t["hip right" if al == (ms if scope[q][2] else ps) else "hip wrong"] += 1
            if ep is not None:
                epar = [0, 0, 0, 0]
                for q, side in em_calls[ep].items():
                    epar[(0 if scope[q][2] else 2) + side] += 1
                ems, eps = (0 if epar[0] >= epar[1] else 1), (0 if epar[2] >= epar[3] else 1)
                if ems != eps:
                    for q, al in out:
                        if q in em_calls[ep]:
                            t["em right" if em_calls[ep][q] == (ems if scope[q][2] else eps) else "em wrong"] += 1
        hc, ec, n = t["hip right"] + t["hip wrong"], t["em right"] + t["em wrong"], max(1, t["reads"])
        print(f"  {g[0] + ', ' + g[1]:36} {t['sites']:6,} {t['reads']:7,}   {hc:7,} ({100 * hc / n:4.0f}%) "
              f"{100 * t['hip wrong'] / max(1, hc):7.1f}% {t['ambig']:6,} {t['fallback']:8,}   "
              f"{ec:6,} ({100 * ec / n:3.0f}%) {100 * t['em wrong'] / max(1, ec):7.1f}%")


if __name__ == "__main__":
    main()
