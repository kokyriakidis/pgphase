#!/usr/bin/env python3
"""Real hets we have no EM site for: where each one is lost (evaluation only).

Usage: missing.py DIAG_DIR [--examples N]

Sites: truth het indels on the arms that the rival's calls (the DeepVariant
VCF HiPhase phased) also call het, with no EM site within max(3, repeat run)
bp. For each, at the nearest position within that distance: the MSA's
verdict (msa dump: category, clusters), the alignment candidates there (aln
dump: initial > final category), their injection decision (inj dump:
decision, category, MSA-verified), the graph row (graph dump: category,
voting) and the EM's view of injected rows (emsites: in EM, or left out).
The chain is summarised as one label per site and counted; examples print
the rows for debugging.
"""
import argparse
import bisect
import collections
import os
import sys

import pysam

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import diaglib as D  # noqa: E402
from islands import repeat_run  # noqa: E402
from hiphase_exact import hiphase_inputs  # noqa: E402


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("diag")
    ap.add_argument("--examples", type=int, default=3)
    a = ap.parse_args()
    contig = os.environ["PGDEV_CONTIG"]
    fa = pysam.FastaFile(os.environ["PGDEV_REF"])
    ref_name = contig if contig in fa.references else contig.split("#")[-1]
    cen = tuple(int(x) for x in os.environ.get("PGDEV_CEN", "0-0").split("-"))
    truth = sorted(p for p, recs in D.vcf_records(os.environ["PGDEV_TRUTH_VCF"], contig, 1, 2 ** 31 - 1).items()
                   if any(len(gt) == 2 and None not in gt and gt[0] != gt[1] for _, _, gt, _, _ in recs))
    # the rival's het calls (the VCF HiPhase phased): (p0, end0, allele0, allele1, het)
    dv = []
    for pos, recs in sorted(D.vcf_records(hiphase_inputs()[2], contig, 1, 2 ** 31 - 1).items()):
        for ref, alts, gt, ps, c in recs:
            if len(gt) != 2 or None in gt or gt[0] == gt[1]:
                continue
            al = [ref] + list(alts)
            i0, i1 = sorted(gt)
            dv.append((pos - 1, pos - 1 + len(ref), al[i0].upper(), al[i1].upper(), True))
    ems = D.load_emsites(a.diag)
    em_in = sorted(p for p, rows in ems.items() if any(r["in_em"] for r in rows))
    graph, inj, aln, msa = D.load_graph(a.diag), D.load_inj(a.diag), D.load_aln(a.diag), D.load_msa(a.diag)
    msa_pos = sorted(msa)

    def near_rows(index, pos, tol):
        """every row within tol of pos (all positions, not only the nearest)"""
        hit = D.near(index, pos, tol)
        return [r for p, rows in hit for r in rows] if hit else []

    tab = collections.Counter()
    examples = collections.defaultdict(list)
    for p0, end0, a0, a1, het in dv:
        pos = p0 + 1
        if not het or (len(a0) == 1 and len(a1) == 1) or cen[0] <= p0 < cen[1]:
            continue
        run = repeat_run(fa, ref_name, pos)
        tol = max(3, run)
        j = bisect.bisect_left(truth, pos - tol)
        if not (j < len(truth) and truth[j] <= pos + tol):
            continue
        k = bisect.bisect_left(em_in, pos - tol)
        if k < len(em_in) and em_in[k] <= pos + tol:
            continue
        i = bisect.bisect_left(msa_pos, pos - tol)
        mrows = [r for m in msa_pos[i:i + 10] if m <= pos + tol for r in msa[m]]
        mtag = ("MSA " + "/".join(sorted({r["cat"] for r in mrows}))) if mrows else "no MSA call"
        irows = near_rows(inj, pos, tol)
        if irows:
            best = min(irows, key=lambda r: D.decision_rank(r["decision"]))
            itag = "inj best: " + best["decision"] + (":" + best["cat"] if best["decision"] == "category" else "") + \
                   (" (MSA-verified)" if best.get("msa") else "")
        else:
            itag = "no inj row"
        arows = near_rows(aln, pos, tol)
        atag = ("aln " + "/".join(sorted({f"{r['init']}>{r['final']}" for r in arows}))) if arows else "no aln row"
        grows = near_rows(graph, pos, tol)
        gtag = ("graph " + "/".join(sorted({r["cat"] + ("+vote" if r["voting"] else "") for r in grows}))) if grows else "no graph row"
        erows = near_rows(ems, pos, tol)
        etag = "emsites: left out" if erows else ""
        label = f"{mtag} | {itag}"
        tab[label] += 1
        if len(examples[label]) < a.examples:
            examples[label].append((pos, run, a0, a1, atag, gtag, etag))
    total = sum(tab.values())
    print(f"{contig}: truth-het indels called het by the rival with no EM site: {total:,}")
    for label, n in tab.most_common(14):
        print(f"  {n:5,}  {label}")
        for e in examples[label]:
            print(f"           e.g. {e[0]:,} run {e[1]} {e[2]}>{e[3]}  {e[4]}  {e[5]}  {e[6]}")


if __name__ == "__main__":
    main()
