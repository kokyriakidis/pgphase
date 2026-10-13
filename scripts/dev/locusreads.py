#!/usr/bin/env python3
"""Every read at a locus against every site, in the alignment solve and the EM (evaluation only).

Usage: locusreads.py DIAG_DIR REGION [--rival NAME] [--missed] [--max-sites N]

Sites (columns): the truth heterozygotes in REGION, the rival's phased
heterozygous records there, and our EM sites (junk ones included: a column
with no truth shows the raw allele, 0 / 1 / 2). The legend lists
each column: position, type, truth GT, our furthest stage (EM kind / learned
error / block, else the injection decision, else the alignment category),
and the rival's GT and phase set.

Reads (rows): scored reads overlapping REGION, by start; --missed keeps the
reads a rival phases correctly and we do not. Each row: read, parent (M/P),
MAPQ, our class (C correct, D discordant, U unphased) with HP and phase-set
letter, the rivals' classes, the alignment solve's HP and phase-set letter;
then two grids over the columns:
  aln   the read's allele in the alignment solve (before injection)
  em    the read's allele in the final EM
with = the allele its truth haplotype carries, x the other allele, . the read
spans the site but has no call, blank the read does not span it. An indel
whose representation differs from truth's may show x on every read.
"""
import argparse
import collections
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import diaglib as D  # noqa: E402
import score  # noqa: E402



def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("diag")
    ap.add_argument("region")
    ap.add_argument("--rival", default="hiphase")
    ap.add_argument("--missed", action="store_true")
    ap.add_argument("--max-sites", type=int, default=90)
    a = ap.parse_args()
    contig, beg, end = D.parse_region(a.region)
    comps = D.competitors()
    scope = score.read_scope()
    tools = {"ours": a.diag, **comps}
    cls = {n: score.read_classes(d) for n, d in tools.items()}
    reads = sorted((q for q, v in scope.items() if v[3] < end and v[4] > beg), key=lambda q: scope[q][3])
    if a.missed:
        reads = [q for q in reads if cls["ours"][q] != "correct" and any(cls[n][q] == "correct" for n in comps)]
    names = set(reads)
    tags = {n: score.read_tags(d, names) for n, d in tools.items()}

    # columns: truth hets (biallelic, phased), plus the rival's non-truth sites
    tp = score.truth_tables()["truth_phased"]
    truth = {}
    for pos, recs in D.vcf_records(os.environ["PGDEV_TRUTH_VCF"], contig, beg, end).items():
        for ref, alts, gt, ps, c in recs:
            key = (pos, ref.upper(), alts[0].upper()) if len(alts) == 1 else None
            if key in tp:
                truth[pos] = (ref, alts[0], tp[key])
    rival = {}
    if a.rival in comps:
        for pos, recs in D.vcf_records(os.path.join(comps[a.rival], "phased.vcf"), contig, beg, end).items():
            for ref, alts, gt, ps, c in recs:
                if ps and len(gt) == 2 and None not in gt and gt[0] != gt[1]:
                    rival[pos] = (ref, alts, gt, ps)
    em_reads, em_sites = D.load_em(a.diag, beg - 100_000, end + 100_000)
    ours_sites = {p for p in em_sites if beg <= p <= end}
    cols = sorted(set(truth) | set(rival) | ours_sites)[: a.max_sites]

    ems = D.load_emsites(a.diag, beg, end)
    inj, aln = D.load_inj(a.diag, beg, end), D.load_aln(a.diag, beg, end)
    alnr = D.load_alnreads(a.diag, beg, end, names)

    def stage(pos):
        e = em_sites.get(pos)
        if e:
            return f"em {e['kind']} {e['err']:.2f} b{e['block']}"
        i = inj.get(pos)
        if i:
            r = min(i, key=lambda r: D.decision_rank(r["decision"]))
            return f"inj {r['decision']} {r['cat']}{' msa' if r['msa'] else ''} conf {r['conf']}"
        al = aln.get(pos)
        if al:
            r = max(al, key=lambda r: r["alt_cov"])
            return f"aln {r['init']}>{r['final']} {r['ref']}/{r['alt_cov']}"
        return "absent"

    # which truth haplotype gt[0] belongs to: the parent whose reads carry it in the alignment solve
    agree = collections.Counter()
    for q in reads:
        for pos, (typ, al) in alnr.get(q, {"calls": {}})["calls"].items():
            if pos in truth and typ == "S" and al in (0, 1):
                agree[(al == truth[pos][2]) == scope[q][2]] += 1
    first_is_mat = agree[True] >= agree[False]

    def expected(q, pos):
        g0 = truth[pos][2]
        return g0 if scope[q][2] == first_is_mat else 1 - g0

    def em_allele(q, pos):
        d = em_reads.get(q)
        if not d:
            return None
        rows = ems.get(pos)
        for o in d["obs"]:
            if o[0] != pos:
                continue
            if rows and rows[0]["in_em"]:
                return rows[0]["x"] if o[2] == 0 else rows[0]["y"]
            return o[2]  # a window: its side (0 / 1)
        return None

    def cell(q, pos, allele):
        if not (scope[q][3] < pos <= scope[q][4]):
            return " "
        if allele is None:
            return "."
        if pos not in truth:
            return str(allele)[:1]
        return "=" if allele == expected(q, pos) else "x"

    letters = {n: {} for n in list(tools) + ["aln"]}

    def letter(n, ps):
        if ps is None:
            return "-"
        d = letters[n]
        if ps not in d:
            d[ps] = chr(ord("A") + len(d) % 26)
        return d[ps]

    print(f"region {contig}:{beg:,}-{end:,}  reads {len(reads)}  columns {len(cols)}  "
          f"(truth GT first haplotype = {'maternal' if first_is_mat else 'paternal'}, from {sum(agree.values())} calls)")
    print("\ncolumns:")
    for i, pos in enumerate(cols):
        t = truth.get(pos)
        r = rival.get(pos)
        tcell = f"truth {t[0][:6]}>{t[1][:6]} " if t else "not truth      "
        rcell = f"{a.rival} {'|'.join(map(str, r[2]))} PS {r[3]}" if r else f"{a.rival} -"
        print(f"  {i:>3} {pos:>11,}  {tcell:22} {stage(pos):44} {rcell}")
    ruler = "".join(str(i // 10 % 10) for i in range(len(cols))) + "\n" + "".join(str(i % 10) for i in range(len(cols)))
    head = f"{'read':12} {'p':1} {'mq':>3} {'ours':8} " + " ".join(f"{n[:4]:4}" for n in comps) + f" {'aln':6} "
    pad = " " * len(head)
    print("\n" + pad + ruler.split("\n")[0] + "   " + ruler.split("\n")[0])
    print(head + ruler.split("\n")[1] + "   " + ruler.split("\n")[1] + "\n" + pad + "aln" + " " * max(0, len(cols) - 3) + "   em")
    for q in reads:
        t = tags["ours"].get(q)
        ours = f"{cls['ours'][q][0].upper()}{t[0] if t else ' '}{letter('ours', t[1] if t else None)}"
        riv = " ".join(f"{cls[n][q][0].upper():4}" for n in comps)
        ar = alnr.get(q)
        acell = f"{ar['hap']}{letter('aln', ar['ps']) if ar and ar['hap'] else '-'}" if ar else "--"
        mq = ar["mapq"] if ar else (em_reads.get(q) or {}).get("mapq", -1)
        agrid = "".join(cell(q, pos, ar["calls"].get(pos, (None, None))[1] if ar else None) for pos in cols)
        egrid = "".join(cell(q, pos, em_allele(q, pos)) for pos in cols)
        print(f"{q.split('/')[1][-12:]:12} {'M' if scope[q][2] else 'P'} {mq:>3} {ours:8} {riv} {acell:6} {agrid}   {egrid}")


if __name__ == "__main__":
    main()
