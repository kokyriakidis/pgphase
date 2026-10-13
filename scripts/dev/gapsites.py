#!/usr/bin/env python3
"""Sites a rival uses to close a gap we leave, and how far each gets with us (evaluation only).

Usage: gapsites.py DIAG_DIR [--rival NAME] [--kind clean|all] [--top N] [--out TSV]

Gaps: arm reads crossing a truth het that we leave unphased and the rival
phases correctly, merged into loci (reads overlapping, 5 kb slack). For each
locus, the rival's phase set holding most of those reads, and that phase
set's heterozygous records inside the locus: the sites the rival closes the
gap with. Each site gets
  kind     snp | indel (non-repeat) | str (period 2-6 tract) | hp (homopolymer
           tract) | complex; tract length in bp for str / hp
  truth    het (a truth heterozygote at the position; indels within 3 bp,
           repeat indels within their tract length) | hom | no
  stage    our furthest stage (funnel.py's names) with its detail: EM kind,
           learned error, block; graph row category; injection decision,
           category, MSA, depth / chunk median, confident reads per allele;
           alignment category initial>final and allele counts; for a site
           short of a reliable EM site, its MSA record: read calls ref,alt,
           missing before > after the refresh (LOST when it dropped some) and
           each consensus's longest stretch without a base (a deletion)
  obs      of the locus reads, how many we have a call at this site for
--kind clean (default) keeps loci where a gap read covers a snp or
non-repeat indel the rival uses, ranked by those reads. A summary counts the rival's sites in those loci by kind
and stage, weighted by the locus reads, so the largest losses come first.
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



def period(s):
    for p in range(1, len(s) + 1):
        if len(s) % p == 0 and s[:p] * (len(s) // p) == s:
            return p
    return len(s)


def site_kind(fa, ref_name, pos, ref, alt):
    """snp | indel | str | hp | complex, and the repeat tract length (0 when not in a repeat)."""
    if len(ref) == 1 and len(alt) == 1:
        return "snp", 0
    if len(ref) > 1 and len(alt) > 1 or ref[0].upper() != alt[0].upper():
        return "complex", 0
    s = (alt[1:] if len(alt) > len(ref) else ref[1:]).upper()
    p = period(s)
    if p > 6:
        return "indel", 0
    u = s[:p]
    w0 = max(0, pos - 200)
    seq = fa.fetch(ref_name, w0, pos + 200 + len(ref)).upper()
    i = pos - w0  # 0-based index just after the VCF anchor base
    lo = i
    hi = i if len(alt) > len(ref) else i + len(s)
    k = 0
    while lo - 1 >= 0 and seq[lo - 1] == u[p - 1 - (k % p)]:
        lo -= 1; k += 1
    k = 0
    while hi < len(seq) and seq[hi] == u[k % p]:
        hi += 1; k += 1
    tract = hi - lo
    # a repeat needs the unit at least twice in the reference, beyond the indel itself
    copies = tract // p - (0 if len(alt) > len(ref) else len(s) // p)
    if copies < 2 and tract < 2 * p:
        return "indel", 0
    return ("hp" if p == 1 else "str"), tract


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("diag")
    ap.add_argument("--rival", default="hiphase")
    ap.add_argument("--kind", default="clean", choices=["clean", "all"])
    ap.add_argument("--top", type=int, default=15)
    ap.add_argument("--out")
    a = ap.parse_args()
    contig = os.environ["PGDEV_CONTIG"]
    fa = pysam.FastaFile(os.environ["PGDEV_REF"])
    ref_name = contig if contig in fa.references else contig.split("#")[-1]
    comps = D.competitors()
    rdir = comps[a.rival]
    scope = {q: v for q, v in score.read_scope().items() if v[1] == "arms_phaseable"}
    ours = score.read_classes(a.diag)
    theirs = score.read_classes(rdir)
    target = sorted((q for q in scope if ours[q] == "unphased" and theirs[q] == "correct"), key=lambda q: scope[q][3])
    loci, cur = [], None
    for q in target:
        s, e = scope[q][3], scope[q][4]
        if cur and s <= cur["end"] + 5000:
            cur["end"] = max(cur["end"], e); cur["reads"].append(q)
        else:
            cur = {"beg": s, "end": e, "reads": [q]}; loci.append(cur)

    tags = score.read_tags(rdir, set(target))
    rvcf = D.vcf_records(os.path.join(rdir, "phased.vcf"), contig, 1, 2 ** 31 - 1)
    by_ps = collections.defaultdict(list)
    for pos, recs in rvcf.items():
        for ref, alts, gt, ps, c in recs:
            if ps and len(gt) == 2 and None not in gt and gt[0] != gt[1]:
                by_ps[ps].append((pos, ref, alts[max(gt) - 1]))
    for ps in by_ps:
        by_ps[ps].sort()
    truth_vcf = D.vcf_records(os.environ["PGDEV_TRUTH_VCF"], contig, 1, 2 ** 31 - 1)
    truth = {}
    for pos, recs in truth_vcf.items():
        for ref, alts, gt, ps, c in recs:
            if len(gt) == 2 and None not in gt:
                truth[pos] = "het" if gt[0] != gt[1] else "hom"
    # EM sites from the small emsites dump; EM calls only for the gap reads
    em_sites = {}
    for pos, rows in D.load_emsites(a.diag).items():
        r = next((r for r in rows if r["in_em"]), None)
        if r is not None:
            em_sites[pos] = dict(kind=r["origin"], err=r["err"], block=r["block"])
    want = set(target)
    em_reads = {}
    for f in D.chunk_files(os.path.join(a.diag, "em"), "em"):
        with open(f) as fh:
            for line in fh:
                q = line[:line.index("\t")]
                if q not in want:
                    continue
                p = line.rstrip("\n").split("\t")
                first = 6 if len(p) > 5 and ":" not in p[5] else 5
                obs = [(int(t.split(":", 1)[0]),) for t in p[first:]]
                if q not in em_reads or len(obs) > len(em_reads[q]["obs"]):
                    em_reads[q] = {"obs": obs}
    graph, aln, inj = D.load_graph(a.diag), D.load_aln(a.diag), D.load_inj(a.diag)
    msa = D.load_msa(a.diag)

    def msa_note(pos, tol):
        m = D.near(msa, pos, tol)
        if not m:
            return ""
        r = m[0][1][0]
        lost = " LOST" if r["after"][0] + r["after"][1] < r["before"][0] + r["before"][1] else ""
        return (f" | msa {','.join(map(str, r['before']))}>{','.join(map(str, r['after']))}{lost}"
                f" gap {r['longest_gap'][0]}/{r['longest_gap'][1]}")

    def truth_of(pos, tol):
        hits = D.near(truth, pos, tol)
        return hits[0][1] if hits else "no"

    def stage_of(pos, tol):
        e = D.near(em_sites, pos, tol)
        if e:
            p, s = e[0]
            st = "em_unverified" if s["kind"] == "U" else "em_unreliable" if s["err"] >= 0.3 else "em_reliable"
            return st, f"{s['kind']} err {s['err']:.2f} b{s['block']}", p
        g = D.near(graph, pos, tol)
        if g:
            rows = g[0][1]
            if any(r["voting"] for r in rows):
                return "graph_voting", "covs " + rows[0]["covs"], g[0][0]
            return "graph_" + rows[0]["cat"], "covs " + rows[0]["covs"], g[0][0]
        i = D.near(inj, pos, tol)
        if i:
            r = min(i[0][1], key=lambda r: D.decision_rank(r["decision"]))
            st = "inj_" + r["decision"] + ((":" + r["cat"] + (":msa" if r["msa"] else "")) if r["decision"] == "category" else "")
            return st, f"{r['cat']}{' msa' if r['msa'] else ''} depth {r['depth']}/{r['median']:.0f} conf {r['conf']}", i[0][0]
        al = D.near(aln, pos, tol)
        if al:
            r = max(al[0][1], key=lambda r: r["alt_cov"])
            return f"aln_{r['init']}>{r['final']}", f"ref {r['ref']} alt {r['alt_cov']} of {r['total']}", al[0][0]
        return "absent", "", None

    rows = []
    summary = collections.Counter()
    summary_sites = collections.Counter()
    for l in loci:
        rs = l["reads"]
        ps = collections.Counter(tags[q][1] for q in rs if q in tags)
        if not ps:
            continue
        rps = str(ps.most_common(1)[0][0])
        lst = by_ps.get(rps, [])
        i = bisect.bisect_left(lst, (l["beg"],))
        sites = []
        for pos, ref, alt in lst[i:]:
            if pos > l["end"]:
                break
            kind, tract = site_kind(fa, ref_name, pos, ref, alt)
            # a SNP is never credited to an adjacent indel's site; a repeat indel
            # may sit anywhere in its tract (truth left-aligns it elsewhere)
            tol = 0 if kind == "snp" else max(3, tract)
            st, detail, our_pos = stage_of(pos, tol)
            if not st.startswith("em_reliable"):
                detail += msa_note(pos, tol)
            obs = 0
            if our_pos is not None and st.startswith("em_"):
                obs = sum(1 for q in rs if q in em_reads and any(abs(o[0] - our_pos) <= 1 for o in em_reads[q]["obs"]))
            n_cover = sum(1 for q in rs if scope[q][3] <= pos <= scope[q][4])
            sites.append(dict(pos=pos, ref=ref, alt=alt, kind=kind, tract=tract, truth=truth_of(pos, tol),
                              stage=st, detail=detail, obs=obs, cover=n_cover))
        clean_reads = len({q for q in rs for s in sites
                           if s["kind"] in ("snp", "indel") and scope[q][3] <= s["pos"] <= scope[q][4]})
        clean = clean_reads > 0
        if a.kind == "clean" and not clean:
            continue
        rows.append((clean_reads if a.kind == "clean" else len(rs), l, rps, sites))
        for s in sites:
            if s["kind"] in ("snp", "indel") and s["cover"]:
                summary[(s["kind"], s["truth"], s["stage"])] += s["cover"]
                summary_sites[(s["kind"], s["truth"], s["stage"])] += 1
    rows.sort(key=lambda r: -r[0])
    total = sum(r[0] for r in rows)
    print(f"{a.rival}: gaps it closes ({'with at least one snp / non-repeat indel' if a.kind == 'clean' else 'all'}): "
          f"{len(rows)} loci, {total:,} reads" + (" covering such a site" if a.kind == "clean" else ""))
    print("\nits snp / non-repeat indel sites in those gaps, by our furthest stage (sites; reads covering them):")
    for (k, t, st), n in sorted(summary.items(), key=lambda kv: -kv[1]):
        print(f"  {k:5} truth {t:3}  {st:40} sites {summary_sites[(k, t, st)]:4}  reads {n:5}")
    for n, l, rps, sites in rows[:a.top]:
        print(f"\n{l['beg']:,}-{l['end']:,} ({l['end'] - l['beg']:,} bp)  gap reads {len(l['reads'])}"
              f"{f', on a clean site {n}' if a.kind == 'clean' else ''}  {a.rival} PS {rps}")
        for s in sites:
            k = s["kind"] + (f"{s['tract']}" if s["tract"] else "")
            obs = f"obs {s['obs']}/{s['cover']}" if s["stage"].startswith("em_") else f"cover {s['cover']}"
            print(f"  {s['pos']:>11,} {s['ref'][:8]:>8}>{s['alt'][:8]:<8} {k:8} truth {s['truth']:3}  {s['stage']:34} {s['detail']:60} {obs}")
    if a.out:
        with open(a.out, "w") as fh:
            fh.write("locus_beg\tlocus_end\treads\tpos\tref\talt\tkind\ttract\ttruth\tstage\tdetail\tobs\tcover\n")
            for n, l, rps, sites in rows:
                for s in sites:
                    fh.write(f"{l['beg']}\t{l['end']}\t{n}\t{s['pos']}\t{s['ref']}\t{s['alt']}\t{s['kind']}\t{s['tract']}\t"
                             f"{s['truth']}\t{s['stage']}\t{s['detail']}\t{s['obs']}\t{s['cover']}\n")


if __name__ == "__main__":
    main()
