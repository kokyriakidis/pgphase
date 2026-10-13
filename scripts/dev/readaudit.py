#!/usr/bin/env python3
"""Every read a rival phases correctly and we do not, each in one primary cause (evaluation only).

Usage: readaudit.py DIAG_DIR RIVAL [--reads islands|FILE] [--list CLASS] [--top N]

Population: arms reads crossing a truth het that RIVAL phases correctly and we
leave unphased or discordant (official classes). For each read, the rival's
evidence is its phased heterozygous records inside the read carrying the
phase set it tagged the read with. Each such site gets
  kind   snp | indel | hp | str | complex (repeat tract from the reference)
  truth  het within 0 bp (SNP), 3 bp (indel) or the tract length (repeat)
  stage  our furthest stage at the site (as funnel.py)
and the read gets its own evidence: MAPQ, and its reliable EM calls.
Primary cause, first match wins:
  low_mapq               read MAPQ < 30
  rival_no_truth_site    none of the rival's sites in the read is a truth het
  good_site_called       a clean (snp / indel) truth-het rival site is a
                         reliable EM site and the read has a call there
  good_site_no_call      ... a reliable EM site, but the read has no call there
  good_site_lost:<stage> a clean truth-het rival site we do not phase reliably
  repeat_only            the read's truth-het rival sites are all repeat indels
Discordant reads are then split by their own reliable EM calls, each judged
by which side each parent's reads carry at that site (all scored reads; real
at parental purity >= 0.9, noise <= 0.65), so no allele representation is
compared:
  labelled by noise sites    most calls sit at noise sites
  wrong calls at real sites  more wrong than right calls at real sites
  right calls (block)        all real-site calls right: the block's
                             orientation there is what is wrong
  no clear site              no call at a site with a clear parental split
Then every rival truth-het site in those reads, by kind, our furthest stage
and the read's own EM call there (right / wrong by the parents' sides, or
none): what the rival uses that we lack, reject or call wrong.
--reads restricts the population: "islands" (reads an island block labelled,
from aln/islands.*.tsv) or a file of read names; it also prints our class
against the rival's for the whole set.
--list CLASS lists the loci (5 kb merge) of a primary cause (prefix match).
"""
import argparse
import bisect
import collections
import glob
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
    if len(ref) == 1 and len(alt) == 1:
        return "snp", 0
    if (len(ref) > 1 and len(alt) > 1) or ref[0].upper() != alt[0].upper():
        return "complex", 0
    s = (alt[1:] if len(alt) > len(ref) else ref[1:]).upper()
    p = period(s)
    if p > 6:
        return "indel", 0
    u = s[:p]
    w0 = max(0, pos - 200)
    seq = fa.fetch(ref_name, w0, pos + 200 + len(ref)).upper()
    i = pos - w0
    lo, hi = i, (i if len(alt) > len(ref) else i + len(s))
    k = 0
    while lo - 1 >= 0 and seq[lo - 1] == u[p - 1 - (k % p)]:
        lo -= 1; k += 1
    k = 0
    while hi < len(seq) and seq[hi] == u[k % p]:
        hi += 1; k += 1
    tract = hi - lo
    if tract < 2 * p + (0 if len(alt) > len(ref) else len(s)):
        return "indel", 0
    return ("hp" if p == 1 else "str"), tract


def select_reads(diag, spec):
    """Read names for --reads: "islands" (applied island reads) or a file of names."""
    if spec == "islands":
        claims = {}
        for f in D.chunk_files(os.path.join(diag, "aln"), "islands"):
            with open(f) as fh:
                for line in fh:
                    p = line.split("\t")
                    if p[9].strip() == "1":
                        claims[p[0]] = int(p[2]) + 1_000_000_000
        if not claims:
            sys.exit("no island rows: make the diag run with PGPHASE_ISLAND_MAX_ERR set")
        tags = score.read_tags(diag, set(claims))
        return {q for q, ps in claims.items() if tags.get(q, (0, 0))[1] == ps}
    with open(spec) as fh:
        return {line.split()[0] for line in fh if line.strip()}


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("diag")
    ap.add_argument("rival", nargs="?", default="hiphase")
    ap.add_argument("--reads", default="")
    ap.add_argument("--list", default="")
    ap.add_argument("--top", type=int, default=20)
    a = ap.parse_args()
    contig = os.environ["PGDEV_CONTIG"]
    fa = pysam.FastaFile(os.environ["PGDEV_REF"])
    ref_name = contig if contig in fa.references else contig.split("#")[-1]
    rdir = D.competitors()[a.rival]
    scope = {q: v for q, v in score.read_scope().items() if v[1] == "arms_phaseable"}
    ours = score.read_classes(a.diag)
    theirs = score.read_classes(rdir)
    target = [q for q in scope if ours[q] != "correct" and theirs[q] == "correct"]
    if a.reads:
        chosen = select_reads(a.diag, a.reads)
        sel = [q for q in scope if q in chosen]
        print(f"read set {a.reads}: {len(chosen):,} reads, arms_phaseable {len(sel):,}; ours x {a.rival}:")
        for (x, y), n in collections.Counter((ours[q], theirs[q]) for q in sel).most_common():
            print(f"  ours {x:10} {a.rival} {y:10} {n:6,}")
        print()
        target = [q for q in target if q in chosen]
    names = set(target)
    rtags = score.read_tags(rdir, names)
    # rival hets by phase set
    by_ps = collections.defaultdict(list)
    for pos, recs in D.vcf_records(os.path.join(rdir, "phased.vcf"), contig, 1, 2 ** 31 - 1).items():
        for ref, alts, gt, ps, c in recs:
            if ps and len(gt) == 2 and None not in gt and gt[0] != gt[1]:
                by_ps[ps].append((pos, ref, alts[max(gt) - 1]))
    for ps in by_ps:
        by_ps[ps].sort()
    truth = {}
    for pos, recs in D.vcf_records(os.environ["PGDEV_TRUTH_VCF"], contig, 1, 2 ** 31 - 1).items():
        for ref, alts, gt, ps, c in recs:
            if len(gt) == 2 and None not in gt:
                truth[pos] = "het" if gt[0] != gt[1] else "hom"
    tpos = sorted(truth)
    # our EM sites (small dump) and target reads' EM evidence
    em_sites = {}
    for pos, rows in D.load_emsites(a.diag).items():
        r = next((r for r in rows if r["in_em"]), None)
        if r is not None:
            em_sites[pos] = (r["origin"], r["err"])
    em_obs, em_mapq = {}, {}
    for f in D.chunk_files(os.path.join(a.diag, "em"), "em"):
        with open(f) as fh:
            for line in fh:
                q = line[:line.index("\t")]
                if q not in names:
                    continue
                p = line.rstrip("\n").split("\t")
                first = 6 if len(p) > 5 and ":" not in p[5] else 5
                obs = [(int(t.split(":")[0]), t.split(":")[1], int(t.split(":")[2]), float(t.split(":")[3])) for t in p[first:]]
                if q not in em_obs or len(obs) > len(em_obs[q]):
                    em_obs[q] = obs
                    em_mapq[q] = int(p[1])
    graph, aln, inj = D.load_graph(a.diag), D.load_aln(a.diag), D.load_inj(a.diag)
    scope_all = score.read_scope()
    parent_sides = collections.defaultdict(lambda: [0, 0, 0, 0])   # (pos, kind) -> mat 0/1, pat 0/1
    for f in D.chunk_files(os.path.join(a.diag, "em"), "em"):
        with open(f) as fh:
            for line in fh:
                p = line.rstrip("\n").split("\t")
                sc = scope_all.get(p[0])
                if sc is None:
                    continue
                first = 6 if len(p) > 5 and ":" not in p[5] else 5
                for t in p[first:]:
                    x = t.split(":")
                    if x[2] in ("0", "1"):
                        parent_sides[(int(x[0]), x[1])][(0 if sc[2] else 2) + int(x[2])] += 1

    def site_truth(key):
        c = parent_sides.get(key)
        if c is None or c[0] + c[1] < 3 or c[2] + c[3] < 3:
            return None, None
        ms, ps = (0 if c[0] >= c[1] else 1), (0 if c[2] >= c[3] else 1)
        if ms == ps:
            return "noise", None
        pur = (c[ms] + c[2 + ps]) / sum(c)
        return ("real" if pur >= 0.9 else "noise" if pur <= 0.65 else "mid"), ms

    def truth_of(pos, tol):
        k = bisect.bisect_left(tpos, pos - tol)
        best = "no"
        while k < len(tpos) and tpos[k] <= pos + tol:
            if truth[tpos[k]] == "het":
                return "het"
            best = "hom"
            k += 1
        return best

    def stage(pos, tol):
        for d in range(tol + 1):
            for p in ((pos,) if d == 0 else (pos - d, pos + d)):
                if p in em_sites:
                    o, e = em_sites[p]
                    return ("em_unverified" if o == "U" else "em_unreliable" if e >= 0.3 else "em_reliable"), p
        g = D.near(graph, pos, tol)
        if g:
            rows = g[0][1]
            return ("graph_voting" if any(r["voting"] for r in rows) else "graph_" + rows[0]["cat"]), g[0][0]
        i = D.near(inj, pos, tol)
        if i:
            r = min(i[0][1], key=lambda r: D.decision_rank(r["decision"]))
            d = r["decision"]
            return ("inj_" + (f"category:{r['cat']}" if d == "category" else d)), i[0][0]
        al = D.near(aln, pos, tol)
        if al:
            r = max(al[0][1], key=lambda r: r["alt_cov"])
            return f"aln_{r['init']}>{r['final']}", al[0][0]
        return "absent", None

    # MAPQ for reads the EM dump lacks
    bam = pysam.AlignmentFile(os.environ["PGDEV_BAM"])
    missing = [q for q in target if q not in em_mapq]
    for q in missing:
        s, e = scope[q][3], scope[q][4]
        for r in bam.fetch(contig, s, s + 1):
            if r.query_name == q and not r.is_secondary and not r.is_supplementary:
                em_mapq[q] = r.mapping_quality
                break
    tab = collections.Counter()
    disc_sub = collections.Counter()
    disc_where = collections.defaultdict(collections.Counter)
    detail = collections.defaultdict(list)
    called_quality = collections.Counter()
    site_tab = collections.Counter()
    for q in target:
        t = rtags.get(q)
        s, e = scope[q][3], scope[q][4]
        lst = by_ps.get(str(t[1]), []) if t else []
        i = bisect.bisect_left(lst, (s,))
        sites = []
        for pos, ref, alt in lst[i:]:
            if pos > e:
                break
            kind, tract = site_kind(fa, ref_name, pos, ref, alt)
            tol = 0 if kind == "snp" else max(3, tract)
            st, our_pos = stage(pos, tol)
            sites.append((pos, kind, tract, truth_of(pos, tol), st, our_pos))
        obs = em_obs.get(q, [])
        obs_pos = {o[0] for o in obs}
        mq = em_mapq.get(q, 60)
        th = [x for x in sites if x[3] == "het"]
        clean = [x for x in th if x[1] in ("snp", "indel")]
        if mq < 30:
            cls = "low_mapq"
        elif not th:
            cls = "rival_no_truth_site"
        elif any(x[4] == "em_reliable" and x[5] in obs_pos for x in clean):
            cls = "good_site_called"
        elif any(x[4] == "em_reliable" for x in clean):
            cls = "good_site_no_call"
        elif clean:
            best = min(clean, key=lambda x: ["em_unreliable", "em_unverified", "graph", "inj", "aln", "absent"].index(
                next((p for p in ["em_unreliable", "em_unverified", "graph", "inj", "aln", "absent"] if x[4].startswith(p)), "absent")))
            cls = "good_site_lost:" + best[4]
        else:
            cls = "repeat_only"
        tab[(ours[q], cls)] += 1
        if ours[q] == "discordant":
            right = wrong = noise = 0
            for pos, kind, side, err in obs:
                if err >= 0.3 or kind == "U":
                    continue
                st_, mat_side = site_truth((pos, kind))
                if st_ == "real":
                    if (side == mat_side) == scope[q][2]:
                        right += 1
                    else:
                        wrong += 1
                elif st_ == "noise":
                    noise += 1
            sub = ("no clear site" if right + wrong + noise == 0 else
                   "wrong calls at real sites" if wrong > right else
                   "labelled by noise sites" if noise > right + wrong else
                   "right calls (block)" if wrong == 0 else "mixed")
            disc_sub[sub] += 1
            disc_where[sub][scope[q][3] // 100000] += 1
        detail[cls].append((q, ours[q], s, e, mq, sites, obs))
        obs_at = {o[0]: o for o in obs}
        for pos, kind, tract, tr, st, our_pos in sites:
            if tr != "het":
                continue
            o = obs_at.get(our_pos) if our_pos is not None else None
            if o is None:
                call = "no call"
            else:
                st_, mat_side = site_truth((o[0], o[1]))
                call = ("call ?" if mat_side is None else
                        "call right" if (o[2] == mat_side) == scope[q][2] else "call wrong")
            stg = st if st.startswith("em_") else st.split("_")[0] + "_*" if "_" in st else st
            site_tab[(ours[q], kind, stg, call)] += 1
        if cls == "good_site_called":
            # the read's reliable EM calls against truth: does it carry wrong ones?
            rel = [o for o in obs if o[3] < 0.3 and o[1] != "U"]
            called_quality[(ours[q], "reliable calls", min(len(rel), 5))] += 1
    print(f"{a.rival}: reads it phases correctly and we do not: {len(target):,} "
          f"(unphased {sum(ours[q] == 'unphased' for q in target):,}, discordant {sum(ours[q] == 'discordant' for q in target):,})")
    classes = sorted({c for _, c in tab}, key=lambda c: -(tab[("unphased", c)] + tab[("discordant", c)]))
    print(f"  {'primary cause':44} {'unphased':>9} {'discordant':>11}")
    for c in classes:
        print(f"  {c:44} {tab[('unphased', c)]:9,} {tab[('discordant', c)]:11,}")
    print("\ndiscordant reads by their own calls (judged by the parents' sides at each site):")
    for k, n in disc_sub.most_common():
        print(f"  {k:28} {n:6,}   " + ", ".join(f"{b / 10:.1f} Mb:{c}" for b, c in disc_where[k].most_common(5)))
    print("\nrival truth-het sites in these reads: kind x our stage x the read's EM call there")
    print(f"  {'kind':8} {'our stage':16} {'read call':12} {'unphased':>9} {'discordant':>11}")
    keys = sorted({k[1:] for k in site_tab}, key=lambda k: -(site_tab[("unphased",) + k] + site_tab[("discordant",) + k]))
    for k in keys:
        print(f"  {k[0]:8} {k[1]:16} {k[2]:12} {site_tab[('unphased',) + k]:9,} {site_tab[('discordant',) + k]:11,}")
    if a.list:
        rows = sorted((r for k, v in detail.items() if k.startswith(a.list) for r in v), key=lambda r: r[2])
        # group into loci (5 kb)
        loci, cur = [], None
        for r in rows:
            if cur and r[2] <= cur["end"] + 5000:
                cur["end"] = max(cur["end"], r[3]); cur["reads"].append(r)
            else:
                cur = {"beg": r[2], "end": r[3], "reads": [r]}; loci.append(cur)
        loci.sort(key=lambda l: -len(l["reads"]))
        print(f"\n{a.list}: {len(rows)} reads in {len(loci)} loci; largest:")
        for l in loci[: a.top]:
            rs = l["reads"]
            cl = collections.Counter(r[1] for r in rs)
            sites = {}
            for r in rs:
                for x in r[5]:
                    if x[3] == "het":
                        sites[x[0]] = x
            print(f"  {l['beg']:,}-{l['end']:,}  reads {len(rs)} ({dict(cl)})  sites: " +
                  "; ".join(f"{p}:{x[1]}{x[2] or ''}:{x[4]}" for p, x in sorted(sites.items())[:6]))


if __name__ == "__main__":
    main()
