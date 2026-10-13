#!/usr/bin/env python3
"""The union EM's island blocks in a diagnostic run (evaluation only).

Usage: islands.py DIAG_DIR [--worst N] [--tsv OUT]

Needs a diag run made with PGPHASE_ISLAND_MAX_ERR set, so aln/islands.*.tsv
lists every read an island claimed. An island read is applied when its final
tag is the island's (PS = island PS + 1e9): the rescue and other chunks' labels
come first. Per island site:
  reads      applied arms_phaseable reads: correct / discordant (official classes)
  truth      a truth het within max(3, repeat run) bp of the site
  purity     parental purity of every read's call at the site (EM dump): the
             share of reads on their parent's majority side; 0.5 = no split
Production-visible features, tallied against truth so a filter can be chosen:
  err        the site's learned error;  callers  unlabelled reads calling it
  linked     most member reads sharing one other verified site
  assoc      among reads the EM labelled, agreement of the site's allele with
             their haplotype, per block (sum of max(same, cross) / n); n = reads
  minor      minor allele fraction among reads calling the site
  run        reference repeat run at the site (bp)
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

OFF = 1_000_000_000


def repeat_run(fa, ref_name, pos):
    """Longest unit-1..6 tandem run covering pos or the base after it (a VCF
    indel is anchored on the base before its tract), in bp."""
    lo = max(0, pos - 80)
    s = fa.fetch(ref_name, lo, pos + 80).upper()
    best = 1
    for c in (pos - 1 - lo, pos - lo):
        for unit in range(1, 7):
            for start in range(max(0, c - 70), c + 1):
                u = s[start:start + unit]
                if len(u) < unit or "N" in u:
                    continue
                n = 1
                while s[start + n * unit:start + (n + 1) * unit] == u:
                    n += 1
                span = n * unit
                if n >= 2 and start <= c < start + span and span > best:
                    best = span
    return best


def purity(par):
    """Parental purity from [mat side 0, mat side 1, pat side 0, pat side 1]."""
    if not (par[0] + par[1] and par[2] + par[3]):
        return None
    ms, pp = (0 if par[0] >= par[1] else 1), (0 if par[2] >= par[3] else 1)
    return 0.5 if ms == pp else (par[ms] + par[2 + pp]) / sum(par)


def load_islands(diag):
    rows = []
    for f in D.chunk_files(os.path.join(diag, "aln"), "islands"):
        with open(f) as fh:
            for line in fh:
                p = line.rstrip("\n").split("\t")
                rows.append(dict(read=p[0], pos=int(p[1]), ps=int(p[2]), side=int(p[3]), hap=int(p[4]),
                                 rerr=float(p[5]), err=float(p[6]), callers=int(p[7]), linked=int(p[8]),
                                 kept=p[9] == "1", kind=p[10] if len(p) > 10 else "E"))
    return rows


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("diag")
    ap.add_argument("--worst", type=int, default=12)
    ap.add_argument("--tsv", default="")
    a = ap.parse_args()
    contig = os.environ["PGDEV_CONTIG"]
    fa = pysam.FastaFile(os.environ["PGDEV_REF"])
    ref_name = contig if contig in fa.references else contig.split("#")[-1]
    claims = load_islands(a.diag)
    if not claims:
        sys.exit("no aln/islands.*.tsv rows: make the diag run with PGPHASE_ISLAND_MAX_ERR set")
    scope_all = score.read_scope()
    cls = score.read_classes(a.diag)
    tags = score.read_tags(a.diag, set(scope_all))
    sites = {}
    for c in claims:
        s = sites.setdefault((c["pos"], c["ps"]), dict(pos=c["pos"], ps=c["ps"], err=c["err"], callers=c["callers"],
                                                       linked=c["linked"], kept=c["kept"], kind=c["kind"], members=[], c=0, d=0, o=0,
                                                       applied=0))
        s["members"].append(c["read"])
        t = tags.get(c["read"])
        if c["kept"] and t and t[1] == c["ps"] + OFF:
            s["applied"] += 1
            sc = scope_all.get(c["read"])
            if sc and sc[1] == "arms_phaseable":
                s["c" if cls[c["read"]] == "correct" else "d"] += 1
            else:
                s["o"] += 1
    want = {s["pos"] for s in sites.values()}
    # every read's call at the island sites, and the EM labels, from the EM dump
    calls = collections.defaultdict(dict)
    unl = collections.defaultdict(dict)  # unlabelled read -> rejected site -> side
    for f in D.chunk_files(os.path.join(a.diag, "em"), "em"):
        with open(f) as fh:
            for line in fh:
                p = line.rstrip("\n").split("\t")
                hap, ps = int(p[2]), int(p[3])
                first = 6 if len(p) > 5 and ":" not in p[5] else 5
                for t in p[first:]:
                    x = t.split(":")
                    if x[2] not in ("0", "1"):
                        continue
                    pos = int(x[0])
                    if pos in want and p[0] not in calls[pos]:
                        calls[pos][p[0]] = (int(x[2]), hap, ps)
                    if hap not in (1, 2) and x[1] not in ("L", "U") and float(x[3]) >= 0.3:
                        unl[p[0]].setdefault(pos, int(x[2]))
    by_site = collections.defaultdict(list)
    for q, d in unl.items():
        for pos in d:
            if pos in want:
                by_site[pos].append(q)
    truth = []
    for pos, recs in D.vcf_records(os.environ["PGDEV_TRUTH_VCF"], contig, 1, 2 ** 31 - 1).items():
        if any(len(gt) == 2 and None not in gt and gt[0] != gt[1] for _, _, gt, _, _ in recs):
            truth.append(pos)
    truth.sort()
    for s in sites.values():
        pos = s["pos"]
        s["run"] = repeat_run(fa, ref_name, pos)
        tol = max(3, s["run"])
        i = bisect.bisect_left(truth, pos - tol)
        s["truth"] = i < len(truth) and truth[i] <= pos + tol
        par = [0, 0, 0, 0]
        par_lab, par_unl = [0, 0, 0, 0], [0, 0, 0, 0]
        lab_right = lab_n = 0
        blocks = collections.defaultdict(lambda: [0, 0])
        n1 = 0
        for q, (side, hap, ps) in calls[pos].items():
            n1 += side
            sc = scope_all.get(q)
            if sc is not None:
                par[(0 if sc[2] else 2) + side] += 1
                (par_lab if hap in (1, 2) else par_unl)[(0 if sc[2] else 2) + side] += 1
                if hap in (1, 2):
                    lab_n += 1
                    lab_right += cls.get(q) == "correct"
            if hap in (1, 2) and ps > 0:
                blocks[ps][(side == 0) == (hap == 1)] += 1
        n = len(calls[pos])
        s["minor"] = min(n1, n - n1) / n if n else 0.0
        s["purity"] = purity(par)
        s["purity_lab"], s["purity_unl"] = purity(par_lab), purity(par_unl)
        # co-segregation with another rejected site among unlabelled reads
        pair = collections.defaultdict(lambda: [0, 0])
        for q in by_site.get(pos, []):
            mine = unl[q][pos]
            for j, sj in unl[q].items():
                if j != pos:
                    pair[j][0 if sj == mine else 1] += 1
        best, best_n = None, 0
        for j, (ag, dg) in pair.items():
            if ag + dg >= 4:
                c = max(ag, dg) / (ag + dg)
                if best is None or c > best or (c == best and ag + dg > best_n):
                    best, best_n = c, ag + dg
        s["pair"], s["pair_n"] = best, best_n
        s["lab_correct"] = lab_right / lab_n if lab_n else None
        s["assoc_n"] = sum(sum(v) for v in blocks.values())
        s["assoc"] = sum(max(v) for v in blocks.values()) / s["assoc_n"] if s["assoc_n"] else None
    msa = D.load_msa(a.diag)
    msa_pos = sorted(msa)
    for s in sites.values():
        tol = max(3, s["run"])
        i = bisect.bisect_left(msa_pos, s["pos"] - tol)
        near = [msa_pos[j] for j in range(i, len(msa_pos)) if msa_pos[j] <= s["pos"] + tol]
        if not near:
            s["msa"] = "no MSA call"
        else:
            best = min(near, key=lambda p: abs(p - s["pos"]))
            r = msa[best][0]
            s["msa"] = ("MSA " + r["cat"] + (" 1 cluster" if min(r["clusters"]) == 0 else " 2 clusters"))
    rows = list(sites.values())
    kept = [s for s in rows if s["kept"]]
    print(f"{contig}: island sites {len(rows):,} (kept {len(kept):,}); reads claimed {len(claims):,}, applied "
          f"{sum(s['applied'] for s in kept):,}; applied arms_phaseable correct {sum(s['c'] for s in kept):,} "
          f"discordant {sum(s['d'] for s in kept):,}")

    def tally(name, pred, pool=kept):
        sub = [s for s in pool if pred(s)]
        c, d = sum(s["c"] for s in sub), sum(s["d"] for s in sub)
        nt = sum(1 for s in sub if s["truth"])
        print(f"  {name:40} sites {len(sub):5,} (truth het {nt:4,})  correct {c:6,}  discordant {d:5,}  "
              f"ratio {c / max(1, d):6.1f}")

    def bins(name, key, edges, pool=kept):
        for lo, hi in zip(edges, edges[1:]):
            tally(f"{name} [{lo}, {hi})", lambda s, lo=lo, hi=hi: s[key] is not None and lo <= s[key] < hi, pool)
        if any(s[key] is None for s in pool):
            tally(f"{name} n/a", lambda s: s[key] is None, pool)

    print("site kind (E: site the EM rejects, R: repeat het indel the EM never sees, M: MSA-verified het injection drops):")
    tally("kind E", lambda s: s["kind"] == "E")
    tally("kind R", lambda s: s["kind"] == "R")
    tally("kind M", lambda s: s["kind"] == "M")
    print("the alignment solve's MSA at the site (kept islands):")
    for m in sorted({s["msa"] for s in kept}):
        tally(m, lambda s, m=m: s["msa"] == m)
    print("truth (kept islands):")
    tally("truth het at the site", lambda s: s["truth"])
    tally("no truth het", lambda s: not s["truth"])
    print("parental purity of every read's call at the site:")
    bins("purity", "purity", [0.0, 0.6, 0.75, 0.9, 1.01])
    print("production-visible features (kept islands):")
    bins("assoc", "assoc", [0.0, 0.6, 0.7, 0.8, 0.9, 1.01])
    bins("assoc reads", "assoc_n", [0, 1, 5, 10, 20, 40, 10 ** 6])
    bins("pair", "pair", [0.0, 0.6, 0.7, 0.8, 0.9, 0.95, 1.01])
    bins("err", "err", [0.3, 0.33, 0.36, 0.4, 0.5, 1.0])
    bins("callers", "callers", [2, 3, 5, 10, 20, 10 ** 6])
    bins("linked", "linked", [0, 1, 3, 5, 10, 10 ** 6])
    bins("minor", "minor", [0.0, 0.2, 0.3, 0.4, 0.51])
    bins("run", "run", [1, 4, 8, 13, 20, 30, 10 ** 6])
    print("truth het vs production features, all island sites (kept or not): share of truth hets per bin")
    for key, edges in (("pair", [0.0, 0.6, 0.7, 0.8, 0.9, 0.95, 1.01]), ("assoc", [0.0, 0.6, 0.7, 0.8, 0.9, 1.01]), ("minor", [0.0, 0.2, 0.3, 0.4, 0.51]),
                       ("run", [1, 4, 8, 13, 20, 30, 10 ** 6]), ("err", [0.3, 0.33, 0.36, 0.4, 0.5, 1.0])):
        cells = []
        for lo, hi in zip(edges, edges[1:]):
            sub = [s for s in rows if s[key] is not None and lo <= s[key] < hi]
            if sub:
                cells.append(f"[{lo},{hi}) {sum(s['truth'] for s in sub)}/{len(sub)}")
        print(f"  {key:8} " + "  ".join(cells))
    print("calls at the site by read group (all island sites): median parental purity, labelled reads' class")
    for name, pred in (("truth het", lambda s: s["truth"]), ("no truth het", lambda s: not s["truth"])):
        sub = [s for s in rows if pred(s)]

        def med(key):
            v = sorted(s[key] for s in sub if s[key] is not None)
            return f"{v[len(v) // 2]:.2f} (n {len(v)})" if v else "n/a"
        print(f"  {name:13} sites {len(sub):4}  purity all {med('purity')}  labelled reads {med('purity_lab')}  "
              f"unlabelled reads {med('purity_unl')}  labelled reads correct {med('lab_correct')}")
    print("worst kept islands:")
    for s in sorted(kept, key=lambda s: -s["d"])[:a.worst]:
        f = lambda v: "n/a" if v is None else f"{v:.2f}"  # noqa: E731
        print(f"  {s['pos']:>11,} run {s['run']:3} truth {'het' if s['truth'] else '-  '} correct {s['c']:3} "
              f"discordant {s['d']:3}  purity {f(s['purity'])}  assoc {f(s['assoc'])} (n {s['assoc_n']})  "
              f"err {s['err']:.3f} callers {s['callers']} linked {s['linked']} minor {s['minor']:.2f}")
    if a.tsv:
        keys = ["pos", "ps", "kind", "kept", "applied", "c", "d", "o", "truth", "purity", "purity_lab", "purity_unl",
                "lab_correct", "pair", "pair_n", "assoc", "assoc_n", "err",
                "callers", "linked", "minor", "run", "msa"]
        with open(a.tsv, "w") as fh:
            fh.write("\t".join(keys) + "\n")
            for s in sorted(rows, key=lambda s: s["pos"]):
                fh.write("\t".join(str(s[k]) for k in keys) + "\n")
        print(f"per-island table: {a.tsv}")


if __name__ == "__main__":
    main()
