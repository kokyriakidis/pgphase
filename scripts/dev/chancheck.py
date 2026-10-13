#!/usr/bin/env python3
"""Per-read calls at catalog SNPs: the EM's call versus the alignment solve's, against truth (evaluation only).

Usage: chancheck.py DIAG_DIR

For simple biallelic catalog SNPs whose ALT is truth's ALT, every scored
read's final EM call (emreads dump) and alignment-solve call (alnreads dump)
are compared with the allele its truth haplotype carries. Reports, per
combination of channels present (both agree / both disagree / EM only /
alignment only), how often each is right, and the same for reads we leave
unphased or discordant.
"""
import collections
import glob
import multiprocessing as mp
import os
import sys

sys.path.insert(0, "scripts/dev")
import diaglib as D  # noqa: E402
import score  # noqa: E402

G = {}


def scan(path_pair):
    em_path, kind = path_pair
    out = {}
    with open(em_path) as fh:
        for line in fh:
            p = line.rstrip("\n").split("\t")
            q = p[0]
            if q not in G["scope"]:
                continue
            calls = {}
            if kind == "em":
                first = 6 if len(p) > 5 and ":" not in p[5] else 5
                for tok in p[first:]:
                    x = tok.split(":")
                    pos = int(x[0])
                    if x[1] == "C" and pos in G["sites"] and pos in G["xy"]:
                        calls[pos] = G["xy"][pos][0] if x[2] == "0" else G["xy"][pos][1]
            else:
                for tok in p[5:]:
                    x = tok.split(":")
                    pos = int(x[0])
                    if pos in G["sites"]:
                        calls[pos] = int(x[-1])
            if calls:
                old = out.get(q)
                if old is None or len(calls) > len(old):
                    out[q] = calls
    return kind, out


def main():
    diag = sys.argv[1]
    scope = score.read_scope()
    tp = score.truth_tables()["truth_phased"]
    # simple catalog SNP rows in the EM: pos -> truth hap-1 allele (0/1)
    import pysam
    em_c = set()
    for pos, rows in D.load_emsites(diag).items():
        if any(r["origin"] == "C" and r["type"] == "S" and r["ref_len"] == 1 for r in rows):
            em_c.add(pos)
    sites = {}
    tf = pysam.TabixFile(os.environ["PGDEV_SITES"])
    for line in tf.fetch(D.tabix_contig(tf, os.environ["PGDEV_CONTIG"])):
        f = line.split("\t", 6)
        pos = int(f[1])
        if pos in em_c and len(f[3]) == 1 and len(f[4]) == 1:
            sites[pos] = f[4].upper()
    truth = {}
    for (pos, ref, alt), g0 in tp.items():
        if pos in sites and len(ref) == 1 and alt == sites[pos]:
            truth[pos] = g0
    G["sites"] = set(truth)
    G["scope"] = {q: v for q, v in scope.items()}
    xy = {}
    for pos, rows in D.load_emsites(diag).items():
        for r in rows:
            if r["origin"] == "C" and r["in_em"]:
                xy[pos] = (r["x"], r["y"])
    G["xy"] = xy
    jobs = [(f, "em") for f in D.chunk_files(os.path.join(diag, "em"), "em")] + \
           [(f, "aln") for f in D.chunk_files(os.path.join(diag, "aln"), "alnreads")]
    with mp.Pool(int(os.environ.get("PGDEV_THREADS", 10))) as pool:
        parts = pool.map(scan, jobs)
    em, al = {}, {}
    for kind, out in parts:
        dst = em if kind == "em" else al
        for q, c in out.items():
            if q not in dst or len(c) > len(dst[q]):
                dst[q] = c
    # which truth haplotype gt[0] is: calibrate on alignment calls
    agree = collections.Counter()
    for q, calls in al.items():
        for pos, a in calls.items():
            if a in (0, 1):
                agree[(a == truth[pos]) == scope[q][2]] += 1
    first_mat = agree[True] >= agree[False]
    cls = score.read_classes(diag)
    aln_sites = set()
    for calls in al.values():
        aln_sites.update(calls)
    tab = collections.defaultdict(collections.Counter)
    for q in set(em) | set(al):
        mat = scope[q][2]
        group = cls.get(q, "?")
        e, a = em.get(q, {}), al.get(q, {})
        for pos in set(e) | set(a):
            exp = truth[pos] if mat == first_mat else 1 - truth[pos]
            ce, ca = e.get(pos), a.get(pos)
            ce = ce if ce in (0, 1) else None
            ca = ca if ca in (0, 1) else None
            if ce is None and ca is None:
                continue
            if ce is not None and ca is not None:
                key = "both agree" if ce == ca else "both disagree"
            elif ce is not None:
                key = "EM only, site in alignment" if pos in aln_sites else "EM only, site not in alignment"
            else:
                key = "alignment only"
            right_e = ce == exp if ce is not None else None
            right_a = ca == exp if ca is not None else None
            tab[group][(key, right_e, right_a)] += 1
    print(f"simple catalog SNPs matching a truth het: {len(truth):,}  (truth first haplotype maternal: {first_mat})")
    for group in ("correct", "discordant", "unphased"):
        t = tab[group]
        print(f"\nreads we call {group}:")
        for key in ("both agree", "both disagree", "EM only, site in alignment", "EM only, site not in alignment", "alignment only"):
            rows = {k: n for k, n in t.items() if k[0] == key}
            tot = sum(rows.values())
            if not tot:
                continue
            em_right = sum(n for k, n in rows.items() if k[1])
            al_right = sum(n for k, n in rows.items() if k[2])
            em_n = sum(n for k, n in rows.items() if k[1] is not None)
            al_n = sum(n for k, n in rows.items() if k[2] is not None)
            print(f"  {key:32} calls {tot:9,}   EM right {100 * em_right / max(1, em_n):6.2f}%   "
                  f"alignment right {100 * al_right / max(1, al_n):6.2f}%")


if __name__ == "__main__":
    main()
