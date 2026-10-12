#!/usr/bin/env python3
"""Why phaseable arm reads stay unphased, and which sites flip (evaluation only).

Usage: gaps.py RUN_DIR EM_DUMP_DIR [--out TSV]
Uses PGDEV_* (score.py) and PGDEV_COMPETITORS. For every arms_phaseable read
the run leaves unphased: its class (why.py's decision tree), MAPQ, and
whether each competitor phases it correctly. For every flip (a truth-matched
het that disagrees with its block while both neighbours agree): the site's
VCF category, kind in the EM (C clean, O other, I injected, L window) and
learned error.
"""
import argparse
import collections
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import score  # noqa: E402
import why  # noqa: E402


def flips(run_dir):
    tp = score.truth_tables()["truth_phased"]
    vcf = os.path.join(run_dir, "phased.vcf")
    by_ps = collections.defaultdict(list)
    for pos, ref, alts, gt, ps, cat in score.vcf_records(vcf):
        if ps is None or len(alts) != 1 or set(gt) != {0, 1}:
            continue
        key = (pos, ref.upper(), alts[0].upper())
        if key in tp and not (score.CEN[0] <= pos - 1 < score.CEN[1]):
            by_ps[ps].append((pos, gt[0] == tp[key], cat, len(ref) == 1 and len(alts[0]) == 1))
    out = []
    for sites in by_ps.values():
        sites.sort()
        for i in range(1, len(sites) - 1):
            if sites[i][1] != sites[i - 1][1] and sites[i + 1][1] == sites[i - 1][1]:
                out.append(sites[i])
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("run")
    ap.add_argument("dump")
    ap.add_argument("--out")
    a = ap.parse_args()
    scope = {q: v for q, v in score.read_scope().items() if v[1] == "arms_phaseable"}
    ours = score.classify_reads(a.run, scope)
    comps = {}
    for pair in os.environ.get("PGDEV_COMPETITORS", "").split():
        name, d = pair.split("=", 1)
        comps[name] = score.classify_reads(d, scope)
    dump = why.load_dump(a.dump)
    unph = [q for q in scope if ours[q] == "unphased"]
    tab, sub = collections.Counter(), collections.Counter()
    comp_ok = collections.defaultdict(collections.Counter)
    rows = []
    for q in unph:
        d = dump.get(q)
        cls, detail = why.classify(d, "unphased", scope[q][2])
        mq = d["mapq"] if d else -1
        tab[cls] += 1
        sub[(cls, detail.split(";")[0] if detail else "", "MAPQ>=20" if mq >= 20 else "MAPQ<20" if mq >= 0 else "-")] += 1
        for name, c in comps.items():
            comp_ok[cls][name] += c[q] == "correct"
        rows.append((q, scope[q][3], scope[q][4], cls, detail, mq) + tuple(comps[n][q] for n in comps))
    print(f"unphased arms_phaseable reads: {len(unph):,} of {len(scope):,} ({100 * len(unph) / len(scope):.2f}%)")
    for cls, n in tab.most_common():
        cc = ", ".join(f"{name} correct {comp_ok[cls][name]:,}" for name in comps)
        print(f"  {cls:16} {n:6,}   ({cc})")
        for (c, det, mq), m in sorted(sub.items(), key=lambda kv: -kv[1]):
            if c == cls and m >= max(20, n // 20):
                print(f"      {m:6,}  {det or '-':28} {mq}")
    bins = collections.Counter(int(r[1]) // 1_000_000 for r in rows)
    print("  unphased by 1 Mb bin (top):", bins.most_common(8))
    fl = flips(a.run)
    print(f"\nflips on the arms: {len(fl)}")
    print("  by category:", collections.Counter(f[2] for f in fl).most_common())
    print("  SNP vs indel:", collections.Counter("SNP" if f[3] else "indel" for f in fl).most_common())
    print("  positions:", [f[0] for f in fl][:12])
    if a.out:
        with open(a.out, "w") as fh:
            fh.write("read\tstart\tend\tclass\tdetail\tmapq\t" + "\t".join(comps) + "\n")
            for r in rows:
                fh.write("\t".join(map(str, r)) + "\n")


if __name__ == "__main__":
    main()
