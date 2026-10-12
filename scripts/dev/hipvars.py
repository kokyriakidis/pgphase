#!/usr/bin/env python3
"""For arm reads the competitor phases correctly and we leave unphased: which
truth variants does the competitor phase inside the read, and what is each one
in our solve? (evaluation only)

Usage: hipvars.py RUN_DIR EM_DUMP_DIR COMPETITOR_DIR
A read takes the best status over the competitor's truth-matching phased hets
in its span, in this order:
  read_call_reliable   EM site, the read has a call, site error < 0.3 (should label)
  no_call_at_site      EM site the read has no call at
  site_unreliable_pure EM site with error >= 0.2 whose calls separate the parents (>= 0.9)
  site_bad_calls       EM site whose calls do not separate the parents
  candidate_not_site   in our candidates.tsv, not an EM site (category shown)
  absent               not in our candidates at all
"""
import bisect
import collections
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import score  # noqa: E402
import why  # noqa: E402

ORDER = ["read_call_reliable", "no_call_at_site", "site_unreliable_pure", "site_bad_calls",
         "candidate_not_site", "absent"]


def main():
    run, dump_dir, comp = sys.argv[1:4]
    scope = {q: v for q, v in score.read_scope().items() if v[1] == "arms_phaseable"}
    ours, theirs = score.classify_reads(run, scope), score.classify_reads(comp, scope)
    reads = {q for q in scope if theirs[q] == "correct" and ours[q] == "unphased"}
    tp = score.truth_tables()["truth_phased"]
    # competitor phased hets that match truth exactly
    vcf = os.path.join(comp, "phased.vcf")
    vcf = vcf if os.path.exists(vcf) else vcf + ".gz"
    hv = []
    for pos, ref, alts, gt, ps, cat in score.vcf_records(vcf):
        if ps and len(gt) == 2 and gt[0] != gt[1] and None not in gt and len(alts) == 1 \
                and (pos, ref.upper(), alts[0].upper()) in tp:
            hv.append((pos, ref.upper(), alts[0].upper()))
    hv.sort()
    hpos = [v[0] for v in hv]
    dump = why.load_dump(dump_dir)
    purity = why.site_purity(dump, score.read_scope())
    cand = collections.defaultdict(set)
    with open(os.path.join(run, "candidates.tsv")) as fh:
        hdr = next(fh).rstrip("\n").split("\t")
        for line in fh:
            d = dict(zip(hdr, line.rstrip("\n").split("\t")))
            cand[int(d["POS"])].add(d["CATEGORY"])

    def site_at(v, positions):
        pos, ref, alt = v
        snp = len(ref) == 1 and len(alt) == 1
        for p in ((pos, pos - 1) if snp else range(pos - 2, pos + 3)):
            if p in positions:
                return p
        return None

    tab, vtype, cats, ex = collections.Counter(), collections.Counter(), collections.Counter(), collections.defaultdict(list)
    per_var = collections.Counter()
    for q in reads:
        s, e = scope[q][3], scope[q][4]
        vs = hv[bisect.bisect_left(hpos, s + 1):bisect.bisect_right(hpos, e)]
        d = dump.get(q)
        obs = {o[0]: o for o in d["obs"]} if d else {}
        best, best_v = None, None
        for v in vs:
            p = site_at(v, purity)
            if p is not None:
                n, pur, err, kind = purity[p]
                o = obs.get(p)
                if o is not None and o[3] < 0.3 and pur >= 0.9:
                    st = "read_call_reliable"
                elif o is None and err < 0.3:
                    st = "no_call_at_site"
                elif pur >= 0.9 and err >= 0.2:
                    st = "site_unreliable_pure"
                elif pur < 0.9:
                    st = "site_bad_calls"
                else:
                    st = "no_call_at_site"
            else:
                c = next((cand[p] for p in range(v[0] - 2, v[0] + 3) if p in cand), None)
                st = "candidate_not_site" if c else "absent"
            if best is None or ORDER.index(st) < ORDER.index(best):
                best, best_v = st, v
        if best is None:
            best = "competitor has no truth het in span"
        tab[best] += 1
        if best_v:
            kind = "SNP" if len(best_v[1]) == 1 and len(best_v[2]) == 1 else \
                   f"{'INS' if len(best_v[2]) > len(best_v[1]) else 'DEL'}{abs(len(best_v[2]) - len(best_v[1]))}"
            vtype[(best, kind if kind == "SNP" else ("indel<=2" if abs(len(best_v[2]) - len(best_v[1])) <= 2 else "indel>2"))] += 1
            per_var[(best, best_v[0])] += 1
            if best == "candidate_not_site":
                c = next((cand[p] for p in range(best_v[0] - 2, best_v[0] + 3) if p in cand), set())
                cats[",".join(sorted(c))] += 1
        if len(ex[best]) < 4 and best_v:
            ex[best].append(f"{best_v[0]}:{best_v[1][:6]}>{best_v[2][:6]}")
    print(f"reads competitor-correct, ours unphased (arms_phaseable): {len(reads):,}")
    for st, n in tab.most_common():
        kinds = ", ".join(f"{k[1]} {v}" for k, v in vtype.items() if k[0] == st)
        loci = len({p for (b, p) in per_var if b == st})
        print(f"  {st:24} {n:5,} reads at {loci:4} loci  ({kinds})  e.g. {' '.join(ex[st])}")
    if cats:
        print("  candidate_not_site categories:", cats.most_common())
    print("  loci carrying most reads:", [(b, p, n) for (b, p), n in per_var.most_common(12)])


if __name__ == "__main__":
    main()
