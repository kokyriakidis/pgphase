#!/usr/bin/env python3
"""Everything known about one region, by position (evaluation only).

Usage: inspect_region.py DIAG_DIR REGION [--all] [--reads N]

REGION is chrN:A-B, A-B (contig from PGDEV_CONTIG) or a single position
(+-5 kb). DIAG_DIR comes from `scripts/dev/pg diag`.

Prints
  reads    arm and centromere reads crossing a truth het that overlap the
           region: correct / discordant / unphased for us and each rival
           (official classes: each phase set oriented by all scored reads),
           and for reads a rival phases correctly that we do not, why.py's
           reason on our side
  noisy    the alignment solve's noisy regions overlapping the region
  sites    one row per position with any of: a truth het, a record in our or
           a rival's VCF, an EM site, an injected or hidden-complement
           alignment call, a clean alignment call that was suppressed, a
           voting graph row (--all: every catalog / graph / alignment row)
    truth   type and GT
    catalog REF length / ALT count / nesting level (first record)
    graph   category, voting flag (v), allele coverages
    aln     initial>final category of the alignment call, REF/ALT reads,
            U when kept from a noisy region too long for MSA
    inj     injection decision (inj, comp = hidden complement, depth,
            conf = too few MAPQ>=30 reads, cat = category not injected)
    em      EM kind (C clean, O other, I injected, U unverified, L window),
            learned error, reads, purity against the parents, block
    alleles the final EM's read counts per allele at the candidate
            (allele:reads), u = reads in range without an allele, and the
            two alleles contrasted (x/y) or "out" (left out of the EM)
    ours / rivals  GT, phase-set letter, and = / x for agreement with the
            truth orientation (blocks are compared by letter)
  msa      every noisy-region MSA variant in the region: cluster sizes,
           consensus coverage (a large deletion on one haplotype shows as a
           long stretch without consensus bases), which consensus carries
           it, and read calls before and after the refresh (LOST CALLS when
           the refresh dropped some)
  joins    EM block-cut decisions in and around the region: every cut, and
           joins carried only by weak sites (--all: every join)
"""
import argparse
import collections
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import diaglib as D  # noqa: E402
import score  # noqa: E402
import why  # noqa: E402


def vtype(ref, alt):
    if len(ref) == 1 and len(alt) == 1:
        return "snp"
    d = len(alt) - len(ref)
    return f"ind{d:+d}" if len(ref) == 1 or len(alt) == 1 else f"cx{d:+d}"


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("diag")
    ap.add_argument("region")
    ap.add_argument("--all", action="store_true")
    ap.add_argument("--reads", type=int, default=0, help="list up to N rival-correct reads we miss")
    a = ap.parse_args()
    contig, beg, end = D.parse_region(a.region)
    comps = D.competitors()
    print(f"region {contig}:{beg:,}-{end:,}  diag {a.diag}")

    # reads
    scope = score.read_scope()
    tools = {"ours": a.diag, **comps}
    cls = {name: score.read_classes(d) for name, d in tools.items()}
    inreg = [q for q, v in scope.items() if v[3] < end and v[4] > beg and v[1].endswith("phaseable")]
    em_reads, em_sites = D.load_em(a.diag, beg - 100_000, end + 100_000)
    print(f"\nreads crossing a truth het ({len(inreg):,}):")
    for name in tools:
        c = collections.Counter(cls[name][q] for q in inreg)
        print(f"  {name:10} correct {c['correct']:5,}  discordant {c['discordant']:5,}  unphased {c['unphased']:5,}")
    miss = [q for q in inreg if cls["ours"][q] != "correct" and any(cls[n][q] == "correct" for n in comps)]
    reasons = collections.Counter()
    for q in miss:
        r, det = why.classify(em_reads.get(q), cls["ours"][q], scope[q][2])
        reasons[f"{cls['ours'][q]}:{r}"] += 1
    print(f"  rival correct, ours not: {len(miss):,}  " + ", ".join(f"{k} {v}" for k, v in reasons.most_common()))
    for q in sorted(miss, key=lambda q: scope[q][3])[:a.reads]:
        d = em_reads.get(q)
        obs = " ".join(f"{o[0]}:{o[1]}{o[2]}:{o[3]:.2f}" for o in (d["obs"] if d else []))
        print(f"    {q.split('/')[1]} {scope[q][3]:,}-{scope[q][4]:,} {'M' if scope[q][2] else 'P'} "
              f"ours={cls['ours'][q]} " + " ".join(f"{n}={cls[n][q]}" for n in comps) + f" | {obs}")

    # noisy regions
    noisy = D.load_noisy(a.diag, beg, end)
    print("\nnoisy regions: " + (", ".join(f"{x:,}-{y:,} ({y - x + 1:,} bp)" for x, y in noisy) or "none"))

    # sites
    truth = D.truth_records(contig, beg, end)
    tp = score.truth_tables()["truth_phased"]
    cat = D.catalog(contig, beg, end)
    graph = D.load_graph(a.diag, beg, end)
    aln = D.load_aln(a.diag, beg, end)
    inj = D.load_inj(a.diag, beg, end)
    ems = D.load_emsites(a.diag, beg, end)
    vcfs = {name: D.vcf_records(os.path.join(d, "phased.vcf"), contig, beg, end) for name, d in tools.items()}
    letters = {name: {} for name in tools}

    def vcf_cell(name, pos):
        out = []
        for ref, alts, gt, ps, c in vcfs[name].get(pos, []):
            if len(gt) != 2 or None in gt or gt[0] == gt[1]:
                continue
            g = f"{gt[0]}{'|' if ps else '/'}{gt[1]}"
            if ps:
                lt = letters[name].setdefault(ps, chr(ord("A") + len(letters[name]) % 26))
                g += lt
                k = (pos, ref.upper(), alts[0].upper()) if len(alts) == 1 else None
                if k in tp:
                    g += "=" if gt[0] == tp[k] else "x"
            out.append(g if len(ref) == len(alts[0]) == 1 else g + vtype(ref, alts[0]))
        return ",".join(out)

    rows = set(truth) | {p for n in vcfs for p in vcfs[n]} | {p for p in em_sites if beg <= p <= end}
    rows |= {p for p, l in inj.items() if any(r["decision"] in ("injected", "complement") for r in l)}
    rows |= {p for p, l in aln.items() if any(r["init"].startswith("CLEAN_HET") and r["final"] == "NON_VAR" for r in l)}
    rows |= {p for p, l in graph.items() if any(r["voting"] for r in l)}
    if a.all:
        rows |= set(cat) | set(graph) | set(aln)
    names = list(tools)
    print(f"\n{'pos':>11} {'truth':10} {'catalog':9} {'graph':16} {'aln':16} {'inj':6} {'em':22} {'alleles':20} " +
          " ".join(f"{n[:9]:10}" for n in names))
    for pos in sorted(rows):
        t = truth.get(pos)
        tcell = f"{vtype(t[0], t[1].split(',')[0])} {t[2][0]}/{t[2][1]}" if t else ""
        c = D.near(cat, pos)
        ccell = f"{c[0][1][0][0]}/{c[0][1][0][1]}L{c[0][1][0][2]}" if c else ""
        g = D.near(graph, pos)
        gcell = ""
        if g:
            r = g[0][1][0]
            gcell = f"{D.SHORT.get(r['cat'], r['cat'])}{'v' if r['voting'] else ''} {r['covs']}"[:16]
        al = aln.get(pos, [])
        acell = ""
        if al:
            r = max(al, key=lambda r: r["alt_cov"])
            acell = (f"{D.SHORT.get(r['init'], '?')}>{D.SHORT.get(r['final'], '?')} {r['ref']}/{r['alt_cov']}"
                     + ("U" if r["unresolved"] else ""))
        ij = inj.get(pos, [])
        icell = ""
        if ij:
            short = {"injected": "inj", "complement": "comp", "depth_gate": "depth", "few_confident_reads": "conf",
                     "category": "cat", "outside_chunk": "out", "not_germline": "ng"}
            best = min(ij, key=lambda r: ["injected", "complement", "depth_gate", "few_confident_reads",
                                          "category", "not_germline", "outside_chunk"].index(r["decision"]))
            icell = short.get(best["decision"], best["decision"])
        e = D.near(em_sites, pos)
        ecell = ""
        if e:
            s = e[0][1]
            pur = f"{s['purity']:.2f}" if s["purity"] is not None else "-"
            ecell = f"{s['kind']} {s['err']:.2f} {s['n']}r p{pur} b{s['block'] % 100000}"
        es = D.near(ems, pos)
        scell = ""
        if es:
            r = es[0][1][0]
            scell = f"{r['counts']} u{r['undecided']} " + (f"{r['x']}/{r['y']}" if r["in_em"] else "out")
        print(f"{pos:>11,} {tcell:10} {ccell:9} {gcell:16} {acell:16} {icell:6} {ecell:22} {scell[:20]:20} " +
              " ".join(f"{vcf_cell(n, pos):10}" for n in names))
    for n in names:
        if letters[n]:
            print(f"  {n} blocks: " + ", ".join(f"{v}={k}" for k, v in letters[n].items()))

    # noisy-region MSA calls
    msa = D.load_msa(a.diag, beg, end)
    if msa:
        print("\nMSA variants (cluster reads; per consensus: reference bases aligned / longest stretch without "
              "a consensus base; which consensus carries the variant; read calls ref,alt,missing before > after "
              "the refresh):")
        for pos in sorted(msa):
            for r in msa[pos]:
                lost = r["after"][0] + r["after"][1] < r["before"][0] + r["before"][1]
                print(f"  {pos:>11,} {r['type']}{r['ref_len']} {r['alt'][:12]:12} {D.SHORT.get(r['cat'], r['cat']):5} "
                      f"region {r['reg_beg']:,}+{r['reg_len']} clusters {r['clusters'][0]}/{r['clusters'][1]} "
                      f"cons {r['aligned'][0]}/{r['longest_gap'][0]} {r['aligned'][1]}/{r['longest_gap'][1]} "
                      f"in {'1' if r['in_cons'][0] else '-'}{'2' if r['in_cons'][1] else '-'} "
                      f"calls {','.join(map(str, r['before']))} > {','.join(map(str, r['after']))}"
                      + ("  LOST CALLS" if lost else ""))

    # EM block-cut decisions
    joins = [j for j in D.load_joins(a.diag, beg - 50_000, end + 50_000) if j["right"] >= beg and j["left"] <= end]
    if joins:
        print("\nEM joins between adjacent sites (margin, margin without weak sites, reads crossing, reads with "
              "strong evidence crossing, decision: J join, C cut, S too few reads, W weak-bridge cut; "
              "* = joined as a weak bridge):")
        for j in joins:
            weak = j["strong"] < 4.0 or (j["strong_reads"] is not None and j["strong_reads"] < 3)
            flag = "*" if j["decision"] == "J" and weak else " "
            if j["decision"] != "J" or flag == "*" or a.all:
                sr = j["strong_reads"] if j["strong_reads"] is not None else "-"
                print(f"  {j['left']:>11,} - {j['right']:<11,} margin {j['delta']:8.1f}  strong {j['strong']:8.1f}  "
                      f"reads {j['spanning']:4}  strong reads {sr:>4}  {j['decision']}{flag}")


if __name__ == "__main__":
    main()
