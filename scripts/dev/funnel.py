#!/usr/bin/env python3
"""Where the pipeline loses the reads a rival phases correctly (evaluation only).

Usage: funnel.py DIAG_DIR [--out TSV] [--top N]

For every arms_phaseable read that HiPhase or LongPhase phases correctly and
we leave unphased or discordant (official classes), each truth heterozygote
inside the read is traced through the pipeline and given the furthest stage
it reached; the read takes its best heterozygote's stage:
  em_reliable+obs      a reliable EM site this read has a call at
  em_reliable-no_obs   a reliable EM site; this read has no call there
  em_unreliable        an EM site at learned error >= 0.3
  em_unverified        a call kept from a noisy region too long for MSA
                       (phased, does not label reads)
  graph_voting         a voting graph row the EM left out (skewed alleles)
  graph_<CAT>          a graph row that does not vote (its category)
  inj_<decision>       an alignment call that was not injected (why; for
                       category: the category, and :msa when MSA-verified)
  aln_<INIT>><FINAL>   an alignment call pruned before injection
                       (CLEAN_HET_SNP>NON_VAR: suppressed by a noisy region)
  absent               neither the catalog nor the alignment holds it
Prints reads per stage for unphased and discordant reads, the loci (1 Mb
bins) holding most of each stage, and with --out one row per read.

With --rival NAME the traced sites are not the truth hets but the rival's
own evidence: its phased heterozygous records inside the read that carry the
phase set the rival tagged the read with (indels matched within 3 bp). Each
such site is also marked truth (a truth het within 1 bp) or not.

With --support N, up to N hets per stage are checked in the input BAM: at
the het (indels: within its repeat, +-10 bp), reads with MAPQ >= 20 carrying
the ALT, split by parent. A het is
  separates   ALT on at least 20% of reads and >= 90% of ALT reads from one parent
  supported   ALT on at least 20% of reads, not separating the parents
  unsupported ALT on under 20% of reads (the reads do not show it)
"""
import argparse
import bisect
import collections
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import diaglib as D  # noqa: E402
import score  # noqa: E402

ORDER = ["em_reliable+obs", "em_reliable-no_obs", "em_unreliable", "em_unverified", "graph_voting"]


def rank(stage):
    if stage in ORDER:
        return ORDER.index(stage)
    for i, p in enumerate(("graph_", "inj_", "aln_", "absent")):
        if stage.startswith(p):
            return len(ORDER) + i
    return 99


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("diag")
    ap.add_argument("--out")
    ap.add_argument("--top", type=int, default=4)
    ap.add_argument("--support", type=int, default=0)
    ap.add_argument("--rival", help="trace the rival's phased sites instead of truth hets")
    a = ap.parse_args()
    contig = os.environ["PGDEV_CONTIG"]
    comps = D.competitors()
    scope = {q: v for q, v in score.read_scope().items() if v[1] == "arms_phaseable"}
    ours = score.read_classes(a.diag)
    theirs = {n: score.read_classes(d) for n, d in comps.items()}
    target = [q for q in scope if ours[q] != "correct" and any(t[q] == "correct" for t in theirs.values())]
    print(f"reads a rival phases correctly and we do not: {len(target):,} "
          f"({sum(ours[q] == 'unphased' for q in target):,} unphased, {sum(ours[q] == 'discordant' for q in target):,} discordant)")

    truth = D.truth_records(contig, 1, 10 ** 10)
    tpos = sorted(truth)
    em_reads, em_sites = D.load_em(a.diag)
    graph, aln, inj = D.load_graph(a.diag), D.load_aln(a.diag), D.load_inj(a.diag)

    def het_stage(q, pos, tol=1):
        e = D.near(em_sites, pos, tol)
        if e:
            p, s = e[0]
            if s["kind"] == "U":
                return "em_unverified"
            if s["err"] >= 0.3:
                return "em_unreliable"
            d = em_reads.get(q)
            seen = d is not None and any(abs(o[0] - pos) <= 1 for o in d["obs"])
            return "em_reliable+obs" if seen else "em_reliable-no_obs"
        g = D.near(graph, pos, tol)
        if g:
            rows = g[0][1]
            if any(r["voting"] for r in rows):
                return "graph_voting"
            return "graph_" + rows[0]["cat"]
        i = D.near(inj, pos, tol)
        if i:
            best = min(i[0][1], key=lambda r: D.decision_rank(r["decision"]))
            if best["decision"] == "category":
                return "inj_category:" + best["cat"] + (":msa" if best["msa"] else "")
            return "inj_" + best["decision"]
        al = D.near(aln, pos, tol)
        if al:
            r = max(al[0][1], key=lambda r: r["alt_cov"])
            return f"aln_{r['init']}>{r['final']}"
        return "absent"

    per_read, het_stages = {}, {}
    rival_sites = None
    if a.rival:
        # the rival's phased hets by phase set, and each target read's phase set
        rdir = comps[a.rival]
        rvcf = D.vcf_records(os.path.join(rdir, "phased.vcf"), contig, 1, 2 ** 31 - 1)
        by_ps = collections.defaultdict(list)
        for pos, recs in rvcf.items():
            for ref, alts, gt, ps, c in recs:
                if ps and len(gt) == 2 and None not in gt and gt[0] != gt[1]:
                    snp = len(ref) == 1 and all(len(x) == 1 for x in alts)
                    by_ps[ps].append((pos, snp))
        for ps in by_ps:
            by_ps[ps].sort()
        tags = score.read_tags(rdir, set(target))
        rival_sites = {}
        for q in target:
            t = tags.get(q)
            if t is None:
                continue
            lst = by_ps.get(str(t[1]), [])
            s, e = scope[q][3], scope[q][4]
            i = bisect.bisect_right(lst, (s, True))
            rival_sites[q] = [x for x in lst[i:] if x[0] <= e]
    for q in target:
        s, e = scope[q][3], scope[q][4]
        stages = []
        if rival_sites is not None:
            for p, snp in rival_sites.get(q, []):
                st = het_stage(q, p, 1 if snp else 3)
                stages.append((p, st))
                het_stages.setdefault(p, st if st != "em_reliable+obs" else "em_reliable")
        else:
            hets = tpos[bisect.bisect_right(tpos, s):bisect.bisect_right(tpos, e)]
            for p in hets:
                st = het_stage(q, p)
                stages.append((p, st))
                het_stages.setdefault(p, st if st != "em_reliable+obs" else "em_reliable")
        best = min((st for _, st in stages), key=rank) if stages else "no_site"
        per_read[q] = (best, stages)
    if rival_sites is not None:
        kinds = {}
        for q, lst in rival_sites.items():
            for p, snp in lst:
                kinds[p] = (snp, any(x in truth for x in (p - 1, p, p + 1)))
        rival_used = collections.Counter()
        for p, st in het_stages.items():
            snp, tr = kinds.get(p, (None, None))
            rival_used[(st, "snp" if snp else "indel", "truth" if tr else "not_truth")] += 1
        print(f"\n{a.rival}: phased sites it used for these reads, by our furthest stage (each site once):")
        for (st, typ, tr), n in sorted(rival_used.items(), key=lambda kv: (rank(kv[0][0]), -kv[1])):
            print(f"  {st:34} {typ:5} {tr:9} {n:6,}")

    for klass in ("unphased", "discordant"):
        reads = [q for q in target if ours[q] == klass]
        tab = collections.Counter(per_read[q][0] for q in reads)
        print(f"\n{klass} ({len(reads):,}): best stage of the read's truth hets")
        for st, n in sorted(tab.items(), key=lambda kv: (rank(kv[0]), -kv[1])):
            bins = collections.Counter((scope[q][3] + scope[q][4]) // 2 // 1_000_000 for q in reads if per_read[q][0] == st)
            loci = ", ".join(f"{b} Mb:{c}" for b, c in bins.most_common(a.top))
            rivals = ", ".join(f"{name} {sum(theirs[name][q] == 'correct' for q in reads if per_read[q][0] == st):,}"
                               for name in comps)
            print(f"  {st:34} {n:6,}   ({rivals})   {loci}")
    lost = collections.Counter(st for st in het_stages.values())
    print("\n" + (f"{a.rival} sites" if a.rival else "truth hets") + " on those reads, by furthest stage reached (each once):")
    for st, n in sorted(lost.items(), key=lambda kv: (rank(kv[0]), -kv[1])):
        print(f"  {st:34} {n:6,}")
    if a.support:
        import pysam
        bam = pysam.AlignmentFile(os.environ["PGDEV_BAM"])
        parent = score.truth_tables()["reads"]
        by_stage = collections.defaultdict(list)
        for p, st in sorted(het_stages.items()):
            by_stage[st].append(p)
        print(f"\nread support of those hets (up to {a.support} per stage; MAPQ >= 20):")
        for st in sorted(by_stage, key=rank):
            hets = by_stage[st][:: max(1, len(by_stage[st]) // a.support)][: a.support]
            tab = collections.Counter()
            for p in hets:
                if p not in truth:  # a rival's site that is not a truth het
                    tab["not_truth"] += 1
                    continue
                ref, alt, gt = truth[p]
                alt = alt.split(",")[0]
                snp = len(ref) == 1 and len(alt) == 1
                carry = collections.Counter()
                total = 0
                for r in bam.fetch(contig, p - 1, p):
                    if r.is_secondary or r.is_supplementary or r.mapping_quality < 20 or r.query_name not in parent:
                        continue
                    total += 1
                    has = False
                    if snp:
                        for qp, rp in r.get_aligned_pairs(matches_only=True):
                            if rp == p - 1:
                                has = r.query_sequence[qp].upper() == alt.upper()
                                break
                    else:
                        want_del = len(ref) > len(alt)
                        rpos = r.reference_start
                        for op, ln in r.cigartuples:
                            if op in (0, 7, 8):
                                rpos += ln
                            elif op == 2:
                                if want_del and abs(rpos - p) <= 10:
                                    has = True
                                rpos += ln
                            elif op == 1:
                                if not want_del and abs(rpos - p) <= 10:
                                    has = True
                            if rpos > p + 12:
                                break
                    if has:
                        carry["M" if parent[r.query_name][2] else "P"] += 1
                n_alt = carry["M"] + carry["P"]
                if total == 0 or n_alt < 0.2 * total:
                    tab["unsupported"] += 1
                elif max(carry["M"], carry["P"]) >= 0.9 * n_alt:
                    tab["separates"] += 1
                else:
                    tab["supported"] += 1
            print(f"  {st:34} checked {len(hets):4}: " + ", ".join(f"{k} {v}" for k, v in tab.most_common()))
    if a.out:
        with open(a.out, "w") as fh:
            fh.write("read\tstart\tend\tours\t" + "\t".join(comps) + "\tbest\thets\n")
            for q, (best, stages) in sorted(per_read.items(), key=lambda kv: scope[kv[0]][3]):
                fh.write(f"{q}\t{scope[q][3]}\t{scope[q][4]}\t{ours[q]}\t" +
                         "\t".join(theirs[n][q] for n in comps) + f"\t{best}\t" +
                         ",".join(f"{p}:{st}" for p, st in stages) + "\n")
        print(f"\nper-read table: {a.out}")


if __name__ == "__main__":
    main()
