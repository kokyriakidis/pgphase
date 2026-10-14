#!/usr/bin/env python3
"""Would whole-read graph genotyping label the reads we miss? (evaluation only)

Usage: graphgeno.py RUN_DIR [--classes unphased,discordant] [--control N]
                    [--workers W] [--min-margin M] [--seed S] [--sites truth|rival]

Each target read (arms_phaseable, our class in --classes; plus --control N
reads we phase correctly) is aligned whole to the graph of every het and
hom-alt call in our own phased VCF that it spans (hiphase_port.py: HiPhase's
global realignment, local fallback). Its allele at every phased het votes for
the haplotype of that het's phase set carrying it; the read is labelled in
the phase set with the largest vote margin when the margin is at least
--min-margin (votes for one side minus the other). A label is right when it
puts the read on its parent's side under the phase set's orientation (the
majority of our labelled reads in it), as score.py orients blocks.
Tables: labels by margin bucket and right / wrong, per class; the largest
loci of right and wrong labels.
--sites truth (every truth het and hom, phased as in the truth) or rival
(the DeepVariant calls HiPhase phases, each het that matches a truth het
phased as in the truth) replaces our VCF: a ceiling for reads that carry
phaseable sequence at sites we do not call. One phase set; its orientation
comes from the control reads' labels.
"""
import argparse
import bisect
import collections
import multiprocessing
import os
import random
import sys

import pysam

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import diaglib as D  # noqa: E402
import hiphase_port as H  # noqa: E402
import score  # noqa: E402

G = {}


def work(job):
    q, s, e = job
    hets, homs, seq, starts = G["hets"], G["homs"], G["seq"], G["starts"]
    bam = G.get("bam_fh")
    if bam is None:
        bam = G["bam_fh"] = pysam.AlignmentFile(G["bam"])
    rec = None
    for r in bam.fetch(G["contig"], s, s + 1):
        if r.query_name == q and not r.is_secondary and not r.is_supplementary:
            rec = r
            break
    if rec is None:
        return q, None, "no read"
    lo = bisect.bisect_left(starts, rec.reference_start - 20)
    hi = bisect.bisect_right(starts, rec.reference_end + 20)
    if hi <= lo:
        return q, {}, "no het"
    sub = hets[lo:hi]
    hpos = [v.position for v in homs]
    hsub = homs[bisect.bisect_left(hpos, rec.reference_start - 20):bisect.bisect_right(hpos, rec.reference_end + 20)]
    mode = "global"
    try:
        alleles, _ = H.global_realignment(rec, sub, hsub, seq, 0)
    except H.MaxEditDistance:
        alleles = H.local_realignment(rec, sub)
        mode = "local"
    votes = collections.defaultdict(lambda: [0, 0])
    for v, a in zip(sub, alleles):
        ph = G["phase"].get(v.key)
        if ph is None or a not in (H.REF, H.ALT):
            continue
        ps, g0, lo_al, hi_al = ph
        carried = lo_al if a == H.REF else hi_al
        votes[ps][0 if carried == g0 else 1] += 1
    return q, dict(votes), mode


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("run")
    ap.add_argument("--classes", default="unphased,discordant")
    ap.add_argument("--control", type=int, default=300)
    ap.add_argument("--workers", type=int, default=12)
    ap.add_argument("--min-margin", type=int, default=2)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--sites", choices=["ours", "truth", "rival"], default="ours")
    a = ap.parse_args()
    contig = os.environ["PGDEV_CONTIG"]
    fa = pysam.FastaFile(os.environ["PGDEV_REF"])
    ref_name = contig if contig in fa.references else contig.split("#")[-1]
    bam_path = os.environ["PGDEV_BAM"]
    with pysam.AlignmentFile(bam_path) as fh:
        bam_contig = contig if contig in fh.references else contig.split("#")[-1]
    seq = fa.fetch(ref_name).upper()
    records, phase_of = [], {}
    if a.sites == "ours":
        for pos, recs in sorted(D.vcf_records(os.path.join(a.run, "phased.vcf"), contig, 1, 2 ** 31 - 1).items()):
            for ref, alts, gt, ps, c in recs:
                if len(gt) != 2 or None in gt:
                    continue
                records.append((pos - 1, ref.upper(), [x.upper() for x in alts], gt))
                if ps is not None and gt[0] != gt[1]:
                    phase_of[(pos - 1, ref.upper(), tuple(x.upper() for x in alts), gt)] = (ps, gt[0], min(gt), max(gt))
    else:
        truth_gt = {}
        for pos, recs in D.vcf_records(os.environ["PGDEV_TRUTH_VCF"], contig, 1, 2 ** 31 - 1).items():
            for ref, alts, gt, ps, c in recs:
                if len(gt) == 2 and None not in gt:
                    truth_gt[(pos, ref.upper(), tuple(x.upper() for x in alts))] = gt
        if a.sites == "truth":
            src = os.environ["PGDEV_TRUTH_VCF"]
        else:
            from hiphase_exact import hiphase_inputs
            src = hiphase_inputs()[2]
        for pos, recs in sorted(D.vcf_records(src, contig, 1, 2 ** 31 - 1).items()):
            for ref, alts, gt, ps, c in recs:
                if len(gt) != 2 or None in gt:
                    continue
                alts_u = tuple(x.upper() for x in alts)
                tg = truth_gt.get((pos, ref.upper(), alts_u))
                if gt[0] != gt[1] and tg is not None and tg[0] != tg[1]:
                    gt = tg  # truth phase (and alleles) for a matching het
                records.append((pos - 1, ref.upper(), list(alts_u), gt))
                if gt[0] != gt[1] and tg is not None and tg[0] != tg[1]:
                    phase_of[(pos - 1, ref.upper(), alts_u, gt)] = ("1", gt[0], min(gt), max(gt))
    hets, homs = H.load_variants(records, lambda x, y: seq[x:y], reference_buffer=15)
    # load_variants keys use the record as given
    phase = {}
    for v in hets:
        k = (v.key[0], v.key[1].upper(), tuple(x.upper() for x in v.key[2]), v.key[3])
        if k in phase_of:
            phase[v.key] = phase_of[k]
    scope = score.read_scope()
    classes = score.read_classes(a.run)
    orient = {}
    if a.sites == "ours":
        tags = score.read_tags(a.run, set(scope))
        votes = collections.defaultdict(lambda: [0, 0])
        for q, (hp, ps) in tags.items():
            votes[str(ps)][((hp == 1) == scope[q][2])] += 1
        orient = {ps: v[1] >= v[0] for ps, v in votes.items()}
    want = set(a.classes.split(","))
    arms = [q for q, v in scope.items() if v[1] == "arms_phaseable"]
    targets = [q for q in arms if classes[q] in want]
    rng = random.Random(a.seed)
    correct = [q for q in arms if classes[q] == "correct"]
    control = rng.sample(correct, min(a.control, len(correct)))
    jobs = [(q, scope[q][3], scope[q][4]) for q in targets + control]
    G.update(hets=hets, homs=homs, seq=seq, starts=[v.position for v in hets], bam=bam_path, contig=bam_contig,
             phase=phase)
    with multiprocessing.get_context("fork").Pool(a.workers) as pool:
        results = list(pool.imap_unordered(work, jobs, chunksize=4))
    if a.sites != "ours":
        # one phase set: orient it by the control reads' graph labels
        ov = [0, 0]
        for q, v, mode in results:
            if q in set(control) and v and "1" in v and abs(v["1"][0] - v["1"][1]) >= a.min_margin:
                ov[((v["1"][0] > v["1"][1]) == scope[q][2])] += 1
        orient = {"1": ov[1] >= ov[0]}
    tab = collections.defaultdict(collections.Counter)
    loci = collections.defaultdict(collections.Counter)
    modes = collections.Counter()
    ctrl = set(control)
    for q, v, mode in results:
        group = "control (our correct)" if q in ctrl else classes[q]
        modes[(group, mode)] += 1
        if not v:
            tab[group]["no phased het called"] += 1
            continue
        ps, (n1, n2) = max(v.items(), key=lambda kv: abs(kv[1][0] - kv[1][1]))
        m = abs(n1 - n2)
        if m < a.min_margin:
            tab[group][f"margin < {a.min_margin}"] += 1
            continue
        hap = 1 if n1 > n2 else 2
        if str(ps) not in orient:
            tab[group]["phase set has no labelled reads"] += 1
            continue
        right = ((hap == 1) == scope[q][2]) == orient[str(ps)]
        mb = "2-3" if m < 4 else "4-9" if m < 10 else ">=10"
        clean = "unanimous" if min(n1, n2) == 0 else "split"
        tab[group][f"{'right' if right else 'wrong'} margin {mb} {clean}"] += 1
        loci[(group, "right" if right else "wrong")][scope[q][3] // 100000] += 1
    print(f"{contig}: graph genotyping on {a.sites} sites ({len(hets):,} hets, {len(homs):,} homs; "
          f"{len(phase):,} phased); targets {len(targets):,}, control {len(control):,}")
    for group in sorted(tab):
        t = tab[group]
        right = sum(v for k, v in t.items() if k.startswith("right"))
        wrong = sum(v for k, v in t.items() if k.startswith("wrong"))
        n = sum(t.values())
        print(f"  {group}: {n:,} reads; labelled right {right:,}, wrong {wrong:,}; "
              f"modes " + ", ".join(f"{m} {c}" for (g, m), c in modes.items() if g == group))
        for k in sorted(t):
            print(f"    {k:40} {t[k]:6,}")
    for (group, verdict), c in sorted(loci.items()):
        if group.startswith("control"):
            continue
        print(f"  {group} {verdict}: " + ", ".join(f"{b / 10:.1f} Mb:{n}" for b, n in c.most_common(8)))


if __name__ == "__main__":
    main()
