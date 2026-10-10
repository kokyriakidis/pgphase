#!/usr/bin/env python3
"""Score a phasing run against parental truth (evaluation only).

  score.py run   DIR [--regions BED] [--json OUT]   score one run directory
  score.py table JSON [JSON ...] [--names a,b,...]  side-by-side table, first is the reference
  score.py bins  DIR [--bin BP]                     per-bin phaseable-read profile of a run
  score.py diff  BASE DIR --out TSV [--regions BED] arms_phaseable reads whose class changed

DIR holds phased.bam (HP/PS tags) and phased.vcf or phased.vcf.gz.

Reads: truth-labelled primary alignments of PGDEV_BAM. Each phase set takes
the majority parental orientation over its scored reads; a read is correct
when it agrees with its phase set's orientation, discordant when it does not,
unphased when it has no HP/PS. A read is centromeric ("cen") when its
alignment midpoint lies in PGDEV_CEN, otherwise on the arms; it is
"phaseable" when its aligned span contains a heterozygous truth record,
otherwise "nohet". With --regions only reads whose midpoint lies in a region
are scored, and only VCF records inside a region count.

Variants: heterozygous biallelic records matched to the phased truth VCF by
(POS, REF, ALT); within a phase set, a single disagreeing site between
agreeing neighbours is a flip, a lasting change a switch.

Truth tables are cached under PGDEV_CACHE keyed by input file identity.
"""
import argparse
import bisect
import collections
import gzip
import hashlib
import json
import os
import pickle
import sys

import pysam

ENV = os.environ
CEN = tuple(int(x) for x in ENV.get("PGDEV_CEN", "26000000-32000000").split("-"))
CACHE = ENV.get("PGDEV_CACHE", ".dev-cache")


def file_id(path):
    st = os.stat(path)
    return f"{os.path.realpath(path)}:{st.st_size}:{int(st.st_mtime)}"


def truth_tables():
    """name -> (start, end, parent); sorted truth het positions; truth phased hets."""
    bam, reads_tsv, vcf = ENV["PGDEV_BAM"], ENV["PGDEV_TRUTH_READS"], ENV["PGDEV_TRUTH_VCF"]
    key = hashlib.sha256("|".join(map(file_id, (bam, reads_tsv, vcf))).encode()).hexdigest()[:16]
    path = os.path.join(CACHE, "truth", f"{key}.pkl")
    if os.path.exists(path):
        with open(path, "rb") as fh:
            return pickle.load(fh)
    parent = {}
    with open(reads_tsv) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) > 1 and f[1] in ("MATERNAL", "PATERNAL"):
                parent[f[0]] = f[1] == "MATERNAL"
    reads = {}
    with pysam.AlignmentFile(bam) as fh:
        for r in fh.fetch(until_eof=True):
            if r.is_unmapped or r.is_secondary or r.is_supplementary or r.is_duplicate:
                continue
            if r.query_name in parent:
                reads[r.query_name] = (r.reference_start, r.reference_end, parent[r.query_name])
    het_pos, phased = [], {}
    for rec in vcf_records(vcf):
        pos, ref, alts, gt, _, _ = rec
        if len(gt) == 2 and gt[0] != gt[1] and None not in gt:
            het_pos.append(pos - 1)
            if len(alts) == 1 and set(gt) == {0, 1}:
                phased[(pos, ref.upper(), alts[0].upper())] = gt[0]
    het_pos.sort()
    tables = {"reads": reads, "het_pos": het_pos, "truth_phased": phased}
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path + ".tmp", "wb") as fh:
        pickle.dump(tables, fh, protocol=pickle.HIGHEST_PROTOCOL)
    os.replace(path + ".tmp", path)
    return tables


def vcf_records(path):
    """(pos, ref, alts, gt tuple, PS or None, CAT) per record."""
    opener = gzip.open if path.endswith(".gz") else open
    with opener(path, "rt") as fh:
        for line in fh:
            if line[0] == "#":
                continue
            f = line.rstrip("\n").split("\t")
            fmt = dict(zip(f[8].split(":"), f[9].split(":")))
            raw = fmt.get("GT", ".")
            sep = "|" if "|" in raw else "/"
            gt = tuple(None if a == "." else int(a) for a in raw.split(sep))
            ps = fmt.get("PS") if sep == "|" else None
            cat = next((kv[4:] for kv in f[7].split(";") if kv.startswith("CAT=")), "?")
            yield int(f[1]), f[3], f[4].split(","), gt, (None if ps in (None, ".", "") else ps), cat


def load_regions(path):
    if not path:
        return None
    out = []
    with open(path) as fh:
        for line in fh:
            if line.strip() and line[0] != "#":
                f = line.split()
                out.append((int(f[1]), int(f[2])))
    return sorted(out)


def in_regions(regions, pos):
    if regions is None:
        return True
    i = bisect.bisect_right(regions, (pos, float("inf"))) - 1
    return i >= 0 and regions[i][0] <= pos < regions[i][1]


def read_tags(run_dir, names):
    """qname -> (HP, PS) for scored names. Aligned BAMs are fetched whole."""
    tags = {}
    with pysam.AlignmentFile(os.path.join(run_dir, "phased.bam"), check_sq=False) as fh:
        for r in fh.fetch(until_eof=True):
            if r.is_secondary or r.is_supplementary or r.query_name not in names:
                continue
            if r.has_tag("HP") and r.has_tag("PS"):
                tags[r.query_name] = (r.get_tag("HP"), r.get_tag("PS"))
    return tags


def n50(lengths):
    total, run = sum(lengths), 0
    for x in sorted(lengths, reverse=True):
        run += x
        if 2 * run >= total:
            return x
    return 0


def read_scope(regions=None):
    """Scored reads: name -> (arm, arm_sub, maternal, start, end)."""
    t = truth_tables()
    het_pos = t["het_pos"]
    scope = {}
    for q, (s, e, mat) in t["reads"].items():
        mid = (s + e) // 2
        if not in_regions(regions, mid):
            continue
        arm = "cen" if CEN[0] <= mid < CEN[1] else "arms"
        crosses = bisect.bisect_left(het_pos, e) > bisect.bisect_left(het_pos, s)
        scope[q] = (arm, f"{arm}_{'phaseable' if crosses else 'nohet'}", mat, s, e)
    return scope


def classify_reads(run_dir, scope):
    """name -> correct / discordant / unphased under per-phase-set majority orientation."""
    tags = read_tags(run_dir, scope)
    votes = collections.defaultdict(lambda: [0, 0])
    for q, (hp, ps) in tags.items():
        votes[ps][((hp == 1) == scope[q][2])] += 1
    orient = {ps: v[1] >= v[0] for ps, v in votes.items()}
    out = {}
    for q, info in scope.items():
        tag = tags.get(q)
        out[q] = "unphased" if tag is None else (
            "correct" if ((tag[0] == 1) == info[2]) == orient[tag[1]] else "discordant")
    return out


def score_run(run_dir, regions=None):
    t = truth_tables()
    scope = read_scope(regions)
    classes = classify_reads(run_dir, scope)
    groups = ("all", "arms", "cen", "arms_phaseable", "arms_nohet", "cen_phaseable", "cen_nohet")
    reads_out = {g: {"correct": 0, "discordant": 0, "unphased": 0} for g in groups}
    for q, (arm, sub, *_) in scope.items():
        for g in ("all", arm, sub):
            reads_out[g][classes[q]] += 1

    vcf = os.path.join(run_dir, "phased.vcf")
    if not os.path.exists(vcf):
        vcf += ".gz"
    truth_phased = t["truth_phased"]
    blocks, by_ps = {}, collections.defaultdict(list)
    for pos, ref, alts, gt, ps, cat in vcf_records(vcf):
        if ps is None or len(gt) != 2 or gt[0] == gt[1] or None in gt or not in_regions(regions, pos - 1):
            continue
        lo, hi = blocks.get(ps, (pos, pos))
        blocks[ps] = (min(lo, pos), max(hi, pos))
        if len(alts) == 1 and set(gt) == {0, 1}:
            key = (pos, ref.upper(), alts[0].upper())
            if key in truth_phased:
                by_ps[ps].append((pos, gt[0] == truth_phased[key]))
    var = collections.Counter()
    for ps, sites in by_ps.items():
        sites.sort()
        var["truth_matched_phased_hets"] += len(sites)
        same = [x[1] for x in sites]
        i = 1
        while i < len(same):
            if same[i] != same[i - 1]:
                where = "cen" if CEN[0] <= sites[i][0] - 1 < CEN[1] else "arms"
                if i + 1 < len(same) and same[i + 1] == same[i - 1]:
                    var[f"flips_{where}"] += 1
                    i += 2
                    continue
                var[f"switches_{where}"] += 1
            i += 1
    spans = [hi - lo for lo, hi in blocks.values()]
    return {"run": os.path.abspath(run_dir), "reads": reads_out,
            "variants": {k: var.get(k, 0) for k in ("truth_matched_phased_hets", "switches_arms",
                                                     "switches_cen", "flips_arms", "flips_cen")},
            "blocks": {"phase_sets": len(blocks), "n50": n50(spans), "largest": max(spans, default=0)}}


ROWS = [  # (label, getter, higher is better)
    ("arms_phaseable correct", lambda s: s["reads"]["arms_phaseable"]["correct"], True),
    ("arms_phaseable discordant", lambda s: s["reads"]["arms_phaseable"]["discordant"], False),
    ("arms_phaseable unphased", lambda s: s["reads"]["arms_phaseable"]["unphased"], False),
    ("arms_phaseable correct %", lambda s: pct(s["reads"]["arms_phaseable"]), True),
    ("arms_nohet discordant", lambda s: s["reads"]["arms_nohet"]["discordant"], False),
    ("cen_phaseable correct", lambda s: s["reads"]["cen_phaseable"]["correct"], True),
    ("cen_phaseable discordant", lambda s: s["reads"]["cen_phaseable"]["discordant"], False),
    ("all correct", lambda s: s["reads"]["all"]["correct"], True),
    ("all discordant", lambda s: s["reads"]["all"]["discordant"], False),
    ("switches arms", lambda s: s["variants"]["switches_arms"], False),
    ("switches cen", lambda s: s["variants"]["switches_cen"], False),
    ("flips arms", lambda s: s["variants"]["flips_arms"], False),
    ("truth hets phased", lambda s: s["variants"]["truth_matched_phased_hets"], True),
    ("phase sets", lambda s: s["blocks"]["phase_sets"], False),
    ("block N50", lambda s: s["blocks"]["n50"], True),
]


def pct(r):
    n = r["correct"] + r["discordant"] + r["unphased"]
    return round(100.0 * r["correct"] / n, 3) if n else 0.0


def table(scores, names):
    """Rows of metrics; columns after the first show their difference from it,
    marked + when better and - when worse."""
    grid = []
    for label, get, up in ROWS:
        ref = get(scores[0])
        cells = []
        for i, s in enumerate(scores):
            v = get(s)
            cell = f"{v:,.3f}" if isinstance(v, float) else f"{v:,}"
            if i > 0 and v != ref:
                d = v - ref
                mark = "+" if (d > 0) == up else "-"
                cell += f" ({d:+,.3f}{mark})" if isinstance(v, float) else f" ({d:+,}{mark})"
            cells.append(cell)
        grid.append((label, cells))
    widths = [max(len(names[i]), *(len(c[i]) for _, c in grid)) + 2 for i in range(len(names))]
    out = ["%-27s" % "" + "".join("%*s" % (w, n) for w, n in zip(widths, names))]
    out += ["%-27s" % label + "".join("%*s" % (w, c) for w, c in zip(widths, cells)) for label, cells in grid]
    return "\n".join(out)


def diff(base_dir, run_dir, regions, out_path, bin_size=1_000_000):
    """Per-read class changes from base_dir to run_dir on arms_phaseable reads."""
    scope = {q: v for q, v in read_scope(regions).items() if v[1] == "arms_phaseable"}
    a, b = classify_reads(base_dir, scope), classify_reads(run_dir, scope)
    moves = collections.Counter((a[q], b[q]) for q in scope if a[q] != b[q])
    print("arms_phaseable reads that changed class (base -> run):")
    for (x, y), n in sorted(moves.items(), key=lambda kv: -kv[1]):
        good = {"correct": 2, "unphased": 1, "discordant": 0}
        print(f"  {x:>10} -> {y:<10} {n:6,}  {'better' if good[y] > good[x] else 'worse'}")
    if not moves:
        print("  none")
        return
    per_bin = collections.defaultdict(collections.Counter)
    with open(out_path, "w") as fh:
        fh.write("read\tstart\tend\tbase\trun\n")
        for q in sorted(scope, key=lambda q: scope[q][3]):
            if a[q] != b[q]:
                s, e = scope[q][3], scope[q][4]
                fh.write(f"{q}\t{s}\t{e}\t{a[q]}\t{b[q]}\n")
                per_bin[(s + e) // 2 // bin_size][f"{a[q]}->{b[q]}"] += 1
    print(f"most changed {bin_size // 1_000_000} Mb bins:")
    for k, c in sorted(per_bin.items(), key=lambda kv: -sum(kv[1].values()))[:8]:
        print(f"  {k * bin_size:>11,}  " + "  ".join(f"{m} {n}" for m, n in c.most_common()))
    print(f"changed reads: {out_path}")


def bins(run_dir, size):
    scope = read_scope()
    classes = classify_reads(run_dir, scope)
    out = collections.defaultdict(collections.Counter)
    for q, (_, sub, _, s, e) in scope.items():
        if sub.endswith("phaseable"):
            out[((s + e) // 2) // size][classes[q]] += 1
    print("bin_start\tbin_end\tcorrect\tdiscordant\tunphased")
    for b in sorted(out):
        c = out[b]
        print(f"{b * size}\t{(b + 1) * size}\t{c['correct']}\t{c['discordant']}\t{c['unphased']}")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    r = sub.add_parser("run")
    r.add_argument("dir")
    r.add_argument("--regions")
    r.add_argument("--json")
    t = sub.add_parser("table")
    t.add_argument("json", nargs="+")
    t.add_argument("--names")
    b = sub.add_parser("bins")
    b.add_argument("dir")
    b.add_argument("--bin", type=int, default=1_000_000)
    d = sub.add_parser("diff")
    d.add_argument("base")
    d.add_argument("dir")
    d.add_argument("--regions")
    d.add_argument("--out", required=True)
    a = ap.parse_args()
    if a.cmd == "run":
        s = score_run(a.dir, load_regions(a.regions))
        text = json.dumps(s, indent=1)
        if a.json:
            with open(a.json, "w") as fh:
                fh.write(text + "\n")
        else:
            print(text)
    elif a.cmd == "table":
        scores = [json.load(open(p)) for p in a.json]
        names = a.names.split(",") if a.names else [os.path.basename(os.path.dirname(p)) for p in a.json]
        print(table(scores, names))
    elif a.cmd == "diff":
        diff(a.base, a.dir, load_regions(a.regions), a.out)
    else:
        bins(a.dir, a.bin)


if __name__ == "__main__":
    sys.exit(main())
