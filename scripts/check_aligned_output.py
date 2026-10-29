#!/usr/bin/env python3
"""collect-graph-variation's aligned output (--out-sam), with --refine-aln against without (stdlib only).

  check_aligned_output.py --refined REFINED.sam --plain PLAIN.sam [--ref REF.fa[.gz]]

Both are sorted and hold the same records (by name, flag and, unmapped, position) with the same SEQ, QUAL and HP/PS;
the refined one rewrites some primaries' alignments and no other record's, each with a CIGAR as long as its SEQ and,
with --ref, NM and MD (where present) those of the new alignment. Exit status 1 on any failure.
"""
import argparse
import collections
import gzip
import re
import sys

CIGAR = re.compile(r"(\d+)([MIDNSHP=X])")


def read_sam(path):
    order, records = {}, []
    for line in open(path):
        if line.startswith("@"):
            if line.startswith("@SQ"):
                name = re.search(r"\tSN:(\S+)", line).group(1)
                order[name] = len(order)
            continue
        c = line.rstrip("\n").split("\t")
        tags = {t[:2]: t[5:] for t in c[11:]}
        records.append({"qname": c[0], "flag": int(c[1]), "rname": c[2], "pos": int(c[3]), "cigar": c[5],
                        "seq": c[9], "qual": c[10], "tags": tags})
    return order, records


def key(r):
    return (r["qname"], r["flag"], r["pos"] if r["flag"] & 4 else 0)


def ops(cigar):
    return [(int(n), op) for n, op in CIGAR.findall(cigar)] if cigar != "*" else []


def query_len(cigar):
    return sum(n for n, op in ops(cigar) if op in "MIS=X")


def read_fasta(path):
    seqs, name, parts = {}, None, []
    with (gzip.open(path, "rt") if path.endswith(".gz") else open(path)) as f:
        for line in f:
            if line.startswith(">"):
                if name is not None:
                    seqs[name] = "".join(parts).upper()
                name, parts = line[1:].split()[0], []
            else:
                parts.append(line.strip())
    if name is not None:
        seqs[name] = "".join(parts).upper()
    return seqs


def nm_md(r, ref):
    rp, qp, nm, run, md = r["pos"] - 1, 0, 0, 0, ""
    seq = r["seq"].upper()
    for n, op in ops(r["cigar"]):
        if op in "M=X":
            for j in range(n):
                if ref[rp + j] == seq[qp + j]:
                    run += 1
                else:
                    nm, md, run = nm + 1, md + f"{run}{ref[rp + j]}", 0
            rp, qp = rp + n, qp + n
        elif op == "I":
            nm, qp = nm + n, qp + n
        elif op == "D":
            nm, md, run, rp = nm + n, md + f"{run}^{ref[rp:rp + n]}", 0, rp + n
        elif op == "N":
            rp += n
        elif op == "S":
            qp += n
    return nm, md + str(run)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--refined", required=True)
    ap.add_argument("--plain", required=True)
    ap.add_argument("--ref")
    a = ap.parse_args()
    fail = []
    sams = {}
    for name, path in (("refined", a.refined), ("plain", a.plain)):
        order, recs = read_sam(path)
        sams[name] = recs
        last = (-1, -1)
        for r in recs:
            at = (order.get(r["rname"], len(order)), r["pos"])
            if r["rname"] != "*" and at < last:
                fail.append(f"{name}: {r['qname']} out of order")
                break
            last = at
        dups = [k for k, n in collections.Counter(key(r) for r in recs).items() if n > 1]
        if dups:
            fail.append(f"{name}: {len(dups)} records written twice, e.g. {dups[0]}")
    refined = {key(r): r for r in sams["refined"]}
    plain = {key(r): r for r in sams["plain"]}
    if refined.keys() != plain.keys():
        fail.append(f"records differ: {len(refined.keys() - plain.keys())} only refined, "
                    f"{len(plain.keys() - refined.keys())} only plain")
    ref = read_fasta(a.ref) if a.ref else None
    rewritten = tagged = 0
    for k in refined.keys() & plain.keys():
        r, p = refined[k], plain[k]
        if r["seq"] != p["seq"] or r["qual"] != p["qual"]:
            fail.append(f"{r['qname']}: SEQ or QUAL changed")
        if (r["tags"].get("HP"), r["tags"].get("PS")) != (p["tags"].get("HP"), p["tags"].get("PS")):
            fail.append(f"{r['qname']}: HP/PS differ with --refine-aln")
        tagged += "HP" in r["tags"]
        if (r["pos"], r["cigar"]) == (p["pos"], p["cigar"]):
            continue
        rewritten += 1
        if r["flag"] & (0x100 | 0x800 | 0x4):
            fail.append(f"{r['qname']}: a non-primary record rewritten")
        if r["seq"] != "*" and query_len(r["cigar"]) != len(r["seq"]):
            fail.append(f"{r['qname']}: CIGAR {query_len(r['cigar'])} long, SEQ {len(r['seq'])}")
        if ref is not None and ("NM" in r["tags"] or "MD" in r["tags"]):
            nm, md = nm_md(r, ref[r["rname"]])
            if ("NM" in r["tags"] and int(r["tags"]["NM"]) != nm) or ("MD" in r["tags"] and r["tags"]["MD"] != md):
                fail.append(f"{r['qname']}: NM/MD not the new alignment's")
    if rewritten == 0:
        fail.append("no alignment refined")
    for f in fail[:20]:
        print(f"  {f}", file=sys.stderr)
    print(f"aligned output: {len(refined)} records, {tagged} with HP/PS, {rewritten} refined"
          + ("" if not fail else f"; {len(fail)} failures"))
    sys.exit(1 if fail else 0)


if __name__ == "__main__":
    main()
