#!/usr/bin/env python3
"""The stage trace's own contract: every line well formed and every reason documented.

  check_trace.py --trace DIR_OR_TRACE_TSV --doc docs/TRACE.md

Checks DIR/trace.tsv (and DIR/reads.tsv when the run wrote one, with --phased-bam-out) against docs/TRACE.md:

  header     the columns docs/TRACE.md names, in order
  columns    every line has them all; chunk and pass are integers, pass 0, 1 or 2; beg and end integers, end >= beg
  reasons    every (stage, event, reason) appears in a row of docs/TRACE.md's tables
  detail     the detail is space-separated key=value tokens or bare flags (lower-case words), nothing else
  reads      reads.tsv: its columns, hap 0, 1 or 2

so a new trace point cannot land undocumented, and the harness's parser (core.parse_detail) reads every detail.
Python 3 standard library only.
"""
import argparse
import collections
import os
import re
import sys

TRACE_COLUMNS = ["contig", "chunk", "pass", "stage", "event", "beg", "end", "key", "reason", "detail"]
READS_COLUMNS = ["contig", "beg", "end", "read", "hap", "ps", "copies", "chunks"]
TOKEN = re.compile(r"^(?:[a-z][a-z0-9_]*=\S*|[a-z][a-z0-9_]*)$")


def documented(doc_path):
    """(stage, event, reason) triples of docs/TRACE.md's table rows: `| stage | event(s) | reasons | ...`; a row's
    events (`a / b`) and backquoted reasons combine every way."""
    out = set()
    for line in open(doc_path):
        cells = [c.strip() for c in line.strip().split("|")]
        if len(cells) < 5 or cells[1] in ("stage", "") or set(cells[1]) <= set("-"):
            continue
        stage = cells[1]
        events = [e.strip() for e in cells[2].split("/")]
        reasons = re.findall(r"`([^`]+)`", cells[3])
        for ev in events:
            for r in reasons:
                out.add((stage, ev, r))
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--trace", required=True)
    ap.add_argument("--doc", required=True)
    ap.add_argument("--examples", type=int, default=3)
    a = ap.parse_args()
    trace = os.path.join(a.trace, "trace.tsv") if os.path.isdir(a.trace) else a.trace
    reads = os.path.join(os.path.dirname(trace), "reads.tsv")
    doc = documented(a.doc)
    bad = collections.defaultdict(list)
    n = 0
    with open(trace) as f:
        header = f.readline().rstrip("\n").split("\t")
        if header != TRACE_COLUMNS:
            bad["header"].append(" ".join(header))
        for line in f:
            n += 1
            c = line.rstrip("\n").split("\t")
            if len(c) != len(TRACE_COLUMNS):
                bad["columns"].append(line.rstrip("\n")[:160])
                continue
            contig, chunk, pass_, stage, event, beg, end, key, reason, detail = c
            try:
                ok = int(pass_) in (0, 1, 2) and int(chunk) >= -1 and int(end) >= int(beg) and contig != ""
            except ValueError:
                ok = False
            if not ok:
                bad["columns"].append(line.rstrip("\n")[:160])
            if (stage, event, reason) not in doc:
                bad["reasons"].append(f"{stage} {event} {reason}")
            wrong = [t for t in detail.split(" ") if t and not TOKEN.match(t)]
            if wrong or "  " in detail or detail != detail.strip():
                bad["detail"].append(f"{stage} {event} {reason}: {detail[:120]}")
    if os.path.exists(reads):
        with open(reads) as f:
            header = f.readline().rstrip("\n").split("\t")
            if header != READS_COLUMNS:
                bad["reads"].append("header " + " ".join(header))
            for line in f:
                c = line.rstrip("\n").split("\t")
                if len(c) != len(READS_COLUMNS) or c[4] not in ("0", "1", "2"):
                    bad["reads"].append(line.rstrip("\n")[:160])
    print(f"trace: {n} lines" + (", reads.tsv" if os.path.exists(reads) else ""))
    failed = False
    for check in ("header", "columns", "reasons", "detail", "reads"):
        items = bad.get(check, [])
        if not items:
            print(f"  {check:8} ok")
            continue
        failed = True
        distinct = collections.Counter(items)
        print(f"  {check:8} {len(items)} lines, {len(distinct)} distinct")
        for item, k in distinct.most_common(a.examples if check != "reasons" else None):
            print(f"    {k:6}  {item}")
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
