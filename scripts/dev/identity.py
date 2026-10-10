#!/usr/bin/env python3
"""Are two pgphase runs output-identical?

Usage: identity.py DIR_A DIR_B
Compares candidates.tsv byte for byte, phased.vcf without its ## header
lines (they carry the command line), and the HP/PS tags of every primary
record in phased.bam. Exits 0 when all three match, 1 otherwise.
"""
import itertools
import os
import sys

import pysam


def vcf_body(path):
    with open(path) as fh:
        return [line for line in fh if not line.startswith("##")]


def first_difference(a, b):
    for i, (x, y) in enumerate(itertools.zip_longest(a, b)):
        if x != y:
            return i, x, y
    return None


def bam_tags(path):
    out = {}
    with pysam.AlignmentFile(path, check_sq=False) as fh:
        for r in fh.fetch(until_eof=True):
            if r.is_secondary or r.is_supplementary:
                continue
            out[r.query_name] = (r.get_tag("HP") if r.has_tag("HP") else None,
                                 r.get_tag("PS") if r.has_tag("PS") else None)
    return out


def main():
    a, b = sys.argv[1:3]
    same = True
    with open(os.path.join(a, "candidates.tsv"), "rb") as fa, open(os.path.join(b, "candidates.tsv"), "rb") as fb:
        ok = fa.read() == fb.read()
    print(f"candidates.tsv  {'identical' if ok else 'DIFFERS'}")
    same &= ok
    va, vb = vcf_body(os.path.join(a, "phased.vcf")), vcf_body(os.path.join(b, "phased.vcf"))
    d = first_difference(va, vb)
    print(f"phased.vcf      {'identical' if d is None else 'DIFFERS'} ({len(va):,} vs {len(vb):,} lines)")
    if d is not None:
        print(f"  first difference at body line {d[0] + 1}:\n  - {(d[1] or '').rstrip()[:160]}\n  + {(d[2] or '').rstrip()[:160]}")
    same &= d is None
    ta, tb = bam_tags(os.path.join(a, "phased.bam")), bam_tags(os.path.join(b, "phased.bam"))
    changed = [q for q in ta.keys() | tb.keys() if ta.get(q) != tb.get(q)]
    print(f"phased.bam tags {'identical' if not changed else 'DIFFERS'} "
          f"({len(ta):,} vs {len(tb):,} reads, {len(changed):,} changed)")
    same &= not changed
    return 0 if same else 1


if __name__ == "__main__":
    sys.exit(main())
