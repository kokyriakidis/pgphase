#!/usr/bin/env python3
"""Truth-free contracts a pgphase phased VCF must keep, checked against the run's own MSA dump and stage trace.

  check_contracts.py --vcf PHASED_VCF --ref FASTA[.gz] [--msa-dump TSV] [--trace TSV]
                     [--known TSV] [--write-known TSV]

Each contract compares two things the pipeline itself produced, so a violation is a bug or a deliberate policy, never
a disagreement with a truth set:

  vcf_order    records in POS order; no record written twice
  ref_match    every REF is the reference's bases at POS
  hap_overlap  no two records whose REF spans overlap put ALT alleles on one haplotype (1/1 beside another ALT; two
               phased records of one phase set on the same haplotype)
  msa_spell    (--msa-dump, PGPHASE_MSA_DUMP) inside each noisy region the MSA resolved, the records spell its
               consensus sequences: where the two differ, exactly that pair; where they agree, no het. With --trace
               (PGPHASE_TRACE), each violation is given a cause from the region's injection decisions: a policy that
               keeps the MSA het out on purpose, no MSA het candidate, or unexplained.

vcf_order and ref_match must hold. The others are compared with --known, the violations accepted so far (one key per
line): a violation not in it fails the check; one in it that no longer occurs is reported as fixed. --write-known
writes this run's violations as the new list. Python 3 standard library only.
"""
import argparse
import bisect
import collections
import difflib
import gzip
import itertools
import sys

Record = collections.namedtuple("Record", "chrom pos ref alleles gt phased ps info")
Region = collections.namedtuple("Region", "chrom beg end branch ref consensus")


def open_text(path):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path)


def read_vcf(path):
    out = []
    with open_text(path) as f:
        for line in f:
            if line.startswith("#"):
                continue
            c = line.rstrip("\n").split("\t")
            sample = dict(zip(c[8].split(":"), c[9].split(":"))) if len(c) > 9 else {}
            gt_text = sample.get("GT", "./.")
            phased = "|" in gt_text
            gt = tuple(None if x == "." else int(x) for x in gt_text.replace("|", "/").split("/"))
            ps = sample.get("PS")
            out.append(Record(c[0], int(c[1]), c[3].upper(), [c[3].upper()] + [a.upper() for a in c[4].split(",")],
                              gt, phased, int(ps) if phased and ps not in (None, ".") else None, c[7]))
    return out


def read_fasta(path, names):
    """The sequences of the contigs in names."""
    out, name, parts = {}, None, []
    with open_text(path) as f:
        for line in f:
            if line.startswith(">"):
                if name in names:
                    out[name] = "".join(parts).upper()
                name, parts = line[1:].split()[0], []
            elif name in names:
                parts.append(line.strip())
    if name in names:
        out[name] = "".join(parts).upper()
    return out


def vcf_contigs(path):
    """The contig names of the VCF header, in order."""
    out = []
    with open_text(path) as f:
        for line in f:
            if not line.startswith("##"):
                break
            if line.startswith("##contig=<ID="):
                out.append(line[len("##contig=<ID="):].split(",")[0].rstrip(">\n"))
    return out


def read_msa_dump(path, contigs):
    """Each region the MSA resolved (the last dump line of a region wins): reference and consensus sequences. The dump
    names a contig by its index in the alignments' header: the VCF's one contig, or its header's contig there."""
    regions = {}
    for line in open_text(path):
        c = line.rstrip("\n").split("\t")
        n = int(c[5])
        chrom = contigs[0] if len(contigs) == 1 else contigs[int(c[0])]
        regions[(chrom, int(c[1]), int(c[2]))] = Region(chrom, int(c[1]), int(c[2]), c[3], c[4].upper(),
                                                        [c[6 + 2 * i].upper() for i in range(n)])
    return [regions[k] for k in sorted(regions)]


def read_injection_decisions(path):
    """Sorted (pos, type, decision, category): each alignment candidate's last injection decision. The trace names
    no contig: a run with --msa-dump and --trace is checked on one contig."""
    last = {}
    with open_text(path) as f:
        next(f)
        for line in f:
            c = line.rstrip("\n").split("\t")
            if len(c) >= 7 and c[0] == "inject" and c[1] == "decision":
                cat = c[6].split()[0].split("=")[1] if c[6].startswith("cat=") else "-1"
                last[(int(c[2]), c[4])] = (int(c[2]), c[4][0], c[5], int(cat))
    return sorted(last.values())


def end(r):
    return r.pos + len(r.ref) - 1


def carried(r):
    """The allele index on each haplotype (0 for REF or no call)."""
    return [a if a is not None else 0 for a in r.gt[:2]] if len(r.gt) >= 2 else [0, 0]


def is_het(r):
    h = carried(r)
    return h[0] != h[1]


def is_hom_alt(r):
    h = carried(r)
    return h[0] == h[1] and h[0] > 0


def alt_text(r):
    return f"{r.chrom}:{r.pos} {r.ref}>{','.join(r.alleles[1:])} {'|'.join(map(str, carried(r)))}"


def spell(beg, seq, records, flips=()):
    """The window's two haplotype sequences with records applied (record k's alleles swapped when k is in flips);
    None when two records overlap on a haplotype."""
    haps = []
    for h in range(2):
        out, at = [], beg
        for k, r in enumerate(records):
            a = carried(r)[1 - h if k in flips else h]
            if a == 0 or r.alleles[a] == "*":
                continue
            if r.pos < at:
                return None
            out.append(seq[at - beg:r.pos - beg])
            out.append(r.alleles[a])
            at = r.pos + len(r.ref)
        out.append(seq[at - beg:])
        haps.append("".join(out))
    return haps


def spells(beg, seq, records, want, max_free=8):
    """(a pair the records spell, whether it is want): phased hets of one phase set swap together, every other het
    alone. (None, False) when the records cannot be spelled at all."""
    records = sorted(records, key=lambda r: r.pos)
    units = collections.OrderedDict()
    for k, r in enumerate(records):
        if is_het(r):
            units.setdefault(("ps", r.ps) if r.ps is not None else ("rec", k), []).append(k)
    free = list(units.values())[1:max_free + 1]
    first = None
    for choice in itertools.product([False, True], repeat=len(free)):
        haps = spell(beg, seq, records, {k for unit, f in zip(free, choice) if f for k in unit})
        if haps is None:
            continue
        if sorted(haps) == sorted(want):
            return haps, True
        first = first or haps
    return first, False


def check_order(records):
    out, seen = [], set()
    for prev, r in zip([None] + records[:-1], records):
        if prev is not None and prev.chrom == r.chrom and r.pos < prev.pos:
            out.append(("vcf_order", f"{r.chrom}:{r.pos} after {prev.pos}"))
        key = (r.chrom, r.pos, r.ref, tuple(r.alleles[1:]), r.gt)
        if key in seen:
            out.append(("vcf_order", f"written twice: {alt_text(r)}"))
        seen.add(key)
    return out


def check_ref(records, fasta):
    out = []
    for r in records:
        want = fasta.get(r.chrom, "")[r.pos - 1:r.pos - 1 + len(r.ref)]
        if want != r.ref:
            out.append(("ref_match", f"{alt_text(r)}: the reference has {want or '-'}"))
    return out


def nested_kind(outer, inner):
    """Whether the later record, inside the earlier one's REF span, repeats a change the earlier one's ALT makes there
    ("redundant") or contradicts it; "partial" when it reaches past it."""
    if outer.pos > inner.pos:
        outer, inner = inner, outer
    if end(inner) > end(outer):
        return "partial"
    alt = outer.alleles[max(carried(outer))]
    inner_alt = inner.alleles[max(carried(inner))]
    off = inner.pos - outer.pos
    for tag, i1, i2, j1, _ in difflib.SequenceMatcher(None, outer.ref, alt, autojunk=False).get_opcodes():
        if i1 <= off < i2 or i1 == i2 == off:
            if tag == "equal":
                return "contradictory" if inner.ref != inner_alt else "redundant"
            span = alt[j1 + off - i1:j1 + off - i1 + len(inner.ref)] if tag == "replace" else ""
            return "redundant" if span == inner_alt else "contradictory"
    return "contradictory"


def check_overlap(records):
    out, active = [], []
    for r in sorted((r for r in records if any(a > 0 for a in carried(r))), key=lambda r: (r.chrom, r.pos)):
        active = [a for a in active if a.chrom == r.chrom and end(a) >= r.pos]
        for a in active:
            clash = None
            if is_hom_alt(a) or is_hom_alt(r):
                clash = "1/1 beside another ALT"
            elif a.ps is not None and a.ps == r.ps:
                if any(carried(a)[h] > 0 and carried(r)[h] > 0 for h in range(2)):
                    clash = "one haplotype of a phase set"
            if clash:
                out.append(("hap_overlap", f"{alt_text(a)} and {alt_text(r)}: {clash}, {nested_kind(a, r)}"))
        active.append(r)
    return out


def check_msa(records, regions, decisions):
    out, checked = [], 0
    by_chrom = collections.defaultdict(list)
    for r in records:
        by_chrom[r.chrom].append(r)
    starts = {c: [r.pos for r in rs] for c, rs in by_chrom.items()}
    dpos = [d[0] for d in decisions]
    for reg in regions:
        if not reg.consensus:
            continue
        rs, ps = by_chrom.get(reg.chrom, []), starts.get(reg.chrom, [])
        # From the base before the region (an insertion anchors there) to its end; a record reaching past either
        # edge leaves the region unchecked.
        near = [r for r in rs[bisect.bisect_left(ps, reg.beg - 2000):bisect.bisect_right(ps, reg.end)]
                if end(r) >= reg.beg - 1]
        if any(r.pos < reg.beg or end(r) > reg.end for r in near):
            continue
        checked += 1
        het_region = len(reg.consensus) == 2 and reg.consensus[0] != reg.consensus[1]
        cands = [d for d in decisions[bisect.bisect_left(dpos, reg.beg):bisect.bisect_right(dpos, reg.end)]
                 if d[3] in (6, 7)]
        het_ds = [d[2] for d in cands if d[3] == 6]
        if not decisions:
            cause = "no trace"
        elif not het_ds:
            cause = "no MSA het candidate"
        elif any(x.startswith("island:") for x in het_ds):
            cause = "policy: island-only het"
        elif any(x in ("depth_gate", "few_confident_reads") for x in het_ds):
            cause = "policy: gated het"
        else:
            cause = "unexplained"
        where = f"{reg.chrom}:{reg.beg}-{reg.end}"
        if het_region:
            got, ok = spells(reg.beg, reg.ref, near, reg.consensus)
            if ok:
                continue
            if got is None:
                kind = "records overlap"
            elif got[0] == got[1]:
                kind = "het written hom" if any(is_het(r) for r in near) else "het unwritten"
            else:
                kind = "another pair"
            out.append(("msa_spell", f"{where} {kind} [{cause}]"))
        elif any(is_het(r) for r in near):
            out.append(("msa_spell", f"{where} het where the MSA reads one haplotype [{cause}]"))
    return out, checked


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--vcf", required=True)
    ap.add_argument("--ref", required=True)
    ap.add_argument("--msa-dump")
    ap.add_argument("--trace")
    ap.add_argument("--known")
    ap.add_argument("--write-known")
    a = ap.parse_args()
    records = read_vcf(a.vcf)
    fasta = read_fasta(a.ref, {r.chrom for r in records})
    found = check_order(records) + check_ref(records, fasta) + check_overlap(records)
    checked = 0
    if a.msa_dump:
        msa, checked = check_msa(records, read_msa_dump(a.msa_dump, vcf_contigs(a.vcf)),
                                 read_injection_decisions(a.trace) if a.trace else [])
        found += msa
    counts = collections.Counter(c for c, _ in found)
    print(f"contracts: {len(records)} records" + (f", {checked} MSA regions" if a.msa_dump else ""))
    for c in ("vcf_order", "ref_match", "hap_overlap", "msa_spell"):
        state = "not checked" if c == "msa_spell" and not a.msa_dump else \
            "ok" if not counts[c] else f"{counts[c]} violations"
        print(f"  {c:12} {state}")
    keys = sorted({f"{c}\t{d}" for c, d in found if c not in ("vcf_order", "ref_match")})
    if a.write_known:
        with open(a.write_known, "w") as f:
            f.write("".join(k + "\n" for k in keys))
        print(f"  accepted violations written: {a.write_known} ({len(keys)})")
    failed = [f"{c}\t{d}" for c, d in found if c in ("vcf_order", "ref_match")]
    if a.known and not a.write_known:
        known = {line.rstrip("\n") for line in open(a.known) if line.strip()}
        failed += [k for k in keys if k not in known]
        fixed = sorted(known - set(keys))
        for k in fixed:
            print(f"  fixed (no longer occurs; refresh {a.known}): {k}")
    for k in failed:
        print(f"  VIOLATION {k}", file=sys.stderr)
    sys.exit(1 if failed else 0)


if __name__ == "__main__":
    main()
