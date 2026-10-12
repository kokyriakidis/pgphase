"""Shared loaders for the diagnostic harness (evaluation only).

A diagnostic directory (`scripts/dev/pg diag`) holds a full-contig run of
./pgphase made with every dump on:
  phased.vcf, phased.bam      the run's output
  em/em.TID.BEG.tsv           final EM evidence per read (PGPHASE_EM_DUMP)
  aln/aln.TID.BEG.tsv         alignment-solve candidates before pruning
  aln/noisy.TID.BEG.tsv       the solve's noisy regions
  aln/inj.TID.BEG.tsv         each solve candidate's injection decision
  aln/graph.TID.BEG.tsv       the graph chunk's rows before the merge
  aln/emsites.TID.BEG.tsv     the final EM's voting candidates: allele counts,
                              undecided reads, contrast, error, phase, block
  aln/joins.TID.BEG.tsv       every EM block-cut decision between adjacent sites
  aln/msa.TID.BEG.tsv         every noisy-region MSA variant: cluster sizes,
                              consensus coverage, read calls before / after
                              the refresh of assigned reads
  aln/alnreads.TID.BEG.tsv    each read's calls in the alignment solve
Positions are VCF-convention (indels at their anchor base) throughout, so
every source lines up with the VCFs; a catalog or EM SNP may sit one base
off, which matching tolerates.
"""
import bisect
import collections
import glob
import gzip
import os
import re
import sys

import pysam

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import score  # noqa: E402
import why  # noqa: E402

ENV = os.environ
CATS = ["LOW_COV", "LOW_AF", "STRAND_BIAS", "CLEAN_HET_SNP", "CLEAN_HET_INDEL", "CLEAN_HOM",
        "NOISY_CAND_HET", "NOISY_CAND_HOM", "NOISY_RESOLVED", "REP_HET_INDEL", "NON_VAR"]
SHORT = {"LOW_COV": "LC", "LOW_AF": "LAF", "STRAND_BIAS": "SB", "CLEAN_HET_SNP": "HS",
         "CLEAN_HET_INDEL": "HI", "CLEAN_HOM": "HOM", "NOISY_CAND_HET": "NH", "NOISY_CAND_HOM": "NHOM",
         "NOISY_RESOLVED": "NR", "REP_HET_INDEL": "RHI", "NON_VAR": "NV"}
CHUNK_SLACK = 2_000_000


def cat(code):
    return CATS[int(code)] if 0 <= int(code) < len(CATS) else "?"


def competitors():
    out = {}
    for pair in ENV.get("PGDEV_COMPETITORS", "").split():
        name, d = pair.split("=", 1)
        out[name] = d
    return out


def parse_region(text):
    """'chrN:A-B', 'A-B' (contig from PGDEV_CONTIG) or 'A' (+-5 kb); 1-based inclusive."""
    contig = ENV["PGDEV_CONTIG"]
    if ":" in text:
        contig, text = text.split(":", 1)
    text = text.replace(",", "")
    if "-" in text:
        a, b = text.split("-", 1)
        return contig, int(a), int(b)
    p = int(text)
    return contig, p - 5000, p + 5000


def chunk_files(directory, prefix, beg=None, end=None):
    """Files PREFIX.TID.BEG.tsv in DIRECTORY whose chunk may overlap [beg, end]."""
    out = []
    for f in glob.glob(os.path.join(directory, f"{prefix}.*.tsv")):
        m = re.search(r"\.(\d+)\.tsv$", f)
        cb = int(m.group(1)) if m else 0
        if beg is None or (cb <= end and cb >= beg - CHUNK_SLACK):
            out.append(f)
    return sorted(out)


def load_rows(directory, prefix, beg=None, end=None):
    """Rows (split lines) of every chunk file, first column a position in [beg, end]."""
    rows = []
    for f in chunk_files(directory, prefix, beg, end):
        with open(f) as fh:
            for line in fh:
                p = line.rstrip("\n").split("\t")
                pos = int(p[0])
                if beg is None or beg <= pos <= end:
                    rows.append(p)
    return rows


def load_aln(diag, beg=None, end=None):
    """pos -> list of alignment candidates: dict(type, ref_len, alt, init, final, ref, alt_cov, total, unresolved)."""
    out = collections.defaultdict(list)
    seen = set()
    for p in load_rows(os.path.join(diag, "aln"), "aln", beg, end):
        key = (p[0], p[2], p[3], p[4])
        if key in seen:  # overlapping chunks
            continue
        seen.add(key)
        out[int(p[0])].append(dict(type=p[2], ref_len=int(p[3]), alt=p[4], init=cat(p[5]), final=cat(p[6]),
                                   ref=int(p[7]), alt_cov=int(p[8]), total=int(p[9]), unresolved=p[10] == "1"))
    return out


def load_inj(diag, beg=None, end=None):
    """pos -> list of injection decisions: dict(type, ref_len, alt, cat, msa, unresolved, depth, median, conf, decision)."""
    out = collections.defaultdict(list)
    rank = {"injected": 0, "complement": 1, "depth_gate": 2, "few_confident_reads": 3,
            "outside_chunk": 9, "category": 5, "not_germline": 6}
    best = {}
    for p in load_rows(os.path.join(diag, "aln"), "inj", beg, end):
        key = (p[0], p[2], p[3], p[4])
        rec = dict(type=p[2], ref_len=int(p[3]), alt=p[4], cat=cat(p[5]), msa=p[6] == "1", unresolved=p[7] == "1",
                   depth=int(p[8]), median=float(p[9]), conf=p[10], decision=p[11])
        if key not in best or rank.get(rec["decision"], 8) < rank.get(best[key]["decision"], 8):
            best[key] = rec
    for key, rec in best.items():
        out[int(key[0])].append(rec)
    return out


def load_graph(diag, beg=None, end=None):
    """pos -> list of graph rows: dict(type, ref_len, alt, cat, voting, id, covs)."""
    out = collections.defaultdict(list)
    seen = set()
    for p in load_rows(os.path.join(diag, "aln"), "graph", beg, end):
        key = (p[0], p[2], p[3], p[4], p[7])
        if key in seen:
            continue
        seen.add(key)
        out[int(p[0])].append(dict(type=p[2], ref_len=int(p[3]), alt=p[4], cat=cat(p[5]), voting=p[6] == "1",
                                   id=p[7], covs=p[8]))
    return out


def load_emsites(diag, beg=None, end=None):
    """pos -> list of dict(type, ref_len, alt, origin, counts, undecided, in_em, x, y, err, hap1, block)."""
    out = collections.defaultdict(list)
    seen = set()
    for p in load_rows(os.path.join(diag, "aln"), "emsites", beg, end):
        key = (p[0], p[1], p[2], p[3])
        if key in seen:
            continue
        seen.add(key)
        em = p[7] == "em"
        out[int(p[0])].append(dict(type=p[1], ref_len=int(p[2]), alt=p[3], origin=p[4], counts=p[5],
                                   undecided=int(p[6]), in_em=em, x=int(p[8]) if em else None,
                                   y=int(p[9]) if em else None, err=float(p[10]) if em else None,
                                   hap1=int(p[11]) if em else None, block=int(p[12]) if em else None))
    return out


def load_noisy(diag, beg, end):
    regs = set()
    for f in chunk_files(os.path.join(diag, "aln"), "noisy", beg, end):
        with open(f) as fh:
            for line in fh:
                a, b = map(int, line.split())
                if b >= beg and a <= end:
                    regs.add((a, b))
    return sorted(regs)


def load_em(diag, beg=None, end=None):
    """(reads, sites): reads as why.load_dump; sites pos -> dict(kind, err, block, weak, n, mat, pat, sides)."""
    em_dir = os.path.join(diag, "em")
    if beg is not None:
        # why.load_dump reads whole directories; restrict to the chunks near the region
        tmp = {}
        for f in chunk_files(em_dir, "em", beg, end):
            tmp.update(_load_em_file(f, tmp))
        reads = tmp
    else:
        reads = why.load_dump(em_dir)
    scope = score.read_scope()
    sites = {}
    for q, d in reads.items():
        for pos, kind, side, err, blk in d["obs"]:
            if beg is not None and not (beg - 1 <= pos <= end):
                continue
            s = sites.setdefault(pos, dict(kind=kind, err=err, block=blk, n=0, agree=0, truth=0))
            s["n"] += 1
            if q in scope:
                s["truth"] += 1
                s["agree"] += (side == 0) == scope[q][2]
    for s in sites.values():
        s["purity"] = max(s["agree"], s["truth"] - s["agree"]) / s["truth"] if s["truth"] else None
    return reads, sites


def _load_em_file(path, prev):
    out = {}
    with open(path) as fh:
        for line in fh:
            p = line.rstrip("\n").split("\t")
            obs = []
            first = 6 if len(p) > 5 and ":" not in p[5] else 5
            for tok in p[first:]:
                f = tok.split(":")
                obs.append((int(f[0]), f[1], int(f[2]), float(f[3]), int(f[4])))
            old = prev.get(p[0]) or out.get(p[0])
            if old is None or len(obs) > len(old["obs"]):
                out[p[0]] = {"mapq": int(p[1]), "hap": int(p[2]), "ps": int(p[3]), "margins": {}, "obs": obs}
    return out


def vcf_records(path, contig, beg, end):
    """pos -> list of (ref, alts, gt, ps, cat) for records in [beg, end]; tabix when indexed."""
    if not os.path.exists(path) and os.path.exists(path + ".gz"):
        path += ".gz"
    out = collections.defaultdict(list)
    if path.endswith(".gz") and (os.path.exists(path + ".tbi") or os.path.exists(path + ".csi")):
        with pysam.TabixFile(path) as tf:
            name = tabix_contig(tf, contig)
            if name is not None:
                for line in tf.fetch(name, max(0, beg - 1), min(end, 2 ** 31 - 1)):
                    _vcf_line(line, out)
    else:
        opener = gzip.open if path.endswith(".gz") else open
        with opener(path, "rt") as fh:
            for line in fh:
                if line[0] == "#":
                    continue
                pos = int(line.split("\t", 2)[1])
                if beg <= pos <= end:
                    _vcf_line(line, out)
    return out


def _vcf_line(line, out):
    f = line.rstrip("\n").split("\t")
    fmt = dict(zip(f[8].split(":"), f[9].split(":"))) if len(f) > 9 else {}
    raw = fmt.get("GT", ".")
    sep = "|" if "|" in raw else "/"
    gt = tuple(None if a == "." else int(a) for a in raw.split(sep))
    ps = fmt.get("PS") if sep == "|" else None
    c = next((kv[4:] for kv in f[7].split(";") if kv.startswith("CAT=")), "")
    out[int(f[1])].append((f[3], f[4].split(","), gt, None if ps in (None, ".", "") else ps, c))


def tabix_contig(tf, contig):
    """The indexed name for contig: as given, or with or without a PanSN prefix (CHM13#0#)."""
    names = set(tf.contigs)
    for c in (contig, contig.split("#")[-1]):
        if c in names:
            return c
    for c in names:
        if c.split("#")[-1] == contig.split("#")[-1]:
            return c
    return None


def catalog(contig, beg, end):
    """pos -> list of (ref_len, n_alts, level)."""
    out = collections.defaultdict(list)
    with pysam.TabixFile(ENV["PGDEV_SITES"]) as tf:
        name = tabix_contig(tf, contig)
        if name is None:
            return out
        for line in tf.fetch(name, max(0, beg - 1), min(end, 2 ** 31 - 1)):
            f = line.split("\t", 8)
            lv = re.search(r"LV=(\d+)", f[7])
            out[int(f[1])].append((len(f[3]), len(f[4].split(",")), int(lv.group(1)) if lv else 0))
    return out


def truth_records(contig, beg, end):
    """pos -> (ref, alt, gt, hap1_allele or None) for heterozygous truth records."""
    out = {}
    for pos, recs in vcf_records(ENV["PGDEV_TRUTH_VCF"], contig, beg, end).items():
        for ref, alts, gt, ps, c in recs:
            if len(gt) == 2 and None not in gt and gt[0] != gt[1]:
                out[pos] = (ref, ",".join(alts), gt)
    return out


def near(index, pos, tol=1):
    """Entries of a pos-keyed dict within tol of pos, closest first."""
    hits = []
    for d in range(tol + 1):
        for p in ((pos,) if d == 0 else (pos - d, pos + d)):
            if p in index:
                hits.append((p, index[p]))
    return hits


def load_joins(diag, beg=None, end=None):
    """Block-cut decisions: list of dict(left, right, delta, strong, spanning, decision, strong_reads)
    sorted by left; strong_reads = reads with strong evidence on both sides."""
    out, seen = [], set()
    for p in load_rows(os.path.join(diag, "aln"), "joins", beg, end):
        key = (p[0], p[1])
        if key in seen:
            continue
        seen.add(key)
        out.append(dict(left=int(p[0]), right=int(p[1]), delta=float(p[2]), strong=float(p[3]),
                        spanning=int(p[4]), decision=p[5], strong_reads=int(p[6]) if len(p) > 6 else None))
    return sorted(out, key=lambda r: r["left"])


def load_msa(diag, beg=None, end=None):
    """pos -> list of MSA variants: dict(reg_beg, reg_len, clusters (n0, n1), aligned, longest_gap, type,
    ref_len, alt, cat, in_cons (first, second), before (ref, alt, missing), after (ref, alt, missing))."""
    out, seen = collections.defaultdict(list), set()
    for f in chunk_files(os.path.join(diag, "aln"), "msa", beg, end):
        with open(f) as fh:
            for line in fh:
                p = line.rstrip("\n").split("\t")
                pos = int(p[8])
                if beg is not None and not (beg <= pos <= end):
                    continue
                key = tuple(p[8:12])
                if key in seen:
                    continue
                seen.add(key)
                out[pos].append(dict(reg_beg=int(p[0]), reg_len=int(p[1]), clusters=(int(p[2]), int(p[3])),
                                     aligned=(int(p[4]), int(p[6])), longest_gap=(int(p[5]), int(p[7])),
                                     type=p[9], ref_len=int(p[10]), alt=p[11], cat=cat(p[12]),
                                     in_cons=(p[13] == "1", p[14] == "1"),
                                     before=tuple(int(x) for x in p[15].split(",")),
                                     after=tuple(int(x) for x in p[16].split(","))))
    return out


def load_alnreads(diag, beg, end, names=None):
    """name -> dict(mapq, skipped, hap, ps, calls {pos: (type, allele)}) from the alignment solve,
    calls restricted to [beg, end]; reads from every chunk file near the region."""
    out = {}
    for f in chunk_files(os.path.join(diag, "aln"), "alnreads", beg, end):
        with open(f) as fh:
            for line in fh:
                p = line.rstrip("\n").split("\t")
                if names is not None and p[0] not in names:
                    continue
                calls = {}
                for tok in p[5:]:
                    pos, typ, al = tok.split(":")
                    pos = int(pos)
                    if beg <= pos <= end:
                        calls[pos] = (typ, int(al))
                old = out.get(p[0])
                if old is None or len(calls) > len(old["calls"]):
                    out[p[0]] = dict(mapq=int(p[1]), skipped=p[2] == "1", hap=int(p[3]), ps=int(p[4]), calls=calls)
    return out
