"""Our unreliable injected sites behind remaining reads vs HiPhase and truth."""
import collections
import gzip
import sys

sys.path.insert(0, ".")
import audit_hybrid as A

run, dump_path, qlist = sys.argv[1:4]
C = "/home/kokyriakidis/Downloads/pgphase-eval-data/results/chr12-18-20-comparison/chr20"
truth = A.load_truth("../derived/chr20_truth_hap.tsv")
tags = A.load_tags(f"{run}/phased.bam")
targets = {q for q in open(qlist).read().split() if q not in tags}
per_site = collections.Counter()
calls = collections.defaultdict(collections.Counter)
for line in open(dump_path):
    f = line.rstrip("\n").split("\t")
    for o in (x.split(":") for x in f[6:]):
        if o[1] == "-1" or o[2] != "1":
            continue
        p = int(o[0])
        if f[0] in targets and float(o[3]) >= 0.3:
            per_site[p] += 1
        if f[0] in truth:
            calls[p][(o[5], truth[f[0]])] += 1


def vcf_near(path, p, w):
    out = []
    op = gzip.open if path.endswith(".gz") else open
    for line in op(path, "rt"):
        if line[0] == "#":
            continue
        f = line.split("\t", 10)
        q = int(f[1])
        if q > p + w:
            break
        if q >= p - w:
            out.append((q, f[3][:14], f[4][:24], f[9].split(":")[0]))
    return out


rows = collections.defaultdict(list)
with open(f"{run}/candidates.tsv") as fh:
    next(fh)
    for line in fh:
        f = line.rstrip("\n").split("\t")
        rows[int(f[1])].append((f[2], f[3][:14], f[4][:16], f[6], f[7], f[14]))
for p, n in per_site.most_common(6):
    c = calls[p]
    cis = c[("1", "MATERNAL")] + c[("2", "PATERNAL")]
    tr = c[("2", "MATERNAL")] + c[("1", "PATERNAL")]
    print(f"== site {p}: {n} target reads; truth split cis {cis} trans {tr}; our rows {rows.get(p) or rows.get(p+1) or rows.get(p-1)}")
    print("   HiPhase:", vcf_near(f"{C}/hiphase/phased.vcf.gz", p, 30))
    print("   truth:  ", vcf_near(f"{C}/truth.vcf.gz", p, 30))
