"""Reads correct in run A that are discordant in run B, by position."""
import collections
import sys

import pysam

sys.path.insert(0, ".")
import audit_hybrid as A

a_dir, b_dir = sys.argv[1:3]
truth = A.load_truth("../derived/chr20_truth_hap.tsv")
ta, tb = A.load_tags(f"{a_dir}/phased.bam"), A.load_tags(f"{b_dir}/phased.bam")
_, oa = A.score_reads(ta, truth, set(truth))
_, ob = A.score_reads(tb, truth, set(truth))


def ok(t, o, q):
    return ((t[q][0] == 1) == (truth[q] == "MATERNAL")) == o[t[q][1]]


changed = collections.Counter()
ps_b = collections.Counter()
bam = pysam.AlignmentFile("../HG002_chr20_hifi_mapped_to_CHM13_chr20_annotated.bam")
pos = {}
for r in bam.fetch(until_eof=True):
    if not r.is_secondary and not r.is_supplementary and r.query_name in truth:
        pos[r.query_name] = r.reference_start
for q in truth:
    if q not in tb or ok(tb, ob, q):
        continue
    was = "unphased" if q not in ta else ("correct" if ok(ta, oa, q) else "discordant")
    if was == "discordant":
        continue
    changed[(was, pos.get(q, -1) // 1_000_000)] += 1
    ps_b[tb[q][1]] += 1
print("newly discordant reads by (previous state, Mb):")
for k, v in sorted(changed.items(), key=lambda kv: -kv[1])[:12]:
    print(v, k)
print("top phase sets in B:", ps_b.most_common(5))
