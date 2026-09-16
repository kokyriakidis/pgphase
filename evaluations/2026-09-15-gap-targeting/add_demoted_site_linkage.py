#!/usr/bin/env python3
"""Re-test each gap's bottleneck using graph sites the pipeline currently excludes.

The first pass measured linkage with the sites the pipeline calls phase-informative
(CLEAN_HET_SNP / CLEAN_HET_INDEL) and, as a second opinion, a pileup caller's het
calls. Both omit repeat-context indels, so a bottleneck can look like a
heterozygosity desert while the snarl catalog holds usable het sites there that
apply_graph_noise_filter demoted (REP_HET_INDEL) or the AF gate dropped (low_af).
This pass adds those classes and records whether they would bridge.

Linkage here is a candidate, not a conclusion: these sites were demoted because
per-read genotypes at homopolymer/STR indels are unreliable, so any that bridge
must be MSA-verified inside the target window before being trusted as anchors.
"""
import argparse
import bisect
import collections
import csv
import re
import subprocess
from pathlib import Path

p = argparse.ArgumentParser(description=__doc__)
p.add_argument('--targets', type=Path, required=True, help='gap_targets.tsv from locate_gap_evidence.py')
p.add_argument('--graph-variants', type=Path, required=True, help='retained sites from collect-graph-variation')
p.add_argument('--graph-filtered', type=Path, required=True, help='filtered-sites diagnostic TSV')
p.add_argument('--bam', type=Path, required=True)
p.add_argument('--contig', default='CHM13#0#chr20')
p.add_argument('--flank', type=int, default=25000)
p.add_argument('--min-mapq', type=int, default=30)
p.add_argument('--min-link-reads', type=int, default=2)
p.add_argument('--output', type=Path, required=True)
a = p.parse_args()

CIGAR = re.compile(r'(\d+)([MIDNSHP=X])')
REF_OPS = frozenset('MDN=X')

by_category = collections.defaultdict(list)
with a.graph_variants.open() as stream:
    for r in csv.DictReader(stream, delimiter='\t'):
        by_category[r['CATEGORY']].append(int(r['POS']))
low_af = []
with a.graph_filtered.open() as stream:
    for r in csv.DictReader(stream, delimiter='\t'):
        if r['REASON'] == 'low_af':
            try:
                low_af.append(int(r['POS']))
            except ValueError:
                pass
clean = sorted(by_category['CLEAN_HET_SNP'] + by_category['CLEAN_HET_INDEL'])
repeat_indels = sorted(by_category['REP_HET_INDEL'])
low_af.sort()


def reads_near(lo, hi):
    out = subprocess.run(['samtools', 'view', str(a.bam), f'{a.contig}:{lo}-{hi}'],
                         capture_output=True, text=True).stdout
    rows = []
    for line in out.splitlines():
        f = line.split('\t', 6)
        if len(f) < 6:
            continue
        pos = int(f[3])
        span = sum(int(n) for n, op in CIGAR.findall(f[5]) if op in REF_OPS)
        rows.append((pos, pos + max(span, 1) - 1, int(f[4])))
    return rows


def window(sorted_sites, lo, hi):
    return sorted_sites[bisect.bisect_left(sorted_sites, lo):bisect.bisect_right(sorted_sites, hi)]


def linking_reads(reads, sites, cut):
    n = 0
    for start, end, mapq in reads:
        if mapq < a.min_mapq or start >= cut or end <= cut:
            continue
        left = bisect.bisect_right(sites, cut - 1) - bisect.bisect_left(sites, start)
        right = bisect.bisect_right(sites, end) - bisect.bisect_left(sites, cut + 1)
        if left and right:
            n += 1
    return n


rows = []
for r in csv.DictReader(a.targets.open(), delimiter='\t'):
    cut = int(r['bottleneck_pos'])
    lo, hi = cut - a.flank, cut + a.flank
    reads = reads_near(lo, hi)
    c, p_rep, p_low = window(clean, lo, hi), window(repeat_indels, lo, hi), window(low_af, lo, hi)
    link_clean = linking_reads(reads, c, cut)
    link_rep = linking_reads(reads, sorted(c + p_rep), cut)
    link_all = linking_reads(reads, sorted(c + p_rep + p_low), cut)
    if link_clean >= a.min_link_reads:
        revised = r['cause']
    elif link_rep >= a.min_link_reads:
        revised = 'repeat_indels_would_bridge'
    elif link_all >= a.min_link_reads:
        revised = 'low_af_sites_would_bridge'
    else:
        revised = 'no_linkage_from_any_site_class'
    out = dict(r)
    out.update({'graph_clean_near': len(c), 'graph_repeat_indels_near': len(p_rep),
                'graph_low_af_near': len(p_low), 'link_clean': link_clean,
                'link_with_repeat_indels': link_rep, 'link_with_low_af': link_all,
                'revised_cause': revised})
    rows.append(out)

with a.output.open('w', newline='') as stream:
    w = csv.DictWriter(stream, fieldnames=list(rows[0]), delimiter='\t')
    w.writeheader()
    w.writerows(rows)

print('gaps re-tested: %d' % len(rows))
for cause, n in collections.Counter(r['revised_cause'] for r in rows).most_common():
    span = sum(int(r['gap_bp']) for r in rows if r['revised_cause'] == cause)
    print('  %-32s %3d gaps  %7.2f Mb' % (cause, n, span / 1e6))
