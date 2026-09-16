#!/usr/bin/env python3
"""Locate, for every unresolved gap, the sub-interval of the BAM to attack.

A gap is bridgeable only where some read crosses it carrying informative sites on
both sides of the crossing. Scanning candidate cut points and counting such reads
gives a linkage profile; its minimum is the bottleneck, and the position of that
minimum -- not the gap interval as a whole -- is the region a targeted recovery
subprocess should be pointed at.

Each gap is then classified by which ingredient the bottleneck lacks, because the
remedy differs: a site desert needs discovery, a read desert cannot be fixed from
this BAM at all, a MAPQ-starved cut needs gated admission, and a cut with both
reads and sites present means the evidence is there and the solver or the stitch
is what failed.
"""
import argparse
import bisect
import collections
import csv
import re
import subprocess
from pathlib import Path

p = argparse.ArgumentParser(description=__doc__)
p.add_argument('--tiers', type=Path, required=True)
p.add_argument('--candidates', type=Path, required=True, help='our called sites (pipeline output)')
p.add_argument('--potential-vcf', type=Path, required=True, help='het calls a caller finds using all reads')
p.add_argument('--bam', type=Path, required=True)
p.add_argument('--contig', default='CHM13#0#chr20')
p.add_argument('--potential-contig', default='chr20')
p.add_argument('--flank', type=int, default=25000, help='read-reach added either side of the gap')
p.add_argument('--min-mapq', type=int, default=30)
p.add_argument('--max-cuts', type=int, default=160)
p.add_argument('--min-link-reads', type=int, default=2, help='reads needed across a cut to call it linkable')
p.add_argument('--output', type=Path, required=True)
a = p.parse_args()

CIGAR = re.compile(r'(\d+)([MIDNSHP=X])')
REF_OPS = frozenset('MDN=X')


def ref_span(pos, cigar):
    span = 0
    for length, op in CIGAR.findall(cigar):
        if op in REF_OPS:
            span += int(length)
    return pos, pos + max(span, 1) - 1


def reads_in(region):
    out = subprocess.run(['samtools', 'view', str(a.bam), region], capture_output=True, text=True).stdout
    rows = []
    for line in out.splitlines():
        f = line.split('\t', 6)
        if len(f) < 6:
            continue
        start, end = ref_span(int(f[3]), f[5])
        rows.append((start, end, int(f[4])))
    return rows


def potential_sites(region):
    out = subprocess.run(['tabix', str(a.potential_vcf), region], capture_output=True, text=True).stdout
    sites = []
    for line in out.splitlines():
        f = line.split('\t')
        if f[9].split(':')[0] in ('0|1', '1|0', '0/1', '1/0'):
            sites.append(int(f[1]))
    return sorted(sites)


INFORMATIVE = frozenset({'CLEAN_HET_SNP', 'CLEAN_HET_INDEL'})
our_sites = []
with a.candidates.open() as stream:
    for r in csv.DictReader(stream, delimiter='\t'):
        if r['CATEGORY'] in INFORMATIVE:
            our_sites.append(int(r['POS']))
our_sites.sort()

RANK = {'joined': 0, 'split': 1, 'partial': 2, 'vetoed': 3, 'rejected': 4, 'open': 5}
best = {}
with a.tiers.open() as stream:
    for r in csv.DictReader(stream, delimiter='\t'):
        key = (int(r['GAP_LEFT']), int(r['GAP_RIGHT']))
        if key not in best or RANK.get(r['STATUS'], 9) < RANK.get(best[key], 9):
            best[key] = r['STATUS']


def count_between(sorted_sites, lo, hi):
    """Informative sites in [lo, hi]."""
    return bisect.bisect_right(sorted_sites, hi) - bisect.bisect_left(sorted_sites, lo)


def profile(cuts, reads, sites, min_mapq):
    """For each cut: reads crossing it that carry >=1 site on each side."""
    out = []
    for cut in cuts:
        n = 0
        for start, end, mapq in reads:
            if mapq < min_mapq or start >= cut or end <= cut:
                continue
            if count_between(sites, start, cut - 1) and count_between(sites, cut + 1, end):
                n += 1
        out.append(n)
    return out


rows = []
for (left, right), status in sorted(best.items()):
    if status == 'joined':
        continue
    lo, hi = left - a.flank, right + a.flank
    reads = reads_in(f'{a.contig}:{lo}-{hi}')
    pot = potential_sites(f'{a.potential_contig}:{lo}-{hi}')
    ours_lo = bisect.bisect_left(our_sites, lo)
    ours_hi = bisect.bisect_right(our_sites, hi)
    ours = our_sites[ours_lo:ours_hi]

    # Cut points: the two junctions plus an even scan across the interior.
    interior = [left, right]
    if right - left > 2:
        step = max(1, (right - left) // a.max_cuts)
        interior.extend(range(left + step, right, step))
    cuts = sorted(set(interior))

    ours_prof = profile(cuts, reads, ours, a.min_mapq)
    pot_prof = profile(cuts, reads, pot, 0)
    i_min = min(range(len(cuts)), key=lambda i: (ours_prof[i], pot_prof[i]))
    bottleneck = cuts[i_min]

    # Which ingredient is missing at the bottleneck?
    win_lo, win_hi = bottleneck - a.flank, bottleneck + a.flank
    reads_here = sum(1 for s, e, q in reads if s < bottleneck < e)
    reads_here_pass = sum(1 for s, e, q in reads if s < bottleneck < e and q >= a.min_mapq)
    ours_here = count_between(ours, win_lo, win_hi)
    pot_here = count_between(pot, win_lo, win_hi)
    if reads_here < a.min_link_reads:
        cause, remedy = 'read_desert', 'unfixable_from_this_bam'
    elif pot_here == 0:
        cause, remedy = 'site_desert', 'msa_discovery_required'
    elif ours_here == 0 or (ours_prof[i_min] < a.min_link_reads <= pot_prof[i_min]):
        if reads_here_pass < a.min_link_reads <= reads_here:
            cause, remedy = 'mapq_starved', 'gated_admission_then_msa'
        else:
            cause, remedy = 'sites_not_called', 'msa_discovery_required'
    elif ours_prof[i_min] >= a.min_link_reads:
        cause, remedy = 'linkage_present', 'solver_or_stitch'
    else:
        cause, remedy = 'thin_linkage', 'gated_admission_then_msa'

    rows.append({
        'gap_left': left, 'gap_right': right, 'gap_bp': right - left, 'status': status,
        'bottleneck_pos': bottleneck,
        'bottleneck_is_junction': bottleneck in (left, right),
        'link_reads_ours': ours_prof[i_min], 'link_reads_potential': pot_prof[i_min],
        'reads_spanning': reads_here, 'reads_spanning_pass_mapq': reads_here_pass,
        'our_sites_near': ours_here, 'potential_sites_near': pot_here,
        'target_beg': max(1, win_lo), 'target_end': win_hi,
        'cause': cause, 'remedy': remedy,
    })

with a.output.open('w', newline='') as stream:
    w = csv.DictWriter(stream, fieldnames=list(rows[0]), delimiter='\t')
    w.writeheader()
    w.writerows(rows)
print('gaps profiled: %d' % len(rows))
for cause, n in collections.Counter(r['cause'] for r in rows).most_common():
    span = sum(r['gap_bp'] for r in rows if r['cause'] == cause)
    print('  %-16s %3d gaps  %8.2f Mb' % (cause, n, span / 1e6))
print('bottleneck at a junction rather than the interior: %d of %d' % (
    sum(1 for r in rows if r['bottleneck_is_junction']), len(rows)))
