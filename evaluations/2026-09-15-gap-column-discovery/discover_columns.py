#!/usr/bin/env python3
"""Find phasing columns in a gap straight from the BAM, judged by mutual consistency.

Per-site classification cannot work here -- repeat-context sites measured 23%
informative, so any accept/reject on one site is near a coin flip. Real sites do
correlate with each other, though, because they reflect the same two haplotypes,
while systematic errors correlate with nothing. So candidate columns are taken
from the raw pileup with no catalog lookup and no repeat screening, and are then
admitted or dropped by how well they agree with the read partition the other
columns imply.

Outputs, per window: the columns admitted, the read partition, the partition's
self-consistency (computed without truth), and -- for evaluation only -- how the
admitted columns and the partition score against the read-level truth.
"""
import argparse
import collections
import csv
import statistics
from pathlib import Path

import pysam

p = argparse.ArgumentParser(description=__doc__)
p.add_argument('--bam', type=Path, required=True)
p.add_argument('--windows', type=Path, required=True, help='TSV with name,beg,end,kind')
p.add_argument('--truth-map', type=Path, required=True)
p.add_argument('--called-sites', type=Path, help='our own candidate TSV, for comparison')
p.add_argument('--reference-vcf', type=Path, help='het calls from a pileup caller, for comparison')
p.add_argument('--contig', default='CHM13#0#chr20')
p.add_argument('--reference-contig', default='chr20')
p.add_argument('--min-mapq', type=int, default=1)
p.add_argument('--min-depth', type=int, default=10)
p.add_argument('--min-allele-reads', type=int, default=3)
p.add_argument('--balance', type=float, default=0.25, help='keep columns with |AF-0.5| <= this')
p.add_argument('--agree-at', type=float, default=0.90, help='column-to-partition agreement to admit')
p.add_argument('--min-shared', type=int, default=10)
p.add_argument('--rounds', type=int, default=3)
p.add_argument('--site-output', type=Path, required=True)
p.add_argument('--window-output', type=Path, required=True)
a = p.parse_args()

truth = {}
with a.truth_map.open() as stream:
    for line in stream:
        name, hap = line.rstrip('\n').split('\t')
        truth[name] = hap
bam = pysam.AlignmentFile(str(a.bam), 'rb')


def candidate_columns(lo, hi):
    """Biallelic-looking pileup columns: no catalog, no repeat screening."""
    out = {}
    for column in bam.pileup(a.contig, lo, hi, truncate=True, min_base_quality=0,
                             max_depth=20000, stepper='nofilter'):
        calls = {}
        for read in column.pileups:
            if read.alignment.mapping_quality < a.min_mapq:
                continue
            if read.is_del or read.is_refskip or read.query_position is None:
                continue
            calls[read.alignment.query_name] = read.alignment.query_sequence[read.query_position].upper()
        if len(calls) < a.min_depth:
            continue
        counts = collections.Counter(calls.values()).most_common()
        if len(counts) < 2 or counts[1][1] < a.min_allele_reads:
            continue
        major, minor = counts[0][0], counts[1][0]
        af = counts[1][1] / (counts[0][1] + counts[1][1])
        if abs(af - 0.5) > a.balance:
            continue
        out[column.reference_pos + 1] = {
            'alleles': (major, minor), 'af': af,
            'calls': {k: v for k, v in calls.items() if v in (major, minor)}}
    return out


def partition_from(columns, admitted):
    """Majority vote: each read gets the side its admitted columns agree on."""
    votes = collections.defaultdict(lambda: [0, 0])
    for pos in admitted:
        major, minor = columns[pos]['alleles']
        for name, base in columns[pos]['calls'].items():
            votes[name][0 if base == major else 1] += 1
    side = {}
    for name, (a_votes, b_votes) in votes.items():
        if a_votes != b_votes:
            side[name] = 0 if a_votes > b_votes else 1
    return side


def agreement(columns, pos, side):
    major, minor = columns[pos]['alleles']
    shared = [(name, base) for name, base in columns[pos]['calls'].items() if name in side]
    if len(shared) < a.min_shared:
        return None, len(shared)
    straight = sum(1 for name, base in shared
                   if (base == major) == (side[name] == 0))
    return max(straight, len(shared) - straight) / len(shared), len(shared)


called = collections.defaultdict(set)
if a.called_sites:
    with a.called_sites.open() as stream:
        for r in csv.DictReader(stream, delimiter='\t'):
            called[r['CATEGORY']].add(int(r['POS']))
reference_het = set()
if a.reference_vcf:
    import subprocess
    done = subprocess.run(['tabix', '-l', str(a.reference_vcf)], capture_output=True, text=True)
    contig = a.reference_contig if a.reference_contig in done.stdout.split() else a.contig
    windows = [l.split('\t') for l in a.windows.open().read().splitlines()[1:]]
    for _, beg, end, _ in windows:
        got = subprocess.run(['tabix', str(a.reference_vcf), f'{contig}:{beg}-{end}'],
                             capture_output=True, text=True)
        for line in got.stdout.splitlines():
            f = line.split('\t')
            if f[9].split(':')[0] in ('0|1', '1|0', '0/1', '1/0'):
                reference_het.add(int(f[1]))

site_rows, window_rows = [], []
for row in csv.DictReader(a.windows.open(), delimiter='\t'):
    lo, hi = int(row['beg']), int(row['end'])
    columns = candidate_columns(lo, hi)
    if not columns:
        window_rows.append({'name': row['name'], 'kind': row['kind'], 'beg': lo, 'end': hi,
                            'candidate_columns': 0, 'admitted': 0, 'consistency': 0.0,
                            'reads_partitioned': 0, 'truth_concordance': 0.0,
                            'admitted_in_our_calls': 0, 'admitted_in_reference': 0,
                            'median_segregation': 0.0})
        continue
    # Seed on the best-balanced, deepest columns, then iterate admit/re-partition.
    seeds = sorted(columns, key=lambda p_: (abs(columns[p_]['af'] - 0.5),
                                            -len(columns[p_]['calls'])))[:max(3, len(columns) // 20)]
    admitted = list(seeds)
    side = partition_from(columns, admitted)
    for _ in range(a.rounds):
        scored = {}
        for pos in columns:
            score, shared = agreement(columns, pos, side)
            if score is not None and score >= a.agree_at:
                scored[pos] = score
        if not scored:
            break
        admitted = sorted(scored)
        side = partition_from(columns, admitted)
    if not admitted or not side:
        continue

    # Self-consistency: mean column-to-partition agreement over admitted columns.
    agreements = [agreement(columns, pos, side)[0] for pos in admitted]
    agreements = [x for x in agreements if x is not None]
    consistency = statistics.mean(agreements) if agreements else 0.0

    pairs = [(side[n], truth[n]) for n in side if n in truth]
    counts = collections.Counter(pairs)
    straight = counts[(0, 'MATERNAL')] + counts[(1, 'PATERNAL')]
    flipped = counts[(0, 'PATERNAL')] + counts[(1, 'MATERNAL')]
    concordance = max(straight, flipped) / len(pairs) if pairs else 0.0

    segs = []
    for pos in admitted:
        major, minor = columns[pos]['alleles']
        tp = [(base == major, truth[n]) for n, base in columns[pos]['calls'].items() if n in truth]
        if len(tp) < 10:
            continue
        c = collections.Counter(tp)
        s = max(c[(True, 'MATERNAL')] + c[(False, 'PATERNAL')],
                c[(True, 'PATERNAL')] + c[(False, 'MATERNAL')]) / len(tp)
        segs.append(s)
        site_rows.append({'window': row['name'], 'kind': row['kind'], 'pos': pos,
                          'af': round(columns[pos]['af'], 4), 'reads': len(columns[pos]['calls']),
                          'segregation': round(s, 4),
                          'in_our_calls': pos in set().union(*called.values()) if called else '',
                          'in_reference_vcf': pos in reference_het if reference_het else ''})
    window_rows.append({
        'name': row['name'], 'kind': row['kind'], 'beg': lo, 'end': hi,
        'candidate_columns': len(columns), 'admitted': len(admitted),
        'consistency': round(consistency, 4), 'reads_partitioned': len(side),
        'truth_concordance': round(concordance, 4),
        'admitted_in_our_calls': sum(1 for pos in admitted
                                     if called and pos in set().union(*called.values())),
        'admitted_in_reference': sum(1 for pos in admitted if pos in reference_het),
        'median_segregation': round(statistics.median(segs), 4) if segs else 0.0,
    })

for path, rows in ((a.site_output, site_rows), (a.window_output, window_rows)):
    with path.open('w', newline='') as stream:
        w = csv.DictWriter(stream, fieldnames=list(rows[0]), delimiter='\t')
        w.writeheader()
        w.writerows(rows)

print('%-22s %-8s %6s %6s %7s %7s %7s %8s %8s %8s' % (
    'window', 'kind', 'cand', 'admit', 'consist', 'reads', 'truth%', 'medSeg', 'ours', 'ref'))
for r in window_rows:
    print('%-22s %-8s %6d %6d %7.3f %7d %6.1f%% %8.3f %8d %8d' % (
        r['name'], r['kind'], r['candidate_columns'], r['admitted'], r['consistency'],
        r['reads_partitioned'], 100 * r['truth_concordance'], r['median_segregation'],
        r['admitted_in_our_calls'], r['admitted_in_reference']))
