#!/usr/bin/env python3
"""Test whether the sub-MAPQ-floor reads in starved gaps carry correct phase.

For each MAPQ-starved gap, take the reads pgphase leaves unphased and ask what a
competitor made of them, scored against the diplinator read truth. This
separates gaps where the discarded reads are genuinely informative from gaps
where they are noise, which the raw starvation count cannot distinguish.
"""
import argparse
import csv
import gzip
import statistics
import subprocess
from pathlib import Path

p = argparse.ArgumentParser(description=__doc__)
p.add_argument('--screen', type=Path, required=True, help='unresolved_gaps.tsv from screen.py')
p.add_argument('--bam', type=Path, required=True, help='surjected input BAM')
p.add_argument('--phased-bam', type=Path, required=True, help='our phased output BAM')
p.add_argument('--eval-root', type=Path, required=True, help='frozen per-tool read evaluations')
p.add_argument('--competitor-vcf', type=Path, required=True)
p.add_argument('--contig', default='CHM13#0#chr20')
p.add_argument('--competitor-contig', default='chr20')
p.add_argument('--min-mapq', type=int, default=30)
p.add_argument('--min-input-depth', type=float, default=20.0)
p.add_argument('--min-depth-lost-pct', type=float, default=50.0)
p.add_argument('--trust-accuracy', type=float, default=95.0)
p.add_argument('--trust-min-reads', type=int, default=20)
p.add_argument('--output', type=Path, required=True)
a = p.parse_args()


def tool_status(tool):
    with gzip.open(a.eval_root / f'{tool}_reads/per_read.tsv.gz', 'rt') as stream:
        return {r['read_name']: r['status'].lower() for r in csv.DictReader(stream, delimiter='\t')}


def input_reads(region):
    out = subprocess.run(['samtools', 'view', str(a.bam), region], capture_output=True, text=True).stdout
    return {line.split('\t', 5)[0]: int(line.split('\t', 5)[4]) for line in out.splitlines()}


def our_phased(region):
    out = subprocess.run(['samtools', 'view', str(a.phased_bam), region], capture_output=True, text=True).stdout
    names = set()
    for line in out.splitlines():
        fields = line.split('\t')
        if any(x.startswith('PS:i:') for x in fields[11:]):
            names.add(fields[0])
    return names


def vaf_deviation(region):
    """Median |VAF-0.5| at competitor het calls: a truth-free noise proxy."""
    out = subprocess.run(['tabix', str(a.competitor_vcf), region], capture_output=True, text=True).stdout
    devs = []
    for line in out.splitlines():
        fields = line.split('\t')
        if fields[9].split(':')[0] not in ('0|1', '1|0'):
            continue
        record = dict(zip(fields[8].split(':'), fields[9].split(':')))
        try:
            devs.append(abs(float(record.get('VAF', 'nan')) - 0.5))
        except ValueError:
            pass
    return (statistics.median(devs) if devs else '', len(devs))


tools = ('longphase', 'hiphase')
status = {t: tool_status(t) for t in tools}
starved = [r for r in csv.DictReader(a.screen.open(), delimiter='\t')
           if float(r['depth_all']) >= a.min_input_depth
           and float(r['depth_lost_pct']) >= a.min_depth_lost_pct]

rows = []
for gap in starved:
    left, right = int(gap['gap_left']), int(gap['gap_right'])
    region = f'{a.contig}:{left}-{right}'
    mapq = input_reads(region)
    phased = our_phased(region)
    unphased = [n for n in mapq if n not in phased]
    dev, sites = vaf_deviation(f'{a.competitor_contig}:{left}-{right}')
    row = {'gap_left': left, 'gap_right': right, 'gap_bp': right - left,
           'status': gap['status'], 'depth_all': gap['depth_all'],
           'depth_pass_mapq': gap['depth_pass_mapq'], 'depth_lost_pct': gap['depth_lost_pct'],
           'reads_in_window': len(mapq), 'unphased_by_us': len(unphased),
           'unphased_below_floor': sum(1 for n in unphased if mapq[n] < a.min_mapq),
           'competitor_het_sites': sites,
           'median_vaf_deviation': round(dev, 4) if dev != '' else ''}
    for tool in tools:
        scored = [n for n in unphased if status[tool].get(n) in ('concordant', 'discordant')]
        concordant = sum(1 for n in scored if status[tool][n] == 'concordant')
        row[f'{tool}_phased'] = len(scored)
        row[f'{tool}_concordant'] = concordant
        row[f'{tool}_accuracy'] = round(100.0 * concordant / len(scored), 1) if scored else ''
    row['verdict'] = ('trustworthy'
                      if row['longphase_phased'] >= a.trust_min_reads
                      and row['longphase_accuracy'] != ''
                      and row['longphase_accuracy'] >= a.trust_accuracy
                      else 'noise')
    rows.append(row)

rows.sort(key=lambda r: -(r['longphase_accuracy'] if r['longphase_accuracy'] != '' else -1))
with a.output.open('w', newline='') as stream:
    writer = csv.DictWriter(stream, fieldnames=list(rows[0]), delimiter='\t')
    writer.writeheader()
    writer.writerows(rows)

trust = [r for r in rows if r['verdict'] == 'trustworthy']
print('starved gaps examined      : %d' % len(rows))
print('trustworthy                : %d gaps, %d kb, %d reads we leave unphased, %d competitor-concordant' % (
    len(trust), sum(r['gap_bp'] for r in trust) // 1000,
    sum(r['unphased_by_us'] for r in trust), sum(r['longphase_concordant'] for r in trust)))
print('noise                      : %d gaps, %d reads we leave unphased' % (
    len(rows) - len(trust), sum(r['unphased_by_us'] for r in rows if r['verdict'] == 'noise')))
print('wrote %s' % a.output)
