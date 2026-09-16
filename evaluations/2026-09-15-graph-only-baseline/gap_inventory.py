#!/usr/bin/env python3
"""Derive the phase-block gap inventory from a graph-only run.

`collect-graph-variation` has no gap finder of its own -- gap recovery lives in
the hybrid path -- so the gaps a follow-up subprocess would have to close are
defined here the same way `find_phase_gaps` defines them: sort the phase blocks
by position and take the space between consecutive blocks. Each gap is reported
with the reads sitting inside it and how many of those the run left untagged,
since that is the work a gap subprocess would be asked to do.
"""
import argparse
import collections
import csv
import re
import subprocess
from pathlib import Path

p = argparse.ArgumentParser(description=__doc__)
p.add_argument('--phase-sites', type=Path, required=True)
p.add_argument('--phase-reads', type=Path, required=True,
               help='per-read assignment TSV; the phased BAM carries no contig header, '
                    'so region queries against it fail and cannot be used here')
p.add_argument('--bam', type=Path, required=True, help='surjected input BAM (read inventory)')
p.add_argument('--contig', default='CHM13#0#chr20')
p.add_argument('--min-block-sites', type=int, default=2)
p.add_argument('--output', type=Path, required=True)
a = p.parse_args()

CIGAR = re.compile(r'(\d+)([MIDNSHP=X])')
REF_OPS = frozenset('MDN=X')

# Phase blocks: span of the sites sharing a phase set.
spans = {}
with a.phase_sites.open() as stream:
    reader = csv.DictReader(stream, delimiter='\t')
    ps_field = 'PHASE_SET' if 'PHASE_SET' in reader.fieldnames else reader.fieldnames[-1]
    pos_field = 'POS' if 'POS' in reader.fieldnames else reader.fieldnames[1]
    for r in reader:
        ps = r[ps_field]
        if ps in ('', '-1', '0', '.'):
            continue
        pos = int(r[pos_field])
        lo, hi, n = spans.get(ps, (pos, pos, 0))
        spans[ps] = (min(lo, pos), max(hi, pos), n + 1)

blocks = sorted((lo, hi, ps, n) for ps, (lo, hi, n) in spans.items()
                if n >= a.min_block_sites)
print('phase blocks with >= %d sites: %d (of %d phase sets seen)'
      % (a.min_block_sites, len(blocks), len(spans)))
total_block_bp = sum(hi - lo for lo, hi, _, _ in blocks)


def read_spans(region):
    done = subprocess.run(['samtools', 'view', str(a.bam), region],
                          capture_output=True, text=True)
    if done.returncode != 0:
        raise SystemExit('samtools view failed on %s: %s' % (region, done.stderr.strip()))
    out = done.stdout
    rows = []
    for line in out.splitlines():
        f = line.split('\t', 6)
        if len(f) < 6:
            continue
        pos = int(f[3])
        span = sum(int(n) for n, op in CIGAR.findall(f[5]) if op in REF_OPS)
        rows.append((f[0], pos, pos + max(span, 1) - 1, int(f[4])))
    return rows


tagged = set()
with a.phase_reads.open() as stream:
    for r in csv.DictReader(stream, delimiter='\t'):
        if r['PHASE_SET'] not in ('', '0', '-1', '.') and r['HAP'] in ('1', '2'):
            tagged.add(r['READ'])
if not tagged:
    raise SystemExit('no phased reads found in %s' % a.phase_reads)
print('reads carrying a haplotype assignment: %d' % len(tagged))


rows = []
# A read can overlap two gaps when the block between them is shorter than a read,
# so cross-gap totals are deduplicated by read name.
untagged_any = set()
untagged_any_pass = set()
for (lo1, hi1, ps1, n1), (lo2, hi2, ps2, n2) in zip(blocks, blocks[1:]):
    if lo2 <= hi1:
        continue  # overlapping blocks: not a gap
    gap_lo, gap_hi = hi1, lo2
    reads = read_spans(f'{a.contig}:{gap_lo}-{gap_hi}')
    inside = [r for r in reads if r[1] >= gap_lo and r[2] <= gap_hi]
    untagged_any.update(r[0] for r in reads if r[0] not in tagged)
    untagged_any_pass.update(r[0] for r in reads if r[0] not in tagged and r[3] >= 30)
    rows.append({
        'gap_left': gap_lo, 'gap_right': gap_hi, 'gap_bp': gap_hi - gap_lo,
        'left_ps': ps1, 'right_ps': ps2,
        'left_block_bp': hi1 - lo1, 'right_block_bp': hi2 - lo2,
        'left_block_sites': n1, 'right_block_sites': n2,
        'reads_overlapping': len(reads),
        'reads_contained': len(inside),
        'reads_untagged': sum(1 for r in reads if r[0] not in tagged),
        'reads_untagged_pass_mapq30': sum(1 for r in reads if r[0] not in tagged and r[3] >= 30),
    })

with a.output.open('w', newline='') as stream:
    w = csv.DictWriter(stream, fieldnames=list(rows[0]), delimiter='\t')
    w.writeheader()
    w.writerows(rows)

span = sum(r['gap_bp'] for r in rows)
print('gaps between consecutive blocks: %d   total gap span %.2f Mb   total block span %.2f Mb'
      % (len(rows), span / 1e6, total_block_bp / 1e6))
buckets = [(0, 1000), (1000, 10000), (10000, 50000), (50000, 200000), (200000, 10 ** 9)]
print('gap size distribution:')
for lo, hi in buckets:
    sel = [r for r in rows if lo <= r['gap_bp'] < hi]
    if sel:
        print('  %8s - %-8s  %4d gaps  %7.2f Mb  %8d untagged reads'
              % (lo, hi, len(sel), sum(r['gap_bp'] for r in sel) / 1e6,
                 sum(r['reads_untagged'] for r in sel)))
print('distinct untagged reads in gap windows: %d (%d at MAPQ >= 30)'
      % (len(untagged_any), len(untagged_any_pass)))
print('  sum over gaps before deduplication: %d (reads spanning two gaps counted twice)'
      % sum(r['reads_untagged'] for r in rows))
