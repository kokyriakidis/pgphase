#!/usr/bin/env python3
"""Audit graph-catalog events omitted or ineligible in frozen gap projections."""
import argparse
import bisect
import collections
import csv
import json
from pathlib import Path

import pysam


def key(pos, ref, alt):
    prefix = 0
    while prefix < min(len(ref), len(alt)) and ref[prefix] == alt[prefix]:
        prefix += 1
    suffix = 0
    if len(ref) == len(alt):
        while suffix + prefix < len(ref) and ref[-suffix-1] == alt[-suffix-1]:
            suffix += 1
        if len(ref) - prefix - suffix == 1:
            return pos + prefix - 1, pos + prefix, alt[prefix], pos + prefix
    beg = pos + prefix - 1
    return beg, pos + len(ref) - 1, alt[prefix:], beg


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('audit', type=Path)
    p.add_argument('catalog', type=Path)
    p.add_argument('output', type=Path)
    a = p.parse_args()
    gaps = []
    for path in sorted(a.audit.glob('*.events.tsv')):
        evidence = path.with_name(path.name.replace('.events.tsv', '.evidence.tsv'))
        sites, skipped = {}, set()
        with evidence.open() as f:
            f.readline()
            gap = f.readline().rstrip().split('\t')
            left, right = int(gap[4]), int(gap[5])
            for line in f:
                fields = line.rstrip('\n').split('\t')
                if fields[0] == 'SITE':
                    beg = int(fields[2]) - 1
                    sites[beg, beg + int(fields[4]), fields[5]] = int(fields[7])
                elif fields[0] == 'READ' and fields[7] == '1':
                    skipped.add((fields[2], fields[3]))
        with path.open() as f:
            events = {int(r['EVENT']): r for r in csv.DictReader(f, delimiter='\t')}
        counts = collections.defaultdict(collections.Counter)
        observations = path.with_name(path.name.replace('.events.tsv', '.observations.tsv'))
        with observations.open() as f:
            for r in csv.DictReader(f, delimiter='\t'):
                if (r['INPUT'], r['READ']) in skipped:
                    continue
                if r['GRAPH_STATUS'] == 'observed':
                    counts[int(r['EVENT'])][int(r['GRAPH_ALLELE'])] += 1
        by_key = {}
        for ei, event in events.items():
            alleles = event['ALLELES'].split(',')
            if len(alleles) < 2:
                continue
            k = int(event['BEGIN_0']), int(event['END_0']), alleles[1]
            by_key[k] = event, counts[ei], sites.get(k)
        gaps.append((left, right, by_key))
    gaps.sort()
    starts = [g[0] for g in gaps]
    rows, seen = [], set()
    catalog = pysam.VariantFile(str(a.catalog))
    for record in catalog:
        for alt in record.alts or ():
            beg, end, alt_seq, sort_pos = key(record.pos, record.ref, alt)
            gi = bisect.bisect_left(starts, sort_pos) - 1
            if gi < 0:
                continue
            left, right, events = gaps[gi]
            if max(sort_pos, end) >= right:
                continue
            k = beg, end, alt_seq
            if (gi, k) in seen:
                continue
            seen.add((gi, k))
            event, counts, flags = events.get(k, ({'ROLE': 'absent'}, {}, None))
            total = sum(counts.values())
            balanced = total >= 5 and sum(n >= 2 and n / total >= .2 for n in counts.values()) >= 2
            role = event['ROLE']
            reason = role if role in ('absent', 'unsupported', 'boundary') else (
                'clean_eligible' if flags is not None and flags & 0x108c else
                'msa_eligible' if flags is not None and flags & 0x300 else
                'category_excluded')
            rows.append([left, right, record.pos, beg, end, alt_seq, reason,
                         '' if flags is None else hex(flags), total, counts.get(0, 0),
                         total - counts.get(0, 0), int(balanced)])
    a.output.mkdir(parents=True, exist_ok=True)
    with (a.output / 'graph_sites.tsv').open('w') as f:
        w = csv.writer(f, delimiter='\t', lineterminator='\n')
        w.writerow(['gap_left', 'gap_right', 'vcf_pos', 'begin_0', 'end_0', 'alt',
                    'reason', 'flags', 'graph_depth', 'graph_ref', 'graph_alt', 'balanced_graph'])
        w.writerows(rows)
    summary = dict(gaps=len(gaps), catalog_events=len(rows),
                   reasons=dict(collections.Counter(r[6] for r in rows)),
                   balanced_reasons=dict(collections.Counter(r[6] for r in rows if r[-1])))
    (a.output / 'summary.json').write_text(json.dumps(summary, indent=2) + '\n')
    print(json.dumps(summary, indent=2))


if __name__ == '__main__':
    main()
