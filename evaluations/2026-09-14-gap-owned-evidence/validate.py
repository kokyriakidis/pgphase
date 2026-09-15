#!/usr/bin/env python3
"""Check gap ownership and reconstruct the frozen pre-recovery read baseline."""
import argparse
import collections
import csv
import json
from pathlib import Path

import pysam


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('run', type=Path)
    parser.add_argument('baseline', type=Path)
    args = parser.parse_args()
    counts = collections.Counter()
    for path in sorted((args.run / 'audit').glob('*.events.tsv')):
        evidence = path.with_name(path.name.replace('.events.tsv', '.evidence.tsv'))
        with evidence.open() as stream:
            assert stream.readline().strip() == 'SCHEMA\t1'
            gap = stream.readline().strip().split('\t')
            left, right = int(gap[4]), int(gap[5])
        with path.open() as stream:
            for row in csv.DictReader(stream, delimiter='\t'):
                counts[row['ROLE']] += 1
                if row['ROLE'] == 'private':
                    beg, end = int(row['BEGIN_0']), int(row['END_0'])
                    # Both indel types sort at their preceding reference base.
                    alleles = row['ALLELES'].split(',')
                    snp = end - beg == 1 and len(alleles[1]) == 1
                    first = beg + 1 if snp else beg
                    last = end if end > beg else first
                    assert first > left and last < right, (path, row, left, right)
        counts['gaps'] += 1
    assert counts['gaps'] > 0
    members = {}
    for path in (args.run / 'audit').glob('*.members.tsv'):
        with path.open() as stream:
            for row in csv.DictReader(stream, delimiter='\t'):
                assert int(row['INPUT']) == 0, 'This validation handles one input BAM'
                name = row['READ']
                assert name not in members
                members[name] = int(row['HP']), int(row['PS'])
    args.baseline.mkdir(parents=True, exist_ok=True)
    output = args.baseline / 'phased.bam'
    assert not output.exists(), output
    seen = set()
    with pysam.AlignmentFile(str(args.run / 'phased.bam'), 'rb') as source:
        with pysam.AlignmentFile(str(output), 'wb', header=source.header) as dest:
            for read in source:
                read.set_tag('HP', None)
                read.set_tag('PS', None)
                if read.query_name in members:
                    hp, ps = members[read.query_name]
                    read.set_tag('HP', hp)
                    read.set_tag('PS', ps)
                    if not read.is_secondary and not read.is_supplementary:
                        seen.add(read.query_name)
                dest.write(read)
    assert seen == set(members), f'{len(set(members) - seen)} original molecules absent'
    pysam.index(str(output))
    counts['original_phased_molecules'] = len(members)
    (args.baseline / 'ownership.json').write_text(json.dumps(counts, indent=2) + '\n')
    print(json.dumps(counts, indent=2))


if __name__ == '__main__':
    main()
