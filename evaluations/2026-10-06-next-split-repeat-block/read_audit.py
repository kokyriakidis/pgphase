#!/usr/bin/env python3
"""Read-level scoring with the committed C++ scorer's deterministic tie rule."""
from collections import Counter, defaultdict
from pathlib import Path
import pysam


def assignments(path, truth):
    tags, votes = {}, defaultdict(Counter)
    with pysam.AlignmentFile(str(path), check_sq=False) as bam:
        for read in bam:
            if read.is_secondary or read.is_supplementary:
                continue
            hp = read.get_tag('HP') if read.has_tag('HP') else 0
            ps = read.get_tag('PS') if read.has_tag('PS') else 0
            tags[read.query_name] = (hp, ps)
            if hp in (1, 2) and ps > 0 and read.query_name in truth:
                votes[ps][(hp == 1) != truth[read.query_name]] += 1
    orientation = {ps: counts[True] >= counts[False] for ps, counts in votes.items()}
    status = {q: 'unphased' if hp not in (1, 2) or ps <= 0 else
              'correct' if ((hp == 1) != truth[q]) == orientation[ps] else 'discordant'
              for q, (hp, ps) in tags.items() if q in truth}
    return tags, status


def variants(path):
    result = {}
    if not path.exists():
        path = path.with_name('native.vcf')
    with open(path) as vcf:
        for line in vcf:
            if line.startswith('#'):
                continue
            fields = line.rstrip().split('\t')
            key = tuple(fields[:5])
            assert key not in result, key
            result[key] = (tuple(fields[5:9]), dict(zip(fields[8].split(':'), fields[9].split(':'))))
    return result


def block_extents(rows):
    positions = defaultdict(list)
    for key, (_, sample) in rows.items():
        gt, ps = sample.get('GT', ''), int(sample.get('PS', '0'))
        if ps > 0 and gt in ('0|1', '1|0'):
            positions[ps].append(int(key[1]))
    return [(min(pos), max(pos)) for pos in positions.values() if len(pos) > 1]
