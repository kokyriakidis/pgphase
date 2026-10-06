#!/usr/bin/env python3
"""Audit a real block union without requiring old PS labels to stay unchanged."""
import argparse
from collections import Counter, defaultdict
import json
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
    orientation = {ps: counts.most_common(1)[0][0] for ps, counts in votes.items()}
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


def audit(before, after, truth):
    old, old_status = assignments(before / 'phased.bam', truth)
    new, new_status = assignments(after / 'phased.bam', truth)
    lost = [q for q, state in old_status.items()
            if state == 'correct' and new_status.get(q) != 'correct']
    assert not lost, lost
    dropped_tags = [q for q, (hp, ps) in old.items()
                    if hp in (1, 2) and ps > 0 and
                    (q not in new or new[q][0] not in (1, 2) or new[q][1] <= 0)]
    assert not dropped_tags, dropped_tags
    changed = [q for q in old if old[q] != new.get(q)]
    gauge_maps = defaultdict(set)
    for q in changed:
        hp, ps = old[q]
        nhp, nps = new[q]
        if hp in (1, 2) and ps > 0:
            gauge_maps[(ps, nps)].add(hp != nhp)
    assert all(len(flips) == 1 for flips in gauge_maps.values()), gauge_maps
    old_vcf, new_vcf = variants(before / 'phased.vcf'), variants(after / 'phased.vcf')
    assert old_vcf.keys() <= new_vcf.keys(), old_vcf.keys() - new_vcf.keys()
    changed_vcf = 0
    for key, (info, sample) in old_vcf.items():
        ninfo, nsample = new_vcf[key]
        assert info == ninfo, key
        assert {k: v for k, v in sample.items() if k not in ('GT', 'PS')} == {
            k: v for k, v in nsample.items() if k not in ('GT', 'PS')}, key
        gt, ngt = sample.get('GT'), nsample.get('GT')
        if gt != ngt:
            assert '|' in gt and ngt == '|'.join(reversed(gt.split('|'))), (key, gt, ngt)
            pair = (int(sample['PS']), int(nsample['PS']))
            assert gauge_maps.get(pair) == {True}, pair
        changed_vcf += sample != nsample
    return {'old_primary_output_reads': len(old), 'new_primary_output_reads': len(new),
            'old_correct': Counter(old_status.values())['correct'],
            'new_correct': Counter(new_status.values())['correct'],
            'old_discordant': Counter(old_status.values())['discordant'],
            'new_discordant': Counter(new_status.values())['discordant'],
            'old_variant_records': len(old_vcf), 'new_variant_records': len(new_vcf),
            'changed_variant_gauges': changed_vcf, 'changed_read_tags': len(changed),
            'lost_correct_assignments': 0, 'lost_phased_assignments': 0,
            'variant_alleles_counts_and_filters_preserved': True,
            'new_assignments': [q for q in new if new[q][0] in (1, 2) and new[q][1] > 0 and
                                (q not in old or old[q][0] not in (1, 2))]}


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--before', type=Path, required=True)
    parser.add_argument('--after', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    truth = {f[0]: f[1] == 'PATERNAL'
             for line in open('test_data/derived/chr20_truth_hap.tsv')
             if len(f := line.rstrip().split('\t')) == 2}
    result = audit(args.before, args.after, truth)
    args.output.write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps(result, indent=2))
