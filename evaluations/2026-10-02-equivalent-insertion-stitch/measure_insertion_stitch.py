#!/usr/bin/env python3
"""Evaluation only: reproduce edit identity and exact output changes for chr20."""
import argparse
from collections import Counter
from pathlib import Path
import json
import pysam

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--before', type=Path, required=True)
parser.add_argument('--after', type=Path, required=True)
parser.add_argument('--matrix', type=Path, required=True)
parser.add_argument('--reference', default='test_data/chm13v2.0.chr20.renamed.fa')
parser.add_argument('--bam', default='test_data/HG002_chr20_hifi_mapped_to_CHM13_chr20_annotated.bam')
parser.add_argument('--output', type=Path, required=True)
args = parser.parse_args()


def calls(directory):
    result = {}
    for line in (directory / 'phased.vcf').open():
        if line.startswith('#'):
            continue
        fields = line.rstrip().split('\t')
        result[(int(fields[1]), fields[3], fields[4])] = dict(
            zip(fields[8].split(':'), fields[9].split(':')))
    return result


def tags(directory):
    with pysam.AlignmentFile(str(directory / 'phased.bam'), check_sq=False) as bam:
        return {read.query_name: (read.get_tag('HP') if read.has_tag('HP') else 0,
                                  read.get_tag('PS') if read.has_tag('PS') else 0)
                for read in bam if not read.is_secondary and not read.is_supplementary}


before, after = calls(args.before), calls(args.after)
left = [key for key in after if key[0] == 1196894]
right = [key for key in after if key[0] == 1196967]
assert len(left) == len(right) == 1
left, right = left[0], right[0]
assert left[1] == right[1] == 'A'
left_alt, right_alt = left[2][1:], right[2][1:]
with pysam.FastaFile(args.reference) as ref:
    background = ref.fetch('CHM13#0#chr20', left[0], right[0]).upper()
a, b = left_alt + background, background + right_alt
bare_mismatches = [i for i, (x, y) in enumerate(zip(a, b)) if x != y]
assert len(left_alt) == len(right_alt) == 155
assert background[1196950 - 1 - left[0]] == 'C'
common = list(background)
common[1196950 - 1 - left[0]] = 'T'
common = ''.join(common)
assert left_alt + common == common + right_alt

variants, observations, graph_mapq, verified = {}, {}, {}, {}
for line in args.matrix.open():
    fields = line.rstrip().split('\t')
    if fields[0] == 'VAR':
        variants[int(fields[1])] = fields
    elif fields[0] == '#META':
        verified[int(fields[1])] = 'msa_verified=1' in fields[2:]
    elif fields[0] == 'READ':
        graph_mapq[fields[1]] = int(fields[4]) if fields[5] == '0' else -1
    elif fields[0] == 'OBS':
        observations.setdefault(fields[1], {})[int(fields[2])] = (
            int(fields[4]), int(fields[5]))  # graph, BAM
left_index = next(i for i, row in variants.items()
                  if int(row[2]) == left[0] + 1 and row[7] == left_alt)
right_index = next(i for i, row in variants.items()
                   if int(row[2]) == right[0] and row[7].startswith('>'))
hom = [row for i, row in variants.items() if int(row[2]) == 1196950 and
       row[7] == 'T' and row[9:11] == ['1', '1'] and verified.get(i)]
assert len(hom) == 1
bam_mapq = {}
with pysam.AlignmentFile(args.bam) as bam:
    for read in bam.fetch('CHM13#0#chr20', left[0] - 1, right[0]):
        if read.flag & (4 | 256 | 2048 | 512 | 1024):
            continue
        bam_mapq[read.query_name] = read.mapping_quality
pairs = Counter()
for name, obs in observations.items():
    if graph_mapq.get(name, -1) < 30 or graph_mapq.get(name) == 255 or \
            bam_mapq.get(name, -1) < 30 or bam_mapq.get(name) == 255:
        continue
    bam_call = obs.get(left_index, (-1, -1))[1]
    graph_call = obs.get(right_index, (-1, -1))[0]
    if bam_call in (0, 1) and graph_call in (0, 1):
        pairs[(bam_call, graph_call)] += 1

bt, at = tags(args.before), tags(args.after)
shared = before.keys() & after.keys()
changed = [key for key in shared if before[key] != after[key]]
result = {
    'reference_edit': {'left_vcf_pos': left[0], 'right_vcf_pos': right[0],
                       'inserted_length': len(left_alt),
                       'bare_reference_mismatch_offsets': bare_mismatches,
                       'common_snp': {'pos': 1196950, 'ref': 'C', 'alt': 'T',
                                      'ref_cov': int(hom[0][14]), 'alt_cov': int(hom[0][15])},
                       'equal_on_common_background': True},
    'paired_calls': {f'{a}/{b}': n for (a, b), n in sorted(pairs.items())},
    'exact_comparison': {
        'reads_before': len(bt), 'reads_after': len(at),
        'read_hap_changes': sum(bt.get(name, (0, 0))[0] != at.get(name, (0, 0))[0]
                                for name in bt.keys() | at.keys()),
        'read_ps_changes': sum(bt.get(name, (0, 0))[1] != at.get(name, (0, 0))[1]
                               for name in bt.keys() | at.keys()),
        'read_only_assignment_changes': sum(
            bt[name] != at.get(name) for name in bt if bt[name][1] >= 1000000000),
        'read_ps_transitions': {
            f'{old}/{new}': n for (old, new), n in Counter(
                (bt[name][1], at[name][1]) for name in bt.keys() & at.keys()
                if bt[name][1] != at[name][1]).items()},
        'lost_vcf_keys': len(before.keys() - after.keys()),
        'gained_vcf_keys': len(after.keys() - before.keys()),
        'changed_vcf_rows': len(changed),
        'changed_gt_depth_or_other_fields': sum(
            {k: v for k, v in before[key].items() if k != 'PS'} !=
            {k: v for k, v in after[key].items() if k != 'PS'} for key in changed)},
    'connection': {'left_gt_ps': [after[left]['GT'], after[left]['PS']],
                   'right_gt_ps': [after[right]['GT'], after[right]['PS']]}}
args.output.write_text(json.dumps(result, indent=2) + '\n')
print(json.dumps(result, indent=2))
