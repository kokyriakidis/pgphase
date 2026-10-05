#!/usr/bin/env python3
"""Audit a gap closure against exact old calls and independent parental truth."""
import argparse
from collections import Counter, defaultdict
import json
from pathlib import Path
import pysam

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--before', type=Path, required=True)
parser.add_argument('--after', type=Path, required=True)
parser.add_argument('--output', type=Path, required=True)
parser.add_argument('--truth', default='test_data/derived/chr20_truth_hap.tsv')
parser.add_argument('--bam', default='test_data/HG002_chr20_hifi_mapped_to_CHM13_chr20_annotated.bam')
args = parser.parse_args()
truth = {f[0]: f[1] == 'PATERNAL' for line in open(args.truth)
         if len(f := line.rstrip().split('\t')) == 2}


def read_tags(folder):
    tags, votes = {}, defaultdict(Counter)
    with pysam.AlignmentFile(str(folder / 'phased.bam'), check_sq=False) as bam:
        for read in bam:
            if read.is_secondary or read.is_supplementary or not read.has_tag('HP') or not read.has_tag('PS'):
                continue
            hp, ps = read.get_tag('HP'), read.get_tag('PS')
            if hp not in (1, 2) or ps <= 0:
                continue
            tags[read.query_name] = (hp, ps)
            if read.query_name in truth:
                votes[ps][(hp == 1) != truth[read.query_name]] += 1
    orientation = {ps: counts.most_common()[0][0] for ps, counts in votes.items()}
    status = {name: ((hp == 1) != truth[name]) == orientation[ps]
              for name, (hp, ps) in tags.items() if name in truth}
    return tags, status


def variants(folder):
    rows, extents = {}, defaultdict(list)
    with pysam.VariantFile(str(folder / 'phased.vcf')) as vcf:
        for record in vcf:
            sample = next(iter(record.samples.values()))
            gt, ps = sample.get('GT'), sample.get('PS')
            key = (record.contig, record.pos, record.ref, record.alts)
            fields = str(record).rstrip().split('\t')
            nonphase = fields[:9] + [v for k, v in zip(fields[8].split(':'), fields[9].split(':')) if k not in ('GT', 'PS')]
            rows[key] = (gt, ps if sample.phased and ps else 0, nonphase)
            if sample.phased and ps and gt and len(set(gt)) > 1:
                extents[ps].append(record.pos)
    return rows, {ps: (min(pos), max(pos)) for ps, pos in extents.items()}


before_tags, before_status = read_tags(args.before)
after_tags, after_status = read_tags(args.after)
gap_fill_ps_offset = 1_000_000_000
changed_output_only_tags = [name for name, tag in before_tags.items()
                            if tag[1] >= gap_fill_ps_offset and after_tags.get(name) != tag]
transitions = Counter((before_status.get(n), after_status.get(n))
                      for n in before_status.keys() | after_status.keys())
before_rescue_groups, after_rescue_groups = defaultdict(set), defaultdict(set)
for name, (hp, ps) in before_tags.items():
    if ps >= gap_fill_ps_offset: before_rescue_groups[ps].add((name, hp))
for name, (hp, ps) in after_tags.items():
    if ps >= gap_fill_ps_offset: after_rescue_groups[ps].add((name, hp))
rescue_groups_preserved = ({frozenset(v) for v in before_rescue_groups.values()} ==
                           {frozenset(v) for v in after_rescue_groups.values()})
before_rows, before_extents = variants(args.before)
after_rows, after_extents = variants(args.after)
transfers = defaultdict(set)
missing_phased, changed_genotypes, changed_nonphase = [], [], []
for key, (gt, ps, nonphase) in before_rows.items():
    if key not in after_rows:
        continue
    new_gt, new_ps, new_nonphase = after_rows[key]
    if sorted(gt) != sorted(new_gt):
        changed_genotypes.append(str(key))
    if nonphase != new_nonphase:
        changed_nonphase.append(str(key))
    if ps and len(set(gt)) > 1:
        if ps == 62637077:
            assert (new_gt != gt) == (key[1] >= 62722021), key
        if not new_ps:
            missing_phased.append(str(key))
        else:
            transfers[ps].add((new_ps, new_gt != gt))


def spans(extents, left, right):
    return any(a <= left and z >= right for a, z in extents.values())


old_gaps = []
previous = None
for left, right in sorted(before_extents.values()):
    if previous is not None and left > previous:
        old_gaps.append((previous, left))
    previous = max(previous or 0, right)
new_connections = [(a, z) for a, z in old_gaps if spans(after_extents, a, z)]
tracked, reopened = [], []
for line in open('evaluations/2026-09-16-test-panel/panel.tsv'):
    if line.startswith('#') or not line.strip() or not line.split('\t')[0].isdigit():
        continue
    left, right = map(int, line.split('\t')[:2])
    tracked.append(spans(after_extents, left, right))
    if spans(before_extents, left, right) and not spans(after_extents, left, right):
        reopened.append((left, right))

left, right = 62623253, 62642316
spanning_ps = [ps for ps, (a, z) in after_extents.items() if a <= left and z >= right]
flanks = [Counter(), Counter()]
local_votes = [defaultdict(Counter), defaultdict(Counter)]
local_reads = set()
with pysam.AlignmentFile(args.bam) as bam:
    for read in bam.fetch('CHM13#0#chr20', left - 10001, right + 10000):
        name = read.query_name
        if read.is_secondary or read.is_supplementary or name not in truth:
            continue
        beg, end = read.reference_start + 1, read.reference_end
        overlapping = beg <= right and end >= left
        if overlapping:
            local_reads.add(name)
        for i, tags in enumerate((before_tags, after_tags)):
            if name not in tags:
                continue
            hp, ps = tags[name]
            bit = (hp == 1) != truth[name]
            if overlapping:
                local_votes[i][ps][bit] += 1
            if i != 1 or ps not in spanning_ps:
                continue
            if end < left and end >= left - 10000:
                flanks[0][bit] += 1
            if beg > right and beg <= right + 10000:
                flanks[1][bit] += 1

result = {
    'truth_transitions': {str(k): v for k, v in transitions.items()},
    'changed_output_only_tags': changed_output_only_tags,
    'output_only_membership_and_hp_preserved': rescue_groups_preserved,
    'newly_discordant_reads': [name for name, correct in before_status.items()
                               if correct and after_status.get(name) is False],
    'lost_keys': len(before_rows.keys() - after_rows.keys()),
    'gained_keys': len(after_rows.keys() - before_rows.keys()),
    'changed_genotype_alleles': changed_genotypes,
    'changed_nonphase_fields': changed_nonphase,
    'old_phased_keys_lost': missing_phased,
    'old_blocks_split_or_mixed_gauge': {str(ps): sorted(changes) for ps, changes in transfers.items() if len(changes) != 1},
    'newly_connected_old_gaps': new_connections,
    'reopened_tracked_gaps': reopened,
    'connected_tracked': [sum(tracked), len(tracked)],
    'target_before_spans': spans(before_extents, left, right),
    'target_after_spans': spans(after_extents, left, right),
    'target_truth_scorable_reads': len(local_reads),
    'target_scores': [{'phased': sum(sum(c.values()) for c in votes.values()),
                       'correct': sum(max(c.values()) for c in votes.values()),
                       'read_phase_sets': len(votes),
                       'dominant_correct': max((max(c.values()) for c in votes.values()), default=0),
                       'separated': max((max(c.values()) for c in votes.values()), default=0) / len(local_reads)} for votes in local_votes],
    'disjoint_flank_parent_votes': [dict(c) for c in flanks],
}
args.output.write_text(json.dumps(result, indent=2) + '\n')
print(json.dumps(result, indent=2))
# This repair corrects an existing internal switch, rather than a uniform
# union of two already correct blocks. Preserve and report the sole remaining
# newly discordant read; its high-quality physical SNP agrees with the repaired
# genotype, so its label cannot be special-cased by parental truth.
newly_discordant = [name for name, correct in before_status.items()
                    if correct and after_status.get(name) is False]
assert transitions[(True, False)] <= 1 and transitions[(True, None)] == 0
assert all(before_tags[name][1] == 62637077 for name in newly_discordant)
assert transitions[(False, True)] >= 329
assert not result['lost_keys'] and not result['changed_genotype_alleles']
assert not changed_nonphase and rescue_groups_preserved
assert not result['old_phased_keys_lost']
assert set(result['old_blocks_split_or_mixed_gauge']) == {'62637077'}
assert new_connections == [(left, right)]
assert not reopened and result['target_after_spans']
assert all(sum(c.values()) >= 5 and max(c.values()) / sum(c.values()) >= 0.95 for c in flanks)
assert flanks[0].most_common()[0][0] == flanks[1].most_common()[0][0]
