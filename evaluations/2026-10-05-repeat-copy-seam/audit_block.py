#!/usr/bin/env python3
"""Audit the last seam in the third-largest HiPhase block on identical original reads."""
from collections import Counter, defaultdict
import hashlib
import importlib.util
import json
from pathlib import Path
import pysam

ROOT = Path(__file__).resolve().parents[2]
OUT = Path(__file__).resolve().parent
BEFORE = ROOT/'test_data/tmp_gap_fix72/frozen_final/0'
PG = ROOT/'test_data/tmp_gap_fix73/frozen_final/0'
HI = ROOT/'test_data/tmp_gap_fix48/competitor/hiphase_dv'
truth = {f[0]: f[1] == 'PATERNAL'
         for line in (ROOT/'test_data/derived/chr20_truth_hap.tsv').open()
         if len(f := line.rstrip().split('\t')) == 2}
spec = importlib.util.spec_from_file_location('audit', OUT/'audit_preservation.py')
audit = importlib.util.module_from_spec(spec)
spec.loader.exec_module(audit)
target = {'left': 24103779, 'right': 25944471, 'length': 1840693, 'hiphase_rank': 3}
left, right = 25855631, 25855633

def blocks(path):
    positions = defaultdict(list)
    with pysam.VariantFile(str(path)) as vcf:
        for row in vcf:
            sample = next(iter(row.samples.values()))
            gt, ps = sample.get('GT'), sample.get('PS')
            if sample.phased and gt and None not in gt and len(set(gt)) > 1 and ps and ps > 0:
                positions[ps].append(row.pos)
    return sorted([{'ps': ps, 'left': min(pos), 'right': max(pos),
                    'length': max(pos)-min(pos)+1, 'rows': len(pos)}
                   for ps, pos in positions.items() if len(pos) > 1],
                  key=lambda b: b['length'], reverse=True)


def geometry(read):
    return (read.reference_start, read.reference_end, read.cigarstring,
            hashlib.sha256(read.query_sequence.encode()).hexdigest())


original, groups, unscorable = {}, defaultdict(set), Counter()
with pysam.AlignmentFile(str(ROOT/'test_data/HG002_chr20_hifi_mapped_to_CHM13_chr20_annotated.bam')) as bam:
    for read in bam.fetch('CHM13#0#chr20', target['left']-1, target['right']):
        if read.is_secondary or read.is_supplementary:
            continue
        in_gap = read.reference_start < right and read.reference_end >= left
        if read.query_name not in truth:
            unscorable['block'] += 1
            unscorable['gap'] += in_gap
            continue
        original[read.query_name] = geometry(read)
        groups['block'].add(read.query_name)
        if in_gap:
            groups['gap'].add(read.query_name)
        if read.reference_start < 24131707 and read.reference_end >= 24121713:
            groups['other_gap'].add(read.query_name)
        if left-50000 < read.reference_end < left:
            groups['left_flank'].add(read.query_name)
        if right <= read.reference_start < right+50000:
            groups['right_flank'].add(read.query_name)
assert not (groups['left_flank'] & groups['right_flank'])
assert not (groups['gap'] & (groups['left_flank'] | groups['right_flank']))

report = {'target': target, 'gap': [left, right], 'gap_boundary_distance_bp': right-left,
          'unscorable_original_reads': dict(unscorable), 'tools': {}}
for tool, path in [('before', BEFORE/'phased.bam'), ('pgphase', PG/'phased.bam'), ('hiphase', HI/'phased.bam')]:
    votes, tags, matched = defaultdict(Counter), {}, set()
    with pysam.AlignmentFile(str(path), check_sq=False) as bam:
        for read in bam:
            if read.is_secondary or read.is_supplementary:
                continue
            hp = read.get_tag('HP') if read.has_tag('HP') else 0
            ps = read.get_tag('PS') if read.has_tag('PS') else 0
            q = read.query_name
            if q in truth and hp in (1, 2) and ps > 0:
                votes[ps][(hp == 1) != truth[q]] += 1
            if q in original:
                tags[q] = hp, ps
                if tool == 'hiphase':
                    assert geometry(read) == original[q], q
                    matched.add(q)
    if tool == 'hiphase':
        assert matched == original.keys()
    orientation = {ps: c.most_common(1)[0][0] for ps, c in votes.items()}
    regions = {}
    for name, names in groups.items():
        counts, cores, local_votes = Counter(), Counter(), defaultdict(Counter)
        for q in names:
            hp, ps = tags.get(q, (0, 0))
            state = 'unphased' if hp not in (1, 2) or ps <= 0 else (
                'correct' if ((hp == 1) != truth[q]) == orientation[ps] else 'discordant')
            counts[state] += 1
            if state == 'correct' and ps < 1_000_000_000:
                cores[ps] += 1
            if hp in (1, 2) and ps > 0:
                local_votes[ps][(hp == 1) != truth[q]] += 1
        regions[name] = {'scorable': len(names), 'counts': dict(counts),
                         'correct_all_reads': counts['correct']/len(names),
                         'dominant_connected_correct_core': max(cores.values(), default=0),
                         'correct_cores': dict(cores),
                         'parental_orientation_votes': {ps: dict(v) for ps, v in local_votes.items()}}
    report['tools'][tool] = {'regions': regions, 'identical_original_alignments': len(matched) if tool == 'hiphase' else None,
                             'phase_set_orientation': {ps: dict(votes[ps]) for ps in {ps for hp, ps in tags.values()} if ps > 0}}
report['preservation'] = audit.compare(BEFORE, PG, truth)
report['metrics'] = {}
for name, path in [('before', BEFORE/'phased.vcf'), ('pgphase', PG/'phased.vcf'), ('hiphase', HI/'phased.vcf.gz')]:
    bs = blocks(path)
    lengths = sorted([b['length'] for b in bs], reverse=True)
    running = 0
    n50 = 0
    for length in lengths:
        running += length
        if 2 * running >= sum(lengths):
            n50 = length
            break
    report['metrics'][name] = {'block_count': len(bs), 'n50_bp': n50,
        'largest_block': bs[0],
        'target_blocks': [b for b in bs if b['left'] <= target['right'] and b['right'] >= target['left']]}
fixed = report['tools']['pgphase']['regions']
competitor = report['tools']['hiphase']['regions']
assert fixed['gap']['correct_all_reads'] >= 0.8
assert fixed['gap']['dominant_connected_correct_core'] >= 88
assert fixed['block']['dominant_connected_correct_core'] >= competitor['block']['dominant_connected_correct_core']
assert fixed['other_gap']['dominant_connected_correct_core'] >= 93
joined = [b for b in report['metrics']['pgphase']['target_blocks']
          if b['left'] <= target['left'] and b['right'] >= target['right']]
assert len(joined) == 1, joined
core = joined[0]['ps']
flanks = [fixed[name]['parental_orientation_votes'].get(core, {}) for name in ('left_flank', 'right_flank')]
assert all(sum(v.values()) >= 5 for v in flanks), flanks
orientation = max(flanks[0], key=flanks[0].get)
assert all(v.get(orientation, 0)/sum(v.values()) >= 0.9 for v in flanks), flanks
(OUT/'results.json').write_text(json.dumps(report, indent=2)+'\n')
print(json.dumps({'metrics': report['metrics'], 'preservation': report['preservation'],
    'gap': {name: result['regions']['gap'] for name, result in report['tools'].items()}}, indent=2))
