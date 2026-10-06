#!/usr/bin/env python3
"""Rank current HiPhase spans and audit the next seam in the next split block."""
from collections import Counter, defaultdict
import hashlib
import json
from pathlib import Path
import pysam

ROOT = Path(__file__).resolve().parents[2]
OUT = Path(__file__).resolve().parent
PG = ROOT/'test_data/tmp_gap_fix75/retired_final/0'
HI = ROOT/'test_data/tmp_gap_fix48/competitor/hiphase_dv'
truth = {f[0]: f[1] == 'PATERNAL'
         for line in (ROOT/'test_data/derived/chr20_truth_hap.tsv').open()
         if len(f := line.rstrip().split('\t')) == 2}


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


pg_blocks, hi_blocks = blocks(PG/'phased.vcf'), blocks(HI/'phased.vcf.gz')
ranking = []
for rank, block in enumerate(hi_blocks, 1):
    overlap = sorted([b for b in pg_blocks if b['left'] <= block['right'] and b['right'] >= block['left']],
                     key=lambda b: b['left'])
    ranking.append({**block, 'hiphase_rank': rank,
                    'exact_span_covered': any(b['left'] <= block['left'] and b['right'] >= block['right'] for b in overlap),
                    'pgphase_blocks': overlap})
(OUT/'ranked-blocks.json').write_text(json.dumps(ranking, indent=2)+'\n')
# Report terminal differences in the ranking; select the largest remaining
# internal split as the next connection target.
target = next(b for b in ranking if not b['exact_span_covered'] and len(b['pgphase_blocks']) > 1)
print('SELECTED',target,flush=True)
gaps = [(a['right'], b['left']) for a, b in zip(target['pgphase_blocks'], target['pgphase_blocks'][1:])]
print('GAPS',gaps,flush=True)


def geometry(read):
    return (read.reference_start, read.reference_end, read.cigarstring,
            hashlib.sha256(read.query_sequence.encode()).hexdigest())


original, groups, unscorable = {}, defaultdict(set), Counter()
with pysam.AlignmentFile(str(ROOT/'test_data/HG002_chr20_hifi_mapped_to_CHM13_chr20_annotated.bam')) as bam:
    for read in bam.fetch('CHM13#0#chr20', target['left']-1, target['right']):
        if read.is_secondary or read.is_supplementary:
            continue
        if read.query_name not in truth:
            unscorable['block'] += 1
            for i, (left, right) in enumerate(gaps, 1):
                unscorable[f'gap_{i}'] += read.reference_start < right and read.reference_end >= left
            continue
        original[read.query_name] = geometry(read)
        groups['block'].add(read.query_name)
        for i, (left, right) in enumerate(gaps, 1):
            if read.reference_start < right and read.reference_end >= left:
                groups[f'gap_{i}'].add(read.query_name)
            if read.reference_start < left and read.reference_end >= right:
                groups[f'gap_{i}_physical_bridges'].add(read.query_name)
            if left-50000 < read.reference_end < left:
                groups[f'gap_{i}_left_flank'].add(read.query_name)
            if right <= read.reference_start < right+50000:
                groups[f'gap_{i}_right_flank'].add(read.query_name)
for i in range(1, len(gaps)+1):
    assert not (groups[f'gap_{i}_left_flank'] & groups[f'gap_{i}_right_flank'])
    assert not (groups[f'gap_{i}'] & (groups[f'gap_{i}_left_flank'] | groups[f'gap_{i}_right_flank']))

report = {'target': target, 'gaps': [{'left': left, 'right': right, 'boundary_distance_bp': right-left} for left, right in gaps],
          'inputs': {'pgphase': str(PG), 'hiphase': str(HI),
                     'production_sha256': hashlib.sha256((ROOT/'pgphase').read_bytes()).hexdigest(),
                     'pgphase_vcf_sha256': hashlib.sha256((PG/'phased.vcf').read_bytes()).hexdigest()},
          'unscorable_original_reads': dict(unscorable), 'tools': {}}
for tool, path in [('pgphase', PG/'phased.bam'), ('hiphase', HI/'phased.bam')]:
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
(OUT/'next-block.json').write_text(json.dumps(report, indent=2)+'\n')
print(json.dumps({'target': target, 'gaps': report['gaps'],
                  'regions': {tool: data['regions'] for tool, data in report['tools'].items()}}, indent=2), flush=True)
