#!/usr/bin/env python3
"""Score terminal counterfactuals on identical original reads, including abstentions."""
from collections import Counter, defaultdict
import importlib.util
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
OUT = Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location(
    'audit_reads', ROOT/'evaluations/2026-10-05-largest-hiphase-block/audit_reads.py')
audit = importlib.util.module_from_spec(spec)
spec.loader.exec_module(audit)
truth = {f[0]: f[1] == 'PATERNAL'
         for line in (ROOT/'test_data/derived/chr20_truth_hap.tsv').open()
         if len(f := line.rstrip().split('\t')) == 2}
reads = json.loads((OUT/'terminal-reads.json').read_text())['reads']
base_tags, base = audit.assignments(
    ROOT/'test_data/tmp_gap_fix70/baseline/66/phased.bam', truth)
result = {}
for stage in ['baseline/66', 'mapq1/66', 'permissive/66']:
    path = ROOT/'test_data/tmp_gap_fix70'/stage
    tags, status = audit.assignments(path/'phased.bam', truth)
    regions = {}
    for name, beg in [('tail', 66194204), ('terminal', 66206480)]:
        names = [r['qname'] for r in reads
                 if r['truth'] and r['end'] >= beg and r['start'] <= 66206480]
        core = Counter(tags[q][1] for q in names
                       if status.get(q) == 'correct' and tags[q][1] < 1_000_000_000)
        regions[name] = {
            'scorable': len(names),
            'statuses': dict(Counter(status.get(q, 'unphased') for q in names)),
            'dominant_correct_core': max(core.values(), default=0),
            'correct_cores': dict(core)}
    blocks = defaultdict(list)
    for line in (path/'phased.vcf').open():
        if line.startswith('#'):
            continue
        f = line.rstrip().split('\t')
        sample = dict(zip(f[8].split(':'), f[9].split(':')))
        gt, ps = sample.get('GT', ''), int(sample.get('PS', '0'))
        if '|' in gt and len(set(gt.split('|'))) > 1 and ps > 0:
            blocks[ps].append(int(f[1]))
    lost = [q for q in base if base[q] == 'correct' and status.get(q) != 'correct']
    result[stage] = {
        'regions': regions,
        'blocks_ending_in_tail': {ps: [min(pos), max(pos), len(pos)]
                                  for ps, pos in blocks.items() if max(pos) > 66194203},
        'preservation': {
            'old_correct': Counter(base.values())['correct'],
            'new_correct': Counter(status.values())['correct'],
            'lost_old_correct': len(lost),
            'lost_old_connected_correct': sum(base_tags[q][1] < 1_000_000_000 for q in lost),
            'newly_correct': sum(state == 'correct' and base.get(q) != 'correct'
                                 for q, state in status.items())}}
assert result['baseline/66']['regions']['tail']['statuses']['correct'] == 25
assert result['mapq1/66']['regions']['tail']['statuses']['correct'] == 32
assert result['mapq1/66']['preservation']['lost_old_correct'] == 74
assert result['permissive/66']['regions']['tail']['statuses']['correct'] == 31
(OUT/'counterfactuals.json').write_text(json.dumps(result, indent=2)+'\n')
print(json.dumps(result, indent=2))
