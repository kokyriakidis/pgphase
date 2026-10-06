#!/usr/bin/env python3
"""Verify the complete chromosome closes only the reviewed VCF interval."""
from collections import defaultdict
import json
from pathlib import Path


def extents(path):
    blocks = defaultdict(list)
    with open(path) as vcf:
        for line in vcf:
            if line.startswith('#'):
                continue
            fields = line.rstrip().split('\t')
            sample = dict(zip(fields[8].split(':'), fields[9].split(':')))
            gt, ps = sample.get('GT', ''), sample.get('PS', '.')
            if '|' not in gt or len(set(gt.split('|'))) != 2 or ps in ('.', '0', '-1'):
                continue
            blocks[int(ps)].append(int(fields[1]))
    return sorted((min(positions), max(positions), ps) for ps, positions in blocks.items())


before = extents('test_data/tmp_gap_fix65/frozen_full/0/phased.vcf')
after = extents('test_data/tmp_gap_fix67/frozen_full/0/phased.vcf')
assert all(any(a <= left and b >= right for a, b, _ in after) for left, right, _ in before)
covered = 0
new_closures = []
for left, right, ps in before:
    if covered and left > covered:
        spans = [nps for a, b, nps in after if a <= covered and b >= left]
        if spans:
            new_closures.append({'left': covered, 'right': left, 'phase_sets': spans})
    covered = max(covered, right)
assert [(gap['left'], gap['right']) for gap in new_closures] == [(62408056, 62432427)], new_closures
result = {'region': 'CHM13#0#chr20:1-66210255',
          'variant_extent_baseline': 'test_data/tmp_gap_fix65/frozen_full/0/phased.vcf',
          'baseline_binary': 'f9efac33e864064418fe09e70c4502b4db43dae48f6825bc4e9556d699e6ed1b',
          'baseline_note': 'The immediately preceding masked-SNP fix added one read tag and changed no VCF records.',
          'final_binary': 'f98aead0dfc6f9e72a2e2cb99de3463234523d5be7437859cecb095a2259e6fd',
          'all_previous_variant_extents_preserved': True, 'new_closures': new_closures}
Path('evaluations/2026-10-05-shared-deletion-bridge/full-closures.json').write_text(json.dumps(result, indent=2) + '\n')
print(json.dumps(result, indent=2))
