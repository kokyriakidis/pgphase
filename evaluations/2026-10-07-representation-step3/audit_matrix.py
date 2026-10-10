#!/usr/bin/env python3
"""Audit molecule reduction against the frozen step-2 owning chunk."""
from collections import Counter
import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
BEFORE = ROOT/'test_data/tmp_representation_step2/matrix5/matrix.chunk0.bam-overlay-output.tsv'
AFTER = ROOT/'test_data/tmp_representation_step3/matrix5/matrix.chunk0.bam-overlay-output.tsv'

def matrix(path):
    result = {kind: {} for kind in ('VAR', 'READ', 'OBS', 'BAMQ')}
    for line in path.read_text().splitlines():
        fields = line.split('\t')
        kind = fields[0]
        if kind == 'VAR':
            result[kind][int(fields[1])] = fields[2:]
        elif kind == 'READ':
            result[kind][fields[1]] = fields[2:]
        elif kind == 'OBS':
            result[kind][(fields[1], int(fields[2]))] = tuple(map(int, fields[3:6]))
        elif kind == 'BAMQ':
            result[kind][(fields[1], int(fields[2]))] = int(fields[3])
    return result

before, after = matrix(BEFORE), matrix(AFTER)
assert before['VAR'] == after['VAR'], 'identity aliases must not alter graph candidates'
assert before['READ'] == after['READ'], 'overlay must not import source read gauges'
changes = Counter()
examples = []
for key in sorted(before['OBS'].keys() | after['OBS'].keys()):
    old = before['OBS'].get(key, (-1, -1, -1))
    new = after['OBS'].get(key, (-1, -1, -1))
    assert old[1] == new[1], 'independent graph evidence changed'
    if old[2] < 0 <= new[2]:
        assert new[1] in (0, 1), 'normalized call borrowed a source REF class without independent graph evidence'
    if old != new:
        changes['changed_observations'] += 1
        changes['changed_effective_calls'] += old[0] != new[0]
        changes['added_bam_calls'] += old[2] < 0 <= new[2]
        changes['added_graph_bam_disagreements'] += old[2] < 0 <= new[2] and new[1] >= 0 and new[1] != new[2]
        changes['withheld_bam_calls'] += new[2] < 0 <= old[2]
        changes['changed_known_bam_calls'] += old[2] >= 0 and new[2] >= 0 and old[2] != new[2]
        if len(examples) < 10:
            examples.append(dict(read=key[0], candidate=key[1], position=after['VAR'][key[1]][0], before=old, after=new))
assert all(before['BAMQ'].get(key, 0) == after['BAMQ'].get(key, 0) for key in before['BAMQ'].keys() | after['BAMQ'].keys()), 'overlay invented or displaced SNP base-quality certificates'
report = dict(before=str(BEFORE.relative_to(ROOT)), after=str(AFTER.relative_to(ROOT)), before_sha256=hashlib.sha256(BEFORE.read_bytes()).hexdigest(), after_sha256=hashlib.sha256(AFTER.read_bytes()).hexdigest(), candidates=len(after['VAR']), reads=len(after['READ']), candidate_state_unchanged=True, read_state_unchanged=True, graph_channel_unchanged=True, snp_certificates_unchanged=True, changes=dict(changes), examples=examples)
(Path(__file__).resolve().parent/'matrix-checks.json').write_text(json.dumps(report, indent=2)+'\n')
print(json.dumps(report, indent=2))
