#!/usr/bin/env python3
"""Replay the existing exact MEC kernel with fixed SNP gauges and gap variables."""
import argparse
import hashlib
import importlib.util
import json
import subprocess
import time
from collections import Counter, defaultdict
from pathlib import Path

import pysam

ROOT = Path(__file__).resolve().parents[2]
OUT = Path(__file__).resolve().parent
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--replays', type=Path, default=ROOT/'test_data/tmp_joint_evidence_investigation')
parser.add_argument('--snp-priority', action='store_true', help='Make one fixed SNP mismatch cost more than all variable mismatches')
args = parser.parse_args()
source = (ROOT/'src/collect_phase.cpp').read_text()
kernel = source[source.index('struct RecoveryMecRow {'):source.index('// Solve one adjacent recovery edge')]
driver = (ROOT/'evaluations/2026-10-02-joint-matrix-investigation/build_replay_solver.py').read_text()
main = driver.split("main=r'''", 1)[1].split("'''", 1)[0]
cpp = args.replays/'mec_replay.cpp'
cpp.write_text('#include <algorithm>\n#include <array>\n#include <climits>\n#include <functional>\n#include <iostream>\n#include <map>\n#include <vector>\n'+kernel+main)
binary = args.replays/'mec_replay'
subprocess.run(['g++', '-O3', '-std=c++17', '-Wall', '-Wextra', str(cpp), '-o', str(binary)], check=True)
spec = importlib.util.spec_from_file_location('read_audit', ROOT/'evaluations/2026-10-06-shared-insertion-source/read_audit.py')
audit = importlib.util.module_from_spec(spec)
spec.loader.exec_module(audit)
truth = {f[0]: f[1] == 'PATERNAL' for line in (ROOT/'test_data/derived/chr20_truth_hap.tsv').open()
         if len(f := line.rstrip().split('\t')) == 2}
hi_tags, hi_status = json.loads((ROOT/'test_data/tmp_gap_fix78/hi-tags.json').read_text())
out = []
for chunk, left, right in [(4, 4766928, 4792960), (5, 5309406, 5345085), (6, 6513891, 6516221)]:
    directory = args.replays/str(chunk)
    tags, status = audit.assignments(directory/'phased.bam', truth)
    with pysam.AlignmentFile(str(ROOT/'test_data/HG002_chr20_hifi_mapped_to_CHM13_chr20_annotated.bam')) as bam:
        names = {r.query_name for r in bam.fetch('CHM13#0#chr20', left-1, right)
                 if not r.is_secondary and not r.is_supplementary and r.query_name in truth}
    core = Counter(tags[q][1] for q in names if q in tags and tags[q][0] in (1, 2) and 0 < tags[q][1] < 1_000_000_000).most_common(1)[0][0]
    core_votes = Counter((hp == 1) != truth[q] for q, (hp, ps) in tags.items() if q in truth and ps == core and hp in (1, 2))
    core_orientation = core_votes[True] >= core_votes[False]
    fixed_positions = {}
    with pysam.VariantFile(str(directory/'phased.vcf')) as vcf:
        for row in vcf:
            s = next(iter(row.samples.values()))
            if len(row.ref) == 1 and row.alts and len(row.alts) == 1 and len(row.alts[0]) == 1 and s.phased and s.get('PS') == core and set(s.get('GT', ())) == {0, 1}:
                fixed_positions[row.pos] = s['GT'][0]
    variants, observations, skipped = {}, defaultdict(dict), set()
    for line in (directory/'matrix.chunk0.bam-overlay-output.tsv').open():
        f = line.rstrip().split('\t')
        if f[0] == 'VAR':
            variants[int(f[1])] = {'pos': int(f[2]), 'type': f[3], 'flags': int(f[4]), 'ps': int(f[8]), 'h1': int(f[9]), 'h2': int(f[10]), 'af': float(f[16]), 'alt': f[7]}
        elif f[0] == 'READ' and int(f[5]):
            skipped.add(f[1])
        elif f[0] == 'OBS':
            observations[f[1]][int(f[2])] = tuple(map(int, f[3:6]))
    fixed = {i: fixed_positions[v['pos']] for i, v in variants.items() if v['type'] == 'X' and v['flags'] == 4 and v['pos'] in fixed_positions}
    beg = max(variants[i]['pos'] for i in fixed if variants[i]['pos'] <= left)
    end = min(variants[i]['pos'] for i in fixed if variants[i]['pos'] >= right)
    for channel in ['effective', 'agreement_union']:
        def call(values):
            a, g, b = values
            if channel == 'effective':
                return a
            return -1 if g in (0, 1) and b in (0, 1) and g != b else g if g in (0, 1) else b
        base_profiles = {q: {i: a for i, values in p.items() if (a := call(values)) in (0, 1)} for q, p in observations.items() if q not in skipped}
        seen = defaultdict(set)
        nonbinary = {i for q, p in observations.items() if q not in skipped
                     for i, values in p.items() if call(values) > 1}
        for p in base_profiles.values():
            for i, a in p.items():
                seen[i].add(a)
        for selection in ['centered', 'all']:
            variables = sorted(i for i, v in variants.items() if i not in fixed and i not in nonbinary and beg <= v['pos'] <= end and seen[i] == {0, 1}
                               and not (v['ps'] > 0 and v['h1'] == v['h2'] and v['h1'] >= 0)
                               and (selection == 'all' or abs(v['af']-.5) <= .12))
            record = {'window': f'{left}-{right}', 'channel': channel, 'selection': selection, 'variable_scope': [beg, end],
                      'variables': [dict(index=i, **variants[i]) for i in variables], 'fixed_core_snps': len(fixed), 'folds': []}
            if len(variables) > 20:
                record['resource_limit'] = True
                out.append(record)
                continue
            variable_index = {i: j for j, i in enumerate(variables)}
            selected = set(fixed) | set(variables)
            profiles = {q: {i: a for i, a in p.items() if i in selected} for q, p in base_profiles.items()}
            snp_weight = 1 + sum(sum(i in variable_index for i in p) for p in profiles.values()) if args.snp_priority else 1
            record['snp_weight'] = snp_weight
            def fold(q):
                value = 14695981039346656037
                for byte in q.encode():
                    value = ((value ^ byte) * 1099511628211) & ((1 << 64)-1)
                return value & 1
            for half in [-1, 0, 1]:
                rows = {q: p for q, p in profiles.items() if p and (half < 0 or fold(q) == half)}
                encoded = [f'{len(rows)} {len(variables)}']
                for p in rows.values():
                    mismatches = snp_weight * sum(a != fixed[i] for i, a in p.items() if i in fixed)
                    observations_count = snp_weight * sum(i in fixed for i in p)
                    costs = [(variable_index[i], a, 1-a) for i, a in p.items() if i in variable_index]
                    encoded.append(f'{mismatches} {observations_count} {len(costs)} '+' '.join(f'{j} {c0} {c1}' for j, c0, c1 in costs))
                started = time.monotonic()
                solved = subprocess.run([str(binary)], input='\n'.join(encoded)+'\n', text=True, capture_output=True, check=True, timeout=30)
                values = list(map(int, solved.stdout.split()))
                calls = dict(fixed)
                calls.update({i: values[j+1] for i, j in variable_index.items()})
                counts, new_counts = Counter(), Counter()
                core_losses = 0
                for q in names:
                    p = profiles.get(q, {})
                    vote = Counter()
                    for i, a in p.items():
                        vote[a == calls[i]] += snp_weight if i in fixed else 1
                    hp = 0 if vote[True] == vote[False] else 1 if vote[True] > vote[False] else 2
                    predicted = 'unphased' if not hp else 'correct' if ((hp == 1) != truth[q]) == core_orientation else 'discordant'
                    counts[predicted] += 1
                    if tags.get(q, (0, 0))[1] != core or status.get(q) != 'correct':
                        new_counts[predicted] += 1
                    elif predicted != 'correct':
                        core_losses += 1
                record['folds'].append({'fold': half, 'mec_score': values[0], 'bits': values[1:], 'rows': len(rows),
                                        'seconds': time.monotonic()-started, 'diagnostic_counts': dict(counts),
                                        'noncorrect_or_noncore_predictions': dict(new_counts), 'old_correct_core_losses': core_losses})
            out.append(record)
            print(record['window'], channel, selection, len(variables), [(r['fold'], r['mec_score'], r['diagnostic_counts'], r['old_correct_core_losses']) for r in record['folds']], flush=True)
report = {'kernel_sha256': hashlib.sha256(kernel.encode()).hexdigest(), 'production_sha256': hashlib.sha256((ROOT/'pgphase').read_bytes()).hexdigest(), 'replays': out}
(OUT/('joint-replay-snp-priority.json' if args.snp_priority else 'joint-replay.json')).write_text(json.dumps(report, indent=2)+'\n')
