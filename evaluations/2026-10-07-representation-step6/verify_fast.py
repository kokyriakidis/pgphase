#!/usr/bin/env python3
"""Measure saved-state scoring and prove focused tests detect production defects."""
import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import time

ROOT = Path(__file__).resolve().parents[2]
OUT = Path(__file__).resolve().parent
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--runs', type=int, default=100)
parser.add_argument('--work', type=Path, default=ROOT/'test_data/tmp_representation_step6/mutants')
args = parser.parse_args()
args.work.mkdir(parents=True, exist_ok=True)
source = (ROOT/'src/allele_context.cpp').read_text()
mutations = {
    'omit_parent_alternatives': ('add_other(sequence);', '(void)sequence;'),
    'omit_reference_decoy': ('add_other(reference);', ''),
    'accept_selected_ties': ('result.distances[nearest] < result.distances[1 - nearest]', 'result.distances[nearest] <= result.distances[1 - nearest]'),
    'ignore_reference_skip': ('if (op == BAM_CREF_SKIP && ref_pos <= context.end && ref_pos + length > context.beg)', 'if (false && op == BAM_CREF_SKIP && ref_pos <= context.end && ref_pos + length > context.beg)'),
    'erase_base_qualities': ('result.qualities.push_back(qualities[i]);', 'result.qualities.push_back(40);'),
    'shift_query_slice': ('for (int i = query_beg; i < query_end; ++i)', 'for (int i = query_beg + 1; i < query_end; ++i)'),
    'use_local_alignment': ('EDLIB_MODE_NW', 'EDLIB_MODE_HW'),
}
test_object = args.work/'test.o'
subprocess.run(['g++', '-O0', '-std=c++17', '-Isrc', '-c', 'src/test_allele_context.cpp', '-o', str(test_object)], cwd=ROOT, check=True)
results = {}
for name, (old, new) in mutations.items():
    assert source.count(old) == 1, name
    path = args.work/(name+'.cpp')
    path.write_text(source.replace(old, new))
    obj, binary = args.work/(name+'.o'), args.work/name
    subprocess.run(['g++', '-O0', '-std=c++17', '-Isrc', '-Ithird_party/edlib/edlib/include', '-c', str(path), '-o', str(obj)], cwd=ROOT, check=True)
    subprocess.run(['g++', '-o', str(binary), str(test_object), str(obj), 'src/allele_identity.o', 'src/edlib.o', '-lhts'], cwd=ROOT, check=True)
    result = subprocess.run([str(binary)], cwd=ROOT, capture_output=True, text=True)
    (args.work/(name+'.log')).write_text(result.stdout+result.stderr)
    assert result.returncode != 0 and 'FAIL:' in result.stderr, name
    results[name] = dict(exit=result.returncode, failed_checks=result.stderr.count('FAIL:'))
(OUT/'mutation-checks.json').write_text(json.dumps(results, indent=2)+'\n')
timings = {}
commands = {'fixtures': ['./test_allele_context']}
for name in ['matrix4-final', 'matrix5-final', 'matrix65-final', 'whole-final']:
    prefix = ROOT/'test_data/tmp_representation_step6/accepted'/name/'matrix.chunk0'
    commands[name] = ['./test_allele_context', '--replay', str(prefix)+'.joint-contexts.tsv', str(prefix)+'.joint-sequences.tsv']
for name, command in commands.items():
    result = subprocess.run(command, cwd=ROOT, capture_output=True, text=True, check=True)
    start = time.monotonic()
    for _ in range(args.runs):
        subprocess.run(command, cwd=ROOT, stdout=subprocess.DEVNULL, check=True)
    timings[name] = dict(executions=args.runs, average_ms=1000*(time.monotonic()-start)/args.runs,
                         result=result.stdout.strip())
report = dict(test_binary_sha256=hashlib.sha256((ROOT/'test_allele_context').read_bytes()).hexdigest(), timings=timings)
(OUT/'fast-checks.json').write_text(json.dumps(report, indent=2)+'\n')
print(json.dumps(dict(mutations=results, fast_checks=report), indent=2))
