#!/usr/bin/env python3
"""Measure the fast fixtures and prove they reject broken molecule reduction."""
import argparse
import json
from pathlib import Path
import subprocess
import time

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--work', type=Path, default=Path('test_data/tmp_representation_step3'))
parser.add_argument('--runs', type=int, default=100)
args = parser.parse_args()
out = Path(__file__).resolve().parent
args.work.mkdir(parents=True, exist_ok=True)
source = Path('src/allele_identity.cpp').read_text()
mutations = {
    'ignore_conflict': ('} else if (allele_ != allele) {', '} else if (allele_ != allele && false) {'),
    'borrow_query_index': ('} else if (query_index_ != query_index) {', '} else if (query_index_ != query_index && false) {'),
    'repeat_provenance': ('previous.query_index == query_index) return;', 'previous.query_index == query_index) break;'),
}
results = {}
for name, (old, new) in mutations.items():
    assert source.count(old) == 1
    path = args.work / (name + '.cpp')
    path.write_text(source.replace(old, new))
    binary = args.work / name
    subprocess.run(['g++', '-O2', '-std=c++17', '-Isrc', 'src/test_allele_identity.cpp',
                    str(path), '-o', str(binary)], check=True)
    run = subprocess.run([str(binary)], text=True, capture_output=True)
    (args.work / (name + '.log')).write_text(run.stdout + run.stderr)
    assert run.returncode != 0, name
    results[name] = {'exit': run.returncode, 'failed_checks': run.stderr.count('FAIL:')}
(out / 'mutation-checks.json').write_text(json.dumps(results, indent=2) + '\n')
start = time.monotonic()
for _ in range(args.runs):
    subprocess.run(['./test_allele_identity'], stdout=subprocess.DEVNULL, check=True)
seconds = time.monotonic() - start
(out / 'fast-checks.json').write_text(json.dumps({
    'executions': args.runs, 'total_seconds': seconds,
    'average_ms': seconds * 1000 / args.runs, 'molecule_permutations_per_execution': 42,
    'complete_haplotype_fixtures_per_execution': 3328,
}, indent=2) + '\n')
print(json.dumps(results, indent=2))
