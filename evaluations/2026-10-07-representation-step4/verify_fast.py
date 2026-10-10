#!/usr/bin/env python3
"""Run focused joint-locus fixtures and reject mutations of the production adapter."""
import argparse
import hashlib
import json
from pathlib import Path
import shlex
import subprocess
import time

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--work', type=Path, default=Path('test_data/tmp_representation_step4/mutants'))
parser.add_argument('--runs', type=int, default=100)
args = parser.parse_args()
args.work.mkdir(parents=True, exist_ok=True)
out = Path(__file__).resolve().parent
source = Path('src/graph_bam_adapter.cpp').read_text()
mutations = {
    'disable_projection': ('if (projected.empty()) {', 'if (true) {'),
    'choose_conflicting_call': ('evidence.allele() < 0 ? -1 : evidence.allele()',
                              'evidence.allele() < 0 ? 0 : evidence.allele()'),
    'lose_original_profiles': ('chunk.read_var_profile = std::move(original_profiles);', ''),
    'ignore_alt_sequence': ('descriptions[*identity].push_back(ci);',
                           'descriptions[{identity->pos, identity->type, identity->ref_len, ""}].push_back(ci);'),
    'nonreference_alt_class': ('(!opts.snarl_allele_phasing || graph_chunk.site_meta[ci].alts.size() == 1)',
                             'true'),
    'ignore_read_footprint': ('if (covered) projected.push_back(&members);',
                            'projected.push_back(&members);'),
    'borrow_parent_reference': ('if (!same_contrast && parent_call != 0 && parent_call != 1) return false;',
                               ''),
}
dry_run = subprocess.check_output(['make', '--dry-run',
    '--assume-new=src/test_graph_bam_adapter.cpp', 'test_graph_bam_adapter'], text=True)
link = shlex.split(next(line for line in reversed(dry_run.splitlines())
                       if ' -o test_graph_bam_adapter ' in line))
test_object = args.work / 'test.o'
subprocess.run(['g++', '-O0', '-std=c++17', '-Isrc', '-c', 'src/test_graph_bam_adapter.cpp',
                '-o', str(test_object)], check=True)
results = {}
for name, (old, new) in mutations.items():
    assert source.count(old) == 1, name
    path = args.work / (name + '.cpp')
    path.write_text(source.replace(old, new))
    obj = args.work / (name + '.o')
    binary = args.work / name
    subprocess.run(['g++', '-O0', '-std=c++17', '-Isrc', '-c', str(path), '-o', str(obj)], check=True)
    command = [str(test_object) if arg == 'src/test_graph_bam_adapter.cpp' else
               str(obj) if arg == 'src/graph_bam_adapter.o' else arg for arg in link]
    command[command.index('-o') + 1] = str(binary)
    subprocess.run(command, check=True)
    run = subprocess.run([str(binary), '--joint-loci'], text=True, capture_output=True)
    (args.work / (name + '.log')).write_text(run.stdout + run.stderr)
    assert run.returncode != 0, name
    results[name] = dict(exit=run.returncode, failed_checks=run.stderr.count('FAIL:'))
(out / 'mutation-checks.json').write_text(json.dumps(results, indent=2) + '\n')
start = time.monotonic()
for _ in range(args.runs):
    subprocess.run(['./test_graph_bam_adapter', '--joint-loci'], stdout=subprocess.DEVNULL, check=True)
seconds = time.monotonic() - start
(out / 'fast-checks.json').write_text(json.dumps(dict(executions=args.runs,
    average_ms=1000 * seconds / args.runs, total_seconds=seconds,
    test_binary_sha256=hashlib.sha256(Path('test_graph_bam_adapter').read_bytes()).hexdigest()), indent=2) + '\n')
print(json.dumps(results, indent=2))
