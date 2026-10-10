#!/usr/bin/env python3
"""Verify complete contrast fixtures and mutations without rebuilding production."""
import argparse
import hashlib
import json
from pathlib import Path
import shlex
import subprocess
import time

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--work', type=Path, default=Path('test_data/tmp_representation_step5/mutants'))
parser.add_argument('--runs', type=int, default=100)
args = parser.parse_args()
args.work.mkdir(parents=True, exist_ok=True)
out = Path(__file__).resolve().parent
sources = {name: Path('src', name + '.cpp').read_text() for name in ('allele_identity', 'graph_bam_adapter')}
mutations = {
    'discard_first_allele': ('allele_identity',
        'if (key.alleles[0] == key.alleles[1]) return std::nullopt;',
        'key.alleles[0].reset();\n    if (key.alleles[0] == key.alleles[1]) return std::nullopt;'),
    'ignore_alt_order': ('graph_bam_adapter',
        'if (contrast->reversed) std::swap(local[0], local[1]);', ''),
    'invent_reference_class': ('graph_bam_adapter',
        ': call == contrast->local_alleles[0] ? 0', ': call == 0 ? 0 : call == contrast->local_alleles[0] ? 0'),
    'prefer_primary': ('graph_bam_adapter',
        'if (profile.graph_alleles.empty())\n            add(profile.alleles, profile.alt_qi);\n        else\n            add(profile.graph_alleles, {});\n        add(profile.bam_alleles, profile.bam_qi);',
        'add(profile.alleles, profile.alt_qi);'),
    'ignore_other_alt_class': ('graph_bam_adapter', 'if (meta.non_selected_alt_class) continue;', ''),
    'erase_existing_conflict': ('allele_identity',
        'if (allele != 0 && allele != 1 && allele != kConflictingBamAllele) return;',
        'if (allele != 0 && allele != 1) return;'),
    'omit_class_provenance': ('graph_bam_adapter', 'new_meta.back().non_selected_alt_class = multi;', ''),
    'ignore_original_walk_mapping': ('graph_bam_adapter',
        'const int source = original[static_cast<size_t>(allele)];', 'const int source = allele;'),
}
dry = subprocess.check_output(['make', '--dry-run', '--assume-new=src/test_graph_bam_adapter.cpp', 'test_graph_bam_adapter'], text=True)
link = shlex.split(next(line for line in reversed(dry.splitlines()) if ' -o test_graph_bam_adapter ' in line))
test_object = args.work / 'test.o'
subprocess.run(['g++', '-O0', '-std=c++17', '-Isrc', '-c', 'src/test_graph_bam_adapter.cpp', '-o', str(test_object)], check=True)
results = {}
for name, (component, old, new) in mutations.items():
    source = sources[component]
    assert source.count(old) == 1, name
    path = args.work / (name + '.cpp')
    path.write_text(source.replace(old, new))
    obj, binary = args.work / (name + '.o'), args.work / name
    subprocess.run(['g++', '-O0', '-std=c++17', '-Isrc', '-c', str(path), '-o', str(obj)], check=True)
    command = [str(test_object) if arg == 'src/test_graph_bam_adapter.cpp' else
               str(obj) if arg == 'src/' + component + '.o' else arg for arg in link]
    command[command.index('-o') + 1] = str(binary)
    subprocess.run(command, check=True)
    run = subprocess.run([str(binary)] + ([] if name == 'omit_class_provenance' else ['--joint-contrasts']), text=True, capture_output=True)
    (args.work / (name + '.log')).write_text(run.stdout + run.stderr)
    assert run.returncode != 0, name
    results[name] = dict(exit=run.returncode, failed_checks=run.stderr.count('FAIL:'))
(out / 'mutation-checks.json').write_text(json.dumps(results, indent=2) + '\n')
timings = {}
for binary, option in [('./test_graph_bam_adapter', '--joint-contrasts'), ('./test_allele_identity', '--contrasts')]:
    start = time.monotonic()
    for _ in range(args.runs):
        subprocess.run([binary, option], stdout=subprocess.DEVNULL, check=True)
    elapsed = time.monotonic() - start
    timings[option] = dict(executions=args.runs, average_ms=1000 * elapsed / args.runs,
        test_binary_sha256=hashlib.sha256(Path(binary).read_bytes()).hexdigest())
(out / 'fast-checks.json').write_text(json.dumps(timings, indent=2) + '\n')
print(json.dumps(results, indent=2))
print(json.dumps(timings, indent=2))
