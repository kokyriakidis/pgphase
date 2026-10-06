#!/usr/bin/env python3
"""Persist successful gap replays across test processes; always rerun assertions."""

import argparse
import fcntl
import hashlib
import json
import os
from pathlib import Path
import shlex
import shutil
import subprocess


OUTPUTS = ('candidates.tsv', 'native.vcf', 'phased.bam', 'stdout.log', 'stderr.log')
INPUT_OPTIONS = ('--ref', '--bam', '--sites', '--gaf')
INDEX_SUFFIXES = ('.fai', '.gzi', '.bai', '.csi', '.tbi')
SHELL_OPERATORS = ('&&', ';', '>', '2>')


def replay_command(command, outdir, replacement, matrix_prefix):
    # Quote arguments again after substitution, retaining the fixed shell
    # operators used by run_arm (including its redirected diagnostic logs).
    tokens = shlex.split(command)
    if '--phase-matrix-dump' in tokens:
        tokens[tokens.index('--phase-matrix-dump') + 1] = matrix_prefix
    return ' '.join(token if token in SHELL_OPERATORS else
                    shlex.quote(token.replace(outdir, replacement))
                    for token in tokens)


def identity(path):
    path = Path(path).resolve(strict=True)
    stat = path.stat()
    return [str(path), stat.st_dev, stat.st_ino, stat.st_size,
            stat.st_mtime_ns, stat.st_ctime_ns]


def signature(command, outdir):
    tokens = shlex.split(command)
    executable = Path(tokens[tokens.index('collect-graph-variation') - 1]).resolve()
    inputs = []
    for option in INPUT_OPTIONS:
        path = Path(tokens[tokens.index(option) + 1])
        inputs.append(identity(path))
        for suffix in INDEX_SUFFIXES:
            index = Path(str(path) + suffix)
            if index.exists():
                inputs.append(identity(index))
    # Dynamic libraries affect execution even if the pgphase binary is unchanged.
    libraries = subprocess.check_output(['ldd', str(executable)], text=True)
    for line in libraries.splitlines():
        fields = line.split()
        path = fields[2] if '=>' in fields and len(fields) > 2 else fields[0]
        if path.startswith('/'):
            inputs.append(identity(path))
    payload = {
        'schema': 1,
        'helper': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        'command': replay_command(command, outdir, '<output>', '<matrix>'),
        'binary': hashlib.sha256(executable.read_bytes()).hexdigest(),
        'inputs': inputs,
        'environment': {key: value for key, value in os.environ.items()
                        if key.startswith('PGPHASE_') and key not in
                        ('PGPHASE_TEST_WORKDIR', 'PGPHASE_TEST_CACHE', 'PGPHASE_BIN')},
        'ld_library_path': os.environ.get('LD_LIBRARY_PATH', ''),
    }
    encoded = json.dumps(payload, sort_keys=True).encode()
    return hashlib.sha256(encoded).hexdigest(), payload


def output_state(directory):
    for name in OUTPUTS:
        identity(directory / name)
    return {path.name: identity(path)[1:] for path in directory.iterdir()
            if path.is_file() and path.name not in ('state.json', 'state.tmp')}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--command', required=True)
    parser.add_argument('--outdir', required=True)
    args = parser.parse_args()
    key, payload = signature(args.command, args.outdir)
    root = Path(os.environ.get('PGPHASE_TEST_CACHE',
                              'test_data/.gap-replay-cache')).resolve()
    root.mkdir(parents=True, exist_ok=True)
    directory = root / key
    # Hold the lock through publication and copying: partial and failed runs
    # never supply a completed state, including when multiple shards collide.
    with (root / (key + '.lock')).open('w') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX)
        state_path = directory / 'state.json'
        try:
            state = json.loads(state_path.read_text())
            hit = state['outputs'] == output_state(directory)
        except (OSError, ValueError, KeyError):
            hit = False
        if not hit:
            directory.mkdir(exist_ok=True)
            state_path.unlink(missing_ok=True)
            for path in directory.iterdir():
                if path.is_file():
                    path.unlink()
            result = subprocess.run(replay_command(args.command, args.outdir,
                                                   str(directory),
                                                   str(directory / 'matrix')),
                                    shell=True, check=False)
            if result.returncode:
                # Keep the case's failure logs available for Catch2 diagnostics.
                Path(args.outdir).mkdir(parents=True, exist_ok=True)
                for name in OUTPUTS:
                    if (directory / name).exists():
                        shutil.copyfile(directory / name, Path(args.outdir) / name)
                return result.returncode
            if any((directory / name).stat().st_size == 0
                   for name in ('native.vcf', 'phased.bam')):
                raise RuntimeError(f'replay produced empty outputs: {directory}')
            state = {'signature': payload, 'outputs': output_state(directory)}
            temporary = directory / 'state.tmp'
            temporary.write_text(json.dumps(state, indent=2) + '\n')
            temporary.replace(state_path)
        target = Path(args.outdir)
        target.mkdir(parents=True, exist_ok=True)
        for name in OUTPUTS:
            destination = target / name
            # A case may still hold hard links from the old in-process cache.
            destination.unlink(missing_ok=True)
            shutil.copyfile(directory / name, destination)
        tokens = shlex.split(args.command)
        if '--phase-matrix-dump' in tokens:
            prefix = tokens[tokens.index('--phase-matrix-dump') + 1]
            Path(prefix).parent.mkdir(parents=True, exist_ok=True)
            for path in directory.glob('matrix*'):
                shutil.copyfile(path, prefix + path.name[len('matrix'):])
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
