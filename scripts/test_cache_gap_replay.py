#!/usr/bin/env python3
"""Check replay invalidation, failed runs, and concurrent shard reuse."""

import os
from pathlib import Path
import shlex
import shutil
import subprocess
import sys
import tempfile
import unittest

import cache_gap_replay


class ReplayCacheTest(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory(prefix="gap-cache-'quote-")
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.binary = self.root / 'pgphase'
        shutil.copyfile('/bin/true', self.binary)
        self.binary.chmod(0o755)
        self.inputs = [self.root / name for name in ('ref', 'bam', 'sites', 'gaf')]
        for path in self.inputs:
            path.write_text('input')
        self.counter = self.root / 'calls'
        self.producer = self.root / 'producer.py'
        self.producer.write_text('''import pathlib, sys, time
out, counter, fail = map(pathlib.Path, sys.argv[1:4])
with counter.open('a') as f: f.write('run\\n')
time.sleep(0.1)
for name in ('candidates.tsv', 'native.vcf', 'phased.bam', 'stdout.log', 'stderr.log'):
 (out / name).write_text('result')
if len(sys.argv) > 4: (out / 'matrix.chunk0.tsv').write_text('observations')
sys.exit(int(str(fail)))
''')
        self.env = dict(os.environ, PGPHASE_TEST_CACHE=str(self.root / 'cache'))

    def command(self, directory, fail=0):
        argv = [str(self.binary), 'collect-graph-variation']
        for option, path in zip(cache_gap_replay.INPUT_OPTIONS, self.inputs):
            argv.extend((option, str(path)))
        extra = []
        if hasattr(self, 'matrix_prefix'):
            argv.extend(('--phase-matrix-dump', str(self.matrix_prefix)))
            extra = ['dump']
        return ('mkdir -p ' + shlex.quote(str(directory)) + ' && ' +
                shlex.join(argv) + ' && ' + shlex.join([
                    sys.executable, str(self.producer), str(directory),
                    str(self.counter), str(fail)] + extra))

    def launch(self, name, fail=0):
        directory = self.root / name
        return subprocess.Popen([
            sys.executable, str(Path(cache_gap_replay.__file__).resolve()),
            '--outdir', str(directory), '--command', self.command(directory, fail)],
            env=self.env, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

    def run_replay(self, name, fail=0):
        process = self.launch(name, fail)
        stdout, stderr = process.communicate()
        self.assertEqual(process.returncode, fail, (stdout, stderr))

    def calls(self):
        return len(self.counter.read_text().splitlines())

    def test_cross_process_reuse_and_invalidation(self):
        self.run_replay('first')
        self.run_replay('another-shard')
        self.assertEqual(self.calls(), 1)
        # Same size and restored mtime must still invalidate through ctime.
        path = self.inputs[1]
        old = path.stat()
        path.write_text('other')
        os.utime(path, ns=(old.st_atime_ns, old.st_mtime_ns))
        self.run_replay('input-changed')
        self.assertEqual(self.calls(), 2)
        index = Path(str(path) + '.bai')
        index.write_text('index')
        self.run_replay('index-added')
        self.assertEqual(self.calls(), 3)
        with self.binary.open('ab') as file:
            file.write(b'changed')
        self.run_replay('binary-changed')
        self.assertEqual(self.calls(), 4)
        # Corrupt cached artifacts are never accepted as a completed run.
        states = list((self.root / 'cache').glob('*/state.json'))
        for state in states:
            (state.parent / 'phased.bam').write_text('corrupt')
        self.run_replay('corrupt-output')
        self.assertEqual(self.calls(), 5)

    def test_failed_replays_never_publish_success(self):
        self.run_replay('failure-one', fail=1)
        self.run_replay('failure-two', fail=1)
        self.assertEqual(self.calls(), 2)
        self.assertEqual(list((self.root / 'cache').glob('*/state.json')), [])

    def test_concurrent_shards_run_once(self):
        processes = [self.launch('shard-' + str(i)) for i in range(3)]
        results = []
        for process in processes:
            stdout, stderr = process.communicate()
            results.append((process.returncode, stdout, stderr))
        for code, stdout, stderr in results:
            self.assertEqual(code, 0, (stdout, stderr))
        self.assertEqual(self.calls(), 1)
        for i in range(3):
            self.assertEqual((self.root / ('shard-' + str(i)) /
                              'phased.bam').read_text(), 'result')

    def test_matrix_state_restores_removed_debug_files(self):
        self.matrix_prefix = self.root / 'first-matrix'
        self.run_replay('matrix-first')
        Path(str(self.matrix_prefix) + '.chunk0.tsv').unlink()
        self.matrix_prefix = self.root / 'second-matrix'
        self.run_replay('matrix-second')
        self.assertEqual(self.calls(), 1)
        self.assertEqual(Path(str(self.matrix_prefix) + '.chunk0.tsv').read_text(),
                         'observations')


if __name__ == '__main__':
    unittest.main()
