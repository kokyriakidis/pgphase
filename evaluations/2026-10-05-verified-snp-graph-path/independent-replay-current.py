#!/usr/bin/env python3
"""Rescore final-binary outputs after an assertion-only test update.

The old and new test sources differ only in one CHECK and its comment.
CLI construction, case ordering and run cache insertion order are identical.
The manifest binds each ordinal to immutable final-binary output hashes.
"""
import hashlib
import json
from pathlib import Path
import shutil
import sys
root = Path('test_data/tmp_gap_fix52')
manifest = json.loads((root / 'replay-manifest.json').read_text())
assert hashlib.sha256(Path(manifest['binary']).read_bytes()).hexdigest() == manifest['binary_sha256']
for name, expected in manifest['input_stat'].items():
    stat = Path(name).stat()
    assert [stat.st_size, stat.st_mtime_ns] == expected, name
record = root / 'replay-requests.jsonl'
rows = record.read_text().splitlines() if record.exists() else []
ordinal = len(rows)
run = manifest['runs'][ordinal]
args = sys.argv[1:]
out = Path(args[args.index('-o') + 1]).parent
out.mkdir(parents=True, exist_ok=True)
assert [x.replace(manifest["replay_workdir"], "<workdir>") for x in args] == run["args"]
source = Path(run['source'])
for name, expected in run['hashes'].items():
    assert hashlib.sha256((source / name).read_bytes()).hexdigest() == expected, name
for flag, name in [('-o', 'candidates.tsv'), ('--phased-vcf-out', 'native.vcf'), ('--phased-bam-out', 'phased.bam')]:
    shutil.copyfile(source / name, args[args.index(flag) + 1])
if run['auxiliary']:
    prefix = args[args.index('--phase-matrix-dump') + 1]
    for auxiliary in run['auxiliary']:
        source_path = Path(auxiliary['source'])
        assert hashlib.sha256(source_path.read_bytes()).hexdigest() == auxiliary['sha256']
        shutil.copyfile(source_path, prefix + auxiliary['suffix'])
with record.open('a') as stream:
    stream.write(json.dumps({'ordinal': ordinal, 'args': args, 'source': str(source)}) + '\n')
sys.stdout.buffer.write((source / 'stdout.log').read_bytes())
sys.stderr.buffer.write((source / 'stderr.log').read_bytes())
