#!/usr/bin/env python3
"""Rescore verified native runs; execute any request absent from the cache.

Evaluation only. Cache keys retain every CLI argument except the output folder.
The saved requests came from the panel's recorded argument vectors; source
files are from the complete fresh final-binary run, not the accepted baseline.
"""
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import sys

manifest = json.loads(Path(__file__).with_name('native-cache.json').read_text())
binary = Path('pgphase').resolve()
assert hashlib.sha256(binary.read_bytes()).hexdigest() == manifest['binary_sha256']
for name, expected in manifest['input_stat'].items():
    stat = Path(name).stat()
    assert [stat.st_size, stat.st_mtime_ns] == expected, name
args = sys.argv[1:]
out = Path(args[args.index('-o') + 1]).parent
key = json.dumps([arg.replace(str(out), '<output>') for arg in args], separators=(',', ':'))
cached = manifest['runs'].get(key)
if cached is None:
    sys.exit(subprocess.run([str(binary), *args]).returncode)
source = Path(cached['source'])
out.mkdir(parents=True, exist_ok=True)
for name, expected in cached['hashes'].items():
    path = source / name
    assert hashlib.sha256(path.read_bytes()).hexdigest() == expected, path
    shutil.copyfile(path, out / name)
