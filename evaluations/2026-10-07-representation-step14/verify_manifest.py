#!/usr/bin/env python3
"""Independently verify every file recorded by the accepted evidence manifest."""
import hashlib
import json
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2]
manifest=json.loads(Path(__file__).with_name('manifest.json').read_text())
for record in manifest['files']:
    path=ROOT/record['path'];assert path.stat().st_size==record['size'],path
    h=hashlib.sha256()
    with path.open('rb') as f:
        while chunk:=f.read(1024*1024):h.update(chunk)
    assert h.hexdigest()==record['sha256'],path
fast=json.loads(Path(__file__).with_name('fast-checks.json').read_text())
assert hashlib.sha256((ROOT/'test_allele_context').read_bytes()).hexdigest()==fast['test_binary_sha256']
print(f"Verified {len(manifest['files'])} accepted file hashes and the timed test binary")
