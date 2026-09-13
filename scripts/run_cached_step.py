#!/usr/bin/env python3
"""Run a command only when its argv, executable, inputs, or outputs changed."""

import argparse
import hashlib
import json
import os
import shlex
import shutil
import subprocess
import sys
import time
from datetime import datetime, timezone
from pathlib import Path


SAMPLE_BYTES = 1024 * 1024


def sampled_sha256(path):
    size = path.stat().st_size
    digest = hashlib.sha256()
    with path.open("rb") as fh:
        offsets = sorted({0, max(0, size // 2 - SAMPLE_BYTES // 2),
                          max(0, size - SAMPLE_BYTES)})
        for offset in offsets:
            fh.seek(offset)
            digest.update(offset.to_bytes(8, "little"))
            digest.update(fh.read(SAMPLE_BYTES))
    return digest.hexdigest()


def fingerprint(path_string):
    path = Path(path_string).resolve()
    stat = path.stat()
    return {
        "path": str(path),
        "size": stat.st_size,
        "sampled_sha256": sampled_sha256(path),
    }


def signature(command, inputs):
    executable = command[0]
    if os.path.sep not in executable:
        resolved = shutil.which(executable)
        if resolved is None:
            raise FileNotFoundError(f"executable not found: {executable}")
        executable = resolved
    payload = {
        "command": command,
        "executable": fingerprint(executable),
        "inputs": [fingerprint(path) for path in inputs],
    }
    encoded = json.dumps(payload, sort_keys=True, separators=(",", ":")).encode()
    return hashlib.sha256(encoded).hexdigest(), payload


def outputs_exist(outputs):
    return all(Path(path).is_file() and Path(path).stat().st_size > 0
               for path in outputs)


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--state", required=True)
    parser.add_argument("--input", action="append", default=[])
    parser.add_argument("--output", action="append", default=[])
    parser.add_argument("--force", action="store_true")
    parser.add_argument("command", nargs=argparse.REMAINDER)
    args = parser.parse_args()
    if args.command and args.command[0] == "--":
        args.command = args.command[1:]
    if not args.command:
        parser.error("a command is required after --")
    if not args.output:
        parser.error("at least one --output is required")
    return args


def main():
    args = parse_args()
    state_path = Path(args.state)
    current_signature, payload = signature(args.command, args.input)
    previous = None
    if state_path.exists():
        previous = json.loads(state_path.read_text())
    if not args.force and previous \
            and previous.get("signature") == current_signature \
            and outputs_exist(args.output):
        print(f"CACHE HIT {state_path}: {args.command[0]}")
        return 0

    state_path.parent.mkdir(parents=True, exist_ok=True)
    print("RUN " + shlex.join(args.command))
    started = time.monotonic()
    completed = subprocess.run(args.command, check=False)
    elapsed = time.monotonic() - started
    if completed.returncode != 0:
        return completed.returncode
    if not outputs_exist(args.output):
        missing = [path for path in args.output
                   if not Path(path).is_file() or Path(path).stat().st_size == 0]
        raise RuntimeError(f"command succeeded without expected outputs: {missing}")

    state = {
        "schema_version": 1,
        "signature": current_signature,
        "completed_at": datetime.now(timezone.utc).isoformat(),
        "elapsed_seconds": round(elapsed, 3),
        **payload,
        "outputs": [fingerprint(path) for path in args.output],
    }
    temporary = state_path.with_suffix(state_path.suffix + ".tmp")
    temporary.write_text(json.dumps(state, indent=2) + "\n")
    temporary.replace(state_path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
