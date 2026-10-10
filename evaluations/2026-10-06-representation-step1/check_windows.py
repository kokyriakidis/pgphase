#!/usr/bin/env python3
"""Run each registered gap check independently against one production binary."""
import argparse
import concurrent.futures
import hashlib
import json
import os
from pathlib import Path
import re
import subprocess
import time


ROOT = Path(__file__).resolve().parents[2]
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--binary", required=True)
parser.add_argument("--out", required=True)
args = parser.parse_args()
binary = (ROOT / args.binary).resolve(strict=True)
output = (ROOT / args.out).resolve()
output.mkdir(parents=True, exist_ok=True)
source = (ROOT / "src/test_gap_windows.cpp").read_text()
registry = source.split("static const GapCheck checks[] = {", 1)[1].split("\n    };", 1)[0]
names = re.findall(r'\{"([^"\n]+)", "\[gap\]', registry)
assert len(names) == len(set(names)) and len(names) >= 80
names.sort(key=lambda name: name != "chr20 gap windows")


def check(item):
    index, name = item
    work = output / str(index)
    work.mkdir(exist_ok=True)
    environment = dict(os.environ, PGPHASE_BIN=args.binary,
                       PGPHASE_GAP_FILTER=name, PGPHASE_TEST_WORKDIR=str(work))
    started = time.monotonic()
    with (work / "test.log").open("w") as log:
        result = subprocess.run(["./test_gap_windows", "all gaps", "--section", name,
                                 "--use-colour", "no"],
                                cwd=ROOT, env=environment, stdout=log, stderr=subprocess.STDOUT)
    text = (work / "test.log").read_text()
    assert "skipped: inputs absent" not in text
    outcome = dict(name=name, index=index, exit=result.returncode,
                   seconds=time.monotonic() - started, log=str(work / "test.log"))
    print(f"{'PASS' if result.returncode == 0 else 'FAIL'} {name}", flush=True)
    return outcome


started = time.monotonic()
results = []
with concurrent.futures.ThreadPoolExecutor(max_workers=3) as executor:
    for result in executor.map(check, enumerate(names)):
        results.append(result)
summary = dict(binary=str(binary), binary_sha256=hashlib.sha256(binary.read_bytes()).hexdigest(),
               registered_checks=len(names), passed=sum(r["exit"] == 0 for r in results),
               seconds=time.monotonic() - started, results=results)
(output / "checks.json").write_text(json.dumps(summary, indent=2) + "\n")
raise SystemExit(any(r["exit"] != 0 for r in results))
