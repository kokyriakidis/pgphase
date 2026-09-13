#!/usr/bin/env python3
"""Focused tests for benchmark artifact locking and step invalidation."""

import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

import benchmark_panel


class BenchmarkFrameworkTest(unittest.TestCase):
    def test_sampled_fingerprint_detects_content_change(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "artifact"
            path.write_bytes(b"before")
            before = benchmark_panel.fingerprint(path)
            path.write_bytes(b"after!")
            after = benchmark_panel.fingerprint(path)
            self.assertNotEqual(before["sampled_sha256"], after["sampled_sha256"])

    def test_cached_step_records_and_invalidates_input(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            source = root / "input.txt"
            output = root / "output.txt"
            state = root / "state.json"
            source.write_text("one\n")
            script = Path(__file__).with_name("run_cached_step.py")

            def command():
                return [
                    sys.executable, str(script), "--state", str(state),
                    "--input", str(source), "--output", str(output), "--",
                    sys.executable, "-c",
                    ("from pathlib import Path; "
                     f"Path({str(output)!r}).write_text(Path({str(source)!r}).read_text())"),
                ]

            first = subprocess.run(command(), check=True, text=True, capture_output=True)
            self.assertIn("RUN ", first.stdout)
            first_state = json.loads(state.read_text())["signature"]
            second = subprocess.run(command(), check=True, text=True, capture_output=True)
            self.assertIn("CACHE HIT", second.stdout)
            source.write_text("two\n")
            subprocess.run(command(), check=True, text=True, capture_output=True)
            second_state = json.loads(state.read_text())["signature"]
            self.assertNotEqual(first_state, second_state)
            self.assertEqual(output.read_text(), "two\n")


if __name__ == "__main__":
    unittest.main()
