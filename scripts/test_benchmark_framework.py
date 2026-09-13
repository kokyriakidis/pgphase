#!/usr/bin/env python3
"""Focused tests for benchmark artifact locking and step invalidation."""

import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

import benchmark_panel
import merge_graph_hybrid_tags


class BenchmarkFrameworkTest(unittest.TestCase):
    def test_chromosome_scoped_competitor(self):
        chromosome = {
            "linear_reference": "/ref", "linear_bam": "/reads",
            "shared_vcf": "/calls", "truth_bam": "/truth",
        }
        manifest = {
            "chromosomes": {"chr12": chromosome, "chr20": chromosome},
            "competitors": {
                "panel": {"version": "1", "commands": ["panel"]},
                "native": {
                    "version": "1", "chromosomes": ["chr20"],
                    "commands": ["native ${chrom}"],
                    "inputs": ["${truth_bam}"],
                    "artifacts": ["${tool_dir}/native.vcf"],
                },
            },
            "competitor_artifacts": ["${tool_dir}/phased.vcf.gz"],
        }
        entries = benchmark_panel.resolved_competitors(manifest)
        self.assertEqual(
            [(entry["chromosome"], entry["tool"]) for entry in entries],
            [("chr12", "panel"), ("chr20", "panel"), ("chr20", "native")],
        )

    def test_graph_bridge_parity_and_distance_gate(self):
        accepted = {
            100: [(1000, 0, 20), (1100, 1, 18), (2000, 0, 17)],
        }
        merged, edges, conflicts, rejected, _ = \
            merge_graph_hybrid_tags.merge_graph_links(accepted, 300)
        self.assertEqual((edges, conflicts, rejected), (1, 0, 1))
        root_1, parity_1 = merged.find(1000)
        root_2, parity_2 = merged.find(1100)
        root_3, parity_3 = merged.find(2000)
        self.assertEqual(root_1, root_2)
        self.assertEqual(parity_1 ^ parity_2, 1)
        self.assertNotEqual(root_1, root_3)
        self.assertEqual(parity_3, 0)

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
