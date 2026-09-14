#!/usr/bin/env python3
"""Regression checks for matching audited VCF endpoints to native events."""

import tempfile
import unittest
from pathlib import Path

from trial_graph_gap_bam import endpoint_blocks


class EndpointTests(unittest.TestCase):
    def test_indel_anchor_coordinates(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "candidates.tsv"
            path.write_text("POS\tTYPE\tHAP_ALT\tHAP_REF\tPHASE_SET\n"
                            "101\tINS\t2\t1\t90\n"
                            "201\tDEL\t1\t2\t90\n"
                            "300\tSNP\t2\t1\t300\n"
                            "401\tINS\t3\t0\t300\n")
            self.assertEqual(endpoint_blocks(path, 100, 200), "joined")
            self.assertEqual(endpoint_blocks(path, 100, 300), "split")
            self.assertEqual(endpoint_blocks(path, 101, 300), "unresolved_endpoint")
            self.assertEqual(endpoint_blocks(path, 300, 400), "unresolved_endpoint")


if __name__ == "__main__":
    unittest.main()
