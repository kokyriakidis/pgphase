#!/usr/bin/env python3
"""Regression checks for matching audited VCF endpoints to native events."""

import tempfile
import unittest
from pathlib import Path

from trial_graph_gap_bam import endpoint_blocks


class EndpointTests(unittest.TestCase):
    def test_indel_anchor_coordinates(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "native.vcf"
            path.write_text("##fileformat=VCFv4.2\n"
                            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n"
                            "chr20\t100\t.\tA\tAT\t60\tPASS\t.\tGT:PS\t0|1:90\n"
                            "chr20\t200\t.\tAT\tA\t60\tPASS\t.\tGT:PS\t1|0:90\n"
                            "chr20\t300\t.\tA\tT\t60\tPASS\t.\tGT:PS\t0|1:300\n"
                            "chr20\t400\t.\tA\tAC,ACC\t60\tPASS\t.\tPS:GT\t300:1|2\n"
                            "chr20\t500\t.\tA\tAC\t60\tPASS\t.\tGT:PS\t1|1:300\n"
                            "chr20\t600\t.\tA\tAC,ACC\t60\tPASS\t.\tGT:PS\t1/2:300\n")
            self.assertEqual(endpoint_blocks(path, 100, 200), "joined")
            self.assertEqual(endpoint_blocks(path, 100, 300), "split")
            self.assertEqual(endpoint_blocks(path, 101, 300), "unresolved_endpoint")
            self.assertEqual(endpoint_blocks(path, 300, 400), "joined")
            self.assertEqual(endpoint_blocks(path, 300, 500), "unresolved_endpoint")
            self.assertEqual(endpoint_blocks(path, 300, 600), "unresolved_endpoint")


if __name__ == "__main__":
    unittest.main()
