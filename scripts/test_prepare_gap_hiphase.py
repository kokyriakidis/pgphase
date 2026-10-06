#!/usr/bin/env python3
"""Exercise competitor measurement reuse and identical-alignment enforcement."""
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
import pysam


class BenchmarkStateTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.root = Path(self.temporary.name)
        (self.root / 'truth.tsv').write_text('m\tMATERNAL\np\tPATERNAL\nu\tMATERNAL\n')
        (self.root / 'panel.tsv').write_text('gap_left\tgap_right\n10\t12\n')
        self.write_bam('input.bam', phased=False)
        self.write_bam('hiphase.bam', phased=True)
        pysam.index(str(self.root / 'input.bam'))

    def tearDown(self):
        self.temporary.cleanup()

    def write_bam(self, name, phased, paternal_hp=2, sequence='AAAA'):
        with pysam.AlignmentFile(self.root / name, 'wb', header={
                'HD': {'VN': '1.6'}, 'SQ': [{'SN': 'CHM13#0#chr20', 'LN': 100}]}) as bam:
            for qname in ('m', 'p', 'u'):
                read = pysam.AlignedSegment()
                read.query_name = qname
                read.query_sequence = sequence
                read.reference_id = 0
                read.reference_start = 9
                read.cigarstring = '4M'
                read.mapping_quality = 60
                if phased and qname != 'u':
                    read.set_tag('HP', paternal_hp if qname == 'p' else 1)
                    read.set_tag('PS', 10)
                bam.write(read)

    def run_measurement(self):
        return subprocess.run([sys.executable, 'scripts/prepare_gap_hiphase.py',
            '--bam', str(self.root / 'input.bam'), '--hiphase', str(self.root / 'hiphase.bam'),
            '--truth', str(self.root / 'truth.tsv'), '--panel', str(self.root / 'panel.tsv'),
            '--out', str(self.root / 'benchmark.tsv')], text=True, capture_output=True)

    def test_reuse_and_competitor_invalidation(self):
        result = self.run_measurement()
        self.assertEqual(result.returncode, 0, result.stderr)
        output = self.root / 'benchmark.tsv'
        self.assertIn('10-12\t3\t2\t2', output.read_text())
        timestamp = output.stat().st_mtime_ns
        result = self.run_measurement()
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn('Reused', result.stdout)
        self.assertEqual(output.stat().st_mtime_ns, timestamp)
        self.write_bam('hiphase.bam', phased=True, paternal_hp=1)
        result = self.run_measurement()
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn('10-12\t3\t1\t1', output.read_text())

    def test_alignment_mismatch_does_not_publish_measurements(self):
        self.write_bam('hiphase.bam', phased=True, sequence='CCCC')
        result = self.run_measurement()
        self.assertNotEqual(result.returncode, 0)
        self.assertIn('alignment differs', result.stderr)
        self.assertFalse((self.root / 'benchmark.tsv').exists())
        self.assertFalse((self.root / 'benchmark.tsv.json').exists())


if __name__ == '__main__':
    unittest.main()
