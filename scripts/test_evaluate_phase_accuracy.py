#!/usr/bin/env python3
"""Regression: parental assembly offsets must not change read diagnostics."""
import json
from pathlib import Path
import subprocess
import tempfile
import unittest


class ReadCoordinateTest(unittest.TestCase):
    def test_parental_offsets_do_not_change_transitions(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            header = '@HD\tVN:1.6\n@SQ\tSN:chr20\tLN:100000\n'
            phased = root / 'input.sam'
            # Five concordant reads followed by three discordant reads.
            statuses = [True] * 5 + [False] * 3
            rows = []
            for i, conc in enumerate(statuses):
                hp = 1 if i % 2 == 0 else 2
                rows.append(f'r{i}\t0\tchr20\t{100+i*100}\t60\t1M\t*\t0\t0\tA\tI\tHP:i:{hp}\tPS:i:100\n')
            phased.write_text(header + ''.join(rows))
            outputs = []
            for offset in (0, 10000):
                truth = root / f'truth{offset}.sam'
                rows = []
                for i, conc in enumerate(statuses):
                    mat = (i % 2 == 0) == conc
                    pos = 100 + i * 100 + (0 if mat else offset)
                    rows.append(f'r{i}\t0\tchr20\t{pos}\t60\t1M\t*\t0\t0\tA\tI\tHO:Z:{"MAT" if mat else "PAT"}\thq:i:60\n')
                truth.write_text(header + ''.join(rows))
                out = root / str(offset)
                out.mkdir()
                subprocess.run(['python3', str(Path(__file__).with_name('evaluate_phase_accuracy.py')),
                                str(phased), str(truth), '0', '0', '1', '', str(out), 'samtools'],
                               check=True, capture_output=True, text=True)
                outputs.append(json.loads((out / 'summary.json').read_text()))
            for result in outputs:
                self.assertEqual(result['discordant_reads'], 3)
                self.assertEqual(result['switch_errors'], 1)
                self.assertEqual(result['flip_errors'], 0)
                self.assertEqual(result['accuracy_metric_unit'], 'read')


if __name__ == '__main__':
    unittest.main()
